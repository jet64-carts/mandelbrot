// carts/mandelbrot/scheduler.c
//
// Ownership (design.md §6.2). Everything static in this file belongs to the main thread. While a
// batch is in flight (from submit to fence) the workers read the frozen contexts and their own
// units' inputs, and write their own units' results and their own tiles' pixels and corners; the
// main thread touches none of those until jet_thread_pool_wait returns, which is the batch's
// happens-before edge in both directions. A reconcile pass is a batch like any other, fenced
// before the settle goes on.
#pragma STDC FP_CONTRACT OFF

#include "scheduler.h"

#include "certify.h"
#include "deep.h"
#include "jet_thread.h"
#include "reuse.h"
#include "tiles.h"
#include "work.h"

// No reconcile this settle.
#define NO_RECONCILE (-1)

// At depth the grid's origin follows the view: once the centre strays this many pixels from it,
// the origin moves (field_move_origin), exactly.
#define ORIGIN_DRIFT_PIXELS 1048576.0

static Field* s_field;
static SchedulerOptions s_options;
static int64_t s_quantum, s_restart_quantum, s_motion_quantum, s_rim_quantum;

static WorkContext s_context;
static ReconcileContext s_reconcile;
static Corner* s_corners;
static uint32_t s_run;       // counts compute runs; a corner record from an earlier run is stale
static bool s_certify_run;   // this run computes corners and certifies cells
static bool s_coarse_first;  // this run plans where its previews are coarsest first (plan)

static WorkUnit s_units[TILE_COUNT];
static ReconcileUnit s_reconcile_units[TILE_COUNT];
static int s_unit_count;
static bool s_in_flight;
static bool s_deep_in_flight;  // the batch in flight carries a slice of the next reference

// Per tile: its place in the schedule, the pixels each step still needs, its pixels by coarse level
// (work.h), and what one of its samples has cost, in lane-iterations.
static int s_position[TILE_COUNT];
static int32_t s_need[TILE_COUNT][TILE_STEPS];
static int32_t s_coarse[TILE_COUNT][WORK_COARSE_LEVELS];
static double s_sample_cost[TILE_COUNT];
static bool s_costs_perturbed;  // the kind of picture s_sample_cost was measured on

static AxisMap s_columns, s_rows;
static SchedulerStats s_stats;

SchedulerOptions scheduler_default_options(void) {
	SchedulerOptions o;
	o.reuse = true;
	o.pipeline = true;
	o.exterior_certificate = true;
	o.interior_certificate = true;
	o.longest_first = true;
	o.deep.bla = true;
	o.deep.nucleus = true;
	o.deep.lookahead = true;
	return o;
}

#define CORNER_BYTES ((uint32_t)(CERT_CORNERS * sizeof(Corner)))

uint32_t scheduler_bytes(void) {
	return CORNER_BYTES + deep_bytes();
}

bool scheduler_init(Field* field, const SchedulerOptions* options) {
	s_field = field;
	s_options = *options;
	const bool web = jet_host_is_web();
	s_quantum = web ? MB_QUANTUM_WEB : MB_QUANTUM_NATIVE;
	s_restart_quantum = web ? MB_RESTART_QUANTUM_WEB : MB_RESTART_QUANTUM_NATIVE;
	s_motion_quantum = web ? MB_MOTION_QUANTUM_WEB : MB_MOTION_QUANTUM_NATIVE;
	s_rim_quantum = web ? MB_RIM_QUANTUM_WEB : MB_RIM_QUANTUM_NATIVE;
	s_corners = (Corner*)jet_memory_reserve(CORNER_BYTES);
	if (!s_corners) return false;
	memset(s_corners, 0, CORNER_BYTES);  // run 0 is never current: every record starts stale
	s_run = 0;
	s_certify_run = false;
	s_coarse_first = false;
	s_unit_count = 0;
	s_in_flight = false;
	tiles_init();
	work_init();
	for (int t = 0; t < TILE_COUNT; t++) {
		s_position[t] = TILE_STEPS;
		memset(s_coarse[t], 0, sizeof s_coarse[t]);
		s_sample_cost[t] = MB_DEFAULT_SAMPLE_COST;
	}
	s_costs_perturbed = false;
	memset(&s_stats, 0, sizeof s_stats);
	s_stats.reconcile = NO_RECONCILE;
	s_stats.columns = &s_columns;
	s_stats.rows = &s_rows;
	if (!deep_init(&s_options.deep)) return false;
	jet_thread_pool_init(MB_WORKER_COUNT);
	return true;
}

// --- Batches -------------------------------------------------------------------------------------
static void submit_units(void) {
	for (int i = 0; i < s_unit_count; i++) jet_thread_pool_add(work_run, &s_units[i]);
	s_in_flight = s_unit_count > 0;
}

// Waits for the compute batch in flight, then folds what each unit did into its tile's progress,
// needs and cost. A unit's cost per sample is exactly counted, so the next plan is as
// deterministic as this one. Returns the lane-iterations the batch spent.
static int64_t fence(void) {
	if (!s_in_flight && !s_deep_in_flight) return 0;
	jet_thread_pool_wait();
	s_deep_in_flight = false;
	int64_t spent = 0;
	for (int i = 0; i < s_unit_count; i++) {
		const WorkUnit* u = &s_units[i];
		s_position[u->tile] = u->last_step;
		memcpy(s_need[u->tile], u->need, sizeof s_need[u->tile]);
		memcpy(s_coarse[u->tile], u->coarse, sizeof s_coarse[u->tile]);
		if (u->samples > 0) s_sample_cost[u->tile] = (double)u->work / (double)u->samples;
		s_stats.certified += u->certified;
		s_stats.samples += u->samples;
		spent += u->work;
	}
	s_unit_count = 0;
	s_in_flight = false;
	return spent;
}

// --- Planning ------------------------------------------------------------------------------------
// The step that enters `level` in a run that certifies: its unit certifies the cells of the level
// before, so it is never skipped as free, and waits for its neighbours (neighbours_ready).
static bool certifies(int step) {
	if (!s_certify_run) return false;
	const int level = tiles_step(step)->level;
	return level >= 1 && level <= CERT_LEVELS && step == tiles_level_start(level);
}

// The first step at or after `step` that has pixels to compute (or cells to certify).
static int skip_free(int tile, int step) {
	while (step < TILE_STEPS && s_need[tile][step] == 0 && !certifies(step)) step++;
	return step;
}

// True when the tiles whose corners tile `t`'s cells share — right, below, below-right — had all
// finished the level below `level` when this batch was planned: what the tile then reads across its
// borders is final, since nothing in the batch still writes it. So a tile certifies every cell,
// never a share of them that happened to be ready, and the converged picture does not depend on
// how the work was cut into batches.
static bool neighbours_ready(int t, int level, const int* base) {
	const int done = tiles_level_start(level);
	const bool right = t % TILE_COLUMNS + 1 < TILE_COLUMNS, below = t / TILE_COLUMNS + 1 < TILE_ROWS;
	if (right && base[t + 1] < done) return false;
	if (below && base[t + TILE_COLUMNS] < done) return false;
	return !(right && below && base[t + TILE_COLUMNS + 1] < done);
}

static WorkUnit* add_unit(int tile, int first_step) {
	WorkUnit* u = &s_units[s_unit_count++];
	u->context = &s_context;
	u->tile = tile;
	u->first_step = first_step;
	u->last_step = first_step;
	u->predicted = 0;
	u->work = 0;
	u->samples = 0;
	u->certified = 0;
	memcpy(u->need, s_need[tile], sizeof u->need);
	memcpy(u->coarse, s_coarse[tile], sizeof u->coarse);
	return u;
}

// A plan in the making: each tile's unit (or -1), where its schedule stood, what is left to spend,
// and the rules it is dealt under.
typedef struct {
	int unit_of[TILE_COUNT];
	int base[TILE_COUNT];
	int64_t budget;
	int forced_levels;
	bool progress;
} Planning;

// Deals no matter how fine the tile's previews already are.
#define DEAL_EVERY_TILE (-1)

// Deals the steps of `level` round-robin, one step per tile per round, nearest the centre first —
// to the tiles holding previews coarser than level `coarser_than`, or with DEAL_EVERY_TILE to every
// tile — until no such tile's next step is of that level or the budget is spent, which ends the plan
// (false). A step below `forced_levels` is dealt whatever it costs, and a batch that must `progress`
// takes its first step however costly.
static bool deal(Planning* p, int level, int coarser_than) {
	const bool forced = level < p->forced_levels;
	const int* order = tiles_order();
	for (bool dealt = true; dealt;) {
		dealt = false;
		for (int i = 0; i < TILE_COUNT; i++) {
			const int t = order[i];
			if (coarser_than != DEAL_EVERY_TILE && work_coarser_than(s_coarse[t], coarser_than) == 0) continue;
			const int step = p->unit_of[t] < 0 ? p->base[t] : s_units[p->unit_of[t]].last_step;
			if (step >= TILE_STEPS || tiles_step(step)->level != level) continue;
			if (certifies(step) && !neighbours_ready(t, level, p->base)) continue;
			const int64_t cost = (int64_t)(s_need[t][step] * s_sample_cost[t]);
			if (!forced && cost > p->budget && !(p->progress && s_unit_count == 0)) return false;
			if (p->unit_of[t] < 0) {
				p->unit_of[t] = s_unit_count;
				add_unit(t, step);
			}
			WorkUnit* u = &s_units[p->unit_of[t]];
			u->last_step = skip_free(t, step + 1);
			u->predicted += cost;
			p->budget -= cost;
			dealt = true;
		}
	}
	return true;
}

// Plans the next batch, priced by what each step still needs, until `quantum` is spent: first, in a
// run that plans so, where the picture is coarsest, then breadth-first by level. Coarsest first
// means level by level: the tiles holding previews coarser than level L — values from farther away
// than its blocks would leave them — are dealt every level up to L before any finer pass begins. So
// while the budget can buy a sample near every pixel, no pixel shows a value from across the screen
// (a zoom that outruns the pool sees its new view coarse rather than its old one smeared), and a
// sustained zoom at depth, whose moving lattice keeps every level 0 to do, still refines the slow
// middle of the screen rather than leaving it to its copies. A picture that holds no previews (a
// snap of a converged one, or a shallow move that reuse carried whole) is planned breadth-first
// either way. With `rim_only`, the plan is that first pass's first level alone: the coarsest level of
// the tiles holding far previews. Every step below level `forced_levels` is planned whatever it costs,
// and a batch that must `progress` takes its first step however costly. Every tile gets at most one
// unit, a contiguous run of its schedule, so no two units of a batch share a pixel.
static void plan(int64_t quantum, int forced_levels, bool progress, bool rim_only) {
	Planning p;
	for (int t = 0; t < TILE_COUNT; t++) {
		s_position[t] = skip_free(t, s_position[t]);  // steps with nothing to do are simply done
		p.base[t] = s_position[t];
		p.unit_of[t] = -1;
	}
	p.budget = quantum;
	p.forced_levels = forced_levels;
	p.progress = progress;
	if (rim_only) {
		deal(&p, 0, 0);
		return;
	}
	for (int coarser_than = 0; s_coarse_first && coarser_than < TILE_LEVELS; coarser_than++) {
		for (int level = 0; level <= coarser_than; level++) {
			if (!deal(&p, level, coarser_than)) return;
		}
	}
	for (int level = 0; level < TILE_LEVELS; level++) {
		if (!deal(&p, level, DEAL_EVERY_TILE)) return;
	}
}

// Orders the planned units by predicted cost, longest first: the classic LPT rule, which keeps
// the costliest unit from being the last one started. An insertion sort, since there are at most
// TILE_COUNT units and they arrive nearly sorted by level.
static void sort_longest_first(void) {
	if (!s_options.longest_first) return;
	for (int i = 1; i < s_unit_count; i++) {
		const WorkUnit u = s_units[i];
		int at = i;
		while (at > 0 && s_units[at - 1].predicted < u.predicted) {
			s_units[at] = s_units[at - 1];
			at--;
		}
		s_units[at] = u;
	}
}

// Whether the field's picture can be computed now, freezing the reference into the batch context if
// it is deep: a batch uses the reference that was current when it was planned, whatever happens
// to the pipeline while it runs.
static bool computable(void) {
	KernelReference reference;
	double origin_re = 0.0, origin_im = 0.0;
	if (!deep_ready(s_field, &reference, &origin_re, &origin_im)) return false;
	if (s_field->spec.perturbed) work_context_reference(&s_context, &reference, origin_re, origin_im);
	return true;
}

// Plans and runs batches synchronously until `quantum` is spent or nothing is left to do. With
// forced levels, the first batch is those levels alone, whatever they cost: it is also how the
// schedule learns what THIS picture's samples cost — a restart's tiles still carry the last
// picture's, which can be wrong by orders of magnitude. Only the first batch must progress: a later
// one whose cheapest step would overrun what is left of the quantum leaves it to the next draw, so
// beyond the first batch's first step the quantum is overrun only by what the plan misprices. With
// `rim_only`, only the rim is planned (plan), and nothing must progress: a deep step that would
// overrun the quantum is the batch's. Returns the lane-iterations spent.
static int64_t compute_now(int64_t quantum, int forced_levels, bool rim_only) {
	int64_t spent = 0;
	if (!computable()) return 0;
	while (spent < quantum) {
		plan(forced_levels > 0 ? 0 : quantum - spent, forced_levels, spent == 0 && !rim_only, rim_only);
		if (s_unit_count == 0) break;
		sort_longest_first();
		submit_units();
		spent += fence();
		forced_levels = 0;
	}
	return spent;
}

// --- Reconciling -----------------------------------------------------------------------------------
// How the settle reaches the new spec: the reconcile's kind, its two flags, and how a preview's
// distance scales from the old grid to the new.
typedef struct {
	int kind;
	bool demote;            // across the crossover
	bool keep;              // a restart at depth, previewing the old picture
	double distance_scale;  // the old pixel spacing over the new
} Plan;

// Brings the field to its newly planned spec: one job per tile rebuilds the pixels (reuse.h) and
// counts what each step needs; then the new grid (and, for a gather, the new pixels) is swapped
// in, every schedule starts again, and a new compute run begins.
static void reconcile(const Plan* plan, int64_t old_limit, bool at_rest) {
	const ReconcileKind kind = (ReconcileKind)plan->kind;
	const bool keep = plan->keep;
	s_reconcile.field = s_field;
	s_reconcile.kind = kind;
	s_reconcile.demote = plan->demote;
	s_reconcile.keep = keep;
	s_reconcile.distance_scale = plan->distance_scale;
	s_reconcile.old_limit = old_limit;
	s_reconcile.new_limit = s_field->spec.max_iterations;
	s_reconcile.columns = &s_columns;
	s_reconcile.rows = &s_rows;
	for (int t = 0; t < TILE_COUNT; t++) {
		s_reconcile_units[t].context = &s_reconcile;
		s_reconcile_units[t].tile = t;
		jet_thread_pool_add(reconcile_run, &s_reconcile_units[t]);
	}
	jet_thread_pool_wait();
	field_swap(s_field, kind == RECONCILE_GATHER || keep);
	for (int t = 0; t < TILE_COUNT; t++) {
		s_position[t] = 0;
		memcpy(s_need[t], s_reconcile_units[t].need, sizeof s_need[t]);
		memcpy(s_coarse[t], s_reconcile_units[t].coarse, sizeof s_coarse[t]);
	}
	// Across the crossover a sample's iterations change price (kernel.h): a first deep batch
	// planned on the shallow picture's costs ran five times its quantum (79 M lane-iterations
	// against 16 M, on a dive's first perturbed tick).
	if (s_field->spec.perturbed != s_costs_perturbed) {
		const double price = (double)KERNEL_PERTURBED_COST_PERCENT / (double)KERNEL_PERCENT;
		for (int t = 0; t < TILE_COUNT; t++) s_sample_cost[t] *= s_field->spec.perturbed ? price : 1.0 / price;
		s_costs_perturbed = s_field->spec.perturbed;
	}

	// A fresh lattice at rest — a snap's, or a restart's that will be kept — computes corners and
	// certifies cells. A moving picture's strips are single pixels, with no cells to certify, and a
	// restart while the view still moves (a zoom too fast to reuse) is shown for a tick and then
	// carried or dropped: its corners, at three times the cost of plain samples, would buy nothing.
	// Its picture is certified when the view comes to rest and snaps.
	s_run++;
	s_certify_run = (kind == RECONCILE_SNAP || (kind == RECONCILE_RESTART && at_rest)) &&
	                (s_options.exterior_certificate || s_options.interior_certificate);
	// A restart's schedule is coarse to fine by construction — every pixel starts empty, or held, or
	// warped from a picture too far off to carry — so it is planned breadth-first, as it always was:
	// ranked by its previews' distances, its plan would only be reordered around the cells its
	// certificates fill (which are not previews), and converge a draw later. Every other run plans
	// where the pictures that motion leaves are coarsest first (plan).
	s_coarse_first = kind != RECONCILE_RESTART;
	CertifyContext cells;
	cells.spacing = view_spacing(&s_field->spec.view);
	cells.max_iterations = s_field->spec.max_iterations;
	cells.exterior_margin = MB_EXTERIOR_MARGIN;
	cells.interior_margin = MB_INTERIOR_MARGIN;
	cells.exterior = s_options.exterior_certificate;
	cells.interior = s_options.interior_certificate;
	cells.run = s_run;
	cells.corners = s_corners;
	work_context_init(&s_context, s_field, s_certify_run, &cells);
	s_stats.kept = keep;
	s_stats.reconcile = kind;
	s_stats.certified = 0;
	s_stats.samples = 0;
}

// Moves the grid's origin to where the spec needs it: the plane's zero in direct mode, the view's
// centre at depth. Rounds the grid's coordinates (the caller demotes or previews the pixels).
static void place_origin(const FieldSpec* spec) {
	const DeepPoint zero = point_make(0.0, 0.0);
	field_set_origin(s_field, spec->perturbed ? &spec->view.centre : &zero);
}

// At depth, follows the view with the origin once it strays ORIGIN_DRIFT_PIXELS on either axis, by
// a double offset: every coordinate then lies within a factor of two of the offset, and the
// subtraction is exact (Sterbenz).
static void follow_with_origin(const FieldSpec* spec) {
	double d_re, d_im;
	point_difference(&spec->view.centre, &s_field->origin, &d_re, &d_im);
	const double drift = ORIGIN_DRIFT_PIXELS * view_spacing(&spec->view);
	const bool re = d_re > drift || d_re < -drift, im = d_im > drift || d_im < -drift;
	if (re || im) field_move_origin(s_field, re ? d_re : 0.0, im ? d_im : 0.0);
}

// Decides how the field reaches `spec`, and plans the new grid for it: nothing, when the field
// already shows it; a limit change, when only the limit moved; a gather, when the view moved and
// enough of the old picture still falls on the new grid; otherwise a restart — which at depth
// keeps a preview of the old picture, since nothing can be computed until a reference is ready.
static Plan classify(const FieldSpec* spec, bool at_rest) {
	Field* f = s_field;
	Plan plan = {NO_RECONCILE, false, false, 1.0};
	const bool jump = !f->has_spec || f->spec.epoch != spec->epoch;
	if (!jump && view_equal(&f->spec.view, &spec->view)) {
		if (f->spec.max_iterations != spec->max_iterations) {
			field_plan_spec(f, spec);
			reuse_keep_grid(f);
			plan.kind = RECONCILE_LIMIT;
		} else if (at_rest && !field_grid_is_ideal(f)) {
			field_plan_spec(f, spec);
			reuse_snap_grid(f);
			plan.kind = RECONCILE_SNAP;
		}
		return plan;
	}
	const bool had_picture = f->has_spec;
	const bool crossing = had_picture && f->spec.perturbed != spec->perturbed;
	if (had_picture) plan.distance_scale = view_spacing(&f->spec.view) / view_spacing(&spec->view);
	if (jump || crossing) {
		place_origin(spec);
	} else if (spec->perturbed) {
		follow_with_origin(spec);
	}
	field_plan_spec(f, spec);
	if (!jump && s_options.reuse) {
		s_stats.matched_share = reuse_match(f, &s_columns, &s_rows);
		if (s_stats.matched_share >= MB_MIN_REUSE_SHARE) {
			plan.kind = RECONCILE_GATHER;
			plan.demote = crossing;
			return plan;
		}
	}
	plan.kind = RECONCILE_RESTART;
	plan.keep = spec->perturbed && had_picture;
	if (plan.keep && jump) {
		reuse_hold_maps(&s_columns, &s_rows);
	} else if (plan.keep) {
		reuse_warp_maps(f, &s_columns, &s_rows);
	}
	reuse_ideal_grid(f);
	return plan;
}

void scheduler_settle(const FieldSpec* spec, bool at_rest) {
	s_stats.batch_work = fence();
	s_stats.settle_work = 0;
	s_stats.reconcile = NO_RECONCILE;
	const int64_t old_limit = s_field->has_spec ? s_field->spec.max_iterations : spec->max_iterations;
	const Plan plan = classify(spec, at_rest);
	if (plan.kind != NO_RECONCILE) reconcile(&plan, old_limit, at_rest);
	// The deep view: promote a finished reference, start the next one if the view wants it.
	if (s_field->has_spec) deep_settle(s_field);
	if (plan.kind == NO_RECONCILE) return;
	// A restart computes its coarsest level whatever it costs, so the frame is complete — and, at
	// rest (a jump the player will look at), finer levels within the restart quantum; while the
	// view still moves the picture is shown for a tick, and gets the motion quantum. A gather or a
	// limit change computes the pixels it could not carry. A snap changes no value, and its
	// refinement is the pipelined batches' work. At depth the draw computes nothing but the RIM:
	// every pixel already shows a preview (a kept restart's old picture, a gather's neighbour), and a
	// deep pixel costs thousands of shallow ones, so the quantum would buy a sliver of the picture
	// while holding the draw for milliseconds; the batch launched beside the frame computes it
	// instead. But a preview a zoom stretched beyond the old picture is FAR — a copy of its edge
	// drawn across the new screen — and a sample in 256 of it, the coarsest level, is cheap: so those
	// tiles' coarsest level is sampled now, within the rim quantum, and the frame shows the new view
	// coarse where the old one has run out.
	if (s_field->spec.perturbed) {
		s_stats.settle_work = compute_now(s_rim_quantum, 0, true);
		return;
	}
	if (plan.kind == RECONCILE_RESTART) {
		s_stats.settle_work = compute_now(at_rest ? s_restart_quantum : s_motion_quantum, 1, false);
	} else if (plan.kind != RECONCILE_SNAP) {
		s_stats.settle_work = compute_now(s_motion_quantum, 0, false);
	}
}

// At rest the whole batch quantum refines the picture the player is looking at. On the move a
// batch's pixels may be carried off the screen by the next gather, and — where ticks run back to
// back — it shares its tick with the next draw's motion compute, so it gets the motion quantum.
// (At depth, where the draw computes nothing, a batch of both motion quanta was measured through the
// door: its p95 tick 7.0 ms and 8-18 late at 10^159, against 4.2-4.4 ms and none; design.md §9.9.)
void scheduler_launch(bool at_rest) {
	if (s_in_flight || s_deep_in_flight) return;
	if ((at_rest || s_options.pipeline) && computable()) {
		plan(at_rest ? s_quantum : s_motion_quantum, 0, true, false);
		sort_longest_first();
		submit_units();
	}
	// The next reference's slice rides with the batch, or alone when there is nothing to compute.
	s_deep_in_flight = deep_launch(!s_in_flight);
}

void scheduler_end(void) {
	fence();
	jet_thread_pool_shutdown();
}

double scheduler_progress(void) {
	int64_t needed = 0;
	for (int t = 0; t < TILE_COUNT; t++) {
		for (int step = 0; step < TILE_STEPS; step++) needed += s_need[t][step];
	}
	return 1.0 - (double)needed / (double)MB_SCREEN_PIXELS;
}

bool scheduler_converged(void) {
	for (int t = 0; t < TILE_COUNT; t++) {
		if (skip_free(t, s_position[t]) < TILE_STEPS) return false;
	}
	return field_grid_is_ideal(s_field);
}

const SchedulerStats* scheduler_stats(void) {
	return &s_stats;
}
