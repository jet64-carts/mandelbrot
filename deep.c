// carts/mandelbrot/deep.c
//
// Ownership (design.md §6.2): everything here belongs to the main thread but the slot being built
// and the pipeline's search state, which the one deep job of a batch advances while the batch is
// in flight; the main thread reads them only after the fence (deep_settle). The current reference
// is read by the batch's units and written by nobody while it is current.
#pragma STDC FP_CONTRACT OFF

#include "deep.h"

#include "config.h"
#include "jet_thread.h"
#include "numerics.h"

#define NO_SLOT (-1)
#define SLOTS 2
#define HALF 0.5

// Perturbation steps per slice of a search: about as long as a slice of bignum orbit.
#define DEEP_SLICE_STEPS 400000
#define DEEP_SLICE_STEPS_WEB 100000

// Slices one job runs when it rides alone — no batch beside it, as after a jump while the first
// reference is built: the pool is otherwise idle, and the tick has room for a longer job. (The
// web's job runs on the main thread at the fence, and is never lengthened.)
#define DEEP_ALONE_SLICES 4

// A nucleus search tries the deepest atom domains first, this many of them.
#define DEEP_NUCLEUS_CANDIDATES 3

// Newton in doubles stops within this fraction of a pixel (the bignum refinement finishes the job),
// or at the double's own resolution of the offset.
#define DEEP_NEWTON_PIXELS 0x1p-8
#define DEEP_NEWTON_RELATIVE 0x1p-48

typedef enum {
	PIPE_IDLE,
	PIPE_PERIOD,  // the target's orbit against the current reference, recording its atom domains
	PIPE_NEWTON,  // Newton for a candidate period against the current reference
	PIPE_BUILD,   // the next reference, in the bignum
} PipeStage;

static Reference s_slots[SLOTS];
static DeepOptions s_options;
static int s_current = NO_SLOT;
static int s_next = NO_SLOT;
static DeepStats s_stats;

// What the next reference is for, and how far its search has got.
static struct {
	PipeStage stage;
	DeepPoint target;     // the view's centre when the pipeline began
	double target_scale;  // its scale_log2
	int frac;             // the new reference's precision
	int32_t limit;        // the iteration limit it must cover
	double max_delta_c;   // its BLA table's reach
	double spacing;       // the view's pixel spacing when it began
	const Reference* base;  // the reference searched against (the current one)
	double dc_re, dc_im;  // the target's offset from it
	PeriodSearch period;
	NucleusSearch nucleus;
	int32_t candidates[REFERENCE_RECORDS];
	int candidate_count, candidate;
	bool from_centre;     // the next reference is a centre reference
	bool refining;        // the search is the current centre reference's own, for its nucleus
	double chosen_scale;  // for DEEP_REREFERENCE_OCTAVES: the scale the current one was chosen at
} s_pipe;

static double s_current_chosen_scale;  // the view's scale when the current reference was started
static int s_job_slices;               // slices the job in flight runs

uint32_t deep_bytes(void) {
	return SLOTS * reference_bytes();
}

bool deep_init(const DeepOptions* options) {
	s_options = *options;
	for (int i = 0; i < SLOTS; i++) {
		if (!reference_reserve(&s_slots[i])) return false;
	}
	s_current = s_next = NO_SLOT;
	s_pipe.stage = PIPE_IDLE;
	memset(&s_stats, 0, sizeof s_stats);
	return true;
}

// --- Serving ----------------------------------------------------------------------------------------
// The view's largest |delta_c| against a reference at c: the centre's distance and the screen's
// half-diagonal.
static double largest_delta_c(const View* view, const DeepPoint* c) {
	double dr, di;
	point_difference(&view->centre, c, &dr, &di);
	const double h = view_spacing(view);
	const double half_w = HALF * MB_SCREEN_WIDTH * h, half_h = HALF * MB_SCREEN_HEIGHT * h;
	return sqrt(dr * dr + di * di) + sqrt(half_w * half_w + half_h * half_h);
}

static bool within_neighbourhood(const View* view, const DeepPoint* c, double pixels) {
	double dr, di;
	point_difference(&view->centre, c, &dr, &di);
	const double reach = pixels * view_spacing(view);
	return dr * dr + di * di <= reach * reach;
}

// True when reference r can compute the field's picture: precision, neighbourhood, table, length.
static bool serves(const Reference* r, const FieldSpec* spec) {
	if (r->stage != REFERENCE_READY) return false;
	if (spec->view.scale_log2 < (double)reference_deepest_log2(r)) return false;
	if (!within_neighbourhood(&spec->view, &r->c, DEEP_NEIGHBOURHOOD_PIXELS)) return false;
	if (r->levels > 0 && largest_delta_c(&spec->view, &r->c) > r->max_delta_c) return false;
	return r->period > 0 || r->escaped || r->computed >= spec->max_iterations;
}

// True when the current reference will soon stop serving, or has been kept long enough that a
// fresher one should be found while it still can be.
static bool wants_successor(const Reference* r, const FieldSpec* spec) {
	const double scale = spec->view.scale_log2;
	if (scale < (double)reference_deepest_log2(r) + DEEP_LOOKAHEAD_OCTAVES) return true;
	if (!within_neighbourhood(&spec->view, &r->c, HALF * DEEP_NEIGHBOURHOOD_PIXELS)) return true;
	// Half the table's reach spent: the view is zooming out beyond what it was built for.
	if (r->levels > 0 && largest_delta_c(&spec->view, &r->c) > HALF * r->max_delta_c) return true;
	if (s_current_chosen_scale - scale > DEEP_REREFERENCE_OCTAVES) return true;
	return r->period == 0 && !r->escaped && r->computed < spec->max_iterations;
}

bool deep_ready(const Field* f, KernelReference* ref, double* origin_re, double* origin_im) {
	if (!f->spec.perturbed) return true;
	if (s_current == NO_SLOT || !serves(&s_slots[s_current], &f->spec)) return false;
	const Reference* r = &s_slots[s_current];
	*ref = reference_kernel(r);
	point_difference(&f->origin, &r->c, origin_re, origin_im);
	return true;
}

// --- The pipeline ------------------------------------------------------------------------------------
static int free_slot(void) {
	return s_current == 0 ? 1 : 0;
}

// Starts building a reference at c in the free slot.
static void begin_build(const DeepPoint* c, int32_t period) {
	s_next = free_slot();
	reference_begin(&s_slots[s_next], c, s_pipe.frac, s_pipe.limit, period, s_pipe.max_delta_c,
	                s_options.bla ? MB_BLA_EPSILON : 0.0);
	s_pipe.stage = PIPE_BUILD;
	s_pipe.from_centre = period == 0;
	s_stats.building = true;
}

// The next candidate period's Newton search, or, when none is left, the centre — unless the search
// was the current centre reference's own (promote), which already is that centre: then the pipeline
// rests, and the view's next need starts a fresh one. (Building the centre again there looped: the
// same point at the same precision, promoted, searched and rebuilt every few ticks, while the view
// dived past what that precision serves — found by the harness's dive to 10^254.)
static void next_candidate(void) {
	while (s_pipe.candidate < s_pipe.candidate_count && s_pipe.candidate < DEEP_NUCLEUS_CANDIDATES) {
		const int32_t period = s_pipe.candidates[s_pipe.candidate_count - 1 - s_pipe.candidate++];
		if (period < 1 || period > REFERENCE_CAPACITY) continue;
		double tolerance = DEEP_NEWTON_PIXELS * s_pipe.spacing;
		const double magnitude = sqrt(s_pipe.dc_re * s_pipe.dc_re + s_pipe.dc_im * s_pipe.dc_im);
		if (tolerance < DEEP_NEWTON_RELATIVE * magnitude) tolerance = DEEP_NEWTON_RELATIVE * magnitude;
		reference_nucleus_begin(&s_pipe.nucleus, period, s_pipe.dc_re, s_pipe.dc_im, tolerance,
		                        HALF * DEEP_NEIGHBOURHOOD_PIXELS * s_pipe.spacing);
		s_pipe.stage = PIPE_NEWTON;
		return;
	}
	if (s_pipe.refining) {
		s_pipe.stage = PIPE_IDLE;
		s_stats.building = false;
		return;
	}
	begin_build(&s_pipe.target, 0);
}

// Starts the candidate searches from a list of atom-domain periods, against `base`.
static void begin_candidates(const Reference* base, const int32_t* periods, int count, bool refining) {
	s_pipe.base = base;
	s_pipe.refining = refining;
	memcpy(s_pipe.candidates, periods, (size_t)count * sizeof(int32_t));
	s_pipe.candidate_count = count;
	s_pipe.candidate = 0;
	point_difference(&s_pipe.target, &base->c, &s_pipe.dc_re, &s_pipe.dc_im);
	next_candidate();
}

// Starts a pipeline for the field's view.
static void begin_pipeline(const FieldSpec* spec) {
	s_pipe.target = spec->view.centre;
	s_pipe.target_scale = spec->view.scale_log2;
	s_pipe.spacing = view_spacing(&spec->view);
	double deepest = spec->view.scale_log2 - DEEP_SPAN_OCTAVES;
	if (deepest < MB_DEPTH_LIMIT_LOG2) deepest = MB_DEPTH_LIMIT_LOG2;
	s_pipe.frac = reference_limbs_for(deepest);
	// The limit it must cover: the view's own, and what the limit will have grown to DEEP_SPAN_OCTAVES
	// deeper (a centre reference's orbit must reach it; a nucleus's covers any).
	View deeper = spec->view;
	deeper.scale_log2 = deepest;
	int64_t limit = view_iteration_limit(&deeper, 0);
	if (limit < spec->max_iterations) limit = spec->max_iterations;
	s_pipe.limit = limit > REFERENCE_CAPACITY ? REFERENCE_CAPACITY : (int32_t)limit;
	s_pipe.max_delta_c = DEEP_TABLE_REACH * largest_delta_c(&spec->view, &spec->view.centre);
	s_pipe.chosen_scale = spec->view.scale_log2;
	const Reference* base = s_current == NO_SLOT ? NULL : &s_slots[s_current];
	if (s_options.nucleus && base && serves(base, spec)) {
		// Search for the view's minibrot against the reference in use.
		double dr, di;
		point_difference(&spec->view.centre, &base->c, &dr, &di);
		reference_period_begin(&s_pipe.period, dr, di, s_pipe.limit);
		s_pipe.base = base;
		s_pipe.stage = PIPE_PERIOD;
		s_stats.building = true;
		return;
	}
	begin_build(&spec->view.centre, 0);
}

static void promote(void) {
	Reference* r = &s_slots[s_next];
	s_current = s_next;
	s_next = NO_SLOT;
	s_current_chosen_scale = s_pipe.chosen_scale;
	s_pipe.stage = PIPE_IDLE;
	s_stats.building = false;
	s_stats.references_built++;
	s_stats.nuclei_found += r->period > 0;
	// A centre reference knows its own atom domains: look for its minibrot's nucleus straight away.
	if (s_options.nucleus && s_pipe.from_centre && r->record_count > 0) {
		s_pipe.target = r->c;
		begin_candidates(r, r->record_at, r->record_count, true);
		s_stats.building = true;
	}
}

void deep_settle(const Field* f) {
	const FieldSpec* spec = &f->spec;
	if (s_pipe.stage == PIPE_BUILD && s_slots[s_next].stage == REFERENCE_READY) promote();
	const Reference* current = s_current == NO_SLOT ? NULL : &s_slots[s_current];
	if (current) {
		s_stats.serving = spec->perturbed && serves(current, spec);
		s_stats.period = current->period;
		s_stats.length = current->computed;
		s_stats.frac = current->frac;
	}
	if (s_pipe.stage != PIPE_IDLE) return;
	const double ahead = s_options.lookahead ? DEEP_CROSSOVER_LOOKAHEAD_OCTAVES : 0.0;
	const bool approaching = spec->view.scale_log2 < view_direct_limit(&spec->view) + ahead;
	if (!spec->perturbed && !approaching) return;
	const bool successor = s_options.lookahead && current && wants_successor(current, spec);
	if (!current || successor || !serves(current, spec)) begin_pipeline(spec);
}

// The deep job: one slice of whichever stage the pipeline is at.
static void deep_job(void* unused) {
	(void)unused;
	const bool web = jet_host_is_web();
	const int64_t steps = (web ? DEEP_SLICE_STEPS_WEB : DEEP_SLICE_STEPS) * s_job_slices;
	switch (s_pipe.stage) {
	case PIPE_PERIOD:
		if (reference_period_advance(&s_pipe.period, s_pipe.base, steps) == SEARCH_FOUND) {
			if (s_pipe.period.record_count > 0) {
				begin_candidates(s_pipe.base, s_pipe.period.record_at, s_pipe.period.record_count, false);
			} else {
				begin_build(&s_pipe.target, 0);
			}
		}
		break;
	case PIPE_NEWTON: {
			const SearchStatus status = reference_nucleus_advance(&s_pipe.nucleus, s_pipe.base, steps);
			if (status == SEARCH_FOUND) {
				DeepPoint nucleus = s_pipe.base->c;
				point_offset(&nucleus, s_pipe.nucleus.cr, s_pipe.nucleus.ci);
				// The table must reach the view from the nucleus, wherever in the neighbourhood it lies.
				s_pipe.max_delta_c += sqrt(s_pipe.nucleus.cr * s_pipe.nucleus.cr + s_pipe.nucleus.ci * s_pipe.nucleus.ci);
				begin_build(&nucleus, s_pipe.nucleus.period);
			} else if (status == SEARCH_FAILED) {
				next_candidate();
			}
			break;
		}
	case PIPE_BUILD:
		reference_advance(&s_slots[s_next],
		                  (web ? REFERENCE_SLICE_PRODUCTS_WEB : REFERENCE_SLICE_PRODUCTS) * s_job_slices,
		                  (web ? REFERENCE_SLICE_ENTRIES_WEB : REFERENCE_SLICE_ENTRIES) * s_job_slices);
		break;
	case PIPE_IDLE:
		break;
	}
}

bool deep_launch(bool alone) {
	if (s_pipe.stage == PIPE_IDLE) return false;
	s_job_slices = alone && !jet_host_is_web() ? DEEP_ALONE_SLICES : 1;
	jet_thread_pool_add(deep_job, NULL);
	return true;
}

const DeepStats* deep_stats(void) {
	return &s_stats;
}
