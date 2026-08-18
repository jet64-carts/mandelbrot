// carts/mandelbrot/work.c
#pragma STDC FP_CONTRACT OFF

#include "work.h"

#include "numerics.h"
#include "tiles.h"

// Samples iterated per kernel call: enough to keep every slot of the kernel busy across several
// sample rows, small enough that the chunk lives comfortably on a worker's stack (about 20 KiB).
#define WORK_CHUNK 256

// How far each level's blocks reach: a block's far corner lies (stride - 1) sqrt 2 from its sample,
// as a preview's distance (field.h).
#define WORK_SQRT2 1.41421356237309504880
#define WORK_LEVEL_REACH(level) \
	((uint16_t)((double)((TILE_COARSEST_STRIDE >> (level)) - 1) * DISTANCE_PIXEL * WORK_SQRT2))
static const uint16_t WORK_REACH[TILE_LEVELS] = {
	WORK_LEVEL_REACH(0), WORK_LEVEL_REACH(1), WORK_LEVEL_REACH(2), WORK_LEVEL_REACH(3), WORK_LEVEL_REACH(4),
};
_Static_assert(TILE_LEVELS == 5, "one reach for each level of the schedule");

// Each pixel's distance from the sample whose block holds it, by its offset from the sample, and the
// coarse level of a preview at that distance: a block is the square below and to the right of its
// sample, the coarsest TILE_COARSEST_STRIDE wide, and every finer one its corner. A sample fills a
// block of 4 to 256 pixels, so the tables spare a square root and a search at each.
static uint16_t s_block_distance[TILE_COARSEST_STRIDE][TILE_COARSEST_STRIDE];
static uint8_t s_block_level[TILE_COARSEST_STRIDE][TILE_COARSEST_STRIDE];

void work_init(void) {
	const uint8_t preview = state_make(RANK_PREVIEW, FATE_ESCAPED);
	for (int dy = 0; dy < TILE_COARSEST_STRIDE; dy++) {
		for (int dx = 0; dx < TILE_COARSEST_STRIDE; dx++) {
			s_block_distance[dy][dx] = distance_of_pixels(sqrt((double)(dx * dx + dy * dy)));
			s_block_level[dy][dx] = (uint8_t)work_coarse_level(preview, s_block_distance[dy][dx]);
		}
	}
}

// Where a sample's value goes: its pixel, the side of the block it fills, and the step of its
// tile's schedule that samples it (whose need it settles).
typedef struct {
	uint16_t x, y;
	uint16_t stride;
	uint16_t step;
} SamplePlace;

// The samples gathered for the next kernel call, all for one kernel variant.
typedef struct {
	KernelPoint points[WORK_CHUNK];
	KernelOrbit orbits[WORK_CHUNK];
	SamplePlace places[WORK_CHUNK];
	int count;
	bool corners;  // iterate with the derivative, and record each sample as a certificate corner
} Chunk;

// A unit's running totals.
typedef struct {
	int64_t work;
	int samples;
	int certified;
} Tally;

void work_context_init(WorkContext* c, Field* field, bool certify, const CertifyContext* cells) {
	const double bailout_sq = (double)MB_BAILOUT_RADIUS * (double)MB_BAILOUT_RADIUS;
	c->field = field;
	c->limits.bailout_sq = bailout_sq;
	c->limits.max_iterations = field->spec.max_iterations;
	c->limits.shape_test = true;
	c->limits.cycle_test = true;
	c->mu_offset = num_log2(num_log2(bailout_sq));
	c->certify = certify;
	c->cells = *cells;
	c->perturbed = field->spec.perturbed;
	c->limits.shape_test = !c->perturbed;  // the closed forms are c's, which a double cannot hold at depth
	memset(&c->reference, 0, sizeof c->reference);
	c->origin_re = c->origin_im = 0.0;
}

void work_context_reference(WorkContext* c, const KernelReference* reference, double origin_re, double origin_im) {
	c->reference = *reference;
	c->origin_re = origin_re;
	c->origin_im = origin_im;
}

// A perturbation kernel's count, at its price in plain lane-iterations.
static int64_t perturbed_price(KernelCount count) {
	return (count.lanes * KERNEL_PERTURBED_COST_PERCENT + count.bla_steps * KERNEL_BLA_STEP_COST_PERCENT) /
	       KERNEL_PERCENT;
}

float work_escape_fraction(const WorkContext* c, double radius_sq) {
	// In exact arithmetic phi lies in (0, 1] — |z_n|^2 lies between the bailout and about its
	// square — and the clamp absorbs the logarithms' last-bit excursions. radius_sq is finite and
	// beyond the bailout (the kernel freezes an orbit at its first escape), so no NaN reaches here.
	const double phi = num_log2(num_log2(radius_sq)) - c->mu_offset;
	return (float)num_clamp(phi, 0.0, 1.0);
}

static Fate fate_of(KernelFate fate) {
	switch (fate) {
	case KERNEL_FATE_INTERIOR_SHAPE: return FATE_SHAPE;
	case KERNEL_FATE_INTERIOR_CYCLE: return FATE_CYCLE;
	case KERNEL_FATE_INTERIOR_LIMIT: return FATE_LIMIT;
	case KERNEL_FATE_ESCAPED:
	default: return FATE_ESCAPED;
	}
}

int work_coarse_level(uint8_t state, uint16_t distance) {
	if (state_rank(state) > RANK_PREVIEW) return TILE_LEVELS;
	int level = 0;
	while (level < TILE_LEVELS && distance <= WORK_REACH[level]) level++;
	return level;
}

// Writes one sample's value: EXACT at its own pixel, and a PREVIEW over
// every pixel of its block whose own value lies farther from it than the sample does — nothing at
// all yet, a coarser sample's (which the finer levels overwrite in turn), or a copy carried from
// farther away; never a nearer copy, nor anything reused or better. The pixel no longer needs work,
// so its step needs one pixel fewer, and the planner's counts of coarse previews follow every pixel
// it writes.
static void store(const WorkContext* c, const SamplePlace* place, int32_t n, float phi, Fate fate, WorkUnit* unit) {
	PixelStore* s = &c->field->front;
	const size_t own = (size_t)place->y * MB_SCREEN_WIDTH + place->x;
	unit->coarse[work_coarse_level(s->state[own], s->distance[own])]--;
	unit->coarse[TILE_LEVELS]++;
	s->n[own] = n;
	s->phi[own] = phi;
	s->state[own] = state_make(RANK_EXACT, fate);
	s->distance[own] = 0;
	unit->need[place->step]--;
	const uint8_t preview = state_make(RANK_PREVIEW, fate);
	for (int dy = 0; dy < place->stride; dy++) {
		const size_t row = (size_t)(place->y + dy) * MB_SCREEN_WIDTH + place->x;
		for (int dx = 0; dx < place->stride; dx++) {
			const uint16_t distance = s_block_distance[dy][dx];
			if (state_rank(s->state[row + dx]) > RANK_PREVIEW || s->distance[row + dx] <= distance) continue;
			unit->coarse[work_coarse_level(s->state[row + dx], s->distance[row + dx])]--;
			unit->coarse[s_block_level[dy][dx]]++;
			s->n[row + dx] = n;
			s->phi[row + dx] = phi;
			s->state[row + dx] = preview;
			s->distance[row + dx] = distance;
		}
	}
}

// Iterates the chunk's samples and stores their values (and, for corners, their certificates).
static void flush(const WorkContext* c, Chunk* chunk, Tally* tally, WorkUnit* unit) {
	if (chunk->count == 0) return;
	// Work is counted in plain lane-iterations, so a quantum buys the same time whichever variant ran.
	if (c->perturbed && chunk->corners) {
		const KernelCount count = kernel_iterate_perturbed_derivative(chunk->points, chunk->orbits, chunk->count,
		                                                              &c->limits, &c->reference);
		tally->work += perturbed_price(count) * KERNEL_DERIVATIVE_COST_PERCENT / KERNEL_PERCENT;
	} else if (c->perturbed) {
		tally->work += perturbed_price(
		    kernel_iterate_perturbed(chunk->points, chunk->orbits, chunk->count, &c->limits, &c->reference));
	} else if (chunk->corners) {
		const int64_t lanes = kernel_iterate_derivative(chunk->points, chunk->orbits, chunk->count, &c->limits);
		tally->work += lanes * KERNEL_DERIVATIVE_COST_PERCENT / KERNEL_PERCENT;
	} else {
		tally->work += kernel_iterate(chunk->points, chunk->orbits, chunk->count, &c->limits);
	}
	tally->samples += chunk->count;
	const int32_t limit = (int32_t)c->limits.max_iterations;
	for (int i = 0; i < chunk->count; i++) {
		const KernelOrbit* o = &chunk->orbits[i];
		const bool escaped = o->fate == KERNEL_FATE_ESCAPED;
		const int32_t n = escaped ? (int32_t)o->iterations : limit;
		const float phi = escaped ? work_escape_fraction(c, o->radius_sq) : 0.0f;
		const SamplePlace* place = &chunk->places[i];
		store(c, place, n, phi, fate_of(o->fate), unit);
		if (chunk->corners) {
			const double mu = (double)n - (double)phi;
			Corner* k = &c->cells.corners[certify_corner_index(place->x, place->y)];
			// Newton for an orbit at the limit (design.md §9.7) costs a corner several times its own
			// orbit, and repays it only where the cell it opens is large: at the coarsest level, whose
			// certified cells spare 255 pixels each. Finer corners at the limit certify nothing.
			const bool at_limit = place->stride == tiles_level_stride(0);
			int64_t steps = 0;
			*k = certify_corner(chunk->points[i].re, chunk->points[i].im, o, mu, c->cells.spacing,
			                    c->perturbed ? &c->reference : NULL, at_limit, &steps);
			tally->work += steps * CERTIFY_STEP_COST_PERCENT / KERNEL_PERCENT;
			k->run = c->cells.run;
		}
	}
	chunk->count = 0;
}

void work_run(void* arg) {
	WorkUnit* unit = (WorkUnit*)arg;
	const WorkContext* c = unit->context;
	Field* f = c->field;
	const TileOrigin origin = tiles_origin(unit->tile);

	Chunk chunk;
	chunk.count = 0;
	chunk.corners = false;
	Tally tally = {0, 0, 0};
	for (int step = unit->first_step; step < unit->last_step; step++) {
		const TileStep* s = tiles_step(step);
		// Entering a level: the level before it is complete here and (the planner saw to it) in the
		// neighbours whose corners its cells share, so its cells can be certified before any of
		// their finer samples is iterated.
		if (c->certify && step == tiles_level_start(s->level) && s->level >= 1 && s->level <= CERT_LEVELS) {
			flush(c, &chunk, &tally, unit);
			tally.certified +=
			    certify_cells(&c->cells, &f->front, origin.x, origin.y, tiles_level_stride(s->level - 1));
		}
		const bool corners = c->certify && s->level < CERT_LEVELS;
		if (corners != chunk.corners) {
			flush(c, &chunk, &tally, unit);
			chunk.corners = corners;
		}
		const int y = origin.y + s->row;
		// A sample's c, or at depth its delta_c: its grid offsets plus the origin's offset from the
		// reference (which is zero in direct mode, where the origin is the plane's).
		const double c_im = c->perturbed ? f->rows.actual[y] + c->origin_im : f->rows.actual[y];
		const uint8_t* state = f->front.state + (size_t)y * MB_SCREEN_WIDTH;
		for (int column = s->first_column; column < MB_TILE_WIDTH; column += s->column_step) {
			const int x = origin.x + column;
			if (!state_needs_work(state[x])) continue;
			if (chunk.count == WORK_CHUNK) flush(c, &chunk, &tally, unit);
			chunk.points[chunk.count].re = c->perturbed ? f->columns.actual[x] + c->origin_re : f->columns.actual[x];
			chunk.points[chunk.count].im = c_im;
			chunk.places[chunk.count].x = (uint16_t)x;
			chunk.places[chunk.count].y = (uint16_t)y;
			chunk.places[chunk.count].stride = s->stride;
			chunk.places[chunk.count].step = (uint16_t)step;
			chunk.count++;
		}
	}
	flush(c, &chunk, &tally, unit);

	// Every sample counted itself off its step's need as it stored its value (a unit at depth may
	// compute a handful of pixels, and recounting its whole tile would cost more than they did);
	// certified cells are counted afresh. Nothing here shades: the look shades the whole frame from
	// the field once the batch is in (look.h).
	if (tally.certified > 0) work_count_needs(&f->front, unit->tile, unit->need, unit->coarse);
	unit->work = tally.work;
	unit->samples = tally.samples;
	unit->certified = tally.certified;
}

void work_count_needs(const PixelStore* store, int tile, int32_t* need, int32_t* coarse) {
	const TileOrigin origin = tiles_origin(tile);
	for (int s = 0; s < TILE_STEPS; s++) need[s] = 0;
	for (int level = 0; level < WORK_COARSE_LEVELS; level++) coarse[level] = 0;
	for (int ty = 0; ty < MB_TILE_HEIGHT; ty++) {
		const size_t row = (size_t)(origin.y + ty) * MB_SCREEN_WIDTH + origin.x;
		const uint8_t* state = store->state + row;
		const uint16_t* distance = store->distance + row;
		const uint8_t* step = tiles_pixel_steps(ty);
		for (int tx = 0; tx < MB_TILE_WIDTH; tx++) {
			if (state_needs_work(state[tx])) need[step[tx]]++;
			coarse[work_coarse_level(state[tx], distance[tx])]++;
		}
	}
}
