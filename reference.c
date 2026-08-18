// carts/mandelbrot/reference.c
//
// The orbit in the bignum, the BLA table in doubles, and Newton's method in perturbation form.
// Every loop here is bounded by a budget or a named maximum, so a slice's cost is known before it
// starts, and every result depends on the inputs alone: the same reference is built to the same
// bits on every runtime and at any slicing.
#pragma STDC FP_CONTRACT OFF

#include "reference.h"

#include "config.h"

#define BAILOUT ((double)MB_BAILOUT_RADIUS)
#define BAILOUT_SQ (BAILOUT * BAILOUT)
#define NEWTON_STEPS 64  // quadratic convergence needs a handful; this many is divergence

// A square of an (f + 1)-limb magnitude takes f + 1 diagonal and (f + 1) f / 2 cross products,
// and an iteration takes three squares: x^2, y^2 and (x + y)^2 (2xy = the last less the others),
// and its fixed work.
#define SQUARES_PER_ITERATION 3

static int64_t iteration_products(int frac) {
	const int64_t n = frac + 1;
	return SQUARES_PER_ITERATION * (n + n * (n - 1) / 2) + REFERENCE_ITERATION_OVERHEAD;
}

uint32_t reference_bytes(void) {
	return REFERENCE_BYTES;
}

bool reference_reserve(Reference* r) {
	memset(r, 0, sizeof *r);
	r->re = (double*)jet_memory_reserve((uint32_t)((REFERENCE_CAPACITY + 1) * sizeof(double)));
	r->im = (double*)jet_memory_reserve((uint32_t)((REFERENCE_CAPACITY + 1) * sizeof(double)));
	r->table = (BlaEntry*)jet_memory_reserve((uint32_t)(REFERENCE_CAPACITY * sizeof(BlaEntry)));
	r->stage = REFERENCE_EMPTY;
	return r->re && r->im && r->table;
}

void reference_begin(Reference* r, const DeepPoint* c, int frac, int32_t target, int32_t period, double max_delta_c,
                     double epsilon) {
	Big cr, ci;
	point_coordinates(c, &cr, &ci);
	big_truncate(&cr, frac);
	big_truncate(&ci, frac);
	r->c = point_from_coordinates(&cr, &ci);
	r->frac = frac;
	r->period = period;
	r->target = period > 0 ? period : target;
	if (r->target > REFERENCE_CAPACITY) r->target = REFERENCE_CAPACITY;
	if (r->target < 1) r->target = 1;
	r->max_delta_c = max_delta_c;
	r->epsilon = epsilon;
	r->stage = REFERENCE_ORBIT;
	r->zr = big_zero();
	r->zi = big_zero();
	r->computed = 0;
	r->escaped = false;
	r->re[0] = r->im[0] = 0.0;
	r->record_count = 0;
	r->record_sq = BAILOUT_SQ;
	r->levels = 0;
	r->table_level = 1;
	r->table_entry = 0;
	r->dzr = r->dzi = 0.0;
	r->refinements = 0;
}

// A nucleus's orbit is done: Newton's step from Z_p and dZ_p/dc, and another pass if it was not
// negligible. Returns true when the orbit must be iterated again from the refined point.
static bool refine_nucleus(Reference* r) {
	if (r->period <= 0 || r->refinements == REFERENCE_REFINEMENTS) return false;
	const double zr = r->re[r->computed], zi = r->im[r->computed];
	const double d2 = r->dzr * r->dzr + r->dzi * r->dzi;
	const double sr = (zr * r->dzr + zi * r->dzi) / d2, si = (zi * r->dzr - zr * r->dzi) / d2;
	if (!(sr == sr && si == si)) return false;
	// Negligible below 2^-REFINE_BITS of the finest spacing the reference serves.
	double negligible = 1.0;
	for (int b = reference_deepest_log2(r) - REFERENCE_REFINE_BITS; b < 0; b++) negligible *= 0.5;
	if (sr * sr + si * si <= negligible * negligible) return false;
	point_offset(&r->c, -sr, -si);
	Big cr, ci;
	point_coordinates(&r->c, &cr, &ci);
	big_truncate(&cr, r->frac);
	big_truncate(&ci, r->frac);
	r->c = point_from_coordinates(&cr, &ci);
	r->refinements++;
	r->zr = big_zero();
	r->zi = big_zero();
	r->computed = 0;
	r->escaped = false;
	r->dzr = r->dzi = 0.0;
	r->record_count = 0;
	r->record_sq = BAILOUT_SQ;
	return true;
}

// |x| >= 1024: the coordinate has escaped, and its square could overflow the integer limb.
static bool beyond_bailout(const Big* x) {
	const int32_t integer = (int32_t)x->limb[BIG_INTEGER];
	return integer >= MB_BAILOUT_RADIUS || integer < -MB_BAILOUT_RADIUS;
}

static void orbit_slice(Reference* r, int64_t products) {
	const int f = r->frac;
	const int64_t cost = iteration_products(f);
	Big cr, ci, xx, yy, sum, r2, bail = big_zero();
	point_coordinates(&r->c, &cr, &ci);  // read afresh: a refinement moves it
	bail.limb[BIG_INTEGER] = (uint32_t)(MB_BAILOUT_RADIUS * MB_BAILOUT_RADIUS);
	while (r->computed < r->target && products >= cost) {
		if (beyond_bailout(&r->zr) || beyond_bailout(&r->zi)) {
			r->escaped = true;
			break;
		}
		big_square(&xx, &r->zr, f);
		big_square(&yy, &r->zi, f);
		big_add(&r2, &xx, &yy, f);
		if (big_compare(&r2, &bail) > 0) {
			r->escaped = true;
			break;
		}
		// dZ/dc for a nucleus's refinement, in doubles: dZ <- 2 Z dZ + 1, from the Z before the step.
		const double pr = r->re[r->computed], pi = r->im[r->computed];
		const double dre = pr * r->dzr - pi * r->dzi, dim = pr * r->dzi + pi * r->dzr;
		r->dzr = dre + dre + 1.0;
		r->dzi = dim + dim;
		big_add(&sum, &r->zr, &r->zi, f);
		big_square(&sum, &sum, f);
		big_sub(&sum, &sum, &xx, f);
		big_sub(&sum, &sum, &yy, f);
		big_add(&r->zi, &sum, &ci, f);
		big_sub(&sum, &xx, &yy, f);
		big_add(&r->zr, &sum, &cr, f);
		const int32_t n = ++r->computed;
		r->re[n] = big_to_double(&r->zr);
		r->im[n] = big_to_double(&r->zi);
		const double d2 = r->re[n] * r->re[n] + r->im[n] * r->im[n];
		if (d2 < r->record_sq) {
			r->record_sq = d2;
			if (r->record_count == REFERENCE_RECORDS) {
				memmove(r->record_at, r->record_at + 1, (REFERENCE_RECORDS - 1) * sizeof(int32_t));
				r->record_count--;
			}
			r->record_at[r->record_count++] = n;
		}
		products -= cost;
	}
	if (r->computed < r->target && !r->escaped) return;
	if (!r->escaped && refine_nucleus(r)) return;
	// The orbit is done: lay out the table's levels, or finish without one.
	int32_t entries = 0;
	r->levels = 0;
	if (r->epsilon > 0.0) {
		for (int level = 1; level <= KERNEL_BLA_LEVELS && (r->computed >> level) > 0; level++) {
			r->level_start[level] = r->table + entries;
			r->level_count[level] = r->computed >> level;
			entries += r->level_count[level];
			r->levels = level;
		}
	}
	r->stage = r->levels > 0 ? REFERENCE_TABLE : REFERENCE_READY;
}

// --- The BLA table ---------------------------------------------------------------------------------
// One step at m: delta <- 2 Z_m delta + delta_c exactly but for the dropped delta^2, which stays
// within epsilon of the kept 2 Z_m delta while |delta| <= epsilon |2 Z_m| (design.md §9.3). A step
// where |Z_m| > 2 is never approximated: the orbit escapes within a few steps of it.
static BlaEntry single_step(const Reference* r, int32_t m) {
	const double zr = r->re[m], zi = r->im[m];
	const double z2 = zr * zr + zi * zi;
	const double radius = 2.0 * r->epsilon;
	BlaEntry e;
	e.ar = zr + zr;
	e.ai = zi + zi;
	e.br = 1.0;
	e.bi = 0.0;
	e.r2 = z2 > BLA_Z_LIMIT_SQ ? 0.0 : radius * radius * z2;
	return e;
}

// x then y: A = A_y A_x, B = A_y B_x + B_y, and the merged block is valid while x is and while x's
// result lies within y's radius for every |delta_c| up to the table's bound:
// r = min(r_x, max(0, (r_y - |B_x| max|delta_c|) / |A_x|)).
static BlaEntry merge(const BlaEntry* x, const BlaEntry* y, double max_delta_c) {
	BlaEntry e;
	e.ar = y->ar * x->ar - y->ai * x->ai;
	e.ai = y->ar * x->ai + y->ai * x->ar;
	e.br = (y->ar * x->br - y->ai * x->bi) + y->br;
	e.bi = (y->ar * x->bi + y->ai * x->br) + y->bi;
	const double rx = sqrt(x->r2), ry = sqrt(y->r2);
	const double ax = sqrt(x->ar * x->ar + x->ai * x->ai);
	const double bx = sqrt(x->br * x->br + x->bi * x->bi);
	double t = (ry - bx * max_delta_c) / ax;
	if (!(t > 0.0)) t = 0.0;  // negative, or 0/0: the block is never valid
	const double r = rx < t ? rx : t;
	e.r2 = r * r;
	return e;
}

static void table_slice(Reference* r, int32_t entries) {
	while (entries > 0 && r->table_level <= r->levels) {
		const int level = r->table_level;
		BlaEntry* out = r->level_start[level];
		while (entries > 0 && r->table_entry < r->level_count[level]) {
			const int32_t j = r->table_entry;
			if (level == 1) {
				const BlaEntry x = single_step(r, 2 * j), y = single_step(r, 2 * j + 1);
				out[j] = merge(&x, &y, r->max_delta_c);
			} else {
				const BlaEntry* below = r->level_start[level - 1];
				out[j] = merge(&below[2 * j], &below[2 * j + 1], r->max_delta_c);
			}
			r->table_entry++;
			entries--;
		}
		if (r->table_entry == r->level_count[level]) {
			r->table_level++;
			r->table_entry = 0;
		}
	}
	if (r->table_level > r->levels) r->stage = REFERENCE_READY;
}

ReferenceStage reference_advance(Reference* r, int64_t products, int32_t entries) {
	if (r->stage == REFERENCE_ORBIT) {
		orbit_slice(r, products);
	} else if (r->stage == REFERENCE_TABLE) {
		table_slice(r, entries);
	}
	return r->stage;
}


KernelReference reference_kernel(const Reference* r) {
	KernelReference k;
	memset(&k, 0, sizeof k);
	k.re = r->re;
	k.im = r->im;
	k.length = r->computed;
	k.levels = r->levels;
	for (int level = 1; level <= r->levels; level++) k.bla[level] = r->level_start[level];
	k.bla_length = r->computed;
	k.epsilon = r->epsilon;
	k.period = r->period;
	return k;
}

int reference_deepest_log2(const Reference* r) {
	return -(r->frac * BIG_LIMB_BITS - REFERENCE_GUARD_BITS);
}

int reference_limbs_for(double deepest_log2) {
	int bits = REFERENCE_GUARD_BITS;
	while (-(double)(bits - REFERENCE_GUARD_BITS) > deepest_log2) bits++;
	return big_limbs_for_bits(bits);
}

// --- Searches --------------------------------------------------------------------------------------
// One perturbation step with rebasing against the reference (kernel_lanes.h's step, scalar).
static void perturb_step(const Reference* r, double dc_re, double dc_im, double* x, double* y, int32_t* m) {
	const double tr = r->re[*m] + r->re[*m] + *x, ti = r->im[*m] + r->im[*m] + *y;
	const double nx = tr * *x - ti * *y + dc_re, ny = tr * *y + ti * *x + dc_im;
	const int32_t next = *m + 1;
	const double zr = r->re[next] + nx, zi = r->im[next] + ny;
	const bool rebase = zr * zr + zi * zi < nx * nx + ny * ny || next == r->computed;
	*x = rebase ? zr : nx;
	*y = rebase ? zi : ny;
	*m = rebase ? 0 : next;
}

void reference_period_begin(PeriodSearch* s, double dc_re, double dc_im, int32_t limit) {
	memset(s, 0, sizeof *s);
	s->dc_re = dc_re;
	s->dc_im = dc_im;
	s->limit = limit;
	s->record_sq = BAILOUT_SQ;
}

SearchStatus reference_period_advance(PeriodSearch* s, const Reference* r, int64_t steps) {
	while (steps-- > 0) {
		if (s->n >= s->limit) return SEARCH_FOUND;
		perturb_step(r, s->dc_re, s->dc_im, &s->x, &s->y, &s->m);
		s->n++;
		const double zr = r->re[s->m] + s->x, zi = r->im[s->m] + s->y;
		const double z2 = zr * zr + zi * zi;
		if (z2 > BAILOUT_SQ) return SEARCH_FOUND;
		if (z2 < s->record_sq) {
			s->record_sq = z2;
			if (s->record_count == REFERENCE_RECORDS) {
				memmove(s->record_at, s->record_at + 1, (REFERENCE_RECORDS - 1) * sizeof(int32_t));
				s->record_count--;
			}
			s->record_at[s->record_count++] = s->n;
		}
	}
	return SEARCH_RUNNING;
}

void reference_nucleus_begin(NucleusSearch* s, int32_t period, double start_re, double start_im, double tolerance,
                             double reach) {
	memset(s, 0, sizeof *s);
	s->period = period;
	s->start_re = s->cr = start_re;
	s->start_im = s->ci = start_im;
	s->tolerance = tolerance;
	s->reach = reach;
}

SearchStatus reference_nucleus_advance(NucleusSearch* s, const Reference* r, int64_t steps) {
	while (steps-- > 0) {
		if (s->n < s->period) {
			// z_n and dz_n/dc at c = reference + estimate, by perturbation with rebasing.
			const double zr = r->re[s->m] + s->x, zi = r->im[s->m] + s->y;
			const double re = zr * s->dx - zi * s->dy, im = zr * s->dy + zi * s->dx;
			s->dx = re + re + 1.0;
			s->dy = im + im;
			perturb_step(r, s->cr, s->ci, &s->x, &s->y, &s->m);
			s->n++;
			s->zr = r->re[s->m] + s->x;
			s->zi = r->im[s->m] + s->y;
			if (s->zr * s->zr + s->zi * s->zi > BAILOUT_SQ) return SEARCH_FAILED;  // not near this period's nucleus
			continue;
		}
		// Newton: c <- c - z_p / (dz_p/dc).
		const double d2 = s->dx * s->dx + s->dy * s->dy;
		const double sr = (s->zr * s->dx + s->zi * s->dy) / d2, si = (s->zi * s->dx - s->zr * s->dy) / d2;
		if (!(sr == sr && si == si)) return SEARCH_FAILED;  // the derivative vanished
		s->cr -= sr;
		s->ci -= si;
		const double wr = s->cr - s->start_re, wi = s->ci - s->start_im;
		if (!(wr * wr + wi * wi <= s->reach * s->reach)) return SEARCH_FAILED;
		if (sr * sr + si * si <= s->tolerance * s->tolerance) return SEARCH_FOUND;
		if (++s->step == NEWTON_STEPS) return SEARCH_FAILED;
		s->x = s->y = s->dx = s->dy = 0.0;
		s->m = 0;
		s->n = 0;
	}
	return SEARCH_RUNNING;
}

bool reference_find_nucleus(const Reference* r, int32_t period, double start_re, double start_im, double tolerance,
                            double reach, double* nucleus_re, double* nucleus_im) {
	NucleusSearch s;
	reference_nucleus_begin(&s, period, start_re, start_im, tolerance, reach);
	SearchStatus status;
	while ((status = reference_nucleus_advance(&s, r, INT64_MAX)) == SEARCH_RUNNING) {
	}
	*nucleus_re = s.cr;
	*nucleus_im = s.ci;
	return status == SEARCH_FOUND;
}
