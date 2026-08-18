// carts/mandelbrot/kernel_lanes.h
//
// THE LANE MACHINERY, shared by the kernel's variants: the lanes' state, the step, speculation and
// replay, parking and refill, and Brent's cycle test (the design is in kernel.h).
//
// ONE BODY, SEVERAL VARIANTS. Every function here is static inline and takes constant
// `derivative` and `perturbed` flags, and each variant is its own translation unit with exactly
// one call to kernel_drive (kernel.c, kernel_derivative.c, kernel_perturbed.c,
// kernel_perturbed_derivative.c). A function with a single call site is always inlined, so each
// unit compiles its own specialisation with the other variants' arithmetic removed — measured, the
// plain kernel runs exactly as fast as it did before the derivative existed (design.md §4.5). One
// call site for two variants in one unit would not be inlined, and the flag would stay a runtime
// branch in the hottest loops, costing 15%. Each unit also names the groups it keeps in flight
// (LANES_SLOTS), since the derivative's doubled state wants fewer (kernel.h); the answer for a
// point is the same at any shape.
//
// The derivative dz/dc is read and written only by its own statements and never feeds back into
// z, so the escape — and with it the whole speculation argument — is untouched: the derivative
// variant's fates, iterations and |z_n|^2 are the plain kernel's bit for bit (the harness asserts
// it). A derivative that overflows can only do so where z has, inside a block that is then
// restored; one that overflowed before escape (a view far deeper than direct iteration reaches)
// comes out non-finite, and every certificate comparison rejects a non-finite value (certify.c).
//
// PERTURBED lanes hold delta in x, y, their reference index in m and z = Z_m + delta in zx, zy; the
// escape test reads z, the cycle test the state (delta, m). A lane's iteration count is the group's
// clock plus what BLA steps have skipped for it alone (skip); in the plain variants skip is zero
// and every count is the group's clock, as before.
//
// Include it only from a kernel variant's unit, below that unit's FP_CONTRACT pragma, after
// defining LANES_SLOTS.
#ifndef MANDELBROT_KERNEL_LANES_H
#define MANDELBROT_KERNEL_LANES_H

#include "kernel.h"

// The unit that includes this header says how many groups it keeps in flight.
#ifndef LANES_SLOTS
#error "define LANES_SLOTS, the variant's groups in flight, before including kernel_lanes.h"
#endif
#define LANES_WIDTH (KERNEL_GROUP * LANES_SLOTS)

// The fixed point a finished lane is parked on: z = 0 under c = 0 maps to itself exactly (and its
// derivative settles at 1, bounded). A perturbed lane parks with delta_c = -Z_1, so c = 0 there too:
// its z comes back to 0 every step, and the rebase holds it at m = 0.
#define PARKED 0.0

// z_0 = 0 does not depend on c: the derivative starts from zero.
#define DERIVATIVE_START 0.0

// A slot with no group in it.
#define NO_GROUP (-1)

// The main cardioid in closed form: with q = (re - 1/4)^2 + im^2, a point lies inside when
// q (q + (re - 1/4)) <= im^2 / 4. The period-2 bulb is the disc of radius 1/4 about -1.
#define CARDIOID_CUSP 0.25
#define CARDIOID_IM_SQ_SCALE 0.25
#define BULB2_CENTRE_RE (-1.0)
#define BULB2_RADIUS_SQ 0.0625

static inline bool in_cardioid_or_bulb(double re, double im) {
	const double im_sq = im * im;
	const double shifted = re - CARDIOID_CUSP;
	const double q = shifted * shifted + im_sq;
	if (q * (q + shifted) <= im_sq * CARDIOID_IM_SQ_SCALE) return true;
	const double bulb = re - BULB2_CENTRE_RE;
	return bulb * bulb + im_sq <= BULB2_RADIUS_SQ;
}

// Every lane of every slot, as fixed-size arrays; slot s owns lanes [s * KERNEL_GROUP,
// (s + 1) * KERNEL_GROUP).
typedef struct {
	double x[LANES_WIDTH], y[LANES_WIDTH];      // z (plain) or delta (perturbed)
	double dx[LANES_WIDTH], dy[LANES_WIDTH];    // dz/dc (the derivative variants only)
	double cx[LANES_WIDTH], cy[LANES_WIDTH];    // c (plain) or delta_c (perturbed)
	double sx[LANES_WIDTH], sy[LANES_WIDTH];    // x, y saved before a speculative block
	double sdx[LANES_WIDTH], sdy[LANES_WIDTH];  // and the derivative
	double px[LANES_WIDTH], py[LANES_WIDTH];    // Brent's reference state, for the cycle test
	int64_t n[LANES_WIDTH];                      // iterations taken, during a checked replay
	bool live[LANES_WIDTH];                      // still iterating (not parked)

	// The perturbed variants only.
	double zx[LANES_WIDTH], zy[LANES_WIDTH];    // z = Z_m + delta
	double szx[LANES_WIDTH], szy[LANES_WIDTH];  // saved before a speculative block
	int32_t m[LANES_WIDTH], sm[LANES_WIDTH], pm[LANES_WIDTH];  // the reference index: now, saved,
	                                                            // and in Brent's reference state
	int64_t skip[LANES_WIDTH];                   // iterations BLA has taken for this lane alone
	int64_t pskip[LANES_WIDTH];                  // and had taken at Brent's reference state

	// Per slot. Each group keeps its own clock, k, so everything that happens to a lane — when
	// it is tested, when its reference moves — depends on that lane's orbit alone and never on
	// its neighbours in the stream.
	int64_t k[LANES_SLOTS];               // iterations every live lane of the group has taken
	int64_t next_reference[LANES_SLOTS];  // Brent: the reference moves at doubling intervals
	int64_t reference_k[LANES_SLOTS];     // the clock when the reference was taken
	int first[LANES_SLOTS];               // stream index of the group's first sample, or NO_GROUP
	int live_count[LANES_SLOTS];
} Lanes;

// z's coordinates, whichever variant.
static inline double lane_zx(const Lanes* s, int l, bool perturbed) {
	return perturbed ? s->zx[l] : s->x[l];
}
static inline double lane_zy(const Lanes* s, int l, bool perturbed) {
	return perturbed ? s->zy[l] : s->y[l];
}

// The whole of the arithmetic, for one lane. Plain: z <- z^2 + c and, in the derivative variant,
// dz <- 2 z dz + 1 from the z before the step. Perturbed: delta <- (2 Z_m + delta) delta + delta_c
// and z = Z_{m+1} + delta, rebased when |z| < |delta| or the reference ends; the derivative the
// same as the plain one's, from z. The unchecked and the checked steps below both call it, which
// is what lets a replayed block reproduce the speculative one bit for bit. (Results go out through
// pointers to locals: the vectoriser keeps a group's lanes in registers across a replay that way,
// and a struct return cost the plain kernel 8% on boundary-dense views.)
static inline void lane_step(const Lanes* s, int l, bool derivative, bool perturbed, const KernelReference* ref,
                             double* x, double* y, double* dx, double* dy, int32_t* m, double* zx, double* zy) {
	if (derivative) {
		const double ox = lane_zx(s, l, perturbed), oy = lane_zy(s, l, perturbed);
		const double re = ox * s->dx[l] - oy * s->dy[l];
		const double im = ox * s->dy[l] + oy * s->dx[l];
		*dx = re + re + 1.0;
		*dy = im + im;
	}
	if (!perturbed) {
		const double xx = s->x[l] * s->x[l];
		const double yy = s->y[l] * s->y[l];
		const double xy = s->x[l] * s->y[l];
		*y = xy + xy + s->cy[l];
		*x = xx - yy + s->cx[l];
		return;
	}
	const int32_t at = s->m[l];
	const double tr = ref->re[at] + ref->re[at] + s->x[l];
	const double ti = ref->im[at] + ref->im[at] + s->y[l];
	const double nx = tr * s->x[l] - ti * s->y[l] + s->cx[l];
	const double ny = tr * s->y[l] + ti * s->x[l] + s->cy[l];
	const int32_t next = at + 1;
	const double zr = ref->re[next] + nx;
	const double zi = ref->im[next] + ny;
	const bool rebase = (zr * zr + zi * zi < nx * nx + ny * ny) | (next == ref->length);
	*x = rebase ? zr : nx;
	*y = rebase ? zi : ny;
	*m = rebase ? 0 : next;
	*zx = zr;
	*zy = zi;
}

// One unchecked step of every lane of every slot.
static inline void step_all(Lanes* s, bool derivative, bool perturbed, const KernelReference* ref) {
	for (int l = 0; l < LANES_WIDTH; l++) {
		double x, y, dx = 0.0, dy = 0.0, zx = 0.0, zy = 0.0;
		int32_t m = 0;
		lane_step(s, l, derivative, perturbed, ref, &x, &y, &dx, &dy, &m, &zx, &zy);
		s->x[l] = x;
		s->y[l] = y;
		if (derivative) {
			s->dx[l] = dx;
			s->dy[l] = dy;
		}
		if (perturbed) {
			s->m[l] = m;
			s->zx[l] = zx;
			s->zy[l] = zy;
		}
	}
}

// One checked step of one group at iteration k: a lane inside the bailout circle steps and
// records that it has now taken k + 1 iterations; a lane outside it stays frozen where it
// escaped, its derivative with it. A perturbed lane counts its own skipped iterations too, and
// stops at the limit on its own count. (The limits come by value: through a pointer the compiler
// would reload them after every store to the lanes.)
static inline void step_checked(Lanes* s, int base, double bailout_sq, int64_t max_iterations, int64_t k,
                                bool derivative, bool perturbed, const KernelReference* ref) {
	for (int l = base; l < base + KERNEL_GROUP; l++) {
		double x, y, dx = 0.0, dy = 0.0, zx = 0.0, zy = 0.0;
		int32_t m = 0;
		if (!perturbed) {
			// Exactly the Phase 2 form, which the vectoriser turns into one two-lane vector run
			// with selects; any other spelling of it (a branch, a combined flag) was measured 8%
			// slower on boundary-dense views.
			const bool inside = s->x[l] * s->x[l] + s->y[l] * s->y[l] <= bailout_sq;
			lane_step(s, l, derivative, false, ref, &x, &y, &dx, &dy, &m, &zx, &zy);
			s->x[l] = inside ? x : s->x[l];
			s->y[l] = inside ? y : s->y[l];
			if (derivative) {
				s->dx[l] = inside ? dx : s->dx[l];
				s->dy[l] = inside ? dy : s->dy[l];
			}
			s->n[l] = inside ? k + 1 : s->n[l];
			continue;
		}
		const bool inside = (s->zx[l] * s->zx[l] + s->zy[l] * s->zy[l] <= bailout_sq) &
		                    (k + s->skip[l] < max_iterations);
		lane_step(s, l, derivative, true, ref, &x, &y, &dx, &dy, &m, &zx, &zy);
		s->x[l] = inside ? x : s->x[l];
		s->y[l] = inside ? y : s->y[l];
		if (derivative) {
			s->dx[l] = inside ? dx : s->dx[l];
			s->dy[l] = inside ? dy : s->dy[l];
		}
		s->m[l] = inside ? m : s->m[l];
		s->zx[l] = inside ? zx : s->zx[l];
		s->zy[l] = inside ? zy : s->zy[l];
		s->n[l] = inside ? k + 1 + s->skip[l] : s->n[l];
	}
}

static inline bool lane_cycled(const Lanes* s, int l, bool perturbed) {
	const bool same = (s->x[l] == s->px[l]) & (s->y[l] == s->py[l]);
	return s->live[l] & same & (!perturbed || s->m[l] == s->pm[l]);
}

// How one group's speculative block ended: whether every lane stayed inside the bailout circle
// (a parked lane sits at z = 0 and always does), and whether any live lane has come back exactly
// to its reference state. Both are read in one pass over the group, which is what makes the
// cycle test nearly free: it costs two compares a lane beside the bailout test the block needs
// anyway (design.md §4.4 has the measurements).
typedef struct {
	bool inside;
	bool cycled;
} BlockVerdict;

static inline BlockVerdict judge_block(const Lanes* s, int base, double bailout_sq, bool cycle_test, bool perturbed) {
	BlockVerdict v = {true, false};
	for (int l = base; l < base + KERNEL_GROUP; l++) {
		const double zx = lane_zx(s, l, perturbed), zy = lane_zy(s, l, perturbed);
		v.inside &= zx * zx + zy * zy <= bailout_sq;
	}
	if (cycle_test) {
		for (int l = base; l < base + KERNEL_GROUP; l++) v.cycled |= lane_cycled(s, l, perturbed);
	}
	return v;
}

static inline void save_all(Lanes* s, bool derivative, bool perturbed) {
	for (int l = 0; l < LANES_WIDTH; l++) {
		s->sx[l] = s->x[l];
		s->sy[l] = s->y[l];
		if (derivative) {
			s->sdx[l] = s->dx[l];
			s->sdy[l] = s->dy[l];
		}
		if (perturbed) {
			s->szx[l] = s->zx[l];
			s->szy[l] = s->zy[l];
			s->sm[l] = s->m[l];
		}
	}
}

static inline void restore_group(Lanes* s, int base, bool derivative, bool perturbed) {
	for (int l = base; l < base + KERNEL_GROUP; l++) {
		s->x[l] = s->sx[l];
		s->y[l] = s->sy[l];
		if (derivative) {
			s->dx[l] = s->sdx[l];
			s->dy[l] = s->sdy[l];
		}
		if (perturbed) {
			s->zx[l] = s->szx[l];
			s->zy[l] = s->szy[l];
			s->m[l] = s->sm[l];
		}
	}
}

// Parks a lane where it stays for ever: z = 0 under c = 0 (see PARKED).
static inline void park(Lanes* s, int l, bool perturbed, const KernelReference* ref) {
	s->x[l] = s->y[l] = PARKED;
	s->cx[l] = perturbed ? -ref->re[1] : PARKED;
	s->cy[l] = perturbed ? -ref->im[1] : PARKED;
	if (perturbed) {
		s->zx[l] = s->zy[l] = PARKED;
		s->m[l] = 0;
	}
}

// Records a lane's fate and parks it on the fixed point, where it no longer holds back the
// unchecked blocks of the lanes still running.
static inline void retire(Lanes* s, KernelOrbit* orbits, int slot, int l, KernelFate fate, int64_t n, double r2,
                          bool perturbed, const KernelReference* ref) {
	KernelOrbit* orbit = &orbits[s->first[slot] + (l - slot * KERNEL_GROUP)];
	orbit->fate = fate;
	orbit->iterations = n;
	orbit->radius_sq = r2;
	orbit->z_re = lane_zx(s, l, perturbed);
	orbit->z_im = lane_zy(s, l, perturbed);
	orbit->dz_re = s->dx[l];
	orbit->dz_im = s->dy[l];
	// A perturbed lane's own count includes what BLA skipped, at the reference state and now.
	const int64_t own = perturbed ? s->skip[l] - s->pskip[l] : 0;
	orbit->cycle_distance = fate == KERNEL_FATE_INTERIOR_CYCLE ? s->k[slot] - s->reference_k[slot] + own : 0;
	orbit->m = perturbed ? s->m[l] : 0;
	orbit->delta_re = s->x[l];
	orbit->delta_im = s->y[l];
	s->live[l] = false;
	park(s, l, perturbed, ref);
	s->live_count[slot]--;
}

// Loads the group of samples starting at stream index `first` into `slot`, or empties the slot
// when the stream is exhausted. A sample the shape test places inside the set retires at once.
static inline void load_group(Lanes* s, const KernelPoint* points, KernelOrbit* orbits, int count,
                              const KernelLimits* limits, int slot, int first, bool perturbed,
                              const KernelReference* ref) {
	const int base = slot * KERNEL_GROUP;
	s->k[slot] = 0;
	s->next_reference[slot] = KERNEL_BLOCK;
	s->reference_k[slot] = 0;
	s->first[slot] = first < count ? first : NO_GROUP;
	s->live_count[slot] = 0;
	for (int l = base; l < base + KERNEL_GROUP; l++) {
		const int i = first + (l - base);
		const bool carries = s->first[slot] != NO_GROUP && i < count;
		park(s, l, perturbed, ref);
		s->px[l] = s->py[l] = PARKED;
		s->dx[l] = s->dy[l] = DERIVATIVE_START;
		if (perturbed) {
			s->pm[l] = 0;
			s->skip[l] = s->pskip[l] = 0;
		}
		if (carries) {
			s->cx[l] = points[i].re;
			s->cy[l] = points[i].im;
		}
		s->live[l] = carries;
		s->live_count[slot] += carries;
		if (carries && limits->shape_test && in_cardioid_or_bulb(s->cx[l], s->cy[l])) {
			retire(s, orbits, slot, l, KERNEL_FATE_INTERIOR_SHAPE, limits->max_iterations, 0.0, perturbed, ref);
		}
	}
}

// Brent's cycle test for one group at a block boundary: retire each live lane that has come back
// exactly to its reference state, then, at doubling intervals of the group's clock, make the
// present state the new reference. A state equal to an EARLIER state is all the proof needs.
static inline void test_cycles(Lanes* s, KernelOrbit* orbits, const KernelLimits* limits, int slot, bool perturbed,
                               const KernelReference* ref) {
	const int base = slot * KERNEL_GROUP;
	for (int l = base; l < base + KERNEL_GROUP; l++) {
		if (lane_cycled(s, l, perturbed)) {
			retire(s, orbits, slot, l, KERNEL_FATE_INTERIOR_CYCLE, limits->max_iterations, 0.0, perturbed, ref);
		}
	}
	if (s->k[slot] >= s->next_reference[slot]) {
		for (int l = base; l < base + KERNEL_GROUP; l++) {
			s->px[l] = s->x[l];
			s->py[l] = s->y[l];
			if (perturbed) {
				s->pm[l] = s->m[l];
				s->pskip[l] = s->skip[l];
			}
		}
		s->reference_k[slot] = s->k[slot];
		s->next_reference[slot] = s->k[slot] + s->k[slot];
	}
}

// A lane's own iteration count at the group's clock.
static inline int64_t lane_count(const Lanes* s, int l, int64_t k, bool perturbed) {
	return k + (perturbed ? s->skip[l] : 0);
}

// Replays `span` iterations of one group checked, then retires what finished. The scalar loop
// tests |z_n| only while n < max_iterations, so a lane that reaches the limit is interior even
// if that last z lies beyond the bailout. Returns the lane-iterations spent.
static inline int64_t replay_group(Lanes* s, KernelOrbit* orbits, const KernelLimits* limits, int slot,
                                   int64_t span, bool derivative, bool perturbed, const KernelReference* ref) {
	const int base = slot * KERNEL_GROUP;
	const int64_t k = s->k[slot];
	for (int l = base; l < base + KERNEL_GROUP; l++) s->n[l] = lane_count(s, l, k, perturbed);
	const double bailout_sq = limits->bailout_sq;
	const int64_t max_iterations = limits->max_iterations;
	for (int64_t u = 0; u < span; u++) {
		step_checked(s, base, bailout_sq, max_iterations, k + u, derivative, perturbed, ref);
	}
	s->k[slot] = k + span;
	for (int l = base; l < base + KERNEL_GROUP; l++) {
		if (!s->live[l]) continue;
		const double zx = lane_zx(s, l, perturbed), zy = lane_zy(s, l, perturbed);
		const double r2 = zx * zx + zy * zy;
		if (s->n[l] < limits->max_iterations && r2 > limits->bailout_sq) {
			retire(s, orbits, slot, l, KERNEL_FATE_ESCAPED, s->n[l], r2, perturbed, ref);
		} else if (lane_count(s, l, s->k[slot], perturbed) >= limits->max_iterations) {
			retire(s, orbits, slot, l, KERNEL_FATE_INTERIOR_LIMIT, limits->max_iterations, 0.0, perturbed, ref);
		}
	}
	return span * KERNEL_GROUP;
}

// --- BLA ------------------------------------------------------------------------------------------
// A block's validity radius is never more than its first step's: a merge takes the smaller of its
// first half's and what its second half allows (reference.c). So a lane whose single step at m is
// invalid can take no block there, and one test spares it the walk down every level — which, at
// an m where Z_m is tiny (every rebase lands at m = 0, where Z_0 = 0), was a load from each level's
// run, every block boundary, for nothing. On views far above a deep nucleus, where BLA seldom
// applies, it halves the level tests and brings the lanes from 1.2-1.3 times their price to
// 1.0-1.2. The test allows BLA_MONOTONE_SLACK over the step's own radius, so a merged radius an ulp
// above it (a square root squared) still gets the full walk: the blocks a lane takes are exactly
// those it took before.
#define BLA_MONOTONE_SLACK (1.0 + 0x1p-40)

// The longest valid block of the table starting at a lane's m, applied, and again, until none is
// valid, the table ends, or `budget` iterations are used. A block of 2^l steps starts only at an
// m aligned to it, so a lane climbs: single steps (level 0, computed from Z_m) to an even m, a
// level-1 block to a multiple of four, and so on. Returns the steps taken.
static inline int64_t bla_advance(Lanes* s, int l, const KernelReference* ref, int64_t budget, bool derivative) {
	int64_t taken = 0;
	const double radius = 2.0 * ref->epsilon;
	for (;;) {
		const int32_t m = s->m[l];
		if (m >= ref->bla_length) return taken;
		const double x = s->x[l], y = s->y[l];
		const double d2 = x * x + y * y;
		const double zr0 = ref->re[m], zi0 = ref->im[m];
		const double z2 = zr0 * zr0 + zi0 * zi0;
		const double step_r2 = z2 > BLA_Z_LIMIT_SQ ? 0.0 : radius * radius * z2;
		if (!(d2 < step_r2 * BLA_MONOTONE_SLACK)) return taken;
		int level = ref->levels;
		if (m) {
			int aligned = 0;
			while (!((m >> aligned) & 1)) aligned++;
			if (aligned < level) level = aligned;
		}
		BlaEntry e = {0.0, 0.0, 0.0, 0.0, 0.0};
		int32_t span = 0;
		for (; level >= 0; level--) {
			const int32_t steps = (int32_t)1 << level;
			if (steps > budget || m + steps > ref->bla_length) continue;
			if (level == 0) {
				e.ar = zr0 + zr0;
				e.ai = zi0 + zi0;
				e.br = 1.0;
				e.bi = 0.0;
				e.r2 = step_r2;
			} else {
				e = ref->bla[level][m >> level];
			}
			if (d2 < e.r2) {
				span = steps;
				break;
			}
		}
		if (!span) return taken;
		const double cx = s->cx[l], cy = s->cy[l];
		const double nx = e.ar * x - e.ai * y + (e.br * cx - e.bi * cy);
		const double ny = e.ar * y + e.ai * x + (e.br * cy + e.bi * cx);
		if (derivative) {
			const double ddx = e.ar * s->dx[l] - e.ai * s->dy[l] + e.br;
			const double ddy = e.ar * s->dy[l] + e.ai * s->dx[l] + e.bi;
			s->dx[l] = ddx;
			s->dy[l] = ddy;
		}
		const int32_t next = m + span;
		const double zr = ref->re[next] + nx, zi = ref->im[next] + ny;
		const bool rebase = (zr * zr + zi * zi < nx * nx + ny * ny) | (next == ref->length);
		s->x[l] = rebase ? zr : nx;
		s->y[l] = rebase ? zi : ny;
		s->m[l] = rebase ? 0 : next;
		s->zx[l] = zr;
		s->zy[l] = zi;
		s->skip[l] += span;
		budget -= span;
		taken++;
	}
}

// The driver every variant shares; `derivative` and `perturbed` are constants at each unit's one
// call site, and `ref` is NULL for the plain variants.
static inline KernelCount kernel_drive(const KernelPoint* points, KernelOrbit* orbits, int count,
                                       const KernelLimits* limits, bool derivative, bool perturbed,
                                       const KernelReference* ref) {
	Lanes s;
	int next = 0;
	KernelCount work = {0, 0};
	const bool bla = perturbed && ref->levels > 0;
	for (int slot = 0; slot < LANES_SLOTS; slot++) {
		load_group(&s, points, orbits, count, limits, slot, next, perturbed, ref);
		next += KERNEL_GROUP;
	}

	for (;;) {
		// Refill every slot whose group has finished; stop when the stream and the slots are
		// both empty.
		bool busy = false;
		for (int slot = 0; slot < LANES_SLOTS; slot++) {
			while (s.first[slot] != NO_GROUP && s.live_count[slot] == 0) {
				load_group(&s, points, orbits, count, limits, slot, next, perturbed, ref);
				next += KERNEL_GROUP;
			}
			busy |= s.first[slot] != NO_GROUP;
		}
		if (!busy) break;

		// BLA first: each live lane skips what the table allows, then the group's furthest lane
		// decides whether a whole block fits before the limit.
		bool speculate[LANES_SLOTS];
		bool any = false;
		for (int slot = 0; slot < LANES_SLOTS; slot++) {
			int64_t furthest = s.k[slot];
			if (perturbed) {
				for (int l = slot * KERNEL_GROUP; l < (slot + 1) * KERNEL_GROUP; l++) {
					if (!s.live[l]) continue;
					if (bla) {
						const int64_t budget = limits->max_iterations - lane_count(&s, l, s.k[slot], true);
						work.bla_steps += bla_advance(&s, l, ref, budget, derivative);
					}
					const int64_t count_l = lane_count(&s, l, s.k[slot], true);
					if (count_l > furthest) furthest = count_l;
				}
			}
			speculate[slot] = s.live_count[slot] > 0 && furthest + KERNEL_BLOCK <= limits->max_iterations;
			any |= speculate[slot];
		}
		if (any) {
			save_all(&s, derivative, perturbed);
			for (int u = 0; u < KERNEL_BLOCK; u++) step_all(&s, derivative, perturbed, ref);
			work.lanes += (int64_t)KERNEL_BLOCK * LANES_WIDTH;
		}

		// Each group commits its block, or restores and replays it: a group whose block failed
		// replays it in full, and a group near the limit finishes what is left checked. Either
		// way the group ends at a block boundary of its own clock, where the cycle test runs.
		for (int slot = 0; slot < LANES_SLOTS; slot++) {
			if (s.live_count[slot] == 0) continue;
			const int base = slot * KERNEL_GROUP;
			BlockVerdict verdict = {false, false};
			if (speculate[slot]) {
				verdict = judge_block(&s, base, limits->bailout_sq, limits->cycle_test, perturbed);
			}
			if (verdict.inside) {
				s.k[slot] += KERNEL_BLOCK;
				if (verdict.cycled || (limits->cycle_test && s.k[slot] >= s.next_reference[slot])) {
					test_cycles(&s, orbits, limits, slot, perturbed, ref);
				}
				continue;
			}
			if (any) restore_group(&s, base, derivative, perturbed);
			// A perturbed group replays a whole block even where one lane is near the limit: that
			// lane stops on its own count, and the others go on.
			const int64_t remaining = limits->max_iterations - s.k[slot];
			const int64_t span = remaining < KERNEL_BLOCK ? remaining : KERNEL_BLOCK;
			work.lanes += replay_group(&s, orbits, limits, slot, span, derivative, perturbed, ref);
			if (limits->cycle_test && s.live_count[slot] > 0) test_cycles(&s, orbits, limits, slot, perturbed, ref);
		}
	}
	return work;
}

#endif // MANDELBROT_KERNEL_LANES_H
