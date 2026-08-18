// carts/mandelbrot/certify.c
//
// Every comparison that grants a certificate is written so that a NaN fails it: a corner whose
// derivative overflowed, or whose cycle could not be resolved, simply certifies nothing.
#pragma STDC FP_CONTRACT OFF

#include "certify.h"

#include "numerics.h"

#define LN_2 0.69314718055994530942
#define INV_LN_2 1.44269504088896340736
#define QUARTER 0.25
#define HALF_DIAGONAL 0.70710678118654752440  // 1/sqrt 2: a cell's farthest point from its corners

// Below this potential, (1 - e^(-2G)) / (2G) is 1 - G to within 1e-12, and the direct formula would
// lose its digits to cancellation.
#define POTENTIAL_SERIES_LIMIT 1e-6

// The least period is where the orbit first comes back to z0 to within this fraction of |z0|
// (2^-30): far above the float cycle's own jitter, far below the spacing of distinct cycle points
// anywhere a certificate could be granted (near a bifurcation, where points crowd, b_int -> 0).
#define PERIOD_TOLERANCE_SQ 0x1p-60

// A certified exterior cell may hold no pixel the limit would call interior. Every pixel lies
// within r = s/sqrt 2 of a corner whose M-free disc has radius R >= m r (m the exterior margin), and
// the escape potential G is positive and harmonic on that disc, so by Harnack's inequality
// G(pixel) >= G(corner) (R - r) / (R + r): mu = const - log2 G rises across the cell by at most
// log2((m + 1) / (m - 1)). An escaped pixel has n <= mu + 1 and n <= limit - 1, so a corner is
// admitted only below limit - ESCAPE_HEADROOM - log2((m + 1) / (m - 1)). (A margin at or below 1
// makes the logarithm a NaN, which fails the comparison: no exterior cell is certified.)
#define ESCAPE_HEADROOM 2.0

typedef struct {
	double re, im;
} Complex;

static inline Complex cx_make(double re, double im) {
	Complex z = {re, im};
	return z;
}
static inline Complex cx_add(Complex a, Complex b) {
	return cx_make(a.re + b.re, a.im + b.im);
}
static inline Complex cx_sub(Complex a, Complex b) {
	return cx_make(a.re - b.re, a.im - b.im);
}
static inline Complex cx_mul(Complex a, Complex b) {
	return cx_make(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static inline Complex cx_twice(Complex a) {
	return cx_make(a.re + a.re, a.im + a.im);
}
static inline double cx_norm(Complex a) {
	return a.re * a.re + a.im * a.im;
}
static inline Complex cx_div(Complex a, Complex b) {
	const double d = cx_norm(b);
	return cx_make((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
// The principal square root, from real square roots (correctly rounded everywhere).
static inline Complex cx_sqrt(Complex a) {
	const double r = sqrt(cx_norm(a));
	const double re = sqrt((r + a.re) * 0.5);
	const double im = sqrt((r - a.re) * 0.5);
	return cx_make(re, a.im < 0.0 ? -im : im);
}

static Corner no_corner(void) {
	Corner k;
	memset(&k, 0, sizeof k);
	k.kind = CORNER_NONE;
	return k;
}

// --- The exterior certificate ------------------------------------------------------------------
static Corner exterior_corner(const KernelOrbit* o, double mu, double spacing) {
	const Complex z = cx_make(o->z_re, o->z_im);
	const Complex dz = cx_make(o->dz_re, o->dz_im);
	const double r2 = cx_norm(z);
	const double ln_r = 0.5 * num_log2(r2) * LN_2;  // ln|z_n|
	const double dz_abs = sqrt(cx_norm(dz));
	const double b = 2.0 * sqrt(r2) * ln_r / dz_abs;  // the distance estimate, in the plane

	// G = 2^-n ln|z_n|; the factor (1 - e^(-2G)) / (2G) turns b/4 into Koebe's exact lower bound.
	const double potential = ln_r * num_exp2(-(double)o->iterations);
	const double factor = potential < POTENTIAL_SERIES_LIMIT
	                          ? 1.0 - potential
	                          : (1.0 - num_exp2(-2.0 * potential * INV_LN_2)) / (2.0 * potential);

	// grad mu = -(Re w, -Im w) / (ln 2 ln|z|) per unit of c, with w = z'/z; per pixel, with the
	// screen's y running down the imaginary axis.
	const Complex w = cx_div(dz, z);
	const double k = -spacing / (LN_2 * ln_r);
	Corner c;
	memset(&c, 0, sizeof c);
	c.kind = CORNER_EXTERIOR;
	c.mu = mu;
	c.gx = (float)(k * w.re);
	c.gy = (float)(k * w.im);
	c.radius = (float)(QUARTER * b * factor / spacing);
	return c;
}

// --- The interior certificate ------------------------------------------------------------------
// From z0 on the attracting cycle of period p: the derivatives of f^p at z0, by the recurrences
// (in this order each step) dcdz <- 2(z dcdz + dz dc); dc <- 2 z dc + 1; dzdz <- 2(z dzdz + dz^2);
// dz <- 2 z dz; z <- z^2 + c; then b_int.
// The certificate from the derivatives of f^p at z0: the multiplier dz, dc, dzdz and dcdz.
static Corner interior_disc(Complex dz, Complex dc, Complex dzdz, Complex dcdz, double spacing) {
	const double multiplier = cx_norm(dz);
	if (!(multiplier < 1.0)) return no_corner();
	const Complex denominator = cx_add(dcdz, cx_div(cx_mul(dzdz, dc), cx_sub(cx_make(1.0, 0.0), dz)));
	const double b = (1.0 - multiplier) / sqrt(cx_norm(denominator));
	Corner k = no_corner();
	k.kind = CORNER_INTERIOR;
	k.radius = (float)(QUARTER * b / spacing);
	return k;
}

static Corner interior_corner(Complex c, Complex z0, int64_t period, double spacing, int64_t* steps) {
	*steps += period;
	Complex z = z0, dz = cx_make(1.0, 0.0), dc = cx_make(0.0, 0.0);
	Complex dzdz = cx_make(0.0, 0.0), dcdz = cx_make(0.0, 0.0);
	for (int64_t i = 0; i < period; i++) {
		dcdz = cx_twice(cx_add(cx_mul(z, dcdz), cx_mul(dz, dc)));
		dc = cx_add(cx_twice(cx_mul(z, dc)), cx_make(1.0, 0.0));
		dzdz = cx_twice(cx_add(cx_mul(z, dzdz), cx_mul(dz, dz)));
		dz = cx_twice(cx_mul(z, dz));
		z = cx_add(cx_mul(z, z), c);
	}
	return interior_disc(dz, dc, dzdz, dcdz, spacing);
}

// The least period of a cycle caught at z0 after `distance` iterations: the float cycle's period
// divides the distance, and the true period divides the float cycle's (a weakly attracting cycle
// can be shadowed by a float cycle of a multiple of its period). The derivative formula needs the
// true one: at a multiple k of it the multiplier is lambda^k, whose map is not univalent, and
// the estimate would be too large near the component's centre.
static int64_t least_period(Complex c, Complex z0, int64_t distance, int64_t* steps) {
	const double tolerance = PERIOD_TOLERANCE_SQ * (1.0 + cx_norm(z0));
	Complex z = z0;
	for (int64_t q = 1; q <= distance; q++) {
		z = cx_add(cx_mul(z, z), c);
		if (cx_norm(cx_sub(z, z0)) <= tolerance) {
			*steps += q;
			return q;
		}
	}
	*steps += distance;
	return 0;
}

// Inside the main cardioid the attracting fixed point is (1 - sqrt(1 - 4c)) / 2; inside the period-2
// bulb the cycle solves z^2 + z + c + 1 = 0, z = (-1 + sqrt(-3 - 4c)) / 2.
#define BULB2_CENTRE_RE (-1.0)
#define BULB2_RADIUS_SQ 0.0625
#define CARDIOID_PERIOD 1
#define BULB2_PERIOD 2

static Corner shape_corner(Complex c, double spacing, int64_t* steps) {
	const double bulb = c.re - BULB2_CENTRE_RE;
	if (bulb * bulb + c.im * c.im <= BULB2_RADIUS_SQ) {
		const Complex root = cx_sqrt(cx_make(-3.0 - 4.0 * c.re, -4.0 * c.im));
		return interior_corner(c, cx_make((root.re - 1.0) * 0.5, root.im * 0.5), BULB2_PERIOD, spacing, steps);
	}
	const Complex root = cx_sqrt(cx_make(1.0 - 4.0 * c.re, -4.0 * c.im));
	return interior_corner(c, cx_make((1.0 - root.re) * 0.5, -root.im * 0.5), CARDIOID_PERIOD, spacing, steps);
}

// --- At depth: the cycle in perturbation form -------------------------------------------------------
// c = C + delta_c with C the reference's point, which no double holds; the orbit is z = Z_m + delta,
// followed as the perturbation kernels follow it, rebased where |z| < |delta| or the reference ends.
// A deep cycle's points are a pixel's width apart at most in delta, not in z, so it is compared in
// delta's terms: two states are the same point when their difference — the references' difference
// plus the deltas' — is within DEEP_PERIOD_TOLERANCE_SQ of the caught state's |delta|. The kernel
// caught the state through BLA blocks, whose landings may sit a relative 3e-10 off the orbit's own
// cycle (design.md §9.4); a return within 2^-24 of |delta| is one, and a near-return of the
// reference orbit itself (a shallower minibrot's period) lies orders of magnitude further off.
#define DEEP_PERIOD_TOLERANCE_SQ 0x1p-48

typedef struct {
	int32_t m;
	Complex delta;
} DeepState;

static Complex reference_at(const KernelReference* r, int32_t m) {
	return cx_make(r->re[m], r->im[m]);
}

// One step, the kernel's arithmetic: delta <- (2 Z_m + delta) delta + delta_c; then z = Z_{m+1} +
// delta, rebased when |z| < |delta| or the reference ends.
static void deep_step(const KernelReference* r, Complex dc, DeepState* s) {
	s->delta = cx_add(cx_mul(cx_add(cx_twice(reference_at(r, s->m)), s->delta), s->delta), dc);
	s->m++;
	const Complex z = cx_add(reference_at(r, s->m), s->delta);
	if (cx_norm(z) < cx_norm(s->delta) || s->m == r->length) {
		s->delta = z;
		s->m = 0;
	}
}

static int64_t deep_least_period(const KernelReference* r, Complex dc, DeepState s0, int64_t distance,
                                 int64_t* steps) {
	const double tolerance = DEEP_PERIOD_TOLERANCE_SQ * cx_norm(s0.delta);
	DeepState s = s0;
	for (int64_t q = 1; q <= distance; q++) {
		deep_step(r, dc, &s);
		const Complex apart = cx_add(cx_sub(reference_at(r, s.m), reference_at(r, s0.m)), cx_sub(s.delta, s0.delta));
		if (cx_norm(apart) <= tolerance) {
			*steps += q;
			return q;
		}
	}
	*steps += distance;
	return 0;
}

// An orbit still bounded at the limit, against a NUCLEUS reference of period P, most likely lies in
// that minibrot's own cardioid, whose cycle has period P: its cycle is found by Newton's method on
// f^P(z) = z from where the orbit stopped (it has been falling toward the cycle all along), each
// step the orbit's P steps from z and their derivative, the correction applied to delta. Newton
// stops within DEEP_PERIOD_TOLERANCE_SQ of |delta|, or fails: then there is no certificate. A cycle
// it finds is a proof whatever made it look: an attracting P-cycle puts c inside the set, and the
// disc about it is interior_disc's. (An orbit in a satellite of period kP has no attracting
// P-cycle: Newton finds none, or one whose multiplier is at least 1, and certifies nothing.)
#define DEEP_NEWTON_STEPS 8

static bool deep_newton(const KernelReference* r, Complex dc, DeepState* s, int32_t period, int64_t* steps) {
	for (int step = 0; step < DEEP_NEWTON_STEPS; step++) {
		*steps += period;
		DeepState end = *s;
		Complex dz = cx_make(1.0, 0.0);
		for (int32_t i = 0; i < period; i++) {
			dz = cx_twice(cx_mul(cx_add(reference_at(r, end.m), end.delta), dz));
			deep_step(r, dc, &end);
		}
		const Complex apart = cx_add(cx_sub(reference_at(r, end.m), reference_at(r, s->m)), cx_sub(end.delta, s->delta));
		const double size = cx_norm(apart);
		if (!(size == size)) return false;
		if (size <= DEEP_PERIOD_TOLERANCE_SQ * cx_norm(s->delta)) return true;
		s->delta = cx_sub(s->delta, cx_div(apart, cx_sub(dz, cx_make(1.0, 0.0))));
	}
	return false;
}

// |a| without squaring it: at depth the derivatives below reach 2^900, whose square no double holds.
static double cx_abs(Complex a) {
	const double x = a.re < 0.0 ? -a.re : a.re, y = a.im < 0.0 ? -a.im : a.im;
	const double m = x > y ? x : y;
	if (!(m > 0.0)) return m;  // zero, or a NaN, which fails every certificate
	const double u = x / m, v = y / m;
	return m * sqrt(u * u + v * v);
}

// interior_disc at depth, with the derivatives in c taken per PIXEL (dc and dcdz times the spacing),
// so the radius comes out in pixels directly. Per unit of c, dz/dc reaches 1/size and d2z/dz2 as
// much, and their product — about 2^1700 at 10^-254 — holds in no double: the first form certified
// no cell deeper than about 10^-150 (measured at p4128: 3 159 interior discs, none a pixel wide).
static Corner deep_interior_disc(Complex dz, Complex dc_px, Complex dzdz, Complex dcdz_px) {
	const double multiplier = cx_norm(dz);
	if (!(multiplier < 1.0)) return no_corner();
	const Complex denominator = cx_add(dcdz_px, cx_div(cx_mul(dzdz, dc_px), cx_sub(cx_make(1.0, 0.0), dz)));
	Corner k = no_corner();
	k.kind = CORNER_INTERIOR;
	k.radius = (float)(QUARTER * (1.0 - multiplier) / cx_abs(denominator));
	return k;
}

// interior_corner's derivatives along the cycle from s0, the orbit stepped in perturbation form and
// the derivatives in c per pixel (deep_interior_disc).
static Corner deep_interior_corner(const KernelReference* r, Complex dc, DeepState s, int64_t period, double spacing,
                                   int64_t* steps) {
	*steps += period;
	Complex dz = cx_make(1.0, 0.0), dcv = cx_make(0.0, 0.0);
	Complex dzdz = cx_make(0.0, 0.0), dcdz = cx_make(0.0, 0.0);
	for (int64_t i = 0; i < period; i++) {
		const Complex z = cx_add(reference_at(r, s.m), s.delta);
		dcdz = cx_twice(cx_add(cx_mul(z, dcdz), cx_mul(dz, dcv)));
		dcv = cx_add(cx_twice(cx_mul(z, dcv)), cx_make(spacing, 0.0));
		dzdz = cx_twice(cx_add(cx_mul(z, dzdz), cx_mul(dz, dz)));
		dz = cx_twice(cx_mul(z, dz));
		deep_step(r, dc, &s);
	}
	return deep_interior_disc(dz, dcv, dzdz, dcdz);
}

Corner certify_corner(double c_re, double c_im, const KernelOrbit* o, double mu, double spacing,
                      const KernelReference* reference, bool at_limit, int64_t* steps) {
	const Complex c = cx_make(c_re, c_im);
	if (reference) {
		if (o->fate == KERNEL_FATE_ESCAPED) return exterior_corner(o, mu, spacing);
		DeepState s0;
		s0.m = o->m;
		s0.delta = cx_make(o->delta_re, o->delta_im);
		if (o->fate == KERNEL_FATE_INTERIOR_CYCLE) {
			const int64_t period = deep_least_period(reference, c, s0, o->cycle_distance, steps);
			return period > 0 ? deep_interior_corner(reference, c, s0, period, spacing, steps) : no_corner();
		}
		if (!at_limit || o->fate != KERNEL_FATE_INTERIOR_LIMIT || reference->period <= 0) return no_corner();
		if (!deep_newton(reference, c, &s0, reference->period, steps)) return no_corner();
		return deep_interior_corner(reference, c, s0, reference->period, spacing, steps);
	}
	switch (o->fate) {
	case KERNEL_FATE_ESCAPED:
		return exterior_corner(o, mu, spacing);
	case KERNEL_FATE_INTERIOR_SHAPE:
		return shape_corner(c, spacing, steps);
	case KERNEL_FATE_INTERIOR_CYCLE: {
			const Complex z0 = cx_make(o->z_re, o->z_im);
			const int64_t period = least_period(c, z0, o->cycle_distance, steps);
			return period > 0 ? interior_corner(c, z0, period, spacing, steps) : no_corner();
		}
	case KERNEL_FATE_INTERIOR_LIMIT:
	default:
		return no_corner();
	}
}

// --- Cells ---------------------------------------------------------------------------------------
// The cubic Hermite basis on [0, 1]: value at 0, value at 1, slope at 0, slope at 1.
typedef struct {
	double v0, v1, s0, s1;
} Hermite;

static Hermite hermite(double t) {
	const double u = 1.0 - t;
	Hermite h;
	h.v0 = (1.0 + 2.0 * t) * u * u;
	h.v1 = t * t * (3.0 - 2.0 * t);
	h.s0 = t * u * u;
	h.s1 = -t * t * u;
	return h;
}

// Bicubic Hermite with zero twist from four corners: continuous with its neighbours across every
// edge (an edge's values and cross-slopes depend on its own two corners alone).
static double interpolate(const Corner* k00, const Corner* k10, const Corner* k01, const Corner* k11,
                          const Hermite* hu, const Hermite* hv, double side) {
	const double values = k00->mu * hu->v0 * hv->v0 + k10->mu * hu->v1 * hv->v0 + k01->mu * hu->v0 * hv->v1 +
	                      k11->mu * hu->v1 * hv->v1;
	const double along = k00->gx * hu->s0 * hv->v0 + k10->gx * hu->s1 * hv->v0 + k01->gx * hu->s0 * hv->v1 +
	                     k11->gx * hu->s1 * hv->v1;
	const double down = k00->gy * hu->v0 * hv->s0 + k10->gy * hu->v1 * hv->s0 + k01->gy * hu->v0 * hv->s1 +
	                    k11->gy * hu->v1 * hv->s1;
	return values + side * (along + down);
}

static bool corner_current(const CertifyContext* c, const Corner* k, CornerKind kind) {
	return k->run == c->run && k->kind == kind;
}

// The kind of certificate a cell earns: all four corners of one kind, each disc wide enough.
// `ceiling` is the highest mu an exterior corner may have (see ESCAPE_HEADROOM).
static CornerKind judge_cell(const CertifyContext* c, const Corner* k[4], int stride, double ceiling) {
	const CornerKind kind = (CornerKind)k[0]->kind;
	const bool wanted = (kind == CORNER_EXTERIOR && c->exterior) || (kind == CORNER_INTERIOR && c->interior);
	if (!wanted) return CORNER_NONE;
	const double margin = kind == CORNER_EXTERIOR ? c->exterior_margin : c->interior_margin;
	const double needed = HALF_DIAGONAL * margin * (double)stride;
	for (int i = 0; i < 4; i++) {
		if (!corner_current(c, k[i], kind)) return CORNER_NONE;
		if (!((double)k[i]->radius >= needed)) return CORNER_NONE;
		if (kind == CORNER_EXTERIOR && !(k[i]->mu < ceiling)) return CORNER_NONE;
	}
	return kind;
}

#define CELL_MAX_STRIDE (CERT_LATTICE_STRIDE << (CERT_LEVELS - 1))

int certify_cells(const CertifyContext* c, PixelStore* store, int tile_x, int tile_y, int stride) {
	Hermite basis[CELL_MAX_STRIDE];
	for (int t = 0; t < stride; t++) basis[t] = hermite((double)t / (double)stride);
	const uint8_t interior_state = state_make(RANK_CERTIFIED, FATE_CERTIFIED);
	const uint8_t exterior_state = state_make(RANK_CERTIFIED, FATE_ESCAPED);
	const int32_t limit = (int32_t)c->max_iterations;
	const double m = c->exterior_margin;
	const double ceiling = (double)c->max_iterations - ESCAPE_HEADROOM - num_log2((m + 1.0) / (m - 1.0));
	int filled = 0;
	for (int y = tile_y; y < tile_y + MB_TILE_HEIGHT; y += stride) {
		if (y + stride >= MB_SCREEN_HEIGHT) break;  // the screen's last cells have no far corners
		for (int x = tile_x; x < tile_x + MB_TILE_WIDTH; x += stride) {
			if (x + stride >= MB_SCREEN_WIDTH) break;
			const Corner* k[4] = {
				&c->corners[certify_corner_index(x, y)],
				&c->corners[certify_corner_index(x + stride, y)],
				&c->corners[certify_corner_index(x, y + stride)],
				&c->corners[certify_corner_index(x + stride, y + stride)],
			};
			const CornerKind kind = judge_cell(c, k, stride, ceiling);
			if (kind == CORNER_NONE) continue;
			for (int v = 0; v < stride; v++) {
				const size_t row = (size_t)(y + v) * MB_SCREEN_WIDTH;
				for (int u = 0; u < stride; u++) {
					const size_t p = row + (size_t)(x + u);
					if (!state_needs_work(store->state[p])) continue;
					store->distance[p] = 0;  // decided at its own coordinate
					if (kind == CORNER_INTERIOR) {
						store->n[p] = limit;
						store->phi[p] = 0.0f;
						store->state[p] = interior_state;
					} else {
						double mu = interpolate(k[0], k[1], k[2], k[3], &basis[u], &basis[v], (double)stride);
						if (!(mu > 0.0)) mu = 0.0;
						const int32_t n = (int32_t)mu + 1;
						store->n[p] = n;
						store->phi[p] = (float)((double)n - mu);
						store->state[p] = exterior_state;
					}
					filled++;
				}
			}
		}
	}
	return filled;
}
