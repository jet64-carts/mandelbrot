// carts/mandelbrot/kernel.h
//
// THE LANE KERNEL — escape-time iteration of z <- z^2 + c over a stream of sample points.
//
// The kernel keeps KERNEL_SLOTS small GROUPS of samples in flight at once, KERNEL_GROUP lanes
// each, and steps all of them together, so every lane loop is a straight run of independent
// SIMD arithmetic that the compiler's vectoriser maps onto NEON (and wasm simd128): plain ISO C,
// no intrinsics, no vector extensions. A group is one 128-bit register of doubles, the SIMD
// width every Jet64 runtime shares; the slots are independent chains that hide each other's
// floating-point latency.
//
// Its speed comes from speculation. A block of KERNEL_BLOCK iterations runs with no escape test
// at all; a group whose lanes are all still inside the bailout circle at the end of the block
// commits it, and a group with an escaped lane restores the state saved before the block and
// replays it checked, freezing each lane where it escapes. A finished lane is PARKED on the
// fixed point z = 0, c = 0 and a finished group is REFILLED with the next samples of the stream,
// so no slot idles while work remains. The groups are independent: one group's escape never
// costs another group a replay.
//
// The answer is bit for bit the plain scalar loop's (tests/mandelbrot/oracle.c):
//
//     n = 0; z = 0
//     while n < max_iterations and |z|^2 <= bailout^2:  z = z^2 + c; n = n + 1
//
// with n = max_iterations meaning interior. Speculation changes how fast that answer is found,
// never what it is. Two interior shortcuts ride along, both exact with respect to that loop:
//   * the SHAPE test: the main cardioid and the period-2 bulb, in closed form;
//   * the CYCLE test: Brent's periodicity detection at block granularity. A lane whose state
//     returns bit for bit to an earlier state is periodic IN THE ARITHMETIC, so the scalar loop
//     could never escape either — which is why its tolerance is exactly zero (design.md §4.4).
//
// AT DEPTH (Phase 3) the same lanes iterate by PERTURBATION against a reference orbit Z_m computed
// in the bignum (reference.h): each lane carries its offset delta from the reference and its own
// index m into it, with z = Z_m + delta, and
//
//     delta <- (2 Z_m + delta) delta + delta_c,   m <- m + 1
//
// REBASED (Zhuoran) whenever |Z_m + delta| < |delta| or m reaches the reference's end: delta <- z,
// m <- 0. Rebasing is what makes one reference serve every pixel without glitches: z never cancels
// catastrophically against the reference, so it keeps the relative precision a double gives. Where
// the reference's BLA table allows, a lane skips whole blocks of steps linearly at a block boundary
// (delta <- A delta + B delta_c, design.md §9.3). The escape test, the speculation and the cycle
// test are the same machinery; the cycle test compares the whole state (delta, m), exactly.
#ifndef MANDELBROT_KERNEL_H
#define MANDELBROT_KERNEL_H

#include "jet64.h"

// The shape of the lockstep, chosen by measurement (design.md §4.2). The guards exist for the
// harness's ablation lane alone, which rebuilds the kernel at other shapes to re-measure the
// choice; the cart itself never overrides them.
#ifndef KERNEL_GROUP
#define KERNEL_GROUP 2  // lanes per group: one 128-bit vector of doubles
#endif
#ifndef KERNEL_SLOTS
#define KERNEL_SLOTS 6  // groups in flight
#endif
#ifndef KERNEL_BLOCK
#define KERNEL_BLOCK 8  // iterations per speculative block
#endif

// The derivative variant carries four doubles a lane where the plain one carries two, and at six
// slots its state no longer fits the register file: measured on one core (kernel_test --bench,
// the seahorse), 0.71 G lane-iterations a second at six slots against 1.47 at four and 1.37 at
// three. So it keeps four groups in flight (design.md §4.5).
#ifndef KERNEL_DERIVATIVE_SLOTS
#define KERNEL_DERIVATIVE_SLOTS 4
#endif

// The perturbation variants' groups in flight: their lanes carry delta, z and delta_c, and the
// derivative's dz too (measured by tests/mandelbrot/deep_test.c --bench, design.md §9.2).
#ifndef KERNEL_PERTURBED_SLOTS
#define KERNEL_PERTURBED_SLOTS 4
#endif
#ifndef KERNEL_PERTURBED_DERIVATIVE_SLOTS
#define KERNEL_PERTURBED_DERIVATIVE_SLOTS 3
#endif

// What a lane-iteration of the derivative variant costs, in percent of a plain one's: the plain
// kernel's 2.32 G lane-iterations a second on the seahorse against the derivative's 1.47 (the
// overview, 2.04 against 1.16). The scheduler budgets in plain lane-iterations, so it counts each
// derivative one at this price.
#define KERNEL_DERIVATIVE_COST_PERCENT 160
#define KERNEL_PERCENT 100

// What a lane-iteration of the perturbation variants costs, in percent of a plain one's, and what
// one BLA step costs, likewise (tests/mandelbrot/deep_test.c --bench, design.md §9.2): the plain
// kernel's 2.17 G lane-iterations a second on the seahorse against perturbation's 0.45 in lockstep
// on p711, p2126 and p3743 (4.75 to 4.85 times), and a BLA step at 22 to 24.5 plain
// lane-iterations — one lane's own scalar walk down the table's levels.
#define KERNEL_PERTURBED_COST_PERCENT 480
#define KERNEL_BLA_STEP_COST_PERCENT 2300

// The deepest level of a BLA table: level l collapses 2^l steps (reference.h sizes the table).
#define KERNEL_BLA_LEVELS 18

// A BLA step holds only where |Z_m| is small enough that the orbit cannot escape within it: beyond
// 2 it escapes in a few steps, and those steps are always taken exactly. The kernel's single steps
// and the table's levels (reference.c) apply the same rule.
#define BLA_Z_LIMIT_SQ 4.0

// How one sample's orbit ended.
typedef enum {
	KERNEL_FATE_ESCAPED = 0,     // |z_n|^2 exceeded the bailout at iteration n
	KERNEL_FATE_INTERIOR_SHAPE,  // inside the main cardioid or the period-2 bulb
	KERNEL_FATE_INTERIOR_CYCLE,  // the orbit returned exactly to an earlier state
	KERNEL_FATE_INTERIOR_LIMIT,  // still bounded after the iteration limit
} KernelFate;

// The rules every orbit in one call is iterated under.
typedef struct {
	double bailout_sq;       // |z|^2 beyond which an orbit has escaped
	int64_t max_iterations;  // the iteration limit; an orbit still bounded here is interior
	bool shape_test;         // test the cardioid and the period-2 bulb before iterating
	bool cycle_test;         // detect exactly periodic orbits while iterating
} KernelLimits;

// One sample point, c = re + i im.
typedef struct {
	double re, im;
} KernelPoint;

// One sample's fate.
typedef struct {
	int64_t iterations;  // n: the index of the first z beyond the bailout (the limit if interior)
	double radius_sq;    // |z_n|^2 for an escaped orbit, else 0
	KernelFate fate;

	// For an escaped orbit: z_n, and — from kernel_iterate_derivative only — its derivative with
	// respect to c, dz_n/dc, which the exterior certificate is built from (certify.h).
	// For a CYCLE orbit: z is the state the cycle was caught at, which lies on the orbit's
	// floating-point cycle, and cycle_distance the number of iterations it took to come back to it
	// (a multiple of the cycle's period), which the interior certificate starts from. At depth z
	// alone cannot name the state to a pixel's precision, so the perturbation variants also give it
	// as perturbation holds it: z = Z_m + delta, against the reference they iterated.
	double z_re, z_im;
	double dz_re, dz_im;
	int64_t cycle_distance;
	int32_t m;
	double delta_re, delta_im;
} KernelOrbit;

// One entry of a bivariate linear approximation: 2^l steps of the perturbation collapsed to
// delta <- A delta + B delta_c, valid while |delta|^2 < r2 (design.md §9.3 derives the radius).
typedef struct {
	double ar, ai;  // A
	double br, bi;  // B
	double r2;      // the validity radius, squared
} BlaEntry;

// A reference orbit as the perturbation kernels read it: Z_0 = 0 .. Z_length as doubles, and its
// BLA table. A lane at m = length rebases, so the orbit may end anywhere: where it escaped, at the
// period of a nucleus (Z_p ~ 0), or at the iteration limit.
typedef struct {
	const double* re;
	const double* im;
	int32_t length;
	// The table: bla[l] for l in 1..levels, bla[l][j] covering steps j 2^l .. (j + 1) 2^l of the
	// first bla_length steps; levels = 0 means no table. Level 0 (single steps) is computed on the
	// fly from Z_m with the same epsilon the table was built with.
	const BlaEntry* bla[KERNEL_BLA_LEVELS + 1];
	int levels;
	int32_t bla_length;
	double epsilon;
	int32_t period;  // a nucleus's period (its orbit repeats), or 0 for a centre reference
} KernelReference;

/**
 * Iterates `count` points to their fates under `limits`, writing orbits[i] for points[i].
 * Returns the lane-iterations executed — every lane of every stepped group, speculation and
 * parked lanes included — which is the kernel's true cost and the unit the scheduler budgets
 * in. A pure function of its arguments: the same points always cost the same.
 */
int64_t kernel_iterate(const KernelPoint* points, KernelOrbit* orbits, int count,
                       const KernelLimits* limits);

/**
 * kernel_iterate, carrying the derivative dz/dc of every orbit beside it (dz <- 2 z dz + 1).
 * The derivative never feeds back into z, so every fate, n and |z_n|^2 is bit for bit
 * kernel_iterate's; it costs about twice as much, and the scheduler asks for it only for the
 * coarse samples that become certificate corners.
 */
int64_t kernel_iterate_derivative(const KernelPoint* points, KernelOrbit* orbits, int count,
                                  const KernelLimits* limits);

// What a kernel did: the lane-iterations it executed (every lane of a group in flight, parked or
// not, as the plain variants count them) and the BLA steps it took (each one lane's alone). The two
// cost very differently, so the caller prices them apart (KERNEL_PERTURBED_COST_PERCENT,
// KERNEL_BLA_STEP_COST_PERCENT).
typedef struct {
	int64_t lanes;
	int64_t bla_steps;
} KernelCount;

/**
 * The perturbation variant: points[i] is delta_c, the sample's offset from the reference's c, and
 * each orbit's z is z_n itself (Z_m + delta). The shape test does not apply (a deep view is nowhere
 * near the cardioid's closed form); the cycle test compares the state (delta, m) exactly, and a
 * cycle's distance counts the iterations BLA took as well as the steps.
 */
KernelCount kernel_iterate_perturbed(const KernelPoint* points, KernelOrbit* orbits, int count,
                                     const KernelLimits* limits, const KernelReference* reference);

/** kernel_iterate_perturbed carrying dz/dc, through BLA steps as dz <- A dz + B. */
KernelCount kernel_iterate_perturbed_derivative(const KernelPoint* points, KernelOrbit* orbits, int count,
                                                const KernelLimits* limits, const KernelReference* reference);

#endif // MANDELBROT_KERNEL_H
