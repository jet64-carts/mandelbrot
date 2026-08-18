// carts/mandelbrot/numerics.h
//
// IEEE-ONLY ELEMENTARY FUNCTIONS — log2, exp2, pow and cbrt from + - * / and bit manipulation.
//
// The standard does not require log, exp, pow or cbrt to be correctly rounded, and the three
// runtimes link three different libms (the host's, TinyCC's, and the web player's), so any of
// them in a path that reaches a pixel would let the same cart draw different pictures on the
// workstation, a desktop export and a browser. Addition, subtraction, multiplication and
// division ARE correctly rounded everywhere, so functions built only from them — under
// FP_CONTRACT OFF, so no compiler may fuse a multiply into an add — give the same bits on every
// runtime. These are those functions. Each is tested against the host libm over its whole
// domain by tests/mandelbrot/numerics_test.c, which asserts the error bound stated below; the
// bounds are measured against that libm, itself within an ulp of the truth.
#ifndef MANDELBROT_NUMERICS_H
#define MANDELBROT_NUMERICS_H

#include "jet64.h"

/**
 * log2(x) for finite x > 0. Arguments below the smallest normal double are treated as that
 * double. Error: within 2 ulp of the true value (asserted by the harness).
 */
double num_log2(double x);

/**
 * 2^x for any finite x; the result is clamped to the normal range, so x below -1022 gives
 * 2^-1022 and x above 1023 gives 2^1023. Error: within 2 ulp (asserted by the harness).
 */
double num_exp2(double x);

/**
 * x^y for finite x > 0, as 2^(y log2 x). Relative error grows with |y log2 x|; within 1e-14
 * over the colour pipeline's arguments (asserted by the harness).
 */
double num_pow(double x, double y);

/**
 * The real cube root of any finite x (negative x gives a negative root, zero gives zero).
 * Relative error within 1e-15 (asserted by the harness).
 */
double num_cbrt(double x);

/** v held to [lo, hi]. */
static inline double num_clamp(double v, double lo, double hi) {
	return v < lo ? lo : (v > hi ? hi : v);
}

#endif // MANDELBROT_NUMERICS_H
