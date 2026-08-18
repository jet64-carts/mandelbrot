// carts/mandelbrot/numerics.c
//
// Each function reduces its argument exactly with bit manipulation, evaluates a fixed series on
// the reduced range by Horner's rule, and scales back exactly. The series are the textbook ones
// (the inverse hyperbolic tangent for the logarithm, the Taylor series for the exponential),
// truncated where the next term falls below half an ulp; their coefficients are written as
// constant expressions the compiler folds with one correctly rounded division each, so the table
// is the same on every toolchain.
#pragma STDC FP_CONTRACT OFF

#include "numerics.h"

// --- The binary64 layout ---------------------------------------------------------------------
#define DOUBLE_MANTISSA_BITS 52
#define DOUBLE_EXPONENT_BIAS 1023
#define DOUBLE_MANTISSA_MASK ((UINT64_C(1) << DOUBLE_MANTISSA_BITS) - 1)
#define DOUBLE_SMALLEST_NORMAL 0x1p-1022
#define DOUBLE_EXPONENT_MIN (-1022)  // the least exponent of a normal double
#define DOUBLE_EXPONENT_MAX 1023     // the greatest

// --- Mathematical constants, to the precision of a double ------------------------------------
#define LN_2 0.69314718055994530942
#define SQRT_2 1.41421356237309504880
#define ONE_THIRD (1.0 / 3.0)

typedef union {
	double value;
	uint64_t bits;
} DoubleBits;

// --- log2 ------------------------------------------------------------------------------------
// With x = m 2^e and m in (1/sqrt 2, sqrt 2], log2 x = e + log(1 + f) / ln 2 with f = m - 1,
// which is exact (Sterbenz). With s = f / (2 + f) (so |s| < 0.172), the inverse hyperbolic
// tangent series gives log(1 + f) = 2s + 2s^3/3 + 2s^5/5 + ..., and since 2s = f - s f,
//
//     log(1 + f) = f - s (f - R),   R = 2s^2/3 + 2s^4/5 + ...
//
// the arrangement of fdlibm's log: the exact f carries the result, and the rounding in the
// series touches only the small correction. Nine terms of R leave a truncation error below
// 2 s^20 / 21 < 5e-17 of the result.
static const double LOG_SERIES[] = {
	2.0 / 3.0,  2.0 / 5.0,  2.0 / 7.0,  2.0 / 9.0,  2.0 / 11.0,
	2.0 / 13.0, 2.0 / 15.0, 2.0 / 17.0, 2.0 / 19.0,
};
#define LOG_TERMS ((int)(sizeof(LOG_SERIES) / sizeof(LOG_SERIES[0])))
#define INV_LN_2 1.44269504088896340736

double num_log2(double x) {
	if (!(x >= DOUBLE_SMALLEST_NORMAL)) x = DOUBLE_SMALLEST_NORMAL;  // NaN lands here too
	DoubleBits u;
	u.value = x;
	int exponent = (int)(u.bits >> DOUBLE_MANTISSA_BITS) - DOUBLE_EXPONENT_BIAS;
	u.bits = (u.bits & DOUBLE_MANTISSA_MASK) | ((uint64_t)DOUBLE_EXPONENT_BIAS << DOUBLE_MANTISSA_BITS);
	double m = u.value;  // in [1, 2)
	if (m > SQRT_2) {
		m *= 0.5;  // exact
		exponent += 1;
	}
	const double f = m - 1.0;
	const double s = f / (2.0 + f);
	const double z = s * s;
	double series = LOG_SERIES[LOG_TERMS - 1];
	for (int k = LOG_TERMS - 2; k >= 0; k--) series = series * z + LOG_SERIES[k];
	const double log_m = f - s * (f - z * series);
	return (double)exponent + log_m * INV_LN_2;
}

// --- exp2 ------------------------------------------------------------------------------------
// With x = n + f, n the nearest integer (so |f| <= 1/2, and f = x - n is exact),
// 2^x = 2^n e^(f ln 2), and the Taylor series of e^y for |y| <= ln 2 / 2 < 0.347 needs fourteen
// terms to fall below y^14 / 14! < 5e-18. The 2^n is built directly in the exponent field.
static const double EXP_SERIES[] = {
	1.0,
	1.0,
	1.0 / 2.0,
	1.0 / 6.0,
	1.0 / 24.0,
	1.0 / 120.0,
	1.0 / 720.0,
	1.0 / 5040.0,
	1.0 / 40320.0,
	1.0 / 362880.0,
	1.0 / 3628800.0,
	1.0 / 39916800.0,
	1.0 / 479001600.0,
	1.0 / 6227020800.0,
};
#define EXP_TERMS ((int)(sizeof(EXP_SERIES) / sizeof(EXP_SERIES[0])))
#define ROUND_HALF 0.5

double num_exp2(double x) {
	if (!(x >= DOUBLE_EXPONENT_MIN)) x = DOUBLE_EXPONENT_MIN;  // NaN lands here too
	if (x > DOUBLE_EXPONENT_MAX) x = DOUBLE_EXPONENT_MAX;
	const int n = (int)(x >= 0.0 ? x + ROUND_HALF : x - ROUND_HALF);  // nearest, ties away
	const double y = (x - (double)n) * LN_2;
	double series = EXP_SERIES[EXP_TERMS - 1];
	for (int k = EXP_TERMS - 2; k >= 0; k--) series = series * y + EXP_SERIES[k];
	// The clamp keeps n inside the normal exponents, and at either end f has the sign that keeps
	// the product normal too: n = -1022 only for f >= 0, n = 1023 only for f <= 0.
	DoubleBits scale;
	scale.bits = (uint64_t)(n + DOUBLE_EXPONENT_BIAS) << DOUBLE_MANTISSA_BITS;
	return series * scale.value;
}

// --- pow and cbrt ----------------------------------------------------------------------------
double num_pow(double x, double y) {
	return num_exp2(y * num_log2(x));
}

double num_cbrt(double x) {
	if (x == 0.0) return 0.0;
	const double a = x < 0.0 ? -x : x;
	// A first root from the logarithm, then one Newton step for r^3 = a, which squares the
	// relative error: r <- (2r + a / r^2) / 3.
	double r = num_exp2(num_log2(a) * ONE_THIRD);
	r = (r + r + a / (r * r)) * ONE_THIRD;
	return x < 0.0 ? -r : r;
}
