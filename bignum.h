// carts/mandelbrot/bignum.h
//
// THE FIXED-POINT BIGNUM — the precision the view's centre and the reference orbit need at depth.
//
// A Big is a two's-complement fixed-point number: one signed 32-bit INTEGER limb and up to
// BIG_FRAC_LIMBS 32-bit FRACTION limbs, least significant first. Only the arithmetic the explorer
// needs is here — add, subtract, negate, multiply, square, shift, compare, exact decimal parsing
// and printing, and exact conversions to and from doubles — in ISO C with 64-bit products of
// 32-bit limbs (no __int128), so TinyCC and wasm32 build it as they build the rest.
//
// PRECISION. Every operation takes the number of fraction limbs to work at, `frac`. The limbs are
// laid out from the TOP: the integer limb is always limb[BIG_FRAC_LIMBS], and a number at
// precision f uses the f fraction limbs beneath it; the limbs below those are zero. So numbers of
// different precisions share one layout, and lowering a number's precision is truncating it. The
// view's centre works at the full BIG_FRAC_LIMBS (a handful of operations a tick); a reference
// orbit at the fewest limbs its depth needs (millions of operations, reference.c).
//
// ROUNDING. Every result is EXACT but a product's, which is truncated toward zero (the magnitude's
// fraction beyond `frac` limbs is dropped); decimal parsing truncates toward zero too, exactly
// (the digit-by-digit division below never loses a bit it keeps); a conversion to a double rounds
// to nearest, ties to even. The integer limb holds |x| < 2^31: the reference iteration squares only
// values inside the bailout circle, whose squares stay below 2^22.
//
// Nothing allocates: a Big is a plain value of BIG_LIMBS limbs.
#ifndef MANDELBROT_BIGNUM_H
#define MANDELBROT_BIGNUM_H

#include "jet64.h"

// The fraction limbs a Big can hold: 30 x 32 = 960 bits, which set the depth limit (config.h: a pixel
// spacing of 2^-880). A pan there moves the centre by a double whose last bit lies about 62 bits
// below the spacing (a sub-pixel step of 2^-10 pixel carries 52 bits beneath its leading one), 942
// bits in all; a reference at the depth limit needs 880 bits and its guard (64, see reference.h),
// 944; the 16 left are margin. reference.h asserts both fit.
#define BIG_FRAC_LIMBS 30
#define BIG_LIMBS (BIG_FRAC_LIMBS + 1)
#define BIG_LIMB_BITS 32
#define BIG_INTEGER BIG_FRAC_LIMBS  // the index of the integer limb

typedef struct {
	uint32_t limb[BIG_LIMBS];
} Big;

/** Zero. */
Big big_zero(void);

/** The double d exactly, if its bits fit the full precision; otherwise truncated toward zero. */
Big big_from_double(double d);

/** x rounded to the nearest double, ties to even (0 when |x| is below the smallest normal). */
double big_to_double(const Big* x);

/**
 * Parses a decimal string, "[+-]digits[.digits]", truncating toward zero to `frac` fraction limbs,
 * exactly. Returns false (and zero) for anything else, or an integer part of 2^31 or more.
 */
bool big_parse(Big* out, const char* text, int frac);

/**
 * Writes x as a decimal string with `digits` digits after the point, truncated toward zero, into
 * `out` (at most `capacity` bytes with the terminator). Returns false if it does not fit. With
 * 32 x frac digits the expansion is exact, and big_parse gives back x itself.
 */
bool big_format(const Big* x, int digits, char* out, int capacity);

/** r = a + b at `frac` fraction limbs (exact; r may alias either). */
void big_add(Big* r, const Big* a, const Big* b, int frac);

/** r = a - b at `frac` fraction limbs (exact; r may alias either). */
void big_sub(Big* r, const Big* a, const Big* b, int frac);

/** r = -a (exact; r may alias a). */
void big_neg(Big* r, const Big* a, int frac);

/** r = a b at `frac` fraction limbs, truncated toward zero (r may alias either). */
void big_mul(Big* r, const Big* a, const Big* b, int frac);

/** r = a^2 at `frac` fraction limbs, truncated toward zero: the reference iteration's hot path,
 * which computes each cross product once. (r may alias a.) */
void big_square(Big* r, const Big* a, int frac);

/** r = a 2^bits (bits > 0 shifts left, < 0 right, truncating toward minus infinity; r may alias a). */
void big_shift(Big* r, const Big* a, int bits, int frac);

/** r = x + d, exactly when d's bits fit `frac` limbs. */
void big_add_double(Big* r, const Big* x, double d, int frac);

/** a - b rounded to the nearest double. */
double big_difference(const Big* a, const Big* b, int frac);

/** -1, 0 or 1 as a < b, a == b, a > b. */
int big_compare(const Big* a, const Big* b);

/** True when x < 0. */
bool big_negative(const Big* x);

/** x truncated to `frac` fraction limbs (the limbs below zeroed). */
void big_truncate(Big* x, int frac);

/** The fraction limbs needed to hold `bits` bits of fraction, at most BIG_FRAC_LIMBS. */
int big_limbs_for_bits(int bits);

#endif // MANDELBROT_BIGNUM_H
