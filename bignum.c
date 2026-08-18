// carts/mandelbrot/bignum.c
//
// A number at precision f occupies limbs [BIG_FRAC_LIMBS - f, BIG_LIMBS): its value is the
// two's-complement integer those limbs spell, times 2^(-32 f). Products work on magnitudes (the
// sign is applied after), so their truncation is toward zero; the double conversions work on bits,
// never on floating-point arithmetic, so they are exact wherever the text above says so.
#pragma STDC FP_CONTRACT OFF

#include "bignum.h"

#define LIMB_MASK 0xFFFFFFFFull
#define PRODUCT_LIMBS (2 * BIG_LIMBS)
#define DECIMAL_BASE 10u
#define INTEGER_LIMIT 0x80000000ull  // the integer limb holds |x| < 2^31

// IEEE 754 binary64.
#define DOUBLE_MANTISSA_BITS 52
#define DOUBLE_EXPONENT_BIAS 1023
#define DOUBLE_EXPONENT_MASK 0x7FFu
#define DOUBLE_SIGN_BIT 63
#define DOUBLE_HIDDEN_BIT (1ull << DOUBLE_MANTISSA_BITS)
#define DOUBLE_MANTISSA_MASK (DOUBLE_HIDDEN_BIT - 1)
#define DOUBLE_SUBNORMAL_EXPONENT (1 - DOUBLE_EXPONENT_BIAS)
#define WORD_BITS 64
#define DECIMAL_WHOLE_DIGITS 16  // |x| < 2^31 has at most ten
#define ROUND_SHIFT (WORD_BITS - 1 - DOUBLE_MANTISSA_BITS)  // 11 bits below a 53-bit mantissa
#define FRACTION_BITS (BIG_FRAC_LIMBS * BIG_LIMB_BITS)       // bit p of the limbs weighs 2^(p - this)

static inline int first_limb(int frac) {
	return BIG_FRAC_LIMBS - frac;
}

static uint64_t double_bits(double d) {
	uint64_t u;
	memcpy(&u, &d, sizeof u);
	return u;
}

static double bits_double(uint64_t u) {
	double d;
	memcpy(&d, &u, sizeof d);
	return d;
}

Big big_zero(void) {
	Big z;
	memset(&z, 0, sizeof z);
	return z;
}

bool big_negative(const Big* x) {
	return (int32_t)x->limb[BIG_INTEGER] < 0;
}

void big_truncate(Big* x, int frac) {
	for (int i = 0; i < first_limb(frac); i++) x->limb[i] = 0;
}

int big_limbs_for_bits(int bits) {
	const int limbs = (bits + BIG_LIMB_BITS - 1) / BIG_LIMB_BITS;
	return limbs < 1 ? 1 : limbs > BIG_FRAC_LIMBS ? BIG_FRAC_LIMBS : limbs;
}

void big_add(Big* r, const Big* a, const Big* b, int frac) {
	uint64_t carry = 0;
	for (int i = first_limb(frac); i < BIG_LIMBS; i++) {
		const uint64_t t = (uint64_t)a->limb[i] + b->limb[i] + carry;
		r->limb[i] = (uint32_t)t;
		carry = t >> BIG_LIMB_BITS;
	}
	big_truncate(r, frac);
}

void big_sub(Big* r, const Big* a, const Big* b, int frac) {
	uint64_t borrow = 0;
	for (int i = first_limb(frac); i < BIG_LIMBS; i++) {
		const uint64_t t = (uint64_t)a->limb[i] - b->limb[i] - borrow;
		r->limb[i] = (uint32_t)t;
		borrow = (t >> BIG_LIMB_BITS) != 0;
	}
	big_truncate(r, frac);
}

void big_neg(Big* r, const Big* a, int frac) {
	const Big zero = big_zero();
	big_sub(r, &zero, a, frac);
}

// |x| at precision f, into f + 1 limbs from the bottom; returns whether x was negative.
static bool magnitude(uint32_t* m, const Big* x, int frac) {
	const bool negative = big_negative(x);
	Big t;
	if (negative) {
		big_neg(&t, x, frac);
		x = &t;
	}
	for (int i = 0; i <= frac; i++) m[i] = x->limb[first_limb(frac) + i];
	return negative;
}

// The fixed-point result of a product of two (f + 1)-limb magnitudes: the limbs f .. 2f of the
// double-width product, signed.
static void place_product(Big* r, const uint32_t* p, bool negative, int frac) {
	for (int i = 0; i < first_limb(frac); i++) r->limb[i] = 0;
	for (int i = 0; i <= frac; i++) r->limb[first_limb(frac) + i] = p[frac + i];
	if (negative) big_neg(r, r, frac);
}

void big_mul(Big* r, const Big* a, const Big* b, int frac) {
	uint32_t ma[BIG_LIMBS], mb[BIG_LIMBS];
	uint32_t p[PRODUCT_LIMBS];
	const int n = frac + 1;
	const bool negative = magnitude(ma, a, frac) != magnitude(mb, b, frac);
	for (int i = 0; i < 2 * n; i++) p[i] = 0;
	for (int i = 0; i < n; i++) {
		uint64_t carry = 0;
		for (int j = 0; j < n; j++) {
			const uint64_t t = (uint64_t)ma[i] * mb[j] + p[i + j] + carry;
			p[i + j] = (uint32_t)t;
			carry = t >> BIG_LIMB_BITS;
		}
		p[i + n] = (uint32_t)carry;
	}
	place_product(r, p, negative, frac);
}

void big_square(Big* r, const Big* a, int frac) {
	uint32_t m[BIG_LIMBS];
	uint32_t p[PRODUCT_LIMBS];
	const int n = frac + 1;
	magnitude(m, a, frac);
	for (int i = 0; i < 2 * n; i++) p[i] = 0;
	// Each cross product a_i a_j (i < j) once ...
	for (int i = 0; i < n; i++) {
		uint64_t carry = 0;
		for (int j = i + 1; j < n; j++) {
			const uint64_t t = (uint64_t)m[i] * m[j] + p[i + j] + carry;
			p[i + j] = (uint32_t)t;
			carry = t >> BIG_LIMB_BITS;
		}
		p[i + n] = (uint32_t)carry;
	}
	// ... doubled ...
	uint32_t top = 0;
	for (int i = 0; i < 2 * n; i++) {
		const uint32_t next = p[i] >> (BIG_LIMB_BITS - 1);
		p[i] = (p[i] << 1) | top;
		top = next;
	}
	// ... and the squares on the diagonal added.
	uint64_t carry = 0;
	for (int i = 0; i < n; i++) {
		const uint64_t sq = (uint64_t)m[i] * m[i];
		uint64_t t = (uint64_t)p[2 * i] + (sq & LIMB_MASK) + carry;
		p[2 * i] = (uint32_t)t;
		t = (uint64_t)p[2 * i + 1] + (sq >> BIG_LIMB_BITS) + (t >> BIG_LIMB_BITS);
		p[2 * i + 1] = (uint32_t)t;
		carry = t >> BIG_LIMB_BITS;
	}
	place_product(r, p, false, frac);
}

void big_shift(Big* r, const Big* a, int bits, int frac) {
	const int lo = first_limb(frac);
	const int n = frac + 1;
	const uint32_t fill = big_negative(a) ? 0xFFFFFFFFu : 0u;  // what an arithmetic shift brings in
	uint32_t w[BIG_LIMBS];
	for (int i = 0; i < n; i++) w[i] = a->limb[lo + i];
	const int limbs = (bits < 0 ? -bits : bits) / BIG_LIMB_BITS;
	const int rest = (bits < 0 ? -bits : bits) % BIG_LIMB_BITS;
	uint32_t out[BIG_LIMBS];
	for (int i = 0; i < n; i++) {
		// The limb that lands at i, and its neighbour whose bits spill across.
		if (bits >= 0) {
			const int src = i - limbs;
			const uint32_t hi = src >= 0 ? w[src] : 0u;
			const uint32_t below = src - 1 >= 0 ? w[src - 1] : 0u;
			out[i] = rest ? (hi << rest) | (below >> (BIG_LIMB_BITS - rest)) : hi;
		} else {
			const int src = i + limbs;
			const uint32_t low = src < n ? w[src] : fill;
			const uint32_t above = src + 1 < n ? w[src + 1] : fill;
			out[i] = rest ? (low >> rest) | (above << (BIG_LIMB_BITS - rest)) : low;
		}
	}
	for (int i = 0; i < lo; i++) r->limb[i] = 0;
	for (int i = 0; i < n; i++) r->limb[lo + i] = out[i];
}

int big_compare(const Big* a, const Big* b) {
	const int32_t ia = (int32_t)a->limb[BIG_INTEGER], ib = (int32_t)b->limb[BIG_INTEGER];
	if (ia != ib) return ia < ib ? -1 : 1;
	for (int i = BIG_INTEGER - 1; i >= 0; i--) {
		if (a->limb[i] != b->limb[i]) return a->limb[i] < b->limb[i] ? -1 : 1;
	}
	return 0;
}

// --- Doubles ---------------------------------------------------------------------------------------
// The number of zero bits above the leading one of a nonzero limb.
static int leading_zeros(uint32_t w) {
	int n = 0;
	while (!(w & (1u << (BIG_LIMB_BITS - 1)))) {
		w <<= 1;
		n++;
	}
	return n;
}

static Big from_double_at(double d, int frac) {
	Big r = big_zero();
	const uint64_t u = double_bits(d);
	const unsigned field = (unsigned)(u >> DOUBLE_MANTISSA_BITS) & DOUBLE_EXPONENT_MASK;
	if (field == DOUBLE_EXPONENT_MASK) return r;  // inf or NaN: not a number the view can hold
	const uint64_t mantissa = (u & DOUBLE_MANTISSA_MASK) | (field ? DOUBLE_HIDDEN_BIT : 0);
	const int exponent = (field ? (int)field : 1) - DOUBLE_EXPONENT_BIAS - DOUBLE_MANTISSA_BITS;
	// |d| = mantissa 2^exponent: its bit b sits at limb bit b + exponent + FRACTION_BITS.
	const int lowest = first_limb(frac) * BIG_LIMB_BITS;
	for (int b = 0; b <= DOUBLE_MANTISSA_BITS; b++) {
		if (!((mantissa >> b) & 1u)) continue;
		const int p = b + exponent + FRACTION_BITS;
		if (p < lowest) continue;  // below the precision: truncated
		if (p >= BIG_LIMBS * BIG_LIMB_BITS - 1) return big_zero();  // beyond the integer limb
		r.limb[p / BIG_LIMB_BITS] |= 1u << (p % BIG_LIMB_BITS);
	}
	if (u >> DOUBLE_SIGN_BIT) big_neg(&r, &r, frac);
	return r;
}

Big big_from_double(double d) {
	return from_double_at(d, BIG_FRAC_LIMBS);
}

double big_to_double(const Big* x) {
	Big m = *x;
	const bool negative = big_negative(x);
	if (negative) big_neg(&m, x, BIG_FRAC_LIMBS);
	int top = BIG_INTEGER;
	while (top >= 0 && m.limb[top] == 0) top--;
	if (top < 0) return 0.0;
	// The 64 bits from the leading one down, from the top three limbs, and whether any bit below
	// them is set.
	const int zeros = leading_zeros(m.limb[top]);
	const uint64_t high = ((uint64_t)m.limb[top] << BIG_LIMB_BITS) | (top >= 1 ? m.limb[top - 1] : 0u);
	const uint32_t next = top >= 2 ? m.limb[top - 2] : 0u;
	const uint64_t word = (high << zeros) | (zeros ? (uint64_t)next >> (BIG_LIMB_BITS - zeros) : 0u);
	bool sticky = (uint32_t)(next << zeros) != 0;
	for (int i = 0; i < top - 2 && !sticky; i++) sticky = m.limb[i] != 0;
	int exponent = top * BIG_LIMB_BITS + (BIG_LIMB_BITS - 1 - zeros) - FRACTION_BITS;  // x = 1.m 2^exponent

	// Round the 64 bits to 53, to nearest, ties to even.
	uint64_t mantissa = word >> ROUND_SHIFT;
	const uint64_t rest = word & ((1ull << ROUND_SHIFT) - 1);
	const uint64_t half = 1ull << (ROUND_SHIFT - 1);
	if (rest > half || (rest == half && (sticky || (mantissa & 1u)))) mantissa++;
	if (mantissa >> (DOUBLE_MANTISSA_BITS + 1)) {
		mantissa >>= 1;
		exponent++;
	}
	if (exponent < DOUBLE_SUBNORMAL_EXPONENT) return 0.0;
	const uint64_t field = (uint64_t)(exponent + DOUBLE_EXPONENT_BIAS);
	const uint64_t u = ((uint64_t)negative << DOUBLE_SIGN_BIT) | (field << DOUBLE_MANTISSA_BITS) |
	                   (mantissa & DOUBLE_MANTISSA_MASK);
	return bits_double(u);
}

void big_add_double(Big* r, const Big* x, double d, int frac) {
	const Big t = from_double_at(d, frac);
	big_add(r, x, &t, frac);
}

double big_difference(const Big* a, const Big* b, int frac) {
	Big t;
	big_sub(&t, a, b, frac);
	return big_to_double(&t);
}

// --- Decimals --------------------------------------------------------------------------------------
static bool is_digit(char c) {
	return c >= '0' && c <= '9';
}

bool big_parse(Big* out, const char* text, int frac) {
	*out = big_zero();
	const char* p = text;
	const bool negative = *p == '-';
	if (*p == '-' || *p == '+') p++;
	uint64_t whole = 0;
	int whole_digits = 0;
	while (is_digit(*p)) {
		whole = whole * DECIMAL_BASE + (uint64_t)(*p - '0');
		if (whole >= INTEGER_LIMIT) return false;
		p++;
		whole_digits++;
	}
	const char* fraction = p;
	int fraction_digits = 0;
	if (*p == '.') {
		fraction = ++p;
		while (is_digit(*p)) {
			p++;
			fraction_digits++;
		}
	}
	if (*p != '\0' || whole_digits + fraction_digits == 0) return false;
	// The fraction from its last digit: x <- (x + d) / 10 keeps x = floor(true x) exactly, since a
	// floored x lies less than one unit below the true one, and dividing by ten cannot carry that
	// shortfall across a unit.
	Big x = big_zero();
	const int lo = first_limb(frac);
	for (int k = fraction_digits - 1; k >= 0; k--) {
		x.limb[BIG_INTEGER] += (uint32_t)(fraction[k] - '0');
		uint64_t remainder = 0;
		for (int i = BIG_INTEGER; i >= lo; i--) {
			const uint64_t current = (remainder << BIG_LIMB_BITS) | x.limb[i];
			x.limb[i] = (uint32_t)(current / DECIMAL_BASE);
			remainder = current % DECIMAL_BASE;
		}
	}
	x.limb[BIG_INTEGER] = (uint32_t)whole;
	if (negative) big_neg(&x, &x, frac);
	*out = x;
	return true;
}

bool big_format(const Big* x, int digits, char* out, int capacity) {
	Big m = *x;
	const bool negative = big_negative(x);
	if (negative) big_neg(&m, x, BIG_FRAC_LIMBS);
	char whole[DECIMAL_WHOLE_DIGITS];
	int w = 0;
	uint32_t integer = m.limb[BIG_INTEGER];
	do {
		whole[w++] = (char)('0' + integer % DECIMAL_BASE);
		integer /= DECIMAL_BASE;
	} while (integer);
	const int needed = (negative ? 1 : 0) + w + (digits > 0 ? 1 + digits : 0) + 1;
	if (needed > capacity) return false;
	int at = 0;
	if (negative) out[at++] = '-';
	while (w > 0) out[at++] = whole[--w];
	if (digits > 0) out[at++] = '.';
	for (int d = 0; d < digits; d++) {
		// Times ten: the digit is what carries out of the top fraction limb.
		uint64_t carry = 0;
		for (int i = 0; i < BIG_INTEGER; i++) {
			const uint64_t t = (uint64_t)m.limb[i] * DECIMAL_BASE + carry;
			m.limb[i] = (uint32_t)t;
			carry = t >> BIG_LIMB_BITS;
		}
		out[at++] = (char)('0' + carry);
	}
	out[at] = '\0';
	return true;
}
