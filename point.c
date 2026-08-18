// carts/mandelbrot/point.c
//
// The DeepPoint as two full-precision bignums. Every verb is exact but the ones that round to a
// double, which round to nearest; so every runtime moves the centre to the same bits.
#pragma STDC FP_CONTRACT OFF

#include "point.h"

DeepPoint point_make(double re, double im) {
	DeepPoint p;
	p.re = big_from_double(re);
	p.im = big_from_double(im);
	return p;
}

bool point_parse(DeepPoint* p, const char* re, const char* im) {
	return big_parse(&p->re, re, BIG_FRAC_LIMBS) && big_parse(&p->im, im, BIG_FRAC_LIMBS);
}

void point_offset(DeepPoint* p, double d_re, double d_im) {
	big_add_double(&p->re, &p->re, d_re, BIG_FRAC_LIMBS);
	big_add_double(&p->im, &p->im, d_im, BIG_FRAC_LIMBS);
}

void point_difference(const DeepPoint* a, const DeepPoint* b, double* d_re, double* d_im) {
	*d_re = big_difference(&a->re, &b->re, BIG_FRAC_LIMBS);
	*d_im = big_difference(&a->im, &b->im, BIG_FRAC_LIMBS);
}

double point_re(const DeepPoint* p) {
	return big_to_double(&p->re);
}

double point_im(const DeepPoint* p) {
	return big_to_double(&p->im);
}

// One coordinate held to [-limit, limit], exactly (the limit is a double, and so a Big exactly).
static void clamp_coordinate(Big* x, double limit) {
	const Big high = big_from_double(limit), low = big_from_double(-limit);
	if (big_compare(x, &high) > 0) *x = high;
	if (big_compare(x, &low) < 0) *x = low;
}

void point_clamp(DeepPoint* p, double limit) {
	clamp_coordinate(&p->re, limit);
	clamp_coordinate(&p->im, limit);
}

bool point_equal(const DeepPoint* a, const DeepPoint* b) {
	return big_compare(&a->re, &b->re) == 0 && big_compare(&a->im, &b->im) == 0;
}

void point_coordinates(const DeepPoint* p, Big* re, Big* im) {
	*re = p->re;
	*im = p->im;
}

DeepPoint point_from_coordinates(const Big* re, const Big* im) {
	DeepPoint p;
	p.re = *re;
	p.im = *im;
	return p;
}
