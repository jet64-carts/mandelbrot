// carts/mandelbrot/point.h
//
// A POINT OF THE COMPLEX PLANE, held at whatever precision the view needs.
//
// The view's centre is the one number in the explorer that must outlive double precision: at a
// magnification of 2^200 a double centre could not even name the pixel under the cursor. So
// nothing outside this module looks inside a DeepPoint. Callers move a point by a small double
// offset, measure the small double offset between two points, and read a point's double
// approximation (which is all direct iteration and the view's limits need). Each coordinate is a
// fixed-point bignum at full precision (bignum.h): a pan adds its double offset EXACTLY, so no
// pan, however long or however deep, ever drifts the view off the grid it has computed (reuse.h),
// and the reference orbits at depth start from the exact point (reference.h). design.md §3.
#ifndef MANDELBROT_POINT_H
#define MANDELBROT_POINT_H

#include "bignum.h"

// The fields are this module's own; code elsewhere uses the functions below.
typedef struct {
	Big re, im;
} DeepPoint;

/** The point re + i im, exactly. */
DeepPoint point_make(double re, double im);

/** The point written in decimal, truncated to full precision; false if either text is malformed. */
bool point_parse(DeepPoint* p, const char* re, const char* im);

/** Moves `p` by the offset (d_re, d_im), exactly. */
void point_offset(DeepPoint* p, double d_re, double d_im);

/** The offset a - b, each coordinate rounded to the nearest double. */
void point_difference(const DeepPoint* a, const DeepPoint* b, double* d_re, double* d_im);

/** The nearest double to each coordinate. */
double point_re(const DeepPoint* p);
double point_im(const DeepPoint* p);

/** Holds each coordinate of `p` within [-limit, limit], exactly. */
void point_clamp(DeepPoint* p, double limit);

/** True when the two points are the same point. */
bool point_equal(const DeepPoint* a, const DeepPoint* b);

/** The coordinates themselves, for the reference orbit's arithmetic (reference.c). */
void point_coordinates(const DeepPoint* p, Big* re, Big* im);

/** The point with these coordinates. */
DeepPoint point_from_coordinates(const Big* re, const Big* im);

#endif // MANDELBROT_POINT_H
