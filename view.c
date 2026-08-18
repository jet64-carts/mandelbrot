// carts/mandelbrot/view.c
#pragma STDC FP_CONTRACT OFF

#include "view.h"

#include "numerics.h"

// log2 of the home view's pixel spacing, computed once by view_init.
static double s_home_scale_log2;

void view_init(void) {
	s_home_scale_log2 = num_log2(MB_HOME_HEIGHT / (double)MB_SCREEN_HEIGHT);
}

View view_home(void) {
	View v;
	v.centre = point_make(MB_HOME_CENTRE_RE, MB_HOME_CENTRE_IM);
	v.scale_log2 = s_home_scale_log2;
	return v;
}

double view_spacing(const View* v) {
	return num_exp2(v->scale_log2);
}

double view_octaves(const View* v) {
	return s_home_scale_log2 - v->scale_log2;
}

void view_zoom_about(View* v, double octaves, double anchor_x, double anchor_y) {
	// With a the anchor's offset from the screen's centre, the plane point under it is
	// centre + a h before and must be centre' + a h' after.
	const double a_re = view_offset_x(anchor_x);
	const double a_im = view_offset_y(anchor_y);
	const double before = view_spacing(v);
	v->scale_log2 -= octaves;
	const double shrink = before - view_spacing(v);
	point_offset(&v->centre, a_re * shrink, a_im * shrink);
}

void view_pan_pixels(View* v, double dx, double dy) {
	const double h = view_spacing(v);
	point_offset(&v->centre, -dx * h, dy * h);
}

double view_depth_limit(const View* v) {
	(void)v;  // the same everywhere: the double's exponent range does not depend on the centre
	return MB_DEPTH_LIMIT_LOG2;
}

double view_direct_limit(const View* v) {
	const double re = point_re(&v->centre);
	const double im = point_im(&v->centre);
	double magnitude = sqrt(re * re + im * im);  // correctly rounded on every runtime
	if (magnitude < MB_DIRECT_MAGNITUDE_FLOOR) magnitude = MB_DIRECT_MAGNITUDE_FLOOR;
	return num_log2(magnitude) - MB_DIRECT_PRECISION_BITS;
}

bool view_perturbed(const View* v) {
	return v->scale_log2 < view_direct_limit(v);
}

double view_outermost_scale(void) {
	return s_home_scale_log2 + MB_ZOOM_OUT_LIMIT_OCTAVES;
}

void view_clamp(View* v) {
	point_clamp(&v->centre, MB_CENTRE_LIMIT);
	const double outermost = view_outermost_scale();
	const double deepest = view_depth_limit(v);
	if (v->scale_log2 > outermost) v->scale_log2 = outermost;
	if (v->scale_log2 < deepest) v->scale_log2 = deepest;
}

int64_t view_iteration_limit(const View* v, int scale_steps) {
	double octaves = view_octaves(v);
	if (octaves < 0.0) octaves = 0.0;
	// Whole octaves only, so the limit (and with it every interior pixel) changes once per
	// doubling rather than on every tick of a zoom.
	const int64_t whole = (int64_t)octaves;
	int64_t limit = MB_ITERATIONS_BASE + MB_ITERATIONS_PER_OCTAVE * whole;
	for (int i = 0; i < scale_steps && limit < MB_ITERATIONS_MAX; i++) limit += limit;
	for (int i = 0; i > scale_steps && limit > MB_ITERATIONS_MIN; i--) limit /= 2;
	if (limit < MB_ITERATIONS_MIN) limit = MB_ITERATIONS_MIN;
	if (limit > MB_ITERATIONS_MAX) limit = MB_ITERATIONS_MAX;
	return limit;
}

bool view_equal(const View* a, const View* b) {
	return point_equal(&a->centre, &b->centre) && a->scale_log2 == b->scale_log2;
}
