// carts/mandelbrot/view.h
//
// THE VIEW — which part of the plane the screen shows, and at what iteration limit.
//
// A view is a centre and a scale. The scale is held as the base-two logarithm of the PIXEL
// SPACING (the distance in the plane between neighbouring pixel centres), so zooming is adding,
// easing is linear, and a depth of 2^-800 is as ordinary a number as 2^-8. Screen point (x, y)
// shows
//
//     c = centre + (view_offset_x(x) h, view_offset_y(y) h),   h = 2^scale_log2
//
// with the imaginary axis pointing up the screen, as it does on paper; pixel (i, j) is the point
// at its own centre, (i + 1/2, j + 1/2). A view deeper than the CROSSOVER is computed by
// perturbation (view_perturbed); none may go deeper than the DEPTH LIMIT (config.h).
#ifndef MANDELBROT_VIEW_H
#define MANDELBROT_VIEW_H

#include "config.h"
#include "point.h"

// Screen coordinates are continuous, in pixels from the top-left corner: pixel (i, j) covers
// [i, i + 1) x [j, j + 1), so its centre sits at (i + VIEW_PIXEL_CENTRE, j + VIEW_PIXEL_CENTRE).
#define VIEW_PIXEL_CENTRE 0.5
#define VIEW_SCREEN_CENTRE_X (MB_SCREEN_WIDTH / 2.0)
#define VIEW_SCREEN_CENTRE_Y (MB_SCREEN_HEIGHT / 2.0)

/** The offset of screen coordinate x from the screen's centre, in pixels: the real axis. */
static inline double view_offset_x(double screen_x) {
	return screen_x - VIEW_SCREEN_CENTRE_X;
}

/** The offset of screen coordinate y from the screen's centre, in pixels, the imaginary axis up. */
static inline double view_offset_y(double screen_y) {
	return VIEW_SCREEN_CENTRE_Y - screen_y;
}

typedef struct {
	DeepPoint centre;
	double scale_log2;  // log2 of the pixel spacing
} View;

/** Computes the home framing. Call once from _init, before any other view function. */
void view_init(void);

/** The home view: the whole set, framed. */
View view_home(void);

/** The pixel spacing h = 2^scale_log2. */
double view_spacing(const View* v);

/** Octaves of magnification from home (negative when zoomed out beyond it). */
double view_octaves(const View* v);

/**
 * Zooms in by `octaves` (out when negative) about the screen point `anchor_x`, `anchor_y` (in
 * screen pixels), which stays where it is on the screen.
 */
void view_zoom_about(View* v, double octaves, double anchor_x, double anchor_y);

/** Slides the picture by (dx, dy) screen pixels, as a hand dragging it would. */
void view_pan_pixels(View* v, double dx, double dy);

/**
 * Holds the view inside its limits: no further out than MB_ZOOM_OUT_LIMIT_OCTAVES beyond home,
 * no deeper than the depth limit (MB_DEPTH_LIMIT_LOG2), and no further from the origin than
 * MB_CENTRE_LIMIT.
 */
void view_clamp(View* v);

/** The deepest scale_log2 a view may take (MB_DEPTH_LIMIT_LOG2). */
double view_depth_limit(const View* v);

/** The deepest scale_log2 direct double iteration resolves at this view's centre: the crossover. */
double view_direct_limit(const View* v);

/** True when the view is deeper than direct iteration resolves, and computes by perturbation. */
bool view_perturbed(const View* v);

/** The outermost scale_log2 a view may take: MB_ZOOM_OUT_LIMIT_OCTAVES beyond home. */
double view_outermost_scale(void);

/**
 * The automatic iteration limit for this view: MB_ITERATIONS_BASE plus MB_ITERATIONS_PER_OCTAVE
 * for each whole octave of zoom, scaled by 2^scale_steps and clamped to
 * [MB_ITERATIONS_MIN, MB_ITERATIONS_MAX].
 */
int64_t view_iteration_limit(const View* v, int scale_steps);

/** True when the two views frame exactly the same picture. */
bool view_equal(const View* a, const View* b);

#endif // MANDELBROT_VIEW_H
