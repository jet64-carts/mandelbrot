// carts/mandelbrot/camera.h
//
// THE CAMERA — the player's hand on the view, advanced once per tick in _update.
//
// A drag moves the picture exactly as far as the pointer moved: panning is never eased, or the
// plane would slide out from under the hand. Zoom is eased: the wheel and the keys move a TARGET
// scale, and each tick the view closes a fixed share of the distance to it, about a point that
// stays still on the screen — the pointer for the wheel, the screen's centre for keys and
// triggers. Keyboard and stick panning ease their velocity, so a tap neither jerks nor
// overshoots. A bookmark is a FLIGHT (flight.h), which any input cancels where it has got to, the
// hand taking over. The TOUR flies through the bookmarks in CAMERA_TOUR's order, dwelling
// CAMERA_TOUR_DWELL_SECONDS at each, until any input stops it. Every rate is per tick of the fixed
// 120 Hz clock, so the camera moves the same on every host.
#ifndef MANDELBROT_CAMERA_H
#define MANDELBROT_CAMERA_H

#include "controls.h"
#include "field.h"
#include "flight.h"
#include "view.h"

// ATTRACT: after this long at rest with no input, the tour starts by itself.
#define CAMERA_ATTRACT_TICKS 7200  // sixty seconds

typedef struct {
	View view;                   // the view shown
	double zoom_target;          // the scale_log2 the zoom is easing toward
	double anchor_x, anchor_y;   // the screen point the zoom holds still
	double velocity_x, velocity_y;  // keyboard and stick panning, screen pixels per tick
	int iteration_steps;         // the player's iteration scale, in powers of two
	bool grabbing;               // the hand is on the picture (a drag is held)
	bool flying;                 // a flight is under way
	Flight flight;
	double readiness;            // the share of the last picture at its own coordinates (camera_set_readiness)
	double flight_rate;          // the flight clock's rate, easing toward what the readiness allows
	bool touring;                // the tour is under way
	int tour_stop;               // the tour's current stop, an index into its order
	int tour_dwell;              // ticks spent at the current stop
	int idle_ticks;              // ticks at rest with no input (ATTRACT)
	uint32_t epoch;              // counts declared jumps: a reset always shows a fresh picture
} Camera;

/** A camera at home. */
void camera_init(Camera* camera);

/** One tick of the camera under the player's controls. */
void camera_update(Camera* camera, const Controls* controls);

/**
 * What the last draw showed: the share of its pixels at their own coordinates (reused or better), in
 * [0, 1], which governs a flight's speed (flight.h). A camera without a picture reads 1.
 */
void camera_set_readiness(Camera* camera, double share);

/** A declared jump to `view`: the picture restarts there (the harness's way to a deep suite view). */
void camera_jump(Camera* camera, const View* view);

/**
 * At rest, a view whose grid is one lattice a fraction of a pixel off its ideal — as a pan leaves
 * it — moves by that fraction, onto the lattice. The picture already holds the exact values there,
 * so nothing is recomputed; snapping the grid to the view instead would recompute every pixel, to
 * move the picture by less than half a pixel. Call after camera_update with what the last draw left:
 * the spec the field `shown`, and its grid's offset (field_lattice_offset). It acts only when that
 * spec is this very view.
 */
void camera_rest_on_grid(Camera* camera, const FieldSpec* shown, double d_re, double d_im);

/** What the field should show for this camera: its view, iteration limit and epoch. */
FieldSpec camera_spec(const Camera* camera);

/**
 * True when the view will not move without a new hand: no flight is under way, the zoom has reached
 * its target, no pan is coasting, and no drag is held (a held drag may move on any tick).
 */
bool camera_at_rest(const Camera* camera);

#endif // MANDELBROT_CAMERA_H
