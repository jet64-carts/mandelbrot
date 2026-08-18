// carts/mandelbrot/controls.h
//
// CONTROLS — one tick of mouse, keyboard and gamepad, read into what the player asked for.
//
// Nothing else in the cart reads an input device. The camera and the explorer act on a Controls
// value, which is what lets the harness drive them with a scripted autopilot instead of a hand.
//
//   mouse      drag with the left button to pan, 1:1; the wheel zooms about the pointer
//   keyboard   arrows or W A S D pan; E zooms in and Q out while held; R returns home; 1 to 9
//              fly to the bookmarks (bookmarks.h); T tours them; - and = halve and double the
//              iteration limit;
//              C changes the style; V the gradient; [ and ] thin and thicken the colour bands; SPACE pauses the
//              colour cycle; H hides the HUD; / opens the help; P saves a screenshot
//   gamepad    left stick pans; right trigger zooms in, left trigger out; START returns home;
//              LB and RB halve and double the iteration limit; A pauses the colour cycle;
//              SELECT hides the HUD; Y opens the help
//
// No control is TAB, a function key or a shifted character (docs/claude.md note #5).
#ifndef MANDELBROT_CONTROLS_H
#define MANDELBROT_CONTROLS_H

#include "jet64.h"

typedef struct {
	// Continuous, this tick.
	double drag_x, drag_y;     // screen pixels the picture was dragged
	double wheel_notches;      // the wheel's turn; positive zooms in
	double zoom_rate;          // octaves per second requested by keys or triggers; positive is in
	double pan_x, pan_y;       // panning direction in [-1, 1] per axis, screen axes (y down)
	double pointer_x, pointer_y;  // the mouse's screen position, the wheel's anchor
	bool grabbing;             // the drag button is held

	// Pressed this tick.
	bool reset;
	int bookmark;  // the bookmark to visit, 1 to BOOKMARK_COUNT, or 0
	bool tour;     // start the tour
	bool iterations_up, iterations_down;
	bool toggle_cycle, toggle_hud, toggle_help;
	bool next_style;
	bool next_gradient;
	int density;  // steps of colour density: +1 thicker bands, -1 thinner, 0 none
	bool screenshot;
} Controls;

/** This tick's controls. Call once per _update. */
Controls controls_read(void);

/** True when the player did anything at all this tick (a pointer that only moved is not an input). */
static inline bool controls_any(const Controls* c) {
	return c->drag_x != 0.0 || c->drag_y != 0.0 || c->wheel_notches != 0.0 || c->zoom_rate != 0.0 ||
	       c->pan_x != 0.0 || c->pan_y != 0.0 || c->grabbing || c->reset || c->bookmark != 0 || c->tour ||
	       c->iterations_up || c->iterations_down || c->toggle_cycle || c->toggle_hud || c->toggle_help ||
	       c->next_style || c->next_gradient || c->density != 0 || c->screenshot;
}

#endif // MANDELBROT_CONTROLS_H
