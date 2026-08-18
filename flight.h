// carts/mandelbrot/flight.h
//
// FLIGHTS — the camera's way to a bookmark: out, across, and in, never a cut.
//
// A flight is three legs in log-scale space. It ZOOMS OUT about the start's centre until both ends
// fit on the screen together (FLIGHT_FIT_SCREENS apart at most) or as far as the start already is;
// PANS across to the target's centre at that scale; and ZOOMS IN about the target's centre to its
// scale. Each leg eases in and out (its speed ramps over FLIGHT_RAMP of the leg at each end), so the
// camera never jerks into motion or stops dead, and a leg that is not needed takes no time.
//
// Built to be exact at any depth: the scale is carried as log2 of the pixel spacing, as the view
// carries it, so a target of 10^-130 costs nothing; a zoom leg holds its centre fixed, so the
// target's bignum centre is the view's centre exactly for the whole descent; and the pan leg's
// centre is the target's minus a double fraction of the offset between the ends, which reaches the
// target exactly at the leg's end. The final view is the target, bit for bit.
//
// Durations are in ticks of the fixed 120 Hz clock: a zoom leg flies FLIGHT_OCTAVES_PER_SECOND,
// and a long flight is planned to fit FLIGHT_MAX_SECONDS. But THE FLIGHT NEVER OUTRUNS THE PICTURE:
// its clock advances each tick by a RATE the caller governs from the picture's readiness (the camera:
// the share of pixels at their own coordinates the last draw showed), so a deep descent slows to
// what the pool can keep detailed and a shallow one flies at the planned speed.
#ifndef MANDELBROT_FLIGHT_H
#define MANDELBROT_FLIGHT_H

#include "view.h"

#define FLIGHT_LEGS 3

typedef struct {
	View from, to;
	double top;                // the pan leg's scale_log2
	double delta_re, delta_im; // to's centre less from's, in the plane
	int ticks[FLIGHT_LEGS];    // each leg's planned length: out, across, in
	double time;               // the flight's clock, in planned ticks
} Flight;

/** Plans a flight from `from` to `to` (both inside the view limits). */
void flight_plan(Flight* flight, const View* from, const View* to);

/** The flight's length in ticks. */
int flight_ticks(const Flight* flight);

/**
 * Flies one tick with the clock advancing `rate` planned ticks (in (0, 1]), and writes the view
 * there; true once the flight has arrived (at `to` exactly).
 */
bool flight_step(Flight* flight, double rate, View* out);

#endif // MANDELBROT_FLIGHT_H
