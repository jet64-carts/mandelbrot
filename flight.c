// carts/mandelbrot/flight.c
#pragma STDC FP_CONTRACT OFF

#include "flight.h"

#include "config.h"
#include "numerics.h"

#define FLIGHT_OCTAVES_PER_SECOND 5.0   // a zoom leg's mean speed, unless the flight is long
#define FLIGHT_MAX_SECONDS 12.0         // ... when its zoom legs are sped up to fit this
#define FLIGHT_MIN_LEG_SECONDS 0.6      // no leg that moves at all is shorter
#define FLIGHT_PAN_SECONDS 1.5          // the pan leg, when the ends are apart
#define FLIGHT_FIT_SCREENS 0.5          // at the pan's scale the ends lie at most this many widths apart
#define FLIGHT_PAN_MIN_SCREENS 1e-3     // ... and closer than this they need no pan
#define FLIGHT_RAMP 0.2                 // the share of a leg over which its speed ramps up, and down
#define FLIGHT_TICKS_PER_SECOND (1.0 / (double)JET_GUEST_TICK_SECONDS)

enum { LEG_OUT, LEG_ACROSS, LEG_IN };

static int seconds_to_ticks(double seconds) {
	return (int)(seconds * FLIGHT_TICKS_PER_SECOND + 0.5);
}

// A zoom leg of `octaves`: at the mean speed, but never shorter than the minimum; none if there is
// nothing to zoom.
static double zoom_seconds(double octaves) {
	if (!(octaves > 0.0)) return 0.0;
	const double seconds = octaves / FLIGHT_OCTAVES_PER_SECOND;
	return seconds > FLIGHT_MIN_LEG_SECONDS ? seconds : FLIGHT_MIN_LEG_SECONDS;
}

void flight_plan(Flight* flight, const View* from, const View* to) {
	flight->from = *from;
	flight->to = *to;
	flight->time = 0.0;
	point_difference(&to->centre, &from->centre, &flight->delta_re, &flight->delta_im);
	const double distance = sqrt(flight->delta_re * flight->delta_re + flight->delta_im * flight->delta_im);

	// The pan's scale: both ends on the screen together, or the start's or the end's if either is
	// wider already; never beyond the outermost view.
	double top = from->scale_log2 > to->scale_log2 ? from->scale_log2 : to->scale_log2;
	if (distance > 0.0) {
		const double fit = num_log2(distance / (FLIGHT_FIT_SCREENS * MB_SCREEN_WIDTH));
		if (fit > top) top = fit;
	}
	if (top > view_outermost_scale()) top = view_outermost_scale();
	flight->top = top;

	const double screens = distance / (num_exp2(top) * MB_SCREEN_WIDTH);
	const double pan = screens > FLIGHT_PAN_MIN_SCREENS ? FLIGHT_PAN_SECONDS : 0.0;
	double out = zoom_seconds(top - from->scale_log2), in = zoom_seconds(top - to->scale_log2);
	const double zooms = FLIGHT_MAX_SECONDS - pan;
	if (out + in > zooms) {  // a long flight: both zooms sped up alike
		const double squeeze = zooms / (out + in);
		out *= squeeze;
		in *= squeeze;
	}
	flight->ticks[LEG_OUT] = seconds_to_ticks(out);
	flight->ticks[LEG_ACROSS] = seconds_to_ticks(pan);
	flight->ticks[LEG_IN] = seconds_to_ticks(in);
}

int flight_ticks(const Flight* flight) {
	return flight->ticks[LEG_OUT] + flight->ticks[LEG_ACROSS] + flight->ticks[LEG_IN];
}

// How far along a leg the camera is at `x` of its time, in [0, 1]: speed ramps up linearly over
// FLIGHT_RAMP, holds, and ramps down over the last FLIGHT_RAMP, so position and speed are continuous
// and the speed is zero at both ends.
static double eased(double x) {
	const double peak = 1.0 / (1.0 - FLIGHT_RAMP);  // the held speed, in legs per leg-time
	if (x <= 0.0) return 0.0;
	if (x >= 1.0) return 1.0;
	if (x < FLIGHT_RAMP) return peak * x * x / (2.0 * FLIGHT_RAMP);
	if (x > 1.0 - FLIGHT_RAMP) return 1.0 - peak * (1.0 - x) * (1.0 - x) / (2.0 * FLIGHT_RAMP);
	return peak * (x - FLIGHT_RAMP / 2.0);
}

bool flight_step(Flight* flight, double rate, View* out) {
	flight->time += rate;
	double t = flight->time;
	if (t < flight->ticks[LEG_OUT]) {
		const double along = eased(t / flight->ticks[LEG_OUT]);
		*out = flight->from;
		out->scale_log2 = flight->from.scale_log2 + (flight->top - flight->from.scale_log2) * along;
		return false;
	}
	t -= flight->ticks[LEG_OUT];
	if (t < flight->ticks[LEG_ACROSS]) {
		// The centre as the target's less the share of the offset still to go: exact at the end.
		const double left = 1.0 - eased(t / flight->ticks[LEG_ACROSS]);
		*out = flight->to;
		point_offset(&out->centre, -flight->delta_re * left, -flight->delta_im * left);
		out->scale_log2 = flight->top;
		return false;
	}
	t -= flight->ticks[LEG_ACROSS];
	if (t < flight->ticks[LEG_IN]) {
		const double along = eased(t / flight->ticks[LEG_IN]);
		*out = flight->to;
		out->scale_log2 = flight->top + (flight->to.scale_log2 - flight->top) * along;
		return false;
	}
	*out = flight->to;
	return true;
}
