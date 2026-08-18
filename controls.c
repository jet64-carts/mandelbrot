// carts/mandelbrot/controls.c
#pragma STDC FP_CONTRACT OFF

#include "controls.h"

#include "bookmarks.h"
#include "config.h"
#include "numerics.h"

// The keys, named once; each letter key is its own ASCII code (jet_input.h).
#define KEY_ZOOM_IN 'E'
#define KEY_ZOOM_OUT 'Q'
#define KEY_PAN_LEFT 'A'
#define KEY_PAN_RIGHT 'D'
#define KEY_PAN_UP 'W'
#define KEY_PAN_DOWN 'S'
#define KEY_RESET 'R'
#define KEY_ITERATIONS_DOWN '-'
#define KEY_ITERATIONS_UP '='
#define KEY_CYCLE JET_KEY_SPACE
#define KEY_HUD 'H'
#define KEY_HELP '/'
#define KEY_SCREENSHOT 'P'
#define KEY_STYLE 'C'
#define KEY_TOUR 'T'
#define KEY_GRADIENT 'V'
#define KEY_DENSITY_DOWN '['
#define KEY_DENSITY_UP ']'
#define KEY_FIRST_BOOKMARK '1'  // and on up the top row: '1' + i visits bookmark i + 1

// A trigger reads 0 at rest and 1 pulled home, with no deadzone of its own; this much travel is
// ignored so a resting trigger never creeps the zoom.
#define TRIGGER_DEADZONE 0.1

static double key_axis(int negative, int alternative_negative, int positive, int alternative_positive) {
	double axis = 0.0;
	if (jet_key_down(negative) || jet_key_down(alternative_negative)) axis -= 1.0;
	if (jet_key_down(positive) || jet_key_down(alternative_positive)) axis += 1.0;
	return axis;
}

// Input is the one door through which a number the cart did not compute can enter it; a NaN or
// an infinity let through would poison the view for good, since every clamp compares false
// against a NaN. So anything non-finite reads as no movement at all.
static double finite_or_zero(double v) {
	return v - v == 0.0 ? v : 0.0;  // false for NaN and for either infinity
}

static double trigger(int pad, int axis) {
	const double pull = num_clamp(finite_or_zero(jet_pad_axis(pad, axis)), 0.0, 1.0);
	if (pull <= TRIGGER_DEADZONE) return 0.0;
	return (pull - TRIGGER_DEADZONE) / (1.0 - TRIGGER_DEADZONE);
}

Controls controls_read(void) {
	Controls c;
	memset(&c, 0, sizeof c);

	// The pointer anchors the wheel's zoom, so it is held to the screen it points into.
	c.pointer_x = num_clamp(finite_or_zero(jet_mouse_x()), 0.0, MB_SCREEN_WIDTH);
	c.pointer_y = num_clamp(finite_or_zero(jet_mouse_y()), 0.0, MB_SCREEN_HEIGHT);
	c.grabbing = jet_mouse_down(MOUSE_BTN_LEFT);
	if (c.grabbing) {
		c.drag_x = finite_or_zero(jet_mouse_dx());
		c.drag_y = finite_or_zero(jet_mouse_dy());
	}
	c.wheel_notches = finite_or_zero(jet_mouse_wheel());

	c.pan_x = key_axis(JET_KEY_LEFT, KEY_PAN_LEFT, JET_KEY_RIGHT, KEY_PAN_RIGHT);
	c.pan_y = key_axis(JET_KEY_UP, KEY_PAN_UP, JET_KEY_DOWN, KEY_PAN_DOWN);
	double zoom = 0.0;
	if (jet_key_down(KEY_ZOOM_IN)) zoom += MB_KEY_ZOOM_OCTAVES_PER_SECOND;
	if (jet_key_down(KEY_ZOOM_OUT)) zoom -= MB_KEY_ZOOM_OCTAVES_PER_SECOND;

	for (int pad = 0; pad < JET_GPAD_MAX_PADS; pad++) {
		if (!jet_pad_connected(pad)) continue;
		float stick_x, stick_y;
		jet_pad_stick(pad, JET_PAD_STICK_LEFT, &stick_x, &stick_y);
		c.pan_x += finite_or_zero(stick_x);
		c.pan_y += finite_or_zero(stick_y);
		zoom += (trigger(pad, JET_PAD_AXIS_RT) - trigger(pad, JET_PAD_AXIS_LT)) *
		        MB_TRIGGER_ZOOM_OCTAVES_PER_SECOND;
	}
	c.pan_x = num_clamp(c.pan_x, -1.0, 1.0);
	c.pan_y = num_clamp(c.pan_y, -1.0, 1.0);
	c.zoom_rate = zoom;

	c.reset = jet_key_pressed(KEY_RESET) || jet_pad_any_pressed(JET_PAD_BTN_START);
	c.iterations_down = jet_key_pressed(KEY_ITERATIONS_DOWN) || jet_key_pressed(JET_KEY_KP_SUBTRACT) ||
	                    jet_pad_any_pressed(JET_PAD_BTN_LB);
	c.iterations_up = jet_key_pressed(KEY_ITERATIONS_UP) || jet_key_pressed(JET_KEY_KP_ADD) ||
	                  jet_pad_any_pressed(JET_PAD_BTN_RB);
	c.toggle_cycle = jet_key_pressed(KEY_CYCLE) || jet_pad_any_pressed(JET_PAD_BTN_A);
	c.toggle_hud = jet_key_pressed(KEY_HUD) || jet_pad_any_pressed(JET_PAD_BTN_SELECT);
	c.toggle_help = jet_key_pressed(KEY_HELP) || jet_pad_any_pressed(JET_PAD_BTN_Y);
	c.screenshot = jet_key_pressed(KEY_SCREENSHOT);
	c.next_style = jet_key_pressed(KEY_STYLE) || jet_pad_any_pressed(JET_PAD_BTN_X);
	c.tour = jet_key_pressed(KEY_TOUR) || jet_pad_any_pressed(JET_PAD_BTN_B);
	c.next_gradient = jet_key_pressed(KEY_GRADIENT) || jet_pad_any_pressed(JET_PAD_BTN_R3);
	c.density = (jet_key_pressed(KEY_DENSITY_UP) ? 1 : 0) - (jet_key_pressed(KEY_DENSITY_DOWN) ? 1 : 0);
	for (int i = 0; i < BOOKMARK_COUNT; i++) {
		if (jet_key_pressed(KEY_FIRST_BOOKMARK + i)) c.bookmark = i + 1;
	}
	return c;
}
