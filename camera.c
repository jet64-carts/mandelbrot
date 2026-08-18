// carts/mandelbrot/camera.c
#pragma STDC FP_CONTRACT OFF

#include "camera.h"

#include "bookmarks.h"
#include "config.h"
#include "numerics.h"

#define TICK_SECONDS ((double)JET_GUEST_TICK_SECONDS)

// The tour's stops, as bookmark indices: out through the valleys and spirals, down to the two deep
// minibrots, and home, then round again.
static const int CAMERA_TOUR[] = {2, 1, 3, 4, 5, 6, 7, 8, 0};
#define CAMERA_TOUR_STOPS ((int)(sizeof(CAMERA_TOUR) / sizeof(CAMERA_TOUR[0])))
#define CAMERA_TOUR_DWELL_SECONDS 6.0  // time to converge, and to look
#define CAMERA_TOUR_DWELL_TICKS ((int)(CAMERA_TOUR_DWELL_SECONDS / TICK_SECONDS))

// THE FLIGHT GOVERNOR. A flight's clock runs at full rate while at least CAMERA_READY_SHARE of the
// last picture sat at its own coordinates, and in proportion below it, never slower than
// CAMERA_RATE_FLOOR (so it cannot stall); the rate eases toward that by CAMERA_RATE_EASE a tick, so
// the motion stays smooth. The readiness is the look's count of the drawn frame (look.h), so the
// governor is as deterministic as the picture.
#define CAMERA_READY_SHARE 0.5
#define CAMERA_RATE_FLOOR 0.02
#define CAMERA_RATE_EASE 0.05

// Everything but the epoch, which only ever counts up.
static void go_home(Camera* camera) {
	camera->view = view_home();
	camera->zoom_target = camera->view.scale_log2;
	camera->anchor_x = VIEW_SCREEN_CENTRE_X;
	camera->anchor_y = VIEW_SCREEN_CENTRE_Y;
	camera->velocity_x = camera->velocity_y = 0.0;
	camera->iteration_steps = 0;
	camera->grabbing = false;
	camera->flying = false;
	camera->touring = false;
}

// Starts the flight to bookmark `index` from wherever the camera is.
static void fly_to(Camera* camera, int index) {
	View visit;
	if (!bookmark_view(index, &visit)) return;
	view_clamp(&visit);
	flight_plan(&camera->flight, &camera->view, &visit);
	camera->flying = true;
	camera->flight_rate = 1.0;
	camera->velocity_x = camera->velocity_y = 0.0;
	camera->anchor_x = VIEW_SCREEN_CENTRE_X;
	camera->anchor_y = VIEW_SCREEN_CENTRE_Y;
}

void camera_init(Camera* camera) {
	go_home(camera);
	camera->epoch = 0;
	camera->readiness = 1.0;
	camera->flight_rate = 1.0;
	camera->idle_ticks = 0;
}

void camera_set_readiness(Camera* camera, double share) {
	camera->readiness = num_clamp(share, 0.0, 1.0);
}

// Eases a pan velocity toward the one the controls ask for. A velocity that has all but died
// away is stopped outright: a view that keeps creeping by a thousandth of a pixel is a view that
// never settles, and so never finishes computing.
static double ease_velocity(double velocity, double wanted) {
	const double eased = velocity + (wanted - velocity) * MB_PAN_EASE_PER_TICK;
	const bool resting = wanted == 0.0 && eased < MB_PAN_SNAP_PIXELS && eased > -MB_PAN_SNAP_PIXELS;
	return resting ? 0.0 : eased;
}

void camera_jump(Camera* camera, const View* view) {
	go_home(camera);
	camera->view = *view;
	view_clamp(&camera->view);
	camera->zoom_target = camera->view.scale_log2;
	camera->epoch++;
}

void camera_update(Camera* camera, const Controls* controls) {
	if (controls->reset) {
		go_home(camera);
		camera->epoch++;  // a declared jump: the picture restarts even if home is close by
		return;
	}
	camera->idle_ticks = controls_any(controls) || !camera_at_rest(camera) ? 0 : camera->idle_ticks + 1;
	if (controls->tour || camera->idle_ticks >= CAMERA_ATTRACT_TICKS) {
		camera->idle_ticks = 0;
		camera->touring = true;
		camera->tour_stop = 0;
		camera->tour_dwell = 0;
		fly_to(camera, CAMERA_TOUR[0]);
	} else if (controls->bookmark > 0 && controls->bookmark <= BOOKMARK_COUNT) {
		camera->touring = false;
		fly_to(camera, controls->bookmark - 1);
	} else if ((camera->flying || camera->touring) && controls_any(controls)) {
		camera->flying = camera->touring = false;  // the hand takes over where the camera has got to
		camera->zoom_target = camera->view.scale_log2;
	} else if (camera->touring && !camera->flying && ++camera->tour_dwell >= CAMERA_TOUR_DWELL_TICKS) {
		camera->tour_stop = (camera->tour_stop + 1) % CAMERA_TOUR_STOPS;
		camera->tour_dwell = 0;
		fly_to(camera, CAMERA_TOUR[camera->tour_stop]);
	}
	if (camera->flying) {
		const double allowed = num_clamp(camera->readiness / CAMERA_READY_SHARE, CAMERA_RATE_FLOOR, 1.0);
		camera->flight_rate += (allowed - camera->flight_rate) * CAMERA_RATE_EASE;
		camera->flying = !flight_step(&camera->flight, camera->flight_rate, &camera->view);
		view_clamp(&camera->view);
		camera->zoom_target = camera->view.scale_log2;
		return;
	}

	// The hand: the picture follows the pointer exactly.
	camera->grabbing = controls->grabbing;
	view_pan_pixels(&camera->view, controls->drag_x, controls->drag_y);

	// Keys and sticks: ease toward the requested velocity, then slide the picture against it,
	// so pressing right travels right.
	const double full_speed = MB_PAN_SCREENS_PER_SECOND * MB_SCREEN_HEIGHT * TICK_SECONDS;
	camera->velocity_x = ease_velocity(camera->velocity_x, controls->pan_x * full_speed);
	camera->velocity_y = ease_velocity(camera->velocity_y, controls->pan_y * full_speed);
	view_pan_pixels(&camera->view, -camera->velocity_x, -camera->velocity_y);

	// The zoom's target: the wheel anchors on the pointer, keys and triggers on the centre.
	if (controls->wheel_notches != 0.0) {
		camera->zoom_target -= controls->wheel_notches * MB_WHEEL_OCTAVES_PER_NOTCH;
		camera->anchor_x = controls->pointer_x;
		camera->anchor_y = controls->pointer_y;
	} else if (controls->zoom_rate != 0.0) {
		camera->zoom_target -= controls->zoom_rate * TICK_SECONDS;
		camera->anchor_x = VIEW_SCREEN_CENTRE_X;
		camera->anchor_y = VIEW_SCREEN_CENTRE_Y;
	}

	// Ease toward the target, closing a fixed share each tick and the last sliver outright — and
	// then landing on the target exactly, since scale + (target - scale) need not round to the
	// target, and a view that misses by an ulp would never settle.
	view_clamp(&camera->view);
	camera->zoom_target =
	    num_clamp(camera->zoom_target, view_depth_limit(&camera->view), view_outermost_scale());
	const double remaining = camera->zoom_target - camera->view.scale_log2;
	if (remaining != 0.0) {
		const bool last = remaining < MB_ZOOM_SNAP_OCTAVES && remaining > -MB_ZOOM_SNAP_OCTAVES;
		const double octaves = last ? -remaining : -remaining * MB_ZOOM_EASE_PER_TICK;
		view_zoom_about(&camera->view, octaves, camera->anchor_x, camera->anchor_y);
		if (last) camera->view.scale_log2 = camera->zoom_target;
		view_clamp(&camera->view);
	}

	if (controls->iterations_up && camera->iteration_steps < MB_ITERATION_SCALE_MAX_STEPS) {
		camera->iteration_steps++;
	}
	if (controls->iterations_down && camera->iteration_steps > -MB_ITERATION_SCALE_MAX_STEPS) {
		camera->iteration_steps--;
	}
}

bool camera_at_rest(const Camera* camera) {
	return !camera->flying && !camera->grabbing && camera->zoom_target == camera->view.scale_log2 &&
	       camera->velocity_x == 0.0 && camera->velocity_y == 0.0;
}

void camera_rest_on_grid(Camera* camera, const FieldSpec* shown, double d_re, double d_im) {
	if (!camera_at_rest(camera) || shown->epoch != camera->epoch || !view_equal(&shown->view, &camera->view)) return;
	point_offset(&camera->view.centre, d_re, d_im);
	view_clamp(&camera->view);
}

FieldSpec camera_spec(const Camera* camera) {
	FieldSpec spec;
	spec.view = camera->view;
	spec.max_iterations = view_iteration_limit(&camera->view, camera->iteration_steps);
	spec.epoch = camera->epoch;
	spec.perturbed = view_perturbed(&camera->view);
	return spec;
}
