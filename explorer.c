// carts/mandelbrot/explorer.c
#pragma STDC FP_CONTRACT OFF

#include "explorer.h"

#include "camera.h"
#include "config.h"
#include "controls.h"
#include "field.h"
#include "hud.h"
#include "look.h"
#include "numerics.h"
#include "palette.h"
#include "scheduler.h"

// The title card shows for this long unless the player starts sooner, then fades over the second
// span.
#define TITLE_TICKS 600
#define TITLE_FADE_TICKS 60
#define TITLE_NOT_FADING (-1)
// The status line holds this long after the last input, then fades over the second span.
#define STATUS_HOLD_TICKS 480
#define STATUS_FADE_TICKS 60

static Field s_field;
static Camera s_camera;
static bool s_show_status;
static bool s_show_help;
static int s_tick;
static int s_title_fade_from;  // the tick the title began to fade, or TITLE_NOT_FADING
static double s_title_share;   // how much of the title card shows, in [0, 1]
static int s_idle_ticks;       // ticks since the player's last input
static double s_status_share;  // how much of the status line shows, in [0, 1]
static uint32_t s_shaded_epoch;

bool explorer_init(void) {
	jet_display_set_resolution(MB_SCREEN_WIDTH, MB_SCREEN_HEIGHT);
	if (!field_init(&s_field)) return false;
	if (!look_init()) return false;
	view_init();
	palette_init();
	camera_init(&s_camera);
	const SchedulerOptions options = scheduler_default_options();
	if (!scheduler_init(&s_field, &options)) return false;
	s_show_status = true;
	s_show_help = false;
	s_tick = 0;
	s_title_fade_from = TITLE_NOT_FADING;
	s_title_share = 1.0;
	s_idle_ticks = 0;
	s_status_share = 1.0;

	// The home picture's restart runs now, during boot, rather than in the first tick: the first
	// frame is presented already refined, and no tick ever pays for a cold start.
	const FieldSpec home = camera_spec(&s_camera);
	scheduler_settle(&home, true);
	s_shaded_epoch = home.epoch;
	printf("MANDELBROT: %u bytes of WRAM, %d workers\n", (unsigned)jet_memory_usage_bytes(), MB_WORKER_COUNT);
	return true;
}

void explorer_update(void) {
	const Controls controls = controls_read();
	camera_set_readiness(&s_camera, look_ready_share());
	const bool was_touring = s_camera.touring;
	camera_update(&s_camera, &controls);
	if (s_camera.touring && !was_touring) palette_set_cycling(true);  // a tour is shown with the colour moving
	double d_re, d_im;
	if (field_lattice_offset(&s_field, &d_re, &d_im)) camera_rest_on_grid(&s_camera, &s_field.spec, d_re, d_im);
	if (controls.toggle_cycle) palette_set_cycling(!palette_cycling());
	if (controls.toggle_hud) s_show_status = !s_show_status;
	if (controls.toggle_help) s_show_help = !s_show_help;
	if (controls.screenshot) jet_host_screenshot();
	if (controls.next_style) look_set_style((LookStyle)((look_style() + 1) % LOOK_STYLE_COUNT));
	if (controls.next_gradient) palette_next_gradient();
	if (controls.density != 0) look_adjust_density(controls.density);
	palette_advance();

	s_tick++;
	if (s_title_fade_from == TITLE_NOT_FADING && (controls_any(&controls) || s_tick >= TITLE_TICKS)) {
		s_title_fade_from = s_tick;
	}
	if (s_title_fade_from != TITLE_NOT_FADING) {
		s_title_share = num_clamp(1.0 - (double)(s_tick - s_title_fade_from) / TITLE_FADE_TICKS, 0.0, 1.0);
	}
	palette_set_title_fade(s_title_share);

	s_idle_ticks = controls_any(&controls) ? 0 : s_idle_ticks + 1;
	s_status_share = num_clamp(1.0 - (double)(s_idle_ticks - STATUS_HOLD_TICKS) / STATUS_FADE_TICKS, 0.0, 1.0);
	palette_set_status_fade(s_status_share);
}

void explorer_draw(void) {
	const FieldSpec spec = camera_spec(&s_camera);
	const bool at_rest = camera_at_rest(&s_camera);
	scheduler_settle(&spec, at_rest);

	// The look reshades the frame when the picture or its scale moved; a jump snaps the scale.
	const SchedulerStats* stats = scheduler_stats();
	const bool changed = stats->reconcile >= 0 || stats->batch_work > 0 || stats->settle_work > 0;
	const bool jump = spec.epoch != s_shaded_epoch;
	s_shaded_epoch = spec.epoch;
	look_update(&s_field, changed, jump);
	palette_set_layout(look_style() == LOOK_SILK ? PAL_LAYOUT_RING : PAL_LAYOUT_PLANE);
	palette_commit();
	memcpy(VRAM, look_frame(), MB_SCREEN_PIXELS);
	HudStatus status;
	status.octaves = view_octaves(&spec.view);
	status.max_iterations = spec.max_iterations;
	status.progress = scheduler_progress();
	status.at_depth_limit = spec.view.scale_log2 <= view_depth_limit(&spec.view);
	hud_draw(&status, s_show_status && s_status_share > 0.0, s_show_help, s_title_share > 0.0);

	scheduler_launch(at_rest);
}

void explorer_end(void) {
	scheduler_end();
}
