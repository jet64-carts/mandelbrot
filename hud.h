// carts/mandelbrot/hud.h
//
// THE HUD — a single quiet line of what the explorer is doing, which fades when the player leaves the
// controls alone and returns at a touch; a help card on request; and a title card at boot that fades
// as the player begins.
//
// Drawn over the picture after it has been copied to the screen, in the HUD's own palette inks,
// which the colour cycle never touches.
#ifndef MANDELBROT_HUD_H
#define MANDELBROT_HUD_H

#include "jet64.h"

// What the status line reports.
typedef struct {
	double octaves;          // the magnification relative to the home view, as a power of two
	int64_t max_iterations;
	double progress;         // the share of pixels computed exactly, in [0, 1]
	bool at_depth_limit;     // the view is at the depth limit (config.h: MB_DEPTH_LIMIT_LOG2)
} HudStatus;

/**
 * Draws the status line (when `show_status`), the help card (when `show_help`) and the title card
 * (when `show_title`, in the PAL_TITLE ink the explorer fades).
 */
void hud_draw(const HudStatus* status, bool show_status, bool show_help, bool show_title);

#endif // MANDELBROT_HUD_H
