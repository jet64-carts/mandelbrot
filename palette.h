// carts/mandelbrot/palette.h
//
// THE PALETTE — how 256 indices are shared out, and the colour cycle that animates them.
//
// Index 0 is the set's interior: with the Copper on, every index-0 pixel shows that scanline's
// Copper colour, so the interior gets a smooth vertical gradient in true colour for nothing. The
// RING, indices 1..PAL_RING_SIZE, is one full turn of a cyclic gradient (or a plane of hues by
// shades, below): shading maps the smooth iteration count onto it (look.h), and the colour cycle
// rewrites it at a
// fractional phase each tick, so colour flows outward from the set without one pixel being
// redrawn. The top entries are the HUD's inks, which the cycle never touches.
//
// The gradient is interpolated in OKLab, a perceptual space, so its steps look even to the eye;
// it is sampled PAL_GRADIENT_SUBSTEPS times per ring entry, which is the resolution the cycle's
// phase moves in.
#ifndef MANDELBROT_PALETTE_H
#define MANDELBROT_PALETTE_H

#include "jet64.h"

#define PAL_INTERIOR 0     // the set itself: the Copper shows through
#define PAL_RING_FIRST 1   // the gradient ring
#define PAL_RING_SIZE 240
#define PAL_HUD_SHADOW 248 // the HUD's inks, darkest to brightest
#define PAL_HUD_PANEL 249
#define PAL_HUD_DIM 250
#define PAL_HUD_TEXT 251
#define PAL_HUD_ACCENT 252
#define PAL_INK 253        // the Silk style's edges
#define PAL_TITLE 254      // the title card's ink, which fades
#define PAL_STATUS 255     // the status line's ink, which fades when the player is idle
#define PAL_ENTRIES 256

// The ring's entries have two LAYOUTS (look.h's styles). As a ring, entry 1 + e is hue e of
// PAL_RING_SIZE. As a PLANE, entry 1 + s x PAL_PLANE_HUES + h is hue h of PAL_PLANE_HUES at shade
// s of PAL_PLANE_SHADES: the gradient's colour scaled toward black in OKLab, to PAL_SHADE_LOW of
// its own lightness at shade 0 and to all of it at the top shade — in steps that close toward the
// top (palette_shade_light), where the relief's gentle slopes live and a coarse step's dither would
// show as texture.
#define PAL_PLANE_HUES 30
#define PAL_PLANE_SHADES 8
#define PAL_SHADE_LOW 0.25

_Static_assert(PAL_RING_FIRST + PAL_RING_SIZE <= PAL_HUD_SHADOW, "the ring must not reach the HUD");
_Static_assert(PAL_PLANE_HUES * PAL_PLANE_SHADES == PAL_RING_SIZE, "the plane fills the ring's entries");
_Static_assert(PAL_STATUS < PAL_ENTRIES, "the inks must fit the palette");

typedef enum {
	PAL_LAYOUT_RING,
	PAL_LAYOUT_PLANE
} PaletteLayout;

/**
 * The lightness of plane shade `level` (continuous in [0, PAL_PLANE_SHADES - 1]), as a share of the
 * gradient's own: 1 - (1 - PAL_SHADE_LOW)(1 - level / (PAL_PLANE_SHADES - 1))^2.
 */
static inline double palette_shade_light(double level) {
	const double below_top = 1.0 - level / (double)(PAL_PLANE_SHADES - 1);
	return 1.0 - (1.0 - PAL_SHADE_LOW) * below_top * below_top;
}

/** Builds the gradient, paints the HUD's inks and the interior's Copper gradient. */
void palette_init(void);

/** One tick of the colour cycle, from _update. Does nothing while cycling is paused. */
void palette_advance(void);

/** Pauses or resumes the colour cycle; the phase holds where it is. */
void palette_set_cycling(bool on);
bool palette_cycling(void);

/** The gradients: Classic, Ember, Glacier (palette.c), all cyclic, all mixed in OKLab. */
#define PAL_GRADIENT_COUNT 3

/** The next gradient from the next commit on. */
void palette_next_gradient(void);

/** Lays the ring's entries out as a ring or as a plane from the next commit on. */
void palette_set_layout(PaletteLayout layout);

/** The title card's ink at `share` of its full brightness, in [0, 1]. */
void palette_set_title_fade(double share);

/** The status line's ink at `share` of its full brightness, in [0, 1]. */
void palette_set_status_fade(double share);

/** Writes the ring at the current phase into the hardware palette, from _draw. */
void palette_commit(void);

#endif // MANDELBROT_PALETTE_H
