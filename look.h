// carts/mandelbrot/look.h
//
// THE LOOK — the picture the player sees, shaded from the field once per changed frame.
//
// The field holds each pixel's smooth iteration count mu at its actual coordinate. The look turns
// that into a palette index with three things a pixel cannot know alone, so it is a pass over the
// whole screen, run on the pool between the settle and the present (design.md §8.6):
//
//   * THE COLOUR SCALE. Hue follows log(mu - base + 1) between robust percentiles of the view's
//     own exterior, spanning LOOK_TURNS turns of the gradient: every view, home or 10^-254, gets
//     the whole gradient across what it shows, and a common offset in mu (about 40 000 at the
//     deepest bookmark) is subtracted before the logarithm, so it costs no precision. The scale
//     eases toward each frame's percentiles, so colours flow while the view moves; a jump snaps it.
//   * DISTANCE-ESTIMATE EDGES. mu = C - log2 G, so the distance to the set in pixels is
//     DE = 1 / (ln 2 |grad mu|), grad mu taken by differences over the grid's actual coordinates.
//     A pixel within LOOK_EDGE_WIDTH of the set darkens by its coverage — which draws every
//     filament as a line where colouring by mu alone is aliased noise. Where the boundary is dense
//     everywhere, the darkening's depth is held back by the share of such pixels, so a view like
//     the Feigenbaum point keeps its colours.
//   * RELIEF. The gradient's direction lights the field as a surface, from the upper left.
//
// Two styles share the pass: RELIEF (the default) indexes a plane of PAL_PLANE_HUES hues by
// PAL_PLANE_SHADES shades; SILK the full ring of PAL_RING_SIZE hues, its edges dithered toward one
// ink. The plane's two axes are dithered by thresholds from two unrelated families, so they never
// lock into one lattice: the shade by Roberts' R2, the hue by interleaved gradient noise.
//
// The pass also counts the frame's READY pixels — those at their own coordinates, reused or better —
// which the camera reads to keep a flight from outrunning the picture (flight.h).
#ifndef MANDELBROT_LOOK_H
#define MANDELBROT_LOOK_H

#include "field.h"

typedef enum {
	LOOK_RELIEF,
	LOOK_SILK,
	LOOK_STYLE_COUNT
} LookStyle;

// The colour scale: hue in turns = (log2(max(mu - base, 0) + 1) - low) * gain; and the edges'
// depth, the lightness an edge pixel keeps at the least.
typedef struct {
	double base;
	double low;
	double gain;
	double edge_floor;
} LookScale;

/** Bytes of WRAM look_init reserves (the shaded frame). */
uint32_t look_bytes(void);

/** Reserves the shaded frame. False if WRAM ran out. */
bool look_init(void);

/** The style the pass shades in; a change reshades the next frame. */
void look_set_style(LookStyle style);
LookStyle look_style(void);

/**
 * The player's colour density, in quarter-octaves of the scale's gain (0 is the automatic scale
 * as measured); a change holds the automatic scale where it is for LOOK_HOLD_TICKS.
 */
void look_adjust_density(int quarter_octaves);

/**
 * The scale a picture's samples ask for: base the least exterior mu sampled, low and the gain from
 * the LOOK_LOW_PERCENTILE and LOOK_HIGH_PERCENTILE of log2(mu - base + 1), and the edge floor from
 * the share of samples within an edge width of the set. False (the scale untouched) when the
 * picture has too few exterior samples to say.
 */
bool look_measure(const Field* f, LookScale* scale);

/**
 * Brings the shaded frame up to date: measures the scale when the field `changed`, eases toward it
 * (or snaps, on a `jump`), and reshades the whole frame on the pool when anything it depends on
 * moved. Call from _draw between the settle and the present, with no batch in flight.
 */
void look_update(const Field* f, bool changed, bool jump);

/** Shades the whole frame of `f` at `scale` into `out`, on the pool (the pass look_update runs). */
void look_shade(const Field* f, const LookScale* scale, LookStyle style, uint8_t* out);

/** The shaded frame. */
const uint8_t* look_frame(void);

/** The share of the shaded frame's pixels at their own coordinates (reused or better), in [0, 1]. */
double look_ready_share(void);

#endif // MANDELBROT_LOOK_H
