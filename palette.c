// carts/mandelbrot/palette.c
#pragma STDC FP_CONTRACT OFF

#include "palette.h"

#include "colour.h"
#include "config.h"
#include "numerics.h"

// The gradient is tabulated at this many samples per ring entry; the cycle's phase moves in
// whole table samples, so a sixteenth of an entry is the finest step colour ever takes.
#define PAL_GRADIENT_SUBSTEPS 16
#define PAL_GRADIENT_SAMPLES (PAL_RING_SIZE * PAL_GRADIENT_SUBSTEPS)

// One stop of a cyclic gradient: where it sits on the turn, in [0, 1), and its colour.
typedef struct {
	double position;
	JetRgb colour;
} GradientStop;

// The GRADIENTS, each a cyclic turn of stops, mixed in OKLab. PAL_GRADIENT_STOPS each.
#define PAL_GRADIENT_STOPS 5
typedef struct {
	GradientStop stops[PAL_GRADIENT_STOPS];
} Gradient;

static const Gradient GRADIENTS[PAL_GRADIENT_COUNT] = {
	// CLASSIC: the gradient generations of explorers grew up on — deep blue, through sky and white,
	// to amber and on into black, before the turn closes back onto the blue.
	{{{0.0, {0, 7, 100}}, {0.16, {32, 107, 203}}, {0.42, {237, 255, 255}}, {0.6425, {255, 170, 0}},
	  {0.8575, {0, 2, 0}}}},
	// EMBER: coals — a deep oxblood, through flame and gold to a white heat, and down through rust.
	{{{0.0, {26, 4, 10}}, {0.24, {150, 24, 14}}, {0.48, {240, 112, 22}}, {0.68, {255, 226, 150}},
	  {0.85, {118, 38, 20}}}},
	// GLACIER: ice — a deep navy, through teal and a pale cyan to ice-white, and back by slate violet.
	{{{0.0, {4, 14, 44}}, {0.28, {22, 96, 146}}, {0.52, {126, 206, 222}}, {0.7, {242, 250, 255}},
	  {0.86, {104, 112, 178}}}},
};
#define TURN 1.0  // one full turn of the gradient

// The HUD's inks, darkest to brightest (their indices are in palette.h).
static const JetRgb HUD_INKS[] = {
	{4, 5, 10},       // PAL_HUD_SHADOW
	{16, 18, 28},     // PAL_HUD_PANEL
	{128, 136, 152},  // PAL_HUD_DIM
	{232, 236, 244},  // PAL_HUD_TEXT
	{255, 196, 110},  // PAL_HUD_ACCENT
};
_Static_assert(sizeof(HUD_INKS) / sizeof(HUD_INKS[0]) == PAL_HUD_ACCENT - PAL_HUD_SHADOW + 1,
               "one ink per HUD index");

// What index 0 shows wherever the Copper is not compositing: the interior's own darkest tone.
static const JetRgb INTERIOR_FALLBACK = {MB_INTERIOR_BOTTOM_R, MB_INTERIOR_BOTTOM_G, MB_INTERIOR_BOTTOM_B};
// The Silk style's edge ink: a deep blue-black, darker than any colour of the gradient but its
// own darkest.
static const JetRgb INK = {8, 10, 24};

// The plane's hues sit this many gradient samples apart.
#define PAL_PLANE_HUE_SAMPLES (PAL_GRADIENT_SAMPLES / PAL_PLANE_HUES)
_Static_assert(PAL_GRADIENT_SAMPLES % PAL_PLANE_HUES == 0, "a plane hue sits on a gradient sample");

// Every gradient is tabulated at boot, so switching costs nothing: 3.8 k ring colours and 31 k plane
// colours each.
static JetRgb s_gradient[PAL_GRADIENT_COUNT][PAL_GRADIENT_SAMPLES];
static JetRgb s_plane[PAL_GRADIENT_COUNT][PAL_PLANE_SHADES][PAL_GRADIENT_SAMPLES];
static int s_current;  // the gradient shown
static PaletteLayout s_layout;
static int s_cycle_tick;  // ticks into the current turn, in [0, MB_CYCLE_PERIOD_TICKS)
static bool s_cycling = true;

// A gradient's colour a fraction `u` of the way round the turn, given its stops in OKLab.
static OkLab gradient_at(const Gradient* g, const OkLab* stops, double u) {
	int k = PAL_GRADIENT_STOPS - 1;
	while (k > 0 && g->stops[k].position > u) k--;
	const int next = (k + 1) % PAL_GRADIENT_STOPS;
	const double end = next == 0 ? TURN + g->stops[0].position : g->stops[next].position;
	const double t = (u - g->stops[k].position) / (end - g->stops[k].position);
	return colour_mix(stops[k], stops[next], t);
}

// A colour darkened toward black along its own ray in OKLab: lightness and chroma scaled together,
// so a shade keeps its hue.
static OkLab darken(OkLab c, double light) {
	c.l *= light;
	c.a *= light;
	c.b *= light;
	return c;
}

// The Copper paints the interior: one true-colour gradient row per scanline, top to bottom.
static void paint_copper(void) {
	const int last = MB_SCREEN_HEIGHT - 1;
	for (int y = 0; y < MB_SCREEN_HEIGHT; y++) {
		const int below = last - y;
		jet_copper_set_rgb((uint32_t)y,
		                   (uint8_t)((MB_INTERIOR_TOP_R * below + MB_INTERIOR_BOTTOM_R * y) / last),
		                   (uint8_t)((MB_INTERIOR_TOP_G * below + MB_INTERIOR_BOTTOM_G * y) / last),
		                   (uint8_t)((MB_INTERIOR_TOP_B * below + MB_INTERIOR_BOTTOM_B * y) / last));
	}
	jet_copper_enable(true);
}

void palette_init(void) {
	for (int g = 0; g < PAL_GRADIENT_COUNT; g++) {
		OkLab stops[PAL_GRADIENT_STOPS];
		for (int k = 0; k < PAL_GRADIENT_STOPS; k++) stops[k] = colour_to_oklab(GRADIENTS[g].stops[k].colour);
		for (int i = 0; i < PAL_GRADIENT_SAMPLES; i++) {
			const OkLab c = gradient_at(&GRADIENTS[g], stops, (double)i / (double)PAL_GRADIENT_SAMPLES);
			s_gradient[g][i] = colour_from_oklab(c);
			for (int shade = 0; shade < PAL_PLANE_SHADES; shade++) {
				s_plane[g][shade][i] = colour_from_oklab(darken(c, palette_shade_light((double)shade)));
			}
		}
	}
	s_current = 0;
	for (int i = 0; i <= PAL_HUD_ACCENT - PAL_HUD_SHADOW; i++) {
		jet_palette_set((uint8_t)(PAL_HUD_SHADOW + i), HUD_INKS[i].r, HUD_INKS[i].g, HUD_INKS[i].b);
	}
	jet_palette_set(PAL_INK, INK.r, INK.g, INK.b);
	jet_palette_set(PAL_INTERIOR, INTERIOR_FALLBACK.r, INTERIOR_FALLBACK.g, INTERIOR_FALLBACK.b);
	paint_copper();
	s_layout = PAL_LAYOUT_PLANE;
	s_cycle_tick = 0;
	s_cycling = true;
	palette_set_title_fade(1.0);
	palette_set_status_fade(1.0);
	palette_commit();
}

void palette_set_layout(PaletteLayout layout) {
	s_layout = layout;
}

void palette_next_gradient(void) {
	s_current = (s_current + 1) % PAL_GRADIENT_COUNT;
}

// A fading ink: the text's colour at `share`, down to the outline's shadow at none, where the glyph
// disappears into its own outline.
static void set_fade(uint8_t ink, double share) {
	const JetRgb full = HUD_INKS[PAL_HUD_TEXT - PAL_HUD_SHADOW];
	const JetRgb dark = HUD_INKS[0];
	const double t = num_clamp(share, 0.0, 1.0);
	jet_palette_set(ink, (uint8_t)(dark.r + (full.r - dark.r) * t + 0.5), (uint8_t)(dark.g + (full.g - dark.g) * t + 0.5),
	                (uint8_t)(dark.b + (full.b - dark.b) * t + 0.5));
}

void palette_set_title_fade(double share) {
	set_fade(PAL_TITLE, share);
}

void palette_set_status_fade(double share) {
	set_fade(PAL_STATUS, share);
}

void palette_advance(void) {
	if (!s_cycling) return;
	s_cycle_tick = (s_cycle_tick + 1) % MB_CYCLE_PERIOD_TICKS;
}

void palette_set_cycling(bool on) {
	s_cycling = on;
}

bool palette_cycling(void) {
	return s_cycling;
}

void palette_commit(void) {
	// The phase, in table samples. Adding it to every entry's own place on the ring moves each
	// colour toward lower entries — toward lower mu, which is outward, away from the set.
	const int phase = (int)((int64_t)s_cycle_tick * PAL_GRADIENT_SAMPLES / MB_CYCLE_PERIOD_TICKS);
	if (s_layout == PAL_LAYOUT_RING) {
		for (int e = 0; e < PAL_RING_SIZE; e++) {
			const JetRgb c = s_gradient[s_current][(e * PAL_GRADIENT_SUBSTEPS + phase) % PAL_GRADIENT_SAMPLES];
			jet_palette_set((uint8_t)(PAL_RING_FIRST + e), c.r, c.g, c.b);
		}
		return;
	}
	for (int shade = 0; shade < PAL_PLANE_SHADES; shade++) {
		for (int h = 0; h < PAL_PLANE_HUES; h++) {
			const JetRgb c = s_plane[s_current][shade][(h * PAL_PLANE_HUE_SAMPLES + phase) % PAL_GRADIENT_SAMPLES];
			jet_palette_set((uint8_t)(PAL_RING_FIRST + shade * PAL_PLANE_HUES + h), c.r, c.g, c.b);
		}
	}
}
