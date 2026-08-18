// carts/mandelbrot/main.c
//
// MANDELBROT — the Mandelbrot set, explored at 120 Hz.
//
// Drag it, wheel into it, fly through it: every frame is computed on all the machine's cores by
// a speculative SIMD kernel that finds each pixel's exact escape, refined coarse-to-fine so a
// new view is on screen in the tick it is asked for, and coloured by a smooth, dithered, cycling
// gradient. design.md in docs/mandelbrot tells the whole story.
//
// The hooks delegate and nothing else.

#include "jet64.h"

#include "explorer.h"

#define FAULT_FONT 0
#define FAULT_INK 15
#define FAULT_SIZE 8
#define FAULT_MARGIN 8

static bool s_ok;

void _init(void) {
	s_ok = explorer_init();
}

void _update(float dt) {
	(void)dt;  // the explorer runs on the machine's fixed tick; see camera.c
	if (s_ok) explorer_update();
}

void _draw(void) {
	if (s_ok) {
		explorer_draw();
		return;
	}
	jet_cls(0);
	TextFormat fmt = jet_make_text_format(FAULT_FONT, FAULT_INK, FAULT_SIZE, 0, 0, 1.0f, 0, 0);
	jet_print_ext("MANDELBROT: WORKING MEMORY RAN OUT", FAULT_MARGIN, FAULT_MARGIN, &fmt);
}

void _end(void) {
	if (s_ok) explorer_end();
}
