// carts/mandelbrot/colour.h
//
// COLOUR SPACES — sRGB, linear light, and OKLab.
//
// Gradients are mixed in OKLab (Bjorn Ottosson, 2020), a perceptual space in which equal steps
// look equal, so a ramp between two stops neither sags through a muddy middle nor rushes through
// its bright end. The conversions use only the IEEE-only functions of numerics.h, so every
// runtime computes the same palette to the bit.
#ifndef MANDELBROT_COLOUR_H
#define MANDELBROT_COLOUR_H

#include "jet64.h"

typedef struct {
	double l, a, b;  // lightness, green-red, blue-yellow
} OkLab;

/** An 8-bit sRGB colour in OKLab. */
OkLab colour_to_oklab(JetRgb srgb);

/** An OKLab colour in 8-bit sRGB, clipped to the sRGB gamut and rounded. */
JetRgb colour_from_oklab(OkLab lab);

/** The colour a fraction `t` of the way from `a` to `b`, in OKLab. */
OkLab colour_mix(OkLab a, OkLab b, double t);

#endif // MANDELBROT_COLOUR_H
