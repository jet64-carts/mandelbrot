// carts/mandelbrot/colour.c
//
// The matrices are Ottosson's published OKLab transform (https://bottosson.github.io/posts/oklab/)
// and the transfer function is IEC 61966-2-1's sRGB curve.
#pragma STDC FP_CONTRACT OFF

#include "colour.h"

#include "numerics.h"

// --- The sRGB transfer function --------------------------------------------------------------
#define CHANNEL_MAX 255.0
#define CHANNEL_ROUND 0.5
#define SRGB_LINEAR_KNEE 0.04045     // encoded values at or below this are on the linear toe
#define SRGB_LINEAR_KNEE_OUT 0.0031308  // the same knee, in linear light
#define SRGB_TOE_SLOPE 12.92
#define SRGB_OFFSET 0.055
#define SRGB_SCALE 1.055
#define SRGB_GAMMA 2.4

static double srgb_to_linear(uint8_t channel) {
	const double c = (double)channel / CHANNEL_MAX;
	if (c <= SRGB_LINEAR_KNEE) return c / SRGB_TOE_SLOPE;
	return num_pow((c + SRGB_OFFSET) / SRGB_SCALE, SRGB_GAMMA);
}

static uint8_t linear_to_srgb(double v) {
	if (!(v > 0.0)) return 0;  // NaN lands here too
	if (v >= 1.0) return (uint8_t)CHANNEL_MAX;
	const double c = v <= SRGB_LINEAR_KNEE_OUT ? v * SRGB_TOE_SLOPE
	                                            : SRGB_SCALE * num_pow(v, 1.0 / SRGB_GAMMA) - SRGB_OFFSET;
	return (uint8_t)(c * CHANNEL_MAX + CHANNEL_ROUND);
}

// --- OKLab -----------------------------------------------------------------------------------
// Linear sRGB -> cone response (LMS), then a cube root and a second matrix to OKLab; the inverse
// runs the same road backwards. Each matrix is written out as rows of coefficients.
#define CHANNELS 3
static const double TO_LMS[CHANNELS][CHANNELS] = {
	{0.4122214708, 0.5363325363, 0.0514459929},
	{0.2119034982, 0.6806995451, 0.1073969566},
	{0.0883024619, 0.2817188376, 0.6299787005},
};
static const double LMS_TO_LAB[CHANNELS][CHANNELS] = {
	{0.2104542553, 0.7936177850, -0.0040720468},
	{1.9779984951, -2.4285922050, 0.4505937099},
	{0.0259040371, 0.7827717662, -0.8086757660},
};
static const double LAB_TO_LMS[CHANNELS][CHANNELS] = {
	{1.0, 0.3963377774, 0.2158037573},
	{1.0, -0.1055613458, -0.0638541728},
	{1.0, -0.0894841775, -1.2914855480},
};
static const double FROM_LMS[CHANNELS][CHANNELS] = {
	{4.0767416621, -3.3077115913, 0.2309699292},
	{-1.2684380046, 2.6097574011, -0.3413193965},
	{-0.0041960863, -0.7034186147, 1.7076147010},
};

static void multiply(const double m[CHANNELS][CHANNELS], const double in[CHANNELS], double out[CHANNELS]) {
	for (int row = 0; row < CHANNELS; row++) {
		out[row] = m[row][0] * in[0] + m[row][1] * in[1] + m[row][2] * in[2];
	}
}

OkLab colour_to_oklab(JetRgb srgb) {
	const double linear[CHANNELS] = {srgb_to_linear(srgb.r), srgb_to_linear(srgb.g), srgb_to_linear(srgb.b)};
	double lms[CHANNELS], cone[CHANNELS], lab[CHANNELS];
	multiply(TO_LMS, linear, lms);
	for (int i = 0; i < CHANNELS; i++) cone[i] = num_cbrt(lms[i]);
	multiply(LMS_TO_LAB, cone, lab);
	OkLab out = {lab[0], lab[1], lab[2]};
	return out;
}

JetRgb colour_from_oklab(OkLab lab) {
	const double in[CHANNELS] = {lab.l, lab.a, lab.b};
	double cone[CHANNELS], lms[CHANNELS], linear[CHANNELS];
	multiply(LAB_TO_LMS, in, cone);
	for (int i = 0; i < CHANNELS; i++) lms[i] = cone[i] * cone[i] * cone[i];
	multiply(FROM_LMS, lms, linear);
	JetRgb out = {linear_to_srgb(linear[0]), linear_to_srgb(linear[1]), linear_to_srgb(linear[2])};
	return out;
}

OkLab colour_mix(OkLab a, OkLab b, double t) {
	OkLab out = {a.l + (b.l - a.l) * t, a.a + (b.a - a.a) * t, a.b + (b.b - a.b) * t};
	return out;
}
