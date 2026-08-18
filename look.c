// carts/mandelbrot/look.c
#pragma STDC FP_CONTRACT OFF

#include "look.h"

#include "config.h"
#include "jet_thread.h"
#include "numerics.h"
#include "palette.h"

// --- The colour scale ------------------------------------------------------------------------------
#define LOOK_TURNS 3.0             // gradient turns across the view's robust range of log mu
#define LOOK_LOW_PERCENTILE 0.02   // the range: these percentiles of the samples' log2(mu - base + 1)
#define LOOK_HIGH_PERCENTILE 0.98
#define LOOK_MIN_SPAN 0.25         // log2 units: a nearly uniform view is not stretched without end
#define LOOK_SAMPLE_STRIDE 8       // the scale samples every eighth pixel of every eighth row
#define LOOK_SAMPLE_CAPACITY ((MB_SCREEN_WIDTH / LOOK_SAMPLE_STRIDE) * (MB_SCREEN_HEIGHT / LOOK_SAMPLE_STRIDE))
#define LOOK_MIN_SAMPLES 64        // fewer exterior samples than this say nothing about a scale
#define LOOK_EASE 0.08             // the share of the way to the measured scale the shown one moves a tick
#define LOOK_SETTLED 1e-4          // closer than this (in turns, or relatively) and the scale arrives
#define LOOK_HOLD_TICKS 90         // after the player sets the density, the measured scale holds
#define LOOK_DENSITY_STEPS 4.0     // density steps an octave of gain
#define LOOK_DENSITY_LIMIT 16      // steps either way

// --- Edges and relief ------------------------------------------------------------------------------
#define LOOK_LN2 0.69314718055994530942
#define LOOK_EDGE_WIDTH 0.7        // pixels: within this distance of the set a pixel darkens by coverage
#define LOOK_FLOOR_CRISP 0.2       // the lightness an edge keeps where edges are few
#define LOOK_FLOOR_DENSE 0.8       // ... and where the boundary fills the view
#define LOOK_DENSE_START 0.45      // the share of edge pixels at which the floor starts to rise
#define LOOK_DENSE_FULL 0.75       // ... and at which it is all the way dense
#define LOOK_MIN_STEP 0.25         // pixels: a neighbour nearer than this is a copy, not a sample
#define LOOK_RELIEF 0.9            // the surface's slope, for lighting
#define LOOK_AMBIENT 0.55          // the lightness of a slope turned fully from the light
#define LOOK_LIGHT_X (-0.5)        // the light, from the upper left at 45 degrees (screen y is down)
#define LOOK_LIGHT_Y 0.5
#define LOOK_LIGHT_Z 0.70710678118654752440

// --- The pass ----------------------------------------------------------------------------------------
#define LOOK_BAND_ROWS 16
#define LOOK_BANDS (MB_SCREEN_HEIGHT / LOOK_BAND_ROWS)
_Static_assert(MB_SCREEN_HEIGHT % LOOK_BAND_ROWS == 0, "the bands tile the screen");

// The two axes' dither thresholds come from two unrelated families, so they never lock into one
// lattice: the SHADE from Roberts' R2 (below); the HUE from interleaved gradient noise (Jimenez,
// 2014), frac(k3 frac(k1 x + k2 y)), whose fine diagonal structure the eye does not assemble into
// strokes — against R2's transpose, measured on a 4x crop of the home view's bright band, it removed
// a brick-like texture that showed where hues change lightness fastest.
//
// A level almost never falls on an entry, so it is DITHERED: a per-pixel threshold t in [0, 1) is
// added before truncation, so a pixel a third of the way from one entry to the next shows the next
// entry with probability one third, and the eye averages the two. R2 is frac(x a1 + y a2), with a1,
// a2 the reciprocals of the plastic number and its square: a low-discrepancy pattern with no visible
// period, computed exactly in integer arithmetic, needing no table and no initialisation. Both
// families depend only on the pixel's position on the screen, so a still picture never shimmers, and
// a moving one is reshaded against the same fixed pattern rather than dragging its dither along.
//
// R2 in 32-bit fixed point: round(2^32 / p) and round(2^32 / p^2), with p the plastic number (the
// real root of x^3 = x + 1). Unsigned arithmetic wraps by definition, which is exactly the
// "fractional part" the sequence needs; LOOK_R2_SCALE maps the fraction onto [0, 1).
#define LOOK_R2_STEP_X 3242174889u
#define LOOK_R2_STEP_Y 2447445414u
#define LOOK_R2_SCALE 0x1p-32
#define LOOK_IGN_X 0.06711056
#define LOOK_IGN_Y 0.00583715
#define LOOK_IGN_SCALE 52.9829189

static double fraction_of(double v) {  // v >= 0
	return v - (double)(int64_t)v;
}

static double hue_threshold_at(int x, int y) {
	return fraction_of(LOOK_IGN_SCALE * fraction_of(LOOK_IGN_X * x + LOOK_IGN_Y * y));
}

// One axis of the grid seen from each pixel: the reciprocal of the distance, in screen pixels, from
// its neighbour before to its neighbour after (central), from itself to the one after, and from the
// one before to itself; 0 where there is no such neighbour, or it sits too near to be a sample.
typedef struct {
	double central[MB_SCREEN_WIDTH];
	double after[MB_SCREEN_WIDTH];
	double before[MB_SCREEN_WIDTH];
} LookRuns;

typedef struct {
	const Field* field;
	LookScale scale;
	LookStyle style;
	uint8_t* out;
	const LookRuns* columns;
	const LookRuns* rows;
} LookJob;

typedef struct {
	int first_row;
} LookBand;

static uint8_t* s_frame;
static LookStyle s_style;
static bool s_style_changed;
static int s_density;       // the player's density, in steps
static int s_hold;          // ticks the measured scale still holds
static LookScale s_target;  // the scale the picture asks for
static LookScale s_shown;   // the scale being shown, easing toward the target
static LookScale s_shaded;  // the scale the frame was last shaded at (the density applied)
static bool s_has_target, s_has_shown, s_has_frame;
static double s_samples[LOOK_SAMPLE_CAPACITY];
static LookJob s_job;
static LookBand s_bands[LOOK_BANDS];
// Each band's three rolling rows of mu — the rows either side of the one whose shading terms are
// being measured, and that row itself — with LOOK_UNPLACED where a pixel is the set's (interior) or
// off the screen, and LOOK_PREVIEW where its value is not of its own coordinate (look.c: placed): a
// pixel's slope reads its four neighbours' mu once each from here rather than from the field's three
// arrays. And its three rolling rows of SHADING TERMS — the share of full lightness each pixel's edge
// keeps, and its relief's light — for the row being shaded and the rows either side, LOOK_UNKNOWN
// where no slope could be measured: such a pixel borrows its neighbours' terms from here.
#define LOOK_UNPLACED (-1.0)  // mu >= 0 for every exterior pixel (n >= 1, phi <= 1)
#define LOOK_PREVIEW (-2.0)
#define LOOK_UNKNOWN (-1.0)   // an edge term is a lightness in [0, 1]
// A placed pixel's slope is trusted while at most this many of its four neighbours are previews.
#define LOOK_MAX_PREVIEW_NEIGHBOURS 1
#define LOOK_RING 3
static double s_rows[LOOK_BANDS][LOOK_RING][MB_SCREEN_WIDTH];
static double s_edges[LOOK_BANDS][LOOK_RING][MB_SCREEN_WIDTH];
static double s_lights[LOOK_BANDS][LOOK_RING][MB_SCREEN_WIDTH];
// And each row's terms summed across each pixel and its two neighbours in the row, with how many
// there are, made at most once a row and only when a row it neighbours has pixels to borrow (so a
// converged picture never makes them): a pixel that borrows reads three of each, one from each row
// of its window, rather than visiting its eight neighbours. In floats, which only a preview's
// borrowed terms are made of, so the sums take half the time and touch half the memory.
static float s_sum_edges[LOOK_BANDS][LOOK_RING][MB_SCREEN_WIDTH];
static float s_sum_lights[LOOK_BANDS][LOOK_RING][MB_SCREEN_WIDTH];
static float s_sum_known[LOOK_BANDS][LOOK_RING][MB_SCREEN_WIDTH];
static int32_t s_band_ready[LOOK_BANDS];  // each band's pixels at their own coordinates
static double s_ready_share = 1.0;
static LookRuns s_column_runs, s_row_runs;

uint32_t look_bytes(void) {
	return MB_SCREEN_PIXELS;
}

bool look_init(void) {
	s_frame = (uint8_t*)jet_memory_reserve(look_bytes());
	if (!s_frame) return false;
	memset(s_frame, PAL_INTERIOR, look_bytes());
	s_style = LOOK_RELIEF;
	s_style_changed = true;
	s_density = 0;
	s_hold = 0;
	s_has_target = s_has_shown = s_has_frame = false;
	for (int b = 0; b < LOOK_BANDS; b++) s_bands[b].first_row = b * LOOK_BAND_ROWS;
	return true;
}

void look_set_style(LookStyle style) {
	s_style_changed = s_style_changed || style != s_style;
	s_style = style;
}

LookStyle look_style(void) {
	return s_style;
}

void look_adjust_density(int steps) {
	s_density += steps;
	if (s_density > LOOK_DENSITY_LIMIT) s_density = LOOK_DENSITY_LIMIT;
	if (s_density < -LOOK_DENSITY_LIMIT) s_density = -LOOK_DENSITY_LIMIT;
	s_hold = LOOK_HOLD_TICKS;
}

const uint8_t* look_frame(void) {
	return s_frame;
}

// --- The gradient of mu --------------------------------------------------------------------------
// A pixel whose value belongs to its own coordinate: exterior, and reused or better (a preview is a
// neighbour's copy, whose mu belongs somewhere else).
static inline bool placed(const PixelStore* s, size_t p) {
	const uint8_t st = s->state[p];
	return !state_interior(st) && state_rank(st) >= RANK_REUSED;
}

// The reciprocal of a distance between two coordinates in screen pixels, or 0 when it is too short.
static double reciprocal_run(double run, double inv_spacing, double sign) {
	const double pixels = sign * run * inv_spacing;
	return pixels > LOOK_MIN_STEP ? 1.0 / pixels : 0.0;
}

// An axis's runs from its actual coordinates; `sign` is +1 when they grow along the screen axis and
// -1 when they shrink (the rows).
static void measure_runs(const double* coordinate, int count, double inv_spacing, double sign, LookRuns* runs) {
	for (int i = 0; i < count; i++) {
		runs->central[i] =
		    i > 0 && i + 1 < count ? reciprocal_run(coordinate[i + 1] - coordinate[i - 1], inv_spacing, sign) : 0.0;
		runs->after[i] = i + 1 < count ? reciprocal_run(coordinate[i + 1] - coordinate[i], inv_spacing, sign) : 0.0;
		runs->before[i] = i > 0 ? reciprocal_run(coordinate[i] - coordinate[i - 1], inv_spacing, sign) : 0.0;
	}
}

// d mu along one axis from the neighbours before and after a pixel: central where both are placed,
// one-sided where one is. False when neither is placed, or they sit too near to be samples.
static bool axis_slope(const PixelStore* s, const LookRuns* runs, int i, int count, size_t p, size_t stride,
                       double mu, double* slope) {
	const bool before = i > 0 && placed(s, p - stride);
	const bool after = i + 1 < count && placed(s, p + stride);
	double rise, per_pixel;
	if (before && after) {
		rise = field_mu(s, p + stride) - field_mu(s, p - stride);
		per_pixel = runs->central[i];
	} else if (after) {
		rise = field_mu(s, p + stride) - mu;
		per_pixel = runs->after[i];
	} else if (before) {
		rise = mu - field_mu(s, p - stride);
		per_pixel = runs->before[i];
	} else {
		return false;
	}
	if (!(per_pixel > 0.0)) return false;
	*slope = rise * per_pixel;
	return true;
}

// grad mu at pixel (x, y) in mu per screen pixel, x rightward and y downward. False when neither
// axis has a slope.
static bool gradient(const Field* f, const LookRuns* columns, const LookRuns* rows, int x, int y, double mu,
                     double* gx, double* gy) {
	const size_t p = (size_t)y * MB_SCREEN_WIDTH + (size_t)x;
	*gx = *gy = 0.0;
	const bool across = axis_slope(&f->front, columns, x, MB_SCREEN_WIDTH, p, 1, mu, gx);
	const bool down = axis_slope(&f->front, rows, y, MB_SCREEN_HEIGHT, p, MB_SCREEN_WIDTH, mu, gy);
	return across || down;
}

// Both axes' runs for the field's grid as it stands.
static void measure_grid(const Field* f) {
	const double inv_spacing = 1.0 / view_spacing(&f->spec.view);
	measure_runs(f->columns.actual, MB_SCREEN_WIDTH, inv_spacing, 1.0, &s_column_runs);
	measure_runs(f->rows.actual, MB_SCREEN_HEIGHT, inv_spacing, -1.0, &s_row_runs);
}

// The share of full lightness a pixel keeps for its distance to the set: 1 beyond an edge width,
// falling smoothly to `floor` at the set itself.
#define LOOK_EDGE_SLOPE (1.0 / (LOOK_LN2 * LOOK_EDGE_WIDTH))  // the slope at a distance of one width
static double edge_light(double slope, double floor) {
	if (!(slope > LOOK_EDGE_SLOPE)) return 1.0;
	const double t = LOOK_EDGE_SLOPE / slope;  // the distance, in edge widths
	if (t >= 1.0) return 1.0;
	return floor + (1.0 - floor) * t * t * (3.0 - 2.0 * t);
}

// THE RELIEF: the surface lit from the upper left, its normal tilted along the gradient's direction
// by a constant slope LOOK_RELIEF; flat ground reads 1, a slope turned from the light down to
// LOOK_AMBIENT (the pass computes it inline, measure_terms).

// log2 x for x >= 1, to within 2e-7: the exponent from the bits, and log2 of the mantissa m in
// [1/sqrt 2, sqrt 2) from the series 2/ln 2 (t + t^3/3 + t^5/5 + t^7/7 + t^9/9), t = (m - 1)/(m + 1),
// |t| < 0.172, whose first omitted term is below 2e-9. Hue needs a thousandth of a turn at most;
// num_log2, to an ulp, was most of the pass's time. IEEE arithmetic only, so every runtime shades
// alike.
#define LOOK_EXPONENT_BIAS 1023
#define LOOK_MANTISSA_BITS 52
#define LOOK_MANTISSA_MASK (((uint64_t)1 << LOOK_MANTISSA_BITS) - 1)
#define LOOK_ONE_BITS ((uint64_t)LOOK_EXPONENT_BIAS << LOOK_MANTISSA_BITS)
#define LOOK_SQRT2 1.41421356237309504880
#define LOOK_TWO_OVER_LN2 2.88539008177792681472
static double look_log2(double x) {
	uint64_t bits;
	memcpy(&bits, &x, sizeof bits);
	int exponent = (int)(bits >> LOOK_MANTISSA_BITS) - LOOK_EXPONENT_BIAS;
	bits = (bits & LOOK_MANTISSA_MASK) | LOOK_ONE_BITS;
	double m;
	memcpy(&m, &bits, sizeof m);
	if (m >= LOOK_SQRT2) {
		m *= 0.5;
		exponent++;
	}
	const double t = (m - 1.0) / (m + 1.0), t2 = t * t;
	const double series = t * (1.0 + t2 * (1.0 / 3.0 + t2 * (1.0 / 5.0 + t2 * (1.0 / 7.0 + t2 * (1.0 / 9.0)))));
	return (double)exponent + LOOK_TWO_OVER_LN2 * series;
}

// The hue of mu on the scale, as a fraction of a turn in [0, 1).
static double hue_of(const LookScale* scale, double mu) {
	const double above = mu - scale->base;
	const double turns = (look_log2((above > 0.0 ? above : 0.0) + 1.0) - scale->low) * scale->gain;
	const double whole = (double)(int64_t)turns;  // the scale keeps turns far inside int64's range
	const double fraction = turns - whole;
	return fraction < 0.0 ? fraction + 1.0 : fraction;
}

// --- The pass ------------------------------------------------------------------------------------
// Row y's mu into `row`: LOOK_PREVIEW where a pixel's value is not of its own coordinate (a preview,
// interior or not), LOOK_UNPLACED where it is the set's; a row off the screen is all unplaced.
static void fill_row(const PixelStore* s, int y, double* row) {
	if (y < 0 || y >= MB_SCREEN_HEIGHT) {
		for (int x = 0; x < MB_SCREEN_WIDTH; x++) row[x] = LOOK_UNPLACED;
		return;
	}
	const size_t start = (size_t)y * MB_SCREEN_WIDTH;
	for (int x = 0; x < MB_SCREEN_WIDTH; x++) {
		const size_t p = start + (size_t)x;
		row[x] = placed(s, p) ? field_mu(s, p) : state_rank(s->state[p]) <= RANK_PREVIEW ? LOOK_PREVIEW : LOOK_UNPLACED;
	}
}

// True when an axis's slope is unknown rather than flat: neither neighbour placed, and one of them a
// preview, whose value says nothing of the slope here. (Beside the set's interior, or the screen's
// edge, the axis is flat, as it always was.)
static bool axis_unknown(double before, double after) {
	return !(before >= 0.0) && !(after >= 0.0) && (before == LOOK_PREVIEW || after == LOOK_PREVIEW);
}

// d mu along one axis from the cached neighbours: central where both are placed, one-sided where one
// is; 0 where neither is, or they sit too near to be samples (their run's reciprocal is 0).
static double cached_slope(double before, double centre, double after, const LookRuns* runs, int i) {
	if (before >= 0.0 && after >= 0.0) return (after - before) * runs->central[i];
	if (after >= 0.0) return (after - centre) * runs->after[i];
	if (before >= 0.0) return (centre - before) * runs->before[i];
	return 0.0;
}

// The relief's light for a unit gradient direction (ux, uy): the surface's normal is
// (-ux R, uy R, 1), whose length is sqrt(R^2 + 1) whatever the direction, so the diffuse term is
// one dot product with these constants.
#define LOOK_NORMAL_LENGTH 1.34536240470737125  // sqrt(LOOK_RELIEF^2 + 1)
#define LOOK_DIFFUSE_X (-LOOK_RELIEF * LOOK_LIGHT_X / LOOK_NORMAL_LENGTH)
#define LOOK_DIFFUSE_Y (LOOK_RELIEF * LOOK_LIGHT_Y / LOOK_NORMAL_LENGTH)
#define LOOK_DIFFUSE_Z (LOOK_LIGHT_Z / LOOK_NORMAL_LENGTH)

// The light of a surface whose slope is unknown: the relief's mean over every direction it might
// face. The diffuse term is DIFFUSE_Z plus a dot product of the direction with (DIFFUSE_X, DIFFUSE_Y),
// which averages to nothing over the circle, and never clips at zero (DIFFUSE_Z, 0.53, exceeds the
// length of that vector, 0.47), so the mean is the light of DIFFUSE_Z alone — 0.88 of full, where
// flat ground reads 1. A preview with no measured neighbour, lit as flat ground, read as a light
// patch among relief-lit samples.
#define LOOK_UNKNOWN_LIGHT (LOOK_AMBIENT + (1.0 - LOOK_AMBIENT) * LOOK_DIFFUSE_Z / LOOK_LIGHT_Z)

// One row's rolling buffers: mu, the two shading terms, and their sums across each pixel's
// neighbourhood in the row (s_sum_*), with whether those are made yet; and how many of its pixels
// have no terms of their own (previews, interior ones among them, and samples whose slope is not
// trusted).
typedef struct {
	double* mu;
	double* edge;
	double* light;
	float* sum_edge;
	float* sum_light;
	float* sum_known;
	bool summed;
	int32_t borrowers;
} LookRow;

// Measures row y's shading terms from the mu of the rows above it, at it and below it: for every
// placed pixel whose slope is TRUSTED, the share of full lightness its edge keeps (falling to `floor`
// at the set itself) and its relief's light — flat ground, 1 and 1, where the slope is zero;
// LOOK_UNKNOWN for every other pixel. A slope is trusted where it is known along both axes and at
// most LOOK_MAX_PREVIEW_NEIGHBOURS of the four neighbours are previews. Known along one axis alone,
// a slope tilts the relief onto that axis; and a sample with a lone neighbour on each axis, among
// copies, extrapolates a slope at right angles from two values. In a moving frame at depth, where
// most of the picture is previews, either lit its pixel (and the copies borrowing from it) by a coin
// toss: a speckle of light and dark across the copies' even ground. A converged picture has no
// previews, so this changes none of its pixels.
// Returns how many of the row's pixels, the set's interior and the screen aside, have no terms.
static int32_t measure_terms(const LookJob* job, const double* above, const double* here, const double* below, int y,
                             double floor, double* edge_out, double* light_out) {
	int32_t borrowers = 0;
	for (int x = 0; x < MB_SCREEN_WIDTH; x++) {
		edge_out[x] = LOOK_UNKNOWN;
		borrowers += here[x] != LOOK_UNPLACED;  // counted back when the pixel has terms after all
		if (!(here[x] >= 0.0)) continue;
		const double left = x > 0 ? here[x - 1] : LOOK_UNPLACED;
		const double right = x + 1 < MB_SCREEN_WIDTH ? here[x + 1] : LOOK_UNPLACED;
		if (axis_unknown(left, right) || axis_unknown(above[x], below[x])) continue;
		const int previews = (left == LOOK_PREVIEW) + (right == LOOK_PREVIEW) + (above[x] == LOOK_PREVIEW) +
		                     (below[x] == LOOK_PREVIEW);
		if (previews > LOOK_MAX_PREVIEW_NEIGHBOURS) continue;
		const double gx = cached_slope(left, here[x], right, job->columns, x);
		const double gy = cached_slope(above[x], here[x], below[x], job->rows, y);
		const double slope = sqrt(gx * gx + gy * gy);
		edge_out[x] = light_out[x] = 1.0;
		borrowers--;
		if (!(slope > 0.0)) continue;
		const double inverse = 1.0 / slope;
		double edge = 1.0;
		if (slope > LOOK_EDGE_SLOPE) {
			const double t = LOOK_EDGE_SLOPE * inverse;  // the distance, in edge widths
			edge = floor + (1.0 - floor) * t * t * (3.0 - 2.0 * t);
		}
		double diffuse = LOOK_DIFFUSE_X * gx * inverse + LOOK_DIFFUSE_Y * gy * inverse + LOOK_DIFFUSE_Z;
		if (diffuse < 0.0) diffuse = 0.0;
		double light = LOOK_AMBIENT + (1.0 - LOOK_AMBIENT) * diffuse / LOOK_LIGHT_Z;
		if (light > 1.0) light = 1.0;
		edge_out[x] = edge;
		light_out[x] = light;
	}
	return borrowers;
}

// A measured row's terms summed across each pixel and its two neighbours in the row, branch-free: a
// pixel's measured terms count, an unknown one's nothing.
static void sum_row(LookRow* row) {
	float edge[MB_SCREEN_WIDTH], light[MB_SCREEN_WIDTH], known[MB_SCREEN_WIDTH];
	for (int x = 0; x < MB_SCREEN_WIDTH; x++) {
		const bool has = row->edge[x] >= 0.0;
		edge[x] = has ? (float)row->edge[x] : 0.0f;
		light[x] = has ? (float)row->light[x] : 0.0f;
		known[x] = has ? 1.0f : 0.0f;
	}
	const int last = MB_SCREEN_WIDTH - 1;
	row->sum_edge[0] = edge[0] + edge[1];
	row->sum_light[0] = light[0] + light[1];
	row->sum_known[0] = known[0] + known[1];
	for (int x = 1; x < last; x++) {
		row->sum_edge[x] = edge[x - 1] + edge[x] + edge[x + 1];
		row->sum_light[x] = light[x - 1] + light[x] + light[x + 1];
		row->sum_known[x] = known[x - 1] + known[x] + known[x + 1];
	}
	row->sum_edge[last] = edge[last - 1] + edge[last];
	row->sum_light[last] = light[last - 1] + light[last];
	row->sum_known[last] = known[last - 1] + known[last];
	row->summed = true;
}

// The shading terms of a pixel that has none of its own — a PREVIEW, or a sample among previews: the
// mean of those of its eight neighbours that have theirs, read from its window's row sums, so a
// copy shades as the pixels it was copied between (and a column or row of copies carries the relief
// and the edges straight across); where no neighbour has any — inside a coarse block, whose own
// sample has no placed neighbour either — no edge, and the light of a slope that could face any way
// (LOOK_UNKNOWN_LIGHT). Without it, every preview in a moving frame read as flat ground beside its
// relief-lit neighbours: bright lines where a zoom left single columns and rows of copies.
static void borrow_terms(const LookRow* rows, int x, double* edge, double* light) {
	const float known = rows[0].sum_known[x] + rows[1].sum_known[x] + rows[2].sum_known[x];
	if (!(known > 0.0f)) {
		*edge = 1.0;
		*light = LOOK_UNKNOWN_LIGHT;
		return;
	}
	const float share = 1.0f / known;
	*edge = (double)((rows[0].sum_edge[x] + rows[1].sum_edge[x] + rows[2].sum_edge[x]) * share);
	*light = (double)((rows[0].sum_light[x] + rows[1].sum_light[x] + rows[2].sum_light[x]) * share);
}

// Shades row y, whose shading terms and its neighbours' are rows[1] and rows[0], rows[2] (their row
// sums made first, when any of its pixels must borrow); returns how many of its pixels sit at their
// own coordinates.
static int32_t shade_row(const LookJob* job, int y, LookRow* rows) {
	const PixelStore* s = &job->field->front;
	const double* here = rows[1].mu;
	if (rows[1].borrowers > 0) {
		for (int r = 0; r < LOOK_RING; r++) {
			if (!rows[r].summed) sum_row(&rows[r]);
		}
	}
	int32_t ready = 0;
	uint32_t shade_bits = (uint32_t)y * LOOK_R2_STEP_Y;
	for (int x = 0; x < MB_SCREEN_WIDTH; x++, shade_bits += LOOK_R2_STEP_X) {
		const size_t p = (size_t)y * MB_SCREEN_WIDTH + (size_t)x;
		const uint8_t st = s->state[p];
		ready += state_rank(st) >= RANK_REUSED;
		if (state_interior(st) || state_rank(st) == RANK_EMPTY) {
			job->out[p] = PAL_INTERIOR;
			continue;
		}
		const double mu = here[x] >= 0.0 ? here[x] : field_mu(s, p);
		const double hue = hue_of(&job->scale, mu);
		const double hue_threshold = hue_threshold_at(x, y);
		const double shade_threshold = (double)shade_bits * LOOK_R2_SCALE;
		double edge, light;
		if (rows[1].edge[x] >= 0.0) {
			edge = rows[1].edge[x];
			light = rows[1].light[x];
		} else {
			borrow_terms(rows, x, &edge, &light);
		}
		if (job->style == LOOK_SILK) {
			int e = (int)(hue * PAL_RING_SIZE + hue_threshold);
			if (e >= PAL_RING_SIZE) e -= PAL_RING_SIZE;
			job->out[p] = shade_threshold >= edge ? PAL_INK : (uint8_t)(PAL_RING_FIRST + e);
			continue;
		}
		// The level whose lightness (palette_shade_light) is the pixel's, inverted.
		const double lit = edge * light;
		const double dark = lit < 1.0 ? (1.0 - lit) / (1.0 - PAL_SHADE_LOW) : 0.0;
		const double level = dark < 1.0 ? (1.0 - sqrt(dark)) * (PAL_PLANE_SHADES - 1) : 0.0;
		int shade = level > 0.0 ? (int)(level + shade_threshold) : 0;
		if (shade > PAL_PLANE_SHADES - 1) shade = PAL_PLANE_SHADES - 1;
		int h = (int)(hue * PAL_PLANE_HUES + hue_threshold);
		if (h >= PAL_PLANE_HUES) h -= PAL_PLANE_HUES;
		job->out[p] = (uint8_t)(PAL_RING_FIRST + shade * PAL_PLANE_HUES + h);
	}
	return ready;
}

// Rolls a ring of rows down one: the first becomes the last, to be refilled.
static void roll(LookRow* rows) {
	const LookRow spent = rows[0];
	rows[0] = rows[1];
	rows[1] = rows[2];
	rows[2] = spent;
}

// One band of the pass. A row's terms need the mu of the rows either side of it, and a row's shading
// needs the terms of the rows either side of it, so the band measures one row past each of its ends.
static void shade_band(void* arg) {
	const LookBand* band = (const LookBand*)arg;
	const LookJob* job = &s_job;
	const PixelStore* s = &job->field->front;
	const int b = band->first_row / LOOK_BAND_ROWS;
	const double floor = job->scale.edge_floor;
	double* mu[LOOK_RING] = {s_rows[b][0], s_rows[b][1], s_rows[b][2]};  // rows r - 1, r, r + 1
	LookRow rows[LOOK_RING];                                               // rows y - 1, y, y + 1
	for (int r = 0; r < LOOK_RING; r++) {
		rows[r].mu = NULL;
		rows[r].edge = s_edges[b][r];
		rows[r].light = s_lights[b][r];
		rows[r].sum_edge = s_sum_edges[b][r];
		rows[r].sum_light = s_sum_lights[b][r];
		rows[r].sum_known = s_sum_known[b][r];
		rows[r].summed = false;
		rows[r].borrowers = 0;
	}
	const int first = band->first_row;
	fill_row(s, first - 2, mu[0]);
	fill_row(s, first - 1, mu[1]);
	fill_row(s, first, mu[2]);
	rows[2].borrowers = measure_terms(job, mu[0], mu[1], mu[2], first - 1, floor, rows[2].edge, rows[2].light);
	rows[2].summed = false;
	rows[2].mu = mu[1];
	int32_t ready = 0;
	for (int y = first - 1; y < first + LOOK_BAND_ROWS; y++) {
		// Bring the mu ring to rows y, y + 1, y + 2, and measure row y + 1's terms.
		double* spent = mu[0];
		mu[0] = mu[1];
		mu[1] = mu[2];
		mu[2] = spent;
		fill_row(s, y + 2, mu[2]);
		roll(rows);
		rows[2].borrowers = measure_terms(job, mu[0], mu[1], mu[2], y + 1, floor, rows[2].edge, rows[2].light);
		rows[2].summed = false;
		rows[2].mu = mu[1];
		if (y >= first) ready += shade_row(job, y, rows);
	}
	s_band_ready[b] = ready;
}

void look_shade(const Field* f, const LookScale* scale, LookStyle style, uint8_t* out) {
	s_job.field = f;
	s_job.scale = *scale;
	s_job.style = style;
	s_job.out = out;
	measure_grid(f);
	s_job.columns = &s_column_runs;
	s_job.rows = &s_row_runs;
	for (int b = 0; b < LOOK_BANDS; b++) jet_thread_pool_add(shade_band, &s_bands[b]);
	jet_thread_pool_wait();
	int64_t ready = 0;
	for (int b = 0; b < LOOK_BANDS; b++) ready += s_band_ready[b];
	s_ready_share = (double)ready / (double)MB_SCREEN_PIXELS;
}

double look_ready_share(void) {
	return s_ready_share;
}

// --- Measuring the scale -------------------------------------------------------------------------
// The k-th smallest of values[0, count), by Hoare's selection with the middle element as pivot:
// linear on average and, the order of the samples being fixed, the same answer on every machine.
// Reorders the values.
static double select_kth(double* values, int count, int k) {
	int lo = 0, hi = count - 1;
	while (lo < hi) {
		const double pivot = values[lo + (hi - lo) / 2];
		int i = lo, j = hi;
		while (i <= j) {
			while (values[i] < pivot) i++;
			while (values[j] > pivot) j--;
			if (i <= j) {
				const double t = values[i];
				values[i] = values[j];
				values[j] = t;
				i++;
				j--;
			}
		}
		if (k <= j) {
			hi = j;
		} else if (k >= i) {
			lo = i;
		} else {
			break;
		}
	}
	return values[k];
}

bool look_measure(const Field* f, LookScale* scale) {
	const PixelStore* s = &f->front;
	int count = 0;
	double base = 0.0;
	for (int y = LOOK_SAMPLE_STRIDE / 2; y < MB_SCREEN_HEIGHT; y += LOOK_SAMPLE_STRIDE) {
		for (int x = LOOK_SAMPLE_STRIDE / 2; x < MB_SCREEN_WIDTH; x += LOOK_SAMPLE_STRIDE) {
			const size_t p = (size_t)y * MB_SCREEN_WIDTH + (size_t)x;
			if (state_interior(s->state[p]) || state_rank(s->state[p]) == RANK_EMPTY) continue;
			const double mu = field_mu(s, p);
			if (count == 0 || mu < base) base = mu;
			count++;
		}
	}
	if (count < LOOK_MIN_SAMPLES) return false;
	measure_grid(f);
	int n = 0, judged = 0, near = 0;
	for (int y = LOOK_SAMPLE_STRIDE / 2; y < MB_SCREEN_HEIGHT; y += LOOK_SAMPLE_STRIDE) {
		for (int x = LOOK_SAMPLE_STRIDE / 2; x < MB_SCREEN_WIDTH; x += LOOK_SAMPLE_STRIDE) {
			const size_t p = (size_t)y * MB_SCREEN_WIDTH + (size_t)x;
			const uint8_t st = s->state[p];
			if (state_interior(st) || state_rank(st) == RANK_EMPTY) continue;
			const double mu = field_mu(s, p);
			s_samples[n++] = look_log2(mu - base + 1.0);
			double gx, gy;
			if (state_rank(st) >= RANK_REUSED && gradient(f, &s_column_runs, &s_row_runs, x, y, mu, &gx, &gy)) {
				judged++;
				near += edge_light(sqrt(gx * gx + gy * gy), 0.0) < 1.0;
			}
		}
	}
	const double low = select_kth(s_samples, n, (int)(LOOK_LOW_PERCENTILE * (n - 1)));
	const double high = select_kth(s_samples, n, (int)(LOOK_HIGH_PERCENTILE * (n - 1)));
	const double span = high - low > LOOK_MIN_SPAN ? high - low : LOOK_MIN_SPAN;
	const double share = judged > 0 ? (double)near / (double)judged : 0.0;
	const double dense = num_clamp((share - LOOK_DENSE_START) / (LOOK_DENSE_FULL - LOOK_DENSE_START), 0.0, 1.0);
	scale->base = base;
	scale->low = low;
	scale->gain = LOOK_TURNS / span;
	scale->edge_floor = LOOK_FLOOR_CRISP + (LOOK_FLOOR_DENSE - LOOK_FLOOR_CRISP) * dense;
	return true;
}

// --- Each frame ----------------------------------------------------------------------------------
static bool scale_equal(const LookScale* a, const LookScale* b) {
	return a->base == b->base && a->low == b->low && a->gain == b->gain && a->edge_floor == b->edge_floor;
}

// One tick of easing toward the target; arrives exactly once the remaining change is invisible.
static void ease(LookScale* shown, const LookScale* target) {
	LookScale next = *shown;
	next.base += (target->base - shown->base) * LOOK_EASE;
	next.low += (target->low - shown->low) * LOOK_EASE;
	next.gain += (target->gain - shown->gain) * LOOK_EASE;
	next.edge_floor += (target->edge_floor - shown->edge_floor) * LOOK_EASE;
	const double base_gap = target->base - next.base, low_gap = target->low - next.low;
	const double gain_gap = target->gain - next.gain, floor_gap = target->edge_floor - next.edge_floor;
	const bool arrived = (base_gap < 0.0 ? -base_gap : base_gap) * target->gain < LOOK_SETTLED &&
	                     (low_gap < 0.0 ? -low_gap : low_gap) * target->gain < LOOK_SETTLED &&
	                     (gain_gap < 0.0 ? -gain_gap : gain_gap) < LOOK_SETTLED * target->gain &&
	                     (floor_gap < 0.0 ? -floor_gap : floor_gap) < LOOK_SETTLED;
	*shown = arrived ? *target : next;
}

void look_update(const Field* f, bool changed, bool jump) {
	if (s_hold > 0) s_hold--;
	const bool measure = changed || !s_has_target;
	if (measure && (s_hold == 0 || jump || !s_has_target)) {
		LookScale measured;
		if (look_measure(f, &measured)) {
			s_target = measured;
			s_has_target = true;
		}
	}
	if (!s_has_target) {  // an all-interior picture: any scale shades it the same
		s_target.base = s_target.low = 0.0;
		s_target.gain = LOOK_TURNS;
		s_target.edge_floor = LOOK_FLOOR_CRISP;
	}
	if (jump || !s_has_shown) {
		s_shown = s_target;
		s_has_shown = true;
	} else {
		ease(&s_shown, &s_target);
	}
	LookScale applied = s_shown;
	applied.gain *= num_exp2((double)s_density / LOOK_DENSITY_STEPS);
	if (!changed && !s_style_changed && s_has_frame && scale_equal(&applied, &s_shaded)) return;
	look_shade(f, &applied, s_style, s_frame);
	s_shaded = applied;
	s_style_changed = false;
	s_has_frame = true;
}
