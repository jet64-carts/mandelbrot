// carts/mandelbrot/field.h
//
// THE FIELD — what the compute knows about every pixel, how good that knowledge is, and the grid
// it was learnt on.
//
// Per pixel the field holds the smooth iteration count as an integer and a fraction,
// mu = n - phi (n the escape iteration, phi in [0, 1]): n is what an iteration-limit change needs
// to decide a pixel without iterating it (design.md §6.7), and phi as a float resolves mu to 6e-8
// of an iteration at any depth — as fine as a double mu would past 2^29, and far finer than
// neighbouring pixels' mu ever differ. Beside them is a state byte — the pixel's RANK and FATE —
// and a preview's DISTANCE from the coordinate its value belongs to. The field holds no colour: the
// look shades the picture from it (look.h).
//
// THE GRID is a tensor grid: pixel (i, j) was computed at (column_re[i], row_im[j]), the ACTUAL
// coordinates. When the view moves, columns and rows are matched to the new view's IDEAL
// coordinates and carried, keeping their actual coordinates (reuse.c). When it comes to rest, a
// grid that is one lattice off its ideal has the view moved onto it (camera_rest_on_grid); any
// other is SNAPPED to the ideal, and the pixels of its moved columns and rows recomputed.
//
// Ownership (design.md §6.2): the grids, the spec and the choice of front buffer belong to the
// main thread and change only while no batch is in flight. A worker writes only the pixels of the
// tile its unit owns — in the back buffer during a gather, in the front buffer otherwise — and
// reads the front buffer and the grids, which nothing writes meanwhile.
#ifndef MANDELBROT_FIELD_H
#define MANDELBROT_FIELD_H

#include "config.h"
#include "view.h"

// A pixel's RANK: how good its value is, best last. MONOTONE QUALITY: within one grid, a pixel's
// value is only ever replaced by one of equal or better rank (design.md §6.6).
typedef enum {
	RANK_EMPTY = 0,      // nothing is known
	RANK_PREVIEW = 1,    // a coarse sample's block fill, or a neighbour's copy
	RANK_REUSED = 2,     // from an earlier spec: within half a pixel of its coordinate (a pan's
	                     // reach; a zoom's is a pixel, reuse.h), or under a lower iteration limit
	RANK_CERTIFIED = 3,  // interpolated or filled at its coordinate under a Koebe certificate
	RANK_EXACT = 4,      // iterated at its coordinate
} Rank;

// How a pixel's orbit ended.
typedef enum {
	FATE_ESCAPED = 0,
	FATE_SHAPE = 1,      // interior: the main cardioid or the period-2 bulb
	FATE_CYCLE = 2,      // interior: the orbit came back exactly to an earlier state
	FATE_LIMIT = 3,      // interior: still bounded at the iteration limit — the one interior fate
	                     // a higher limit could overturn
	FATE_CERTIFIED = 4,  // interior: filled under an interior certificate
} Fate;

// The state byte: the rank, the fate, and the OFFSET bit — set when a snap moved the pixel's
// column or row, so its value belongs to a coordinate up to half a pixel from the grid's. Carried
// once more by a moving view, such a value could be a whole pixel out, so the carry demotes it.
#define STATE_RANK_MASK 0x07u
#define STATE_FATE_SHIFT 3
#define STATE_FATE_MASK 0x38u
#define STATE_OFFSET 0x40u

static inline uint8_t state_make(Rank rank, Fate fate) {
	return (uint8_t)((unsigned)rank | ((unsigned)fate << STATE_FATE_SHIFT));
}
static inline Rank state_rank(uint8_t state) {
	return (Rank)(state & STATE_RANK_MASK);
}
static inline Fate state_fate(uint8_t state) {
	return (Fate)((state & STATE_FATE_MASK) >> STATE_FATE_SHIFT);
}
static inline bool state_interior(uint8_t state) {
	return state_fate(state) != FATE_ESCAPED;
}
static inline uint8_t state_with_rank(uint8_t state, Rank rank) {
	return (uint8_t)((state & ~STATE_RANK_MASK) | (unsigned)rank);
}

// A pixel needs computing when its value is not already exact or certified at its coordinate.
static inline bool state_needs_work(uint8_t state) {
	return state_rank(state) < RANK_CERTIFIED;
}

// A PREVIEW's DISTANCE: how far the coordinate its value belongs to may lie from the pixel's own,
// counted in 1/DISTANCE_PIXEL of a pixel of the store's grid. Previews are ranked among themselves
// by it — within PREVIEW, a value is only ever replaced by a nearer one — so a coarse sample a few
// pixels away replaces a copy carried in from across the screen, and never a neighbour's copy a
// pixel away (design.md §6.6). A value at its own coordinate (reused or better) has distance 0;
// DISTANCE_FAR means nothing near is known: an empty pixel, a picture held after a jump, and any
// copy farther than the type counts.
#define DISTANCE_PIXEL 256
#define DISTANCE_FAR UINT16_MAX

/** `pixels` as a distance, to the nearest unit: DISTANCE_FAR beyond the type's range (or for a NaN). */
static inline uint16_t distance_of_pixels(double pixels) {
	const double units = pixels * DISTANCE_PIXEL + 0.5;
	return units < (double)DISTANCE_FAR ? (uint16_t)units : DISTANCE_FAR;
}

/**
 * A distance measured on one grid, seen on the next: `scale` times as many pixels (the old spacing
 * over the new), and farther by the way from the pixel to the one whose value it shows, `gap_x`
 * across and `gap_y` down. Far stays far.
 */
static inline uint16_t distance_carried(uint16_t distance, double scale, uint16_t gap_x, uint16_t gap_y) {
	if (distance == DISTANCE_FAR || gap_x == DISTANCE_FAR || gap_y == DISTANCE_FAR) return DISTANCE_FAR;
	const double across = (double)gap_x, down = (double)gap_y;
	// A gap along one axis — every carried pixel's, and most copies' — needs no square root.
	const double gap = gap_y == 0 ? across : gap_x == 0 ? down : sqrt(across * across + down * down);
	const double units = (double)distance * scale + gap + 0.5;
	return units < (double)DISTANCE_FAR ? (uint16_t)units : DISTANCE_FAR;
}

// Everything that decides what the field's pixels should hold. Two equal specs describe the same
// picture exactly. The epoch counts declared jumps (a reset): a new epoch always restarts, even
// where the old picture could have been reused.
typedef struct {
	View view;
	int64_t max_iterations;
	uint32_t epoch;
	bool perturbed;  // deeper than the crossover: computed by perturbation (view_perturbed)
} FieldSpec;

// One of the two pixel buffers.
typedef struct {
	int32_t* n;          // the escape iteration; the limit for an interior pixel
	float* phi;          // mu = n - phi, phi in [0, 1]; 0 for an interior pixel
	uint8_t* state;      // Rank | Fate << STATE_FATE_SHIFT | STATE_OFFSET
	uint16_t* distance;  // a preview's distance (above); 0 for a value at its own coordinate
} PixelStore;

// One axis of the grid: a coordinate per column (or row), the ideal coordinate the spec asks
// for, and whether the two are the same (to within field_ideal_tolerance); and the LATTICE a pan
// extends — a pan computes each new column at lattice + n x lattice_spacing, the lattice point
// nearest its ideal, so the axis stays the one lattice it was (reuse.h).
typedef struct {
	double* actual;
	double* ideal;
	uint8_t* at_ideal;       // actual == ideal
	double lattice;          // a coordinate on the axis's lattice
	double lattice_spacing;  // the pixel spacing it was set at: a pan keeps it, a zoom sets another
} GridAxis;

typedef struct {
	PixelStore front;  // the picture
	PixelStore back;   // where a reconcile pass writes the next picture; swapped in after it

	GridAxis columns;  // MB_SCREEN_WIDTH; the real axis
	GridAxis rows;     // MB_SCREEN_HEIGHT; the imaginary axis, decreasing down the screen
	GridAxis next_columns, next_rows;  // a reconcile's new grid, swapped in with the pixels

	FieldSpec spec;  // what the grid and the pixels describe
	bool has_spec;   // false until the first spec is set

	// The point the grid's coordinates are offsets from. Direct iteration keeps it at zero, so a
	// coordinate is c itself; at depth it sits near the view (a double cannot hold c there, but it
	// holds c's offset from a point nearby), and moves with the view, rarely and exactly
	// (field_move_origin).
	DeepPoint origin;
} Field;

/** Bytes of WRAM field_init reserves. */
uint32_t field_bytes(void);

/**
 * Reserves the field's buffers from WRAM. They are not cleared: the first spec's restart writes
 * every pixel before anything reads one. False if WRAM ran out.
 */
bool field_init(Field* f);

/** The smooth iteration count of pixel p of a store. */
static inline double field_mu(const PixelStore* s, size_t p) {
	return (double)s->n[p] - (double)s->phi[p];
}

/** The rank a pixel shows: an exact or certified value off the ideal grid is only reused. */
static inline Rank field_effective_rank(const Field* f, int x, int y) {
	const Rank r = state_rank(f->front.state[(size_t)y * MB_SCREEN_WIDTH + (size_t)x]);
	const bool ideal = f->columns.at_ideal[x] && f->rows.at_ideal[y];
	return (r >= RANK_CERTIFIED && !ideal) ? RANK_REUSED : r;
}

/** Computes the ideal grid of `spec` into the next grid's ideal arrays and records the spec. */
void field_plan_spec(Field* f, const FieldSpec* spec);

/** The ideal coordinates of `spec`, as offsets from `origin`, into `columns` and `rows`. */
void field_ideal_grid(const FieldSpec* spec, const DeepPoint* origin, double* columns, double* rows);

/**
 * Moves the origin by (d_re, d_im), exactly, and re-expresses the grid's coordinates by the same
 * offset. A move of many pixels subtracts exactly (Sterbenz: every coordinate lies within a factor
 * of two of the offset); a move to an unrelated point rounds each coordinate to the nearest double,
 * and the caller treats the pixels as values from another grid.
 */
void field_move_origin(Field* f, double d_re, double d_im);

/** Makes `origin` the origin, re-expressing the grid's coordinates (rounding, as above). */
void field_set_origin(Field* f, const DeepPoint* origin);

/** Swaps in the next grid (and, when `pixels`, the back pixel store) after a reconcile pass. */
void field_swap(Field* f, bool pixels);

/** True when the field already describes `spec`. */
bool field_matches(const Field* f, const FieldSpec* spec);

/** True when every column and row sits at its ideal coordinate. */
bool field_grid_is_ideal(const Field* f);

/**
 * How far a grid coordinate may lie from its ideal and still count as ideal, in the plane: a
 * millionth of a pixel, or sixteen units in the last place of the centre's own magnitude, whichever
 * is the larger. Re-expressing a coordinate about a new centre rounds it by an ulp or two — which at
 * the depth limit of direct iteration (view.h) is a few thousandths of a pixel — and re-iterating a
 * column for that would change nothing visible.
 */
double field_ideal_tolerance(const Field* f);

/**
 * When the grid is ONE LATTICE off its ideal — every column the same distance from its ideal
 * coordinate, and every row, to within the tolerance, as a pan leaves it — writes that distance for
 * each axis, in the plane, and returns true. False when the grid is ideal already, or is not one
 * lattice (a zoom's grid is not).
 */
bool field_lattice_offset(const Field* f, double* d_re, double* d_im);

#endif // MANDELBROT_FIELD_H
