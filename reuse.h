// carts/mandelbrot/reuse.h
//
// PIXEL REUSE — carrying a picture from one view to the next, the way XaoS does.
//
// The field is a tensor grid, so a view change is two one-dimensional problems. Each new IDEAL
// column is MATCHED to an old column whose actual coordinate lies near it; a matched column keeps
// its old actual coordinate (its pixels are exact there), and an unmatched one is MISSING, to be
// computed. Rows alike. How near depends on the move:
//
//   PAN   (the spacing unchanged) within half a new pixel. The windows of half a pixel either side
//         of consecutive ideal coordinates are disjoint and in order, so choosing the nearest old
//         column inside each window is optimal — it matches every window that can be matched,
//         each with the least error — and needs one sweep of two pointers.
//   ZOOM  within ZOOM_REACH_PIXELS, a whole pixel, by an optimal matching (the most columns
//         carried, then the least displacement). A zoom scales every column about its anchor, so
//         columns computed on different ticks drift through each other's phase: under a half-pixel
//         reach two of them keep landing in one window, one is dropped, and the column it leaves
//         missing is computed again — 10% of each axis every tick at one octave a second (measured
//         on the seahorse, design.md §6.7), where the zoom itself uncovers 0.6%. A pixel's reach
//         lets the drifting columns pass each other: 1.2% (simulated, and 1.0% measured), at a mean
//         displacement of 0.31 pixels against 0.25. At rest the snap puts every column on its ideal.
//
// Either way the result is monotone and one-to-one.
//
// A RECONCILE pass then rebuilds the field for the new spec, one job per tile on the pool:
//   RESTART  every pixel EMPTY (a declared jump, or a change reuse cannot serve) — or, at depth,
//            where nothing can be computed until a reference is ready, every pixel a PREVIEW of the
//            old picture: WARPED onto the new view (each column showing the old column nearest its
//            coordinate) when the view moved too fast to carry, HELD in place after a jump;
//   GATHER   every pixel carried through the column and row maps; a missing pixel previews
//            whichever of its nearest carried neighbours — before it and after it along each axis —
//            holds the value nearest it (field.h: a preview's distance), so a band of missing rows
//            is shown from both its sides rather than copied on and on from one; an OFFSET pixel is
//            demoted to PREVIEW (§ field.h); across
//            the crossover (direct iteration to perturbation, or back) every exact or certified
//            value is DEMOTED to REUSED: it is a value from the other arithmetic, recomputed at rest;
//   LIMIT    the iteration limit changed and nothing else: escaped pixels at or beyond a lower
//            limit become interior without iterating; interior-by-limit pixels under a higher one
//            are kept but REUSED, to be iterated again (a gather applies the same rules);
//   SNAP     the view is at rest: every column and row moves to its ideal coordinate, and each
//            exact or certified pixel in a moved one becomes REUSED | OFFSET — its value shown
//            until it is recomputed, never replaced by anything worse.
// Each job also counts, step by step of its tile's schedule, the pixels that still need work, which
// is what the planner prices the next batches with. Nothing here shades: the look reads the field
// (look.h).
#ifndef MANDELBROT_REUSE_H
#define MANDELBROT_REUSE_H

#include "field.h"
#include "tiles.h"
#include "work.h"

typedef enum {
	RECONCILE_RESTART,
	RECONCILE_GATHER,
	RECONCILE_LIMIT,
	RECONCILE_SNAP,
} ReconcileKind;

// A zoom's reach, in new pixels, and the most old columns one ideal chooses among (a window of two
// pixels holds two or three).
#define ZOOM_REACH_PIXELS 1.0
#define ZOOM_CHOICES 8
#define ZOOM_STATES (ZOOM_CHOICES + 1)

// The two sides a missing column may take its preview from: its nearest carried neighbour before
// it, and after it.
enum { REUSE_BEFORE, REUSE_AFTER, REUSE_SIDES };

// One axis's map from the old grid to the new.
typedef struct {
	int16_t source[MB_SCREEN_WIDTH];               // the old column (row) a new one carries, or
	                                               // REUSE_MISSING
	int16_t fill[REUSE_SIDES][MB_SCREEN_WIDTH];    // the old columns whose pixels it may show: its
	                                               // source on both sides, or when missing its
	                                               // nearest carried neighbour's on each (REUSE_MISSING
	                                               // where that side has none)
	uint16_t gap[REUSE_SIDES][MB_SCREEN_WIDTH];    // how far each one's coordinate lies from the new
	                                               // one's, as a preview's distance (field.h): 0
	                                               // where it carries its own
	int matched;
} AxisMap;
#define REUSE_MISSING (-1)

// What every reconcile job of one pass shares, frozen for the pass.
typedef struct {
	Field* field;
	ReconcileKind kind;
	int64_t old_limit, new_limit;
	const AxisMap* columns;  // GATHER, and a RESTART that keeps previews
	const AxisMap* rows;
	bool demote;             // GATHER across the crossover: carried exact and certified values are REUSED
	bool keep;               // RESTART at depth: every pixel a PREVIEW through the maps, not EMPTY
	double distance_scale;   // a preview's distance on the old grid times this is its distance on the
	                         // new one: the old pixel spacing over the new
} ReconcileContext;

typedef struct {
	const ReconcileContext* context;
	int tile;
	int32_t need[TILE_STEPS];  // written by the job: pixels needing work, per schedule step
	int32_t coarse[WORK_COARSE_LEVELS];  // and its pixels by coarse level (work.h)
} ReconcileUnit;

/**
 * Matches the current grid to the ideal grid the field has planned (field_plan_spec) and writes
 * the new actual grid into the field's next grid. Returns the share of columns and of rows
 * matched, the smaller of the two, in [0, 1]. Main thread, no batch in flight.
 */
double reuse_match(Field* f, AxisMap* columns, AxisMap* rows);

/** Sets the next grid to the planned ideal grid (a restart). */
void reuse_ideal_grid(Field* f);

/**
 * Maps for a restart that HOLDS the old picture in place: every pixel previews itself — a picture of
 * somewhere else, so every preview is FAR.
 */
void reuse_hold_maps(AxisMap* columns, AxisMap* rows);

/**
 * Maps for a restart that WARPS the old picture onto the planned grid: every new column previews
 * one of the two old columns whose coordinates lie either side of its ideal one (rows alike),
 * however far, whichever holds the nearer value, at the distance between them.
 */
void reuse_warp_maps(const Field* f, AxisMap* columns, AxisMap* rows);

/** Moves every column and row that is off its ideal coordinate onto it (a snap). */
void reuse_snap_grid(Field* f);

/** Keeps the current grid as the next (a limit change). */
void reuse_keep_grid(Field* f);

/** Runs one ReconcileUnit (passed as the job's argument). The pool's job function. */
void reconcile_run(void* unit);

#endif // MANDELBROT_REUSE_H
