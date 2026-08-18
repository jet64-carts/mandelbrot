// carts/mandelbrot/scheduler.h
//
// THE SCHEDULER — deterministic, budgeted, pipelined compute on the thread pool.
//
// Each draw spends QUANTA of work measured in plain lane-iterations (the kernel's own count of what
// it executed, a derivative one at its price), never in time: so which pixels a frame shows depends
// only on the inputs, and two runs of the same scenario draw the same frames on any machine. The
// draw is split around the present:
//
//   scheduler_settle   FENCES the batch the previous draw launched (its results land in the grid
//                      it was planned for), then brings the field to the new spec:
//                        a declared jump, or a view reuse cannot serve   RESTART: every tile's
//                            coarsest level at once, then finer levels, in the restart quantum
//                            at rest and the motion quantum on the move;
//                        a moving view                                  GATHER the old picture
//                            through the column and row maps (reuse.h), then compute the new
//                            strips in the motion quantum;
//                        only the iteration limit                       LIMIT: decide what can be
//                            decided without iterating, reuse the rest;
//                        a view come to rest                            SNAP the grid to the ideal
//                            and refine in place, certifying whole cells (certify.h).
//                      Every reconcile runs on the pool, one job per tile.
//   (present)          the caller shades the field into the look's frame and presents it (look.h).
//   scheduler_launch   plans the next batch within the quantum, starts it, and returns: the
//                      workers compute while the host presents and runs the next _update. It
//                      launches every draw, moving or not — a batch that lands after the view has
//                      moved is carried by the next gather like any other pixel.
//
// Planning is breadth-first — every tile's level L before any tile's level L + 1 — dealing
// sample rows round-robin, nearest the centre first, priced by the pixels each step still needs
// times what that tile's samples have cost, and submitted longest first (LPT). A tile certifies a
// level's cells only once the tiles its cells reach into — right, below, below-right — finished
// that level in an earlier batch, so every corner it reads is final and the converged picture does
// not depend on how the work was cut into batches.
// On the web player there are no threads: a batch runs inside the next fence, on the main thread,
// so its quanta are one core's share.
#ifndef MANDELBROT_SCHEDULER_H
#define MANDELBROT_SCHEDULER_H

#include "deep.h"
#include "field.h"
#include "reuse.h"

// The techniques Phase 2 adds, each switchable so the harness can measure what it is worth.
typedef struct {
	bool reuse;                 // carry pixels across view changes (otherwise every change restarts)
	bool pipeline;              // launch batches while the view moves (otherwise only at rest)
	bool exterior_certificate;  // interpolate certified exterior cells
	bool interior_certificate;  // fill certified interior cells
	bool longest_first;         // submit a batch's units costliest first (otherwise in plan order)
	DeepOptions deep;           // the deep view's techniques (deep.h)
} SchedulerOptions;

/** Everything on: the cart's own configuration. */
SchedulerOptions scheduler_default_options(void);

/** Builds the schedule and starts the pool. Call once from _init, after field_init. */
bool scheduler_init(Field* field, const SchedulerOptions* options);

/**
 * The fence and the reconcile. Call first in _draw, with the spec the frame should show and
 * whether the camera is at rest (a view at rest is snapped and refined).
 */
void scheduler_settle(const FieldSpec* spec, bool at_rest);

/** Plans and starts the next batch. Call last in _draw, after the field has been presented. */
void scheduler_launch(bool at_rest);

/** Waits for any batch in flight and stops the pool. Call from _end. */
void scheduler_end(void);

/**
 * The share of the picture's pixels that are exact or certified at their own grid coordinates, in
 * [0, 1]: what is left is what the schedule will still compute. A moving view keeps its carried
 * pixels, so its share stays high; a restart starts it from nothing.
 */
double scheduler_progress(void);

/** True once every tile's schedule has run to its end. */
bool scheduler_converged(void);

// What the last settle and launch did, for the harness's instruments and the HUD.
typedef struct {
	int reconcile;          // the last reconcile's kind (ReconcileKind), or -1 if the spec held
	double matched_share;   // the smaller share of columns and rows the last gather carried
	int64_t settle_work;    // work the last settle computed synchronously, in plain lane-iterations
	int64_t batch_work;     // work the last fenced batch computed, likewise
	int64_t certified;      // pixels filled by certificates since the last reconcile
	int64_t samples;        // samples iterated since the last reconcile
	const AxisMap* columns; // the last gather's maps
	const AxisMap* rows;
	bool kept;              // the last restart kept a preview of the old picture (at depth)
} SchedulerStats;

const SchedulerStats* scheduler_stats(void);

/** Bytes of WRAM scheduler_init reserves (the certificate corners). */
uint32_t scheduler_bytes(void);

#endif // MANDELBROT_SCHEDULER_H
