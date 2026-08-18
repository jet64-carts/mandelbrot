// carts/mandelbrot/deep.h
//
// THE DEEP VIEW — which reference orbit serves the picture, and which is being built to serve the
// next one.
//
// Two reference slots (reference.h): the CURRENT one the batches compute against, and the NEXT one,
// built a slice per batch as one more job on the pool — never a long call, so no tick ever waits
// for a reference (on the web, a slice runs at a fence on the main thread, a quarter the size).
//
// THE POLICY. The current reference SERVES a view while
//   * its precision reaches the view's pixel spacing (reference.h);
//   * it lies within DEEP_NEIGHBOURHOOD_PIXELS of the view's centre, and its BLA table was built for
//     the view's largest |delta_c|;
//   * its orbit covers the iteration limit (a nucleus's periodic orbit covers any limit);
// and the next one is started AHEAD of need, while the current one still serves: when the view
// comes within DEEP_LOOKAHEAD_OCTAVES of its precision, strays half its neighbourhood, or has
// zoomed DEEP_REREFERENCE_OCTAVES since it was chosen. So a continuous dive never waits. A reference
// is built for DEEP_SPAN_OCTAVES beyond the view it was started at.
//
// CHOOSING WELL. The next reference is the NUCLEUS of the minibrot that dominates the view where
// one can be found, and the view's centre otherwise:
//   1. the view's centre is iterated by perturbation against the current reference, recording the
//      iterations where |z_n| reaches a new low — the atom domains it lies in, whose periods name
//      the minibrots near it;
//   2. for the deepest of those periods first, Newton's method on z_p(c) = 0 runs in perturbation
//      form against the current reference, in doubles (reference_find_nucleus);
//   3. the nucleus is built in the bignum and refined there by the same Newton step (reference.c),
//      so its orbit is periodic to the reference's precision.
// With no current reference (the first dive past the crossover, or a jump) the view's centre is
// built first, and the nucleus found from it next. Changing the reference never restarts the
// picture: pixels already computed stay, and only batches planned afterwards use the new one.
#ifndef MANDELBROT_DEEP_H
#define MANDELBROT_DEEP_H

#include "field.h"
#include "kernel.h"
#include "reference.h"

// The reference must lie within this many pixels of the view's centre: delta_c's rounding is then
// 2^-33 of a pixel, and the BLA table's radii still leave the view room.
#define DEEP_NEIGHBOURHOOD_PIXELS 4096.0
// Start the next reference this many octaves before the current one's precision runs out ...
#define DEEP_LOOKAHEAD_OCTAVES 12.0
// ... and after this many octaves of zoom in any case: Newton against the current reference finds
// a nucleus to a double's precision of its offset, so the two may not drift too far apart.
#define DEEP_REREFERENCE_OCTAVES 16.0
// A reference is built to serve this many octaves deeper than the view it was started at.
#define DEEP_SPAN_OCTAVES 24.0
// Its BLA table is valid for delta_c up to this many times the view's own largest, so the view may
// zoom out a little before the table must be rebuilt.
#define DEEP_TABLE_REACH 8.0
// A pipeline starts this many octaves above the crossover while the view dives toward it.
#define DEEP_CROSSOVER_LOOKAHEAD_OCTAVES 3.0

// The techniques the policy may use; the cart uses them all, and the harness's ablations switch
// each off to measure what it buys.
typedef struct {
	bool bla;        // build every reference's BLA table (otherwise plain perturbation)
	bool nucleus;    // search for a minibrot's nucleus (otherwise every reference is a centre)
	bool lookahead;  // start the next reference ahead of need (otherwise only once the current
	                 // one stops serving, or the view has crossed over)
} DeepOptions;

/** Bytes of WRAM deep_init reserves (two references). */
uint32_t deep_bytes(void);

/** Reserves the two reference slots. Call once from _init. */
bool deep_init(const DeepOptions* options);

/**
 * After the settle's reconcile, with no batch in flight: promotes a finished reference, and starts
 * building the next when the field's view needs one or soon will.
 */
void deep_settle(const Field* f);

/**
 * Whether the field's picture can be computed now: always in direct mode; at depth, when a
 * reference serves the view. If so, writes the reference and the offset of the field's origin
 * from the reference's point (delta_c = that offset + a pixel's grid offset).
 */
bool deep_ready(const Field* f, KernelReference* ref, double* origin_re, double* origin_im);

/**
 * Adds a slice of the pipeline's work to the batch just submitted — or, `alone`, with no batch
 * beside it, DEEP_ALONE_SLICES of them (deep.c). False if there is nothing to do.
 */
bool deep_launch(bool alone);

// What the deep view is doing, for the harness and the HUD.
typedef struct {
	bool serving;         // a reference serves the view
	int32_t period;       // the current reference's period, 0 for a centre reference
	int32_t length;       // its orbit's length
	int frac;             // its precision, in limbs
	bool building;        // a reference is being built
	int references_built; // since boot
	int nuclei_found;
} DeepStats;

const DeepStats* deep_stats(void);

#endif // MANDELBROT_DEEP_H
