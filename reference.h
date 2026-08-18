// carts/mandelbrot/reference.h
//
// THE REFERENCE ORBIT — one point iterated in the bignum, so that every pixel at depth can follow
// it in doubles (kernel.h's perturbation variants).
//
// A reference is built in DETERMINISTIC SLICES: each slice is a bounded amount of bignum work
// (REFERENCE_SLICE_PRODUCTS limb products, an iteration's fixed work counted as
// REFERENCE_ITERATION_OVERHEAD more, so a slice costs about the same at any precision), never one
// long call, so it can run as one job of a batch on the pool — or, on the web, at a fence on the
// main thread — without ever holding a tick. The orbit is stored as doubles, Z_0 = 0 ..
// Z_length, then its BLA table is built from it, level by level, in slices too.
//
// Two kinds:
//   * a CENTRE reference iterates the view's centre until it escapes or reaches the iteration
//     limit; on the way it keeps the RECORD LOWS of |Z_n| — the atom domains the centre lies in —
//     whose periods name the minibrots near it;
//   * a NUCLEUS reference iterates the nucleus of such a minibrot for exactly one period p: its
//     orbit is periodic (Z_p = 0 but for rounding), so a lane rebases at m = p and follows the same
//     short orbit again, the BLA table covers every period with a handful of blocks, and interior
//     pixels of that minibrot cost a few steps a period instead of p.
// The nucleus is found without bignum division (reference_find_nucleus): Newton's method on
// z_p(c) = 0 in PERTURBATION form against a centre reference, in doubles, lands within double
// precision of the offset; the nucleus is then centre + offset, exactly, in the bignum.
//
// PRECISION. A reference serves pixel spacings down to 2^-(bits - REFERENCE_GUARD_BITS), where
// bits is its fraction precision (design.md §9.1 has the guard's measurement).
#ifndef MANDELBROT_REFERENCE_H
#define MANDELBROT_REFERENCE_H

#include "config.h"
#include "kernel.h"
#include "point.h"

// The longest orbit a reference stores: 2^18 iterations. Two references live at once (the one in
// use and the one being built ahead of it), each with its orbit (16 bytes an iteration) and its
// BLA table (at most one 40-byte entry an iteration), 29.4 MB of WRAM between them.
#define REFERENCE_CAPACITY (1 << 18)
_Static_assert(REFERENCE_CAPACITY >= MB_ITERATIONS_MAX, "a centre reference's orbit must reach any iteration limit");
#define REFERENCE_BYTES (2u * (REFERENCE_CAPACITY + 1u) * (unsigned)sizeof(double) + \
                         REFERENCE_CAPACITY * (unsigned)sizeof(BlaEntry))
// The budget: the two references take less than half of WRAM, and the field (19.4 MB) and the
// certificate corners the rest (design.md §2.3; the soak asserts the whole).
_Static_assert(2u * REFERENCE_BYTES <= JET_WRAM_CAPACITY_BYTES / 2u, "the references must fit half of WRAM");

// Bits of precision beyond the pixel spacing: the reference's own rounding, amplified along its
// orbit, must stay far below a pixel (design.md §9.1).
#ifndef REFERENCE_GUARD_BITS
#define REFERENCE_GUARD_BITS 64
#endif

// The depth limit is where the bignum runs out (config.h, bignum.h): a reference serving it, guard
// and all, and the view's centre a pan's double below its spacing, both fit.
#define REFERENCE_PAN_BITS 62
_Static_assert(MB_DEPTH_LIMIT_BITS + REFERENCE_GUARD_BITS <= BIG_FRAC_LIMBS * BIG_LIMB_BITS,
               "a reference at the depth limit must fit the bignum");
_Static_assert(MB_DEPTH_LIMIT_BITS + REFERENCE_PAN_BITS <= BIG_FRAC_LIMBS * BIG_LIMB_BITS,
               "the centre at the depth limit must fit the bignum");

// Limb products per slice of the orbit and BLA entries per slice of the table: about a millisecond
// of one core each (an iteration costs 0.10 us plus 0.66 ns a limb product, a table entry about
// 30 ns: tests/mandelbrot/deep_test.c --bench). The web's slices, run on the main thread, a quarter
// of those.
#define REFERENCE_SLICE_PRODUCTS 1500000
#define REFERENCE_SLICE_PRODUCTS_WEB 375000
#define REFERENCE_SLICE_ENTRIES 32768
#define REFERENCE_SLICE_ENTRIES_WEB 8192

// An iteration's fixed work beyond its limb products — the additions, the two conversions to
// double, the escape test and the records — in limb products: 0.10 us at 0.66 ns a product.
#define REFERENCE_ITERATION_OVERHEAD 150

typedef enum {
	REFERENCE_EMPTY,  // nothing is being built
	REFERENCE_ORBIT,  // iterating
	REFERENCE_TABLE,  // building the BLA table
	REFERENCE_READY,  // usable
} ReferenceStage;

// The atom domains the centre's orbit passed through, deepest last.
#define REFERENCE_RECORDS 16

// A nucleus reference is refined in the bignum after its orbit: Newton's step from its own Z_p and
// dZ_p/dc (the latter in doubles), until the step is below REFERENCE_REFINE_BITS of the spacing it
// serves, at most this many times.
#define REFERENCE_REFINEMENTS 4
#define REFERENCE_REFINE_BITS 16

typedef struct {
	// What is being built.
	DeepPoint c;          // the reference point, at the reference's precision
	int frac;             // fraction limbs of precision
	int32_t target;       // iterations to compute: the limit, or the period of a nucleus
	int32_t period;       // a nucleus reference's period, or 0
	double max_delta_c;   // the table is valid for |delta_c| up to this
	double epsilon;       // the table's validity tolerance (kernel.h's KernelReference.epsilon)

	// Progress.
	ReferenceStage stage;
	Big zr, zi;           // Z_computed, in the bignum
	int32_t computed;     // Z_0 .. Z_computed are stored
	bool escaped;         // Z_computed lies beyond the bailout
	int table_level;      // the BLA level being built (1 ..), and the next entry of it
	int32_t table_entry;

	// What the centre's orbit found: the iterations of its record lows of |Z_n| (n >= 1).
	int32_t record_at[REFERENCE_RECORDS];
	int record_count;
	double record_sq;     // the lowest |Z_n|^2 so far

	// A nucleus's refinement: dZ_n/dc along the orbit, and the passes made.
	double dzr, dzi;
	int refinements;

	// Storage, reserved once (reference_reserve).
	double* re;
	double* im;
	BlaEntry* table;
	BlaEntry* level_start[KERNEL_BLA_LEVELS + 1];
	int32_t level_count[KERNEL_BLA_LEVELS + 1];
	int levels;
} Reference;

/** Bytes of WRAM one reference reserves. */
uint32_t reference_bytes(void);

/** Reserves one reference's storage from WRAM; false if it ran out. */
bool reference_reserve(Reference* r);

/**
 * Starts building a reference at `c` with `frac` fraction limbs: `target` iterations (at most
 * REFERENCE_CAPACITY), `period` > 0 for a nucleus (then target = period), a BLA table valid for
 * |delta_c| <= max_delta_c at tolerance `epsilon` (0 for no table).
 */
void reference_begin(Reference* r, const DeepPoint* c, int frac, int32_t target, int32_t period, double max_delta_c,
                     double epsilon);

/**
 * One slice: iterates or builds the table for at most `products` limb products or `entries` table
 * entries, deterministically. Returns the stage reached. deep.c's job runs it, a slice a batch.
 */
ReferenceStage reference_advance(Reference* r, int64_t products, int32_t entries);

/** The kernel's view of a READY reference. */
KernelReference reference_kernel(const Reference* r);

/** The pixel spacing's log2 below which the reference's precision no longer suffices. */
int reference_deepest_log2(const Reference* r);

/** The fraction limbs a reference needs to serve pixel spacings down to 2^deepest_log2. */
int reference_limbs_for(double deepest_log2);

// --- Searching against a READY reference, in resumable slices ------------------------------------------
typedef enum {
	SEARCH_RUNNING,
	SEARCH_FOUND,
	SEARCH_FAILED,
} SearchStatus;

// The record lows of |z_n| along the orbit of reference + (dc_re, dc_im), by perturbation.
typedef struct {
	double dc_re, dc_im;
	double x, y;          // delta
	int32_t m;            // the reference index
	int32_t n;            // iterations taken
	int32_t limit;
	double record_sq;
	int32_t record_at[REFERENCE_RECORDS];
	int record_count;
} PeriodSearch;

/** Starts a period search for the point at (dc_re, dc_im) from the reference, to `limit` iterations. */
void reference_period_begin(PeriodSearch* s, double dc_re, double dc_im, int32_t limit);

/** One slice of at most `steps` perturbation steps; FOUND when the orbit escaped or hit the limit. */
SearchStatus reference_period_advance(PeriodSearch* s, const Reference* r, int64_t steps);

// Newton's method for the nucleus of one period, in perturbation form and doubles.
typedef struct {
	int32_t period;
	double start_re, start_im;  // where it began, as an offset from the reference
	double cr, ci;              // the estimate
	double tolerance, reach;
	int step;                   // Newton steps taken
	// the orbit under way at the estimate
	double x, y, dx, dy, zr, zi;
	int32_t m, n;
} NucleusSearch;

/**
 * Starts Newton's method for the nucleus of period `period` near reference + start; it succeeds
 * when a step is below `tolerance`, and fails beyond `reach` from the start, on escape, on a
 * vanishing derivative, or after a bounded number of steps.
 */
void reference_nucleus_begin(NucleusSearch* s, int32_t period, double start_re, double start_im, double tolerance,
                             double reach);

/** One slice of at most `steps` perturbation steps. FOUND leaves the nucleus's offset in (cr, ci). */
SearchStatus reference_nucleus_advance(NucleusSearch* s, const Reference* r, int64_t steps);

/**
 * The whole search at once, for the harness: the nucleus of period `period` near r's c + start.
 * Writes its offset from r's c and returns true on success.
 */
bool reference_find_nucleus(const Reference* r, int32_t period, double start_re, double start_im, double tolerance,
                            double reach, double* nucleus_re, double* nucleus_im);

#endif // MANDELBROT_REFERENCE_H
