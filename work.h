// carts/mandelbrot/work.h
//
// A UNIT OF WORK — a run of one tile's schedule steps: its samples iterated, their values stored
// and block-filled into the field, and the tile's cells certified at each level they complete. A
// unit shades nothing: the look turns the field into the picture (look.h).
//
// It works IN PLACE and MONOTONELY: a sample whose pixel is already exact or certified is skipped;
// an iterated sample becomes EXACT; a block fill writes PREVIEW only over a pixel that is EMPTY, or a
// PREVIEW whose value lies farther from it than the sample does (field.h: a preview's distance); a
// certified cell writes CERTIFIED only over pixels below it. So one schedule serves a picture that
// starts from nothing (a restart), from a snapped picture (every pixel REUSED), and from a moving one
// (only the new strips need work, and a zoom that outruns the pool shows coarse samples of the new
// view rather than copies stretched across it), and never paints a coarse block over a good pixel or
// a nearer preview.
//
// work_run is the job the thread pool runs. Everything it reads is either the WorkContext, frozen
// for the batch, or the field's front store and grid, of which it writes only its own tile's
// pixels and the corner records of its own lattice points (the planner gives each tile at most one
// unit per batch, and only lets a tile certify a level once its neighbours finished that level in
// an earlier batch).
#ifndef MANDELBROT_WORK_H
#define MANDELBROT_WORK_H

#include "certify.h"
#include "field.h"
#include "kernel.h"
#include "tiles.h"

// A pixel's coarse level is one of the schedule's levels, or TILE_LEVELS for a pixel no block
// improves (work_coarse_level).
#define WORK_COARSE_LEVELS (TILE_LEVELS + 1)

// What every unit of one batch shares, frozen for the batch's lifetime.
typedef struct {
	Field* field;
	KernelLimits limits;
	double mu_offset;  // log2(log2(bailout^2)): phi = log2(log2 |z_n|^2) - mu_offset
	bool certify;      // this run computes corners and certifies cells
	CertifyContext cells;

	// At depth: the reference the samples perturb against, and the field's origin less its point,
	// so a sample's delta_c is that plus its grid offsets (deep.h).
	bool perturbed;
	KernelReference reference;
	double origin_re, origin_im;
} WorkContext;

typedef struct {
	const WorkContext* context;
	int tile;
	int first_step, last_step;  // the steps [first_step, last_step) of the tile's schedule
	int64_t predicted;          // the planner's estimate of its cost, in lane-iterations

	// Written by the job.
	int64_t work;              // what the kernel spent, in plain lane-iterations (kernel.h)
	int samples;               // samples iterated
	int certified;             // pixels filled by certificates
	int32_t need[TILE_STEPS];  // what each step of the tile needs: before, as the planner saw it;
	                           // afterwards, what is left (certificates settle pixels the planner
	                           // had priced as work)
	int32_t coarse[WORK_COARSE_LEVELS];  // the tile's pixels by coarse level (work_coarse_level):
	                                     // before, as the planner saw them; afterwards, as they are
} WorkUnit;

/** Builds the table of block distances. Call once from _init, before any worker starts. */
void work_init(void);

/** Freezes a batch's shared context for the field's current spec. */
void work_context_init(WorkContext* c, Field* field, bool certify, const CertifyContext* cells);

/** Sets the reference a perturbed field's samples iterate against (deep_ready's answer). */
void work_context_reference(WorkContext* c, const KernelReference* reference, double origin_re, double origin_im);

/**
 * The fraction phi of an orbit that escaped with |z_n|^2 = radius_sq: mu = n - phi, continuous as
 * c crosses from one n to the next. Clamped to [0, 1].
 */
float work_escape_fraction(const WorkContext* c, double radius_sq);

/** Runs one WorkUnit (passed as the job's argument). The pool's job function. */
void work_run(void* unit);

/**
 * Counts, step by step of the tile's schedule, the pixels of `store` that still need work; and the
 * tile's pixels by coarse level (work_coarse_level).
 */
void work_count_needs(const PixelStore* store, int tile, int32_t* need, int32_t* coarse);

/**
 * A pixel's COARSE LEVEL: the coarsest level of the schedule whose block would bring it a nearer
 * value (a block's far corner lies (stride - 1) sqrt 2 from its sample) — 0 for an empty pixel, or a
 * copy from across the screen; TILE_LEVELS for one that no block improves, its value already at its
 * own coordinate. A tile's previews coarser than level L are those whose coarse level is L or less.
 */
int work_coarse_level(uint8_t state, uint16_t distance);

/** A tile's previews coarser than `level`, from its pixels by coarse level. */
static inline int32_t work_coarser_than(const int32_t* coarse, int level) {
	int32_t coarser = 0;
	for (int l = 0; l <= level; l++) coarser += coarse[l];
	return coarser;
}

#endif // MANDELBROT_WORK_H
