// carts/mandelbrot/tiles.h
//
// TILES AND THE REFINEMENT SCHEDULE — the order in which a picture is computed.
//
// The screen is cut into MB_TILE_WIDTH x MB_TILE_HEIGHT tiles, and each tile is computed
// coarse-to-fine in TILE_LEVELS levels of halving STRIDE: level 0 samples every 16th pixel in
// both directions and fills each 16 x 16 block with its sample's value, level 1 samples the new
// pixels of the stride-8 lattice and fills 8 x 8 blocks, and so on down to stride 1, where every
// pixel is exact. A level never recomputes a pixel an earlier level sampled: on a row of the
// coarser lattice it samples only the odd multiples of its stride, and on a row between them it
// samples every multiple. So the whole schedule computes each pixel exactly once.
//
// A tile's schedule is a fixed sequence of STEPS, one per sample row of each level in order:
// the three rows of level 0, then the six of level 1, and so on. A tile's progress is simply how
// many of its steps are done, and a unit of work is a contiguous run of one tile's steps. The
// step table is the same for every tile, built once by tiles_init.
#ifndef MANDELBROT_TILES_H
#define MANDELBROT_TILES_H

#include "config.h"
#include "jet64.h"

#define TILE_LEVELS 5
#define TILE_COARSEST_STRIDE (1 << (TILE_LEVELS - 1))  // 16
#define TILE_COLUMNS (MB_SCREEN_WIDTH / MB_TILE_WIDTH)
#define TILE_ROWS (MB_SCREEN_HEIGHT / MB_TILE_HEIGHT)
#define TILE_COUNT (TILE_COLUMNS * TILE_ROWS)

// Steps in a tile's schedule: MB_TILE_HEIGHT / stride summed over the strides 16, 8, 4, 2, 1,
// which is MB_TILE_HEIGHT (2^LEVELS - 1) / 2^(LEVELS - 1).
#define TILE_STEPS (MB_TILE_HEIGHT * ((1 << TILE_LEVELS) - 1) / TILE_COARSEST_STRIDE)
// The steps of level 0, which a restart computes before the first frame is shown.
#define TILE_COARSE_STEPS (MB_TILE_HEIGHT / TILE_COARSEST_STRIDE)

_Static_assert(MB_TILE_WIDTH % TILE_COARSEST_STRIDE == 0, "a tile must hold whole coarse blocks");
_Static_assert(MB_TILE_HEIGHT % TILE_COARSEST_STRIDE == 0, "a tile must hold whole coarse blocks");
_Static_assert(TILE_STEPS <= UINT8_MAX + 1, "a step index must fit the pixel-to-step table's bytes");

// One step: one sample row of one level.
typedef struct {
	uint8_t level;         // 0 (coarsest) .. TILE_LEVELS - 1 (every pixel)
	uint8_t stride;        // the level's lattice spacing, which is also its fill block's side
	uint16_t row;          // the sample row's pixel row within the tile
	uint16_t first_column; // the first sampled pixel column within the tile
	uint16_t column_step;  // the distance between samples along the row
	uint16_t samples;      // how many samples the row holds
} TileStep;

// The pixel rectangle of one tile.
typedef struct {
	int x, y;  // top-left pixel
} TileOrigin;

/** Builds the step table. Call once from _init, before any worker starts. */
void tiles_init(void);

/** Step `step` of every tile's schedule, 0 <= step < TILE_STEPS. */
const TileStep* tiles_step(int step);

/** Samples in the first `steps` steps of a tile's schedule: its progress in pixels. */
int tiles_samples_before(int steps);

/** The first step of `level`; level TILE_LEVELS gives TILE_STEPS, the schedule's end. */
int tiles_level_start(int level);

/** The lattice spacing of `level`: 16, 8, 4, 2, 1. */
int tiles_level_stride(int level);

/**
 * The step that samples each pixel of tile-local row y, indexed by column: every pixel belongs to
 * exactly one step. A row at a time, since the counts of what a tile still needs (work.c) read
 * every pixel of it after every job.
 */
const uint8_t* tiles_pixel_steps(int y);

/** Where tile `tile` sits on the screen. */
TileOrigin tiles_origin(int tile);

/**
 * The tiles in the order the planner visits them: nearest the screen's centre first, so a
 * picture sharpens where the eye already is.
 */
const int* tiles_order(void);

#endif // MANDELBROT_TILES_H
