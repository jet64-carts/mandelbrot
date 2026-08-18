// carts/mandelbrot/tiles.c
#include "tiles.h"

static TileStep s_steps[TILE_STEPS];
static int s_samples_before[TILE_STEPS + 1];
static int s_level_start[TILE_LEVELS + 1];
static uint8_t s_pixel_step[MB_TILE_HEIGHT][MB_TILE_WIDTH];
static int s_order[TILE_COUNT];

// Squared distance from a tile's centre to the screen's, in doubled pixels so it stays integral.
static int64_t distance_from_centre_sq(int tile) {
	const TileOrigin o = tiles_origin(tile);
	const int64_t dx = 2 * (int64_t)o.x + MB_TILE_WIDTH - MB_SCREEN_WIDTH;
	const int64_t dy = 2 * (int64_t)o.y + MB_TILE_HEIGHT - MB_SCREEN_HEIGHT;
	return dx * dx + dy * dy;
}

void tiles_init(void) {
	int step = 0;
	for (int level = 0; level < TILE_LEVELS; level++) {
		s_level_start[level] = step;
		const int stride = tiles_level_stride(level);
		for (int row = 0; row < MB_TILE_HEIGHT; row += stride) {
			// Level 0 samples its whole lattice. Every later level samples only what is new: on a
			// row of the previous (twice as coarse) lattice, the odd multiples of its stride; on a
			// row between them, every multiple.
			const bool coarse_row = level > 0 && row % (2 * stride) == 0;
			TileStep* s = &s_steps[step];
			s->level = (uint8_t)level;
			s->stride = (uint8_t)stride;
			s->row = (uint16_t)row;
			s->first_column = (uint16_t)(coarse_row ? stride : 0);
			s->column_step = (uint16_t)(coarse_row ? 2 * stride : stride);
			s->samples = (uint16_t)((MB_TILE_WIDTH - s->first_column + s->column_step - 1) / s->column_step);
			for (int column = s->first_column; column < MB_TILE_WIDTH; column += s->column_step) {
				s_pixel_step[row][column] = (uint8_t)step;
			}
			s_samples_before[step + 1] = s_samples_before[step] + s->samples;
			step++;
		}
	}
	s_level_start[TILE_LEVELS] = step;

	// Nearest the centre first; equal distances keep their row-major order, so the order is a
	// pure function of the geometry.
	for (int t = 0; t < TILE_COUNT; t++) {
		int at = t;
		const int64_t d = distance_from_centre_sq(t);
		while (at > 0 && distance_from_centre_sq(s_order[at - 1]) > d) {
			s_order[at] = s_order[at - 1];
			at--;
		}
		s_order[at] = t;
	}
}

const TileStep* tiles_step(int step) {
	return &s_steps[step];
}

int tiles_samples_before(int steps) {
	return s_samples_before[steps];
}

int tiles_level_start(int level) {
	return s_level_start[level];
}

int tiles_level_stride(int level) {
	return TILE_COARSEST_STRIDE >> level;
}

const uint8_t* tiles_pixel_steps(int y) {
	return s_pixel_step[y];
}

TileOrigin tiles_origin(int tile) {
	TileOrigin o;
	o.x = (tile % TILE_COLUMNS) * MB_TILE_WIDTH;
	o.y = (tile / TILE_COLUMNS) * MB_TILE_HEIGHT;
	return o;
}

const int* tiles_order(void) {
	return s_order;
}
