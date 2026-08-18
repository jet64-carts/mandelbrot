// carts/mandelbrot/reuse.c
#pragma STDC FP_CONTRACT OFF

#include "reuse.h"

#include "work.h"

_Static_assert(MB_SCREEN_HEIGHT <= MB_SCREEN_WIDTH, "an AxisMap holds a row map in its first entries");

// The whole number nearest x, as a double. |x| is a count of pixels between two points of the view's
// plane, which the view limits keep far inside an int64_t.
static double nearest_whole(double x) {
	return (double)(int64_t)(x < 0.0 ? x - 0.5 : x + 0.5);
}

// The distance between two coordinates.
static double distance_between(double a, double b) {
	return a > b ? a - b : b - a;
}

// A PAN's matching. A window reaches half a pixel either side, and the grid's tolerance beyond
// (field.h): a drag by a whole number of pixels and a half puts every old column exactly half a
// pixel from two new ones, and which side of the half the rounding lands must not decide whether the
// column is carried. So neighbouring windows may share a boundary coordinate; the sweep consumes
// every old column it matches, and the first window to reach it — its lower neighbour — takes it.
static void match_pan(const GridAxis* from, const GridAxis* to, int count, double sign, double reach, AxisMap* map) {
	int k = 0;
	for (int i = 0; i < count; i++) {
		const double t = sign * to->ideal[i];
		// The old coordinates are increasing too; skip those below this window for good, then take
		// the nearest inside it (the lower of two equally near).
		while (k < count && sign * from->actual[k] < t - reach) k++;
		int best = REUSE_MISSING;
		double best_distance = reach;
		for (int q = k; q < count && sign * from->actual[q] <= t + reach; q++) {
			const double distance = distance_between(sign * from->actual[q], t);
			if (distance < best_distance || (best == REUSE_MISSING && distance <= reach)) {
				best = q;
				best_distance = distance;
			}
		}
		map->source[i] = (int16_t)best;
		if (best != REUSE_MISSING) k = best + 1;
	}
}

// A ZOOM's matching, by dynamic programming over the ideal columns in order: the most old columns
// carried, and among those matchings the least total displacement. After ideal i the state is the
// last old column the matching has used: FREE (below ideal i's window, so any column of the next
// window may follow it) or the k-th column of i's window. Each match gains 1 less its displacement
// over (count + 1) reaches, so one more match always outweighs any saving in displacement.
static int16_t s_window_first[MB_SCREEN_WIDTH];
static int8_t s_window_size[MB_SCREEN_WIDTH];
static int8_t s_previous[MB_SCREEN_WIDTH][ZOOM_STATES];  // each state's state after the ideal before
static bool s_took[MB_SCREEN_WIDTH][ZOOM_STATES];        // whether reaching it matched this ideal
#define ZOOM_FREE 0
#define ZOOM_UNREACHED (-1.0)

static void match_zoom(const GridAxis* from, const GridAxis* to, int count, double sign, double reach, AxisMap* map) {
	double score[ZOOM_STATES], next[ZOOM_STATES];
	for (int k = 0; k < ZOOM_STATES; k++) score[k] = ZOOM_UNREACHED;
	score[ZOOM_FREE] = 0.0;
	int previous_first = 0, previous_size = 0;
	int first = 0;
	const double per_match = 1.0 / ((double)(count + 1) * reach);
	for (int i = 0; i < count; i++) {
		const double t = sign * to->ideal[i];
		while (first < count && sign * from->actual[first] < t - reach) first++;
		int size = 0;
		while (size < ZOOM_CHOICES && first + size < count && sign * from->actual[first + size] <= t + reach) size++;
		s_window_first[i] = (int16_t)first;
		s_window_size[i] = (int8_t)size;
		for (int k = 0; k <= size; k++) next[k] = ZOOM_UNREACHED;
		for (int p = 0; p <= previous_size; p++) {
			if (score[p] < 0.0) continue;
			// The last column used, as a state of this window: FREE if it lies below it.
			const int last = p == ZOOM_FREE ? -1 : previous_first + p - 1;
			const int k_last = last < first ? ZOOM_FREE : last - first + 1;
			if (score[p] > next[k_last]) {
				next[k_last] = score[p];
				s_previous[i][k_last] = (int8_t)p;
				s_took[i][k_last] = false;
			}
			for (int k = k_last + 1; k <= size; k++) {
				const double gain = 1.0 - distance_between(sign * from->actual[first + k - 1], t) * per_match;
				if (score[p] + gain > next[k]) {
					next[k] = score[p] + gain;
					s_previous[i][k] = (int8_t)p;
					s_took[i][k] = true;
				}
			}
		}
		for (int k = 0; k <= size; k++) score[k] = next[k];
		previous_first = first;
		previous_size = size;
	}
	int state = ZOOM_FREE;
	for (int k = 1; k <= previous_size; k++) {
		if (score[k] > score[state]) state = k;
	}
	for (int i = count - 1; i >= 0; i--) {
		map->source[i] = s_took[i][state] ? (int16_t)(s_window_first[i] + state - 1) : REUSE_MISSING;
		state = s_previous[i][state];
	}
}

// Matches one axis. `sign` is +1 where coordinates increase along the axis (columns) and -1 where
// they decrease (rows), so both matchings run over increasing values. A matched column keeps its
// old actual coordinate; the result is monotone and one-to-one.
//
// A missing column shows one of its nearest carried neighbours until it is computed — the one before
// it or the one after it, pixel by pixel, whichever holds the nearer value (preview_pixel) — each
// side's gap the distance between their coordinates (a preview's distance, field.h). Where it is
// computed
// depends on the move. A PAN (the spacing unchanged) computes it at the lattice point nearest its
// ideal, less than half a pixel away, so the axis stays one lattice and a view at rest can move onto
// it (camera_rest_on_grid) instead of recomputing it; and since every lattice point is reckoned
// afresh from one anchor, no rounding gathers however long the pan. Placed at its ideal instead, a
// strip computed while the grid stood half a pixel off would seam against it, and every seam,
// carried on across the screen, would lose a column each time the pan came back to that phase. A
// ZOOM has no lattice to keep: a missing column is computed at its ideal, and the first carried
// column anchors the lattice the next pan extends.
static int match_axis(const GridAxis* from, GridAxis* to, int count, double sign, double spacing, double tolerance,
                      AxisMap* map) {
	const bool pan = from->lattice_spacing == spacing;
	if (pan) {
		match_pan(from, to, count, sign, 0.5 * spacing + tolerance, map);
	} else {
		match_zoom(from, to, count, sign, ZOOM_REACH_PIXELS * spacing + tolerance, map);
	}
	int matched = 0, first = REUSE_MISSING;
	for (int i = 0; i < count; i++) {
		if (map->source[i] == REUSE_MISSING) continue;
		to->actual[i] = from->actual[map->source[i]];
		if (first == REUSE_MISSING) first = i;
		matched++;
	}

	to->lattice = pan ? from->lattice : first != REUSE_MISSING ? to->actual[first] : to->ideal[0];
	to->lattice_spacing = spacing;
	int last = REUSE_MISSING;
	for (int i = 0; i < count; i++) {
		if (map->source[i] != REUSE_MISSING) last = i;
		map->fill[REUSE_BEFORE][i] = (int16_t)last;  // the nearest carried new column at or before i, for
		                                             // now: the sweep back turns it into its old one
	}
	int next = REUSE_MISSING;
	for (int i = count - 1; i >= 0; i--) {
		if (map->source[i] != REUSE_MISSING) next = i;
		if (map->source[i] == REUSE_MISSING) {
			to->actual[i] = pan ? to->lattice + nearest_whole((to->ideal[i] - to->lattice) / spacing) * spacing
			                    : to->ideal[i];
		}
		const int nearest[REUSE_SIDES] = {map->fill[REUSE_BEFORE][i], next};
		for (int side = 0; side < REUSE_SIDES; side++) {
			const int n = nearest[side];
			map->fill[side][i] = n == REUSE_MISSING ? REUSE_MISSING : map->source[n];
			map->gap[side][i] = n == REUSE_MISSING
			                        ? DISTANCE_FAR
			                        : distance_of_pixels(distance_between(to->actual[i], to->actual[n]) / spacing);
		}
		const double off = to->actual[i] - to->ideal[i];
		to->at_ideal[i] = off <= tolerance && off >= -tolerance;
	}
	map->matched = matched;
	return matched;
}

double reuse_match(Field* f, AxisMap* columns, AxisMap* rows) {
	const double spacing = view_spacing(&f->spec.view);
	const double tolerance = field_ideal_tolerance(f);
	const int c = match_axis(&f->columns, &f->next_columns, MB_SCREEN_WIDTH, 1.0, spacing, tolerance, columns);
	const int r = match_axis(&f->rows, &f->next_rows, MB_SCREEN_HEIGHT, -1.0, spacing, tolerance, rows);
	const double share_c = (double)c / MB_SCREEN_WIDTH;
	const double share_r = (double)r / MB_SCREEN_HEIGHT;
	return share_c < share_r ? share_c : share_r;
}

// An axis on its ideal lattice, anchored at its first coordinate.
static void anchor_ideal(GridAxis* next, double spacing) {
	next->lattice = next->ideal[0];
	next->lattice_spacing = spacing;
}

static void ideal_axis(GridAxis* next, int count, double spacing) {
	for (int i = 0; i < count; i++) {
		next->actual[i] = next->ideal[i];
		next->at_ideal[i] = true;
	}
	anchor_ideal(next, spacing);
}

void reuse_hold_maps(AxisMap* columns, AxisMap* rows) {
	for (int i = 0; i < MB_SCREEN_WIDTH; i++) {
		columns->source[i] = rows->source[i] = REUSE_MISSING;
		for (int side = 0; side < REUSE_SIDES; side++) {
			columns->fill[side][i] = rows->fill[side][i] = (int16_t)i;
			columns->gap[side][i] = rows->gap[side][i] = DISTANCE_FAR;
		}
	}
	columns->matched = rows->matched = 0;
}

// The old actual coordinates either side of each new ideal one, by one sweep: both run monotonically.
// (Beyond either end of the old axis, both are its end column.)
static void nearest_axis(const GridAxis* from, const GridAxis* to, int count, double sign, double spacing,
                         AxisMap* map) {
	int k = 0;
	for (int i = 0; i < count; i++) {
		const double t = sign * to->ideal[i];
		while (k + 1 < count && sign * from->actual[k + 1] <= t) k++;
		const int around[REUSE_SIDES] = {k, k + 1 < count && sign * from->actual[k] <= t ? k + 1 : k};
		map->source[i] = REUSE_MISSING;
		for (int side = 0; side < REUSE_SIDES; side++) {
			map->fill[side][i] = (int16_t)around[side];
			map->gap[side][i] = distance_of_pixels(distance_between(to->ideal[i], from->actual[around[side]]) / spacing);
		}
	}
	map->matched = 0;
}

void reuse_warp_maps(const Field* f, AxisMap* columns, AxisMap* rows) {
	const double spacing = view_spacing(&f->spec.view);
	nearest_axis(&f->columns, &f->next_columns, MB_SCREEN_WIDTH, 1.0, spacing, columns);
	nearest_axis(&f->rows, &f->next_rows, MB_SCREEN_HEIGHT, -1.0, spacing, rows);
}

void reuse_ideal_grid(Field* f) {
	const double spacing = view_spacing(&f->spec.view);
	ideal_axis(&f->next_columns, MB_SCREEN_WIDTH, spacing);
	ideal_axis(&f->next_rows, MB_SCREEN_HEIGHT, spacing);
}

// A snap moves every column that is off its ideal coordinate onto it, and leaves the rest where
// they are: a column within the tolerance is ideal already, its pixels exact at its own coordinate,
// and moving it by a few ulps would make every one of them a value for somewhere else.
static void snap_axis(const GridAxis* from, GridAxis* to, int count, double spacing) {
	for (int i = 0; i < count; i++) {
		to->actual[i] = from->at_ideal[i] ? from->actual[i] : to->ideal[i];
		to->at_ideal[i] = true;
	}
	anchor_ideal(to, spacing);
}

void reuse_snap_grid(Field* f) {
	const double spacing = view_spacing(&f->spec.view);
	snap_axis(&f->columns, &f->next_columns, MB_SCREEN_WIDTH, spacing);
	snap_axis(&f->rows, &f->next_rows, MB_SCREEN_HEIGHT, spacing);
}

static void copy_axis(const GridAxis* from, GridAxis* to, int count) {
	for (int i = 0; i < count; i++) {
		to->actual[i] = from->actual[i];
		to->at_ideal[i] = from->at_ideal[i];
	}
	to->lattice = from->lattice;
	to->lattice_spacing = from->lattice_spacing;
}

void reuse_keep_grid(Field* f) {
	copy_axis(&f->columns, &f->next_columns, MB_SCREEN_WIDTH);
	copy_axis(&f->rows, &f->next_rows, MB_SCREEN_HEIGHT);
}

// The iteration limit's rules, applied to one pixel. An interior pixel's n is the limit, whatever
// proved it interior. Escaped at n with a limit now at or below n, the oracle would call it
// interior: it becomes so, exactly, with no iteration. Interior only because the old limit ran out,
// a higher one might let it escape: its value stays but its rank drops to REUSED, so it is iterated
// again (a shape, a cycle or an interior certificate is interior under any limit). A certified
// exterior value is interpolated, its n uncertain by one, so near a lower limit it is only a
// preview.
static void apply_limit(int32_t* n, float* phi, uint8_t* state, int64_t old_limit, int64_t new_limit) {
	if (new_limit == old_limit) return;
	const Fate fate = state_fate(*state);
	const Rank rank = state_rank(*state);
	if (fate != FATE_ESCAPED) {
		*n = (int32_t)new_limit;
		if (fate == FATE_LIMIT && new_limit > old_limit && rank > RANK_REUSED) {
			*state = state_with_rank(*state, RANK_REUSED);
		}
		return;
	}
	if (new_limit > old_limit) return;
	if (rank == RANK_CERTIFIED) {
		if (*n >= new_limit - 1) *state = state_with_rank(*state, RANK_PREVIEW);
		return;
	}
	if (*n >= new_limit) {
		*n = (int32_t)new_limit;
		*phi = 0.0f;
		*state = (uint8_t)(state_make(rank, FATE_LIMIT) | (*state & STATE_OFFSET));
	}
}

// How far an old pixel's value lies from its own coordinate: a preview's distance; a pixel's width
// for a value a snap had moved (OFFSET) — up to half a pixel from its column already, and up to a
// whole one once carried again; nothing for any other.
static uint16_t own_distance(const PixelStore* in, size_t q) {
	const uint8_t state = in->state[q];
	if (state_rank(state) <= RANK_PREVIEW) return in->distance[q];
	return (state & STATE_OFFSET) ? DISTANCE_PIXEL : 0;
}

// One pixel of a gather, carried: old pixel q's value at its own coordinate, and its state — but a
// value a snap had moved (OFFSET) goes back to being a preview, and across the crossover an exact or
// certified one is only reused; its distance seen on the new grid. Most of every gather, so kept to
// a few loads and stores.
static inline void carry_pixel(const ReconcileContext* c, size_t q, PixelStore* out, size_t p) {
	const PixelStore* in = &c->field->front;
	uint8_t state = in->state[q];
	const Rank rank = state_rank(state);
	uint16_t distance = 0;
	if (rank <= RANK_PREVIEW) {
		distance = distance_carried(in->distance[q], c->distance_scale, 0, 0);
	} else if (state & STATE_OFFSET) {
		distance = distance_carried(DISTANCE_PIXEL, c->distance_scale, 0, 0);
		state = (uint8_t)(state_with_rank(state, RANK_PREVIEW) & ~STATE_OFFSET);
	} else if (c->demote && rank > RANK_REUSED) {
		state = state_with_rank(state, RANK_REUSED);
	}
	out->n[p] = in->n[q];
	out->phi[p] = in->phi[q];
	out->state[p] = state;
	out->distance[p] = distance;
}

// One pixel of a gather that its maps do not carry, or of a kept restart: a PREVIEW of whichever of
// its nearest carried neighbours (before and after it, along each axis) holds the value nearest it —
// a preview's distance being its source's, seen on the new grid, and the gap between the two pixels.
// A carried column (row) offers one candidate, the same on both sides.
static void preview_pixel(const ReconcileContext* c, int x, int y, PixelStore* out, size_t p) {
	const PixelStore* in = &c->field->front;
	const int sides_x = c->columns->fill[REUSE_AFTER][x] == c->columns->fill[REUSE_BEFORE][x] ? 1 : REUSE_SIDES;
	const int sides_y = c->rows->fill[REUSE_AFTER][y] == c->rows->fill[REUSE_BEFORE][y] ? 1 : REUSE_SIDES;
	bool found = false;
	size_t q = 0;
	uint16_t distance = DISTANCE_FAR;
	for (int side_x = 0; side_x < sides_x; side_x++) {
		const int sx = c->columns->fill[side_x][x];
		if (sx == REUSE_MISSING) continue;
		for (int side_y = 0; side_y < sides_y; side_y++) {
			const int sy = c->rows->fill[side_y][y];
			if (sy == REUSE_MISSING) continue;
			const size_t candidate = (size_t)sy * MB_SCREEN_WIDTH + (size_t)sx;
			const uint16_t d = distance_carried(own_distance(in, candidate), c->distance_scale,
			                                    c->columns->gap[side_x][x], c->rows->gap[side_y][y]);
			if (!found || d < distance) {
				q = candidate;
				distance = d;
				found = true;
			}
		}
	}
	if (!found) {
		out->state[p] = state_make(RANK_EMPTY, FATE_ESCAPED);
		out->distance[p] = DISTANCE_FAR;
		return;
	}
	const uint8_t state = in->state[q];
	const Rank rank = state_rank(state);
	out->n[p] = in->n[q];
	out->phi[p] = in->phi[q];
	out->state[p] = (uint8_t)(state_with_rank(state, rank < RANK_PREVIEW ? rank : RANK_PREVIEW) & ~STATE_OFFSET);
	out->distance[p] = distance;
}

// A snap: a pixel whose column or row moves to its ideal coordinate keeps its value, now up to
// half a pixel from where it belongs, as REUSED | OFFSET — unless it was below REUSED anyway.
static void snap_pixel(const Field* f, int x, int y, size_t p) {
	if (f->columns.at_ideal[x] && f->rows.at_ideal[y]) return;
	uint8_t* state = &f->front.state[p];
	if (state_rank(*state) >= RANK_CERTIFIED) *state = (uint8_t)(state_with_rank(*state, RANK_REUSED) | STATE_OFFSET);
}

void reconcile_run(void* arg) {
	ReconcileUnit* unit = (ReconcileUnit*)arg;
	const ReconcileContext* c = unit->context;
	Field* f = c->field;
	const TileOrigin origin = tiles_origin(unit->tile);
	const bool through_maps = c->kind == RECONCILE_GATHER || (c->kind == RECONCILE_RESTART && c->keep);
	PixelStore* out = through_maps ? &f->back : &f->front;

	for (int ty = 0; ty < MB_TILE_HEIGHT; ty++) {
		const int y = origin.y + ty;
		const int row_source = through_maps ? c->rows->source[y] : REUSE_MISSING;
		for (int tx = 0; tx < MB_TILE_WIDTH; tx++) {
			const int x = origin.x + tx;
			const size_t p = (size_t)y * MB_SCREEN_WIDTH + (size_t)x;
			switch (c->kind) {
			case RECONCILE_RESTART:
				if (c->keep) {
					preview_pixel(c, x, y, out, p);  // every map entry is missing
				} else {
					out->state[p] = state_make(RANK_EMPTY, FATE_ESCAPED);
					out->distance[p] = DISTANCE_FAR;
				}
				break;
			case RECONCILE_GATHER: {
					const int column_source = c->columns->source[x];
					if (row_source != REUSE_MISSING && column_source != REUSE_MISSING) {
						carry_pixel(c, (size_t)row_source * MB_SCREEN_WIDTH + (size_t)column_source, out, p);
					} else {
						preview_pixel(c, x, y, out, p);
					}
					apply_limit(&out->n[p], &out->phi[p], &out->state[p], c->old_limit, c->new_limit);
					break;
				}
			case RECONCILE_LIMIT:
				apply_limit(&out->n[p], &out->phi[p], &out->state[p], c->old_limit, c->new_limit);
				break;
			case RECONCILE_SNAP:
				snap_pixel(f, x, y, p);
				break;
			}
		}
	}
	work_count_needs(out, unit->tile, unit->need, unit->coarse);
}
