// carts/mandelbrot/field.c
#pragma STDC FP_CONTRACT OFF

#include "field.h"

#define PIXEL_STORE_BYTES \
	((uint32_t)(MB_SCREEN_PIXELS * (sizeof(int32_t) + sizeof(float) + sizeof(uint8_t) + sizeof(uint16_t))))
#define AXIS_BYTES(count) ((uint32_t)((count) * (2 * sizeof(double) + sizeof(uint8_t))))
#define GRID_BYTES (AXIS_BYTES(MB_SCREEN_WIDTH) + AXIS_BYTES(MB_SCREEN_HEIGHT))
#define STORES 2  // front and back
#define GRIDS 2   // the grid and the next grid

#define IDEAL_TOLERANCE_PIXELS 0x1p-20
#define IDEAL_TOLERANCE_RELATIVE 0x1p-48  // sixteen ulps of a double's 2^-52
#define SMALLEST_MAGNITUDE 0.25           // no coordinate's ulp is judged finer than a quarter's

// The bump allocator rounds every block up to 8 bytes, so each reservation below is 8-byte
// aligned whatever came before it; the wider types are reserved first only for tidiness.
static bool reserve_store(PixelStore* s) {
	s->n = (int32_t*)jet_memory_reserve((uint32_t)(MB_SCREEN_PIXELS * sizeof(int32_t)));
	s->phi = (float*)jet_memory_reserve((uint32_t)(MB_SCREEN_PIXELS * sizeof(float)));
	s->distance = (uint16_t*)jet_memory_reserve((uint32_t)(MB_SCREEN_PIXELS * sizeof(uint16_t)));
	s->state = (uint8_t*)jet_memory_reserve((uint32_t)MB_SCREEN_PIXELS);
	return s->n && s->phi && s->distance && s->state;
}

static bool reserve_axis(GridAxis* a, int count) {
	a->actual = (double*)jet_memory_reserve((uint32_t)(count * sizeof(double)));
	a->ideal = (double*)jet_memory_reserve((uint32_t)(count * sizeof(double)));
	a->at_ideal = (uint8_t*)jet_memory_reserve((uint32_t)count);
	a->lattice = a->lattice_spacing = 0.0;  // no spacing is zero: the first reconcile sets both
	return a->actual && a->ideal && a->at_ideal;
}

uint32_t field_bytes(void) {
	return STORES * PIXEL_STORE_BYTES + GRIDS * GRID_BYTES;
}

bool field_init(Field* f) {
	// Nothing is cleared: the first spec's restart writes every pixel of the front store before
	// anything reads them (tests/mandelbrot/render_test.c proves it).
	f->has_spec = false;
	f->origin = point_make(0.0, 0.0);
	return reserve_store(&f->front) && reserve_store(&f->back) &&
	       reserve_axis(&f->columns, MB_SCREEN_WIDTH) && reserve_axis(&f->rows, MB_SCREEN_HEIGHT) &&
	       reserve_axis(&f->next_columns, MB_SCREEN_WIDTH) && reserve_axis(&f->next_rows, MB_SCREEN_HEIGHT);
}

void field_ideal_grid(const FieldSpec* spec, const DeepPoint* origin, double* columns, double* rows) {
	// Each pixel samples the plane at its own centre (view.h).
	const double h = view_spacing(&spec->view);
	double centre_re, centre_im;
	point_difference(&spec->view.centre, origin, &centre_re, &centre_im);
	for (int i = 0; i < MB_SCREEN_WIDTH; i++) {
		columns[i] = centre_re + view_offset_x((double)i + VIEW_PIXEL_CENTRE) * h;
	}
	for (int j = 0; j < MB_SCREEN_HEIGHT; j++) {
		rows[j] = centre_im + view_offset_y((double)j + VIEW_PIXEL_CENTRE) * h;
	}
}

void field_plan_spec(Field* f, const FieldSpec* spec) {
	field_ideal_grid(spec, &f->origin, f->next_columns.ideal, f->next_rows.ideal);
	f->spec = *spec;
	f->has_spec = true;
}

static void swap_axis(GridAxis* a, GridAxis* b) {
	const GridAxis t = *a;
	*a = *b;
	*b = t;
}

void field_swap(Field* f, bool pixels) {
	swap_axis(&f->columns, &f->next_columns);
	swap_axis(&f->rows, &f->next_rows);
	if (pixels) {
		const PixelStore t = f->front;
		f->front = f->back;
		f->back = t;
	}
}

bool field_matches(const Field* f, const FieldSpec* spec) {
	return f->has_spec && f->spec.max_iterations == spec->max_iterations && f->spec.epoch == spec->epoch &&
	       view_equal(&f->spec.view, &spec->view);
}

// Shifts one axis's coordinates (actual, ideal and its lattice) by -d.
static void shift_axis(GridAxis* a, int count, double d) {
	for (int i = 0; i < count; i++) {
		a->actual[i] -= d;
		a->ideal[i] -= d;
	}
	a->lattice -= d;
}

void field_move_origin(Field* f, double d_re, double d_im) {
	point_offset(&f->origin, d_re, d_im);
	shift_axis(&f->columns, MB_SCREEN_WIDTH, d_re);
	shift_axis(&f->rows, MB_SCREEN_HEIGHT, d_im);
}

void field_set_origin(Field* f, const DeepPoint* origin) {
	double d_re, d_im;
	point_difference(origin, &f->origin, &d_re, &d_im);
	f->origin = *origin;
	shift_axis(&f->columns, MB_SCREEN_WIDTH, d_re);
	shift_axis(&f->rows, MB_SCREEN_HEIGHT, d_im);
}

double field_ideal_tolerance(const Field* f) {
	const FieldSpec* spec = &f->spec;
	// The magnitude the coordinates are held at: the centre's offset from the origin, and in direct
	// mode (origin zero, coordinates absolute) never finer than a quarter's.
	double re, im;
	point_difference(&spec->view.centre, &f->origin, &re, &im);
	const double floor = spec->perturbed ? 0.0 : SMALLEST_MAGNITUDE;
	const double magnitude = (re < 0.0 ? -re : re) + (im < 0.0 ? -im : im) + floor;
	const double pixels = IDEAL_TOLERANCE_PIXELS * view_spacing(&spec->view);
	const double rounding = IDEAL_TOLERANCE_RELATIVE * magnitude;
	return pixels > rounding ? pixels : rounding;
}

// The common offset of one axis, when it has one: every coordinate within half the tolerance of the
// first one's offset, so that after the view moves by that offset each lies within the tolerance.
static bool axis_offset(const GridAxis* a, int count, double tolerance, double* offset) {
	const double first = a->actual[0] - a->ideal[0];
	for (int i = 1; i < count; i++) {
		const double spread = a->actual[i] - a->ideal[i] - first;
		if (!(spread <= 0.5 * tolerance && spread >= -0.5 * tolerance)) return false;
	}
	*offset = first;
	return true;
}

bool field_lattice_offset(const Field* f, double* d_re, double* d_im) {
	if (!f->has_spec || field_grid_is_ideal(f)) return false;
	const double tolerance = field_ideal_tolerance(f);
	return axis_offset(&f->columns, MB_SCREEN_WIDTH, tolerance, d_re) &&
	       axis_offset(&f->rows, MB_SCREEN_HEIGHT, tolerance, d_im);
}

bool field_grid_is_ideal(const Field* f) {
	for (int i = 0; i < MB_SCREEN_WIDTH; i++) {
		if (!f->columns.at_ideal[i]) return false;
	}
	for (int j = 0; j < MB_SCREEN_HEIGHT; j++) {
		if (!f->rows.at_ideal[j]) return false;
	}
	return true;
}
