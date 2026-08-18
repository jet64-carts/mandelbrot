// carts/mandelbrot/certify.h
//
// THE CERTIFIED SAMPLER — whole cells of pixels decided without iterating a single one of them.
//
// Two theorems of Koebe type give every sample c a disc about it that provably contains no point
// of the set's boundary:
//
//   EXTERIOR. With G the escape potential, the distance d from c to the set satisfies
//       sinh(G) / (2 e^G |G'|)  <  d  <  2 sinh(G) / |G'|,
//   and from an escaped orbit, G/|G'| = |z_n| ln|z_n| / |z_n'| to within O(1/|z_n|^2). So c owns an
//   M-free disc of radius (b/4) (1 - e^(-2G)) / (2G), with b = 2 |z_n| ln|z_n| / |z_n'| the usual
//   distance estimate. (The factor tends to 1 near the set, and matters in the far field.)
//   INTERIOR. For c inside a hyperbolic component of period p, z0 on its attracting cycle and the
//   derivatives of f^p there, b_int = (1 - |dz|^2) / |dcdz + dzdz dc / (1 - dz)|, and the distance
//   to the component's boundary lies in (b_int / 4, b_int]: c owns an interior disc of b_int / 4.
//
// Every point of an s x s cell lies within s/sqrt 2 of its nearest corner, so when all four corners
// own discs of radius at least s/sqrt 2 of the same kind, the discs cover the cell and the whole
// cell is exterior (or interior): no filament, however fine, can hide inside it. An interior cell
// is filled as interior; an exterior cell's smooth count is interpolated by bicubic Hermite
// interpolation from its corners' mu and gradient, which is smooth there because the escape
// potential is harmonic off the set (Harnack). The margins over the theorem are measured
// (design.md §7): the interior needs none; the exterior's buys interpolation accuracy.
//
// The corners are the coarse levels' samples (strides 16, 8 and 4), computed with the derivative
// kernel; the cells are the refinement schedule's own lattice squares, so the schedule's quadtree
// is the certificate's quadtree: a cell that passes skips every finer sample inside it.
#ifndef MANDELBROT_CERTIFY_H
#define MANDELBROT_CERTIFY_H

#include "field.h"
#include "kernel.h"

// Corners live on the finest certified lattice; the levels coarser than it supply them.
#define CERT_LATTICE_STRIDE 4
#define CERT_LEVELS 3  // levels 0, 1, 2: strides 16, 8, 4
#define CERT_CORNER_COLUMNS (MB_SCREEN_WIDTH / CERT_LATTICE_STRIDE)
#define CERT_CORNER_ROWS (MB_SCREEN_HEIGHT / CERT_LATTICE_STRIDE)
#define CERT_CORNERS (CERT_CORNER_COLUMNS * CERT_CORNER_ROWS)

typedef enum {
	CORNER_NONE = 0,  // no certificate: a cell touching it is computed
	CORNER_EXTERIOR,
	CORNER_INTERIOR,
} CornerKind;

// What one corner sample knows.
typedef struct {
	double mu;         // the smooth iteration count (exterior)
	float gx, gy;      // its gradient, per pixel, in screen axes (exterior)
	float radius;      // the certified disc's radius, in pixels
	uint8_t kind;      // CornerKind
	uint32_t run;      // the compute run that wrote it; a record from another run is stale
} Corner;

// What a certification pass needs: the grid's scale, the limit, the margins, and the run.
typedef struct {
	double spacing;          // the pixel spacing h
	int64_t max_iterations;
	double exterior_margin;  // multiples of the theorem's minimum radius, s / sqrt 2
	double interior_margin;
	bool exterior;           // which certificates are in use
	bool interior;
	uint32_t run;
	Corner* corners;         // CERT_CORNERS records, indexed by certify_corner_index
} CertifyContext;

/** The record index of lattice point (x, y), both multiples of CERT_LATTICE_STRIDE. */
static inline int certify_corner_index(int x, int y) {
	return (y / CERT_LATTICE_STRIDE) * CERT_CORNER_COLUMNS + x / CERT_LATTICE_STRIDE;
}

/**
 * The corner record of a sample at c whose orbit `o` came from a derivative kernel, with mu its
 * smooth count. Handles every fate: escaped (the exterior certificate), cycle and shape (the
 * interior certificate), and limit (no certificate). At depth (`reference` not NULL) (c_re, c_im)
 * is the sample's delta_c, and its cycle is followed against the reference in perturbation form,
 * from the state the kernel caught it at (design.md §9.7), and — when `at_limit` — an orbit that
 * reached the limit against a nucleus reference is given Newton's method for its cycle. Adds to
 * `steps` the orbit steps the interior certificate took (the caller prices them:
 * CERTIFY_STEP_COST_PERCENT).
 */
Corner certify_corner(double c_re, double c_im, const KernelOrbit* o, double mu, double spacing,
                      const KernelReference* reference, bool at_limit, int64_t* steps);

// What one of those steps costs, in percent of a plain lane-iteration (the scheduler's unit): a
// scalar step of the orbit and its derivatives, 3.8 ns against the plain kernel's 0.47 ns a
// lane-iteration (tests/mandelbrot/deep_test.c --bench). At depth, against a nucleus of period P,
// a corner that reached the limit takes about 3.2 P of them — at p4128 eighteen times its own orbit
// through the BLA table — so a certifying batch unpriced for them ran twenty times its quantum
// (design.md §9.7).
#define CERTIFY_STEP_COST_PERCENT 820

/**
 * Tests every cell of side `stride` in the tile at (tile_x, tile_y) and fills the ones that pass:
 * each pixel of a certified cell whose value is not yet certified or exact becomes CERTIFIED.
 * Writes only the tile's own pixels (a cell fills [x, x + s) x [y, y + s)); reads corner records,
 * which may belong to the right, lower and lower-right tiles — the caller guarantees those are
 * final. Returns the pixels filled.
 */
int certify_cells(const CertifyContext* c, PixelStore* store, int tile_x, int tile_y, int stride);

#endif // MANDELBROT_CERTIFY_H
