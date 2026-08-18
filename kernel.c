// carts/mandelbrot/kernel.c
//
// The plain lane kernel: kernel_lanes.h's machinery specialised without the derivative. The escape
// test is `<=` against the bailout on purpose: a NaN or an infinity, which an unchecked block
// produces once a lane has escaped and gone on squaring, compares false — so a block that ran into
// one simply fails and is replayed from its saved state.
#pragma STDC FP_CONTRACT OFF

#include "kernel.h"

#define LANES_SLOTS KERNEL_SLOTS
#include "kernel_lanes.h"

int64_t kernel_iterate(const KernelPoint* points, KernelOrbit* orbits, int count, const KernelLimits* limits) {
	return kernel_drive(points, orbits, count, limits, false, false, NULL).lanes;
}
