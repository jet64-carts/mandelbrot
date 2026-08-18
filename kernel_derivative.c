// carts/mandelbrot/kernel_derivative.c
//
// The derivative lane kernel: kernel_lanes.h's machinery specialised to carry dz/dc beside every
// orbit, for the coarse samples that become certificate corners (certify.h). Its state is twice the
// plain kernel's, so it keeps fewer groups in flight (kernel.h).
#pragma STDC FP_CONTRACT OFF

#include "kernel.h"

#define LANES_SLOTS KERNEL_DERIVATIVE_SLOTS
#include "kernel_lanes.h"

int64_t kernel_iterate_derivative(const KernelPoint* points, KernelOrbit* orbits, int count,
                                  const KernelLimits* limits) {
	return kernel_drive(points, orbits, count, limits, true, false, NULL).lanes;
}
