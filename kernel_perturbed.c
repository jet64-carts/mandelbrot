// carts/mandelbrot/kernel_perturbed.c
//
// The perturbation lane kernel: kernel_lanes.h's machinery specialised to iterate each lane's
// offset from a reference orbit, rebasing it, with BLA steps where the reference's table allows.
// Its lanes hold three pairs of doubles where the plain kernel's hold two, so it keeps its own
// number of groups in flight (KERNEL_PERTURBED_SLOTS, kernel.h).
#pragma STDC FP_CONTRACT OFF

#include "kernel.h"

#define LANES_SLOTS KERNEL_PERTURBED_SLOTS
#include "kernel_lanes.h"

KernelCount kernel_iterate_perturbed(const KernelPoint* points, KernelOrbit* orbits, int count,
                                 const KernelLimits* limits, const KernelReference* reference) {
	return kernel_drive(points, orbits, count, limits, false, true, reference);
}
