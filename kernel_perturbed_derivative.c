// carts/mandelbrot/kernel_perturbed_derivative.c
//
// The perturbation lane kernel carrying dz/dc beside every orbit, for the certificate corners at
// depth (certify.h): kernel_lanes.h's machinery with both flags set.
#pragma STDC FP_CONTRACT OFF

#include "kernel.h"

#define LANES_SLOTS KERNEL_PERTURBED_DERIVATIVE_SLOTS
#include "kernel_lanes.h"

KernelCount kernel_iterate_perturbed_derivative(const KernelPoint* points, KernelOrbit* orbits, int count,
                                            const KernelLimits* limits, const KernelReference* reference) {
	return kernel_drive(points, orbits, count, limits, true, true, reference);
}
