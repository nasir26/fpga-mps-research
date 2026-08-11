#ifndef HYDROGEN_RADIAL_KERNEL_H
#define HYDROGEN_RADIAL_KERNEL_H

extern "C" void hydrogen_radial_kernel(
        int    n,
        int    l,
        double r_max,
        int    npts,
        double* r_grid,
        double* u_wave,
        double* diagnostics);

#endif
