/*
 * Author: Nasir Ali
 * Organization: C-DAC Noida
 * Date: August 2026
 *
 * Software C-simulation (csim) testbench for hydrogen_radial_kernel.cpp.
 * Calls the exact HLS-synthesizable kernel function as a plain C++ call
 * (no XRT/hardware needed) and checks its output against the closed-form
 * hydrogen-atom formulas:
 *
 *   E_nl   = -1 / (2 n^2)                     Hartree
 *   <r>_nl = (1/2) [ 3 n^2 - l(l+1) ]          Bohr radii
 *
 * Also dumps r,u(r),R(r),|R(r)|^2*r^2 per state to CSV for plot_waveforms.py.
 *
 * Build:  g++ -O2 -o hydrogen_testbench hydrogen_radial_kernel.cpp host_testbench.cpp
 * Run:    ./hydrogen_testbench
 */

#include "hydrogen_radial_kernel.h"
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <sys/stat.h>

struct OrbitalSpec { int n; int l; const char* name; };

int main() {
    const std::vector<OrbitalSpec> states = {
        {1, 0, "1s"}, {2, 0, "2s"}, {2, 1, "2p"},
        {3, 0, "3s"}, {3, 1, "3p"}, {3, 2, "3d"},
        {4, 0, "4s"}, {4, 2, "4d"}, {4, 3, "4f"},
    };

    const int NPTS = 6000;
    mkdir("data", 0755);

    printf("%-6s %-14s %-14s %-10s   %-14s %-14s %-10s\n",
           "state", "E_kernel(Ha)", "E_exact(Ha)", "rel.err",
           "<r>_kernel(a0)", "<r>_exact(a0)", "rel.err");
    printf("--------------------------------------------------------------------------------------------------\n");

    bool all_pass = true;

    for (const auto& s : states) {
        double r_max = 8.0 * s.n * s.n + 20.0;   // past classical turning
                                                   // point + exp-decay tail
        std::vector<double> r_grid(NPTS), u_wave(NPTS);
        double diagnostics[2];

        hydrogen_radial_kernel(s.n, s.l, r_max, NPTS,
                                r_grid.data(), u_wave.data(), diagnostics);

        double E_kernel = diagnostics[0];
        double r_kernel = diagnostics[1];

        double E_exact = -1.0 / (2.0 * s.n * s.n);
        double r_exact = 0.5 * (3.0 * s.n * s.n - s.l * (s.l + 1));

        double E_rel_err = std::fabs((E_kernel - E_exact) / E_exact);
        double r_rel_err = std::fabs((r_kernel - r_exact) / r_exact);

        printf("%-6s %-14.8f %-14.8f %-10.2e   %-14.6f %-14.6f %-10.2e\n",
               s.name, E_kernel, E_exact, E_rel_err,
               r_kernel, r_exact, r_rel_err);

        if (E_rel_err > 1e-4 || r_rel_err > 1e-3) all_pass = false;

        char path[256];
        std::snprintf(path, sizeof(path), "data/%s.csv", s.name);
        FILE* f = std::fopen(path, "w");
        std::fprintf(f, "r,u,R,radial_prob_density\n");
        for (int i = 0; i < NPTS; i++) {
            double r = r_grid[i];
            double u = u_wave[i];
            double R = u / r;                 // R(r) = u(r)/r
            double P = u * u;                 // r^2 |R|^2 = u^2, already
                                                // the radial probability
                                                // density r^2|R(r)|^2 dr
            std::fprintf(f, "%.10e,%.10e,%.10e,%.10e\n", r, u, R, P);
        }
        std::fclose(f);
    }

    printf("--------------------------------------------------------------------------------------------------\n");
    printf(all_pass
           ? "PASS: all states within tolerance (E rel.err < 1e-4, <r> rel.err < 1e-3)\n"
           : "FAIL: one or more states outside tolerance\n");

    return all_pass ? 0 : 1;
}
