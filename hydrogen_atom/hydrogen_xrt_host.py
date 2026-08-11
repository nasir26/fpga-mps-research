"""
Author: Nasir Ali
Organization: C-DAC Noida
Date: August 2026

XRT host runner for hydrogen_radial_kernel.xclbin on the Alveo U55C.
Mirrors the pyxrt buffer-object / kernel-launch pattern used by
fpga_simulator_v07.py's QuantumFPGACircuit for the statevector kernel.

Falls back to nothing automatically -- pyxrt + a built xclbin are both
required. Build the xclbin first with build_hydrogen/run_build.sh on a
machine with Vitis 2023.2 + XRT installed.

Usage:
    python3 hydrogen_xrt_host.py --n 2 --l 1 --xclbin build_hydrogen/hydrogen_radial_kernel.xclbin
"""
from __future__ import annotations
import argparse
import os
import numpy as np

try:
    import pyxrt as xrt
    XRT_AVAILABLE = True
except ImportError:
    XRT_AVAILABLE = False

KERNEL_NAME = "hydrogen_radial_kernel"


class HydrogenFPGA:
    def __init__(self, xclbin_path: str, device_index: int = 0):
        if not XRT_AVAILABLE:
            raise RuntimeError("pyxrt unavailable -- run on a machine with XRT installed.")
        if not os.path.exists(xclbin_path):
            raise FileNotFoundError(f"XCLBIN not found: {xclbin_path}")

        self.device = xrt.device(device_index)
        self.uuid = self.device.load_xclbin(xclbin_path)
        self.kernel = xrt.kernel(self.device, self.uuid, KERNEL_NAME)

    def solve(self, n: int, l: int, r_max: float, npts: int):
        if not (0 <= l < n):
            raise ValueError(f"require 0 <= l < n, got n={n}, l={l}")

        nbytes = npts * 8
        # Kernel arg order: n, l, r_max, npts, r_grid, u_wave, diagnostics
        # -> group_id(4)=r_grid, group_id(5)=u_wave, group_id(6)=diagnostics
        r_bo = xrt.bo(self.device, nbytes, xrt.bo.flags.normal, self.kernel.group_id(4))
        u_bo = xrt.bo(self.device, nbytes, xrt.bo.flags.normal, self.kernel.group_id(5))
        d_bo = xrt.bo(self.device, 2 * 8, xrt.bo.flags.normal, self.kernel.group_id(6))

        run = self.kernel(n, l, r_max, npts, r_bo, u_bo, d_bo)
        run.wait()

        r_bo.sync(xrt.xclBOSyncDirection.XCL_BO_SYNC_BO_FROM_DEVICE, nbytes, 0)
        u_bo.sync(xrt.xclBOSyncDirection.XCL_BO_SYNC_BO_FROM_DEVICE, nbytes, 0)
        d_bo.sync(xrt.xclBOSyncDirection.XCL_BO_SYNC_BO_FROM_DEVICE, 2 * 8, 0)

        r_grid = np.frombuffer(r_bo.read(nbytes, 0), dtype=np.float64).copy()
        u_wave = np.frombuffer(u_bo.read(nbytes, 0), dtype=np.float64).copy()
        diagnostics = np.frombuffer(d_bo.read(2 * 8, 0), dtype=np.float64).copy()

        return r_grid, u_wave, diagnostics[0], diagnostics[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1)
    ap.add_argument("--l", type=int, default=0)
    ap.add_argument("--npts", type=int, default=6000)
    ap.add_argument("--xclbin", type=str, default="build_hydrogen/hydrogen_radial_kernel.xclbin")
    args = ap.parse_args()

    r_max = 8.0 * args.n * args.n + 20.0
    fpga = HydrogenFPGA(args.xclbin)
    r_grid, u_wave, E, r_expect = fpga.solve(args.n, args.l, r_max, args.npts)

    E_exact = -1.0 / (2.0 * args.n ** 2)
    r_exact = 0.5 * (3.0 * args.n ** 2 - args.l * (args.l + 1))
    print(f"n={args.n} l={args.l}  E_fpga={E:.8f} Ha  (exact {E_exact:.8f})")
    print(f"                <r>_fpga={r_expect:.6f} a0  (exact {r_exact:.6f})")


if __name__ == "__main__":
    main()
