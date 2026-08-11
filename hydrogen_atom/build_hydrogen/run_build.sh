#!/bin/bash
set -e
set -o pipefail
cd "$(dirname "$0")"
source /tools/Xilinx/Vitis/2023.2/settings64.sh
if [ -f /opt/xilinx/xrt/setup.sh ]; then source /opt/xilinx/xrt/setup.sh; fi

PLATFORM=/opt/xilinx/platforms/xilinx_u55c_gen3x16_xdma_3_202210_1/xilinx_u55c_gen3x16_xdma_3_202210_1.xpfm
KNAME=hydrogen_radial_kernel
SRC=../hydrogen_radial_kernel.cpp

echo "[$(date)] === COMPILE (v++ -c) ==="
v++ -c -t hw \
    --platform "$PLATFORM" \
    -k "$KNAME" \
    --config hydrogen_kernel.cfg \
    -o "${KNAME}.xo" \
    "$SRC"

echo "[$(date)] === LINK (v++ -l) ==="
v++ -l -t hw \
    --platform "$PLATFORM" \
    --config hydrogen_kernel.cfg \
    -o "${KNAME}.xclbin" \
    "${KNAME}.xo"

echo "[$(date)] === BUILD COMPLETE ==="
ls -la "${KNAME}.xclbin"
