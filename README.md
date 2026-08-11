# FPGA MPS Simulator — Research

Research work, paper drafts, and reference material accompanying the FPGA
quantum circuit simulator code in the sibling repo:
https://github.com/nasir26/fpga-mps-simulator

## `paper/`

Draft manuscript: **"Hardware-in-the-Loop Verification of Distance-3
Surface-Code Quantum Memory on a 28-Qubit HBM-Accelerated FPGA Statevector
Simulator"** (`main.tex`, compiled to `main.pdf`). Reports a distance-3
rotated surface-code quantum-memory experiment run on real Alveo U55C
hardware, cross-validated against `qiskit.quantum_info` and an independently
tracked classical Pauli frame. Results (150 trials / 600 FPGA-executed
syndrome rounds) are reproduced by `surface_code_qec.py` in the code repo;
`figures/results.csv` and `figures/fig_logical_error_rate.png` hold the
plotted data.

## `hydrogen_atom/`

Hydrogen-atom radial time-independent Schrodinger equation solver
implemented as an FPGA HLS kernel (`hydrogen_radial_kernel.cpp/.h`,
Vitis 2023.2, Alveo U55C), with a host testbench (`host_testbench.cpp`,
`hydrogen_xrt_host.py`) and reference wavefunction data/plots
(`data/`, `graphs/`) for radial states 1s through 4f. `build_hydrogen/`
holds the Vitis build config. The compiled `hydrogen_testbench` binary is
not tracked — rebuild it from the `.cpp` sources.

## `qiskit_nature/`

Qiskit Nature tutorial notebooks (01–12, official Qiskit Nature examples)
used as the reference source for the electronic-structure circuits that
`fpga_validate_nature_circuits.py` (in the code repo) reconstructs and runs
on the real FPGA MPS simulator, cross-checking against the exact/Aer values
each notebook reports.

## Reference documents

- `12.pdf` — Dennis, Kitaev, Landahl, Preskill, "Topological quantum
  memory" (arXiv:quant-ph/0110143) — source for the surface-code stabilizer
  formalism and error model used in `paper/` and `surface_code_qec.py`.
- `113.pdf` — supporting reference paper.
- `Bilateral_Grant_Proposal_CDAC_L7TR_22207_RH.docx` — CDAC grant proposal
  associated with this project (kept private).
