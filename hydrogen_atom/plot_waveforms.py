"""
Author: Nasir Ali
Organization: C-DAC Noida
Date: August 2026

Plots the radial "waveforms" produced by hydrogen_testbench (csim of
hydrogen_radial_kernel.cpp) for a set of hydrogen orbitals:
  - u(r) = r*R(r), the quantity the Numerov kernel actually propagates
  - R(r), the physical radial wavefunction
  - r^2 |R(r)|^2, the radial probability density (probability per dr
    of finding the electron at radius r, after integrating out angles)

Run hydrogen_testbench first to (re)generate data/*.csv.

Usage: python3 plot_waveforms.py
"""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
OUT_DIR = os.path.join(os.path.dirname(__file__), "graphs")

STATES = ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4d", "4f"]
COLORS = plt.cm.viridis(np.linspace(0.05, 0.9, len(STATES)))


def load(state):
    path = os.path.join(DATA_DIR, f"{state}.csv")
    with open(path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    r = np.array([float(row["r"]) for row in rows])
    u = np.array([float(row["u"]) for row in rows])
    R = np.array([float(row["R"]) for row in rows])
    P = np.array([float(row["radial_prob_density"]) for row in rows])
    return r, u, R, P


def plot_overview():
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    for ax, ylabel, key, title in zip(
        axes,
        ["u(r) = r R(r)  (Bohr radii^-1/2)", "R(r)  (Bohr radii^-3/2)", "r^2 |R(r)|^2  (probability / a0)"],
        ["u", "R", "P"],
        ["Numerov kernel output u(r)", "Radial wavefunction R(r)", "Radial probability density"],
    ):
        for state, color in zip(STATES, COLORS):
            r, u, R, P = load(state)
            y = {"u": u, "R": R, "P": P}[key]
            r_cut = min(len(r), np.searchsorted(r, 8.0 * (int(state[0]) ** 2) / 2.0 + 25))
            ax.plot(r[:r_cut], y[:r_cut], label=state, color=color, linewidth=1.8)
        ax.set_xlabel("r (Bohr radii, a0)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.axhline(0, color="0.6", linewidth=0.6)
        ax.set_xlim(left=0)
    axes[0].legend(ncol=3, fontsize=9, framealpha=0.9)
    fig.suptitle("Hydrogen atom radial Schrodinger solutions (Numerov shooting kernel, csim)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, "hydrogen_radial_overview.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def plot_per_state():
    os.makedirs(OUT_DIR, exist_ok=True)
    for state in STATES:
        r, u, R, P = load(state)
        n = int(state[0])
        r_max_plot = min(r[-1], 8.0 * n * n)
        cut = np.searchsorted(r, r_max_plot)

        fig, ax1 = plt.subplots(figsize=(7, 4.5))
        ax1.plot(r[:cut], u[:cut], color="tab:blue", linewidth=2, label="u(r) = r R(r)  [kernel output]")
        ax1.axhline(0, color="0.7", linewidth=0.6)
        ax1.set_xlabel("r (Bohr radii, a0)")
        ax1.set_ylabel("u(r)", color="tab:blue")
        ax1.tick_params(axis="y", labelcolor="tab:blue")

        ax2 = ax1.twinx()
        ax2.fill_between(r[:cut], P[:cut], color="tab:orange", alpha=0.25)
        ax2.plot(r[:cut], P[:cut], color="tab:orange", linewidth=1.6, label="r^2|R(r)|^2  [radial density]")
        ax2.set_ylabel("radial probability density", color="tab:orange")
        ax2.tick_params(axis="y", labelcolor="tab:orange")
        ax2.set_ylim(bottom=0)

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=8)

        ax1.set_title(f"Hydrogen {state} orbital -- radial waveform")
        fig.tight_layout()
        out_path = os.path.join(OUT_DIR, f"hydrogen_{state}.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"wrote {out_path}")


if __name__ == "__main__":
    plot_overview()
    plot_per_state()
