#!/usr/bin/env python3
"""
Plots the local divergence of the velocity field before (u*) and after (u^{n+1})
the Poisson/FFT projection, using the CSVs written by case_B.cpp.

Usage:
    python plot_divergence.py [dir_with_csvs] [Lx] [Ly]
Defaults: dir = results/case_B_results, Lx = Ly = 1.0
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def load(path):
    # each CSV row = one mesh row (Nx values), N rows in total
    return np.loadtxt(path, delimiter=",", ndmin=2)


def stats(div, dx, dy):
    return {
        "total": div.sum(),                      # sum of the local divergences
        "integral": div.sum() * dx * dy,         # ∫ div dA
        "max": np.abs(div).max(),                # L-infinity norm
        "l2": np.sqrt(np.mean(div ** 2)),        # RMS value
    }


def draw(ax, div, st, title, Lx, Ly):
    vmax = np.abs(div).max()
    vmax = vmax if vmax > 0 else 1.0
    im = ax.imshow(div, origin="lower", extent=[0, Lx, 0, Ly],
                   cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="equal")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(r"$\nabla\cdot\mathbf{u}$")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    txt = (f"Total divergence (Σ) = {st['total']:.3e}\n"
           f"∫ div dA = {st['integral']:.3e}\n"
           f"max |div| = {st['max']:.3e}\n"
           f"RMS div = {st['l2']:.3e}")
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))


def main():
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/case_B_run")
    out_folder = Path("results/case_B_results")
    Lx = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    Ly = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0

    div_up = load(folder / "upFieldDivergence.csv")   # arbitrary initial field
    div_un = load(folder / "unFieldDivergence.csv")   # field after projection

    Ny, Nx = div_up.shape
    dx, dy = Lx / Nx, Ly / Ny

    st_up = stats(div_up, dx, dy)
    st_un = stats(div_un, dx, dy)

    fig, axs = plt.subplots(1, 2, figsize=(13, 5.5))
    draw(axs[0], div_up, st_up, "Initial field (before projection)", Lx, Ly)
    draw(axs[1], div_un, st_un, r"Field at $n+1$ (after projection)", Lx, Ly)
    fig.suptitle(f"Local divergence of the velocity field  ({Nx}×{Ny} mesh)")
    fig.tight_layout()
    fig.savefig(out_folder / "divergence_check.png", dpi=200)

    print(f"Before : total = {st_up['total']:.3e}, max|div| = {st_up['max']:.3e}")
    print(f"After  : total = {st_un['total']:.3e}, max|div| = {st_un['max']:.3e}")


if __name__ == "__main__":
    main()
