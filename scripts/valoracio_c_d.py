#!/usr/bin/env python3
"""
Convergence check for the convective / diffusive terms of the 2D periodic
Navier-Stokes FVM solver (staggered grid, Taylor-Green field).

For each mesh size N it:
  1. runs the C++ binary with N as command-line argument,
  2. reads the CSV it writes (columns: convU, convV, difU, difV),
  3. divides the integrated terms by the cell volume h^2,
  4. compares with the analytic terms at the *staggered* positions,
  5. plots the max error vs h (log-log) with an O(h^2) reference.

Assumed interface of the C++ code (adjust the constants below if different):
    ./NS <N>      ->  writes  datos<N>.csv   (or datos.csv) in the working dir
    CSV rows written with i as outer loop and j as inner loop (row = i*N + j)

Usage:
    python convergence_check.py --exe ./NS --sizes 10 50 100 150 200
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent   # project root

NU = 1.51e-5
LX = 1.0

# ---------------------------------------------------------------- analytic --
x, y = sp.symbols("x y")
ua = sp.cos(2 * sp.pi * x) * sp.sin(2 * sp.pi * y)
va = -sp.sin(2 * sp.pi * x) * sp.cos(2 * sp.pi * y)

conv_x = sp.diff(ua*ua, x) + sp.diff(ua*va, y)
conv_y = sp.diff(va*ua, x) + sp.diff(va*va, y)
dif_x = NU * (sp.diff(ua, x, 2) + sp.diff(ua, y, 2))
dif_y = NU * (sp.diff(va, x, 2) + sp.diff(va, y, 2))

f_convU = sp.lambdify((x, y), conv_x, "numpy")
f_convV = sp.lambdify((x, y), conv_y, "numpy")
f_difU = sp.lambdify((x, y), dif_x, "numpy")
f_difV = sp.lambdify((x, y), dif_y, "numpy")


# --------------------------------------------------------------- C++ runner --
def run_case(exe: Path, n: int, workdir: Path) -> np.ndarray:
    """Run the solver for mesh size n and return the CSV as an (n*n, 4) array."""
    out_dir = ROOT / "results" / "case_A_run" 
    (out_dir / f"datos{n}.csv").unlink(missing_ok=True)

    res = subprocess.run([str(exe.resolve()), str(n)], cwd=ROOT,
                         capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"Solver failed for N={n}:\n{res.stderr}")
    path = out_dir / f"datos{n}.csv"
    if not path.exists():
        sys.exit(f"No CSV found in {out_dir} for N={n}")
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    if data.shape != (n * n, 4):
        sys.exit(f"{path}: expected shape {(n * n, 4)}, got {data.shape}")
    return data
    
# ------------------------------------------------------------------- errors --
def errors_for(n: int, data: np.ndarray) -> list:
    h = LX / n
    conv_u, conv_v, dif_u, dif_v = (data[:, c] / h**2 for c in range(4))

    # CSV row = i*N + j  (i outer, j inner)  ->  reshape to [i, j]
    conv_u, conv_v, dif_u, dif_v = (a.reshape(n, n).T for a in (conv_u, conv_v, dif_u, dif_v))

    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")  # 0-based
    # Staggered positions: u(i,j) at ((i+1)h, (j+1/2)h); v(i,j) at ((i+1/2)h, (j+1)h)
    xu, yu = (i + 1.0) * h, (j + 0.5) * h
    xv, yv = (i + 0.5) * h, (j + 1.0) * h

    return [
        np.abs(conv_u - f_convU(xu, yu)).max() * 2,
        np.abs(conv_v - f_convV(xv, yv)).max(),
        np.abs(dif_u - f_difU(xu, yu)).max()*2,
        np.abs(dif_v - f_difV(xv, yv)).max(),
    ]


# --------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="../case_A", help="path to the compiled C++ binary")
    ap.add_argument("--sizes", type=int, nargs="+", default=[4, 8, 16, 32, 64, 128, 256, 512, 1024])
    ap.add_argument("--workdir", default="../results/c_d_validation", help="directory where the binary is run")
    ap.add_argument("--out", default="../results/case_A_results/c_d_convergence.png")
    args = ap.parse_args()

    exe, workdir = Path(args.exe), Path(args.workdir)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    sizes = sorted(args.sizes)
    h = np.array([LX / n for n in sizes])
    err = np.zeros((len(sizes), 4))

    for k, n in enumerate(sizes):
        print(f"N = {n} ...", flush=True)
        err[k] = errors_for(n, run_case(exe, n, workdir))

    labels = ["convU", "convV", "difU", "difV"]
    print("\n   N         h   " + "  ".join(f"{l:>10}" for l in labels))
    for n, hh, e in zip(sizes, h, err):
        print(f"{n:4d}  {hh:8.4f}   " + "  ".join(f"{v:10.3e}" for v in e))
    if len(sizes) > 1:
        slopes = [np.polyfit(np.log(h), np.log(err[:, c]), 1)[0] for c in range(4)]
        print("\nFitted order (slope in log-log): " +
              ", ".join(f"{l}={s:.2f}" for l, s in zip(labels, slopes)))

    # --------------------------------------------------------------- plot --
    fig, ax = plt.subplots()
    for c, (l, m) in enumerate(zip(labels, "os^d")):
        ax.loglog(h, err[:, c], m + "-", lw=1.5, label=l)
    ref = err[0, 0] / h[0] ** 2 * h**2
    ax.loglog(h, ref, "--", lw=1.5, label=r"$O(h^2)$")
    ax.set_xlabel("h")
    ax.set_ylabel("Max error")
    ax.set_title("Convergence of the errors")
    ax.grid(True, which="both")
    ax.legend(loc="best")
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to {args.out}")


if __name__ == "__main__":
    main()
