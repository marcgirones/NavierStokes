#!/usr/bin/env python3
"""Validate the time-dependent Taylor-Green vortex against its exact solution.

Runs the C++ solver for several N values at the same final time, reads u_final.csv
and v_final.csv, compares them at their staggered locations, and saves field
comparisons, divergence checks, and a mesh-convergence plot.

Expected executable interface for NS_modified.cpp: NS N tFinal nu dt
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

L = 1.0
RHO_DEFAULT = 1.0
NU_DEFAULT = 0.01  # Re = Umax*L/nu ~= 100 for Umax=L=1


def exact_fields(n, time, nu, rho=RHO_DEFAULT):
    """Return exact u,v,p sampled on the staggered U/V and cell-centred P grids."""
    h = L / n
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="xy")
    # U nodes: x=(i+1)h, y=(j+1/2)h; V nodes: x=(i+1/2)h, y=(j+1)h.
    xu, yu = (i + 1.0) * h, (j + 0.5) * h
    xv, yv = (i + 0.5) * h, (j + 1.0) * h
    xp, yp = (i + 0.5) * h, (j + 0.5) * h
    f = np.exp(-8.0 * np.pi**2 * nu * time)
    u = f * np.cos(2.0 * np.pi * xu) * np.sin(2.0 * np.pi * yu)
    v = -f * np.sin(2.0 * np.pi * xv) * np.cos(2.0 * np.pi * yv)
    p = -rho * f**2 * (0.5 * np.cos(2.0 * np.pi * xp)**2
                         + 0.5 * np.cos(2.0 * np.pi * yp)**2)
    return u, v, p


def run_case(exe, n, nsteps, dt, nu, rho, workdir):
    case_dir = workdir / f"N{n}"
    case_dir.mkdir(parents=True, exist_ok=True)

    # NS_modified.cpp expects: N, tFinal, nu, dt.
    # The older script interface (N, nSteps, dt, nu, rho) was incompatible.
    t_final = nsteps * dt
    cmd = [
        str(exe.resolve()),
        str(n),
        f"{t_final:.17g}",
        f"{nu:.17g}",
        f"{dt:.17g}",
    ]
    result = subprocess.run(cmd, cwd=case_dir, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(
            f"Solver failed for N={n}.\\nCommand: {' '.join(cmd)}\\n"
            f"stdout:\\n{result.stdout}\\nstderr:\\n{result.stderr}"
        )

    fields = {}
    for name in ("u_final.csv", "v_final.csv", "pressure_final.csv"):
        path = case_dir / name
        if not path.exists():
            raise FileNotFoundError(f"Solver did not create {path}")
        arr = np.loadtxt(path, delimiter=",")
        if arr.shape != (n, n):
            raise ValueError(f"{path}: expected {(n, n)}, got {arr.shape}")
        fields[name] = arr
    return case_dir, fields, result.stdout


def divergence(u, v, h):
    # Periodic finite-volume divergence at pressure-cell centres.
    return (u - np.roll(u, 1, axis=1)) / h + (v - np.roll(v, 1, axis=0)) / h


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exe", default="./NS", help="path to compiled C++ solver")
    ap.add_argument("--sizes", type=int, nargs="+", default=[16, 32, 64],
                    help="mesh sizes N (square N x N grids)")
    ap.add_argument("--time", type=float, default=0.01, help="target final time")
    ap.add_argument("--dt", type=float, default=1.0e-4, help="fixed time step")
    ap.add_argument("--nu", type=float, default=NU_DEFAULT, help="kinematic viscosity")
    ap.add_argument("--rho", type=float, default=RHO_DEFAULT, help="density")
    ap.add_argument("--workdir", default="validation_runs", help="output directory")
    ap.add_argument("--out", default="validation_results", help="plots/output folder")
    args = ap.parse_args()

    if args.time < 0 or args.dt <= 0 or args.nu < 0 or args.rho <= 0:
        ap.error("Require time>=0, dt>0, nu>=0 and rho>0")
    sizes = sorted(set(args.sizes))
    if any(n < 2 for n in sizes):
        ap.error("All N values must be >= 2")
    nsteps = int(round(args.time / args.dt))
    if nsteps < 0:
        ap.error("Invalid number of time steps")
    actual_time = nsteps * args.dt
    if not np.isclose(actual_time, args.time, rtol=0.0, atol=1e-10 * max(1.0, args.time)):
        print(f"Note: target time {args.time:g} adjusted to {actual_time:g} = nSteps*dt")

    exe = Path(args.exe).expanduser().resolve()
    if not exe.exists():
        ap.error(f"Executable not found: {exe}. Compile the C++ solver first.")
    workdir = Path(args.workdir).expanduser().resolve()
    outdir = Path(args.out).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    errors = []
    div_results = []
    results = []
    for n in sizes:
        print(f"Running N={n}, steps={nsteps}, t={actual_time:g} ...", flush=True)
        case_dir, fields, _ = run_case(exe, n, nsteps, args.dt, args.nu, args.rho, workdir)
        un = fields["u_final.csv"]
        vn = fields["v_final.csv"]
        pn = fields["pressure_final.csv"]
        ue, ve, pe = exact_fields(n, actual_time, args.nu, args.rho)
        eu = float(np.max(np.abs(un - ue)))
        ev = float(np.max(np.abs(vn - ve)))
        div = divergence(un, vn, L / n)
        divmax = float(np.max(np.abs(div)))
        divrms = float(np.sqrt(np.mean(div**2)))
        # Pressure is only defined up to a constant for periodic flow.
        dpn = pn - pn[0, 0]
        dpe = pe - pe[0, 0]
        ep = float(np.max(np.abs(dpn - dpe)))
        errors.append((eu, ev, ep))
        div_results.append((divmax, divrms))
        results.append((n, case_dir, un, vn, ue, ve))
        print(f"  max|u-uexact|={eu:.4e}; max|v-vexact|={ev:.4e}; "
              f"max|p-pref - exact|={ep:.4e}; maxDiv={divmax:.4e}; RMSDiv={divrms:.4e}")

        # For each mesh: exact and numerical velocity fields, plus absolute error.
        fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
        panels = [
            (axes[0, 0], ue, "Exact u"), (axes[0, 1], un, "Numerical u"),
            (axes[0, 2], np.abs(un - ue), "|u error|"),
            (axes[1, 0], ve, "Exact v"), (axes[1, 1], vn, "Numerical v"),
            (axes[1, 2], np.abs(vn - ve), "|v error|"),
        ]
        for ax, field, title in panels:
            im = ax.imshow(field, origin="lower", extent=(0, L, 0, L),
                           aspect="equal", interpolation="nearest")
            ax.set_title(title)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            fig.colorbar(im, ax=ax, shrink=0.82)
        fig.suptitle(f"Taylor–Green vortex, N={n}, t={actual_time:g}, nu={args.nu:g}")
        fig.savefig(outdir / f"fields_N{n}.png", dpi=160)
        plt.close(fig)

    err = np.asarray(errors)
    divs = np.asarray(div_results)
    print("\n N       h          max error u    max error v    max error p*    max div       RMS div")
    for k, n in enumerate(sizes):
        print(f"{n:3d}  {L/n:9.3e}  {err[k,0]:12.4e}  {err[k,1]:12.4e}  "
              f"{err[k,2]:12.4e}  {divs[k,0]:12.4e}  {divs[k,1]:12.4e}")

    hvals = np.array([L / n for n in sizes])
    fig, ax = plt.subplots(figsize=(7, 5))
    for col, label, marker in [(0, "u max error", "o"), (1, "v max error", "s")]:
        ax.loglog(hvals, err[:, col], marker + "-", label=label)
    if len(sizes) > 1:
        ref = err[0, 0] * (hvals / hvals[0])**2
        ax.loglog(hvals, ref, "--", label=r"$O(h^2)$ reference")
    ax.set_xlabel("h = L/N")
    ax.set_ylabel("Maximum absolute error")
    ax.set_title(f"Velocity mesh convergence at t={actual_time:g}")
    ax.grid(True, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "mesh_convergence.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(sizes, divs[:, 0], "o-", label="maximum |div u|")
    ax.loglog(sizes, divs[:, 1], "s-", label="RMS div u")
    ax.set_xlabel("N")
    ax.set_ylabel("Divergence")
    ax.set_title("Incompressibility check after projection")
    ax.grid(True, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "divergence_check.png", dpi=160)
    plt.close(fig)

    summary = outdir / "validation_summary.csv"
    np.savetxt(summary, np.column_stack((sizes, hvals, err, divs)), delimiter=",")
    print(f"\nPlots and summary saved in: {outdir}")
    print("Files: fields_N*.png, mesh_convergence.png, divergence_check.png, validation_summary.csv")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, FileNotFoundError, ValueError) as exc:
        sys.exit(str(exc))
