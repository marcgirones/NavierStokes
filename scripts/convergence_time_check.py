#!/usr/bin/env python3
"""Validate the time-dependent Taylor-Green vortex against its exact solution.

Runs the C++ solver once per N. The solver (case_C.cpp) writes EVERY saved
time step to history_u.bin / history_v.bin / history_p.bin (+ history_times.csv),
so no re-runs are needed to build time histories. The script compares each
snapshot with the exact solution at its staggered locations and saves:
  - field comparisons at the final time, divergence check, mesh convergence;
  1. Numerical vs exact u and v at a chosen node, vs time (one figure per N).
  2. Max absolute error of u and v vs time at every step (one line per N).
  3. Max absolute error of the relative pressure dp = p - p[0,0] vs time
     (one line per N).

Default viscosity gives Re = U*L/nu = 100 (U = 1, L = 1  ->  nu = 0.01).

Executable interface for case_C: case_C N tFinal nu [saveEvery]
"""
import argparse
import re
import subprocess
import time as clock
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

L = 1.0
RHO_DEFAULT = 1.0
RE_DEFAULT = 100.0
NU_DEFAULT = 1.0 * L / RE_DEFAULT   # Re = U L / nu = 100


def exact_fields(n, t, nu, rho=RHO_DEFAULT):
    """Return exact u,v,p sampled on the staggered U/V and cell-centred P grids."""
    h = L / n
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="xy")
    # U nodes: x=(i+1)h, y=(j+1/2)h; V nodes: x=(i+1/2)h, y=(j+1)h.
    xu, yu = (i + 1.0) * h, (j + 0.5) * h
    xv, yv = (i + 0.5) * h, (j + 1.0) * h
    xp, yp = (i + 0.5) * h, (j + 0.5) * h
    f = np.exp(-8.0 * np.pi**2 * nu * t)
    u = f * np.cos(2.0 * np.pi * xu) * np.sin(2.0 * np.pi * yu)
    v = -f * np.sin(2.0 * np.pi * xv) * np.cos(2.0 * np.pi * yv)
    p = -rho * f**2 * (0.5 * np.cos(2.0 * np.pi * xp)**2
                       + 0.5 * np.cos(2.0 * np.pi * yp)**2)
    return u, v, p


def run_case(exe, n, t_final, nu, case_dir, save_every=1):
    case_dir.mkdir(parents=True, exist_ok=True)

    # case_C expects: N, tFinal, nu, [saveEvery].
    cmd = [str(exe.resolve()), str(n), f"{t_final:.17g}", f"{nu:.17g}", str(save_every)]
    result = subprocess.run(cmd, cwd=case_dir, capture_output=True, text=True)

    if result.returncode:
        raise RuntimeError(
            f"Solver failed for N={n}.\nCommand: {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )

    # Snapshot times written by the solver: columns step, t, dt.
    times_path = case_dir / "history_times.csv"
    if not times_path.exists():
        raise FileNotFoundError(f"Solver did not create {times_path}")
    meta = np.loadtxt(times_path, delimiter=",", skiprows=1, ndmin=2)
    t = meta[:, 1]
    ns = len(t)

    # Raw float64 snapshots, each n*n values (row = j, column = i).
    hist = {}
    for key, name in (("u", "history_u.bin"), ("v", "history_v.bin"),
                      ("p", "history_p.bin")):
        path = case_dir / name
        if not path.exists():
            raise FileNotFoundError(f"Solver did not create {path}")
        if path.stat().st_size != ns * n * n * 8:
            raise ValueError(f"{path}: size does not match {ns} snapshots of {n}x{n}")
        hist[key] = np.memmap(path, dtype=np.float64, mode="r", shape=(ns, n, n))

    return t, hist, result.stdout


def divergence(u, v, h):
    # Periodic finite-volume divergence at pressure-cell centres.
    return (u - np.roll(u, 1, axis=1)) / h + (v - np.roll(v, 1, axis=0)) / h


def node_index(n, x, y, stagger):
    """Array index (row=j, col=i) of the grid node closest to (x, y).

    stagger = 'u': x=(i+1)h,   y=(j+1/2)h
    stagger = 'v': x=(i+1/2)h, y=(j+1)h
    """
    h = L / n
    if stagger == "u":
        i, j = round(x / h - 1.0), round(y / h - 0.5)
    else:
        i, j = round(x / h - 0.5), round(y / h - 1.0)
    return int(j) % n, int(i) % n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", default="case_C", help="path to compiled C++ solver")
    ap.add_argument("--sizes", type=int, nargs="+", default=[8, 16, 32, 64, 128],
                    help="mesh sizes N (square N x N grids)")
    ap.add_argument("--time", type=float, default=12.0, help="target final time")
    ap.add_argument("--nu", type=float, default=NU_DEFAULT, help="kinematic viscosity")
    ap.add_argument("--rho", type=float, default=RHO_DEFAULT, help="density")
    ap.add_argument("--workdir", default="results/validation_runs", help="output directory")
    ap.add_argument("--out", default="results/validation_results", help="plots/output folder")
    ap.add_argument("--save-every", type=int, default=1,
                    help="solver saves a snapshot every this many steps (1 = every step; "
                         "raise it for large N to save disk space)")
    ap.add_argument("--u-node", type=float, nargs=2, default=[0.5, 0.25],
                    metavar=("X", "Y"), help="physical point (0..1) whose nearest U node is tracked")
    ap.add_argument("--v-node", type=float, nargs=2, default=[0.25, 0.5],
                    metavar=("X", "Y"), help="physical point (0..1) whose nearest V node is tracked")
    ap.add_argument("--node-index", type=int, nargs=2, default=None, metavar=("I", "J"),
                    help="track U and V at the same 0-based array index (i=x, j=y); "
                         "overrides --u-node/--v-node")
    args = ap.parse_args()

    if args.time < 0 or args.nu < 0 or args.rho <= 0:
        ap.error("Require time>=0, nu>=0 and rho>0")
    if args.save_every < 1:
        ap.error("--save-every must be >= 1")

    sizes = sorted(set(args.sizes))
    if any(n < 2 for n in sizes):
        ap.error("All N values must be >= 2")

    exe = Path(args.exe).expanduser().resolve()
    if not exe.exists():
        ap.error(f"Executable not found: {exe}. Compile the C++ solver first.")

    workdir = Path(args.workdir).expanduser().resolve()
    outdir = Path(args.out).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    errors = []
    div_results = []
    history = {}  # n -> dict of time series

    for n in sizes:
        print(f"Running N={n}, target t={args.time:g} ...", flush=True)

        if args.node_index is not None:
            iu = iv = args.node_index[0] % n
            ju = jv = args.node_index[1] % n
        else:
            ju, iu = node_index(n, *args.u_node, "u")
            jv, iv = node_index(n, *args.v_node, "v")

        start = clock.perf_counter()
        t, snap, _ = run_case(exe, n, args.time, args.nu,
                              workdir / f"N{n}", args.save_every)
        elapsed = clock.perf_counter() - start
        ns = len(t)
        print(f"  Time computing: {elapsed:.2f} s; {ns} snapshots saved")

        hist = {k: np.empty(ns) for k in ("un", "ue", "vn", "ve", "eu", "ev", "ep")}
        for k in range(ns):
            un, vn, pn = np.asarray(snap["u"][k]), np.asarray(snap["v"][k]), np.asarray(snap["p"][k])
            ue, ve, pe = exact_fields(n, t[k], args.nu, args.rho)

            # Pressure is only defined up to a constant: dp_i = p_i - p_1.
            dpn = pn - pn[0, 0]
            dpe = pe - pe[0, 0]

            hist["un"][k], hist["ue"][k] = un[ju, iu], ue[ju, iu]
            hist["vn"][k], hist["ve"][k] = vn[jv, iv], ve[jv, iv]
            hist["eu"][k] = np.max(np.abs(un - ue))
            hist["ev"][k] = np.max(np.abs(vn - ve))
            hist["ep"][k] = np.max(np.abs(dpn - dpe))
        hist["t"] = t
        history[n] = hist

        # ---- Final-time diagnostics (last snapshot = final state) ----
        eu, ev, ep = hist["eu"][-1], hist["ev"][-1], hist["ep"][-1]
        actual_time = t[-1]
        div = divergence(un, vn, L / n)
        divmax = float(np.max(np.abs(div)))
        divrms = float(np.sqrt(np.mean(div**2)))

        errors.append((eu, ev, ep))
        div_results.append((divmax, divrms))

        print(f"  Actual final time: {actual_time:.12g}")
        print(f"  max|u-uexact|={eu:.4e}; max|v-vexact|={ev:.4e}; "
              f"max|dp - dp_exact|={ep:.4e}; maxDiv={divmax:.4e}; "
              f"RMSDiv={divrms:.4e}")

        # For each mesh: exact and numerical velocity fields, plus absolute error.
        fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
        panels = [
            (axes[0, 0], ue, "Exact u"),
            (axes[0, 1], un, "Numerical u"),
            (axes[0, 2], np.abs(un - ue), "|u error|"),
            (axes[1, 0], ve, "Exact v"),
            (axes[1, 1], vn, "Numerical v"),
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

        # ---- Feature 1: node time series (numerical vs exact u and v) ----
        h = L / n
        xu_n, yu_n = (iu + 1.0) * h, (ju + 0.5) * h
        xv_n, yv_n = (iv + 0.5) * h, (jv + 1.0) * h
        H = history[n]
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(H["t"], H["un"], color="red", label="$u_n$")
        ax.plot(H["t"], H["vn"], color="lime", label="$v_n$")
        ax.plot(H["t"], H["ue"], color="blue", label="$u_a$")
        ax.plot(H["t"], H["ve"], color="black", label="$v_a$")
        ax.set_xlabel("time (s)")
        ax.set_ylabel("velocity at node (m/s)")
        ax.set_title(f"N={n}: u at (x={xu_n:.4f}, y={yu_n:.4f}), "
                     f"v at (x={xv_n:.4f}, y={yv_n:.4f})", fontsize=10)
        ax.set_xlim(0, H["t"][-1])
        ax.grid(True)
        ax.legend(loc="upper right")
        fig.tight_layout()
        fig.savefig(outdir / f"node_timeseries_N{n}.png", dpi=160)
        plt.close(fig)

    err = np.asarray(errors)
    divs = np.asarray(div_results)

    print("\n N       h          max error u    max error v    max error dp    "
          "max div       RMS div")
    for k, n in enumerate(sizes):
        print(f"{n:3d}  {L/n:9.3e}  {err[k,0]:12.4e}  "
              f"{err[k,1]:12.4e}  {err[k,2]:12.4e}  "
              f"{divs[k,0]:12.4e}  {divs[k,1]:12.4e}")

    hvals = np.array([L / n for n in sizes])

    # ---- Mesh convergence (final time) ----
    fig, ax = plt.subplots(figsize=(7, 5))
    for col, label, marker in [(0, "u max error", "o"), (1, "v max error", "s")]:
        ax.loglog(hvals, err[:, col], marker + "-", label=label)
    if len(sizes) > 1:
        ref = err[0, 0] * (hvals / hvals[0])**2
        ax.loglog(hvals, ref, "--", label=r"$O(h^2)$ reference")
    ax.set_xlabel("h = L/N")
    ax.set_ylabel("Maximum absolute error")
    ax.set_title(f"Velocity mesh convergence at t={args.time:g}")
    ax.grid(True, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "mesh_convergence.png", dpi=160)
    plt.close(fig)

    # ---- Divergence check ----
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

    # ---- Feature 2: velocity absolute error vs time, u and v in the same plot ----
    fig, ax = plt.subplots(figsize=(10, 6.5))
    for n in sizes:
        H = history[n]
        # t=0 skipped (log axis; error is zero there).
        ax.semilogy(H["t"][1:], H["eu"][1:], "-", lw=3, label=f"u error, N = {n}")
        ax.semilogy(H["t"][1:], H["ev"][1:], "--", lw=3, label=f"v error, N = {n}")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Abs error")
    ax.set_title("Error vs time for different N", fontweight="bold")
    ax.set_xlim(0, args.time)
    ax.grid(True, which="both", alpha=0.4)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(outdir / "velocity_error_vs_time.png", dpi=160)
    plt.close(fig)

    # ---- Feature 3: relative-pressure absolute error vs time, one line per N ----
    fig, ax = plt.subplots(figsize=(7, 5))
    for n in sizes:
        H = history[n]
        # t=0 skipped: p is only the initial guess (0) there, not a solution.
        ax.semilogy(H["t"][1:], H["ep"][1:], "-", label=f"N={n}")
    ax.set_xlabel("t")
    ax.set_ylabel(r"max $|\Delta p - \Delta p_{exact}|$,  $\Delta p_i = p_i - p_1$")
    ax.set_title("Relative pressure error vs time")
    ax.grid(True, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "pressure_error_vs_time.png", dpi=160)
    plt.close(fig)

    # ---- Summaries ----
    summary = outdir / "validation_summary.csv"
    np.savetxt(summary, np.column_stack((sizes, hvals, err, divs)), delimiter=",")

    for n in sizes:
        H = history[n]
        np.savetxt(
            outdir / f"time_history_N{n}.csv",
            np.column_stack((H["t"], H["un"], H["ue"], H["vn"], H["ve"],
                             H["eu"], H["ev"], H["ep"])),
            delimiter=",",
            header="t,u_num,u_exact,v_num,v_exact,max_err_u,max_err_v,max_err_dp",
            comments="",
        )

    print(f"\nPlots and summary saved in: {outdir}")
    print("Files: fields_N*.png, node_timeseries_N*.png, velocity_error_vs_time.png, "
          "pressure_error_vs_time.png, mesh_convergence.png, divergence_check.png, "
          "validation_summary.csv, time_history_N*.csv")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}")
