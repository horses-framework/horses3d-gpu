#!/usr/bin/env python3
"""
Validate CylinderProbeFileASCII output.

Checks:
  1. Classic inline probe and all 4 file-probe ASCII outputs exist.
  2. All 4 file probes recorded the same number of time steps.
  3. Physical sanity: the inline-probe 'u' values are finite and non-trivial.
  4. Agreement: every file-probe 'u' column matches the inline-probe 'u'
     column to within a tight tolerance (same point, same variable).
  5. Physical sanity on the remaining file-probe variables (rho > 0, p > 0).

Usage:
  python3 validate.py [RESULTS_DIR]
  RESULTS_DIR defaults to ./RESULTS
"""

import sys
import os
import numpy as np


def parse_probe_file(path):
    """Return (header_vars, data) where data is a 2-D float array (nsteps, nvars).

    The .probe file written by HORSES3D looks like:
      Monitor name:      <name>
      x, y, z coordinates: ...
      <blank>
      Iteration   Time   var1   var2 ...
      <data rows>
    """
    with open(path) as fh:
        lines = fh.readlines()

    # Find the header line (contains "Iteration")
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Iteration"):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(f"No 'Iteration' header found in {path}")

    header_vars = lines[header_idx].split()[2:]  # skip "Iteration" and "Time"
    data_lines = [l for l in lines[header_idx + 1:] if l.strip()]
    if not data_lines:
        raise ValueError(f"No data rows found in {path}")

    data = np.array([[float(x) for x in l.split()] for l in data_lines])
    # columns: iter, time, var1, var2 ...
    return header_vars, data


def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else "RESULTS"
    probes_dir = os.path.join(results_dir, "probes")
    # solution file name is CylinderProbeFileASCII.hsol, so probe files
    # carry that stem: CylinderProbeFileASCII.hsol.<name>.probe
    stem = "CylinderProbeFileASCII.hsol"

    # File paths
    inline_path = os.path.join(probes_dir, f"{stem}.wake_u.probe")
    # 1 classic probe + 4 file probes → file probes are probe_2 … probe_5
    fp_paths = [os.path.join(probes_dir, f"{stem}.probe_{i}.probe") for i in range(2, 6)]

    errors = []

    # --- 1. Existence ---
    for p in [inline_path] + fp_paths:
        if not os.path.isfile(p):
            errors.append(f"Missing file: {p}")

    if errors:
        for e in errors:
            print(f"ERROR: {e}", flush=True)
        sys.exit(1)

    # --- Parse ---
    inline_vars, inline_data = parse_probe_file(inline_path)
    fp_data = []
    for p in fp_paths:
        _, d = parse_probe_file(p)
        fp_data.append(d)

    n_steps = inline_data.shape[0]

    # --- 2. Step count consistency ---
    for i, d in enumerate(fp_data):
        if d.shape[0] != n_steps:
            errors.append(
                f"probe_{i+2} has {d.shape[0]} steps but inline probe has {n_steps}"
            )

    # --- 3. Physical sanity on inline probe (u column, index 2) ---
    u_inline = inline_data[:, 2]
    if not np.all(np.isfinite(u_inline)):
        errors.append("Inline probe 'u' contains non-finite values")
    if np.all(u_inline == 0.0):
        errors.append("Inline probe 'u' is identically zero (probe not found?)")

    # --- 4. Agreement between file probes and inline probe ---
    tol = 1.0e-10
    for i, d in enumerate(fp_data):
        u_fp = d[:, 2]  # column 2 is 'u' (same variables= u v w rho pressure mach k velocity)
        max_diff = np.max(np.abs(u_fp - u_inline))
        if max_diff > tol:
            errors.append(
                f"probe_{i+2} 'u' differs from inline probe by {max_diff:.3e} (tol {tol:.0e})"
            )

    # --- 5. Physical sanity on other variables (rho > 0, pressure > 0) ---
    # File-probe variable order: u v w rho pressure mach k velocity
    # columns in data: 0=iter, 1=time, 2=u, 3=v, 4=w, 5=rho, 6=pressure, 7=mach, 8=k, 9=velocity
    rho_col = 5
    p_col = 6
    for i, d in enumerate(fp_data):
        if d.shape[1] <= p_col:
            errors.append(f"probe_{i+2} has fewer columns than expected")
            continue
        rho = d[:, rho_col]
        pres = d[:, p_col]
        if not np.all(rho > 0.0):
            errors.append(f"probe_{i+2}: rho <= 0 detected (min={rho.min():.4e})")
        if not np.all(pres > 0.0):
            errors.append(f"probe_{i+2}: pressure <= 0 detected (min={pres.min():.4e})")

    if errors:
        print("VALIDATION FAILED:", flush=True)
        for e in errors:
            print(f"  ERROR: {e}", flush=True)
        sys.exit(1)

    print(
        f"OK: {n_steps} steps, 4 file probes match inline probe "
        f"(max |du| < {tol:.0e}), rho>0, p>0",
        flush=True,
    )


if __name__ == "__main__":
    main()
