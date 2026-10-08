#!/usr/bin/env python3
"""
Validate CylinderProbeFileHDF5 output.

Checks:
  1. The .probes.h5 file and the classic inline-probe ASCII file exist.
  2. The HDF5 file contains the expected datasets (/coordinates, /time,
     /iteration, /u, /v, /w, /rho, /pressure, /mach, /k, /velocity).
  3. Physical sanity: rho > 0, pressure > 0, mach >= 0, velocity >= 0.
  4. Agreement: the 'u' values for all 4 file probes match the inline-probe
     'u' values to within a tight tolerance.

The script reads HDF5 data using h5dump (subprocess) + numpy so that it
does not require h5py and avoids MPI library conflicts with pvpython.

Usage:
  python3 validate.py [RESULTS_DIR]
  RESULTS_DIR defaults to ./RESULTS
"""

import sys
import os
import subprocess
import tempfile
import numpy as np


# ---------------------------------------------------------------------------
# Minimal HDF5 reader via h5dump + h5ls (no h5py)
# ---------------------------------------------------------------------------

def _h5ls_shapes(h5_file):
    """Return {dataset_name: shape_tuple} by parsing h5ls output."""
    result = subprocess.run(
        ["h5ls", h5_file],
        capture_output=True, text=True, check=True,
    )
    shapes = {}
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        name = parts[0]
        brace_start = line.find("{")
        brace_end = line.find("}")
        if brace_start == -1:
            continue
        dims_str = line[brace_start + 1:brace_end]
        try:
            dims = tuple(int(d.strip()) for d in dims_str.split(","))
        except ValueError:
            continue
        shapes[name] = dims
    return shapes


def _h5dump_array(h5_file, dset, start, count, dtype=np.float64):
    """Extract a hyperslab from an HDF5 dataset using h5dump."""
    with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as tmp:
        tmppath = tmp.name

    try:
        cmd = [
            "h5dump",
            "-d", dset,
            "-s", ",".join(str(s) for s in start),
            "-c", ",".join(str(c) for c in count),
            "-b", "LE",
            "-o", tmppath,
            h5_file,
        ]
        subprocess.run(cmd, capture_output=True, check=True)
        arr = np.fromfile(tmppath, dtype=dtype)
    finally:
        os.unlink(tmppath)

    return arr.reshape(count)


# ---------------------------------------------------------------------------
# Classic probe parser
# ---------------------------------------------------------------------------

def parse_probe_file(path):
    with open(path) as fh:
        lines = fh.readlines()

    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Iteration"):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(f"No 'Iteration' header found in {path}")

    data_lines = [l for l in lines[header_idx + 1:] if l.strip()]
    if not data_lines:
        raise ValueError(f"No data rows found in {path}")

    data = np.array([[float(x) for x in l.split()] for l in data_lines])
    return data


# ---------------------------------------------------------------------------
# Main validation
# ---------------------------------------------------------------------------

def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else "RESULTS"
    probes_dir = os.path.join(results_dir, "probes")
    # solution file name is CylinderProbeFileHDF5.hsol, so probe files
    # carry that stem: CylinderProbeFileHDF5.hsol.<name>.probe / .probes.h5
    stem = "CylinderProbeFileHDF5.hsol"

    h5_path = os.path.join(probes_dir, f"{stem}.probes.h5")
    inline_path = os.path.join(probes_dir, f"{stem}.wake_u.probe")

    errors = []

    # --- 1. Existence ---
    for p in [h5_path, inline_path]:
        if not os.path.isfile(p):
            errors.append(f"Missing file: {p}")

    if errors:
        for e in errors:
            print(f"ERROR: {e}", flush=True)
        sys.exit(1)

    # --- 2. Dataset presence ---
    expected_dsets = {"coordinates", "time", "iteration",
                      "u", "v", "w", "rho", "pressure", "mach", "k", "velocity"}
    shapes = _h5ls_shapes(h5_path)
    missing = expected_dsets - set(shapes.keys())
    if missing:
        errors.append(f"Missing HDF5 datasets: {sorted(missing)}")
        for e in errors:
            print(f"ERROR: {e}", flush=True)
        sys.exit(1)

    n_probes = shapes["u"][0]   # (nProbes, nSteps) — Fortran-transposed
    n_steps  = shapes["u"][1]

    if n_probes != 4:
        errors.append(f"Expected 4 file probes in HDF5, got {n_probes}")

    # --- Read HDF5 'u' for all 4 probes: shape (n_probes, n_steps) → transpose ---
    u_h5 = _h5dump_array(h5_path, "u", start=(0, 0), count=(n_probes, n_steps)).T
    # u_h5 is now (n_steps, n_probes)

    # --- Read inline probe ---
    inline_data = parse_probe_file(inline_path)
    u_inline = inline_data[:, 2]  # col 0=iter, 1=time, 2=u

    if len(u_inline) != n_steps:
        errors.append(
            f"Inline probe has {len(u_inline)} steps but HDF5 has {n_steps}"
        )
    else:
        tol = 1.0e-10
        for i in range(n_probes):
            max_diff = np.max(np.abs(u_h5[:, i] - u_inline))
            if max_diff > tol:
                errors.append(
                    f"HDF5 probe {i} 'u' differs from inline by {max_diff:.3e} "
                    f"(tol {tol:.0e})"
                )

    # --- 5. Physical sanity (rho > 0, pressure > 0) ---
    for var, check_positive in [("rho", True), ("pressure", True),
                                  ("mach", False), ("velocity", False)]:
        arr = _h5dump_array(h5_path, var, start=(0, 0), count=(n_probes, n_steps))
        if not np.all(np.isfinite(arr)):
            errors.append(f"HDF5 '{var}' contains non-finite values")
        if check_positive and not np.all(arr > 0.0):
            errors.append(f"HDF5 '{var}' <= 0 detected (min={arr.min():.4e})")

    if errors:
        print("VALIDATION FAILED:", flush=True)
        for e in errors:
            print(f"  ERROR: {e}", flush=True)
        sys.exit(1)

    print(
        f"OK: {n_probes} file probes, {n_steps} steps; "
        f"HDF5 'u' matches inline probe (max |du| < {tol:.0e}); "
        f"rho>0, p>0",
        flush=True,
    )


if __name__ == "__main__":
    main()
