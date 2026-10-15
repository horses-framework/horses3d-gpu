#!/usr/bin/env python3
"""
Validate the probe resume of CylinderProbeFileHDF5.

CylinderProbeFileHDF5.control runs iterations 0..100 and writes the probe file;
CylinderProbeFileHDF5_resume.control restarts from iteration 50 with
"resume = .true." and runs to iteration 100 again. A copy of the probe file
taken before the restart is the reference. After the restart the file must:

  1. keep the same shape (no steps lost, none duplicated),
  2. have strictly increasing time and iteration,
  3. contain the same data as the reference in every dataset (the restart
     recomputes iterations 50..100 and must overwrite them, not append).

Gradient probes at the restart iteration are reported but not compared: in the
uninterrupted run the stored gradients come from the last Runge-Kutta stage,
whereas after a restart the first record uses gradients of the restored state.

Usage:
  python3 validate_resume.py REFERENCE.probes.h5 RESUMED.probes.h5 [--restart-iteration N]
"""

import os
import re
import subprocess
import sys
import tempfile

import numpy as np


def h5ls_shapes(h5_file):
    result = subprocess.run(["h5ls", h5_file], capture_output=True, text=True, check=True)
    shapes = {}
    for line in result.stdout.splitlines():
        b0, b1 = line.find("{"), line.find("}")
        if b0 == -1 or len(line.split()) < 3:
            continue
        try:
            shapes[line.split()[0]] = tuple(
                int(d.strip().split("/")[0]) for d in line[b0 + 1:b1].split(","))
        except ValueError:
            continue
    return shapes


def h5dump_array(h5_file, dset, shape, dtype):
    with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as tmp:
        path = tmp.name
    try:
        subprocess.run(
            ["h5dump", "-d", dset, "-s", ",".join("0" for _ in shape),
             "-c", ",".join(str(n) for n in shape), "-b", "LE", "-o", path, h5_file],
            capture_output=True, check=True)
        arr = np.fromfile(path, dtype=dtype)
    finally:
        os.unlink(path)
    return arr.reshape(shape)


def main():
    args = sys.argv[1:]
    restart_iter = None
    if "--restart-iteration" in args:
        i = args.index("--restart-iteration")
        restart_iter = int(args[i + 1])
        del args[i:i + 2]
    if len(args) != 2:
        print(__doc__)
        sys.exit(2)
    ref_path, new_path = args
    errors = []

    for p in (ref_path, new_path):
        if not os.path.isfile(p):
            errors.append(f"Missing file: {p}")
    if errors:
        for e in errors:
            print(f"ERROR: {e}", flush=True)
        sys.exit(1)

    ref_shapes, new_shapes = h5ls_shapes(ref_path), h5ls_shapes(new_path)
    if ref_shapes != new_shapes:
        errors.append(f"Dataset shapes differ after resume: reference {ref_shapes} vs resumed {new_shapes}")
        for e in errors:
            print(f"ERROR: {e}", flush=True)
        sys.exit(1)

    n_steps = ref_shapes["time"][0]
    time_new = h5dump_array(new_path, "time", (n_steps,), np.float64)
    iter_new = h5dump_array(new_path, "iteration", (n_steps,), np.int32)

    if not np.all(np.diff(time_new) > 0.0):
        errors.append("time is not strictly increasing (duplicated or unordered steps)")
    if not np.all(np.diff(iter_new) > 0):
        errors.append("iteration is not strictly increasing (duplicated or unordered steps)")
    if iter_new[0] != 0:
        errors.append(f"first saved iteration is {iter_new[0]}, expected 0 (steps before the restart were lost)")

    grad_re = re.compile(r"^(u|v|w|p|rho)_[xyz]$")
    worst = {}
    for name, shape in ref_shapes.items():
        dtype = np.int32 if name == "iteration" else np.float64
        a = h5dump_array(ref_path, name, shape, dtype)
        b = h5dump_array(new_path, name, shape, dtype)
        if not np.all(np.isfinite(b)):
            errors.append(f"'{name}' has non-finite values after resume")
            continue
        rel = np.abs(a - b) / (1.0 + np.abs(a))
        if rel.ndim == 2:
            err_step = rel.max(axis=0)          # per saved step
        elif rel.ndim == 1 and name in ("time", "iteration"):
            err_step = rel
        else:
            err_step = np.array([rel.max()])
        skip = np.zeros(err_step.shape, dtype=bool)
        if restart_iter is not None and grad_re.match(name) and err_step.size == n_steps:
            skip = (iter_new == restart_iter)
            if skip.any():
                print(f"  info: '{name}' at restart iteration {restart_iter}: "
                      f"relative diff {err_step[skip].max():.3e} (not compared)", flush=True)
        checked = np.where(skip, 0.0, err_step)
        k = int(np.argmax(checked))
        worst[name] = float(checked[k])
        if checked[k] > 1.0e-9:
            where = f" (worst at iteration {iter_new[k]})" if err_step.size == n_steps else ""
            errors.append(f"'{name}' differs from the reference by {checked[k]:.3e} "
                          f"(relative, tol 1e-9){where}")

    if errors:
        print("VALIDATION FAILED:", flush=True)
        for e in errors:
            print(f"  ERROR: {e}", flush=True)
        sys.exit(1)

    print(f"OK: resumed probe file has {n_steps} steps (iterations {iter_new[0]}..{iter_new[-1]}), "
          f"no duplicates; {len(worst)} datasets match the reference "
          f"(max relative diff {max(worst.values()):.2e})", flush=True)


if __name__ == "__main__":
    main()
