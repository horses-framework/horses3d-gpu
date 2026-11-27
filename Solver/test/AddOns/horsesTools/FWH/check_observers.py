#!/usr/bin/env python3
"""
Checks the FWH observers of the horses.tools CI test.

Usage: ./check_observers.py ON_THE_FLY_PREFIX POST_PREFIX [--rtol-post R] [--rtol-ref R]

   e.g. ./check_observers.py RESULTS/CylinderFWH RESULTS/FWHPost

For every observer in REFERENCE, reads <prefix>.<observer>.observer for both runs and checks:
  1. On the fly vs post-processing: same iterations and source times, and the same acoustic
     pressure (P'T, P'L, P'). The post-processing reads Q and Qdot from the saved surface
     files, so it must reproduce the on the fly integrals up to round-off.
  2. Reference values: P'L and P' of the on the fly run at some iterations. The surface files
     are written from the same face storage used on the fly, so a wrong prolongation of Q or Qdot
     to the faces changes both runs alike: only this check catches it.
  3. The signal is not trivially zero: P' must be non-zero after the retarded time start up.

Tolerances are relative to the maximum |P'| of each observer, to be robust at zero crossings.
"""
import argparse
import math
import sys

# Observer name -> {iteration: (P'L, P')}, from CylinderFWH.control (solid surface, so P'T = 0)
REFERENCE = {
    "side": {
        25: (-7.3683079330274410e+01, -7.3683079330274410e+01),
        40: (1.4931395382570478e+00, 1.4931395382570478e+00),
        50: (2.1027754917901174e+00, 2.1027754917901174e+00),
    },
    "upstream": {
        25: (-3.6312103110338882e+02, -3.6312103110338882e+02),
        40: (4.3878832052623249e+00, 4.3878832052623249e+00),
        50: (1.3724032359002623e+01, 1.3724032359002623e+01),
    },
}
# First iteration with a non-zero signal (before it, the observer time cannot be interpolated)
FIRST_NONZERO = 20


def read_observer(fileName):
    """Returns {iteration: (time, observer_time, P'T, P'L, P')}"""
    data = {}
    try:
        with open(fileName) as f:
            for line in f:
                tokens = line.split()
                if len(tokens) != 6:
                    continue
                try:
                    it = int(tokens[0])
                    values = tuple(float(v) for v in tokens[1:])
                except ValueError:
                    continue  # header
                data[it] = values
    except OSError as e:
        sys.exit(f"FAILED: cannot read {fileName}: {e}")
    if not data:
        sys.exit(f"FAILED: no data in {fileName}")
    for it, v in data.items():
        if not all(math.isfinite(x) for x in v):
            sys.exit(f"FAILED [{fileName}]: non finite value at iteration {it}")
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fly")
    parser.add_argument("post")
    parser.add_argument("--rtol-post", type=float, default=1.0e-10)
    parser.add_argument("--rtol-ref", type=float, default=1.0e-6)
    args = parser.parse_args()

    errors = []
    for obs, ref in REFERENCE.items():
        flyFile = f"{args.fly}.{obs}.observer"
        postFile = f"{args.post}.{obs}.observer"
        fly = read_observer(flyFile)
        post = read_observer(postFile)
        scale = max(abs(v[4]) for v in fly.values())
        print(f"Observer {obs}: {len(fly)} lines on the fly, {len(post)} post-processed, max |P'| = {scale:.6e}")
        if scale == 0.0:
            errors.append(f"{obs}: the acoustic pressure is zero everywhere")
            continue

        # 1) on the fly vs post-processing
        if sorted(fly) != sorted(post):
            errors.append(f"{obs}: different iterations on the fly and post-processed")
        else:
            maxDiff = 0.0
            for it in fly:
                f, p = fly[it], post[it]
                if abs(f[0] - p[0]) > 1.0e-12 * max(1.0, abs(f[0])):
                    errors.append(f"{obs}: different source time at iteration {it}: {f[0]} vs {p[0]}")
                    break
                maxDiff = max(maxDiff, max(abs(a - b) for a, b in zip(f[2:], p[2:])) / scale)
            print(f"   on the fly vs post: max relative difference = {maxDiff:.3e} (tol {args.rtol_post:.1e})")
            if maxDiff > args.rtol_post:
                errors.append(f"{obs}: on the fly and post-processed acoustic pressure differ ({maxDiff:.3e})")

        # 2) reference values
        for it, (pl, p) in ref.items():
            if it not in fly:
                errors.append(f"{obs}: iteration {it} not found")
                continue
            got = fly[it]
            for name, r, v in (("P'L", pl, got[3]), ("P'", p, got[4])):
                diff = abs(v - r) / scale
                status = "ok" if diff <= args.rtol_ref else "FAILED"
                print(f"   iter {it:3d} {name:3s} = {v: .10e}  ref = {r: .10e}  rel diff = {diff:.2e}  {status}")
                if diff > args.rtol_ref:
                    errors.append(f"{obs}: {name} at iteration {it} = {v}, reference {r}")

        # 3) non-zero signal after the start up
        late = [v[4] for it, v in fly.items() if it >= FIRST_NONZERO]
        if not late or any(x == 0.0 for x in late):
            errors.append(f"{obs}: zero acoustic pressure after iteration {FIRST_NONZERO}")

    if errors:
        print("\nFAILED:")
        for e in errors:
            print("  * " + e)
        sys.exit(1)
    print("\nFWH observers OK")


if __name__ == "__main__":
    main()
