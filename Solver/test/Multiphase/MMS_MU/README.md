# MMS_MU — Method of Manufactured Solutions convergence suite (Multiphase, Navier–Stokes–Cahn–Hilliard)

This suite verifies the spatial (and temporal) convergence of the multiphase
solver (`horses3d.mu`), which solves the incompressible Navier–Stokes–Cahn–Hilliard
system (a diffuse-interface, two-fluid model). It forces a chosen analytical
("manufactured") solution through a source term, sets it as the initial
condition, and measures the discrete `L2` error against it. For a smooth
manufactured solution the error decreases spectrally with the polynomial order
`P`, the signature of a correct high-order discretisation.

This is the multiphase analogue of `test/NavierStokes/MMS_NS`; see that README
for the shared background. The two suites share the same cube meshes.

## Background: how the method works

1. **Pick** a smooth field `Q_ex(x,t)` compatible with the boundary conditions.
2. **Substitute** it into the governing equations — it leaves a residual
   `S = ∂Q_ex/∂t + ∇·F(Q_ex)` (advective + viscous/Cahn–Hilliard fluxes).
3. **Add** `S` as a source term, so `Q_ex` is an exact solution of the forced system.
4. **Solve** numerically from `Q_ex` at `t=0` and compare to `Q_ex`; the
   difference is pure discretisation error.

Because the manufactured fields are analytic, the error falls **exponentially
(spectrally) with `P`** at a fixed mesh, and algebraically, `O(h^{P+1})`, under
mesh refinement at fixed `P`.

## ⚠️ Works only WITHOUT OpenACC (gfortran / CPU, serial)

The manufactured forcing is applied through `UserDefinedSourceTermNS`, a host
routine in the problemfile library. On the GPU build (`nvfortran`/OpenACC) it is
compiled out — see the `#ifndef _OPENACC` guard in
`Solver/src/MultiphaseSolver/SpatialDiscretization.f90`
(`TimeDerivative_ComputeQDot`). There the forcing is **not** applied and the test
is meaningless (error flat in `P` and `dt`).

**Build and run this suite with `gfortran`, `ENABLE_THREADS=NO`, `COMM=SEQUENTIAL`.**

## The multiphase equations and variables

The `NCONS = 5` conserved variables are
`Q = [c, √ρ·u, √ρ·v, √ρ·w, p]` (`IMC, IMSQRHOU, IMSQRHOV, IMSQRHOW, IMP`), with a
`√ρ` transformation on the momentum. The mixture density is `ρ = c·ρ₁ + (1−c)·ρ₂`.
The system is: a Cahn–Hilliard equation for the concentration `c`, three momentum
equations, and an artificial-compressibility pressure equation.

## The default manufactured solution

Defined in `SETUP/generate_problemfile_MU.py`, periodic on `[-1,1]^3`:

```
c = 1/2 (1 + cos(πx) cos(πy) cos(πz) sin(t))    (concentration, in [0,1])
u =  sin(πx) cos(πy) cos(πz) sin(t)
v =  cos(πx) sin(πy) cos(πz) sin(t)
w = -2 cos(πx) cos(πy) sin(πz) sin(t)            (⇒ ∇·u = 0)
p =  sin(πx) sin(πy) sin(πz) cos(t)
```

All fields are periodic; the velocity is divergence-free.

## Contents

| File | Role |
|------|------|
| `SETUP/generate_problemfile_MU.py` | Generates `SETUP/ProblemFile.f90` (initial condition, source terms, `L2` error) from the symbolic manufactured solution. **Edit this, not the generated `.f90`.** |
| `control_template.control` | Control-file template; the runner fills the `{...}` tokens per case. |
| `run_convergence_MU.py` | Runs the `(mesh, P)` sweep, collects results into `errors.csv`. |
| `plot_convergence_MU.py` | Plots the convergence from `errors.csv`. |
| `../../TestMeshes/MMS_cube{N}.h5` | Shared cube meshes: `N` elements per direction, periodic, on `[-1,1]^3` (also used by `test/NavierStokes/MMS_NS`). |

Generated / transient (not tracked): `SETUP/ProblemFile.f90`, `SETUP/*.o/.mod/.so`,
`CONTROL/`, `RESULTS/`, `errors.csv`, `mms_l2_error.dat`, `convergence.*`.

## Prerequisites

- `gfortran` and GNU `make`.
- HDF5 with Fortran bindings (to read the `.h5` meshes). Use the system serial
  HDF5 (HDF5 1.10, matches horses3d's API):
  `export HDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/serial`.
- Python 3 with `sympy` and `numpy` (generator), and `matplotlib` (plotting).

## Build & run

All paths below are relative to the `Solver/` directory.

1. **Build the solver** (once), with HDF5 and no threads/OpenACC:
   ```bash
   export HDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/serial
   make mu COMPILER=gfortran ENABLE_THREADS=NO WITH_HDF5=YES COMM=SEQUENTIAL
   ```
   Run `make clean` first when switching compiler/HDF5 (stale-module errors otherwise).

2. **Generate and compile the problemfile** (in the test's `SETUP/`):
   ```bash
   cd test/Multiphase/MMS_MU/SETUP
   python3 generate_problemfile_MU.py            # writes ProblemFile.f90
   make mu COMPILER=gfortran ENABLE_THREADS=NO MODE=RELEASE   # builds libproblemfile_mu.so
   cd ..
   ```
   The solver picks up `./SETUP/libproblemfile_mu.so` at runtime via its RUNPATH, so
   run from the `MMS_MU/` directory (the runner does this). Confirm with
   `ldd ./horses3d.mu | grep problemfile`.

3. **Run the study** (in `MMS_MU/`):
   ```bash
   export LD_LIBRARY_PATH=$HDF5_ROOT/lib:$LD_LIBRARY_PATH
   python3 run_convergence_MU.py                 # add --dry-run to preview cases
   ```
   Results are appended to `errors.csv` (`nelems, P, NDOF, L2_error, t_final`).
   Cases already in `errors.csv` are skipped (delete the row, or the file, to re-run).

4. **Plot** the convergence:
   ```bash
   python3 plot_convergence_MU.py                # writes convergence.pdf / .png
   ```

## Interpreting the results

`errors.csv` has one row per case: `nelems, P, NDOF, L2_error, t_final`. The
`L2_error` is the true `L2` norm over the domain, `sqrt(∫ ‖Q_ex − Q_h‖² dV)`.

A correct run shows, at each fixed mesh, an error that drops by roughly an order
of magnitude per `+1` in `P` until it hits the floating-point / time-integration
floor. For example (one mesh, `P = 2…6`):

```
P=2  ~6e-7
P=3  ~8e-8
P=4  ~8e-9      each step ≈ ×0.1  → spectral convergence
P=5  ~6e-10
P=6  ~4e-11
```

If instead the error is **flat in both `P` and `dt`**, the forcing is not
entering — you built with OpenACC/nvfortran (where the source is guarded out) or
the default problemfile is loaded instead of `./SETUP`'s.

## Editing / customising

Three files must stay mutually consistent; the parameters that must agree are the
fluid densities `ρ₁, ρ₂`, the quadrature node type, and the maximum order.

- **Manufactured solution and physics** — `SETUP/generate_problemfile_MU.py`:
  - `SECTION 1`: `RHO1`, `RHO2` (documentation only; the source reads
    `dimensionless_ % rho(1,2)` at runtime — they must match the control file).
  - `SECTION 2`: `c, u, v, w, p` as SymPy expressions, **periodic on `[-1,1]^3`**;
    keep `c ∈ [0,1]` (the solver clips density/viscosity/sound-speed blends to `[0,1]`).
  - `SECTION 3`: `NODE_TYPE` must match `Discretization nodes`; `MAX_ORDER` is the
    highest `P` for which quadrature weights are tabulated; `SIMPLIFY_STRATEGY`
    (`"none"` is fastest and fine — simplification is only cosmetic).
  - The source matches what the solver evaluates, including the **variable
    (c-dependent) viscosity** `μ(c)=c·μ₁+(1−c)μ₂` and the **variable speed of sound**
    `invMa2=(cs₁c+cs₂(1−c))²ρ` (all read at runtime), so different per-fluid
    viscosities / sound speeds won't break it.
  - **Re-run the generator and rebuild the problemfile** after any change here.

- **Study parameters** — `run_convergence_MU.py`, `SECTION 2`: `MESH_SIZES`,
  `P_VALUES` (≤ `MAX_ORDER`), `DT`, `T_FINAL`.

- **Solver / numerical settings** — `control_template.control`: fluid densities and
  viscosities, per-fluid sound speeds, interface width, surface tension, mobility,
  `Discretization nodes`, `Explicit method` (e.g. `RK3`). The per-fluid
  `fluid 1/2 sound speed square (m/s)` keys enable the variable speed-of-sound
  branch the source term assumes — keep them set.

## Troubleshooting

- **`Fatal Error: Mismatch in components of derived type 'c_h5o_info_t'`** (or other
  `.mod` mismatch): stale build — `make clean` and rebuild.
- **`Cannot open module file 'hdf5.mod'`** when building `mu`: the MU Makefile
  historically dropped the HDF5 include (used an undefined `EXTLIB_MU`); the
  main-program compile rule in `src/MultiphaseSolver/Makefile` must use
  `$(EXTLIB_INC)` like the NS solver.
- **`libhdf5_serial_fortran.so … not found`** at run time:
  `export LD_LIBRARY_PATH=$HDF5_ROOT/lib:$LD_LIBRARY_PATH`.
- **`MMS: order N not tabulated`**: `P` exceeds `MAX_ORDER` — raise it, regenerate.
- **Error flat in `P` and `dt`**: source not applied — build with gfortran (not
  nvfortran/OpenACC) and check `./SETUP/libproblemfile_mu.so` is loaded (`ldd`).
- **`NaN` / clipped fields**: `DT`/`T_FINAL` too large, or the manufactured `c` left
  `[0,1]` (density/viscosity/sound-speed blends are clipped there).

## How it works (details)

`generate_problemfile_MU.py` forms the strong-form residual
`S = ∂Q/∂t + ∇·F(Q_ex)` for the manufactured fields, matching the solver's
formulation:

- **concentration (Cahn–Hilliard):** `∂c/∂t + ∇·(c u) − M₀ ∇²μ`, with chemical
  potential `μ = 48(σ/ε) c(c−1)(c−½) − 1.5 σε ∇²c` (explicit path, no `M₀` factor).
- **momentum** (`√ρ·u` variable): `√ρ ∂(√ρ u)/∂t + ½ρ (u·∇)u + c ∇μ + ∇p − ∇·τ`,
  with variable viscosity `τ = μ(c)(∇u+∇uᵀ)`, `μ(c)=c·μ₁+(1−c)μ₂`.
- **pressure** (artificial compressibility): `∂p/∂t + invMa2 ∇·u`,
  `invMa2 = (cs₁c + cs₂(1−c))² ρ`, `cs_i = √(invMa2(i)/ρ(i))`.

Because `μ`, `ρ` and `cs` vary in space, the symbolic divergence automatically
includes the gradient contributions the DG operator also produces. The source is
added to `QDot` with the `√ρ` momentum scaling: `QDot += S / [1, √ρ, √ρ, √ρ, 1]`
(`Q(IMSQRHOU)=√ρ u`, so `QDot = ∂(√ρ u)/∂t`). The generator emits
`UserDefinedInitialCondition` (sets `Q` to the exact solution),
`UserDefinedSourceTermNS` (the forcing, reading `dimensionless_/multiphase_` at
runtime), and `UserDefinedFinalize` (the `L2` error → `mms_l2_error.dat`).
