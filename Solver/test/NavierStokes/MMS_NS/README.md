# MMS_NS — Method of Manufactured Solutions convergence suite (compressible Navier–Stokes)

This suite verifies the spatial (and temporal) convergence of the compressible
Navier–Stokes solver (`horses3d.ns`). It forces a chosen analytical
("manufactured") solution through a source term, sets that solution as the
initial condition, and measures the discrete `L2` error against it. For a smooth
manufactured solution the error decreases spectrally with the polynomial order
`P`, which is the signature of a correct high-order discretisation.

## Background: how the method works

The Method of Manufactured Solutions is a way to obtain an *exact* error for a
discretisation, even when the PDE has no analytical solution:

1. **Pick** a smooth field `Q_ex(x,t)` — it need not be physical, only smooth and
   compatible with the boundary conditions.
2. **Substitute** it into the governing equations. It will not satisfy them, so it
   leaves a residual `S = ∂Q_ex/∂t + ∇·F(Q_ex)` (inviscid + viscous fluxes).
3. **Add** `S` as a source term. Now `Q_ex` is, by construction, an exact solution
   of the *forced* equations.
4. **Solve** the forced equations numerically starting from `Q_ex` at `t=0`, and
   compare the numerical result to `Q_ex` at the final time. The difference is
   pure discretisation error.

Because the manufactured solution here is analytic (products of sines/cosines),
the error is expected to fall **exponentially (spectrally) with `P`** at a fixed
mesh — roughly an order of magnitude per order — and algebraically,
`O(h^{P+1})`, under mesh refinement at fixed `P`. With the short default final
time the spatial error dominates (temporal Runge–Kutta error is negligible), so
the `P`-sweep on a fixed mesh is the cleanest convergence signal.

## ⚠️ Works only WITHOUT OpenACC (gfortran / CPU, serial)

The manufactured forcing is applied through `UserDefinedSourceTermNS`, which
lives in the problemfile shared library and runs **on the host**. On the GPU
build (`nvfortran`/OpenACC) that source term is compiled out — see the
`#ifndef _OPENACC` guard in
`Solver/src/NavierStokesSolver/SpatialDiscretization.f90`
(`TimeDerivative_ComputeQDot`). There the forcing is **not** applied, so the
manufactured solution is not a solution of the discrete system and the test is
meaningless (the error will look flat in `P` and `dt`).

**Build and run this suite with `gfortran`, `ENABLE_THREADS=NO`, `COMM=SEQUENTIAL`.**
Do the convergence analysis in serial; if you later want to check the GPU, run a
few *relevant* (non-MMS) cases and confirm they match the serial results.

## The default manufactured solution

Defined in `SETUP/generate_problemfile_NS.py`, periodic on `[-1,1]^3`:

```
ρ = 1                                             (constant)
u =  sin(πx) cos(πy) cos(πz) sin(t)
v =  cos(πx) sin(πy) cos(πz) sin(t)
w = -2 cos(πx) cos(πy) sin(πz) sin(t)             (⇒ ∇·u = 0)
p = 1 + 0.1 sin(πx) sin(πy) sin(πz) cos(t)        (stays positive)
```

The velocity is divergence-free and every field is periodic, so it is compatible
with the periodic meshes. You can replace any of these — see *Editing* below.

## Contents

| File | Role |
|------|------|
| `SETUP/generate_problemfile_NS.py` | Generates `SETUP/ProblemFile.f90` (initial condition, source terms, `L2` error) from a symbolic manufactured solution. **This is the source of truth — edit this, not the generated `.f90`.** |
| `control_template.control` | Control-file template; the runner fills the `{...}` tokens per case. |
| `run_convergence_NS.py` | Runs the `(mesh, P)` sweep, collects results into `errors.csv`. |
| `plot_convergence_NS.py` | Plots the convergence from `errors.csv`. |
| `../../TestMeshes/MMS_cube{N}.h5` | Shared cube meshes: `N` elements per direction, periodic, on `[-1,1]^3` (also used by the Multiphase suite `test/Multiphase/MMS_MU`). |

Generated / transient (not tracked): `SETUP/ProblemFile.f90`, `SETUP/*.o/.mod/.so`,
`CONTROL/`, `RESULTS/`, `errors.csv`, `mms_l2_error.dat`, `convergence.*`.

## Prerequisites

- `gfortran` and GNU `make`.
- HDF5 with Fortran bindings (to read the `.h5` meshes), built with a compatible
  gfortran. Point `HDF5_DIR` at an install whose layout is
  `$HDF5_DIR/lib/libhdf5_fortran.*` and `$HDF5_DIR/include/*.mod`.
- Python 3 with `sympy` and `numpy` (generator), and `matplotlib` (plotting).

## Build & run

All paths below are relative to the `Solver/` directory.

1. **Build the solver** (once), with HDF5 and no threads/OpenACC:
   ```bash
   export HDF5_DIR=/path/to/hdf5
   make ns COMPILER=gfortran ENABLE_THREADS=NO WITH_HDF5=YES COMM=SEQUENTIAL HDF5_DIR=$HDF5_DIR
   ```
   If you are switching from a previous build (e.g. a different compiler or HDF5),
   run `make clean` first to avoid stale-module errors (see *Troubleshooting*).

2. **Generate and compile the problemfile** (in the test's `SETUP/`):
   ```bash
   cd test/NavierStokes/MMS_NS/SETUP
   python3 generate_problemfile_NS.py            # writes ProblemFile.f90
   make ns COMPILER=gfortran ENABLE_THREADS=NO MODE=RELEASE   # builds libproblemfile_ns.so
   cd ..
   ```
   The solver picks up `./SETUP/libproblemfile_ns.so` at runtime via its RUNPATH,
   so it must be run from the `MMS_NS/` directory (the runner does this). You can
   confirm the right library is loaded with `ldd ./horses3d.ns | grep problemfile`.

3. **Run the study** (in `MMS_NS/`). Make the HDF5 runtime library visible so the
   solver can read the meshes:
   ```bash
   export LD_LIBRARY_PATH=$HDF5_DIR/lib:$LD_LIBRARY_PATH
   python3 run_convergence_NS.py                 # add --dry-run to preview cases
   ```
   Results are appended to `errors.csv` (`nelems, P, NDOF, L2_error, t_final`).
   The run is resumable: cases already in `errors.csv` are skipped (delete the row,
   or the file, to re-run).

4. **Plot** the convergence:
   ```bash
   python3 plot_convergence_NS.py
   ```

## Interpreting the results

`errors.csv` has one row per case: `nelems, P, NDOF, L2_error, t_final`. The
`L2_error` is the true `L2` norm over the domain, `sqrt(∫ ‖Q_ex − Q_h‖² dV)`.

A correct run shows, at each fixed mesh, an error that drops by roughly an order
of magnitude per `+1` in `P` until it hits the floating-point / time-integration
floor (~`1e-12`–`1e-13`). For example (one mesh, `P = 1…5`):

```
P=1  ~5e-8
P=2  ~1e-8
P=3  ~2e-9      each step ≈ ×0.1  → spectral convergence
P=4  ~3e-10
P=5  ~3e-11
```

If instead the error is **flat in both `P` and `dt`**, the forcing is not
entering — you almost certainly built with OpenACC/nvfortran (where the source is
guarded out) or the default problemfile is being loaded instead of `./SETUP`'s.

## Editing / customising

Three files must stay mutually consistent. The parameters that must agree are the
ratio of specific heats `GAMMA`, the quadrature node type, and the maximum order.

- **Manufactured solution and physics** — `SETUP/generate_problemfile_NS.py`:
  - `SECTION 1`: `RHO0`, `GAMMA` (must match the control file).
  - `SECTION 2`: `rho, u, v, w, p` as SymPy expressions. They must be **periodic on
    `[-1,1]^3`** (the meshes use periodic BCs) and keep pressure positive. A
    divergence-free velocity is convenient but not required.
  - `SECTION 3`: `NODE_TYPE` must match `Discretization nodes` in the control file;
    `MAX_ORDER` is the highest `P` for which quadrature weights are tabulated (raise
    it if you want higher orders); `SIMPLIFY_STRATEGY` (`"none"` is fastest and is
    fine — simplification is only cosmetic, not required for correctness).
  - The source term is derived to match what the solver actually evaluates,
    including the **temperature-dependent (Sutherland) viscosity** and thermal
    conductivity, so it stays consistent with `horses3d`. If you change the
    Reynolds/Mach/Prandtl numbers or the Sutherland settings in the control file,
    no generator change is needed — those are read at runtime — but a change to the
    *fields* or `GAMMA`/node type does require regenerating.
  - **Re-run the generator and rebuild the problemfile** after any change here.

- **Study parameters** — `run_convergence_NS.py`, `SECTION 2`: `MESH_SIZES`,
  `P_VALUES` (≤ `MAX_ORDER`), `DT`, `T_FINAL`. Keep `DT` small enough that the
  time-integration error stays below the spatial error you want to observe.

- **Solver / numerical settings** — `control_template.control`: Mach, Reynolds,
  Prandtl, `Discretization nodes`, Riemann/viscous scheme, time integration. These
  must stay consistent with the generator (`GAMMA`, node type) and the runner. Note
  the GPU build does not support all keywords (e.g. `viscous discretization = BR2`
  / `IP`, `Gradient Variables = Energy`); the template uses BR1, which is fine.

- **Meshes** — the runner expects `test/TestMeshes/MMS_cube{N}.h5` with `N` elements per
  direction on `[-1,1]^3`, all-periodic (`MESH_DIR = ../../TestMeshes`). Add more by
  generating additional `MMS_cube{N}.h5`
  (e.g. with HOPR) and listing `N` in `MESH_SIZES`.

## Troubleshooting

- **`Fatal Error: Mismatch in components of derived type 'c_h5o_info_t'`** (or other
  `.mod` mismatch): a stale build mixing objects/modules. Run `make clean` (or
  `make allclean`) and rebuild.
- **`libhdf5_fortran.so.200 => not found`** at run time: the HDF5 runtime library
  is not on the loader path. `export LD_LIBRARY_PATH=$HDF5_DIR/lib:$LD_LIBRARY_PATH`.
- **`MMS: order N not tabulated`** (`error stop`): `P` exceeds `MAX_ORDER` in the
  generator. Raise `MAX_ORDER`, regenerate, rebuild the problemfile.
- **Error is flat in `P` and `dt` / no convergence**: the source term is not being
  applied. Check that you built with gfortran (not nvfortran/OpenACC) and that the
  local `./SETUP/libproblemfile_ns.so` is loaded (`ldd ./horses3d.ns`).
- **`NaN` / negative pressure**: `DT` (or `T_FINAL`) too large for the finest mesh /
  highest `P`, or a manufactured pressure that goes non-positive. Reduce `DT` or
  adjust the pressure field.
- **Wrong / non-converging error with correct-looking runs**: `NODE_TYPE` in the
  generator disagrees with `Discretization nodes` in the control file, so the error
  quadrature uses the wrong weights.

## How it works (details)

`generate_problemfile_NS.py` forms the strong-form NS residual
`S = ∂Q/∂t + ∇·F_inv − ∇·F_visc` for the chosen manufactured fields, using the
**same nondimensional conventions as `horses3d`**:

- viscous stress `τ = μ(∇u + ∇uᵀ − ⅔(∇·u)I)` with `μ = suther(T)/Re` (Sutherland's
  law, `T = γM² p/ρ`), matching `get_laminar_mu_kappa`/`ViscousFlux_STATE`;
- heat flux `q = −κ ∇T` with `κ = μ γ / ((γ−1) Pr)`, which reproduces the solver's
  `κ · ∇T = dimensionless%kappa · γM² ∇(p/ρ)`;
- total energy `ρe = p/(γ−1) + ½ρ‖u‖²`.

Because `μ` and `κ` vary in space, the symbolic divergence automatically includes
the `∇μ`, `∇κ` contributions the DG operator also produces. The sign convention
matches the solver: `S` is added to `QDot` (`QDot += S_NS`), and `S = ∂Q/∂t + ∇·F`.

The generator emits three routines into `ProblemFile.f90`:

- `UserDefinedInitialCondition` — sets `Q` to the exact solution at `t = 0`;
- `UserDefinedSourceTermNS` — the forcing `S` (with viscosity/conductivity read at
  runtime from the solver's `dimensionless`/`refValues`);
- `UserDefinedFinalize` — the Gauss-quadrature `L2` error, written to
  `mms_l2_error.dat` as `nelems, P, NDOF, L2_error, t_final`.

The runner fills the control template for each `(mesh, P)`, runs the solver, reads
`mms_l2_error.dat`, and appends the row to `errors.csv`.
