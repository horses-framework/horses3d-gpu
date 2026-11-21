# MMS convergence tests

Method of manufactured solutions for the compressible Navier-Stokes (`NS`) and
multiphase (`MU`) solvers. A smooth periodic solution is imposed with a source term, and
the L2 error at the final time gives the order of convergence under h-refinement.

The user source term is a host routine, so the tests run on CPU builds (gfortran,
single-threaded); with OpenACC the solvers do not call it. The scripts are the same as
in horses3d, where the other equation sets have cases too.

## Run on GitHub

Actions > MMS convergence > Run workflow. Choose the case (or all), the polynomial
orders and the meshes (empty for the defaults of each case in `sweeps.json`) and,
optionally, changes to the control file, which apply to every selected case (e.g.
`Reynolds number = 100; Riemann solver = Roe` with case NS: MU has no Roe solver). Each
case and order runs in its own job, in parallel; a job has a 6-hour limit, which mesh 32
at high P may exceed. A final job per case merges them, writes the table of errors and
orders to the job summary and uploads `errors.csv`, `convergence.png` and the control
files and logs of the runs. The default sweep, P = 3 to 6 on meshes 4, 8 and 16 for both
cases, takes about 40 minutes, set by NS at P = 6 on one thread.

GitHub shows the Run workflow button only for workflows on the default branch. On
another branch, once the workflow has run in the repository, use the GitHub CLI:

```
gh workflow run mms.yml --ref <branch> -f case=NS -f P="2 3 4" -f meshes="4 8 16" -f set="Reynolds number = 100"
```

## Run locally

Build the solver (`./configure && make ns COMPILER=gfortran ENABLE_THREADS=NO` in `Solver`), then:

```
pip install sympy gmsh matplotlib
python3 run.py NS                                   # default sweep of NS (sweeps.json)
python3 run.py NS --P 3 4 --meshes 4 8              # quick run
python3 run.py NS --set "Reynolds number = 100"     # change the control file
python3 run.py MU --P 3 4 --meshes 4 8 --check      # compare with MU/reference.csv
python3 plot.py NS                                  # NS/convergence.png
```

The scripts need Linux, since the problem file is built as a shared library (`.so`).
The gmsh module also needs the system libraries that the workflow installs
(`apt install libglu1-mesa libgl1 libxcursor1 libxft2 libxinerama1 libfontconfig1`).

`run.py` generates `<case>/SETUP/ProblemFile.f90`, builds it with the HORSES3D build
settings (`COMPILER`, `ENABLE_THREADS`, ... from the environment and `Solver/make.inc`),
and runs every mesh and P with `<case>/control_template.control`, filling in the mesh,
the solution file and the order. Everything else, including `dt` and the number of
steps, comes from the control file; edit it or override keys with `--set`. Keep a fixed
time step: with `cfl` the solver ignores `dt`, so the runs stop at another time and are
reported as failed. Results go
to `<case>/errors.csv`, and pairs already run with the same solver, problem file and
control file are skipped. `run.py` exits with an error if a run fails and, with
`--check` and the committed control file, if an error differs from
`<case>/reference.csv` by more than 1e-3.

The meshes are periodic cubes [-1,1]^3 with N^3 elements, made by the cartesian-box
generator of [horses3d-meshes](https://github.com/horses-framework/horses3d-meshes)
(Gmsh Python API). `run.py` downloads the script at a fixed commit into `MESH/`, or uses
`$HORSES3D_MESHES_DIR/cartesian-box/` if that variable points to a clone, and runs it
with `$GMSH_PYTHON` if set. No mesh files are kept in this repository.

## Files

* `generate.py`: manufactured solutions and equations of each case. The solutions are at
  full amplitude from t = 0, with a variable density (concentration in MU) and a velocity
  that is neither divergence-free nor irrotational, so that the terms of the equations
  are active from the first step. The source term is the residual of the solution,
  emitted without symbolic simplification. The physical parameters are read from the
  solver at run time.
* `ProblemFile.f90.in`: problem file template. It sets the initial condition, adds the
  source term and writes `mms_l2_error.dat` with the total and per-variable L2 errors,
  normalised by the mesh volume (integral of 1).
* `<case>/control_template.control`: control file of each case.
* `sweeps.json`: default orders and meshes of each case, used by `run.py` and the workflow.
* `<case>/reference.csv`: errors of the default sweep of the committed case, the same as
  in horses3d (gfortran, RELEASE). After changing a case on purpose, copy its `errors.csv`
  here.

## Results

Orders of convergence between meshes 8 and 16, gfortran, the same as in horses3d. Each
cell gives the order of the total error and, in brackets, the lowest order of any single
variable. The job summary of each run gives the order of every variable between
consecutive meshes, and the lowest of them.

| P | NS | MU |
|---|---|---|
| 3 | 3.18 (3.08) | 3.03 (2.47) |
| 4 | 4.31 (3.92) | 3.96 (2.52) |
| 5 | 5.18 (5.04) | 5.03 (4.95) |
| 6 | 6.75 (6.16) | 6.17 (6.14) |

* These are the orders of short runs between two meshes, not asymptotic rates. The runs
  end at t = 0.02 in NS (200 steps), and at t = 2e-4 in MU, whose explicit time step is
  limited by Cahn-Hilliard. The initial condition is exact at the nodes, where the error
  is measured, so the whole error is made during the run by the truncation error of the
  spatial operator at the nodes, which is the residual printed at iteration 0 of each log
  (maximum norm). Between meshes 8 and 16 its order is between P - 1.6 and P + 0.2 in NS
  (lowest at even P, from BR1), and lower in MU (on meshes 12 and 16, about 2 at P = 4 and
  4 at P = 6). The dissipation of the scheme damps part of this error during the run,
  faster on finer meshes, so the orders change with the final time and with the meshes.
* Every variable converges at order P - 0.1 or higher between meshes 8 and 16, except in
  MU: the momentum at P = 3 and 4, and c at P = 4. In MU the total is mostly the pressure
  error (at least 92% of the squared total on mesh 16), so the order of the total is the
  order of the pressure.
* In horses3d, one-off runs on finer meshes with the same final time give these orders:
  NS between meshes 16 and 32, 3.25 (3.15), 4.75 (4.61), 5.11 (5.00), 7.23 (6.88) for P = 3
  to 6; MU between meshes 16 and 24, 3.05 (2.60), 4.19 (4.18), 5.08 (4.95), 6.30 (6.24),
  and between 24 and 32, 3.06 (2.73) and 4.45 (3.62) for P = 3 and 4. Runs five times
  longer give NS 3.51, 5.04, 5.12, 7.81 between meshes 8 and 16. Halving the time step
  changes the errors at P = 6 on mesh 16 by at most 2.6e-8 (NS) and 1.1e-6 (MU), so they
  are errors of the spatial discretization. To run longer, raise `Number of time steps`.
* NS runs at Re = 10, so that the viscous terms are about a fifth of the momentum source.
  The two MU fluids have different densities, viscosities and sound speeds, so that the
  mixture laws, which weight each property by the concentration, are exercised. The MU
  errors of this repository agree with those of horses3d within the run-to-run noise of
  its multi-threaded MU runs (up to about 1e-4 relative).
* Not tested: Riemann solvers other than matrix dissipation in NS and the exact one in
  MU; two-point fluxes other than Chandrasekar; BR2 and IP; Gauss nodes; time integrators
  other than RK3; non-periodic boundary conditions; curved, non-affine or non-conforming
  meshes; gravity; MPI; the OpenACC paths of this repository.

## Assumptions

* NS uses the default Sutherland constants, so do not set `sutherland temperature` or
  `sutherland reference temperature`. Gravity is not included.
* The time steps of the control files are stable up to P = 6 on mesh 16, the finest
  pair of the default sweeps, where twice the step is still stable; higher orders or
  finer meshes may need a smaller `dt`, with more steps to keep the final time. In MU the
  explicit step is limited by Cahn-Hilliard: the limit scales like h^4, falls steeply
  with P and grows with the chemical characteristic time. With 20000 s (a small
  mobility), dt = 1e-6 is stable; at P = 6 on mesh 16, 2.5 times this step is stable and
  5 times diverges (measured in horses3d).
* The default sweeps start at P = 3: at P = 2 the truncation error of the MU momentum
  does not converge.
* The error sums are rank-local, so run without MPI.
