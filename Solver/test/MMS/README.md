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
files and logs of the runs. The default sweep, P = 3 to 6 on meshes 4, 8 and 16 for NS
and P = 2 to 6 on meshes 4, 6, 8, 12 and 16 for MU, takes about 45 minutes, set by MU
at P = 6 on one thread.

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

Orders of the total error between meshes 8 and 16, gfortran, the same as in horses3d:

| P | NS | MU |
|---|---|---|
| 2 |      | 1.99 |
| 3 | 3.18 | 3.00 |
| 4 | 4.31 | 4.00 |
| 5 | 5.18 | 5.01 |
| 6 | 6.75 | 6.03 |

In MU the total error is that of the pressure, one to three orders of magnitude above
the others. Orders of each MU variable between meshes 8 and 16 (the job summary lists
every variable):

| P | c | sqrt(rho) u | p |
|---|---|---|---|
| 2 | -0.69 | 0.69 | 1.99 |
| 3 | 1.88 | 2.72 | 3.00 |
| 4 | 1.73 | 2.86 | 4.00 |
| 5 | 5.25 | 5.10 | 5.01 |
| 6 | 6.10 | 6.08 | 6.03 |

* The runs are short: t = 0.02 in NS (200 steps), and t = 2e-5 in MU, whose time step is
  limited by Cahn-Hilliard. The initial condition is exact at the nodes, where the error
  is measured, so the error comes from the time evolution alone. At first it grows like
  t times the truncation error of the spatial operator at the nodes, and it approaches
  the O(h^(P+1)) error of DG only later, as the scheme damps it. For the first
  derivatives that truncation error converges at order P, so the orders of the total
  error are about P: P in MU, the shortest run, whose error grows linearly in time, and
  between P and P + 0.8 in NS.
* The truncation error of higher derivatives converges more slowly, and in MU the run is
  too short for the scheme to damp it: c, through the fourth derivatives of Cahn-Hilliard,
  and the momentum, through the viscous and capillary terms, converge at lower orders up
  to P = 4, and at P = 2 the error of c grows on the finer meshes.
* In horses3d, runs five times longer give orders of 3.51 in NS and 3.04 in MU at P = 3
  (2.09 for c in MU), at five times the cost, and halving the time step changes the
  errors at P = 6 on mesh 16 by about 1e-5 at most, so they are errors of the spatial
  discretization. To run longer, raise `Number of time steps`.
* NS runs at Re = 10, so that the viscous terms are about a fifth of the momentum source.

## Assumptions

* NS uses the default Sutherland constants, so do not set `sutherland temperature` or
  `sutherland reference temperature`. Gravity is not included.
* The time steps of the control files are stable up to P = 6 on mesh 16, the finest
  pair of the default sweeps; higher orders or finer meshes may need a smaller `dt`, with
  more steps to keep the final time. In MU the explicit step is limited by Cahn-Hilliard,
  which scales like h^4/P^8; the chemical characteristic time of 2000 s (a small
  mobility) keeps dt = 1e-7 stable.
* The error sums are rank-local, so run without MPI.
