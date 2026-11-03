# MMS convergence tests

Method of manufactured solutions for the compressible Navier-Stokes (`NS`) and
multiphase (`MU`) solvers. A smooth periodic solution is imposed with a source term, and
the L2 error at the final time gives the order of convergence under h-refinement.

The user source term is a host routine, so the tests run on CPU builds (gfortran,
single-threaded); with OpenACC the solvers do not call it. The scripts are the same as
in horses3d, where the other equation sets have cases too.

## Run on GitHub

Actions > MMS convergence > Run workflow. Choose the case (or all), the polynomial
orders, the meshes and, optionally, changes to the control file such as
`Reynolds number = 100; Riemann solver = Roe`. Each case builds its solver, writes the
table of errors and orders to the job summary, and uploads `errors.csv`,
`convergence.png` and the control files and logs of the runs.

The workflow also runs, with the default inputs, on pushes that change the tests.
GitHub shows the Run workflow button only for workflows on the default branch; on
another branch, use the GitHub CLI:

```
gh workflow run mms.yml --ref <branch> -f case=NS -f P="2 3 4" -f meshes="4 8 16" -f set="Reynolds number = 100"
```

## Run locally

Build the solver (`./configure && make ns COMPILER=gfortran ENABLE_THREADS=NO` in `Solver`), then:

```
pip install sympy gmsh matplotlib
python3 run.py NS                                   # P = 2 3 on meshes 4 8
python3 run.py NS --P 2 3 4 5 --meshes 4 8 16       # full study
python3 run.py NS --set "Reynolds number = 100"     # change the control file
python3 run.py MU --check                           # compare with MU/reference.csv
python3 plot.py NS                                  # NS/convergence.png
```

`run.py` generates `<case>/SETUP/ProblemFile.f90`, builds it with the HORSES3D build
settings (`COMPILER`, `ENABLE_THREADS`, ... from the environment and `Solver/make.inc`),
and runs every mesh and P with `<case>/control_template.control`, filling in the mesh,
the solution file and the order. Everything else, including `dt` and the number of
steps, comes from the control file; edit it or override keys with `--set`. Results go
to `<case>/errors.csv`, and pairs already run with the same solver, problem file and
control file are skipped. `--check` fails if a run fails or, with the committed
control file, an error differs from `<case>/reference.csv` by more than 1e-3.

The meshes are periodic cubes [-1,1]^3 with N^3 elements, made by the cartesian-box
generator of [horses3d-meshes](https://github.com/horses-framework/horses3d-meshes)
(Gmsh Python API). `run.py` downloads the script at a fixed commit into `MESH/`, or uses
`$HORSES3D_MESHES/cartesian-box/` if that variable points to a clone. No mesh files
are kept in this repository.

## Files

* `generate.py`: manufactured solutions and equations of each case. The source term is
  the residual of the solution, emitted without symbolic simplification. The physical
  parameters are read from the solver at run time.
* `ProblemFile.f90.in`: problem file template. It sets the initial condition, adds the
  source term and writes `mms_l2_error.dat` with the total and per-variable L2 errors,
  normalised by the mesh volume (integral of 1).
* `<case>/control_template.control`: control file of each case.
* `<case>/reference.csv`: errors of the committed case, the same as in horses3d (gfortran,
  RELEASE). After changing a case on purpose, copy its `errors.csv` here.

## Results

Orders of the total error between meshes 8 and 16 (theory P+1), gfortran:

| P | NS | MU |
|---|---|---|
| 2 | 2.65 | 2.00 |
| 3 | 3.90 | 3.04 |
| 4 | 4.87 | 4.06 |
| 5 | 5.69 | |

* The meshes up to 16 are not fully in the asymptotic range: between 16 and 32 elements
  per side NS gives 2.95, 3.89, 4.72 and 5.54.
* MU converges closer to P because the artificial compressibility multiplies the
  divergence error by c0^2. At P = 2 the momentum and the concentration converge more
  slowly than the pressure, which dominates the total.
* NS runs at Re = 1e5, where the viscous terms are a small part of the source term.

## Assumptions

* NS uses the default Sutherland constants, so do not set `sutherland temperature` or
  `sutherland reference temperature`. MU assumes per-fluid sound speeds, not
  `artificial sound speed square`. Gravity is not included.
* Keep the final time: with u ~ sin(t), as in NS, a short run leaves the flow nearly at
  rest, the source is dominated by du/dt and the measured order drops.
* The error sums are rank-local, so run without MPI.
