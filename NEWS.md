# Changelog

This file documents notable changes in the repository.

## Release Notes — v1.0
 
This release brings major updates to CI coverage, solver capabilities, MPI/mesh workflows, and overall robustness.
 
### 🚀 Added
- Code coverage reporting integration (Codecov) for PR validation.
- Ability to write/read mesh partitioning files for MPI runs, avoiding repartitioning on each execution.
- New and expanded regression cases (including additional channel/cylinder configurations and limiter-focused tests).
- Additional Navier–Stokes flux/Riemann-solver configurations with test coverage.
- GPU-enabled and validated explicit time-integration paths for Euler, RK5, SSPRK33, and SSPRK43.

### 🔧 Improved
- The horses3d and horses3d-gpu repositories are now hosted in a GitHub organization
  https://github.com/horses-framework/
- CI workflows for CPU/GPU, serial/parallel, and MPI scenarios, with broader coverage and better reliability.
- Parallel preprocessing flow for HOPR/SFC partitioning-related workflows. Modified the preprocessing of HOPR HDF5 meshes. The read-in is now MPI-parallel when choosing `partitioning = SFC`
- Handling of anisotropic polynomial-order scenarios.
- Solver documentation and run-configuration defaults alignment.
 
### 🐛 Fixed
- Wall-function connectivity issues when interfaces cross MPI partitions.
- Channel-forcing behavior in parallel runs.
- GPU/IBM and MPI execution/workflow issues affecting CI and runtime stability.
- PointLinkedList memory-safety issues (including double-free and freed-memory access cases).
- Multiple correctness fixes and typo-level cleanups across solver and docs.
 
### 📝 Documentation
- Updated README and user manual sections (including partitioning-from-files workflow and default flag clarifications).
- Added/updated NEWS and project references for organization migration.
 
### Summary
`main` is updated with all current `develop` advancements: wider automated validation, new runtime capabilities, and a substantial set of stability/correctness fixes across GPU, MPI, and solver components.

