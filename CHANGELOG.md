# QMOMKit-Couplings Changelog

All notable changes to this project will be documented in this file.

Each release has the same version number as the
[QMOMKit](https://github.com/polycfd/QMOMKit) release it has been tested
against; see the [README](README.md#versioning).
The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- `hyqmomTransport1D` in `HyQMOM-FV-1D/`, a one-dimensional coupling of the
  HyQMOM of QMOMKit: first-order finite-volume transport of the velocity
  moments, with the fluxes built from the weighted points, the characteristic
  speeds and the closure moment that QMOMKit returns for each cell. The flux
  is the kinetic flux of the weighted points, the kinetic flux of the
  characteristic speeds, or the HLL flux. It needs the characteristic speed
  weights and the closure moment that QMOMKit added in version 1.1.
- Verification cases of `hyqmomTransport1D` with exact solutions: the cases
  of Li & Marchisio (2022, Sections 3.1 and 3.2) with one, two or three
  velocities and with crossing groups of particles, a variant with normally
  distributed velocities, and the Riemann problem of Fox & Laurent (2022,
  Section 6). A plotting script, `plot_hyqmom.py`, and a script that runs
  every case and plots the results, `run_all.sh`.

## [1.0.0] - 2026-11-07

First public release, tested against [QMOMKit v1.0.0](https://github.com/polycfd/QMOMKit/releases/tag/v1.0.0).

### Added

- `qmomAdvectionDG`, a coupling of QMOMKit to deal.II 9.8: discontinuous
  Galerkin transport of the moments with `FE_DGQ(1)` and SSP-RK3, Lie-split
  with the QMOMKit phase-space step. A Zhang–Shu limiter with a single factor
  shared by all moments, a realizability safeguard, optional adaptive
  inversion and repair. VTK output and a plotting script, `plot_vtk.py`.
- `qmomAdvectionFoam`, a coupling of QMOMKit to OpenFOAM-14: finite-volume
  transport of the moments with an implicit upwind scheme, Lie-split with the
  QMOMKit phase-space step. Optional adaptive inversion and repair of cells
  whose moments cannot be inverted, set in `constant/qmomProperties`.
- Verification problems with analytical solutions, starting from a Gaussian
  phase-space distribution: growth only (A), advection only (B), advection
  with constant growth (C) and a 2-D Taylor–Green vortex with growth (D), in
  both couplings. The deal.II coupling adds growth with a multiplication
  source (E) and with binary fragmentation (F, no analytical solution).
- L2 errors of every moment against the analytical solution, reported by
  both solvers.

[Unreleased]: https://github.com/polycfd/QMOMKit-Couplings/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/polycfd/QMOMKit-Couplings/releases/tag/v1.0.0
