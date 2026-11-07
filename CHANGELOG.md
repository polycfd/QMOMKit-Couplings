# QMOMKit-Couplings Changelog

All notable changes to this project will be documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Each release has the same version number as the
[QMOMKit](https://github.com/polycfd/QMOMKit) release it has been tested
against; see the [README](README.md#versioning).

## [Unreleased]

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
