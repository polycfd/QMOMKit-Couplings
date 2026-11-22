# QMOMKit-Couplings

Reference couplings of [QMOMKit](https://github.com/polycfd/QMOMKit) to
physical-space transport solvers. QMOMKit handles the phase-space part of a
moment method (moment inversion, closures and evolution of the quadrature). The
host solver transports the moments through physical space. 

## Versioning

QMOMKit-Couplings is released alongside [QMOMKit](https://github.com/polycfd/QMOMKit):
each release of QMOMKit gets a matching release of this repository with the
same version number. Release `vX.Y.Z` of QMOMKit-Couplings has been tested
against release `vX.Y.Z` of QMOMKit, so check out the tag that matches your
QMOMKit version.

## Couplings

| Directory                                     | Solver              | Moment method | Physical-space discretization                                                                       |
| --------------------------------------------- | ------------------- | ------------- | --------------------------------------------------------------------------------------------------- |
| [DG-dealii-v9-8](DG-dealii-v9-8/README.md) | `qmomAdvectionDG`   | QMOM          | discontinuous Galerkin `FE_DGQ(1)` ([deal.II](https://dealii.org/) 9.8), SSP-RK3, Zhang–Shu limiter |
| [HyQMOM-FV-1D](HyQMOM-FV-1D/README.md)     | `hyqmomTransport1D` | HyQMOM        | finite volume, first order, kinetic or HLL fluxes                                                   |
| [OpenFOAM-v14](OpenFOAM-v14/README.md)     | `qmomAdvectionFoam` | QMOM          | finite volume (OpenFOAM-14), implicit upwind                                                        |

The OpenFOAM and deal.II couplings solve a population balance for an internal
coordinate, such as the particle size, which a physical velocity carries through
space. The HyQMOM coupling transports the moments of the particle velocity
itself, as in a kinetic equation. The README of each coupling describes its
model, solution strategy, cases, requirements, building and running.

## Authors

- [Bruno Blais](https://www.polymtl.ca/expertises/en/blais-bruno)
- Bruna Campos
- [Fabian Denner](https://www.polymtl.ca/expertises/en/denner-fabian)

## License

Copyright 2026 The QMOMKit-Couplings developers.

QMOMKit-Couplings is released under the Apache License 2.0. See [LICENSE](LICENSE) for details.
