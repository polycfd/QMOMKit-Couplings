# QMOMKit-Couplings

Reference couplings of [QMOMKit](https://github.com/polycfd/QMOMKit) to
physical-space transport solvers. QMOMKit handles the phase-space part of a
population balance (moment inversion and evolution of the quadrature). The host
solver transports the moments through physical space. Each coupling includes
verification cases with analytical solutions.

## Equations

A number density `f(x, ξ, t)` over physical space `x` and an internal
coordinate `ξ` (e.g. particle size) obeys the population balance equation

```
  ∂f/∂t + ∇ₓ·( u f ) + ∂/∂ξ ( G(ξ) f ) = 0
```

with `u(x)` the physical velocity and `G(ξ) = dξ/dt` the phase-space velocity
(growth law). Instead of solving for `f`, the solvers transport its raw moments
`M_k = ∫ ξᵏ f dξ`, `k = 0 … 2n−1`:

```
  ∂M_k/∂t + ∇ₓ·( u M_k ) = S_k,      S_k = k ∫ ξ^(k−1) G(ξ) f dξ
```

`S_k` cannot be computed from the moments alone. The Quadrature Method of
Moments (QMOM) closes it with an `n`-node quadrature `f ≈ Σᵢ wᵢ δ(ξ − ξᵢ)`,
whose weights and abscissas QMOMKit recovers from the `2n` moments with the
Wheeler algorithm, giving

```
  S_k = Σᵢ wᵢ k G(ξᵢ) ξᵢ^(k−1)
```

The deal.II coupling also supports additional moment source terms, such as
multiplication and binary fragmentation.

## Solution strategy

- **Operator splitting (Lie).** Each time step first transports every moment
  through physical space, then advances the phase space cell by cell.
- **Phase space (QMOMKit).** Instead of integrating `S_k`, QMOMKit moves the
  quadrature nodes along `dξᵢ/dt = G(ξᵢ)` with fixed weights, then rebuilds the
  moments from the advanced quadrature. This step is realizable by
  construction; only the inversion of the incoming moments can fail.
- **Physical space (host solver).** All moments use the same bounded scheme,
  so the transported moment vector stays in the (convex) realizable set.

## Couplings

| Directory | Solver | Physical-space discretization | Cases |
|---|---|---|---|
| [`OpenFOAM-v14/`](OpenFOAM-v14/README.md) | `qmomAdvectionFoam` | finite volume (OpenFOAM-14), implicit upwind | A–D |
| [`DG-dealii-v9-8/`](DG-dealii-v9-8/README.md) | `qmomAdvectionDG` | discontinuous Galerkin `FE_DGQ(1)` ([deal.II](https://dealii.org/) 9.8), SSP-RK3, Zhang–Shu limiter | A–F |

## Verification problems

All start from a Gaussian phase-space distribution `ξ ~ N(10, 1)`.

- **A**: growth only (`G = ξ`, `u = 0`); exact solution `M_k(t) = M_k(0) eᵏᵗ`. Tests the QMOM closure on its own.
- **B**: pure advection of a Gaussian blob (`u = 1`, no growth).
- **C**: advection plus growth; comparing with B separates the two error sources.
- **D**: 2-D Taylor–Green vortex plus growth; reference solution obtained by tracing characteristics backwards.
- **E, F** (deal.II only): growth plus a multiplication (E) or binary-fragmentation (F) source term; F has no analytical solution.

## Requirements

- An installed QMOMKit (static C++17 library), passed to each build through
  `QMOMKIT_DIR`.
- OpenFOAM-14 and/or deal.II ≥ 9.8.

The build, run and configuration instructions are in each coupling's README.

## License

Apache License 2.0; see [`LICENSE`](LICENSE).
