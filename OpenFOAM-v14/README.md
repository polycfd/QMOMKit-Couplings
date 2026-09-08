# QMOM with advection in OpenFOAM

Finite-volume transport in physical space (OpenFOAM-14) coupled to the
Quadrature Method of Moments in phase space (QMOMKit). Ships as the solver
`qmomAdvectionFoam` plus four validation cases under `run/`.

## Model

A number density `f(x, ξ, t)` obeys

```
  ∂f/∂t + ∇ₓ·( u f ) + ∂/∂ξ ( G(ξ) f ) = 0
```

The raw moments `M_k = ∫ ξᵏ f dξ` are transported instead:

```
  ∂M_k/∂t + ∇ₓ·( u M_k ) = S_k
```

QMOMKit closes `S_k` by Wheeler inversion of the first `2n` moments. Rather
than integrate `S_k`, the phase-space half advances the quadrature itself: the
nodes move along their trajectories and the moments are reconstructed from the
advanced quadrature. Two calls in `phaseSpaceStep.H`:

```cpp
qmomkit::advance_qmom_continuous(qmomCase.initial_state(), t, dt, growthModel, options)
```

Advecting nodes instead of moments is what makes the step realizable: the
weights are carried unchanged and the abscissas only move, so the result cannot
leave the realizable cone. Only the inversion of the incoming moments can fail.

Moments are dimensionless (`[0 0 0 0 0 0 0]`). Operators are Lie-split:
implicit `ddt`+`div` advection, then node advection (SSP-RK3) in phase space —
first-order in time overall.

## Problems

One case per problem under `run/`. All start from `X ~ N(10, 1)`; B, C, D also
carry a Gaussian blob in space (`σ_space = 0.25`). Nodes `n = N/2`. Every
problem has an analytical solution, so the solver reports a true L2 error, both
absolute and relative — the relative values are the meaningful ones, since the
high moments are enormous in absolute terms.

| | Physical space | Phase space | Mesh | N | Δt | t_end |
|---|---|---|---|---|---|---|
| **A** | none | growth `G=ξ` | 5 cells, periodic | 10 | 0.0025 | 2 |
| **B** | advection `u=1` | none | 200 cells, periodic | 10 | 0.0025 | 2 |
| **C** | advection `u=1` | growth `G=ξ` | 200 cells, periodic | 10 | 0.0025 | 2 |
| **D** | Taylor–Green vortex | growth `G=ξ` | 50×50 | 6 | 0.005 | 20 |

A validates QMOM alone: with `u = 0` the advection operator vanishes and the
closure is exact (`M_k(t) = M_k(0) e^{kt}`). B is the converse, pure transport.
C is both, so comparing it against B separates the two error contributions.
D's exact solution is the initial field at the foot of the characteristic,
traced backwards with RK4 in `exactSolution::traceBackToOrigin()`.

## Boundedness

All moments must use the **same** `div` scheme, and it should give face values
that are convex combinations of cell values: `upwind` does (donor-cell moment
vector, and the realizable set is convex), so it is the default everywhere. A
per-moment limiter distorts ratios between moments and can leave the realizable
cone.

| `div(phi,M)` | realizability |
|---|---|
| `Gauss upwind` | safe by construction |
| `Gauss limitedLinear 1` | held on A–D, but limits per moment |
| `Gauss vanLeer` | held on A–D, but limits per moment |
| `Gauss linear` | unbounded |

The higher-order schemes are noticeably more accurate but carry the caveat
above — switch with eyes open, and watch for the repair message. What `upwind`
buys is boundedness for free: no over- or undershoot in `M₀`, and `∫M₀ dV`
conserved to round-off on a periodic box. The price is peak decay from the
numerical diffusion.

Two fallbacks in `constant/qmomProperties` cover the inversion, both off by
default (on for D):

- `qmomOptions/adaptive` — Wheeler drops nodes instead of failing. Also fires
  on merely ill-conditioned sequences.
- `repair/repairBrokenCells` — rebuilds a non-invertible cell from the nominal
  shape rescaled to its `m₀` (`m₀` conserved, higher moments overwritten), then
  retries once. `repair/m0Floor` additionally clamps `m₀`. A run leaning on this
  is under-resolved and says so.

Without either, a failed inversion is a `FatalError` naming the cell.

## Building

```bash
. $HOME/openfoam/OpenFOAM-14/etc/bashrc
QMOMKIT_DIR=$HOME/QMOMkit/inst ./Allwmake
```

`QMOMKIT_DIR` must point at a QMOMKit install prefix. The executable lands in
`$FOAM_USER_APPBIN`. `Make/options` appends `-std=c++17` (QMOMKit needs it,
OpenFOAM-14 defaults to C++14) and links `-lqmomkit` statically.

## Running

```bash
cd run/problemC && ./Allrun          # or problemA, problemB, problemD
```

`Allrun` runs `blockMesh` then `qmomAdvectionFoam`; `Allclean` resets one case
and `./Allwclean` at the top level resets all of them and removes the build
artefacts. Output is ordinary time directories with `M0 … M<N-1>`, `U`, `phi` —
open `case.foam` in ParaView. `errorAnalysis/reportEvery` throttles the error
report (D reports every 200 steps, since its exact solution costs a backward
trace per cell).

## Configuring a case

Everything is in `constant/qmomProperties`:

```
nMoments        10;            // must be even; nNodes defaults to nMoments/2

phaseSpace      { growthModel linear; growthRate 1; }  // linear | constant | none

qmomOptions     { adaptive false; rmin 1e-10; eabs 1e-10; cutoff 1e-30; }

initialDistribution
{
    mu 10; sigma 1;
    spatialShape gaussianBlob;  // uniform | gaussianBlob
    sigmaSpace 0.25; centre (0 0 0);
}

velocityModel   { type constantX; U0 1; }   // noFlow | constantX | taylorGreen

errorAnalysis   { type advectedGaussian; }  // none | constantSpace
                                            // | advectedGaussian | swirlGaussian

repair          { repairBrokenCells false; m0Floor 0; }
```

The moment fields are not stored one file per moment: the solver builds
`M0 … M<N-1>` from the template field `0/M` and overwrites the internal field
from `initialDistribution`. `0/U` likewise only fixes boundary types; values
come from `velocityModel`. `phi` is built from the analytical velocity at face
centres, which makes D's boundary flux exactly zero (so `zeroGradient` needs no
boundary work) and keeps the discrete divergence at round-off.

## Code map

| Item | Location |
|---|---|
| Growth laws, velocity fields, exact solutions | `qmomModels.H` |
| Field and dictionary setup | `createFields.H` |
| Physical-space step | `advectionStep.H` |
| Phase-space step (QMOMKit call, repair) | `phaseSpaceStep.H` |
| Error measurement | `momentErrors.H` |
| Time loop and splitting | `qmomAdvectionFoam.C` |
