# QMOM with advection in deal.II

Finite-element Discontinuous Galerkin (DG) transport in physical space ([deal.II](https://dealii.org/) version 9.8) coupled to the Quadrature Method of Moments in phase space ([QMOMKit](https://github.com/polycfd/QMOMKit)). Ships as the solver `qmomAdvectionDG` plus five validation cases.

## Model

A number density $n(x, \xi, t)$ obeys

```math
  \frac{\partial n}{\partial t} + \nabla_x \cdot (\mathbf{u}\, n) + \frac{\partial}{\partial \xi} [G(\xi)\, n] = 0
```

The raw moments $M_k = \int \xi^k \, n d\xi$ are transported instead:

```math
  \frac{\partial M_k}{\partial t} + \nabla_x \cdot (\mathbf{u} \, M_k) = S_k
```

QMOMKit closes $S_k$ by Wheeler inversion of the first $2n$ moments. Rather
than integrate $S_k$, the phase-space half advances the quadrature itself: the
nodes move along their trajectories and the moments are reconstructed from the
advanced quadrature. The phase-space evolution is done by the following call:

```cpp
qmomkit::advance_qmom_phase_space(qmomCase.initial_state(), time, dt, phase_space_velocity, source, options)
```

Advecting nodes instead of moments is what makes the step realizable: the
weights are carried unchanged and the abscissas only move, so the result cannot
leave the realizable cone. Only the inversion of the incoming moments can fail.

## Problems

Six problem types (A-F) are available. All start from $X ~ N(10, 1)$; problems B, C, and D also carry a Gaussian blob in space ($σ_space = 0.25$). Nodes $n = n_moments/2$. Problems E and F also include source terms. Problems A to E have an analytical solution, and in these cases the solver reports a true absolute L2 error. Problems A, D, and F have a linear growth velocity in phase space, while problem C has a constant growth velocity.

|       | Physical space      | Phase advection | Phase moment source                                | Mesh               | $N$ | $\delta t$ | t_end |
| ----- | ------------------- | --------------- | -------------------------------------------------- | ------------------ | --- | ---------- | ----- |
| **A** | none                | growth $G=ξ$    | none                                               | 1 cell, periodic   | 10  | 0.002      | 2     |
| **B** | advection $u=1$     | none            | none                                               | 50 cells, periodic | 10  | 0.002      | 2     |
| **C** | advection $u=1$     | growth $G=0.1$  | none                                               | 50 cells, periodic | 10  | 0.002      | 2     |
| **D** | Taylor–Green vortex | growth $G=ξ$    | none                                               | 50×50              | 6   | 0.005      | 20    |
| **E** | none                | growth $G=ξ$    | multiplication $S=\alpha \xi^k$                    | 1 cell, periodic   | 10  | 0.002      | 2     |
| **F** | none                | growth $G=ξ$    | binary fragmentation $S=\alpha (2^(k-1) - 1)\xi^k$ | 1 cell, periodic   | 10  | 0.002      | 2     |

Problem A validates QMOM alone: with $\mathbf{u} = 0$ the advection operator vanishes and the closure is exact ($M_k(t) = M_k(0) \exp^{kt}$). Problem B is the inverse, pure transport. Problem C is both, so comparing it against B separates the two error contributions. Problem D's exact solution is the initial field at the foot of the characteristic, traced backwards with RK4 in `exactSolution::traceBackToOrigin()`. Problem E validates QMOM including a source term, and F showcases a QMOM-alone case with binary fragmentation (no analytical solution).

## Discretisation

**Space.** `FE_DGQ(1)` in an `FESystem`, one component per moment. `QGauss(2)`.
Standard upwind DG: volume term $\mathbf{u} \cdot \nabla \psi_i M$ plus an upwind interior-face flux. No
boundary work — A/B/C/E/F are periodic, D has zero normal velocity. DG mass matrices are block-diagonal, so no global system is assembled; each cell's local mass matrix is inverted in the `mesh_loop` copier.

**Time.** Explicit SSP-RK3, with the **Lie-split** operator: the three
physical-space stages run first, then `QMOMKit` is responsible for evolving the phase space. First-order in time overall — Strang splitting would be the natural upgrade.

## Boundedness

Unlimited DG overshoots at steep gradients (problem B: +2.1% above and −59%
below the initial range after 800 steps), and moments that leave their range
tend to leave the _realizable_ cone, at which point the inversion has no
solution and the run aborts. Four knobs, all members of `QMOM1DPhaseSpace`:

| knob                    | default                | what it does                                     |
| ----------------------- | ---------------------- | ------------------------------------------------ |
| `use_limiter`           | `true`                 | Zhang–Shu limiter after every advection RK stage |
| `enforce_realizability` | `true`                 | limiter also guarantees the state is invertible  |
| `qmom_options.adaptive` | `false` (`true` for D) | Wheeler drops nodes instead of failing           |
| `repair_broken_cells`   | `false`                | last-resort rebuild of non-invertible cells      |

**Limiter.** Each cell polynomial is squeezed towards its own cell average,
$\tilde{u}_K = \tilde{u}_K + \theta (u_K − \tilde{u}_K)$, with $\theta$ the largest factor keeping it inside the
pointwise range the solution had entering the timestep. The average is
untouched, so the limiter is exactly conservative (verified to 6e-16), and
$\theta = 1$ in smooth regions. Result on B: zero overshoot and undershoot, at a
modest accuracy cost (moment 0 at $t=2$: 0.0099 unlimited → 0.0118 limited).

A **single $\theta$ is shared by all moments**, deliberately — a per-moment factor
distorts the ratios between moments and breaks realizability even when each
moment is individually in range. With one $\theta$ the result is a convex
combination, and the realizable set is convex.

**Realizability safeguard.** Halves $\theta$ until the state can actually be
inverted, falling back on the cell average. Prints `realizability safeguard
active in N cells` when it engages; silence is the normal case. The gate is
`compute_qmom(...).ok()`, _not_ `check_realizability()` — the latter rejects
even the exact analytical ten-moment Gaussian sequence, and gating on it
flattens the solution onto its cell averages (error 0.0118 → 0.209).

**Adaptive inversion.** A negative recurrence coefficient means the sequence has
lost realizability and no $n$-node quadrature exists. `qmom_options.adaptive`
makes Wheeler reduce the node count until one does, instead of failing. Not
free: it also fires on merely ill-conditioned sequences, inflating problem A's
high-moment error by 12% ($m_2$) to 68% ($m_9$). Off by default; on for D.

**Repair.** `repair_broken_cells` rewrites a cell whose moments cannot be
inverted with the nominal phase-space shape rescaled to the $m_{\theta}$ it already
holds — realizable by construction, since the realizable set is a cone. $m_{\theta}$ is
conserved exactly; the higher moments are overwritten. $m_{\theta, floor}$ additionally
clamps $m_{\theta}$. This cannot fail, which is why it is the last line of defence, but
a run leaning on it is under-resolved.

One failure mode none of this fixes: if $m_{\theta}$ decays below QMOMKit's `cutoff`
($10 \times 10^{-30}$) in the far field, even the cell average is unusable. With
$\sigma_{space} = 0.25$ the smallest $m_{\theta}$ on the mesh is ~$1.3 \times 10^{-9}$; at $\sigma_{space} = 0.1$
the initial condition is already unrepresentable at $2 \times 10^{-56}$.

## Building

Needs deal.II ≥ 9.8.0 and an installed QMOMKit (a dependency-free static C++17 library). Neither path is hard-coded:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DDEAL_II_DIR=/path/to/dealii/inst \
  -DQMOMKIT_DIR=/path/to/qmomkit/inst
cmake --build build -j
```

`QMOMKIT_DIR` is the install prefix. Linking happens after
`deal_ii_invoke_autopilot()`, which is what creates the target, and must use the
keyword-less `target_link_libraries` signature — see the comment in
`CMakeLists.txt`.

## Running and output

```bash
cd build && ./qmomAdvectionDG A     # or B, C, D, E, F — defaults to A
```

Each step prints the L2 error of every moment. Note the high
moments are enormous in absolute terms ($m_9 \approx 9 \times 10^{16}$ at $t=2$), so read those
errors relative to the moment values.

The `output/` folder contains a file with the computed L2 error and a file with the moments, both at each time step. VTK files `solution<P>-<n>.vtk` are also stored there, containing one scalar field $u_k$ per moment $M_k$. The variable `output_frequency` controls how often; the initial condition and the final step are always written. The `plot_vtk` script will plot the moment distribution over the physical space at the final simulation time.

```bash
python plot_vtk.py <directory> <prefix>
python plot_vtk.py build/output solutionB
```

## Code map

| Item                                      | Location                                                                    |
| ----------------------------------------- | --------------------------------------------------------------------------- |
| Growth laws                               | `linear_growth_velocity()`, `no_growth_velocity()`                          |
| Moment source terms                       | `binary_fragmentation_source()`, `multiplication_source()`                  |
| Initial conditions / analytical solutions | `ConstantPhysicalGaussianPhase` and siblings                                |
| Backward characteristic trace (D)         | `SwirlGaussianPhysicalGaussianPhase::trace_back_to_origin()`                |
| Velocity fields                           | `NoFlow`, `ConstantXFlow`, `TaylorGreenVortex`                              |
| Solver, all the knobs                     | `QMOM1DPhaseSpace`                                                          |
| Problem setup                             | constructor + `setup_mesh()` + `setup_system()`                             |
| **Phase-space step — the QMOMKit call**   | `rk_step_phase_space()`                                                     |
| Physical-space DG step                    | `rk_step_physical_space()`                                                  |
| Limiter, safeguard, repair                | `apply_zhang_shu_limiter()`, `moments_are_invertible()`, `repair_moments()` |
| Time loop and splitting                   | `solve()`                                                                   |
| Error measurement                         | `calculate_L2_error()`                                                      |

## Known limitations

- Lie splitting, so first-order in time regardless of the SSP-RK3 stages.
- No CFL control: $\delta t$ is fixed per problem.
- The limiter is written for `FE_DGQ(1)`, whose nodal points coincide with the
  Gauss–Lobatto points the Zhang–Shu argument needs. A higher degree would need
  a larger test-point set.
- D's reference costs a backward trace per evaluation point, and its cost grows
  linearly with $t$, so late steps are markedly slower than early ones.
- `print_solution_vector()` exists as a debugging aid; its call site in
  `solve()` is commented out.
