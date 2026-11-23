# HyQMOM transport in one dimension

Transport of the velocity moments of a number density in one space dimension, with the fluxes in physical space built from the HyQMOM of
[QMOMKit](https://github.com/polycfd/QMOMKit), with first-order finite volumes and forward Euler. The code shows how a host obtains the closure of the
flux, the characteristic speeds and the weighted points from QMOMKit, and how it builds three different fluxes from them. `hyqmomTransport1D.cc` holds
only `main()`, with the time loop and every call to QMOMKit, and the other files hold the parts that do not involve QMOMKit. The cases are those of Li
& Marchisio (2022, Sections 3.1 and 3.2), a variant of them with normally distributed velocities, and the Riemann problem of Fox & Laurent (2022,
Section 6). Every case has an exact solution.

The coupling needs the characteristic speed weights and the closure moment of `compute_hyqmom_characteristic_speeds(...)`, which QMOMKit added in
version 1.1.

## Model

The phase-space coordinate $\xi$ is the velocity of the particles. Without forces and collisions, the number density $n(x, \xi, t)$ obeys

```math
  \frac{\partial n}{\partial t} + \xi \frac{\partial n}{\partial x} = 0 ,
```

so every particle keeps its velocity, and $n(x, \xi, t) = n(x - \xi t, \xi, 0)$. The code transports the raw moments $m_k = \int \xi^k n \, d\xi$ for
$k = 0, \ldots, 2N$:

```math
  \frac{\partial m_k}{\partial t} + \frac{\partial m_{k+1}}{\partial x} = 0 .
```

The flux of $m_{2N}$ is $m_{2N+1}$, which the transported moments do not contain. This is the closure problem of the transport, and HyQMOM [Fox &
Laurent, 2022] solves it: from $m_0, \ldots, m_{2N}$, QMOMKit reconstructs $N + 1$ weighted points, and their moment $m_{2N+1}$ closes the system. The
closure makes the system globally hyperbolic. Its $2N + 1$ characteristic speeds are the eigenvalues of two Jacobi matrices, and QMOMKit returns them
together with a weight for each speed and the closure moment.

Nothing changes the velocities, so unlike the other couplings in this repository, this coupling has no phase-space step and no call to
`advance_qmom_phase_space(...)`.

## How the code uses QMOMKit

The case is built once:

```cpp
qmomkit::QmomCase qmom_case;
qmom_case.node_count = N;
qmom_case.family = qmomkit::QmomFamily::hyperbolic;
```

In every step, `main()` does the following for each cell:

1. `qmom_case.state_candidate(qmomkit::MomentSet(moments))` creates a state from the moments of the cell. A cell whose $m_0$ is below
   `QmomOptions::cutoff` is `empty`, has no state and sends no flux. Any other status stops the run with QMOMKit's message.
2. The moments of the state, `state.moment_values()`, replace the moments of the cell, as described below.
3. `compute_hyqmom_characteristic_speeds(state)` returns the speeds in ascending order, their weights and `closure_moment`.
4. `state.quadrature()` holds the abscissas and weights of the weighted points.

The default options of the case are adaptive. When a cell holds only $k \le N$ different velocities, its moments lie on the boundary of the realizable
set, and QMOMKit returns the $k$ weighted points of a classical quadrature instead of $N + 1$.

**The moments and the weighted points form one unit.** The moments of the state, $m_k = \sum_j w_j \xi_j^k$, replace those of the cell, as QMOMKit's
coupling guide prescribes for every host that changes moments itself. The weighted points carry the fluxes, so the transport then changes the moments
of the same weighted points that it moves.

## Fluxes

The face between a left cell L and a right cell R gets one of three fluxes, selected with `--flux`. All three are conservative, because each face flux
leaves one cell and enters the other.

- **`nodes`**, the kinetic flux of Desjardins, Fox & Villedieu (2008, eqs. 33 to 36). The weighted points of L that move right and those of R that
  move left carry their moments through the face:

  ```math
    F_k = \sum_{j \in L} w_j \max(\xi_j, 0)\, \xi_j^{k} + \sum_{j \in R} w_j \min(\xi_j, 0)\, \xi_j^{k} .
  ```

- **`speeds`**, the same kinetic flux with the $2N + 1$ characteristic speeds and their weights in place of the weighted points (Patel, Desjardins &
  Fox, 2019, eqs. 13 and 14).
- **`hll`**, the HLL flux of Fox & Laurent (2022, Section 6), with the physical flux $(m_1, \ldots, m_{2N}, m_{2N+1})$ of each cell and
  `closure_moment` as $m_{2N+1}$. The wave speeds are the slowest and the fastest speed of L and R with `--hll-speeds local`, the default, or of the
  whole domain with `--hll-speeds global`, as in Fox & Laurent.

## Cases

| `--case` | Setup | Domain and boundaries | Cells | Time step | Output times |
|---|---|---|---|---|---|
| `lm-mono` | a block on $0.3 \le x \le 0.5$ with $m_0 = 1$, all particles at velocity 0.6 | $[0, 1]$, periodic | 600 | $10^{-4}$ | 0.05, 0.5, 1, 1.5, 2, 2.5 |
| `lm-bi` | the same block, velocities 0.6 and 0.4 with number fractions 0.7 and 0.3 | as above | 600 | $10^{-4}$ | as above |
| `lm-tri` | the same block, velocities 0.6, 0.4 and 0.2 with number fractions 0.5, 0.3 and 0.2 | as above | 600 | $10^{-4}$ | as above |
| `lm-tri-gauss` | as `lm-tri`, with a normal distribution of standard deviation `--velocity-width` (0.02 by default) around each velocity | as above | 600 | $10^{-4}$ | as above |
| `lm-crossing-a` | a block on $[0.25, 0.35]$ at velocity 0.5 and a block at rest on $[0.45, 0.55]$, each with $m_0 = 1$ | as above | 600 | $10^{-4}$ | 0.1 to 0.6 |
| `lm-crossing-b` | blocks on $[0.25, 0.35]$ at velocity 0.5 and on $[0.65, 0.75]$ at velocity -0.5 | as above | 600 | $10^{-4}$ | 0.1 to 0.8 |
| `lm-crossing-c` | blocks on $[0.2, 0.3]$ at velocity 0.5, at rest on $[0.45, 0.55]$ and on $[0.7, 0.8]$ at velocity -0.5 | as above | 600 | $10^{-4}$ | 0.1 to 1.0 |
| `fl-riemann` | normally distributed velocities with $m_0 = 1$ and variance $1/3$ everywhere, mean velocity 1 for $x < 0$ and -1 for $x \ge 0$ | $[-0.5, 0.5]$, outflow | 2000 | Courant number 0.9 | 0.05, 0.1 |

The `lm-` cases follow Li & Marchisio (2022). Section 3.2 gives no numbers, so the parameters of the crossing cases are read off their Figs. 5 and 6.
The three nodes of their HyQMOM correspond to $N = 2$ here. `fl-riemann` follows Fox & Laurent (2022, Section 6).

`--cells` and `--courant` override the number of cells and the time step of a case. With `--courant`, the time step follows from the Courant number
and the largest characteristic speed in the domain in each step. With the 600 cells and the time step of $10^{-4}$ of the `lm-` cases, the Courant
number of the particles is at most 0.036. A run stops when the Courant number of the fastest characteristic speed exceeds one, because the fluxes then
no longer keep the moments realizable.

## Exact solutions and the reference

The exact cell averages follow from moving the particles of each block along their velocities, $n(x, \xi, t) = n(x - \xi t, \xi, 0)$.

When every block has a single velocity, the code also transports each block on its own with first-order upwinding and the same time step, and adds the
moments of the blocks. This sum is what the flux `nodes` gives when the weighted points of every cell are the velocities of the blocks, which is the
case when $N \ge K$ for $K$ velocities. The difference between the coupled run and the sum is then a check of the code, and for $N < K$ it is the
error of the closure.

## Building

Build and install QMOMKit as described in its installation guide, then:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DQMOMKIT_DIR=/path/to/qmomkit/install
cmake --build build -j
```

`QMOMKIT_DIR` is the QMOMKit install prefix and can also be set in the environment.

## Running and output

`run_all.sh` runs every case of the Results section and makes every plot, which takes about 20 minutes:

```bash
./run_all.sh
```

Single runs, from the directory `build`:

```bash
cd build
./hyqmomTransport1D --case lm-tri --nodes 3 --flux nodes
./hyqmomTransport1D --case fl-riemann --nodes 4 --flux hll --hll-speeds global
```

Each run writes into `output/<case>/`, which `--output` changes:

- `N<N>_<flux>_t<time>.txt` at $t = 0$ and at every output time, with the columns `x`, the computed moments `m0`, ..., the exact cell averages
  `exact_m0`, ..., and, for the discrete-velocity cases, the reference `reference_m0`, ...;
- `N<N>_<flux>_summary.txt`, the errors and other lines that the run also prints.

`<flux>` is `nodes`, `speeds`, `hll-local` or `hll-global`. A run that stops prints the reason, writes it into the summary file and returns 1.

`plot_hyqmom.py` writes the plots as PDF files into `output/figures/`, named after the case. For the two runs above:

```bash
python3 ../plot_hyqmom.py profiles output lm-tri --flux nodes   # m0 over x at every output time
python3 ../plot_hyqmom.py riemann  output --flux hll-global     # Figs. 3 and 4 of Fox & Laurent
```

Two more commands compare several runs: `fluxes` compares the fluxes for one node count, and `convergence` the errors of `fl-riemann` against the node
count. The Results section shows both.

## Results

### `lm-mono`

```bash
for n in 1 2 3 4; do ./hyqmomTransport1D --case lm-mono --nodes $n --flux nodes; done
python3 ../plot_hyqmom.py profiles output lm-mono --flux nodes
```

All particles have the velocity 0.6, so QMOMKit reconstructs a single weighted point in every cell, and the runs for $N = 1$ to $4$ are identical. At
$t = 2.5$, the relative L1 error of $m_0$ is 0.392, and the relative L1 difference to the reference is below $10^{-13}$, so the error is the numerical
diffusion of first-order upwinding.

### `lm-bi`

```bash
for n in 1 2 3; do ./hyqmomTransport1D --case lm-bi --nodes $n --flux nodes; done
python3 ../plot_hyqmom.py profiles output lm-bi --flux nodes
```

At $t = 2.5$, with the flux `nodes`:

| $N$ | relative L1 error of $m_0$ | relative L1 difference of $m_0$ to the reference |
|---|---|---|
| 1 | 1.15 | 0.838 |
| 2 | 0.371 | $< 10^{-7}$ |
| 3 | 0.371 | $< 10^{-7}$ |

Two velocities need $N = 2$. From $N = 2$ on, the weighted points of every cell are the two velocities, the runs agree with the reference, and the
error is the numerical diffusion, in agreement with Fig. 4b of Li & Marchisio. With $N = 1$, the two weighted points cannot represent the two
velocities, and the blocks do not separate.

### `lm-tri`

```bash
for n in 2 3 4; do ./hyqmomTransport1D --case lm-tri --nodes $n --flux nodes; done
python3 ../plot_hyqmom.py profiles output lm-tri --flux nodes
```

At $t = 2.5$, with the flux `nodes`:

| $N$ | relative L1 error of $m_0$ | relative L1 difference of $m_0$ to the reference |
|---|---|---|
| 2 | 0.949 | 0.627 |
| 3 | 0.338 | $< 10^{-7}$ |
| 4 | 0.338 | $< 10^{-7}$ |

Three velocities need $N = 3$. From $N = 3$ on, the runs agree with the reference, and the error is the numerical diffusion. With $N = 2$, as with the
three nodes of Li & Marchisio in their Fig. 4c, the blocks do not separate. $N = 1$ is left out.

### `lm-tri-gauss`

```bash
for n in 2 3 4 5 6; do ./hyqmomTransport1D --case lm-tri-gauss --nodes $n --flux nodes; done
python3 ../plot_hyqmom.py profiles output lm-tri-gauss --flux nodes
```

The velocities of each block of `lm-tri` are normally distributed with a standard deviation of 0.02, so the moments lie inside the realizable set and
the closure of HyQMOM matters. With the flux `nodes`:

| $N$ | relative L1 error of $m_0$ at $t = 1$ | relative L1 error of $m_0$ at $t = 2.5$ |
|---|---|---|
| 2 | 0.155 | 0.597 |
| 3 | 0.082 | 0.258 |
| 4 | 0.074 | 0.162 |
| 5 | 0.068 | 0.116 |
| 6 | 0.068 | 0.115 |

Unlike in `lm-tri`, the error falls with every additional weighted point up to $N = 5$. The groups have no single velocity, so there is no reference.
$N = 1$ is left out.

### `lm-crossing-a`

```bash
for n in 2 3 4; do ./hyqmomTransport1D --case lm-crossing-a --nodes $n --flux nodes; done
python3 ../plot_hyqmom.py profiles output lm-crossing-a --flux nodes
```

A block moving at 0.5 passes through a block at rest and leaves it unchanged. Two velocities need $N = 2$, and the runs for $N = 2$, 3 and 4 are
identical. At $t = 0.6$, the relative L1 error of $m_0$ is 0.176, and the relative L1 difference to the reference is below $10^{-7}$. $N = 1$ is left
out.

### `lm-crossing-b`

```bash
for n in 2 3 4; do ./hyqmomTransport1D --case lm-crossing-b --nodes $n --flux nodes; done
python3 ../plot_hyqmom.py profiles output lm-crossing-b --flux nodes
```

Two blocks moving at 0.5 and -0.5 pass through each other and separate again. Two velocities need $N = 2$, and the runs for $N = 2$, 3 and 4 are
identical. At $t = 0.8$, the relative L1 error of $m_0$ is 0.406, and the relative L1 difference to the reference is below $10^{-7}$. $N = 1$ is left
out.

### `lm-crossing-c`

```bash
./hyqmomTransport1D --case lm-crossing-c --nodes 2 --flux nodes --courant 0.5
for n in 3 4; do ./hyqmomTransport1D --case lm-crossing-c --nodes $n --flux nodes; done
python3 ../plot_hyqmom.py profiles output lm-crossing-c --flux nodes
```

Blocks moving at 0.5 and -0.5 pass through a block at rest. At $t = 1$, with the flux `nodes`:

| $N$ | time step | relative L1 error of $m_0$ | relative L1 difference of $m_0$ to the reference |
|---|---|---|---|
| 2 | Courant number 0.5 | 0.999 | 0.749 |
| 3 | $10^{-4}$ | 0.302 | $< 10^{-5}$ |
| 4 | $10^{-4}$ | 0.302 | $< 10^{-7}$ |

Three velocities need $N = 3$, and with $N = 3$ and $N = 4$ the blocks pass through each other as in the exact solution. With $N = 2$, as with the
three nodes of Li & Marchisio in their Fig. 6c, the closure does not represent the overlapping blocks, and particles accumulate where they meet. The
closure also creates large characteristic speeds in nearly empty cells, so with the time step of $10^{-4}$ the Courant number exceeds one and the run
stops. $N = 2$ therefore runs with `--courant 0.5`. $N = 1$ is left out.

### `fl-riemann`

```bash
for n in 1 2 3 4 5 6 7 8 9 10; do ./hyqmomTransport1D --case fl-riemann --nodes $n --flux hll --hll-speeds global; done
python3 ../plot_hyqmom.py riemann output --flux hll-global
python3 ../plot_hyqmom.py convergence output --flux hll-global
```

This is the case that uses all of HyQMOM: the moments stay inside the realizable set, and the closure moment and the characteristic speeds differ from
those of a classical quadrature. With 2000 cells and the flux `hll-global`:

| $N$ | relative L2 error of $m_0$ at $t = 0.1$ | largest characteristic speed |
|---|---|---|
| 1 | 0.104 | 2.34 |
| 2 | 0.065 | 3.35 |
| 3 | 0.052 | 3.80 |
| 4 | 0.035 | 4.18 |
| 5 | 0.027 | 4.44 |
| 6 | 0.020 | 4.66 |
| 7 | 0.016 | 4.84 |
| 8 | 0.012 | 5.00 |
| 9 | 0.0097 | 5.14 |
| 10 | 0.0080 | 5.28 |

The error falls steadily with $N$, as in Fig. 5 of Fox & Laurent. The largest speed of the run agrees with the values that Fox & Laurent report for
their final solution, 3.34 for $N = 2$ and 5.32 for $N = 10$.

The three fluxes for $N = 4$, with local wave speeds for HLL, together with the run of `hll-global` above:

```bash
for flux in nodes speeds hll; do ./hyqmomTransport1D --case fl-riemann --nodes 4 --flux $flux; done
python3 ../plot_hyqmom.py fluxes output fl-riemann --nodes 4
```

## Code map

| File | Content |
|---|---|
| `hyqmomTransport1D.cc` | `main()`: the QMOMKit case, the state of every cell, the replacement of its moments, the characteristic speeds, the fluxes through the faces and the time step |
| `settings.hpp`, `settings.cc` | the command line |
| `cases.hpp`, `cases.cc` | the cases and the mesh, with the cells on each side of a face |
| `exact_solution.hpp`, `exact_solution.cc` | the exact cell averages of the moments |
| `fluxes.hpp`, `fluxes.cc` | the kinetic and the HLL flux, built from what QMOMKit returns for each cell |
| `reference.hpp`, `reference.cc` | the groups of particles transported separately |
| `output.hpp`, `output.cc` | the profiles and the summary |
| `plot_hyqmom.py` | the plots |
| `run_all.sh` | every case of the Results section, with its plots |

Only `hyqmomTransport1D.cc` calls QMOMKit, and everything from QMOMKit carries the prefix `qmomkit::`. The other files use the namespace
`hyqmom_transport`.

## Known limitations

- The scheme is first order in space and time. With 600 cells and the small Courant number of the `lm-` cases, numerical diffusion spreads each block
  considerably, and the profiles are smoother than in Li & Marchisio, who used the quasi-high-order scheme of Vikas et al. (2011).

## References

- O. Desjardins, R. O. Fox, P. Villedieu, A quadrature-based moment method for dilute fluid-particle flows, J. Comput. Phys. 227 (2008) 2514-2539,
  [doi:10.1016/j.jcp.2007.10.026](https://doi.org/10.1016/j.jcp.2007.10.026).
- R. O. Fox, F. Laurent, Hyperbolic quadrature method of moments for the one-dimensional kinetic equation, SIAM J. Appl. Math. 82 (2022) 750-771,
  [doi:10.1137/21M1406143](https://doi.org/10.1137/21M1406143).
- D. Li, D. Marchisio, Implementation of CHyQMOM in OpenFOAM for the simulation of non-equilibrium gas-particle flows under one-way and two-way
  coupling, Powder Technol. 396 (2022) 765-784, [doi:10.1016/j.powtec.2021.11.004](https://doi.org/10.1016/j.powtec.2021.11.004).
- R. G. Patel, O. Desjardins, R. O. Fox, Three-dimensional conditional hyperbolic quadrature method of moments, J. Comput. Phys. X 1 (2019) 100006,
  [doi:10.1016/j.jcpx.2019.100006](https://doi.org/10.1016/j.jcpx.2019.100006).
- V. Vikas, Z. J. Wang, A. Passalacqua, R. O. Fox, Realizable high-order finite-volume schemes for quadrature-based moment methods, J. Comput. Phys.
  230 (2011) 5328-5352, [doi:10.1016/j.jcp.2011.03.038](https://doi.org/10.1016/j.jcp.2011.03.038).
