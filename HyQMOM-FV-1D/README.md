# HyQMOM transport in one dimension

Transport of the velocity moments of a number density in one space dimension,
with the fluxes in physical space built from the HyQMOM of
[QMOMKit](https://github.com/polycfd/QMOMKit), with first-order finite volumes
and forward Euler. The code shows how a host obtains the closure of the flux,
the characteristic speeds and the weighted points from QMOMKit, and how it
builds three different fluxes from them. `hyqmomTransport1D.cc` holds only
`main()`, with the time loop and every call to QMOMKit, and the other files
hold the parts that do not involve QMOMKit. The cases are those of Li & Marchisio (2022, Sections 3.1 and 3.2), a
variant of them with normally distributed velocities, and the Riemann problem
of Fox & Laurent (2022, Section 6). Every case has an exact solution.

The coupling needs the characteristic speed weights and the closure moment of
`compute_hyqmom_characteristic_speeds(...)`, which QMOMKit added after
version 1.0.0.

## Model

The phase-space coordinate $\xi$ is the velocity of the particles. Without
forces and collisions, the number density $n(x, \xi, t)$ obeys

```math
  \frac{\partial n}{\partial t} + \xi \frac{\partial n}{\partial x} = 0 ,
```

so every particle keeps its velocity, and $n(x, \xi, t) = n(x - \xi t, \xi, 0)$.
The code transports the raw moments $m_k = \int \xi^k n \, d\xi$ for
$k = 0, \ldots, 2N$:

```math
  \frac{\partial m_k}{\partial t} + \frac{\partial m_{k+1}}{\partial x} = 0 .
```

The flux of $m_{2N}$ is $m_{2N+1}$, which the transported moments do not
contain. This is the closure problem of the transport, and HyQMOM
[Fox & Laurent, 2022] solves it: from $m_0, \ldots, m_{2N}$, QMOMKit
reconstructs $N + 1$ weighted points, and their moment $m_{2N+1}$ closes the
system. The closure makes the system globally hyperbolic. Its $2N + 1$
characteristic speeds are the eigenvalues of two Jacobi matrices, and QMOMKit
returns them together with a weight for each speed and the closure moment.

This coupling differs from the other couplings in this repository. There, the
phase-space coordinate is an internal coordinate, such as the size of the
particles, which a physical velocity carries through space, and QMOMKit
advances the phase space in a separate step. Here, the phase-space coordinate
is the velocity itself, so the moments carry each other through space. Nothing
changes the velocities, so there is no operator splitting and no call to
`advance_qmom_phase_space(...)`.

## How the code uses QMOMKit

The case is built once:

```cpp
qmomkit::QmomCase qmom_case;
qmom_case.node_count = N;
qmom_case.family = qmomkit::QmomFamily::hyperbolic;
```

In every step, `main()` does the following for each cell:

1. `qmom_case.state_candidate(qmomkit::MomentSet(moments))` creates a state
   from the moments of the cell. A cell whose $m_0$ is below
   `QmomOptions::cutoff` is `empty`, has no state and sends no flux. Any other
   status stops the run with QMOMKit's message.
2. The moments of the state, `state.moment_values()`, replace the moments of
   the cell, as described below.
3. `compute_hyqmom_characteristic_speeds(state)` returns the speeds in
   ascending order, their weights and `closure_moment`.
4. `state.quadrature()` holds the abscissas and weights of the weighted points.

The default options of the case are adaptive. When a cell holds only $k \le N$
different velocities, its moments lie on the boundary of the realizable set,
and QMOMKit returns the $k$ weighted points of a classical quadrature instead of
$N + 1$. The discrete-velocity cases below consist almost entirely of such
cells.

**The moments and the weighted points form one unit.** Once QMOMKit has
reconstructed the weighted points of a cell, the moments of the cell are the
moments of these weighted points, $m_k = \sum_j w_j \xi_j^k$, which the state
stores. QMOMKit's coupling guide prescribes this for every host that changes
moments itself, and the other couplings in this repository do the same after
their phase-space step. The weighted points also carry the fluxes, so the
transport changes the moments of the same weighted points that it moves. For a
cell inside the realizable set, the $N + 1$ weighted points reproduce all
moments, and the replacement changes them only by rounding. For a cell on the
boundary, the $k$ weighted points reproduce $m_0, \ldots, m_{2k-1}$ and
determine $m_{2k}, \ldots, m_{2N}$. Moments kept apart from the weighted points
would there accumulate rounding errors that no flux removes, and in the tail
that first-order upwinding leaves behind the particles, where the content of a
cell decays exponentially, these errors would soon dominate the moments. The
adaptive reconstruction drops weighted points of negligible weight, so the
summed $m_0$ is conserved to about $10^{-12}$ rather than to rounding.

## Fluxes

The face between a left cell L and a right cell R gets one of three fluxes,
selected with `--flux`:

- **`nodes`**, the kinetic flux of Desjardins, Fox & Villedieu (2008,
  eqs. 33 to 36). The weighted points of L that move right and those of R that
  move left carry their moments through the face:

  ```math
    F_k = \sum_{j \in L} w_j \max(\xi_j, 0)\, \xi_j^{k} + \sum_{j \in R} w_j \min(\xi_j, 0)\, \xi_j^{k} .
  ```

  With a Courant number below one for every weighted point, the updated
  moments are those of a number density that is nowhere negative (Desjardins
  et al., 2008, Proposition 3.1).
- **`speeds`**, the same kinetic flux with the $2N + 1$ characteristic speeds
  and their weights in place of the weighted points (Patel, Desjardins &
  Fox, 2019, eqs. 13 and 14). With their weights, the speeds have the moments
  $m_0, \ldots, m_{2N+1}$ of the weighted points, and inside the realizable set
  the weights are positive. For a classical quadrature, the two copies of each
  velocity share its weight and the mean has weight zero, so this flux equals
  the flux `nodes`.
- **`hll`**, the HLL flux used by Fox & Laurent (2022, Section 6). The
  physical flux of a cell is $(m_1, \ldots, m_{2N}, m_{2N+1})$, with the
  moments of the weighted points and `closure_moment` as $m_{2N+1}$. With
  `--hll-speeds local`, the default, the wave speeds are the slowest and the
  fastest speed of L and R. With `--hll-speeds global`, they are the slowest
  and the fastest speed over the whole domain in each step, as in Fox &
  Laurent. An empty cell has no speeds, so only the other cell sets the local
  wave speeds, and the flux between two empty cells is zero. Fox & Laurent did
  not prove that this scheme keeps the moments realizable. They checked it in
  every step, and this code stops when a state cannot be created.

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

The `lm-` cases follow Li & Marchisio (2022). Section 3.1 gives the block, the
velocities, the number fractions, the 600 cells and the time step of
$10^{-4}$ s. The domain $[0, 1]$ and the output times come from their Fig. 4.
Section 3.2 gives no numbers, so the blocks, velocities and times of the
crossing cases are read off their Figs. 5 and 6. Li & Marchisio use the
three-node HyQMOM of Fox, Laurent & Vié (2018), with its middle node fixed at
the mean velocity. QMOMKit implements the closure of Fox & Laurent (2022),
which has no such restriction and no limit on the node count. In the notation
of this code, the three nodes of Li & Marchisio correspond to $N = 2$, so the
comparison with their figures is qualitative.

`fl-riemann` follows Fox & Laurent (2022, Section 6). They used 16000 cells and
the HLL flux with global wave speeds. The boundaries copy the cell next to them
into the ghost cell, which is exact here because the states outside the domain
stay uniform until $t = 0.1$.

`--cells` and `--courant` override the number of cells and the time step of a
case. With `--courant`, the time step follows from the Courant number and the
largest characteristic speed in the domain in each step.

## Exact solutions and the reference

Every particle keeps its velocity, so the particles of a block with a single
velocity $u$ have moved by $u t$ at time $t$. The exact cell averages of such a
block follow from the overlap of the moved block with each cell, and in a
periodic domain the copies of the block shifted by whole domain lengths add
their particles. For a block with normally distributed velocities, the
particles at $x$ are those whose velocities lie in
$[(x - x_\mathrm{right}) / t, (x - x_\mathrm{left}) / t]$, so the moments at
$x$ are the moments of the normal distribution over this interval. They follow
from a recursion for the moments of a truncated normal distribution, and a
composite Gauss-Legendre rule averages them over each cell. The recursion
agrees with direct numerical integration to within $5 \times 10^{-14}$ of the
moment scale up to order 21.

When every block has a single velocity, the code also transports each block on
its own with first-order upwinding and the same time step, and adds the moments
of the blocks. This sum is what the flux `nodes` gives when the weighted points
of every cell are the velocities of the blocks, which is the case when
$N \ge K$ for $K$ velocities. The difference between the coupled run and the
sum is then a check of the code, and for $N < K$ it is the error of the
closure.

## Building

Build and install QMOMKit as described in its installation guide, then:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DQMOMKIT_DIR=/path/to/qmomkit/install
cmake --build build -j
```

`QMOMKIT_DIR` is the QMOMKit install prefix and can also be set in the
environment.

## Running and output

```bash
cd build
./hyqmomTransport1D --case lm-tri --nodes 3 --flux nodes
./hyqmomTransport1D --case fl-riemann --nodes 4 --flux hll --hll-speeds global
```

Each run writes into `output/<case>/`, which `--output` changes:

- `N<N>_<flux>_t<time>.txt` at $t = 0$ and at every output time, with the
  columns `x`, the computed moments `m0`, ..., the exact cell averages
  `exact_m0`, ..., and, for the discrete-velocity cases, the sum of the
  separately transported blocks `reference_m0`, ...;
- `N<N>_<flux>_summary.txt`, the lines the run also prints. For every output
  time, they give the relative L1 error of $m_0$ and $m_1$, the relative L1
  difference to the reference, the number of occupied cells whose quadrature
  holds fewer than $N + 1$ weighted points, and the relative change of the
  summed $m_0$. At the final time, they give the relative L2 error of every
  moment, normalized as in Fig. 5 of Fox & Laurent, the number of steps, the
  largest characteristic speed and the largest Courant number computed from
  the speeds.

`<flux>` is `nodes`, `speeds`, `hll-local` or `hll-global`. A run that stops
prints QMOMKit's message, with the cell, its position and the time, and returns
1.

`plot_hyqmom.py` plots the output:

```bash
python plot_hyqmom.py profiles    output lm-tri --flux nodes   # m0 over x at every output time, all N
python plot_hyqmom.py fluxes      output lm-crossing-c --nodes 3 # the fluxes at the final time
python plot_hyqmom.py riemann     output --flux hll-global      # Figs. 3 and 4 of Fox & Laurent
python plot_hyqmom.py convergence output --flux hll-global      # Fig. 5 of Fox & Laurent
```

## Results

### `lm-mono`

```bash
for n in 1 2 3 4; do ./hyqmomTransport1D --case lm-mono --nodes $n --flux nodes; done
python ../plot_hyqmom.py profiles output lm-mono --flux nodes
```

All particles have the velocity 0.6, so every occupied cell holds a single
velocity, and QMOMKit reconstructs a single weighted point for every $N$. The
runs for $N = 1$ to $4$ are identical. At $t = 2.5$, the relative L1 error of
$m_0$ is 0.392, the relative L1 difference to the reference is
$7 \times 10^{-14}$, the largest characteristic speed is 0.6, and the summed
$m_0$ changes by $4 \times 10^{-15}$. The error is entirely the numerical
diffusion of first-order upwinding.

### `lm-bi`

```bash
for n in 1 2 3; do ./hyqmomTransport1D --case lm-bi --nodes $n --flux nodes; done
python ../plot_hyqmom.py profiles output lm-bi --flux nodes
```

At $t = 2.5$, with the flux `nodes`:

| $N$ | relative L1 error of $m_0$ | relative L1 difference of $m_0$ to the reference | occupied cells with fewer than $N + 1$ weighted points | largest characteristic speed |
|---|---|---|---|---|
| 1 | 1.15 | 0.838 | 0 of 600 | 2.17 |
| 2 | 0.371 | $1.6 \times 10^{-8}$ | 600 of 600 | 0.600 |
| 3 | 0.371 | $1.6 \times 10^{-8}$ | 600 of 600 | 0.600 |

Two velocities need $N = 2$. From $N = 2$ on, the weighted points of every
cell are the two velocities of the particles, the runs are identical and agree
with the reference to $1.6 \times 10^{-8}$, and the error of 0.371 is the
numerical diffusion of first-order upwinding. $N = 4$ gives the same result.
This agrees with Li & Marchisio, whose three weighted points, $N = 2$ here,
separate the two blocks in their Fig. 4b. With $N = 1$, the two weighted points
lie symmetrically around the mean velocity and carry equal weights, so they
cannot represent number fractions of 0.7 and 0.3, and the blocks do not
separate. The summed $m_0$ changes by at most $4.6 \times 10^{-12}$.

### `lm-tri`

```bash
for n in 2 3 4; do ./hyqmomTransport1D --case lm-tri --nodes $n --flux nodes; done
python ../plot_hyqmom.py profiles output lm-tri --flux nodes
```

At $t = 2.5$, with the flux `nodes`:

| $N$ | relative L1 error of $m_0$ | relative L1 difference of $m_0$ to the reference | occupied cells with fewer than $N + 1$ weighted points | largest characteristic speed |
|---|---|---|---|---|
| 2 | 0.949 | 0.627 | 0 of 600 | 1.97 |
| 3 | 0.338 | $5.6 \times 10^{-8}$ | 541 of 600 | 0.698 |
| 4 | 0.338 | $5.5 \times 10^{-8}$ | 600 of 600 | 0.641 |

For $N = 3$ and $N = 4$, the weighted points of every cell are the three
velocities of the particles, and the runs agree with the reference, the three
blocks transported separately, to $6 \times 10^{-8}$. The difference comes from
the adaptive reconstruction, which drops weighted points of negligible weight.
The whole error of 0.338 is therefore the numerical diffusion of first-order
upwinding: the three blocks separate as in the exact solution, but their edges
are smeared over many cells. $N = 2$ corresponds to the three weighted points
of Li & Marchisio, and as in their Fig. 4c, the blocks do not separate, because
three weighted points reconstructed from $m_0, \ldots, m_4$ cannot recover
three velocities. The summed $m_0$ changes by at most $2.4 \times 10^{-12}$.
$N = 1$ is left out, because two weighted points cannot represent three
velocities.

### `lm-tri-gauss`

```bash
for n in 2 3 4 5 6; do ./hyqmomTransport1D --case lm-tri-gauss --nodes $n --flux nodes; done
python ../plot_hyqmom.py profiles output lm-tri-gauss --flux nodes
```

The velocities of each block of `lm-tri` are normally distributed with a
standard deviation of 0.02 around 0.6, 0.4 and 0.2. Each block therefore
spreads as it moves, and no cell holds a finite number of velocities, so the
moments lie inside the realizable set and the closure moment of HyQMOM matters
in every cell. With the flux `nodes`:

| $N$ | relative L1 error of $m_0$ at $t = 1$ | relative L1 error of $m_0$ at $t = 2.5$ | occupied cells with fewer than $N + 1$ weighted points at $t = 2.5$ | largest characteristic speed |
|---|---|---|---|---|
| 2 | 0.155 | 0.597 | 0 of 600 | 2.01 |
| 3 | 0.082 | 0.258 | 0 of 600 | 1.57 |
| 4 | 0.074 | 0.162 | 0 of 600 | 1.69 |
| 5 | 0.068 | 0.116 | 0 of 600 | 1.37 |
| 6 | 0.068 | 0.115 | 36 of 600 | 1.23 |

Unlike in `lm-tri`, where the error drops in one step once $N$ reaches the
number of velocities, the error here falls with every additional weighted
point up to $N = 5$. At $N = 5$ and $N = 6$, the error no longer changes. The
blocks are then probably resolved as well as the numerical diffusion of the
600 cells allows, which has not been checked on a finer mesh. The groups have
no single velocity, so the code computes no reference for this case. The
summed $m_0$ changes by at most $1.1 \times 10^{-11}$. $N = 1$ is left out.

### `lm-crossing-a`

```bash
for n in 2 3 4; do ./hyqmomTransport1D --case lm-crossing-a --nodes $n --flux nodes; done
python ../plot_hyqmom.py profiles output lm-crossing-a --flux nodes
```

A block moving at 0.5 passes through a block at rest. Two velocities need
$N = 2$, and the runs for $N = 2$, 3 and 4 are identical. At $t = 0.6$, the
relative L1 error of $m_0$ is 0.176, the relative L1 difference to the
reference is $8.9 \times 10^{-9}$, every occupied cell holds two or fewer
weighted points, the largest characteristic speed is 0.5, and the summed $m_0$
changes by at most $1.2 \times 10^{-15}$. The moving block passes through the
block at rest and leaves it unchanged, as in the exact solution. The block at
rest keeps its sharp edges, because particles at rest cross no faces. The peaks
and dips where the two blocks meet come from the smeared edges of the moving
block next to the sharp edges of the block at rest, and the reference shows
them as well. $N = 1$ is left out, because two velocities need $N = 2$.

### `lm-crossing-b`

```bash
for n in 2 3 4; do ./hyqmomTransport1D --case lm-crossing-b --nodes $n --flux nodes; done
python ../plot_hyqmom.py profiles output lm-crossing-b --flux nodes
```

Two blocks with the velocities 0.5 and -0.5 pass through each other. Two
velocities need $N = 2$, and the runs for $N = 2$, 3 and 4 are identical. At
$t = 0.8$, the relative L1 error of $m_0$ is 0.406, the relative L1 difference
to the reference is $7.9 \times 10^{-9}$, the largest characteristic speed is
0.500005, and the summed $m_0$ changes by $1.8 \times 10^{-14}$. The blocks
cross and separate again with their original shapes, apart from the numerical
diffusion, and no particles accumulate where they overlap. With $N = 2$, the
fluxes `speeds` and `hll` agree with the flux `nodes` to a relative L1
difference of $10^{-14}$ and $1.2 \times 10^{-11}$: every particle moves at
0.5 or -0.5, and with wave speeds of -0.5 and 0.5, the HLL flux reduces to the
upwind flux of each block. $N = 1$ is left out, because two velocities need
$N = 2$.

### `lm-crossing-c`

```bash
./hyqmomTransport1D --case lm-crossing-c --nodes 2 --flux nodes --courant 0.5
for n in 3 4; do ./hyqmomTransport1D --case lm-crossing-c --nodes $n --flux nodes; done
python ../plot_hyqmom.py profiles output lm-crossing-c --flux nodes
```

Blocks with the velocities 0.5 and -0.5 pass through a block at rest. At
$t = 1$, with the flux `nodes`:

| $N$ | time step | relative L1 error of $m_0$ | relative L1 difference of $m_0$ to the reference | occupied cells with fewer than $N + 1$ weighted points | largest characteristic speed |
|---|---|---|---|---|---|
| 2 | Courant number 0.5 | 0.999 | 0.749 | 0 of 600 | 127 |
| 3 | $10^{-4}$ | 0.302 | $1.6 \times 10^{-6}$ | 500 of 600 | 1.16 |
| 4 | $10^{-4}$ | 0.302 | $5.4 \times 10^{-8}$ | 558 of 600 | 0.906 |

Three velocities need $N = 3$. With $N = 3$ and $N = 4$, the three blocks pass
through each other and separate with their original shapes, apart from the
numerical diffusion, and the runs agree with the reference. With $N = 3$, 100
occupied cells get four weighted points instead of at most three, and these
cells are the likely source of the larger difference to the reference and of
the largest speed of 1.16. In exact arithmetic, the moments of a cell with
three velocities lie on the boundary of the realizable set, so rounding errors
that place them just inside are a probable cause, which has not been checked.

$N = 2$ corresponds to the three weighted points of Li & Marchisio, and as in
their Fig. 6c, the closure does not represent three blocks that overlap. Where
the blocks meet, the three weighted points reconstructed from
$m_0, \ldots, m_4$ differ from the three velocities, particles accumulate at
$x = 0.5$, with a peak of $m_0 = 5.9$ at $t = 0.5$ where the exact value is 3,
and the outer blocks arrive deformed. The closure also creates characteristic
speeds that grow in nearly empty cells, here up to 127. With the time step of
$10^{-4}$ of the other runs, the Courant number of these speeds exceeds one,
$m_0$ becomes negative, and the run stops at $t = 0.58$. With
`--courant 0.5`, the time step follows the largest speed, and the run
completes in 7023 steps. Its time step is $1.4 \times 10^{-4}$ on average, so
its numerical diffusion differs from that of $N = 3$ and $N = 4$. The summed
$m_0$ changes by at most $6 \times 10^{-13}$. $N = 1$ is left out.

### `fl-riemann`

```bash
for n in 1 2 3 4 5 6 7 8 9 10; do ./hyqmomTransport1D --case fl-riemann --nodes $n --flux hll --hll-speeds global; done
python ../plot_hyqmom.py riemann output --flux hll-global
python ../plot_hyqmom.py convergence output --flux hll-global
```

This is the case that uses all of HyQMOM: the moments stay inside the
realizable set, the closure moment and the characteristic speeds differ from
those of a classical quadrature, and the three fluxes differ. The relative L2
error of $m_0$ at $t = 0.1$ on 2000 cells, and the largest characteristic speed
of the runs with `hll-global`:

| $N$ | `nodes` | `speeds` | `hll-local` | `hll-global` | largest speed |
|---|---|---|---|---|---|
| 1 | 0.104 | 0.110 | 0.104 | 0.104 | 2.34 |
| 2 | 0.076 | 0.075 | 0.070 | 0.065 | 3.35 |
| 3 | 0.066 | 0.064 | 0.058 | 0.052 | 3.80 |
| 4 | 0.050 | 0.049 | 0.039 | 0.035 | 4.18 |
| 5 | 0.044 | 0.042 | 0.032 | 0.027 | 4.44 |
| 6 | 0.037 | 0.036 | 0.024 | 0.020 | 4.66 |
| 7 | 0.032 | 0.031 | 0.019 | 0.016 | 4.84 |
| 8 | 0.028 | 0.028 | 0.015 | 0.012 | 5.00 |
| 9 | 0.025 | 0.024 | 0.012 | 0.0097 | 5.14 |
| 10 | 0.023 | 0.022 | 0.0095 | 0.0080 | 5.28 |

The error of every moment falls steadily with $N$, as in Fig. 5 of Fox &
Laurent, with the odd moments above the even ones. The profiles show the
plateaus and steps of Figs. 3 and 4 of Fox & Laurent, which come from the
finite number of characteristic speeds. The largest speed agrees with the
values that Fox & Laurent report for their final solution, 3.34 for $N = 2$ and
5.32 for $N = 10$.

On finer meshes, the error grows: with `hll-global`, from 0.065 on 2000 cells
to 0.081 on 16000 cells for $N = 2$, from 0.035 to 0.053 for $N = 4$, and from
0.0080 to 0.013 on 4000 cells for $N = 10$. The solution of the closed system
has steps that the exact solution does not have, and the numerical diffusion of
a coarse mesh smooths them, which probably hides part of the closure error. The
errors on 2000 cells are therefore lower than the error of the closure itself,
and a flux with a smaller error on this mesh is not necessarily more accurate.

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

Only `hyqmomTransport1D.cc` calls QMOMKit, and everything from QMOMKit carries
the prefix `qmomkit::`. The other files use the namespace `hyqmom_transport`.

## Known limitations

- The scheme is first order in space and time. With 600 cells and the small
  Courant number of the `lm-` cases, numerical diffusion spreads each block
  considerably, and the profiles are smoother than in Li & Marchisio, who used
  the quasi-high-order scheme of Vikas et al. (2011).
- The HLL flux is first order only. A higher-order HLL flux that reconstructs
  the moments would not keep them realizable.
- The parameters of the crossing cases are read off the figures of Li &
  Marchisio.
- The resolution study of `fl-riemann` stops at 16000 cells for $N \le 4$ and
  at 4000 cells for $N = 10$, because a step of $N = 10$ on 16000 cells is
  expensive.

## References

- O. Desjardins, R. O. Fox, P. Villedieu, A quadrature-based moment method for
  dilute fluid-particle flows, J. Comput. Phys. 227 (2008) 2514-2539,
  [doi:10.1016/j.jcp.2007.10.026](https://doi.org/10.1016/j.jcp.2007.10.026).
- R. O. Fox, F. Laurent, Hyperbolic quadrature method of moments for the
  one-dimensional kinetic equation, SIAM J. Appl. Math. 82 (2022) 750-771,
  [doi:10.1137/21M1406143](https://doi.org/10.1137/21M1406143).
- R. O. Fox, F. Laurent, A. Vié, Conditional hyperbolic quadrature method of
  moments for kinetic equations, J. Comput. Phys. 365 (2018) 269-293.
- D. Li, D. Marchisio, Implementation of CHyQMOM in OpenFOAM for the simulation
  of non-equilibrium gas-particle flows under one-way and two-way coupling,
  Powder Technol. 396 (2022) 765-784,
  [doi:10.1016/j.powtec.2021.11.004](https://doi.org/10.1016/j.powtec.2021.11.004).
- R. G. Patel, O. Desjardins, R. O. Fox, Three-dimensional conditional
  hyperbolic quadrature method of moments, J. Comput. Phys. X 1 (2019) 100006,
  [doi:10.1016/j.jcpx.2019.100006](https://doi.org/10.1016/j.jcpx.2019.100006).
- V. Vikas, Z. J. Wang, A. Passalacqua, R. O. Fox, Realizable high-order
  finite-volume schemes for quadrature-based moment methods, J. Comput. Phys.
  230 (2011) 5328-5352,
  [doi:10.1016/j.jcp.2011.03.038](https://doi.org/10.1016/j.jcp.2011.03.038).
