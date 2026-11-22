#!/bin/sh
# Runs every case of the README and plots the results, after building as described in the README. The output is written into build/output/<case>/
# and the plots into build/output/figures/.

set -e
cd "$(dirname "$0")/build"
plot() { python3 ../plot_hyqmom.py "$@"; }
rm -rf output

for n in 1 2 3 4; do ./hyqmomTransport1D --case lm-mono --nodes $n --flux nodes; done
plot profiles output lm-mono --flux nodes

for n in 1 2 3; do ./hyqmomTransport1D --case lm-bi --nodes $n --flux nodes; done
plot profiles output lm-bi --flux nodes

for n in 2 3 4; do ./hyqmomTransport1D --case lm-tri --nodes $n --flux nodes; done
plot profiles output lm-tri --flux nodes

for n in 2 3 4 5 6; do ./hyqmomTransport1D --case lm-tri-gauss --nodes $n --flux nodes; done
plot profiles output lm-tri-gauss --flux nodes

for n in 2 3 4; do ./hyqmomTransport1D --case lm-crossing-a --nodes $n --flux nodes; done
plot profiles output lm-crossing-a --flux nodes

for n in 2 3 4; do ./hyqmomTransport1D --case lm-crossing-b --nodes $n --flux nodes; done
plot profiles output lm-crossing-b --flux nodes

./hyqmomTransport1D --case lm-crossing-c --nodes 2 --flux nodes --courant 0.5
for n in 3 4; do ./hyqmomTransport1D --case lm-crossing-c --nodes $n --flux nodes; done
plot profiles output lm-crossing-c --flux nodes

for n in 1 2 3 4 5 6 7 8 9 10; do ./hyqmomTransport1D --case fl-riemann --nodes $n --flux hll --hll-speeds global; done
plot riemann output --flux hll-global
plot convergence output --flux hll-global

for flux in nodes speeds hll; do ./hyqmomTransport1D --case fl-riemann --nodes 4 --flux $flux; done
plot fluxes output fl-riemann --nodes 4
