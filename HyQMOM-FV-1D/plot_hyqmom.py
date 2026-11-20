"""Plots of the output of hyqmomTransport1D.

Usage:
    python plot_hyqmom.py profiles    <output> <case> [--flux nodes] [--nodes 2 3 4]
    python plot_hyqmom.py fluxes      <output> <case> --nodes N
    python plot_hyqmom.py riemann     <output> [--flux hll-global] [--nodes 2 3 4]
    python plot_hyqmom.py convergence <output> [--flux hll-global]

`profiles` plots m0 over x at every output time for all node counts, or for
those given with --nodes, as in
Figs. 4 and 6 of Li & Marchisio (2022). `fluxes` compares the three fluxes at
the final time for one node count. `riemann` plots the moments and the
standardized moments of fl-riemann at t = 0.1, as in Figs. 3 and 4 of Fox &
Laurent (2022), and `convergence` the relative L2 error of every moment
against the node count, as in their Fig. 5. Each command writes a PNG file
into the directory of the case.
"""

import argparse
import glob
import math
import os
import re

import matplotlib.pyplot as plt
import numpy as np

# One colour per node count or flux, always in this order, and the exact
# solution in dark grey.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
EXACT = "#2b2b2b"


def read_profile(path):
    """The columns of one output file, by name."""
    with open(path) as file:
        names = file.readline().lstrip("#").split()
    data = np.loadtxt(path, ndmin=2)
    return {name: data[:, i] for i, name in enumerate(names)}


def profile_files(directory, flux):
    """The output files of one flux, by node count and output time."""
    files = {}
    for path in glob.glob(os.path.join(directory, f"N*_{flux}_t*.txt")):
        match = re.search(r"N(\d+)_.*_t([0-9.]+)\.txt$", os.path.basename(path))
        files.setdefault(int(match.group(1)), {})[float(match.group(2))] = path
    return files


def style(axis):
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(color="#e4e3df", linewidth=0.6)
    axis.set_axisbelow(True)


def plot_profiles(directory, case, flux, selected=None):
    files = profile_files(directory, flux)
    node_counts = sorted(n for n in files if selected is None or n in selected)
    times = sorted({t for n in node_counts for t in files[n] if t > 0.0})
    figure, axes = plt.subplots(len(times), 1, figsize=(7.0, 1.4 * len(times) + 0.6), sharex=True)
    for axis, time in zip(np.atleast_1d(axes), times):
        exact = read_profile(next(files[n][time] for n in node_counts if time in files[n]))
        axis.plot(exact["x"], exact["exact_m0"], color=EXACT, linewidth=1.2, label="exact")
        for slot, n in enumerate(node_counts):
            if time in files[n]:
                profile = read_profile(files[n][time])
                # The colour belongs to the node count. Each line is narrower
                # than the one before and every second line is dashed, so that
                # runs giving the same result stay visible.
                axis.plot(profile["x"], profile["m0"], color=SERIES[(n - 1) % len(SERIES)],
                          linewidth=max(2.6 - 0.8 * slot, 0.8), linestyle="--" if slot % 2 else "-", label=f"N = {n}")
        axis.set_ylabel("$m_0$")
        axis.text(0.01, 0.85, f"t = {time:g}", transform=axis.transAxes)
        style(axis)
    np.atleast_1d(axes)[0].legend(ncol=len(node_counts) + 1, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    np.atleast_1d(axes)[-1].set_xlabel("x")
    figure.suptitle(f"{case}, flux {flux}", y=1.0)
    figure.tight_layout()
    path = os.path.join(directory, f"profiles_{flux}.png")
    figure.savefig(path, dpi=150, bbox_inches="tight")
    print(path)


def plot_fluxes(directory, case, node_count):
    fluxes = ["nodes", "speeds", "hll-local", "hll-global"]
    figure, axis = plt.subplots(figsize=(7.0, 3.0))
    plotted_exact = False
    for slot, flux in enumerate(fluxes):
        files = profile_files(directory, flux).get(node_count)
        if not files:
            continue
        profile = read_profile(files[max(files)])
        if not plotted_exact:
            axis.plot(profile["x"], profile["exact_m0"], color=EXACT, linewidth=1.2, label="exact")
            plotted_exact = True
        axis.plot(profile["x"], profile["m0"], color=SERIES[slot], linewidth=1.5, label=flux)
        final_time = max(files)
    axis.set_xlabel("x")
    axis.set_ylabel("$m_0$")
    axis.legend(frameon=False)
    axis.set_title(f"{case}, N = {node_count}, t = {final_time:g}")
    style(axis)
    figure.tight_layout()
    path = os.path.join(directory, f"fluxes_N{node_count}.png")
    figure.savefig(path, dpi=150, bbox_inches="tight")
    print(path)


def standardized(moments, order):
    """The mean velocity, the variance and the standardized moments S_k up to
    `order`, from the raw moments m_0, ..., m_order."""
    m0 = moments[0]
    mean = moments[1] / m0
    central = []
    for k in range(order + 1):
        central.append(sum(math.comb(k, i) * (-mean) ** (k - i) * moments[i] / m0 for i in range(k + 1)))
    variance = central[2]
    return mean, variance, [central[k] / variance ** (k / 2) for k in range(order + 1)]


def plot_riemann(directory, flux, node_counts):
    files = profile_files(directory, flux)
    final_time = max(files[node_counts[0]])
    largest = max(n for n in files if final_time in files[n])
    exact_profile = read_profile(files[largest][final_time])
    order = min(8, 2 * largest)
    exact_moments = [exact_profile[f"exact_m{k}"] for k in range(order + 1)]
    exact_mean, exact_variance, exact_standardized = standardized(exact_moments, order)

    for kind in ["moments", "standardized"]:
        figure, axes = plt.subplots(3, 3, figsize=(10.0, 7.5))
        for panel, axis in enumerate(axes.flat):
            if kind == "moments":
                label, exact = f"$M_{panel}$", (exact_moments[panel] if panel <= order else None)
            else:
                names = ["$M_0$", "$u$", "$C_2$"] + [f"$S_{k}$" for k in range(3, 9)]
                values = [exact_moments[0], exact_mean, exact_variance] + exact_standardized[3:]
                label, exact = names[panel], (values[panel] if panel < len(values) else None)
            if exact is None:
                axis.set_visible(False)
                continue
            axis.plot(exact_profile["x"], exact, color=EXACT, linewidth=1.2, label="exact")
            for slot, n in enumerate(node_counts):
                if final_time not in files.get(n, {}):
                    continue
                profile = read_profile(files[n][final_time])
                available = 2 * n
                moments = [profile[f"m{k}"] for k in range(min(order, available) + 1)]
                if kind == "moments":
                    if panel > available:
                        continue
                    value = moments[panel]
                else:
                    mean, variance, values_n = standardized(moments, min(order, available))
                    candidates = [moments[0], mean, variance] + values_n[3:]
                    if panel >= len(candidates):
                        continue
                    value = candidates[panel]
                axis.plot(profile["x"], value, color=SERIES[(n - 1) % len(SERIES)], linewidth=1.5, label=f"N = {n}")
            axis.set_ylabel(label)
            axis.set_xlim(-0.4, 0.4)
            style(axis)
        axes.flat[0].legend(frameon=False, fontsize=8)
        figure.suptitle(f"fl-riemann, flux {flux}, t = {final_time:g}")
        figure.tight_layout()
        path = os.path.join(directory, f"riemann_{kind}_{flux}.png")
        figure.savefig(path, dpi=150, bbox_inches="tight")
        print(path)


def plot_convergence(directory, flux):
    errors = {}
    for path in glob.glob(os.path.join(directory, f"N*_{flux}_summary.txt")):
        n = int(re.search(r"N(\d+)_", os.path.basename(path)).group(1))
        with open(path) as file:
            line = [l for l in file if l.startswith("relative L2 error")][-1]
        values = re.findall(r"m(\d+) ([0-9.eE+-]+)", line)
        errors[n] = {int(k): float(v) for k, v in values}
    node_counts = sorted(errors)
    orders = sorted({k for n in node_counts for k in errors[n]})
    # The moment order is a magnitude, so one hue from light to dark encodes
    # it, with even orders solid and odd orders dashed as in Fox & Laurent.
    colours = plt.cm.Blues(np.linspace(0.35, 1.0, len(orders)))
    figure, axis = plt.subplots(figsize=(7.0, 4.5))
    for colour, k in zip(colours, orders):
        ns = [n for n in node_counts if k in errors[n]]
        axis.semilogy(ns, [errors[n][k] for n in ns], color=colour, linewidth=1.2, linestyle="-" if k % 2 == 0 else "--",
                      marker="o", markersize=3, label=f"$m_{{{k}}}$" if k <= 3 or k == orders[-1] else None)
    axis.set_xlabel("N")
    axis.set_ylabel("relative L2 error")
    axis.set_xticks(node_counts)
    axis.legend(frameon=False, fontsize=8, title="darker: higher order", title_fontsize=8)
    axis.set_title(f"fl-riemann, flux {flux}")
    style(axis)
    figure.tight_layout()
    path = os.path.join(directory, f"convergence_{flux}.png")
    figure.savefig(path, dpi=150, bbox_inches="tight")
    print(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["profiles", "fluxes", "riemann", "convergence"])
    parser.add_argument("output", help="the output directory of hyqmomTransport1D")
    parser.add_argument("case", nargs="?", default="fl-riemann")
    parser.add_argument("--flux", default=None)
    parser.add_argument("--nodes", type=int, nargs="+", default=None)
    arguments = parser.parse_args()
    directory = os.path.join(arguments.output, arguments.case)
    if arguments.command == "profiles":
        plot_profiles(directory, arguments.case, arguments.flux or "nodes", arguments.nodes)
    elif arguments.command == "fluxes":
        plot_fluxes(directory, arguments.case, (arguments.nodes or [3])[0])
    elif arguments.command == "riemann":
        plot_riemann(directory, arguments.flux or "hll-global", arguments.nodes or [2, 3, 4])
    else:
        plot_convergence(directory, arguments.flux or "hll-global")


if __name__ == "__main__":
    main()
