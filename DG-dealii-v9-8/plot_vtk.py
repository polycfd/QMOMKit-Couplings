#!/usr/bin/env python3
"""
Read VTK files produced by deal.II and plot moments vs x at the final time step.

Usage:
    python plot_vtk.py <directory> <prefix> [--ymin YMIN] [--ymax YMAX]

Examples:
    python plot_vtk.py build solutionB 
    python plot_vtk.py build solutionB --ymin 0 --ymax 1.2
"""

import argparse
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (17,8)
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

def parse_vtk(filepath):
    moments = []
    """Parse a deal.II legacy VTK file and return a data structure."""
    with open(filepath, "r") as f:
        lines = f.readlines()

    # read points
    i = 0
    n_points = 0
    while i < len(lines):
        if lines[i].startswith("POINTS"):
            n_points = int(lines[i].split()[1])
            i += 1
            break
        i += 1

    coords = []
    while len(coords) < n_points:
        parts = lines[i].split()
        # each line is "x y z" — we only need x
        coords.append(float(parts[0]))
        i += 1

    x = np.array(coords)

    # read moments
    while i < len(lines):
        if lines[i].startswith("SCALARS"):
            # skip LOOKUP_TABLE line
            i += 2
            moments.append(np.array([float(v) for v in lines[i].split()]))
        i += 1

    moments = np.array(moments)
    return x, moments


def discover_file(directory, prefix):
    """The last file with the matching prefix-<N>.vtk."""
    # Data extraction
    max_number = -1
    pattern = re.compile(rf"^{re.escape(prefix)}-(\d+)\.vtk$")

    # Find the maximum number of files and the last file
    for filename in os.listdir(directory):
        match = pattern.match(filename)
        if match:
            number = int(match.group(1))
            if number > max_number:
                max_number = number
     
    file = os.path.join(directory, prefix + f"-{max_number}.vtk")

    print(f"Last file found: {file}")
    return file


def main():
    parser = argparse.ArgumentParser(description="Plot moments over physical space.")
    parser.add_argument("directory", help="Directory containing the VTK files")
    parser.add_argument("prefix", help="File prefix (e.g. solutionB)")
    parser.add_argument("--ymin", type=float, default=None, help="Y-axis minimum")
    parser.add_argument("--ymax", type=float, default=None, help="Y-axis maximum")

    args = parser.parse_args()

    # Find file corresponding to last time step
    file = discover_file(args.directory, args.prefix)
    if not file:
        sys.exit(f"No files matching {args.prefix}-*.vtk in {args.directory}")

    # Parse vtk file
    x, moments = parse_vtk(file)
    # Initialize figure
    fig, ax = plt.subplots(2, 5)
    ax = ax.flatten()

    # Plot each moment
    for i, moment in enumerate(moments):
        ax[i].plot(x, moment, "k-")
        ax[i].set_xlabel("x")
        ax[i].set_ylabel(f"Moment {i}")
        if args.ymin is not None or args.ymax is not None:
            ax[i].set_ylim(args.ymin, args.ymax)

    fig.tight_layout()
    fig.savefig(f"{args.prefix}_moments.png")
    plt.show()


if __name__ == "__main__":
    main()
