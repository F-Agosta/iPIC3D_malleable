#!/usr/bin/env python3
"""
Read and reassemble iPIC3D shdf5 output (proc*.hdf / restart_x*.hdf + settings.hdf)
and plot or compare grid fields.

Each proc file stores its subdomain of a quantity under
    /<group>/<comp>/cycle_<N>        e.g. /fields/Ez/cycle_20
plus /topology/cartesian_coord with the proc's coordinates in the domain
decomposition. Reading ALL proc files and placing each subdomain at
coord*(cells_per_proc) reassembles the full global grid, independent of how
many processes the run used -- this is what makes runs with different MPI
topologies (e.g. after a repartitioning restart) directly comparable.

Usage:
  plot a field over a cycle range (default: all cycles present):
    read_hdf.py --directory runs/weibel/save/baseline-2 --group fields --comp Ez
  compare the same field at the same cycle between two runs:
    read_hdf.py --directory runs/weibel/save/restart-same-2 \
                --compare runs/weibel/save/restart-2-to-4 \
                --cycle 40 --group fields --comp Ez --tol 1e-8 \
                --out-prefix verify_2to4_
  inspect a restart checkpoint directory instead of a save directory:
    read_hdf.py --directory runs/weibel/restart/baseline-2 --glob 'restart_*.hdf' ...

Compare mode exits with code 1 if max|A-B|/max|A| exceeds --tol.
"""

import argparse
import glob
import os
import sys

import numpy as np
import tables
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


# ---------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------

def read_settings(directory):
    setting_filename = os.path.join(directory, "settings.hdf")
    with tables.open_file(setting_filename, mode="r", title="Setting_file") as h5:
        s = {
            "Nprocs": int(h5.root.topology.Nprocs.read()[0]),
            "XLEN": int(h5.root.topology.XLEN.read()[0]),
            "YLEN": int(h5.root.topology.YLEN.read()[0]),
            "ZLEN": int(h5.root.topology.ZLEN.read()[0]),
            "dx": float(h5.root.collective.Dx.read()[0]),
            "dy": float(h5.root.collective.Dy.read()[0]),
            "dz": float(h5.root.collective.Dz.read()[0]),
            "Nxc": int(h5.root.collective.Nxc.read()[0]),
            "Nyc": int(h5.root.collective.Nyc.read()[0]),
            "Nzc": int(h5.root.collective.Nzc.read()[0]),
            "Ns": int(h5.root.collective.Ns.read()[0]),
        }
    return s


def proc_files(directory, proc_glob):
    files = sorted(glob.glob(os.path.join(directory, proc_glob)))
    if not files:
        raise FileNotFoundError(
            f"No files matching '{proc_glob}' in '{directory}'"
        )
    return files


def cycle_of_node(node):
    # dataset names look like "cycle_20"
    return int(node.name.split("_")[-1])


def available_cycles(directory, group, comp, proc_glob):
    files = proc_files(directory, proc_glob)
    cycles = set()
    with tables.open_file(files[0], mode="r") as h5:
        for node in h5.walk_nodes(f"/{group}/{comp}"):
            if isinstance(node, tables.Array):
                cycles.add(cycle_of_node(node))
    return sorted(cycles)


def read_field(directory, group, comp, cycle, proc_glob):
    """
    Reassemble the full global grid of /<group>/<comp>/cycle_<cycle>
    from all proc files in `directory`. Returns a float32 array of
    shape (Nxc+1, Nyc+1, Nzc+1).
    """
    s = read_settings(directory)
    shape = (s["Nxc"] // s["XLEN"] + 1,
             s["Nyc"] // s["YLEN"] + 1,
             s["Nzc"] // s["ZLEN"] + 1)
    field = np.zeros((s["Nxc"] + 1, s["Nyc"] + 1, s["Nzc"] + 1), dtype="float32")

    found = False
    for proc_filename in proc_files(directory, proc_glob):
        with tables.open_file(proc_filename, mode="r", title="Proc_file") as h5:
            coord = h5.root.topology.cartesian_coord.read()
            x0 = int(coord[0]) * (shape[0] - 1)
            y0 = int(coord[1]) * (shape[1] - 1)
            z0 = int(coord[2]) * (shape[2] - 1)
            for node in h5.walk_nodes(f"/{group}/{comp}"):
                if isinstance(node, tables.Array) and cycle_of_node(node) == cycle:
                    field[x0:x0 + shape[0], y0:y0 + shape[1], z0:z0 + shape[2]] = node.read()
                    found = True
    if not found:
        avail = available_cycles(directory, group, comp, proc_glob)
        raise ValueError(
            f"Cycle {cycle} not found for /{group}/{comp} in '{directory}'. "
            f"Available: {avail}"
        )
    return field


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def plot_field(field, comp, cycle, zindex, out_prefix, vmin, vmax):
    slice2d = np.transpose(field[:, :, zindex])  # x on the horizontal axis
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cc = ax.imshow(slice2d, cmap="seismic", vmin=vmin, vmax=vmax,
                   interpolation="antialiased", origin="lower")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    cbar = fig.colorbar(cc)
    cbar.set_label(comp)
    out = f"{out_prefix}{comp}_{str(cycle).zfill(6)}.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved {out}")


def plot_range(directory, args):
    cycles = available_cycles(directory, args.group, args.comp, args.glob)
    if args.first_cycle is not None or args.last_cycle is not None:
        lo = args.first_cycle if args.first_cycle is not None else cycles[0]
        hi = args.last_cycle if args.last_cycle is not None else cycles[-1]
        cycles = [c for c in cycles if lo <= c <= hi]
    step = args.step or 1
    cycles = cycles[::step]
    if not cycles:
        print(f"No cycles to plot for /{args.group}/{args.comp} in {directory}")
        return
    for cycle in cycles:
        field = read_field(directory, args.group, args.comp, cycle, args.glob)
        plot_field(field, args.comp, cycle, args.zindex, args.out_prefix,
                   args.vmin, args.vmax)


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------

def compare_directories(dir_a, dir_b, args):
    fa = read_field(dir_a, args.group, args.comp, args.cycle, args.glob)
    fb = read_field(dir_b, args.group, args.comp, args.cycle, args.glob)
    if fa.shape != fb.shape:
        print(f"FAIL: grid shapes differ: {dir_a} {fa.shape} vs {dir_b} {fb.shape}",
              file=sys.stderr)
        return 1
    diff = np.abs(fa - fb)
    max_abs = float(diff.max())
    denom = float(np.abs(fa).max())
    rel = max_abs / denom if denom > 0 else max_abs
    ok = rel <= args.tol
    print(f"Compare /{args.group}/{args.comp} at cycle {args.cycle}")
    print(f"  A (reference): {dir_a}  shape={fa.shape}")
    print(f"  B (test)     : {dir_b}  shape={fb.shape}")
    print(f"  max|A-B|        = {max_abs:.3e}")
    print(f"  max|A-B|/max|A| = {rel:.3e}   (tol={args.tol:g})")
    print(f"  {'OK' if ok else 'FAIL'}")

    # side-by-side figure: reference, test, difference
    z = args.zindex
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, data, title in (
        (axes[0], fa[:, :, z].T, f"A {args.comp}"),
        (axes[1], fb[:, :, z].T, f"B {args.comp}"),
        (axes[2], diff[:, :, z].T, "|A-B|"),
    ):
        cc = ax.imshow(data, cmap="seismic", origin="lower",
                       interpolation="antialiased")
        ax.set_title(title)
        fig.colorbar(cc, ax=ax)
    out = f"{args.out_prefix}compare_{args.comp}_{str(args.cycle).zfill(6)}.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved {out}")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--directory", default="data/",
                   help="directory with settings.hdf and proc/restart hdf files")
    p.add_argument("--compare", metavar="DIR_B",
                   help="second run directory to compare against --directory")
    p.add_argument("--glob", default="proc*.hdf",
                   help="proc file pattern (use 'restart_*.hdf' for checkpoint dirs)")
    p.add_argument("--group", default="fields", help="hdf group, e.g. fields, moments, moments/species_0")
    p.add_argument("--comp", default="Ez", help="component, e.g. Ez, Bz, Jx, rho, pXX")
    p.add_argument("--cycle", type=int, default=None,
                   help="cycle for --compare (required) or to restrict plotting to")
    p.add_argument("--first-cycle", type=int, default=None)
    p.add_argument("--last-cycle", type=int, default=None)
    p.add_argument("--step", type=int, default=None)
    p.add_argument("--zindex", type=int, default=0, help="z slice to plot")
    p.add_argument("--tol", type=float, default=1e-8, help="max relative difference in compare mode")
    p.add_argument("--out-prefix", default="", help="prefix for output png files")
    p.add_argument("--vmin", type=float, default=None)
    p.add_argument("--vmax", type=float, default=None)
    args = p.parse_args()

    if args.compare:
        if args.cycle is None:
            print("--compare requires --cycle", file=sys.stderr)
            sys.exit(2)
        sys.exit(compare_directories(args.directory, args.compare, args))

    if args.cycle is not None:
        field = read_field(args.directory, args.group, args.comp, args.cycle, args.glob)
        plot_field(field, args.comp, args.cycle, args.zindex, args.out_prefix,
                   args.vmin, args.vmax)
    else:
        plot_range(args.directory, args)


if __name__ == "__main__":
    main()
