#!/usr/bin/env python3
"""
Plot / compare ConservedQuantities.txt files produced by iPIC3D.

File format (written by c_Solver::WriteConserved in main/iPIC3Dlib.cpp, one row
every DiagnosticsOutputCycle, tab separated; an optional text header is written
at cycle 0 and is skipped here):

    cycle  <empty>  <empty>  Total_Energy  Momentum  Eenergy  Benergy  Kenergy  Ke[species 0..ns-1]  BulkEnergy[species 0..ns-1]

Whitespace-collapsed, the columns are therefore (1-based):

    1  Cycle
    2  Total energy (E + B + K)
    3  Total momentum
    4  Electric energy
    5  Magnetic energy
    6  Kinetic energy (total)
    7+i        Kinetic energy, species i
    7+ns+i     Bulk energy, species i

Usage:
  interactive (no arguments): asks for folders and a column, as before
  plot several runs:
    plot_conserved.py --files run1 run2 ... --column 2 --out comparison.pdf
  compare a test run against a reference at the cycles they have in common
  (exit code 1 if any requested column deviates more than --tol):
    plot_conserved.py --compare REF_DIR TEST_DIR --columns 2,3,4,5,6 --tol 1e-8

Each entry of --files / --compare may be a run directory (containing
ConservedQuantities.txt or data/ConservedQuantities.txt) or the .txt file itself.
"""

import os
import sys
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------

def resolve_file(path):
    """Accept a run directory or a direct ConservedQuantities.txt path."""
    if os.path.isfile(path):
        return path
    for cand in (os.path.join(path, "ConservedQuantities.txt"),
                 os.path.join(path, "data", "ConservedQuantities.txt")):
        if os.path.isfile(cand):
            return cand
    raise FileNotFoundError(
        f"No ConservedQuantities.txt found in '{path}' "
        "(looked for <dir>/ConservedQuantities.txt and <dir>/data/ConservedQuantities.txt)"
    )


def load_table(filename):
    """
    Load a ConservedQuantities.txt file.

    Returns (columns_dict, data) where data is a 2D float array with
    column 0 = cycle. The optional text header line is skipped.
    """
    rows = []
    with open(filename) as f:
        for line in f:
            fields = line.split()
            if not fields:
                continue
            try:
                int(fields[0])
            except ValueError:
                continue  # header line ("Total_Energy ...")
            rows.append([float(x) for x in fields])
    if not rows:
        raise ValueError(f"No numeric rows found in {filename}")
    data = np.array(rows)
    ncols = data.shape[1]
    ns = (ncols - 6) // 2  # cycle + 5 scalars + ns Ke + ns Bulk
    if ns < 0 or 6 + 2 * ns != ncols:
        raise ValueError(
            f"{filename}: unexpected number of columns ({ncols}); "
            "expected cycle + 5 scalars + 2*ns species columns"
        )
    columns = {
        1: "Cycle",
        2: "Total energy (E+B+K)",
        3: "Total momentum",
        4: "Electric energy",
        5: "Magnetic energy",
        6: "Kinetic energy",
    }
    for i in range(ns):
        columns[7 + i] = f"Kinetic energy species {i}"
    for i in range(ns):
        columns[7 + ns + i] = f"Bulk energy species {i}"
    return columns, data


def load_column(filename, column):
    """Return (cycles, values) for a 1-based column number."""
    columns, data = load_table(filename)
    if column not in columns or column < 1 or column > data.shape[1]:
        raise ValueError(
            f"{filename}: invalid column {column} "
            f"(file has {data.shape[1]} columns)"
        )
    return data[:, 0], data[:, column - 1]


def label_for(path):
    if os.path.isdir(path):
        return os.path.basename(os.path.normpath(path))
    f = resolve_file(path)
    d = os.path.dirname(f)
    if os.path.basename(d) == "data":
        d = os.path.dirname(d)
    return os.path.basename(d)


def print_columns(columns):
    print("\nAvailable columns:")
    for key, value in columns.items():
        print(f"{key:2d} -> {value}")


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def plot_files(files, column, out, yscale):
    columns = None
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_axes([0.2, 0.2, 0.55, 0.75])
    linestyles = ["-", "--", ":", "-."]
    for i, path in enumerate(files):
        filename = resolve_file(path)
        columns, data = load_table(filename)
        if column not in columns:
            raise ValueError(f"{filename}: no column {column}")
        ax.plot(data[:, 0], data[:, column - 1],
                linestyle=linestyles[i % len(linestyles)],
                label=label_for(path))
    ax.set_yscale(yscale)
    ax.set_xlabel(columns[1], fontsize=14)
    ax.set_ylabel(columns[column], fontsize=14)
    ax.legend(fontsize=10)
    fig.savefig(out)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------

def compare_pair(ref_path, test_path, columns, tol, out, yscale):
    ref_file = resolve_file(ref_path)
    test_file = resolve_file(test_path)
    columns_dict, ref = load_table(ref_file)
    _, test = load_table(test_file)

    ref_cycles = ref[:, 0].astype(int)
    test_cycles = test[:, 0].astype(int)
    common = np.intersect1d(ref_cycles, test_cycles)
    if common.size == 0:
        print(f"FAIL: no common cycles between {ref_file} and {test_file}",
              file=sys.stderr)
        return 1
    ref_idx = {c: i for i, c in enumerate(ref_cycles)}
    test_idx = {c: i for i, c in enumerate(test_cycles)}

    print(f"Reference : {ref_file}")
    print(f"Test      : {test_file}")
    print(f"Common cycles: {common.size} "
          f"({common[0]} .. {common[-1]})")
    print(f"{'column':>6}  {'description':<28}  {'max rel err':>12}  status")

    # one multi-panel plot of the compared columns
    n = len(columns)
    fig, axes = plt.subplots(n, 1, figsize=(7, 3 * n), squeeze=False)
    failed = False
    for row, col in enumerate(columns):
        if col not in columns_dict:
            print(f"FAIL: column {col} missing in {ref_file}", file=sys.stderr)
            return 1
        r = np.array([ref[ref_idx[c], col - 1] for c in common])
        t = np.array([test[test_idx[c], col - 1] for c in common])
        denom = np.max(np.abs(r))
        if denom == 0.0:
            err = np.max(np.abs(r - t))  # absolute if reference is identically zero
            rel = err
        else:
            rel = np.max(np.abs(r - t)) / denom
        ok = rel <= tol
        failed = failed or not ok
        print(f"{col:>6}  {columns_dict[col]:<28}  {rel:>12.3e}  "
              f"{'OK' if ok else 'FAIL (tol=%g)' % tol}")
        ax = axes[row, 0]
        ax.plot(common, r, "-", label=label_for(ref_path))
        ax.plot(common, t, "--", label=label_for(test_path))
        ax.set_yscale(yscale)
        ax.set_ylabel(columns_dict[col], fontsize=9)
        ax.legend(fontsize=8)
        if row == n - 1:
            ax.set_xlabel(columns_dict[1], fontsize=9)
    fig.savefig(out)
    print(f"Saved {out}")
    if failed:
        print("COMPARE: FAILED", file=sys.stderr)
        return 1
    print("COMPARE: OK")
    return 0


# ---------------------------------------------------------------------------
# interactive mode (original behaviour)
# ---------------------------------------------------------------------------

def interactive():
    n = int(input("\nHow many files do you want to plot? "))
    selected = []
    for i in range(n):
        while True:
            folder = input(f"Enter folder path of file {i + 1}: ").strip()
            try:
                selected.append(resolve_file(folder))
                break
            except FileNotFoundError as e:
                print(e)
    columns, _ = load_table(selected[0])
    print_columns(columns)
    column = int(input("\nWhich column do you want to plot? "))
    if column not in columns:
        print("Invalid column.")
        return
    plot_files(selected, column, "comparison.pdf", "log")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--files", nargs="+",
                        help="run directories or ConservedQuantities.txt files to plot")
    parser.add_argument("--column", type=int, default=2,
                        help="1-based column to plot (default: 2 = total energy)")
    parser.add_argument("--compare", nargs=2, metavar=("REF", "TEST"),
                        help="compare TEST against REF at common cycles")
    parser.add_argument("--columns", default="2,3,4,5,6",
                        help="comma-separated columns for --compare (default: 2,3,4,5,6)")
    parser.add_argument("--tol", type=float, default=1e-8,
                        help="max relative error accepted in --compare mode")
    parser.add_argument("--out", default="comparison.pdf",
                        help="output figure (plot mode)")
    parser.add_argument("--yscale", default="log", choices=["log", "linear", "symlog"],
                        help="y scale of the plots (default: log)")
    args = parser.parse_args()

    if args.compare:
        ref, test = args.compare
        columns = [int(c) for c in args.columns.split(",")]
        out = args.out if args.out != "comparison.pdf" else "comparison.pdf"
        sys.exit(compare_pair(ref, test, columns, args.tol, out, args.yscale))

    if args.files:
        plot_files(args.files, args.column, args.out, args.yscale)
        return

    interactive()


if __name__ == "__main__":
    main()
