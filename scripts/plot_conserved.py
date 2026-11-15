#!/usr/bin/env python3

import os
import glob
import numpy as np
import matplotlib.pyplot as plt


# Column descriptions
columns = {
    1: "Cycle",
    2: "Total internal energy",
    3: "Variation of total internal energy",
    4: "Electric energy",
    5: "Local magnetic energy",
    6: "Kinetic energy",
    7: "Momentum",
    8: "Total magnetic energy",
    9: "Internal magnetic energy",
    10: "Kinetic energy removed",
    11: "Electric energy removed",
    12: "Number of particles of species 0",
    13: "Total charge of species 0",
    14: "Kinetic energy of species 0",
    15: "Number of particles of species 1",
    16: "Total charge of species 1",
    17: "Kinetic energy of species 1",
}


def print_columns():
    print("\nAvailable columns:")
    for key, value in columns.items():
        print(f"{key:2d} -> {value}")


import os

def select_files():
    n = int(input("\nHow many files do you want to plot? "))

    selected = []

    for i in range(n):
        while True:
            folder = input(f"Enter folder path of file {i + 1}: ").strip()

            filename = os.path.join(folder, "data", "ConservedQuantities.txt")

            if os.path.isfile(filename):
                selected.append(filename)
                break
            else:
                print(f"File not found: {filename}. Try again.")

    return selected

def load_column(filename, column):
    """
    Load selected column.
    Column indexing starts at 1.
    """

    try:
        data = np.loadtxt(filename)

    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None, None

    if column > data.shape[1]:
        raise ValueError(
            f"{filename} has only {data.shape[1]} columns."
        )

    cycle = data[:, 0]
    values = data[:, column - 1]
    print(values[0:10])

    return cycle, values


def main():

    files = select_files()

    print_columns()

    column = int(input("\nWhich column do you want to plot? "))

    if column not in columns:
        print("Invalid column.")
        return

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_axes([0.2, 0.2, 0.55, 0.75])
    
    linestyles = ["-", "--", ":", "-."]
    #linecolor = ["k","r"]    

    for i, filename in enumerate(files):

        print(filename)
        cycle, values = load_column(filename, column)

        if cycle is not None:
            ax.plot(
                cycle,
                values,
                linestyle=linestyles[i % len(linestyles)],
                label = os.path.basename(os.path.dirname(os.path.dirname(filename)))
            )

    ax.set_yscale("log")
    ax.set_xlabel(columns[1], fontname='Helvetica', fontsize=28)
    ax.set_ylabel(columns[column], fontname='Helvetica', fontsize=28)
    
    ax.tick_params(labelsize=28)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
       label.set_fontname('Helvetica')
    #plt.title(columns[column], fontname='Helvetica', fontsize=28)

    #plt.grid(True)
    #plt.legend(fontsize=28)

    #plt.tight_layout(rect=[0, 0, 1, 1])
    #plt.show()
    plt.savefig("comparison.pdf", format="pdf")


if __name__ == "__main__":
    main()
