import matplotlib.pyplot as plt
import numpy as np
import sys

import logging


def main():
    if not len(sys.argv) == 2:
        logging.info(f"Usage: {__file__}\n" "   Expected 1 arguments:")
        logging.info(
            "shape_file (str): `.shapes` file to plot data from"
        )
        sys.exit()
    else:
        shape_file = sys.argv[1]

    name = shape_file.replace(".shapes", "")
    output_file = shape_file.replace(".shapes", ".pdf")

    fig, ax = plt.subplots(figsize=(8, 5))

    xmin = 0
    xmax = 10
    width = 0.2
    X_bins = np.arange(xmin, xmax + width * 2, width)

    xdata = []
    with open(shape_file, "r") as f:
        for line in f.readlines():
            xdata.append(float(line.strip()))

    hist, bin_edges = np.histogram(
        a=xdata,
        bins=X_bins,
        density=True,
    )
    ax.plot(
        bin_edges[:-1],
        hist,
        lw=3,
        label=name,
    )

    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.set_xlabel("HP-6", fontsize=16)
    ax.set_ylabel("density", fontsize=16)
    ax.legend(fontsize=16)

    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()
    fig.savefig(output_file, dpi=720, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
