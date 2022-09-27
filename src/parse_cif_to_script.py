from ase.io import read
import numpy as np
import stk
import rdkit.Chem.AllChem as rdkit
import sys


import logging


def main():
    if not len(sys.argv) == 2:
        logging.info(f"Usage: {__file__}\n" "   Expected 1 arguments:")
        logging.info(
            "cif_file (str): `.cif` file to write topology from"
        )
        sys.exit()
    else:
        cif_file = sys.argv[1]


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
