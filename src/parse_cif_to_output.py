from ase.io import read
import numpy as np
import stk
import rdkit.Chem.AllChem as rdkit
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import sys

import shape_module

import logging


def get_bonded_pairs(atom_ids, dists):
    bonded_pairs = set()
    for atom_id in atom_ids:
        dds = dists[atom_id]
        logging.warning(
            "THIS NEEDS WORK, TO MODIFY LENGTHS FOR SOME SYSTEMS!"
        )
        long_bonded_ids = [
            i
            for i in np.argwhere((dds < 26) & (dds > 22))
            if i in atom_ids
        ]
        short_bonded_ids = [
            i
            for i in np.argwhere((dds < 22) & (dds > 20))
            if i in atom_ids
        ]
        # long_bonded_ids = [
        #     i
        #     for i in np.argwhere((dds < 31) & (dds > 28))
        #     if i in atom_ids
        # ]
        # short_bonded_ids = [
        #     i
        #     for i in np.argwhere((dds < 25) & (dds > 20))
        #     if i in atom_ids
        # ]
        if len(long_bonded_ids) + len(short_bonded_ids) != 3:
            continue
        all_bonded_ids = short_bonded_ids + long_bonded_ids
        for i in all_bonded_ids:
            bonded_pairs.add(tuple(sorted((atom_id, int(i)))))

    return bonded_pairs


def get_stk_molecule(ase_system, id_map, bonded_pairs, name):
    stk_mol = stk.BuildingBlock.init(
        atoms=tuple(
            stk.C(id=id_map[i])
            for i, aa in enumerate(ase_system)
            if aa.symbol == "Ti"
        ),
        bonds=tuple(
            stk.Bond(stk.C(id_map[id1]), stk.C(id_map[id2]), 1)
            for id1, id2 in bonded_pairs
        ),
        position_matrix=np.array(
            tuple(i.position for i in ase_system if i.symbol == "Ti")
        ),
    )
    stk_mol.write(f"{name}_just_ti.mol")
    return stk_mol


def yield_hexagons(stk_mol):
    logging.warning("THIS NEEDS WORK, IT IS NOT FINDING ALL HEXAGONS!")
    smarts = "[#6]1[#6][#6][#6][#6][#6]1"
    rdkit_mol = stk_mol.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)

    yield from rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts(smarts),
        uniquify=True,
        maxMatches=10000000,
    )


def plot_hexagons(
    name,
    all_atom_coordinates,
    all_atom_coordinates_w_shape,
):

    # Matplotlib thing.
    fig, axs = plt.subplots(
        ncols=2,
        figsize=(9, 8),
        gridspec_kw={"width_ratios": [12, 1]},
    )
    ax = axs[0]

    ax.scatter(
        [i[0] for i in all_atom_coordinates],
        [i[1] for i in all_atom_coordinates],
        c="k",
        s=20,
        edgecolor="none",
        alpha=1.0,
    )

    cmap = plt.get_cmap("Reds")
    shape_min = 0
    shape_max = 10
    norm = mpl.colors.Normalize(vmin=shape_min, vmax=shape_max)

    for hexid in all_atom_coordinates_w_shape:
        xs = []
        ys = []
        for aic in all_atom_coordinates_w_shape[hexid]:
            x, y, _, c = aic
            xs.append(x)
            ys.append(y)
            # Normalize.
            col = (c - shape_min) / (shape_max - shape_min)
        ax.fill(
            xs,
            ys,
            facecolor=cmap(col),
            edgecolor="k",
            linewidth=1,
        )

    # ax.tick_params(axis='both', which='major', labelsize=16)
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    # ax.set_xlabel('X [A]', fontsize=16)
    # ax.set_ylabel('Y [A]', fontsize=16)

    cb1 = mpl.colorbar.ColorbarBase(
        axs[1],
        cmap=cmap,
        norm=norm,
    )
    cb1.set_label("HP-6", fontsize=16)
    cb1.ax.tick_params(labelsize=16)

    fig.tight_layout()
    fig.savefig(f"{name}_hexes.pdf", dpi=720, bbox_inches="tight")
    plt.close()


def main():
    if not len(sys.argv) == 2:
        logging.info(f"Usage: {__file__}\n" "   Expected 1 arguments:")
        logging.info("cif_file (str): `.cif` file to plot data from")
        sys.exit()
    else:
        cif_file = sys.argv[1]

    ase_system = read(cif_file)
    name = cif_file.replace(".cif", "")

    logging.info(
        f"getting Ti atoms in {cif_file}, and defining bonding..."
    )
    ti_ids = []
    id_map = {}
    for i, aa in enumerate(ase_system):
        if aa.symbol == "Ti":
            ti_ids.append(i)
    id_map = {ti_ids[i]: i for i in range(len(ti_ids))}

    dists = ase_system.get_all_distances()
    bonded_pairs = get_bonded_pairs(atom_ids=ti_ids, dists=dists)

    logging.info(f"defining an stk molecule to {name}_just_ti.mol...")
    stk_mol = get_stk_molecule(ase_system, id_map, bonded_pairs, name)

    all_shapes = []
    centroid_strings = []
    all_atom_strings = []
    all_centroid_coordinates = []
    all_atom_coordinates_w_shape = {}
    for hex_id, atom_ids in enumerate(yield_hexagons(stk_mol)):
        stk_mol.write("temp.mol", atom_ids=atom_ids)
        new_mol = stk.BuildingBlock.init_from_file("temp.mol")

        logging.info("calculating a shape...")
        shape = shape_module.calculate_hp6(new_mol, f"s_{hex_id}")
        all_shapes.append(shape)

        # This is the coordinate of the hexagon centroid.
        cent = new_mol.get_centroid()
        all_centroid_coordinates.append(cent)

        # Added another column with the hexagon id, so you can see the
        # same ones.
        centroid_strings.append(
            f"H {round(cent[0], 2)} {round(cent[1], 2)} "
            f"{round(cent[2], 2)} {hex_id} {shape}\n"
        )

        # This part writes all hexagon atoms to a bigger xyz file.
        coords = new_mol.get_position_matrix()
        for i, atom_id in enumerate(range(new_mol.get_num_atoms()), 1):
            x, y, z = (i for i in coords[atom_id])
            # Added another column with the hexagon id, so you can see
            # the same ones.
            all_atom_strings.append(
                f"C {round(x, 2)} {round(y, 2)} {round(z, 2)} {hex_id}"
                f" {shape}\n"
            )
            if hex_id not in all_atom_coordinates_w_shape:
                all_atom_coordinates_w_shape[hex_id] = []
            all_atom_coordinates_w_shape[hex_id].append(
                (round(x, 2), round(y, 2), round(z, 2), shape)
            )

    plot_hexagons(
        name=name,
        all_atom_coordinates=stk_mol.get_position_matrix(),
        all_atom_coordinates_w_shape=all_atom_coordinates_w_shape,
    )

    with open(f"{name}_hexagons.xyz", "w") as f:
        f.write(f"{len(centroid_strings)}\nnAll hexagon centroids\n")
        for line in centroid_strings:
            f.write(line)

    with open(f"{name}_allatoms.xyz", "w") as f:
        f.write(f"{len(all_atom_strings)}\nAll hexagon atoms\n")
        for line in all_atom_strings:
            f.write(line)

    with open(f"{name}.shapes", "w") as f:
        for sh in all_shapes:
            f.write(f"{sh}\n")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
