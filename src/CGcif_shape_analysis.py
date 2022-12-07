from ase.io import read
import stk
import sys
import logging


from utilities import (
    get_bonded_pairs,
    yield_hexagons,
    get_stk_molecule,
    plot_hexagons,
)
import shape_module


def main():
    if not len(sys.argv) == 4:
        logging.info(f"Usage: {__file__}\n" "   Expected 3 arguments:")
        logging.info("xyz_file (str): `.xyz` file to plot data from")
        logging.info(
            "small_bb_name (str): select bb from `bb1`, `bb2`, `bb3`"
        )
        logging.info(
            "large_bb_name (str): select bb from `bb1`, `bb2`, `bb3`"
        )
        sys.exit()
    else:
        xyz_file = sys.argv[1]
        small_bb_name = sys.argv[2]
        large_bb_name = sys.argv[3]

    ase_system = read(xyz_file)
    name = xyz_file.replace(".xyz", "")

    logging.info(
        f"getting Ti atoms in {xyz_file}, and defining bonding..."
    )
    ti_ids = []
    id_map = {}
    for i, aa in enumerate(ase_system):
        if aa.symbol == "Ti":
            ti_ids.append(i)
    id_map = {ti_ids[i]: i for i in range(len(ti_ids))}

    dists = ase_system.get_all_distances()
    bonded_pairs = get_bonded_pairs(
        atom_ids=ti_ids,
        dists=dists,
        large_bb_name=large_bb_name,
        small_bb_name=small_bb_name,
    )

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
