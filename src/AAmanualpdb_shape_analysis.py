import stk
import rdkit.Chem.AllChem as rdkit
from ase import Atoms
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
        logging.info("pdb_file (str): `.pdb` file to plot data from")
        logging.info(
            "small_bb_name (str): select bb from `bb1`, `bb2`, `bb3`"
        )
        logging.info(
            "large_bb_name (str): select bb from `bb1`, `bb2`, `bb3`"
        )
        sys.exit()
    else:
        pdb_file = sys.argv[1]
        small_bb_name = sys.argv[2]
        large_bb_name = sys.argv[3]

    prefix = pdb_file.replace(".pdb", "")
    stk_mol = stk.BuildingBlock.init_from_file(pdb_file)

    smarts = "[#6]~[#6]1~[#6]~[#6](~[#6]~[#6](~[#6]1)~[#6])~[#6]"
    rdkit_mol = stk_mol.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    ti_atoms = []
    for atom_ids in rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts(smarts),
    ):
        centroid = stk_mol.get_centroid(atom_ids=atom_ids)
        ti_atoms.append(centroid)

    with open(f"{prefix}_justti.xyz", "w") as f:
        f.write(f"{len(ti_atoms)}\n\n")
        for a in ti_atoms:
            f.write(
                f"C {round(a[0], 2)} {round(a[1], 2)} "
                f"{round(a[2], 2)}\n"
            )

    ase_system = Atoms(
        f"Ti{len(ti_atoms)}", positions=[i for i in ti_atoms]
    )

    dists = ase_system.get_all_distances()
    bonded_pairs = get_bonded_pairs(
        atom_ids=[i for i in range(len(ti_atoms))],
        dists=dists,
        large_bb_name=large_bb_name,
        small_bb_name=small_bb_name,
    )

    logging.info(f"defining an stk molecule to {prefix}_just_ti.mol...")
    id_map = {i: i for i in range(len(ti_atoms))}
    stk_mol = get_stk_molecule(ase_system, id_map, bonded_pairs, prefix)
    stk_mol.write(f"{prefix}_just_ti.mol")

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
        name=prefix,
        all_atom_coordinates=stk_mol.get_position_matrix(),
        all_atom_coordinates_w_shape=all_atom_coordinates_w_shape,
    )

    with open(f"{prefix}_hexagons.xyz", "w") as f:
        f.write(f"{len(centroid_strings)}\nnAll hexagon centroids\n")
        for line in centroid_strings:
            f.write(line)

    with open(f"{prefix}_allatoms.xyz", "w") as f:
        f.write(f"{len(all_atom_strings)}\nAll hexagon atoms\n")
        for line in all_atom_strings:
            f.write(line)

    with open(f"{prefix}.shapes", "w") as f:
        for sh in all_shapes:
            f.write(f"{sh}\n")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
