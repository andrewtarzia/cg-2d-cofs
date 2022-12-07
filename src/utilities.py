import numpy as np
import stk
import rdkit.Chem.AllChem as rdkit
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker


def bb_name_to_dist_ranges(bb_name):
    return {
        "bb1": (19, 22),
        "bb2": (24, 26),
        "bb3": (29, 31),
    }[bb_name]


def get_bonded_pairs(atom_ids, dists, small_bb_name, large_bb_name):

    short_range = bb_name_to_dist_ranges(small_bb_name)
    large_range = bb_name_to_dist_ranges(large_bb_name)

    bonded_pairs = set()
    for atom_id in atom_ids:
        dds = dists[atom_id]
        long_bonded_ids = [
            i
            for i in np.argwhere(
                (dds < large_range[1]) & (dds > large_range[0])
            )
            if i in atom_ids
        ]
        short_bonded_ids = [
            i
            for i in np.argwhere(
                (dds < short_range[1]) & (dds > short_range[0])
            )
            if i in atom_ids
        ]
        all_bonded_ids = short_bonded_ids + long_bonded_ids
        if len(all_bonded_ids) != 3:
            continue

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
        figsize=(14, 8),
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
