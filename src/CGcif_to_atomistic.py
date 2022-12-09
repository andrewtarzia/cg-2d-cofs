import numpy as np
import stk
import stko
import sys
import ase
import logging

from cif_writer import CifWriter


class CustomPeriodicTopology:
    def __init__(
        self,
        building_blocks,
        lattice_size,
        lattice_matrix,
        vertex_prototypes,
        edge_prototypes,
        vertex_alignments=None,
        reaction_factory=stk.GenericReactionFactory(),
        num_processes=1,
        optimizer=stk.NullOptimizer(),
    ):
        class InternalTopology(stk.cof.Cof):
            def _get_vertices(self, lattice):

                xdim, ydim, zdim = self._lattice_size
                num_vertices = (
                    xdim * ydim * zdim * len(self._vertex_prototypes)
                )
                vertices = [None for i in range(num_vertices)]
                for vertex in self._vertex_prototypes:
                    vertices[vertex.get_id()] = vertex

                return tuple(vertices)

            def _get_scale(self, building_block_vertices):
                return 1

            _lattice_constants = _a, _b, _c = (
                np.array(lattice_matrix[0]),
                np.array(lattice_matrix[1]),
                np.array(lattice_matrix[2]),
            )

            _vertex_prototypes = vertex_prototypes
            _edge_prototypes = edge_prototypes

        self._topology_graph = InternalTopology(
            building_blocks=building_blocks,
            lattice_size=lattice_size,
            periodic=True,
            vertex_alignments=vertex_alignments,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
            optimizer=optimizer,
        )

    def construct(self):
        return self._topology_graph.construct()


def parse_cif(path, supercell_size):
    n = int(int(supercell_size) / 2)
    coords = []
    coord_type = []

    with open(path, "r") as f:
        lines = f.readlines()
        lines = [i.split() for i in lines]

    cell_lengths = (
        float(lines[8][1]),
        float(lines[9][1]),
        float(lines[10][1]),
    )

    cell_angles = (
        float(lines[11][1]),
        float(lines[12][1]),
        float(lines[13][1]),
    )

    lattice_matrix = ase.geometry.cellpar_to_cell(
        cellpar=[
            cell_lengths[0],
            cell_lengths[1],
            cell_lengths[2],
            cell_angles[0],
            cell_angles[1],
            cell_angles[2],
        ]
    )

    for i in range(20 * n * n):
        x, y, z = lines[i + 21][1:4]
        coords.append(
            [float(i) * j for i, j in zip([x, y, z], cell_lengths)]
        )
        coord_type.append(lines[i + 21][0][:2])

    # Transform coordinates into array of shape n*3.
    coords = np.asarray(coords, dtype=np.float64, order="C")

    return {
        "coords": coords,
        "coord_type": coord_type,
        "cell_lengths": cell_lengths,
        "cell_angles": cell_angles,
        "lattice_matrix": lattice_matrix,
    }


def write_output_txt(path, lattice_matrix):
    with open(path, "w") as f:
        f.write(
            "    _lattice_constants = _a, _b, _c = (\n"
            f"       np.array([{lattice_matrix[0]}]),\n"
            f"       np.array([{lattice_matrix[1]}]),\n"
            f"       np.array([{lattice_matrix[2]}]),\n"
            ")\n\n"
        )


def get_vertex_prototypes(cif_data):

    nonlinear_vertices = {}
    linear_vertices = {}
    gamma = cif_data["cell_angles"][2] * np.pi / 180
    for i, ctype in enumerate(cif_data["coord_type"]):
        coordinates = cif_data["coords"][i]
        if ctype == "Ti":
            vertex = stk.cof.NonLinearVertex(
                id=len(nonlinear_vertices),
                position=(
                    coordinates[0] + coordinates[1] * np.cos(gamma),
                    coordinates[1] * np.sin(gamma),
                    coordinates[2],
                ),
            )
            nonlinear_vertices[i] = vertex
        else:
            vertex = stk.cof.LinearVertex(
                id=len(linear_vertices),
                position=(
                    coordinates[0] + coordinates[1] * np.cos(gamma),
                    coordinates[1] * np.sin(gamma),
                    coordinates[2],
                ),
            )
            linear_vertices[i] = vertex

    linear_vertices = {
        i: stk.cof.LinearVertex(
            id=len(nonlinear_vertices) + vert.get_id(),
            position=vert.get_position(),
        )
        for i, vert in linear_vertices.items()
    }

    return nonlinear_vertices, linear_vertices


def get_bb_dictionary(
    linear_vertices,
    nonlinear_vertices,
    coord_type,
    small_bb_name,
    large_bb_name,
):

    # Define atomistic bbs.
    bbs = {
        "bb1": stk.BuildingBlock(
            smiles="Br/N=C/c1ccc(/C=N/Br)cc1",
            functional_groups=(stk.BromoFactory(),),
        ),
        "bb2": stk.BuildingBlock(
            smiles=r"Br/N=C/c2ccc(c1ccc(/C=N\Br)cc1)cc2",
            functional_groups=(stk.BromoFactory(),),
        ),
        "bb3": stk.BuildingBlock(
            smiles=(
                r"C1=CC(C2=CC=C(C3=CC=C(/C=N/Br)C=C3)C=C2)=CC=C1/C=N/Br"
            ),
            functional_groups=(stk.BromoFactory(),),
        ),
    }
    tritopic = stk.BuildingBlock(
        smiles="Brc4ccc(c3cc(c1ccc(Br)cc1)cc(c2ccc(Br)cc2)c3)cc4",
        functional_groups=(stk.BromoFactory(),),
    )

    small_bb_in_model = bbs[small_bb_name]
    large_bb_in_model = bbs[large_bb_name]
    logging.info(f"selected small bb: {small_bb_in_model}")
    logging.info(f"selected large bb: {large_bb_in_model}")

    bb1_values = []
    bb2_values = []
    for i in linear_vertices:
        vert = linear_vertices[i]
        id_ = vert.get_id()
        if coord_type[i] == "Mn":
            bb2_values.append(id_)
        # elif coord_type[i] == "Pb":
        #     bb2_values.append(id_)
        elif coord_type[i] == "Pb":
            bb1_values.append(id_)

    tritopic_values = []
    for i in nonlinear_vertices:
        vert = nonlinear_vertices[i]
        tritopic_values.append(vert.get_id())

    return {
        small_bb_in_model: tuple(bb1_values),
        large_bb_in_model: tuple(bb2_values),
        tritopic: tuple(tritopic_values),
    }


def run_construction(
    linear_vertices,
    nonlinear_vertices,
    edge_prototypes,
    cif_data,
    small_bb_name,
    large_bb_name,
):
    coord_type = cif_data["coord_type"]
    lattice_matrix = cif_data["lattice_matrix"]

    bb_dict = get_bb_dictionary(
        linear_vertices=linear_vertices,
        nonlinear_vertices=nonlinear_vertices,
        coord_type=coord_type,
        small_bb_name=small_bb_name,
        large_bb_name=large_bb_name,
    )

    lattice_size = (1, 1, 1)

    vertex_prototypes = tuple(
        nonlinear_vertices[i] for i in nonlinear_vertices
    ) + tuple(linear_vertices[i] for i in linear_vertices)

    topology_graph = CustomPeriodicTopology(
        building_blocks=bb_dict,
        lattice_size=lattice_size,
        lattice_matrix=lattice_matrix,
        vertex_prototypes=vertex_prototypes,
        edge_prototypes=edge_prototypes,
    )

    construction_result = topology_graph.construct()
    cof = stk.ConstructedMolecule.init_from_construction_result(
        construction_result=construction_result,
    )
    unit_cell = stko.UnitCell(
        vector_1=lattice_matrix[0],
        vector_2=lattice_matrix[1],
        vector_3=lattice_matrix[2],
    )

    return cof, unit_cell


def get_periodicity(
    nlc,
    lc,
    lattice_constant_a,
    lattice_constant_b,
    lattice_constant_c,
    gamma,
):

    vec1 = np.zeros((2))
    vec1[0] = (nlc[0] + nlc[1] * np.cos(gamma)) - (
        lc[0] + lc[1] * np.cos(gamma)
    )
    vec1[1] = (nlc[1] - lc[1]) * np.sin(gamma)

    # +x
    vec2 = np.zeros((2))
    vec2[0] = (nlc[0] + lattice_constant_a + nlc[1] * np.cos(gamma)) - (
        lc[0] + lc[1] * np.cos(gamma)
    )
    vec2[1] = (nlc[1] - lc[1]) * np.sin(gamma)

    # +x+y
    vec3 = np.zeros((2))
    vec3[0] = (
        nlc[0]
        + lattice_constant_a
        + (nlc[1] - lattice_constant_b) * np.cos(gamma)
    ) + (lc[0] - lc[1] * np.cos(gamma))
    vec3[1] = (nlc[1] - lattice_constant_b - lc[1]) * np.sin(gamma)

    # +y
    vec4 = np.zeros((2))
    vec4[0] = (
        nlc[0] + (nlc[1] - lattice_constant_b) * np.cos(gamma)
    ) - (lc[0] + lc[1] * np.cos(gamma))
    vec4[1] = (nlc[1] - lattice_constant_b - lc[1]) * np.sin(gamma)

    # -x
    vec5 = np.zeros((2))
    vec5[0] = (nlc[0] + nlc[1] * np.cos(gamma)) - (
        lc[0] + lattice_constant_a + lc[1] * np.cos(gamma)
    )
    vec5[1] = (nlc[1] - lc[1]) * np.sin(gamma)

    # -x-y
    vec6 = np.zeros((2))
    vec6[0] = (nlc[0] + (nlc[1]) * np.cos(gamma)) - (
        lc[0]
        + lattice_constant_a
        + (lc[1] - lattice_constant_b) * np.cos(gamma)
    )
    vec6[1] = (nlc[1] - (lc[1] - lattice_constant_b)) * np.sin(gamma)

    # -y
    vec7 = np.zeros((2))
    vec7[0] = (nlc[0] + nlc[1] * np.cos(gamma)) - (
        lc[0] + (lc[1] - lattice_constant_b) * np.cos(gamma)
    )
    vec7[1] = (nlc[1] - (lc[1] - lattice_constant_b)) * np.sin(gamma)

    # x-y
    vec8 = np.zeros((2))
    vec8[0] = (nlc[0] + lattice_constant_a + nlc[1] * np.cos(gamma)) - (
        lc[0] + (lc[1] - lattice_constant_b) * np.cos(gamma)
    )
    vec8[1] = (nlc[1] - (lc[1] - lattice_constant_b)) * np.sin(gamma)

    # -x+y
    vec9 = np.zeros((2))
    vec9[0] = (
        nlc[0] + (nlc[1] - lattice_constant_b) * np.cos(gamma)
    ) - (lc[0] + lattice_constant_a + (lc[1]) * np.cos(gamma))
    vec9[1] = (nlc[1] - lattice_constant_b - (lc[1])) * np.sin(gamma)

    vec_mag_1 = np.sqrt(vec1[1] ** 2 + vec1[0] ** 2)
    vec_mag_2 = np.sqrt(vec2[1] ** 2 + vec2[0] ** 2)
    vec_mag_3 = np.sqrt(vec3[1] ** 2 + vec3[0] ** 2)
    vec_mag_4 = np.sqrt(vec4[1] ** 2 + vec4[0] ** 2)
    vec_mag_5 = np.sqrt(vec5[1] ** 2 + vec5[0] ** 2)
    vec_mag_6 = np.sqrt(vec6[1] ** 2 + vec6[0] ** 2)
    vec_mag_7 = np.sqrt(vec7[1] ** 2 + vec7[0] ** 2)
    vec_mag_8 = np.sqrt(vec8[1] ** 2 + vec8[0] ** 2)
    vec_mag_9 = np.sqrt(vec9[1] ** 2 + vec9[0] ** 2)

    no_edge = False
    if vec_mag_1 < 20:
        periodicity = None
    elif vec_mag_2 < 20:
        periodicity = (-1, 0, 0)
    elif vec_mag_3 < 20:
        periodicity = (-1, 1, 0)
    elif vec_mag_4 < 20:
        periodicity = (0, 1, 0)
    elif vec_mag_5 < 20:
        periodicity = (1, 0, 0)
    elif vec_mag_6 < 20:
        periodicity = (1, -1, 0)
    elif vec_mag_7 < 20:
        periodicity = (0, -1, 0)
    elif vec_mag_8 < 20:
        periodicity = (-1, -1, 0)
    elif vec_mag_9 < 20:
        periodicity = (1, 1, 0)
    else:
        no_edge = True
        periodicity = None

    return periodicity, no_edge


def get_edge_prototypes(
    cif_data,
    nonlinear_vertices,
    linear_vertices,
):

    coordinates = cif_data["coords"]
    cell_lengths = cif_data["cell_lengths"]
    cell_angles = cif_data["cell_angles"]

    edge_prototypes = []
    for nl_c_id, nlv in nonlinear_vertices.items():
        for l_c_id, lv in linear_vertices.items():
            edge_id = len(edge_prototypes)

            # If statements to get periodicity.
            periodicity, no_edge = get_periodicity(
                nlc=coordinates[nl_c_id],
                lc=coordinates[l_c_id],
                lattice_constant_a=cell_lengths[0],
                lattice_constant_b=cell_lengths[1],
                lattice_constant_c=cell_lengths[2],
                gamma=cell_angles[2] * np.pi / 180.0,
            )

            if no_edge:
                continue

            if periodicity is None:
                edge = stk.Edge(
                    id=edge_id,
                    vertex1=nlv,
                    vertex2=lv,
                )
            else:
                edge = stk.Edge(
                    id=edge_id,
                    vertex1=nlv,
                    vertex2=lv,
                    periodicity=periodicity,
                )

            edge_prototypes.append(edge)

    return edge_prototypes


def main():
    if not len(sys.argv) == 5:
        logging.info(f"Usage: {__file__}\n" "   Expected 4 arguments:")
        logging.info(
            "cif_file (str): `.cif` file to python script from"
        )
        logging.info(
            "run_optimisation (str): `t` if true, anything otherwise"
        )
        logging.info(
            "small_bb_name (str): select bb from `bb1`, `bb2`, `bb3`"
        )
        logging.info(
            "large_bb_name (str): select bb from `bb1`, `bb2`, `bb3`"
        )
        sys.exit()
    else:
        cif_file = sys.argv[1]
        run_optimisation = True if sys.argv[2] == "t" else False
        small_bb_name = sys.argv[3]
        large_bb_name = sys.argv[4]

    name = cif_file.replace(".cif", "")
    _, t_on_j, supercell_size = name.split("_")
    sim_name = f"cof_n_{supercell_size}_T_{t_on_j}"
    logging.info(f"variables: n={supercell_size}/2, T/J={t_on_j}")
    cif_data = parse_cif(cif_file, supercell_size=int(supercell_size))

    write_output_txt(
        path=f"output_{t_on_j}.txt",
        lattice_matrix=cif_data["lattice_matrix"],
    )

    # Define vertex prototypes.
    nl_vertices, l_vertices = get_vertex_prototypes(cif_data)
    edge_prototypes = get_edge_prototypes(
        cif_data=cif_data,
        nonlinear_vertices=nl_vertices,
        linear_vertices=l_vertices,
    )

    # Run construction.
    cof, unit_cell = run_construction(
        linear_vertices=l_vertices,
        nonlinear_vertices=nl_vertices,
        edge_prototypes=edge_prototypes,
        cif_data=cif_data,
        small_bb_name=small_bb_name,
        large_bb_name=large_bb_name,
    )

    CifWriter().write(
        molecule=cof,
        path=f"{sim_name}.cif",
        periodic_info=unit_cell,
    )
    stk.XyzWriter().write(
        molecule=cof,
        path=f"{sim_name}.xyz",
    )

    if run_optimisation:
        logging.info("optimising the cof with Gulp...")
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path="/home/atarzia/software/gulp-6.1/Src/gulp",
            output_dir=f"{name}_gulpopt",
            maxcyc=50,
            conjugate_gradient=True,
        )
        gulp_opt.assign_FF(cof)
        cof, unit_cell = gulp_opt.p_optimize(
            mol=cof,
            unit_cell=unit_cell,
        )
        CifWriter().write(
            molecule=cof,
            path=f"{sim_name}_opt.cif",
            periodic_info=unit_cell,
        )
        stk.XyzWriter().write(
            molecule=cof,
            path=f"{sim_name}_opt.xyz",
        )
        stk.MolWriter().write(
            molecule=cof,
            path=f"{sim_name}_opt.mol",
        )
        stk.PdbWriter().write(
            molecule=cof,
            path=f"{sim_name}_opt.pdb",
            periodic_info=unit_cell,
        )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
