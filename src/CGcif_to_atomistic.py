import numpy as np
import stk
import stko
import sys
import logging

from cif_writer import CifWriter


class NewNonLinearVertex(stk.cof.NonLinearVertex):
    def place_building_block(self, building_block, edges):
        logging.info("check why this is needed.")
        assert building_block.get_num_functional_groups() > 2, (
            f"{building_block} needs to have more than 2 functional "
            "groups but has "
            f"{building_block.get_num_functional_groups()}."
        )
        edges = sorted(edges, key=lambda edge: edge.get_parent_id())
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        normal = building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = stk.get_acute_vector(
            reference=core_centroid - self._position,
            vector=normal,
        )
        building_block = building_block.with_rotation_between_vectors(
            start=normal,
            target=[0, 0, 1],
            origin=self._position,
        )
        (fg,) = building_block.get_functional_groups(0)
        fg_centroid = building_block.get_centroid(fg.get_placer_ids())
        edge_position = edges[self._aligner_edge].get_position()
        return building_block.with_rotation_to_minimize_angle(
            start=fg_centroid - self._position,
            target=edge_position - self._position,
            axis=np.array([0, 0, 1], dtype=np.float64),
            origin=self._position,
        ).get_position_matrix()


class CustomPeriodicTopology:
    def __init__(
        self,
        building_blocks,
        lattice_size,
        lattice_constant_a,
        lattice_constant_b,
        lattice_constant_c,
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
                np.array([lattice_constant_a, 0.0, 0.0]),
                np.array(
                    [
                        -lattice_constant_a * 0.5,
                        -lattice_constant_b * np.sqrt(3) / 2,
                        0.0,
                    ]
                ),
                np.array([0.0, 0.0, lattice_constant_c]),
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


def parse_cif(path, n):
    logging.info("check this variable")

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

    for i in range(20 * n * n):
        x, y, z = lines[i + 21][1:4]
        coords.append(
            [float(i) * j for i, j in zip([x, y, z], cell_lengths)]
        )
        coord_type.append(lines[i + 21][0][:2])

    # Transform coordinates into array of shape n*3.
    coords = np.asarray(coords, dtype=np.float64, order="C")
    # coords = coords.reshape(20 * n * n, 3)
    # print(coords.shape)
    # for i in range(20 * n * n):
    #     coords[i][0] = coords[i][0] * float(cell_lengths[0])
    #     coords[i][1] = coords[i][1] * float(cell_lengths[1])
    #     coords[i][2] = coords[i][2] * float(cell_lengths[2])
    # print(coords.shape)
    # raise SystemExit()
    return {
        "coords": coords,
        "coord_type": coord_type,
        "cell_lengths": cell_lengths,
    }


def write_output_txt(
    path,
    lattice_constant_a,
    lattice_constant_b,
    lattice_constant_c,
):
    with open(path, "w") as f:
        f.write(
            "    _lattice_constants = _a, _b, _c = (\n"
            f"       np.array([{lattice_constant_a}, 0., 0.]),\n"
            f"       np.array([{-lattice_constant_a*.5}, "
            f"{-lattice_constant_b * np.sqrt(3) / 2}, 0.]),\n"
            f"       np.array([0., 0., {lattice_constant_c}])\n    "
            ")\n\n"
        )


def get_vertex_prototypes(cif_data):

    nonlinear_vertices = {}
    linear_vertices = {}
    for i, ctype in enumerate(cif_data["coord_type"]):
        coordinates = cif_data["coords"][i]
        if ctype == "Ti":
            vertex = NewNonLinearVertex(
                id=len(nonlinear_vertices),
                position=(
                    coordinates[0] - coordinates[1] * 0.5,
                    -coordinates[1] * np.sqrt(3) / 2,
                    coordinates[2],
                ),
            )
            nonlinear_vertices[i] = vertex
        else:
            vertex = stk.cof.LinearVertex(
                id=len(linear_vertices),
                position=(
                    coordinates[0] - coordinates[1] * 0.5,
                    -coordinates[1] * np.sqrt(3) / 2,
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


def get_bb_dictionary(linear_vertices, nonlinear_vertices, coord_type):
    # Define atomistic bbs.
    bb1 = stk.BuildingBlock(
        smiles="Br/N=C/c2ccc(c1ccc(/C=N\Br)cc1)cc2",
        functional_groups=(stk.BromoFactory(),),
    )
    bb2 = stk.BuildingBlock(
        smiles="Br/N=C/c1ccc(/C=N/Br)cc1",
        functional_groups=(stk.BromoFactory(),),
    )
    bb3 = stk.BuildingBlock(
        smiles="Brc4ccc(c3cc(c1ccc(Br)cc1)cc(c2ccc(Br)cc2)c3)cc4",
        functional_groups=(stk.BromoFactory(),),
    )

    bb1_values = []
    bb2_values = []
    for i in linear_vertices:
        vert = linear_vertices[i]
        id_ = vert.get_id()
        if coord_type[i] == "Mn":
            bb1_values.append(id_)
        elif coord_type[i] == "Pb":
            bb2_values.append(id_)

    bb3_values = []
    for i in nonlinear_vertices:
        vert = nonlinear_vertices[i]
        bb3_values.append(vert.get_id())

    return {
        bb1: tuple(bb1_values),
        bb2: tuple(bb2_values),
        bb3: tuple(bb3_values),
    }


def run_construction(
    linear_vertices,
    nonlinear_vertices,
    edge_prototypes,
    cif_data,
):
    coord_type = cif_data["coord_type"]
    cell_lengths = cif_data["cell_lengths"]
    bb_dict = get_bb_dictionary(
        linear_vertices=linear_vertices,
        nonlinear_vertices=nonlinear_vertices,
        coord_type=coord_type,
    )

    lattice_size = (1, 1, 1)

    vertex_prototypes = tuple(
        nonlinear_vertices[i] for i in nonlinear_vertices
    ) + tuple(linear_vertices[i] for i in linear_vertices)

    topology_graph = CustomPeriodicTopology(
        building_blocks=bb_dict,
        lattice_size=lattice_size,
        lattice_constant_a=cell_lengths[0],
        lattice_constant_b=cell_lengths[1],
        lattice_constant_c=5.0,
        vertex_prototypes=vertex_prototypes,
        edge_prototypes=edge_prototypes,
    )

    construction_result = topology_graph.construct()
    cof = stk.ConstructedMolecule.init_from_construction_result(
        construction_result=construction_result,
    )
    periodic_info = construction_result.get_periodic_info()
    unit_cell = stko.UnitCell(
        vector_1=periodic_info.get_vector_1(),
        vector_2=periodic_info.get_vector_2(),
        vector_3=periodic_info.get_vector_3(),
    )

    return cof, unit_cell


def get_periodicity(
    nlc,
    lc,
    lattice_constant_a,
    lattice_constant_b,
    lattice_constant_c,
):
    vec1 = np.zeros((2))
    vec1[0] = (nlc[0] - nlc[1] * 0.5) - (lc[0] - lc[1] * 0.5)
    vec1[1] = -(nlc[1] - lc[1]) * np.sqrt(3) / 2

    # +x
    vec2 = np.zeros((2))
    vec2[0] = (nlc[0] + lattice_constant_a - nlc[1] / 2) - (
        lc[0] - lc[1] / 2
    )
    vec2[1] = -(nlc[1] - lc[1]) * np.sqrt(3) / 2

    # +x+y
    vec3 = np.zeros((2))
    vec3[0] = (
        nlc[0] + lattice_constant_a - (nlc[1] - lattice_constant_b) / 2
    ) - (lc[0] - lc[1] / 2)
    vec3[1] = -(nlc[1] - lattice_constant_b - lc[1]) * np.sqrt(3) / 2

    # +y
    vec4 = np.zeros((2))
    vec4[0] = (nlc[0] - (nlc[1] - lattice_constant_b) / 2) - (
        lc[0] - lc[1] / 2
    )
    vec4[1] = -(nlc[1] - lattice_constant_b - lc[1]) * np.sqrt(3) / 2

    # -x
    vec5 = np.zeros((2))
    vec5[0] = (nlc[0] - nlc[1] / 2) - (
        lc[0] + lattice_constant_a - lc[1] / 2
    )
    vec5[1] = -(nlc[1] - lc[1]) * np.sqrt(3) / 2

    # -x-y
    vec6 = np.zeros((2))
    vec6[0] = (nlc[0] - (nlc[1]) / 2) - (
        lc[0] + lattice_constant_a - (lc[1] - lattice_constant_b) / 2
    )
    vec6[1] = -(nlc[1] - (lc[1] - lattice_constant_b)) * np.sqrt(3) / 2

    # -y
    vec7 = np.zeros((2))
    vec7[0] = (nlc[0] - nlc[1] / 2) - (
        lc[0] - (lc[1] - lattice_constant_b) / 2
    )
    vec7[1] = -(nlc[1] - (lc[1] - lattice_constant_b)) * np.sqrt(3) / 2

    # x-y
    vec8 = np.zeros((2))
    vec8[0] = (nlc[0] + lattice_constant_a - nlc[1] / 2) - (
        lc[0] - (lc[1] - lattice_constant_b) / 2
    )
    vec8[1] = -(nlc[1] - (lc[1] - lattice_constant_b)) * np.sqrt(3) / 2

    # -x+y
    vec9 = np.zeros((2))
    vec9[0] = (nlc[0] - (nlc[1] - lattice_constant_b) / 2) - (
        lc[0] + lattice_constant_a - (lc[1]) / 2
    )
    vec9[1] = -(nlc[1] - lattice_constant_b - (lc[1])) * np.sqrt(3) / 2

    vec_mag_1 = np.sqrt(vec1[1] ** 2 + vec1[0] ** 2)
    vec_mag_2 = np.sqrt(vec2[1] ** 2 + vec2[0] ** 2)
    vec_mag_3 = np.sqrt(vec3[1] ** 2 + vec3[0] ** 2)
    vec_mag_4 = np.sqrt(vec4[1] ** 2 + vec4[0] ** 2)
    vec_mag_5 = np.sqrt(vec5[1] ** 2 + vec5[0] ** 2)
    vec_mag_6 = np.sqrt(vec6[1] ** 2 + vec6[0] ** 2)
    vec_mag_7 = np.sqrt(vec7[1] ** 2 + vec7[0] ** 2)
    vec_mag_8 = np.sqrt(vec8[1] ** 2 + vec8[0] ** 2)
    vec_mag_9 = np.sqrt(vec9[1] ** 2 + vec9[0] ** 2)

    # print(
    #     vec_mag_1,
    #     vec_mag_2,
    #     vec_mag_3,
    #     vec_mag_4,
    #     vec_mag_5,
    #     vec_mag_6,
    #     vec_mag_7,
    #     vec_mag_8,
    #     vec_mag_9,
    # )

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

    logging.info("need to check why I do not get none in vecmag")
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
    if not len(sys.argv) == 2:
        logging.info(f"Usage: {__file__}\n" "   Expected 1 arguments:")
        logging.info(
            "cif_file (str): `.cif` file to python script from"
        )
        sys.exit()
    else:
        cif_file = sys.argv[1]

    name = cif_file.replace(".cif", "")
    _, temp, n = name.split("_")
    logging.info(f"ask about this variable, n: {n}, temp: {temp}")
    cif_data = parse_cif(cif_file, n=int(n))

    write_output_txt(
        path=f"output_{temp}.txt",
        lattice_constant_a=cif_data["cell_lengths"][0],
        lattice_constant_b=cif_data["cell_lengths"][1],
        lattice_constant_c=cif_data["cell_lengths"][2],
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
    )

    CifWriter().write(
        molecule=cof,
        path=f"cof_n_{n}_T_{temp}.cif",
        periodic_info=unit_cell,
    )
    stk.XyzWriter().write(molecule=cof, path="cof.xyz")
    stk.MolWriter().write(molecule=cof, path="cof.mol")

    logging.info('optimising the cof with Gulp...')
    gulp_opt = stko.GulpUFFOptimizer(
        gulp_path="/home/atarzia/software/gulp-6.1/Src/gulp",
        output_dir=f'{name}_gulpopt',
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
        path=f"cof_n_{n}_T_{temp}_opt.cif",
        periodic_info=unit_cell,
    )
    stk.XyzWriter().write(molecule=cof, path="cof_opt.xyz")
    stk.MolWriter().write(molecule=cof, path="cof_opt.mol")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
