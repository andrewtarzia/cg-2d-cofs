#!/usr/bin/env python

# -*- coding: utf-8 -*-

# Distributed under the terms of the MIT License.

"""
Module for calculating Oh shape.

Author: Andrew Tarzia

"""

import subprocess as sp
import os
import logging
import pathlib


def shape_path():

    path = pathlib.Path(
        "/home/atarzia/software/shape_2.1_linux_64/SHAPE_2.1_linux_64/"
        "shape_2.1_linux64"
    )
    if not os.path.exists(path):
        raise FileNotFoundError(
            "you need to install SHAPE 2.1 or update this function in "
            f"{__file__} with the right directory"
        )

    return path


def calculate_hp6(molecule, name):
    return calculate_shape(molecule, name)["HP-6"]


def calculate_shape(molecule, name):
    """
    Calculate the shape of a 6 atom molecule.

    Shape: http://www.ee.ub.edu/index.php?option=com_content&view=
    article&id=575:shape-available&catid=80:news&Itemid=466

    """

    if molecule.get_num_atoms() != 6:
        raise ValueError("Molecule does not have 6 atoms.")

    shape_dicts = (ref_shape_dict()["hp6"],)
    n_verts = list(set([i["vertices"] for i in shape_dicts]))
    if len(n_verts) != 1:
        raise ValueError("Different vertex shapes selected.")

    input_file = "temp_shp.dat"
    std_out = "temp_shp.out"
    output_file = "temp_shp.tab"
    write_shape_input_file(
        input_file=input_file,
        name=name,
        structure=molecule,
        num_vertices=n_verts[0],
        central_atom_id=0,
        ref_shapes=[i["code"] for i in shape_dicts],
    )

    run_shape(input_file, shape_path(), std_out)
    shapes = collect_all_shape_values(output_file)
    return shapes


def ref_shape_dict():
    return {
        "hp6": {
            "vertices": "6",
            "label": "HP-6",
            "code": "1",
        },
    }


def write_shape_input_file(
    input_file,
    name,
    structure,
    num_vertices,
    central_atom_id,
    ref_shapes,
):
    """
    Write input file for shape.

    """

    title = "$shape run by Andrew Tarzia.\n"
    size_of_poly = f"{num_vertices} {central_atom_id}\n"
    codes = " ".join(ref_shapes) + "\n"

    structure_string = f"{name}\n"
    pos_mat = structure.get_position_matrix()
    for atom in structure.get_atoms():
        ele = atom.__class__.__name__
        x, y, z = pos_mat[atom.get_id()]
        structure_string += f"{ele} {x} {y} {z}\n"

    string = title + size_of_poly + codes + structure_string

    with open(input_file, "w") as f:
        f.write(string)


def run_shape(input_file, shape_path, std_out):
    """
    Run input file for shape.

    """

    cmd = f"{shape_path} {input_file}"
    logging.info(f"running {cmd}")
    with open(std_out, "w") as f:
        # Note that sp.call will hold the program until completion
        # of the calculation.
        sp.call(
            cmd,
            stdin=sp.PIPE,
            stdout=f,
            stderr=sp.PIPE,
            # Shell is required to run complex arguments.
            shell=True,
        )


def collect_all_shape_values(output_file):
    """
    Collect shape values from output.

    """

    with open(output_file, "r") as f:
        lines = f.readlines()

    label_idx_map = {}
    for line in reversed(lines):
        if "Structure" in line:
            line = [
                i.strip()
                for i in line.rstrip().split("]")[1].split(" ")
                if i.strip()
            ]
            for idx, symb in enumerate(line):
                label_idx_map[symb] = idx
            break
        line = [i.strip() for i in line.rstrip().split(",")]
        values = line

    shapes = {
        i: float(values[1 + label_idx_map[i]]) for i in label_idx_map
    }

    return shapes
