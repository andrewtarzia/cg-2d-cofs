# cg-2d-cofs
Scripts to convert CIFs to `stk` molecules for analysis

Associated with the publication at DOI:

The code can be used by downloading this repository and creating a new conda environment based on `environment.yml`.


`AAmanualpdb_shape_analysis.py`

    Takes an atomistic `.pdb` model, places Ti atoms at the centre of the tritopic building blocks, and performs shape calculations based on these atoms. Here, bonding between Ti atoms is defined based on hard-coded distance criteria.

`CGcif_shape_analysis.py`

    Takes a coarse-grained `.cif` model, finds the Ti atoms (centre of the tritopic building block), and performs the shape analysis on them. Here, bonding between Ti atoms is defined based on hard-coded distance criteria.

`CGcif_to_atomistic.py`

    Takes a `.cif` of the coarse-grained models and uses `stk` to place the atomistic building blocks. Optimisation of the atomistic model with GULP/UFF is optional using `stko`.

`plot_shape_map.py`:

    Plot the shapes and hexagons on xy plane.