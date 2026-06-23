# Crystallography

Python implementation of fundamental concepts in crystallography.
Designed to be a module for various other projects.

## What this library does

Its main functionality is to define Crystal = Lattice + SpaceGroup + Atoms
and investigate their properties, while taking symmetry constraints very seriously. 
It implements symmetry operations, space groups and magnetic space groups as well as atoms and lattices.


## Module map

```
bragg/
├── crystal/
│   ├── atom.py          # Atom: position, moment, spin, g-tensor
│   ├── lattice.py       # Lattice, LatticeOriented: coordinate transforms
│   └── crystal.py       # Crystal(Lattice): atoms + MSG, generates full unit cell
├── symmetry/
│   ├── group.py             # SymOp (ABC), Group[T] (generic, Cayley BFS)
│   ├── magnetic_symmetry.py # mSymOp(SymOp), MSG(Group[mSymOp])
│   └── crystall_space_group.py  # cSymOp(SymOp), SG(Group[cSymOp])
├── plotting/
│   ├── supercell_plotter.py     # SupercellPlotter (ABC): balls, lines, arrows, ellipsoids
│   ├── supercell_plotter_mpl.py # MPLSupercellPlotter: matplotlib backend
│   └── crystal_plotter_mixin.py # CrystalPlotterMixin: stub for Crystal plotting methods
├── databases/
│   ├── database.py          # db_entry (ABC), Database, load_database_from_txt
│   ├── implement_txt.py     # atom_data, isotope_data, color_data, magion_data, xrayion_data
│   └── implement_spglib.py  # SG_entry, MSG_entry (TODO: complete spglib integration)
├── utils/
│   ├── linalg.py            # Rotation matrices, DMI, dipolar, vector ops
│   └── arrays.py            # @ensure_shape decorator, create_mesh
└── data_tables/             # Raw database text files
    ├── atom.dat, isotope.dat, magion.dat, color.dat, xrayion.dat
    ├── spglib_SGnames.dat, spglib_MSGnames.dat
    └── SG_generators.txt    # TODO: generate from spglib; MSG_generators.txt planned
```


## Core class hierarchy

```
SymOp (ABC)  <- symmetry/group.py
├── cSymOp   <- symmetry/crystall_space_group.py  (matrix + translation)
└── mSymOp   <- symmetry/magnetic_symmetry.py     (matrix + translation + time_reversal)

Group[T: SymOp]  <- symmetry/group.py  (generic, Cayley-BFS construction)
├── SG   = Group[cSymOp]   <- crystall_space_group.py
└── MSG  = Group[mSymOp]   <- magnetic_symmetry.py

Lattice              <- crystal/lattice.py
└── LatticeOriented  <- crystal/lattice.py  (Lattice + orientation matrix U)

Crystal(Lattice)     <- crystal/crystal.py  (atoms + MSG)

SupercellPlotter (ABC)       <- plotting/supercell_plotter.py
└── MPLSupercellPlotter      <- plotting/supercell_plotter_mpl.py

db_entry (ABC, dataclass)  <- databases/database.py
└── atom_entry, isotope_entry, magion_entry, color_entry, xrayion_entry  <- implement_txt.py
└── SG_entry, MSG_entry    <- implement_spglib.py (TODO: complete)

Database  <- databases/database.py
```


# License
Please cite:
https://github.com/mstekiel/crystallography

# DEV

## Bugs

 - [ ] There is some warning in `database` module that causes a lot of empty lines printout when warning levels are higher.

## TODO

- Core
  - [ ] Implement magnetic field
    - [x] implement
    - [ ] test -> sign convention unclear
    - [ ] SpinW -> only obscuretutorials found
  - [ ] Implement primitive/reduced cell such that the calculation is on smaller cell=matrix=faster, while the setup is in the coordinates of the main cell, which are easier to interpret
  - [ ] standard paths of lattices doi.org/10.1016/j.commatsci.2010.05.010
  - [x] MSG should inherit from Group
  - [x] crystallographic SG should also inherit from Group

- Future
  - [ ] Lot of core functionalities rely on hashing the objects for unique identifiers:
    -  `SymOps` do it from string for unique elements finding
    -  `Atom` same, for unique position finding
    Is there a sturdy way to implement hashing?

- Usage
  - [ ] Plotting of unit cell edges.
  - [x] coordinate systems.
  - [x] lighting.

- Make documentation:
  For list of modules and their descriptions see documentation at: 
  https://mstekiel.github.io/mikibox/build/html/index.html

- Make GUI (PyQt6) with inline editeor (https://qscintilla.com), 3D viewer (https://vispy.org). The GUI can be a direct copy of vspy/.examples/../sandbox
  - VISPY has a lot of beautiful examples for Ctrl+CV coding.
    - whole GUI as in vispy\examples\basics\scene\modular_shaders\sandbox.py
    - tubes vispy\examples\basics\visuals\tube.py
    - check out `Canvas.measure_fps()`
    - demo\gloo\boids.py for animatoins with particles (neutrons)
    - demo\gloo\camera.py with gestures recognition neural network to control the viewing.
  - https://mediapipe-studio.webapps.google.com/studio/demo/face_landmarker for gesture recognition and 3D navigation


## Managing project with uv

Moved the project to uv. 

Now, if you want to use `spinwaves` as library in editable mode include these in the `pyproject.toml`:

```toml
[project]
dependencies = [
    "spinwaves",
]

[tool.uv.sources]
spinwaves = { path = "path_to_spinwaves_project", editable=true }
```
### Testing
Full testing

```python
# run pytest
uv run pytest
 
# generate htmlcov/ report, open htmlcov/index.html
uv run pytest --cov=bragg --cov-report=html

# generate htmlcov/ report, open htmlcov/index.html
uv run pytest --cov=bragg --cov-report=term-missing --cov-fail-under=80

# full matrix — slow, run deliberately
uvx nox -f tests/noxfile.py -s version_matrix
```

## Test structure

Tests mirror the bragg library file structure:

```
tests/
├── test_crystal.py              # Crystal construction, atom indexing, DMI validation
├── test_linalg.py               # Rotation matrices (Rx/Ry/Rz), RtoZ, Rodrigues
├── test_utils_arrays.py         # create_mesh, ensure_shape
├── symmetry/
│   ├── __init__.py              # make_perm_symop(n) factory + S3/C4 generators
│   ├── test_group.py            # Group axioms, Cayley graph, known orders
│   ├── test_magnetic_symmetry.py
│   └── test_symmetrization.py
└── databases/
    └── test_txt_databases.py    # All 5 txt databases: types, values, search, iteration
```


## General notes
- I tried to follow the structure of https://github.com/pypa/sampleproject for the development of this project. Following descriptions from https://python-packaging-tutorial.readthedocs.io/en/latest/setup_py.html