# `PyBigDFT/BigDFT/Fragments.py`

## Overview

The `Fragments.py` module provides classes and functions for working with molecular fragments and systems composed of multiple fragments. It is central to BigDFT's capabilities for fragment-based calculations, allowing for the definition, manipulation, and analysis of subsystems. Key functionalities include reading/writing XYZ files for fragments, performing geometric transformations (rotations, translations), and handling multipole data associated with fragments.

## Key Components

### `XYZfile` (class)
Manages the creation and writing of XYZ files.
-   **`__init__(self, filename: str = None, units: str = 'atomic')`**: Initializes an XYZ file writer.
    -   `filename`: Name of the file to write. If `None`, output goes to `sys.stdout`.
    -   `units`: Units for coordinates ('atomic' or 'angstroem').
-   **`append(self, array: list[list[float]], basename: str = '', names: list[str] = None, attributes: list[dict] = None)`**: Adds atomic coordinates to the file.
    - `array`: List of atomic positions.
    - `basename`: Base string for atom names.
    - `names`: List of specific atom names/symbols.
    - `attributes`: List of dictionaries containing additional data per atom to be written in the comment section of its line.
-   **`dump(self, position: str = 'w')`**: Writes the accumulated atomic data to the file or standard output.

### `Lattice` (class)
Defines lattice vectors for periodic systems.
-   **`__init__(self, vectors: list[list[float]])`**: Initializes with a list of lattice vectors.
-   **`grid(self, origin: list[float] = [0.0, 0.0, 0.0], extremes: list[list[int]] = None, radius: float = None) -> list[list[float]]`**: Generates a grid of translation vectors based on the lattice vectors.
    - `origin`: The starting point for the grid.
    - `extremes`: List of `[min_i, max_i]` for each lattice vector, defining the range of integer multiples.
    - `radius`: If provided, only includes translation vectors whose norm is less than this radius.

### `RotoTranslation` (class)
Represents a rigid-body rototranslation.
-   **`__init__(self, pos1: numpy.ndarray, pos2: numpy.ndarray)`**: Calculates the optimal rotation (`R`) and translation (`t`) to align `pos1` (a set of points) onto `pos2` using `wahba.rigid_transform_3D`. The alignment error is stored in `J`.
-   **`dot(self, pos: numpy.ndarray) -> numpy.ndarray`**: Applies the stored rototranslation (R*pos + t) to a set of positions.
-   **`invert(self)`**: Inverts the transformation (R becomes R.T, t becomes -R.T*t).

#### `Translation(RotoTranslation)` (subclass)
Specialized `RotoTranslation` for pure translations.
-   **`__init__(self, t: list[float])`**: Initializes with a translation vector `t`. `R` is None.

#### `Rotation(RotoTranslation)` (subclass)
Specialized `RotoTranslation` for pure rotations.
-   **`__init__(self, R: numpy.ndarray)`**: Initializes with a rotation matrix `R`. `t` is None.

### `Fragment` (class)
Represents a molecular fragment, a collection of atoms with associated properties.
-   **`__init__(self, atomlist: list[dict] = None, id: str = 'Unknown', units: str = 'AU')`**:
    -   `atomlist`: A list of atom dictionaries (BigDFT format, e.g., `[{'C': [0,0,0]}, {'H': [1,0,0]}]`).
    -   `id`: A string identifier for the fragment.
    -   `units`: Units for input coordinates ('AU' or 'A'). Internal storage is in AU.
-   **`append(self, atom: dict = None, sym: str = None, positions: list[float] = None)`**: Adds an atom to the fragment.
-   **`xyz(self, filename: str = None, units: str = 'atomic')`**: Writes the fragment to an XYZ file. `units` here refers to output units.
-   **`get_posinp(self, units: str = "AU") -> dict`**: Returns a dictionary suitable for `posinp` in a BigDFT input file (e.g., `{'positions': [{'C': [0,0,0]}...], 'units': 'AU'}`).
-   **`dict(self) -> list[dict]`**: Returns a list of atom dictionaries, potentially including multipole information, suitable for external potential definitions. Coordinates are in AU.
-   **`centroid(self) -> numpy.ndarray`**: Calculates the geometric center of the fragment (in AU).
-   **`transform(self, Rt: RotoTranslation)`**: Applies a `RotoTranslation` to the fragment's atoms. Also transforms atomic dipoles if present.
-   **`line_up(self)`**: Aligns the fragment's principal axes of inertia with the coordinate axes and centers its centroid at the origin.
-   **`Q(self, atom: list = None, order: int = 0) -> float or numpy.ndarray or None`**:
    - If `set_fragment_multipoles` was called, returns the pre-stored fragment monopole (order 0), dipole (order 1), or quadrupole (order 2).
    - If `order=0` and fragment multipoles are not set, it computes the sum of atomic charges (`q0` from atom dictionaries). `atom` argument is deprecated.
-   **`set_fragment_multipoles(self, q0, q1, q2)`**: Stores the pre-calculated multipoles (monopole, dipole, quadrupole) for the entire fragment.
-   **`d0(self, center: list[float] = None) -> numpy.ndarray or None`**: Calculates the dipole moment from atomic charges (`q0`) relative to `center` (or fragment centroid if `center` is None). Result in AU.
-   **`d1(self, center: list[float] = None) -> numpy.ndarray or None`**: Calculates the total dipole moment including atomic dipoles (`q1`) and charge-based dipole (`d0`). Result in AU.

### `System` (class)
Represents a collection of `Fragment` objects.
-   **`__init__(self, mp_dict: list[dict] = None, xyz: str = None, nat_reference: int = None, units: str = 'AU', transformations: list[dict] = None, reference_fragments: list[Fragment] = None, posinp_dict: dict = None, frag_partition: list[int] = None, fragmentation: list[list] = None)`**:
    Initializes a system from various sources.
    -   `mp_dict`: List of atom dictionaries with multipole information.
    -   `xyz`: Path to an XYZ file. `nat_reference` or `fragmentation` is needed to partition atoms into fragments.
    -   `nat_reference`: If `xyz` is given, number of atoms per fragment for uniform fragmentation.
    -   `units`: Default units for data read from `xyz` or `posinp_dict`.
    -   `transformations`, `reference_fragments`: Used by `recompose`.
    -   `posinp_dict`: A BigDFT `posinp` dictionary.
    -   `frag_partition`: List of indices to partition a flat list of atoms (e.g., from `mp_dict`).
    -   `fragmentation`: A list like `[['label1', [atom_indices_1_based]], ['label2', [atom_indices_1_based]]]` to define fragments from a flat list of atoms (e.g., from an XYZ file).
-   **`append(self, frag: Fragment)`**: Adds a `Fragment` to the system.
-   **`fill_from_xyz(self, file: str, nat_reference: int, fragmentation: list)`**: Populates the system by reading an XYZ file and partitioning it into fragments.
-   **`fill_from_mp_dict(self, mpd: list[dict], nat_reference: int = None, frag_partition: list[int] = None)`**: Populates from a list of atom dictionaries with multipole information.
-   **`fill_from_posinp_dict(self, dct: dict)`**: Populates from a BigDFT `posinp` dictionary, respecting `frag` field in atom entries.
-   **`xyz(self, filename: str = None, units: str = 'atomic')`**: Writes the entire system to an XYZ file, with fragment information (name and index) in comments for each atom.
-   **`write_fragfile(self, filename: str = "log.yaml")`**: Creates a fragment list file (e.g., `frag.yaml`) defining the 1-based atom indices for each fragment in the system, suitable for BigDFT fragment calculations.
-   **`dict(self, filename: str = None) -> dict`**: Returns a dictionary representation of the system (e.g., `{'units': 'AU', 'values': [atom_dicts...], 'global monopole': ...}`), suitable for some BigDFT inputs or analysis.
-   **`decompose(self, reference_fragments: list[Fragment])`**: Decomposes the system's fragments by finding the best rototranslation to match them against a list of `reference_fragments`. Stores transformations in `self.decomposition`.
-   **`recompose(self, transformations: list[dict] = None, reference_fragments: list[Fragment] = None)`**: Rebuilds the system from a list of transformations and reference fragments.
-   **`Q(self) -> float or None`**: Sum of monopoles of all fragments that have `Q()` defined (either from `set_fragment_multipoles` or summed atomic charges).

### Utility Functions
-   **`GetSymbol(atom: dict) -> str`**: Extracts the chemical symbol from an atom dictionary (e.g., if key is 'C', returns 'C'; if key is 'r' and 'sym' exists, returns `atom['sym']`).
-   **`SetFragId(name: str, fragid: int) -> str`**: Creates a fragment identifier string (e.g., "MOL:1").
-   **`GetFragTuple(id: str) -> tuple[str, int]`**: Parses a fragment identifier string (e.g., "MOL:1") into a name and ID tuple (e.g., `('MOL', 1)`).
-   **`fragmentation_dict(pattern: list[int], repeat: int, labels: list[str] = None) -> list`**: Creates a `fragmentation` list (for `System` class) based on a repeating pattern of fragment sizes.
-   **`CreateFragDict(start: dict) -> dict`**: Processes a `posinp` dictionary (with `frag` field in atoms) to group atoms by fragment name and ID, returning a nested dictionary.
-   **`MergeFragmentsTogether(frag_dict: dict, merge_list: list[list[tuple]]) -> tuple[dict, list]`**: Merges specified fragments within a `frag_dict` (as created by `CreateFragDict`).
-   **`CreateFragList(frag_dict: dict, merge_list: list = None) -> list`**: Creates a flat list of fragments (tuples of `(fragment_label, atom_indices_1_based)`) from a `frag_dict`.
-   **`prepare_fragment_inputs(...)`**: Helper to set up input files for fragment calculations (not fully detailed here, involves creating YAML and XYZ files).
-   **`find_reference_fragment(...)`**: Identifies which fragments in a list correspond to a set of reference templates by object identity.
-   **`frag_average(...)`**: Averages attributes (like multipoles `q0`, `q1`, `q2` and their sigmas) over a list of similar fragments, using the first fragment as a positional reference.
-   **`distance(i: Fragment, j: Fragment, cell: list[float] = None) -> float`**: Calculates the distance between the centroids of two fragments, considering periodic boundary conditions if `cell` (list of 3 cell lengths) is provided.

## Important Variables/Constants

-   **`AU_to_A: float`**: Conversion factor from Bohr (Atomic Units of length) to Angstroem. Value: `0.52917721092`.
-   **`Debye_to_AU: float`**: Conversion factor from Debye to Atomic Units for dipoles. Value: `0.393430307`.
-   **`MULTIPOLE_ANALYSIS_KEYS: list[str]`**: Keys relevant for multipole data: `['q0', 'q1', 'q2', 'sigma']`.
-   **`PROTECTED_KEYS: list[str]`**: Keys with special meaning in atom dictionaries that are not element symbols: `MULTIPOLE_ANALYSIS_KEYS + ["frag"]`.

## Usage Examples

```python
from PyBigDFT.BigDFT import Fragments
import numpy as np
import os # For file operations

# 1. Create a Fragment
atom_c = {"C": [0.0, 0.0, 0.0]}
atom_h1 = {"H": [0.0, 0.0, 1.08]} # Approx C-H bond in Angstroem
atom_h2 = {"H": [1.02, 0.0, -0.36]}
atom_h3 = {"H": [-0.51, -0.88, -0.36]}
atom_h4 = {"H": [-0.51, 0.88, -0.36]}
methane_frag = Fragments.Fragment([atom_c, atom_h1, atom_h2, atom_h3, atom_h4],
                                  id="Methane_Mol", units="A")

# Get fragment centroid (result in AU, convert to Angstroem for printing)
print(f"Methane centroid (Angstroem): {methane_frag.centroid() * Fragments.AU_to_A}")

# Write fragment to XYZ file
methane_frag.xyz("methane.xyz", units="angstroem")
print("Wrote methane.xyz")

# 2. Create a System from two fragments
water1_atoms = [{"O": [0.0, 0.0, 0.0]}, {"H": [0.96, 0.0, 0.0]}, {"H": [-0.24, 0.93, 0.0]}]
water1 = Fragments.Fragment(water1_atoms, id="Water1", units="A")

water2_atoms = [{"O": [3.0, 0.0, 0.0]}, {"H": [2.04, 0.0, 0.0]}, {"H": [3.24, 0.93, 0.0]}]
water2 = Fragments.Fragment(water2_atoms, id="Water2", units="A")

water_system = Fragments.System(units="A") # System units also default to AU if not specified
water_system.append(water1)
water_system.append(water2)

# Write system to XYZ (includes fragment info in comments)
water_system.xyz("water_dimer.xyz", units="angstroem")
print("Wrote water_dimer.xyz")

# Write a fragment definition file (for BigDFT fragment calculations)
# This defines which atoms belong to which fragment.
water_system.write_fragfile("water_dimer.frags.yaml")
print("Wrote water_dimer.frags.yaml")

# 3. Using RotoTranslation
# Create a translation vector (in Atomic Units)
trans_vec_au = np.array([1.0, 0.5, 0.0])
translation = Fragments.Translation(trans_vec_au)

# Apply to water2 (its internal coordinates are in AU)
print(f"Water2 original centroid (AU): {water2.centroid()}")
water2.transform(translation)
print(f"Water2 translated centroid (AU): {water2.centroid()}")

# 4. Example of reading a system from XYZ and defining fragmentation
# Create a dummy XYZ file first for this example
xyz_content = """4 angstroem
Two H2 molecules, coordinates in Angstroem
H 0.0 0.0 0.0
H 0.0 0.0 0.7
H 2.0 0.0 0.0
H 2.0 0.0 0.7
"""
with open("dummy_system.xyz", "w") as f:
    f.write(xyz_content)

# Define fragmentation: two fragments, each with 2 atoms.
# Atom indices in the 'fragmentation' list are 1-based.
frag_definition = [["H2_frag_A", [1, 2]], ["H2_frag_B", [3, 4]]]

system_from_xyz = Fragments.System(xyz="dummy_system.xyz",
                                   fragmentation=frag_definition,
                                   units="A") # Specify units of the XYZ file
print(f"\nNumber of fragments in system_from_xyz: {len(system_from_xyz.fragments)}")
for frag in system_from_xyz.fragments:
    print(f"Fragment ID: {frag.id}, Number of atoms: {len(frag)}")
    # Verify centroid calculation (coordinates are converted to AU internally)
    print(f"  Centroid (AU): {frag.centroid()}")

# Clean up example files
files_to_remove = ["methane.xyz", "water_dimer.xyz",
                   "water_dimer.frags.yaml", "dummy_system.xyz"]
for f in files_to_remove:
    if os.path.exists(f):
        os.remove(f)
        print(f"Removed {f}")

```

## Dependencies and Interactions

-   **`numpy`**: Extensively used for all numerical operations, especially coordinate manipulations, matrix operations for transformations.
-   **`futile.Utils.write`**: Used for safe printing.
-   **`yaml`**: Used for writing fragment definition files (implicitly via `write_fragfile` which calls `yaml.dump`).
-   **`BigDFT.wahba`**: The `RotoTranslation` class uses `wahba.rigid_transform_3D` to find optimal alignments between sets of points.
-   This module is fundamental for any calculations involving the "fragment" capabilities of BigDFT, including linear scaling calculations (where fragments might represent localized orbitals or basis functions) and analyses based on molecular decomposition (e.g., QMMM, subsystem DFT).
```
