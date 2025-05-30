# `PyBigDFT/BigDFT/Atom.py`

## Overview

The `Atom.py` module provides the `Atom` class, which serves as a fundamental data structure for representing individual atoms within the PyBigDFT ecosystem. It allows for storing various atomic properties such as symbol, position, and multipole moments, and handles necessary unit conversions, primarily between Bohr and Angstroem for coordinates.

## Key Components

### `nzion(sym: str) -> float`
A utility function that returns the charge of the nucleus associated with a given atomic symbol.
-   `sym`: The atomic symbol (e.g., "H", "He").
-   Returns: The nuclear charge as a float.

### `class Atom(MutableMapping)`
A dictionary-like class to store and manage data for a single atom. It allows dynamic association of various properties while providing standardized methods for common operations.

#### Key Methods and Properties:

-   **`__init__(self, *args, **kwargs)`**: Initializes an `Atom` object. Can accept a dictionary or keyword arguments to populate atomic data.
-   **`dict(self) -> dict`**: Returns the internal dictionary store of the atom.
-   **`get_external_potential(self, units="bohr") -> dict`**: Transforms the atom's data into a dictionary format suitable for defining an external potential in BigDFT calculations.
    -   `units`: Desired units for the position in the output dictionary (default: "bohr").
-   **`q0` (property)**: Returns the monopole (charge) of the atom, if defined.
-   **`q1` (property)**: Returns the dipole of the atom, if defined (converted to a standard orientation).
-   **`set_multipole(self, mp: dict or Atom)`**: Sets multipole-related values (charge, dipole, quadrupole, sigma) for this atom from another atom or dictionary.
-   **`sym` (property)**: Returns the atomic symbol (e.g., "H", "C").
-   **`get_position(self, units="bohr") -> list[float]`**: Retrieves the atom's Cartesian coordinates in the specified units.
    -   `units`: Desired units ("bohr" or "angstroem").
-   **`set_position(self, new_pos: list[float], units="bohr")`**: Updates the atom's position.
    -   `new_pos`: A list of three floats representing the new x, y, z coordinates.
    -   `units`: The units of `new_pos` (default: "bohr").
-   **`__eq__(self, other: dict or Atom) -> bool`**: Compares two atoms for equality based on their symbol and position (within a small tolerance).

## Important Variables/Constants

-   **`AU_to_A: float`**: Conversion factor from Atomic Units (Bohr) to Angstroem. Value: `0.52917721092`.
-   **`MULTIPOLE_ANALYSIS_KEYS: list[str]`**: A list of keys used for multipole analysis data: `['q0', 'q1', 'q2', 'sigma', 'multipole character']`.
-   **`PROTECTED_KEYS: list[str]`**: A list of keys that have special meaning or are protected from general deletion/modification: `MULTIPOLE_ANALYSIS_KEYS + ["frag"] + ["r", "units"]`.

## Usage Examples

```python
from PyBigDFT.BigDFT import Atom # Assuming Atom class is directly available under this import
from futile.Utils import write as safe_print # Used in the original example

# Create an Atom object using dictionary initialization
test_atom_dict = Atom.Atom({'r': [1.0, 0.0, 0.0], 'sym': "He", 'units': 'bohr'})
safe_print("Atom from dict:", test_atom_dict.dict())
safe_print("Symbol:", test_atom_dict.sym)
safe_print("Position (bohr):", test_atom_dict.get_position())
safe_print("Position (angstroem):", test_atom_dict.get_position('angstroem'))
safe_print("-" * 20)

# Create an Atom object using keyword arguments (common alternative style)
test_atom_kwargs = Atom.Atom(He=[1.0, 0.0, 0.0], units='bohr')
safe_print("Atom from kwargs:", test_atom_kwargs.dict())
safe_print("Symbol:", test_atom_kwargs.sym)
safe_print("Position (bohr):", test_atom_kwargs.get_position())
safe_print("-" * 20)

# Compare two atoms
# Create a new atom with position in Angstroem that matches test_atom_dict
pos_angstroem = test_atom_dict.get_position('angstroem')
new_atom_ang = Atom.Atom({'r': pos_angstroem, 'sym': "He", 'units': 'angstroem'})
safe_print(f"Are test_atom_dict and new_atom_ang equal? {new_atom_ang == test_atom_dict}")
safe_print("-" * 20)

# Add custom properties (Atom behaves like a dictionary)
test_atom_dict["custom_frag"] = "MoleculeFragment1"
test_atom_dict["energy"] = -79.0
safe_print("Atom with custom properties:", test_atom_dict.dict())
safe_print(f"Custom fragment: {test_atom_dict['custom_frag']}")
safe_print("-" * 20)

# Set a new position
safe_print(f"Original position (bohr): {test_atom_kwargs.get_position()}")
test_atom_kwargs.set_position([0.5, 0.5, 0.5], units='angstroem')
safe_print(f"New position (bohr): {test_atom_kwargs.get_position()}")
safe_print(f"New position (angstroem): {test_atom_kwargs.get_position(units='angstroem')}")
safe_print("-" * 20)

# Example of changing symbol (if 'sym' and 'r' representation is used)
atom_sym_r = Atom.Atom({'r': [0.0, 0.0, 0.0], 'sym': "C", 'units': 'bohr'})
safe_print("Original atom (sym/r):", atom_sym_r.dict())
atom_sym_r["sym"] = "N" # Change symbol
safe_print("Atom with changed symbol:", atom_sym_r.dict())
safe_print(f"New symbol: {atom_sym_r.sym}")
```

## Dependencies and Interactions

-   **`collections.abc.MutableMapping`**: The `Atom` class inherits from `MutableMapping` to provide dictionary-like behavior.
-   **`numpy`**: Used internally for array operations, especially for positions and multipole moments (though not explicitly imported in every method, its usage is implied for array-like objects).
-   **`futile.Utils.write`**: Used for safe printing in the module's `if __name__ == "__main__":` block (demonstrated in usage examples).
```
