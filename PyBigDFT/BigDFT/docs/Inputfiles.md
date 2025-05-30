# `PyBigDFT/BigDFT/Inputfiles.py`

## Overview

The `Inputfiles.py` module defines the `Inputfile` class, which is the primary object for representing and manipulating BigDFT input parameters in Python. BigDFT input files use the YAML syntax, which maps directly to Python dictionaries. The `Inputfile` class leverages this by inheriting from the standard Python `dict` class.

A key feature of the `Inputfile` class is its dynamic integration with the functions defined in the `BigDFT.InputActions` module. At runtime, when an `Inputfile` object is instantiated, the action functions from `InputActions` (like `set_xc`, `optimize_geometry`, etc.) are bound to the instance as methods. This provides a convenient and expressive object-oriented interface for building up a BigDFT input.

## Key Components

### `class Inputfile(dict)`
Represents a BigDFT input file. It behaves like a standard Python dictionary but with added methods for common BigDFT input configurations.

-   **`__init__(self, *args, **kwargs)`**: Initializes the `Inputfile` object. It can be initialized in the same ways as a standard Python dictionary (e.g., empty, from another dictionary, or from keyword arguments).
    During initialization, it iterates through functions in `BigDFT.InputActions` and binds them as methods to the instance using `functools.partial`. This means an action like `InputActions.set_xc(inp, xc='PBE')` can be called as `inp_obj.set_xc(xc='PBE')` on an `Inputfile` instance `inp_obj`.

-   **Dictionary Behavior**: Since it inherits from `dict`, all standard dictionary methods are available (e.g., `inp['dft'] = {'ixc': 'LDA'}`, `inp.get('dft')`, `del inp['posinp']`). The `Inputfile` object itself *is* the dictionary that will be written as a YAML input file.

-   **Dynamically Added Methods**: All functions from `BigDFT.InputActions` become methods of an `Inputfile` instance. For example, if `set_hgrid(inp, hgrids=0.4)` is a function in `InputActions`, an `Inputfile` instance `inp_obj` can call `inp_obj.set_hgrid(hgrids=0.4)`. Refer to the documentation for `PyBigDFT/BigDFT/InputActions.md` for a list and description of these action methods.
    The `remove(action_function)` method, also dynamically bound from `InputActions`, allows previously set actions to be undone (e.g., `inp_obj.remove(InputActions.set_xc)`).

## Important Variables/Constants

There are no module-level variables or constants in `Inputfiles.py` intended for direct external configuration. The primary content is the `Inputfile` class itself.

## Usage Examples

```python
from PyBigDFT.BigDFT import Inputfiles
from PyBigDFT.BigDFT import InputActions # Needed for the 'remove' method's argument
import yaml # For pretty printing the dictionary

# Create an empty Inputfile object
inp = Inputfiles.Inputfile()

# Use dynamically added methods (from InputActions) to set parameters
# These calls modify the 'inp' dictionary directly.
inp.set_xc(xc='PBE')
inp.set_hgrid(hgrids=0.35)
inp.add_empty_SCF_orbitals(norbs=20)
inp.optimize_geometry(method='FIRE', nsteps=100)

# Add positions using dictionary-like assignment (standard dict behavior)
inp['posinp'] = {
    'units': 'angstroem',
    'positions': [
        {'H': [0.0, 0.0, 0.0]},
        {'H': [0.0, 0.0, 0.74]}
    ]
}

# Print the resulting input dictionary using yaml.dump for readability
print("#-- Inputfile content --")
print(yaml.dump(inp, default_flow_style=False))
# Expected output (order of top-level keys might vary):
# dft:
#   hgrids: 0.35
#   ixc: PBE
# geopt:
#   method: FIRE
#   ncount_cluster_x: 100
# mix:
#   norbsempty: 20
# posinp:
#   positions:
#   - H:
#     - 0.0
#     - 0.0
#     - 0.0
#   - H:
#     - 0.0
#     - 0.0
#     - 0.74
#   units: angstroem

# Remove one of the actions
# Note: You pass the original function from InputActions module to the remove method
inp.remove(InputActions.optimize_geometry)

print("\n#-- Inputfile after removing optimize_geometry --")
print(yaml.dump(inp, default_flow_style=False))
# Expected output (geopt section removed):
# dft:
#   hgrids: 0.35
#   ixc: PBE
# mix:
#   norbsempty: 20
# posinp:
#   positions:
#   - H:
#     - 0.0
#     - 0.0
#     - 0.0
#   - H:
#     - 0.0
#     - 0.0
#     - 0.74
#   units: angstroem

# Initialize from an existing dictionary
existing_dict = {"dft": {"ixc": "LDA"}, "output": {"orbitals": "text"}}
inp_from_dict = Inputfiles.Inputfile(existing_dict)
# Now use action methods to add more settings
inp_from_dict.set_rmult(coarse=4.5, fine=7.5)
print("\n#-- Inputfile initialized from dict and modified --")
print(yaml.dump(inp_from_dict, default_flow_style=False))
# Expected output:
# dft:
#   ixc: LDA
#   rmult:
#   - 4.5
#   - 7.5
# output:
#   orbitals: text
```

## Dependencies and Interactions

-   **`BigDFT.InputActions`**: This is a critical dependency. The `Inputfile` class dynamically imports all functions from `InputActions` (except private ones) and makes them available as methods on its instances.
-   **`functools.partial`**: Used internally during initialization to bind the action functions from `InputActions` as methods to the `Inputfile` instance, pre-filling the first `inp` argument.
-   **Standard Python `dict`**: The `Inputfile` class inherits from `dict`, so it relies on the standard dictionary implementation for its base behavior and methods.
```
