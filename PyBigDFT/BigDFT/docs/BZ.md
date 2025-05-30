# `PyBigDFT/BigDFT/BZ.py`

## Overview

The `BZ.py` module is dedicated to handling calculations and representations related to the Brillouin Zone (BZ) in electronic structure calculations. It provides classes to manage band data, define k-point paths, and interpolate band structures within the BZ. This module is essential for analyzing the electronic band structure of periodic systems.

## Key Components

### `get_ev(ev: dict, keys: list[str] = None, ikpt: int = 1) -> list or bool`
A utility function to extract eigenvalue energies from a dictionary, typically representing a single k-point's data from a log file. It can search for energies under various common keys.
-   `ev`: Dictionary containing eigenvalue data.
-   `keys`: Optional list of keys to search for energy values.
-   `ikpt`: The k-point index to match if k-point specific data is present.
-   Returns: A list of energies (possibly `[None, energy_spin_down]` for spin-polarized cases) or `False` if not found.

### `astruct_to_cell(astruct: dict) -> tuple`
Converts atomic structure information (as parsed by `Logfiles` module) into the cell format required by the `spglib` library.
-   `astruct`: A dictionary containing 'cell' vectors and 'positions' (list of atom dicts).
-   Returns: A tuple `(lattice, positions, atom_types)` compatible with `spglib`.

### `class BandArray(numpy.ndarray)`
A specialized `numpy.ndarray` subclass designed to store band energies for one or more k-points. It handles spin-polarized data by having a shape of (1 or 2, num_bands).
-   **`__new__(cls, *args, **kwargs)`**: Constructor that can take raw data (`data=...`) or parse eigenvalues from `logdata` for a specific `ikpt`.
-   **`__init__(self, *args, **kwargs)`**: Initializes k-point information (`ikpt`, `kpt` coordinates, `kwgt` weight).
-   **`set_kpt(self, ikpt: int, kpt: tuple[float, float, float], kwgt: float = 1.0)`**: Sets or updates the k-point identifier, coordinates, and weight.
-   **`info` (attribute)**: Stores the number of orbitals for each spin channel (e.g., `(num_spin_up, num_spin_down)`).

### `class BZPath()`
Represents a path in the k-space (Brillouin Zone) defined by a series of special k-points or explicit coordinates.
-   **`__init__(self, lattice: list[list[float]], path: list[str or dict], special_points: dict, npts: int = 50)`**:
    -   `lattice`: The reciprocal lattice vectors.
    -   `path`: A list defining the path. Elements can be strings (keys from `special_points`) or dictionaries like `{'label': [kx, ky, kz]}`.
    -   `special_points`: A dictionary mapping special k-point labels (e.g., "G", "X") to their reciprocal coordinates.
    -   `npts`: Number of points to interpolate along each segment of the path.
-   **Attributes**:
    -   `path`: Interpolated k-point coordinates along the path.
    -   `xaxis`: Corresponding x-axis values for plotting.
    -   `xlabel`: Positions of the special points on the x-axis for labels.
    -   `symbols`: Labels for the special points.

### `class BrillouinZone()`
The main class for representing the full Brillouin Zone and performing operations like band structure interpolation.
-   **`__init__(self, astruct: dict, mesh: list[int], evals: list[BandArray], fermi_energy: float)`**:
    -   `astruct`: Atomic structure dictionary.
    -   `mesh`: The k-point mesh dimensions (e.g., `[nkx, nky, nkz]`).
    -   `evals`: A list of `BandArray` objects, one for each irreducible k-point.
    -   `fermi_energy`: The Fermi energy in Hartrees.
-   **`plot(self, path: BZPath = None, npts: int = 50)`**: Plots the interpolated band structure along the given `BZPath` (or a default path if `None`).
-   **`interpolator` (attribute)**: A `scipy.interpolate.interpnd.LinearNDInterpolator` object used for band interpolation.

## Important Variables/Constants

None of the module-level variables appear to be critical constants intended for direct external use beyond their role in internal logic.

## Usage Examples

```python
# This module is typically used with data from a Logfile
from PyBigDFT.Logfiles import Logfile
from PyBigDFT.BigDFT import BZ # Assuming BZ classes are available under this import

# Assuming 'log.yaml' is a BigDFT logfile from a periodic calculation
# For this example to run, you would need a real 'log.yaml' file.
# log = Logfile('log.yaml')

# Placeholder for logfile data for demonstration if a real logfile isn't available
# In a real scenario, these would be populated from 'log'.
class MockLogfile:
    def __init__(self):
        self.kpts = [([0.0, 0.0, 0.0], 0.25), ([0.5, 0.0, 0.0], 0.25)] # Example k-points
        self.nkpt = len(self.kpts)
        self.astruct = {
            'cell': [10.0, 10.0, 10.0], # Example cubic cell
            'positions': [{'Si': [0.0, 0.0, 0.0]}, {'Si': [0.25, 0.25, 0.25]}]
        }
        self.kpt_mesh = [4, 4, 4] # Example k-point mesh
        # Example BandArray objects (simplified structure for this mock)
        # In reality, BandArray instances would be properly initialized with energy data
        band_data1 = BZ.BandArray(data=[[-5.0, -4.0, 1.0, 2.0]], ikpt=1, kpt=(0.0,0.0,0.0), kwgt=0.125)
        band_data2 = BZ.BandArray(data=[[-4.5, -3.5, 1.5, 2.5]], ikpt=2, kpt=(0.5,0.0,0.0), kwgt=0.125)
        # ... more BandArray objects for all irreducible k-points
        self.evals = [band_data1, band_data2] # Needs to cover all k-points in the mesh after symmetry reduction
        self.fermi_level = 0.0 # Example Fermi level in Hartrees

log = MockLogfile() # Using the mock for this example snippet

if hasattr(log, 'kpts') and log.nkpt > 1:
    astruct = log.astruct
    mesh = log.kpt_mesh if isinstance(log.kpt_mesh, list) else [log.kpt_mesh]*3
    if astruct['cell'][1] == float('inf'): mesh[1] = 1
    if astruct['cell'][2] == float('inf'): mesh[2] = 1

    evals_list = log.evals
    fermi_energy_hartree = log.fermi_level

    # Create a BrillouinZone object
    # Note: For a real run, evals_list must contain BandArray objects for *all*
    # irreducible k-points corresponding to the mesh and symmetries.
    # The MockLogfile provides a very simplified 'evals' for structure demonstration.
    # A real 'evals' list from a Logfile object would be more complex.
    # The following line might raise errors with the mock data due to spglib processing
    # if the mock data isn't perfectly consistent for symmetry analysis.
    try:
        bz = BZ.BrillouinZone(astruct, mesh, evals_list, fermi_energy_hartree)

        # Define a custom path (example for a cubic lattice, adapt as needed)
        # path_definition = ['G', 'X', 'M', 'G', 'R', 'X']
        # custom_bz_path = BZ.BZPath(bz.lattice, path_definition, bz.special_points, npts=50)

        # Plot the band structure
        # bz.plot(path=custom_bz_path) # For a custom path
        # For a default path (if defined for the detected lattice type):
        print("Attempting to plot Brillouin Zone (requires matplotlib and full k-point data)...")
        # bz.plot() # This line would generate a plot.
        print("Plotting example skipped in non-graphical environment or with mock data.")
    except Exception as e:
        print(f"Could not initialize or plot BrillouinZone, possibly due to mock data limitations or missing dependencies: {e}")
else:
    print("Logfile does not contain k-point information or has only one k-point.")

```
*Note: Running the example requires a valid BigDFT log file from a periodic calculation with k-points and necessary libraries like `matplotlib` and `spglib`. The example code includes a placeholder `MockLogfile` to illustrate data structures, but real data is needed for full functionality.*

## Dependencies and Interactions

-   **`numpy`**: Heavily used for numerical operations and array manipulations.
-   **`matplotlib.pyplot`**: Used for plotting band structures in `BrillouinZone.plot()`.
-   **`scipy.interpolate.interpnd`**: Used by `BrillouinZone` for band interpolation.
-   **`ase.dft.kpoints`**: Used by `BZPath` and `BrillouinZone` for getting special k-points and defining k-point paths.
-   **`spglib`**: Used by `BrillouinZone` for space group analysis and handling irreducible k-points.
-   **`futile.Utils.write`**: Used for safe printing within the module.
-   **`PyBigDFT.Logfiles`**: `BrillouinZone` is often initialized with data extracted via the `Logfile` class from this sibling module.
```
