# `PyBigDFT/BigDFT/DoS.py`

## Overview

The `DoS.py` module is responsible for calculating and visualizing the Density of States (DoS) from a set of energy eigenvalues. It can handle data from single calculations or k-point sampled calculations, apply smearing (Gaussian by default), and plot the resulting DoS curves. It also has capabilities to incorporate spatially resolved DoS (SDOS) if such data is available.

## Key Components

### `class DiracSuperposition()`
Represents a DoS as a sum of Dirac delta functions (which are then typically broadened).
-   **`__init__(self, dos: numpy.ndarray, wgts: list[float] or float = [1.0])`**:
    -   `dos`: A NumPy array of energy eigenvalues. For k-point data, this might be a list of arrays (or an array of arrays).
    -   `wgts`: Weights associated with each set of energies (e.g., k-point weights or spin polarization factors). Can be a single float or a list matching `dos`.
-   **`curve(self, xs: numpy.ndarray, sigma: float, wgts: list[float] = None) -> tuple[numpy.ndarray, numpy.ndarray]`**: Calculates the broadened DoS curve.
    -   `xs`: Array of energy values at which to evaluate the DoS.
    -   `sigma`: Smearing width (e.g., Gaussian sigma).
    -   `wgts`: Optional further weighting for sub-components of the DoS (e.g., for SDOS contributions).
    -   Returns: A tuple `(xs, broadened_dos_values)`.
-   **`peak(self, omega: float, e: float, sigma: float) -> float`**: Calculates the value of a single broadened peak (currently Gaussian).
-   **`peaks(self, xs: numpy.ndarray, dos: numpy.ndarray, norms: numpy.ndarray, sigma: float) -> numpy.ndarray`**: Calculates the sum of broadened peaks for a set of energies, applying given normalizations.

### `class DoS()`
The main class for managing and plotting Density of States.
-   **`__init__(self, bandarrays: list[BandArray] = None, energies: numpy.ndarray = None, evals: list[dict] = None, units: str = 'eV', label: str = '1', sigma: float = 0.2/AU_eV, npts: int = 2500, fermi_level: float = None, norm: float or list[float] = 1.0, sdos: list[dict] = None)`**:
    -   Can be initialized from `bandarrays` (list of `BigDFT.BZ.BandArray` objects), raw `energies` (a NumPy array or list of arrays), or `evals` (list of eigenvalue dictionaries from a log file).
    -   `units`: Units of input energies and `fermi_level` ('eV' or 'AU'). Output plots are always in eV.
    -   `label`: Default label for the first dataset.
    -   `sigma`: Default smearing width in Hartrees (converted to eV internally if input units are AU).
    -   `npts`: Number of points for the energy grid of the plot.
    -   `fermi_level`: The Fermi energy.
    -   `norm`: Default normalization factor for the first dataset.
    -   `sdos`: Spatially resolved DoS data.
-   **`append_from_bandarray(self, bandarrays: list[BandArray], label: str)`**: Adds DoS data from a list of `BandArray` objects (typically for k-points).
-   **`append_from_dict(self, evals: list[dict], label: str)`**: Adds DoS data from a list of eigenvalue dictionaries.
-   **`append(self, energies: numpy.ndarray, label: str = None, units: str = 'eV', norm: float or list[float] = 1.0)`**: Appends a new set of energies to the DoS object. `energies` can be a single array or a list of arrays (e.g., for different k-points or spin channels if handled manually).
-   **`fermi_level(self, fermi_level: float, units: str = 'eV')`**: Sets or updates the Fermi level.
-   **`dump(self, sigma: float = None)`**: Prints the DoS curve data (energy and DoS value for each dataset) to standard output, suitable for tools like Gnuplot.
-   **`plot(self, sigma: float = None, legend: bool = True, xlmin: float = None, xlmax: float = None, ylmin: float = None, ylmax: float = None)`**: Plots the DoS using `matplotlib`. Includes an interactive slider for adjusting the smearing (`sigma`). If SDOS data is present, additional controls for visualizing SDOS might appear.
-   **`update(self, val: float)`**: Callback function for the smearing slider to update the plot.

### `_bandarray_to_data(jspin: int, bandarrays: list[BandArray]) -> tuple`
Internal helper function to extract energy and weight data for a specific spin channel from a list of `BigDFT.BZ.BandArray` objects.

## Important Variables/Constants

-   **`AU_eV: float`**: Conversion factor from Hartrees (Atomic Units of energy) to electronVolts. Value: `27.21138386`.

## Usage Examples

```python
from PyBigDFT.BigDFT import DoS
from PyBigDFT.Logfiles import Logfile # Typically used with Logfile data
import numpy as np
import matplotlib.pyplot as plt # For displaying the plot

# Example 1: From a simple list of energies (adapted from DoS.py __main__)
raw_energies_au = np.array([
    -0.815, -0.803, -0.780, -0.750, -0.723, -0.714, -0.710, -0.687,
    -0.672, -0.659, -0.625, -0.608, -0.565, -0.561, -0.551, -0.541,
    -0.532, -0.515, -0.474, -0.473, -0.465, -0.445, -0.433, -0.416,
    -0.407, -0.406, -0.403, -0.389, -0.380, -0.375, -0.375, -0.367,
    -0.367, -0.359, -0.358, -0.354, -0.334, -0.332, -0.315, -0.308,
    -0.298, -0.294, -0.292, -0.285, -0.284, -0.267, -0.259, -0.239,
    -0.224, -0.204, -0.164, -0.117, -0.071, -0.052, -0.034, -0.016,
    -0.013, -0.010, 0.007,  0.009,  0.010,  0.021,  0.023,  0.041,
    0.042,  0.045,  0.051
]) # in Hartrees

# Initialize DoS object with the first set of energies
# Note: sigma in DoS constructor is in Hartrees (AU_eV is used for conversion if units='eV')
# Here, units='AU', so sigma is directly taken as Hartrees.
dos_from_energies = DoS.DoS(energies=raw_energies_au, units='AU', fermi_level=-0.1, sigma=0.01)

# Append another set of energies (e.g., shifted)
dos_from_energies.append(0.2 + raw_energies_au, label="Shifted (+0.2 Ha)", units='AU')

# To plot (interactive window will open if matplotlib backend supports it)
# print("Plotting DoS from raw energies...")
# dos_from_energies.plot()
# plt.show() # Call this to keep the plot window open if not running interactively
# print("Plotting example finished or window closed.")

# Example of dumping data for gnuplot (using a smearing of 0.01 Ha)
# print("\nDoS data for gnuplot (sigma=0.01 Ha):")
# dos_from_energies.dump(sigma=0.01)


# Example 2: From a BigDFT Logfile (conceptual)
# Assuming 'log.yaml' is a BigDFT log file from a calculation.
try:
    # Replace 'log.yaml' with the actual path to your logfile
    log = Logfile('log.yaml')

    if hasattr(log, 'evals') and log.evals is not None and \
       hasattr(log, 'fermi_level') and log.fermi_level is not None:

        # log.evals is a list of BandArray objects.
        # The DoS class handles units conversion if 'units' is specified.
        # Fermi level from logfile is typically in Hartrees.
        dos_from_log = DoS.DoS(bandarrays=log.evals, units='AU',
                               fermi_level=log.fermi_level, label=log.logname) # Use logname as label

        print(f"\nDoS object created from logfile '{log.logname}'.")
        # To plot (interactive window will open)
        # print("Plotting DoS from logfile...")
        # dos_from_log.plot()
        # plt.show()
        # print("Logfile DoS plotting example finished or window closed.")

    else:
        print("\nLogfile 'log.yaml' does not have 'evals' or 'fermi_level' attributes, or they are None.")

except FileNotFoundError:
    print("\nLogfile 'log.yaml' not found. Skipping Logfile-based DoS example.")
except Exception as e:
    print(f"\nAn error occurred with Logfile-based DoS example: {e}")

```
*Note: Running the plotting examples will open interactive matplotlib windows. The Logfile example requires a valid BigDFT log file.*

## Dependencies and Interactions

-   **`numpy`**: Essential for numerical operations, especially handling arrays of energies and performing calculations for broadening.
-   **`matplotlib.pyplot`** and **`matplotlib.widgets.Slider`**: Used for interactive plotting of DoS curves.
-   **`futile.Utils.write`**: Used for safe printing (e.g., in the `dump` method).
-   **`BigDFT.Logfiles.Logfile`**: The `DoS` class is often initialized using data extracted by the `Logfile` class (specifically `evals` which are lists of `BandArray` objects, and `fermi_level`).
-   **`BigDFT.BZ.BandArray`**: `DoS` can directly process `BandArray` objects, which store eigenvalue data, often originating from Brillouin zone calculations.
```
