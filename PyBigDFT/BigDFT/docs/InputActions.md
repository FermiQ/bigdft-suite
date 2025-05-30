# `PyBigDFT/BigDFT/InputActions.py`

## Overview

The `InputActions.py` module provides a collection of utility functions designed to modify BigDFT input dictionaries. These functions encapsulate common input parameter settings, making it easier to configure calculations. A key feature of this module is that its functions are dynamically added as methods to the `Inputfiles.Inputfile` class at runtime (this binding happens in the `Inputfiles` module itself). This allows for a fluent and intuitive API when working with `Inputfile` objects.

Each function in this module typically takes an input dictionary (`inp`) as its first argument, followed by parameters specific to the action. When called as a method of an `Inputfile` object, the `inp` argument is implicitly the object's internal dictionary.

## Key Components (Action Functions)

This module consists primarily of functions that manipulate the input dictionary. Below are some of the most commonly used actions:

### General Setup
-   **`set_hgrid(inp: dict, hgrids: float or list[float] = 0.4)`**: Sets the wavelet grid spacing(s) in Bohr.
-   **`set_rmult(inp: dict, rmult: list[float] = None, coarse: float = 5.0, fine: float = 8.0)`**: Sets the coarse and fine wavelet grid extension radii multipliers. `rmult` takes precedence if provided.
-   **`set_xc(inp: dict, xc: str = 'PBE')`**: Sets the exchange-correlation functional (e.g., 'LDA', 'PBE', 'HF', or LibXC codes).
-   **`set_mesh_sizes(inp: dict, ngrids: int or list[int] = 64)`**: Constrains the number of grid points in each direction. Useful for maintaining constant degrees of freedom in variable cell calculations.
-   **`set_symmetry(inp: dict, yes: bool = True)`**: Enables (`True`) or disables (`False`) symmetry detection.
-   **`change_data_directory(inp: dict, name: str = '')`**: Modifies the name of the `data-` directory (sets the `radical` key). Useful for linking to data from a previous run. If `name` is empty, it might remove the `radical` key.
-   **`connect_run_data(inp: dict, log: Logfile = None)`**: Sets the data directory (`radical`) of the current input to match that of a provided `BigDFT.Logfiles.Logfile` object.

### Self-Consistent Field (SCF) Cycle
-   **`set_wavefunction_convergence(inp: dict, gnrm: float = 1.0e-04)`**: Sets the convergence tolerance for the wavefunction residue (`gnrm_cv`).
-   **`set_SCF_method(inp: dict, method: str = 'dirmin', mixing_on: str = 'density', mixing_scheme: str = 'Pulay')`**: Configures the SCF algorithm.
    -   `method`: 'dirmin', 'mixing' (cubic), 'foe', 'pexsi', 'diag' (linear scaling).
    -   `mixing_on`: 'density' or 'potential' (for `method='mixing'`).
    -   `mixing_scheme`: 'Pulay', 'Simple', 'Anderson', 'Anderson2', 'CG' (for `method='mixing'`).
-   **`set_wavefunction_iterations(inp: dict, nit: list[int] or int = [50,1])`**: Sets the maximum number of iterations for the SCF cycle (`itermax`) and, if `nit` is a list of two, subspace iterations (`nrepmax`).
-   **`set_electronic_temperature(inp: dict, kT: float = 1.e-3, T: float = 0)`**: Defines the electronic temperature for Fermi-Dirac smearing. `kT` is in Hartrees, `T` is in Kelvin (T takes precedence if non-zero).
-   **`spin_polarize(inp: dict, mpol: int = 1)`**: Enables spin-polarized calculation (`nspin=2`) with a given total magnetic moment (`mpol`).
-   **`charge(inp: dict, charge: float = -1)`**: Sets the total charge of the system (`qcharge`).
-   **`charge_and_polarize(inp: dict)`**: Convenience function to set charge to +1 and mpol to 1 (typically for removing an electron from a closed-shell system).

### Orbitals and Density
-   **`set_random_inputguess(inp: dict)`**: Initializes input orbitals with random coefficients (`inputpsiid=-2`).
-   **`add_empty_SCF_orbitals(inp: dict, norbs: int = 10)`**: Includes a specified number of empty orbitals (`norbsempty` or `extra_states`) in the SCF procedure.
-   **`extract_virtual_states(inp: dict, nvirt: int, davidson: bool = False)`**: Calculates a number of virtual (empty) states after the SCF cycle. Uses Davidson diagonalization if `davidson=True`, otherwise Trace Minimization.
-   **`write_orbitals_on_disk(inp: dict, format: str = 'binary')`**: Instructs BigDFT to write converged orbitals to disk. Formats: 'binary', 'text', 'etsf'.
-   **`read_orbitals_from_disk(inp: dict)`**: Instructs BigDFT to read orbitals from the data directory (restart, sets `inputpsiid` to 2 or 102).
-   **`write_density_on_disk(inp: dict)`**: Writes the converged charge density to disk (`output_denspot=21`).
-   **`calculate_dipole(inp: dict)`**: Requests the calculation of the electric dipole moment (sets `lin_general.calc_dipole=True`, useful for linear-scaling).
-   **`write_cubefiles_around_fermi_level(inp: dict, nplot: int = 1)`**: Writes a specified number of orbitals (`nplot`) around the Fermi level in CUBE format.

### Geometry Optimization
-   **`optimize_geometry(inp: dict, method: str = 'FIRE', nsteps: int = 50)`**: Configures a geometry optimization run.
    -   `method`: Optimizer algorithm like 'FIRE', 'LBFGS', 'SDCG', 'VSSD', 'BFGS', 'DIIS', 'SQNM'.
    -   `nsteps`: Maximum number of optimization steps (`ncount_cluster_x`).

### External Fields and Advanced Features
-   **`apply_electric_field(inp: dict, elecfield: list[float] = [0,0,1.e-3])`**: Applies a static external electric field.
-   **`use_gpu_acceleration(inp: dict)`**: Enables GPU acceleration (`perf.accel='OCLGPU'`, `psolver.setup.accel='CUDA'`) if available.
-   **`calculate_tddft_coupling_matrix(inp: dict, tda: bool = False, rpa: bool = True, fxc: bool = True)`**: Configures a Time-Dependent DFT (TDDFT) Casida coupling matrix calculation.
    - `tda`: Use Tamm-Dancoff approximation if `True`.
    - `rpa`, `fxc`: Control calculation of RPA and fxc terms.
-   **`set_linear_scaling(inp: dict)`**: Activates linear scaling calculation mode (sets `inputpsiid` to 'linear' or 102 if restarting).
-   **`set_atomic_positions(inp: dict, posinp: dict = None)`**: Directly embeds atomic positions (a `posinp` dictionary, e.g., `{"positions": [{"H": [0,0,0]}], "units": "angstroem"}`) into the main input dictionary under the 'posinp' key.

### Removing Actions
-   **`remove(inp: dict, action: callable)`**: Reverts a previously applied action by attempting to delete the keys that `action` would have set. The `action` argument must be one of the functions from this module (e.g., `InputActions.set_xc`). This works by temporarily redirecting the internal setter to a remover function and then calling the specified action with its default arguments.

## Important Variables/Constants

-   **`__set__`**: An internal module variable that points to `futile.Utils.dict_set` by default. It's temporarily changed by the `remove` function to an internal `__undo__` function to revert settings. This is an implementation detail and not for external use.

## Usage Examples

```python
from PyBigDFT.BigDFT import InputActions, Inputfiles # Assuming Inputfiles is where Inputfile class is defined
from PyBigDFT.Logfiles import Logfile # For connect_run_data example

# Example 1: Using functions directly on a dictionary
inp_dict = {}
InputActions.set_xc(inp_dict, xc='LDA')
InputActions.set_hgrid(inp_dict, hgrids=0.35)
InputActions.spin_polarize(inp_dict, mpol=2)
print("inp_dict after actions:", inp_dict)

InputActions.remove(inp_dict, InputActions.spin_polarize) # Attempt to remove spin_polarize settings
print("inp_dict after removing spin_polarize:", inp_dict) # Should revert 'dft.nspin' and 'dft.mpol'

# Example 2: Using as methods of an Inputfile object
# This is the more common and recommended usage.
inp_obj = Inputfiles.Inputfile() # Create an Inputfile instance
inp_obj.set_xc(xc='PBE') # 'inp' argument is implicit (it's inp_obj.inp)
inp_obj.optimize_geometry(method='LBFGS', nsteps=100)
inp_obj.write_orbitals_on_disk(format='text')
print("\nInputfile object's dictionary after actions:", inp_obj.inp)

# To remove an action when using Inputfile object:
# The remove method of Inputfile also takes the action function from InputActions
inp_obj.remove(InputActions.optimize_geometry)
print("Inputfile object's dictionary after removing optimize_geometry:", inp_obj.inp)

# Example with connect_run_data (conceptual, needs a real Logfile object)
# mock_log = Logfile() # In reality, load a Logfile from a file
# mock_log.log = {'radical': 'previous_run_data_dir'} # Simulate a loaded logfile's radical
# inp_for_restart = Inputfiles.Inputfile()
# inp_for_restart.connect_run_data(log=mock_log)
# print("\nInput for restart after connect_run_data:", inp_for_restart.inp)

```

## Dependencies and Interactions

-   **`futile.Utils.dict_set`**: This function is the primary mechanism used by most actions to set values within the nested input dictionary.
-   **`BigDFT.Inputfiles.Inputfile`**: The action functions in this module are designed to be dynamically bound as methods to `Inputfile` objects (this binding occurs in `Inputfiles.__init__.py`). This provides a convenient, object-oriented way to construct input files.
-   **`BigDFT.Logfiles.Logfile`**: The `connect_run_data` action takes a `Logfile` object as an argument to link data directories.
```
