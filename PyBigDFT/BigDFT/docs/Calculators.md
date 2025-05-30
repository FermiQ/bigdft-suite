# `PyBigDFT/BigDFT/Calculators.py`

## Overview

The `Calculators.py` module provides the primary means for executing BigDFT calculations from Python. It offers two main calculator classes: `GIBinding` for interacting with BigDFT via GObject Introspection bindings, and `SystemCalculator` for running BigDFT as an external system command. Both inherit from a base `Runner` class that standardizes the execution flow.

This module aims to provide a light interface, somewhat compatible with the concepts in ASE (Atomic Simulation Environment).

## Key Components

### `class Runner()`
A base class that defines a generic framework for running processes. It handles global and local options for a run.
-   **`__init__(self, **kwargs)`**: Initializes the runner with global options.
-   **`global_options(self) -> dict`**: Returns the current global options.
-   **`get_global_option(self, key: str)`**: Retrieves a specific global option.
-   **`update_global_options(self, **kwargs)`**: Updates the global options.
-   **`run(self, **kwargs)`**: The main execution method. It calls `pre_processing`, `process_run`, and `post_processing` in sequence.
-   **`pre_processing(self) -> dict`**: Prepares arguments for `process_run`. Returns a dictionary of arguments.
-   **`process_run(self, **kwargs) -> dict`**: The core method that executes the main task. Returns a dictionary of results.
-   **`post_processing(self, **kwargs) -> any`**: Processes the results from `process_run` and returns the final output of the `run` method.

### `class GIBinding(Runner)`
A calculator that uses GObject Introspection bindings to interact with the BigDFT library directly. This is suitable for more tightly coupled integrations.
-   **`__init__(self)`**: Initializes MPI and the BigDFT library bindings.
-   **`update(self, inputfile: dict)`**: Updates the BigDFT run object with new input parameters, typically for a restart or continued calculation.
-   **`run(self) -> BigDFT.Run`**: Executes the calculation using the bindings. Returns a BigDFT Run object (specific to the GI bindings, often referred to as `out` in examples, with attributes like `eKS`).
-   **`set(self, inputfile: dict = None)`**: (Re)initializes the BigDFT run object with a given input file dictionary.
-   **`__del__(self)`**: Finalizes MPI and cleans up resources.

### `class SystemCalculator(Runner)`
A calculator that runs BigDFT by making system calls to the `bigdft` executable. This is suitable for workstation-based runs or for scripts submitted to batch schedulers.
-   **`__init__(self, omp: str = <OMP_NUM_THREADS_env_var_or_1>, mpi_run: str = <BIGDFT_MPIRUN_env_var_or_''>, dry_run: bool = False, skip: bool = False, verbose: bool = True)`**:
    -   `omp`: Number of OpenMP threads.
    -   `mpi_run`: MPI command (e.g., `mpirun -np 4`).
    -   `dry_run`: If `True`, checks input and estimates memory without running.
    -   `skip`: If `True`, skips calculation if the log file already exists.
    -   `verbose`: If `True`, prints information about operations.
-   **`pre_processing(self) -> dict`**: Creates the YAML input file from the `input` and `posinp` arguments provided to `run`, and prepares the command line.
-   **`process_run(self, command: str) -> dict`**: Executes the `bigdft` command.
-   **`post_processing(self, timedbg: float, logname: str, command: str) -> Logfile or None`**: Checks for errors and returns a `BigDFT.Logfiles.Logfile` instance. Returns `None` if an error occurred or the log file is problematic.
-   **Key `run` arguments for `SystemCalculator`** (passed as `**kwargs` to `run` and accessed via `self.run_options.get()`):
    -   `name: str`: Basename for input (`<name>.yaml`) and output (`log-<name>.yaml`) files.
    -   `run_dir: str`: Directory where `bigdft` will be executed (defaults to current directory).
    -   `input: dict`: Dictionary representing the BigDFT input parameters.
    -   `posinp: str or dict`: Path to the atomic position file (e.g., XYZ) or a dictionary representing positions.
    -   `outdir: str`: Output directory for BigDFT data (passed as `-d` to `bigdft`).
    -   `run_name: str`: File containing a list of run IDs for task-farmed runs (passed as `-r`).
    -   `taskgroup_size: str`: Number of MPI processes per task group (passed as `-t`).

## Important Variables/Constants

No module-level variables are intended for direct external configuration. The behavior is primarily controlled by class initialization parameters and environment variables (`OMP_NUM_THREADS`, `BIGDFT_MPIRUN`, `BIGDFT_ROOT`).

## Usage Examples

### `SystemCalculator` Example:

```python
from PyBigDFT.BigDFT import Calculators, Inputfiles

# Define a simple input
inp = Inputfiles.Inputfile() # Helper to build input dictionaries
inp.set_xc('LDA')
inp.add_posinp({"H": [0.0, 0.0, 0.0]}, units="angstroem") # Add positions directly
inp.set_dft({"hgrids": 0.4})

# Initialize the calculator (ensure BIGDFT_ROOT is set in your environment)
# For this example, we'll assume it runs sequentially if mpi_run is not set.
try:
    calc = Calculators.SystemCalculator(omp="1", verbose=True)

    # Run the calculation
    # The 'input' argument to run() takes the dictionary representation
    log = calc.run(name="H_atom_test", input=inp.get_inputdict())

    if log and log.get_energy() is not None:
        print("SystemCalculator Calculation successful.")
        print(f"Energy: {log.energy}")
    elif log:
        print("SystemCalculator Calculation finished, but energy not found or error in log.")
    else:
        print("SystemCalculator Calculation failed or log file was not created/parsed.")

except AssertionError as e:
    print(f"SystemCalculator requirement failed: {e}")
    print("Please ensure BIGDFT_ROOT is set and 'bigdft' executable is found.")
except ValueError as e:
    print(f"ValueError during SystemCalculator run: {e}")
except Exception as e:
    print(f"An unexpected error occurred during SystemCalculator example: {e}")
finally:
    # Clean up created files for the example
    import os
    import shutil
    files_to_remove = ["H_atom_test.yaml", "log-H_atom_test.yaml", "time-H_atom_test.yaml"]
    dirs_to_remove = ["data-H_atom_test", "debug"]
    for f in files_to_remove:
        if os.path.exists(f): os.remove(f)
    for d in dirs_to_remove:
        if os.path.exists(d): shutil.rmtree(d)

```

### `GIBinding` Example:

*Note: Running the `GIBinding` example requires that the GObject Introspection bindings for BigDFT are correctly installed and accessible in your Python environment. This often involves specific compilation and installation steps for BigDFT.*

```python
from PyBigDFT.BigDFT import Calculators
import yaml

# Example input from the __main__ block of Calculators.py
basic_input_str = """
#mode: {method: lj}
logfile: No
dft: { ixc: HF, nspin: 2} # Using HF for example, LDA/PBE more common for DFT
posinp:
   positions:
   - {Be : [0.0, 0.0, 0.0]}
   - {Be : [0.0, 0.0, 1.0]}
   units: bohr # Assuming Bohr from context, adjust if Angstroem preferred
#   properties: {format: yaml}
ig_occupation:
   Atom 1: {2s: {up: 1.0, down: 0.9}, 2p: {up: 0.0, down: 0.2} }
   Atom 2: {2s: {up: 0.9, down: 1.0}, 2p: {up: 0.2, down: 0.0} }

psppar.Be: {Pseudopotential XC: 11} # Example, ensure this matches available PSPs
"""
inp_dict = yaml.safe_load(basic_input_str)

try:
    # Initialize the GIBinding calculator
    study = Calculators.GIBinding()

    # Set the input
    study.set(inp_dict)

    # Perform the first calculation
    out = study.run() # 'out' is a BigDFT.Run GObject

    if out and hasattr(out, 'eKS'):
        print(f"GIBinding initial calculation successful. Energy (EKS): {out.eKS}")

        # Example: Perform a dissociation curve (simplified from __main__)
        energy_curve = [out.eKS]
        positions_z = [inp_dict['posinp']['positions'][-1]['Be'][2]]

        for i in range(3): # Reduced loop for brevity
            shift = float(i + 1) * 0.05 # Small shifts
            inp_dict['posinp']['positions'][-1]['Be'][2] += shift
            study.update(inp_dict) # Update the calculator with new positions
            out_updated = study.run()
            if out_updated and hasattr(out_updated, 'eKS'):
                energy_curve.append(out_updated.eKS)
                positions_z.append(inp_dict['posinp']['positions'][-1]['Be'][2])
                if study.iproc == 0: # Print output only from the master process
                    print(f"Iteration {i}, Shift {shift:.2f}, New Z: {positions_z[-1]:.2f}, Energy: {out_updated.eKS}")
            else:
                print("Failed to get energy in update step.")
                break

        if study.iproc == 0 and len(energy_curve) > 1:
            print("Dissociation curve data points (Z vs Energy):")
            for z, e in zip(positions_z, energy_curve):
                print(f"{z:.2f}, {e}")
            # Plotting example (requires matplotlib)
            # import matplotlib.pyplot as plt
            # plt.plot(positions_z, energy_curve)
            # plt.xlabel("Z coordinate of second Be atom (bohr)")
            # plt.ylabel("Total Energy (Hartree)")
            # plt.title("Be Dimer Dissociation (GIBinding)")
            # plt.show()
            print("Plotting skipped in this example text.")

    else:
        print("GIBinding initial calculation failed or eKS attribute not found.")

except ImportError:
    print("GIBinding example skipped: GObject Introspection bindings for BigDFT not found or gi repository not available.")
except Exception as e:
    print(f"An error occurred during GIBinding example: {e}")

```

## Dependencies and Interactions

-   **`os`**: Used by `SystemCalculator` for environment variables and system calls.
-   **`shutil`**: Used by `SystemCalculator` (e.g., for directory operations if `run_dir` is created, and in cleanup for examples).
-   **`futile.Utils`**: For utility functions like `write` (safe_print) and `dict_merge`.
-   **`futile.YamlIO`**: Used by `SystemCalculator` for reading/writing YAML input files.
-   **`BigDFT.Logfiles`**: `SystemCalculator` returns `Logfile` objects.
-   **`BigDFT.Inputfiles`**: Used in the `SystemCalculator` example to prepare input.
-   **`gi.repository.BigDFT`**: Required by `GIBinding` for GObject Introspection.
-   **`yaml`**: Used in examples for loading string-based input definitions.
-   The `bigdft` executable and related scripts/tools from the BigDFT installation are crucial for `SystemCalculator`.
-   Environment variables `OMP_NUM_THREADS`, `BIGDFT_MPIRUN`, and `BIGDFT_ROOT` significantly influence the behavior of `SystemCalculator`.
```
