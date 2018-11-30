"""A module for post processing BigDFT calculations.

"""


def _system_command(command, options):
    """
    Run the command as ``os.system (command + options)``

    Todo:
       remove outfile, this should be controlled by the script.
    Args:
       command (str): the actual command to run.
       options (str): the options to pass to the command.
    """
    from subprocess import call

    command_str = command + " " + options
    call(command_str, shell=True)


def _get_datadir(log):
    """
    Given a log file, this returns the path to the data directory.

    Args:
      log (Logfile): logfile from a BigDFT run.
    Returns:
      str: path to the data directory.
    """
    from os.path import join
    return join(log.srcdir, log.data_directory)


class BigDFTool(object):
    """

    This module defines the actions that can be performed using the
    ``bigdft-tool`` python script. Such a script sometimes invokes the
    ``memguess`` executable, which is executed serially. Sometimes the
    ``utilities`` main program is necessary to take care of more time
    consuming operations which require parallelism.

    Args:
      omp (int): number of OpenMP threads.
        It defaults to the $OMP_NUM_THREADS variable in the environment,
        if present, otherwise it fixes the run to 1 thread.
      mpi_run (str): define the MPI command to be used.
        It defaults to the value $BIGDFT_MPIRUN of the environment, if present.
        When using this calculator into a job submission script, the value of
        $BIGDFT_MPIRUN variable may be set appropriately to launch the bigdft
        executable.
    """

    def __init__(self, omp=1, mpi_run=""):
        from futile.YamlArgparse import get_python_function
        from futile import YamlIO
        from os import environ
        from os.path import join
        from functools import partial
        from copy import deepcopy

        # Executables
        bigdftroot = environ['BIGDFT_ROOT']
        self.bigdft_tool_command = join(bigdftroot, "bigdft-tool")
        self.utilities_command = mpi_run + join(bigdftroot, "utilities")
        self.memguess_command = join(bigdftroot, "memguess")
        self.bigpoly_command = join(bigdftroot, "BigPoly")
        environ['OMP_NUM_THREADS'] = str(omp)

        # Load the dictionary that defines all the operations.
        # This path should be easier to specify if you use the python
        # egg setup.
        with open(join(environ['PYTHONPATH'], "BigDFT", "Database",
                       "postprocess.yaml")) as ifile:
            db = YamlIO.load(stream=ifile)[0]

        # Create the subroutines of the class from the dictionary
        for action, spec in db.items():
            naction = action.replace("_", "-")
            nspec = deepcopy(spec)
            nspec["args"]["mpirun"] = {"default": "" + mpi_run + ""}
            nspec["args"]["action"] = {"default": "" + naction + ""}
            my_action = get_python_function(partial(self._invoke_command,
                                                    self.bigdft_tool_command),
                                            action, nspec)
            setattr(self, action, my_action)

    def _invoke_command(self, command, **kwargs):
        from futile.Utils import option_line_generator
        _system_command(command, option_line_generator(**kwargs))

    def _run_fragment_multipoles(self, log, system=None, orbitals=None):
        """
        Performs the actual run of the fragment multipoles. This means we
        process the default parameters, override the parameters if they are
        specified, write all the necessary input files, and then run.

        Returns:
          (str): the name of the file containing the multipoles
        """
        from os.path import join, isfile
        from os import remove
        from inspect import getargspec
        from yaml import dump
        from futile.YamlIO import load

        # Convert the arguments of the function to a dictionary
        args, vargs, keywords, default = getargspec(self.fragment_multipoles)
        options = {a: d for a, d in zip(args, default)}

        # Use the logfile to determine the right matrix format
        format = log.log["lin_general"]["output_mat"]
        if format == 4:
            options["matrix_format"] = "parallel_mpi-native"
            for key in options:
                if ".txt" in options[key]:
                    options[key] = options[key].replace(".txt", ".mpi")

        # Replace the default directory with the appropriate one if it is
        # available
        if log.log["radical"]:
            data_dir = _get_datadir(log)
            for a, d in options.items():
                if a == "mpirun" or a == "action" or a == "matrix_format":
                    continue
                options[a] = join(data_dir, d)

        # Create the frag.yaml file from the provided system.
        if system:
            system.write_fragfile(options["fragment_file"], log)

        # Create the orbitals.yaml file from the provided orbitals.
        if orbitals is None:
            with open(options["orbital_file"], "w") as ofile:
                ofile.write(dump([-1]))

        if isfile(options["log_file"]):
            remove(options["log_file"])
        self.fragment_multipoles(**options)


        return load(options["log_file"], doc_lists=False)

    def set_fragment_multipoles(self, system, log):
        """
        Set the fragment multipoles of a system based on a run.

        The multipoles and purity values of the system are
        updated by this call.

        Args:
          system (System): instance of a System class, which defines the
          fragments we will use.
          log (Logfile): instance of a Logfile class
        """
        from os.path import join
        from os import remove

        mp_data = self._run_fragment_multipoles(log, system)

        # Update multipoles and purity values.
        mp_dict = mp_data["Orbital occupation"][0]["Fragment multipoles"]

        for frag, fdata in zip(system.values(), mp_dict):
            frag.purity_indicator = float(fdata["Purity indicator"])
            frag.q0 = [float(x) for x in fdata["q0"]]
            frag.q1 = [float(x) for x in fdata["q1"]]
            frag.q2 = [float(x) for x in fdata["q2"]]

    def compute_spillage(self, system, log, targetid):
        """
        Compute a measure of the spillage interaction between fragments.

        Args:
          system (System): instance of a System class, which defines the
            fragments we will use.
          log (Logfile): instance of a Logfile class
          targetid (str): which fragment to compute the spillage of all other
            fragments with. You can either set this or target.
        Result:
          (dict): for each fragment id, what the spillage value is.
        """
        from os.path import join, isfile
        from Spillage import MatrixMetadata, compute_spillbase, compute_spillage

        # Define the input files.
        spillage_array = []
        data_dir = _get_datadir(log)
        sfile = join(data_dir, "overlap_sparse.txt")
        hfile = join(data_dir, "hamiltonian_sparse.txt")

        # Convert to text format if necessary
        if log.log["lin_general"]["output_mat"] == 4:
            for fname in ["overlap_sparse.txt", "hamiltonian_sparse.txt"]:
                infile = join(data_dir, fname.replace(".txt", ".mpi"))
                outfile = join(data_dir, fname)
                self.convert_matrix_format(conversion="binary_to_bigdft",
                                           infile=infile, outfile=outfile)

        # Get the metadata
        metadatafile = join(data_dir, "sparsematrix_metadata.dat")
        metadata = MatrixMetadata(metadatafile)
        frag_indices = metadata.get_frag_indices(system)

        # Check whether to use python or bigpoly version
        # And then perform the computational of the inverse etc
        if isfile(self.bigpoly_command):
            sinvxh, sinvxh2 = compute_spillbase(sfile, hfile, bigpoly=True)
        else:
            # First convert to binary matrix format
            soutfile = join(data_dir, "overlap_sparse.ccs")
            houtfile = join(data_dir, "hamiltonian_sparse.ccs")
            self.convert_matrix_format(conversion="bigdft_to_ccs", infile=sfile,
                                       outfile=soutfile)
            self.convert_matrix_format(conversion="bigdft_to_ccs", infile=hfile,
                                       outfile=houtfile)
            sinvxh, sinvxh2 = compute_spillbase(soutfile, houtfile,
                                                bigpoly=False)

        # Compute the spillage array
        spillage_array = compute_spillage(
            sinvxh, sinvxh2, frag_indices, targetid)

        return spillage_array

    def plot_spillage(self, axs, spillvals):
        """
        Plot the spillage values.

        Args:
          axs: the axs we we should plot on.
          spillvals (dict): a dictionary mapping fragments to spillage values.
        """
        from Fragments import GetFragTuple
        axs.set_xlabel("Fragment", fontsize=12)
        axs.set_ylabel("Spillage Values", fontsize=12)
        axs.set_yscale("log")

        svals = []
        labels = []
        for id, val in spillvals.items():
            svals.append(val)
            labels.append(id)

        axs.set_xticks(range(len(labels)))
        axs.set_xticklabels(sorted(labels, key=lambda x: GetFragTuple(x)[1]),
                            rotation=90)
        axs.plot([v for _, v in sorted(zip(labels, svals),
                                       key=lambda x: GetFragTuple(x[0])[1])], 'x--')
