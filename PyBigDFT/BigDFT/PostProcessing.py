"""A module for post processing BigDFT calculations.

"""


def _system_command(command, options, outfile=None):
    """
    Run the command as ``os.system (command + options)``

    Todo:
       remove outfile, this should be controlled by the script.
    Args:
       command (str):
       options (str):
    """
    from subprocess import call

    command_str = command + " " + options

    if outfile:
        with open(outfile, "w") as ofile:
            call(command_str, shell=True, stdout=ofile)
    else:
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
    return join(log.srcdir, "data-" + log.log["radical"])


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
            nspec["args"] = [
                {"mpirun": {"default": "" + mpi_run + ""}}] + nspec["args"]
            nspec["args"] = [
                {"action": {"default": "" + naction + ""}}] + nspec["args"]
            my_action = get_python_function(partial(self.__invoke_command,
                                                    self.bigdft_tool_command),
                                            action, nspec)
            setattr(self, action, my_action)

    def __invoke_command(self, command, **kwargs):
        from futile.Utils import option_line_generator
        _system_command(command, option_line_generator(**kwargs), self.outfile)

    def _run_fragment_multipoles(self, log, system=None, orbitals=None):
        from os.path import join
        from inspect import getargspec
        import yaml

        # Convert the arguments of the function to a dictionary
        args, vargs, keywords, default = getargspec(self.fragment_multipoles)
        options = {a: d for a, d in zip(args, default)}

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
            system.write_fragfile(filename=options["fragment_file"])

        # Create the orbitals.yaml file from the provided orbitals.
        if orbitals is None:
            with open(options["orbital_file"], "w") as ofile:
                ofile.write(yaml.dump([-1]))

        self.fragment_multipoles(**options)

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
        from futile.YamlIO import load

        self.outfile = join(log.srcdir, "mp.yaml")
        self._run_fragment_multipoles(log, system)

        # Update multipoles and purity values.
        mp_data = load(self.outfile, doc_lists=False)
        mp_dict = mp_data["Orbital occupation"][0]["Fragment multipoles"]

        for frag, fdata in zip(system.fragments, mp_dict):
            frag.purity_indicator = float(fdata["Purity indicator"])
            q0 = [float(x) for x in fdata["q0"]]
            q1 = [float(x) for x in fdata["q1"]]
            q2 = [float(x) for x in fdata["q2"]]
            frag.set_fragment_multipoles(q0, q1, q2)

    def compute_spillage(self, system, log, target):
        """
        Compute a measure of the spillage interaction between fragments.

        Args:
          system (System): instance of a System class, which defines the
          fragments we will use.
          log (Logfile): instance of a Logfile class
          target (int): which fragment to compute the spillage of all other
            fragments with.
        """
        from os.path import join
        from Spillage import compute_spillbase, process_metadata, \
            compute_spillage

        # Define the input files.
        self.outfile = join(log.srcdir, "spillage.yaml")
        spillage_array = []
        data_dir = _get_datadir(log)
        sfile = join(data_dir, "overlap_sparse.txt")
        soutfile = join(data_dir, "overlap_sparse.ccs")
        hfile = join(data_dir, "hamiltonian_sparse.txt")
        houtfile = join(data_dir, "hamiltonian_sparse.ccs")

        # First convert to CCS matrix format
        self.convert_matrix_format(conversion="bigdft_to_ccs", infile=sfile,
                                   outfile=soutfile)
        self.convert_matrix_format(conversion="bigdft_to_ccs", infile=hfile,
                                   outfile=houtfile)

        # Read in the matrices from file
        sinvxh, sinvxh2 = compute_spillbase(soutfile, houtfile)

        # Next we compute the indices associated with the target fragment.
        metadatafile = join(data_dir, "sparsematrix_metadata.dat")
        frag_indices = process_metadata(metadatafile, system)

        spillage_array = compute_spillage(
            sinvxh, sinvxh2, frag_indices, target)

        return spillage_array
