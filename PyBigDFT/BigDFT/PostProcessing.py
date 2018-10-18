"""A module for post processing.

"""

database = """
fragment_multipoles:
  help: |
    Compute the multipoles and purity value of a system given its
    partitioning into fragments.
  args:
  - matrix_format:
      shorthelp: format of the sparse matrices.
      default: serial_text
  - metadata_file:
      shorthelp: input file with the sparse matrix metadata
      default: data/sparsematrix_metadata.dat
  - fragment_file:
      shorthelp: file with the fragment to be analyzed
      default: data/fragment.dat
  - overlap_file:
      shorthelp: input file with the overlap matrix
      default: data/overlap_sparse.txt
  - kernel_file:
      shorthelp: input file with the coefficients
      default: data/density_kernel_sparse.txt
  - kernel_matmul_file:
      shorthelp: kernel matrix multiplication sparsity pattern file.
      default: data/density_kernel_sparse_matmul.txt
  - multipole_matrix_0_0:
      default: data/mpmat_0_0.txt
  - multipole_matrix_1_0:
      default: data/mpmat_1_-1.txt
  - multipole_matrix_1_1:
      default: data/mpmat_1_0.txt
  - multipole_matrix_1_2:
      default: data/mpmat_1_1.txt
  - multipole_matrix_2_0:
      default: data/mpmat_2_-2.txt
  - multipole_matrix_2_1:
      default: data/mpmat_2_-1.txt
  - multipole_matrix_2_2:
      default: data/mpmat_2_0.txt
  - multipole_matrix_2_3:
      default: data/mpmat_2_1.txt
  - multipole_matrix_2_4:
      default: data/mpmat_2_2.txt
  - orbital_file:
      shorthelp: input file specifying which orbitals to include
      default: data/orbitals.yaml
  - coeff_file:
      shorthelp: input file with the coefficients
      default: data/coeff.txt
"""
# from Futile import YamlIO as _Y
# db = _Y.load(stream=database)


# def function_generator(spec):
#     """
#     Method for generating a function according to the arguments provided
#     by the ``spec`` dictionary.
#
#     Args:
#       spec (dict): specification of the function provided as a dictionary
#         of the keyword arguments that the function will have. The
#         key ``help`` of the spec dictionary is passed as a docstring of
#         the function. The other keys are considered as arguments and
#         should contain the default values and help information provided
#         in the ``default`` and the ``shorthelp`` keys respectively.
#
#         Example:
#           {"help": "Compute the multipoles and purity value.",
#            "matrix_format": {
#            "shorthelp": "format of the sparse matrices.",
#            "default": "parallel_mpi-native"},
#            "metadata_file": {
#            "shorthelp": "input file with the sparse matrix metadata",
#            "default": "sparsematrix_metadata.dat"}}
#
#     """
#     docstring = spec.pop("help")
#     if "shorthelp" in spec:
#         spec.pop("shorthelp")
#     fun_args = {key: spec[key]["default"] for key in spec}
#
#     def fun(**fun_args):
#         command = ""
#         for key in fun_args:
#             command += "--" + key + "=" + fun_args[key] + " "
#         return command
#     fun.__doc__ = docstring
#     return fun

def _system_command(command, options):
    """
    Run the command as ``os.system (command + options)``

    Args:
       command (str):
       options (str):
    """
    from os import system

    command_str = command + " " + options
    system(command_str + " > check.txt")


class BigDFTool(object):
    """

    This module defines the actions that can be performed using the
    ``bigdft-tool`` python script. Such a script sometimes invokes the
    ``memguess`` executable, which is executed serially. Sometimes the
    ``utilities`` main program is necessary to take care of more time
    consuming operations which require parallelism.

    Args:
      omp (int): number of OpenMP threads.
        It defaults to the $OMP_NUM_THREADS variable in the environment, if present, otherwise it fixes the run to 1 thread.
      mpi_run (str): define the MPI command to be used.
        It defaults to the value $BIGDFT_MPIRUN of the environment, if present.
        When using this calculator into a job submission script, the value of
        $BIGDFT_MPIRUN variable may be set appropriately to launch the bigdft executable.
    """

    def __init__(self, omp=1, mpi_run=""):
        from futile.YamlArgparse import get_python_function
        from futile import YamlIO
        import os
        from functools import partial
        from copy import deepcopy

        self.bigdft_tool_command = "$BIGDFT_ROOT/bigdft-tool"
        self.utilities_command = mpi_run + "$BIGDFT_ROOT/utilities"
        self.memguess_command = "$BIGDFT_ROOT/memguess"
        os.environ['OMP_NUM_THREADS'] = str(omp)

        db = YamlIO.load(stream=database, doc_lists=False)
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
        _system_command(command, option_line_generator(**kwargs))

    def _run_fragment_multipoles(self, logfile):
        from os.path import join
        from inspect import getargspec

        args, vargs, keywords, default = getargspec(self.fragment_multipoles)
        options = {a:d for a,d in zip(args, default)}

        self.fragment_multipoles(**options)

    def get_fragment_multipoles(self, logfile):
        """
        Extract the fragment multipoles coming from a run specified by
        """
        self._run_fragment_multipoles(logfile)
