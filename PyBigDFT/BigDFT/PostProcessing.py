"""A module for post processing BigDFT calculations.

"""


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
        self.bigdft_tool_command = join("$BIGDFT_ROOT", "bigdft-tool")
        self.utilities_command = mpi_run + join("$BIGDFT_ROOT", "utilities")
        self.memguess_command = join("$BIGDFT_ROOT", "memguess")
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
        _system_command(command, option_line_generator(**kwargs))

    def _run_fragment_multipoles(self, logfile):
        from os.path import join
        from inspect import getargspec

        # Convert the arguments of the function to a dictionary
        args, vargs, keywords, default = getargspec(self.fragment_multipoles)
        options = {a: d for a, d in zip(args, default)}

        # Replace the default directory with the appropriate one if it is
        # available
        if logfile.log["radical"]:
            data_dir = join(logfile.srcdir, "data-" + logfile.log["radical"])
            for a, d in options.items():
                if a == "mpirun" or a == "action" or a == "matrix_format":
                    continue
                options[a] = join(data_dir, d.lstrip("data/"))
                
        self.fragment_multipoles(**options)

    def get_fragment_multipoles(self, logfile):
        """
        Extract the fragment multipoles coming from a run specified by the
        passed logfile.
        """
        self._run_fragment_multipoles(logfile)
