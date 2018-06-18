"""
This module defines some classes to perform a calculation using BigDFT
using binding (GIBinding) or using system call (SystemCalculator).

The goal is to have a light Calculator almost compatible with ASE (Atomic Simulation environment, see https://gitlab.com/ase/ase)
.. todo::
    In a future we add our method to ASE which is at a higher level (workflow of simulations).

For the class GIBinding using the Gobject Introspection bindings, two methods set and update are added.
"""

import os
import shutil
import yaml
from futile.Utils import write as safe_print
import BigDFT.Logfiles as Lf

class GIBinding():
    """
    Calculator for BigDFT from Gobject Introspection bindings.
    """

    def __init__(self):
        #Import bindings about BigDFT (if the bindings are not generated, do not work at all)
        from gi.repository import BigDFT
        self.runObj = -1
        # MPI initialisation
        (ierr, self.iproc, self.nproc, igroup, ngroup) = BigDFT.lib_init(0)
        self.runObj = None

    def update(self, inputfile):
        # If the inputpsiid is not present in the inputfile
        # assumes that the user wants to do a restart
        if "dft" in inputfile and "inputpsiid" in inputfile["dft"]:
            var = inputfile
        else:
            var = inputfile.copy()
            var.update({'dft': {'inputpsiid': 1}})
        self.runObj.update(BigDFT.Dict(var))

    def run(self):
        self.out = self.runObj.calculate(self.iproc, self.nproc)
        return self.out

    def set(self, inputfile=None):
        if inputfile is None:
            var = {}
        else:
            var = inputfile
        # Free memory first
        self.out = None
        self.runObj = None
        self.runObj = BigDFT.Run.new_from_dict(BigDFT.Dict(var))

    def __del__(self):
        if self.runObj == -1:
            return
        # MPI finalisation.
        self.out = None
        self.runObj = None
        BigDFT.lib_finalize()


class SystemCalculator():
    """
    Define a BigDFT calculator.

    :param int omp: number of OpenMP threads (if none set to $OMP_NUM_THREADS)
    :param int mpi: number of MPI processes (if none set to 1)
    :param str mpi_run: define the MPI command to be used
    :param bool skip: if True, do not run the calculation if the corresponding logfile exists.
    :param bool verbose: default verbosity (2 levels true or False)

    Check is the environment variable $BIGDFT_ROOT is defined.
    Use also two environment variables:
        * OMP_NUM_THREADS to set the number of OMP_NUM_THREADS
        * BIGDFT_ROOT_MPIRUN to define the MPI execution command.

    :Example:
        >>> dico = { 'dft': { 'ixc': 'LDA' }}
        >>> study = SystemCalculator(omp=1)
        >>> logf = study.run(name="test",input=dico)
    """

    def __init__(self, omp=None, mpi=1,mpi_run=None,skip=False,verbose=True):
        """
        Initialize the SystemCalculator defining the Unix command, the number of openMP threads and MPI processes.
        """
        # Save variables for future use
        if omp == None:
            self.omp = os.environ.setdefault('OMP_NUM_THREADS','1')
        else:
            self.omp = str(omp)
        self.mpi = str(mpi)
        #Default value for all runs
        self.skip = skip
        #Default values for all runs
        self.verbose = verbose
        # Verify if $BIGDFT_ROOT is in the environment
        assert 'BIGDFT_ROOT' in os.environ
        if mpi_run == None:
            # Verify if $BIGDFT_MPIRUN is in the environment
            mpi_run = os.environ.setdefault('BIGDFT_MPIRUN','')
        if mpi_run != '':
            mpi_run = mpi_run + ' ' + self.mpi + ' '
        #Build the command setting the number of omp threads
        self.command =  mpi_run + os.environ['BIGDFT_ROOT']+'/bigdft'
        safe_print('Initialize a Calculator with OMP_NUM_THREADS=%s and command %s' % (self.omp,self.command) )

    def run(self, name='', outdir='', run_name='', input=None, posinp=None, skip=None,verbose=None):
        """
        Run a calculation building the input file from a dictionary.

        :param str name: naming schme of the run i.e. <name>.yaml is the input file and log-<name>.yaml the output one.
        :param str outdir: specify the output directory
        :param str run_name: File containing the list of the run ids which have to be launched independently 
                             (list in yaml format). The option runs-file is not compatible with the name option.
        :param input: give the input parameters
        :type input: dict or filename
        :param bool skip: avoid to rerun the calculation. Check if the file log-<name> already exists.
        :param bool verbose: verbosity
        :param posinp: indicate the posinp file (atomic position file)
        :type posinp: filename
        :return: a Logfile instance is returned.
        :rtype: Logfile
        """
        # Set the number of omp threads
        os.environ['OMP_NUM_THREADS'] = self.omp
        # Creating the yaml input file from a dictionary or another file
        if skip == None: skip = self.skip
        if verbose == None: verbose = self.verbose
        if len(name) > 0:
            input_file = "%s.yaml" % name
            logname = "log-%s.yaml" % name
        else:
            input_file = "input.yaml" #default name
            logname = "log.yaml"
        # Copying the posinp input file if needed
        if posinp != None:
            #Check if the file does exist
            assert os.path.exists(posinp)
            if isinstance(input,dict):
                #Add into the dictionary a posinp key
                input['posinp'] = posinp
                if verbose: safe_print('Add a "posinp" key with the value "%s" into the dictionary input' % posinp)
            else:
                if len(name) > 0:
                    fn,fn_ext = os.path.splitext(posinp)
                    posinp_file = name+fn_ext
                if posinp != posinp_file:
                    shutil.copyfile(posinp,posinp_file)
                    if verbose: safe_print('Copying from "%s" the posinp file "%s"' % (posinp,posinp_file))
        if isinstance(input,dict):
            #Creating the yaml input file
            yaml.dump(input,stream=open(input_file,"w"),default_flow_style=False)
            if verbose: safe_print('Creating from a dictionary the yaml input file "%s"' % input_file)
        elif isinstance(input,str):
             #Check if the file input does exist
            assert os.path.exists(input)
            if input != input_file:
                shutil.copyfile(input,input_file)
                if verbose: safe_print('Copying from "%s" the yaml input file "%s"' % (input,input_file))
        # Adjust the command line with options
        command = self.command
        if len(name) > 0:
            command += ' -n '+name
        if len(run_name) > 0:
            command += ' -r '+run_name
        if len(outdir) > 0:
            command += ' -d '+outdir
        if skip:
            command += ' -s Yes'
        safe_print('Executing command: ', command)
        os.system(command)
        #Check the existence and the log file and return an instance logfile
        if os.path.exists(logname):
            return Lf.Logfile(logname)
        else:
            return None


# Test the calculators
if __name__ == '__main__':

    import yaml
    import matplotlib.pyplot as plt

    basicinput = """
#mode: {method: lj}
logfile: No
dft: { ixc: HF, nspin: 2}
posinp:
   positions:
   - {Be : [0.0, 0.0, 0.0]}#, IGSpin: -1}
   - {Be : [0.0, 0.0, 1.0]}#, IGSpin: 1}
#   properties: {format: yaml}
ig_occupation:
   Atom 1: {2s: {up: 1.0, down: 0.9}, 2p: {up: 0.0, down: 0.2} }
   Atom 2: {2s: {up: 0.9, down: 1.0}, 2p: {up: 0.2, down: 0.0} }

psppar.Be: {Pseudopotential XC: 11}
"""

    # Initialize the calculator
    study = GIBinding()
    inp = yaml.load(basicinput)
    study.set(inp)

    # Perform the first calculation
    out = study.run()
    safe_print('starting energy', out.eKS)
    energy = [out.eKS]
    pos = [1.0]

    # Perform a dissociation curve
    for i in range(10):
        sh = float(i+1) * 0.02
        inp['posinp']['positions'][-1]['Be'][2] += sh
        study.update(inp)
        out = study.run()
        energy.append(out.eKS)
        pos.append(pos[-1]+sh)
        if study.iproc == 0:
            safe_print('iter', i, 'shift', sh, 'energy', out.eKS)
    out = None
    safe_print('End of the calculations')

    # Plot the dissociation curve
    if study.iproc == 0:
        plt.plot(pos, energy)
        plt.show()
