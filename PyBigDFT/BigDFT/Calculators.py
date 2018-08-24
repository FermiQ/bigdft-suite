"""
This module defines some classes to perform a calculation using BigDFT
using binding (GIBinding) or using system call (SystemCalculator).

"""

##In our case for the class SystemCalculator which uses system calls:
##* We define posinp (equivalent of Atoms)
##* We have a python dictionary for the parameter
##* We define a calculator (equivalent of BFGS which is an Optimizer (a method to optimize))
##Then we perform the method run.
##
##For the class GIBinding using the Gobject Introspection bindings, two methods set and update are added.
##
##The goal is to have a light Calculator almost compatible with ASE (Atomic Simulation environment, see https://gitlab.com/ase/ase)
##.. todo::
##    In a future we add our method to ASE which is at a higher level (workflow of simulations).
##:Example:
##   >>> from ase import Atoms
##   >>> from ase.optimize import BFGS
##   >>>  from ase.calculators.nwchem import NWChem
##   >>>  from ase.io import write
##   >>>  h2 = Atoms('H2',
##   >>>             positions=[[0, 0, 0],
##   >>>                        [0, 0, 0.7]])
##   >>>  h2.calc = NWChem(xc='PBE')
##   >>>  opt = BFGS(h2, trajectory='h2.traj')
##   >>>  opt.run(fmax=0.02)
##   >>>  BFGS:   0  19:10:49    -31.435229     2.2691
##   >>>  BFGS:   1  19:10:50    -31.490773     0.3740
##   >>>  BFGS:   2  19:10:50    -31.492791     0.0630
##   >>>  BFGS:   3  19:10:51    -31.492848     0.0023
##   >>>  write('H2.xyz', h2)
##   >>>  h2.get_potential_energy()  # ASE's units are eV and Ang
##   >>>  -31.492847800329216
##

import os
import shutil
import yaml
import copy
from futile.Utils import write as safe_print
import BigDFT.Logfiles as Lf

import collections

def dict_merge(dct, merge_dct):
    """ Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    
    From [angstwad/dict-merge.py](https://gist.github.com/angstwad/bf22d1822c38a92ec0a9)
    """
    for k, v in merge_dct.iteritems():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]


def get_debugfile_date():
    """
    Get the information about the debug time of the last file in the current directory
    """
    from futile.Utils import file_time
    return file_time(os.path.join('debug','bigdft-err-0.yaml'))


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
            dict_merge(var,{'dft': {'inputpsiid': 1}})
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


class Runner():
    """
    Define a generic class parent of SystemCalculator and Workflow which
    defines a run method and functions to handle some global options.
    The method run can be used many times.
    All arguments in the __init__ call is stored as global options.
    For each run, these global options are used updated by the arguments of the run call
    all stored in a run_options dictionary which can be used by the run call.
    """
    def __init__(self,**kwargs):
        """All arguments for the runs are saved in a private dictionary of global options"""
        import copy
        self._global_options=copy.deepcopy(kwargs)

    def global_options(self):
        """Get all global options"""
        return self._global_options

    def get_global_option(self,key):
        """Get one key in global options"""
        return self._global_options[key]

    def update_global_options(self,**kwargs):
        """Update the global options"""
        self._global_options.update(kwargs)

    def unset_global_option(self,key):
        """Remove a given global option"""
        self._global_option.pop(key)

    def _run_options(self,**kwargs):
        """Create a local dictionary for a specific run."""
        import copy
        #First deepcopy from global_options and update from kwargs (warning: a dictionary is not update)
        self.run_options=copy.deepcopy(self._global_options)
        self.run_options.update(kwargs)

    def run(self,**kwargs):
        """Implement a run method by default (do nothing except updating run options)"""
        self._run_options(**kwargs)
        #Do nothing
        return self.post_processing()
    
    def post_processing(self):
        """Post-processing"""
        return None

class SystemCalculator(Runner):
    """
    Define a BigDFT calculator.

    :param int omp: number of OpenMP threads
      It defaults to the $OMP_NUM_THREADS variable in the environment (if present), otherwise it fixes the run to 1 thread
    :param str mpi_run: define the MPI command to be used.
      It defaults to the value $BIGDFT_MPIRUN of the environment, if present
      When using this calculator into a job submission script, the value of
      $BIGDFT_MPIRUN variable may be set appropriately to launch the bigdft executable.
    :param bool skip: if True, do not run the calculation if the corresponding logfile exists.
    :param bool verbose: default verbosity (2 levels true or False)      
    :param bool dry_run: check the input, estimate the memory but do not perform the calculation
    :param int dry_mpi: Number of MPI processes for the estimation of the memory when dry_run is True

    Check if the environment variable $BIGDFT_ROOT is defined.
    This is a signal that the environment has been properly set prior to the evaluation of the python command.
    Use also two environment variables:
        * OMP_NUM_THREADS to set the number of OMP_NUM_THREADS
        * BIGDFT_MPIRUN to define the MPI execution command.

    :Example:
        >>> inpdict = { 'dft': { 'ixc': 'LDA' }} #a simple input file
        >>> study = SystemCalculator(omp=1)
        >>> logf = study.run(name="test",input=inpdict)
        Executing command:  $BIGDFT_MPIRUN <path_to_$BIGDFT_ROOT>/bigdft test
    """

    import os,shutil
    def __init__(self, 
                 omp=os.environ.get('OMP_NUM_THREADS','1'),
                 mpi_run=os.environ.get('BIGDFT_MPIRUN',''),
                 dry_run=False,skip=False,verbose=True):
        """
        Initialize the SystemCalculator defining the Unix command, the number of openMP threads and MPI processes.
        """
        #Use the initialization from the Runner class (so all options inside __global_options)
        Runner.__init__(self,omp=omp,mpi_run=mpi_run,dry_run=dry_run,skip=skip,verbose=verbose)
        assert 'BIGDFT_ROOT' in os.environ      # Verify if $BIGDFT_ROOT is in the environment
        #Build the command setting the number of omp threads
        self.command = (self._global_options['mpi_run'] + ' ' + os.environ['BIGDFT_ROOT']+'/bigdft').strip()
        safe_print('Initialize a Calculator with OMP_NUM_THREADS=%s and command %s' % (self._global_options['omp'],self.command) )

    def run(self, **kwargs):
    #def run(self, name='', outdir='', run_name='', input={}, posinp=None,**kwargs):
        """
        Run a calculation building the input file from a dictionary.

        :param str name: naming scheme of the run i.e. <name>.yaml is the input file and log-<name>.yaml the output one.
           Data will then be written in the directory `data-<name>.yaml`, unless the "radical" keyword is specified in the input dictionary.
        :param str run_dir: specify the directory where bigdft will be executed (the input and log file will be created in it)
                            it must be a simple 
        :param str outdir: specify the output directory for all data coming from bigdft (parameter of bigdft)
        :param str run_name: File containing the list of the run ids which have to be launched independently 
                             (list in yaml format). The option runs-file is not compatible with the name option.
        :param input: give the input parameters (a dictionary or a list of dictionary)
        :type input: dict
        :param posinp: indicate the posinp file (atomic position file). 
        :type posinp: filename
        :return: a Logfile instance is returned. It returns None if an error occurred
        :rtype: Logfile

        .. todo::
           
           Set the return value of run in the case of a run_file. It should be a list of Logfile classes

        """
        #Create a temporary dictionary of options
        self._run_options(**kwargs)
        #Give the default values in case of unset_global_options are called
        name = self.run_options.get('name','')
        outdir = self.run_options.get('outdir','')
        run_dir = self.run_options.get('run_dir','.')
        run_name = self.run_options.get('run_name','')
        posinp = self.run_options.get('posinp',None)
        inp = self.run_options.get('input',{})
        verbose = self.run_options['verbose']
        dry_run = self.run_options['dry_run']
        #Restrict run_dir to a sub-directory
        if  ("/" in run_dir or run_dir == ".."):
            raise ValueError("run_dir '%s' where bigdft is executed must be a sub-directory" % run_dir)
        #Create the run_dir if not exist
        if not os.path.exists(run_dir):
            #Creation of the sub-directory run_dir
            os.mkdir(run_dir)
            if verbose: safe_print("Create the sub-directory '%s'" % run_dir)
        #Create the input file (deepcopy because we modify it!)
        local_input=copy.deepcopy(inp)
        #Check posinp file
        if posinp != None:
            #Check if the file does exist
            if not os.path.isfile(posinp):
                raise ValueError("posinp: The atomic position file '%s' does not exist" % posinp)
            #Add into the dictionary a posinp key
            if run_dir == '.':
                local_input['posinp'] = posinp
            else:
                #Copy the posinp if not identical
                cp_posinp = "%s/%s" % (run_dir,posinp)
                if not (os.path.isfile(cp_posinp) and os.stat(posinp) != os.stat(cp_posinp)):
                    shutil.copy2(posinp,run_dir)
                    if verbose: safe_print("Copy the posinp file '%s' into '%s'" % (posinp,run_dir))
                local_input['posinp'] = os.path.basename(posinp)
        if len(name) > 0:
            input_file = "%s/%s.yaml" % (run_dir,name)
            logname = "%s/log-%s.yaml" % (run_dir,name)
        else:
            input_file = "input.yaml" #default name
            logname = "log.yaml"
        #check if the debug file will be updated (case of erroneous run)
        timedbg=get_debugfile_date()
        #Creating the yaml input file
        open(input_file,"w").write(yaml.dump(local_input,default_flow_style=None))
        if verbose: safe_print('Creating the yaml input file "%s"' % input_file)
        #Check if it is a dry run
        if dry_run:
            #Use bigdft-tool
            command = self.run_options['mpi_run'] + ' ' + os.environ['BIGDFT_ROOT']+'/bigdft-tool -l'
            if dry_run:
                command += ' -n ' + str(dry_run)
            if len(name) > 0:
                command += ' --name='+name
        else:   
            # Adjust the command line with options
            command = self.command
            if len(name) > 0:
                command += ' -n ' + name
            if len(run_name) > 0:
                command += ' -r ' + run_name
            if len(outdir) > 0:
                command += ' -d ' + outdir
                #Change logname
                logname = outdir+"/"+logname
            if self.run_options['skip']:
                command += ' -s Yes'
        if verbose: safe_print('Executing command: ', command)
        # Set the number of omp threads
        os.environ['OMP_NUM_THREADS'] = self.run_options['omp']
        #Run the command
        os.system("cd "+run_dir+"; "+command)
        #verify that no debug file has been created
        if get_debugfile_date() > timedbg :
            if verbose: 
                safe_print("ERROR: some problem occured during the execution of the command, check the 'debug/' directory and the logfile")
            return None
        return self.post_processing(logname)
    
    def post_processing(self,logname):
        """
        Check the existence and the log file and return an instance logfile
        valid only without run_name
        """
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
