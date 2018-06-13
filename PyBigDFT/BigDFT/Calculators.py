"""file for BigDFT calculators"""

import os
import shutil
import yaml
from futile.Utils import write as safe_print
import BigDFT.Logfiles as Lf

class GIBinding():
    """Calculator for BigDFT from Gobject Introspection bindings"""

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
    """Define the calculator from the system"""

    def __init__(self, omp=None, mpi=1,mpi_run=None):
        """Initialize the SystemCalculator defining the Unix command, the number of openMP threads and MPI processes."""
        # Save variables for future use
        if omp == None:
            self.omp = os.environ.setdefault('OMP_NUM_THREADS','1')
        self.mpi = str(mpi)
        # Verify if $BIGDFT_ROOT is in the environment
        assert 'BIGDFT_ROOT' in os.environ
        if mpi_run == None:
            # Verify if $MPI_RUN is in the environment
            mpi_run = os.environ.setdefault('MPI_RUN','')
        if mpi_run != '':
            mpi_run = mpi_run + ' ' + self.mpi + ' '
        #Build the command setting the number of omp threads
        self.command =  'env OMP_NUM_THREADS='+self.omp+' '+mpi_run + os.environ['BIGDFT_ROOT']+'/bigdft'

    def run(self, name='', outdir='', run_name='', input=None, posinp=None, skip=False):
        """Run a calculation building the input file from a dictionary
           a Logfile instance is returned."""
        # Set the number of omp threads
        #os.environ['OMP_NUM_THREADS'] = self.omp
        # Creating the yaml input file from a dictionary or another file
        if len(name) > 0:
            input_file = "%s.yaml" % name
            logname = "log-%s.yaml" % name
        else:
            input_file = "input.yaml" #default name
            logname = "log.yaml"
        if isinstance(input,dict):
            #Creating the yaml input file
            yaml.dump(input,stream=open(input_file,"w"),default_flow_style=False)
            safe_print('Creating from a dictionary the yaml input file "%s"' % input_file)
        elif isinstance(input,str):
             #Check if the file input does exist
            assert os.path.exists(input)
            if input != input_file:
                shutil.copyfile(input,input_file)
                safe_print('Copying from "%s" the yaml input file "%s"' % (input,input_file))
        # Copying the posinp input file if needed
        if posinp != None:
            #Check if the file does exist
            assert os.path.exists(posinp)
            if len(name) > 0:
                fn,fn_ext = os.path.splitext(posinp)
                posinp_file = name+fn_ext
                if posinp != posinp_file:
                    shutil.copyfile(posinp,posinp_file)
                    safe_print('Copying from "%s" the posinp file "%s"' % (posinp,posinp_file))
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
