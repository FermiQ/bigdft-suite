"""file for BigDFT calculators"""

import os
from gi.repository import BigDFT
from futile.Utils import write as safe_print


class GIBinding():
    """Calculator for BigDFT from Gobject Introspection bindings"""

    def __init__(self):
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

    def __init__(self, omp, mpi):
        # Save variables for future use
        self.omp = str(omp)
        self.mpi = str(mpi)
        # Verify if $BIGDFT_ROOT is in the environment
        assert 'BIGDFT_ROOT' in os.environ
        self.command = 'mpirun -np ' + self.mpi + ' $BIGDFT_ROOT/bigdft'

    def run(self, name='', outdir='', run_name='', skip=False):
        # Set the number of omp threads
        os.environ['OMP_NUM_THREADS'] = self.omp
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
