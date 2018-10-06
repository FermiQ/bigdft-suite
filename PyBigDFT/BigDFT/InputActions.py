"""Actions to define on the Input parameters.

This module defines some of the most common actions that a BigDFT user might like to 
perform on the input file. Such module therefore set some of the keys of the input 
dictionary to the values needed to perform the operations.
Users might also inspire to the actions performed in order to customize the runs in a different way.
All the functions of this module have as first argument ``inp``, the dictionary of the input parameters.

Many other actions are available in BigDFT code. This module only regroups the most common.

.. autosummary::

   set_xc
   set_atomic_positions
   optimize_geometry
   spin_polarize
   charge
   set_random_inputguess
   write_orbitals_on_disk
   read_orbitals_from_disk
   write_density_on_disk
   change_data_directory
   add_empty_SCF_orbitals
   set_electronic_temperature


"""

from futile.Utils import dict_set

__set__ = dict_set
"""func: Action function.

This is the pointer to the set function, useful to modify the action with the undo method

"""

def undo(inp,*subfields):
    """
    Eliminate the last item of the subfields as provided to dict_set
    """
    from futile.Utils import push_path
    keys=subfields[:-1]
    tmp,k=push_path(inp,*keys)
    tmp.pop(k)


def remove(action):
    """Remove input action.
    
    Remove an action from the input file, thereby restoring the **default** value, as if the action were not specified.

    Args:
       action (func): one of the actions of this module. It does not need to be specified before, in which case it produces no effect.
    
    Example:
       >>> from Calculators import SystemCalculator as C
       >>> code=C()
       >>> inp={}
       >>> set_xc(inp,'PBE')
       >>> write_orbitals_on_disk(inp)
       >>> log=code.run(input=inp) # perform calculations
       >>> remove(write_orbitals_on_disk) #remove the action
       >>> read_orbitals_from_disk(inp)
       >>> log2=code.run(input=inp) #this will restart from the previous one
    """
    global __set__
    __set__ = undo
    action(inp)
    __set__ = dict_set

def spin_polarize(inp,mpol=1):
    """
    Add a collinear spin polarization to the system.

    Arguments:
       mpol (int): spin polarization in Bohr magneton units.
    """
    __set__(inp,'dft','nspin',2)
    __set__(inp,'dft','mpol',mpol)

def charge(inp,charge=-1):
    """
    Charge the system

    Arguments:
        charge (int,float): value of the charge in units of *e* (the electron has charge -1). Also accept floating point numbers.
    """
    __set__(inp,'dft','charge',charge)

def make_cation(inp):
    """
    Charge the system by removing one electron. Assume that the original system is closed shell, thus polarize.
    """
    charge(inp,charge=1)
    spin_polarize(inp,mpol=1)

def add_empty_SCF_orbitals(inp,norbs=10):
    """
    Insert ``norbs`` empty orbitals in the SCF procedure

    Args:
       norbs (int): Number of empty orbitals
    """
    __set__(inp,'mix','norbsempty',norbs)

def write_orbitals_on_disk(inp,format='binary'):
    """
    Set the code to write the orbitals on disk in the provided format

    Args:
      format (str): The format to write the orbitals with. Accepts the strings:

         * 'binary'
         * 'text'
         * 'text_with_densities'
         * 'text_with_cube'
         * 'etsf' (requires etsf-io enabled)
    """
    __set__(inp,'output','orbitals',format)

def set_atomic_positions(inp,posinp=None):
    """
    Insert the atomic positions as a part of the input dictionary
    """
    __set__(inp,'posinp',posinp)

def read_orbitals_from_disk(inp):
    """
    Read the orbitals from data directory, if available
    """
    __set__(inp,'dft','inputpsiid',2)

def set_random_inputguess(inp):
    """
    Input orbitals are initialized as random coefficients
    """
    __set__(inp,'dft','inputpsiid',-2)

def set_electronic_temperature(inp,kT=1.e-3,T=0):
    """
    Define the electronic temperature, in AU (``kT``) or K (``T``)
    """
    TtokT=8.617343e-5/27.21138505
    tel= TtoKT*T if T != 0 else kT
    __set__(inp,'mix','tel',tel)
    
def optimize_geometry(inp,method='FIRE',nsteps=50):
    """
    Optimize the geometry of the system

    Args:
       nsteps (int): maximum number of atomic steps.
       method (str): Geometry optimizer. Available keys:
          * SDCG:   A combination of Steepest Descent and Conjugate Gradient
          * VSSD:   Variable Stepsize Steepest Descent method
          * LBFGS:  Limited-memory BFGS
          * BFGS:   Broyden-Fletcher-Goldfarb-Shanno
          * PBFGS:  Same as BFGS with an initial Hessian obtained from a force field
          * DIIS:   Direct inversion of iterative subspace
          * FIRE:   Fast Inertial Relaxation Engine as described by Bitzek et al.
          * SBFGS:  SQNM minimizer, keyword deprecated, will be replaced by SQNM in future release
          * SQNM:   Stabilized quasi-Newton minimzer
    """
    __set__(inp,'geopt','method',method)
    if nsteps !=50: __set__(inp,'geopt','ncount_cluster_x',nsteps)

def set_xc(inp,xc='PBE'):
    """
    Set the exchange and correlation approximation
    
    Args:
       xc (str): the Acronym of the XC approximation

    Todo:
       Insert the XC codes corresponding to ``libXC`` conventions
    """
    __set__(inp,'dft','ixc',xc)

def write_density_on_disk(inp):
    """
    Write the charge density on the disk after the last SCF convergence
    """
    __set__(inp,'dft','output_denspot',21)


def change_data_directory(inp,name=''):
    """
    Modify the name of the ``data-`` directory.
    Useful to grab the orbitals from another directory than the run name
    """
    __set__(inp,'radical',name)


