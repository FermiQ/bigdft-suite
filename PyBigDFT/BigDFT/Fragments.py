"""
This module is related to the usage of BigDFT with Fragment-related Quantities.
Input as well as Logfiles might be processed with the classes and methods provided by it.

"""

from futile.Utils import write as safe_print
try:
    from collections.abc import MutableMapping, MutableSequence
except:
    from collections import MutableMapping, MutableSequence

#: Conversion between Atomic Units and Bohr
AU_to_A = 0.52917721092
#: Conversion between Debye and Atomic units
Debye_to_AU = 0.393430307

MULTIPOLE_ANALYSIS_KEYS = ['q0', 'q1', 'q2', 'sigma']
PROTECTED_KEYS = MULTIPOLE_ANALYSIS_KEYS + ["frag"]

class Fragment(MutableSequence):
    def __init__(self, atomlist=None, xyzfile=None):
        from Atom import Atom
        self.atoms = []

        # insert atoms.
        if atomlist and isinstance(atomlist, list):
            for atom in atomlist:
                self.append(Atom(atom))
        elif xyzfile and isinstance(xyzfile, XYZReader):
            xyzfile.open()
            for line in xyzfile:
                self.append(Atom(line))
            xyzfile.close()

        # Values
        self.purity_indicator = None
        self.q0_value = None
        self.q1_value = None
        self.q2_value = None
        self.qcharge = None

    def __len__(self):
        return len(self.atoms)

    def __delitem__(self, index):
        self.atoms.__delitem__(index)

    def insert(self, index, value):
        from Atom import Atom
        self.atoms.insert(index, Atom(value))

    def __setitem__(self, index, value):
        from Atom import Atom
        self.atoms.__setitem__(index, Atom(value))

    def __getitem__(self, index):
        return self.atoms.__getitem__(index)

    def d0(self, center=None):
        pass

    def d1(self, center=None):
        pass

if __name__ == "__main__":
    from XYZ import XYZReader, XYZWriter
    from os.path import join
    from os import system
    from copy import deepcopy

    safe_print("Read in an xyz file and build from a list.")
    atom_list = []
    with XYZReader(join("Database", "XYZs", "SiO.xyz")) as reader:
        for at in reader:
            atom_list.append(at)
    frag1 = Fragment(atomlist=atom_list)
    for at in frag1:
        safe_print(at.sym, at.get_position())
    safe_print()

    safe_print("Build from an xyz file directory.")
    reader = XYZReader(join("Database", "XYZs", "Si4.xyz"))
    frag2 = Fragment(xyzfile=reader)
    for at in frag2:
        safe_print(at.sym, at.get_position())
    safe_print()

    safe_print("We can combine two fragments with +=")
    frag3 = deepcopy(frag1)
    frag3 += frag2
    for at in frag3:
        print(at.sym, at.get_position())
    safe_print("Length of frag3", len(frag3))
    safe_print()

    safe_print("Since we can iterate easily, we can also write easily.")
    with XYZWriter("test.xyz", len(frag3), "angstroem") as writer:
        for at in frag3:
            writer.write(at)
    system("cat test.xyz")
