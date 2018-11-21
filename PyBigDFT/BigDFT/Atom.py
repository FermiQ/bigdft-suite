'''
This module defines the atom class, which is a class which contains very
general descriptions of a singel atom.
'''
from collections import MutableMapping
from futile.Utils import write as safe_print

AU_to_A = 0.52917721092
MULTIPOLE_ANALYSIS_KEYS = ['q0', 'q1', 'q2', 'sigma']
PROTECTED_KEYS = MULTIPOLE_ANALYSIS_KEYS + ["frag"] + ["r"]


def _GetSymbol(atom):
    """
    Provides the key which defines the element of the of atom.

    Arguments:
      atom (dict): a dictionary describing the atom.
    Returns:
      (str): the symbol the atom.
    """
    ks = atom.keys()
    if 'sym' in ks:
        return atom['sym']

    for k in ks:
        if k not in PROTECTED_KEYS and type(atom[k]) == type([]):
            if len(atom[k]) == 3:
                return k

    raise ValueError


class Atom(MutableMapping):
    """
    Defines a wrapper for atoms.

    An atom may have many quantities associated with it. These quantities
    are get and set in a dictionary like fashion, allow this atom to
    dynamically hold whatever data you need. But we still wrap it in a class
    so that we can have some common operations for it, as well as so we
    can maintain suitable units.

    It is this class's responsibility to extract the main properties of an
    atom (position, symbol) from the dictionary.

    Args:
      data (dict):
        A dictionary of miscellaneous values to associate with this atom.
    """

    def __init__(self, data):
        self.store = dict()
        self._extract_symbol(data)
        self.update(data)
        self._extract_position()

    def _extract_symbol(self, data):
        # Setup symbol
        self.sym = _GetSymbol(data)
        if self.sym == 'r':
            self.sym = data['sym']

    def _extract_position(self):
        # Setup position
        if 'r' in self.store:
            self._set_position(self.store['r'], self.store["units"])
        else:
            self._set_position(self.store[self.sym], self.store["units"])

    def dict(self):
        """
        Convert to a dictionary.
        """
        return self.store

    def _set_position(self, new_pos, units):
        """
        Set the position of the atom.

        Args:
          new_pos (list): a list of positions.
          units (str): the units the position is being passed in.
        """
        from numpy import array
        self._position = array([float(x) for x in new_pos])
        if units == 'angstroem' or units == 'angstroemd0':
            self._position /= AU_to_A

    def get_position(self, units="bohr"):
        """
        Returns the position of the atom in the desired units.

        Args:
          units (str): the units to return the position in. Default is bohr.

        Returns:
          An array of position values.
        """
        if units == "angstroem" or units == "angstroemd0":
            return self._position * AU_to_A
        return self._position

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value
        if "units" in self and (key == "r" or key == self.sym):
            self._extract_position()
        if key == 'sym':
            self._extract_symbol(self.store)

    def __delitem__(self, key):
        if key == self.sym:
            raise ValueError("You can't delete the symbol from an atom.")
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key

    def __eq__(self, other):
        """
        Compare two atoms. They are equal if they have the same position and
        symbol.

        other (dict, Atom): the atom (or something that can be cast to one)
          to compare with.
        """
        from numpy.linalg import norm
        if not isinstance(other, Atom):
            othercomp = Atom(other)
        else:
            othercomp = other
        return norm(othercomp._position - self._position) < 1e-10 and \
            othercomp.sym == self.sym


if __name__ == "__main__":
    """Test the atom module"""
    test_atom = Atom({'r': [1.0, 0.0, 0.0], 'sym': "He", 'units': 'bohr'})
    # Access the full data
    safe_print(dict(test_atom))
    # Access the derived data
    safe_print(test_atom.sym)
    safe_print(test_atom.get_position())
    safe_print(test_atom.get_position('angstroem'))

    # Create a new atom with different units
    safe_print()
    new_atom = Atom({
        'r': [float(x) for x in test_atom.get_position('angstroem')],
        'sym': test_atom.sym, 'units': 'angstroem'})
    # Are these atoms equal?
    safe_print(new_atom == test_atom)

    # Now other times we get an array that looks more like this
    safe_print()
    test_atom = Atom({'He': [1.0, 0.0, 0.0], 'units': 'bohr'})
    # We still preserve the old style of dictionary
    safe_print(dict(test_atom))
    # But everything else works as expected
    safe_print(test_atom.sym)
    safe_print(test_atom.get_position())
    safe_print(new_atom == test_atom)

    # The atom can be used as a dictionary for adding new properties.
    safe_print()
    test_atom["frag"] = "ANA"
    for key, value in test_atom.items():
        safe_print(key, value)
    # And if we update the dictionary position or symbol, everything else
    # reacts with suitable caution.
    test_atom["He"] = [-1.0, 0.0, 0.0]
    safe_print(dict(test_atom))
    safe_print(test_atom.get_position('angstroem'))

    # There is a protection to prevent you from deleting the symbol of an
    # atom from the dictionary.
    safe_print()
    try:
        del test_atom["frag"]
        del test_atom["He"]
    except ValueError as v:
        safe_print(v)
    safe_print(dict(test_atom))

    # But you can change the symbol if you are working with the other
    # representation.
    safe_print()
    safe_print(dict(new_atom))
    new_atom["sym"] = "Na"
    safe_print(new_atom.sym)
    safe_print(dict(new_atom))

    # One final check of the atom comparison
    safe_print()
    new_atom["units"] = "bohr"
    new_atom["r"] = [-1.0, 0.0, 0.0]
    safe_print(new_atom.sym, new_atom.get_position())
    safe_print(test_atom.sym, test_atom.get_position())
    safe_print(new_atom == test_atom)
