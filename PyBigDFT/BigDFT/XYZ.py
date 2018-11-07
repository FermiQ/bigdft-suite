"""
This module defines some tools used for writing XYZ files.

"""

class XYZfile():
    """
    .. |filename_docs| replace::
         The file which will be created. If None, the file will be
         eventually dumped in :class:~`sys.stdout`.

    A class associated to a xyz input file as processed by BigDFT

    :param filename: |filename_docs|
    :type filename: string
    :param units: The units of measure of the positions. Allowed avlues are
      'atomic' or 'angstroem'
    :type units: string
    """

    def __init__(self, filename=None, units='atomic'):
        self.filename = filename
        self.lines = []
        self.units = units
        self.fac = 1.0
        if units == 'angstroem':
            self.fac = AU_to_A

    def append(self, array, basename='', names=None, attributes=None):
        """
        Add lines to the file position list

        :param array: list of the atomic positions
        :type array: list of  float triples
        :param basename: base for the name of the atoms
        :type basename: string
        :param names: list of atom names. Will be appended to `basename` if the latter is present
        :type names: list of strings
        :param attributes: list of further attributes to be associated to each of the atoms.
            Will be serialized close to each of the atomic positions
        :type attributes: list of dictionaries
        """
        nm = basename
        for i, r in enumerate(array):
            if names is not None:
                nm = basename + names[i]
            line = str(nm)
            for t in r:
                line += ' ' + str(self.fac * t)
            if attributes is not None:
                line += ' ' + str(attributes[i])
            self.lines.append(line + '\n')

    def dump(self, position='w'):
        """
        Dump the file on the file system if filename has been provided,
        otherwise dump on sys.stdout.

        :param position: filename position statement. Only menaingful for a file dumping.
        :type position: char
        """
        import sys
        f = sys.stdout
        if self.filename is not None:
            f = open(self.filename, position)
        f.write(str(len(self.lines)) + ' ' + str(self.units) + '\n')
        f.write('# xyz dump \n')
        # then the positions
        for l in self.lines:
            f.write(l)
        if self.filename is not None:
            f.close()


def open_xyz(filename, nat, unit, comment, position='a'):
    import sys
    f = sys.stdout
    if filename is not None:
        f = open(filename, position)
    if (position != 'a'):
        f.write(str(nat) + ' ' + str(unit) + '\n')
        f.write(comment + '\n')
    return f


def close_xyz(f, filename):
    if filename is not None:
        f.close()


def dump_xyz_positions(f, array, basename='', names=None):
    nm = basename
    for i, r in enumerate(array):
        if names is not None:
            nm = basename + names[i]
        f.write(str(nm) + ' ' + str(r[0]) + ' ' +
                str(r[1]) + ' ' + str(r[2]) + '\n')


def xyz_bc_spec(cell):
    """
    Defines the specification for expressing the Boundary Conditions starting from a cell vector.

    :param cell: array of the (orthorhombic) cell. Should be 0.0 on directions with free BC.
       If None is given, the BC are assumed to be Free.
    :type cell: triple of floats or None
    :returns: comment Line of the xyz file specifying the bc
    :rtype: string
    """
    if cell is None:
        return ""
    elif cell[1] == 0.0 and cell[2] != 0.0:
        return "surface " + str(cell[0]) + " 0.0 " + str(cell[2]) + " "
    elif cell[1] == 0.0 and cell[2] == 0.0:
        return "wire 0.0 0.0 " + cell[2] + " "
    else:
        return "periodic " + str(cell[0]) + " " + str(cell[1]) + " " + str(cell[2]) + " "


def dump_xyz(array, basename='', units='atomic', names=None, filename=None, position='a', comment=None, cell=None):
    """
    Create a BigDFT xyz filename. Duplicates the meaning of the :class:`XYZfile` class.

    :param filename: |filename_docs|
    :type filename: string

    .. todo::
       Remove the duplication of operations in favour or the class.
       Move the related operation into a lower-level module.
    """
    cmt = xyz_bc_spec(cell)
    cmt += comment if comment is not None else '# xyz dump with basename "' + basename + '"'
    f = open_xyz(filename, len(array), units, cmt, position)
    dump_xyz_positions(f, array, basename=basename, names=names)
    close_xyz(f, filename)
