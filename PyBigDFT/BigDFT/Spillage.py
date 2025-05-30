"""
A module which contains the routines needed for computing the spillage values.
"""

#: Conversion between Atomic Units and Bohr
AU_to_A = 0.52917721092


class MatrixMetadata:
    def __init__(self, filename):
        """
        This class contains the information stored in the sparse matrix
        metadata file.

        Args:
          filename (str): the name of the metadata file to initialize from.

        Attributes:
          atoms (list): a list of atoms with their positions and associated
            basis functions.
          matdim (int): the dimension of the matrix.
        """
        from Atom import Atom
        self.atoms = []
        symbol_lookup = []
        with open(filename, "r") as ifile:
            # Read the first line
            matinfo = next(ifile).split()
            self.matdim, natoms, ntypes = [int(x) for x in matinfo[:3]]
            # Units
            units = next(ifile).split()[0]
            # skip geocode
            line = next(ifile)
            # skip shift?
            line = next(ifile)
            # get the symbol lookup information
            for i in range(0, ntypes):
                line = next(ifile).split()
                symbol_lookup.append(line[2])
            # Get the atom positions
            for i in range(0, natoms):
                line = next(ifile).split()
                adict = {}
                adict["sym"] = symbol_lookup[int(line[0]) - 1]
                adict["r"] = [float(x) for x in line[1:4]]
                # THIS IS BECAUSE THE METADATA FILE CURRENTLY PRINTS THE
                # WRONG UNITS!
                adict["units"] = "bohr"
                adict["indices"] = []
                self.atoms.append(Atom(adict))
            # Get the indices
            for i in range(0, self.matdim):
                line = next(ifile).split()
                self.atoms[int(line[0]) - 1]["indices"].append(i)

    def get_frag_indices(self, system):
        """
        Retrive the indices of the matrix associated with each fragment
        of a system.

        Args:
          system (Fragments.System): the set of fragments to get the indices
            of.
        """

        from numpy import array
        from numpy.linalg import norm
        # Distance threshold to consider two atoms equal
        thresh = 1e-5

        # Now we associate this information with fragments.
        frag_indices = {}
        for fragid in system:
            frag_indices[fragid] = []
            temp_ind = []
            # Find the atom that matches
            for aj in system[fragid].atoms:
                for ai in self.atoms:
                    if aj == ai:
                        found = True
                        temp_ind += ai["indices"]
                        break
            frag_indices[fragid] = sorted(temp_ind)

        return frag_indices


def compute_spillage_values(sinvxh, sinvxh2, frag_indices, target):
    """
    Computes the actual spillage values.

    Args:
      sinvxh (scipy.sparse.csc): S^-1 * H
      sinvxh2 (scipy.sparse.csc): (S^-1 * H)^2
      frag_indices (list): list of indices associated with each fragment.
      target (int): which fragment is the target fragment.

    Returns:
      list: interaction spillage between the target and each fragment.
    """
    from numpy import trace
    indices_f = frag_indices[target]
    spillage_values = {}

    # Compute the denominator tr(HRfHRf)
    denom = sinvxh[:, indices_f]
    denom = denom[indices_f, :]
    denom = denom.dot(denom)
    denom_t = trace(denom.todense())

    # Compute the left side tr(HS-1HRf)
    H2T = sinvxh2[:, indices_f]
    H2T = H2T[indices_f, :]
    left_t = trace(H2T.todense())

    # Compute the right side values tr(HRgHRf)
    for id_g in frag_indices:
        indices_g = frag_indices[id_g]
        TFH = sinvxh[indices_f, :]
        TFHTG = TFH[:, indices_g]

        TGH = sinvxh[indices_g, :]
        TGHTF = TGH[:, indices_f]

        right_mat = (TFHTG.dot(TGHTF))
        right_t = trace(right_mat.todense())
        spillage_values[id_g] = right_t / denom_t

    return spillage_values


def serial_compute_spillbase(sfile, hfile):
    """
    This routine computes the matrix (S^-1 * H)^2 using python.

    You will need this for computing the spillage values.

    Args:
      sfile (str): the file name of the overlap matrix. Must be in ccs format
        if ``bigpoly=False``, bigdft format otherwise.
      hfile (str): the file name of the hamiltonian. Must be in ccs format
        if ``bigpoly=False``, bigdft format otherwise.

     Returns:
       (scipy.sparse.csc, scipy.sparse.csc): (S^-1 * H), (S^-1 * H)^2
    """
    from scipy.sparse.linalg import inv
    from scipy.sparse import csc_matrix
    from scipy.io import mmread

    # Read from file
    smat = _read_ccs(sfile)
    hmat = _read_ccs(hfile)

    # Compute the matrix (S^-1H)^2
    sinv = inv(smat)
    sinvxh = sinv.dot(hmat)
    sinvxh2 = sinvxh.dot(sinvxh)

    return sinvxh, sinvxh2


def _read_ccs(fname):
    """
    Read a CCS matrix from file into a scipy csc matrix.

    This routine uses the CCS matrix files that CheSS generates.
    Args:
      fname (str): name of the file.

    Result:
      scipy.sparse.csc_matrix: the matrix in the file.
    """
    from scipy.sparse import csc_matrix

    data = []
    indices = []
    indptr = []

    with open(fname, "r") as ifile:
        # Read the meta data
        line = next(ifile)
        split = line.split()
        matdim = int(split[0])
        nel = int(split[2])
        split = next(ifile).split()
        indptr = [int(x) - 1 for x in split]

        # Read the indices
        added = 0
        while(added < nel):
            split = next(ifile).split()
            values = [int(x) - 1 for x in split]
            indices.extend(values)
            added += len(values)

        # Read the data
        for line in ifile:
            data.append(float(line))

    matrix = csc_matrix((data, indices, indptr), shape=(matdim, matdim))
    return matrix
