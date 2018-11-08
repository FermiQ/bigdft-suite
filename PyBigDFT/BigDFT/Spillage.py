"""
A module which contains the routines needed for computing the spillage values.
"""

#: Conversion between Atomic Units and Bohr
AU_to_A = 0.52917721092


def compute_spillage(sinvxh, sinvxh2, frag_indices, target):
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
    spillage_values = []

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
    for indices_g in frag_indices:
        TFH = sinvxh[indices_f, :]
        TFHTG = TFH[:, indices_g]

        TGH = sinvxh[indices_g, :]
        TGHTF = TGH[:, indices_f]

        right_mat = (TFHTG.dot(TGHTF))
        right_t = trace(right_mat.todense())
        spillage_values.append(right_t / denom_t)

    return spillage_values


def process_metadata(mfile, system):
    """
    This routine processes the sparse matrix meta data file.

    Args:
      mfile (str): the filename of the metadata file.
      system (list): the system associated with the file.
    Returns:
      list: a list that has for each fragment the indices of the matrix
        associated with it.
    """
    from numpy import array
    from numpy.linalg import norm
    # Distance threshold to consider two atoms equal
    thresh = 1e-5

    # Get information about to each atom
    atoms = _getatoms_from_metadata(mfile)

    # Now we associate this information with fragments.
    frag_indices = []
    for fragment in system.fragments:
        temp_ind = []
        # Find the atom that matches
        for aj in fragment.atoms:
            for ai in atoms:
                if norm(array(aj["r"]) - array(ai["r"])) < thresh:
                    found = True
                    temp_ind += ai["indices"]
                    break
        frag_indices.append(sorted(temp_ind))

    return frag_indices


def _getatoms_from_metadata(mfile):
    """
    This routine gets per atom information from a metadata file.

    Args:
      mfile (str): the filename of the metadata file.
    Returns:
      list: a list that has for each atom the information provided in the
        metadata file in dictionary form.
    """
    atoms = []
    symbol_lookup = []
    with open(mfile, "r") as ifile:
        # Read the first line
        matinfo = next(ifile).split()
        dim, natoms, ntypes = [int(x) for x in matinfo[:3]]
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
            positions = [float(x) for x in line[1:4]]
            if units == "angstroem":
                positions = [x * AU_to_A for x in positions]
            adict["r"] = positions
            adict["indices"] = []
            atoms.append(adict)
        # Get the indices
        for i in range(0, dim):
            line = next(ifile).split()
            atoms[int(line[0]) - 1]["indices"].append(i)

        return atoms


def compute_spillbase(sfile, hfile):
    """
    This routine computes the matrix (S^-1 * H)^2.

    You will need this for computing the spillage values.

    Todo:
      Compute this matrix with a Fortran routine.

    Args:
      sfile (str): the file name of a CCS representation of the overlap matrix.
      hfile (str): the file name of a CCS representation of the hamiltonian.

     Returns:
       scipy.sparse.csc: (S^-1 * H)^2
    """
    from scipy.sparse.linalg import inv

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
