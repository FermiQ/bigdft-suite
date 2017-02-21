//> @file
// Routines wraping OpenBabel read routines.
// @author
//    Copyright (C) 2017 BigDFT group
//    This file is distributed under the terms of the
//    GNU General Public License, see ~/COPYING file
//    or http://www.gnu.org/copyleft/gpl.txt .
//    For the list of contributors, see ~/AUTHORS

#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/math/vector3.h>

#include <config.h>
#include <stdlib.h>

extern "C" void FC_FUNC_(astruct_set_cell, ASTRUCT_SET_CELL)(void *, const double [3]);
extern "C" void FC_FUNC_(astruct_add_atom, ASTRUCT_ADD_ATOM)(void *, const double [3], const char*, int *, unsigned int);

extern "C" void FC_FUNC_(openbabel_load, OPENBABEL_LOAD)(void *dict_posinp,
                                                       const char *filename, unsigned int flen)
{
  char *fname = (char*)malloc(sizeof(char) * (flen + 1));
  memcpy(fname, filename, sizeof(char) * flen);
  fname[flen] = '\0';
  std::ifstream fin(fname);
  OpenBabel::OBConversion conv(&fin, NULL);

  OpenBabel::OBFormat *pFormat;
  pFormat = conv.FormatFromExt(fname);

  free(fname);

  if (!pFormat || (pFormat->Flags() & NOTREADABLE))
    return;

  conv.SetInFormat(pFormat);

  OpenBabel::OBMol mol;

  if (!conv.Read(&mol))
    return;

  /* Store if the file is periodic or not. */
  double vect[3], cell[3];
  OpenBabel::OBUnitCell *uc = (OpenBabel::OBUnitCell*)mol.GetData(OpenBabel::OBGenericDataType::UnitCell);
  if (uc)
    {
      double rprimdFull[9];
      uc->GetCellMatrix().GetArray(rprimdFull);
      if (rprimdFull[1] > 1e-12 || rprimdFull[1] < -1e-12 ||
          rprimdFull[2] > 1e-12 || rprimdFull[2] < -1e-12 ||
          rprimdFull[3] > 1e-12 || rprimdFull[3] < -1e-12 ||
          rprimdFull[5] > 1e-12 || rprimdFull[5] < -1e-12 ||
          rprimdFull[6] > 1e-12 || rprimdFull[6] < -1e-12 ||
          rprimdFull[7] > 1e-12 || rprimdFull[7] < -1e-12)
        return;
      cell[0] = rprimdFull[0];
      cell[1] = rprimdFull[4];
      cell[2] = rprimdFull[8];
      uc->GetOffset().Get(vect);
      uc->FillUnitCell(&mol);
      FC_FUNC_(astruct_set_cell, ASTRUCT_SET_CELL)(dict_posinp, cell);
    }
  else
    {
      vect[0] = 0.;
      vect[1] = 0.;
      vect[2] = 0.;
    }

  /* Stores coordinates. */
  FOR_ATOMS_OF_MOL(a, mol)
    {
      double xyz[3];
      xyz[0] = a->x() + vect[0];
      xyz[1] = a->y() + vect[1];
      xyz[2] = a->z() + vect[2];
      if (uc && (xyz[0] / cell[0] > 1 - 1e-6 ||
                 xyz[1] / cell[1] > 1 - 1e-6 ||
                 xyz[2] / cell[2] > 1 - 1e-6))
        continue;
      const char *symbol = OpenBabel::etab.GetSymbol(a->GetAtomicNum());
      int len = strlen(symbol);
      FC_FUNC_(astruct_add_atom, ASTRUCT_ADD_ATOM)(dict_posinp, xyz,
                                                   symbol, &len, len);
    }
}
