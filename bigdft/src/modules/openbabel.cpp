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
extern "C" {
#include <futile.h>
}
//extern "C" void FC_FUNC_(astruct_set_cell, ASTRUCT_SET_CELL)(void *, const double [3]);
//extern "C" void FC_FUNC_(astruct_add_atom, ASTRUCT_ADD_ATOM)(void *, const double [3], const char*, int *, unsigned int);
extern "C" void FC_FUNC_(astruct_get_types_dict, ASTRUCT_GET_TYPES_DICT)(f90_dictionary_pointer* dict, f90_dictionary_pointer* types);

extern "C" void FC_FUNC_(openbabel_load, OPENBABEL_LOAD)(f90_dictionary_pointer *dict_posinp,
                                                       const char *filename, unsigned int *flen)
{
  char *fname = (char*)malloc(sizeof(char) * (*flen + 1));
  memcpy(fname, filename, sizeof(char) * *flen);
  fname[*flen] = '\0';
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
      //FC_FUNC_(astruct_set_cell, ASTRUCT_SET_CELL)(dict_posinp, cell);
      dict_set_double_array(dict_posinp, "cell", cell, 3);
    }
  else
    {
      vect[0] = 0.;
      vect[1] = 0.;
      vect[2] = 0.;
    }

  /* retrieve positions */
  f90_dictionary_pointer dict_positions;
  dict_init(&dict_positions);

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
      //const char *symbol = OpenBabel::etab.GetSymbol(a->GetAtomicNum());
      //int len = strlen(symbol);
      //FC_FUNC_(astruct_add_atom, ASTRUCT_ADD_ATOM)(dict_posinp, xyz,
      //                                             symbol, &len, len);
      f90_dictionary_pointer atom;
      dict_init(&atom);
      dict_set_double_array(&atom, OpenBabel::etab.GetSymbol(a->GetAtomicNum()), xyz, 3);
      dict_add_dict(&dict_positions, &atom);
    }

  dict_set_dict(dict_posinp, "positions", &dict_positions);
}

extern "C" void FC_FUNC_(openbabel_dump, OPENBABEL_DUMP)(f90_dictionary_pointer *dict_posinp,
							 const char *filename, unsigned int *flen)
{

  /* Ensure output file */
  char *fname = (char*)malloc(sizeof(char) * (*flen + 1));
  memcpy(fname, filename, sizeof(char) * *flen);
  fname[*flen] = '\0';
  std::ofstream fout(fname);

  OpenBabel::OBConversion conv(NULL, &fout);

  OpenBabel::OBFormat *pFormat;
  pFormat = conv.FormatFromExt(fname);

  free(fname);

  if (!pFormat || (pFormat->Flags() & NOTWRITABLE))
    return;

  conv.SetOutFormat(pFormat);


  /* Create a new OpenBabel object. */
  OpenBabel::OBMol mol;
  OpenBabel::OBUnitCell *cell;

  /* if cell key exists */
  double acell[3];
  if (dict_get_double_array(dict_posinp, "cell", acell, 3))
    {
      OpenBabel::vector3 a(0.,0.,0.), b(0.,0.,0.), c(0.,0.,0.);
      double rprimdFull[3][3];
      memset(rprimdFull, 0, sizeof(double)*9);
      rprimdFull[0][0]=acell[0];
      rprimdFull[1][1]=acell[1];
      rprimdFull[2][2]=acell[2];
      a.Set(rprimdFull[0]);
      b.Set(rprimdFull[1]);
      c.Set(rprimdFull[2]);
      cell=new OpenBabel::OBUnitCell;
      cell->SetData(a, b, c);
      mol.SetData(cell);    
    }

  /* get types dict. */
  f90_dictionary_pointer types;
  FC_FUNC_(astruct_get_types_dict,ASTRUCT_GET_TYPES_DICT)(dict_posinp,&types);

  /*iterate on atoms */
  f90_dictionary_pointer dict_positions;
  if (dict_get_dict(dict_posinp, "positions", &dict_positions))
    {
      f90_dictionary_iterator it_atom;
      dict_iter_new(&it_atom, &dict_positions);
      while (iterate(&it_atom))
	{
	  f90_dictionary_pointer coord;
	  char *symbol;
	  double xyz[3]={-123456789.0,-123456789.0,-123456789.0};
	  f90_dictionary_iterator it_type;
	  dict_iter_new(&it_type, &types);
	  while (iterate(&it_type))
	    {
	      symbol = it_type.key;	      
	      if (dict_get_double_array(&it_atom.dict, symbol, xyz, 3)) //dict_get_dict(&it_atom.dict, symbol, &coord))
		break;
	    }
	  OpenBabel::OBAtom *atom;
	  atom = mol.NewAtom();
	  atom->SetAtomicNum(OpenBabel::etab.GetAtomicNum(symbol));
	  atom->SetVector(xyz[0], xyz[1], xyz[2]);
	}
    }
  dict_free(&types);
  if (!conv.Write(&mol))
    return;
  
  fout.close();
 
}
