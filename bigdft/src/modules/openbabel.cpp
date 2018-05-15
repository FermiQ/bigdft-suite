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
#include <openbabel/oberror.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/math/vector3.h>

#include <config.h>
#include <stdlib.h>
#include <iostream>
extern "C" {
#include <futile.h>
}
extern "C" void FC_FUNC_(astruct_get_types_dict, ASTRUCT_GET_TYPES_DICT)(f90_dictionary_pointer* dict, f90_dictionary_pointer* types);

extern "C" void FC_FUNC_(openbabel_formats,OPENBABEL_FORMATS)(f90_dictionary_pointer* input_formats, f90_dictionary_pointer* output_formats)
{
  OpenBabel::OBConversion conv(NULL, NULL);
  std::vector<std::string> informats, outformats;
  std::vector<std::string>::iterator it;

  informats = conv.GetSupportedInputFormat();
  outformats = conv.GetSupportedOutputFormat();

  dict_init(input_formats);
  int i=0;
  for(it = informats.begin(); it != informats.end(); it++,i++)    {
    dict_add_string(input_formats, it->c_str());
    }

  dict_init(output_formats);
  i=0;
  for(it = outformats.begin(); it != outformats.end(); it++,i++)    {
    dict_add_string(output_formats, it->c_str());
  }

}

extern "C" void FC_FUNC_(openbabel_load, OPENBABEL_LOAD)(f90_dictionary_pointer *dict_posinp,
                                                         const char *filename) //, unsigned int *flen)
{
  OpenBabel::OBMol mol;
  {
    /*char *fname = (char*)malloc(sizeof(char) * (*flen + 1));
    memcpy(fname, filename, sizeof(char) * *flen);
    fname[*flen] = '\0';
    std::ifstream fin(fname); */
    std::ifstream fin(filename);
    OpenBabel::OBConversion conv(&fin, NULL);

    OpenBabel::OBFormat *pFormat;
    pFormat = conv.FormatFromExt(filename);
    //pFormat = conv.FormatFromExt(fname);

    //free(fname);

    if (!pFormat || (pFormat->Flags() & NOTREADABLE))
      {
        err_throw_by_name("Unknown format for OpenBabel.",
                          "BIGDFT_INPUT_VARIABLES_ERROR");
        return;
      }

    std::ofstream ferr("OpenBabel.err");
    OpenBabel::obErrorLog.SetOutputStream(&ferr);
    conv.SetInFormat(pFormat);
    conv.Read(&mol);
  }
  if (mol.Empty()) {
    std::vector<std::string> messages(OpenBabel::obErrorLog.GetMessagesOfLevel(OpenBabel::obError));
    if (messages.empty())
      messages = OpenBabel::obErrorLog.GetMessagesOfLevel(OpenBabel::obWarning);
    if (!messages.empty())
      err_throw_by_name(messages[0].c_str(),
                        "BIGDFT_INPUT_VARIABLES_ERROR");
    else
      err_throw_by_name("Error while reading OpenBabel format.",
                        "BIGDFT_INPUT_VARIABLES_ERROR");
    return;
  }

  /* Store if the file is periodic or not. */
  double vect[3], alphabetagamma[3], cell[3];
  OpenBabel::OBUnitCell *uc(static_cast<OpenBabel::OBUnitCell*>(mol.GetData(OpenBabel::OBGenericDataType::UnitCell)));
  if (uc)
    {
      double rprimdFull[9];
      uc->GetCellMatrix().GetArray(rprimdFull);
      alphabetagamma[0] = uc->GetAlpha();
      alphabetagamma[1] = uc->GetBeta();
      alphabetagamma[2] = uc->GetGamma();
      cell[0] = uc->GetA();
      cell[1] = uc->GetB();
      cell[2] = uc->GetC();
      uc->GetOffset().Get(vect);
      uc->FillUnitCell(&mol);
      dict_set_double_array(dict_posinp, "cell", cell, 3);
      dict_set_double(dict_posinp, "alpha", alphabetagamma[0]);
      dict_set_double(dict_posinp, "beta", alphabetagamma[1]);
      dict_set_double(dict_posinp, "gamma", alphabetagamma[2]);
      dict_set_double_matrix(dict_posinp, "abc", rprimdFull, 3, 3);
    }
  else
    {
      vect[0] = 0.;
      vect[1] = 0.;
      vect[2] = 0.;
    } 
  /* the units for openbabel formats is supposed to be angstroem unless otherwise specified */
  dict_set_string(dict_posinp, "units", "angstroem");
  /* retrieve positions */
  f90_dictionary_pointer dict_positions;
  dict_init(&dict_positions);

  /* Stores coordinates. */
  FOR_ATOMS_OF_MOL(a, mol)
    {
      double xyz[3],xred[3];
      xyz[0] = a->x() + vect[0];
      xyz[1] = a->y() + vect[1];
      xyz[2] = a->z() + vect[2];
      xred[0] = xyz[0] / cell[0];
      xred[1] = xyz[1] / cell[1];
      xred[2] = xyz[2] / cell[2];
      //very dangerous condition, should be removed
      /*if (uc && (xred[0] > 1 - 1e-6 && xred [0] < 1 ||
                 xred[1] > 1 - 1e-6 && xred [1] < 1 ||
		 xred[2] > 1 - 1e-6 && xred [2] < 1 ))
		 continue;*/
      f90_dictionary_pointer atom;
      dict_init(&atom);
      dict_set_double_array(&atom, OpenBabel::etab.GetSymbol(a->GetAtomicNum()), xyz, 3);
      OpenBabel::OBResidue *residue = a->GetResidue();
      if (residue)
        {
          f90_dictionary_pointer frag;
          dict_init(&frag);
          dict_add_string(&frag, residue->GetName().c_str());
          dict_add_int(&frag, residue->GetNum());
          dict_set_dict(&atom, "frag", &frag);
        }
      dict_add_dict(&dict_positions, &atom);
    }

  dict_set_dict(dict_posinp, "positions", &dict_positions);
}

extern "C" void FC_FUNC_(openbabel_dump, OPENBABEL_DUMP)(f90_dictionary_pointer *dict_posinp,
                                                         f90_dictionary_pointer *dict_types,
							 const char *filename) //, unsigned int *flen)
{

  /* Ensure output file */
  /*char *fname = (char*)malloc(sizeof(char) * (*flen + 1));
  memcpy(fname, filename, sizeof(char) * *flen);
  fname[*flen] = '\0'; 
  std::ofstream fout(fname); */

  std::ofstream fout(filename);

  OpenBabel::OBConversion conv(NULL, &fout);

  OpenBabel::OBFormat *pFormat;
  pFormat = conv.FormatFromExt(filename);
  //pFormat = conv.FormatFromExt(fname);

  //free(fname);

  if (!pFormat || (pFormat->Flags() & NOTWRITABLE))
    {
      err_throw_by_name("Unknown format for OpenBabel.",
                        "BIGDFT_RUNTIME_ERROR");
      return;
    }

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

  /*iterate on atoms */
  f90_dictionary_pointer dict_positions;
  if (dict_get_dict(dict_posinp, "positions", &dict_positions))
    {
      f90_dictionary_iterator it_atom;
      dict_iter_new(&it_atom, &dict_positions);
      while (iterate(&it_atom))
	{
          #define UNSET -123456789.0
	  f90_dictionary_pointer coord;
	  double xyz[3] = {UNSET, UNSET, UNSET};
	  f90_dictionary_iterator it_type;
	  dict_iter_new(&it_type, dict_types);
	  while (iterate(&it_type))
	    {
	      if (dict_get_double_array(&it_atom.dict, it_type.key, xyz, 3))
		break;
	    }
          if (xyz[0] != UNSET && xyz[1] != UNSET && xyz[2] != UNSET)
            {
              OpenBabel::OBAtom *atom;
              atom = mol.NewAtom();
              atom->SetAtomicNum(OpenBabel::etab.GetAtomicNum(it_type.value));
              atom->SetVector(xyz[0], xyz[1], xyz[2]);
            }
	}
    }
  if (!conv.Write(&mol))
    {
      err_throw_by_name("Error while writing OpenBabel format.",
                        "BIGDFT_RUNTIME_ERROR");
      return;
    }
  
  fout.close();
 
}
