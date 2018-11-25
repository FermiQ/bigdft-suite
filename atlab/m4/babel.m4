# -*- Autoconf -*-
#
# Copyright (c) 2018 BigDFT Group (D. Caliste, L.Genovese)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#
dnl we might think of inserting this option in the suitepkg.m4 file
AC_DEFUN([AX_OPENBABEL],
[dnl Test for PSPIO
ax_have_OPENBABEL_search="yes"
AC_ARG_ENABLE(openbabel, AS_HELP_STRING([--enable-openbabel], [Enable detection of openbabel compilation (enabled by default).]),
			 ax_have_OPENBABEL_search=$enableval, ax_have_OPENBABEL_search="yes")
if test x"$ax_have_OPENBABEL_search" = "xyes" ; then
AX_PACKAGE([OPENBABEL],[2.0],[-lopenbabel],[],[],
           [
#include <openbabel/obconversion.h>

int main(int argc, char **argv)
{
  std::ifstream fin("test.xyz");
  std::istream* pIn = &fin;
  OpenBabel::OBConversion conv(pIn, NULL);

  return 0;
}
],
           [OpenBabel::OBConversion conv],
           [C++], [openbabel-2.0],[#include <openbabel/obconversion.h>])
else
  ax_have_OPENBABEL="no"
fi
AM_CONDITIONAL(HAVE_OPENBABEL, test $ax_have_OPENBABEL = "yes")
])
