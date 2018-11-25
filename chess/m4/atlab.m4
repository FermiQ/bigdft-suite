# -*- Autoconf -*-
#
# Copyright (c) 2018 BigDFT Group (Luigi Genovese, Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
AC_DEFUN([AX_ATLAB],
[dnl Test for PSolver
AC_REQUIRE([AX_FLIB])
AC_REQUIRE([AX_LINALG])
AX_PACKAGE([ATLAB],[1.0],[-latlab-1],[$LIB_FUTILE_LIBS $LINALG_LIBS],[$LIB_FUTILE_CFLAGS],
[program main
  use f_harmonics
  type(f_multipoles) :: mp
  call f_multipoles_create(mp,2)
  end program],
  [call ISF_family()])
])
