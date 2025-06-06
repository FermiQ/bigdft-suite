# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Luigi Genovese, Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
AC_DEFUN([AX_CHESS],
[dnl Test for PSolver
dnl AC_REQUIRE([AX_FLIB])
dnl AC_REQUIRE([AX_LINALG])
AC_REQUIRE([AX_ATLAB])
AC_REQUIRE([AX_MPI])
AX_PACKAGE([CHESS],[0.1.3],[-lCheSS-1],[$LIB_ATLAB_LIBS],[$LIB_ATLAB_CFLAGS],
[program main
    use sparsematrix_base
    use sparsematrix_highlevel
    implicit none
    type(sparse_matrix) :: smat
    type(matrices) :: mat
    call matrices_init(smat, mat)
  end program],
  [     use sparsematrix_base
  use sparsematrix_highlevel
  implicit none
  type(sparse_matrix) :: smat
  type(matrices) :: mat
  call matrices_init(smat, mat)
])
])
