# -*- Autoconf -*-
#
# Copyright (c) 2018 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#
AC_DEFUN([AX_PSPIO],
[dnl Test for PSPIO
  AC_REQUIRE([AX_LINALG])
  AX_PACKAGE([PSPIO],[0.2.0],[-lpspiof -lpspio],[-lgsl $LINALG_LIBS],[],
             [program main
    use pspiof_m
  
    write(*,*) PSPIO_SUCCESS
    end program],
         [use pspiof_m, only: pspiof_pspdata_alloc, pspiof_pspdata_t
    type(pspiof_pspdata_t) :: pspio
    integer :: ierr
  
    ierr = pspiof_pspdata_alloc(pspio)])
])
