!> @file  
!! Test the usage of the localization regions
!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program locreg_test
   use module_base
   use locregs
   use locreg_operations
   implicit none
   logical :: correct
   type(locreg_descriptors) :: lr,Glr
   integer, dimension(2,3) :: nbox
   real(f_double), dimension(3) :: hgrids

   call f_lib_initialize()

   call lr_box(lr,Glr,hgrids,nbox,correct)

   call f_lib_finalize()

end program locreg_test
