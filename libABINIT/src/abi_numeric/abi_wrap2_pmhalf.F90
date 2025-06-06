!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_wrap2_pmhalf
!! NAME
!! abi_wrap2_pmhalf
!!
!! FUNCTION
!! Transforms a real number (num) in its corresponding reduced number
!! (red) in the interval ]-1/2,1/2] where -1/2 is not included (tol12)
!! num=red+shift
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! num=real number
!!
!! OUTPUT
!! red=reduced number of num in the interval ]-1/2,1/2] where -1/2 is not included
!! shift=num-red
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_wrap2_pmhalf(num,red,shift)

 use abi_defs_basis

 implicit none

!Arguments -------------------------------
!scalars
 real(dp),intent(in) :: num
 real(dp),intent(out) :: red,shift

!Local variables-------------------------------

! *********************************************************************

 if (num>zero) then
   red=mod((num+half-tol12),one)-half+tol12
 else
   red=-mod(-(num-half-tol12),one)+half+tol12
 end if
 if(abs(red)<tol12)red=0.0d0
 shift=num-red

end subroutine abi_wrap2_pmhalf
!!***
