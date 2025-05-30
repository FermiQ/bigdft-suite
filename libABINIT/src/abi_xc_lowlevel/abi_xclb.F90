!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_xclb
!! NAME
!! abi_xclb
!!
!! FUNCTION
!! Computes the GGA like part (vx_lb) of the Leeuwen-Baerends
!! exchange-correlation potential (vxc_lb) and adds it to the
!! lda exchange-correlation potential (vxc_lda) which
!! must be provided as input,
!!            vxci  <--  vxc_lb =: vxc_lda + vx_lb
!!
!! [R van Leeuwen and EJ Baerends, Phys Rev A 49, 2421 (1994)]
!!
!! With respect to spin, the van Leeuwen-Baerends
!! potential is "exchange-like" : separate contributions from
!! spin up and spin down.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (MF,LG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npts= number of points to be computed
!!  nspden=1 for unpolarized, 2 for spin-polarized
!!  grho2_updn(npts,2*nspden-1)=square of the gradient of the spin-up,
!!     and, if nspden==2, spin-down, and total density (Hartree/Bohr**2)
!!  rho_updn(npts,nspden)=spin-up and spin-down density (Hartree/bohr**3)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  vxci(npts,nspden)=input xc potential to which Leeuwen-Baerends correction
!!   is added at output.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_xclb(grho2_updn,npts,nspden,rho_updn,vxci)

 use abi_defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts,nspden
!arrays
 real(dp),intent(in) :: grho2_updn(npts,2*nspden-1),rho_updn(npts,nspden)
 real(dp),intent(inout) :: vxci(npts,nspden)

!Local variables-------------------------------
!scalars
 integer :: ipts,ispden
 real(dp),parameter :: beta=0.05_dp
 real(dp) :: density,density_gradient,density_t13,s_g_sq,scaled_gradient
 real(dp) :: scaling_factor,vx_lb

! *************************************************************************

!scale the spin densities for evaluating spin up or down exchange
 scaling_factor=one
 if(nspden == 2) scaling_factor=two

 do ispden=1,nspden

   do ipts=1,npts

     density= scaling_factor * rho_updn(ipts,ispden)
     density_gradient= scaling_factor * sqrt(grho2_updn(ipts,ispden))

     density_t13= density**third
     scaled_gradient= density_gradient/max(density*density_t13,1.e-12_dp)

     s_g_sq= scaled_gradient*scaled_gradient

     vx_lb= -beta*density_t13 * s_g_sq/ &
&     (one+3.d0*beta* scaled_gradient*log(scaled_gradient+sqrt(one+s_g_sq*s_g_sq)))

     vxci(ipts,ispden)=vxci(ipts,ispden)+vx_lb
   end do

 end do

end subroutine abi_xclb
!!***
