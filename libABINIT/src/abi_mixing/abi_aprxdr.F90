!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_aprxdr
!! NAME
!! abi_aprxdr
!!
!! FUNCTION
!! Compute the approximative derivatives of the energy at different
!! points along the line search, thanks to a finite-difference formula.
!! This formula is the projection along the line search of the
!! Eq.(11) in PRB54, 4383 (1996).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2011 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex: if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!! choice= if==3, compute dedv_new, dedv_old, and dedv_mix,
!! if/=3, compute only dedv_new and dedv_old.
!! i_vresid and i_rhor, see the next lines.
!! f_fftgr(nfft,nspden,n_fftgr)=different functions defined on the fft grid :
!! The last residual potential is in f_fftgr(:,:,i_vresid(1)).
!! The old  residual potential is in f_fftgr(:,:,i_vresid(2)).
!! The previous old residual potential is in f_fftgr(:,:,i_vresid(3)).
!! (needed only when choice==3)
!! The old  density is in f_fftgr(:,:,i_rhor2).
!! f_atm(3,natom,n_fftgr)=different functions defined for each atom :
!! The last HF force is in f_atm(:,:,i_vresid(1)).
!! The old  HF force is in f_fftgr(:,:,i_vresid(2)).
!! The previous old HF force is in f_fftgr(:,:,i_vresid(3)).
!! (needed only when choice==3)
!! The old atomic positions are in f_atm(:,:,i_rhor2)
!! moved_atm_inside: if==1, the atoms are allowed to move.
!! mpi_comm=the mpi communicator used for the summation
!! mpi_sumarize=set it to .true. if parallelisation is done over FFT
!! natom=number of atoms in unit cell
!! nfft=(effective) number of FFT grid points (for this processor)
!! nfftot=total number of FFT grid points
!! nspden=number of spin-density components
!! rhor(nfft,nspden)=actual density
!! xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!! dedv_mix=approximate derivative from previous old residual
!! dedv_new=approximate derivative from new residual
!! dedv_old=approximate derivative from old residual (output only when choice==3)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_aprxdr(cplex,choice,dedv_mix,dedv_new,dedv_old,&
&  f_atm,f_fftgr,i_rhor2,i_vresid,moved_atm_inside,&
&  natom,nfft,nfftot,nspden,n_fftgr,rhor,ucvol,xred,fdot,user_data)

 use abi_defs_basis
 use abi_interfaces_mixing, except_this_one => abi_aprxdr

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,i_rhor2,moved_atm_inside,n_fftgr,natom,nfft
 integer,intent(in) :: nfftot,nspden
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: dedv_mix,dedv_new,dedv_old
!arrays
 integer,intent(in) :: i_vresid(3)
 integer, intent(in) :: user_data(:)
 real(dp),intent(in) :: f_atm(3,natom,n_fftgr)
 real(dp),intent(in) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
 real(dp),intent(in) :: rhor(cplex*nfft,nspden),xred(3,natom)

 interface
    function fdot(x,y,cplex,nfft,nspden,opt_denpot,user_data)
      integer, intent(in) :: cplex,nfft,nspden,opt_denpot
      double precision, intent(in) :: x(*), y(*)
      integer, intent(in) :: user_data(:)
      
      double precision :: fdot
    end function fdot
 end interface

!Local variables-------------------------------
!scalars
 integer :: iatom,idir
!arrays
 real(dp) :: dedv_temp(1)
 real(dp),allocatable :: ddens(:,:,:)

! *************************************************************************

 allocate(ddens(cplex*nfft,nspden,1))

!Compute approximative derivative of the energy
!with respect to change of potential

 ddens(:,:,1)=rhor(:,:)-f_fftgr(:,:,i_rhor2)

!call dotprod_vn(cplex,1,ddens,dedv_old,nfft,nfftot,nspden,1,vresid,ucvol)
!Dot product ddens(:,:,1) f_fftgr(:,:,i_vresid(2))
!!$ call abi_dotprodm_vn(cplex,1,ddens,dedv_temp,1,i_vresid(2),mpi_comm,mpi_summarize,1,1,1,&
!!$& nfft,n_fftgr,nspden,f_fftgr)
 dedv_temp(1) = fdot(ddens(1,1,1), f_fftgr(1,1,i_vresid(2)), &
      & cplex, nfft, nspden, 2, user_data) * ucvol / dble(nfftot)
 dedv_old = dedv_temp(1)

!Dot product ddens(:,:,1) f_fftgr(:,:,i_vresid(1))
!!$ call abi_dotprodm_vn(cplex,1,ddens,dedv_temp,1,i_vresid(1),mpi_comm,mpi_summarize,1,1,1,&
!!$& nfft,nfftot,n_fftgr,nspden,f_fftgr,ucvol)
 dedv_temp(1) = fdot(ddens(1,1,1), f_fftgr(1,1,i_vresid(1)), &
      & cplex, nfft, nspden, 2, user_data) * ucvol / dble(nfftot)
 dedv_new= dedv_temp(1)

 if(choice==3)then
!  Dot product ddens(:,:,1) f_fftgr(:,:,i_vresid(3))
!!$   call abi_dotprodm_vn(cplex,1,ddens,dedv_temp,1,i_vresid(3),mpi_comm,mpi_summarize,1,1,1,&
!!$&   nfft,nfftot,n_fftgr,nspden,f_fftgr,ucvol)
   dedv_temp(1) = fdot(ddens(1,1,1), f_fftgr(1,1,i_vresid(3)), &
        & cplex, nfft, nspden, 2, user_data) * ucvol / dble(nfftot)
   dedv_mix = dedv_temp(1)
 end if

 deallocate(ddens)

!-------------------------------------------------------

!Now, take care of eventual atomic displacements

 if(moved_atm_inside==1)then
   do idir=1,3
     do iatom=1,natom
       dedv_new=dedv_new+&
&       f_atm(idir,iatom,i_vresid(1))*(xred(idir,iatom)-f_atm(idir,iatom,i_rhor2))
       dedv_old=dedv_old+&
&       f_atm(idir,iatom,i_vresid(2))*(xred(idir,iatom)-f_atm(idir,iatom,i_rhor2))
       if(choice==3) dedv_mix=dedv_mix+&
&       f_atm(idir,iatom,i_vresid(3))*(xred(idir,iatom)-f_atm(idir,iatom,i_rhor2))
     end do
   end do
 end if

end subroutine abi_aprxdr
!!***
