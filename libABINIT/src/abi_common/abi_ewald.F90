!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_ewald
!! NAME
!! abi_ewald
!!
!! FUNCTION
!! Compute Ewald energy and derivatives with respect to dimensionless
!!  reduced atom coordinates xred.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (bohr^-2)
!! natom=number of atoms in unit cell
!! ntypat=numbe of type of atoms
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume (bohr^3)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! eew=final ewald energy in hartrees
!! grewtn(3,natom)=grads of eew wrt xred(3,natom), hartrees.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_ewald(iproc,nproc,mpi_comm,eew,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)
 use abi_defs_basis
 use abi_interfaces_lowlevel
 use abi_interfaces_numeric

 implicit none
 !SM there are probably better ways than this...
 include 'mpif.h'

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iproc,nproc,mpi_comm,natom,ntypat
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: eew
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: gmet(3,3),rmet(3,3),xred(3,natom),zion(ntypat)
 real(dp),intent(out) :: grewtn(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ia,ib,ig1,ig2,ig3,ir1,ir2,ir3,newg,newr,ng,nr
 real(dp) :: arg,c1i,ch,chsq,derfc_arg,direct,drdta1,drdta2,drdta3,eta,fac
 real(dp) :: fraca1,fraca2,fraca3,fracb1,fracb2,fracb3,gsq,gsum,phi,phr,r1
 real(dp) :: r1a1d,r2,r2a2d,r3,r3a3d,recip,reta,rmagn,rsq,sumg,summi,summr,sumr
 real(dp) :: t1,term
 character(len=500) :: message
 real(dp) :: tt
 integer :: natp, isat, ii, iia, ierr
 real(dp),dimension(2) :: sumarr
 real(dp),dimension(:,:,:),allocatable :: grewtn_tmp

! *************************************************************************

!SM: MPI parallelization over the atoms
 tt = real(natom,kind=dp)/nproc
 natp = floor(tt) !number of atoms per proc
 isat = iproc*natp !offset for each proc
 ii = natom-nproc*natp !remaining atoms
 if (iproc<ii) then
     natp = natp+1 !one more atom for this proc
     isat = isat+iproc !offset increases by the number of additional atoms up to iproc
 else
     isat = isat+ii !offset increases by the number of additional atoms
 end if
 ! Check
 ii = natp
 if (nproc>1) then
     !call mpiallred(ii, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
     iia=0
     call mpi_allreduce(ii, iia, 1, mpi_integer, mpi_sum, mpi_comm, ierr)
     ii = iia
 end if
 !if (ii/=natom) call f_err_throw('ii/=natom',err_name='BIGDFT_RUNTIME_ERROR')
 if (ii/=natom) stop 'ii/=natom'
 allocate(grewtn_tmp(3,natom,2))

!Add up total charge and sum of $charge^2$ in cell
 chsq=0._dp
 ch=0._dp
 !$omp parallel default(none) shared(natom,zion, typat,ch,chsq) private(ia)
 !$omp do reduction(+:ch,chsq) schedule(static)
 do ia=1,natom
   ch=ch+zion(typat(ia))
   chsq=chsq+zion(typat(ia))**2
 end do
 !$omp end do
 !$omp end parallel

!Compute eta, the Ewald summation convergence parameter,
!for approximately optimized summations:
 direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
 recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
!A bias is introduced, because G-space summation scales
!better than r space summation ! Note : debugging is the most
!easier at fixed eta.
 eta=pi*200.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)

!Conduct reciprocal space summations
 fac=pi**2/eta
 gsum=0._dp
 grewtn(:,:)=0.0_dp

!Sum over G space, done shell after shell until all
!contributions are too small.
 ng=0
 do
   ng=ng+1
   newg=0

   do ig3=-ng,ng
     do ig2=-ng,ng
       do ig1=-ng,ng

!        Exclude shells previously summed over
         if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng&
&         .or. ng==1 ) then

!          gsq is G dot G = |G|^2
           gsq=gmet(1,1)*dble(ig1*ig1)+gmet(2,2)*dble(ig2*ig2)+&
&           gmet(3,3)*dble(ig3*ig3)+2._dp*(gmet(2,1)*dble(ig1*ig2)+&
&           gmet(3,1)*dble(ig1*ig3)+gmet(3,2)*dble(ig3*ig2))

!          Skip g=0:
           if (gsq>1.0d-20) then
             arg=fac*gsq

!            Larger arg gives 0 contribution because of exp(-arg)
             if (arg <= 80._dp) then
!              When any term contributes then include next shell
               newg=1
               term=exp(-arg)/gsq
               summr = 0.0_dp
               summi = 0.0_dp
!              Note that if reduced atomic coordinates xred drift outside
!              of unit cell (outside [0,1)) it is irrelevant in the following
!              term, which only computes a phase.
!              OCL SCALAR ! by MM for Fujitsu
               !$omp parallel default(none) &
               !$omp shared(natom,ig1,ig2,ig3,xred,zion,typat,summr,summi) private(ia,arg)
               !$omp do reduction(+:summr,summi) schedule(static)
               do ia=1,natom
                 arg=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
!                Sum real and imaginary parts (avoid complex variables)
                 summr=summr+zion(typat(ia))*cos(arg)
                 summi=summi+zion(typat(ia))*sin(arg)
               end do
               !$omp end do
               !$omp end parallel

!              The following two checks avoid an annoying
!              underflow error message
               if (abs(summr)<1.d-16) summr=0.0_dp
               if (abs(summi)<1.d-16) summi=0.0_dp

!              The product of term and summr**2 or summi**2 below
!              can underflow if not for checks above
               t1=term*(summr*summr+summi*summi)
               gsum=gsum+t1

!              OCL SCALAR ! by MM for Fujitsu
               !$omp parallel default(none) &
               !$omp shared(natom,ig1,ig2,ig3,xred,summr,summi,term,zion,typat,grewtn) &
               !$omp private(ia,arg,phr,phi,c1i)
               !$omp do schedule(static)
               do ia=1,natom
!                Again only phase is computed so xred may fall outside [0,1).
                 arg=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
                 phr= cos(arg)
                 phi=-sin(arg)
!                (note: do not need real part, commented out)
!                c1r=(phr*summr-phi*summi)*(term*zion(typat(ia)))
                 c1i=(phi*summr+phr*summi)*(term*zion(typat(ia)))
!                compute coordinate gradients
                 grewtn(1,ia)=grewtn(1,ia)-c1i*ig1
                 grewtn(2,ia)=grewtn(2,ia)-c1i*ig2
                 grewtn(3,ia)=grewtn(3,ia)-c1i*ig3
               end do
               !$omp end do
               !$omp end parallel

!              End condition of not larger than 80.0
             end if

!            End skip g=0
           end if

!          End triple loop over G s and associated new shell condition
         end if
       end do
     end do
   end do

!  Check if new shell must be calculated
   if (newg==0) exit

!  End the loop on ng (new shells). Note that there is one exit
!  from this loop.
 end do
!
 sumg=gsum/(two_pi*ucvol)

!Stress tensor is now computed elsewhere (abi_ewald2) hence do not need
!length scale gradients (used to compute them here).

!normalize coordinate gradients by unit cell volume ucvol
 term=-2._dp/ucvol
 grewtn(:,:)=grewtn(:,:)*term
!call DSCAL(3*natom,term,grewtn,1)

!Conduct real space summations
 reta=sqrt(eta)
 fac=2._dp*sqrt(eta/pi)
 sumr=0.0_dp

 grewtn_tmp(:,:,:)=0.0_dp

!In the following a summation is being conducted over all
!unit cells (ir1, ir2, ir3) so it is appropriate to map all
!reduced coordinates xred back into [0,1).
!
!Loop on shells in r-space as was done in g-space
 nr=0
 do
   nr=nr+1
   newr=0
!  
   do ir3=-nr,nr
     do ir2=-nr,nr
       do ir1=-nr,nr
         if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr&
&         .or. nr==1 )then

           do ia=1,natp!natom
             iia=isat+ia
!            Map reduced coordinate xred(mu,ia) into [0,1)
             fraca1=xred(1,iia)-aint(xred(1,iia))+0.5_dp-sign(0.5_dp,xred(1,iia))
             fraca2=xred(2,iia)-aint(xred(2,iia))+0.5_dp-sign(0.5_dp,xred(2,iia))
             fraca3=xred(3,iia)-aint(xred(3,iia))+0.5_dp-sign(0.5_dp,xred(3,iia))
             drdta1=0.0_dp
             drdta2=0.0_dp
             drdta3=0.0_dp
!            OCL SCALAR ! by MM for Fujitsu
             !$omp parallel default(none) &
             !$omp shared(natom,xred,fraca1,fraca2,fraca3,rmet,reta,zion,typat) &
             !$omp shared(iia,newr,sumr,drdta1,drdta2,drdta3,ir1,ir2,ir3,eta,fac) &
             !$omp private(ib,fracb1,fracb2,fracb3,r1,r2,r3,rsq,term,rmagn,arg,derfc_arg,r1a1d,r2a2d,r3a3d)
             !$omp do reduction(+:newr,sumr,drdta1,drdta2,drdta3)
             do ib=1,natom
               fracb1=xred(1,ib)-aint(xred(1,ib))+0.5_dp-sign(0.5_dp,xred(1,ib))
               fracb2=xred(2,ib)-aint(xred(2,ib))+0.5_dp-sign(0.5_dp,xred(2,ib))
               fracb3=xred(3,ib)-aint(xred(3,ib))+0.5_dp-sign(0.5_dp,xred(3,ib))
               r1=dble(ir1)+fracb1-fraca1
               r2=dble(ir2)+fracb2-fraca2
               r3=dble(ir3)+fracb3-fraca3
               rsq=rmet(1,1)*r1*r1+rmet(2,2)*r2*r2+rmet(3,3)*r3*r3+&
&               2.0_dp*(rmet(2,1)*r2*r1+rmet(3,2)*r3*r2+rmet(3,1)*r1*r3)

!              Avoid zero denominators in 'term':
               if (rsq>=1.0d-24) then

!                Note: erfc(8) is about 1.1e-29,
!                so do not bother with larger arg.
!                Also: exp(-64) is about 1.6e-28,
!                so do not bother with larger arg**2 in exp.
                 term=0._dp
                 if (eta*rsq<64.0_dp) then
                   newr=newr+1
                   rmagn=sqrt(rsq)
                   arg=reta*rmagn
!                  derfc is the real(dp) complementary error function
                   call abi_derfcf(derfc_arg,arg)
                   term=derfc_arg/rmagn
                   sumr=sumr+zion(typat(iia))*zion(typat(ib))*term
                   term=zion(typat(iia))*zion(typat(ib))*&
&                   (term+fac*exp(-eta*rsq))/rsq
!                  Length scale grads now handled with stress tensor in abi_ewald2
                   r1a1d=rmet(1,1)*r1+rmet(1,2)*r2+rmet(1,3)*r3
                   r2a2d=rmet(2,1)*r1+rmet(2,2)*r2+rmet(2,3)*r3
                   r3a3d=rmet(3,1)*r1+rmet(3,2)*r2+rmet(3,3)*r3
!                  Compute terms related to coordinate gradients
                   drdta1=drdta1+term*r1a1d
                   drdta2=drdta2+term*r2a2d
                   drdta3=drdta3+term*r3a3d
                 end if

!                End avoid zero denominators in'term'
               end if

!              end loop over ib:
             end do
             !$omp end do
             !$omp end parallel

             grewtn_tmp(1,iia,1)=grewtn_tmp(1,iia,1)+drdta1
             grewtn_tmp(2,iia,1)=grewtn_tmp(2,iia,1)+drdta2
             grewtn_tmp(3,iia,1)=grewtn_tmp(3,iia,1)+drdta3

!            end loop over ia:
           end do

!          end triple loop over real space points and associated condition of new shell
         end if
       end do
     end do
   end do

   if (nproc>1) then
     !call mpiallred(newr, mpi_sum, bigdft_mpi%mpi_comm)
     ii=0
     call mpi_allreduce(newr, ii, 1, &
          mpi_integer, mpi_sum, mpi_comm, ierr)
     newr=ii
   end if

!  Check if new shell must be calculated
   if(newr==0) exit

!  End loop on nr (new shells). Note that there is an exit within the loop
 end do

 if (nproc > 1 .and. natom > 0) then
   !call mpiallred(grewtn, mpi_sum, bigdft_mpi%mpi_comm)
   call mpi_allreduce(grewtn_tmp(1,1,1), grewtn_tmp(1,1,2), 3*natom, &
        mpi_double_precision, mpi_sum, mpi_comm, ierr)
   tt=0.0_dp
   call mpi_allreduce(sumr, tt, 1, &
        mpi_double_precision, mpi_sum, mpi_comm, ierr)
   sumr=tt
 else
   grewtn_tmp(:,:,2)=grewtn_tmp(:,:,1)
 end if
 grewtn(:,:)=grewtn(:,:)+grewtn_tmp(:,:,2)
!
 sumr=0.5_dp*sumr
 fac=pi*ch**2/(2.0_dp*eta*ucvol)

!Finally assemble Ewald energy, eew
 eew=sumg+sumr-chsq*reta/sqrt(pi)-fac

!Length scale grads handled with stress tensor, abi_ewald2

!Output the final values of ng and nr
 write(message, '(a,a,i4,a,i4)' )ch10,&
& ' abi_ewald : nr and ng are ',nr,' and ',ng
 call abi_wrtout(std_out,message,'COLL')

end subroutine abi_ewald
!!***
