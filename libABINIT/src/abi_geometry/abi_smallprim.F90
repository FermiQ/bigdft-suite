!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_smallprim
!! NAME
!! abi_smallprim
!!
!! FUNCTION
!! Find the smallest possible primitive vectors for an input lattice
!! This algorithm is not as restrictive as the conditions mentioned at p.740
!! of the international tables for crystallography (1983).
!! The final vectors form a right-handed basis, while their
!! sign and ordering is chosen such as to maximize the overlap
!! with the original vectors in order.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  rprimd(3,3)=primitive vectors
!!
!! OUTPUT
!!  metmin(3,3)=metric for the new (minimal) primitive vectors
!!  minim(3,3)=minimal primitive translations
!!
!! NOTES
!! The routine might as well be defined without
!! metmin as argument, but it is more convenient to have it
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_smallprim(metmin,minim,rprimd)

 use abi_defs_basis
 use abi_interfaces_lowlevel
 use abi_interfaces_geometry, except_this_one => abi_smallprim

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: metmin(3,3),minim(3,3)

!Local variables-------------------------------
!scalars
 integer :: ia,ib,ii,itrial,minimal
 real(dp) :: determinant,length2,metsum,ucvol
 character(len=500) :: message
!arrays
 integer :: nvecta(3),nvectb(3)
 real(dp) :: gprimd(3,3),gmet(3,3),rmet(3,3),scprod(3),tmpvect(3)

!**************************************************************************

 call abi_metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 nvecta(1)=2 ; nvectb(1)=3
 nvecta(2)=1 ; nvectb(2)=3
 nvecta(3)=1 ; nvectb(3)=2

 minim(:,:)=rprimd(:,:)
 metmin(:,:)=rmet(:,:)

!Note this loop without index
 do

!  Will exit if minimal=1 is still valid after a trial
!  to reduce the vectors of each of the three pairs
   minimal=1

   do itrial=1,3

     ia=nvecta(itrial) ; ib=nvectb(itrial)
!    Make sure the scalar product is negative
     if(metmin(ia,ib)>tol8)then
       minim(:,ia)=-minim(:,ia)
       metmin(ia,ib)=-metmin(ia,ib) ; metmin(ib,ia)=-metmin(ib,ia)
       metmin(ia,itrial)=-metmin(ia,itrial)
       metmin(itrial,ia)=-metmin(itrial,ia)
     end if
!    Compute the length of the sum vector
     length2=metmin(ia,ia)+2*metmin(ia,ib)+metmin(ib,ib)
!    Replace the first vector by the sum vector if the latter is smaller
     if(length2<metmin(ia,ia)-tol8)then
       minim(:,ia)=minim(:,ia)+minim(:,ib)
       metmin(ia,ia)=length2
       metmin(ia,ib)=metmin(ia,ib)+metmin(ib,ib)
       metmin(ia,itrial)=metmin(ia,itrial)+metmin(ib,itrial)
       metmin(ib,ia)=metmin(ia,ib)
       metmin(itrial,ia)=metmin(ia,itrial)
       minimal=0
!      Replace the second vector by the sum vector if the latter is smaller
     else if(length2<metmin(ib,ib)-tol8)then
       minim(:,ib)=minim(:,ia)+minim(:,ib)
       metmin(ib,ib)=length2
       metmin(ia,ib)=metmin(ia,ib)+metmin(ia,ia)
       metmin(itrial,ib)=metmin(itrial,ib)+metmin(itrial,ia)
       metmin(ib,ia)=metmin(ia,ib)
       metmin(ib,itrial)=metmin(itrial,ib)
       minimal=0
     end if

   end do

   if(minimal==1)exit

 end do

!At this stage, the three vectors have angles between each other that are
!comprised between 90 and 120 degrees. It might still be that minus the vector
!that is the sum of the three vectors is smaller than the longest of these vectors
 do
   
!  Will exit if minimal=1 is still valid after a trial
!  to replace each of the three vectors by minus the summ of the three vectors
   minimal=1
   metsum=sum(metmin(:,:))
   do itrial=1,3
     ia=nvecta(itrial) ; ib=nvectb(itrial)
     if(metmin(ia,ia)>metsum+tol8)then
       minim(:,ia)=-minim(:,1)-minim(:,2)-minim(:,3)
       metmin(ia,ib)=-sum(metmin(:,ib))
       metmin(ia,itrial)=-sum(metmin(:,itrial))
       metmin(ia,ia)=metsum
       metmin(ib,ia)=metmin(ia,ib)
       metmin(itrial,ia)=metmin(ia,itrial)
       minimal=0      
     end if
   end do

   if(minimal==1)exit
   
 end do

!Choose the first vector
!Compute the scalar product of the three minimal vectors
!with the first original vector
 scprod(:)=zero
 do ii=1,3
   scprod(:)=scprod(:)+minim(ii,:)*rprimd(ii,1)
 end do
!Determine the vector with the maximal absolute overlap
 itrial=1
 if(abs(scprod(2))>abs(scprod(1))+tol8)itrial=2
 if(abs(scprod(3))>abs(scprod(itrial))+tol8)itrial=3
!Switch the vectors if needed
 if(itrial/=1)then
   tmpvect(:)=minim(:,1)
   minim(:,1)=minim(:,itrial)
   minim(:,itrial)=tmpvect(:)
 end if
!Choose the sign
 if(scprod(itrial)<tol8)minim(:,1)=-minim(:,1)

!Choose the second vector
!Compute the scalar product of the second and third minimal vectors
!with the second original vector
 scprod(2:3)=zero
 do ii=1,3
   scprod(2:3)=scprod(2:3)+minim(ii,2:3)*rprimd(ii,2)
 end do
!Determine the vector with the maximal absolute overlap
 itrial=2
 if(abs(scprod(3))>abs(scprod(2))+tol8)itrial=3
!Switch the vectors if needed
 if(itrial/=2)then
   tmpvect(:)=minim(:,2)
   minim(:,2)=minim(:,itrial)
   minim(:,itrial)=tmpvect(:)
 end if
!Choose the sign
 if(scprod(itrial)<tol8)minim(:,2)=-minim(:,2)

!Change sign of the third vector if not right-handed basis
 determinant=minim(1,1)*(minim(2,2)*minim(3,3)-minim(3,2)*minim(2,3))+&
& minim(2,1)*(minim(3,2)*minim(1,3)-minim(1,2)*minim(3,3))+&
& minim(3,1)*(minim(1,2)*minim(2,3)-minim(2,2)*minim(1,3))
 if(determinant<-tol8)minim(:,3)=-minim(:,3)
 if(abs(determinant)<tol8)then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' abi_smallprim : BUG -',ch10,&
&   '  minim gives vanishing unit cell volume.'
   call abi_wrtout(std_out,message,'COLL')
   call abi_leave_new('COLL')
 end if

!Final computation of metmin
 do ii=1,3
   metmin(ii,:)=minim(1,ii)*minim(1,:)+&
&   minim(2,ii)*minim(2,:)+&
&   minim(3,ii)*minim(3,:)
 end do

end subroutine abi_smallprim
!!***
