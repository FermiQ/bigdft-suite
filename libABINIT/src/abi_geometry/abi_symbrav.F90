!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_symbrav
!! NAME
!! abi_symbrav
!!
!! FUNCTION
!! From the list of symmetry operations, and the lattice vectors,
!! determine the Bravais information (including the holohedry, the centering,
!! the coordinate of the primitive vectors in the conventional vectors), 
!! as well as the point group.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! msym=dimension of symrel
!! nsym=actual number of symmetries
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! symrel(3,3,msym)=symmetry operations in real space in terms
!!                  of primitive translations
!! tolsym=tolerance for the symmetries
!!
!! OUTPUT
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprimd in the axes
!!              of the conventional bravais lattice (*2 if center/=0)
!! ptgroup=symmetry point group
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_symbrav(bravais,msym,nsym,ptgroup,rprimd,symrel,tolsym)

 use abi_defs_basis
 use abi_interfaces_lowlevel
 use abi_interfaces_numeric
 use abi_interfaces_geometry, except_this_one => abi_symbrav

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,nsym
 real(dp),intent(in) :: tolsym
 character(len=5),intent(out) :: ptgroup
!arrays
 integer,intent(in) :: symrel(3,3,msym)
 integer,intent(out) :: bravais(11)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: iaxis,ii,bravais1now,ideform,iholohedry,invariant,isym
 integer :: jaxis,next_stage,nptsym,problem
 real(dp) :: norm,scprod
 character(len=500) :: message
!arrays
 integer :: identity(3,3),axis_trial(3),hexa_axes(3,7),ortho_axes(3,13)
 integer,allocatable :: ptsymrel(:,:,:),symrelconv(:,:,:)
 real(dp) :: axes(3,3),axis_cart(3),axis_red(3)
 real(dp) :: rprimdconv(3,3),rprimdtry(3,3),rprimdnow(3,3)
 real(dp) :: rprimdconv_invt(3,3)

!**************************************************************************

 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1

 ortho_axes(:,:)=0
 ortho_axes(1,1)=1
 ortho_axes(2,2)=1
 ortho_axes(3,3)=1
 ortho_axes(:,4)=(/0,1,1/)
 ortho_axes(:,5)=(/1,0,1/)
 ortho_axes(:,6)=(/1,1,0/)
 ortho_axes(:,7)=(/0,1,-1/)
 ortho_axes(:,8)=(/-1,0,1/)
 ortho_axes(:,9)=(/1,-1,0/)
 ortho_axes(:,10)=(/1,1,1/)
 ortho_axes(:,11)=(/-1,1,1/)
 ortho_axes(:,12)=(/1,-1,1/)
 ortho_axes(:,13)=(/1,1,-1/)

 hexa_axes(:,:)=0
 hexa_axes(1,1)=1
 hexa_axes(2,2)=1
 hexa_axes(3,3)=1
 hexa_axes(:,4)=(/1,-1,0/)
 hexa_axes(:,5)=(/2,1,0/)
 hexa_axes(:,6)=(/1,1,0/)
 hexa_axes(:,7)=(/1,2,0/)

!Determine the point group from the list of symmetry operations.
!Also determine the holohedry, up to one undeterminacy : hR versus hP
 call abi_symptgroup(iholohedry,nsym,ptgroup,symrel)

!Loop over trial deformations
!This is needed in case the Bravais lattice determination from the lattice vectors
!has a higher holohedry than the real one, in which the symmetry
!operations for the atoms (or electric field, etc) are taken into account
 iaxis=0
 invariant=0
 next_stage=0
 rprimdnow(:,:)=rprimd(:,:)
 rprimdtry(:,:)=rprimd(:,:)
 allocate(symrelconv(3,3,nsym))

!At most will have to try 65 deformations (13 axes, five stages)
 do ideform=1,65

   allocate(ptsymrel(3,3,msym))
   call abi_symlatt(bravais,msym,nptsym,ptsymrel,rprimdtry,tolsym)
   deallocate(ptsymrel)

!  Examine the agreement with bravais(1)
!  Warning : might change Bravais lattice hR to hP, if hexagonal axes
   problem=0
   select case (bravais(1))
     case(7)
       if(iholohedry<6)problem=1
       if(iholohedry==6)problem=2
     case(6)
       if(iholohedry<4)problem=1
       if(iholohedry==7 .or. iholohedry==4)problem=2
!      Here, change hR into hP
       if(iholohedry==5)iholohedry=6
     case(5)
       if(iholohedry<4)problem=1
       if(iholohedry==7 .or. iholohedry==6 .or. iholohedry==4)problem=2
     case(4)
       if(iholohedry<4)problem=1
       if(iholohedry>4)problem=2
     case(3)
       if(iholohedry<3)problem=1
       if(iholohedry>3)problem=2
     case(2)
       if(iholohedry<2)problem=1
       if(iholohedry>2)problem=2
     case(1)
       if(iholohedry>1)problem=2
   end select

!  This is the usual situation, in which the lattice belong to the same holohedry 
!  as the lattice+atoms (+electric field + ...)
   if(problem==0)exit

  if(problem==2)then
     if(iaxis==0)then
       write(message, '(6a,i3,3a,i3,7a)' )ch10,&
&       ' abi_symbrav : BUG -',ch10,&
&       '  The Bravais lattice determined only from the primitive',ch10,&
&       '  vectors (rprim or angdeg), bravais(1)=',bravais(1),', is not compatible',ch10,&
&       '  with the real one, iholohedry=',iholohedry,', obtained by taking into',ch10,&
&       '  account the symmetry operations. This might be due to an insufficient',ch10,&
&       '  number of digits in the specification of rprim (at least 10),',ch10,&
&       '  or to an erroneous rprim or angdeg. If this is not the case, then ...'
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     if(iaxis==1)then
       write(message, '(6a,3i3,2a,i3,2a,i3)' )ch10,&
&       ' abi_symbrav : BUG -',ch10,&
&       '  Could not succeed to determine the bravais lattice',ch10,&
&       '  problem,iaxis,invariant=',problem,iaxis,invariant,ch10,&
&       '  bravais(1)=',bravais(1),ch10,&
&       '  iholohedry=',iholohedry
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
   end if

   if(problem==1)then
!    One is left with the problem=1 case, basically iholohedry is lower than bravais(1)
     if(iaxis==0)then
       write(message, '(a,a,a,a,a,a,i3,a,a,a,i3,a,a,a)' )ch10,&
&       ' abi_symbrav : COMMENT -',ch10,&
&       '  The Bravais lattice determined only from the primitive',ch10,&
&       '  vectors, bravais(1)=',bravais(1),', is more symmetric',ch10,&
&       '  than the real one, iholohedry=',iholohedry,', obtained by taking into',ch10,&
&       '  account the atomic positions. Start deforming the primitive vector set.'
       call abi_wrtout(std_out,message,'COLL')
       next_stage=1
     else if(iaxis/=0)then
       if(bravais(1)<bravais1now)then
         write(message, '(6a,i3,3a,i3,2a)' )ch10,&
&         ' abi_symbrav : COMMENT -',ch10,&
&         '  The Bravais lattice determined from modified primitive',ch10,&
&         '  vectors, bravais(1)=',bravais(1),', has a lower symmetry than before,',ch10,&
&         '  but is still more symmetric than the real one, iholohedry=',iholohedry,ch10,&
&         '  obtained by taking into account the atomic positions.'
         call abi_wrtout(std_out,message,'COLL')
         next_stage=1
       else if(iaxis==1)then
         write(message, '(6a,3i3,2a,i3,2a,i3)' )ch10,&
&         ' abi_symbrav : BUG -',ch10,&
&         '  Could not succeed to determine the bravais lattice',ch10,&
&         '  problem,iaxis,invariant=',problem,iaxis,invariant,ch10,&
&         '  bravais(1)=',bravais(1),ch10,&
&         '  iholohedry=',iholohedry
         call abi_wrtout(std_out,message,'COLL')
         call abi_leave_new('COLL')
       end if
     end if
   end if ! problem==1

   if(next_stage==1)then
     bravais1now=bravais(1)
     rprimdnow(:,:)=rprimdtry(:,:)
!    Generate the symmetry operations in the conventional vector coordinates
     rprimdconv(:,1)=bravais(3:5)
     rprimdconv(:,2)=bravais(6:8)
     rprimdconv(:,3)=bravais(9:11)
     axes(:,:)=zero
     axes(1,1)=one ; axes(2,2)=one ; axes(3,3)=one
     symrelconv(:,:,1:nsym)=symrel(:,:,1:nsym)
     call abi_symrelrot(nsym,rprimdconv,axes,symrelconv,tolsym)
     if(bravais(1)/=6)then
       iaxis=14
     else
       iaxis=8
     end if
     next_stage=0
   end if

   iaxis=iaxis-1
   do jaxis=iaxis,1,-1
     if(bravais(1)/=6)then
       axis_trial(:)=ortho_axes(:,jaxis)
     else
       axis_trial(:)=hexa_axes(:,jaxis)
     end if
     invariant=1
!    Examine whether all symmetry operations leave the axis invariant (might be reversed, though)
     do isym=1,nsym
       if(sum(abs(matmul(symrelconv(:,:,isym),axis_trial)-axis_trial(:)))/=0 .and. &
&       sum(abs(matmul(symrelconv(:,:,isym),axis_trial)+axis_trial(:)))/=0 )invariant=0
     end do
     if(invariant==1)then
       iaxis=jaxis
       write(message, '(2a,i3)' )ch10,&
&       ' abi_symbrav : found invariant axis, jaxis=',iaxis
       call abi_wrtout(std_out,message,'COLL')
       exit
     end if
   end do

   if(invariant==0)then
!    Not a single axis was invariant with respect to all operations ?!
     write(message, '(6a,3i3,2a,i3,2a,i3)' )ch10,&
&     ' abi_symbrav : BUG -',ch10,&
&     '  Could not succeed to determine the bravais lattice (not a single invariant)',ch10,&
&     '  problem,iaxis,invariant=',problem,iaxis,invariant,ch10,&
&     '  bravais(1)=',bravais(1),ch10,&
&     '  iholohedry=',iholohedry
     call abi_wrtout(std_out,message,'COLL')
     call abi_leave_new('COLL')
   end if

   call abi_matr3inv(rprimdconv,rprimdconv_invt)
   axis_red(:)=axis_trial(1)*rprimdconv_invt(1,:)+ &
&   axis_trial(2)*rprimdconv_invt(2,:)+ &
&   axis_trial(3)*rprimdconv_invt(3,:)
   axis_cart(:)=axis_red(1)*rprimdnow(:,1)+ &
&   axis_red(2)*rprimdnow(:,2)+ &
&   axis_red(3)*rprimdnow(:,3)
   norm=sum(axis_cart(:)**2)
!  Expand by a uniform, quite arbitrary, dilatation, along the invariant axis
!  Note : make these dilatation different, according to ideform 
   do ii=1,3
     scprod=axis_cart(1)*rprimdnow(1,ii)+axis_cart(2)*rprimdnow(2,ii)+axis_cart(3)*rprimdnow(3,ii) 
     rprimdtry(:,ii)=rprimdnow(:,ii)+ideform*(tol3-tol6)*scprod/norm*axis_cart(:)
   end do

 end do ! ideform

 if(bravais(1)/=iholohedry)then
   write(message, '(6a,3i3,2a,i3,2a,i3)' )ch10,&
&   ' abi_symbrav : BUG -',ch10,&
&   '  Despite efforts, Could not succeed to determine the bravais lattice :',ch10,&
&   '  bravais(1)=',bravais(1),ch10,&
&   '  iholohedry=',iholohedry
   call abi_wrtout(std_out,message,'COLL')
   call abi_leave_new('COLL')
 end if

 deallocate(symrelconv)

end subroutine abi_symbrav
!!***
