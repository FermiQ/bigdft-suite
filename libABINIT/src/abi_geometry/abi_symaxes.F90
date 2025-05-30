!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_symaxes
!! NAME
!! abi_symaxes
!!
!! FUNCTION
!! Determines the type of symmetry operation, for
!! the proper symmetries 2,2_1,3,3_1,3_2,4,4_1,4_2,4_3,6,6_1,...6_5
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (RC, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! center=type of bravais lattice centering
!!   	  center=0        no centering
!!        center=-1       body-centered
!!        center=-3       face-centered
!!        center=1        A-face centered
!!        center=2        B-face centered
!!        center=3        C-face centered
!! iholohedry=type of holohedry
!!            iholohedry=1   triclinic      1bar
!!            iholohedry=2   monoclinic     2/m
!!            iholohedry=3   orthorhombic   mmm
!!            iholohedry=4   tetragonal     4/mmm
!!            iholohedry=5   trigonal       3bar m  (rhombohedral Bravais latt)
!!            iholohedry=6   hexagonal      6/mmm
!!            iholohedry=7   cubic          m3bar m
!! isym=number of the symmetry operation that is currently analyzed
!! isymrelconv=symrel matrix for the particular operation, in conv. axes
!! ordersym=order of the symmetry operation
!! tnons_order=order of the screw translation
!! trialt(3)=screw translation associated with the symmetry operation
!!           in conventional axes (all components in the range ]-1/2,1/2] )
!!
!! OUTPUT
!! type_axis=type of the symmetry operation
!!
!! NOTES
!! It is assumed that the symmetry operations will be entered in the
!! symrel tnonsconv arrays, for the CONVENTIONAL cell.
!! For proper symmetries (rotations), the
!! associated translation is determined.
!!
!! There is a subtlety with translations associated with rotations :
!! all the rotations with axis
!! parallel to the one analysed do not all have the
!! same translation characteristics. This is clearly seen
!! in the extended Hermann-Mauguin symbols, see the international
!! table for crystallography, chapter 4.
!! In the treatment that we adopt, one will distinguish
!! the cases of primitive Bravais lattices, and centered
!! bravais lattices. In the latter case, in the present routine,
!! at the exception of the trigonal axis for the
!! cubic system, we explicitely generate the correct ratio of different
!! translations, so that their type can be explicitely assigned,
!! without confusion. By contrast, for primitive lattices,
!! the "tnons" that has been transmitted to the present routine
!! might be one of the few possible translations vectors,
!! nearly at random. We deal with this case by the explicit
!! examination of the system classes, and the identification
!! of such a possibility. In particular:
!! (1) for the trigonal axis in the rhombohedral Bravais lattice,
!! or in the cubic system, there is an equal number of 3, 3_1,
!! and 3_2 axes parallel to each other, in a cell that
!! is primitive (as well as conventional). In this particular case,
!! in the present
!! routine, all 3, 3_1 and 3_2 axes are assigned to be 3 axes,
!! independently of the centering.
!! (2) for the 4- or 6- axes, no confusion is possible :
!! in the primitive cell, there is only one possible translation,
!! while in the centered cells, the correct ratio of translation
!! vectors will be generated
!! (3) for the binary axes, there is no problem when the cell
!! is centered, but there are problems
!! (3a) for the tP Bravais lattice, for an axis in a tertiary direction,
!! (see the description of the lattice symmetry directions
!!  table 2.4.1 of the international tables for crystallography),
!!  where the family of axes is made equally of 2 and 2_1 axis.
!!  In this case, we attribute the binary axis to the specific class
!!  of "tertiary 2-axis". We keep track of the 2 or 2_1
!!  characteristics of all other binary axes
!! (3b) for the tI Bravais lattice, in all the directions,
!!  there is an equal number of 2 and 2_1 axes. We distinguish
!!  the primary and secondary family from the tertiary family.
!! (3c) for the hP Bravais lattice, each binary axis can present
!!  no translation or be a screw axis (in the same direction).
!!  For primary axes, one need the "2" and "2_1" classification,
!!  while for secondary and tertiary axes, the associated
!!  translation vector will have not importance.
!!  However, one will need to distinguish secondary from
!!  tertiary, and these from primary axes.
!!  So, this is the most complicated case, for binary axes,
!!  with the following sets of binary axes : "2", "2_1",
!!  "secondary 2" and "tertiary 2".
!! (3d) for the hR Bravais lattice, each binary axis can present
!!  no translation or be a screw axis (in the same direction).
!!  There is no distinction between tertiary axes and other, so that
!!  we simply assign a binary axis to "2-axis"
!! (3e) for the cP lattice, the binary axes along tertiary directions
!!  can also have different translation vectors, while for the primary
!!  direction, there is no such ambiguity. So, we will attribute
!!  tertiary 2 axis to the "tertiary 2-axis" set (there are always 6),
!!  and attribute 2 and 2_1 primary axes to the corresponding sets.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_symaxes(center,iholohedry,&
&isym,isymrelconv,ordersym,tnons_order,trialt,type_axis)

 use abi_defs_basis
 use abi_interfaces_lowlevel

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: center,iholohedry,isym,ordersym,tnons_order
 integer,intent(out) :: type_axis
!arrays
 integer,intent(in) :: isymrelconv(3,3)
 real(dp),intent(in) :: trialt(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: message
 integer :: direction,directiontype
 real(dp),parameter :: nzero=1.0d-6

!**************************************************************************

 select case(ordersym)

   case(2)                       ! point symmetry 2
!    Must characterize directiontype for cP, tP, tI, and hP Bravais lattices
     directiontype=1
     if( iholohedry==4 .or. iholohedry==7) then ! tP or cP Bravais lattices
       if(abs(isymrelconv(1,1))+ &
&       abs(isymrelconv(2,2))+ &
&       abs(isymrelconv(3,3))  ==1) directiontype=3
     else if(iholohedry==6)then   ! hP Bravais lattice
       if(sum(isymrelconv(:,:))/=-1 )directiontype=2
       if(sum(isymrelconv(:,:))==0 .or. sum(isymrelconv(:,:))==-3 )&
&       directiontype=3
!      directiontype=1 corresponds to a primary axis
!      directiontype=2 corresponds to a tertiary axis
!      directiontype=3 corresponds to a secondary axis
     end if

!    Now, classify the 2 axes
     if(directiontype==2)then
       type_axis=4                 ! secondary 2  (only in the hP Bravais latt case)
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a secondary 2-axis '

     else if(directiontype==3 .and. iholohedry==4)then
       type_axis=21                ! tertiary 2
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a tertiary 2-axis '
     else if(directiontype==3 .and. &
&       center==0 .and. (iholohedry==6.or.iholohedry==7) )then
       type_axis=21                ! tertiary 2
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a tertiary 2-axis '
     else if(tnons_order==1 .or. (iholohedry==4 .and. center==-1) .or. &
&       iholohedry==5)then
       type_axis=9                 ! 2
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 2-axis '
     else
       type_axis=20                ! 2_1
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 2_1-axis '
     end if
     call abi_wrtout(6,message,'COLL')

   case(3)                       ! point symmetry 3
     if(tnons_order==1)then
       type_axis=10                ! 3
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 3-axis '
     else if(iholohedry==5 .or. iholohedry==7)then
!      This is a special situation : in the same family of parallel 3-axis,
!      one will have an equal number of 3, 3_1 and 3_2 axes, so that
!      it is non-sense to try to classify one of them.
       type_axis=10                ! 3, 3_1 or 3_2, undistinguishable
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 3, 3_1 or 3_2 axis '
     else
!      Must recognize 3_1 or 3_2
       if(isymrelconv(1,1)==0)then  ! 3+
         if(abs(trialt(3)-third)<nzero)type_axis=22   ! 3_1
         if(abs(trialt(3)+third)<nzero)type_axis=23   ! 3_2
       else if(isymrelconv(1,1)==-1)then  ! 3-
         if(abs(trialt(3)-third)<nzero)type_axis=23   ! 3_2
         if(abs(trialt(3)+third)<nzero)type_axis=22   ! 3_1
       end if
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 3_1 or 3_2-axis '
     end if
     call abi_wrtout(6,message,'COLL')

   case(4)                       ! point symmetry 4
     if(tnons_order==1)then
       type_axis=12                ! 4
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 4-axis '
     else if(tnons_order==2)then
       type_axis=25                ! 4_2
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 4_2-axis '
     else if(center/=0)then
       type_axis=24                ! 4_1 or 4_3
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 4_1 or 4_3-axis '
     else
!      Must recognize 4_1 or 4_3, along the three primary directions
       do direction=1,3
         if(isymrelconv(direction,direction)==1)then  !
           if( (direction==1 .and. isymrelconv(2,3)==-1) .or. &
&           (direction==2 .and. isymrelconv(3,1)==-1) .or. &
&           (direction==3 .and. isymrelconv(1,2)==-1)       )then ! 4+
             if(abs(trialt(direction)-quarter)<nzero)type_axis=24    ! 4_1
             if(abs(trialt(direction)+quarter)<nzero)type_axis=26    ! 4_3
           else if( (direction==1 .and. isymrelconv(2,3)==1) .or. &
&             (direction==2 .and. isymrelconv(3,1)==1) .or. &
&             (direction==3 .and. isymrelconv(1,2)==1)       )then ! 4-
             if(abs(trialt(direction)-quarter)<nzero)type_axis=26    ! 4_3
             if(abs(trialt(direction)+quarter)<nzero)type_axis=24    ! 4_1
           end if
         end if
       end do
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 4_1 or 4_3-axis '
     end if
     call abi_wrtout(6,message,'COLL')

   case(6)                       ! point symmetry 6
     if(tnons_order==1)then
       type_axis=14                ! 6
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 6-axis '
     else if(tnons_order==2)then
       type_axis=29                ! 6_3
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 6_3-axis '
     else if(tnons_order==3)then
!      Must recognize 6_2 or 6_4
       if(isymrelconv(1,1)==1)then  ! 6+
         if(abs(trialt(3)-third)<nzero)type_axis=28   ! 6_2
         if(abs(trialt(3)+third)<nzero)type_axis=30   ! 6_4
       else if(isymrelconv(1,1)==0)then  ! 6-
         if(abs(trialt(3)-third)<nzero)type_axis=30   ! 6_4
         if(abs(trialt(3)+third)<nzero)type_axis=28   ! 6_2
       end if
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 6_2 or 6_4-axis '
     else
!      Must recognize 6_1 or 6_5
       if(isymrelconv(1,1)==1)then  ! 6+
         if(abs(trialt(3)-sixth)<nzero)type_axis=27   ! 6_1
         if(abs(trialt(3)+sixth)<nzero)type_axis=31   ! 6_5
       else if(isymrelconv(1,1)==0)then  ! 6-
         if(abs(trialt(3)-sixth)<nzero)type_axis=31   ! 6_5
         if(abs(trialt(3)+sixth)<nzero)type_axis=27   ! 6_1
       end if
       write(message,'(a,i3,a)') &
&       ' abi_symaxes : the symmetry operation no. ',isym,' is a 6_1 or 6_5-axis '
     end if
     call abi_wrtout(6,message,'COLL')

 end select

end subroutine abi_symaxes
!!***
