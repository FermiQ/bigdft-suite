!> @file
!!  Common convolutions
!! @author
!!    Copyright (C) 2018-2018 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!>  A analysis wavelet transformation where the size of the data is forced to shrink
!!  The input array y is overwritten
subroutine analyse_wire_self(n1,n2,n3,y,x)
  implicit none
  integer,intent(in)::n1,n2,n3
  real(kind=8),dimension((2*n1+16)*(2*n2+16)*(2*n3+2))::x,y
  integer nt

  ! i1,I2,i3 -> I2,i3,i1
  nt=(2*n2+16)*(2*n3+2)
  call  ana_rot_shrink(n1,nt,y,x)
  ! I2,i3,i1 -> i3,i1,i2
  nt=(2*n3+2)*(2*n1+16)
  call  ana_rot_shrink(n2,nt,x,y)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_per(n3,nt,y,x)

END SUBROUTINE analyse_wire_self


!>   A synthesis wavelet transformation where the size of the data is allowed to grow
!!   The input array x is overwritten
subroutine synthese_wire_self(n1,n2,n3,x,y)
  implicit none
  integer,intent(in)::n1,n2,n3
  real(kind=8),dimension((2*n1+16)*(2*n2+16)*(2*n3+2))::x,y
  integer nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_grow(n1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,I2
  nt=(2*n3+2)*(2*n1+16)
  call  syn_rot_grow(n2,nt,y,x)
  ! i3,i1,I2  -> i1,I2,i3
  nt=(2*n1+16)*(2*n2+16)
  call  syn_rot_per(n3,nt,x,y)

END SUBROUTINE synthese_wire_self


!>   Applies the magic filter matrix in wire BC ( no transposition)
!!   The input array x is overwritten
!!   this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_wire_self(n1,n2,n3,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(-7:n1+8,-7:n2+8,0:n3), intent(inout) :: x
  real(wp), dimension(-7:n1+8,-7:n2+8,0:n3), intent(inout) :: y
  !local variables
  integer :: ndat

  !  (i1,i2*i3) -> (i2*i3,i1)
  ndat=(n2+1)*(n3+1)
  call convrot_grow(n1,ndat,x,y)
  !  (i2,i3*i1) -> (i3*i1,I2)
  ndat=(n3+1)*(n1+16)
  call convrot_grow(n2,ndat,y,x)
  !  (i3,i1*I2) -> (i1*I2,i3)
  ndat=(n1+16)*(n2+16)
  call convrot_n_per(n3,ndat,x,y)

END SUBROUTINE convolut_magic_n_wire_self


!>   Applies the magic filter matrix in periodic BC ( no transposition)
!!   The input array x is not overwritten
!!   this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_wire(n1,n2,n3,x,y,ww)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(-7:n1+8,-7:n2+8,0:n3), intent(inout) :: y
  !local variables
  integer :: ndat
  real(wp), dimension(-7:n1+8,-7:n2+8,0:n3):: ww ! work array

  !  (i1,i2*i3) -> (i2*i3,i1)
  ndat=(n2+1)*(n3+1)
  call convrot_grow(n1,ndat,x,y)
  !  (i2,i3*i1) -> (i3*i1,I2)
  ndat=(n3+1)*(n1+16)
  call convrot_grow(n2,ndat,y,ww)
  !  (i3,i1*I2) -> (i1*I2,i3)
  ndat=(n1+16)*(n2+16)
  call convrot_n_per(n3,ndat,ww,y)

END SUBROUTINE convolut_magic_n_wire


!>   Applies the magic filter matrix transposed in periodic BC 
!!   The input array x is overwritten
!!   this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_t_wire_self(n1,n2,n3,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(-7:n1+8,-7:n2+8,0:n3), intent(inout) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer :: ndat

  !  (i1,I2*i3) -> (I2*i3,i1)
  ndat=(n2+16)*(n3+1)
  call convrot_shrink(n1,ndat,x,y)
  !  (I2,i3*i1) -> (i3*i1,i2)
  ndat=(n3+1)*(n1+16)
  call convrot_shrink(n2,ndat,y,x)
  !  (i3,i1*i2) -> (i1*i2,i3)
  ndat=(n1+1)*(n2+1)
  call convrot_t_per(n3,ndat,x,y)


END SUBROUTINE convolut_magic_t_wire_self
