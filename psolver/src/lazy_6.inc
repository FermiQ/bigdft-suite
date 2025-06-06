!> @file
!!   Filters for interpolating scaling functions (order 6)
!! @author
!!    Copyright (C) 2015-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

integer, parameter :: m=8
real(kind=8), dimension(-m:m) ::  ch,cg,cht,cgt

!******** coefficients for wavelet transform *********************
do i=-m,m
   ch(i)=0.d0
   cht(i)=0.d0
   cg(i)=0.d0
   cgt(i)=0.d0
enddo

! The normalization is chosen such that a constant function remains the same constant 
! on each level of the transform

ch(-5)=3.d0/256.d0
ch(-4)=0.d0
ch(-3)=-25.d0/256.d0
ch(-2)=0.d0
ch(-1)=75.d0/128.d0
ch( 0)=1.d0
ch( 1)=75.d0/128.d0
ch( 2)=0.d0
ch( 3)=-25.d0/256.d0
ch( 4)=0.d0
ch( 5)=3.d0/256.d0
! 
cht( 0)=1.d0

! g coefficients from h coefficients
do i=-m,m-1
   cg(i+1)=cht(-i)*(-1.d0)**(i+1)
   cgt(i+1)=ch(-i)*(-1.d0)**(i+1)
enddo

