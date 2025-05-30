!> @file
!!   Filters for interpolating scaling functions (order 8)
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

integer, parameter :: m=10
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

ch(-7)=-5.d0/2048.d0
ch(-6)=0.d0
ch(-5)=49.d0/2048.d0
ch(-4)=0.d0
ch(-3)=-245.d0/2048.d0
ch(-2)=0.d0
ch(-1)=1225.d0/2048.d0
ch( 0)=1.d0
ch( 1)=1225.d0/2048.d0
ch( 2)=0.d0
ch( 3)=-245.d0/2048.d0
ch( 4)=0.d0
ch( 5)=49.d0/2048.d0
ch( 6)=0.d0
ch( 7)=-5.d0/2048.d0
! 
!cht( 0)=1.d0
cht(-8)=5.d0/8192.d0
cht(-7)=0.d0
cht(-6)=-11.d0/2048.d0
cht(-5)=0.d0
cht(-4)=49.d0/2048.d0
cht(-3)=0.d0
cht(-2)=-245.d0/2048.d0
cht(-1)=1.d0/4.d0
cht( 0)=2871.d0/4096.d0
cht(1)=1.d0/4.d0
cht(2)=-245.d0/2048.d0
cht(3)=0.d0
cht(4)=49.d0/2048.d0
cht(5)=0.d0
cht(6)=-11.d0/2048.d0
cht(7)=0.d0
cht(8)=5.d0/8192.d0

! g coefficients from h coefficients
do i=-m,m-1
   cg(i+1)=cht(-i)*(-1.d0)**(i+1)
   cgt(i+1)=ch(-i)*(-1.d0)**(i+1)
enddo

