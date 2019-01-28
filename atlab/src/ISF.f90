!> @file
!!  Routines which define and use scaling functions
!! @author
!! Copyright (C) 2002-2017 BigDFT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS

!> Calculate the values of a scaling function in a real uniform grid
subroutine ISF_family(itype,nmoms,nd,nrange,a,x)
  use f_utils, only: f_open_file,f_close,f_zero
  use yaml_strings, only: yaml_toa
  use dictionaries, only: f_err_throw
  use dynamic_memory

  implicit none
  !Arguments
  integer, intent(in) :: itype  !< Type of interpolating functions
  !> Number of moments of the lifted dual function
  !! to be preserved. If this value is different from
  !! 0, then the dual scaling function is given as output
  integer, intent(in) :: nmoms
  integer, intent(in) :: nd      !< Number of points: must be a power of 2!!
  integer, intent(out) :: nrange !< Range of the ISF [-nrange/2+1,nrange/2+1]
  real(kind=8), dimension(0:nd), intent(out) :: a !< Abscissae
  real(kind=8), dimension(0:nd), intent(out) :: x !< Values of the ISF
  !Local variables
  character(len=*), parameter :: subname='scaling_function'
  integer, parameter :: m_max=200 !< Maximum length of the filters
  integer :: i,nt,ni,m !, m, j
  !real(kind=8) :: mom
  real(kind=8), dimension(:), allocatable :: y
  double precision, dimension(-m_max:m_max) :: ch,cg !< Filters

  !Give the range of the scaling function
  !from -itype to itype for the direct case, more general for the dual
  call get_isf_family(m_max,itype,nmoms,m,ni,ch,cg)
  nrange = ni

  y = f_malloc(0.to.nd,id='y')

  ! plot scaling function
  call f_zero(x)
  call f_zero(y)
  nt=ni
  x(nt/2-1)=1.d0
  !x(nt+nt/2-1)=1.d0 !delta function: to obtain the wavelet
  !Then we iterate up to the corresponding number of points
  loop1: do
     nt=2*nt
     call back_trans(m,ch(-m),cg(-m),nd,nt,x,y)
     do i=0,nt-1
        x(i)=y(i)
     end do
     if (nt == nd) exit loop1
  end do loop1

  !Plot the interpolating scaling function if needed
!!$  unt=1
!!$  call f_open_file(unt,file='scfunction')
!!$  !open (unit=1,file='scfunction',status='unknown')
  do i=0,nd
     a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
!!$     write(unt,*) a(i),x(i)
  end do
!!$  call f_close(unt)
!!$  !now calculate the moments
!!$  call f_open_file(unt,file='scf_moments'//&
!!$       trim(adjustl(yaml_toa(itype)))//'-'//&
!!$       trim(adjustl(yaml_toa(nmoms))))
!!$  do j=0,itype
!!$     mom=0.d0
!!$     do i=0,nd
!!$        a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
!!$        mom=mom+a(i)**j/real(ni,kind=8)*x(i)
!!$     end do
!!$     write(unt,*)j,mom!*real(nd,kind=8)
!!$  end do
!!$  call f_close(unt)

  call f_free(y)

end subroutine ISF_family

subroutine get_isf_family(m_max,itype,nmoms,m,nrange,ch,cg)
  use dictionaries, only: f_err_throw
  use yaml_strings, only: yaml_toa
  implicit none
  integer, intent(in) :: itype,m_max,nmoms
  integer, intent(out) :: m,nrange
  double precision, dimension(-m_max:m_max), intent(out) :: ch,cg
  !local variables
  integer :: i,dsflb,dsfrb,sflb,sfrb,unsflb,unsfrb,i1,ifac!,sfl,unsfl
  double precision, dimension(-m_max:m_max) :: cht,chu

  !Only itype=2,8,14,16,20,24,30,40,50,60,100
  select case(itype)
  case(2,4,6,8,14,16,20,24,30,40,50,60,100)
     !O.K.
  case default
     call f_err_throw('"Only interpolating functions 2, 4, 6, 8, 14, 16, 20, 24, 30, 40, 50, 60, 100, used:' &
          & // trim(yaml_toa(itype))//'"', err_name='BIGDFT_RUNTIME_ERROR')
  end select

  !case of the direct function
  ch=0.d0
  cg=0.d0
  cht=0.d0
  m=itype+nmoms+2
  ifac=1
  if (nmoms==0) ifac=2

  call scal_filters(m,itype,ch(-m),ifac)
  nrange=2*(itype+nmoms)

  if (nmoms > 0) then
     !lifted case
     call scal_filters(m_max,nmoms,chu,ifac)
     ! prepare
     chu(0)=0.d0
     chu=-2*chu
     ch(0)=-ch(0)
     ! dualscal/wave - lifted
     dsflb=-itype-nmoms+2
     dsfrb=-dsflb
!  sfl=2*itype-1
!  unsfl=2*nmoms-1
     sflb=-itype+1
     sfrb=itype-1
     unsflb=-nmoms+1
     unsfrb=nmoms-1
!sfrb-(sfl+unsfl-1)+i1-dsflb+1=-nmoms+1
!sflb+i1-1-dsflb+1=nmoms-1+i1
!unsfrb-i1+1+dsflb-1=1-itype-i1=sflb-i1
!unsflb+(sfl+unsfl-1)-i1+dsflb-1=itype-1-i1=sfrb-i1
   ! dualscal/wave - lifted
     !print *,'ch',ch(sflb:sfrb)
     !print *,'chu',chu(unsflb:unsfrb)
     do i1 = dsflb, dsfrb
        if( i1==0) then
           cht(i1) = 1.d0
        else
           cht(i1) = 0.d0
        end if
        cht(i1) = cht(i1) + &
             sum( &
             ch( max( sflb, unsflb+i1):min( sfrb, unsfrb+i1 ) ) * &
             chu( max( unsflb, sflb-i1 ):min( unsfrb, sfrb-i1 )))
     end do
     ch(0) = -ch(0)
     !then we have to invert ch and cht
     !chu is not needed anymore
     chu=cht
     cht=2*ch
     ch=2*chu
!!$     cht=2*cht
!!$     ch=2*ch
  else
     cht( 0)=1.d0
  end if

  ! g coefficients from h coefficients
  do i=-m,m-1
     cg(i+1)=cht(-i)*(-1.d0)**(i+1)
  enddo

!!$  do i=-m,m
!!$     print *,'i',i,ch(i),cg(i)
!!$  end do
!!$  stop

end subroutine get_isf_family

subroutine scal_filters(m,itype,ch,fact)
  implicit none
  integer, intent(in) :: m,itype
  integer, intent(in) :: fact
  double precision, dimension(-m:m), intent(out) :: ch
  !local variables
  integer :: i
  double precision :: pref

  ch=0.d0
  ch(0)=0.5d0*fact
  pref=fact*prefactor(itype/2)
  do i=1,itype/2
     ch(2*i-1)=h_odd(itype/2,i)*pref
     ch(-2*i+1)=ch(2*i-1)
  end do

contains
  !>calculate the filters for the interpolating family of order k
  !! has to be multiplied by prefactor*2
  pure function h_odd(p,k)
    implicit none
    integer, intent(in) :: p
    integer, intent(in) :: k
    double precision :: h_odd

    h_odd=(-1)**(k+1)*dble(p+k)/dble(2*k-1)*binomial(2*p,p-k)

  end function h_odd

  !> calculate the binomial coefficient n over k
  pure function binomial(n,k)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: k
    double precision :: binomial
    !local variables
    integer :: i

    binomial=1.d0
    do i=1,k
       binomial=binomial*dble(n+1-i)/dble(i)
    end do
  end function binomial

  !> depends of the isf family
  pure function prefactor(p)
    implicit none
    integer, intent(in) :: p
    double precision :: prefactor
    !local variables
    integer :: i

    prefactor=1.d0
    do i=1,p
       prefactor=prefactor*dble(2*p+1-i)/dble(16*i)
    end do

  end function prefactor

end subroutine scal_filters

!> generic routine for the wavelet transformation
subroutine back_trans(m,ch,cg,nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: m !<size of the filters
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  !> wavelet filters, low-pass and high-pass
  real(kind=8), dimension(-m:m) ::  ch,cg
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0

     do j=-m/2,m/2-1

        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do

  end do
end subroutine back_trans
