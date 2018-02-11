!> @file
!!  Definition of the Spherical Harmonics, Multipoles and related operations
!! @author
!!    Copyright (C) 2016-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_harmonics
  use numerics
  use f_precisions, only: dp => f_double
  use f_enums
  use f_arrays
  implicit none

  integer, parameter :: S_=0,X_=1,Y_=2,Z_=3
  integer, parameter :: XY_=4,YZ_=5,XZ_=6
  integer, parameter :: X2_=7,Y2_=8,Z2_=9

  type(f_enumerator) :: SOLID_HARMONIC_ENUM=f_enumerator('SOLID_HARMONIC',1,null())

  !> Multipole of a scalar field, can be given in different formats and units
  type, public :: f_multipoles
     type(f_enumerator) :: fmt !< specifies the format of the multipole
     integer :: lmax=-1 !<maximum value to construct the multipole
     integer :: nmonomials=-1
     real(dp), dimension(3) :: rxyz=0.0_dp !< center of the multipole
     !type(f_vector), dimension(:), pointer :: Q=>null() !,data of the multipole
     real(dp), dimension(0:9) :: monomials=0.0_dp !expression in of reductible monomials for which the multipoles can be expressed
  end type f_multipoles

  private

  public :: solid_harmonic,f_multipoles_create,f_multipoles_free
  public :: field_multipoles,get_monopole,get_dipole

  contains

    pure subroutine nullify_f_multipoles(mp)
      use f_utils
      implicit none
      type(f_multipoles), intent(out) :: mp
      call nullify_f_enum(mp%fmt)
      mp%lmax=-1
      mp%nmonomials=-1
      mp%rxyz=0.0_dp
      !nullify(mp%Q)
      mp%monomials=0.0_dp
    end subroutine nullify_f_multipoles

    subroutine f_multipoles_create(mp,lmax,center)
      use dynamic_memory
      use yaml_strings
      implicit none
      integer, intent(in) :: lmax
      type(f_multipoles), intent(out) :: mp
      real(dp), dimension(3), intent(in), optional :: center
      !local variables
      !integer :: l
      mp%fmt=SOLID_HARMONIC_ENUM
      if (present(center)) then
         mp%rxyz=center
      else
         mp%rxyz=0.0_dp
      end if
      mp%lmax=lmax
      mp%nmonomials=9 !(3**(lmax+1)-1)/2-1
      !mp%Q=f_malloc_ptr(0.to.lmax,id='multipoles')
      !do l=0,lmax
      !   mp%Q(l)=f_malloc0_ptr(-l.to.l,id='ql'+yaml_toa(l))
      !end do
    end subroutine f_multipoles_create

    subroutine f_multipoles_free(mp)
      use dynamic_memory
      implicit none
      type(f_multipoles), intent(inout) :: mp
      !call f_array_ptr_free(mp%Q)
      call nullify_f_multipoles(mp)
    end subroutine f_multipoles_free

    !>Calculate a set of multipoles from a given box iterator.
    !! Defined up to quadrupoles
    subroutine field_multipoles(bit,field,nfield,mp)
      use box
      implicit none
      integer, intent(in) :: nfield !<number of components of the field to reduce on (useful for spin)
      type(box_iterator), intent(inout) :: bit
      real(dp), dimension(bit%mesh%ndims(1)*bit%mesh%ndims(2)*(bit%i3e-bit%i3s+1),nfield), intent(in) :: field
      type(f_multipoles), intent(inout) :: mp !<multipole structure

      !local variables
      integer :: ispin
      real(dp) :: q

      !the multipoles have to be created and initialized before.
      do ispin=1,nfield
         do while(box_next_point(bit))
            q= field(bit%ind,ispin)*bit%mesh%volume_element
            bit%tmp=closest_r(bit%mesh,bit%rxyz,mp%rxyz)
            call accumulate_multipoles(bit%tmp,q,mp%nmonomials,mp%monomials)
         end do
      end do

    end subroutine field_multipoles   

    subroutine vector_multipoles(nat,rxyz,mp,origin,charges,lookup)
      implicit none
      integer, intent(in) :: nat
      real(dp), dimension(3,nat), intent(in) :: rxyz
      type(f_multipoles), intent(inout) :: mp
      real(dp), dimension(3), intent(in), optional :: origin
      real(dp), dimension(*), intent(in), optional :: charges
      integer, dimension(nat), intent(in), optional :: lookup
      !local variables
      integer :: iat,jat
      real(dp), dimension(3) :: oxyz,tmp
      oxyz=0.0_dp
      if (present(origin)) oxyz=origin
      do iat=1,nat
         tmp=rxyz(:,iat)-oxyz
         jat=iat
         if (present(lookup)) jat=lookup(iat)
         call accumulate_multipoles(tmp,charges(jat),&
              mp%nmonomials,mp%monomials)
      end do
    end subroutine vector_multipoles

    pure subroutine accumulate_multipoles(rxyz,density,n,monomials)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: density
      real(dp), dimension(3), intent(in) :: rxyz
      real(dp), dimension(0:n), intent(inout) :: monomials

      monomials(S_)=monomials(S_)+density

      monomials(X_)=monomials(X_)+density*rxyz(X_)
      monomials(Y_)=monomials(Y_)+density*rxyz(Y_)
      monomials(Z_)=monomials(Z_)+density*rxyz(Z_)

      monomials(XY_)=monomials(XY_)+density*rxyz(X_)*rxyz(Y_)
      monomials(XZ_)=monomials(XZ_)+density*rxyz(X_)*rxyz(Z_)
      monomials(YZ_)=monomials(YZ_)+density*rxyz(Y_)*rxyz(Z_)

      monomials(X2_)=monomials(X2_)+density*rxyz(X_)*rxyz(X_)
      monomials(Y2_)=monomials(Y2_)+density*rxyz(Y_)*rxyz(Y_)
      monomials(Z2_)=monomials(Z2_)+density*rxyz(Z_)*rxyz(Z_)

    end subroutine accumulate_multipoles

!!$    pure subroutine f_multipoles_accumulate(Q,lmax,rxyz,density)
!!$      implicit none
!!$      integer, intent(in) :: lmax
!!$      real(dp), intent(in) :: density
!!$      real(dp), dimension(3), intent(in) :: rxyz
!!$      type(f_vector), dimension(0:lmax), intent(inout) :: Q
!!$      !local variables
!!$      integer :: l,m
!!$      real(dp) :: tt,factor
!!$
!!$      Q%monomials(0)=density !monopole does not need pows
!!$      do ipow=1,nmonomials
!!$         tt=rxyz_monomials(rxyz)
!!$         Q%monomials(ipow)=tt*density
!!$      end do
!!$      !do l=0,lmax
!!$      !   factor=sqrt(fourpi/real(2*l+1,dp))
!!$      !   do m=-l,l
!!$      !      tt = solid_harmonic(0, l, m,rxyz(1),rxyz(2),rxyz(3))
!!$      !      tt = tt*factor
!!$      !      Q(l)%ptr(m)=Q(l)%ptr(m)+tt*density
!!$      !   end do
!!$      !end do
!!$    end subroutine f_multipoles_accumulate

    pure function get_monopole(mp) result(q)
      implicit none
      type(f_multipoles), intent(in) :: mp
      real(dp) :: q
      q=mp%monomials(0)
      !q=mp%Q(0)%ptr(0)
    end function get_monopole

    pure function get_dipole(mp) result(d)
      implicit none
      type(f_multipoles), intent(in) :: mp
      real(dp), dimension(3) :: d
      !d=mp%Q(1)%ptr
      d=mp%monomials(1:3)
    end function get_dipole


    !> Calculates the solid harmonic S_lm (possibly multplied by a power or r) for given values of l, m, x, y, z.
    !! They are normalized such that the integral over the angle gives r^2, i.e.
    !! \int d\Omega S_{lm}*S_{l'm'}/r^{2l} = r^2 \delta_{ll'}\delta_{mm'}
    !! r_exponent indicates how the function is multiplied by r: The final result is given by S_lm*r^(r_exponent*l), with the
    !! definition of the S_lm given above.
    !! rmin gives the minimal radius that is used for the multiplication by r^(r_exponent*l) (can be used to avoid the divergence
    !! around r=0)
    pure function solid_harmonic(r_exponent, l, m, x, y, z) result(sh)

      implicit none
      ! Calling arguments
      integer,intent(in) :: r_exponent
      integer,intent(in) :: l, m
      real(dp),intent(in) ::  x, y, z !<given in cartesian, orthorhombic form
      real(dp) :: sh

      ! Local variables
      integer,parameter :: l_max=2
      real(dp) :: r, r2

!!$      if (l<0) call f_err_throw('l must be non-negative',err_name='BIGDFT_RUNTIME_ERROR')
!!$      if (l>l_max) call f_err_throw('solid harmonics only implemented up to l='//trim(yaml_toa(l_max)),&
!!$           err_name='BIGDFT_RUNTIME_ERROR')
!!$      if (abs(m)>l) call f_err_throw('abs of m ('//trim(yaml_toa(m))//') must not be larger than l ('//trim(yaml_toa(l))//')', &
!!$           err_name='BIGDFT_RUNTIME_ERROR')

      sh=0.0_dp
      select case (l)
      case (0)
         ! No need for r, as l=0
         sh = sqrt(oneofourpi)
      case (1)
         r2 = x**2+y**2+z**2
         r = sqrt(r2)
         select case (m)
         case (-1)
            sh = sqrt(3.0_dp/(fourpi))*y
         case (0)
            sh = sqrt(3.0_dp/(fourpi))*z
         case (1)
            sh = sqrt(3.0_dp/(fourpi))*x
         end select
         ! Multiply by r^{r_exp*l}
         sh = sh*r**r_exponent
      case (2)
         r2 = x**2+y**2+z**2
         select case (m)
         case (-2)
            sh = sqrt(15.d0/(4.d0*pi))*x*y
         case (-1)
            sh = sqrt(15.d0/(4.d0*pi))*y*z
         case (0)
            sh = sqrt(5.d0/(16.d0*pi))*(-x**2-y**2+2.d0*z**2)
         case (1)
            sh = sqrt(15.d0/(4.d0*pi))*z*x
         case (2)
            sh = sqrt(15.d0/(16.d0*pi))*(x**2-y**2)
         end select
         ! Multiply by r^{r_exp*l}
         sh = sh*r2**r_exponent
      end select

    end function solid_harmonic

end module f_harmonics
