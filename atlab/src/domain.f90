!> @file
!!    Modulefile for handling fundamental data structed and methods of the simulation domain
!!
!! @author
!!    G. Fisicaro, L. Genovese (June 2018)
!!    Copyright (C) 2018-2018 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module at_domain
  use f_precisions, gp=>f_double
  use numerics, only: onehalf,pi,Radian_Degree
  use f_enums
  use dictionaries
  use yaml_strings
  implicit none
  private

  !>parameter for the definition of the bc
  integer, parameter :: NULL_PAR=-100
  integer, parameter :: INT_BC_IS_PERIODIC=-1000
  integer, parameter :: FREE=0
  integer, parameter :: PERIODIC=1
  integer, parameter :: AB_=3,BC_=1,AC_=2,A_=1,B_=2,C_=3,X_=1,Y_=2,Z_=3
  integer, parameter :: BOHR=0
  integer, parameter :: ANGSTROEM=1
  integer, parameter :: NANOMETER=3

  character(len=*), parameter, public :: DOMAIN_UNITS = 'units'
  character(len=*), parameter, public :: DOMAIN_CELL = 'cell'
  character(len=*), parameter, public :: DOMAIN_ABC = 'abc'
  character(len=*), parameter, public :: DOMAIN_ALPHA = 'alpha'
  character(len=*), parameter, public :: DOMAIN_BETA = 'beta'
  character(len=*), parameter, public :: DOMAIN_GAMMA = 'gamma'


  type(f_enumerator), parameter,  public :: FREE_BC=f_enumerator('FREE_BC',FREE,null())
  type(f_enumerator), parameter,  public :: PERIODIC_BC=f_enumerator('PERIODIC_BC',PERIODIC,null())
  type(f_enumerator), parameter,  public :: ATOMIC_UNITS=f_enumerator('atomic',BOHR,null())
  type(f_enumerator), parameter,  public :: ANGSTROEM_UNITS=f_enumerator('angstroem',ANGSTROEM,null())
  type(f_enumerator), parameter,  public :: NANOMETER_UNITS=f_enumerator('nanometer',NANOMETER,null())

  !as it can be used in enum_attr calls
  type(f_enumerator), public :: DOUBLE_PRECISION_ENUM=f_enumerator('double precision',487,null())

!!$  !>data type for the simulation domain
!!$  type, public :: domain
!!$     integer :: units=NULL_PAR
!!$     logical :: orthorhombic=.true. !<true if the cell is orthorhombic
!!$     integer, dimension(3) :: bc=NULL_PAR !< boundary conditions on each direction
!!$     real(gp), dimension(3) :: angrad=0.0_gp !<angles between the dimensions in radians (alpha_bc,beta_ac,gamma_bc)
!!$     real(gp), dimension(3,3) :: abc=0.0_gp !<matrix of the normalized translation vectors direction
!!$     real(gp), dimension(3,3) :: uabc=0.0_gp !<matrix of the normalized translation vectors direction
!!$     real(gp), dimension(3) :: acell=0.0_gp !lengths of the primitive vectors abc
!!$     real(gp), dimension(3,3) :: gd=0.0_gp !<covariant metric needed for non-orthorhombic operations
!!$     real(gp), dimension(3,3) :: gu=0.0_gp !<controvariant metric needed for non-orthorhombic operations
!!$     real(gp) :: detgd=0.0_gp !<determinant of the covariant matrix
!!$  end type domain

  !>data type for the simulation domain
  type, public :: domain
     integer :: units
     logical :: orthorhombic !<true if the cell is orthorhombic
     integer, dimension(3) :: bc !< boundary conditions on each direction
     real(gp), dimension(3) :: angrad !<angles between the dimensions in radians (alpha_bc,beta_ac,gamma_bc)
     real(gp), dimension(3,3) :: abc !<matrix of the normalized translation vectors direction
     real(gp), dimension(3,3) :: uabc !<matrix of the normalized translation vectors direction
     real(gp), dimension(3) :: acell !lengths of the primitive vectors abc
     real(gp), dimension(3,3) :: gd !<covariant metric needed for non-orthorhombic operations
     real(gp), dimension(3,3) :: gu !<controvariant metric needed for non-orthorhombic operations
     real(gp) :: detgd !<determinant of the covariant matrix
  end type domain

  interface dotp_gu
     module procedure dotp_gu,dotp_gu_add1,dotp_gu_add2
  end interface dotp_gu

  interface square_gu
     module procedure square,square_add
  end interface square_gu

  interface dotp_gd
     module procedure dotp_gd,dotp_gd_add1,dotp_gd_add2,dotp_gd_add12
  end interface dotp_gd

  interface square_gd
     module procedure square_gd,square_gd_add
  end interface square_gd

  public :: domain_null,domain_set_from_dict,units_enum_from_str,domain_new
  public :: rxyz_ortho,distance,closest_r
  public :: dotp_gu,dotp_gd,square_gu,square_gd
  public :: domain_geocode,bc_periodic_dims,geocode_to_bc,domain_periodic_dims
  public :: geocode_to_bc_enum,change_domain_BC

contains

!!$  function domain_null()
!!$    type(domain) :: domain_null
!!$  end function domain_null

  pure function domain_null() result(dom)
    implicit none
    type(domain) :: dom
    call nullify_domain(dom)
  end function domain_null

  pure subroutine nullify_domain(dom)
    implicit none
    type(domain), intent(out) :: dom
    dom%units=NULL_PAR
    dom%orthorhombic=.true.
    dom%bc=NULL_PAR
    dom%angrad=0.0_gp
    dom%abc=0.0_gp
    dom%uabc=0.0_gp
    dom%acell=0.0_gp
    dom%gd=0.0_gp
    dom%gu=0.0_gp
    dom%detgd=0.0_gp
  end subroutine nullify_domain

  pure elemental function bc_is_periodic(bc) result(yes)
    implicit none
    type(f_enumerator), intent(in) :: bc
    logical :: yes
    yes = bc==PERIODIC_BC
  end function bc_is_periodic

  function domain_new(units,bc,abc,alpha_bc,beta_ac,gamma_ab,acell) result(dom)
    use f_utils
    implicit none
    type(f_enumerator), intent(in) :: units
    type(f_enumerator), dimension(3), intent(in) :: bc !boundary conditions of the system
    !length of the domain
    real(gp), dimension(3), intent(in), optional :: acell
    !> arrays of the unit vectors of the cell. Normalized, in fortran order a_i=abc(i,1), b_i=abc(i,2)
    real(gp), dimension(3,3), intent(in), optional :: abc
    real(gp), intent(in), optional :: alpha_bc,beta_ac,gamma_ab
    type(domain) :: dom
    !local variables
    integer :: i,j

    dom%units=toi(units)
    dom%bc=toi(bc)

    !default orthorhombic
    dom%orthorhombic=.true.
    dom%angrad=onehalf*pi
    !Set the covariant metric
    dom%gd(1,1) = 1.0_gp
    dom%gd(1,2) = 0.0_gp
    dom%gd(1,3) = 0.0_gp
    dom%gd(2,2) = 1.0_gp
    dom%gd(2,3) = 0.0_gp
    dom%gd(3,3) = 1.0_gp
    dom%gd(2,1) = dom%gd(1,2)
    dom%gd(3,1) = dom%gd(1,3)
    dom%gd(3,2) = dom%gd(2,3)

    dom%detgd = 1.0_gp
    !Set the contravariant metric
    dom%gu(1,1) = 1.0_gp
    dom%gu(1,2) = 0.0_gp
    dom%gu(1,3) = 0.0_gp
    dom%gu(2,2) = 1.0_gp
    dom%gu(2,3) = 0.0_gp
    dom%gu(3,3) = 1.0_gp
    dom%gu(2,1) = dom%gu(1,2)
    dom%gu(3,1) = dom%gu(1,3)
    dom%gu(3,2) = dom%gu(2,3)

    dom%uabc=dom%gd
    dom%abc=0.0_gp
    dom%acell=0.0_gp

    !!this means free BC
    !if (.not. any(bc_is_periodic(bc))) return
    !this means only Periodic BC go ahead
    if ( .not. any(bc_is_periodic(bc))) return

    !check of the availablilty of the options
    if (present(acell) .eqv. present(abc)) call f_err_throw('Domain inconsistency: only acell or abc may be given')

    if (any([present(alpha_bc),present(beta_ac),present(gamma_ab)]) .and. present(abc)) &
         call f_err_throw('Domain inconsistency: angles cannot be present if unit cell vectors are given')

    if (present(abc)) then
       dom%abc=abc
       call abc_to_angrad_and_acell(abc,dom%angrad,dom%acell)
    else
       dom%acell=acell
       if (present(alpha_bc)) dom%angrad(BC_)=alpha_bc
       if (present(beta_ac)) dom%angrad(AC_)=beta_ac
       if (present(gamma_ab)) dom%angrad(AB_)=gamma_ab
       call angrad_and_acell_to_abc(dom%angrad,dom%acell,dom%abc)
    end if
    do i=1,3
       if (dom%acell(i) /= 0.0_gp) dom%uabc(:,i)=dom%abc(:,i)/dom%acell(i)
    end do

    call f_assert(all(dom%angrad > 0.0_gp),'Domain inconsistency: some of the angles are not positive')
    if (dom%angrad(AB_) /= onehalf*pi .and. .not.(bc_is_periodic(bc(A_)) .and. bc_is_periodic(bc(B_)))) &
         call f_err_throw('Domain inconsistency: cannot provide angle gamma_ab if the dimensions a and b are not periodic')
    if (dom%angrad(BC_) /= onehalf*pi .and. .not.(bc_is_periodic(bc(B_)) .and. bc_is_periodic(bc(C_)))) &
         call f_err_throw('Domain inconsistency: cannot provide angle alpha_bc if the dimensions b and c are not periodic')
    if (dom%angrad(AC_) /= onehalf*pi .and. .not.(bc_is_periodic(bc(A_)) .and. bc_is_periodic(bc(C_)))) &
         call f_err_throw('Domain inconsistency: cannot provide angle beta_ac if the dimensions a and c are not periodic')

    call f_assert(get_abc_angrad_acell_consistency(dom%abc,dom%angrad,dom%acell)< 1.e-10,&
          'Domain inconsisntency: bad definition of angles angrad or acell wrt abc')

    dom%orthorhombic=all(dom%angrad==onehalf*pi)

    if ((domain_geocode(dom) == 'F' .or. domain_geocode(dom) == 'W') .and. (.not. dom%orthorhombic)) &
         call f_err_throw('For geocode="F","W" the cell must be orthorhombic')

    if (.not. dom%orthorhombic) then
       !some consistency check on the angles should be performed
       !1) sum(angrad) < twopi
       if (all(dom%angrad==dom%angrad(1))) then
          !Treat the case of equal angles (except all right angles) :
          !generates trigonal symmetry wrt third axis
          !Set the covariant metric
          dom%gd(1,1) = 1.0_gp
          dom%gd(1,2) = cos(dom%angrad(3)) !gamma_ab
          dom%gd(1,3) = cos(dom%angrad(2)) !beta_ac
          dom%gd(2,2) = 1.0_gp
          dom%gd(2,3) = cos(dom%angrad(1)) !alpha_bc
          dom%gd(3,3) = 1.0_gp
          !Set the determinant of the covariant metric
          dom%detgd = 1.0_gp - cos(dom%angrad(1))**2 - cos(dom%angrad(2))**2 - cos(dom%angrad(3))**2 +&
               2.0_gp*cos(dom%angrad(1))*cos(dom%angrad(2))*cos(dom%angrad(3))
          !Set the contravariant metric
          dom%gu(1,1) = (sin(dom%angrad(1))**2)/dom%detgd
          dom%gu(1,2) = (cos(dom%angrad(2))*cos(dom%angrad(1))-cos(dom%angrad(3)))/dom%detgd
          dom%gu(1,3) = (cos(dom%angrad(3))*cos(dom%angrad(1))-cos(dom%angrad(2)))/dom%detgd
          dom%gu(2,2) = (sin(dom%angrad(2))**2)/dom%detgd
          dom%gu(2,3) = (cos(dom%angrad(3))*cos(dom%angrad(2))-cos(dom%angrad(1)))/dom%detgd
          dom%gu(3,3) = (sin(dom%angrad(3))**2)/dom%detgd
       else if (domain_geocode(dom) == 'P') then
          !Set the covariant metric
          dom%gd(1,1) = 1.0_gp
          dom%gd(1,2) = cos(dom%angrad(3)) !gamma_ab
          dom%gd(1,3) = cos(dom%angrad(2)) !beta_ac
          dom%gd(2,2) = 1.0_gp
          dom%gd(2,3) = cos(dom%angrad(1)) !alpha_bc
          dom%gd(3,3) = 1.0_gp
          !Set the determinant of the covariant metric
          dom%detgd = 1.0_gp - cos(dom%angrad(1))**2 - cos(dom%angrad(2))**2 - cos(dom%angrad(3))**2 +&
               2.0_gp*cos(dom%angrad(1))*cos(dom%angrad(2))*cos(dom%angrad(3))
          !Set the contravariant metric
          dom%gu(1,1) = (sin(dom%angrad(1))**2)/dom%detgd
          dom%gu(1,2) = (cos(dom%angrad(2))*cos(dom%angrad(1))-cos(dom%angrad(3)))/dom%detgd
          dom%gu(1,3) = (cos(dom%angrad(3))*cos(dom%angrad(1))-cos(dom%angrad(2)))/dom%detgd
          dom%gu(2,2) = (sin(dom%angrad(2))**2)/dom%detgd
          dom%gu(2,3) = (cos(dom%angrad(3))*cos(dom%angrad(2))-cos(dom%angrad(1)))/dom%detgd
          dom%gu(3,3) = (sin(dom%angrad(3))**2)/dom%detgd
       else !only Surfaces is possible here
          !Set the covariant metric
          dom%gd=0.0_gp
          dom%gd(1,1) = 1.0_gp
          dom%gd(1,3) = cos(dom%angrad(2)) !beta_ac
          dom%gd(2,2) = 1.0_gp
          dom%gd(3,3) = 1.0_gp
          !Set the determinant of the covariant metric
          dom%detgd = sin(dom%angrad(2))**2
          !Set the contravariant metric
          dom%gu=0.0_gp
          dom%gu(1,1) = 1.0_gp/dom%detgd
          dom%gu(1,3) = -cos(dom%angrad(2))/dom%detgd
          dom%gu(2,2) = 1.0_gp!/mesh%detgd
          dom%gu(3,3) = 1.0_gp/dom%detgd
       end if
    else
       dom%angrad=onehalf*pi
       !Set the covariant metric
       dom%gd(1,1) = 1.0_gp
       dom%gd(1,2) = 0.0_gp
       dom%gd(1,3) = 0.0_gp
       dom%gd(2,2) = 1.0_gp
       dom%gd(2,3) = 0.0_gp
       dom%gd(3,3) = 1.0_gp
       dom%detgd = 1.0_gp
       !Set the contravariant metric
       dom%gu(1,1) = 1.0_gp
       dom%gu(1,2) = 0.0_gp
       dom%gu(1,3) = 0.0_gp
       dom%gu(2,2) = 1.0_gp
       dom%gu(2,3) = 0.0_gp
       dom%gu(3,3) = 1.0_gp
    end if
    dom%gd(2,1) = dom%gd(1,2)
    dom%gd(3,1) = dom%gd(1,3)
    dom%gd(3,2) = dom%gd(2,3)

    dom%gu(2,1) = dom%gu(1,2)
    dom%gu(3,1) = dom%gu(1,3)
    dom%gu(3,2) = dom%gu(2,3)
    do i=1,3
       do j=1,3
          if (abs(dom%gd(i,j)).lt.1.0d-15) dom%gd(i,j)=0.0_gp
          if (abs(dom%gu(i,j)).lt.1.0d-15) dom%gu(i,j)=0.0_gp
       end do
    end do

    !here we should verify that the the inverse metric times the metric is the identity

  end function domain_new

  function change_domain_BC(dom_in,geocode) result(dom)
    implicit none
    type(domain), intent(in) :: dom_in
    character(len=1), intent(in) :: geocode
    type(domain) :: dom

    dom=dom_in

    select case(geocode)
    case('P')
        dom%bc=[PERIODIC,PERIODIC,PERIODIC]
    case('S')
        dom%bc=[PERIODIC,FREE,PERIODIC]
    case('W')
        dom%bc=[FREE,FREE,PERIODIC]
    case('F')
        dom%bc=[FREE,FREE,FREE]
    case default
        dom%bc=[FREE,FREE,FREE]
    end select

    where (dom%bc==FREE) dom%acell=0.0_gp

  end function change_domain_BC

  subroutine abc_to_angrad_and_acell(abc,angrad,acell)
    implicit none
    real(gp), dimension(3,3), intent(in) :: abc
    real(gp), dimension(3), intent(out) :: angrad,acell

    !construct reduced cell vectors (if they do not exist)
    acell(A_) = sqrt(abc(X_,A_)**2 + abc(Y_,A_)**2 + abc(Z_,A_)**2)
    acell(B_) = sqrt(abc(X_,B_)**2 + abc(Y_,B_)**2 + abc(Z_,B_)**2)
    acell(C_) = sqrt(abc(X_,C_)**2 + abc(Y_,C_)**2 + abc(Z_,C_)**2)

    angrad(AB_) = acos(dot_product(abc(:,A_),abc(:,B_))/(acell(A_)*acell(B_)))
    angrad(BC_) = acos(dot_product(abc(:,B_),abc(:,C_))/(acell(B_)*acell(C_)))
    angrad(AC_) = acos(dot_product(abc(:,A_),abc(:,C_))/(acell(A_)*acell(C_)))

  end subroutine abc_to_angrad_and_acell

  subroutine angrad_and_acell_to_abc(angrad,acell,abc)
    implicit none
    real(gp), dimension(3), intent(in) :: angrad,acell
    real(gp), dimension(3,3), intent(out) :: abc
    !local variables
    integer :: i
    real(gp) :: aa,cc,a2,cosang


    if (all(angrad== onehalf*pi)) then
       ! orthorhombic case
       abc=0.0_gp
       do i=1,3
          abc(i,i)=1.0_gp
       end do
    else if(all(angrad == angrad(BC_))) then
       !Treat the case of equal angles (except all right angles) :
       !generates trigonal symmetry wrt third axis
       cosang=cos(angrad(BC_))
       a2=2.0_gp/3.0_gp*(1.0_gp-cosang)
       aa=sqrt(a2)
       cc=sqrt(1.0_gp-a2)

       abc(X_,A_)=aa
       abc(Y_,A_)=0.0_gp
       abc(Z_,A_)=cc

       abc(X_,B_)=-0.5_gp*aa
       abc(Y_,B_)=sqrt(3.0_gp)*0.5_gp*aa
       abc(Z_,B_)=cc

       abc(X_,C_)=-0.5_gp*aa
       abc(Y_,C_)=-sqrt(3.0_gp)*0.5_gp*aa
       abc(Z_,C_)=cc
    else if (angrad(AC_)/= onehalf*pi .and. angrad(AB_) == onehalf*pi .and. angrad(BC_)==onehalf*pi ) then
       abc=0.0_gp
       abc(X_,A_)=1.0_gp
       abc(Y_,B_)=1.0_gp

       abc(X_,C_)=cos(angrad(AC_))
       !abc(Z_,C_)=sin(angrad(AC_))
       abc(Z_,C_)=sqrt(1.0_gp-abc(X_,C_)**2)
    else if (angrad(BC_)/= onehalf*pi .and. angrad(AB_) == onehalf*pi .and. angrad(AC_)==onehalf*pi ) then
       abc=0.0_gp
       abc(X_,A_)=1.0_gp
       abc(Y_,B_)=1.0_gp
       abc(Y_,C_)=cos(angrad(BC_))
       !abc(Z_,C_)=sin(angrad(BC_))
       abc(Z_,C_)=sqrt(1.0_gp-abc(Y_,C_)**2)
    else if (angrad(AB_)/= onehalf*pi .and. angrad(AC_) == onehalf*pi .and. angrad(BC_)==onehalf*pi ) then
       abc=0.0_gp
       abc(X_,A_)=1.0_gp
       abc(X_,B_)=cos(angrad(AB_))
       abc(Y_,B_)=sin(angrad(AB_))
       abc(Z_,C_)=1.0_gp
    else if (angrad(AC_)/= onehalf*pi .and. angrad(AB_) /= onehalf*pi .and. angrad(BC_)==onehalf*pi ) then
       abc=0.0_gp
       abc(X_,A_)=1.0_gp
       abc(X_,B_)=cos(angrad(AB_))
       abc(Y_,B_)=sin(angrad(AB_))
       abc(X_,C_)=cos(angrad(AC_))
       abc(Y_,C_)=(-abc(X_,B_)*abc(X_,C_))/abc(Y_,B_)
       abc(Z_,C_)=sqrt(1.0_gp-abc(X_,C_)**2-abc(Y_,C_)**2)
    else if (angrad(AC_)/= onehalf*pi .and. angrad(AB_) == onehalf*pi .and. angrad(BC_)/=onehalf*pi ) then
       abc=0.0_gp
       abc(X_,A_)=1.0_gp
       abc(Y_,B_)=1.0_gp
       abc(X_,C_)=cos(angrad(AC_))
       abc(Y_,C_)=cos(angrad(BC_))
       abc(Z_,C_)=sqrt(1.0_gp-abc(X_,C_)**2-abc(Y_,C_)**2)
    else if (angrad(AC_)== onehalf*pi .and. angrad(AB_) /= onehalf*pi .and. angrad(BC_)/=onehalf*pi ) then
       abc=0.0_gp
       abc(X_,A_)=1.0_gp
       abc(X_,B_)=cos(angrad(AB_))
       abc(Y_,B_)=sin(angrad(AB_))
       abc(Y_,C_)=cos(angrad(BC_))/abc(Y_,B_)
       abc(Z_,C_)=sqrt(1.0_gp-abc(Y_,C_)**2)
    else if (all((angrad) /= onehalf*pi)) then
       abc=0.0_gp
       abc(X_,A_)=1.0_gp
       abc(X_,B_)=cos(angrad(AB_))
       abc(Y_,B_)=sin(angrad(AB_))
       abc(X_,C_)=cos(angrad(AC_))
       abc(Y_,C_)=(cos(angrad(BC_))-abc(X_,B_)*abc(X_,C_))/abc(Y_,B_)
       abc(Z_,C_)=sqrt(1.0_gp-abc(X_,C_)**2-abc(Y_,C_)**2)
    end if
    do i=1,3
       abc(:,i)=acell(i)*abc(:,i)
    end do

  end subroutine angrad_and_acell_to_abc

  function get_abc_angrad_acell_consistency(abc,angrad,acell) result(tol)
    implicit none
    real(gp), dimension(3,3), intent(in) :: abc
    real(gp), dimension(3), intent(in) :: angrad,acell
    real(gp) :: tol
    !local variables
    integer, parameter :: ANGLE_=1,AXIS1_=2,AXIS2_=3
    integer :: i,angle,ax1,ax2
    real(gp), dimension(3) :: vect_norm,cang
    integer, dimension(3,3) :: icontrol

    tol=0.0_gp
    do i=1,3
       vect_norm(i) = sqrt(abc(1,i)**2 + abc(2,i)**2 + abc(3,i)**2)
       if (vect_norm(i) > 0.0_gp) tol=max(tol,abs(vect_norm(i) - acell(i)))
    end do

    icontrol(ANGLE_,1)=AB_
    icontrol(AXIS1_,1)=A_
    icontrol(AXIS2_,1)=B_

    icontrol(ANGLE_,2)=BC_
    icontrol(AXIS1_,2)=B_
    icontrol(AXIS2_,2)=C_

    icontrol(ANGLE_,3)=AC_
    icontrol(AXIS1_,3)=A_
    icontrol(AXIS2_,3)=C_

    do i=1,3
       angle=icontrol(ANGLE_,i)
       ax1=icontrol(AXIS1_,i)
       ax2=icontrol(AXIS2_,i)
       if (vect_norm(ax1)*vect_norm(ax2) == 0.0_gp) cycle
       cang(angle) = &
            dot_product(abc(:,ax1),abc(:,ax2))/(vect_norm(ax1)*vect_norm(ax2))
       tol=max(tol,abs(cang(angle) - cos(angrad(angle))))
    end do

  end function get_abc_angrad_acell_consistency

  subroutine bc_and_cell_from_dict(dict,bc,cell)
    implicit none
    type(dictionary), pointer :: dict
    type(f_enumerator), dimension(3), intent(out) :: bc !boundary conditions of the system
    real(gp), dimension(3), intent(out) :: cell
    !local variables
    character(len = max_field_length), dimension(3) :: str
    integer :: i

    !get the specific cells
    bc=PERIODIC_BC
    do i=0,2
       str(i+1)=dict//i
    end do

    where (str == YAML_INFINITY) bc=FREE_BC

    !for the periodic bc retrieve the cell size
    cell=0.0_gp
    do i=0,2
       if (bc_is_periodic(bc(i+1))) cell(i+1)=dict//i
    end do

  end subroutine bc_and_cell_from_dict

  subroutine bc_and_cell_to_dict(dict,bc,cell)
    implicit none
    type(dictionary), pointer :: dict
    type(f_enumerator), dimension(3), intent(in) :: bc !boundary conditions of the system
    real(gp), dimension(3), intent(in) :: cell
    !local variables
    integer :: i

    do i=1,3
       if (bc_is_periodic(bc(i))) then
          call set(dict//(i-1),cell(i))
       else
          call set(dict//(i-1),YAML_INFINITY)
       end if
    end do

  end subroutine bc_and_cell_to_dict

  !that function cannot be pure
  function bc_enum_from_int(ibc) result(bc)
    implicit none
    integer, intent(in) :: ibc
    type(f_enumerator) :: bc

    select case(ibc)
    case(PERIODIC)
       bc=PERIODIC_BC
    case(FREE)
       bc=FREE_BC
    case default
       bc=FREE_BC
    end select
  end function bc_enum_from_int

  function bc_enums_from_int(ibc) result(bc)
    implicit none
    integer, dimension(3), intent(in) :: ibc
    type(f_enumerator), dimension(3) :: bc
    !local variables
    integer :: idim

    do idim=1,3
       select case(ibc(idim))
       case(PERIODIC)
          bc(idim)=PERIODIC_BC
       case(FREE)
          bc(idim)=FREE_BC
       case default
          bc(idim)=FREE_BC
       end select
    end do

  end function bc_enums_from_int

  function units_enum_from_str(str) result(units)
    implicit none
    character(len=*), intent(in) :: str
    type(f_enumerator) :: units

    Units_case: select case(trim(str))
    case('angstroem','angstroemd0')
       units=ANGSTROEM_UNITS
    case('atomic','atomicd0','bohr','bohrd0')
       units=ATOMIC_UNITS
    case('nanometer')
       units=NANOMETER_UNITS
    case default
       units=ANGSTROEM_UNITS
    end select Units_case

    !double precision attribute, to be used from files for backward compatibility
    select case(trim(str))
    case('angstroemd0','atomicd0','bohrd0')
       call f_enum_attr(dest=units,attr=DOUBLE_PRECISION_ENUM)
    end select

  end function units_enum_from_str

  function units_enum_from_int(iunit) result(units)
    implicit none
    integer, intent(in) :: iunit
    type(f_enumerator) :: units

    Units_case: select case(iunit)
    case(ANGSTROEM)
       units=ANGSTROEM_UNITS
    case(BOHR)
       units=ATOMIC_UNITS
    case(NANOMETER)
       units=NANOMETER_UNITS
    case default
       units=ANGSTROEM_UNITS
    end select Units_case

  end function units_enum_from_int

  subroutine domain_set_from_dict(dict,dom)
    use yaml_output
    use f_utils
    implicit none
    type(dictionary), pointer :: dict
    type(domain), intent(out) :: dom
    !local variables
    integer :: i
    real(gp), parameter :: ths2 = 1.e-13_gp
    real(gp) :: alpha,beta,gamma,tol
    real(gp), dimension(3) :: cell,angrad
    real(gp), dimension(3,3) :: abc
    character(len = max_field_length) :: str
    type(f_enumerator) :: units
    type(f_enumerator), dimension(3) :: bc

    ! The units
    call f_strcpy(src='bohr',dest=str)
    str=dict .get. DOMAIN_UNITS
    units=units_enum_from_str(str)

    ! The cell and the boundary conditions
    bc=FREE_BC
    cell=0.0_gp
    if (DOMAIN_CELL .in. dict) call bc_and_cell_from_dict(dict//DOMAIN_CELL,bc,cell)

    !take the default values if the key does not exists in the dictionary
    alpha=90.0_gp
    beta=90.0_gp
    gamma=90.0_gp
    alpha= dict .get. DOMAIN_ALPHA
    beta= dict .get. DOMAIN_BETA
    gamma= dict .get. DOMAIN_GAMMA

    do i=1,3
       if (abs(alpha-90.0_gp).lt.ths2) alpha = 90.0_gp
       if (abs(beta-90.0_gp).lt.ths2) beta = 90.0_gp
       if (abs(gamma-90.0_gp).lt.ths2) gamma = 90.0_gp
    end do

    angrad(BC_)=alpha/Radian_Degree
    angrad(AC_)=beta/Radian_Degree
    angrad(AB_)=gamma/Radian_Degree

    !if (all(dom%bc==PERIODIC_BC) .and. (DOMAIN_ABC .in. dict)) then
    if (DOMAIN_ABC .in. dict) then
       bc=PERIODIC_BC
       call f_zero(abc)
       abc = dict//DOMAIN_ABC
       if (DOMAIN_CELL .in. dict) then
          tol=get_abc_angrad_acell_consistency(abc,angrad,cell)
          if (tol > 1.e-6) call yaml_warning('Inconsitency of tol="'//trim(yaml_toa(tol,fmt='1pe12.5'))//&
               'for the provided domain dictionary, assuming abc is correct')
       end if
       dom=domain_new(units,bc,abc=abc)
    else
       !in this case abc is not provided
       dom=domain_new(units,bc,acell=cell,alpha_bc=angrad(BC_),beta_ac=angrad(AC_),gamma_ab=angrad(AB_))
    end if

  end subroutine domain_set_from_dict

  subroutine domain_merge_to_dict(dict,dom)
    implicit none
    type(dictionary), pointer :: dict
    type(domain), intent(in) :: dom

    call set(dict//DOMAIN_UNITS,toa(units_enum_from_int(dom%units)))

    if (all(dom%bc == PERIODIC)) then
       call set(dict//DOMAIN_ABC,dom%abc)
    else if (any(dom%bc == PERIODIC)) then
       call bc_and_cell_to_dict(dict//DOMAIN_CELL,bc_enums_from_int(dom%bc),dom%acell)
       if (dom%angrad(BC_) /= onehalf*pi) &
            call set(dict//DOMAIN_ALPHA,dom%angrad(BC_))
       if (dom%angrad(AC_) /= onehalf*pi) &
            call set(dict//DOMAIN_BETA,dom%angrad(AC_))
       if (dom%angrad(AB_) /= onehalf*pi) &
            call set(dict//DOMAIN_GAMMA,dom%angrad(AB_))
    end if

  end subroutine domain_merge_to_dict

  ! determine the array of BC from a geometry code
  function geocode_to_bc_enum(geocode) result(bc)
    use dictionaries, only: f_err_throw
    implicit none
    character(len=1), intent(in) :: geocode
    type(f_enumerator), dimension(3) :: bc

    select case(geocode)
    case('P')
       bc=PERIODIC_BC
    case('S')
       bc=PERIODIC_BC
       bc(2)=FREE_BC
    case('F')
       bc=FREE_BC
    case('W')
       bc=FREE_BC
       bc(3)=PERIODIC_BC
    case default
       bc=FREE_BC
    end select
  end function geocode_to_bc_enum

  ! determine the array of BC from a geometry code
  pure function geocode_to_bc(geocode) result(bc)
    use dictionaries, only: f_err_throw
    implicit none
    character(len=1), intent(in) :: geocode
    integer, dimension(3) :: bc
    select case(geocode)
    case('P')
       bc=PERIODIC
    case('S')
       bc=PERIODIC
       bc(2)=FREE
    case('F')
       bc=FREE
    case('W')
       bc=FREE
       bc(3)=PERIODIC
    case default
       bc=NULL_PAR
    end select
  end function geocode_to_bc

  !>give the associated geocode, 'X' for unknown
  pure function domain_geocode(dom) result(dom_geocode)
    implicit none
    type(domain), intent(in) :: dom
    character(len=1) :: dom_geocode
    !local variables
    logical, dimension(3) :: peri

    peri= dom%bc == PERIODIC
    if (all(peri)) then
       dom_geocode='P'
    else if (.not. any(peri)) then
       dom_geocode='F'
    else if (peri(1) .and. .not. peri(2) .and. peri(3)) then
       dom_geocode='S'
    else if (.not. peri(1) .and. .not. peri(2) .and. peri(3)) then
       dom_geocode='W'
    else
       dom_geocode='X'
    end if

  end function domain_geocode

  !> returns a logical array of size 3 which is .true. for all the periodic dimensions
  pure function bc_periodic_dims(bc) result(peri)
    implicit none
    integer, dimension(3), intent(in) :: bc
    logical, dimension(3) :: peri
    peri= bc == PERIODIC
  end function bc_periodic_dims

  !> returns a logical array of size 3 which is .true. for all the periodic dimensions
  pure function domain_periodic_dims(dom) result(peri)
    implicit none
    type(domain), intent(in) :: dom
    logical, dimension(3) :: peri
    !local variables

    peri=bc_periodic_dims(dom%bc)

  end function domain_periodic_dims

  !> Calculates the square of the vector r in the domain defined by dom
  !! Takes into account the non-orthorhombicity of the box
  !! with the controvariant metric (dom%gu)
  pure function square(dom,v)
    implicit none
    !> array of coordinate in the domain reference frame
    real(gp), dimension(3), intent(in) :: v
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: square

    if (dom%orthorhombic) then
       square=v(1)**2+v(2)**2+v(3)**2
    else
       square=dotp_gu(dom,v,v)
    end if

  end function square

  function square_add(dom,v_add) result(square)
    implicit none
    !> array of coordinate in the domain reference frame
    real(gp) :: v_add
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: square

    if (dom%orthorhombic) then
       call dotp_external_ortho(v_add,v_add,square)
    else
       call dotp_external_nonortho(dom%gu,v_add,v_add,square)
    end if

  end function square_add

  !> Calculates the square of the vector r in the domain defined by dom
  !! Takes into account the non-orthorhombicity of the box
  !! with the covariant metric (dom%gd)
  pure function square_gd(dom,v)
    implicit none
    !> array of coordinate in the domain reference frame
    real(gp), dimension(3), intent(in) :: v
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: square_gd

    if (dom%orthorhombic) then
       square_gd=v(1)**2+v(2)**2+v(3)**2
    else
       square_gd=dotp_gd(dom,v,v)
    end if

  end function square_gd

  function square_gd_add(dom,v_add) result(square)
    implicit none
    !> array of coordinate in the domain reference frame
    real(gp) :: v_add
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: square

    if (dom%orthorhombic) then
       call dotp_external_ortho(v_add,v_add,square)
    else
       call dotp_external_nonortho(dom%gd,v_add,v_add,square)
    end if

  end function square_gd_add

  pure function dotp_gu(dom,v1,v2)
    implicit none
    real(gp), dimension(3), intent(in) :: v1,v2
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: dotp_gu
    !local variables
    integer :: i,j

    if (dom%orthorhombic) then
       dotp_gu=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
    else
       dotp_gu=0.0_gp
       do i=1,3
          do j=1,3
             dotp_gu=dotp_gu+dom%gu(i,j)*v1(i)*v2(j)
          end do
       end do
    end if

  end function dotp_gu

  function dotp_gu_add2(dom,v1,v2_add) result(dotp)
    implicit none
    real(gp), dimension(3), intent(in) :: v1
    real(gp) :: v2_add !<intent in, cannot be declared as such
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: dotp

    if (dom%orthorhombic) then
       call dotp_external_ortho(v1,v2_add,dotp)
    else
       call dotp_external_nonortho(dom%gu,v1,v2_add,dotp)
    end if

  end function dotp_gu_add2

  function dotp_gu_add1(dom,v1_add,v2) result(dotp)
    implicit none
    real(gp), dimension(3), intent(in) :: v2
    real(gp) :: v1_add !<intent in, cannot be declared as such
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: dotp

    if (dom%orthorhombic) then
       call dotp_external_ortho(v1_add,v2,dotp)
    else
       call dotp_external_nonortho(dom%gu,v1_add,v2,dotp)
    end if

  end function dotp_gu_add1

  pure function dotp_gd(dom,v1,v2)
    implicit none
    real(gp), dimension(3), intent(in) :: v1,v2
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: dotp_gd
    !local variables
    integer :: i,j

    if (dom%orthorhombic) then
       dotp_gd=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
    else
       dotp_gd=0.0_gp
       do i=1,3
          do j=1,3
             dotp_gd=dotp_gd+dom%gd(i,j)*v1(i)*v2(j)
          end do
       end do
    end if

  end function dotp_gd

  function dotp_gd_add2(dom,v1,v2_add) result(dotp)
    implicit none
    real(gp), dimension(3), intent(in) :: v1
    real(gp) :: v2_add !<intent in, cannot be declared as such
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: dotp

    if (dom%orthorhombic) then
       call dotp_external_ortho(v1,v2_add,dotp)
    else
       call dotp_external_nonortho(dom%gd,v1,v2_add,dotp)
    end if

  end function dotp_gd_add2

  function dotp_gd_add1(dom,v1_add,v2) result(dotp)
    implicit none
    real(gp), dimension(3), intent(in) :: v2
    real(gp) :: v1_add !<intent in, cannot be declared as such
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: dotp

    if (dom%orthorhombic) then
       call dotp_external_ortho(v1_add,v2,dotp)
    else
       call dotp_external_nonortho(dom%gd,v1_add,v2,dotp)
    end if

  end function dotp_gd_add1

  function dotp_gd_add12(dom,v1_add,v2_add) result(dotp)
    implicit none
    real(gp) :: v1_add,v2_add !<intent in, cannot be declared as such
    type(domain), intent(in) :: dom !<definition of the domain
    real(gp) :: dotp

    if (dom%orthorhombic) then
       call dotp_external_ortho(v1_add,v2_add,dotp)
    else
       call dotp_external_nonortho(dom%gd,v1_add,v2_add,dotp)
    end if

  end function dotp_gd_add12

  !>gives the value of the coordinates for an orthorhombic reference system
  !! from their value wrt a nonorthorhombic system
  pure function rxyz_ortho(dom,rxyz)
    implicit none
    type(domain), intent(in) :: dom
    real(gp), dimension(3), intent(in) :: rxyz
    real(gp), dimension(3) :: rxyz_ortho
    ! local variables
    integer :: i,j

    if (dom%orthorhombic) then
     rxyz_ortho(1:3)=rxyz(1:3)
    else
     do i=1,3
      rxyz_ortho(i)=0.0_gp
      do j=1,3
       rxyz_ortho(i)=rxyz_ortho(i)+dom%uabc(i,j)*rxyz(j)
      end do
     end do
    end if

  end function rxyz_ortho

  !>gives the value of the coordinates for a nonorthorhombic reference system
  !! from their value wrt an orthorhombic system
  pure function rxyz_nonortho(dom,rxyz,mtmp)
    implicit none
    type(domain), intent(in) :: dom
    real(gp), dimension(3), intent(in) :: rxyz
    real(gp), dimension(3,3), intent(in) :: mtmp
    real(gp), dimension(3) :: rxyz_nonortho
    ! local variables
    integer :: i,j

    if (dom%orthorhombic) then
     rxyz_nonortho(1:3)=rxyz(1:3)
    else
     do i=1,3
      rxyz_nonortho(i)=0.0_gp
      do j=1,3
       rxyz_nonortho(i)=rxyz_nonortho(i)+mtmp(i,j)*rxyz(j)
      end do
     end do
    end if

  end function rxyz_nonortho

  pure function distance(dom,r,c) result(d)
    use dictionaries, only: f_err_throw
    implicit none
    real(gp), dimension(3), intent(in) :: r,c
    type(domain), intent(in) :: dom
    real(gp) :: d
    !local variables
    integer :: i !,j,k,ii
    real(gp) :: d2!,dold
    real(gp), dimension(3) :: rt!,ri,ci

!!$    rt=closest_r(mesh,r,c)
!!$    d2=square_gd(mesh,rt)
!!$    d=sqrt(d2)

    d=0.0_gp
    if (dom%orthorhombic) then
       d2=0.0_gp
       do i=1,3
          d2=d2+r_wrap(dom%bc(i),dom%acell(i),r(i),c(i))**2
       end do
       d=sqrt(d2)
    else
       rt=closest_r(dom,r,c)
       d2=square_gd(dom,rt)
       d=sqrt(d2)
!       dold=1.0d100 !huge_number
!       do ii=1,3
!        if (mesh%bc(ii)==PERIODIC) then
!          ri(ii)=mod(r(ii),mesh%ndims(ii)*mesh%hgrids(ii))
!          ci(ii)=mod(c(ii),mesh%ndims(ii)*mesh%hgrids(ii))
!        else
!          ri(ii)=r(ii)
!          ci(ii)=c(ii)
!        end if
!       end do
!       do i=-mesh%bc(1),mesh%bc(1)
!        do j=-mesh%bc(2),mesh%bc(2)
!         do k=-mesh%bc(3),mesh%bc(3)
!            rt(1)=ri(1)+real(i,kind=8)*mesh%ndims(1)*mesh%hgrids(1)
!            rt(2)=ri(2)+real(j,kind=8)*mesh%ndims(2)*mesh%hgrids(2)
!            rt(3)=ri(3)+real(k,kind=8)*mesh%ndims(3)*mesh%hgrids(3)
!            d2=square_gd(mesh,rt-ci)
!            d=sqrt(d2)
!            if (d.lt.dold) then
!               dold=d
!            end if
!         end do
!        end do
!       end do
!       d=dold
    end if

  end function distance

  !> Calculates the minimum difference between two coordinates
  !!@warning: this is only valid if the coordinates wrap once.
  pure function r_wrap(bc,alat,ri,ci)
    implicit none
    integer, intent(in) :: bc
    real(gp), intent(in) :: ri,ci,alat
    real(gp) :: r_wrap
    ! local variables
    real(gp) :: r,c

    !for periodic BC calculate mindist only if the center of mass can be defined without the modulo
    r=ri
    c=ci
    r_wrap=r-c
    if (bc==PERIODIC) then
      r=mod(ri,alat)
      c=mod(ci,alat)
      r_wrap=r-c
       if (abs(r_wrap) > 0.5_gp*alat) then
          if (r < 0.5_gp*alat) then
             r_wrap=r+alat-c
          else
             r_wrap=r-alat-c
          end if
       end if
    end if

  end function r_wrap

!!$  pure function min_dist(bc,alat,r,r_old)
!!$    implicit none
!!$    integer, intent(in) :: bc
!!$    real(gp), intent(in) :: r,r_old,alat
!!$    real(gp) :: min_dist
!!$
!!$    !for periodic BC calculate mindist only if the center of mass can be defined without the modulo
!!$    min_dist=abs(r-r_old)
!!$    if (bc==PERIODIC) then
!!$       if (min_dist > 0.5_gp*alat) then
!!$          if (r < 0.5_gp*alat) then
!!$             min_dist=abs(r+alat-r_old)
!!$          else
!!$             min_dist=abs(r-alat-r_old)
!!$          end if
!!$       end if
!!$    end if
!!$
!!$  end function min_dist

  !>find the closest center according to the periodiciy of the
  !! box and provide the vector
  pure function closest_r(dom,v,center) result(r)
    implicit none
    real(gp), dimension(3), intent(in) :: v,center
    type(domain), intent(in) :: dom
    real(gp), dimension(3) :: r
    !local variables
    integer :: i,j,k,ii,icurr,jcurr,kcurr
    real(gp) :: d,d2,dold
    real(gp), dimension(3) :: rt,ri,ci!,c_ortho,r_ortho

    if (dom%orthorhombic) then
       do i=1,3
          r(i)=r_wrap(dom%bc(i),dom%acell(i),v(i),center(i))
       end do
    else
       dold=1.0d100 !huge_number
       icurr=0
       jcurr=0
       kcurr=0
       do ii=1,3
        if (dom%bc(ii)==PERIODIC) then
          ri(ii)=mod(v(ii),dom%acell(ii))
          ci(ii)=mod(center(ii),dom%acell(ii))
        else
          ri(ii)=v(ii)
          ci(ii)=center(ii)
        end if
       end do
       ri=ri-ci
       do i=-dom%bc(1),dom%bc(1)
        do j=-dom%bc(2),dom%bc(2)
         do k=-dom%bc(3),dom%bc(3)
            rt(1)=ri(1)+real(i,gp)*dom%acell(1)
            rt(2)=ri(2)+real(j,gp)*dom%acell(2)
            rt(3)=ri(3)+real(k,gp)*dom%acell(3)
            d2=square_gd(dom,rt)!-ci)
            d=sqrt(d2)
            if (d.lt.dold) then
               dold=d
               icurr=i
               jcurr=j
               kcurr=k
            end if
         end do
        end do
       end do
       d=dold
       r(1)=ri(1)+real(icurr,gp)*dom%acell(1)! - ci(1)
       r(2)=ri(2)+real(jcurr,gp)*dom%acell(2)! - ci(2)
       r(3)=ri(3)+real(kcurr,gp)*dom%acell(3)! - ci(3)
    end if

  end function closest_r

!!$  subroutine parse_xyz_header(ios,dom_dict)
!!$    use f_iostream
!!$    implicit none
!!$    type(io_stream), intent(inout) :: ios
!!$    type(dictionary), pointer :: dom_dict
!!$    !local variables
!!$    logical :: eol
!!$    character(len=256) :: line
!!$
!!$    !first line of the header
!!$    eol=f_iostream_line_read(ios,&
!!$         '[nat,units_str,energy,comment]',&
!!$         dom_dict)
!!$    !second line of the header
!!$    eol=f_iostream_line_read(ios,&
!!$         '[boundary_conditions,alat(3),properties(dict)]',&
!!$         dom_dict)
!!$
!!$    !check eol, should be always true
!!$
!!$  end subroutine parse_xyz_header
!!$
!!$  subroutine parse_ascii_header(ios,dom_dict)
!!$    use f_iostream
!!$    implicit none
!!$    type(io_stream), intent(inout) :: ios
!!$    type(dictionary), pointer :: dom_dict
!!$    !local variables
!!$    logical :: eof
!!$    character(len=256) :: line
!!$
!!$    !second line of the header
!!$    eol=f_iostream_line_read(ios,&
!!$         '[boundary_conditions,alat(3),properties(dict)]',&
!!$         dom_dict)
!!$
!!$  end subroutine parse_ascii_header



!!$  !direct metric
!!$  function gd(bc,angrad)
!!$  end function gd
!!$
!!$  !direct metric
!!$  function gu(bc,angrad)
!!$  end function gu
!!$
!!$

end module at_domain

subroutine dotp_external_ortho(v1,v2,dotp)
  use f_precisions, only: gp=>f_double
  implicit none
  real(gp), dimension(3), intent(in) :: v1,v2
  real(gp), intent(out) :: dotp

  dotp=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
end subroutine dotp_external_ortho

subroutine dotp_external_nonortho(g,v1,v2,dotp)
  use f_precisions, only: gp=>f_double
  implicit none
  real(gp), dimension(3,3), intent(in) :: g
  real(gp), dimension(3), intent(in) :: v1,v2
  real(gp), intent(out) :: dotp
  !local variables
  integer :: i,j

       dotp=0.0_gp
       do i=1,3
          do j=1,3
             dotp=dotp+g(i,j)*v1(i)*v2(j)
          end do
       end do
end subroutine dotp_external_nonortho
