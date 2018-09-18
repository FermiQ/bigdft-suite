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

  !>data type for the simulation domain
  type, public :: domain
     integer :: units=NULL_PAR
     integer, dimension(3) :: bc=NULL_PAR !< boundary conditions on each direction
     real(gp), dimension(3) :: angrad=0.0_gp !<angles between the dimensions in radians (alpha_bc,beta_ac,gamma_bc)
     real(gp), dimension(3,3) :: abc=0.0_gp !<matrix of the normalized translation vectors direction
     real(gp), dimension(3) :: acell=0.0_gp !lengths of the primitive vectors abc
  end type domain

  public :: domain_null,domain_set_from_dict,domain_geocode,units_enum_from_str

contains

  function domain_null()
    type(domain) :: domain_null
  end function domain_null

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

    dom%units=toi(units)
    dom%bc=toi(bc)

    !this means free BC
    if (.not. any(bc_is_periodic(bc))) return

    !check of the availablilty of the options
    if (present(acell) .eqv. present(abc)) call f_err_throw('Domain inconsistency: only acell or abc may be given')

    if (any([present(alpha_bc),present(beta_ac),present(gamma_ab)]) .and. present(abc)) &
         call f_err_throw('Domain inconsistency: angles cannot be present if unit cell vectors are given')
    
    if (present(abc)) then
       dom%abc=abc
       call abc_to_angrad_and_acell(abc,dom%angrad,dom%acell)
    else
       dom%acell=acell
       !default orthorhombic
       dom%angrad=onehalf*pi
       if (present(alpha_bc)) dom%angrad(BC_)=alpha_bc
       if (present(beta_ac)) dom%angrad(AC_)=beta_ac
       if (present(gamma_ab)) dom%angrad(AB_)=gamma_ab      
       call angrad_and_acell_to_abc(dom%angrad,dom%acell,dom%abc)
    end if

    call f_assert(all(dom%angrad > 0.0_gp),'Domain inconsistency: some of the angles are not positive')
    if (dom%angrad(AB_) /= onehalf*pi .and. .not.(bc_is_periodic(bc(A_)) .and. bc_is_periodic(bc(B_)))) &
         call f_err_throw('Domain inconsistency: cannot provide angle gamma_ab if the dimensions a and b are not periodic')
    if (dom%angrad(BC_) /= onehalf*pi .and. .not.(bc_is_periodic(bc(B_)) .and. bc_is_periodic(bc(C_)))) &
         call f_err_throw('Domain inconsistency: cannot provide angle alpha_bc if the dimensions b and c are not periodic')
    if (dom%angrad(AC_) /= onehalf*pi .and. .not.(bc_is_periodic(bc(A_)) .and. bc_is_periodic(bc(C_)))) &
         call f_err_throw('Domain inconsistency: cannot provide angle beta_ac if the dimensions a and c are not periodic')

     call f_assert(get_abc_angrad_acell_consistency(dom%abc,dom%angrad,dom%acell)< 1.e-10,&
          'Domain inconsisntency: bad definition of angles angrad or acell wrt abc')

  end function domain_new

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
       abc(Z_,C_)=sin(angrad(AC_))
    else if (angrad(BC_)/= onehalf*pi .and. angrad(AB_) == onehalf*pi .and. angrad(AC_)==onehalf*pi ) then
       abc=0.0_gp      
       abc(X_,A_)=1.0_gp
       abc(Y_,B_)=1.0_gp
       abc(Y_,C_)=cos(angrad(BC_))
       abc(Z_,C_)=sin(angrad(BC_))
    else if (angrad(AB_)/= onehalf*pi .and. angrad(AC_) == onehalf*pi .and. angrad(BC_)==onehalf*pi ) then
       abc=0.0_gp      
       abc(X_,A_)=1.0_gp
       abc(X_,B_)=cos(angrad(AB_))
       abc(Y_,B_)=sin(angrad(AB_))
       abc(Z_,C_)=1.0_gp
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

  pure elemental function bc_enum_from_int(ibc) result(bc)
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

  pure function units_enum_from_int(iunit) result(units)
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

    if (all(dom%bc==PERIODIC_BC) .and. (DOMAIN_ABC .in. dict)) then
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
       call bc_and_cell_to_dict(dict//DOMAIN_CELL,bc_enum_from_int(dom%bc),dom%acell)
       if (dom%angrad(BC_) /= onehalf*pi) call set(dict//DOMAIN_ALPHA,dom%angrad(BC_))
       if (dom%angrad(AC_) /= onehalf*pi) call set(dict//DOMAIN_BETA,dom%angrad(AC_))
       if (dom%angrad(AB_) /= onehalf*pi) call set(dict//DOMAIN_GAMMA,dom%angrad(AB_))
    end if

  end subroutine domain_merge_to_dict

  ! determine the array of BC from a geometry code
  pure function geocode_to_bc(geocode) result(bc)
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
  end function geocode_to_bc
 
  !>give the associated geocode, 'X' for unknown
  pure function domain_geocode(dom) result(cell_geocode)
    implicit none
    type(domain), intent(in) :: dom
    character(len=1) :: cell_geocode
    !local variables
    logical, dimension(3) :: peri

    peri= dom%bc == PERIODIC
    if (all(peri)) then
       cell_geocode='P'
    else if (.not. any(peri)) then
       cell_geocode='F'
    else if (peri(1) .and. .not. peri(2) .and. peri(3)) then
       cell_geocode='S'
    else if (.not. peri(1) .and. .not. peri(2) .and. peri(3)) then
       cell_geocode='W'
    else
       cell_geocode='X'
    end if

  end function domain_geocode 

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
