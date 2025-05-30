!> @file
!!   Routines used to find the Fermi level during FOE
!! @author
!!   Copyright (C) 2016 CheSS developers
!!
!!   This file is part of CheSS.
!!   
!!   CheSS is free software: you can redistribute it and/or modify
!!   it under the terms of the GNU Lesser General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.
!!   
!!   CheSS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU Lesser General Public License for more details.
!!   
!!   You should have received a copy of the GNU Lesser General Public License
!!   along with CheSS.  If not, see <http://www.gnu.org/licenses/>.


!> Determination of the Fermi level for the density matrix
module fermi_level
  use dictionaries, only: f_err_throw
  use sparsematrix_base
  use wrapper_linalg
  implicit none

  private

  ! Public parameters
  !> Function to determine the occupation numbers
  integer, parameter, public :: SMEARING_DIST_ERF   = 1  !< Tends to 0 and 1 faster \f$1/2\left[1-erf\left(\frac{E-\mu}{\delta E}\right)\right]\f$
  integer, parameter, public :: SMEARING_DIST_FERMI = 2  !< Normal Fermi distribution i.e.\f$\frac{1}{1+e^{E-\mu}/k_BT}\f$
  integer, parameter, public :: SMEARING_DIST_COLD1 = 3  !< Marzari's cold smearing with a=-.5634 (bumb minimization)
  integer, parameter, public :: SMEARING_DIST_COLD2 = 4  !< Marzari's cold smearing with a=-.8165 (monotonic tail)
  integer, parameter, public :: SMEARING_DIST_METPX = 5  !< Methfessel and Paxton (same as COLD with a=0)

  ! Public routines
  public :: init_fermi_level
  public :: determine_fermi_level
  public :: fermilevel_get_real
  public :: fermilevel_get_logical
  public :: eval_to_occ
  !!public :: get_roots_of_cubic_polynomial
  !!public :: determinant

  ! Auxiliary structure that holds the required data
  type,public :: fermi_aux
    logical :: adjust_lower_bound, adjust_upper_bound
    real(kind=mp) :: target_charge, bisection_shift, ef_old, sumn_old
    real(kind=mp) :: ef_interpol_chargediff, ef_interpol_det
    real(kind=mp),dimension(2) :: sumnarr, efarr
    logical,dimension(2) :: bisection_bounds_ok
    real(kind=mp),dimension(4,4) :: interpol_matrix
    real(kind=mp),dimension(4) :: interpol_vector
    integer :: it, it_solver, verbosity
  end type fermi_aux


  contains
    
    !> Initialize the internal variables
    subroutine init_fermi_level(target_charge, ef, f, bisection_shift, ef_interpol_chargediff, ef_interpol_det, verbosity)
      use dynamic_memory
      implicit none

      ! Calling arguments
      real(kind=mp),intent(in) :: target_charge                   !< total charge of the system
      real(kind=mp),intent(in) :: ef                              !< initial guess for the fermi level
      type(fermi_aux),intent(out) :: f                           !< type that holds the internal data
      real(kind=mp),intent(in),optional :: bisection_shift        !< shift to be used for the determination of the bisection bounds
      real(kind=mp),intent(in),optional :: ef_interpol_chargediff !< charge difference below which the cubic interpolation is allowed
      real(kind=mp),intent(in),optional :: ef_interpol_det        !< determinant of the interpolation matrix above which the cubic interpolation is allowed
      integer,intent(in),optional :: verbosity                   !< verbosity of the output: 0 for no output, 1 for more detailed output

      call f_routine(id='init_fermi_level')

      f%adjust_lower_bound = .true.
      f%adjust_upper_bound = .true.
      f%target_charge = target_charge
      f%sumnarr(1) = 0.d0
      f%sumnarr(2) = 1.d100
      if (present(bisection_shift)) then
          f%bisection_shift = bisection_shift
      else
          f%bisection_shift = 0.1d0
      end if
      f%ef_old = 0.d0
      f%sumn_old = 0.d0
      f%efarr(1) = ef - f%bisection_shift
      f%efarr(2) = ef + f%bisection_shift
      if (present(ef_interpol_chargediff)) then
          f%ef_interpol_chargediff = ef_interpol_chargediff
      else
          f%ef_interpol_chargediff = 1.d0
      end if
      if (present(ef_interpol_det)) then
          f%ef_interpol_det = ef_interpol_det
      else
          f%ef_interpol_det = 1.d-20
      end if
      f%bisection_bounds_ok(1) = .false.
      f%bisection_bounds_ok(2) = .false.
      f%interpol_matrix(:,:) = 0.d0
      f%interpol_vector(:) = 0.d0
      f%it = 0
      f%it_solver = 0
      if (present(verbosity)) then
          f%verbosity = verbosity
      else
          f%verbosity = 0
      end if

      call f_release_routine()

    end subroutine init_fermi_level



    subroutine determine_fermi_level(iproc, f, sumn, ef, info)
      use dynamic_memory
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc          !< task ID
      type(fermi_aux),intent(inout) :: f   !< type that holds the internal data
      real(kind=mp),intent(in) :: sumn      !< charge of the system (which should be equal to the target charge once the correct Fermi level is found),
                                           !    obtained with the current value of ef
      real(kind=mp),intent(inout) :: ef     !< on input: current value of the Fermi level
                                           !  on output: new guess for the Fermi level, depending on the value of info
      integer,intent(out),optional :: info !< info parameter: * -1: adjusting the lower bisection bound, ef not meaningful
                                           !                  * -2: adjusting the upper bisection bound, ef not meaningful
                                           !                  *  0: searching the correct fermi level, ef meaningful
      ! Local variables
      real(kind=mp) :: charge_diff
      logical :: interpolation_possible
      integer :: internal_info


!      integer :: iproc,ierr

      call f_routine(id='determine_fermi_level')

      ! Make sure that the bounds for the bisection are negative and positive
      charge_diff = sumn-f%target_charge
!iproc =mpirank(bigdft_mpi%mpi_comm)
!      call mpi_comm_rank(mpi_comm_world, iproc, ierr)
      if (f%adjust_lower_bound) then
          if (charge_diff <= 0.d0) then
              ! Lower bound okay
              f%adjust_lower_bound = .false.
              f%bisection_shift = f%bisection_shift*0.9d0
              f%sumnarr(1) = sumn
              f%bisection_bounds_ok(1) = .true.
          else
              f%efarr(1) = f%efarr(1)-f%bisection_shift
              f%bisection_shift = f%bisection_shift*1.1d0
              f%bisection_bounds_ok(1) = .false.
          end if
      else if (f%adjust_upper_bound) then
          if (charge_diff >= 0.d0) then
              ! Upper bound okay
              f%adjust_upper_bound = .false.
              f%bisection_shift = f%bisection_shift*0.9d0
              f%sumnarr(2) = sumn
              f%bisection_bounds_ok(2) = .true.
          else
              f%efarr(2) = f%efarr(2)+f%bisection_shift
              f%bisection_shift = f%bisection_shift*1.1d0
              f%bisection_bounds_ok(2)=.false.
          end if
      end if


      internal_info = 0
      if (f%adjust_lower_bound) then
          ef = f%efarr(1)
          internal_info = -1
      else if (f%adjust_upper_bound) then
          ef = f%efarr(2)
          internal_info = -2
      end if
      if (present(info)) then
          info = internal_info
      end if

      if (internal_info < 0) then
          call f_release_routine()
          if (f%verbosity>=1 .and. iproc==0) call yaml_map('new eF','bisec bounds')
          return ! no need to proceed further
      end if

      ! If we made it here, the bounds are ok (i.e. the lower bound gives a negative charge difference and the upper one a positive one).


      ! Adjust the bounds for the bisection, i.e. make the search interval more narrow
      if (charge_diff < 0.d0) then
          f%efarr(1) = ef
          f%sumnarr(1) = sumn
      else if (charge_diff >= 0.d0) then
          f%efarr(2) = ef
          f%sumnarr(2) = sumn
      end if


      f%it_solver = f%it_solver+1

      ! Check whether the system behaves reasonably.
      interpolation_possible=.true.
      if (f%it_solver > 1) then
          if (ef > f%ef_old .and. sumn < f%sumn_old) then
              interpolation_possible = .false.
          else if (ef < f%ef_old .and. sumn > f%sumn_old) then
              interpolation_possible = .false.
          end if
          if (abs(sumn-f%sumn_old)<1.d-10) then
              interpolation_possible = .false.
          end if
          if (f%verbosity >= 2 .and. iproc==0) then
              call yaml_newline()
              call yaml_mapping_open('interpol check',flow=.true.)
                 call yaml_map('D eF',ef-f%ef_old,fmt='(es13.6)')
                 call yaml_map('D Tr',sumn-f%sumn_old,fmt='(es13.6)')
                 call yaml_map('interpol possible',interpolation_possible)
              call yaml_mapping_close()
              call yaml_newline()
           end if
      end if
      if (.not.interpolation_possible) then
          ! Set the history for the interpolation to zero.
          f%it_solver=0
      end if

      f%ef_old = ef
      f%sumn_old = sumn



      call determine_new_fermi_level()

      call f_release_routine()

      contains

        subroutine determine_new_fermi_level()
          implicit none
          integer :: info, i, ii
          real(kind=mp) :: m, b, ef_interpol, det
          real(kind=mp),dimension(4,4) :: tmp_matrix
          real(kind=mp),dimension(4) :: interpol_solution
          integer,dimension(4) :: ipiv
          logical :: interpolation_nonsense, cubicinterpol_possible

          !call yaml_map('sumn',sumn)

          ! Shift up the old results.
          if (f%it_solver>4) then
              do i=1,4
                  f%interpol_matrix(1,i)=f%interpol_matrix(2,i)
                  f%interpol_matrix(2,i)=f%interpol_matrix(3,i)
                  f%interpol_matrix(3,i)=f%interpol_matrix(4,i)
              end do
              f%interpol_vector(1)=f%interpol_vector(2)
              f%interpol_vector(2)=f%interpol_vector(3)
              f%interpol_vector(3)=f%interpol_vector(4)
          end if
          !LG: if f%it_solver==0 this index comes out of bounds!
          ii=max(min(f%it_solver,4),1)
          f%interpol_matrix(ii,1)=ef**3
          f%interpol_matrix(ii,2)=ef**2
          f%interpol_matrix(ii,3)=ef
          f%interpol_matrix(ii,4)=1.d0
          f%interpol_vector(ii)=sumn-f%target_charge
        
          ! Solve the linear system f%interpol_matrix*interpol_solution=f%interpol_vector
          if (f%it_solver>=4) then
              ! Calculate the determinant of the matrix used for the interpolation
              det=determinant(iproc,4,f%interpol_matrix)
              if (abs(det) > f%ef_interpol_det) then
                  cubicinterpol_possible = .true.
                  do i=1,ii
                      interpol_solution(i)=f%interpol_vector(i)
                      tmp_matrix(i,1)=f%interpol_matrix(i,1)
                      tmp_matrix(i,2)=f%interpol_matrix(i,2)
                      tmp_matrix(i,3)=f%interpol_matrix(i,3)
                      tmp_matrix(i,4)=f%interpol_matrix(i,4)
                  end do
                  !if (bigdft_mpi%iproc==0) then
                     !call yaml_map('matrix',tmp_matrix,fmt='(es10.3)')
                     !call yaml_map('interpol_vector',f%interpol_vector,fmt='(es12.5)')
                     !call yaml_newline()
                     !call yaml_map('solution',interpol_solution,fmt='(es10.3)')
                     !call yaml_map('determinant',determinant(bigdft_mpi%iproc,4,f%interpol_matrix),fmt='(es10.3)')
                  call dgesv(ii, 1, tmp_matrix, 4, ipiv, interpol_solution, 4, info)
                  if (info/=0) then
                     if (iproc==0) write(*,'(1x,a,i0)') 'ERROR in dgesv (FOE), info=',info
                  end if
        
                  !if (bigdft_mpi%iproc==0) call yaml_map('a x^3+b x^2 + c x + d',interpol_solution,fmt='(es10.3)')
                  call get_roots_of_cubic_polynomial(interpol_solution(1), interpol_solution(2), &
                       interpol_solution(3), interpol_solution(4), ef, ef_interpol)
                  !if (bigdft_mpi%iproc==0) then
                  !    call yaml_newline()
                  !    call yaml_map('zero of cubic polynomial',ef_interpol,fmt='(es10.3)')
                  !end if
                  ! Sanity check: If the charge was too small, then new new guess
                  ! for the Fermi energy must be larger than the actual value, and
                  ! analogously if the charge was too large.
                  interpolation_nonsense = .false.
                  if (f%interpol_vector(ii)<0) then
                      ! Charge too small, the new guess must be larger
                      if (ef_interpol<ef) interpolation_nonsense = .true.
                  else if (f%interpol_vector(ii)>0) then
                      ! Charge too large, the new guess must be smaller
                      if (ef_interpol>ef) interpolation_nonsense = .true.
                  end if
              else
                  cubicinterpol_possible = .false.
              end if
          end if
        
          ! Calculate the new Fermi energy.
          if (f%verbosity>=2 .and. iproc==0) then
              call yaml_newline()
              call yaml_mapping_open('Search new eF',flow=.true.)
          end if
          if (f%it_solver>=4 .and.  &
              abs(sumn-f%target_charge) < f%ef_interpol_chargediff) then! .and. &
              !.not.interpolation_nonsense) then
              !det=determinant(bigdft_mpi%iproc,4,f%interpol_matrix)
              if (f%verbosity >= 2 .and. iproc==0) then
                  call yaml_map('det',det,fmt='(es10.3)')
                  call yaml_map('limit',f%ef_interpol_det,fmt='(es10.3)')
              end if
              !if(abs(det) > f%ef_interpol_det) then
              if(cubicinterpol_possible .and. .not.interpolation_nonsense) then
                  ef = ef_interpol
                  if (f%verbosity>=1 .and. iproc==0) call yaml_map('new eF','cubic interpol')
              else
                  ! linear interpolation
                  m = (f%interpol_vector(4)-f%interpol_vector(3))/(f%interpol_matrix(4,3)-f%interpol_matrix(3,3))
                  b = f%interpol_vector(4)-m*f%interpol_matrix(4,3)
                  ef = -b/m
                  if (f%verbosity>=1 .and. iproc==0) call yaml_map('new eF','linear interpol')
              end if
          else
              ! Use mean value of bisection and secant method if possible,
              ! otherwise only the bisection.
              ! Bisection solution
              ef = 0.5d0*(f%efarr(1)+f%efarr(2))
              !write(*,'(a,i6,5es16.7)') 'iproc, f%efarr, f%sumnarr, f%target_charge', &
              !     bigdft_mpi%iproc, f%efarr, f%sumnarr, f%target_charge
              if (abs(f%sumnarr(2)-f%sumnarr(1))>1.d-6) then !otherwise secant method numerically unstable
                  ! Add Secant method solution
                  ef = ef + f%efarr(2)-(f%sumnarr(2)-f%target_charge)*(f%efarr(2)-f%efarr(1))/(f%sumnarr(2)-f%sumnarr(1))
                  ! Take the mean value
                  ef = 0.5d0*ef
                  if (f%verbosity>=1 .and. iproc==0) call yaml_map('new eF','bisection/secant')
              else
                  if (f%verbosity>=1 .and. iproc==0) call yaml_map('new eF','bisection')
              end if
          end if
          if (f%verbosity>=2 .and. iproc==0) then
              !call yaml_map('guess for new ef',ef,fmt='(es15.8)')
              call yaml_mapping_close()
          end if

        end subroutine determine_new_fermi_level


    end subroutine determine_fermi_level

    function fermilevel_get_real(f, fieldname) result(val)
        ! Calling arguments
        type(fermi_aux),intent(in) :: f      !< type that holds the internal data
        character(len=*),intent(in) :: fieldname
        real(kind=mp) :: val

        !LG: wouldn't be better to use explicit meaning of the 
        !variables instead of their name in the structure?
        select case (trim(fieldname))
        case ("efarr(1)")
            val = f%efarr(1)
        case ("efarr(2)")
            val = f%efarr(2)
        case ("bisection_shift")
            val = f%bisection_shift
        case default
           val=0.d0
            call f_err_throw("wrong argument for "//trim(fieldname))
        end select
    end function fermilevel_get_real

    function fermilevel_get_logical(f, fieldname) result(val)
        ! Calling arguments
        type(fermi_aux),intent(in) :: f      !< type that holds the internal data
        character(len=*),intent(in) :: fieldname
        logical :: val

        select case (trim(fieldname))
        case ("bisection_bounds_ok(1)")
            val = f%bisection_bounds_ok(1)
        case ("bisection_bounds_ok(2)")
            val = f%bisection_bounds_ok(2)
        case default
           val=.false.
            call f_err_throw("wrong argument for "//trim(fieldname))
        end select
    end function fermilevel_get_logical


    ! Finds the real root of the equation ax**3 + bx**2 + cx + d which is closest to target_solution
    subroutine get_roots_of_cubic_polynomial(a, b, c, d, target_solution, solution)
      implicit none
    
      ! Calling arguments
      real(kind=mp),intent(in) :: a, b, c, d
      real(kind=mp),intent(in) :: target_solution
      real(kind=mp),intent(out) :: solution
    
      ! Local variables
      complex(kind=mp) :: a_c, b_c, c_c, d_c, Q_c, S_c, ttp_c, ttm_c
      complex(kind=mp),dimension(3) :: sol_c
      !double complex :: test
      real(kind=mp) :: ttmin, tt
      integer :: i
    
      a_c=cmplx(a,0.d0,kind=mp)
      b_c=cmplx(b,0.d0,kind=mp)
      c_c=cmplx(c,0.d0,kind=mp)
      d_c=cmplx(d,0.d0,kind=mp)
    
      Q_c = sqrt( (2*b_c**3-9*a_c*b_c*c_c+27*a_c**2*d_c)**2 - 4*(b_c**2-3*a_c*c_c)**3 )
      S_c = ( .5d0*(Q_c+2*b_c**3-9*a_c*b_c*c_c+27*a_c**2*d_c) )**(1.d0/3.d0)
      ttp_c = cmplx(1.d0,sqrt(3.d0),kind=mp)
      ttm_c = cmplx(1.d0,-sqrt(3.d0),kind=mp)
    
      sol_c(1) = -b_c/(3*a_c) &
           - S_c/(3*a_c) &
           - (b_c**2-3*a_c*c_c)/(3*a_c*S_c)
      sol_c(2) = -b_c/(3*a_c) + (S_c*ttp_c)/(6*a_c) + ttm_c*(b_c**2-3*a_c*c_c)/(6*a_c*S_c)
      sol_c(3) = -b_c/(3*a_c) + (S_c*ttm_c)/(6*a_c) + ttp_c*(b_c**2-3*a_c*c_c)/(6*a_c*S_c)
      !!if (iproc==0) then
      !!    write(*,*) 'sol 1', sol_c(1)
      !!    write(*,*) 'sol 2', sol_c(2)
      !!    write(*,*) 'sol 3', sol_c(3)
      !!end if
    
      ! Select the real solution that is closest to target_solution
      ttmin=1.d100
      do i=1,3
          if (abs(aimag(sol_c(i)))>1.d-14) cycle !complex solution
          tt=abs(real(sol_c(i),kind=mp)-target_solution)
          if (tt<ttmin) then
              ttmin=tt
              solution=real(sol_c(i),kind=mp)
          end if
      end do
    
    end subroutine get_roots_of_cubic_polynomial

    !> such function should be moved in the linalg wrappers
    real(kind=mp) function determinant(iproc, n, mat)
        implicit none
    
        ! Calling arguments
        integer,intent(in) :: iproc, n
        real(kind=mp),dimension(n,n),intent(in) :: mat
    
        ! Local variables
        integer :: i, info
        integer,dimension(n) :: ipiv
        real(kind=mp),dimension(n,n) :: mat_tmp
        real(kind=mp) :: sgn
    
        call vcopy(n**2, mat(1,1), 1, mat_tmp(1,1), 1)
    
        call dgetrf(n, n, mat_tmp, n, ipiv, info)
        if (info/=0) then
            if (iproc==0) write(*,'(a,i0,a)') 'ERROR in dgetrf, info=',info,'. Set determinant to zero.'
            determinant=0
            return
        end if
    
        determinant=1.d0
        do i=1,n
            determinant=determinant*mat_tmp(i,i)
        end do
    
        sgn=1.d0
        do i=1,n
            if(ipiv(i)/=i) then
                sgn=-sgn
            end if
        end do
    
        determinant=sgn*determinant   
    
    end function determinant

    pure subroutine get_entropy(occopt,e,ef,kT,s)
      use numerics
      implicit none
      integer, intent(in) :: occopt !< it should become an enumerator
      real(mp), intent(in) :: e !< energy of the state
      real(mp), intent(in) :: ef !< fermi level
      real(mp), intent(in) :: kT !< width of the fermi function
      real(mp), intent(out) :: s
      !local variables
      real(mp) :: arg,x,a,f,df

      if (kT == 0.0_mp) then
         s=0.0_mp
         return
      end if

      arg=(e-ef)/kT
      select case(occopt)
      case(SMEARING_DIST_ERF)
         s=0.5_mp*kT*oneosqrtpi*safe_exp(-(arg)**2)
      case(SMEARING_DIST_FERMI)
         call get_occupation(occopt,e,ef,kT,f,df)
         s=-kT*(f*safe_log(f) + (1.0_mp-f)*safe_log(1._mp-f))
      case(SMEARING_DIST_COLD1)
         a=-.5634d0
         x= -arg
         s=0.5_mp*kT*(a*x**3-x**2+0.5_mp)*safe_exp(-x**2)*oneosqrtpi
      case(SMEARING_DIST_COLD2)
         a=-.8165d0
         x= -arg
         s=0.5_mp*kT*(1.0_mp-sqrt(2.0_mp)*x)*&
              safe_exp(-(x-1.0_mp/sqrt(2.0_mp))**2)*oneosqrtpi
      case(SMEARING_DIST_METPX)
         x= -arg
         s=-0.5_mp*kT*(x**2-0.5_mp)*safe_exp(-x**2)*oneosqrtpi
      case default
         f  = 0.d0
         df = 0.d0
         s = 0.d0
      end select

    end subroutine get_entropy

    pure subroutine get_occupation(occopt,e,ef,kT,f,df)
      use numerics
      implicit none
      integer, intent(in) :: occopt !< it should become an enumerator
      real(mp), intent(in) :: e !< energy of the state
      real(mp), intent(in) :: ef !< fermi level
      real(mp), intent(in) :: kT !< width of the fermi function
      real(mp), intent(out) :: f !< occupation # for band i 
      real(mp), intent(out) :: df  !< df/darg
      !local variables
      real(mp) :: arg,a,res,x

      arg=(e-ef)/kT
      select case(occopt)
      case(SMEARING_DIST_ERF)
         res=safe_erf(arg)
         f =.5d0*(1.d0-res)
         df=-safe_exp(-arg**2)*oneosqrtpi
      case(SMEARING_DIST_FERMI)
         f =1.d0/(1.d0+safe_exp(arg))
         df=-1.d0/(2.d0+safe_exp(arg)+safe_exp(-arg))
      case(SMEARING_DIST_COLD1)
         a=-.5634d0
         x= -arg
         res=safe_erf(x)
         f =.5d0*(1.d0+res +safe_exp(-x**2)*(-a*x**2 + .5d0*a+x)*oneosqrtpi)
         df=-safe_exp(-x**2) * (a*x**3 -x**2 -1.5d0*a*x +1.5d0)*oneosqrtpi   ! df:=df/darg=-df/dx
         !POSHM=(2.D0-qe_erfc(X))+(-2.d0*a*x*x+2*x+a)*EXP(-X*X)/SQRT(PI)/2.d0
      case(SMEARING_DIST_COLD2)
         a=-.8165d0
         x= -arg-1.0_mp/sqrt(2.0_mp)
         res=safe_erf(x)
         !POSHM2=(2.D0-qe_erfc(X-1.d0/sqrt(2.d0)))+ &
         !&  sqrt(2.d0)*exp(-x*x+sqrt(2.d0)*x-0.5d0)/sqrt(pi)
         !f =.5d0*(1.d0+res +safe_exp(-x**2)*(-a*x**2 + .5d0*a+x)*oneosqrtpi)
         !f=0.5_mp*(1.0_mp+res+sqrt(2.0_mp)*safe_exp(-x**2+sqrt(2.0_mp)*x-0.5_mp)*&
         f=0.5_mp*(1.0_mp+res+sqrt(2.0_mp)*safe_exp(-x**2)*oneosqrtpi)
         !@todo calculate this derivative (a is meaningless here)
         df=-safe_exp(-x**2) * (a*x**3 -x**2 -1.5d0*a*x +1.5d0)*oneosqrtpi 
      case(SMEARING_DIST_METPX)
         x= -arg
         res=safe_erf(x)
         f =.5d0*(1.d0+res +safe_exp(-x**2)*x*oneosqrtpi)
         df=-safe_exp(-x**2) * ( -x**2 +1.5d0)*oneosqrtpi 
      case default
         f  = 0.d0
         df = 0.d0
      end select
    end subroutine get_occupation


    !> Finds the fermi level ef for an error function distribution with a width wf
    !! eval are the Kohn Sham eigenvalues and melec is the total number of electrons
    !assume that the eigenvalues are ordered
    subroutine eval_to_occ(iproc, norbu, norbd, norb, nkpts, kwgts, &
               eval, occup, not_initialized, kT, occopt, efermi, eTS, &
               norbu_res, norbd_res)
      !use module_base
      use futile
      use numerics
      !use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level
      !use public_enums
      !use abi_interfaces_numeric, only: abi_derf_ab
      implicit none
      integer,intent(in) :: nkpts, norbu, norbd, norb
      real(mp),dimension(nkpts),intent(in) :: kwgts
      real(mp),dimension(nkpts*norb),intent(in) :: eval
      real(mp),dimension(nkpts*norb),intent(inout) :: occup
      logical, intent(in) :: not_initialized
      integer, intent(in) :: iproc
      integer, intent(in) :: occopt
      real(mp), intent(in) :: kT   ! width of Fermi function, i.e. k*T
      !type(orbitals_data), intent(inout) :: orbs
      real(mp),intent(inout) :: efermi, eTS
      integer, intent(in) :: norbu_res, norbd_res !<restricted values of norbu and norbd where the fermi level has to be found
      !local variables
      !   real(gp), parameter :: pi=3.1415926535897932d0
      !real(mp), dimension(1,1,1) :: fakepsi
      integer :: ikpt,iorb,newnorbu,newnorbd !,info_fermi
      real(mp) :: charge,wf0,ev
      real(mp) :: ef,cutoffu,cutoffd,full
      real(mp) :: f, df, s
      !integer :: ierr
      !type(fermi_aux) :: ft

      if (iproc==0) then
          call yaml_mapping_open('Determine Fermi level and occupation numbers')
          call yaml_map('Smearing method',occopt)
          call yaml_map('Electronic temperature',kT)
          call yaml_mapping_close()
      end if

      !exitfermi=.false.
      !if (iproc.lt.1)  write(1000+iproc,*)  'ENTER Fermilevel',norbu,norbd,occopt

      eTS=0.0_mp

      !this might be misleading: norbd==0 might also mean that we are with only one electron or with two electrons and mpol==2
      if (norbd==0) then
         full=2.d0   ! maximum occupation for closed shell  orbital
      else
         full=1.d0   ! maximum occupation for spin polarized orbital
      endif
     
      newnorbu=min(norbu_res,norbu)
      newnorbd=min(norbd_res,norbd)

      charge=0.0_mp
      do ikpt=1,nkpts
         !number of zero orbitals for the given k-point
         !overall charge of the system
         do iorb=1,norb
            charge=charge+occup(iorb+(ikpt-1)*norb) * kwgts(ikpt)
         end do
      end do

      !!! Send all eigenvalues to all procs (presumably not necessary)
      !!call broadcast_kpt_objects(nproc, nkpts, norb, &
      !!     &   eval, orbs%ikptproc)

      if (not_initialized) then
         !last value as a guess
         efermi = maxval(eval)
         ! Take initial value at gamma point.
         do ikpt=1,nkpts
            do iorb = 1, norb
               if (occup(iorb+(ikpt-1)*norb) <= onehalf*full) then
                  efermi = min(eval(iorb+(ikpt-1)*norb),efermi)
                  exit
               end if
            end do
         end do
      end if
      ef=efermi

      if (kT == 0.0_mp) then
         wf0 = 1.e-5
      else
         wf0 = kT
      end if

!!$      if (wf0 > 0.0_mp) then

      call find_fermi_level(iproc==0,occopt,nkpts,norbu,norbd,newnorbu,newnorbd,&
           eval,kwgts,wf0,charge/full,ef)
      
      !update the occupation number
      do ikpt=1,nkpts
         do iorb=1,norbu + norbd
            call get_occupation(occopt,eval((ikpt-1)*norb+iorb),ef,wf0,&
                 f,df)
            occup((ikpt-1)*norb+iorb)=full* f
            !if(iproc==0) print*,  eval((ikpt-1)*norb+iorb), occup((ikpt-1)*norb+iorb)
         end do
      end do

      if (kT==0.0_mp) call occup_nosmearing(norb*nkpts,occup,1.e-5_mp)

      !update electronic entropy S; eTS=T_ele*S is the electronic entropy term the negative of which is added to energy: Free energy = energy-T*S
      eTS=0.0_mp
      do ikpt=1,nkpts
         do iorb=1,norbu + norbd
            call get_entropy(occopt,eval((ikpt-1)*norb+iorb),ef,kT,s)
            eTS=eTS+full*s*kwgts(ikpt)
         end do
      end do

      cutoffu=0.0_mp
      cutoffd=0.0_mp
      do ikpt=1,nkpts
         ev=eval((ikpt-1)*norb+norbu)
         call get_occupation(occopt,ev,ef,wf0,&
              f,df)
         if (ev > ef) cutoffu=max(cutoffu,f)
         ev=eval((ikpt-1)*norb+norbu+norbd)
         call get_occupation(occopt,ev,ef,wf0,&
              f,df)
         if (ev > ef) cutoffd=max(cutoffd,f)
      enddo

      if ((cutoffu > 1.d-12 .or. cutoffd > 1.d-12) .and. iproc == 0) then
         call yaml_warning('Occupation numbers do not fill all available levels' // &
              ' lastu=' // trim(yaml_toa(cutoffu,fmt='(1pe8.1)')) // &
              ' lastd=' // trim(yaml_toa(cutoffd,fmt='(1pe8.1)')))
      end if
      !if (iproc.lt.1) write(1000+iproc,'(1x,a,1pe21.14,2(1x,e8.1))') 'Fermi level, Fermi distribution cut off at:  ',ef,cutoffu,cutoffd
      !      if (iproc.lt.1) flush(1000+iproc)
      
      efermi=ef

!!$      else if(full==1.0_mp) then
!!$         !call eFermi_nosmearing(iproc,orbs)
!!$         call eFermi_nosmearing(iproc, nkpts, norbu, norbd, norb, eval, occup, efermi)
!!$         ! no entropic term when electronc temperature is zero
!!$      end if

    END SUBROUTINE eval_to_occ

    !only allow for 1,0 or onehalf values
    pure subroutine occup_nosmearing(ntot,occup,tol)
      implicit none
      integer, intent(in) :: ntot
      real(mp), intent(in) :: tol
      real(mp), dimension(ntot), intent(inout) :: occup
      !local variables
      integer :: iorb
      real(mp) :: occ
      do iorb=1,ntot
         occ=occup(iorb)
         if (abs(occ-0.5_mp) < tol) then
            occ=0.5_mp
         else
            occ=real(nint(occ),mp)
         end if
      end do
    end subroutine occup_nosmearing
      

    subroutine eFermi_nosmearing(iproc, nkpts, norbu, norbd, norb, eval, occup, efermi)
       use yaml_output
       implicit none
       integer,intent(in) :: iproc, nkpts, norbu, norbd, norb
       real(mp),dimension(norb,nkpts),intent(in) :: eval
       real(mp),dimension(norb,nkpts),intent(inout) :: occup
       !type(orbitals_data), intent(inout) :: orbs
       real(mp),intent(out) :: efermi
       !local variables
       integer :: iu,id,n,nzeroorbs,ikpt,iorb
       real(mp) :: charge
       real(mp) :: eF
    
       !SM: I think iu and id should be initialized to these values, in case the
       ! large if will not be executed.
       iu=norbu
       id=norbd
       eF = 0._mp
       do ikpt=1,nkpts
          !number of zero orbitals for the given k-point
          nzeroorbs=0
          !overall charge of the system
          charge=0.0_mp
          do iorb=1,norb
             if (occup(iorb,ikpt) == 0.0_mp) then
                nzeroorbs=nzeroorbs+1
             else
                charge=charge+occup(iorb,ikpt)
             end if
          end do
          if (nzeroorbs /= 0 .and. norbd .gt.0) then
             do iorb=1,norbu-1
                if (eval(iorb,ikpt) > eval(iorb+1,ikpt)) &
                   &   write(*,*) 'wrong ordering of up EVs',iorb,iorb+1
             end do
             do iorb=1,norbd-1
                if (eval(iorb+norbu,ikpt) > eval(iorb+1+norbu,ikpt))&
                   &   write(*,*) 'wrong ordering of dw EVs',iorb+norbu,iorb+1+norbu
             enddo
    
             iu=0
             id=0
             n=0
             do while (real(n,mp) < charge)
                if (eval(iu+1,ikpt) <= eval(id+1+norbu,ikpt)) then
                   iu=iu+1
                   eF=eval(iu+1,ikpt)
                else
                   id=id+1
                   eF=eval(id+1+norbu,ikpt)
                endif
                n=n+1
             enddo
             if (iproc==0) then
                !write(*,'(1x,a,1pe21.14,a,i4)') 'Suggested Homo energy level',eF,', Spin polarization',iu-id
                call yaml_map('Suggested Fermi Level',ef,fmt='(1pe21.14)')
                call yaml_map('Suggested Spin pol.',iu-id,fmt='(i4)')
             end if
             !write(*,*) 'up,down, up-down',iu,id,iu-id
          end if
       end do
       efermi=eF
       !assign the values for the occupation numbers
       do iorb=1,iu
          occup(iorb,1)=1.0_mp
       end do
       do iorb=iu+1,norbu
          occup(iorb,1)=0.0_mp
       end do
       do iorb=1,id
          occup(iorb+norbu,1)=1.0_mp
       end do
       do iorb=id+1,norbd
          occup(iorb+norbu,1)=0.0_mp
       end do
    
    END SUBROUTINE eFermi_nosmearing

    pure subroutine correct_fermi(correction_allowed,electrons,dlectrons,charge,kT,wf,ef,found)
      implicit none
      logical, intent(in) :: correction_allowed
      real(mp), intent(in) :: electrons,dlectrons,kT,charge
      real(mp), intent(inout) :: wf,ef
      logical, intent(out) :: found
      !local variables
      real(mp) :: diff,corr

      diff=-charge+electrons
      !if (iproc.lt.1) write(1000+iproc,*) diff,full,melec,real(melec,gp)
      !         if (iproc.lt.1) flush(1000+iproc)
      !if (iproc.lt.1) write(1000+iproc,*) diff,1.d-11*sqrt(electrons),wf
      !if (iproc.lt.1) flush(1000+iproc)
      !Exit criterion satiesfied, Nevertheles do one mor update of fermi level
      found =  (abs(diff) < 1.d-11*sqrt(electrons) .and. wf == kT ) ! Assume noise grows as sqrt(electrons)

      !alternative solution to avoid division by so high value
      !if (dlectrons == 0.d0) dlectrons=1.d-100  !line to be added
      if (dlectrons == 0.d0) then
         !always enter into first case below
         corr=0.d0
         if (diff > 0.d0) corr=1.d0*wf
         if (diff < 0.d0) corr=-1.d0*wf
         if (correction_allowed .and. wf < 0.1d0) wf=2.d0*wf  ! speed up search of approximate Fermi level by using higher Temperature
      else
         corr=diff/abs(dlectrons) ! for case of no-monotonic func. abs is needed
         if (abs(corr) > wf) then   !for such a large correction the linear approximation is not any more valid
            if (corr > 0.d0) corr=1.d0*wf
            if (corr < 0.d0*wf) corr=-1.d0*wf
            if (correction_allowed .and. wf < 0.1d0) wf=2.d0*wf  ! speed up search of approximate Fermi level by using higher Temperature
         else
            wf=max(kT,.5d0*wf)
         endif
      end if
      ef=ef-corr  ! Ef=Ef_guess+corr.
    end subroutine correct_fermi

    pure subroutine electrons_and_delectrons(occopt,nkpts,norbu,norbd,nu_restricted,nd_restricted,eval,kwgts,ef,kT,&
         electrons,dlectrons)
      implicit none
      integer, intent(in) :: occopt,nkpts,norbu,norbd,nu_restricted,nd_restricted
      real(mp), intent(in) :: ef,kT
      real(mp), dimension(nkpts), intent(in) :: kwgts
      real(mp), dimension(norbu+norbd,nkpts), intent(in) :: eval
      real(mp), intent(out) :: electrons,dlectrons
      !local variables
      integer :: ikpt,iorb
      real(mp) :: f,df

      electrons=0.d0
      dlectrons=0.d0
      do ikpt=1,nkpts
         do iorb=1,norbd+norbu
            call get_occupation(occopt,eval(iorb,ikpt),ef,kT,f,df)
            if (iorb > norbu+nd_restricted .or. (iorb <= norbu .and. iorb > nu_restricted)) then
               f  = 0.d0
               df = 0.d0
            end if
            !call yaml_map('arg,f,kwgts(ikpt)',(/arg,f,kwgts(ikpt)/))
            electrons=electrons+ f  * kwgts(ikpt)  ! electrons := N_e(Ef+corr.)
            dlectrons=dlectrons+ df * kwgts(ikpt)  ! delectrons:= dN_e/darg ( Well! later we need dN_e/dEf=-1/kT*dN_e/darg
            !if(iproc==0) write(1000,*) iorb,arg,   f , df,dlectrons
         enddo
      enddo
      dlectrons=dlectrons/(-kT)  ! df/dEf=df/darg * -1/kT
    end subroutine electrons_and_delectrons

    subroutine find_fermi_level(allow_for_warning,occopt,nkpts,norbu,norbd,newnorbu,newnorbd,&
                 eval,kwgts,kT,charge,ef)
      use yaml_output
      implicit none
      logical, intent(in) :: allow_for_warning
      integer, intent(in) :: occopt,nkpts,norbu,norbd,newnorbu,newnorbd
      real(mp), intent(in) :: kT,charge
      real(mp), dimension(nkpts), intent(in) :: kwgts
      real(mp), dimension(norbu+norbd,nkpts), intent(in) :: eval
      real(mp), intent(inout) :: ef
      !local variables
      logical :: exitfermi
      integer :: ii
      real(mp) :: wf,electrons,dlectrons

      ! electrons is N_electons = sum f_i * Weight_i
      ! dlectrons is dN_electrons/dEf =dN_electrons/darg * darg/dEf= sum df_i/darg /(-wf) , darg/dEf=-1/wf
      !  f:= occupation # for band i ,  df:=df/darg
      wf=kT
      loop_fermi: do ii=1,100
         !write(1000+iproc,*) 'iteration',ii,' -------------------------------- '
         if (ii == 100 .and. allow_for_warning) &
              call yaml_warning('Fermilevel could not have been adjusted in the available iterations')
!!$            electrons=0.d0
!!$            dlectrons=0.d0
!!$            do ikpt=1,nkpts
!!$               do iorb=1,norbd+norbu
!!$                  call get_occupation(occopt,eval((ikpt-1)*norb+iorb),ef,wf,f,df)
!!$                  if (iorb > norbu+newnorbd .or. (iorb <= norbu .and. iorb > newnorbu)) then
!!$                     f  = 0.d0
!!$                     df = 0.d0
!!$                  end if
!!$                  !call yaml_map('arg,f,kwgts(ikpt)',(/arg,f,kwgts(ikpt)/))
!!$                  electrons=electrons+ f  * kwgts(ikpt)  ! electrons := N_e(Ef+corr.)
!!$                  dlectrons=dlectrons+ df * kwgts(ikpt)  ! delectrons:= dN_e/darg ( Well! later we need dN_e/dEf=-1/wf*dN_e/darg
!!$                  !if(iproc==0) write(1000,*) iorb,arg,   f , df,dlectrons
!!$               enddo
!!$            enddo
!!$
!!$            !call yaml_map('ef',ef)
!!$            !call yaml_map('electrons',electrons)
!!$
!!$            dlectrons=dlectrons/(-wf)  ! df/dEf=df/darg * -1/wf

         call electrons_and_delectrons(occopt,nkpts,norbu,norbd,newnorbu,newnorbd,&
              eval,kwgts,ef,wf,&
              electrons,dlectrons)
         call correct_fermi(ii <= 50,electrons,dlectrons,charge,kT,wf,ef,exitfermi)

!!$            dlectrons=dlectrons/(-wf)  ! df/dEf=df/darg * -1/wf
!!$            diff=-charge/full+electrons
!!$            !if (iproc.lt.1) write(1000+iproc,*) diff,full,melec,real(melec,gp)
!!$            !         if (iproc.lt.1) flush(1000+iproc)
!!$            !if (iproc.lt.1) write(1000+iproc,*) diff,1.d-11*sqrt(electrons),wf
!!$            !if (iproc.lt.1) flush(1000+iproc)
!!$            !Exit criterion satiesfied, Nevertheles do one mor update of fermi level
!!$            if (abs(diff) < 1.d-11*sqrt(electrons) .and. wf == wf0 ) exitfermi=.true.     ! Assume noise grows as sqrt(electrons)
!!$
!!$            !alternative solution to avoid division by so high value
!!$            !if (dlectrons == 0.d0) dlectrons=1.d-100  !line to be added
!!$            if (dlectrons == 0.d0) then
!!$               !always enter into first case below
!!$               corr=0.d0
!!$               if (diff > 0.d0) corr=1.d0*wf
!!$               if (diff < 0.d0) corr=-1.d0*wf
!!$               if (ii <= 50 .and. wf < 0.1d0) wf=2.d0*wf  ! speed up search of approximate Fermi level by using higher Temperature
!!$            else
!!$               corr=diff/abs(dlectrons) ! for case of no-monotonic func. abs is needed
!!$               if (abs(corr).gt.wf) then   !for such a large correction the linear approximation is not any more valid
!!$                  if (corr > 0.d0) corr=1.d0*wf
!!$                  if (corr < 0.d0*wf) corr=-1.d0*wf
!!$                  if (ii <= 50 .and. wf < 0.1d0) wf=2.d0*wf  ! speed up search of approximate Fermi level by using higher Temperature
!!$               else
!!$                  wf=max(wf0,.5d0*wf)
!!$               endif
!!$            end if
!!$            ef=ef-corr  ! Ef=Ef_guess+corr.
!!$            !if (iproc.lt.1) write(1000+iproc,'(i5,5(1pe17.8))') ii,electrons,ef,dlectrons,abs(dlectrons),corr
!!$            !if (iproc.lt.1) flush(1000+iproc)
!!$            !call determine_fermi_level(ft, electrons, ef,info_fermi)
!!$            !if (info_fermi /= 0) then
!!$            !   call f_err_throw('Difficulties in guessing the new Fermi energy, info='//trim(yaml_toa(info_fermi)),&
!!$            !        err_name='BIGDFT_RUNTIME_ERROR')
!!$            !end if

            !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr) !debug
            if (exitfermi) exit loop_fermi
         end do loop_fermi
       end subroutine find_fermi_level



end module fermi_level
