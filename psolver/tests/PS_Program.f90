!> @file
!!  Program test for Poisson
!!  Laplacian V = 4pi rho
!!  May work either in parallel or in serial case
!!  And for different geometries
!! @author
!!    Copyright (C) 2006-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test program for the Poisson Solver
program PSolver_Program
  use Poisson_Solver
  use wrapper_mpi
  use futile
  use numerics
  use PSbox
  use box
  use f_harmonics
  use f_blas
  use IObox
  use at_domain
  implicit none
  !Order of interpolating scaling function
  !integer, parameter :: itype_scf=8
  character(len=*), parameter :: subname='Poisson_Solver'
  real(kind=8), parameter :: a_gauss = 1.0d0, a2 = a_gauss**2
  !Length of the box
  real(kind=8), parameter :: acell = 10.0d0
  real(kind=8), parameter :: EulerGamma = 0.5772156649015328d0

  !Type of function
  integer, parameter :: FUNC_CONSTANT = 1
  integer, parameter :: FUNC_GAUSSIAN = 2
  integer, parameter :: FUNC_GAUSSIAN_SHRINKED = 3
  integer, parameter :: FUNC_COSINE = 4
  integer, parameter :: FUNC_EXP_COSINE = 5
  integer, parameter :: FUNC_SHRINK_GAUSSIAN = 6
  integer, parameter :: FUNC_SINE = 7
  integer, parameter :: FUNC_ATAN = 8

  character(len=*), parameter :: inputs=&
       "- {name: ndim, shortname: n, default: 30,"//&
       "  help_string: Size of the simulation domain,"//&
       "  help_dict: {Allowed values: list of integers}}"//f_cr//&

       "- {name: geocode,shortname: g, default: P,"//&
       "  help_string: Boundary conditions,"//&
       "  help_dict: {Usage: set the boundary conditions of the run,"//&
       "              Allowed values: [F, S , W, P]}}"//f_cr//&

       "- {name: angdeg, shortname: d, default: 90.0,"//&
       "  help_string: Degrees of the angles between the directions,"//&
       "  help_dict: {Allowed values: arrays of floats}}"//f_cr//&

       "- {name: input, shortname: i, default: None,"//&
       "  help_string: Inpufile of Poisson Solver,"//&
       "  help_dict: {Allowed values: dictionary in yaml format (mapping)}}"


  character(len=1) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
  character(len=30) :: mode,info
  real(kind=8), dimension(:,:,:), allocatable :: density,rhopot,potential,pot_ion
  type(coulomb_operator) :: karray
  integer, dimension(3) :: ndims
  real(f_double), dimension(3) :: angdeg
  type(dictionary), pointer :: dict
  real(kind=8) :: hx,hy,hz,hgrid,offset
  real(kind=8) :: ehartree,eexcu,diff_par,e1,ehartree_exp
  integer :: n01,n02,n03
  integer :: i1,i2,i3,j1,j2,j3,i1_max,i2_max,i3_max,iproc,nproc,i3sd,ncomp
  integer :: n_cell,ixc,i3s,unit
  !integer :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
  !triclinic lattice
  real(kind=8) :: alpha,beta,gamma
  real(kind=8), dimension(:,:,:,:), pointer :: rhocore_fake
  type(dictionary), pointer :: options,input
  type(f_multipoles) :: multipoles
  type(box_iterator) :: bit
  external :: gather_timings
  real(dp), parameter :: tol=1.e-8
  logical, parameter :: wrtfiles=.true. !.false.
  real(kind=8), dimension(3) :: hgrids
  type(domain) :: dom
  nullify(rhocore_fake)

  !mode = "charged_thin_wire"
  !mode="cylindrical_capacitor"
  !mode="monopolar"
  mode="zigzag_model_wire"


  call f_lib_initialize()

  call mpiinit()
  iproc=mpirank()
  nproc=mpisize()
  call f_malloc_set_status(iproc=iproc)

  nullify(dict,input)
  call yaml_argparse(options,inputs)
  if (iproc==0) then
     call yaml_new_document()
     call yaml_map('Commandline options provided',options)
  end if
  ndims=options//'ndim'
  n01=ndims(1)
  n02=ndims(2)
  n03=ndims(3)
  ixc=0 !not needed anymore
  geocode=options//'geocode'
  angdeg=options//'angdeg'
  input=options .get. 'input'
  call dict_copy(dict,input) !null if absent

  call dict_free(options)

  alpha = angdeg(1)/180.0_f_double*pi!2.0_dp*datan(1.0_dp) !to be modified
  beta  = angdeg(2)/180.0_f_double*pi!2.0_dp*datan(1.0_dp)
  gamma = angdeg(3)/180.0_f_double*pi!2.0_dp*datan(1.0_dp)
  !alpha = 1.0_dp*datan(1.0_dp)
  !beta  = 1.0_dp*datan(1.0_dp)
  !gamma = 1.0_dp*datan(1.0_dp)

  !perform also the comparison with the serial case
  !code for the Poisson Solver in the parallel case
  select case(geocode)
  case('P')
     info='periodic BC'
  case('S')
     info='surface BC'
  case('F')
     info='free BC'
  case('W')
     info='wires BC'
  end select

  if (iproc==0) then
     call yaml_map('PSolver, '//trim(info),ndims)
     call yaml_map('processes',nproc)
  end if

  !initialize memory counting
  !call memocc(0,iproc,'count','start')

  !Step size
  n_cell = max(n01,n02,n03)
  hx=acell/real(n01,kind=8)
  hy=acell/real(n02,kind=8)
  hz=acell/real(n03,kind=8)


  !grid for the free BC case
  hgrid=max(hx,hy,hz)
  !hgrid=hx
  hgrids=(/hx,hy,hz/) 

  dom=domain_new(units=ATOMIC_UNITS,bc=geocode_to_bc_enum(geocode),&
            alpha_bc=alpha,beta_ac=beta,gamma_ab=gamma,acell=ndims*hgrids)

  !we must choose properly a test case with a positive density
  !itype_scf=16

  !write(*,'(a12,i4)') ' itype_scf = ', itype_scf
  !call yaml_map('itype_scf',itype_scf)

  call f_timing_reset(filename='time.yaml',master=iproc==0)
  !call timing(nproc,'time.prc','IN')

!  karray=pkernel_init(iproc,nproc,dict,&
!       geocode,ndims,(/hx,hy,hz/),&
!       alpha_bc=alpha,beta_ac=beta,gamma_ab=gamma)

  karray=pkernel_init(iproc,nproc,dict,&
       dom,ndims,(/hx,hy,hz/),&
       alpha_bc=alpha,beta_ac=beta,gamma_ab=gamma)

  call pkernel_set(karray,verbose=.true.)

  !Allocations
  !Density
  density = f_malloc(ndims,id='density')
  !Density then potential
  rhopot = f_malloc(ndims,id='rhopot')
  potential = f_malloc(ndims,id='potential')
  !ionic potential
  pot_ion = f_malloc(ndims,id='pot_ion')

  if (iproc==0) call yaml_map('ixc',ixc)

  call test_functions_new2(karray%mesh,acell,a_gauss,karray%mu,density,potential)
  !calculate the expected hartree energy
  ehartree_exp=0.5_dp*f_dot(potential,density)*karray%mesh%volume_element

  call f_multipoles_create(multipoles,lmax=2)
  bit=box_iter(karray%mesh,centered=.true.)

  call field_multipoles(bit,density,1,multipoles)

  if (iproc==0) then
     call yaml_map('Expected Eh',ehartree_exp)
     call yaml_map('monopole',get_monopole(multipoles),fmt='(1pe15.7)')
     call yaml_map('dipole',get_dipole(multipoles),fmt='(1pe15.7)')
  end if
  call f_multipoles_release(multipoles)

  call f_multipoles_create(multipoles,lmax=0)
  call field_multipoles(bit,potential,1,multipoles)
  offset=get_monopole(multipoles)
  if (iproc==0) then
     call yaml_map('Offset (potential monopole)',offset)
  end if
  call f_multipoles_release(multipoles)

  call f_memcpy(src=density,dest=rhopot)

  if (wrtfiles) then
   i2=n02/2
   do i3=1,n03
     do i1=1,n01
        j1=n01/2+1-abs(n01/2+1-i1)
        j2=n02/2+1-abs(n02/2+1-i2)
        j3=n03/2+1-abs(n03/2+1-i3)
        write(110,'(2(1x,I8),2(1x,e22.15))')i1,i3,rhopot(i1,i2,i3),potential(i1,i2,i3)
     end do
   end do
  end if

!!$  !offset, used only for the periodic solver case
!!$  offset=0.0_gp!potential(1,1,1)!-pot_ion(1,1,1)
!!$  do i3=1,n03
!!$     do i2=1,n02
!!$        do i1=1,n01
!!$           offset=offset+potential(i1,i2,i3)
!!$        end do
!!$     end do
!!$  end do
!!$  offset=offset*hx*hy*hz*sqrt(detg) ! /// to be fixed ///
!!$  !write(*,*) 'offset = ',offset
!!$  if (iproc==0) call yaml_map('Offset (potential monopole)',offset)
  call PS_set_options(karray,potential_integral=offset)

  ncomp=n03
  if (karray%opt%datacode=='G') then
     i3sd=1
  else
     i3sd=karray%grid%istart+1
     !ncomp=karray%grid%n3p
  end if
  i3s=i3sd

  !apply the Poisson Solver (case with distributed potential)
  call Electrostatic_Solver(karray,density(1,1,i3sd),ehartree=ehartree)

  call PS_gather(density,karray)
  !this has to be corrected with the volume element of mesh
  eexcu=sum(density)*karray%mesh%volume_element
  if (iproc==0) call yaml_map('potential integral',eexcu)

  if (wrtfiles) then
   i3=n03/2
   do i2=1,n02
     do i1=1,n01
        !j1=n01/2+1-abs(n01/2+1-i1)
        !j2=n02/2+1-abs(n02/2+1-i2)
        !j3=n03/2+1-abs(n03/2+1-i3)
        write(111,'(2(1x,i6),3(1x,1pe25.16e3))') i1,i2,rhopot(i1,i2,i3),potential(i1,i2,i3),density(i1,i2,i3)
        !write(111,*) i1*hx+hy*i2*dcos(alpha)+i3*hz*dcos(beta), &
        !     i2*hy*dsin(alpha)+i3*hz*(-dcos(alpha)*dcos(beta)+dcos(gamma))/dsin(alpha), &
        !     rhopot(i1,i2,i3),potential(i1,i2,i3), &
        !     density(i1,i2,i3)
     end do
   end do

   i2=n02/2
   do i3=1,n03
     !do i2=1,n02
     do i1=1,n01
        !j1=n01/2+1-abs(n01/2+1-i1)
        !j2=n02/2+1-abs(n02/2+1-i2)
        !j3=n03/2+1-abs(n03/2+1-i3)
        write(112,'(2(1x,i6),3(1x,1pe25.16e3))')i1,i3,rhopot(i1,i2,i3),potential(i1,i2,i3),density(i1,i2,i3)
     end do
     !end do
   end do
  end if

  call f_timing_stop(mpi_comm=karray%mpi_env%mpi_comm,nproc=karray%mpi_env%nproc,gather_routine=gather_timings)
  call pkernel_free(karray)

  call compare(n01,n02,n03,potential,density,i1_max,i2_max,i3_max,diff_par)

  if (iproc == 0) then
     call yaml_mapping_open('Report on comparison')
     call yaml_map('Ehartree',ehartree)
     call yaml_map('Ehartree diff',ehartree-ehartree_exp)
     call yaml_map('Max diff at',[i1_max,i2_max,i3_max])
     call yaml_get_default_stream(unit)
!     if (f_tty(unit) .and. abs(diff_par) > tol) then
!        call yaml_map('Max diff',yaml_blink(trim(yaml_toa(diff_par))))
!     else
        call yaml_map('Max diff',diff_par)
!     end if
     call yaml_map('Result',density(i1_max,i2_max,i3_max))
     call yaml_map('Original',potential(i1_max,i2_max,i3_max))
     call yaml_mapping_close()
  end if

  call dict_free(dict)
  call f_free(density)
  call f_free(rhopot)
  call f_free(potential)
  call f_free(pot_ion)

  call mpifinalize()
  call f_lib_finalize()

contains

!> This subroutine builds some analytic functions that can be used for
!! testing the poisson solver.
!! The default choice is already well-tuned for comparison.
!! WARNING: not all the test functions can be used for all the boundary conditions of
!! the poisson solver, in order to have a reliable analytic comparison.
!! The parameters of the functions must be adjusted in order to have a sufficiently localized
!! function in the isolated direction and an explicitly periodic function in the periodic ones.
!! Beware of the high-frequency components that may falsify the results when hgrid is too high.
subroutine test_functions(geocode,ixc,n01,n02,n03,acell,a_gauss,hx,hy,hz,&
     density,potential,rhopot,pot_ion,mu0,alpha,beta,gamma)
  use yaml_output
  use f_utils
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
  integer, intent(in) :: n01,n02,n03,ixc
  real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz,mu0
  !triclinic lattice
  real(kind=8), intent(in) :: alpha, beta, gamma
  real(kind=8), dimension(n01,n02,n03), intent(out) :: density,potential,rhopot,pot_ion

  !local variables
  integer :: i1,i2,i3,ifx,ify,ifz,unit
  real(kind=8) :: x,x1,x2,x3,y,z,length,denval,a2,derf,factor,r,r2,r0
  real(kind=8) :: fx,fx1,fx2,fy,fy1,fy2,fz,fz1,fz2,a,ax,ay,az,bx,by,bz,tt,potion_fac
  real(kind=8) :: monopole
  real(kind=8), dimension(3) :: dipole

  !non-orthorhombic lattice
  real(kind=8), dimension(3,3) :: gu,gd
  real(kind=8) :: detg

  !triclinic cell
  !covariant metric
  gd(1,1) = 1.0_dp
  gd(1,2) = dcos(alpha)
  gd(1,3) = dcos(beta)
  gd(2,2) = 1.0_dp
  gd(2,3) = dcos(gamma)
  gd(3,3) = 1.0_dp

  gd(2,1) = gd(1,2)
  gd(3,1) = gd(1,3)
  gd(3,2) = gd(2,3)
  !
  detg = 1.0_dp - dcos(alpha)**2 - dcos(beta)**2 - dcos(gamma)**2 + 2.0_dp*dcos(alpha)*dcos(beta)*dcos(gamma)

  !write(*,*) 'detg =', detg
  if (iproc==0) call yaml_map('detg',detg)
  !
  !contravariant metric
  gu(1,1) = (dsin(gamma)**2)/detg
  gu(1,2) = (dcos(beta)*dcos(gamma)-dcos(alpha))/detg
  gu(1,3) = (dcos(alpha)*dcos(gamma)-dcos(beta))/detg
  gu(2,2) = (dsin(beta)**2)/detg
  gu(2,3) = (dcos(alpha)*dcos(beta)-dcos(gamma))/detg
  gu(3,3) = (dsin(alpha)**2)/detg
  !
  gu(2,1) = gu(1,2)
  gu(3,1) = gu(1,3)
  gu(3,2) = gu(2,3)

  !gu=gd !test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (iproc==0) then
     call yaml_map('Angles',[alpha,beta,gamma]*180.0_dp*oneopi)
     call yaml_map('Contravariant Metric',gu)
     call yaml_map('Covariant Metric',gd)
     call yaml_map('Product of the two',matmul(gu,gd))
  end if

  unit=200
  call f_open_file(unit=unit,file='references.dat')

  if (ixc==0) denval=0.d0

  if (trim(geocode) == 'P') then
     !parameters for the test functions
     length=acell
     a=0.5d0/a_gauss**2
     !test functions in the three directions
     ifx=FUNC_SHRINK_GAUSSIAN
     ify=FUNC_SHRINK_GAUSSIAN
     ifz=FUNC_SHRINK_GAUSSIAN
     !parameters of the test functions
     ax=length
     ay=length
     az=length

     !the following b's are not used, actually
     bx=2.d0!real(nu,kind=8)
     by=2.d0!real(nu,kind=8)
     bz=2.d0


     ! !original version

     if (wrtfiles) then
      !plot of the functions used
      do i1=1,n03
        x = hx*real(i1-n01/2-1,kind=8)!valid if hy=hz
        y = hz*real(i1-n03/2-1,kind=8)
        call functions(x,ax,bx,fx,fx1,fx2,ifx)
        call functions(y,az,bz,fz,fz1,fz2,ifz)
        write(20,'(1x,I8,4(1x,e22.15))') i1,fx,fx2,fz,fz2
      end do
     end if
     !Initialization of density and potential
     denval=0.d0 !value for keeping the density positive
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,bz,fz,fz1,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,by,fy,fy1,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n01/2-1,kind=8)
              call functions(x1,ax,bx,fx,fx1,fx2,ifx)
              potential(i1,i2,i3) =  -16.d0*datan(1.d0)*fx*fy*fz
              !density(i1,i2,i3) = fx2*fy*fz+fx*fy2*fz+fx*fy*fz2-mu0**2*fx*fy*fz
              !triclinic lattice
              density(i1,i2,i3) = -mu0**2*fx*fy*fz
              density(i1,i2,i3) = density(i1,i2,i3) + gu(1,1)*fx2*fy*fz+gu(2,2)*fx*fy2*fz+gu(3,3)*fx*fy*fz2
              density(i1,i2,i3) = density(i1,i2,i3) + 2.0_dp*(gu(1,2)*fx1*fy1*fz+gu(1,3)*fx1*fy*fz1+gu(2,3)*fx*fy1*fz1)
              denval=max(denval,-density(i1,i2,i3))
           end do
        end do
     end do


     ! !tweaked version: for debugging the solver for non-orthorhombic cells

     ! pi = 4.d0*atan(1.d0)
     ! a2 = a_gauss**2/2
     ! !mu0 = 1.e0_dp

     ! !Normalization
     ! !factor = a_gauss*sqrt(pi)/2.0_dp
     ! factor = 2.0_dp
     ! !gaussian function
     ! do i3=1,n03
     !    !x3 = hz*real(i3-n03/2,kind=8)
     !    do i2=1,n02
     !       !x2 = hy*real(i2-n02/2,kind=8)
     !       do i1=1,n01
     !          x1 = hx*real(i1-n01/2,kind=8)+hy*real(i2-n02/2,kind=8)*dcos(alpha)+hz*real(i3-n03/2,kind=8)*dcos(beta)
     !          x2 = hy*real(i2-n02/2,kind=8)*dsin(alpha) + &
     !               & hz*real(i3-n03/2,kind=8)*(-dcos(alpha)*dcos(beta)+dcos(gamma))/dsin(alpha)
     !          x3 = hz*real(i3-n03/2,kind=8)*sqrt(detg)/dsin(alpha)
     !          !r2 = x1*x1+x2*x2+x3*x3
     !          !triclinic lattice:
     !          !r2 = gd(1,1)*x1*x1+gd(2,2)*x2*x2+gd(3,3)*x3*x3+2.0_dp*(gd(1,2)*x1*x2+gd(1,3)*x1*x3+gd(2,3)*x2*x3)
     !          r2 = x1*x1+x2*x2+x3*x3
     !          !density(i1,i2,i3) = factor*exp(-r2/a2)
     !          r = sqrt(r2)
     !          !Potential from a gaussian
     !          potential(i1,i2,i3) = dexp(-r2/a2)
     !          density(i1,i2,i3) = 4.0_dp*r2/a2**2-6.0_dp/a2
     !          density(i1,i2,i3) = -potential(i1,i2,i3)*density(i1,i2,i3)/16.0_dp/datan(1.0_dp)
     !       end do
     !    end do
     ! end do



    if (wrtfiles) then
     i2=n02/2
     do i3=1,n03
        do i1=1,n01
           !j1=n01/2+1-abs(n01/2+1-i1)
           !j2=n02/2+1-abs(n02/2+1-i2)
           !j3=n03/2+1-abs(n03/2+1-i3)
           write(unit,'(2(1x,I8),2(1x,e22.15))') i1,i3,density(i1,i2,i3),potential(i1,i2,i3)
        end do
     end do
    end if
!plane capacitor oriented along the y direction
!!     do i2=1,n02
!!        if (i2==n02/4) then
!!           do i3=1,n03
!!              do i1=1,n01
!!                 density(i1,i2,i3)=1.d0!real(i2,kind=8)
!!              end do
!!           end do
!!        else if (i2==3*n02/4) then
!!           do i3=1,n03
!!              do i1=1,n01
!!                 density(i1,i2,i3)=-1.d0!real(i2,kind=8)
!!              end do
!!           end do
!!        else
!!           do i3=1,n03
!!              do i1=1,n01
!!                 density(i1,i2,i3)=0.d0
!!              end do
!!           end do
!!        end if
!!     end do
!!     denval=0.d0

     if (ixc==0) denval=0.d0



  else if (trim(geocode) == 'S') then
     !parameters for the test functions
     length=acell
     a=0.5d0/a_gauss**2
     !test functions in the three directions
     ifx=FUNC_EXP_COSINE
     ifz=FUNC_EXP_COSINE !FUNC_CONSTANT
     !non-periodic dimension
     ify=FUNC_SHRINK_GAUSSIAN
     !parameters of the test functions
     ax=length
     ay=length
     az=length
     !b's are not used, actually
     bx=2.d0!real(nu,kind=8)
     by=2.d0!real(nu,kind=8)
     bz=2.d0!real(nu,kind=8)

     call f_assert(alpha-onehalf*pi,id='Alpha angle invalid')
     call f_assert(gamma-onehalf*pi,id='Gamma angle invalid for S BC')

     !non-periodic dimension
     !ay=length!4.d0*a

     density(:,:,:) = 0.d0!1d-20 !added

     if (wrtfiles) then
      !plot of the functions used
      do i1=1,n02
        x = hx*real(i1-n02/2-1,kind=8)!valid if hy=hz
        y = hy*real(i1-n02/2-1,kind=8)
        z = hz*real(i1-n02/2-1,kind=8)
        call functions(x,ax,bx,fx,fx1,fx2,ifx)
        call functions(y,ay,by,fy,fy1,fy2,ify)
        call functions(z,az,bz,fz,fz1,fz2,ifz)
        write(20,'(1x,I8,6(1x,e22.15))') i1,fx,fx2,fy,fy2,fz,fz2
      end do
     end if
     !Initialisation of density and potential
     !Normalisation
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,bz,fz,fz1,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,by,fy,fy1,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n01/2-1,kind=8)
              call functions(x1,ax,bx,fx,fx1,fx2,ifx)
              potential(i1,i2,i3) =  -fourpi*fx*fy*fz
              density(i1,i2,i3) = -mu0**2*fx*fy*fz
              density(i1,i2,i3) = density(i1,i2,i3) + gu(1,1)*fx2*fy*fz+gu(2,2)*fx*fy2*fz+gu(3,3)*fx*fy*fz2
              density(i1,i2,i3) = density(i1,i2,i3) + 2.0_dp*(gu(1,2)*fx1*fy1*fz+gu(1,3)*fx1*fy*fz1+gu(2,3)*fx*fy1*fz1)
              !old:
              !density(i1,i2,i3) = fx2*fy*fz+fx*fy2*fz+fx*fy*fz2 - mu0**2*fx*fy*fz
              denval=max(denval,-density(i1,i2,i3))
           end do
        end do
     end do


     ! !plane capacitor oriented along the y direction
     ! do i2=1,n02
     !    if (i2==n02/4) then
     !       do i3=1,n03
     !          do i1=1,n01
     !             density(i1,i2,i3)=1.d0!real(i2,kind=8)
     !          end do
     !       end do
     !    else if (i2==3*n02/4) then
     !       do i3=1,n03
     !          do i1=1,n01
     !             density(i1,i2,i3)=-1.d0!real(i2,kind=8)
     !          end do
     !       end do
     !    else
     !       do i3=1,n03
     !          do i1=1,n01
     !             density(i1,i2,i3)=0.d0
     !          end do
     !       end do
     !    end if
     ! end do



    if (wrtfiles) then
     i2=n02/2
     do i3=1,n03
        do i1=1,n01
           !j1=n01/2+1-abs(n01/2+1-i1)
           !j2=n02/2+1-abs(n02/2+1-i2)
           !j3=n03/2+1-abs(n03/2+1-i3)
           z=real(i3,dp)*sin(beta)
           write(unit,'(2(1x,i6),3(1x,1pe26.14e3))') &
                i1,i3,z,density(i1,i2,i3),potential(i1,i2,i3)
        end do
     end do
    end if

     if (ixc==0) denval=0.d0


   else if (trim(geocode) == 'F') then

      !grid for the free BC case
      !hgrid=max(hx,hy,hz)

!      pi = 4.d0*atan(1.d0)
      a2 = a_gauss**2

      !Normalization
      factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
      !gaussian function
      do i3=1,n03
         x3 = hz*real(i3-n03/2,kind=8)
         do i2=1,n02
            x2 = hy*real(i2-n02/2,kind=8)
            do i1=1,n01
               x1 = hx*real(i1-n01/2,kind=8)
               r2 = x1*x1+x2*x2+x3*x3
               density(i1,i2,i3) = factor*exp(-r2/a2)
               r = sqrt(r2)
               !Potential from a gaussian
               if (r == 0.d0) then
                  potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
               else
                  potential(i1,i2,i3) = derf(r/a_gauss)/r
               end if
            end do
         end do
      end do

      if (wrtfiles) then
       i2=n02/2
       do i3=1,n03
         do i1=1,n01
            !j1=n01/2+1-abs(n01/2+1-i1)
            !j2=n02/2+1-abs(n02/2+1-i2)
            !j3=n03/2+1-abs(n03/2+1-i3)
            write(200,*) i1,i3,density(i1,i2,i3),potential(i1,i2,i3)
         end do
       end do
      end if

! !plane capacitor oriented along the y direction
! !!     do i2=1,n02
! !!        if (i2==n02/4) then
! !!           do i3=1,n03
! !!              do i1=1,n01
! !!                 density(i1,i2,i3)=1.d0!real(i2,kind=8)
! !!              end do
! !!           end do
! !!        else if (i2==3*n02/4) then
! !!           do i3=1,n03
! !!              do i1=1,n01
! !!                 density(i1,i2,i3)=-1.d0!real(i2,kind=8)
! !!              end do
! !!           end do
! !!        else
! !!           do i3=1,n03
! !!              do i1=1,n01
! !!                 density(i1,i2,i3)=0.d0
! !!              end do
! !!           end do
! !!        end if
! !!     end do

      denval=0.d0


!!!>  else if (trim(geocode) == 'H' .or. trim(geocode) == 'F') then
!!!>
!!!>     !hgrid=max(hx,hy,hz)
!!!>
!!!>     a2 = a_gauss**2
!!!>     !mu0 = 1.e0_dp
!!!>
!!!>     !Normalization
!!!>     !factor = a_gauss*sqrt(pi)/2.0_dp
!!!>     factor = 2.0_dp-a_gauss*dexp(a2*mu0**2/4.0_dp)*sqrt(pi)*mu0*derfc(mu0*a_gauss/2.0_dp)
!!!>     !gaussian function
!!!>     do i3=1,n03
!!!>        x3 = hz*real(i3-n03/2,kind=8)
!!!>        do i2=1,n02
!!!>           x2 = hy*real(i2-n02/2,kind=8)
!!!>           do i1=1,n01
!!!>              x1 = hx*real(i1-n01/2,kind=8)
!!!>              !r2 = x1*x1+x2*x2+x3*x3
!!!>              !triclinic lattice:
!!!>              r2 = x1*x1+x2*x2+x3*x3
!!!>              !density(i1,i2,i3) = factor*exp(-r2/a2)
!!!>              r = sqrt(r2)
!!!>              !Potential from a gaussian
!!!>              if (r == 0.d0) then
!!!>                 potential(i1,i2,i3) = 1.0_dp
!!!>              else
!!!>                 call derf_local(erf_yy,r/a_gauss-a_gauss*mu0/2.0_dp)
!!!>                 erfc_yy=1.0_dp-erf_yy
!!!>                 potential(i1,i2,i3) = -2.0_dp+erfc_yy
!!!>                 !potential(i1,i2,i3) = -2+derfc(r/a_gauss-a_gauss*mu0/2.0_dp)
!!!>                 call derf_local(erf_yy,r/a_gauss+a_gauss*mu0/2.0_dp)
!!!>                 erfc_yy=1.0_dp-erf_yy
!!!>                 potential(i1,i2,i3) = potential(i1,i2,i3)+dexp(2.0_dp*r*mu0)*erfc_yy
!!!>                 !potential(i1,i2,i3) = potential(i1,i2,i3)+dexp(2.0_dp*r*mu0)*derfc(r/a_gauss+a_gauss*mu0/2.0_dp)
!!!>                 potential(i1,i2,i3) = potential(i1,i2,i3)*a_gauss*dexp(-mu0*r)*sqrt(pi)/(-2*r*factor*dexp(-a2*mu0**2/4.0_dp))
!!!>              end if
!!!>              !density(i1,i2,i3) = exp(-r2/a2)/4.0_dp/factor**2 + 0.1_dp**2/(4*pi)*potential(i1,i2,i3)
!!!>              density(i1,i2,i3) = safe_exp(-r2/a2)/factor/(a2*pi)
!!!>           end do
!!!>        end do
!!!>     end do
!!!>
!!!>     i2=n02/2
!!!>     do i3=1,n03
!!!>        do i1=1,n01
!!!>           !j1=n01/2+1-abs(n01/2+1-i1)
!!!>           !j2=n02/2+1-abs(n02/2+1-i2)
!!!>           !j3=n03/2+1-abs(n03/2+1-i3)
!!!>           write(unit,*) i1,i3,density(i1,i2,i3),potential(i1,i2,i3)
!!!>        end do
!!!>     end do
!!!>
!!!>
!!!>     denval=0.d0
!!!>

  else if (trim(geocode) == 'W') then
     !parameters for the test functions
     length=acell
     !a=0.5d0/a_gauss**2
     !test functions in the three directions
     !isolated directions
     ifx=FUNC_SHRINK_GAUSSIAN
     ify=FUNC_SHRINK_GAUSSIAN
     !periodic direction
     ifz=5
     !parameters of the test functions

     ax = length
     ay = length
     az = length

     bx = 2.d0
     by = 2.d0
     bz = 2.d0


     density(:,:,:) = 0.d0!1d-20 !added
     factor = 2.0d0


     if (wrtfiles) then
      !plot of the functions used
      do i1=1,min(n01,n03)
        x = hx*real(i1-n01/2-1,kind=8)!isolated
        z = hz*real(i1-n03/2-1,kind=8)!periodic
        call functions(x,ax,bx,fx,fx1,fx2,ifx)
        call functions(z,az,bz,fz,fz1,fz2,ifz)
        write(20,*) i1,fx,fx2,fz,fz2
      end do
     end if

     !Initialization of density and potential


     select case(mode)
     case ("monopolar")
        ! Gaussian Density Distribution in (x,y)
        ! this is the configuration yielding non-zero monopole
        do i3=1,n03
           x3 = hz*real(i3-n03/2-1,kind=8)
           call functions(x3,az,bz,fz,fz1,fz2,1)
           do i2=1,n02
              x2 = hy*real(i2-n02/2-1,kind=8)
              !call functions(x2,ay,by,fy,fy1,fy2,ify)
              do i1=1,n01
                 x1 = hx*real(i1-n01/2-1,kind=8)
                 r2 = x1*x1+x2*x2
                 r = sqrt(r2)
                 !call functions(x1,ax,bx,fx,fx1,fx2,ifx)
                 if  (r == 0.d0) then
                    !EulerGamma = 0.5772156649015328d
                    density(i1,i2,i3) = dexp(-factor*r2)
                    potential(i1,i2,i3) = (-EulerGamma - dlog(factor))/(4.0d0*factor)
                 else
                    call e1xb(factor*r2,e1)
                    density(i1,i2,i3) = dexp(-factor*r2)
                    potential(i1,i2,i3) = (e1+dlog(r2))/(4.0d0*factor)
                 end if
                 !note that in this case we cannot account for the screening in the following way,
                 !density(i1,i2,i3) = density(i1,i2,i3) + mu0**2*potential(i1,i2,i3)
                 !because the extra-term, proportional to the potential, is not localized
                 !in the non-periodic directions
                 potential(i1,i2,i3) = -16.0d0*datan(1.0d0)*potential(i1,i2,i3)
              end do
           end do
        end do

     case("zigzag_model_wire")

        density = 0.d0
        potential = 0.d0

!!$        density(n01/4,n02/2,n03/4) = -1.0d0
!!$        density(3*n01/4,n02/2,3*n03/4) = 1.0d0
        density(n01/2,n02/2,n03/2) = 1.0d0

        !factor=(16.d0/acell)**2
        !r0 = acell/4.d0
        !the following is evaluated analytically by imposing that
        !\int_0^\infty r*(-exp(-factor(r-r0)^2)+denval*exp(-factor*r^2)) = 0
        !denval=sqrt(4.d0*datan(1.d0)*factor)*r0*(1.d0+derf(sqrt(factor)*r0))
        !do i3=1,n03
        ! x3 = hz*real(i3-n03/2-1,kind=8)
        ! do i2=1,n02
        ! x2 = hy*real(i2-n02/2-1,kind=8)
        ! do i1=1,n01
        ! x1 = hx*real(i1-n01/2-1,kind=8)
        ! !r2 = x1*x1+x2*x2+x3*x3
        ! !r = sqrt(r2)
        ! !in this configuration denval is used so as to achieve zero monopole
        ! density(i1,i2,i3) = -1.d0*dexp(-factor*(x1-r0)**2)*dexp(-factor*x2**2)*dexp(-factor*(x3-r0)**2) &
        ! + 1.0d0*dexp(-factor*(x1+r0)**2)*dexp(-factor*x2**2)*dexp(-factor*(x3+r0)**2)
        ! density(i1,i2,i3) = density(i1,i2,i3)*(factor/4.d0/datan(1.d0))**(3.d0/2.d0)
        ! end do
        ! end do
        !end do

     case ("charged_thin_wire")
        do i3=1,n03
           do i2=1,n02
              do i1=1,n01
                 if (i1 == n01/2+1 .and. i2 == n02/2+1) density(i1,i2,i3) = 1.0d0
              end do
           end do
        end do
     case("cylindrical_capacitor")
        !mimicked by two Gaussian charge distributions,
        !one localized around r = 0,
        !the other around r0 =acell/4
        factor=3.d0*acell
        r0 = acell/4.d0
        !the following is evaluated analytically by imposing that
        !\int_0^\infty r*(-exp(-factor(r-r0)^2)+denval*exp(-factor*r^2)) = 0
        denval=sqrt(4.d0*datan(1.d0)*factor)*r0*(1.d0+derf(sqrt(factor)*r0))
        do i3=1,n03
           x3 = hz*real(i3-n03/2-1,kind=8)
           do i2=1,n02
              x2 = hy*real(i2-n02/2-1,kind=8)
              do i1=1,n01
                 x1 = hx*real(i1-n01/2-1,kind=8)
                 r2 = x1*x1+x2*x2
                 r = sqrt(r2)
                 !in this configuration denval is used so as to achieve zero monopole
                 density(i1,i2,i3) = density(i1,i2,i3) + denval*dexp(-factor*r2) - dexp(-factor*(r-r0)**2)
              end do
           end do
        end do
     case default
        denval=0.d0 !value for keeping the density positive
        do i3=1,n03
           x3 = hz*real(i3-n03/2-1,kind=8)
           call functions(x3,az,bz,fz,fz1,fz2,ifz)
           do i2=1,n02
              x2 = hy*real(i2-n02/2-1,kind=8)
              call functions(x2,ay,by,fy,fy1,fy2,ify)
              do i1=1,n01
                 x1 = hx*real(i1-n01/2-1,kind=8)
                 call functions(x1,ax,bx,fx,fx1,fx2,ifx)
                 potential(i1,i2,i3) = -fx*fy*fz
                 density(i1,i2,i3) = (fx2*fy*fz+fx*fy2*fz+fx*fy*fz2+mu0**2*potential(i1,i2,i3))/(16.d0*datan(1.d0))
                 denval=max(denval,-density(i1,i2,i3))
              end do
           end do
        end do
     end select


     ! !acerioni: r = sqrt(x**2+y**2); V(x,y,z) = ArcTan(a*r)*f(z)/(a*r)
     ! do i3=1,n03
     !    x3 = hz*real(i3-n03/2-1,kind=8)
     !    call functions(x3,az,bz,fz,fz,1fz2,ifz)
     !    do i2=1,n02
     !       x2 = hy*real(i2-n02/2-1,kind=8)
     !       !call functions(x2,ay,by,fy,fy1,fy2,ify)
     !       do i1=1,n01
     !          x1 = hx*real(i1-n01/2-1,kind=8)
     !          r2 = x1*x1+x2*x2
     !          r = sqrt(r2)
     !          fxy = datan(factor*r)/(factor*r)
     !          !call functions(x1,ax,bx,fx,fx1,fx2,ifx)
     !          if (r == 0.d0) then
     !             potential(i1,i2,i3) = potential(i1,i2,i3) + 1.d0*fz
     !             density(i1,i2,i3) = density(i1,i2,i3) - fz*4.d0/3.d0*factor**2
     !             density(i1,i2,i3) = density(i1,i2,i3) + 1.d0*fz2
     !          else
     !             density(i1,i2,i3) = density(i1,i2,i3) + &
     !                  fz*(-3.d0*factor**2/(1+factor**2*r2)**2 - 1.d0/r2/(1+factor**2*r2)**2 + fxy/r2)
     !             density(i1,i2,i3) = density(i1,i2,i3) + fxy*fz2
     !             !denval=max(denval,-density(i1,i2,i3))
     !             potential(i1,i2,i3) = potential(i1,i2,i3) + fxy*fz
     !          end if
     !          density(i1,i2,i3) = density(i1,i2,i3) / (-16.d0*datan(1.d0))
     !          !density(i1,i2,i3) = -density(i1,i2,i3)
     !       end do
     !    end do
     ! end do
     ! !acerioni

     ! !acerioni: density = delta(x,y)*"constant = 1 along z"
     ! do i3=1,n03
     !    x3 = hz*real(i3-n03/2-1,kind=8)
     !    call functions(x3,az,bz,fz,fz1,fz2,ifz)
     !    do i2=1,n02
     !       x2 = hy*real(i2-n02/2-1,kind=8)
     !       !call functions(x2,ay,by,fy,fy1,fy2,ify)
     !       do i1=1,n01
     !          x1 = hx*real(i1-n01/2-1,kind=8)
     !          r2 = x1*x1+x2*x2
     !          r = sqrt(r2)
     !          fxy = datan(factor*r)/(factor*r)
     !          !call functions(x1,ax,bx,fx,fx1,fx2,ifx)
     !          if (r == 0.d0) then
     !             potential(i1,i2,i3) = 0.d0*fz
     !             density(i1,i2,i3) = 1.d0/(hx*hy)
     !          else
     !             potential(i1,i2,i3) = - 2.d0*log(r)
     !          end if
     !          !density(i1,i2,i3) = density(i1,i2,i3) / (-16.d0*datan(1.d0))
     !       end do
     !    end do
     ! end do
     ! !acerioni



     if (wrtfiles) then
     i2=n02/2
     do i3=1,n03
        do i1=1,n01
           !j1=n01/2+1-abs(n01/2+1-i1)
           !j2=n02/2+1-abs(n02/2+1-i2)
           !j3=n03/2+1-abs(n03/2+1-i3)
           write(unit,*) i1,i3,density(i1,i2,i3),potential(i1,i2,i3)
        end do
     end do
     end if


     if (ixc==0) denval=0.d0

  else

     !print *,'geometry code not admitted',geocode
     !stop
     call f_err_throw('geometry code not admitted "'//geocode//'"')

  end if



  !!! evaluation of the monopolar contribution !!!
  monopole = 0.d0
  dipole=0.0d0
  do i3 = 1, n03
     do i2 = 1, n02
        do i1 = 1, n01
           monopole = monopole + density(i1,i2,i3)
           dipole(1) = dipole(1) + density(i1,i2,i3)*real(i1-n01/2*hx,kind=8)
           dipole(2) = dipole(2) + density(i1,i2,i3)*real(i2-n02/2*hy,kind=8)
           dipole(3) = dipole(3) + density(i1,i2,i3)*real(i3-n03/2*hz,kind=8)
        end do
     end do
  end do
  !write(*,*) 'monopole = ', monopole
  !write(*,*) 'dipole = ', dipole
  if (iproc==0) then
     call yaml_map('monopole',monopole)
     call yaml_map('dipole',dipole)
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call f_close(unit)
  ! For ixc/=0 the XC potential is added to the solution, and an analytic comparison is no more
  ! possible. In that case the only possible comparison is between the serial and the parallel case
  ! To ease the comparison between the serial and the parallel case we add a random pot_ion
  ! to the potential.


  if (ixc==0) then
     potion_fac=0.d0
  else
     potion_fac=1.d0
  end if

  rhopot(:,:,:) = density(:,:,:) + denval
     do i3=1,n03
        do i2=1,n02
           do i1=1,n01
              call random_number(tt)
              !tt=0.d0!1.d0
              pot_ion(i1,i2,i3)=tt
              potential(i1,i2,i3)=potential(i1,i2,i3)+potion_fac*tt
!!              !for the ixc/=0 case
!!              call random_number(tt)
!!              rhopot(i1,i2,i3)=abs(tt)
           end do
        end do
     end do
     if (denval /= 0.d0) density=rhopot

end subroutine test_functions


!> Purpose: Compute exponential integral E1(x)
subroutine e1xb(x,e1)
  implicit none
  !Arguments
  real(kind=8), intent(in) :: x   !< x  Argument of E1(x)
  real(kind=8), intent(out) :: e1 !< E1 --- E1(x)  ( x > 0 )
  !Local variables
  real(kind=8), parameter :: ga=0.5772156649015328d0 !< EulerGamma
  real(kind=8) :: r,t0,t
  integer :: k,m

  if (x.eq.0.0) then
     e1=1.0d+300
  else if (x.le.1.0) then
     e1=1.0d0
     r=1.0d0
     do k=1,25
        r=-r*k*x/(k+1.0d0)**2
        e1=e1+r
        if (abs(r) <= abs(e1)*1.0d-15) then
           exit
        end if
     end do
      e1=-ga-dlog(x)+x*e1
   else
        m=20+int(80.0/x)
        t0=0.0d0
        do k=m,1,-1
           t0=k/(1.0d0+k/(x+t0))
        end do
           t=1.0d0/(x+t0)
           e1=dexp(-x)*t
        endif

end subroutine e1xb

end program PSolver_Program


!> Define the test functions
subroutine functions(x,a,b,f,f1,f2,whichone)
  use futile, dp => f_double
  use numerics
  implicit none
  integer, intent(in) :: whichone   !< Choose the function
  real(kind=8), intent(in) :: x     !< Argument of the function
  real(kind=8), intent(in) :: a,b   !< Parameter of the functions
  real(kind=8), intent(out) :: f    !< The value of the function
  real(kind=8), intent(out) :: f1   !< The value of the first derivative
  real(kind=8), intent(out) :: f2   !< The value of the second derivative
  !local variables
  !Type of function
  integer, parameter :: FUNC_CONSTANT = 1
  integer, parameter :: FUNC_GAUSSIAN = 2
  integer, parameter :: FUNC_GAUSSIAN_SHRINKED = 3
  integer, parameter :: FUNC_COSINE = 4
  integer, parameter :: FUNC_EXP_COSINE = 5
  integer, parameter :: FUNC_SHRINK_GAUSSIAN = 6
  integer, parameter :: FUNC_SINE = 7
  integer, parameter :: FUNC_ATAN = 8
  integer, parameter :: FUNC_ERF = 9

  real(kind=8) :: r,r2,y,yp,ys,factor,g,h,g1,g2,h1,h2
  real(kind=8) :: length,frequency,nu,sigma,agauss,derf

  select case(whichone)
  case(FUNC_CONSTANT)
     !constant
     f=1.d0
     f1=0.d0
     f2=0.d0
  case(FUNC_GAUSSIAN)
     !gaussian of sigma s.t. a=1/(2*sigma^2)
     r2=a*x**2
     f=dexp(-r2) !<checed
     f1=-2.d0*a*x*f !<checked
     f2=(-2.d0*a+4.d0*a*r2)*f !<checked
  case(FUNC_GAUSSIAN_SHRINKED)
     !gaussian "shrinked" with a=length of the system
     length=a
     r=pi*x/length
     y=tan(r)
!!$     yp=pi/length*1.d0/(dcos(r))**2
!!$     ys=2.d0*pi/length*y*yp
!!$     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
!!$     !!!!here we still need the first derivative
!!$     f2=factor*dexp(-y**2)
     f=dexp(-y**2) !<checked
     f1=-2.d0*pi*f*y/(length*cos(r)**2) !<checked
     f2=2.d0*pi**2*(2.d0*y**6 + y**4 - 2.d0*y**2 - 1.d0)/length**2*f !<checked
  case(FUNC_COSINE)
     !cosine with a=length, b=frequency
     length=a
     frequency=b
     r=frequency*pi*x/length
     f=dcos(r) !<checked
     f1=-dsin(r)*frequency*pi/length !<checked
     f2=-(frequency*pi/length)**2*dcos(r) !<checked
  case(FUNC_EXP_COSINE)
     !exp of a cosine, a=length
     nu=2.d0
     r=pi*nu/a*x
     y=cos(r)
     yp=-sin(r)
     f=exp(y) !<checked /dexp(1.0_dp) !<to be checked
     factor=(pi*nu/a)**2*(-y+yp**2)
     f1 = f*pi*nu/a*yp !<checked
     f2 = factor*f !<checked
  case(FUNC_SHRINK_GAUSSIAN)
     !gaussian times "shrinked" gaussian, sigma=length/10
     length=1.d0*a
     r=pi*x/length
     y=dtan(r)
     yp=pi/length*1.d0/(dcos(r))**2
     ys=2.d0*pi/length*y*yp
     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
     g=dexp(-y**2) !<checked
     g1=-2.d0*y*yp*g !<checked
     !g2=factor*dexp(-y**2)
     g2=2.d0*pi**2*(2.d0*y**6 + y**4 - 2.d0*y**2 - 1.d0)/length**2*g !<check
     sigma=length/10.0d0
     agauss=0.5d0/sigma**2
     r2=agauss*x**2
     h=dexp(-r2) !<checked
     h1=-2.d0*agauss*x*h !<checked
     h2=(-2.d0*agauss+4.d0*agauss*r2)*h !<checked
     f=g*h !<checked
     f1=g1*h+g*h1 !<checked
     f2=g2*h+g*h2+2.d0*g1*h1 !<checked
  case(FUNC_SINE)
     !sine with a=length, b=frequency
     length=a
     frequency=b
     r=frequency*pi*x/length
     f=dsin(r) !<checked
     f1=frequency*pi*cos(r)/length !<checked
     f2=-(frequency*pi/length)**2*sin(r) !<checked
  case(FUNC_ATAN)
     r=a*x
     factor=r**2+1.d0
     f=atan(r) !<checked
     f1=a/factor !<checked
     f2=-2.d0*r*a**2/factor**2 !<checked
!!$     !atan with a=length, b=frequency
!!$     length=a
!!$     nu = length
!!$     f=(datan(nu*x/length))**2 !<checked
!!$     !!here first derivative is lacking
!!$     f2=2.0d0*nu**2*length*(length-2.0d0*nu*x*f)/(length**2+nu**2*x**2)**2
  case(FUNC_ERF)
     !error function with a=sigma
     factor=sqrt(2.d0/pi)/a
     r=x
     y=x/(sqrt(2.d0)*a)
     if (abs(x)<=1.d-15) then
        f=factor
        f1=0.d0 !<checked
        f2=-sqrt(2.d0/pi)/(3.d0*a**3) !<checked
     else
        f=derf(y)/r
        y=x*x
        y=y/(2.d0*a**2)
        g=dexp(-y)
        h=1.d0/a**2+2.d0/x**2
        f1=-f/x+factor*g/x !<checked
        f2=-factor*g*h+2.d0*f/x**2  !<checked
     end if
  case default
     !print *,"Unknow function:",whichone
     !stop
     call f_err_throw('Unknown function '//trim(yaml_toa(whichone)))
  end select

end subroutine functions

subroutine compare(n01,n02,n03,potential,density,i1_max,i2_max,i3_max,max_diff)
  implicit none
  integer, intent(in) :: n01,n02,n03
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,density
  integer, intent(out) :: i1_max,i2_max,i3_max
  real(kind=8), intent(out) :: max_diff

  !local variables
  integer :: i1,i2,i3
  real(kind=8) :: factor
  max_diff = 0.d0
  i1_max = 1
  i2_max = 1
  i3_max = 1
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           factor=abs(potential(i1,i2,i3)-density(i1,i2,i3))
           if (max_diff < factor) then
              max_diff = factor
              i1_max = i1
              i2_max = i2
              i3_max = i3
           end if
        end do
     end do
  end do
end subroutine compare
