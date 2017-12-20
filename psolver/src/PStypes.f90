!> @file
!!    Modulefile for the definition of the basic structures
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)<br/>
!!    Copyright (C) 2002-2017 BigDFT group<br/>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Define the fortran types and the related routines to handle it.
module PStypes
  use f_enums
  use wrapper_MPI
  use PSbase
  use psolver_environment, only: cavity_data,cavity_default
  use dynamic_memory
  use f_input_file, only: ATTRS
  use box
  use dictionaries

  implicit none

  private

  character(len=*), parameter :: KERNEL_VARIABLES        = 'kernel'
  character(len=*), parameter :: SCREENING               = 'screening' 
  character(len=*), parameter :: ISF_ORDER               = 'isf_order' 
  character(len=*), parameter :: STRESS_TENSOR           = 'stress_tensor' 
  character(len=*), parameter :: ENVIRONMENT_VARIABLES   = 'environment'
  character(len=*), parameter :: CAVITY_KEY              = 'cavity' 
  character(len=*), parameter :: EPSILON_KEY             = 'epsilon' 
  character(len=*), parameter :: EDENSMAXMIN             = 'edensmaxmin' 
  character(len=*), parameter :: DELTA_KEY               = 'delta' 
  character(len=*), parameter :: FACT_RIGID              = 'fact_rigid' 
  character(len=*), parameter :: CAVITATION              = 'cavitation' 
  character(len=*), parameter :: GAMMAS_KEY              = 'gammaS' 
  character(len=*), parameter :: ALPHAS_KEY              = 'alphaS' 
  character(len=*), parameter :: BETAV_KEY               = 'betaV' 
  character(len=*), parameter :: GPS_ALGORITHM           = 'gps_algorithm'
  character(len=*), parameter :: RADII_SET               = 'radii_set'
  character(len=*), parameter :: ATOMIC_RADII            = 'atomic_radii'
  character(len=*), parameter :: PI_ETA                  = 'pi_eta' 
  character(len=*), parameter :: INPUT_GUESS             = 'input_guess' 
  character(len=*), parameter :: FD_ORDER                = 'fd_order' 
  character(len=*), parameter :: ITERMAX                 = 'itermax' 
  character(len=*), parameter :: MINRES                  = 'minres' 
  character(len=*), parameter :: PB_METHOD               = 'pb_method' 
  character(len=*), parameter :: PB_MINRES               = 'pb_minres' 
  character(len=*), parameter :: PB_ITERMAX              = 'pb_itermax' 
  character(len=*), parameter :: PB_INPUT_GUESS          = 'pb_input_guess' 
  character(len=*), parameter :: PB_ETA                  = 'pb_eta'
  character(len=*), parameter, public :: SETUP_VARIABLES = 'setup'
  character(len=*), parameter :: ACCEL                   = 'accel' 
  character(len=*), parameter :: KEEP_GPU_MEMORY         = 'keep_gpu_memory' 
  character(len=*), parameter :: USE_GPU_DIRECT          = 'use_gpu_direct' 
  character(len=*), parameter :: TASKGROUP_SIZE_KEY      = 'taskgroup_size' 
  character(len=*), parameter :: GLOBAL_DATA             = 'global_data' 
  character(len=*), parameter, public :: VERBOSITY               = 'verbose' 
  character(len=*), parameter :: OUTPUT                  = 'output' 
  character(len=*), parameter :: DICT_COMPLETED          = '__dict_has_been_checked__'//ATTRS

  integer, parameter :: RADII_PAULING_ID = 1
  integer, parameter :: RADII_BONDI_ID = 2
  integer, parameter :: RADII_UFF_ID = 3
  

  !> Defines the internal information for application of the FFT between the kernel and the density
  type, public :: FFT_metadata
     integer :: m1,m2,m3 !<original real dimension, with m2 in z direction and and m3 in y direction
     integer :: n1,n2,n3 !<dimension of the FFT operation, taking into account the zero padding if needed
     integer :: md1,md2,md3 !< Dimension of the real unpadded space, 
     !!md2 is further enlarged to be a multiple of number of processes
     integer :: nd1,nd2,nd3 !<fourier dimensions for which the kernel is injective,
     !!                formally 1/8 of the fourier grid. Here the dimension nd3 is
     !!                enlarged to be a multiple of nproc
     integer :: istart,iend,n3p !<start, endpoints and number of planes of the given processor
     real(dp) :: scal !<factor to rescale the solver such that the FFT is unitary, divided by 4pi
  end type FFT_metadata


  !> Define the work arrays needed for the treatment of the generalized (nonvacuum) or
  !! standard (vacuum) Poisson Equation. Not all of them are allocated, the actual memory usage 
  !! depends on the treatment.
  type, public :: PS_workarrays
     integer :: nat !< dimensions of the atomic based cavity. Zero if unused
     !>positions of the atoms of the cavity in the simulation box
     real(dp), dimension(:,:), pointer :: rxyz !< atomic coordinates
     real(dp), dimension(:), pointer :: radii !<radii of the cavity per atom
     !> dielectric function epsilon, continuous and differentiable in the whole
     !domain, to be used in the case of Preconditioned Conjugate Gradient (PCG) method
     real(dp), dimension(:,:), pointer :: eps
     !> logaritmic derivative of the dielectric function epsilon,
     !! to be used in the case of Polarization Iteration (SC) method
     real(dp), dimension(:,:,:,:), pointer :: dlogeps
     !> inverse of the dielectric function
     !! in the case of Polarization Iteration (SC) method
     !! inverse of the square root of epsilon
     !! in the case of the Preconditioned Conjugate Gradient (PCG)
     real(dp), dimension(:,:), pointer :: oneoeps
     !> correction term, given in terms of the multiplicative factor of nabla*eps*nabla
     !! to be used for Preconditioned Conjugate Gradient 
     real(dp), dimension(:,:), pointer :: corr
     !> ionic density, in the case of a Poisson-Boltzmann (PB) approach
     real(dp), dimension(:,:), pointer :: rho_ions
     !> inner rigid cavity to be integrated in the sccs method to avoit inner
     !! cavity discontinuity due to near-zero edens near atoms
     real(dp), dimension(:,:), pointer :: epsinnersccs
     !>work array needed to store the zero-padded part of the 
     !!density and the potential. Cannot be used as-is, it must be
     !!copied back into a distributed array with the 
     !!finalize_hartree_results routine.
     real(dp), dimension(:,:,:), pointer :: zf
     !> input guess vectors to be preserved for future use
     !!or work arrays, might be of variable dimension
     !!(either full of distributed)
     real(dp), dimension(:,:), pointer :: pot
     real(dp), dimension(:,:), pointer :: rho,rho_pb
     !> Polarization charge vector for print purpose only.
     real(dp), dimension(:,:), pointer :: rho_pol
     !> arrays for the execution of the PCG algorithm
     real(dp), dimension(:,:), pointer :: res,z,p,q

     integer(f_address) :: work1_GPU,work2_GPU,rho_GPU,pot_ion_GPU,k_GPU !<addresses for the GPU memory 
     integer(f_address) :: p_GPU,q_GPU,r_GPU,x_GPU,z_GPU,oneoeps_GPU,corr_GPU!<addresses for the GPU memory 
     !> GPU scalars. Event if they are scalars of course their address is needed
     integer(f_address) :: alpha_GPU, beta_GPU, kappa_GPU, beta0_GPU, eexctX_GPU, reduc_GPU, ehart_GPU
  end type PS_workarrays


  !> Datatype defining the mode for the running of the electrostatic solver
  type, public :: PSolver_options
     !> @copydoc poisson_solver::doc::datacode
     character(len=1) :: datacode
     !> integer variable setting the verbosity, from silent (0) to high
     integer :: verbosity_level
     !> if .true., and the cavity is set to 'sccs' attribute, then the epsilon is updated according to rhov
     logical :: update_cavity
     !> if .true. calculate the stress tensor components.
     logical :: calculate_strten
     !> Use the input guess procedure for the solution of the generalized Poisson equation (PCG or SC).
     !! This value is ignored in the case of a vacuum solver
     logical :: use_input_guess
     !> For an implicit solvation run, trigger the energy calculation to the only electrostatic contribution.
     !! If .true., the code only calculates the electrostatic contribution
     !! and the non electrostatic terms to the total solvation energy (cavitation, 
     !! repulsion and dispersion) are neglected.
     logical :: only_electrostatic
     !> extract the polarization charge and the dielectric function, to be used for plotting purposes
     logical :: cavity_info
     !> Use the input guess procedure for the solution in the case of the Poisson Boltzmann Equation.
     !! This option is ignored in the case of a neutral solution
     logical :: use_pb_input_guess
     !> Total integral on the supercell of the final potential on output
     !! clearly meaningful only for Fully periodic BC, ignored in the other cases.
     !> Prepare the information in the rigid cavity case which is needed
     !! to calculate the forces (nabla2pot times epsilon minus one)
     logical :: final_call
     real(gp) :: potential_integral
  end type PSolver_options


  !> Defines the fundamental structure for the kernel
  type, public :: coulomb_operator
     !variables with physical meaning
     integer :: itype_scf             !< Order of the ISF family to be used
     real(gp) :: mu                   !< Inverse screening length for the Helmholtz Eq. (Poisson Eq. -> mu=0)
!!$     !> geocode is used in all the code to specify the boundary conditions (BC) the problem:
!!$     !!          - 'F' free BC, isolated systems.
!!$     !!                The program calculates the solution as if the given density is
!!$     !!                "alone" in R^3 space.
!!$     !!          - 'S' surface BC, isolated in y direction, periodic in xz plane                
!!$     !!                The given density is supposed to be periodic in the xz plane,
!!$     !!                so the dimensions in these direction mus be compatible with the FFT
!!$     !!                Beware of the fact that the isolated direction is y!
!!$     !!          - 'P' periodic BC.
!!$     !!                The density is supposed to be periodic in all the three directions,
!!$     !!                then all the dimensions must be compatible with the FFT.
!!$     !!                No need for setting up the kernel (in principle for Plane Waves)
!!$     !!          - 'W' Wires BC.
!!$     !!                The density is supposed to be periodic in z direction, 
!!$     !!                which has to be compatible with the FFT.
!!$     !!          - 'H' Helmholtz Equation Solver
!!$!     character(len=1) :: geocode
!!$     !> method of embedding in the environment
!!$     !!          - 'VAC' Poisson Equation in vacuum. Default case.
!!$     !!          - 'PCG' Generalized Poisson Equation, Preconditioned Conjugate Gradient
!!$     !!          - 'PI'  Generalized Poisson Equation, Polarization Iteration method
!!$     !character(len=3) :: method 
!!$     !! this represents the information for the equation and the algorithm to be solved
!!$     !! this enumerator contains the algorithm and has the attribute associated to the 
!!$     !! type of cavity to be used
     type(f_enumerator) :: method
     type(cell) :: mesh !< structure which includes all cell informations 
!!$     integer, dimension(3) :: ndims   !< dimension of the box of the density
!!$     real(gp), dimension(3) :: hgrids !<grid spacings in each direction
!!$     real(gp), dimension(3) :: angrad !< angles in radiants between each of the axis
     type(cavity_data) :: cavity !< description of the cavity for the dielectric medium
     type(PSolver_options) :: opt !<Datatype controlling the operations of the solver
     real(dp), dimension(:), pointer :: kernel !< kernel of the Poisson Solver
     integer, dimension(5) :: plan
     integer, dimension(3) :: geo
     !>workarrays for the application of the Solver. Might have different
     !!memory footprints dependently of the treatment.
     type(PS_workarrays) :: w
     !variables with computational meaning
     type(mpi_environment) :: mpi_env !< complete environment for the Poisson Solver
     type(mpi_environment) :: inplane_mpi,part_mpi !<mpi_environment for internal ini-plane parallelization
     type(FFT_metadata) :: grid !<dimensions of the FFT grid associated to this kernel
     logical :: use_gpu_direct 
     integer :: igpu !< control the usage of the GPU
     integer :: gpuPCGRed !< control if GPU can be used for PCG reductions
     integer :: initCufftPlan
     integer :: keepGPUmemory
     integer :: stay_on_gpu
     integer :: keepzf
     !parameters for the iterative methods
     !> Order of accuracy for the finite difference nabla, used for nonvacuum calculations (only in SC).
     !! Usead also for the additional contribution to the KS potential for an implicit solvation run.
     !! Represents the total number of points at left and right of the x0 where we want to calculate the derivative.
     integer :: nord
     !> default set of radii for the rigid cavity
     integer :: radii_set
     !> dictionary of the atom species radii defined by the user
     type(dictionary), pointer :: radii_dict
     integer :: max_iter    !< Maximum number of convergence iterations for PCG or SC.
     real(dp) :: minres     !< Convergence criterion for the PCG or SC loop.
     real(dp) :: PI_eta     !< Mixing parameter for the update of SC iteration
     integer :: max_iter_PB !< Max conv iterations for PB treatment
     real(dp) :: minres_PB  !< Convergence criterion for PB residue
     real(dp) :: PB_eta     !< Mixing scheme for PB
     real(dp) :: IntVol     !< Volume integral needed for the non-electrostatic energy contributions
     real(dp) :: IntSur     !< Surface integral needed for the non-electrostatic energy contributions
     
     integer, dimension(:), pointer :: counts    !< Array needed to gather the information of the Poisson solver
     integer, dimension(:), pointer :: displs    !< Array needed to gather the information of the Poisson solver
     integer, dimension(:), pointer :: rhocounts !< Array needed to gather the information of the Poisson solver on multiple gpus
     integer, dimension(:), pointer :: rhodispls !< Array needed to gather the information of the Poisson solver on multiple gpus
  end type coulomb_operator


  !> Define the energy terms for the Poisson and Generalized poisson operator applications
  type, public :: PSolver_energies
     !> hartree energy, defined as the @f$\int \rho(\mathbf{r}) V(\mathbf{r}) \mathrm d r @f$, with @f$\rho@f$ being the 
     !! input density and @f$V@f$ the potential defined by this density
     !! in the case when rho_ion is passed, the electrostatic contribution is only filled
     real(gp) :: hartree 
     !> electrostatic energy, defined as the hartree energy but with @f$\rho@f$ and @f$V@f$ coming from @f$\rho + \rho_{ion}@f$
     !! the hartree energy can be obtained by subtraction with the potential energy terms
     real(gp) :: elec 
     !> Energetic term coming from the @f$\int \rho V_extra@f$, in the case of a @f$\rho@f$-dependent cavity.
     !! Clearly this term is calculated only if the potential is corrected. When the cavity is fixed, the eVextra is zero.
     real(gp) :: eVextra
     !> Nonelectrostatic contribution, function of the cavity surface and volume (only for implicit solvation).
     real(gp) :: cavitation
     !> stress tensor, to be calculated when calculate_strten is .true.
     real(gp), dimension(6) :: strten
  end type PSolver_energies


  public :: pkernel_null,PSolver_energies_null,pkernel_free,pkernel_allocate_cavity
  public :: pkernel_set_epsilon,PS_allocate_cavity_workarrays,build_cavity_from_rho
  public :: ps_allocate_lowlevel_workarrays,PSolver_options_null,PS_input_dict
  public :: release_PS_potential,PS_release_lowlevel_workarrays,PS_set_options,pkernel_init
  public :: ps_soft_PCM_forces,pkernel_get_radius,epsilon_inner_cavity

  !To specify properly to doxygen (test)
  private :: free_PS_workarrays, PS_fill_variables, PS_input_fill


contains


  pure function PSolver_energies_null() result(e)
    implicit none
    type(PSolver_energies) :: e
    e%hartree    =0.0_gp
    e%elec       =0.0_gp
    e%eVextra    =0.0_gp
    e%cavitation =0.0_gp
    e%strten     =0.0_gp
  end function PSolver_energies_null

  pure function PSolver_options_null() result(o)
    implicit none
    type(PSolver_options) :: o

    o%datacode           ='X'
    o%verbosity_level    =0
    o%update_cavity      =.false.
    o%calculate_strten   =.false.
    o%use_input_guess    =.false.
    o%use_pb_input_guess =.false.
    o%cavity_info        =.false.
    o%only_electrostatic =.true.
    o%final_call         =.false.
    o%potential_integral =0.0_gp
   end function PSolver_options_null

  pure function FFT_metadata_null() result(d)
    implicit none
    type(FFT_metadata) :: d
    d%m1=0
    d%m2=0
    d%m3=0
    d%n1=0
    d%n2=0
    d%n3=0
    d%md1=0
    d%md2=0
    d%md3=0
    d%nd1=0
    d%nd2=0
    d%nd3=0
    d%istart=0
    d%iend=0
    d%n3p=0
    d%scal=0.0_dp
  end function FFT_metadata_null

  pure subroutine nullify_work_arrays(w)
    use f_utils, only: f_zero
    implicit none
    type(PS_workarrays), intent(out) :: w
    call f_zero(w%nat)
    nullify(w%radii)
    nullify(w%rxyz)
    nullify(w%eps)
    nullify(w%dlogeps)
    nullify(w%oneoeps)
    nullify(w%corr)
    nullify(w%epsinnersccs)
    nullify(w%rho_pol)
    nullify(w%rho)
    nullify(w%rho_pb)
    nullify(w%pot)
    nullify(w%res)
    nullify(w%zf)
    nullify(w%z)
    nullify(w%p)
    nullify(w%q)
    nullify(w%eps)
    nullify(w%rho_ions)
    call f_zero(w%work1_GPU)
    call f_zero(w%work2_GPU)
    call f_zero(w%rho_GPU)
    call f_zero(w%pot_ion_GPU)
    call f_zero(w%k_GPU)
    call f_zero(w%p_GPU)
    call f_zero(w%q_GPU)
    call f_zero(w%r_GPU)
    call f_zero(w%x_GPU)
    call f_zero(w%z_GPU)
    call f_zero(w%oneoeps_GPU)
    call f_zero(w%corr_GPU)
    call f_zero(w%alpha_GPU)
    call f_zero(w%beta_GPU)
    call f_zero(w%kappa_GPU)
    call f_zero(w%beta0_GPU)
    call f_zero(w%eexctX_GPU)
    call f_zero(w%ehart_GPU)
    call f_zero(w%reduc_GPU)
  end subroutine nullify_work_arrays

  pure function pkernel_null() result(k)
    use psolver_environment, only : PS_VAC_ENUM
    implicit none
    type(coulomb_operator) :: k
    k%itype_scf=0
!    k%geocode='F'
    call nullify_f_enum(k%method)
    k%cavity=cavity_default()
    k%opt=PSolver_options_null()
    k%mu=0.0_gp
    k%mesh=cell_null()
!    k%ndims=(/0,0,0/)
!    k%hgrids=(/0.0_gp,0.0_gp,0.0_gp/)
    nullify(k%kernel)
    k%plan=(/0,0,0,0,0/)
    k%geo=(/0,0,0/)
    call nullify_work_arrays(k%w)
    call nullify_mpi_environment(k%mpi_env)
    call nullify_mpi_environment(k%inplane_mpi)
    call nullify_mpi_environment(k%part_mpi)
    k%grid=FFT_metadata_null()
    k%igpu=0
    k%initCufftPlan=0
    k%keepGPUmemory=1
    k%use_gpu_direct=.false.
    k%keepzf=1
    k%nord=0
    k%max_iter=0
    k%PI_eta=0.0_dp
    k%minres=0.0_dp
    k%minres_PB=0.0_dp
    k%radii_set=0
    k%IntVol=0.0_dp
    k%IntSur=0.0_dp
    nullify(k%radii_dict)
    nullify(k%counts)
    nullify(k%displs)
  end function pkernel_null

  subroutine free_PS_workarrays(iproc,igpu,keepzf,gpuPCGred,keepGPUmemory,w)
    use dictionaries, only: f_err_throw
    use f_utils, only: f_zero
    implicit none
    integer, intent(in) :: keepzf,gpuPCGred,keepGPUmemory,igpu,iproc
    integer :: i_stat
    type(PS_workarrays), intent(inout) :: w
    call f_zero(w%nat)
    call f_free_ptr(w%radii)
    call f_free_ptr(w%rxyz)
    call f_free_ptr(w%eps)
    call f_free_ptr(w%dlogeps)
    call f_free_ptr(w%oneoeps)
    call f_free_ptr(w%corr)
    call f_free_ptr(w%epsinnersccs)
    call f_free_ptr(w%rho_pol)
    call f_free_ptr(w%pot)
    call f_free_ptr(w%rho)
    call f_free_ptr(w%rho_ions)
    call f_free_ptr(w%res)
    call f_free_ptr(w%z)
    call f_free_ptr(w%p)
    call f_free_ptr(w%q)
    if(keepzf == 1) call f_free_ptr(w%zf)
    if (gpuPCGRed == 1) then
       if (keepGPUmemory == 1) then
          call cudafree(w%z_GPU)
          call cudafree(w%r_GPU)
          call cudafree(w%oneoeps_GPU)
          call cudafree(w%p_GPU)
          call cudafree(w%q_GPU)
          call cudafree(w%x_GPU)
          call cudafree(w%corr_GPU)
          call cudafree(w%alpha_GPU)
          call cudafree(w%beta_GPU)
          call cudafree(w%beta0_GPU)
          call cudafree(w%kappa_GPU)
          call cudafree(w%ehart_GPU)
          call cudafree(w%eexctX_GPU)
          call cudafree(w%reduc_GPU)
       end if
    end if
    if (igpu == 1) then
       if (iproc == 0) then
         call cudadestroystream(i_stat)
         if (i_stat /= 0) call f_err_throw('error freeing stream ')
         call cudadestroycublashandle()
          if (keepGPUmemory == 1) then
             call cudafree(w%work1_GPU)
             call cudafree(w%work2_GPU)
             call cudafree(w%rho_GPU)
             call cudafree(w%pot_ion_GPU)
          endif
          call cudafree(w%k_GPU)
       endif
    end if

  end subroutine free_PS_workarrays

  subroutine release_PS_potential(keepzf,w,use_input_guess)
    implicit none
    integer, intent(in) :: keepzf
    type(PS_workarrays), intent(inout) :: w
    logical, intent(in) :: use_input_guess
    if(keepzf /= 1) call f_free_ptr(w%zf)
!    call f_free_ptr(w%pot)
    if (.not. use_input_guess) call f_free_ptr(w%pot)
  end subroutine release_PS_potential

  !> Free memory used by the kernel operation and all workarrays.
  subroutine pkernel_free(kernel)
    use dictionaries, only: dict_free
    implicit none
    type(coulomb_operator), intent(inout) :: kernel
    call dict_free(kernel%radii_dict)
    call f_free_ptr(kernel%kernel)
    call f_free_ptr(kernel%counts)
    call f_free_ptr(kernel%displs)
    call free_PS_workarrays(kernel%mpi_env%iproc,kernel%igpu,&
         kernel%keepzf,kernel%gpuPCGred,kernel%keepGPUmemory,kernel%w)
    !free GPU data
    if (kernel%igpu == 1) then
       if (kernel%mpi_env%iproc == 0) then
          call f_free_ptr(kernel%rhocounts)
          call f_free_ptr(kernel%rhodispls)
          if (kernel%initCufftPlan == 1) then
             call cufftDestroy(kernel%plan(1))
             call cufftDestroy(kernel%plan(2))
             call cufftDestroy(kernel%plan(3))
             call cufftDestroy(kernel%plan(4))
          endif
       endif
    end if
    call release_mpi_environment(kernel%inplane_mpi)
    call release_mpi_environment(kernel%part_mpi)
    call release_mpi_environment(kernel%mpi_env)
  end subroutine pkernel_free

  !> Initialization of the Poisson kernel starting from the box and boundary
  !! conditions data.
  function pkernel_init(iproc,nproc,dict,geocode,ndims,hgrids,alpha_bc,beta_ac,gamma_ab,mpi_env) result(kernel)
    use yaml_output
    use dictionaries
    use numerics
    use wrapper_MPI
    use f_enums
    use psolver_environment
    use box, only: cell_new
    use f_input_file, only: input_file_dump
    implicit none
    integer, intent(in) :: iproc      !< Proc Id
    integer, intent(in) :: nproc      !< Number of processes
    type(dictionary), pointer :: dict !< dictionary of the input variables
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, dimension(3), intent(in) :: ndims
    real(gp), dimension(3), intent(in) :: hgrids
    real(gp), intent(in), optional :: alpha_bc,beta_ac,gamma_ab
    type(mpi_environment), intent(in), optional :: mpi_env
    type(coulomb_operator) :: kernel
    !local variables
    integer :: nthreads,group_size,taskgroup_size
    !$ integer :: omp_get_max_threads

    !nullification
    kernel=pkernel_null()

    !mesh initialization
    kernel%mesh=cell_new(geocode,ndims,hgrids,alpha_bc,beta_ac,gamma_ab)

    !new treatment for the kernel input variables
    kernel%method=PS_VAC_ENUM
    if (DICT_COMPLETED .notin. dict) call PS_input_dict(dict) !complete the dictionary

    call PS_fill_variables(kernel,kernel%opt,dict) !fill the structure with basic results

    kernel%keepzf=1

    !import the mpi_environment if present
    if (present(mpi_env)) then
       call copy_mpi_environment(src=mpi_env,dest=kernel%mpi_env)
    else

       !specialized treatment
       taskgroup_size=dict//SETUP_VARIABLES//TASKGROUP_SIZE_KEY

       group_size=nproc
       !if the taskgroup size is not a divisor of nproc do not create taskgroups
       if (nproc >1 .and. taskgroup_size > 0 .and. taskgroup_size < nproc .and.&
            mod(nproc,taskgroup_size)==0) then
          group_size=taskgroup_size
       end if
       call mpi_environment_set(kernel%mpi_env,iproc,nproc,mpiworld(),group_size)
    end if

    !gpu can be used only for one nproc
    if (nproc > 1) kernel%igpu=0

    !-------------------
    nthreads=0
    if (kernel%mpi_env%iproc == 0 .and. kernel%mpi_env%igroup == 0 .and. kernel%opt%verbosity_level==1) then
       if (kernel%mu==0.0_gp) then 
          call yaml_comment('Kernel Initialization',hfill='-')
          call yaml_mapping_open('Poisson Kernel Initialization')
       else
          call yaml_mapping_open('Helmholtz Kernel Initialization')
          call yaml_map('Screening Length (AU)',1.0_gp/kernel%mu,fmt='(g25.17)')
       end if
       !we might also perform an input_file dump if needed
       call input_file_dump(dict)
       !$ nthreads = omp_get_max_threads()
       call yaml_map('MPI tasks',kernel%mpi_env%nproc)
       if (nthreads /=0) call yaml_map('OpenMP threads per MPI task',nthreads)
       if (kernel%igpu==1) call yaml_map('Kernel copied on GPU',.true.)
       if (kernel%method /= 'VAC') call yaml_map('Iterative method for Generalised Equation',str(kernel%method))
       if (kernel%method .hasattr. PS_RIGID_ENUM) call yaml_map('Cavity determination','rigid')
       if (kernel%method .hasattr. PS_SCCS_ENUM) call yaml_map('Cavity determination','sccs')
       call yaml_mapping_close() !kernel
    end if

  end function pkernel_init

  !>modifies the options of the poisson solver to switch certain options
  subroutine PS_set_options(kernel,global_data,calculate_strten,verbose,&
       update_cavity,use_input_guess,cavity_info,cavitation_terms,&
       potential_integral,final_call)
    implicit none
    type(coulomb_operator), intent(inout) :: kernel
    logical, intent(in), optional :: global_data,calculate_strten,verbose
    logical, intent(in), optional :: update_cavity,use_input_guess
    logical, intent(in), optional :: cavity_info,cavitation_terms
    logical, intent(in), optional :: final_call
    real(gp), intent(in), optional :: potential_integral

    if (present(global_data)) then
       if (global_data) then
          kernel%opt%datacode='G'
       else
          kernel%opt%datacode='D'
       end if
    end if
    if (present(calculate_strten)) kernel%opt%calculate_strten=calculate_strten
    if (present(verbose)) then
       if (verbose) then
          kernel%opt%verbosity_level=1
       else
          kernel%opt%verbosity_level=0
       end if
    end if
    if (present(update_cavity   )) kernel%opt%update_cavity=update_cavity
    if (present(use_input_guess )) kernel%opt%use_input_guess=use_input_guess
    if (present(cavity_info     )) kernel%opt%cavity_info=cavity_info
    if (present(cavitation_terms)) kernel%opt%only_electrostatic=.not. cavitation_terms
    if (present(potential_integral)) kernel%opt%potential_integral =potential_integral
    if (present(final_call)) kernel%opt%final_call=final_call

  end subroutine PS_set_options


  !>routine to fill the input variables of the kernel
  subroutine PS_input_dict(dict,dict_minimal)
    use dictionaries
    use f_input_file
    use yaml_parse
    implicit none
    !>input dictionary, a copy of the user input, to be filled
    !!with all the variables on exit
    type(dictionary), pointer :: dict
    type(dictionary), pointer, optional :: dict_minimal
    !local variables
    !integer(f_integer) :: params_size
    !integer(kind = 8) :: cbuf_add !< address of c buffer
    !character, dimension(:), allocatable :: params
    type(dictionary), pointer :: parameters
    type(dictionary), pointer :: parsed_parameters
    type(dictionary), pointer :: profiles
    type(dictionary), pointer :: nested,asis
    external :: get_ps_inputvars

    call f_routine(id='PS_input_dict')

    nullify(parameters,parsed_parameters,profiles)

!!$    !alternative filling of parameters from hard-coded source file
!!$    call getpsinputdefsize(params_size)
!!$    !allocate array
!!$    params=f_malloc_str(1,params_size,id='params')
!!$    !fill it and parse dictionary
!!$    call getpsinputdef(params)
!!$    call yaml_parse_from_char_array(parsed_parameters,params)
!!$    call f_free_str(1,params)

    !new filing method, uses database parsing
    call yaml_parse_database(parsed_parameters,get_ps_inputvars)
    !for each of the documents in the input variables specifications
    parameters=>parsed_parameters//0
    profiles => parsed_parameters//1


    call input_file_complete(parameters,dict,imports=profiles)

    if (present(dict_minimal)) then
       nullify(nested,asis)
       call input_file_minimal(parameters,dict,dict_minimal,nested,asis)
    end if

    if (associated(parsed_parameters)) then
       call dict_free(parsed_parameters)
       nullify(parameters)
       nullify(profiles)
    else
       call dict_free(parameters)
    end if

    !write in the dictionary that it has been completed
    call set(dict//DICT_COMPLETED,.true.)

    call f_release_routine()

  end subroutine PS_input_dict

  subroutine PS_fill_variables(k,opt,dict)
    use dictionaries
    implicit none
    type(coulomb_operator), intent(inout) :: k
    type(PSolver_options), intent(inout) :: opt
    type(dictionary), pointer :: dict
    !local variables
    type(dictionary), pointer :: lvl,var

    ! Transfer dict values into input_variables structure.
    lvl => dict_iter(dict)
    do while(associated(lvl))
       var => dict_iter(lvl)
       do while(associated(var))
          call PS_input_fill(k,opt,dict_key(lvl),var)
          var => dict_next(var)
       end do
       lvl => dict_next(lvl)
    end do

  end subroutine PS_fill_variables

  !> Set the dictionary from the input variables
  subroutine PS_input_fill(k,opt, level, val)
    use PSbase
    use psolver_environment
    use yaml_output, only: yaml_warning
    use dictionaries
    use numerics
    implicit none
    type(coulomb_operator), intent(inout) :: k
    type(PSolver_options), intent(inout) :: opt
    type(dictionary), pointer :: val
    character(len = *), intent(in) :: level
    !local variables
    logical :: dummy_l
    real(gp) :: dummy_d
    !integer, dimension(2) :: dummy_int !<to use as filling for input variables
    real(gp), dimension(2) :: dummy_gp !< to fill the input variables
    !logical, dimension(2) :: dummy_log !< to fill the input variables
    !character(len=256) :: dummy_char
    character(len = max_field_length) :: strn

    if (index(dict_key(val), "_attributes") > 0) return

    select case(trim(level))
    case(KERNEL_VARIABLES)
       select case (trim(dict_key(val)))
       case(SCREENING)
          k%mu=val
       case(ISF_ORDER)
          k%itype_scf=val
       case(STRESS_TENSOR)
          opt%calculate_strten=val
       case DEFAULT
          if (k%mpi_env%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (ENVIRONMENT_VARIABLES)
       select case (trim(dict_key(val)))
       case (CAVITY_KEY)
          strn=val

          select case(trim(strn))
          case('vacuum')
             call f_enum_attr(k%method,PS_NONE_ENUM)
          case('soft-sphere')
             call f_enum_attr(k%method,PS_RIGID_ENUM)
          case('sccs')   
             call f_enum_attr(k%method,PS_SCCS_ENUM)
          end select
       case (EPSILON_KEY)
          k%cavity%epsilon0=val
       case (EDENSMAXMIN)
          dummy_gp=val
          k%cavity%edensmin=dummy_gp(1)
          k%cavity%edensmax=dummy_gp(2)
       case (DELTA_KEY)
          dummy_d=val
          ! Divided by 4 because both rigid cavities are 4*delta spread 
          k%cavity%delta=dummy_d  !0.25_gp*dummy_d
       case (FACT_RIGID)
          k%cavity%fact_rigid=val
       case (CAVITATION)
          dummy_l=val
          opt%only_electrostatic=.not. dummy_l
       case (GAMMAS_KEY)
          dummy_d=val
          k%cavity%gammaS=dummy_d*SurfAU
       case (ALPHAS_KEY)
          dummy_d=val
          k%cavity%alphaS=dummy_d*SurfAU
       case (BETAV_KEY)
          dummy_d=val
          k%cavity%betaV=dummy_d/AU_GPa
       case (GPS_ALGORITHM)
          strn=val
          select case(trim(strn))
          case('SC')
             call f_enum_update(dest=k%method,src=PS_PI_ENUM)
          case('PCG')
             call f_enum_update(dest=k%method,src=PS_PCG_ENUM)
          end select
       case(RADII_SET)
          strn=val
          select case(trim(strn))
          case('UFF')
             k%radii_set=RADII_UFF_ID
          case('Pauling')
             k%radii_set=RADII_PAULING_ID
          case('Bondi')
             k%radii_set=RADII_BONDI_ID
          end select
       case(ATOMIC_RADII)
          if (dict_size(val) > 0) call dict_copy(k%radii_dict,val)
       case (PI_ETA)
          k%PI_eta=val
       case (INPUT_GUESS)
          opt%use_input_guess=val
       case (FD_ORDER)
          k%nord=val
       case (ITERMAX)
          k%max_iter=val
       case (MINRES)
          k%minres=val
       case (PB_METHOD)
          strn=val
          select case(trim(strn))
          case('none')
             call f_enum_attr(k%method,PS_PB_NONE_ENUM)
          case('linear')
             call f_enum_attr(k%method,PS_PB_LINEAR_ENUM)
          case('standard')
             call f_enum_attr(k%method,PS_PB_STANDARD_ENUM)
          case('modified')
             call f_enum_attr(k%method,PS_PB_MODIFIED_ENUM)
          end select
       case (PB_MINRES)
          k%minres_PB=val
       case (PB_ITERMAX)
          k%max_iter_PB=val
       case (PB_INPUT_GUESS)
          opt%use_pb_input_guess=val
       case (PB_ETA)
          k%PB_eta=val
       case DEFAULT
          if (k%mpi_env%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (SETUP_VARIABLES)
       select case (trim(dict_key(val)))       
       case (ACCEL)
          strn=val
          select case(trim(strn))
          case('CUDA')
             k%igpu=1
          case('none')
             k%igpu=0
          end select
       case (KEEP_GPU_MEMORY)
          dummy_l=val
          if (dummy_l) then
             k%keepGPUmemory=1
          else
             k%keepGPUmemory=0
          end if
       case (USE_GPU_DIRECT)
          k%use_gpu_direct=val
       case (TASKGROUP_SIZE_KEY)

       case (GLOBAL_DATA)
          dummy_l=val
          if (dummy_l) then
             opt%datacode='G'
          else
             opt%datacode='D'
          end if
       case (VERBOSITY)
          dummy_l=val
          if (dummy_l) then
             opt%verbosity_level=1
          else
             opt%verbosity_level=0
          end if
       case (OUTPUT)
          !for the moment no treatment, to be added
       case DEFAULT
          if (k%mpi_env%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case DEFAULT
    end select
  END SUBROUTINE PS_input_fill


  !> allocate the workarrays needed to perform a nonvacuum GPe operation
  subroutine PS_allocate_lowlevel_workarrays(poisson_boltzmann,cudasolver,rho,kernel)
    use f_utils, only: f_zero
    use wrapper_linalg, only: axpy
    implicit none
    logical, intent(in) :: cudasolver,poisson_boltzmann
    type(coulomb_operator), intent(inout) :: kernel
    real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(in) :: rho !< initial rho, needed for PCG
    !local variables
    integer :: n1,n23

    !we need to reallocate the zf array with the right size when called with stress_tensor and gpu
    if(kernel%keepzf == 1) then
       if(kernel%igpu==1 .and. .not. cudasolver) then !LG: what means that?
          call f_free_ptr(kernel%w%zf)
          kernel%w%zf = f_malloc_ptr([kernel%grid%md1, kernel%grid%md3, &
               2*kernel%grid%md2/kernel%mpi_env%nproc],id='zf')
       end if
    else
       kernel%w%zf = f_malloc_ptr([kernel%grid%md1, kernel%grid%md3, &
            2*kernel%grid%md2/kernel%mpi_env%nproc],id='zf')
    end if
    n23=kernel%grid%m3*kernel%grid%n3p
    n1=kernel%grid%m1

    select case(trim(str(kernel%method)))
    case('PCG')
!!$       if (use_input_guess .and. &
!!$            associated(kernel%w%pot)) then
!!$       else
!!$          !allocate if it is the first time
!!$          if (associated(kernel%w%pot)) then
!!$             !call f_zero(kernel%w%pot)
!!$          else
!!$             kernel%w%pot=f_malloc0_ptr([n1,n23],id='pot')
!!$          end if
!!$       end if
       if (.not. associated(kernel%w%pot)) kernel%w%pot=f_malloc0_ptr([n1,n23],id='pot')
       kernel%w%res=f_malloc_ptr([n1,n23],id='res')
       call f_memcpy(src=rho,dest=kernel%w%res)

       kernel%w%q=f_malloc0_ptr([n1,n23],id='q')
       kernel%w%p=f_malloc0_ptr([n1,n23],id='p')
       kernel%w%z=f_malloc_ptr([n1,n23],id='z')
       kernel%w%rho_pol=f_malloc_ptr([n1,n23],id='rho_pol')
    case('PI')
!!$       if (use_input_guess .and. &
!!$            associated(kernel%w%pot)) then
!!$       else
!!$          !allocate if it is the first time
!!$          if (associated(kernel%w%pot)) then
!!$             !call f_zero(kernel%w%pot)
!!$          else
!!$             kernel%w%pot=f_malloc_ptr([kernel%ndims(1),kernel%ndims(2)*kernel%ndims(3)],id='pot')
!!$          end if
!!$       end if

       if (.not. associated(kernel%w%pot))&
            kernel%w%pot=f_malloc0_ptr([kernel%mesh%ndims(1),kernel%mesh%ndims(2)*kernel%mesh%ndims(3)],id='pot')

       !kernel%w%pot=f_malloc_ptr([kernel%mesh%ndims(1),kernel%mesh%ndims(2)*kernel%mesh%ndims(3)],id='pot')
       kernel%w%rho=f_malloc0_ptr([kernel%mesh%ndims(1),kernel%mesh%ndims(2)*kernel%mesh%ndims(3)],id='rho')
       kernel%w%rho_pol=f_malloc_ptr([n1,n23],id='rho_pol')
    end select

    if (poisson_boltzmann) then
       kernel%w%rho_pb=f_malloc_ptr([n1,n23],id='rho_pb')
       call f_memcpy(src=rho,dest=kernel%w%rho_pb)
    end if

  end subroutine PS_allocate_lowlevel_workarrays

  !> this is useful to deallocate useless space and to 
  !! also perform extra treatment for the inputguess
  subroutine PS_release_lowlevel_workarrays(kernel,keep_rhopol)
    use wrapper_linalg, only: axpy
    implicit none
    type(coulomb_operator), intent(inout) :: kernel
    logical, intent(in) :: keep_rhopol

    select case(trim(str(kernel%method)))
    case('PCG')
       call f_free_ptr(kernel%w%res)
       call f_free_ptr(kernel%w%q)
       call f_free_ptr(kernel%w%p)
       call f_free_ptr(kernel%w%z)
       if (.not. keep_rhopol) call f_free_ptr(kernel%w%rho_pol) 
    case('PI')
       call f_free_ptr(kernel%w%rho)
       if (.not. keep_rhopol) call f_free_ptr(kernel%w%rho_pol)
    end select

    call f_free_ptr(kernel%w%rho_pb)

  end subroutine PS_release_lowlevel_workarrays

  !> allocate the workarrays for the first initialization
  !! their allocation depends on the treatment which we are going to
  !! apply
  subroutine PS_allocate_cavity_workarrays(n1,n23,ndims,method,w)
    use psolver_environment, only: PS_PB_NONE_ENUM
    use dynamic_memory
    implicit none
    integer, intent(in) :: n1,n23
    integer, dimension(3), intent(in) :: ndims
    type(f_enumerator), intent(in) :: method
    type(PS_workarrays), intent(inout) :: w

    select case(trim(str(method)))
    case('PCG')
       !w%rho_pol=f_malloc_ptr([n1,n23],id='rho_pol') !>>>>>>>>>>here the switch
       w%eps=f_malloc_ptr([n1,n23],id='eps')
       w%corr=f_malloc_ptr([n1,n23],id='corr')
       w%oneoeps=f_malloc_ptr([n1,n23],id='oneosqrteps')
       w%epsinnersccs=f_malloc_ptr([n1,n23],id='epsinnersccs')
    case('PI')
       !w%rho_pol=f_malloc_ptr([n1,n23],id='rho_pol') !>>>>>>>>>>here the switch
       w%eps=f_malloc_ptr([n1,n23],id='eps')
       w%dlogeps=f_malloc_ptr([3,ndims(1),ndims(2),ndims(3)],id='dlogeps')
       w%oneoeps=f_malloc_ptr([n1,n23],id='oneoeps')
       w%epsinnersccs=f_malloc_ptr([n1,n23],id='epsinnersccs')
       !w%epsinnersccs=f_malloc_ptr([n1,ndims(2)*ndims(3)],id='epsinnersccs')
    end select

   if (.not. (method .hasattr. PS_PB_NONE_ENUM))&
        w%rho_ions=f_malloc0_ptr([n1,n23],id='rho_ions')

   end subroutine PS_allocate_cavity_workarrays

  !> create the memory space needed to store the arrays for the 
  !! description of the cavity
  subroutine pkernel_allocate_cavity(kernel,vacuum)
    use psolver_environment, only: PS_SCCS_ENUM,vacuum_eps
    use f_utils, only: f_zero
    implicit none
    type(coulomb_operator), intent(inout) :: kernel
    logical, intent(in), optional :: vacuum !<if .true. the cavity is allocated as no cavity exists, i.e. only vacuum
    !local variables
    integer :: n1,n23,i1,i23

    n1=kernel%mesh%ndims(1)
    n23=kernel%mesh%ndims(2)*kernel%grid%n3p
!!$    call PS_allocate_cavity_workarrays(n1,n23,kernel%ndims,&
!!$         kernel%method,kernel%w)
    if (present(vacuum)) then
       if (vacuum) then
          select case(trim(str(kernel%method)))
          case('PCG')
             call f_zero(kernel%w%corr)
             do i23=1,n23
                do i1=1,n1
                   kernel%w%eps(i1,i23)=vacuum_eps
                   kernel%w%oneoeps(i1,i23)=1.0_dp/sqrt(vacuum_eps)
                end do
             end do
          case('PI')
             call f_zero(kernel%w%dlogeps)
             do i23=1,n23
                do i1=1,n1
                   kernel%w%eps(i1,i23)=vacuum_eps
                   kernel%w%oneoeps(i1,i23)=1.0_dp/vacuum_eps
                end do
             end do
          end select
          if (kernel%method .hasattr. PS_SCCS_ENUM) &
               call f_zero(kernel%w%epsinnersccs)
       end if
    end if

  end subroutine pkernel_allocate_cavity


  !> set the array radii on the basis of the information provided by the user
  function pkernel_get_radius(kernel,atname) result(radius)
    use numerics, only: Bohr_Ang
    use psolver_environment
    implicit none
    !> Poisson Solver kernel
    type(coulomb_operator), intent(inout) :: kernel
    !> name of the atom 
    character(len=*), intent(in) :: atname
    !> radii of each of the atom types, calculated on the basis of the input values
    real(dp)  :: radius

    if (atname .in. kernel%radii_dict) then
       radius=kernel%radii_dict//atname
    else
       select case (kernel%radii_set)
       case(RADII_PAULING_ID)
          radius = radii_Pau(atname)
       case(RADII_BONDI_ID)
          radius = radii_Bondi(atname)
       case(RADII_UFF_ID)
          radius = radii_UFF(atname)
       end select
    end if

  end function pkernel_get_radius

  !> Set the dielectric cavity epsilon in the pkernel structure for a nonvacuum
  !! treatment. It generate the dielectric cavity epsilon(r) and all working
  !! arrays needed by the PCG or SC solver of the generalized Poisson equation.
  !! There are several methods to do that:
  !! 1. In case of the soft-sphere model we need just to pass some informations 
  !!    of our atomistic system like the total number of atoms (nat), their positions 
  !!    (rxyz(1:3,nat)) and their radii.
  !!      call pkernel_set_epsilon(pkernel,nat=nat,rxyz=rxyz,radii=radii)
  !! 2. If you have a given cavity epsilon(r) (on the same real space grid of the input
  !!    charge) build up with other methods, you can directly set up all working arrays:
  !!      call pkernel_set_epsilon(pkernel,eps=eps)
  !! 3. If you have both the cavity and the working arrays you can pass all of them to
  !!    the solver. In the case of PCG solver you need these three vectors:
  !!      call pkernel_set_epsilon(pkernel,eps=eps,oneosqrteps=oneosqrteps,corr=corr)
  !!    Here eps represent the dielectric cavity epsilon(r), oneosqrteps the inverse of
  !!    its square root, and corr the vector q(r) of equation (16) in Fisicaro J. Chem.
  !!    Phys. 144, 014103 (2016). All these vectors have to be on the same real space
  !!    grid of the input charge.
  !!    In the case of SC solver you need to pass these vectors:
  !!      call pkernel_set_epsilon(pkernel,eps=eps,oneoeps=oneoeps,dlogeps=dlogeps)
  !!    where oneoeps is the inverse of epsilon and dlogeps the derivative of its
  !!    natural logarithm.
  subroutine pkernel_set_epsilon(kernel,eps,dlogeps,oneoeps,oneosqrteps,corr,nat,rxyz,radii)
    use yaml_strings
    use dynamic_memory
    use FDder
    use dictionaries, only: f_err_throw
    use numerics, only: pi
    use psolver_environment, only: rigid_cavity_arrays,vacuum_eps
    use f_utils, only: f_zero
    implicit none
    !> Poisson Solver kernel
    type(coulomb_operator), intent(inout) :: kernel
    !> dielectric function. Needed for non VAC methods, given in full dimensions
    real(dp), dimension(:,:,:), intent(in), optional :: eps
    !> logarithmic derivative of epsilon. Needed for SC method.
    !! if absent, it will be calculated from the array of epsilon
    real(dp), dimension(:,:,:,:), intent(in), optional :: dlogeps
    !> inverse of epsilon. Needed for SC method.
    !! if absent, it will be calculated from the array of epsilon
    real(dp), dimension(:,:,:), intent(in), optional :: oneoeps
    !> inverse square root of epsilon. Needed for PCG method.
    !! if absent, it will be calculated from the array of epsilon
    real(dp), dimension(:,:,:), intent(in), optional :: oneosqrteps
    !> correction term for PCG
    !! if absent, it will be calculated from the array of epsilon
    real(dp), dimension(:,:,:), intent(in), optional :: corr
    !> number of atoms, can be inserted to define the rigid cavity
    integer, intent(in), optional :: nat
    !> if nat is present, also the positions of the atoms have to be defined
    !! (dimension 3,nat). The rxyz values are defined accordingly to the
    !! position of the atoms in the grid starting from the point [1,1,1] to ndims(:)
    real(dp), dimension(*), intent(in), optional :: rxyz
    !> and the radii aroud each atoms also have to be defined (dimension nat)
    real(dp), dimension(*), intent(in), optional :: radii
    !local variables
    logical, dimension(3) :: prst
    integer :: n1,n23,i3s,i23,i3,i2,i1
    real(dp) :: cc,ep,depsr,epsm1,hh,kk
    real(dp), dimension(3) :: v,dleps
    real(dp), dimension(:,:,:), allocatable :: de2,ddeps
    real(dp), dimension(:,:,:,:), allocatable :: deps
    type(cell) :: mesh

    if (present(corr)) then
       !check the dimensions (for the moment no parallelism)
       if (any(shape(corr) /= kernel%mesh%ndims)) &
            call f_err_throw('Error in the dimensions of the array corr,'//&
            trim(yaml_toa(shape(corr))))
    end if
    if (present(eps)) then
       !check the dimensions (for the moment no parallelism)
       if (any(shape(eps) /= kernel%mesh%ndims)) &
            call f_err_throw('Error in the dimensions of the array epsilon,'//&
            trim(yaml_toa(shape(eps))))
    end if
    if (present(oneoeps)) then
       !check the dimensions (for the moment no parallelism)
       if (any(shape(oneoeps) /= kernel%mesh%ndims)) &
            call f_err_throw('Error in the dimensions of the array oneoeps,'//&
            trim(yaml_toa(shape(oneoeps))))
    end if
    if (present(oneosqrteps)) then
       !check the dimensions (for the moment no parallelism)
       if (any(shape(oneosqrteps) /= kernel%mesh%ndims)) &
            call f_err_throw('Error in the dimensions of the array oneosqrteps,'//&
            trim(yaml_toa(shape(oneosqrteps))))
    end if
    if (present(dlogeps)) then
       !check the dimensions (for the moment no parallelism)
       if (any(shape(dlogeps) /= &
            [3,kernel%mesh%ndims(1),kernel%mesh%ndims(2),kernel%mesh%ndims(3)])) &
            call f_err_throw('Error in the dimensions of the array dlogeps,'//&
            trim(yaml_toa(shape(dlogeps))))
    end if

    prst=[present(nat),present(rxyz),present(radii)]
    if (.not. all(prst) .and. any(prst)) &
         call f_err_throw('All rxyz, radii and nat have to be present in the '//&
         'pkernel_set_epsilon routine')

    if (all(prst) .and. .not. (kernel%method .hasattr. 'rigid')) then 
         call f_err_throw('Where rxyz, radii and nat are present in the '//&
         'pkernel_set_epsilon routine the cavity has to be set to "rigid"') 
    else if (all(prst)) then
       kernel%w%nat=nat
       kernel%w%rxyz=f_malloc_ptr([3,nat],id='rxyz')
       if (nat >0) call f_memcpy(n=3*nat,src=rxyz(1),dest=kernel%w%rxyz(1,1))
       kernel%w%radii=f_malloc_ptr(nat,id='radii')
       if (nat >0) call f_memcpy(n=nat,src=radii(1),dest=kernel%w%radii(1))
    end if

    !store the arrays needed for the method
    !the stored arrays are of rank two to collapse indices for
    !omp parallelism
    n1=kernel%mesh%ndims(1)
    n23=kernel%mesh%ndims(2)*kernel%grid%n3p
    !starting point in third direction
    i3s=kernel%grid%istart+1
    if (kernel%grid%n3p==0) i3s=1

    select case(trim(str(kernel%method)))
    case('PCG')
       !check the dimensions of the associated arrays
       if (all([associated(kernel%w%corr),associated(kernel%w%oneoeps)])) then
          !then check the shapes
          if (any(shape(kernel%w%oneoeps) /= [n1,n23])) &
               call f_err_throw('Incorrect shape of oneoeps')
          if (any(shape(kernel%w%corr) /= [n1,n23])) &
               call f_err_throw('Incorrect shape of corr')
          if (present(corr)) then
             call f_memcpy(n=n1*n23,src=corr(1,1,i3s),dest=kernel%w%corr)
          else if (present(eps)) then
        !allocate work arrays
             deps=f_malloc([kernel%mesh%ndims(1),kernel%mesh%ndims(2),kernel%mesh%ndims(3),3],id='deps')
             de2 =f_malloc(kernel%mesh%ndims,id='de2')
             ddeps=f_malloc(kernel%mesh%ndims,id='ddeps')

             call nabla_u_and_square(kernel%mesh,eps,deps,de2,kernel%nord)

             call div_u_i(kernel%mesh,deps,ddeps,kernel%nord)
             i23=1
             do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
                do i2=1,kernel%mesh%ndims(2)
                   do i1=1,kernel%mesh%ndims(1)
                      kernel%w%corr(i1,i23)=(-0.125d0/pi)*&
                           (0.5d0*de2(i1,i2,i3)/eps(i1,i2,i3)-ddeps(i1,i2,i3))
                   end do
                   i23=i23+1
                end do
             end do
             call f_free(deps)
             call f_free(ddeps)
             call f_free(de2)
          else if (all(prst)) then
             mesh=kernel%mesh !cell_new(kernel%geocode,kernel%mesh%ndims,kernel%mesh%hgrids)
             epsm1=(kernel%cavity%epsilon0-vacuum_eps)
             call f_zero(kernel%IntSur)
             call f_zero(kernel%IntVol)
             hh=mesh%volume_element
             do i3=i3s,kernel%grid%n3p+i3s-1
                v(3)=cell_r(mesh,i3,dim=3)
                do i2=1,mesh%ndims(2)
                   v(2)=cell_r(mesh,i2,dim=2)
                   i23=i2+(i3-i3s)*mesh%ndims(2)
                   do i1=1,mesh%ndims(1)
                      v(1)=cell_r(mesh,i1,dim=1)
                      call rigid_cavity_arrays(kernel%cavity,mesh,v,&
                           kernel%w%nat,kernel%w%rxyz,kernel%w%radii,ep,depsr,dleps,cc,kk)
                      kernel%w%eps(i1,i23)=ep
                      kernel%w%corr(i1,i23)=cc
                      kernel%IntVol=kernel%IntVol+(kernel%cavity%epsilon0-ep)
                      kernel%IntSur=kernel%IntSur+depsr
                   end do
                end do
             end do
             kernel%IntVol=kernel%IntVol*hh/epsm1
             kernel%IntSur=kernel%IntSur*hh/epsm1
          else
             call f_err_throw('For method "PCG" the arrays corr or epsilon should be present')   
          end if
          if (present(oneosqrteps)) then
             call f_memcpy(n=n1*n23,src=oneosqrteps(1,1,i3s),&
                  dest=kernel%w%oneoeps)
          else if (present(eps)) then
             i23=1
             do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
                do i2=1,kernel%mesh%ndims(2)
                   do i1=1,kernel%mesh%ndims(1)
                      kernel%w%oneoeps(i1,i23)=1.0_dp/sqrt(eps(i1,i2,i3))
                   end do
                   i23=i23+1
                end do
             end do
          else if (all(prst)) then
             do i23=1,n23
                do i1=1,n1
                   ep=kernel%w%eps(i1,i23)
                   kernel%w%oneoeps(i1,i23)=1.d0/sqrt(ep) !here square root is correct
                end do
             end do
          else
             call f_err_throw('For method "PCG" the arrays oneosqrteps or epsilon should be present')
          end if
          if (present(eps)) then
             call f_memcpy(n=n1*n23,src=eps(1,1,i3s),&
                  dest=kernel%w%eps)
          else if (.not. all(prst)) then
             call f_err_throw('For method "PCG" the arrays eps should be present')
          end if
       else
          call f_err_throw('For method "PCG" the arrays oneosqrteps'//&
               ' and corr have to be associated, call PS_allocate_cavity_workarrays')
       end if
    case('PI')
       if (all([associated(kernel%w%dlogeps),associated(kernel%w%oneoeps)])) then
          !then check the shapes
          if (any(shape(kernel%w%oneoeps) /= [n1,n23])) &
               call f_err_throw('Incorrect shape of oneoeps')
          if (any(shape(kernel%w%dlogeps) /= [3,kernel%mesh%ndims(1),kernel%mesh%ndims(2),kernel%mesh%ndims(3)])) &
               call f_err_throw('Incorrect shape of dlogeps')
          if (present(dlogeps)) then
             call f_memcpy(src=dlogeps,dest=kernel%w%dlogeps)
          else if (present(eps)) then
             !allocate arrays
             deps=f_malloc([kernel%mesh%ndims(1),kernel%mesh%ndims(2),kernel%mesh%ndims(3),3],id='deps')
             call nabla_u(kernel%mesh,eps,deps,kernel%nord)
             do i3=1,kernel%mesh%ndims(3)
                do i2=1,kernel%mesh%ndims(2)
                   do i1=1,kernel%mesh%ndims(1)
                      !switch and create the logarithmic derivative of epsilon
                      kernel%w%dlogeps(1,i1,i2,i3)=deps(i1,i2,i3,1)/eps(i1,i2,i3)
                      kernel%w%dlogeps(2,i1,i2,i3)=deps(i1,i2,i3,2)/eps(i1,i2,i3)
                      kernel%w%dlogeps(3,i1,i2,i3)=deps(i1,i2,i3,3)/eps(i1,i2,i3)
                   end do
                end do
             end do
             call f_free(deps)
          else if (all(prst)) then
             mesh=kernel%mesh !cell_new(kernel%geocode,kernel%mesh%ndims,kernel%mesh%hgrids)
             epsm1=(kernel%cavity%epsilon0-vacuum_eps)
             call f_zero(kernel%IntSur)
             call f_zero(kernel%IntVol)
             hh=mesh%volume_element

             do i3=1,kernel%mesh%ndims(3)
                v(3)=cell_r(mesh,i3,dim=3)
                do i2=1,mesh%ndims(2)
                   v(2)=cell_r(mesh,i2,dim=2)
                   i23=i2+(i3-i3s)*mesh%ndims(2)
                   do i1=1,mesh%ndims(1)
                      v(1)=cell_r(mesh,i1,dim=1)
                      call rigid_cavity_arrays(kernel%cavity,mesh,v,kernel%w%nat,&
                           kernel%w%rxyz,kernel%w%radii,ep,depsr,dleps,cc,kk)
                      if (i23 <= n23 .and. i23 >=1) then
                         kernel%w%eps(i1,i23)=ep
                         kernel%IntVol=kernel%IntVol+(kernel%cavity%epsilon0-ep)
                         kernel%IntSur=kernel%IntSur+depsr
                      end if
                      kernel%w%dlogeps(:,i1,i2,i3)=dleps
                   end do
                end do
             end do
             kernel%IntVol=kernel%IntVol*hh/epsm1
             kernel%IntSur=kernel%IntSur*hh/epsm1
          else
             call f_err_throw('For method "PI" the arrays dlogeps or epsilon should be present')
          end if
          if (present(oneoeps)) then
             call f_memcpy(n=n1*n23,src=oneoeps(1,1,i3s),&
                  dest=kernel%w%oneoeps)
          else if (present(eps)) then
             i23=1
             do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
                do i2=1,kernel%mesh%ndims(2)
                   do i1=1,kernel%mesh%ndims(1)
                      kernel%w%oneoeps(i1,i23)=1.0_dp/eps(i1,i2,i3)
                   end do
                   i23=i23+1
                end do
             end do
             else if (all(prst)) then
                do i23=1,n23
                   do i1=1,n1
                      ep=kernel%w%eps(i1,i23)
                      kernel%w%oneoeps(i1,i23)=1.d0/ep
                   end do
                end do
          else
             call f_err_throw('For method "PI" the arrays oneoeps or epsilon should be present')
          end if
          if (present(eps)) then
             call f_memcpy(n=n1*n23,src=eps(1,1,i3s),&
                  dest=kernel%w%eps)
          else if (.not. all(prst)) then
             call f_err_throw('For method "PCG" the arrays eps should be present')
          end if
       else
          call f_err_throw('For method "PI" the arrays oneoeps '//&
               'and dlogeps have to be associated, call PS_allocate_cavity_workarrays')
       end if
    end select

  end subroutine pkernel_set_epsilon

  !>calculate the extra contribution to the forces and the cavitation terms
  !!to be given at the end of the calculation
  subroutine ps_soft_PCM_forces(kernel,fpcm)
    use psolver_environment
    use f_utils, only: f_zero
    implicit none
    type(coulomb_operator), intent(in) :: kernel
    real(dp), dimension(3,kernel%w%nat), intent(inout) :: fpcm
    !local variables
    real(dp), parameter :: thr=1.e-12
    integer :: i1,i2,i3,i23
    real(dp) :: cc,epr,depsr,tt,kk
    type(cell) :: mesh
    real(dp), dimension(3) :: v,dleps,deps
    mesh=kernel%mesh !cell_new(kernel%geocode,kernel%mesh%ndims,kernel%mesh%hgrids)

    do i3=1,kernel%grid%n3p
       v(3)=cell_r(mesh,i3+kernel%grid%istart,dim=3)
       do i2=1,kernel%mesh%ndims(2)
          v(2)=cell_r(mesh,i2,dim=2)
          i23=i2+(i3-1)*kernel%mesh%ndims(2)
          do i1=1,kernel%mesh%ndims(1)
             tt=kernel%w%oneoeps(i1,i23) !nablapot2(r)
             v(1)=cell_r(mesh,i1,dim=1)
             !this is done to obtain the depsilon
             call rigid_cavity_arrays(kernel%cavity,mesh,v,kernel%w%nat,&
                  kernel%w%rxyz,kernel%w%radii,epr,depsr,dleps,cc,kk)
             if (abs(epr-vacuum_eps) < thr) cycle
             deps=dleps*epr
             call rigid_cavity_forces(kernel%opt%only_electrostatic,kernel%cavity,mesh,v,&
                  kernel%w%nat,kernel%w%rxyz,kernel%w%radii,epr,tt,fpcm,deps,kk)
          end do
       end do
    end do

  end subroutine ps_soft_PCM_forces


  subroutine build_cavity_from_rho(rho,nabla2_rho,delta_rho,cc_rho,kernel,&
       depsdrho,dsurfdrho,IntSur,IntVol)
    use psolver_environment
    implicit none
    type(coulomb_operator), intent(inout) :: kernel
    real(dp), dimension(kernel%mesh%ndims(1),kernel%mesh%ndims(2)*kernel%grid%n3p), intent(in) :: rho,nabla2_rho,delta_rho,cc_rho
    !> functional derivative of the sc epsilon with respect to 
    !! the electronic density, in distributed memory
    real(dp), dimension(kernel%mesh%ndims(1),kernel%mesh%ndims(2)*kernel%grid%n3p), intent(out) :: depsdrho
    !> functional derivative of the surface integral with respect to 
    !! the electronic density, in distributed memory
    real(dp), dimension(kernel%mesh%ndims(1),kernel%mesh%ndims(2)*kernel%grid%n3p), intent(out) :: dsurfdrho
    real(dp), intent(out) :: IntSur,IntVol
    !local variables
    real(dp), parameter :: innervalue = 0.9d0 !to be defined differently
    integer :: n01,n02,n03,i3s,i1,i2,i3,i23
    real(dp) :: rh,d2,d,dd,de,epsm1

    IntSur=0.d0
    IntVol=0.d0

    n01=kernel%mesh%ndims(1)
    n02=kernel%mesh%ndims(2)
    n03=kernel%mesh%ndims(3)
    !starting point in third direction
    i3s=kernel%grid%istart+1
    epsm1=(kernel%cavity%epsilon0-vacuum_eps)
    !now fill the pkernel arrays according the the chosen method
    select case(trim(str(kernel%method)))
    case('PCG')
       !in PCG we only need corr, oneosqrtepsilon
       i23=1
       do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
          !do i3=1,n03
          do i2=1,n02
             do i1=1,n01
                if (kernel%w%epsinnersccs(i1,i23).gt.innervalue) then 
                   kernel%w%eps(i1,i23)=1.d0 !eps(i1,i2,i3)
                   kernel%w%oneoeps(i1,i23)=1.d0 !oneosqrteps(i1,i2,i3)
                   kernel%w%corr(i1,i23)=0.d0 !corr(i1,i2,i3)
                   depsdrho(i1,i23)=0.d0
                   dsurfdrho(i1,i23)=0.d0
                else
                   rh=rho(i1,i23)
                   d2=nabla2_rho(i1,i23)
                   d=sqrt(d2)
                   dd = delta_rho(i1,i23)
                   de=epsprime(rh,kernel%cavity)
                   depsdrho(i1,i23)=de
                   kernel%w%eps(i1,i23)=eps(rh,kernel%cavity)
                   kernel%w%oneoeps(i1,i23)=oneosqrteps(rh,kernel%cavity)
                   kernel%w%corr(i1,i23)=corr_term(rh,d2,dd,kernel%cavity)
                   dsurfdrho(i1,i23)=-surf_term(rh,d2,dd,cc_rho(i1,i23),kernel%cavity)/epsm1
                   !evaluate surfaces and volume integrals
                   IntSur=IntSur - de*d
                   IntVol=IntVol + (kernel%cavity%epsilon0-eps(rh,kernel%cavity))
                end if
             end do
             i23=i23+1
          end do
       end do
    case('PI')
       !for PI we need  dlogeps,oneoeps
       !first oneovereps
       i23=1
       do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
          do i2=1,n02
             do i1=1,n01
                if (kernel%w%epsinnersccs(i1,i23).gt.innervalue) then ! Check for inner sccs cavity value to fix as vacuum
                   kernel%w%eps(i1,i23)=1.d0 !eps(i1,i2,i3)
                   kernel%w%oneoeps(i1,i23)=1.d0 !oneoeps(i1,i2,i3)
                   depsdrho(i1,i23)=0.d0
                   dsurfdrho(i1,i23)=0.d0
                else
                   rh=rho(i1,i23)
                   d2=nabla2_rho(i1,i23)
                   d=sqrt(d2)
                   dd = delta_rho(i1,i23)
                   de=epsprime(rh,kernel%cavity)
                   depsdrho(i1,i23)=de
                   kernel%w%eps(i1,i23)=eps(rh,kernel%cavity)
                   kernel%w%oneoeps(i1,i23)=oneoeps(rh,kernel%cavity) 
                   dsurfdrho(i1,i23)=-surf_term(rh,d2,dd,cc_rho(i1,i23),kernel%cavity)/epsm1

                   !evaluate surfaces and volume integrals
                   IntSur=IntSur - de*d
                   IntVol=IntVol + (kernel%cavity%epsilon0-eps(rh,kernel%cavity))
                end if
             end do
             i23=i23+1
          end do
       end do
    end select

    !IntSur=IntSur*product(kernel%mesh%hgrids)/epsm1
    !IntVol=IntVol*product(kernel%mesh%hgrids)/epsm1
    IntSur=IntSur*kernel%mesh%volume_element/epsm1
    IntVol=IntVol*kernel%mesh%volume_element/epsm1
     

  end subroutine build_cavity_from_rho


!!!  !> calculates the inner cavity vector epsinnersccs for sccs run 
!!!  !! given a set of centres. Based on error function.
!!!  !! Need the radius of the cavit and its smoothness
!!!  subroutine epsinnersccs_rigid_cavity_error_multiatoms_bc(geocode,ndims,hgrids,natreal,rxyzreal,radiireal,delta,eps)
!!!    use f_utils
!!!    use numerics, only : Bohr_Ang
!!!    use module_base, only: bigdft_mpi
!!!    use f_enums
!!!    use yaml_output
!!!    use dynamic_memory
!!!    use bounds, only: ext_buffers
!!!    use psolver_environment, only: epsl,d1eps,epsle0
!!!    implicit none
!!!!    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
!!!    integer, intent(in) :: natreal !< number of centres defining the cavity
!!!    real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
!!!!    integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
!!!!    real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
!!!    real(kind=8), dimension(natreal), intent(in) :: radiireal !< radii of each of the atoms
!!!    !> position of all the atoms in the grid coordinates
!!!    real(kind=8), dimension(3,natreal), intent(in) :: rxyzreal
!!!    real(kind=8), dimension(bit%mesh%ndim), intent(out) :: eps !< dielectric inner cavity for sccs run
!!!
!!!    !local variables
!!!    logical :: perx,pery,perz
!!!    integer :: i,i1,i2,i3,iat,jat,ii,j,k,l,px,py,pz,unt
!!!    integer :: nbl1,nbl2,nbl3,nbr1,nbr2,nbr3,imin
!!!    real(kind=8) :: r2,x,y2,z2,d,d2,d12,y,z,eps_min,eps1,pi,de2,dde,d1,oneod,h,coeff,dmin,dmax
!!!    real(kind=8) :: value,valuemin
!!!    real(kind=8), dimension(3) :: v,rv,shift,sh
!!!
!!!    real(kind=8), parameter :: valuebc=1.d0
!!!    real(kind=8), dimension(6) :: plandist
!!!    integer, dimension(6) :: ba
!!!    real(kind=8), dimension(:), allocatable :: ep
!!!    real(kind=8), dimension(3,27*natreal) :: rxyztot
!!!    real(kind=8), dimension(27*natreal) :: radiitot
!!!    real(kind=8), dimension(:), allocatable :: radii
!!!    real(kind=8), dimension(:,:), allocatable :: rxyz
!!!    logical, parameter :: dumpeps=.false.
!!!
!!!    !buffers associated to the geocode
!!!    !conditions for periodicity in the three directions
!!!    perx=(geocode /= 'F')
!!!    pery=(geocode == 'P')
!!!    perz=(geocode /= 'F')
!!!
!!!    call ext_buffers(perx,nbl1,nbr1)
!!!    call ext_buffers(pery,nbl2,nbr2)
!!!    call ext_buffers(perz,nbl3,nbr3)
!!!
!!!    pi = 4.d0*datan(1.d0)
!!!
!!!    shift(1)=hgrids(1)*ndims(1)
!!!    shift(2)=hgrids(2)*ndims(2)
!!!    shift(3)=hgrids(3)*ndims(3)
!!!
!!!    !------------------------------------------------------------------------------------------------------
!!!    ! Depending of Free, Periodic or Surface bc, image atoms are or not included.
!!!
!!!    !  if (bigdft_mpi%iproc==0) then
!!!    !   do iat=1,natreal
!!!    !    call yaml_map('real input atoms',iat)
!!!    !    call yaml_map('radii',radiireal(iat))
!!!    !    call yaml_map('rxyz',rxyzreal(:,iat))
!!!    !   end do
!!!    !  end if
!!!
!!!    px=0
!!!    py=0
!!!    pz=0
!!!    if (perx) px=1
!!!    if (pery) py=1
!!!    if (perz) pz=1
!!!
!!!    rxyztot(:,:)=0.d0
!!!
!!!    i=0
!!!    do iat=1,natreal
!!!       ba(1:6)=0
!!!       ! checking what are the image atoms to include in the calculation of the
!!!       ! cavity.
!!!       rv(1:3)=rxyzreal(1:3,iat)
!!!       plandist(1)=dabs(rv(1))
!!!       plandist(2)=dabs(shift(1)-rv(1))
!!!       plandist(3)=dabs(rv(2))
!!!       plandist(4)=dabs(shift(2)-rv(2))
!!!       plandist(5)=dabs(rv(3))
!!!       plandist(6)=dabs(shift(3)-rv(3))
!!!       do ii=1,6
!!!          valuemin=1.d0
!!!          d=plandist(ii)
!!!          value=epsl(d,radiireal(iat),delta)
!!!          if (value.lt.valuebc) then ! valuebc is the value to check on the box border to accept or refuse an image atom.
!!!             if (abs(value).lt.valuemin) then
!!!                valuemin=abs(value)
!!!                imin=ii
!!!             end if
!!!             select case(ii)
!!!             case (1)
!!!                ba(1)=1*px
!!!             case (2)
!!!                ba(2)=-1*px
!!!             case (3)
!!!                ba(3)=1*py
!!!             case (4)
!!!                ba(4)=-1*py
!!!             case (5)
!!!                ba(5)=1*pz
!!!             case (6)
!!!                ba(6)=-1*pz
!!!             end select
!!!          end if
!!!       end do
!!!
!!!       do j=ba(6),ba(5)
!!!          sh(3)=real(j,kind=8)*shift(3)
!!!          do k=ba(4),ba(3)
!!!             sh(2)=real(k,kind=8)*shift(2)
!!!             do l=ba(2),ba(1)
!!!                sh(1)=real(l,kind=8)*shift(1)
!!!                rv(1:3)=rxyzreal(1:3,iat) + sh(1:3)
!!!                i=i+1
!!!                rxyztot(1:3,i)=rv(1:3)
!!!                radiitot(i)=radiireal(iat)
!!!             end do
!!!          end do
!!!       end do
!!!
!!!    end do
!!!
!!!    nat=i
!!!
!!!    ep=f_malloc(nat,id='ep')
!!!    rxyz=f_malloc([3,nat],id='rxyz')
!!!    radii=f_malloc(nat,id='radii')
!!!
!!!    rxyz(1:3,1:nat)=rxyztot(1:3,1:nat)
!!!    radii(1:nat)=radiitot(1:nat)
!!!
!!!    !   if (bigdft_mpi%iproc==0) then
!!!    !    write(*,*)plandist
!!!    !    write(*,'(1x,a,1x,e14.7,1x,a,1x,i4)')'Value min =',valuemin,'at bc side',imin
!!!    !    call yaml_map('nat',nat)
!!!    !    do iat=1,nat
!!!    !     call yaml_map('atom',iat)
!!!    !     call yaml_map('radii',radii(iat))
!!!    !     call yaml_map('rxyz',rxyz(:,iat))
!!!    !    end do
!!!    !   end if
!!!
!!!    !------------------------------------------------------------------------------------------------------
!!!    ! Starting the cavity building for rxyztot atoms=real+image atoms (total natcurr) for periodic
!!!    ! and surface boundary conditions or atoms=real for free bc.
!!!
!!!    do i3=1,ndims(3)
!!!       z=hgrids(3)*(i3-1-nbl3)
!!!       v(3)=z
!!!       do i2=1,ndims(2)
!!!          y=hgrids(2)*(i2-1-nbl2)
!!!          v(2)=y
!!!          do i1=1,ndims(1)
!!!             x=hgrids(1)*(i1-1-nbl1)
!!!             v(1)=x
!!!
!!!             do iat=1,nat
!!!                d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
!!!                d=dsqrt(d2)
!!!
!!!                if (d2.eq.0.d0) then
!!!                   d2=1.0d-30
!!!                   ep(iat)=epsl(d,radii(iat),delta)
!!!                else
!!!                   ep(iat)=epsl(d,radii(iat),delta)
!!!                end if
!!!             end do
!!!
!!!             eps(i1,i2,i3)= 1.d0 - product(ep)
!!!
!!!          end do
!!!       end do
!!!    end do
!!!
!!!
!!!!!$    if (dumpeps) then
!!!!!$
!!!!!$       unt=f_get_free_unit(20)
!!!!!$       call f_open_file(unt,file='epsinnersccs.dat')
!!!!!$       i1=1!n03/2
!!!!!$       do i2=1,ndims(2)
!!!!!$          do i3=1,ndims(3)
!!!!!$             write(unt,'(2(1x,I4),2(1x,e14.7))')i2,i3,eps(i1,i2,i3),eps(ndims(1)/2,i2,i3)
!!!!!$          end do
!!!!!$       end do
!!!!!$       call f_close(unt)
!!!!!$
!!!!!$       unt=f_get_free_unit(22)
!!!!!$       call f_open_file(unt,file='epsinnersccs_line_y.dat')
!!!!!$       do i2=1,ndims(2)
!!!!!$          write(unt,'(1x,I8,1(1x,e22.15))')i2,eps(ndims(1)/2,i2,ndims(3)/2)
!!!!!$       end do
!!!!!$       call f_close(unt)
!!!!!$
!!!!!$       unt=f_get_free_unit(23)
!!!!!$       call f_open_file(unt,file='epsinnersccs_line_z.dat')
!!!!!$       do i3=1,ndims(3)
!!!!!$          write(unt,'(1x,I8,1(1x,e22.15))')i3,eps(ndims(1)/2,ndims(2)/2,i3)
!!!!!$       end do
!!!!!$       call f_close(unt)
!!!!!$
!!!!!$    end if
!!!
!!!    call f_free(ep)
!!!    call f_free(rxyz)
!!!    call f_free(radii)
!!!
!!!  end subroutine epsinnersccs_rigid_cavity_error_multiatoms_bc
!!!!!!


  subroutine epsilon_inner_cavity(pkernel,nat,rxyz,radii,delta,oxyz)
    use box
    use psolver_environment
    implicit none
    integer, intent(in) :: nat
    real(dp), intent(in) :: delta
    type(coulomb_operator), intent(inout) :: pkernel
    real(dp), dimension(nat), intent(in) :: radii
    real(dp), dimension(3,nat), intent(in) :: rxyz,oxyz
    !local variables
    integer :: iat
    real(dp) :: tt,d
    type(box_iterator) :: bit

    bit=box_iter(pkernel%mesh,i3s=pkernel%grid%istart+1,n3p=pkernel%grid%n3p,origin=oxyz)
    do while(box_next_point(bit))
       !loop on atoms
       tt=1.0_dp
       do iat=1,nat
          d=distance(bit%mesh,rxyz(1,iat),bit%rxyz)
          tt=tt*epsl(d,radii(iat),delta)
       end do
       pkernel%w%epsinnersccs(bit%i,bit%i23+1)=1.0_dp-tt
    end do

  end subroutine epsilon_inner_cavity


end module PStypes
