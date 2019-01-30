!> This file aims to show the main features/routines of psolver
!! and how them can be called from external programs.
!! @author
!!    Copyright (C) 2002-2017 BigDFT group  (Giuseppe Fisicaro)<br/>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

program PSolver_examples

   use wrapper_mpi
   use Poisson_Solver
   use PSbox
   use yaml_output
   use dynamic_memory
   use dictionaries, dict_set => set
   use time_profiling
   use f_utils
   use yaml_strings
   use box
   use PSbase
   use PStypes, only: build_cavity_from_rho,PS_allocate_cavity_workarrays
   use numerics
   use f_blas, only: f_dot
   use at_domain
   implicit none
   type(coulomb_operator) :: pkernel !< @copydoc poisson_solver::coulomb_operator
   type(PSolver_energies) :: energies !< @copydoc poisson_solver::coulomb_operator::energies
   type(dictionary), pointer :: options,dict_input
   character(len=2) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
   type(cell) :: mesh !< @copydoc poisson_solver::coulomb_operator::mesh
   character(len=2), parameter :: datacode = 'G'  !< @copydoc poisson_solver::coulomb_operator::datacode
   integer, dimension(3) :: ndims !< @copydoc poisson_solver::coulomb_operator::ndims
   real(8), dimension(3) :: hgrids !< @copydoc poisson_solver::coulomb_operator::hgrids
   real(8), dimension(3) :: angdeg,angrad ! angles for non-othorhombic cell 
   !> Set the method to solve the generalized Poisson equation. Allowed values:
   !! VAC -> standart vacuum solver of bigdft
   !! PCG -> preconditioned conjugate gradient
   !! PI  -> self-consistent approach
   !! L. Genovese, T. Deutsch, A. Neelov, S. Goedecker, and G. Beylkin, J. Chem.
   !! Phys. 125, 074105 (2006)
   !! see G. Fisicaro, L. Genovese, O. Andreussi, N. Marzari, and S. Goedecker, 
   !! J. Chem. Phys. 144, 014103 (2016)
   character(len=4) :: PSol
   !> Set the approach to build up the dielectric cavity. Allowed values:
   !! soft-sphere -> soft-sphere cavity
   !! sccs -> charge-dependent cavity
   !! see G. Fisicaro, L. Genovese, O. Andreussi, S. Mandal, N. N. Nair, N.
   !! Marzari, and S. Goedecker, J. Chem. Theory Comput. 13, 3829 (2017) 
   character(len=11) :: cav
   logical :: usegpu
   !> Dummy box dimension takes equal in each direction
   real(kind=8), parameter :: acell = 10.d0 
   !> Dummy number of atoms to build the soft-sphere cavity
   integer :: nat = 1
   !> Dummy van der Waals radius for the atoms, only for the soft-sphere cavity.
   real(kind=8), parameter :: rad_cav = 2.7d0 
   !> Dummy extension for the transition region, only for the soft-sphere cavity.
   real(kind=8) :: delta 
   !> Parameters for analytical input gaussian charge density
   real(kind=8), parameter :: a_gauss = 1.0d0,a2 = a_gauss**2
   !> local variables
   integer :: iproc,nproc,ixc,n01,n02,n03,n1,n23,i1,i2,i3
   real(kind=8) :: hx,hy,hz,einit,IntSur,IntVol,offset,alpha,beta,gamma
   logical :: logyes
   integer :: igpu
   real(kind=8), dimension(:,:), allocatable :: rxyz
   real(kind=8), dimension(:), allocatable :: radii
   real(kind=8), dimension(:,:,:), allocatable :: density,potential,rho,rhopot
   real(kind=8), dimension(:,:,:), allocatable :: delta_rho,cc_rho,nabla2_rhopot
   real(kind=8), dimension(:,:,:,:), allocatable :: nabla_rho
   real(kind=8), dimension(:,:), allocatable :: depsdrho,dsurfdrho
   type(domain) :: dom
 
!GPEx------------------------------------------------------------------------
! How to run:
! Please type "./PSolver_examples --help" to see all possible entries
! Command line example to run:
! vacuum run  -> "./PSolver_examples -n 100 -m VAC -l no"
! soft-sphere with PCG -> "./PSolver_examples -n 100 -m PCG -l no -c soft-sphere"
! soft-sphere with PI -> "./PSolver_examples -n 100 -m PI -l no -c soft-sphere"
! sccs with PCG -> "./PSolver_examples -n 100 -m PCG -l no -c sccs"
! sccs with PI -> "./PSolver_examples -n 100 -m PI -l no -c sccs"
! Try all these run to compare number of iterations needed to solve the
! generalized Poisson equation with different methods and cavities.

!GPEx------------------------------------------------------------------------
! Reading of the input variables from the command line
   call f_lib_initialize()

   !read command line
   call PS_Check_command_line_options(options)

   call f_zero(PSol)
   call f_zero(cav)
   PSol=options .get. 'method'
   cav=options .get. 'cav'
   if (len_trim(PSol)==0) call f_strcpy(src='VAC',dest=PSol)
   if (len_trim(cav)==0) call f_strcpy(src='VAC',dest=cav)
   ndims=options // 'ndim'
   geocode=options//'geocode'
!   SetEps =options//'seteps'
   usegpu = options // 'accel'
   logyes= options // 'logfile'
   angdeg=options // 'angdeg'
   delta= options .get. 'deltacav'
   call f_zero(einit)

   call dict_init(dict_input)

   if ('input' .in. options) &
        call dict_copy(dest=dict_input,src=options//'input')

   call dict_free(options)

   igpu=0
   if (usegpu) igpu=1
   ! number of grid points in each direction, taken from the input ndims
   n01=ndims(1)
   n02=ndims(2)
   n03=ndims(3)
   ! real space mesh in each direction
   hx=acell/real(n01,kind=8)
   hy=acell/real(n02,kind=8)
   hz=acell/real(n03,kind=8)
   hgrids=(/hx,hy,hz/)

   ! Set the angles in radiant, taken from the input angdeg 
   alpha = angdeg(1)/180.0_f_double*pi
   beta  = angdeg(2)/180.0_f_double*pi
   gamma = angdeg(3)/180.0_f_double*pi
   angrad(1) = angdeg(1)/180.0_f_double*pi
   angrad(2) = angdeg(2)/180.0_f_double*pi
   angrad(3) = angdeg(3)/180.0_f_double*pi
   ! Set the mesh type, which contains all informations of the simulation cell
   ! (i.e. number of grid points "ndims", box mesh "hgrids", metric for
   ! non-orthorhombic, etc ...).  
   !mesh=cell_new(geocode,ndims,hgrids,alpha_bc=angrad(1),beta_ac=angrad(2),gamma_ab=angrad(3)) 
   dom=domain_new(units=ATOMIC_UNITS,bc=geocode_to_bc_enum(geocode),&
            alpha_bc=angrad(1),beta_ac=angrad(2),gamma_ab=angrad(3),acell=ndims*hgrids)
   mesh=cell_new(dom,ndims,hgrids)

   call mpiinit()
   iproc=mpirank()
   nproc=mpisize()

   !control memory profiling
   call f_malloc_set_status(iproc=iproc)
   if (iproc ==0) then
        if (logyes) then
         call yaml_set_stream(record_length=92,tabbing=30,unit=70,filename='log.yaml',position='rewind')
        else
         call yaml_set_stream(record_length=92,tabbing=30)
        end if
      call yaml_new_document()
   end if

   density=f_malloc(ndims,id='density')
   rhopot=f_malloc(ndims,id='rhopot')
   rxyz   =f_malloc([3,nat],id='rxyz')
   radii   =f_malloc([nat],id='radii')
   potential=f_malloc(ndims,id='potential')


!GPEx------------------------------------------------------------------------
! Here we define our system: number of atoms (nat), their coordinates (rxyz)
    if (nat.eq.1) then
     rxyz(1,1) = hx*real((ndims(1)-1)/2,kind=8)
     rxyz(2,1) = hy*real((ndims(2)-1)/2,kind=8)
     rxyz(3,1) = hz*real((ndims(3)-1)/2,kind=8)
    else if (nat.eq.2) then
     rxyz(1,1) = hx*real((ndims(1)-1)/2,kind=8)
     rxyz(2,1) = hy*real((ndims(2)-1)/2,kind=8)
     rxyz(3,1) = hz*real((ndims(3)-1)/2,kind=8)
     rxyz(1,2) = hx*real((ndims(1)-1)/2,kind=8)
     rxyz(2,2) = hy*real((ndims(2)-1)/3,kind=8)
     rxyz(3,2) = hz*real((ndims(3)-1)/2,kind=8)
    else if (nat.eq.3) then
     rxyz(1,1) = hx*real((ndims(1)-1)/2,kind=8)
     rxyz(2,1) = hy*real((ndims(2)-1)/2,kind=8)
     rxyz(3,1) = hz*real((ndims(3)-1)/2,kind=8)
     rxyz(1,2) = hx*real((ndims(1)-1)/3,kind=8)
     rxyz(2,2) = hy*real((ndims(2)-1)/2,kind=8)
     rxyz(3,2) = hz*real((ndims(3)-1)/2,kind=8)
     rxyz(1,3) = hx*real((ndims(1)-1)/2,kind=8)
     rxyz(2,3) = hy*real((ndims(2)-1)/3,kind=8)
     rxyz(3,3) = hz*real((ndims(3)-1)/2,kind=8)
    end if

!GPEx------------------------------------------------------------------------
! Here we define the input parameters of the dielectric cavity:
!      soft-sphere cavity: van der Waals radii for atoms, delta for the transition region
!      sccs charge-dependent cavity: rhomin, rhomax
    if (nat.eq.1) then
     radii(1)=rad_cav
    else if (nat.eq.2) then
     radii=[rad_cav,rad_cav]
    else if (nat.eq.3) then
     radii=[rad_cav,rad_cav,rad_cav]
    end if
    delta=2.0d0

    if (usegpu) call dict_set(dict_input//'setup'//'accel','CUDA')
    call dict_set(dict_input//'environment'//'delta',delta)
    if (trim(PSol) /= 'VAC') then
     select case(trim(cav))
      case('soft-sphere')
       call dict_set(dict_input//'environment'//'cavity','rigid')
      case('sccs')
       call dict_set(dict_input//'environment'//'cavity','sccs')
      end select     
      call dict_set(dict_input//'environment'//'gps_algorithm',PSol)
    end if
! If you want to modify input parameter for the PCG or PI methods, please charge
! these entries on the pkernel type:
!     !> Order of accuracy for derivatives into ApplyLaplace subroutine = Total number 
!     !! of points at left and right of the x0 where we want to calculate the derivative.
!     integer :: nord
!     !> default set of radii for the rigid cavity
!     integer :: radii_set
!     !> dictionary of the atom species radii defined by the user
!     type(dictionary), pointer :: radii_dict
!     integer :: max_iter    !< maximum number of convergence iterations
!     real(dp) :: minres     !< Convergence criterion for the iteration
!     real(dp) :: PI_eta     !< Parameter for the update of PI iteration
!     integer :: max_iter_PB !< Max conv iterations for PB treatment
!     real(dp) :: minres_PB  !< Convergence criterion for PB residue
!     real(dp) :: PB_eta     !< Mixing scheme for PB
!     real(dp) :: IntVol     !< Volume integral needed for the non-electrostatic energy contributions
!     real(dp) :: IntSur     !< Surface integral needed for the non-electrostatic energy contributions

! If you want to modify input parameter for the dielectric cavity (soft-sphere
! or sccs) please charge these entries on the cavity_data type, which is within
! pkernel
!  !> define the cavity type
!  type, public :: cavity_data
!     real(gp) :: epsilon0 !< dielectriv constant of the medium
!     real(gp) :: edensmax !<maximum value of the density for the cavity
!     real(gp) :: edensmin !<minimum  value of the density for the cavity
!     real(gp) :: delta !<parameter for the PCM cavity in the case of rigid
!     real(gp) :: fact_rigid !<multiplying factor for the whole PCM cavity in the case of rigid
!     real(dp) :: gammaS !< surface tension of the solvent [dyn/cm]
!     real(dp) :: alphaS !< proportionality factor for the repulsion free energy in term of the surface integral [dyn/cm]
!     real(dp) :: betaV !<proportionality factor for the dispersion free energy in term of the volume integral [GPa]
!  end type cavity_data


!GPEx------------------------------------------------------------------------
! Here we build up a gaussian input charge density (density). 

   call test_functions_new2(mesh,acell,a_gauss,pkernel%mu,density,potential)

   !offset, used only for the periodic solver case
   call set_offset()

!GPEx------------------------------------------------------------------------
! Here we set up the poisson kernel (pkernel) for:
!      1. vacuum calculations; 
!      2. neutral solvent described by the soft-sphere model;
!      3. neutral solvent described by the sccs charge-dependent model.

   pkernel=pkernel_init(iproc,nproc,dict_input,geocode,ndims,hgrids,&
           alpha_bc=alpha,beta_ac=beta,gamma_ab=gamma)
   call dict_free(dict_input)

   call pkernel_set(pkernel,verbose=.true.)

! allocate cavity vectors if needed for nonvacuum treatments
   if ( trim(PSol)/='VAC') then
    select case(trim(cav))
     case('soft-sphere')
      call pkernel_set_epsilon(pkernel,nat=nat,rxyz=rxyz,radii=radii)
     case('sccs')

!     atoms%astruct%nat=nat
!     call epsinnersccs_cavity(atoms,rxyz,denspot%pkernel)
     n1=pkernel%grid%m1
     n23=pkernel%grid%m3*pkernel%grid%n3p
     rho=f_malloc(ndims,id='rho')
     nabla_rho=f_malloc([ndims(1),ndims(2),ndims(3),3],id='nabla_rho')
     nabla2_rhopot=f_malloc([ndims(1),ndims(2),ndims(3)],id='nabla2_rhopot')
     delta_rho=f_malloc(ndims,id='delta_rho')
     cc_rho=f_malloc(ndims,id='cc_rho')
     depsdrho=f_malloc([n1,n23],id='depsdrho')
     dsurfdrho=f_malloc([n1,n23],id='dsurfdrho')

     call f_memcpy(n=product(pkernel%mesh%ndims),src=density(1,1,1),dest=rho(1,1,1))

     call rebuild_cavity_from_rho(rho,nabla_rho,nabla2_rhopot,delta_rho,cc_rho,depsdrho,dsurfdrho,&
          pkernel,IntSur,IntVol)    

     call f_free(rho)
     call f_free(nabla_rho)
     call f_free(nabla2_rhopot)
     call f_free(delta_rho)
     call f_free(cc_rho)
     call f_free(depsdrho)
     call f_free(dsurfdrho)
    end select  
!    Use this line in case you want to gather the cavity vector pkernel%w%eps in
!    a vector eps(i1,i2,i3) to be ploted (eps to be declared)
!     call PS_gather(kernel=pkernel,src=pkernel%w%eps,dest=eps)
!     i1=n01/2
!     i3=n03/2
!     do i2=1,n02
!      write(10,'(1x,I8,3(1x,1pe26.14e3))') i2,potential(i1,i2,i3),density(i1,i2,i3),&
!            eps(i1,i2,i3)
!     end do
   end if

   einit=0.5_dp*f_dot(density,potential)*pkernel%mesh%volume_element

   rhopot(:,:,:)=density(:,:,:)
!GPEx------------------------------------------------------------------------
! Here we call the electrostatic solver of bigdft: rhopot represents in
! input the charge density rho, in output the potential coming from the solution
! of the vacuum Poisson equation $\Delta \phi(\textbf{r}) = -4 \pi \rho(\textbf{r}) $
! or the Generalized Poisson equation $\nabla \cdot \epsilon( \textbf{r}) \nabla
! \phi(\textbf{r}) = -4 \pi \rho(\textbf{r})$

!   call H_potential('D',pkernel,density(1,1,pkernel%grid%istart+1),density(1,1,pkernel%grid%istart+1),&
!        ehartree,offset,.false.)
   call Electrostatic_Solver(pkernel,rhopot,energies)
   call PS_gather(src=density,kernel=pkernel)

!GPEx------------------------------------------------------------------------
! Here we compare the initial analitical potential (potential) and the output
! potential from the bigdft electrostatic solver (rhopot).
  if (iproc==0) then
     call writeroutinePot(n01,n02,n03,rhopot,potential)
     call yaml_map('Expected hartree energy',einit)
     call yaml_map('Computed Hartree energy',energies%hartree)
     call yaml_map('Diff of expected-computed Hartree energy',einit-energies%hartree)
  end if

 call pkernel_free(pkernel)


  call f_free(density)
  call f_free(rhopot)
  call f_free(rxyz)
  call f_free(radii)
  call f_free(potential)

  call mpifinalize()
  call f_lib_finalize()
  
contains

  !>identify the options from command line
  !! and write the result in options dict
  subroutine PS_Check_command_line_options(options)
    use yaml_parse
    use dictionaries
    implicit none
    !> dictionary of the options of the run
    !! on entry, it contains the options for initializing
    !! on exit, it contains in the key "BigDFT", a list of the 
    !! dictionaries of each of the run that the local instance of BigDFT
    !! code has to execute
    type(dictionary), pointer :: options
    !local variables
    type(yaml_cl_parse) :: parser !< command line parser

    !define command-line options
    parser=yaml_cl_parse_null()
    !between these lines, for another executable using BigDFT as a blackbox,
    !other command line options can be specified
    !then the bigdft options can be specified
    call PS_check_options(parser)
    !parse command line, and retrieve arguments
    call yaml_cl_parse_cmd_line(parser,args=options)
    !free command line parser information
    call yaml_cl_parse_free(parser)

  end subroutine PS_Check_command_line_options


  subroutine set_offset()
   ixc=0
   if (ixc==0) then
      offset=0.0d0
      do i3=1,ndims(3)
         do i2=1,ndims(2)
            do i1=1,ndims(1)
               offset=offset+potential(i1,i2,i3)
            end do
         end do
      end do
      offset=offset*hx*hy*hz*sqrt(mesh%dom%detgd) ! /// to be fixed ///
      if (iproc==0) call yaml_map('offset',offset)
   end if
   pkernel%opt%potential_integral=offset
  end subroutine set_offset

end program PSolver_examples

subroutine writeroutinePot(n01,n02,n03,ri,potential)
  use yaml_output
  use dynamic_memory
  use f_utils
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  real(kind=8), dimension(n01,n02,n03,1), intent(in) :: ri
  real(kind=8), dimension(n01,n02,n03),intent(in) :: potential
  !automatic array, to be check is stack poses problem
  real(kind=8), dimension(:,:,:,:), allocatable :: re
  integer :: i1,i2,i3,i1_max,i2_max,i3_max
  real(kind=8) :: max_val,fact
  re=f_malloc([n01,n02,n03,1],id='re')
  
      max_val = 0.d0
      i1_max = 1
      i2_max = 1
      i3_max = 1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               re(i1,i2,i3,1) = ri(i1,i2,i3,1) - potential(i1,i2,i3)
               fact=abs(re(i1,i2,i3,1))
               if (max_val < fact) then
                  max_val = fact
                  i1_max = i1
                  i2_max = i2
                  i3_max = i3
               end if
            end do
         end do
      end do
      if (max_val == 0.d0) then
         call yaml_map('Inf. Norm difference with reference',0.d0)
      else
         call yaml_mapping_open('Inf. Norm difference with reference')
         call yaml_map('Value',max_val,fmt='(1pe22.15)')
         call yaml_map('Point',[i1_max,i2_max,i3_max],fmt='(i4)')
         call yaml_map('Some values',[re(n01/2,n02/2,n03/2,1),re(2,n02/2,n03/2,1),re(10,n02/2,n03/2,1)],&
              fmt='(1pe22.15)')
         call yaml_mapping_close()
      end if

      call f_free(re)
end subroutine writeroutinePot
