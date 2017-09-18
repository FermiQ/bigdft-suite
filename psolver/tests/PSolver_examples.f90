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
   use psolver_environment, only: rigid_cavity_arrays,rigid_cavity_forces,vacuum_eps,epsprime
   use FDder
   use f_blas, only: f_dot
   implicit none
   
   character(len=4) :: PSol 
   integer :: SetEps
   logical :: usegpu
   real(kind=8), parameter :: acell = 10.d0 
   real(kind=8), parameter :: rad_cav = 2.7d0 
   integer :: nat = 1 ! Number of atoms to build rigid cavity with nat=1.
   character(len=2) :: geocode
   type(cell) :: mesh
   character(len=2), parameter :: datacode = 'G'
   integer, parameter :: nord = 16
   integer, dimension(3) :: ndims
   real(8), dimension(3) :: hgrids,angdeg,angrad
   real(kind=8), parameter :: a_gauss = 1.0d0,a2 = a_gauss**2
   real(kind=8), parameter :: thr = 1.0d-10
   !integer :: m1,m2,m3,md1,md2,md3,nd1,nd2,nd3,n1,n2,n3,
   integer :: itype_scf,n_cell,iproc,nproc,ixc,n01,n02,n03,ifx,ify,ifz,n1,n23
   real(kind=8) :: hx,hy,hz,hgrid,delta
   real(kind=8) :: einit
   real(kind=8) :: ehartree,offset
   real(kind=8), dimension(:,:,:), allocatable :: density,rhopot,rvApp,eps

   logical :: logyes
   logical :: wrtfiles=.true.
   integer :: i_check,unt,igpu
   real(kind=8), dimension(:,:), allocatable :: rxyz
   real(kind=8), dimension(:), allocatable :: radii
   type(coulomb_operator) :: pkernel
   real(kind=8), dimension(:,:,:), allocatable :: potential,rho
   integer :: i1,i2,i3,unit,unit2
   type(dictionary), pointer :: options,dict_input
   real(kind=8) :: alpha,beta,gamma
 
   call f_lib_initialize()

   !read command line
   call PS_Check_command_line_options(options)

   call f_zero(PSol)
   PSol=options .get. 'method'
   if (len_trim(PSol)==0) call f_strcpy(src='VAC',dest=PSol)
   ndims=options // 'ndim'
   geocode=options//'geocode'
   SetEps =options//'seteps'
   usegpu = options // 'accel'
   logyes= options // 'logfile'
   angdeg=options // 'angdeg'
   delta=0.3d0
   delta= options .get. 'deltacav'
   call f_zero(einit)

   call dict_init(dict_input)

   if ('input' .in. options) &
        call dict_copy(dest=dict_input,src=options//'input')

   call dict_free(options)

   igpu=0
   if (usegpu) igpu=1

   n01=ndims(1)
   n02=ndims(2)
   n03=ndims(3)

   hx=acell/real(n01,kind=8)
   hy=acell/real(n02,kind=8)
   hz=acell/real(n03,kind=8)
   hgrids=(/hx,hy,hz/)

   ! Set the angles in radiant. 
   alpha = angdeg(1)/180.0_f_double*pi
   beta  = angdeg(2)/180.0_f_double*pi
   gamma = angdeg(3)/180.0_f_double*pi
   angrad(1) = angdeg(1)/180.0_f_double*pi
   angrad(2) = angdeg(2)/180.0_f_double*pi
   angrad(3) = angdeg(3)/180.0_f_double*pi
  
   mesh=cell_new(geocode,ndims,hgrids,alpha_bc=angrad(1),beta_ac=angrad(2),gamma_ab=angrad(3)) 

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
   rhopot =f_malloc(ndims,id='rhopot')
   rvApp  =f_malloc(ndims,id='rvApp')
   rxyz   =f_malloc([3,nat],id='rxyz')
   radii   =f_malloc([nat],id='radii')
   potential=f_malloc(ndims,id='potential')
   eps=f_malloc(ndims,id='eps')

   n_cell = max(n01,n02,n03)
   hgrid=max(hx,hy,hz)
   ixc=0
   itype_scf=16
   ehartree=0.0

!GPEx------------------------------------------------------------------------
! Here we define our system: number of atoms (nat), their coordinates (rxyz)
    if (nat.eq.1) then
     rxyz(1,1) = hx*real((ndims(1)-1)/2,kind=8)
     rxyz(2,1) = hy*real((ndims(2)-1)/2,kind=8)
     rxyz(3,1) = hz*real((ndims(3)-1)/2,kind=8)
    else if (nat.eq.2) then
     rxyz(1,1) = hx*real(20,kind=8)
     rxyz(2,1) = hy*real(10,kind=8)
     rxyz(3,1) = hz*real(20,kind=8)
     rxyz(1,2) = hx*real(10,kind=8)
     rxyz(2,2) = hy*real(10,kind=8)
     rxyz(3,2) = hz*real(10,kind=8)
    else if (nat.eq.3) then
     rxyz(1:3,1)=[9.300000d0, 9.300337d0, 9.243250d0]!-[2.30d0,2.85d0,3.7d0]
     rxyz(1:3,2)=[9.300000d0, 8.415319d0, 8.700265d0]!-[2.30d0,2.85d0,3.7d0]
     rxyz(1:3,3)=[9.300000d0, 7.299663d0, 9.156750d0]!-[2.30d0,2.85d0,3.7d0]
    end if

!GPEx------------------------------------------------------------------------
! Here we define the input parameters of the dielectric cavity:
!      soft-sphere cavity:van der Waals radii for atoms, delta for the transition region
!      sccs charge-dependent cavity: rhomin, rhomax
    if (nat.eq.1) then
     radii(1)=rad_cav!*1.5d0/0.52917721092d0
    else if (nat.eq.2) then
     radii(1)=rad_cav!*1.5d0/0.52917721092d0
     radii(2)=rad_cav!*1.5d0/0.52917721092d0
    else if (nat.eq.3) then
     radii=[1.4d0,1.0d0,1.0d0]
    end if
    delta=0.3d0


!GPEx------------------------------------------------------------------------
! Here we build up out analytical potential and the corresponding charge density
! rho.

   call test_functions_new2(mesh,acell,a_gauss,pkernel%mu,density,potential)

     !offset, used only for the periodic solver case
     if (ixc==0) then
        offset=0.0d0!potential(1,1,1)!-pot_ion(1,1,1)
        do i3=1,ndims(3)
           do i2=1,ndims(2)
              do i1=1,ndims(1)
                 offset=offset+potential(i1,i2,i3)
              end do
           end do
        end do
        offset=offset*hx*hy*hz*sqrt(mesh%detgd) ! /// to be fixed ///
        !write(*,*) 'offset = ',offset
        if (iproc==0) call yaml_map('offset',offset)
     end if
!------------------------------------------------------------------------

!GPEx------------------------------------------------------------------------
! Here we apply the operator to the potential in order to check
! that input potential and rho are consistents.

   ! Calculate the charge starting from the potential applying the proper Laplace operator.
   call ApplyLaplace(mesh,geocode,n01,n02,n03,hx,hy,hz,potential,rvApp,nord)

  if (iproc==0) then
   call yaml_mapping_open('Comparison between Generalized Poisson operator and analytical density')
   call writeroutinePot(n01,n02,n03,density,0,rvApp)
   call yaml_mapping_close()
  end if

  i_check=1
  rhopot(:,:,:) = density(:,:,:)

!GPEx------------------------------------------------------------------------
! Here we set up the poisson kernel (pkernel) for:
!      1. vacuum calculations; 
!      2. neutral solvent described by the soft-sphere model;
!      3. neutral solvent described by the sccs charge-dependent model.


   if (usegpu) call dict_set(dict_input//'setup'//'accel','CUDA')
   call dict_set(dict_input//'environment'//'delta',delta)
   if (trim(PSol) /= 'VAC') then
      call dict_set(dict_input//'environment'//'cavity','rigid')
      call dict_set(dict_input//'environment'//'gps_algorithm',PSol)
!      if (mPB) call dict_set(dict_input//'environment'//'pb_method','modified')
   end if

      pkernel=pkernel_init(iproc,nproc,dict_input,geocode,ndims,hgrids,&
           alpha_bc=alpha,beta_ac=beta,gamma_ab=gamma)

! allocate cavity vectors if needed for nonvacuum treatments
   if ( trim(PSol)/='VAC') then
    n1=pkernel%ndims(1)
    n23=pkernel%ndims(2)*pkernel%grid%n3p
    call PS_allocate_cavity_workarrays(n1,n23,pkernel%ndims,&
      pkernel%method,pkernel%w)
    call pkernel_set_epsilon(pkernel,nat=nat,rxyz=rxyz,radii=radii)

!    call PS_gather(kernel=pkernel,src=pkernel%w%eps,dest=eps)

    if (wrtfiles) then
     unit=200
     unit2=201
     call f_open_file(unit=unit,file='references_xz_cav.dat')
     call f_open_file(unit=unit2,file='references_line_y_cav.dat')
     i2=n02/2
     do i3=1,n03
        do i1=1,n01
           write(unit,'(2(1x,I8),3(1x,1pe26.14e3))') i1,i3,potential(i1,i2,i3),density(i1,i2,i3),&
             eps(i1,i2,i3)
        end do
     end do
     i1=n01/2
     i3=n03/2
     do i2=1,n02
      write(unit2,'(1x,I8,3(1x,1pe26.14e3))') i2,potential(i1,i2,i3),density(i1,i2,i3),&
             eps(i1,i2,i3)
     end do
     call f_close(unit)
     call f_close(unit2)
     end if
    end if
   einit=0.5_dp*f_dot(rhopot,potential)*pkernel%mesh%volume_element

   call dict_free(dict_input)
   call pkernel_set(pkernel,verbose=.true.)

!GPEx------------------------------------------------------------------------
! Here we call the electrostatic solver of bigdft: rhopot represents in
! input the charge density rho, in output the potential coming from the solution
! of the vacuum Poisson equation $\Delta \phi(\textbf{r}) = -4 \pi \rho(\textbf{r}) $
! or the Generalized Poisson equation $\nabla \cdot \epsilon( \textbf{r}) \nabla
! \phi(\textbf{r}) = -4 \pi \rho(\textbf{r})$

   call H_potential('D',pkernel,rhopot(1,1,pkernel%grid%istart+1),rhopot(1,1,pkernel%grid%istart+1),&
        ehartree,offset,.false.)
   call PS_gather(src=rhopot,kernel=pkernel)

!GPEx------------------------------------------------------------------------
! Here we compare the initial analitical potential (potential) and the output
! potential from the bigdft electrostatic solver (rhopot).
  if (iproc==0) then
     call writeroutinePot(n01,n02,n03,rhopot,pkernel%max_iter,potential)
     call yaml_map('Expected hartree energy',einit)
     call yaml_map('Computed Hartree energy',ehartree)
     call yaml_map('Diff of expected-computed Hartree energy',einit-ehartree)
  end if

 call pkernel_free(pkernel)


  call f_free(density)
  call f_free(rhopot)
  call f_free(rxyz)
  call f_free(radii)
  call f_free(rvApp)
  call f_free(potential)
  call f_free(eps)

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

end program PSolver_examples

subroutine ApplyLaplace(mesh,geocode,n01,n02,n03,hx,hy,hz,x,y,nord)
  use dynamic_memory
  use box
  implicit none
  type(cell), intent(in) :: mesh
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), dimension(n01,n02,n03), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(out) :: y

  ! Local variables.
  real(kind=8), dimension(:,:,:), allocatable :: ddx
  real(kind=8), dimension(:,:,:,:), allocatable :: dx
  real(kind=8) :: pi
  integer :: i1,i2,i3,isp,i

  pi = 4.d0*datan(1.d0)   

  ddx=f_malloc([n01,n02,n03],id='ddx')
  dx=f_malloc([n01,n02,n03,3],id='dx')

   call fssnord3DmatNabla_nonortho(mesh,x,dx,nord)
 
!       isp=1
!       do i3=1,n03
!        do i2=1,n02
!         do i1=1,n01
!          do i=1,3
!           dx(i1,i2,i3,isp,i)=eps(i1,i2,i3)*dx(i1,i2,i3,isp,i)
!          end do
!         end do
!        end do
!       end do
 
   call fssnord3DmatDiv(geocode,n01,n02,n03,hx,hy,hz,dx,y,nord)
 
    y(:,:,:)=-y(:,:,:)/(4.d0*pi)
 
   call f_free(ddx)
   call f_free(dx)

end subroutine ApplyLaplace


subroutine fssnord3DmatNabla_nonortho(mesh,u,du,nord)
      use box
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      type(cell), intent(in) :: mesh
      real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3),1) :: u
      real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3),1,3) :: du
      integer, intent(in) :: nord

!c..local variables
      integer :: n,m,n_cell,n01,n02,n03
      integer :: i,j,i1,i2,i3,isp,ii
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D,c1DF
      real(kind=8) :: hx,hy,hz
      logical :: perx,pery,perz
      real(kind=8), dimension(3) :: der
      !>parameter for the definition of the bc
      integer, parameter :: FREE=0
      integer, parameter :: PERIODIC=1

      n01=mesh%ndims(1)
      n02=mesh%ndims(2)
      n03=mesh%ndims(3)
      hx=mesh%hgrids(1)
      hy=mesh%hgrids(2)
      hz=mesh%hgrids(3)
      n = nord+1
      m = nord/2
      n_cell = maxval(mesh%ndims)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      perx=(mesh%bc(1) /= FREE)
      pery=(mesh%bc(2) == PERIODIC)
      perz=(mesh%bc(3) /= FREE)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
        c1DF(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      isp=1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3,isp,1) = 0.0d0

             if (i1.le.m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,0)*u(ii,i2,i3,isp)!/hx
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,i1-m-1)*u(j+m+1,i2,i3,isp)!/hx
               end do
              end if
             else if (i1.gt.n01-m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j - 1, n01 ) + 1
                du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,0)*u(ii,i2,i3,isp)!/hx
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3,isp)!/hx
               end do
              end if
             else
              do j=-m,m
               du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,0)*u(i1 + j,i2,i3,isp)!/hx
              end do
             end if
             du(i1,i2,i3,isp,1)=du(i1,i2,i3,isp,1)/hx

             du(i1,i2,i3,isp,2) = 0.0d0

             if (i2.le.m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,0)*u(i1,ii,i3,isp)!/hy
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,i2-m-1)*u(i1,j+m+1,i3,isp)!/hy
               end do 
              end if
             else if (i2.gt.n02-m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j - 1, n02 ) + 1
                du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,0)*u(i1,ii,i3,isp)!/hy
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3,isp)!/hy
               end do
              end if
             else
              do j=-m,m
               du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,0)*u(i1,i2 + j,i3,isp)!/hy
              end do
             end if
              du(i1,i2,i3,isp,2)=du(i1,i2,i3,isp,2)/hy

             du(i1,i2,i3,isp,3) = 0.0d0

             if (i3.le.m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,0)*u(i1,i2,ii,isp)!/hx
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,i3-m-1)*u(i1,i2,j+m+1,isp)!/hz
               end do
              end if
             else if (i3.gt.n03-m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j - 1, n03 ) + 1
                du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,0)*u(i1,i2,ii,isp)!/hx
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m,isp)!/hz
               end do
              end if
             else
              do j=-m,m
               du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,0)*u(i1,i2,i3 + j,isp)!/hz
              end do
             end if
              du(i1,i2,i3,isp,3)=du(i1,i2,i3,isp,3)/hz

              der(1:3)=0.d0
              do i=1,3
               do j=1,3
                der(i) =  der(i) + mesh%gu(i,j)*du(i1,i2,i3,isp,j)
               end do
              end do
              du(i1,i2,i3,isp,1:3) = der(1:3)

            end do
         end do
      end do

end subroutine fssnord3DmatNabla_nonortho


subroutine fssnord3DmatDiv(geocode,n01,n02,n03,hx,hy,hz,u,du,nord)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      character(len=2), intent(in) :: geocode
      integer, intent(in) :: n01,n02,n03,nord
      real(kind=8), intent(in) :: hx,hy,hz
      real(kind=8), dimension(n01,n02,n03,1,3) :: u
      real(kind=8), dimension(n01,n02,n03,1) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,i1,i2,i3,isp,ii
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D
      real(kind=8) :: d1,d2,d3
      real(kind=8), parameter :: zero = 0.d0! 1.0d-11
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      n_cell = max(n01,n02,n03)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      isp=1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3,isp) = 0.0d0

             d1 = 0.d0
             if (i1.le.m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                d1 = d1 + c1D(j,0)*u(ii,i2,i3,isp,1)!/hx
               end do
              else
               do j=-m,m
                d1 = d1 + c1D(j,i1-m-1)*u(j+m+1,i2,i3,isp,1)!/hx
               end do
              end if
             else if (i1.gt.n01-m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j - 1, n01 ) + 1
                d1 = d1 + c1D(j,0)*u(ii,i2,i3,isp,1)!/hx
               end do
              else
               do j=-m,m
                d1 = d1 + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3,isp,1)!/hx
               end do
              end if
             else
              do j=-m,m
               d1 = d1 + c1D(j,0)*u(i1 + j,i2,i3,isp,1)!/hx
              end do
             end if
              d1=d1/hx

             d2 = 0.d0
             if (i2.le.m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                d2 = d2 + c1D(j,0)*u(i1,ii,i3,isp,2)!/hy
               end do
              else
               do j=-m,m
                d2 = d2 + c1D(j,i2-m-1)*u(i1,j+m+1,i3,isp,2)!/hy
               end do 
              end if
             else if (i2.gt.n02-m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j - 1, n02 ) + 1
                d2 = d2 + c1D(j,0)*u(i1,ii,i3,isp,2)!/hy
               end do
              else
               do j=-m,m
                d2 = d2 + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3,isp,2)!/hy
               end do
              end if
             else
              do j=-m,m
               d2 = d2 + c1D(j,0)*u(i1,i2 + j,i3,isp,2)!/hy
              end do
             end if
              d2=d2/hy

             d3 = 0.d0
             if (i3.le.m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                d3 = d3 + c1D(j,0)*u(i1,i2,ii,isp,3)!/hz
               end do
              else
               do j=-m,m
                d3 = d3 + c1D(j,i3-m-1)*u(i1,i2,j+m+1,isp,3)!/hz
               end do
              end if
             else if (i3.gt.n03-m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j - 1, n03 ) + 1
                d3 = d3 + c1D(j,0)*u(i1,i2,ii,isp,3)!/hz
               end do
              else
               do j=-m,m
                d3 = d3 + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m,isp,3)!/hz
               end do
              end if
             else
              do j=-m,m
               d3 = d3 + c1D(j,0)*u(i1,i2,i3 + j,isp,3)!/hz
              end do
             end if
              d3=d3/hz

             du(i1,i2,i3,isp) = d1+d2+d3

            end do
         end do
      end do

end subroutine fssnord3DmatDiv

subroutine writeroutinePot(n01,n02,n03,ri,i,potential)
  use yaml_output
  use dynamic_memory
  use f_utils
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: i
  real(kind=8), dimension(n01,n02,n03,1), intent(in) :: ri
  real(kind=8), dimension(n01,n02,n03),intent(in) :: potential
  !automatic array, to be check is stack poses problem
  real(kind=8), dimension(:,:,:,:), allocatable :: re
  integer :: i1,i2,i3,i1_max,i2_max,i3_max,unt
  real(kind=8) :: max_val,fact
  logical :: wrtfiles=.false.
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
      if (wrtfiles) then      
      write(38,'(4(1x,I4),2(1x,e22.15))')i,i1_max,i2_max,i3_max,max_val,&
           re(n01/2,n02/2,n03/2,1)
!!$      write(38,'(4(1x,I4),4(1x,e22.15))')i,i1_max,i2_max,i3_max,max_val,&
!!$           re(n01/2,n02/2,n03/2,1),re(2,n02/2,n03/2,1),re(10,n02/2,n03/2,1)
      !write(*,'(4(1x,I4),4(1x,e22.15))')i,i1_max,i2_max,i3_max,max_val,&
      !     re(n01/2,n02/2,n03/2,1),re(2,n02/2,n03/2,1),re(10,n02/2,n03/2,1)
      end if
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

      if (wrtfiles) then      
       unt=f_get_free_unit(21)
       call f_open_file(unt,file='final.dat')
       !i1=n01/2
       i1=(n01-1)/2+1
       do i2=1,n02
          do i3=1,n03
             write(unt,'(2(1x,I4),2(1x,e14.7))')i2,i3,ri(i1,i2,i3,1),potential(i1,i2,i3)
          end do
       end do
       call f_close(unt)
 
       unt=f_get_free_unit(22)
       i1=(n01-1)/2+1
       i3=(n03-1)/2+1
       call f_open_file(unt,file='final_line.dat')
       !write (str, *) i
       !str = adjustl(str)
       !call f_open_file(unt,file='final_line'//trim(str)//'.dat')
       do i2=1,n02
        write(unt,'(1x,I8,3(1x,e22.15))') i2,ri(i1,i2,i3,1),potential(i1,i2,i3),re(i1,i2,i3,1)
       end do
       call f_close(unt)
      end if

      call f_free(re)
end subroutine writeroutinePot
