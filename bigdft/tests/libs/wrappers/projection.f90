!> test program for the function projection in wavelets
program projection
  use module_defs, only: UNINITIALIZED
  use futile
  use box
  use f_functions
  use locregs
  use gaussians
  use locreg_operations
  use f_trees
  use BigDFT_API, only: bigdft_init_errors,bigdft_init_timing_categories
  use numerics
  use compression, only: wnrm2
  use multipole_preserving
  implicit none
  real(f_double) :: crmult,frmult,maxdiff_RS,maxdiff_MP,maxdiff_RS_MP,sigma
  type(locreg_descriptors) :: lr
  real(f_double), dimension(3) :: kpoint,oxyz,angrad,hgrids
  type(f_tree) :: dict_posinp
  type(workarrays_projectors) :: wp
  real(f_double), dimension(:,:,:), allocatable :: psi,RSpsi,MPpsi
  type(dictionary), pointer :: options
  real(f_double), dimension(3) :: rxyz
  integer :: n !<principal quantum number
  integer :: l !<angular momentum of the shell
  integer :: m !<magnetic momentum of the shell
  integer, parameter :: nmax=1 !<maximum principal quantum number
  integer, parameter :: lmax=4 !<maximum angular momentum of the shell
  integer, parameter :: ider=0 !<direction in which to perform the derivative (0 if any)
  integer, parameter :: nterm_max=20
  integer, parameter :: ncplx_g=1 !< 1 or 2 if the gaussian factor is real or complex respectively
  integer, parameter :: ncplx_p=1 !< 2 if the projector is supposed to be complex, 1 otherwise
  real(f_double), dimension(ncplx_g) :: expo !<exponents (1/2sigma^2 for the first element) of the gaussian (real and imaginary part)
  real(f_double), dimension(ncplx_g) :: coeff !<prefactor of the gaussian
  integer :: iproc,nn,qn,ql,nstart,nend,lstart,lend
  !type(workarr_sumrho) :: w
  !real(f_double), dimension(:), allocatable :: projector_real,gaussian
  !real(f_double), dimension(:), allocatable :: ps,RSps,MPps 
  !integer :: ni,nr

  ! To run the test:
  ! ./projection -> test all projectors

  qn=0
  ql=0
  iproc=0
  call f_lib_initialize()
 
  call bigdft_init_errors()
  call bigdft_init_timing_categories()

  call yaml_argparse(options,&
       '- {name: hgrid, shortname: g, default: 0.333, help_string: hgrid}'//f_cr//&
       '- {name: qn, shortname: n, default: 0, help_string: n}'//f_cr//&
       '- {name: ql, shortname: l, default: 0, help_string: l}'//f_cr//&
       '- {name: sigma, shortname: s, default: 0.3  , help_string: sigma}')

  !hgrids=0.5_f_double
  hgrids=options//'hgrid'
  sigma=options//'sigma'
  qn=options//'qn'
  ql=options//'ql'

  dict_posinp=f_tree_load('{positions: [{ C: [0.0, 0.0, 0.0]}], cell: [10,15,11]}')
  crmult=10.0_f_double
  frmult=10.0_f_double
  angrad=onehalf*pi
  oxyz=5.0_f_double
  kpoint=0.0_f_double

  coeff=[1.0_f_double]
  expo=[0.5_f_double/sigma**2]
  rxyz=[5.0_f_double,5.0_f_double,5.0_f_double]

  call dict_free(options)

  
  call define_lr(lr,dict_posinp,crmult,frmult,hgrids)

  call f_tree_free(dict_posinp)

  if (qn == 0 .or. ql == 0) then
     nstart=1
     nend=nmax
     lstart=1
     lend=lmax
  else
     nstart=qn
     nend=qn
     lstart=ql
     lend=ql
  end if

  do n=nstart,nend
     do l=lstart,lend
   
        if (iproc==0) then
          call yaml_comment('Projectors check',hfill='-')
          call yaml_mapping_open('Projector type')
          call yaml_map('principal quantum number n', n)
          call yaml_map('angular momentum of the shell l', l)
        end if
        call allocate_workarrays_projectors(lr%d%n1, lr%d%n2, lr%d%n3, wp)
        nn=lr%wfd%nvctr_c+7*lr%wfd%nvctr_f 
        psi=f_malloc0([nn,ncplx_p,2*l-1],id='psi')
        RSpsi=f_malloc0([nn,ncplx_p,2*l-1],id='RSpsi')
        MPpsi=f_malloc0([nn,ncplx_p,2*l-1],id='MPpsi')
      
        call project(psi,PROJECTION_1D_SEPARABLE)

        call initialize_real_space_conversion(isf_m=16)

        ! Use of the Collocation-based separable Projection
        call project(RSpsi,PROJECTION_RS_COLLOCATION)

        ! Use of the Multipole-preserving-based separable Projection
        call project(MPpsi,PROJECTION_MP_COLLOCATION)

        do m=1,2*l-1
      
           !calculate norm and difference of the two arrays
           call f_diff(f_size(psi(1:nn,ncplx_p,m)),psi(1:nn,ncplx_p,m),RSpsi(1:nn,ncplx_p,m),maxdiff_RS)
           call f_diff(f_size(psi(1:nn,ncplx_p,m)),psi(1:nn,ncplx_p,m),MPpsi(1:nn,ncplx_p,m),maxdiff_MP)
           call f_diff(f_size(RSpsi(1:nn,ncplx_p,m)),RSpsi(1:nn,ncplx_p,m),MPpsi(1:nn,ncplx_p,m),maxdiff_RS_MP)

           if (iproc==0) then
              call yaml_mapping_open('Shell projectors')
              call yaml_map('magnetic numer m', m)
              call yaml_mapping_open('Checks')

              call yaml_mapping_open('Traditional separable Projection')
              call yaml_map('Norm of the calculated projector',wnrm2(1,lr%wfd,psi(1:nn,ncplx_p,m)))
              call yaml_mapping_close()
      
              call yaml_mapping_open('Collocation-based separable Projection')
              call yaml_map('Norm of the calculated projector',wnrm2(1,lr%wfd,RSpsi(1:nn,ncplx_p,m)))
              call yaml_map('Maximum difference with Traditional',maxdiff_RS)
              call yaml_mapping_close()
      
              call yaml_mapping_open('Multipole-preserving-based separable Projection')
              call yaml_map('Norm of the calculated projector',wnrm2(1,lr%wfd,MPpsi(1:nn,ncplx_p,m)))
              call yaml_map('Maximum difference with Traditional',maxdiff_MP)
              call yaml_mapping_close()

              call yaml_map('Maximum difference between C-b and M-p',maxdiff_RS_MP)

              call yaml_mapping_close()
              call yaml_mapping_close()
           end if

        end do ! loop on m quantum number

      !---------------------------------------------------------------------------------------
      !  ! printing of the projectors for debug
      !  call initialize_work_arrays_sumrho(lr,.true.,w)
      !  nr=lr%mesh%ndim
      !  ps=f_malloc(nr,id='ps') 
      !  call daub_to_isf(lr,w,psi,ps)

      !  RSps=f_malloc(nr,id='RSps') 
      !  call daub_to_isf(lr,w,RSpsi,RSps)

      !  MPps=f_malloc(nr,id='MPps') 
      !  call daub_to_isf(lr,w,MPpsi,MPps)

      !  ni=nr**(1.0_f_double/3.0_f_double)
      !  call print_vect(ni,ni,ni,16,1,ps)
      !  call print_vect(ni,ni,ni,17,1,RSps)
      !  call print_vect(ni,ni,ni,18,1,MPps)

      !  call f_free(ps)
      !  call f_free(RSps)
      !  call f_free(MPps)
      !  call deallocate_work_arrays_sumrho(w)
      !---------------------------------------------------------------------------------------

        call f_free(psi)
        call f_free(RSpsi)
        call f_free(MPpsi)
      
      !---------------------------------------------------------------------------------------
      !  ! compare input analytical gaussian and the roundtrip one to the daubechies
      !  projector_real=f_malloc(lr%mesh%ndim,id='projector_real') 
      !  call initialize_work_arrays_sumrho(lr,.true.,w)
      !  call daub_to_isf(lr,w,tpsi,projector_real)
      !
      !  !build up the input gaussian as done in gaussian_to_wavelets_locreg
      !  gaussian=f_malloc(lr%mesh%ndim,id='gaussian')
      !  oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
      !  call real_space_gaussian(ncplx_g,n,l,ider,nterm_max,coeff,expo,rxyz,lr,gaussian)
      !
      !  !calculate the difference of the two arrays
      !  call f_diff(f_size(gaussian),gaussian,projector_real,maxdiff)
      !
      !  call yaml_map('Maximum difference of the in/out gaussian',maxdiff)
      
      !  call deallocate_work_arrays_sumrho(w)
      !  call f_free(projector_real)
      !  call f_free(gaussian)
      !---------------------------------------------------------------------------------------
      
        call finalize_real_space_conversion()
        call deallocate_workarrays_projectors(wp)
      
        if (iproc==0) call yaml_mapping_close()
   
     end do ! loop on l quantum number
  end do ! loop on n quantum number

  call deallocate_locreg_descriptors(lr)

  call f_lib_finalize()

  contains

    subroutine project(psi,method)
      use f_enums
      implicit none
      type(f_enumerator), intent(in) :: method
      real(f_double), dimension(*) :: psi
      call gaussian_to_wavelets_locreg(lr%mesh_coarse,0,&
           ncplx_g,coeff,expo,UNINITIALIZED(1.0_f_double),&
           n,l,rxyz,kpoint,&
           ncplx_p,lr,wp,psi,method=method)
    end subroutine project    

end program projection

    subroutine real_space_gaussian(ncplx_g,n,l,ider,nterm_max,coeff,expo,rxyz,lr,gaussian)
      use gaussians
      use futile
      use locregs
      use box
      use f_functions
      implicit none
      integer, intent(in) :: ncplx_g !< 1 or 2 if the gaussian factor is real or complex respectively
      integer, intent(in) :: n !<principal quantum number
      integer, intent(in) :: l !<angular momentum of the shell
      integer, intent(in) :: ider !<direction in which to perform the derivative (0 if any)
      integer, intent(in) :: nterm_max !if GTH nterm_max=4 (this value should go in a module)
      real(f_double), dimension(ncplx_g), intent(in) :: coeff !<prefactor of the gaussian
      real(f_double), dimension(ncplx_g), intent(in) :: expo 
      real(f_double), dimension(3), intent(in) :: rxyz!,oxyz 
      type(locreg_descriptors), intent(in) :: lr
      real(f_double), dimension(lr%mesh%ndim), intent(out) :: gaussian
      ! Local variables
      integer, dimension(2*l-1) :: nterms
      integer, dimension(nterm_max,3,2*l-1) :: lxyz
      real(f_double), dimension(ncplx_g) :: sigma_and_expo
      real(f_double), dimension(ncplx_g,nterm_max,2*l-1) :: factors
      !integer :: m,i
      type(box_iterator) :: bit
      !type(f_function), dimension(3) :: funcs
      type(gaussian_real_space) :: g
      !real(f_double), dimension(3) :: noxyz 

      call get_projector_coeffs(ncplx_g,l,n,ider,nterm_max,coeff,expo,&
           nterms,lxyz,sigma_and_expo,factors)

      !call gaussian_real_space_set(g,sqrt(onehalf/expo(1)),1,factors,lxyz)
      call gaussian_real_space_set(g,sigma_and_expo(1),1,factors,lxyz(1,:,1),[0],16)
      call f_zero(gaussian)
      bit=box_iter(lr%mesh)
      call three_dimensional_density(bit,g,sqrt(lr%mesh%volume_element),rxyz,gaussian)

!!$      !for the moment only with s projectors (l=0,n=1)
!!$      noxyz=rxyz-oxyz
!!$
!!$      bit=box_iter(lr%mesh,origin=noxyz) !use here the real space mesh of the projector locreg
!!$
!!$      
!!$
!!$      do m=1,2*l-1
!!$         do i=1,3
!!$            funcs(i)=f_function_new(f_gaussian,exponent=expo(1))
!!$         end do
!!$
!!$         !here we do not consider the lxyz terms yet
!!$         !take the reference functions
!!$         !print *,size(gaussian),'real',lr%mesh%ndims,&
!!$         !     lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3],&
!!$         !     lr%mesh_coarse%hgrids*[lr%ns1,lr%ns2,lr%ns3],rxyz,noxyz
!!$         call separable_3d_function(bit,funcs,factors(1,1,m)*sqrt(lr%mesh%volume_element),gaussian)
!!$      end do !not correctly written, it should be used to define the functions
    end subroutine real_space_gaussian

subroutine print_vect(n01,n02,n03,uni,dir,psi)
  use futile
  implicit none
  integer, intent(in) :: n01,n02,n03,uni,dir
  real(f_double), dimension(n01,n02,n03), intent(in) :: psi
  integer :: i1,i2,i3

  select case(dir)
  case(1) 
     i3=n03/2
     i2=n02/2
     do i1=1,n01
     write(uni,'(1(1x,I4),1x,e14.7)')i1,psi(i1,i2,i3)
     end do
  case(2) 
     i3=n03/2
     i1=n01/2
     do i2=1,n02
     write(uni,'(1(1x,I4),1x,e14.7)')i2,psi(i1,i2,i3)
     end do
  case(3) 
     i1=n01/2
     i2=n02/2
     do i3=1,n03
     write(uni,'(1(1x,I4),1x,e14.7)')i3,psi(i1,i2,i3)
     end do
  end select 

end subroutine print_vect
