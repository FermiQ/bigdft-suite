!> test program for the function projection in wavelets
program projection
  use module_defs, only: UNINITIALIZED
  use futile
  use box
  use f_functions
  use locregs
  use locreg_operations
  use f_trees
  use BigDFT_API, only: bigdft_init_errors,bigdft_init_timing_categories
  use numerics
  use compression, only: wnrm2
  implicit none
  real(f_double) :: crmult,frmult,maxdiff,sigma
  type(locreg_descriptors) :: lr
  real(f_double), dimension(3) :: kpoint,oxyz,angrad,hgrids
  type(f_tree) :: dict_posinp
  type(workarrays_projectors) :: wp
  type(workarr_sumrho) :: w
  real(f_double), dimension(:), allocatable :: psi,tpsi
  type(dictionary), pointer :: options
  real(f_double), dimension(:), allocatable :: projector_real,gaussian
  real(f_double), dimension(3) :: rxyz
  integer, parameter :: n=1 !<principal quantum number
  integer, parameter :: l=1 !<angular momentum of the shell
  integer, parameter :: ider=0 !<direction in which to perform the derivative (0 if any)
  integer, parameter :: nterm_max=20
  integer, parameter :: ncplx_g=1
  real(f_double), dimension(ncplx_g) :: expo 
  real(f_double), dimension(ncplx_g) :: coeff !<prefactor of the gaussian

  call f_lib_initialize()
 
  call bigdft_init_errors()
  call bigdft_init_timing_categories()


  call yaml_argparse(options,&
       '- {name: hgrid, shortname: g, default: 0.333, help_string: hgrid}'//f_cr//&
       '- {name: sigma, shortname: s, default: 0.3  , help_string: sigma}')

  !hgrids=0.5_f_double
  hgrids=options//'hgrid'
  sigma=options//'sigma'
  dict_posinp=f_tree_load('{positions: [{ C: [0.0, 0.0, 0.0]}], cell: [10,15,11]}')
  crmult=10.0_f_double
  frmult=10.0_f_double
  angrad=onehalf*pi
  oxyz=5.0_f_double
  kpoint=0.0_f_double

  coeff=[1.0_f_double]
  expo=[0.5_f_double/sigma]
  rxyz=[5.0_f_double,5.0_f_double,5.0_f_double]

  call dict_free(options)
  
  call define_lr(lr,dict_posinp,crmult,frmult,hgrids)

  call f_tree_free(dict_posinp)
  call allocate_workarrays_projectors(lr%d%n1, lr%d%n2, lr%d%n3, wp)
  psi=f_malloc0(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,id='psi')
  tpsi=f_malloc0(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,id='tpsi')

  call project(psi,PROJECTION_1D_SEPARABLE)

  call yaml_mapping_open('Traditional separable Projection')
  call yaml_map('Norm of the calculated projector',wnrm2(1,lr%wfd,psi))
  call yaml_mapping_close()

  call project(tpsi,PROJECTION_RS_COLLOCATION)

  call yaml_mapping_open('Collocation-based separable Projection')
  call yaml_map('Norm of the calculated projector',wnrm2(1,lr%wfd,tpsi))
  call yaml_mapping_close()

  !calculate the difference of the two arrays
  call f_diff(f_size(psi),psi,tpsi,maxdiff)
  call yaml_map('Maximum difference of the two arrays',maxdiff)

  ! compare input analytical gaussian and the roundtrip one to the daubechies
  projector_real=f_malloc(lr%mesh%ndim,id='projector_real') 
  call initialize_work_arrays_sumrho(lr,.true.,w)
  call daub_to_isf(lr,w,tpsi,projector_real)

  !build up the input gaussian as done in gaussian_to_wavelets_locreg
  gaussian=f_malloc(lr%mesh%ndim,id='gaussian')
  oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
  call real_space_gaussian(ncplx_g,n,l,ider,nterm_max,coeff,expo,rxyz,oxyz,lr,gaussian)

  !calculate the difference of the two arrays
  call f_diff(f_size(gaussian),gaussian,projector_real,maxdiff)
  call yaml_map('Maximum difference of the in/out gaussian',maxdiff)

  call deallocate_workarrays_projectors(wp)
  call deallocate_work_arrays_sumrho(w)
  call deallocate_locreg_descriptors(lr)
  call f_free(projector_real)
  call f_free(gaussian)

  call f_free(psi,tpsi)
  call f_lib_finalize()

  contains

    subroutine project(psi,method)
      use f_enums
      implicit none
      type(f_enumerator), intent(in) :: method
      real(f_double), dimension(*) :: psi
      call gaussian_to_wavelets_locreg(lr%mesh_coarse,0,&
           1,coeff,expo,UNINITIALIZED(1.0_f_double),&
           1,1,rxyz,kpoint,&
           1,lr,wp,psi,method=method)
    end subroutine project

    subroutine real_space_gaussian(ncplx_g,n,l,ider,nterm_max,coeff,expo,rxyz,oxyz,lr,gaussian)
      implicit none
      integer, intent(in) :: ncplx_g !< 1 or 2 if the gaussian factor is real or complex respectively
      integer, intent(in) :: n !<principal quantum number
      integer, intent(in) :: l !<angular momentum of the shell
      integer, intent(in) :: ider !<direction in which to perform the derivative (0 if any)
      integer, intent(in) :: nterm_max !if GTH nterm_max=4 (this value should go in a module)
      real(f_double), dimension(ncplx_g), intent(in) :: coeff !<prefactor of the gaussian
      real(f_double), dimension(ncplx_g), intent(in) :: expo 
      real(f_double), dimension(3), intent(in) :: rxyz,oxyz 
      type(locreg_descriptors), intent(in) :: lr
      real(f_double), dimension(lr%mesh%ndim), intent(out) :: gaussian
      ! Local variables
      integer, dimension(2*l-1) :: nterms
      integer, dimension(nterm_max,3,2*l-1) :: lxyz
      real(f_double), dimension(ncplx_g) :: sigma_and_expo
      real(f_double), dimension(ncplx_g,nterm_max,2*l-1) :: factors
      integer :: m,i
      type(box_iterator) :: bit
      type(f_function), dimension(3) :: funcs
      real(f_double), dimension(3) :: noxyz 

      call get_projector_coeffs(ncplx_g,l,n,ider,nterm_max,coeff,expo,&
           nterms,lxyz,sigma_and_expo,factors)
      !for the moment only with s projectors (l=0,n=1)
      noxyz=rxyz-oxyz
      bit=box_iter(lr%mesh,origin=noxyz) !use here the real space mesh of the projector locreg
      do m=1,2*l-1
         do i=1,3
            funcs(i)=f_function_new(f_gaussian,exponent=expo(1))
         end do
         !here we do not consider the lxyz terms yet
         !take the reference functions
         print *,size(gaussian),'real',lr%mesh%ndims,&
              lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3],&
              lr%mesh_coarse%hgrids*[lr%ns1,lr%ns2,lr%ns3],rxyz,noxyz
         call separable_3d_function(bit,funcs,factors(1,1,m)*sqrt(lr%mesh%volume_element),gaussian)
      end do !not correctly written, it should be used to define the functions
    end subroutine real_space_gaussian

end program projection
