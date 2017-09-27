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
  real(f_double), dimension(:), allocatable :: psi,tpsi
  type(dictionary), pointer :: options
  
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

  call deallocate_workarrays_projectors(wp)

  call deallocate_locreg_descriptors(lr)

  call f_free(psi,tpsi)
  call f_lib_finalize()


  contains

    subroutine project(psi,method)
      use f_enums
      implicit none
      type(f_enumerator), intent(in) :: method
      real(f_double), dimension(*) :: psi
      call gaussian_to_wavelets_locreg(lr%mesh_coarse,0,&
           1,[1.0_f_double],[0.5_f_double/sigma],UNINITIALIZED(1.0_f_double),&
           1,1, [5.0_f_double,5.0_f_double,5.0_f_double],kpoint,&
           1,lr,wp,psi,method=method)
    end subroutine project
end program projection
