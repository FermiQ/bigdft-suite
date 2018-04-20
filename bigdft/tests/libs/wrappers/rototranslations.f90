!> test program for the function projection in wavelets
program all_you_can_rotate_and_translate
  use module_defs, only: UNINITIALIZED
  use futile
  use box
  use f_functions
  use locregs
  use gaussians
  use locreg_operations
  use f_trees
  use BigDFT_API, only: bigdft_init_errors,bigdft_init_timing_categories
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

  call f_lib_initialize()

  call bigdft_init_errors()
  call bigdft_init_timing_categories()

  call yaml_argparse(options,&
       '- {name: hgrid, shortname: g, default: 0.333, help_string: hgrid}'//f_cr//&
       '- {name: sigma, shortname: s, default: 0.3  , help_string: sigma}')

  !hgrids=0.5_f_double
  hgrids=options//'hgrid'
  sigma=options//'sigma'
  dict_src=f_tree_load('{positions: [{ C: [0.0, 0.0, 0.0]}], cell: [10,15,11]}')
  
  dict_dest=f_tree_load('{positions: [{ C: [0.0, 0.0, 0.0]}], cell: [10,15,11]}')
  
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

  call f_lib_finalize()
end program all_you_can_rotate_and_translate
