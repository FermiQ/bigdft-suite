!> iniitalize the global localization region of the system given the posinp tree
subroutine define_lr(lr,tree_posinp,crmult,frmult,hgrids)
  use module_atoms
  use locregs
  use locregs_init
  use futile
  use f_trees
  use pseudopotentials
  implicit none
  type(f_tree), intent(in) :: tree_posinp
  type(locreg_descriptors), intent(out) :: lr
  real(f_double), intent(in) :: crmult,frmult
  real(f_double), dimension(3), intent(inout) :: hgrids
  !local variables
  type(atoms_data) :: atoms
  type(dictionary), pointer :: dict_targets,dict_dm,types,var
  real(f_double), dimension(:,:), allocatable :: rxyz

  call dict_init(dict_dm)
  call dict_init(dict_targets)
  call nullify_atoms_data(atoms)

  !> fill the atomic structure datatype
  call astruct_set(atoms%astruct,tree_posinp%d,0.0_f_double,.true.,1.e-8_f_double,&
       [0.0_f_double,0.0_f_double,0.0_f_double],1,.true.)

  ! Generate the dict of types for later use.
  call astruct_dict_get_types(tree_posinp%d, types)

  nullify(var)
  do while(iterating(var,on=types))
     call psp_dict_fill_all(dict_dm, trim(dict_key(var)), 1, 8.0_f_double,crmult,frmult)
  end do
  call dict_free(types)

  call atoms_fill(atoms,dict_dm,frmult,1,.false.,16,0,0.0_f_double)

  rxyz=f_malloc(src=atoms%astruct%rxyz,id='rxyz')

  call lr_set(lr,0,.false.,.true.,crmult,frmult,hgrids,rxyz,atoms,.true.,.false.)

  call f_free(rxyz)

  call deallocate_atoms_data(atoms)
  call dict_free(dict_dm,dict_targets)

end subroutine define_lr
