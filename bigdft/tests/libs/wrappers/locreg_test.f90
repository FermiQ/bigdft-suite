!> @file
!! Test the usage of the localization regions
!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program locreg_test
  use futile, dict_set => set
  use locregs
  use locreg_operations
  use box
  use f_trees
  use BigDFT_API, only: bigdft_init_errors,bigdft_init_timing_categories
  implicit none
  !logical :: correct
  real(f_double) :: cutoff,crmult,frmult
  type(locreg_descriptors) :: lr,Glr
  type(f_tree) :: dict_posinp
  character(len=max_field_length) :: cellstr
  !integer, dimension(2,3) :: nbox
  real(f_double), dimension(3) :: hgrids,oxyz
  real(f_double), dimension(:), allocatable :: psi
  type(dictionary), pointer :: options
  

  call f_lib_initialize()

  call bigdft_init_errors()
  call bigdft_init_timing_categories()

  call yaml_argparse(options,&
       '- {name: hgrid, shortname: g, default: 0.333, help_string: hgrid}'//f_cr//&
       '- {name: cell, shortname: c, default: "[10,15,11]", help_string: cell dimensions}'//f_cr//&
       '- {name: center, shortname: o, default: 0.0, help_string: center}'//f_cr//&
       '- {name: radius, shortname: r, default: 1.0, help_string: radius}')


  hgrids=options//'hgrid'
  crmult=options//'radius'
  frmult=10.0_f_double
  oxyz=options//'center'
  cutoff=options//'radius'
  cellstr=options//'cell'
  dict_posinp=f_tree_load('{positions: [{ C: '//yaml_toa(oxyz)//&
       '}], cell: '//trim(cellstr)//'}')
  !call yaml_map('Dict posinp',dict_posinp%d)

  !define the global localization region
  call define_lr(Glr,dict_posinp,crmult,frmult,hgrids)

  call dict_free(options)
  call f_tree_free(dict_posinp)

  !define a internal localization region starting from glr
  lr=locreg_null()
  call lr_box(lr,Glr,hgrids,&
       nbox=box_nbox_from_cutoff(Glr%mesh_coarse,oxyz,cutoff))

  !let us now define a dummy wavefunction in the internal locreg
  !define also a method to retrieve the number of components of a locreg
  !to be discussed if this method should be called as lr_nvctr(lr) or
  !by accessing to a ultimate component e.g. lr%wfd%nvctr
  psi=f_malloc0(lr_nvctr(lr),id='psidummy')

  !!!! HERE FOLLOWS API EXAMPLES
!!$
!!$  !>>>> example of the writing API for a locreg (support function)
!!$
!!$  !then we might create a file with the information
!!$  !the idea might be to create a lr_file structure that
!!$  !will contain the information about the header
!!$  lr_file=f_locreg_file(filename='lrfile',lr=lr)
!!$  !then the user will define the type of the file he is going to
!!$  !dump with this structure (by using a public enumerator of the module)
!!$  !example of a binary file
!!$  call lr_file_set_type(lr_file,type=SUPPORT_FUNCTION_TYPE,binary=.true.)
!!$  !then there might be an accessor of this structure that 
!!$  !would add further information in the header of the file
!!$  !the dict_info file is free of format for the sake of 
!!$  !the locreg module but it will be normalized for the orbitals viewpoint
!!$  call lr_file_set_info(lr_file,dict_info)
!!$
!!$  !then the file can be be written
!!$  call lr_file_dump(lr_file,psi)
!!$  !and closed
!!$  call lr_file_close(lr_file)
!!$
!!$  !>>>> example of the loading of a locreg
!!$  lr_file=lr_file_load(filename='lrfile')
!!$
!!$  !with this approach the loaded wavefunction is allocated internally
!!$  !the accessors would be
!!$  call lr_set_from_file(lr,lr_file)
!!$  psi_ptr=>lr_file_get_data_ptr(lr_file) !this will give us the pointer to the data, for storage or reallocation
!!$  dict_info=>lr_file_get_info(lr_file) !gives back the dictionary, in a separate copy
!!$
!!$  call lr_file_close(lr_file) !do what one thinks it does.
  

  

  call deallocate_locreg_descriptors(lr)
  call deallocate_locreg_descriptors(Glr)

  call f_free(psi)
  call f_lib_finalize()

  contains

    !in this section of the mai program we might put tentative methods on
    !the locregs that eventually - after agreement - will be shifted in 
    !the appropriate modules of the locreg directory

    function lr_nvctr(lr)
      !get the number of components associated to a locreg
      !usefule for allocating support functions associated to a locreg
      implicit none
      type(locreg_descriptors), intent (in) :: lr !localisation region
      integer(f_long) :: lr_nvctr !number of components (long integer for Glr)
      lr_nvctr=int(lr%wfd%nvctr_c,f_long)+7*int(lr%wfd%nvctr_f,f_long)
    end function lr_nvctr

    
  
end program locreg_test
