program babel
  use futile
  use module_atoms
  character(len=*), parameter :: input1=&
       "  {name: input, shortname: i, default: None,"//&
       "  help_string: Input file,"//&
       "  help_dict: {Allowed values: filepath}}"
  character(len=*), parameter :: input2=&
       "  {name: output,shortname: o, default: outfile.xyz,"//&
       "  help_string: Output file,"//&
       "  help_dict: {Allowed values: filepath}}"

  character(len=*), parameter :: inputs=&
       '-'//input1//f_cr//&
       '-'//input2

  character(len=64) :: fin, fout
  type(dictionary), pointer :: dict,types,options,iter
  type(atomic_structure) :: astruct

  interface
     subroutine openbabel_load(d, f)
       use dictionaries
       implicit none
       type(dictionary), pointer :: d
       character, dimension(*), intent(in) :: f
       !character(len = *), intent(in) :: f
       !integer, intent(in) :: ln
     end subroutine openbabel_load
     subroutine openbabel_dump(d, t, f)
       use dictionaries
       implicit none
       type(dictionary), pointer :: d, t
       character, dimension(*), intent(in) :: f
       !character(len = *), intent(in) :: f
       !integer, intent(in) :: ln
     end subroutine openbabel_dump
  end interface

  call f_lib_initialize()

  call yaml_argparse(options,inputs)

  call yaml_comment('Welcome to the F90 openbabel wrapper',hfill='-')

  fin=options//'input'
  fout=options//'output'
  
  call yaml_mapping_open('Reading positions')

  !dict=>dict_new()
  call dict_init(dict)
  call openbabel_load(dict,f_char_ptr(trim(fin)))!,len_trim(fin))

  call yaml_map(fin,dict)
  call yaml_mapping_close()

  call dict_to_frags(dict//'positions')

  call astruct_dict_get_types(dict, types)
  nullify(iter)
  do while (iterating(iter, on = types))
     call set(iter, dict_key(iter))
  end do
  
  call openbabel_dump(dict,types,f_char_ptr(trim(fout)))! fout,len_trim(fout))
  call yaml_map('Positions dumped into file',fout)

  call dict_free(options,dict, types)

  ! Reload the generated file.
  astruct = atomic_structure_null()
  call set_astruct_from_file("outfile.xyz", 0, astruct)

  call dict_init(dict)
  call astruct_merge_to_dict(dict, astruct, astruct%rxyz)
  call deallocate_atomic_structure(astruct)

  call yaml_map("outfile.xyz", dict)

  call dict_free(dict)

  call f_lib_finalize()

  contains

    subroutine dict_to_frags(dict)
      !convert dictionaries to fragments to be passed to the chess-toolbox routine
      implicit none
      type(dictionary), pointer :: dict
      !local variables
      integer :: iat,id,ifrag_max,fileunit
      character(len=32) :: fragname
      type(dictionary), pointer :: frag,iter,atom,frag_list

      frag=>dict_new()
      iter => null()
      iat=0
      ifrag_max=0
      !create the dict of the atoms which belong to each fragment
      do while(iterating(iter,on=dict))
         call f_increment(iat)
         atom = iter .get. 'frag'
         if (.not. associated(atom)) then
            call yaml_warning('The atom "'+iat+'" is not assigned to a fragment')
            cycle
         end if
         id=atom//1
         ifrag_max=max(ifrag_max,id)
         fragname=atom//0
         !now store the data in the dictionary
         call set(frag//trim(yaml_toa(id))//'name',fragname)
         call add(frag//trim(yaml_toa(id))//'atoms',iat)
      end do
      fileunit=12
      call f_open_file(fileunit,'frag.yaml')
      call dict_init(frag_list)
      call yaml_set_stream(unit=fileunit)
      do id=1,ifrag_max
         iter=>frag//trim(yaml_toa(id))//'atoms'
         call dict_copy(dest=frag_list//id-1,src=iter)
      end do
      call yaml_dict_dump(frag_list,flow=.true.)
      call yaml_close_stream(unit=fileunit)
      call f_close(fileunit)

      call dict_free(frag_list,frag)
    end subroutine dict_to_frags

end program babel

