program babel
  use futile
  use at_babel
!  use module_atoms
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
  type(dictionary), pointer :: dict,types,options,iter,it
!  type(atomic_structure) :: astruct

  call f_lib_initialize()

  call yaml_argparse(options,inputs)

  call yaml_comment('Welcome to the F90 openbabel wrapper',hfill='-')

  fin=options//'input'
  fout=options//'output'
  call dict_free(options)
  
  call yaml_mapping_open('Read positions')

  !dict=>dict_new()
  call dict_init(dict)
  call load_dict_from_openbabel(dict,trim(fin))

  call yaml_mapping_close()

  !call dict_to_frags(dict//'positions')

!!$  call astruct_dict_get_types(dict, types)
  call dict_init(types)
  nullify(iter)
  do while (iterating(iter, on = dict // "positions"))
     nullify(it)
     do while (iterating(it, on = iter))
        if (ichar(dict_key(it)) >= ichar('A') .and. &
             & ichar(dict_key(it)) <= ichar('Z') .and. &
             & trim(dict_key(it)) /= 'X') then
           call set(types // dict_key(it), dict_key(it))
        end if
     end do
  end do

  call yaml_map('Read types', types)

  call dump_dict_with_openbabel(dict,types,trim(fout))! fout,len_trim(fout))
  call yaml_map('Positions dumped into file',fout)

  call dict_free(dict, types)

  ! Reload the generated file.
  call dict_init(dict)
  call yaml_mapping_open('Written positions')
  call load_dict_from_openbabel(dict,trim(fout))
  call yaml_mapping_close()
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
         print *,'here',iat,ifrag_max
      end do
      call yaml_map('Fragment dict',frag)
      fileunit=12
      call f_open_file(fileunit,'frag.yaml')
      call dict_init(frag_list)
      call yaml_set_stream(unit=fileunit)
      do id=1,ifrag_max
         atom = frag .get. trim(yaml_toa(id))
         if (.not. associated(atom)) then
            call yaml_warning('The fragment"'+id+'" is not assigned to a fragment')
            cycle
         end if

         iter=>frag//trim(yaml_toa(id))//'atoms'
         call dict_copy(dest=frag_list//(id-1),src=iter)
      end do
      call yaml_dict_dump(frag_list,flow=.true.)
      call yaml_close_stream(unit=fileunit)
      call f_close(fileunit)

      call dict_free(frag_list,frag)
    end subroutine dict_to_frags

end program babel

