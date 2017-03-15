program babel
  use futile
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
  type(dictionary), pointer :: dict,options

  interface
     subroutine openbabel_load(d, f,ln)
       use dictionaries
       implicit none
       type(dictionary), pointer :: d
       character(len = *), intent(in) :: f
       integer, intent(in) :: ln
     end subroutine openbabel_load
     subroutine openbabel_dump(d, f,ln)
       use dictionaries
       implicit none
       type(dictionary), pointer :: d
       character(len = *), intent(in) :: f
       integer, intent(in) :: ln
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
  call openbabel_load(dict,fin,len_trim(fin))

  call yaml_map(fin,dict)
  call yaml_mapping_close()

  call openbabel_dump(dict,fout,len_trim(fout))
  call yaml_map('Positions dumped into file',fout)
  
  call dict_free(options,dict)
  call f_lib_finalize()
end program babel
