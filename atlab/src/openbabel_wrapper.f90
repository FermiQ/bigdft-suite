!> @file
!!    Modulefile for handling the conversion from a openbabel supported file
!!
!! @author
!!    D. Caliste, L. Genovese (September 2015)
!!    Copyright (C) 2015-2019 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU Lesser General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module at_babel

  implicit none


contains

  subroutine load_dict_from_openbabel(dict,obfile)
    use dictionaries
    use yaml_strings, only: f_char_ptr
    use yaml_output
    use f_utils
    implicit none
    character(len = *), intent(in) :: obfile
    type(dictionary), pointer :: dict,indict,outdict

    interface
       subroutine openbabel_load(d, f)!,ln)
         use dictionaries
         implicit none
         type(dictionary), pointer :: d
         !character(len = *), intent(in) :: f
         character, dimension(*), intent(in) :: f
         !integer, intent(in) :: ln
       end subroutine openbabel_load
       subroutine openbabel_formats(indict,outdict)
         use dictionaries, only: dictionary
         implicit none
         type(dictionary), pointer :: indict,outdict
       end subroutine openbabel_formats
    end interface
    call dict_init(dict)
    call openbabel_load(dict,f_char_ptr(trim(obfile)))!,len(obfile))
    call openbabel_formats(indict,outdict)

    call yaml_map('Supported input formats',indict)
    call yaml_map('Supported output formats',outdict)
    call dict_free(indict,outdict)

    call yaml_map('Parsed posinp dict from openbabel',dict)


  end subroutine load_dict_from_openbabel


  subroutine dump_dict_with_openbabel(dict,dict_types,fout)
    use dictionaries
    use yaml_strings, only: f_char_ptr
    implicit none
    type(dictionary), pointer :: dict,dict_types
    character(len=*), intent(in) :: fout
    !local variables

    interface
       subroutine openbabel_dump(d, dt, f)!,ln)
         use dictionaries
         implicit none
         type(dictionary), pointer :: d,dt
         !character(len = *), intent(in) :: f
         character, dimension(*), intent(in) :: f
         !integer, intent(in) :: ln
       end subroutine openbabel_dump
    end interface

    call openbabel_dump(dict,dict_types,f_char_ptr(trim(fout)))! fout,len_trim(fout))

  end subroutine dump_dict_with_openbabel

end module
