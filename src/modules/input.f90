!> @file
!!  Module to handle input variables
!! @author
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


module module_input

  use module_base

  implicit none
  private

  integer, parameter :: nmax_lines=500,max_length=100
  character(len=max_length) :: input_file,line_being_processed
  logical :: output,lmpinit
  integer :: iline_parsed,iline_written,iargument,ipos,nlines_total
  character(len=max_length), dimension(:), allocatable :: inout_lines


  interface input_var
     module procedure var_character, var_logical, var_integer, &
          & var_integer_array, var_double, var_keyword, var_ids,&
          & var_double_compulsory,var_real_compulsory,var_int_compulsory,var_logical_compulsory,&
          & var_char_compulsory
  end interface

  public :: input_set_file
  public :: input_var
  public :: input_free
  public :: case_insensitive_equiv
  public :: read_fraction_string
  public :: read_fraction_string_old

contains

  subroutine input_set_file(iproc, filename, exists,comment_file_usage)
    integer, intent(in) :: iproc
    character(len = *), intent(in) :: filename,comment_file_usage
    logical, intent(out) :: exists

    character(len=max_length), dimension(nmax_lines) :: lines
    integer :: i,nlines,ierror,ierr=0 !for MPIfake BCAST

    !no line parsed if the file not exists
    iline_parsed=0
    !no line has been written in the output
    iline_written=0
    !no argument has been read yet
    iargument=1
    !the present line is empty
    line_being_processed=repeat(' ',max_length)
    !the starting position for the lines is zero
    ipos=0
    !there are no lines for the moment
    nlines_total=0

    !verify whether MPI has been initialized
    lmpinit=.false.
    call MPI_INITIALIZED(lmpinit,ierr)

    write(input_file, "(A)") trim(filename)
           

    !check if the file is present
    inquire(file=trim(input_file),exist=exists)
    if (exists) then
       !only the root processor parse the file
       if (iproc==0) then
          open(unit = 1, file = trim(filename), status = 'old')
          i = 1
          parse_file: do 
             lines(i)=repeat(' ',max_length) !initialize lines
             read(1, fmt = '(a)', iostat = ierror) lines(i)
             !eliminate leading blanks from the line
             lines(i)=adjustl(lines(i))
             if (ierror /= 0) exit parse_file
             i = i + 1
          end do parse_file
          close(1)
          nlines=i-1
       end if
       !broadcast the number of lines
       if (lmpinit) call MPI_BCAST(nlines,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       if (ierr /=0) stop 'input_file BCAST (1) '
       nlines_total=nlines
       !broadcast all the lines
       if (lmpinit) call MPI_BCAST(lines,nmax_lines*nlines,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
       if (ierr /=0) stop 'input_file BCAST (2) '

!!$    write(0,*) "Setup input file '", trim(filename), "' with ", i - 1, "lines."

       allocate(inout_lines(0:nlines)) !the 0-th line is for the description of the file
       do i=1,nlines
          inout_lines(i)=lines(i)
       end do
       !start parsing from the first line
       iline_parsed=1
    else
       !in this case the array constitute the output results
       allocate(inout_lines(0:nmax_lines))
       do i=1,nmax_lines
          inout_lines(i)=repeat(' ',max_length) !initialize lines
       end do
       nlines_total=nmax_lines
    end if

    !write the first line in the output
    if (exists) then
       write(inout_lines(iline_written),'(1x,5a)')&
            '|--- (file:', trim(filename),')',&
            repeat('-',83-(len(trim(filename)//trim(comment_file_usage))+11)),&
            trim(comment_file_usage)
    else
       write(inout_lines(iline_written),'(1x,5a)')&
            '|--- (file:',trim(filename),'-- not present)',&
            repeat('-',83-(len(trim(filename)//trim(comment_file_usage))+26)),&
            trim(comment_file_usage)
    end if
    iline_written=iline_written+1

    output = (iproc == 0)
    !dump the 0-th line on the screen
    if (iproc == 0) then
       write(*,'(1x,a)') inout_lines(0)
    end if
!!$    output = (iproc == 0)
!!$    ! Output
!!$    if (iproc == 0) then
!!$       write(*,*)
!!$       if (exists) then
!!$          write(*,'(1x,3a)') '--- (file: ', trim(filename), &
!!$               & ') -----------------------------------------'//&
!!$               trim(comment_file_usage)
!!$       else
!!$          write(*,'(1x,a)')&
!!$               '--- (file:'//trim(filename)//'-- not present) --------------------------'//&
!!$               trim(comment_file_usage)
!!$       end if
!!$    end if

  END SUBROUTINE input_set_file

  subroutine input_free(iproc)
    implicit none
    integer, intent(in), optional :: iproc
    !Local variables
    integer, parameter :: iunit=11
    integer :: ierr,iline

    if (present(iproc)) then !case for compulsory variables
       !if (iline_written==1) iline_written=2
       if (iproc ==0) then
          if (iline_parsed==0) then !the file does not exist
             !add the writing of the file in the given unit
             open(unit=iunit,file=trim(input_file)//'_default', status ='unknown')
             do iline=1,iline_written-1
                write(iunit,*) inout_lines(iline)
             end do
             close(unit=iunit)
          end if
          !dump the file on the screen
          do iline=1,iline_written-1
          write(*,fmt='(1x,a,a)') '|',inout_lines(iline)
          end do
       end if
    end if

    if (allocated(inout_lines)) deallocate(inout_lines)
    if (lmpinit) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  END SUBROUTINE input_free

  subroutine leave()
    implicit none
    !local variables
    integer :: ierr
    if (output) then
       write(*,'(1x,a,a,2(a,i3))')'Error while reading the file "', &
            & trim(input_file), '", line=', iline_written,' argument=', iargument
       if (iline_written <= nlines_total) write(*,*)inout_lines(iline_written),line_being_processed
       !to be called only if mpi is initialized
    end if
    if (lmpinit) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    stop
  end subroutine leave

  subroutine check(ierror)
    implicit none
    !Arguments
    integer, intent(in) :: ierror

    if (ierror/=0) then
       call leave()
    end if
    !increment the argument at each passed check
    iargument=iargument+1

  END SUBROUTINE check

  !> Process the line needed with the default value in mind
  subroutine process_line(default,line_comment)
    implicit none
    character(len=*), intent(in), optional :: default,line_comment
    !Local variables
    integer :: i,iblank,istart,nchars

    if (iargument==1) then
       ipos=0
       line_being_processed=repeat(' ',max_length)
    end if

    if (present(default) .and. iline_parsed==0) then
       !case without file, write default and continue the line
       do i=1,len_trim(default)
          inout_lines(iline_written)(i+ipos:i+ipos)=default(i:i)
       end do
       ipos=ipos+len_trim(default)+1
       inout_lines(iline_written)(ipos:ipos)=' '
    end if
    if (present(line_comment) .and. iline_parsed==0) then
       !case without file, close the line. Start the comment at column 16 if possible
       istart=max(ipos+1,16)
       nchars=min(len(line_comment),max_length-istart)
       inout_lines(iline_written)(istart:istart+nchars)=line_comment(1:nchars)
       iline_written=iline_written+1
       iargument=0
    else if (.not. present(default) .and. .not. present(line_comment).and. iline_parsed/=0) then
       !traditional case, the argument should be parsed one after another
       !start with the entire line
       if (iargument==1) then 
          !print *,'prsed',iline_parsed,inout_lines(iline_parsed)
          line_being_processed=inout_lines(iline_parsed)
       else
          !search in the line the first blank
          iblank=scan(line_being_processed,' ')
          do i=1,max_length-iblank
             line_being_processed(i:i)=line_being_processed(i+iblank:i+iblank)
          end do
          do i=max_length-iblank+1,max_length
             line_being_processed(i:i)=' '
          end do
       end if
       !adjust the line to eliminate further blanks
       line_being_processed=adjustl(line_being_processed)
    else if (.not. present(default) .and. present(line_comment) .and. iline_parsed/=0) then
       !traditional case, close the line and skip to the next one
       iargument=1
       iline_parsed=iline_parsed+1
       iline_written=iline_written+1
    end if
       
  end subroutine process_line

  subroutine find(name, iline, ii)
    character(len = *), intent(in) :: name
    integer, intent(out) :: iline, ii
    
    integer :: k
    !change the allocate condition, since input lines is always used now
    !if (allocated(inout_lines)) then
    if (iline_parsed /= 0) then
       do iline = 1, size(inout_lines)-1 !there is also the zero now
          k = 1
          do ii = 1, len(inout_lines(iline)), 1
             if (ichar(inout_lines(iline)(ii:ii)) == ichar(name(k:k)) .or. &
                  & ichar(inout_lines(iline)(ii:ii)) == ichar(name(k:k)) + 32 .or. &
                  & ichar(inout_lines(iline)(ii:ii)) == ichar(name(k:k)) - 32) then
                k = k + 1
             else
                k = 1
             end if
             if (k == len(name) + 1) then
                return
             end if
          end do
       end do
    end if
    iline = 0
  END SUBROUTINE find


  !> Read a real or real/real, real:real 
  !! Here the fraction is indicated by the ':' or '/'
  !! The problem is that / is a separator for Fortran
  subroutine read_fraction_string(string,var,ierror)
     use module_base
     implicit none
     !Arguments
     character(len=*), intent(in) :: string
     real(gp), intent(out) :: var
     integer, intent(out) :: ierror
     !Local variables
     character(len=200) :: tmp
     integer :: num,den,pfr,psp

     !First look at the first blank after trim
     tmp=trim(string)
     psp = scan(tmp,' ')

     !see whether there is a fraction in the string
     if(psp==0) psp=len(tmp)
     pfr = scan(tmp(1:psp),':')
     if (pfr == 0) pfr = scan(tmp(1:psp),'/')
     !It is not a fraction
     if (pfr == 0) then
        read(tmp(1:psp),*,iostat=ierror) var
     else 
        read(tmp(1:pfr-1),*,iostat=ierror) num
        read(tmp(pfr+1:psp),*,iostat=ierror) den
        if (ierror == 0) var=real(num,gp)/real(den,gp)
     end if
     !Value by defaut
     if (ierror /= 0) var = huge(1_gp)
  END SUBROUTINE read_fraction_string


  !>  Here the fraction is indicated by the :
  subroutine read_fraction_string_old(l,string,occ)
     use module_base
     implicit none
     integer, intent(in) :: l
     character(len=*), intent(in) :: string
     real(gp), intent(out) :: occ
     !local variables
     integer :: num,den,pfr

     !see whether there is a fraction in the string
     if (l>3) then
        pfr=3
     else
        pfr=2
     end if
     if (string(pfr:pfr) == ':') then
        read(string(1:pfr-1),*)num
        read(string(pfr+1:2*pfr-1),*)den
        occ=real(num,gp)/real(den,gp)
     else
        read(string,*)occ
     end if
  END SUBROUTINE read_fraction_string_old


  !> Compare two strings (case-insensitive). Blanks are relevant!
  function case_insensitive_equiv(stra,strb)
    implicit none
    character(len=*), intent(in) :: stra,strb
    logical :: case_insensitive_equiv
    !Local variables
    integer :: i,ica,icb,ila,ilb,ilength
    ila=len(stra)
    ilb=len(strb)
    ilength=min(ila,ilb)
    ica=ichar(stra(1:1))
    icb=ichar(strb(1:1))
    case_insensitive_equiv=(modulo(ica-icb,32) == 0) .and. (ila==ilb)
    do i=2,ilength
       ica=ichar(stra(i:i))
       icb=ichar(strb(i:i))
       case_insensitive_equiv=case_insensitive_equiv .and. &
            (modulo(ica-icb,32) == 0)
       if (.not. case_insensitive_equiv) exit
    end do

  end function case_insensitive_equiv


!--routines for compulsory file

  subroutine var_double_compulsory(var,default,ranges,exclusive,comment,input_iostat)
    implicit none
    character(len=*), intent(in) :: default
    real(kind=8), intent(out) :: var
    character(len=*), intent(in), optional :: comment
    integer, intent(out), optional :: input_iostat
    real(kind=8), dimension(2), intent(in), optional :: ranges
    real(kind=8), dimension(:), intent(in), optional :: exclusive
    !Local variables
    logical :: found
    integer :: ierror,ilist

    if (present(input_iostat)) then
       !first, check if the line is correct
       if (iline_written>nlines_total ) then
          input_iostat=-1
          return
       else
          input_iostat=0 !no error for the moment
       end if
    else
       if (iline_written>nlines_total) then
          call leave()
       end if
    end if

    !if the file has not been opened, use the default variable 
    !then write in the output lines the default
    if (iline_parsed==0) then
       !finalize the line if the comment is present
       if (present(comment)) then
          call process_line(default=default,line_comment=comment)
       else
          call process_line(default=default)
       end if
       call read_fraction_string(default,var,ierror)
       call check(ierror)
    !otherwise read the corresponding argument and check its validity
    else
       !read the argument
       call process_line()
       !print *,line_being_processed
       call read_fraction_string(line_being_processed,var,ierror)
       call check(ierror)

       !check the validity of the variable
       if (present(ranges)) then
          if (var < ranges(1) .or. var > ranges(2)) then
             if (output) then
                write(*,'(1x,a,i0,a,i0)') &
                'ERROR in parsing file '//trim(input_file)//', line=', iline_written,' argument=', iargument-1
                write(*,*)'      values should be in range: [',ranges(1),'-',ranges(2),']'
             end if
             if (present(input_iostat)) then
                input_iostat=1
                return
             else
                call leave()
             end if
          end if
       else if (present(exclusive)) then
          found=.false.
          found_loop: do ilist=1,size(exclusive)
             if (var == exclusive(ilist)) then
                found=.true.
                exit found_loop
             end if
          end do found_loop
          if (.not. found) then
             if (output) then
                write(*,'(1x,a,i0,a,i0)') &
                     'ERROR in parsing file '//trim(input_file)//', line=', iline_written,' argument=', iargument-1
                write(*,*)'      values should be in list: ',exclusive(:)
             end if
             if (present(input_iostat)) then
                input_iostat=1
                return
             else
                call leave()
             end if
          end if
       end if 
      
       !increment the line if comment is present, do not touch the input file
       if (present(comment)) then
          call process_line(line_comment=comment)
       end if
    end if   
  END SUBROUTINE var_double_compulsory

  subroutine var_real_compulsory(var,default,ranges,exclusive,comment,input_iostat)
    implicit none
    character(len=*), intent(in) :: default
    real(kind=4), intent(out) :: var
    character(len=*), intent(in), optional :: comment
    integer, intent(out), optional :: input_iostat
    real(kind=4), dimension(2), intent(in), optional :: ranges
    real(kind=4), dimension(:), intent(in), optional :: exclusive
    !Local variables
    real(gp) :: double_var
    logical :: found
    integer :: ierror,ilist

    if (present(input_iostat)) then
       !first, check if the line is correct
       if (iline_written>nlines_total ) then
          input_iostat=-1
          return
       else
          input_iostat=0 !no error for the moment
       end if
    else
       if (iline_written>nlines_total) then
          call leave()
       end if
    end if

    !if the file has not been opened, use the default variable 
    !then write in the output lines the default
    if (iline_parsed==0) then
       !finalize the line if the comment is present
       if (present(comment)) then
          call process_line(default=default,line_comment=comment)
       else
          call process_line(default=default)
       end if
       call read_fraction_string(default,double_var,ierror)
       call check(ierror)
       var=real(double_var,kind=4)
    !otherwise read the corresponding argument and check its validity
    else
       !read the argument
       call process_line()
       call read_fraction_string(line_being_processed,double_var,ierror)
       call check(ierror)
       var=real(double_var,kind=4)

       !check the validity of the variable
       if (present(ranges)) then
          if (var < ranges(1) .or. var > ranges(2)) then
             if (output) then
                write(*,'(1x,a,i0,a,i0)') &
                     'ERROR in parsing file '//trim(input_file)//', line=', iline_written,' argument=', iargument-1
                write(*,*)'      values should be in range: [',ranges(1),'-',ranges(2),']'
             end if
             if (present(input_iostat)) then
                input_iostat=1
                return
             else
                call leave()
             end if
          end if
       else if (present(exclusive)) then
          found=.false.
          found_loop: do ilist=1,size(exclusive)
             if (var == exclusive(ilist)) then
                found=.true.
                exit found_loop
             end if
          end do found_loop
          if (.not. found) then
             if (output) then 
                write(*,'(1x,a,i0,a,i0)') &
                     'ERROR in parsing file '//trim(input_file)//', line=', iline_written,' argument=', iargument-1
                write(*,*)'      values should be in list: ',exclusive(:)
             end if
             if (present(input_iostat)) then
                input_iostat=1
                return
             else
                call leave()
             end if
          end if
       end if 
      
       !increment the line if comment is present, do not touch the input file
       if (present(comment)) then
          call process_line(line_comment=comment)
       end if
    end if   
  END SUBROUTINE var_real_compulsory

  subroutine var_int_compulsory(var,default,ranges,exclusive,comment,input_iostat)
    implicit none
    character(len=*), intent(in) :: default
    integer, intent(out) :: var
    character(len=*), intent(in), optional :: comment
    integer, intent(out), optional :: input_iostat
    integer, dimension(2), intent(in), optional :: ranges
    integer, dimension(:), intent(in), optional :: exclusive
    !Local variables
    logical :: found
    integer :: ierror,ilist

    if (present(input_iostat)) then
       !first, check if the line is correct
       if (iline_written>nlines_total) then
          input_iostat=-1
          return
       else
          input_iostat=0 !no error for the moment
       end if
    else
       if (iline_written>nlines_total) then
          call leave()
       end if
    end if

    !if the file has not been opened, use the default variable 
    !then write in the output lines the default
    if (iline_parsed==0) then
       !finalize the line if the comment is present
       if (present(comment)) then
          call process_line(default=default,line_comment=comment)
       else
          call process_line(default=default)
       end if
       read(default,*,iostat=ierror)var
       call check(ierror)
    !otherwise read the corresponding argument and check its validity
    else
       !read the argument
       call process_line()

       read(line_being_processed,fmt=*,iostat=ierror) var
       call check(ierror)

       !check the validity of the variable
       if (present(ranges)) then
          if (var < ranges(1) .or. var > ranges(2)) then
             if (output) then
                write(*,'(1x,a,i0,a,i0)') &
                     'ERROR in parsing file '//trim(input_file)//', line=', iline_written,' argument=', iargument-1
                write(*,*)'      values should be in range: [',ranges(1),'-',ranges(2),']'
             end if
             if (present(input_iostat)) then
                input_iostat=1
                return
             else
                call leave()
             end if
          end if
       else if (present(exclusive)) then
          found=.false.
          found_loop: do ilist=1,size(exclusive)
             if (var == exclusive(ilist)) then
                found=.true.
                exit found_loop
             end if
          end do found_loop
          if (.not. found) then
             if (output) then
                write(*,'(1x,a,i0,a,i0)') &
                     'ERROR in parsing file '//trim(input_file)//', line=', iline_written,' argument=', iargument-1
                write(*,*)'      values should be in list: ',exclusive(:)
             end if
             if (present(input_iostat)) then
                input_iostat=1
                return
             else
                call leave()
             end if
          end if
       end if 
      
       !increment the line if comment is present, do not touch the input file
       if (present(comment)) then
          call process_line(line_comment=comment)
       end if
    end if   
  END SUBROUTINE var_int_compulsory

  subroutine var_char_compulsory(var,default,exclusive,comment,input_iostat)
    implicit none
    character(len=*), intent(in) :: default
    character(len=*), intent(out) :: var
    character(len=*), intent(in), optional :: comment
    integer, intent(out), optional :: input_iostat
    character(len=*), dimension(:), intent(in), optional :: exclusive
    !Local variables
    logical :: found
    integer :: ierror,ilist

    if (present(input_iostat)) then
       !first, check if the line is correct (or if it is an optional line)
       if (iline_written>nlines_total) then
          input_iostat=-1
          return
       else
          input_iostat=0 !no error for the moment
       end if
    else
       if (iline_written>nlines_total) then
          call leave()
       end if
    end if

    !if the file has not been opened, use the default variable 
    !then write in the output lines the default
    if (iline_parsed==0) then
       !finalize the line if the comment is present
       if (present(comment)) then
          call process_line(default=default,line_comment=comment)
       else
          call process_line(default=default)
       end if
       read(default,*,iostat=ierror)var
       call check(ierror)
    !otherwise read the corresponding argument and check its validity
    else
       !read the argument
       call process_line()
       read(line_being_processed,fmt=*,iostat=ierror) var
       call check(ierror)

       if (present(exclusive)) then
          found=.false.
          found_loop: do ilist=1,size(exclusive)
             if (case_insensitive_equiv(trim(var),trim(exclusive(ilist)))) then
                found=.true.
                exit found_loop
             end if
          end do found_loop
          if (.not. found) then
             if (output) then 
                write(*,'(1x,a,i0,a,i0)') &
                     'ERROR in parsing file '//trim(input_file)//', line=', iline_written,' argument=', iargument-1
                write(*,'(6x,a,30(1x,a))')&
                     'values should be in list: ',exclusive(:)
             end if
             if (present(input_iostat)) then
                input_iostat=1
                return
             else
                call leave()
             end if
          end if
       end if 
      
       !increment the line if comment is present, do not touch the input file
       if (present(comment)) then
          call process_line(line_comment=comment)
       end if
    end if   
  END SUBROUTINE var_char_compulsory

  subroutine var_logical_compulsory(var,default,comment)
    implicit none
    character(len=*), intent(in) :: default
    logical, intent(out) :: var
    character(len=*), intent(in), optional :: comment
    !Local variables
    integer :: ierror

    !if the file has not been opened, use the default variable 
    !then write in the output lines the default
    if (iline_parsed==0) then
       !finalize the line if the comment is present
       if (present(comment)) then
          call process_line(default=default,line_comment=comment)
       else
          call process_line(default=default)
       end if
       read(default,*,iostat=ierror)var
       call check(ierror)
    !otherwise read the corresponding argument and check its validity
    else
       !read the argument
       call process_line()

       read(line_being_processed,fmt=*,iostat=ierror) var
       call check(ierror)
      
       !increment the line if comment is present, do not touch the input file
       if (present(comment)) then
          call process_line(line_comment=comment)
       end if
    end if   
  END SUBROUTINE var_logical_compulsory

!-routines for non-compulsory file (input.perf)
  subroutine var_character(name, default, description, var)
    character(len = *), intent(in) :: name
    character(len = *), intent(in) :: default
    character(len = *), intent(in) :: description
    character(len = *), intent(out) :: var

    integer :: i, j, ierror, ierr

    write(var, "(A)") default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,a,t30,2a)") "|", name, var, '!', description
  END SUBROUTINE var_character

  subroutine var_logical(name, default, description, var)
    character(len = *), intent(in) :: name
    logical, intent(in) :: default
    character(len = *), intent(in) :: description
    logical, intent(out) :: var

    integer :: i, j, ierror

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror /= 0) then
          var = .true.
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,l1,t30,2a)") "|", name, var, '!', description
  END SUBROUTINE var_logical

  subroutine var_integer(name, default, description, var)
    character(len = *), intent(in) :: name
    integer, intent(in) :: default
    character(len = *), intent(in) :: description
    integer, intent(out) :: var

    integer :: i, j, ierror, ierr

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,I0,t30,2a)") "|", name, var, '!', description
  END SUBROUTINE var_integer

  subroutine var_integer_array(name, default, description, var)
    character(len = *), intent(in) :: name
    integer, intent(in) :: default(:)
    character(len = *), intent(in) :: description
    integer, intent(out) :: var(:)

    integer :: i, j, ierror, ierr

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) then
       write(*,"(1x,a,3x,a,1x)", advance = "NO") "|", name
       do i = 1, size(var), 1
          write(*,"(1x,I0)", advance = "NO") var(i)
       end do
       write(*,"(t7,2a)") '!', description
    end if
  END SUBROUTINE var_integer_array

  subroutine var_double(name, default, description, var)
    character(len = *), intent(in) :: name
    double precision, intent(in) :: default
    character(len = *), intent(in) :: description
    double precision, intent(out) :: var

    integer :: i, j, ierror, ierr

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,es9.2,t30,2a)") "|", name, var, '!', description
  END SUBROUTINE var_double

  subroutine var_keyword(name, length, default, list, description, var)
    character(len = *), intent(in) :: name
    integer, intent(in) :: length
    character(len = length), intent(in) :: default
    character(len = length), intent(in) :: list(:)
    character(len = *), intent(in) :: description
    integer, intent(out) :: var

    integer :: i, j, ierror, ierr
    character(len = length) :: buf

    ! Set the default value to var.
    do i = 1, size(list), 1
       if (trim(default) == trim(list(i))) exit
    end do
    var = i - 1
    ! Find the keyword name in the file.
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) buf
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
       ! Look for buf in list.
       do j = 1, size(list), 1
          if (trim(buf) == trim(list(j))) exit
       end do
       if (j > size(list)) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
       var = j - 1
    end if
    if (output) then
       write(*,"(1x,a,3x,a,1x,a,t30,3a)", advance = "NO") &
            & "|", name, list(var + 1), '!', description, " ("
       write(*,"(A)", advance = "NO") trim(list(1))
       do i = 2, size(list), 1
          write(*,"(2A)", advance = "NO") ", ", trim(list(i))
       end do
       write(*,"(A)") ")"
    end if
  END SUBROUTINE var_keyword

  subroutine var_ids(name, default, list, description, var)
    character(len = *), intent(in) :: name
    integer, intent(in) :: default
    integer, intent(in) :: list(:)
    character(len = *), intent(in) :: description
    integer, intent(out) :: var

    integer :: i, j, ierror, ierr
    
    var = default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
       do j = 1, size(list), 1
          if (var == list(j)) exit
       end do
       if (j > size(list)) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,I0,t30,2a)") "|", name, var, '!', description
  END SUBROUTINE var_ids

end module module_input
