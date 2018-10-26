!!$subroutine read_wave_to_isf(lstat, filename, ln, iorbp, hx, hy, hz, &
!!$     & n1, n2, n3, nspinor, psiscf)
!!$  use module_base
!!$  use module_types
!!$  use module_interfaces, only: readwavetoisf, readwavetoisf_etsf
!!$  use public_enums
!!$  use module_input_keys
!!$  implicit none
!!$
!!$  integer, intent(in) :: ln
!!$  character, intent(in) :: filename(ln)
!!$  integer, intent(in) :: iorbp
!!$  integer, intent(out) :: n1, n2, n3, nspinor
!!$  real(gp), intent(out) :: hx, hy, hz
!!$  real(wp), dimension(:,:,:,:), pointer :: psiscf
!!$  logical, intent(out) :: lstat
!!$
!!$  character(len = 1024) :: filename_
!!$  integer :: iformat, i
!!$
!!$  write(filename_, "(A)") " "
!!$  do i = 1, ln, 1
!!$     filename_(i:i) = filename(i)
!!$  end do
!!$
!!$  ! Find format from name.
!!$  iformat = wave_format_from_filename(1, trim(filename_))
!!$
!!$  ! Call proper wraping routine.
!!$  if (iformat == WF_FORMAT_ETSF) then
!!$     call readwavetoisf_etsf(lstat, trim(filename_), iorbp, hx, hy, hz, &
!!$          & n1, n2, n3, nspinor, psiscf)
!!$  else
!!$     call readwavetoisf(lstat, trim(filename_), (iformat == WF_FORMAT_PLAIN), hx, hy, hz, &
!!$          & n1, n2, n3, nspinor, psiscf)
!!$  end if
!!$END SUBROUTINE read_wave_to_isf

!!$subroutine readwavetoisf(lstat, filename, formatted, hx, hy, hz, &
!!$     & n1, n2, n3, nspinor, psiscf)
!!$  use module_base
!!$!  use module_types
!!$  use io, only: io_open, io_read_descr, io_warning, read_psi_compress, io_gcoordtolocreg
!!$  use locregs
!!$  use locreg_operations
!!$  implicit none
!!$
!!$  character(len = *), intent(in) :: filename
!!$  logical, intent(in) :: formatted
!!$  integer, intent(out) :: n1, n2, n3, nspinor
!!$  real(gp), intent(out) :: hx, hy, hz
!!$  real(wp), dimension(:,:,:,:), pointer :: psiscf
!!$  logical, intent(out) :: lstat
!!$
!!$  character(len = *), parameter :: subname = "readwavetoisf"
!!$  integer :: unitwf, iorb, ispinor, ispin, ikpt
!!$  integer, dimension(:,:), allocatable :: gcoord_c, gcoord_f
!!$  real(wp) :: eval
!!$  real(wp), dimension(:), allocatable :: psi
!!$  type(locreg_descriptors) :: lr
!!$  character(len = 256) :: error
!!$  type(workarr_sumrho) :: w
!!$  character(len = 1024) :: fileRI
!!$  integer :: n1_old, n2_old, n3_old, nvctr_c_old, nvctr_f_old
!!$  real(gp) :: hx_old, hy_old, hz_old
!!$
!!$  call f_routine(id=subname)
!!$
!!$
!!$  ! We open the Fortran file
!!$  call io_open(unitwf, filename, formatted)
!!$  if (unitwf < 0) then
!!$     call f_release_routine()
!!$     return
!!$  end if
!!$
!!$  ! We read the basis set description and the atomic definition.
!!$  call io_read_descr(unitwf, formatted, iorb, eval, n1, n2, n3, &
!!$       & hx, hy, hz, lstat, error,  nvctr_c_old, nvctr_f_old)
!!$  if (.not. lstat) then
!!$     call io_warning(trim(error))
!!$     call f_release_routine()
!!$     return
!!$  end if
!!$  ! Do a magic here with the filenames...
!!$  call readwavedescr(lstat, filename, iorb, ispin, ikpt, ispinor, nspinor, fileRI)
!!$  if (.not. lstat) then
!!$     call io_warning("cannot read wave ids from filename.")
!!$     call f_release_routine()
!!$     return
!!$  end if
!!$
!!$  ! Initial allocations.
!!$  gcoord_c = f_malloc((/ 3, nvctr_c_old  /),id='gcoord_c')
!!$  gcoord_f = f_malloc((/ 3, nvctr_f_old  /),id='gcoord_f')
!!$  psi = f_malloc( nvctr_c_old + 7 * nvctr_f_old ,id='psi')
!!$  ! Read psi and the basis-set
!!$  call read_psi_compress(unitwf, formatted, nvctr_c_old, nvctr_f_old, psi, lstat, error, gcoord_c, gcoord_f)
!!$  if (.not. lstat) then
!!$     call io_warning(trim(error))
!!$     call deallocate_local()
!!$     return
!!$  end if
!!$  call io_gcoordToLocreg(n1, n2, n3, nvctr_c_old, nvctr_f_old, &
!!$       & gcoord_c, gcoord_f, lr)
!!$
!!$  psiscf = f_malloc_ptr((/ lr%d%n1i, lr%d%n2i, lr%d%n3i, nspinor  /),id='psiscf')
!!$
!!$  call initialize_work_arrays_sumrho(lr,.true.,w)
!!$
!!$  ! Magic-filter to isf
!!$  call daub_to_isf(lr, w, psi, psiscf(1,1,1,ispinor))
!!$
!!$  ! Read the other psi part, if any
!!$  if (nspinor > 1) then
!!$     close(unitwf)
!!$     n1_old = n1
!!$     n2_old = n2
!!$     n3_old = n3
!!$     hx_old = hx
!!$     hy_old = hy
!!$     hz_old = hz
!!$     nvctr_c_old = lr%wfd%nvctr_c
!!$     nvctr_f_old = lr%wfd%nvctr_f
!!$
!!$     ispinor = modulo(ispinor, 2) + 1
!!$     call io_open(unitwf, trim(fileRI), formatted)
!!$     if (unitwf < 0) then
!!$        call io_warning("cannot read other spinor part from '" // trim(fileRI) // "'.")
!!$        call deallocate_local()
!!$        return
!!$     end if
!!$     ! We read the basis set description and the atomic definition.
!!$     call io_read_descr(unitwf, formatted, iorb, eval, n1, n2, n3, &
!!$          & hx, hy, hz, lstat, error, lr%wfd%nvctr_c, lr%wfd%nvctr_f)
!!$     if (.not. lstat) then
!!$        call io_warning(trim(error))
!!$        call deallocate_local()
!!$        return
!!$     end if
!!$
!!$     ! Check consistency of the basis-set.
!!$     if (n1_old == n1 .and. n2_old == n2 .and. n3_old == n3 .and. &
!!$          & hx_old == hx .and. hy_old == hy .and. hz_old == hz .and. &
!!$          & nvctr_c_old == lr%wfd%nvctr_c .and. nvctr_f_old == lr%wfd%nvctr_f) then
!!$        call read_psi_compress(unitwf, formatted, lr%wfd%nvctr_c, lr%wfd%nvctr_f, psi, lstat, error)
!!$        if (.not. lstat) then
!!$           call io_warning(trim(error))
!!$           call deallocate_local()
!!$           return
!!$        end if
!!$        call daub_to_isf(lr, w, psi, psiscf(1,1,1,ispinor))
!!$     else
!!$        call io_warning("It exists a file with the same naming convention" // &
!!$             & " but with a different basis-set.")
!!$        hx = hx_old
!!$        hy = hy_old
!!$        hz = hz_old
!!$        psiscf(:,:,:,ispinor) = real(0, wp)
!!$     end if
!!$  end if
!!$
!!$  ! We update the size values to match the allocation of psiscf.
!!$  n1 = lr%d%n1i
!!$  n2 = lr%d%n2i
!!$  n3 = lr%d%n3i
!!$  hx = hx * 0.5d0
!!$  hy = hy * 0.5d0
!!$  hz = hz * 0.5d0
!!$
!!$  call deallocate_local()
!!$  lstat = .true.
!!$
!!$contains
!!$
!!$  subroutine deallocate_local()
!!$    implicit none
!!$
!!$    ! We close the file.
!!$    close(unit=unitwf)
!!$
!!$    !allocation status of a allocatable array is undefined, cannot do that
!!$    !as we do not have any way to deallocate an array
!!$    call f_free(psi)
!!$    call f_free(gcoord_c)
!!$    call f_free(gcoord_f)
!!$
!!$    if (associated(w%x_c)) then
!!$       call deallocate_work_arrays_sumrho(w)
!!$    end if
!!$    call deallocate_locreg_descriptors(lr)
!!$    !call deallocate_convolutions_bounds(lr%bounds)
!!$    !call deallocate_wfd(lr%wfd)
!!$
!!$    call f_release_routine()
!!$  END SUBROUTINE deallocate_local
!!$
!!$END SUBROUTINE readwavetoisf

!!$    subroutine io_gcoordToLocreg(n1, n2, n3, nvctr_c, nvctr_f, gcoord_c, gcoord_f, lr)
!!$      use module_base
!!$      use locregs
!!$      implicit none
!!$      !Arguments
!!$      integer, intent(in) :: n1, n2, n3, nvctr_c, nvctr_f
!!$      integer, dimension(3, nvctr_c), intent(in) :: gcoord_c
!!$      integer, dimension(3, nvctr_f), intent(in) :: gcoord_f
!!$      type(locreg_descriptors), intent(out) :: lr
!!$      !Local variables
!!$      character(len = *), parameter :: subname = "io_gcoordToLocreg"
!!$      integer :: i
!!$      logical, dimension(:,:,:), allocatable :: logrid_c, logrid_f
!!$  
!!$      call f_routine(id=subname)
!!$  
!!$      call nullify_locreg_descriptors(lr)
!!$  
!!$      lr%geocode = "P"
!!$      lr%hybrid_on = .false.
!!$  
!!$      lr%ns1 = 0
!!$      lr%ns2 = 0
!!$      lr%ns3 = 0
!!$  
!!$      lr%d%n1 = n1
!!$      lr%d%n2 = n2
!!$      lr%d%n3 = n3
!!$  
!!$      lr%d%n1i = 2 * n1 + 2
!!$      lr%d%n2i = 2 * n2 + 2
!!$      lr%d%n3i = 2 * n3 + 2
!!$  
!!$      logrid_c = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_c')
!!$      logrid_f = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_f')
!!$  
!!$      lr%d%nfl1 = n1
!!$      lr%d%nfl2 = n2
!!$      lr%d%nfl3 = n3
!!$      lr%d%nfu1 = 0
!!$      lr%d%nfu2 = 0
!!$      lr%d%nfu3 = 0
!!$  
!!$      logrid_c(:,:,:) = .false.
!!$      do i = 1, nvctr_c, 1
!!$         logrid_c(gcoord_c(1, i), gcoord_c(2, i), gcoord_c(3, i)) = .true.
!!$      end do
!!$      logrid_f(:,:,:) = .false.
!!$      do i = 1, nvctr_f, 1
!!$         logrid_f(gcoord_f(1, i), gcoord_f(2, i), gcoord_f(3, i)) = .true.
!!$         lr%d%nfl1 = min(lr%d%nfl1, gcoord_f(1, i))
!!$         lr%d%nfl2 = min(lr%d%nfl2, gcoord_f(2, i))
!!$         lr%d%nfl3 = min(lr%d%nfl3, gcoord_f(3, i))
!!$         lr%d%nfu1 = max(lr%d%nfu1, gcoord_f(1, i))
!!$         lr%d%nfu2 = max(lr%d%nfu2, gcoord_f(2, i))
!!$         lr%d%nfu3 = max(lr%d%nfu3, gcoord_f(3, i))
!!$      end do
!!$  
!!$      !correct the values of the delimiter if there are no wavelets
!!$      if (lr%d%nfl1 == n1 .and. lr%d%nfu1 == 0) then
!!$         lr%d%nfl1 = n1 / 2
!!$         lr%d%nfu1 = n1 / 2
!!$      end if
!!$      if (lr%d%nfl2 == n2 .and. lr%d%nfu2 == 0) then
!!$         lr%d%nfl2 = n2 / 2
!!$         lr%d%nfu2 = n2 / 2
!!$      end if
!!$      if (lr%d%nfl3 == n3 .and. lr%d%nfu3 == 0) then
!!$         lr%d%nfl3 = n3 / 2
!!$         lr%d%nfu3 = n3 / 2
!!$      end if
!!$  
!!$      call wfd_from_grids(logrid_c, logrid_f, .true., lr)
!!$  
!!$      call f_free(logrid_c)
!!$      call f_free(logrid_f)
!!$  
!!$      call f_release_routine()
!!$  
!!$    END SUBROUTINE io_gcoordToLocreg
