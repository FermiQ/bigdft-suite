!> @file
!! Datatypes and associated methods relative to the localization regions (mesh grid)
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Datatypes for localization regions descriptors
module locregs
  use module_base
  use box, only: cell,box_iterator
  use compression
  use bounds, only: convolutions_bounds
  implicit none

  private

  integer, parameter :: SIZE_=1,LR_=2

  !> Grid dimensions in all different wavelet basis
  type, public :: grid_dimensions
     integer :: n1   = 0
     integer :: n2   = 0
     integer :: n3   = 0 !< Coarse grid dimensions
     integer :: nfl1 = 0
     integer :: nfu1 = 0
     integer :: nfl2 = 0
     integer :: nfu2 = 0
     integer :: nfl3 = 0
     integer :: nfu3 = 0 !< Lower and upper indices of fine grid in 3D
     integer :: n1i  = 0
     integer :: n2i  = 0
     integer :: n3i  = 0 !< ISF grid dimension (roughly 2*n+buffer)
  end type grid_dimensions

  !> Contains the information needed for describing completely a wavefunction localisation region
  type, public :: locreg_descriptors
!!$     character(len=1) :: geocode            !< @copydoc poisson_solver::doc::geocode
     logical :: hybrid_on = .false.                  !< Interesting for global, periodic, localisation regions
     integer :: ns1 = 0
     integer :: ns2 = 0
     integer :: ns3 = 0                     !< Starting point of the localisation region in global coordinates
     integer :: nsi1 = 0
     integer :: nsi2 = 0
     integer :: nsi3 = 0                    !< Starting point of locreg for interpolating grid
     integer :: Localnorb = 0                  !< Number of orbitals contained in locreg
     integer, dimension(3) :: outofzone = 0     !< Vector of points outside of the zone outside Glr for periodic systems
     real(gp), dimension(3) :: locregCenter = 0.0_gp !< Center of the locreg
     real(gp) :: locrad = 0.0_gp                    !< Cutoff radius of the localization region
     real(gp) :: locrad_kernel = 0.0_gp             !< Cutoff radius of the localization region (kernel)
     real(gp) :: locrad_mult = 0.0_gp               !< Cutoff radius of the localization region for the sparse matrix multiplications
     type(grid_dimensions) :: d               !< Grid dimensions in old different wavelet basis
     type(wavefunctions_descriptors) :: wfd
     type(convolutions_bounds) :: bounds
     type(cell) :: mesh !<defines the cell of the system
                        !! (should replace the other geometrical informations)
     !> grid in the fine scaling functions box
     type(cell) :: mesh_fine
     type(cell) :: mesh_coarse !<discretization of the coarse egrees of freedom
     !>iterator over the mesh degrees of freedom
     type(box_iterator) :: bit
  end type locreg_descriptors

  !>storage of localization regions, to be used for
  !!communicating the descriptors
  type, public :: locregs_ptr
     type(locreg_descriptors), pointer :: lr=>null()
     logical :: owner = .false.
  end type locregs_ptr

  type, public :: locreg_storage
     type(locregs_ptr), dimension(:), pointer  :: lrs_ptr=>null()
     integer, dimension(:), pointer :: encode_buffer => null()
     integer, dimension(:,:), pointer :: lr_full_sizes => null()
     integer :: lr_size = 0
  end type locreg_storage

  interface assignment(=)
     module procedure allocate_locregs_ptr
  end interface assignment(=)

  public :: nullify_locreg_descriptors,locreg_null
  public :: deallocate_locreg_descriptors,deallocate_wfd
  public :: allocate_wfd,copy_locreg_descriptors,copy_grid_dimensions
  public :: check_overlap,check_overlap_cubic_periodic,check_overlap_from_descriptors_periodic,lr_box
  public :: init_lr,reset_lr,extract_lr,store_lr,steal_lr,gather_locreg_storage
  public :: lr_storage_init, lr_storage_free
  public :: communicate_locreg_descriptors_basics
  public :: get_isf_offset,ensure_locreg_bounds

contains


  !> Create a grid_dimensions instance with all dimensions equal to 0.
  pure function grid_null() result(g)
    type(grid_dimensions) :: g
  end function grid_null


  pure function locreg_null() result(lr)
    implicit none
    type(locreg_descriptors) :: lr
    call nullify_locreg_descriptors(lr)
  end function locreg_null


  pure subroutine nullify_locreg_descriptors(lr)
    use box
    use bounds, only: nullify_convolutions_bounds
    implicit none
    type(locreg_descriptors), intent(out) :: lr
    lr%hybrid_on=.false.
    lr%ns1=0
    lr%ns2=0
    lr%ns3=0
    lr%nsi1=0
    lr%nsi2=0
    lr%nsi3=0
    lr%Localnorb=0
    lr%outofzone=(/0,0,0/)
         lr%d=grid_null()
    call nullify_wfd(lr%wfd)
    call nullify_convolutions_bounds(lr%bounds)
    lr%locregCenter=(/0.0_gp,0.0_gp,0.0_gp/)
    lr%locrad_kernel = 0.0_gp
    lr%locrad_mult = 0.0_gp
    lr%locrad=0.0_gp
    lr%mesh=cell_null()
    lr%mesh_fine=cell_null()
    lr%mesh_coarse=cell_null()
    call nullify_box_iterator(lr%bit)
  end subroutine nullify_locreg_descriptors

  pure subroutine nullify_lr_storage(lr_storage)
    implicit none
    type(locreg_storage), intent(out) :: lr_storage
  end subroutine nullify_lr_storage

  !> Destructors
  subroutine deallocate_locreg_descriptors(lr)
    use bounds
    implicit none
    ! Calling arguments
    type(locreg_descriptors),intent(inout):: lr

    call deallocate_wfd(lr%wfd)
    call deallocate_convolutions_bounds(lr%bounds)

  end subroutine deallocate_locreg_descriptors


  subroutine nullify_lr_pointers(lr)
    use bounds
    use compression
    implicit none
    type(locreg_descriptors), intent(inout) :: lr

    !nullify pointers internal to the structure to avoid fake deallocation
    call nullify_wfd_pointers(lr%wfd)
    call nullify_convolutions_bounds(lr%bounds)
  end subroutine nullify_lr_pointers

  function lr_ptr_sizeof(array)
    implicit none
    type(locregs_ptr), dimension(:), intent(in) :: array
    integer :: lr_ptr_sizeof
    !local variables
    type(locregs_ptr), dimension(2) :: lrs_ptr

    if (size(array) > 1) then
       lr_ptr_sizeof=int(f_loc(array(2))-f_loc(array(1)))
    else
       lr_ptr_sizeof=int(f_loc(lrs_ptr(2))-f_loc(lrs_ptr(1)))
    end if
  end function lr_ptr_sizeof

  subroutine allocate_locregs_ptr(array,m)
    use dynamic_memory
    implicit none
    type(locregs_ptr), dimension(:), pointer, intent(inout) :: array
    type(malloc_information_ptr), intent(in) :: m
    !local variables
    integer :: ierror

    call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)

    allocate(array(m%lbounds(1):m%ubounds(1)),stat=ierror)

    if (.not. malloc_validate(ierror,size(shape(array)),m)) return

    !here the database for the allocation might be updated
    call update_allocation_database(f_loc(array),&
         product(int(m%shape(1:m%rank),f_long)),lr_ptr_sizeof(array),m)

    call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS

  end subroutine allocate_locregs_ptr

  subroutine locregs_ptr_free(array)
    implicit none
    type(locregs_ptr), dimension(:), pointer, intent(inout) :: array
    !local variables
    integer :: ierror, i

    ! Deallocate the locregs we are owner of.
    do i = 1, size(array)
       if (array(i)%owner) then
          call deallocate_locreg_descriptors(array(i)%lr)
          deallocate(array(i)%lr)
          nullify(array(i)%lr)
       end if
    end do

    ! Deallocate container.
    call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)

    call f_purge_database(product(int(shape(array),f_long)),lr_ptr_sizeof(array),&
         f_loc(array))

    deallocate(array,stat=ierror)

    if (.not. free_validate(ierror)) return
    nullify(array)
    call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
  end subroutine locregs_ptr_free

  subroutine lr_storage_init(lr_storage,nlr)
    implicit none
    type(locreg_storage), intent(out) :: lr_storage
    integer, intent(in) :: nlr

    lr_storage%lrs_ptr=f_malloc_ptr(nlr,id='lrs_ptr')
    lr_storage%lr_size=get_locreg_encode_size()
  end subroutine lr_storage_init

  subroutine lr_storage_free(lr_storage)
    implicit none
    type(locreg_storage), intent(inout) :: lr_storage
    call locregs_ptr_free(lr_storage%lrs_ptr)
    call f_free_ptr(lr_storage%encode_buffer)
    call f_free_ptr(lr_storage%lr_full_sizes)
    call nullify_lr_storage(lr_storage)
  end subroutine lr_storage_free

  !might be used in favour of copying of the datatypes
  function get_locreg_encode_size() result(s)
    implicit none
    integer :: s
    !local variable
    type(locreg_descriptors) :: lr
    integer, dimension(3) :: dummy_char_mold

    s=size(transfer(lr,dummy_char_mold))

  end function get_locreg_encode_size

  function get_locreg_full_encode_size(lr) result(s)
    implicit none
    type(locreg_descriptors), intent(in) :: lr
    integer :: s

    s=get_locreg_encode_size()
    if (associated(lr%wfd%buffer)) s=s+size(lr%wfd%buffer)

  end function get_locreg_full_encode_size

  subroutine locreg_encode(lr,lr_size,dest)
    implicit none
    type(locreg_descriptors), intent(in) :: lr
    integer, intent(in) :: lr_size !<obtained from locreg_encode_size
    !array of dimension at least equal to locreg_encode_size
    integer, dimension(lr_size), intent(out) :: dest

    dest=transfer(lr,dest)
  end subroutine locreg_encode

  subroutine locreg_full_encode(lr,lr_size,lr_full_size,dest)
    implicit none
    type(locreg_descriptors), intent(in) :: lr
    integer, intent(in) :: lr_size,lr_full_size !<obtained from locreg_encode_size
    !array of dimension at least equal to locreg_full_encode_size
    integer, dimension(lr_full_size), intent(out) :: dest

    call locreg_encode(lr,lr_size,dest)
    if (lr_full_size-lr_size > 0) &
         & call f_memcpy(n=lr_full_size-lr_size,src=lr%wfd%buffer,dest=dest(lr_size+1))
  end subroutine locreg_full_encode

  subroutine locregs_encode(llr,nlr,lr_size,dest_arr,mask)
    implicit none
    integer, intent(in) :: lr_size,nlr
    type(locreg_descriptors), dimension(nlr), intent(in) :: llr
    !>array of minimal second dimension count(mask(nlr)), if present
    integer, dimension(lr_size,*), intent(out) :: dest_arr
    !> array of logical telling the locregs to mask
    logical, dimension(nlr), intent(in), optional :: mask
    !local variables
    integer :: ilr,iilr

    iilr=0
    do ilr=1,nlr
       iilr=iilr+1
       if (present(mask)) then
          if (.not. mask(ilr)) then
             iilr=iilr-1
             cycle
          end if
       end if
       call locreg_encode(llr(ilr),lr_size,dest_arr(1,iilr))
    end do

  end subroutine locregs_encode


  !> Decode a locreg from a src array
  !! warning, the status of the pointers of lr is nullified after communication
  subroutine locreg_decode(src,lr_size,lr)
    implicit none
    integer, intent(in) :: lr_size
    type(locreg_descriptors), target, intent(out) :: lr
    integer, dimension(lr_size), intent(in) :: src

    lr=transfer(src,lr)
    call nullify_lr_pointers(lr)
    lr%bit%mesh => lr%mesh
  end subroutine locreg_decode


  !> Decode a group of locregs given a src array
  subroutine locreg_full_decode(src,lr_size,lr_full_size,lr,bounds)
    use compression
    implicit none
    integer, intent(in) :: lr_full_size,lr_size
    type(locreg_descriptors), intent(out) :: lr
    integer, dimension(lr_full_size), intent(in) :: src
    logical, intent(in), optional :: bounds
    !local variables
    logical :: bounds_
    integer, dimension(:), pointer :: buffer

    call locreg_decode(src,lr_size,lr)
    if (lr_size == lr_full_size) return
    buffer => f_subptr(src,from=lr_size+1,size=lr_full_size-lr_size)
    call wfd_keys_from_buffer(lr%wfd,buffer)

    if (f_get_option(default=.true.,opt=bounds)) call ensure_locreg_bounds(lr)

  end subroutine locreg_full_decode

  subroutine locregs_decode(src_arr,lr_size,nlr,llr,ipiv)
    implicit none
    integer, intent(in) :: lr_size,nlr
    type(locreg_descriptors), dimension(nlr), intent(inout) :: llr
    integer, dimension(lr_size,nlr), intent(in) :: src_arr
    integer, dimension(nlr), optional  :: ipiv !<array expressing the order of the lrs in the src_arr.
                                               !!When its values are put to zero the update is not performed
    !local variables
    integer :: ilr,iilr

    do ilr=1,nlr
       iilr=ilr
       if (present(ipiv)) iilr=ipiv(ilr)
       !only decode locregs which were not present already
       if (iilr /= 0) call locreg_decode(src_arr(1,ilr),lr_size,llr(iilr))
    end do

  end subroutine locregs_decode


  !> Locreg communication
  subroutine communicate_locreg_descriptors_basics(iproc, nproc, nlr, rootarr, llr)
    implicit none

    ! Calling arguments
    integer,intent(in) :: iproc, nproc, nlr
    integer,dimension(nlr),intent(in) :: rootarr
    type(locreg_descriptors), dimension(nlr), intent(inout) :: llr

    ! Local variables
    character(len=*),parameter :: subname='communicate_locreg_descriptors_basics'
    integer :: ilr,iilr,jlr,jproc,lr_size,nlrp
    !type(locreg_descriptors) :: lr_tmp
    logical, dimension(:), allocatable :: mask
    integer, dimension(:), allocatable :: ipiv,recvcounts
    integer, dimension(:,:), allocatable :: encoded_send,encoded_recv

    call f_routine(id=subname)

    !first encode and communicate
    mask=f_malloc(nlr,id='mask')
    recvcounts=f_malloc0(0.to.nproc-1,id='recvcounts')

    do ilr=1,nlr
       !count order the locregs per process
       jproc=rootarr(ilr)
       recvcounts(jproc)=recvcounts(jproc)+1
       !mask the number of locreg that are associated to the
       !present mpi process
       mask(ilr) = iproc == jproc
    end do

    nlrp=recvcounts(iproc)

    lr_size=get_locreg_encode_size()
    encoded_send=f_malloc([lr_size,nlrp],id='encoded_send')

!!$    iilr=0
!!$    do ilr=1,nlr
!!$       if (mask(ilr)) then
!!$          iilr=iilr+1
!!$          print *,'ilr,iproc',iproc,ilr,llr(ilr)%wfd%nseg_c
!!$          call locreg_encode(llr(ilr),lr_size,encoded_send(1,iilr))
!!$          call locreg_decode(encoded_send(1,iilr),lr_size,lr_tmp)
!!$          print *,'ilr,iproc,after',iproc,ilr,lr_tmp%wfd%nseg_c
!!$       end if
!!$    end do

    call locregs_encode(llr,nlr,lr_size,encoded_send,mask)

    encoded_recv=f_malloc([lr_size,nlr],id='encoded_recv')

    recvcounts=lr_size*recvcounts

    call fmpi_allgather(sendbuf=encoded_send,&
         recvbuf=encoded_recv,&
         recvcounts=recvcounts,comm=bigdft_mpi%mpi_comm)

    call f_free(mask)
    call f_free(recvcounts)
    call f_free(encoded_send)

    !then decode
    ipiv=f_malloc(nlr,id='ipiv')
    iilr=0
    do jproc=0,nproc-1
       do ilr=1,nlr
          jlr=ilr
          if (jproc == iproc) jlr=0 !do not decode locregs already present on the process
          if (rootarr(ilr) == jproc) then
             iilr=iilr+1
             ipiv(iilr)=jlr
          end if
       end do
    end do
    call locregs_decode(encoded_recv,lr_size,nlr,llr,ipiv)
    call f_free(ipiv)
    call f_free(encoded_recv)

    call f_release_routine()

  end subroutine communicate_locreg_descriptors_basics

  !put the localization region lr in the lrs_ptr array
  subroutine store_lr(lr_storage,ilr,lr)
    implicit none
    type(locreg_storage), intent(inout) :: lr_storage
    integer, intent(in) :: ilr
    type(locreg_descriptors), intent(in), target :: lr

    lr_storage%lrs_ptr(ilr)%lr => lr
    lr_storage%lrs_ptr(ilr)%owner = .false. ! Caller is responsible to free lr later
  end subroutine store_lr

  !take ownership of the localization region lr in the lrs_ptr array
  subroutine steal_lr(lr_storage,ilr,lr)
    implicit none
    type(locreg_storage), intent(inout) :: lr_storage
    integer, intent(in) :: ilr
    type(locreg_descriptors), intent(in) :: lr

    allocate(lr_storage%lrs_ptr(ilr)%lr)
    lr_storage%lrs_ptr(ilr)%lr = lr
    !call copy_locreg_descriptors(lr,lr_storage%lrs_ptr(ilr)%lr)
    lr_storage%lrs_ptr(ilr)%owner = .true. ! Caller should not touch lr anymore
  end subroutine steal_lr

  pure function lr_is_stored(lr_storage,ilr) result(ok)
    implicit none
    type(locreg_storage), intent(in) :: lr_storage
    integer, intent(in) :: ilr
    logical :: ok
    ok=associated(lr_storage%lrs_ptr(ilr)%lr) .and. .not. lr_storage%lrs_ptr(ilr)%owner
  end function lr_is_stored

  pure function lr_is_stolen(lr_storage,ilr) result(ok)
    implicit none
    type(locreg_storage), intent(in) :: lr_storage
    integer, intent(in) :: ilr
    logical :: ok
    ok=associated(lr_storage%lrs_ptr(ilr)%lr) .and. lr_storage%lrs_ptr(ilr)%owner
  end function lr_is_stolen


  !> Communicate the total locreg quantities given in a pointer of localisation regions
  !! WARNING: we assume that the association of the storage is performed in a mutually exclusive way, that if a locreg is associated on one proc it is not in the others.
  subroutine gather_locreg_storage(lr_storage)
    implicit none
    type(locreg_storage), intent(inout) :: lr_storage
    !local variables
    integer :: ilr,iilr,nproc,encoding_buffer_idx,nlr
    integer :: lr_full_size,nlrp,encoding_buffer_size
    integer(f_long) :: full_encoding_buffer_size
    type(fmpi_win) :: win_counts
    integer, dimension(:), allocatable :: encoding_buffer
    integer, dimension(:,:), pointer :: lr_sizes


    nlr=size(lr_storage%lrs_ptr)
    lr_storage%lr_full_sizes=f_malloc0_ptr([2,nlr],id='lr_full_sizes')
    encoding_buffer_size=0
    iilr=0
    do ilr=1,nlr
       if ( .not. associated(lr_storage%lrs_ptr(ilr)%lr)) cycle
       iilr=iilr+1
       lr_full_size=get_locreg_full_encode_size(lr_storage%lrs_ptr(ilr)%lr)
       lr_storage%lr_full_sizes(SIZE_,iilr)=lr_full_size
       lr_storage%lr_full_sizes(LR_,iilr)=ilr
       encoding_buffer_size=encoding_buffer_size+lr_full_size
    end do

    if (bigdft_mpi%nproc > 1) then
       lr_sizes=f_malloc_ptr([2,iilr],id='lr_sizes')
       if (iilr > 0) &
            & call f_memcpy(n=2*iilr,src=lr_storage%lr_full_sizes,dest=lr_sizes(1,1))
       call fmpi_allgather(sendbuf=lr_sizes,recvbuf=lr_storage%lr_full_sizes,&
            comm=bigdft_mpi%mpi_comm,win=win_counts)
    else
       lr_sizes => lr_storage%lr_full_sizes
    end if

    encoding_buffer=f_malloc(encoding_buffer_size,id='encoding_buffer')
    !redo the loop to fill the encoding buffer
    encoding_buffer_idx=1
    iilr=0
    do ilr=1,nlr
       if ( .not. associated(lr_storage%lrs_ptr(ilr)%lr)) cycle
       iilr=iilr+1
       lr_full_size=lr_sizes(SIZE_,iilr) !storage%lr_full_sizes(SIZE_,ilr)
       call locreg_full_encode(lr_storage%lrs_ptr(ilr)%lr,lr_storage%lr_size, &
            & lr_full_size,encoding_buffer(encoding_buffer_idx))
       encoding_buffer_idx=encoding_buffer_idx+lr_full_size
    end do

    !close the previous communication window
    if (bigdft_mpi%nproc > 1) then
       call fmpi_win_shut(win_counts)
       call f_free_ptr(lr_sizes)
    else
       nullify(lr_sizes)
    end if

    !allocate the full sized array
    full_encoding_buffer_size=0
    do ilr=1,nlr
       full_encoding_buffer_size=full_encoding_buffer_size+lr_storage%lr_full_sizes(SIZE_,ilr)
    end do
    lr_storage%encode_buffer=f_malloc_ptr(full_encoding_buffer_size,id='lr_storage%encode_buffer')

    !gather the array in the full encoding buffer
    if (bigdft_mpi%nproc > 1) then
       call fmpi_allgather(sendbuf=encoding_buffer,recvbuf=lr_storage%encode_buffer,comm=bigdft_mpi%mpi_comm)
    else
       call f_memcpy(n = encoding_buffer_size, src = encoding_buffer, dest = lr_storage%encode_buffer(1))
    end if

    call f_free(encoding_buffer)

  end subroutine gather_locreg_storage

  subroutine extract_lr(lr_storage,ilr,lr,bounds)
    implicit none
    type(locreg_storage), intent(in) :: lr_storage
    integer, intent(in) :: ilr
    type(locreg_descriptors), intent(out) :: lr
    logical, intent(in), optional :: bounds
    !local variables
    integer :: iilr,nlr
    integer(f_long) :: encoding_buffer_idx
    integer, dimension(:), pointer :: src_buf

    nlr=size(lr_storage%lr_full_sizes,dim=2)
    !assume that lr_storage is ready for extraction
    encoding_buffer_idx=1
    do iilr=1,nlr
       if (lr_storage%lr_full_sizes(LR_,iilr) == ilr) exit
       encoding_buffer_idx=encoding_buffer_idx+lr_storage%lr_full_sizes(SIZE_,iilr)
    end do
    if (iilr == nlr+1) call f_err_throw('The locreg has not been found in the lr_storage',&
         err_name='BIGDFT_RUNTIME_ERROR')

    !we should generalize the API for the from optional variable (it should accept also f_long)
    src_buf => f_subptr(lr_storage%encode_buffer,&
         from=int(encoding_buffer_idx),size=lr_storage%lr_full_sizes(SIZE_,iilr))

    call locreg_full_decode(src_buf,&
         lr_storage%lr_size,lr_storage%lr_full_sizes(SIZE_,iilr),lr,bounds)

  end subroutine extract_lr

  !> Methods for copying the structures, can be needed to avoid recalculating them
  !! should be better by defining a f_malloc inheriting the shapes and the structure from other array
  !! of the type dest=f_malloc(src=source,id='dest')
  subroutine copy_locreg_descriptors(glrin, glrout)
    use bounds
    implicit none
    ! Calling arguments
    type(locreg_descriptors), intent(in) :: glrin !<input locreg. Unchanged on exit.
    type(locreg_descriptors), intent(out), target :: glrout !<output locreg. Must be freed on input.

    !we should here use encode and decode for safety

!!$    glrout%geocode = glrin%geocode
    glrout%hybrid_on = glrin%hybrid_on
    glrout%ns1 = glrin%ns1
    glrout%ns2 = glrin%ns2
    glrout%ns3 = glrin%ns3
    glrout%nsi1 = glrin%nsi1
    glrout%nsi2 = glrin%nsi2
    glrout%nsi3 = glrin%nsi3
    glrout%Localnorb = glrin%Localnorb
    glrout%locrad=glrin%locrad
    glrout%locrad_kernel=glrin%locrad_kernel
    glrout%locrad_mult=glrin%locrad_mult
    glrout%locregCenter=glrin%locregCenter
    glrout%outofzone= glrin%outofzone

    call copy_grid_dimensions(glrin%d, glrout%d)
    call copy_wavefunctions_descriptors(glrin%wfd, glrout%wfd)
    call copy_convolutions_bounds(glrin%bounds, glrout%bounds)

    glrout%mesh=glrin%mesh
    glrout%mesh_fine=glrin%mesh_fine
    glrout%mesh_coarse=glrin%mesh_coarse
    glrout%bit=glrin%bit
    glrout%bit%mesh => glrout%mesh
  end subroutine copy_locreg_descriptors


  pure subroutine copy_grid_dimensions(din, dout)
    implicit none
    ! Calling arguments
    type(grid_dimensions),intent(in):: din
    type(grid_dimensions),intent(out):: dout

    dout%n1 = din%n1
    dout%n2 = din%n2
    dout%n3 = din%n3
    dout%nfl1 = din%nfl1
    dout%nfu1 = din%nfu1
    dout%nfl2 = din%nfl2
    dout%nfu2 = din%nfu2
    dout%nfl3 = din%nfl3
    dout%nfu3 = din%nfu3
    dout%n1i = din%n1i
    dout%n2i = din%n2i
    dout%n3i = din%n3i

  end subroutine copy_grid_dimensions


  subroutine ensure_locreg_bounds(lr)
    use bounds, only: locreg_bounds
    implicit none
    type(locreg_descriptors), intent(inout) :: lr

    !take this as exemple of already associated bounds
    if (associated(lr%bounds%kb%ibyz_c)) return

    call locreg_bounds(lr%d%n1,lr%d%n2,lr%d%n3,&
         lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,&
         lr%d%nfl3,lr%d%nfu3,lr%wfd,lr%bounds)

  end subroutine ensure_locreg_bounds


  !> Almost degenerate with get_number_of_overlap_region
  !! should merge the two... prefering this one since argument list is better
  subroutine check_overlap_cubic_periodic(Glr,Ilr,Jlr,isoverlap)
    use module_base
    use bounds, only: check_whether_bounds_overlap
    implicit none
    type(locreg_descriptors), intent(in) :: Glr
    type(locreg_descriptors), intent(in) :: Ilr
    type(locreg_descriptors), intent(in) :: Jlr
    logical, intent(out) :: isoverlap
    !Local variables
    integer :: is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3
    logical :: overlap1, overlap2, overlap3
  !!  integer :: azones,bzones,ii,izones,jzones !, i_stat, i_all
  !!  logical :: go1, go2, go3
  !!  integer,dimension(3,8) :: astart,bstart,aend,bend

  !!  azones = 1
  !!  bzones = 1
  !!! Calculate the number of regions to cut alr and blr
  !!  do ii=1,3
  !!     if(Ilr%outofzone(ii) > 0) azones = azones * 2
  !!     if(Jlr%outofzone(ii) > 0) bzones = bzones * 2
  !!  end do
  !!
  !!!FRACTURE THE FIRST LOCALIZATION REGION
  !!  call fracture_periodic_zone(azones,Glr,Ilr,Ilr%outofzone,astart,aend)
  !!
  !!!FRACTURE SECOND LOCREG
  !!  call fracture_periodic_zone(bzones,Glr,Jlr,Jlr%outofzone,bstart,bend)
  !!
  !!! Now check if they overlap
  !!  isoverlap = .false.
  !!  loop_izones: do izones=1,azones
  !!    do jzones=1,bzones
  !!      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones))
  !!      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones))
  !!      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones))
  !!      if(go1 .and. go2 .and. go3) then
  !!        isoverlap = .true.
  !!        exit loop_izones
  !!      end if
  !!    end do
  !!  end do loop_izones


    !@ NEW VERSION #########################################
    ! Shift all the indices into the periodic cell. This can result is starting
    ! indices being larger than ending indices
    isoverlap = .false.
    is3 = modulo(ilr%ns3,glr%d%n3+1)
    ie3 = modulo(ilr%ns3+ilr%d%n3,glr%d%n3+1)
    js3 = modulo(jlr%ns3,glr%d%n3+1)
    je3 = modulo(jlr%ns3+jlr%d%n3,glr%d%n3+1)
    overlap3 = check_whether_bounds_overlap(is3, ie3, js3, je3)
    if (overlap3) then
        is2 = modulo(ilr%ns2,glr%d%n2+1)
        ie2 = modulo(ilr%ns2+ilr%d%n2,glr%d%n2+1)
        js2 = modulo(jlr%ns2,glr%d%n2+1)
        je2 = modulo(jlr%ns2+jlr%d%n2,glr%d%n2+1)
        overlap2 = check_whether_bounds_overlap(is2, ie2, js2, je2)
        if (overlap2) then
            is1 = modulo(ilr%ns1,glr%d%n1+1)
            ie1 = modulo(ilr%ns1+ilr%d%n1,glr%d%n1+1)
            js1 = modulo(jlr%ns1,glr%d%n1+1)
            je1 = modulo(jlr%ns1+jlr%d%n1,glr%d%n1+1)
            overlap1 = check_whether_bounds_overlap(is1, ie1, js1, je1)
            if (overlap1) then
                ! If we are here, all three overlaps are true
                isoverlap = .true.
            end if
        end if
    end if

    !!if (overlap1 .and. overlap2 .and. overlap3) then
    !!    isoverlap = .true.
    !!else
    !!    isoverlap = .false.
    !!end if

    !@ END NEW VERSION #####################################

    !!!debug
    !!isoverlap=.true.

  end subroutine check_overlap_cubic_periodic

    subroutine check_overlap(Llr_i, Llr_j, Glr, overlap)
      implicit none

      ! Calling arguments
      type(locreg_descriptors),intent(in) :: Llr_i, Llr_j, Glr
      logical, intent(out) :: overlap

      ! Local variables
      integer :: onseg

      call check_overlap_cubic_periodic(Glr,Llr_i,Llr_j,overlap)
      if(overlap) then
         call check_overlap_from_descriptors_periodic(Llr_i%wfd%nseg_c, Llr_j%wfd%nseg_c,&
              Llr_i%wfd%keyglob, Llr_j%wfd%keyglob, overlap, onseg)
      end if
    end subroutine check_overlap

    ! check if Llrs overlap from there descriptors
    ! The periodicity is hidden in the fact that we are using the keyglobs
    ! which are correctly defined.
    subroutine check_overlap_from_descriptors_periodic(nseg_i, nseg_j, keyg_i, keyg_j,  &
         isoverlap, onseg)
      implicit none
      ! Calling arguments
      integer :: nseg_i, nseg_j
      integer,dimension(2,nseg_i),intent(in) :: keyg_i
      integer,dimension(2,nseg_j),intent(in) :: keyg_j
      logical,intent(out) :: isoverlap
      integer, intent(out) :: onseg
      ! Local variables
      integer :: iseg, jseg, istart, jstart, kstartg
      integer :: iend, jend, kendg, nseg_k


      ! Initialize some counters
      iseg=1
      jseg=1
      nseg_k=0
      isoverlap = .false.
      onseg = 0  ! in case they don't overlap
      ! Check whether all segments of both localization regions have been processed.
      if ((iseg>=nseg_i .and. jseg>=nseg_j) .or. nseg_i==0 .or. nseg_j==0) return

      segment_loop: do

         ! Starting point already in global coordinates
         istart=keyg_i(1,iseg)
         jstart=keyg_j(1,jseg)

         ! Ending point already in global coordinates
         iend=keyg_i(2,iseg)
         jend=keyg_j(2,jseg)
         ! Determine starting and ending point of the common segment in global coordinates.
         kstartg=max(istart,jstart)
         kendg=min(iend,jend)

         ! Check whether this common segment has a non-zero length
         if(kendg-kstartg+1>0) then
            isoverlap = .true.
            nseg_k=nseg_k+1
         end if

         ! Check whether all segments of both localization regions have been processed.
         if(iseg>=nseg_i .and. jseg>=nseg_j) exit segment_loop

         ! Increase the segment index
         if((iend<=jend .and. iseg<nseg_i) .or. jseg==nseg_j) then
            iseg=iseg+1
         else if(jseg<nseg_j) then
            jseg=jseg+1
         end if

      end do segment_loop

      if(isoverlap) then
         onseg = nseg_k
      end if

    end subroutine check_overlap_from_descriptors_periodic

    pure function grid_init(peri,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,&
         ns1,ns2,ns3) result(g)
      implicit none
      integer, intent(in) :: n1,n2,n3
      integer, intent(in) :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
      integer, intent(in) :: ns1,ns2,ns3
      logical, dimension(3), intent(in) :: peri !<periodic dimensions
      type(grid_dimensions) :: g
      !local variables
      integer, parameter :: ISF_GROW_BUFFER=31

      g%n1=n1-ns1
      g%n2=n2-ns2
      g%n3=n3-ns3

      !dimensions of the fine grid inside the localisation region
      g%nfl1=max(ns1,nfl1)-ns1
      g%nfl2=max(ns2,nfl2)-ns2
      g%nfl3=max(ns3,nfl3)-ns3

      !NOTE: This will not work with symmetries (must change it)
      g%nfu1=min(n1,nfu1)-ns1
      g%nfu2=min(n2,nfu2)-ns2
      g%nfu3=min(n3,nfu3)-ns3

      if (peri(1)) then
         g%n1i=2*g%n1+2
      else
         g%n1i=2*g%n1+ISF_GROW_BUFFER
      end if
      if (peri(2)) then
         g%n2i=2*g%n2+2
      else
         g%n2i=2*g%n2+ISF_GROW_BUFFER
      end if
      if (peri(3)) then
         g%n3i=2*g%n3+2
      else
         g%n3i=2*g%n3+ISF_GROW_BUFFER
      end if

    end function grid_init

    !> Create the localisation region information for cubic code
    subroutine init_lr(lr,geocode,hgridsh,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,&
         hybrid_flag,isx,isy,isz,global_geocode,wfd,bnds)
      use compression
      use bounds
      use box
      implicit none
      logical, intent(in) :: hybrid_flag
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
      real(gp), dimension(3), intent(in) :: hgridsh
      character(len=1), intent(in), optional :: global_geocode
      !>have to be present with global_geocode
      integer, intent(in), optional :: isx,isy,isz
      type(wavefunctions_descriptors), intent(in), optional :: wfd
      type(convolutions_bounds), intent(in), optional :: bnds
      type(locreg_descriptors), intent(inout) :: lr
      !local variables
      integer, parameter :: S0_GROW_BUFFER=14
      integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
      integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
      logical, dimension(3) :: peri,peri_glob
      integer, dimension(3) :: ndims
      real(gp), dimension(3) :: oxyz,hgrids

      call f_routine(id='init_lr')

!!$      lr%geocode=geocode
      if (present(isx)) then
         lr%ns1=isx
      else
         lr%ns1=0
      end if
      if (present(isy)) then
         lr%ns2=isy
      else
         lr%ns2=0
      end if
      if (present(isz)) then
         lr%ns3=isz
      else
         lr%ns3=0
      end if

!!$      peri(1)=geocode /= 'F'
!!$      peri(2)=geocode == 'P'
!!$      peri(3)=geocode /= 'F'

      peri=bc_periodic_dims(geocode_to_bc(geocode))

      lr%d=grid_init(peri,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,&
         lr%ns1,lr%ns2,lr%ns3)
      ndims(1)=lr%d%n1i
      ndims(2)=lr%d%n2i
      ndims(3)=lr%d%n3i

      !this is the mesh in real space (ISF basis)
      lr%mesh=cell_new(geocode,ndims,hgridsh)

      call ext_buffers_coarse(peri(1),Lnbl1)
      call ext_buffers_coarse(peri(2),Lnbl2)
      call ext_buffers_coarse(peri(3),Lnbl3)

      ndims(1)=2*lr%d%n1+2+2*Lnbl1
      ndims(2)=2*lr%d%n2+2+2*Lnbl2
      ndims(3)=2*lr%d%n3+2+2*Lnbl3
      !this is the mesh in the fine scaling function
      lr%mesh_fine=cell_new(geocode,ndims,hgridsh)

      ndims(1)=lr%d%n1+1
      ndims(2)=lr%d%n2+1
      ndims(3)=lr%d%n3+1
      hgrids=2.0_gp*hgridsh
      lr%mesh_coarse=cell_new(geocode,ndims,hgrids) !we should write the number of points here

      Gnbl1=0
      Gnbl2=0
      Gnbl3=0
      if (present(global_geocode)) then
         peri_glob(1)=(global_geocode /= 'F')
         peri_glob(2)=(global_geocode == 'P')
         peri_glob(3)=(global_geocode /= 'F')
         call ext_buffers(peri_glob(1),Gnbl1,Gnbr1)
         call ext_buffers(peri_glob(2),Gnbl2,Gnbr2)
         call ext_buffers(peri_glob(3),Gnbl3,Gnbr3)
         call ext_buffers(peri(1),Lnbl1,Lnbr1)
         call ext_buffers(peri(2),Lnbl2,Lnbr2)
         call ext_buffers(peri(3),Lnbl3,Lnbr3)
      end if
      lr%nsi1=0
      lr%nsi2=0
      lr%nsi3=0
      if (present(isx)) lr%nsi1= 2 * lr%ns1 - (Lnbl1 - Gnbl1)
      if (present(isy)) lr%nsi2= 2 * lr%ns2 - (Lnbl2 - Gnbl2)
      if (present(isz)) lr%nsi3= 2 * lr%ns3 - (Lnbl3 - Gnbl3)

      lr%hybrid_on = hybrid_flag
      lr%hybrid_on=lr%hybrid_on .and. (nfu1-nfl1+S0_GROW_BUFFER < n1+1)
      lr%hybrid_on=lr%hybrid_on .and. (nfu2-nfl2+S0_GROW_BUFFER < n2+1)
      lr%hybrid_on=lr%hybrid_on .and. (nfu3-nfl3+S0_GROW_BUFFER < n3+1)

      if (present(wfd)) lr%wfd=wfd !it just associates the pointers

      !this is a point where the geocode is stull used
      if (geocode == 'F' .and. present(bnds)) lr%bounds=bnds

      !here we have to put the modifications of the origin for the
      !iterator of the lr. get_isf_offset should be used as
      !soon as global_geocode is replaced
      oxyz=locreg_mesh_origin(lr%mesh)
      lr%bit=box_iter(lr%mesh,origin=oxyz)

      call f_release_routine()

    END SUBROUTINE init_lr

    subroutine correct_dimensions(is,ie,ns,n,ln,outofzone,correct,periodic)
      implicit none
      logical, intent(in) :: correct
      integer, intent(in) :: ns,n,ln
      integer, intent(inout) :: outofzone
      logical, intent(inout) :: periodic
      integer, intent(inout) :: is,ie

      if (ie - is >= n) then
         is=ns
         ie=ns + n
         periodic = .true.
      else
         if (correct) then
            is=modulo(is,n+1) + ns
            ie= ln + is
         end if
         if (ie > ns+n) then
            outofzone=modulo(ie,n+1)
         end if
      end if
    end subroutine correct_dimensions


    subroutine correct_lr_extremes(lr,Glr,geocode,correct,nbox_lr,nbox)
      use box, only: cell_geocode,cell_periodic_dims
      implicit none
      logical, intent(in) :: correct
      type(locreg_descriptors), intent(inout) :: lr
      type(locreg_descriptors), intent(in) :: Glr
      character(len=1), intent(out) :: geocode
      integer, dimension(2,3), intent(out) :: nbox_lr
      integer, dimension(2,3), intent(in), optional :: nbox
      !local variables
      logical :: xperiodic,yperiodic,zperiodic
      integer :: isx,iex,isy,iey,isz,iez
      integer :: ln1,ln2,ln3
      logical, dimension(3) :: peri


      !initialize out of zone
      lr%outofzone (:) = 0

      if (present(nbox)) then
         ! Localization regions should have free boundary conditions by default
         isx=nbox(1,1)
         iex=nbox(2,1)
         isy=nbox(1,2)
         iey=nbox(2,2)
         isz=nbox(1,3)
         iez=nbox(2,3)
      else !otherwise get box from previously initialized locreg
         isx=lr%ns1
         iex=lr%ns1+lr%d%n1
         isy=lr%ns2
         iey=lr%ns2+lr%d%n2
         isz=lr%ns3
         iez=lr%ns3+lr%d%n3
      end if
      ln1 = iex-isx
      ln2 = iey-isy
      ln3 = iez-isz

      geocode='F'

      xperiodic = .false.
      yperiodic = .false.
      zperiodic = .false.

      peri=cell_periodic_dims(Glr%mesh)
      if (peri(1)) then
         call correct_dimensions(isx,iex,Glr%ns1,Glr%d%n1,ln1,lr%outofzone(1),correct,xperiodic)
      else
         isx=max(isx,Glr%ns1)
         iex=min(iex,Glr%ns1+Glr%d%n1)
         lr%outofzone(1) = 0
      end if
      if (peri(2)) then
         call correct_dimensions(isy,iey,Glr%ns2,Glr%d%n2,ln2,lr%outofzone(2),correct,yperiodic)
      else
         isy=max(isy,Glr%ns2)
         iey=min(iey,Glr%ns2+Glr%d%n2)
         lr%outofzone(2) = 0
      end if
      if (peri(3)) then
         call correct_dimensions(isz,iez,Glr%ns3,Glr%d%n3,ln3,lr%outofzone(3),correct,zperiodic)
      else
         isz=max(isz,Glr%ns3)
         iez=min(iez,Glr%ns3+Glr%d%n3)
         lr%outofzone(3) = 0
      end if
      if (zperiodic) geocode = 'W'
      if (xperiodic .and. zperiodic) geocode = 'S'
      if (xperiodic .and. yperiodic .and. zperiodic) geocode = 'P'

!!$      !assign the starting/ending points and outofzone for the different
!!$      ! geometries
!!$      !!!select case(Glr%geocode)
!!$      select case(cell_geocode(Glr%mesh))
!!$      case('F')
!!$         isx=max(isx,Glr%ns1)
!!$         isy=max(isy,Glr%ns2)
!!$         isz=max(isz,Glr%ns3)
!!$
!!$         iex=min(iex,Glr%ns1+Glr%d%n1)
!!$         iey=min(iey,Glr%ns2+Glr%d%n2)
!!$         iez=min(iez,Glr%ns3+Glr%d%n3)
!!$      case('S')
!!$         ! Get starting and ending for x direction
!!$         if (iex - isx >= Glr%d%n1) then
!!$            isx=Glr%ns1
!!$            iex=Glr%ns1 + Glr%d%n1
!!$            xperiodic = .true.
!!$         else
!!$            if (correct) then
!!$               isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
!!$               iex= ln1 + isx
!!$            end if
!!$            if (iex > Glr%ns1+Glr%d%n1) then
!!$               lr%outofzone(1)=modulo(iex,Glr%d%n1+1)
!!$            end if
!!$         end if
!!$
!!$         ! Get starting and ending for y direction (perpendicular to surface)
!!$         isy=max(isy,Glr%ns2)
!!$         iey=min(iey,Glr%ns2 + Glr%d%n2)
!!$         lr%outofzone(2) = 0
!!$
!!$         !Get starting and ending for z direction
!!$         if (iez - isz >= Glr%d%n3) then
!!$            isz=Glr%ns3
!!$            iez=Glr%ns3 + Glr%d%n3
!!$            zperiodic = .true.
!!$         else
!!$            if (correct) then
!!$               isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
!!$               iez= ln3 + isz
!!$            end if
!!$            if (iez > Glr%ns3+Glr%d%n3) then
!!$               lr%outofzone(3)=modulo(iez,Glr%d%n3+1)
!!$            end if
!!$         end if
!!$         if(xperiodic .and. zperiodic) then
!!$            geocode = 'S'
!!$         end if
!!$
!!$      case('P')
!!$         ! Get starting and ending for x direction
!!$         if (iex - isx >= Glr%d%n1) then
!!$            isx=Glr%ns1
!!$            iex=Glr%ns1 + Glr%d%n1
!!$            xperiodic = .true.
!!$         else
!!$            if (correct) then
!!$               isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
!!$               iex= ln1 + isx
!!$            end if
!!$            if (iex > Glr%ns1+Glr%d%n1) then
!!$               lr%outofzone(1)=modulo(iex,Glr%d%n1+1)
!!$            end if
!!$         end if
!!$
!!$         ! Get starting and ending for y direction (perpendicular to surface)
!!$         if (iey - isy >= Glr%d%n2) then
!!$            isy=Glr%ns2
!!$            iey=Glr%ns2 + Glr%d%n2
!!$            yperiodic = .true.
!!$         else
!!$            if (correct) then
!!$               isy=modulo(isy,Glr%d%n2+1) + Glr%ns2
!!$               iey= ln2 + isy
!!$            end if
!!$            if (iey > Glr%ns2+Glr%d%n2) then
!!$               lr%outofzone(2)=modulo(iey,Glr%d%n2+1)
!!$            end if
!!$         end if
!!$
!!$         !Get starting and ending for z direction
!!$         if (iez - isz >= Glr%d%n3) then
!!$            isz=Glr%ns3
!!$            iez=Glr%ns3 + Glr%d%n3
!!$            zperiodic = .true.
!!$         else
!!$            if (correct) then
!!$               isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
!!$               iez= ln3 + isz
!!$            end if
!!$            if (iez > Glr%ns3+Glr%d%n3) then
!!$               lr%outofzone(3)=modulo(iez,Glr%d%n3+1)
!!$            end if
!!$         end if
!!$         if(xperiodic .and. yperiodic .and. zperiodic ) then
!!$            geocode = 'P'
!!$         end if
!!$      case('W')
!!$         call f_err_throw("Wires bc has to be implemented here", &
!!$              err_name='BIGDFT_RUNTIME_ERROR')
!!$      end select
      ! Make sure that the localization regions are not periodic
      if (xperiodic .or. yperiodic .or. zperiodic) then
         call f_err_throw('The localization region '//&
              ' is supposed to be fully free BC.'//&
              ' Reduce the localization radii or use the cubic version',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end if

      nbox_lr(1,1)=isx
      nbox_lr(1,2)=isy
      nbox_lr(1,3)=isz

      nbox_lr(2,1)=iex
      nbox_lr(2,2)=iey
      nbox_lr(2,3)=iez

    end subroutine correct_lr_extremes


    !initialize the parts of lr that can be initialized from a load from a
    !traditional restart. The output files have to be redefined such that
    !this routine would not need to be called
    subroutine reset_lr(lr,geocode,hgrids,nbox,global_geocode)
      implicit none
      character(len=1), intent(in) :: geocode,global_geocode
      real(gp), dimension(3), intent(in) :: hgrids
      integer, dimension(2,3), intent(in) :: nbox !fine box of the environmental region
      type(locreg_descriptors), intent(inout) :: lr

      call init_lr(lr,geocode,0.5_gp*hgrids,&
           lr%ns1+lr%d%n1,lr%ns2+lr%d%n2,lr%ns3+lr%d%n3,&
           nbox(1,1),nbox(1,2),nbox(1,3),&
           nbox(2,1),nbox(2,2),nbox(2,3),&
           .false.,lr%ns1,lr%ns2,lr%ns3,global_geocode)

    end subroutine reset_lr

    !> initalize the box-related components of the localization regions
    subroutine lr_box(lr,Glr,hgrids,nbox,correction)
      use bounds, only: ext_buffers
      use box, only: cell_geocode
      implicit none
      !> Sub-box to iterate over the points (ex. around atoms)
      !! start and end points for each direction
      real(gp), dimension(3), intent(in) :: hgrids
      type(locreg_descriptors), intent(in) :: Glr
      type(locreg_descriptors), intent(inout) :: lr
      integer, dimension(2,3), intent(in), optional :: nbox
      logical, intent(in), optional :: correction
      !local variables
      character(len=1) :: geocode
      logical :: correct
      integer, dimension(2,3) :: nbox_lr
      real(gp), dimension(3) :: hgridsh

      call f_routine(id='lr_box')

      correct=.false.
      if (present(correction)) correct=correction

      call correct_lr_extremes(lr,Glr,geocode,correct,nbox_lr,nbox)

      hgridsh=0.5_gp*hgrids
      call init_lr(lr,geocode,hgridsh,nbox_lr(2,1),nbox_lr(2,2),nbox_lr(2,3),&
           Glr%d%nfl1,Glr%d%nfl2,Glr%d%nfl3,&
           Glr%d%nfu1,Glr%d%nfu2,Glr%d%nfu3,&
!!$           .false.,nbox_lr(1,1),nbox_lr(1,2),nbox_lr(1,3),Glr%geocode)
           .false.,nbox_lr(1,1),nbox_lr(1,2),nbox_lr(1,3),cell_geocode(Glr%mesh))

      ! Make sure that the extent of the interpolating functions grid for the
      ! locreg is not larger than the that of the global box.
      if (lr%d%n1i>Glr%d%n1i) then
         call f_err_throw('The interpolating functions grid in x dimension for locreg '&
              &//&!trim(yaml_toa(ilr,fmt='(i0)'))//&
              '('//trim(yaml_toa(lr%d%n1i,fmt='(i0)'))//')&
              & is larger than that of the global region('//trim(yaml_toa(Glr%d%n1i,fmt='(i0)'))//').&
              & Reduce the localization radii or use the cubic scaling version',&
              & err_name='BIGDFT_RUNTIME_ERROR')
      end if
      if (lr%d%n2i>Glr%d%n2i) then
         call f_err_throw('The interpolating functions grid in y dimension for locreg '&
              !&//trim(yaml_toa(ilr,fmt='(i0)'))&
              //'('//trim(yaml_toa(lr%d%n2i,fmt='(i0)'))//')&
              & is larger than that of the global region('//trim(yaml_toa(Glr%d%n2i,fmt='(i0)'))//').&
              & Reduce the localization radii or use the cubic scaling version',&
              & err_name='BIGDFT_RUNTIME_ERROR')
      end if
      if (lr%d%n3i>Glr%d%n3i) then
         call f_err_throw('The interpolating functions grid in z dimension for locreg '&
              !&//trim(yaml_toa(ilr,fmt='(i0)'))&
              //'('//trim(yaml_toa(lr%d%n3i,fmt='(i0)'))//')&
              & is larger than that of the global region('//trim(yaml_toa(Glr%d%n3i,fmt='(i0)'))//').&
              & Reduce the localization radii or use the cubic scaling version',&
              & err_name='BIGDFT_RUNTIME_ERROR')
      end if

      call f_release_routine()

    end subroutine lr_box

    !>get the offset of the isf description of the support function
    pure function get_isf_offset(lr,mesh_global) result(ioffset)
      use box, only: cell,cell_periodic_dims
      use bounds, only: isf_box_buffers
      implicit none
      type(locreg_descriptors), intent(in) :: lr
      type(cell), intent(in) :: mesh_global
      integer, dimension(3) :: ioffset
      !local variables
      logical, dimension(3) :: peri_local,peri_global
      integer, dimension(3) :: nli
      !integer :: nl1, nl2, nl3, nr1, nr2, nr3

      !geocode_buffers
      !conditions for periodicity in the three directions
      peri_local=cell_periodic_dims(lr%mesh)
      peri_global=cell_periodic_dims(mesh_global)

      nli=isf_box_buffers(peri_local,peri_global)

      ! offset
      ioffset(1) = lr%nsi1 - nli(1) - 1
      ioffset(2) = lr%nsi2 - nli(2) - 1
      ioffset(3) = lr%nsi3 - nli(3) - 1

    end function get_isf_offset


end module locregs
