!> @file
!! Wrapper for mpi_alltoall flavours
!! Use error handling
!! @author
!!    Copyright (C) 2012-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_alltoall
  use f_enums
  use f_precisions
  use fmpi_types
  use dictionaries, only: f_err_throw
  use f_allreduce
  use f_onesided
  implicit none

  private
  
  integer, parameter :: AUTOMATIC=0
  integer, parameter :: NOT_VARIABLE=10
  integer, parameter :: VARIABLE=11
  integer, parameter :: VARIABLE_ONE_SIDED_GET=12

  type(f_enumerator), public :: AUTOMATIC_ENUM =f_enumerator('AUTOMATIC',AUTOMATIC,null())
  type(f_enumerator), public :: ALLTOALL_ENUM =f_enumerator('ALLTOALL',NOT_VARIABLE,null())
  type(f_enumerator), public :: ALLTOALLV_ENUM =f_enumerator('ALLTOALLV',VARIABLE,null())
  type(f_enumerator), public :: ALLTOALL_GET_ENUM =f_enumerator('ALLTOALL_GET',VARIABLE_ONE_SIDED_GET,null())


  interface fmpi_alltoall
     module procedure mpialltoallw_d11
  end interface fmpi_alltoall
  
!!$  interface mpialltoallv
!!$     module procedure mpialltoallv_int11, mpialltoallv_long11, mpialltoallv_double11
!!$     !module procedure mpialltoallv_double61
!!$  end interface mpialltoallv
!!$
!!$  interface mpiialltoallv
!!$     module procedure mpiialltoallv_double
!!$  end interface mpiialltoallv
!!$
!!$  interface mpi_get_alltoallv
!!$     module procedure mpi_get_alltoallv_i, mpi_get_alltoallv_l, mpi_get_alltoallv_d
!!$  end interface mpi_get_alltoallv

  public :: fmpi_alltoall
  
  contains

    !>choose the correct algorithm to be applied as a function of the sendrecv arrays
    function alltoall_algorithm(nproc,sizets,sizetr,&
         sendcount,sendcounts,sdispls,&
         recvcount,recvcounts,rdispls,algorithm,comm) result(algo)
      use yaml_strings
      implicit none
      integer, intent(in) :: nproc
      integer(f_long), intent(in) :: sizets,sizetr
      integer, intent(in), optional :: sendcount,recvcount
      integer, dimension(0:), intent(in), optional :: sendcounts,sdispls
      integer, dimension(0:), intent(in), optional :: recvcounts,rdispls
      type(f_enumerator), intent(in), optional :: algorithm
      integer, intent(in), optional :: comm
      integer :: algo
      !local variables
      integer, parameter :: max_size = 100000000 !maximal number of elements sent and/or received on a proc using the standard MPI function
      logical :: large
      integer(f_long) :: sizecomms,sizecommr

      !priority to variable version
      if (present(sendcounts)) then
         sizecomms=sum(int(sendcounts,f_long))
      else if (present(sendcount)) then
         sizecomms=int(sendcount,f_long)*nproc
      end if
      sizecommr=sizecomms

      if (present(recvcounts)) then
         sizecommr=sum(int(recvcounts,f_long))
      else if (present(recvcount)) then
         sizecommr=int(recvcount,f_long)*nproc
      end if

      !check the validity of the buffers
      if (sizets < sizecomms .and. sizets /=0) call f_err_throw(&
           'Size of send buffer and of data are not consistent ('&
           //trim(yaml_toa([sizets,sizecomms]))//')',&
           err_name='ERR_MPI_WRAPPERS')

      !check the validity of the buffers
      if (sizetr < sizecommr .and. sizetr /=0) call f_err_throw(&
           'Size of recv buffer and of data are not consistent ('&
           //trim(yaml_toa([sizetr,sizecommr]))//')',&
           err_name='ERR_MPI_WRAPPERS')

      if (present(sendcounts) .and. present(sdispls) .and. present(recvcounts) .and. present(rdispls)) then
         algo=toi(ALLTOALLV_ENUM)
      else if(present(sendcount) .and. present(recvcount)) then
         algo=toi(ALLTOALL_ENUM)
      else
         call f_err_throw('Illegal arguments for alltoall',err_id=ERR_MPI_WRAPPERS)
      end if

      if (present(algorithm)) then
         !if the algorithm has been already provided check if it is possible
         select case(toi(algorithm))
         case(AUTOMATIC)
            !nothing to do here
         case(NOT_VARIABLE)
            !check if the proposition is compatible
            if (algo == toi(ALLTOALLV_ENUM)) &
                 call f_err_throw('Algorithm not compatible with arguments',err_id=ERR_MPI_WRAPPERS)
            return !ok
         case(VARIABLE,VARIABLE_ONE_SIDED_GET)
            !check if the proposition is compatible
            if (algo == toi(ALLTOALL_ENUM)) &
                 call f_err_throw('Algorithm (variable) not compatible with arguments',err_id=ERR_MPI_WRAPPERS)
            if (algorithm == ALLTOALL_GET_ENUM) algo=VARIABLE_ONE_SIDED_GET
            return !ok
         case default
            call f_err_throw('Illegal value for algorithm',err_id=ERR_MPI_WRAPPERS)
         end select
      end if

      if (algo==VARIABLE) then !here we are in the automatic case
         ! Check whether we are having a "big" case. If so, use the hand-made
         ! version using mpi_get, otherwise the standard function.
         !large = (sizets>max_size .or. sizetr>max_size) .and. associated(sdispls_) .and. associated(rdispls_)
         large = (sizecommr>max_size .or. sizecomms>max_size)
         call fmpi_allreduce(large, 1,FMPI_LOR, comm=comm)
         if (large) algo=VARIABLE_ONE_SIDED_GET
      end if
     
    end function alltoall_algorithm

    !recursive
    subroutine mpialltoallw_d11(sendbuf, recvbuf, &
         count, sendcounts, sdispls, recvcounts, rdispls, comm, request, win, algorithm)
      use dynamic_memory
      use yaml_output
      use f_utils, only: f_size
      use dictionaries
      !use iso_c_binding
      implicit none
      real(f_double), dimension(:),intent(in)  :: sendbuf
      real(f_double), dimension(:),intent(out) :: recvbuf
      integer, intent(in), optional :: count
      integer, dimension(0:), intent(in), optional :: sendcounts,sdispls
      integer, dimension(0:), intent(in), optional :: recvcounts,rdispls
      integer(fmpi_integer), intent(out), optional :: request
      type(fmpi_win), intent(out), optional :: win
      integer,intent(in), optional :: comm
      type(f_enumerator), intent(in), optional :: algorithm
      ! Local variables
      integer :: algo,nproc,jproc
      integer(fmpi_integer) :: ierr,cnt
      integer(f_long) :: sizets,sizetr
      integer, dimension(:), allocatable :: nsenddspls_remote
      type(dictionary), pointer :: dict_info
      type(fmpi_win) :: window
      external :: MPI_ALLTOALL,MPI_ALLTOALLV
      
      nproc=mpisize(comm)
      sizets=f_size(sendbuf)
      sizetr=f_size(recvbuf)
      if (present(count)) then
         cnt=count
      else
         !division of the components per nproc
         cnt=min(sizets,sizetr)/nproc
      end if
      !otherwise determine algorithm
      algo=alltoall_algorithm(nproc,sizets,sizetr,&
           cnt,sendcounts,sdispls,&
           cnt,recvcounts,rdispls,algorithm,comm)
      if (present(request)) request=FMPI_REQUEST_NULL
      if (present(win)) algo=VARIABLE_ONE_SIDED_GET
      select case(algo)
      case(NOT_VARIABLE)
         call MPI_ALLTOALL(sendbuf,cnt,mpitype(sendbuf), &
              recvbuf,cnt,mpitype(recvbuf),comm,ierr)
      case(VARIABLE)
         if (present(request)) then
            call MPI_IALLTOALLV(sendbuf,sendcounts,sdispls,mpitype(sendbuf), &
                 recvbuf,recvcounts,rdispls,mpitype(recvbuf),comm,request,ierr)
         else
            call MPI_ALLTOALLV(sendbuf,sendcounts,sdispls,mpitype(sendbuf), &
                 recvbuf,recvcounts,rdispls,mpitype(recvbuf),comm,ierr)
         end if
      case(VARIABLE_ONE_SIDED_GET)
         nsenddspls_remote = f_malloc(0.to.nproc-1,id='nsenddspls_remote')
         !call fmpi_alltoall(sendbuf=sdispls,recvbuf=nsenddspls_remote,count=1,comm=comm) !TO BE ADDED ASAP
         !info=mpiinfo("no_locks", "true")
         dict_info=>dict_new('nolocks' .is. 'true')
         !window = mpiwindow(size(sendbuf), sendbuf, comm)
         call fmpi_win_create(window,sendbuf(1),size=sizets,dict_info=dict_info,comm=comm)
         call fmpi_win_fence(window,FMPI_WIN_OPEN)

         call dict_free(dict_info)
         do jproc=0,nproc-1
            if (recvcounts(jproc)>0) then
               call fmpi_get(origin_addr=recvbuf(1),origin_displ=rdispls(jproc),target_rank=jproc,&
                    target_disp=int(nsenddspls_remote(jproc),fmpi_address),count=recvcounts(jproc),win=window)
            end if
         end do
         if (present(win)) then
            win=window
            !there should be no need to nullify window
         else
            call fmpi_win_fence(window,FMPI_WIN_CLOSE)
            call fmpi_win_free(window)
         end if
         !call mpiinfofree(info)
         
         call f_free(nsenddspls_remote)
         ierr=FMPI_SUCCESS
      end select
      if (ierr/=FMPI_SUCCESS) then
         call f_err_throw('An error in calling to FMPI_ALLTOALL occured',&
              err_id=ERR_MPI_WRAPPERS)
         return
      end if

    end subroutine mpialltoallw_d11

!!$    subroutine mpialltoallv_long11(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, comm, algorithm)
!!$      use dictionaries, only: f_err_throw,f_err_define
!!$      use dynamic_memory
!!$      use yaml_output
!!$      !use iso_c_binding
!!$      implicit none
!!$      integer(f_long),dimension(:),intent(in),target :: sendbuf
!!$      integer(f_long),dimension(:),intent(out),target :: recvbuf
!!$      integer(f_long),dimension(:),pointer :: sendbuf_1d
!!$      integer(f_long),dimension(:),pointer :: recvbuf_1d
!!$      include 'alltoallv-inc.f90'
!!$    end subroutine mpialltoallv_long11
!!$
!!$    subroutine mpialltoallv_double11(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, comm, algorithm)
!!$      use dictionaries, only: f_err_throw,f_err_define
!!$      use dynamic_memory
!!$      use yaml_output
!!$      implicit none
!!$      double precision,dimension(:),intent(in),target :: sendbuf
!!$      double precision,dimension(:),intent(out),target :: recvbuf
!!$      double precision,dimension(:),pointer :: sendbuf_1d
!!$      double precision,dimension(:),pointer :: recvbuf_1d
!!$      include 'alltoallv-inc.f90'
!!$    end subroutine mpialltoallv_double11
!!$
!!$    !!subroutine mpialltoallv_double61(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, comm)
!!$    !!  use dictionaries, only: f_err_throw,f_err_define
!!$    !!  use dynamic_memory
!!$    !!  use yaml_output
!!$    !!  use iso_c_binding
!!$    !!  implicit none
!!$    !!  double precision,dimension(:,:,:,:,:,:),intent(in),target :: sendbuf
!!$    !!  double precision,dimension(:),intent(out),target :: recvbuf
!!$    !!  double precision,dimension(:),pointer :: sendbuf_1d
!!$    !!  double precision,dimension(:),pointer :: recvbuf_1d
!!$    !!  include 'alltoallv-inc.f90'
!!$    !!end subroutine mpialltoallv_double61
!!$
!!$    subroutine mpiialltoallv_double(sendbuf, sendcounts, senddspls, sendtype, &
!!$         recvbuf, recvcounts, recvdspls, recvtype, comm, request)
!!$      use dictionaries, only: f_err_throw,f_err_define
!!$      implicit none
!!$      ! Calling arguments
!!$      integer,intent(in) :: sendcounts, senddspls, sendtype, recvcounts, recvdspls, recvtype, comm
!!$      double precision,intent(in) :: sendbuf
!!$      double precision,intent(out) :: recvbuf
!!$      integer,intent(out) :: request
!!$      ! Local variables
!!$      integer :: ierr
!!$
!!$#ifdef HAVE_MPI3
!!$      call mpi_ialltoallv(sendbuf, sendcounts, senddspls, sendtype, &
!!$           recvbuf, recvcounts, recvdspls, recvtype, comm, request, ierr)
!!$      if (ierr/=0) then
!!$         call f_err_throw('An error in calling to MPI_IALLTOALLV occured',&
!!$              err_id=ERR_MPI_WRAPPERS)
!!$         return
!!$      end if
!!$#else
!!$      call mpi_alltoallv(sendbuf, sendcounts, senddspls, sendtype, &
!!$           recvbuf, recvcounts, recvdspls, recvtype, comm, ierr)
!!$      if (ierr/=0) then
!!$         call f_err_throw('An error in calling to MPI_IALLTOALLV occured',&
!!$              err_id=ERR_MPI_WRAPPERS)
!!$         return
!!$      end if
!!$      request = MPI_REQUEST_NULL
!!$#endif
!!$
!!$    end subroutine mpiialltoallv_double
!!$
!!$    subroutine mpi_get_alltoallv_i(iproc, nproc, comm, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, sendbuf, recvbuf)
!!$      use dynamic_memory
!!$      implicit none
!!$      integer(f_integer),dimension(:),intent(in) :: sendbuf
!!$      integer(f_integer),dimension(:),intent(in) :: recvbuf
!!$
!!$      include 'mpi_get_alltoallv-inc.f90'
!!$
!!$    end subroutine mpi_get_alltoallv_i
!!$
!!$    subroutine mpi_get_alltoallv_l(iproc, nproc, comm, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, sendbuf, recvbuf)
!!$      use dynamic_memory
!!$      implicit none
!!$      integer(f_long),dimension(:),intent(in) :: sendbuf
!!$      integer(f_long),dimension(:),intent(in) :: recvbuf
!!$
!!$      include 'mpi_get_alltoallv-inc.f90'
!!$
!!$    end subroutine mpi_get_alltoallv_l
!!$
!!$    subroutine mpi_get_alltoallv_d(iproc, nproc, comm, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, sendbuf, recvbuf)
!!$      use dynamic_memory
!!$      implicit none
!!$      double precision,dimension(:),intent(in) :: sendbuf
!!$      double precision,dimension(:),intent(in) :: recvbuf
!!$
!!$      include 'mpi_get_alltoallv-inc.f90'
!!$
!!$    end subroutine mpi_get_alltoallv_d
!!$
!!$
end module f_alltoall
