!> @file
!! Wrapper for mpi_allreduce flavours
!! Use error handling
!! @author
!!    Copyright (C) 2012-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_allreduce
  use time_profiling, only: TIMING_UNINITIALIZED
  use f_precisions
  use f_enums
  use fmpi_types
  use f_utils, only: f_size
  use yaml_strings
  use dictionaries, only: f_err_throw
  use dynamic_memory
  implicit none
  private
  
  !> Interface for MPI_ALLREDUCE routine, to be updated little by little
  interface fmpi_allreduce
     module procedure mpiallred_int,mpiallred_real, &
          & mpiallred_double,&!,mpiallred_double_1,mpiallred_double_2,&
          & mpiallred_log
     module procedure mpiallred_long
     module procedure mpiallred_r1,mpiallred_r2,mpiallred_r3,mpiallred_r4
     module procedure mpiallred_d1,mpiallred_d2,mpiallred_d3,mpiallred_d4,mpiallred_d5
     module procedure mpiallred_i1,mpiallred_i2,mpiallred_i3
     module procedure mpiallred_l3
  end interface fmpi_allreduce

!!$  interface mpiiallred
!!$     module procedure mpiiallred_double
!!$  end interface mpiiallred

  public :: fmpi_allreduce

  logical, parameter :: have_mpi2=.true.

  integer, public, save :: TCAT_ALLRED_SMALL = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_ALLRED_LARGE = TIMING_UNINITIALIZED
  integer, parameter, public :: smallsize=5 !< limit for a communication with small size

  contains
   
!!$    subroutine mpiiallred_double(sendbuf, recvbuf, ncount, op, comm, request)
!!$      use dictionaries, only: f_err_throw,f_err_define
!!$      implicit none
!!$      ! Calling arguments
!!$      integer,intent(in) :: ncount, op, comm
!!$      double precision,intent(in) :: sendbuf
!!$      double precision,intent(out) :: recvbuf
!!$      integer,intent(out) :: request
!!$      ! Local variables
!!$      integer :: ierr
!!$
!!$#ifdef HAVE_MPI3
!!$      call mpi_iallreduce(sendbuf, recvbuf, ncount, mpitype(sendbuf), op, comm, request, ierr)
!!$      if (ierr/=0) then
!!$         call f_err_throw('An error in calling to MPI_IALLREDUCE occured',&
!!$              err_id=ERR_MPI_WRAPPERS)
!!$         return
!!$      end if
!!$#else
!!$      call mpi_allreduce(sendbuf, recvbuf, ncount, mpitype(sendbuf), op, comm, ierr)
!!$      if (ierr/=0) then
!!$         call f_err_throw('An error in calling to MPI_ALLREDUCE occured',&
!!$              err_id=ERR_MPI_WRAPPERS)
!!$         return
!!$      end if
!!$      request = MPI_REQUEST_NULL
!!$#endif
!!$
!!$    end subroutine mpiiallred_double
    
    !> Interface for MPI_ALLREDUCE operations
    subroutine mpiallred_int(sendbuf,count,op,comm,recvbuf,request)
      use dictionaries, only: f_err_throw,f_err_define
      use dynamic_memory
      implicit none
      integer(f_integer) :: sendbuf
      integer(f_integer), intent(inout), optional :: recvbuf
      integer(f_integer), dimension(:), allocatable :: copybuf
      include 'allreduce-inc.f90'
    end subroutine mpiallred_int

    subroutine mpiallred_long(sendbuf,count,op,comm,recvbuf,request)
      use dictionaries, only: f_err_throw,f_err_define
      use dynamic_memory
      implicit none
      integer(f_long) :: sendbuf
      integer(f_long), intent(inout), optional :: recvbuf
      integer(f_long), dimension(:), allocatable :: copybuf
      include 'allreduce-inc.f90'
    end subroutine mpiallred_long

    !> Interface for MPI_ALLREDUCE operations
    subroutine mpiallred_real(sendbuf,count,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw,f_err_define
      implicit none
      real(f_simple) :: sendbuf
      real(f_simple), intent(inout), optional :: recvbuf
      real(f_simple), dimension(:), allocatable :: copybuf
      include 'allreduce-inc.f90'
    end subroutine mpiallred_real

    subroutine mpiallred_double(sendbuf,count,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw,f_err_define
      implicit none
      real(f_double) :: sendbuf
      real(f_double), intent(inout), optional :: recvbuf
      real(f_double), dimension(:), allocatable :: copybuf
      include 'allreduce-inc.f90'
    end subroutine mpiallred_double

    subroutine mpiallred_log(sendbuf,count,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw,f_err_define
      implicit none
      logical :: sendbuf
      logical, intent(inout), optional :: recvbuf
      logical, dimension(:), allocatable :: copybuf
      include 'allreduce-inc.f90'
    end subroutine mpiallred_log

    subroutine mpiallred_i1(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      integer, dimension(:), intent(inout) :: sendbuf
      integer, dimension(:), intent(inout), optional :: recvbuf
      integer, dimension(:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_i1

    subroutine mpiallred_i2(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      integer, dimension(:,:), intent(inout) :: sendbuf
      integer, dimension(:,:), intent(inout), optional :: recvbuf
      integer, dimension(:,:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_i2

    subroutine mpiallred_i3(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      integer, dimension(:,:,:), intent(inout) :: sendbuf
      integer, dimension(:,:,:), intent(inout), optional :: recvbuf
      integer, dimension(:,:,:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_i3

    subroutine mpiallred_l3(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      logical, dimension(:,:,:), intent(inout) :: sendbuf
      logical, dimension(:,:,:), intent(inout), optional :: recvbuf
      logical, dimension(:,:,:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_l3


    subroutine mpiallred_r1(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_refine
      !    use yaml_strings, only: yaml_toa
      implicit none
      real, dimension(:), intent(inout) :: sendbuf
      real, dimension(:), intent(inout), optional :: recvbuf
      real, dimension(:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_r1

    subroutine mpiallred_r2(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      real, dimension(:,:), intent(inout) :: sendbuf
      real, dimension(:,:), intent(inout), optional :: recvbuf
      real, dimension(:,:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_r2

    subroutine mpiallred_r3(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      real, dimension(:,:,:), intent(inout) :: sendbuf
      real, dimension(:,:,:), intent(inout), optional :: recvbuf
      real, dimension(:,:,:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_r3

    subroutine mpiallred_r4(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      real, dimension(:,:,:,:), intent(inout) :: sendbuf
      real, dimension(:,:,:,:), intent(inout), optional :: recvbuf
      real, dimension(:,:,:,:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_r4

    subroutine mpiallred_d1(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      real(f_double), dimension(:), intent(inout) :: sendbuf
      real(f_double), dimension(:), intent(inout), optional :: recvbuf
      real(f_double), dimension(:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_d1

    subroutine mpiallred_d2(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      real(f_double), dimension(:,:), intent(inout) :: sendbuf
      real(f_double), dimension(:,:), intent(inout), optional :: recvbuf
      real(f_double), dimension(:,:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_d2

    subroutine mpiallred_d3(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      real(f_double), dimension(:,:,:), intent(inout) :: sendbuf
      real(f_double), dimension(:,:,:), intent(inout), optional :: recvbuf
      real(f_double), dimension(:,:,:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_d3

    subroutine mpiallred_d4(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      real(f_double), dimension(:,:,:,:), intent(inout) :: sendbuf
      real(f_double), dimension(:,:,:,:), intent(inout), optional :: recvbuf
      real(f_double), dimension(:,:,:,:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_d4

    subroutine mpiallred_d5(sendbuf,op,comm,recvbuf,request)
      use dynamic_memory
      use dictionaries, only: f_err_throw!,f_err_define
      !    use yaml_strings, only: yaml_toa
      implicit none
      real(f_double), dimension(:,:,:,:,:), intent(inout) :: sendbuf
      real(f_double), dimension(:,:,:,:,:), intent(inout), optional :: recvbuf
      real(f_double), dimension(:,:,:,:,:), allocatable :: copybuf
      include 'allreduce-arr-inc.f90'
    end subroutine mpiallred_d5

end module f_allreduce
