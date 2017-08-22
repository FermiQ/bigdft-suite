
  integer,dimension(0:),intent(in) :: sendcounts, sdispls
  integer,dimension(0:),intent(in) :: recvcounts, rdispls
  integer,intent(in) :: comm
  ! Local variables
  integer :: ierr, iproc, nproc
  logical :: large
  integer,parameter :: max_size = 100000000 !maximal number of elements sent and/or received on a proc using the standard MPI function

  ! Check whether the sizes are correct
  if (size(sendbuf)<sum(sendcounts)) then
      call f_err_throw("insufficient size of array 'sendbuf': "//trim(yaml_toa(size(sendbuf)))//&
           &" instead of "//trim(yaml_toa(sum(sendcounts))))
  end if
  if (size(recvbuf)<sum(recvcounts)) then
      call f_err_throw("insufficient size of array 'recvbuf': "//trim(yaml_toa(size(recvbuf)))//&
           &" instead of "//trim(yaml_toa(sum(recvcounts))))
  end if

  ! Check whether we are having a "big" case. If so, use the hand-made
  ! version using mpi_get, otherwise the standard function.
  large = .false.
  if (sum(sendcounts)>max_size .or. sum(recvcounts)>max_size) then
      large = .true.
  end if
  call mpiallred(large, 1, mpi_lor, comm=comm)

  ! Now call the corresponding MPI function
  if (large) then

      iproc = mpirank(comm)
      nproc = mpisize(comm)
      if (iproc==0) then
          call yaml_warning('Large arrays in this call to mpi_alltoallv, &
               &using mpi_get_alltoallv instead')
      end if
      !call mpi_get_alltoallv(iproc, nproc, comm, sendcounts, sdispls, &
      !     recvcounts, rdispls, sendbuf, recvbuf)
      call c_f_pointer(c_loc(sendbuf), sendbuf_1d, [size(sendbuf)])
      call c_f_pointer(c_loc(recvbuf), recvbuf_1d, [size(recvbuf)])
      call mpi_get_alltoallv(iproc, nproc, comm, sendcounts, sdispls, &
           recvcounts, rdispls, sendbuf_1d, recvbuf_1d)
  else
  
      call mpi_alltoallv(sendbuf, sendcounts, sdispls, mpitype(sendbuf), &
           recvbuf, recvcounts, rdispls, mpitype(recvbuf), comm, ierr)
      
      if (ierr/=0) then
         call f_err_throw('An error in calling to MPI_ALLTOALLV occured',&
              err_id=ERR_MPI_WRAPPERS)
         return
      end if

  end if
