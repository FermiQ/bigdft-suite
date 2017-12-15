  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(comms_cubic), intent(in) :: comms
  !local variables
  external :: switch_waves_v,psitransspi

  call timing(iproc,'Un-TransSwitch','ON')

!!$  if (nproc > 1 .and. .not. present(workbuf)) then
!!$     call f_err_throw('Workbuf must be present when the number of processor is high',&
!!$          err_name='BIGDFT_RUNTIME_ERROR')
!!$  end if

  if (nproc > 1) then
     call switch_waves_v(nproc,orbs,&
          wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par(0,1),sendbuf,workbuf)

     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     if (present(recvbuf)) then
        call fmpi_alltoall(sendbuf=workbuf,sendcounts=comms%ncntd,sdispls=comms%ndspld, &
             recvbuf=recvbuf,recvcounts=comms%ncntt,rdispls=comms%ndsplt,comm=bigdft_mpi%mpi_comm)
     else
        call fmpi_alltoall(sendbuf=workbuf,sendcounts=comms%ncntd,sdispls=comms%ndspld, &
             recvbuf=sendbuf,recvcounts=comms%ncntt,rdispls=comms%ndsplt,comm=bigdft_mpi%mpi_comm)
     end if
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
  else
     if(orbs%nspinor /= 1) then
        !for only one processor there is no need to transform this
        call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,sendbuf,.true.)
     end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')
