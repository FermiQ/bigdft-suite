subroutine post_p2p_communication(iproc, nproc, nsendbuf, sendbuf, nrecvbuf, recvbuf, comm)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
  real(8),dimension(nsendbuf),intent(in):: sendbuf
  real(8),dimension(nrecvbuf),intent(out):: recvbuf
  type(p2pComms),intent(inout):: comm
  
  ! Local variables
  integer:: jproc, joverlap, nsends, nreceives, mpisource, istsource, ncount, mpidest, istdest, tag, ierr
  
  nreceives=0
  nsends=0
  ! First only post receives
  do jproc=0,nproc-1
      do joverlap=1,comm%noverlaps(jproc)
          mpisource=comm%comarr(1,joverlap,jproc)
          istsource=comm%comarr(2,joverlap,jproc)
          ncount=comm%comarr(3,joverlap,jproc)
          mpidest=comm%comarr(4,joverlap,jproc)
          istdest=comm%comarr(5,joverlap,jproc)
          tag=comm%comarr(6,joverlap,jproc)
          if(nproc>1) then
              if(iproc==mpidest) then
                  nreceives=nreceives+1
                  call mpi_irecv(recvbuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
                       comm%requests(nreceives,2), ierr)
              end if
          else
              nsends=nsends+1
              nreceives=nreceives+1
              call dcopy(ncount, sendbuf(istsource), 1, recvbuf(istdest), 1)
          end if
      end do
  end do
  
  ! Now the sends
  do jproc=0,nproc-1
      do joverlap=1,comm%noverlaps(jproc)
          mpisource=comm%comarr(1,joverlap,jproc)
          istsource=comm%comarr(2,joverlap,jproc)
          ncount=comm%comarr(3,joverlap,jproc)
          mpidest=comm%comarr(4,joverlap,jproc)
          istdest=comm%comarr(5,joverlap,jproc)
          tag=comm%comarr(6,joverlap,jproc)
          if(nproc>1) then
              if(iproc==mpisource) then
                  nsends=nsends+1
                  call mpi_isend(sendbuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
                       comm%requests(nsends,1), ierr)
              end if
          end if
      end do
  end do
  
  !!if(iproc==0) write(*,'(a)') 'done.'
  
  comm%nsend=nsends
  comm%nrecv=nreceives
  
  if(nreceives/=comm%noverlaps(iproc)) then
      write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comm%noverlaps(iproc)',&
           nreceives, comm%noverlaps(iproc)
    stop
  end if
  
  ! Flag indicating whether the communication is complete or not
  if(nproc>1) then
      comm%communication_complete=.false.
  else
      comm%communication_complete=.true.
  end if
  write(*,*) 'after post: comm%communication_complete:',comm%communication_complete


end subroutine post_p2p_communication


subroutine wait_p2p_communication(iproc, nproc, comm)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(p2pComms),intent(inout):: comm
  
  ! Local variables
  integer:: ierr, ind, i, nsend, nrecv
  
  
  write(*,*) 'comm%communication_complete',comm%communication_complete
  if(.not.comm%communication_complete) then

      ! Wait for the sends to complete.
      nsend=0
      if(comm%nsend>0) then
          wait_sends: do
             call mpi_waitany(comm%nsend-nsend, comm%requests(1,1), ind, mpi_status_ignore, ierr)
             nsend=nsend+1
             do i=ind,comm%nsend-nsend
                comm%requests(i,1)=comm%requests(i+1,1)
             end do
             if(nsend==comm%nsend) exit wait_sends
          end do wait_sends
      end if
 
 
      ! Wait for the receives to complete.
      nrecv=0
      if(comm%nrecv>0) then
          wait_recvs: do
             call mpi_waitany(comm%nrecv-nrecv, comm%requests(1,2), ind, mpi_status_ignore, ierr)
             nrecv=nrecv+1
             do i=ind,comm%nrecv-nrecv
                comm%requests(i,2)=comm%requests(i+1,2)
             end do
             if(nrecv==comm%nrecv) exit wait_recvs
          end do wait_recvs
      end if

  end if

  ! Flag indicating that the communication is complete
  comm%communication_complete=.true.

end subroutine wait_p2p_communication
