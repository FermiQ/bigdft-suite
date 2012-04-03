subroutine init_collective_comms(iproc, nproc, orbs, lzd, collcom)
use module_base
use module_types
use module_interfaces, except_this_one => init_collective_comms
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
type(collective_comms),intent(out):: collcom

! Local variables
integer:: ii, istat, iorb, iiorb, ilr, iall, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
real(8),dimension(:,:,:),allocatable:: weight_c, weight_c_temp, weight_f, weight_f_temp
real(8):: weight_c_tot, weight_f_tot, weightp_c, weightp_f, tt, ierr, t1, t2
integer,dimension(:,:),allocatable:: istartend_c, istartend_f
integer,dimension(:,:,:),allocatable:: index_in_global_c, index_in_global_f
character(len=*),parameter:: subname='init_collective_comms'

allocate(weight_c(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
call memocc(istat, weight_c, 'weight_c', subname)
allocate(weight_f(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
call memocc(istat, weight_f, 'weight_f', subname)
allocate(index_in_global_c(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
call memocc(istat, index_in_global_c, 'index_in_global_c', subname)
allocate(index_in_global_f(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
call memocc(istat, index_in_global_f, 'index_in_global_f', subname)

call get_weights(iproc, nproc, orbs, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot)


!!
!!  tt=weight_tot
!!  call mpi_alzd%llreduce(tt, weight_tot, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!  weight_temp=weight
!!  ii=(lzd%glr%ie1-lzd%glr%is1+1)*(lzd%glr%ie2-lzd%glr%is2+1)*(lzd%glr%ie3-lzd%glr%is3+1)
!!  call mpi_alzd%llreduce(weight_temp, weight, ii, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!
  !!if(iproc==0) write(*,'(a,2es14.5)') 'total weights (coarse / fine):',weight_c_tot, weight_f_tot

!!  ! Assign the grid points to the processes such that the work is equally dsitributed
  allocate(istartend_c(2,0:nproc-1), stat=istat)
  call memocc(istat, istartend_c, 'istartend_c', subname)
  allocate(istartend_f(2,0:nproc-1), stat=istat)
  call memocc(istat, istartend_f, 'istartend_f', subname)
call assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot, &
     istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
     weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f)

  iall=-product(shape(weight_c))*kind(weight_c)
  deallocate(weight_c, stat=istat)
  call memocc(istat, iall, 'weight_c', subname)
  iall=-product(shape(weight_f))*kind(weight_f)
  deallocate(weight_f, stat=istat)
  call memocc(istat, iall, 'weight_f', subname)

!!  call assign_weight_to_process(iproc, nproc, lzd%glr, weight, weight_tot, istartend, weightp, collcom%nptsp)
  !!write(*,'(a,i0,a,i0,a,es14.5)') 'process ',iproc,' has ',collcom%nptsp_c,' coarse points and weight ',weightp_c
  !!write(*,'(a,i0,a,i0,a,es14.5)') 'process ',iproc,' has ',collcom%nptsp_f,' fine points and weight ',weightp_f
!!!!  if(istartend(2,nproc-1)/=(lzd%glr%ie1-lzd%glr%is1+1)*(lzd%glr%ie2-lzd%glr%is2+1)*(lzd%glr%ie3-lzd%glr%is3+1)) stop 'istartend(2,nproc-1)/=(lzd%glr%ie1-lzd%glr%is1+1)*(lzd%glr%ie2-lzd%glr%is2+1)*(lzd%glr%ie3-lzd%glr%is3+1)'
!!

  ! some checks
  if(nproc>1) then
      call mpi_allreduce(weightp_c, tt, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  else
      tt=weightp_c
  end if
  if(tt/=weight_c_tot) stop 'wrong partition of coarse weights'
  if(nproc>1) then
      call mpi_allreduce(weightp_f, tt, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  else
      tt=weightp_f
  end if     
  if(tt/=weight_f_tot) stop 'wrong partition of fine weights'
  if(nproc>1) then
      call mpi_allreduce(collcom%nptsp_c, ii, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  else
      ii=collcom%nptsp_c
  end if
  if(ii/=lzd%glr%wfd%nvctr_c) stop 'wrong partition of coarse grid points'
  if(nproc>1) then
      call mpi_allreduce(collcom%nptsp_f, ii, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  else
      ii=collcom%nptsp_f
  end if
  if(ii/=lzd%glr%wfd%nvctr_f) stop 'init_collective_comms: wrong partition of fine grid points'

!!  ! Allocate the keys
  allocate(collcom%norb_per_gridpoint_c(collcom%nptsp_c), stat=istat)
  call memocc(istat, collcom%norb_per_gridpoint_c, 'collcom%norb_per_gridpoint_c', subname)
  allocate(collcom%norb_per_gridpoint_f(collcom%nptsp_f), stat=istat)
  call memocc(istat, collcom%norb_per_gridpoint_f, 'collcom%norb_per_gridpoint_f', subname)
  call mpi_barrier(mpi_comm_world, ierr)
  !!t1=mpi_wtime()
  call determine_num_orbs_per_gridpoint(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
       istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
       weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, &
       collcom%norb_per_gridpoint_c, collcom%norb_per_gridpoint_f)
  !!t2=mpi_wtime()
  !!write(*,'(a,i5,es15.5)') 'iproc, time determine_num_orbs_per_gridpoint:',iproc, t2-t1

  ! Determine the index of a grid point i1,i2,i3 in the compressed array
  !!call mpi_barrier(mpi_comm_world, ierr)
  !!t1=mpi_wtime()
  call get_index_in_global2(lzd%glr, index_in_global_c, index_in_global_f)
  !!t2=mpi_wtime()
  !!write(*,'(a,i5,es15.5)') 'iproc, time get_index_in_global2:',iproc, t2-t1


  ! Determine values for mpi_alltoallv
  allocate(collcom%nsendcounts_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsendcounts_c, 'collcom%nsendcounts_c', subname)
  allocate(collcom%nsenddspls_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsenddspls_c, 'collcom%nsenddspls_c', subname)
  allocate(collcom%nrecvcounts_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvcounts_c, 'collcom%nrecvcounts_c', subname)
  allocate(collcom%nrecvdspls_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvdspls_c, 'collcom%nrecvdspls_c', subname)
  allocate(collcom%nsendcounts_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsendcounts_f, 'collcom%nsendcounts_f', subname)
  allocate(collcom%nsenddspls_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsenddspls_f, 'collcom%nsenddspls_f', subname)
  allocate(collcom%nrecvcounts_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvcounts_f, 'collcom%nrecvcounts_f', subname)
  allocate(collcom%nrecvdspls_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvdspls_f, 'collcom%nrecvdspls_f', subname)
  !!call mpi_barrier(mpi_comm_world, ierr)
  !!t1=mpi_wtime()
  call determine_communication_arrays(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
       index_in_global_c, index_in_global_f, weightp_c, weightp_f, &
       collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, &
       collcom%nsendcounts_f, collcom%nsenddspls_f, collcom%nrecvcounts_f, collcom%nrecvdspls_f)
  !!t2=mpi_wtime()
  !!write(*,'(a,i5,es15.5)') 'iproc, time determine_communication_arrays:',iproc, t2-t1
!!
!!
  ! Now rearrange the data on the process to communicate them
  collcom%ndimpsi_c=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      collcom%ndimpsi_c=collcom%ndimpsi_c+lzd%llr(ilr)%wfd%nvctr_c
  end do
  allocate(collcom%irecvbuf_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, collcom%irecvbuf_c, 'collcom%irecvbuf_c', subname)
  allocate(collcom%indexrecvorbital_c(sum(collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, collcom%indexrecvorbital_c, 'collcom%indexrecvorbital_c', subname)
  allocate(collcom%iextract_c(sum(collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, collcom%iextract_c, 'collcom%iextract_c', subname)
  allocate(collcom%iexpand_c(sum(collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, collcom%iexpand_c, 'collcom%iexpand_c', subname)
  allocate(collcom%isendbuf_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, collcom%isendbuf_c, 'collcom%isendbuf_c', subname)

  collcom%ndimpsi_f=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      collcom%ndimpsi_f=collcom%ndimpsi_f+lzd%llr(ilr)%wfd%nvctr_f
  end do
  allocate(collcom%irecvbuf_f(collcom%ndimpsi_f), stat=istat)
  call memocc(istat, collcom%irecvbuf_f, 'collcom%irecvbuf_f', subname)
  allocate(collcom%indexrecvorbital_f(sum(collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, collcom%indexrecvorbital_f, 'collcom%indexrecvorbital_f', subname)
  allocate(collcom%iextract_f(sum(collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, collcom%iextract_f, 'collcom%iextract_f', subname)
  allocate(collcom%iexpand_f(sum(collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, collcom%iexpand_f, 'collcom%iexpand_f', subname)
  allocate(collcom%isendbuf_f(collcom%ndimpsi_f), stat=istat)
  call memocc(istat, collcom%isendbuf_f, 'collcom%isendbuf_f', subname)

  !!call mpi_barrier(mpi_comm_world, ierr)
  !!t1=mpi_wtime()
  call get_switch_indices(iproc, nproc, orbs, lzd, collcom%ndimpsi_c, collcom%ndimpsi_f, istartend_c, istartend_f, &
       collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, &
       collcom%nsendcounts_f, collcom%nsenddspls_f, collcom%nrecvcounts_f, collcom%nrecvdspls_f, &
       index_in_global_c, index_in_global_f, &
       weightp_c, weightp_f, collcom%isendbuf_c, collcom%irecvbuf_c, collcom%isendbuf_f, collcom%irecvbuf_f, &
       collcom%indexrecvorbital_c, collcom%iextract_c, collcom%iexpand_c, &
       collcom%indexrecvorbital_f, collcom%iextract_f, collcom%iexpand_f)
  !!t2=mpi_wtime()
  !!write(*,'(a,i5,es15.5)') 'iproc, time get_switch_indices:',iproc, t2-t1

  iall=-product(shape(istartend_c))*kind(istartend_c)
  deallocate(istartend_c, stat=istat)
  call memocc(istat, iall, 'istartend_c', subname)

  iall=-product(shape(istartend_f))*kind(istartend_f)
  deallocate(istartend_f, stat=istat)
  call memocc(istat, iall, 'istartend_f', subname)

  iall=-product(shape(index_in_global_c))*kind(index_in_global_c)
  deallocate(index_in_global_c, stat=istat)
  call memocc(istat, iall, 'index_in_global_c', subname)

  iall=-product(shape(index_in_global_f))*kind(index_in_global_f)
  deallocate(index_in_global_f, stat=istat)
  call memocc(istat, iall, 'index_in_global_f', subname)
  
  
end subroutine init_collective_comms

subroutine get_weights(iproc, nproc, orbs, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(out):: weight_c, weight_f
real(8),intent(out):: weight_c_tot, weight_f_tot

! Local variables
integer:: iorb, iiorb, i0, i1, i2, i3, ii, jj, iseg, ierr, ilr, istart, iend, i, j0, j1, ii1, ii2, ii3


  weight_c=0.d0
  weight_f=0.d0
  weight_c_tot=0.d0
  weight_f_tot=0.d0

  ! coarse part
  do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inwhichlocreg(iiorb)
    do iseg=1,lzd%llr(ilr)%wfd%nseg_c
       jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
       j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
       j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
       ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
       i2=ii/(lzd%llr(ilr)%d%n1+1)
       i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
       i1=i0+j1-j0
       !write(*,'(a,8i8)') 'jj, ii, j0, j1, i0, i1, i2, i3',jj,ii,j0,j1,i0,i1,i2,i3
       do i=i0,i1
          ii1=i+lzd%llr(ilr)%ns1
          ii2=i2+lzd%llr(ilr)%ns2
          ii3=i3+lzd%llr(ilr)%ns3
          weight_c(ii1,ii2,ii3)=weight_c(ii1,ii2,ii3)+1.d0
          weight_c_tot=weight_c_tot+1.d0
       enddo
    enddo
  
    ! fine part
    istart=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
    iend=istart+lzd%llr(ilr)%wfd%nseg_f-1
    do iseg=istart,iend
       jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
       j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
       j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
       ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
       i2=ii/(lzd%llr(ilr)%d%n1+1)
       i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
          ii1=i+lzd%llr(ilr)%ns1
          ii2=i2+lzd%llr(ilr)%ns2
          ii3=i3+lzd%llr(ilr)%ns3
          weight_f(ii1,ii2,ii3)=weight_f(ii1,ii2,ii3)+1.d0
          weight_f_tot=weight_f_tot+1.d0
       enddo
    enddo
  end do


  if(nproc>1) then
      call mpiallred(weight_c_tot, 1, mpi_sum, mpi_comm_world, ierr)
      call mpiallred(weight_f_tot, 1, mpi_sum, mpi_comm_world, ierr)
      ii=(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)
      call mpiallred(weight_c(0,0,0), ii,  mpi_sum, mpi_comm_world, ierr)
      call mpiallred(weight_f(0,0,0), ii,  mpi_sum, mpi_comm_world, ierr)
  end if


end subroutine get_weights



subroutine assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_tot_c, weight_tot_f, &
           istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, nptsp_c, nptsp_f)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: weight_c, weight_f
real(8),intent(in):: weight_tot_c, weight_tot_f
integer,dimension(2,0:nproc-1),intent(out):: istartend_c, istartend_f
integer,intent(out):: istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
real(8),intent(out):: weightp_c, weightp_f
integer,intent(out):: nptsp_c, nptsp_f

! Local variables
integer:: jproc, i1, i2, i3, ii, istartp_c, iendp_c, ii2, istartp_f, iendp_f, istart, iend, jj, j0, j1
integer:: i, iseg, i0, iitot, ierr, iiseg
real(8):: tt, tt2, weight_c_ideal, weight_f_ideal

  weight_c_ideal=weight_tot_c/dble(nproc)
  weight_f_ideal=weight_tot_f/dble(nproc)
  !!if(iproc==0) write(*,'(a,2es14.5)') 'ideal weight per process (coarse / fine):',weight_c_ideal,weight_f_ideal

  jproc=0
  tt=0.d0
  tt2=0.d0
  iitot=0
  ii2=0
  iiseg=1
  weightp_c=0.d0
    do iseg=1,lzd%glr%wfd%nseg_c
       jj=lzd%glr%wfd%keyvloc(iseg)
       j0=lzd%glr%wfd%keygloc(1,iseg)
       j1=lzd%glr%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
       ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
       i2=ii/(lzd%glr%d%n1+1)
       i0=ii-i2*(lzd%glr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
           tt=tt+weight_c(i,i2,i3)
           iitot=iitot+1
           if(tt>weight_c_ideal) then
               if(iproc==jproc) then
                   weightp_c=tt
                   nptsp_c=iitot
                   istartp_c=ii2+1
                   iendp_c=istartp_c+iitot-1
                   istartp_seg_c=iiseg
                   iendp_seg_c=iseg
               end if
               istartend_c(1,jproc)=ii2+1
               istartend_c(2,jproc)=istartend_c(1,jproc)+iitot-1
               tt2=tt2+tt
               tt=0.d0
               ii2=ii2+iitot
               iitot=0
               jproc=jproc+1
               iiseg=iseg
           end if
       end do
   end do
  if(iproc==nproc-1) then
      ! Take the rest
      istartp_c=ii2+1
      iendp_c=istartp_c+iitot-1
      weightp_c=weight_tot_c-tt2
      nptsp_c=lzd%glr%wfd%nvctr_c-ii2
      istartp_seg_c=iiseg
      iendp_seg_c=lzd%glr%wfd%nseg_c
  end if
  istartend_c(1,nproc-1)=ii2+1
  istartend_c(2,nproc-1)=istartend_c(1,nproc-1)+iitot-1

  ! some check
  ii=istartend_c(2,iproc)-istartend_c(1,iproc)+1
  if(nproc>1) call mpiallred(ii, 1, mpi_sum, mpi_comm_world, ierr)
  if(ii/=lzd%glr%wfd%nvctr_c) stop 'ii/=lzd%glr%wfd%nvctr_c'


  jproc=0
  tt=0.d0
  tt2=0.d0
  iitot=0
  ii2=0
  weightp_f=0.d0
  istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
  iend=istart+lzd%glr%wfd%nseg_f-1
  iiseg=istart
    do iseg=istart,iend
       jj=lzd%glr%wfd%keyvloc(iseg)
       j0=lzd%glr%wfd%keygloc(1,iseg)
       j1=lzd%glr%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
       ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
       i2=ii/(lzd%glr%d%n1+1)
       i0=ii-i2*(lzd%glr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
           tt=tt+weight_f(i,i2,i3)
           iitot=iitot+1
           if(tt>weight_f_ideal) then
               if(iproc==jproc) then
                   weightp_f=tt
                   nptsp_f=iitot
                   istartp_f=ii2+1
                   iendp_f=istartp_f+iitot-1
                   istartp_seg_f=iiseg
                   iendp_seg_f=iseg
               end if
               istartend_f(1,jproc)=ii2+1
               istartend_f(2,jproc)=istartend_f(1,jproc)+iitot-1
               tt2=tt2+tt
               tt=0.d0
               ii2=ii2+iitot
               iitot=0
               jproc=jproc+1
               iiseg=iseg
           end if
       end do
   end do
  if(iproc==nproc-1) then
      ! Take the rest
      istartp_f=ii2+1
      iendp_f=istartp_f+iitot-1
      weightp_f=weight_tot_f-tt2
      nptsp_f=lzd%glr%wfd%nvctr_f-ii2
      istartp_seg_f=iiseg
      iendp_seg_f=iend
  end if
  istartend_f(1,nproc-1)=ii2+1
  istartend_f(2,nproc-1)=istartend_f(1,nproc-1)+iitot-1

  ! some check
  ii=istartend_f(2,iproc)-istartend_f(1,iproc)+1
  if(nproc>1) call mpiallred(ii, 1, mpi_sum, mpi_comm_world, ierr)
  if(ii/=lzd%glr%wfd%nvctr_f) stop 'assign_weight_to_process: ii/=lzd%glr%wfd%nvctr_f'



end subroutine assign_weight_to_process



subroutine determine_num_orbs_per_gridpoint(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
           istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, nptsp_c, nptsp_f, &
           norb_per_gridpoint_c, norb_per_gridpoint_f)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nptsp_c, nptsp_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
real(8),intent(in):: weightp_c, weightp_f
integer,dimension(nptsp_c),intent(out):: norb_per_gridpoint_c
integer,dimension(nptsp_f),intent(out):: norb_per_gridpoint_f

! Local variables
integer:: ii, iiorb, i1, i2, i3, iipt, iorb, iii, npgp, iseg, jj, j0, j1, iitot, ilr, i, istart, iend, i0, istat, iall
logical:: found, overlap_possible
integer,dimension(:),allocatable:: iseg_start_c, iseg_start_f
character(len=*),parameter:: subname='determine_num_orbs_per_gridpoint'
!!real(8):: t1, t2, t1tot, t2tot, t_check_gridpoint

  allocate(iseg_start_c(lzd%nlr), stat=istat)
  call memocc(istat, iseg_start_c, 'iseg_start_c', subname)
  allocate(iseg_start_f(lzd%nlr), stat=istat)
  call memocc(istat, iseg_start_f, 'iseg_start_f', subname)

  iseg_start_c=1
  iseg_start_f=1

  iitot=0
  iiorb=0
  iipt=0
!!t_check_gridpoint=0.d0
!!t1tot=mpi_wtime()
  !write(*,*) 'iproc, istartp_seg_c,iendp_seg_c', iproc, istartp_seg_c,iendp_seg_c
    !do iseg=1,lzd%glr%wfd%nseg_c
    do iseg=istartp_seg_c,iendp_seg_c
       jj=lzd%glr%wfd%keyvloc(iseg)
       j0=lzd%glr%wfd%keygloc(1,iseg)
       j1=lzd%glr%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
       ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
       i2=ii/(lzd%glr%d%n1+1)
       i0=ii-i2*(lzd%glr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
           !iitot=iitot+1
           iitot=jj+i-i0
           if(iitot>=istartend_c(1,iproc) .and. iitot<=istartend_c(2,iproc)) then
               !write(200+iproc,'(5i10)') iitot, iseg, iitot, jj, jj+i-i0
               iipt=iipt+1
               npgp=0
               do iorb=1,orbs%norb
                   ilr=orbs%inwhichlocreg(iorb)
                   ! Check whether this orbitals extends here
                   call check_grid_point_from_boxes(i, i2, i3, lzd%llr(ilr), overlap_possible)
                   if(.not. overlap_possible) then
                       found=.false.
                   else
                       !!t1=mpi_wtime()
                       call check_gridpoint(lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
                            lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, lzd%llr(ilr)%wfd%keygloc, &
                            i, i2, i3, iseg_start_c(ilr), found)
                       !!t2=mpi_wtime()
                       !!t_check_gridpoint=t_check_gridpoint+t2-t1
                   end if
                   if(found) then
                       npgp=npgp+1
                       iiorb=iiorb+1
                   end if
               end do
               norb_per_gridpoint_c(iipt)=npgp
           end if
      end do
  end do

  if(iipt/=nptsp_c) stop 'iipt/=nptsp_c'
  if(iiorb/=nint(weightp_c)) stop 'iiorb/=weightp_c'



  iitot=0
  iiorb=0
  iipt=0
    istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
    iend=istart+lzd%glr%wfd%nseg_f-1
    !do iseg=istart,iend
    do iseg=istartp_seg_f,iendp_seg_f
       jj=lzd%glr%wfd%keyvloc(iseg)
       j0=lzd%glr%wfd%keygloc(1,iseg)
       j1=lzd%glr%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
       ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
       i2=ii/(lzd%glr%d%n1+1)
       i0=ii-i2*(lzd%glr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
           !iitot=iitot+1
           iitot=jj+i-i0
           if(iitot>=istartend_f(1,iproc) .and. iitot<=istartend_f(2,iproc)) then
               iipt=iipt+1
               npgp=0
               do iorb=1,orbs%norb
                   ilr=orbs%inwhichlocreg(iorb)
                   ! Check whether this orbitals extends here
                   call check_grid_point_from_boxes(i, i2, i3, lzd%llr(ilr), overlap_possible)
                   if(.not. overlap_possible) then
                       found=.false.
                   else
                       iii=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
                       !!t1=mpi_wtime()
                       call check_gridpoint(lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
                            lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
                            lzd%llr(ilr)%wfd%keygloc(1,iii), &
                            i, i2, i3, iseg_start_f(ilr), found)
                       !!t2=mpi_wtime()
                       !!t_check_gridpoint=t_check_gridpoint+t2-t1
                   end if
                   if(found) then
                       npgp=npgp+1
                       iiorb=iiorb+1
                   end if
               end do
               norb_per_gridpoint_f(iipt)=npgp
           end if
      end do
  end do

  if(iipt/=nptsp_f) stop 'iipt/=nptsp_f'
  !!write(*,*) 'iiorb, weightp_f', iiorb, weightp_f
  if(iiorb/=nint(weightp_f)) stop 'iiorb/=weightp_f'


  iall=-product(shape(iseg_start_c))*kind(iseg_start_c)
  deallocate(iseg_start_c, stat=istat)
  call memocc(istat, iall, 'iseg_start_c', subname)
  iall=-product(shape(iseg_start_f))*kind(iseg_start_f)
  deallocate(iseg_start_f, stat=istat)
  call memocc(istat, iall, 'iseg_start_f', subname)

!!t2tot=mpi_wtime()
!!write(*,'(a,es14.5)') 'in sub determine_num_orbs_per_gridpoint: iproc, total time', t2tot-t1tot
!!write(*,'(a,es14.5)') 'in sub determine_num_orbs_per_gridpoint: iproc, time for check_gridpoint', t_check_gridpoint

end subroutine determine_num_orbs_per_gridpoint



subroutine determine_communication_arrays(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
           index_in_global_c, index_in_global_f, &
           weightp_c, weightp_f,  nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c, &
           nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
integer,dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: index_in_global_c, index_in_global_f
real(8),intent(in):: weightp_c, weightp_f
integer,dimension(0:nproc-1),intent(out):: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
integer,dimension(0:nproc-1),intent(out):: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f

! Local variables
integer:: iorb, iiorb, i1, i2, i3, ii, jproc, jproctarget, ierr, jj, ilr, j0, j1, i0, i, ind
integer:: istat, ii1, ii2, ii3, iseg, istart, iend, iall
integer,dimension(:),allocatable:: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp
character(len=*),parameter:: subname='determine_communication_arrays'

  ! Determine values for mpi_alltoallv
  ! first nsendcounts
  nsendcounts_c=0
  do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inwhichlocreg(iiorb)
    do iseg=1,lzd%llr(ilr)%wfd%nseg_c
       jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
       j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
       j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
       ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
       i2=ii/(lzd%llr(ilr)%d%n1+1)
       i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
       i1=i0+j1-j0
       ii2=i2+lzd%llr(ilr)%ns2
       ii3=i3+lzd%llr(ilr)%ns3
       do i=i0,i1
          ii1=i+lzd%llr(ilr)%ns1
          !call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'c', ind)
          ind=index_in_global_c(ii1,ii2,ii3)
          jproctarget=-1
          do jproc=0,nproc-1
              if(ind>=istartend_c(1,jproc) .and. ind<=istartend_c(2,jproc)) then
                  jproctarget=jproc
                  exit
              end if
          end do
          !!if(jproctarget==-1) write(*,*) 'ind, lzd%glr%wfd%nvctr_c',ind, lzd%glr%wfd%nvctr_c
          nsendcounts_c(jproctarget)=nsendcounts_c(jproctarget)+1
        end do
     end do
   end do

   !!write(*,'(a,i3,3x,100i8)') 'iproc, istartend_f(2,:)', iproc, istartend_f(2,:)

  nsendcounts_f=0
  do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inwhichlocreg(iiorb)
    istart=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
    iend=istart+lzd%llr(ilr)%wfd%nseg_f-1
    do iseg=istart,iend
       jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
       j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
       j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
       ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
       i2=ii/(lzd%llr(ilr)%d%n1+1)
       i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
       i1=i0+j1-j0
       ii2=i2+lzd%llr(ilr)%ns2
       ii3=i3+lzd%llr(ilr)%ns3
       do i=i0,i1
          ii1=i+lzd%llr(ilr)%ns1
          !call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'f', ind)
          ind=index_in_global_f(ii1,ii2,ii3)
          do jproc=0,nproc-1
              if(ind>=istartend_f(1,jproc) .and. ind<=istartend_f(2,jproc)) then
                  jproctarget=jproc
                  exit
              end if
          end do
          nsendcounts_f(jproctarget)=nsendcounts_f(jproctarget)+1
      end do
    end do
   end do



  ! The first check is to make sure that there is no stop in case this process has no orbitals (in which case
  ! orbs%npsidim_orbs is 1 and not 0 as assumed by the check)
  if(orbs%npsidim_orbs>1 .and. sum(nsendcounts_c)+7*sum(nsendcounts_f)/=orbs%npsidim_orbs) &
      stop 'sum(nsendcounts_c)+sum(nsendcounts_f)/=orbs%npsidim_orbs'

  
  ! now nsenddspls
  nsenddspls_c(0)=0
  do jproc=1,nproc-1
      nsenddspls_c(jproc)=nsenddspls_c(jproc-1)+nsendcounts_c(jproc-1)
  end do
  nsenddspls_f(0)=0
  do jproc=1,nproc-1
      nsenddspls_f(jproc)=nsenddspls_f(jproc-1)+nsendcounts_f(jproc-1)
  end do



  ! now nrecvcounts
  ! use an mpi_alltoallv to gather the data
  allocate(nsendcounts_tmp(0:nproc-1), stat=istat)
  call memocc(istat, nsendcounts_tmp, 'nsendcounts_tmp', subname)
  allocate(nsenddspls_tmp(0:nproc-1), stat=istat)
  call memocc(istat, nsenddspls_tmp, 'nsenddspls_tmp', subname)
  allocate(nrecvcounts_tmp(0:nproc-1), stat=istat)
  call memocc(istat, nrecvcounts_tmp, 'nrecvcounts_tmp', subname)
  allocate(nrecvdspls_tmp(0:nproc-1), stat=istat)
  call memocc(istat, nrecvdspls_tmp, 'nrecvdspls_tmp', subname)
  nsendcounts_tmp=1
  nrecvcounts_tmp=1
  do jproc=0,nproc-1
      nsenddspls_tmp(jproc)=jproc
      nrecvdspls_tmp(jproc)=jproc
  end do
  if(nproc>1) then
      call mpi_alltoallv(nsendcounts_c, nsendcounts_tmp, nsenddspls_tmp, mpi_integer, nrecvcounts_c, &
           nrecvcounts_tmp, nrecvdspls_tmp, mpi_integer, mpi_comm_world, ierr)
      call mpi_alltoallv(nsendcounts_f, nsendcounts_tmp, nsenddspls_tmp, mpi_integer, nrecvcounts_f, &
           nrecvcounts_tmp, nrecvdspls_tmp, mpi_integer, mpi_comm_world, ierr)
  else
      nrecvcounts_c=nsendcounts_c
      nrecvcounts_f=nsendcounts_f
  end if
  iall=-product(shape(nsendcounts_tmp))*kind(nsendcounts_tmp)
  deallocate(nsendcounts_tmp, stat=istat)
  call memocc(istat, iall, 'nsendcounts_tmp', subname)
  iall=-product(shape(nsenddspls_tmp))*kind(nsenddspls_tmp)
  deallocate(nsenddspls_tmp, stat=istat)
  call memocc(istat, iall, 'nsenddspls_tmp', subname)
  iall=-product(shape(nrecvcounts_tmp))*kind(nrecvcounts_tmp)
  deallocate(nrecvcounts_tmp, stat=istat)
  call memocc(istat, iall, 'nrecvcounts_tmp', subname)
  iall=-product(shape(nrecvdspls_tmp))*kind(nrecvdspls_tmp)
  deallocate(nrecvdspls_tmp, stat=istat)
  call memocc(istat, iall, 'nrecvdspls_tmp', subname)

  ! now recvdspls
  nrecvdspls_c(0)=0
  do jproc=1,nproc-1
      nrecvdspls_c(jproc)=nrecvdspls_c(jproc-1)+nrecvcounts_c(jproc-1)
  end do
  nrecvdspls_f(0)=0
  do jproc=1,nproc-1
      nrecvdspls_f(jproc)=nrecvdspls_f(jproc-1)+nrecvcounts_f(jproc-1)
  end do

  !write(*,*) 'sum(nrecvcounts_c), nint(weightp_c)', sum(nrecvcounts_c), nint(weightp_c)
  !!write(*,*) 'sum(nrecvcounts_f), nint(weightp_f)', sum(nrecvcounts_f), nint(weightp_f)
  if(sum(nrecvcounts_c)/=nint(weightp_c)) stop 'sum(nrecvcounts_c)/=nint(nweightp_c)'
  if(sum(nrecvcounts_f)/=nint(weightp_f)) stop 'sum(nrecvcounts_f)/=nint(nweightp_f)'

end subroutine determine_communication_arrays


subroutine get_switch_indices(iproc, nproc, orbs, lzd, ndimpsi_c, ndimpsi_f, istartend_c, istartend_f, &
           nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c, &
           nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f, &
           index_in_global_c, index_in_global_f, &
           weightp_c, weightp_f,  isendbuf_c, irecvbuf_c, isendbuf_f, irecvbuf_f, &
           indexrecvorbital_c, iextract_c, iexpand_c, indexrecvorbital_f, iextract_f, iexpand_f)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, ndimpsi_c, ndimpsi_f
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
integer,dimension(0:nproc-1),intent(in):: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
integer,dimension(0:nproc-1),intent(in):: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
integer,dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: index_in_global_c, index_in_global_f
real(8),intent(in):: weightp_c, weightp_f
integer,dimension(ndimpsi_c),intent(out):: isendbuf_c, irecvbuf_c
integer,dimension(ndimpsi_f),intent(out):: isendbuf_f, irecvbuf_f
integer,dimension(sum(nrecvcounts_c)),intent(out):: indexrecvorbital_c, iextract_c, iexpand_c
integer,dimension(sum(nrecvcounts_f)),intent(out):: indexrecvorbital_f, iextract_f, iexpand_f

! Local variables
integer:: i, j, iorb, iiorb, i1, i2, i3, ind, jproc, jproctarget, ii, ierr, jj, iseg, iitot, ilr
integer:: istart, iend, indglob, ii1, ii2, ii3, jo, j1, i0, j0, istat, iall
integer,dimension(:),allocatable:: nsend, indexsendorbital2, gridpoint_start_c, gridpoint_start_f, indexrecvorbital2
real(8),dimension(:,:,:),allocatable:: weight_c, weight_f
integer,dimension(:),allocatable:: indexsendorbital_c, indexsendbuf_c, indexrecvbuf_c
integer,dimension(:),allocatable:: indexsendorbital_f, indexsendbuf_f, indexrecvbuf_f
character(len=*),parameter:: subname='get_switch_indices'
!real(8):: t1, t2, t1tot, t2tot, t_reverse

!!t_reverse=0.d0
!!t1tot=mpi_wtime()

allocate(indexsendorbital_c(ndimpsi_c), stat=istat)
call memocc(istat, indexsendorbital_c, 'indexsendorbital_c', subname)
allocate(indexsendbuf_c(ndimpsi_c), stat=istat)
call memocc(istat, indexsendbuf_c, 'indexsendbuf_c', subname)
allocate(indexrecvbuf_c(sum(nrecvcounts_c)), stat=istat)
call memocc(istat, indexrecvbuf_c, 'indexrecvbuf_c', subname)

allocate(indexsendorbital_f(ndimpsi_f), stat=istat)
call memocc(istat, indexsendorbital_f, 'indexsendorbital_f', subname)
allocate(indexsendbuf_f(ndimpsi_f), stat=istat)
call memocc(istat, indexsendbuf_f, 'indexsendbuf_f', subname)
allocate(indexrecvbuf_f(sum(nrecvcounts_f)), stat=istat)
call memocc(istat, indexrecvbuf_f, 'indexrecvbuf_f', subname)

allocate(weight_c(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
call memocc(istat, weight_c, 'weight_c', subname)
allocate(weight_f(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
call memocc(istat, weight_f, 'weight_f', subname)
allocate(gridpoint_start_c((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)), stat=istat)
call memocc(istat, gridpoint_start_c, 'gridpoint_start_c', subname)
allocate(gridpoint_start_f((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)), stat=istat)
call memocc(istat, gridpoint_start_f, 'gridpoint_start_f', subname)
gridpoint_start_c=-1
gridpoint_start_f=-1

!!write(*,*) 'ndimpsi_f, sum(nrecvcounts_f)', ndimpsi_f, sum(nrecvcounts_f)

  allocate(nsend(0:nproc-1), stat=istat)
  call memocc(istat, nsend, 'nsend', subname)

  iitot=0
  nsend=0
  do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inwhichlocreg(iiorb)
    do iseg=1,lzd%llr(ilr)%wfd%nseg_c
       jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
       j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
       j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
       ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
       i2=ii/(lzd%llr(ilr)%d%n1+1)
       i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
       i1=i0+j1-j0
       !write(*,'(a,8i8)') 'jj, ii, j0, j1, i0, i1, i2, i3',jj,ii,j0,j1,i0,i1,i2,i3
       do i=i0,i1
          ii1=i+lzd%llr(ilr)%ns1
          ii2=i2+lzd%llr(ilr)%ns2
          ii3=i3+lzd%llr(ilr)%ns3
          !!call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'c', indglob)
          indglob=index_in_global_c(ii1,ii2,ii3)
              iitot=iitot+1
              do jproc=0,nproc-1
                  if(indglob>=istartend_c(1,jproc) .and. indglob<=istartend_c(2,jproc)) then
                      jproctarget=jproc
                      exit
                  end if
              end do
              !!write(600+iproc,'(a,2(i0,1x),i0,a,i0)') 'point ',ii1,ii2,ii3,' goes to process ',jproctarget
              nsend(jproctarget)=nsend(jproctarget)+1
              ind=nsenddspls_c(jproctarget)+nsend(jproctarget)
              isendbuf_c(iitot)=ind
              indexsendbuf_c(ind)=indglob
              indexsendorbital_c(iitot)=iiorb
              !indexsendorbital(ind)=iiorb
          end do
      end do
  end do

  if(iitot/=ndimpsi_c) stop 'iitot/=ndimpsi_c'

  !check
  do jproc=0,nproc-1
      if(nsend(jproc)/=nsendcounts_c(jproc)) stop 'nsend(jproc)/=nsendcounts_c(jproc)'
  end do


  ! fine part
  iitot=0
  nsend=0
  do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inwhichlocreg(iiorb)
    istart=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
    iend=istart+lzd%llr(ilr)%wfd%nseg_f-1
    do iseg=istart,iend
       jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
       j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
       j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
       ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
       i2=ii/(lzd%llr(ilr)%d%n1+1)
       i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
       i1=i0+j1-j0
       !write(*,'(a,8i8)') 'jj, ii, j0, j1, i0, i1, i2, i3',jj,ii,j0,j1,i0,i1,i2,i3
       do i=i0,i1
          ii1=i+lzd%llr(ilr)%ns1
          ii2=i2+lzd%llr(ilr)%ns2
          ii3=i3+lzd%llr(ilr)%ns3
          !!call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'f', indglob)
          indglob=index_in_global_f(ii1,ii2,ii3)
                  iitot=iitot+1
                  do jproc=0,nproc-1
                      if(indglob>=istartend_f(1,jproc) .and. indglob<=istartend_f(2,jproc)) then
                          jproctarget=jproc
                          exit
                      end if
                  end do
                  nsend(jproctarget)=nsend(jproctarget)+1
                  ind=nsenddspls_f(jproctarget)+nsend(jproctarget)
                  isendbuf_f(iitot)=ind
                  indexsendbuf_f(ind)=indglob
                  indexsendorbital_f(iitot)=iiorb
                  !indexsendorbital(ind)=iiorb
          end do
      end do
  end do

  if(iitot/=ndimpsi_f) stop 'iitot/=ndimpsi_f'

  !check
  do jproc=0,nproc-1
      !!write(*,*) 'nsend(jproc), nsendcounts_f(jproc)', nsend(jproc), nsendcounts_f(jproc)
      if(nsend(jproc)/=nsendcounts_f(jproc)) stop 'nsend(jproc)/=nsendcounts_f(jproc)'
  end do






  allocate(indexsendorbital2(ndimpsi_c), stat=istat)
  call memocc(istat, indexsendorbital2, 'indexsendorbital2', subname)
  indexsendorbital2=indexsendorbital_c
  do i=1,ndimpsi_c
      ind=isendbuf_c(i)
      indexsendorbital_c(ind)=indexsendorbital2(i)
  end do
  ! Inverse of isendbuf
!!t1=mpi_wtime()
  call get_reverse_indices(ndimpsi_c, isendbuf_c, irecvbuf_c)
  !!do i=1,ndimpsi_c
  !!    do j=1,ndimpsi_c
  !!        if(isendbuf_c(j)==i) then
  !!            irecvbuf_c(i)=j
  !!        end if
  !!    end do
  !!end do
!!t2=mpi_wtime()
!!t_reverse=t_reverse+t2-t1
  iall=-product(shape(indexsendorbital2))*kind(indexsendorbital2)
  deallocate(indexsendorbital2, stat=istat)
  call memocc(istat, iall, 'indexsendorbital2', subname)


  allocate(indexsendorbital2(ndimpsi_f), stat=istat)
  call memocc(istat, indexsendorbital2, 'indexsendorbital2', subname)
  indexsendorbital2=indexsendorbital_f
  do i=1,ndimpsi_f
      ind=isendbuf_f(i)
      indexsendorbital_f(ind)=indexsendorbital2(i)
  end do
  ! Inverse of isendbuf
!!t1=mpi_wtime()
  call get_reverse_indices(ndimpsi_f, isendbuf_f, irecvbuf_f)
  !!do i=1,ndimpsi_f
  !!    do j=1,ndimpsi_f
  !!        if(isendbuf_f(j)==i) then
  !!            irecvbuf_f(i)=j
  !!        end if
  !!    end do
  !!end do
!!t2=mpi_wtime()
!!t_reverse=t_reverse+t2-t1
  iall=-product(shape(indexsendorbital2))*kind(indexsendorbital2)
  deallocate(indexsendorbital2, stat=istat)
  call memocc(istat, iall, 'indexsendorbital2', subname)




  if(nproc>1) then
      ! Communicate indexsendbuf
      call mpi_alltoallv(indexsendbuf_c, nsendcounts_c, nsenddspls_c, mpi_integer, indexrecvbuf_c, &
           nrecvcounts_c, nrecvdspls_c, mpi_integer, mpi_comm_world, ierr)
      ! Communicate indexsendorbitals
      call mpi_alltoallv(indexsendorbital_c, nsendcounts_c, nsenddspls_c, mpi_integer, indexrecvorbital_c, &
           nrecvcounts_c, nrecvdspls_c, mpi_integer, mpi_comm_world, ierr)

      ! Communicate indexsendbuf
      call mpi_alltoallv(indexsendbuf_f, nsendcounts_f, nsenddspls_f, mpi_integer, indexrecvbuf_f, &
           nrecvcounts_f, nrecvdspls_f, mpi_integer, mpi_comm_world, ierr)
      ! Communicate indexsendorbitals
      call mpi_alltoallv(indexsendorbital_f, nsendcounts_f, nsenddspls_f, mpi_integer, indexrecvorbital_f, &
           nrecvcounts_f, nrecvdspls_f, mpi_integer, mpi_comm_world, ierr)
   else
       indexrecvbuf_c=indexsendbuf_c
       indexrecvorbital_c=indexsendorbital_c
       indexrecvbuf_f=indexsendbuf_f
       indexrecvorbital_f=indexsendorbital_f
   end if



  !!call get_gridpoint_start(iproc, nproc, norb, glr, llr, nrecvcounts, indexrecvbuf, weight, gridpoint_start)
  call get_gridpoint_start(iproc, nproc, lzd, nrecvcounts_c, nrecvcounts_f, indexrecvbuf_c, indexrecvbuf_f, &
            weight_c, weight_f, gridpoint_start_c, gridpoint_start_f)



  if(maxval(gridpoint_start_c)>sum(nrecvcounts_c)) stop '1: maxval(gridpoint_start_c)>sum(nrecvcounts_c)'
  if(maxval(gridpoint_start_f)>sum(nrecvcounts_f)) stop '1: maxval(gridpoint_start_f)>sum(nrecvcounts_f)'
  ! Rearrange the communicated data
  do i=1,sum(nrecvcounts_c)
      ii=indexrecvbuf_c(i)
      jj=ii-1
      i3=jj/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
      jj=jj-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
      i2=jj/(lzd%glr%d%n1+1)
      i1=jj-i2*(lzd%glr%d%n1+1)
      !!if(weight_c(i1,i2,i3)==0.d0) stop 'coarse: weight is zero!'
      !!if(weight_c(i1,i2,i3)>0.d0) then
      !!    if(gridpoint_start_c(ii)==0) then
      !!        write(*,'(a,5i8)') 'DEBUG: iproc, jj, i1, i2, i3', iproc, jj, i1, i2, i3
      !!        stop 'coarse: weight>0, but gridpoint_start(ii)==0'
      !!    end if
      !!end if

      ind=gridpoint_start_c(ii)
      !!if(ind==0) stop 'ind is zero!'
      iextract_c(i)=ind
      gridpoint_start_c(ii)=gridpoint_start_c(ii)+1  
  end do
  !!write(*,'(a,2i12)') 'sum(iextract_c), nint(weightp_c*(weightp_c+1.d0)*.5d0)', sum(iextract_c), nint(weightp_c*(weightp_c+1.d0)*.5d0)
  !!if(sum(iextract_c)/=nint(weightp_c*(weightp_c+1.d0)*.5d0)) stop 'sum(iextract_c)/=nint(weightp_c*(weightp_c+1.d0)*.5d0)'
  if(maxval(iextract_c)>sum(nrecvcounts_c)) stop 'maxval(iextract_c)>sum(nrecvcounts_c)'
  if(minval(iextract_c)<1) stop 'minval(iextract_c)<1'

  ! Rearrange the communicated data
  do i=1,sum(nrecvcounts_f)
      ii=indexrecvbuf_f(i)
      jj=ii-1
      i3=jj/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
      jj=jj-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
      i2=jj/(lzd%glr%d%n1+1)
      i1=jj-i2*(lzd%glr%d%n1+1)
      !!if(weight_f(i1,i2,i3)==0.d0) stop 'fine: weight is zero!'
      !!if(weight_f(i1,i2,i3)>0.d0) then
      !!    if(gridpoint_start_f(ii)==0) then
      !!        write(*,'(a,5i8)') 'DEBUG: iproc, jj, i1, i2, i3', iproc, jj, i1, i2, i3
      !!        stop 'fine: weight>0, but gridpoint_start(ii)==0'
      !!    end if
      !!end if

      ind=gridpoint_start_f(ii)
      !!if(ind==0) stop 'ind is zero!'
      iextract_f(i)=ind
      gridpoint_start_f(ii)=gridpoint_start_f(ii)+1  
  end do
  if(sum(iextract_f)/=nint(weightp_f*(weightp_f+1.d0)*.5d0)) stop 'sum(iextract_f)/=nint(weightp_f*(weightp_f+1.d0)*.5d0)'
  if(maxval(iextract_f)>sum(nrecvcounts_f)) stop 'maxval(iextract_f)>sum(nrecvcounts_f)'
  if(minval(iextract_f)<1) stop 'minval(iextract_f)<1'




  ! Get the array to transfrom back the data
!!t1=mpi_wtime()
  call get_reverse_indices(sum(nrecvcounts_c), iextract_c, iexpand_c)
  !!do i=1,sum(nrecvcounts_c)
  !!    do j=1,sum(nrecvcounts_c)
  !!        if(iextract_c(j)==i) then
  !!            iexpand_c(i)=j
  !!        end if
  !!    end do
  !!end do
!!t2=mpi_wtime()
!!t_reverse=t_reverse+t2-t1

!!t1=mpi_wtime()
  call get_reverse_indices(sum(nrecvcounts_f), iextract_f, iexpand_f)
  !!do i=1,sum(nrecvcounts_f)
  !!    do j=1,sum(nrecvcounts_f)
  !!        if(iextract_f(j)==i) then
  !!            iexpand_f(i)=j
  !!        end if
  !!    end do
  !!end do
!!t2=mpi_wtime()
!!t_reverse=t_reverse+t2-t1
  



  allocate(indexrecvorbital2(sum(nrecvcounts_c)), stat=istat)
  call memocc(istat, indexrecvorbital2, 'indexrecvorbital2', subname)
  indexrecvorbital2=indexrecvorbital_c
  do i=1,sum(nrecvcounts_c)
      ind=iextract_c(i)
      indexrecvorbital_c(ind)=indexrecvorbital2(i)
  end do
  iall=-product(shape(indexrecvorbital2))*kind(indexrecvorbital2)
  deallocate(indexrecvorbital2, stat=istat)
  call memocc(istat, iall, 'indexrecvorbital2', subname)

  allocate(indexrecvorbital2(sum(nrecvcounts_f)), stat=istat)
  call memocc(istat, indexrecvorbital2, 'indexrecvorbital2', subname)
  indexrecvorbital2=indexrecvorbital_f
  do i=1,sum(nrecvcounts_f)
      ind=iextract_f(i)
      indexrecvorbital_f(ind)=indexrecvorbital2(i)
  end do
  iall=-product(shape(indexrecvorbital2))*kind(indexrecvorbital2)
  deallocate(indexrecvorbital2, stat=istat)
  call memocc(istat, iall, 'indexrecvorbital2', subname)


  if(minval(indexrecvorbital_c)<1) stop 'minval(indexrecvorbital_c)<1'
  if(maxval(indexrecvorbital_c)>orbs%norb) stop 'maxval(indexrecvorbital_c)>orbs%norb'
  if(minval(indexrecvorbital_f)<1) stop 'minval(indexrecvorbital_f)<1'
  if(maxval(indexrecvorbital_f)>orbs%norb) stop 'maxval(indexrecvorbital_f)>orbs%norb'



  iall=-product(shape(indexsendorbital_c))*kind(indexsendorbital_c)
  deallocate(indexsendorbital_c, stat=istat)
  call memocc(istat, iall, 'indexsendorbital_c', subname)
  iall=-product(shape(indexsendbuf_c))*kind(indexsendbuf_c)
  deallocate(indexsendbuf_c, stat=istat)
  call memocc(istat, iall, 'indexsendbuf_c', subname)
  iall=-product(shape(indexrecvbuf_c))*kind(indexrecvbuf_c)
  deallocate(indexrecvbuf_c, stat=istat)
  call memocc(istat, iall, 'indexrecvbuf_c', subname)

  iall=-product(shape(indexsendorbital_f))*kind(indexsendorbital_f)
  deallocate(indexsendorbital_f, stat=istat)
  call memocc(istat, iall, 'indexsendorbital_f', subname)
  iall=-product(shape(indexsendbuf_f))*kind(indexsendbuf_f)
  deallocate(indexsendbuf_f, stat=istat)
  call memocc(istat, iall, 'indexsendbuf_f', subname)
  iall=-product(shape(indexrecvbuf_f))*kind(indexrecvbuf_f)
  deallocate(indexrecvbuf_f, stat=istat)
  call memocc(istat, iall, 'indexrecvbuf_f', subname)

  iall=-product(shape(weight_c))*kind(weight_c)
  deallocate(weight_c, stat=istat)
  call memocc(istat, iall, 'weight_c', subname)
  iall=-product(shape(weight_f))*kind(weight_f)
  deallocate(weight_f, stat=istat)
  call memocc(istat, iall, 'weight_f', subname)

  iall=-product(shape(gridpoint_start_c))*kind(gridpoint_start_c)
  deallocate(gridpoint_start_c, stat=istat)
  call memocc(istat, iall, 'gridpoint_start_c', subname)
  iall=-product(shape(gridpoint_start_f))*kind(gridpoint_start_f)
  deallocate(gridpoint_start_f, stat=istat)
  call memocc(istat, iall, 'gridpoint_start_f', subname)

  iall=-product(shape(nsend))*kind(nsend)
  deallocate(nsend, stat=istat)
  call memocc(istat, iall, 'nsend', subname)

!!t2tot=mpi_wtime()
!!write(*,'(a,i6,es15.5)') 'in sub get_switch_indices: iproc, time reverse',iproc, t_reverse
!!write(*,'(a,i6,es15.5)') 'in sub get_switch_indices: iproc, total time',iproc, t2tot-t1tot

end subroutine get_switch_indices




subroutine get_gridpoint_start(iproc, nproc, lzd, nrecvcounts_c, nrecvcounts_f, indexrecvbuf_c, indexrecvbuf_f, &
            weight_c, weight_f, gridpoint_start_c, gridpoint_start_f)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzd
integer,dimension(0:nproc-1),intent(in):: nrecvcounts_c, nrecvcounts_f
integer,dimension(sum(nrecvcounts_c)),intent(in):: indexrecvbuf_c
integer,dimension(sum(nrecvcounts_f)),intent(in):: indexrecvbuf_f
real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(out):: weight_c, weight_f
integer,dimension((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)),intent(out):: gridpoint_start_c, gridpoint_start_f

! Local variables
integer:: i, ii, jj, i1, i2, i3


  weight_c=0.d0
  do i=1,sum(nrecvcounts_c)
      ii=indexrecvbuf_c(i)
      !!write(650+iproc,*) i, ii
      jj=ii-1
      i3=jj/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
      jj=jj-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
      i2=jj/(lzd%glr%d%n1+1)
      i1=jj-i2*(lzd%glr%d%n1+1)
      weight_c(i1,i2,i3)=weight_c(i1,i2,i3)+1.d0
  end do

  !!write(*,*) 'in get_gridpoint_start: maxval(weight_c)', maxval(weight_c)

  ii=1
  i=0
  gridpoint_start_c=0
  do i3=0,lzd%glr%d%n3
      do i2=0,lzd%glr%d%n2
          do i1=0,lzd%glr%d%n1
              i=i+1
              if(weight_c(i1,i2,i3)>0.d0) then
                  gridpoint_start_c(i)=ii
                  ii=ii+nint(weight_c(i1,i2,i3))
              end if
          end do
      end do
  end do

  !! CHECK
  i=0
  do i3=0,lzd%glr%d%n3
      do i2=0,lzd%glr%d%n2
          do i1=0,lzd%glr%d%n1
              i=i+1
              if(weight_c(i1,i2,i3)>0.d0) then
                  if(gridpoint_start_c(i)==0) stop 'FIRST CHECK: ERROR'
              end if
          end do
      end do
  end do



  ! fine part
  weight_f=0.d0
  do i=1,sum(nrecvcounts_f)
      ii=indexrecvbuf_f(i)
      jj=ii-1
      i3=jj/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
      jj=jj-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
      i2=jj/(lzd%glr%d%n1+1)
      i1=jj-i2*(lzd%glr%d%n1+1)
      weight_f(i1,i2,i3)=weight_f(i1,i2,i3)+1.d0
  end do

  ii=1
  i=0
  gridpoint_start_f=0
  do i3=0,lzd%glr%d%n3
      do i2=0,lzd%glr%d%n2
          do i1=0,lzd%glr%d%n1
              i=i+1
              if(weight_f(i1,i2,i3)>0.d0) then
                  gridpoint_start_f(i)=ii
                  ii=ii+nint(weight_f(i1,i2,i3))
              end if
          end do
      end do
  end do

  !! CHECK
  i=0
  do i3=0,lzd%glr%d%n3
      do i2=0,lzd%glr%d%n2
          do i1=0,lzd%glr%d%n1
              i=i+1
              if(weight_f(i1,i2,i3)>0.d0) then
                  if(gridpoint_start_f(i)==0) stop 'FIRST CHECK: ERROR'
              end if
          end do
      end do
  end do


end subroutine get_gridpoint_start




subroutine check_gridpoint(nseg, n1, n2, noffset1, noffset2, noffset3, keyg, itarget1, itarget2, itarget3, iseg_start, found)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: nseg, n1, n2, noffset1, noffset2, noffset3, itarget1, itarget2, itarget3
integer,dimension(2,nseg),intent(in):: keyg
integer,intent(inout):: iseg_start
logical,intent(out):: found

! Local variables
integer:: j0, j1, ii, i1, i2, i3, i0, ii1, ii2, ii3, iseg, i
logical:: equal_possible, larger_possible, smaller_possible
integer:: iproc

call mpi_comm_rank(mpi_comm_world, iproc, i)

  found=.false.
  !!write(300+iproc,*) '---start---'
  loop_segments: do iseg=iseg_start,nseg
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     ii2=i2+noffset2
     ii3=i3+noffset3
     equal_possible = (ii2==itarget2 .and. ii3==itarget3)
     larger_possible = (ii2>itarget2 .and. ii3>itarget3)
     smaller_possible = (ii2<itarget2 .and. ii3<itarget3)
     ! check whether quick exit is possible since there is no chance to find the point anymore...
     if(ii3>itarget3) then
            exit loop_segments
     end if
     if(ii3>=itarget3 .and. ii2>itarget2) then
         exit loop_segments
     end if
     larger_possible = (ii3>=itarget3 .and. ii2>=itarget2)

     do i=i0,i1
        ii1=i+noffset1
        if(equal_possible .and. ii1==itarget1) then
            found=.true.
            ! no need to search in smaller segments from now on, since the itargets will never decrease any more...
            iseg_start=iseg
            exit loop_segments
        end if
        if(larger_possible .and. ii1>itarget1) then
            ! there is no chance to find the point anymore...
            exit loop_segments
        end if
        if(smaller_possible .and. ii1<itarget1) then
            ! no need to search in these segments from now on, since the itargets will never decrease any more...
            iseg_start=iseg
        end if
     end do
  end do loop_segments
  !!write(300+iproc,*) 'new iseg_start:',iseg_start
  !!write(300+iproc,*) '--- end ---'


end subroutine check_gridpoint




subroutine get_index_in_global(lr, itarget1, itarget2, itarget3, region, ind)
use module_base
use module_types
implicit none

! Calling arguments
type(locreg_descriptors),intent(in):: lr
integer,intent(in):: itarget1, itarget2, itarget3
character(len=1),intent(in):: region
integer,intent(out):: ind

! Local variables
integer:: iitot, iseg, j0, j1, ii, i1, i2, i3, i0, i, istart, iend, ii1, ii2, ii3


 if(region=='c') then
    iitot=0
    loop_segments_c: do iseg=1,lr%wfd%nseg_c
       j0=lr%wfd%keygloc(1,iseg)
       j1=lr%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lr%d%n1+1)*(lr%d%n2+1))
       ii=ii-i3*(lr%d%n1+1)*(lr%d%n2+1)
       i2=ii/(lr%d%n1+1)
       i0=ii-i2*(lr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
          iitot=iitot+1
          ii1=i+lr%ns1
          ii2=i2+lr%ns2
          ii3=i3+lr%ns3
          if(ii1==itarget1 .and. ii2==itarget2 .and. ii3==itarget3) then
              ind=iitot
              exit loop_segments_c
          end if
       end do
    end do loop_segments_c

  else if(region=='f') then

    iitot=0
    istart=lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)
    iend=istart+lr%wfd%nseg_f-1
    loop_segments_f: do iseg=istart,iend
       j0=lr%wfd%keygloc(1,iseg)
       j1=lr%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lr%d%n1+1)*(lr%d%n2+1))
       ii=ii-i3*(lr%d%n1+1)*(lr%d%n2+1)
       i2=ii/(lr%d%n1+1)
       i0=ii-i2*(lr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
          ii1=i+lr%ns1
          ii2=i2+lr%ns2
          ii3=i3+lr%ns3
          iitot=iitot+1
          if(ii1==itarget1 .and. ii2==itarget2 .and. ii3==itarget3) then
              ind=iitot
              exit loop_segments_f
          end if
       end do
    end do loop_segments_f

else
    stop 'wrong region'
end if



end subroutine get_index_in_global





subroutine get_index_in_global2(lr, index_in_global_c, index_in_global_f)
use module_base
use module_types
implicit none

! Calling arguments
type(locreg_descriptors),intent(in):: lr
integer,dimension(0:lr%d%n1,0:lr%d%n2,0:lr%d%n3),intent(out):: index_in_global_c, index_in_global_f

! Local variables
integer:: iitot, iseg, j0, j1, ii, i1, i2, i3, i0, i, istart, iend


    iitot=0
    do iseg=1,lr%wfd%nseg_c
       j0=lr%wfd%keygloc(1,iseg)
       j1=lr%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lr%d%n1+1)*(lr%d%n2+1))
       ii=ii-i3*(lr%d%n1+1)*(lr%d%n2+1)
       i2=ii/(lr%d%n1+1)
       i0=ii-i2*(lr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
          iitot=iitot+1
          index_in_global_c(i,i2,i3)=iitot
       end do
    end do 


    iitot=0
    istart=lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)
    iend=istart+lr%wfd%nseg_f-1
    do iseg=istart,iend
       j0=lr%wfd%keygloc(1,iseg)
       j1=lr%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((lr%d%n1+1)*(lr%d%n2+1))
       ii=ii-i3*(lr%d%n1+1)*(lr%d%n2+1)
       i2=ii/(lr%d%n1+1)
       i0=ii-i2*(lr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
          iitot=iitot+1
          index_in_global_f(i,i2,i3)=iitot
       end do
    end do



end subroutine get_index_in_global2





subroutine transpose_switch_psi(orbs, lzd, collcom, psi, psiwork_c, psiwork_f)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(orbitals_Data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(collective_comms),intent(in):: collcom
  real(8),dimension(orbs%npsidim_orbs),intent(in):: psi
  real(8),dimension(collcom%ndimpsi_c),intent(out):: psiwork_c
  real(8),dimension(7*collcom%ndimpsi_f),intent(out):: psiwork_f
  
  ! Local variables
  integer:: i_tot, i_c, i_f, iorb, iiorb, ilr, i, ind, istat, iall
  real(8),dimension(:),allocatable:: psi_c, psi_f
  character(len=*),parameter:: subname='transpose_switch_psi'
  
  
  ! split up psi into coarse and fine part
  allocate(psi_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, psi_c, 'psi_c', subname)
  allocate(psi_f(7*collcom%ndimpsi_f), stat=istat)
  call memocc(istat, psi_f, 'psi_f', subname)
  i_tot=0
  i_c=0
  i_f=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      do i=1,lzd%llr(ilr)%wfd%nvctr_c
          i_c=i_c+1
          i_tot=i_tot+1
          psi_c(i_c)=psi(i_tot)
      end do
      do i=1,7*lzd%llr(ilr)%wfd%nvctr_f
          i_f=i_f+1
          i_tot=i_tot+1
          psi_f(i_f)=psi(i_tot)
      end do
  end do
  
  
  ! coarse part
  do i=1,collcom%ndimpsi_c
      ind=collcom%isendbuf_c(i)
      psiwork_c(ind)=psi_c(i)
  end do
  
  ! fine part
  do i=1,collcom%ndimpsi_f
      ind=collcom%isendbuf_f(i)
      psiwork_f(7*ind-6)=psi_f(7*i-6)
      psiwork_f(7*ind-5)=psi_f(7*i-5)
      psiwork_f(7*ind-4)=psi_f(7*i-4)
      psiwork_f(7*ind-3)=psi_f(7*i-3)
      psiwork_f(7*ind-2)=psi_f(7*i-2)
      psiwork_f(7*ind-1)=psi_f(7*i-1)
      psiwork_f(7*ind-0)=psi_f(7*i-0)
  end do

  iall=-product(shape(psi_c))*kind(psi_c)
  deallocate(psi_c, stat=istat)
  call memocc(istat, iall, 'psi_c', subname)
  iall=-product(shape(psi_f))*kind(psi_f)
  deallocate(psi_f, stat=istat)
  call memocc(istat, iall, 'psi_f', subname)
  
end subroutine transpose_switch_psi





subroutine transpose_communicate_psi(collcom, psiwork_c, psiwork_f, psitwork_c, psitwork_f)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(collective_comms),intent(in):: collcom
  real(8),dimension(collcom%ndimpsi_c),intent(in):: psiwork_c
  real(8),dimension(7*collcom%ndimpsi_f),intent(in):: psiwork_f
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(out):: psitwork_c
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(out):: psitwork_f
  
  ! Local variables
  integer:: ierr
  
  ! coarse part
  call mpi_alltoallv(psiwork_c, collcom%nsendcounts_c, collcom%nsenddspls_c, mpi_double_precision, psitwork_c, &
       collcom%nrecvcounts_c, collcom%nrecvdspls_c, mpi_double_precision, mpi_comm_world, ierr)
  
  ! fine part
  call mpi_alltoallv(psiwork_f, 7*collcom%nsendcounts_f, 7*collcom%nsenddspls_f, mpi_double_precision, psitwork_f, &
       7*collcom%nrecvcounts_f, 7*collcom%nrecvdspls_f, mpi_double_precision, mpi_comm_world, ierr)

end subroutine transpose_communicate_psi



subroutine transpose_unswitch_psit(collcom, psitwork_c, psitwork_f, psit_c, psit_f)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(collective_comms),intent(in):: collcom
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(in):: psitwork_c
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(in):: psitwork_f
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(out):: psit_c
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(out):: psit_f
  
  ! Local variables
  integer:: i, ind
  
  ! coarse part
  do i=1,sum(collcom%nrecvcounts_c)
      ind=collcom%iextract_c(i)
      psit_c(ind)=psitwork_c(i)
  end do

  ! fine part
  do i=1,sum(collcom%nrecvcounts_f)
      ind=collcom%iextract_f(i)
      psit_f(7*ind-6)=psitwork_f(7*i-6)
      psit_f(7*ind-5)=psitwork_f(7*i-5)
      psit_f(7*ind-4)=psitwork_f(7*i-4)
      psit_f(7*ind-3)=psitwork_f(7*i-3)
      psit_f(7*ind-2)=psitwork_f(7*i-2)
      psit_f(7*ind-1)=psitwork_f(7*i-1)
      psit_f(7*ind-0)=psitwork_f(7*i-0)
  end do

end subroutine transpose_unswitch_psit





subroutine transpose_switch_psit(collcom, psit_c, psit_f, psitwork_c, psitwork_f)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(collective_comms),intent(in):: collcom
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(in):: psit_c
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(in):: psit_f
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(out):: psitwork_c
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(out):: psitwork_f
  
  ! Local variables
  integer:: i, ind

 
  ! coarse part
  do i=1,sum(collcom%nrecvcounts_c)
      ind=collcom%iexpand_c(i)
      psitwork_c(ind)=psit_c(i)
  end do

  ! fine part
  do i=1,sum(collcom%nrecvcounts_f)
      ind=collcom%iexpand_f(i)
      psitwork_f(7*ind-6)=psit_f(7*i-6)
      psitwork_f(7*ind-5)=psit_f(7*i-5)
      psitwork_f(7*ind-4)=psit_f(7*i-4)
      psitwork_f(7*ind-3)=psit_f(7*i-3)
      psitwork_f(7*ind-2)=psit_f(7*i-2)
      psitwork_f(7*ind-1)=psit_f(7*i-1)
      psitwork_f(7*ind-0)=psit_f(7*i-0)
  end do

end subroutine transpose_switch_psit


subroutine transpose_communicate_psit(collcom, psitwork_c, psitwork_f, psiwork_c, psiwork_f)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(collective_comms),intent(in):: collcom
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(in):: psitwork_c
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(in):: psitwork_f
  real(8),dimension(collcom%ndimpsi_c),intent(out):: psiwork_c
  real(8),dimension(7*collcom%ndimpsi_f),intent(out):: psiwork_f
  
  ! Local variables
  integer:: ierr

 ! coarse part
  call mpi_alltoallv(psitwork_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, mpi_double_precision, psiwork_c, &
       collcom%nsendcounts_c, collcom%nsenddspls_c, mpi_double_precision, mpi_comm_world, ierr)

 ! fine part
  call mpi_alltoallv(psitwork_f, 7*collcom%nrecvcounts_f, 7*collcom%nrecvdspls_f, mpi_double_precision, psiwork_f, &
       7*collcom%nsendcounts_f, 7*collcom%nsenddspls_f, mpi_double_precision, mpi_comm_world, ierr)

end subroutine transpose_communicate_psit



subroutine transpose_unswitch_psi(orbs, lzd, collcom, psiwork_c, psiwork_f, psi)
  use module_base
  use module_types
  implicit none
  
  ! Caling arguments
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(collective_comms),intent(in):: collcom
  real(8),dimension(collcom%ndimpsi_c),intent(in):: psiwork_c
  real(8),dimension(7*collcom%ndimpsi_f),intent(in):: psiwork_f
  real(8),dimension(orbs%npsidim_orbs),intent(out):: psi
  
  ! Local variables
  integer:: i, ind, iorb, iiorb, ilr, i_tot, i_c, i_f, istat, iall
  real(8),dimension(:),allocatable:: psi_c, psi_f
  character(len=*),parameter:: subname='transpose_unswitch_psi'
  
  
  allocate(psi_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, psi_c, 'psi_c', subname)
  allocate(psi_f(7*collcom%ndimpsi_f), stat=istat)
  call memocc(istat, psi_f, 'psi_f', subname)
  
  
    ! coarse part
    do i=1,collcom%ndimpsi_c
        ind=collcom%irecvbuf_c(i)
        psi_c(ind)=psiwork_c(i)
    end do
  
    ! fine part
    do i=1,collcom%ndimpsi_f
        ind=collcom%irecvbuf_f(i)
        psi_f(7*ind-6)=psiwork_f(7*i-6)
        psi_f(7*ind-5)=psiwork_f(7*i-5)
        psi_f(7*ind-4)=psiwork_f(7*i-4)
        psi_f(7*ind-3)=psiwork_f(7*i-3)
        psi_f(7*ind-2)=psiwork_f(7*i-2)
        psi_f(7*ind-1)=psiwork_f(7*i-1)
        psi_f(7*ind-0)=psiwork_f(7*i-0)
    end do
  
    ! glue together coarse and fine part
    i_tot=0
    i_c=0
    i_f=0
    do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        ilr=orbs%inwhichlocreg(iiorb)
        do i=1,lzd%llr(ilr)%wfd%nvctr_c
            i_c=i_c+1
            i_tot=i_tot+1
            psi(i_tot)=psi_c(i_c)
        end do
        do i=1,7*lzd%llr(ilr)%wfd%nvctr_f
            i_f=i_f+1
            i_tot=i_tot+1
            psi(i_tot)=psi_f(i_f)
        end do
    end do
  
  iall=-product(shape(psi_c))*kind(psi_c)
  deallocate(psi_c, stat=istat)
  call memocc(istat, iall, 'psi_c', subname)
  iall=-product(shape(psi_f))*kind(psi_f)
  deallocate(psi_f, stat=istat)
  call memocc(istat, iall, 'psi_f', subname)

end subroutine transpose_unswitch_psi




subroutine transpose_localized(iproc, nproc, orbs, lzd, collcom, psi, psit_c, psit_f)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(collective_comms),intent(in):: collcom
  real(8),dimension(orbs%npsidim_orbs),intent(in):: psi
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(out):: psit_c
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(out):: psit_f
  
  ! Local variables
  real(8),dimension(:),allocatable:: psiwork_c, psiwork_f, psitwork_c, psitwork_f
  integer:: istat, iall
  character(len=*),parameter:: subname='transpose_localized'
  
  allocate(psiwork_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, psiwork_c, 'psiwork_c', subname)
  allocate(psiwork_f(7*collcom%ndimpsi_f), stat=istat)
  call memocc(istat, psiwork_f, 'psiwork_f', subname)
  allocate(psitwork_c(sum(collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, psitwork_c, 'psitwork_c', subname)
  allocate(psitwork_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, psitwork_f, 'psitwork_f', subname)
  
  call transpose_switch_psi(orbs, lzd, collcom, psi, psiwork_c, psiwork_f)
  if(nproc>1) then
      call transpose_communicate_psi(collcom, psiwork_c, psiwork_f, psitwork_c, psitwork_f)
  else
      psitwork_c=psiwork_c
      psitwork_f=psiwork_f
  end if
  call transpose_unswitch_psit(collcom, psitwork_c, psitwork_f, psit_c, psit_f)
  
  iall=-product(shape(psiwork_c))*kind(psiwork_c)
  deallocate(psiwork_c, stat=istat)
  call memocc(istat, iall, 'psiwork_c', subname)
  iall=-product(shape(psiwork_f))*kind(psiwork_f)
  deallocate(psiwork_f, stat=istat)
  call memocc(istat, iall, 'psiwork_f', subname)
  iall=-product(shape(psitwork_c))*kind(psitwork_c)
  deallocate(psitwork_c, stat=istat)
  call memocc(istat, iall, 'psitwork_c', subname)
  iall=-product(shape(psitwork_f))*kind(psitwork_f)
  deallocate(psitwork_f, stat=istat)
  call memocc(istat, iall, 'psitwork_f', subname)
  
end subroutine transpose_localized



subroutine untranspose_localized(iproc, nproc, orbs, lzd, collcom, psit_c, psit_f, psi)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(collective_comms),intent(in):: collcom
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(in):: psit_c
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(in):: psit_f
  real(8),dimension(orbs%npsidim_orbs),intent(out):: psi
  
  ! Local variables
  real(8),dimension(:),allocatable:: psiwork_c, psiwork_f, psitwork_c, psitwork_f
  integer:: istat, iall
  character(len=*),parameter:: subname='untranspose_localized'
  
  allocate(psiwork_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, psiwork_c, 'psiwork_c', subname)
  allocate(psiwork_f(7*collcom%ndimpsi_f), stat=istat)
  call memocc(istat, psiwork_f, 'psiwork_f', subname)
  allocate(psitwork_c(sum(collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, psitwork_c, 'psitwork_c', subname)
  allocate(psitwork_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, psitwork_f, 'psitwork_f', subname)

  call transpose_switch_psit(collcom, psit_c, psit_f, psitwork_c, psitwork_f)
  if(nproc>1) then
      call transpose_communicate_psit(collcom, psitwork_c, psitwork_f, psiwork_c, psiwork_f)
  else
      psiwork_c=psitwork_c
      psiwork_f=psitwork_f
  end if
  call transpose_unswitch_psi(orbs, lzd, collcom, psiwork_c, psiwork_f, psi)
  
  iall=-product(shape(psiwork_c))*kind(psiwork_c)
  deallocate(psiwork_c, stat=istat)
  call memocc(istat, iall, 'psiwork_c', subname)
  iall=-product(shape(psiwork_f))*kind(psiwork_f)
  deallocate(psiwork_f, stat=istat)
  call memocc(istat, iall, 'psiwork_f', subname)
  iall=-product(shape(psitwork_c))*kind(psitwork_c)
  deallocate(psitwork_c, stat=istat)
  call memocc(istat, iall, 'psitwork_c', subname)
  iall=-product(shape(psitwork_f))*kind(psitwork_f)
  deallocate(psitwork_f, stat=istat)
  call memocc(istat, iall, 'psitwork_f', subname)
  
end subroutine untranspose_localized


subroutine calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c1, psit_c2, psit_f1, psit_f2, ovrlp)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(collective_comms),intent(in):: collcom
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(in):: psit_c1, psit_c2
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(in):: psit_f1, psit_f2
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
  
  ! Local variables
  integer:: i0, ipt, ii, iiorb, j, jjorb, i, ierr

  ovrlp=0.d0

  i0=0
  do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          do j=1,ii
              jjorb=collcom%indexrecvorbital_c(i0+j)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_c1(i0+i)*psit_c2(i0+j)
          end do
      end do
      i0=i0+ii
  end do

  i0=0
  do ipt=1,collcom%nptsp_f 
      ii=collcom%norb_per_gridpoint_f(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_f(i0+i)
          do j=1,ii
              jjorb=collcom%indexrecvorbital_f(i0+j)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-6)*psit_f2(7*(i0+j)-6)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-5)*psit_f2(7*(i0+j)-5)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-4)*psit_f2(7*(i0+j)-4)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-3)*psit_f2(7*(i0+j)-3)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-2)*psit_f2(7*(i0+j)-2)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-1)*psit_f2(7*(i0+j)-1)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-0)*psit_f2(7*(i0+j)-0)
          end do
      end do
      i0=i0+ii
  end do

  if(nproc>1) call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)

end subroutine calculate_overlap_transposed



subroutine build_linear_combination_transposed(norb, matrix, collcom, psitwork_c, psitwork_f, reset, psit_c, psit_f)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: norb
  real(8),dimension(norb,norb),intent(in):: matrix
  type(collective_comms),intent(in):: collcom
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(in):: psitwork_c
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(in):: psitwork_f
  logical,intent(in):: reset
  real(8),dimension(sum(collcom%nrecvcounts_c)),intent(out):: psit_c
  real(8),dimension(7*sum(collcom%nrecvcounts_f)),intent(out):: psit_f

  ! Local variables
  integer:: i0, ipt, ii, j, iiorb, jjorb, i

  if(reset) then
      psit_c=0.d0
      psit_f=0.d0
  end if

  i0=0
  do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          do j=1,ii
              jjorb=collcom%indexrecvorbital_c(i0+j)
              psit_c(i0+i)=psit_c(i0+i)+matrix(jjorb,iiorb)*psitwork_c(i0+j)
          end do
      end do
      i0=i0+ii
  end do

  i0=0
  do ipt=1,collcom%nptsp_f 
      ii=collcom%norb_per_gridpoint_f(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_f(i0+i)
          do j=1,ii
              jjorb=collcom%indexrecvorbital_f(i0+j)
              psit_f(7*(i0+i)-6) = psit_f(7*(i0+i)-6) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-6)
              psit_f(7*(i0+i)-5) = psit_f(7*(i0+i)-5) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-5)
              psit_f(7*(i0+i)-4) = psit_f(7*(i0+i)-4) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-4)
              psit_f(7*(i0+i)-3) = psit_f(7*(i0+i)-3) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-3)
              psit_f(7*(i0+i)-2) = psit_f(7*(i0+i)-2) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-2)
              psit_f(7*(i0+i)-1) = psit_f(7*(i0+i)-1) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-1)
              psit_f(7*(i0+i)-0) = psit_f(7*(i0+i)-0) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-0)
          end do
      end do
      i0=i0+ii
  end do

end subroutine build_linear_combination_transposed




subroutine check_grid_point_from_boxes(i1, i2, i3, lr, overlap_possible)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: i1, i2, i3
  type(locreg_descriptors),intent(in):: lr  
  logical,intent(out):: overlap_possible

  ! Local variables
  logical:: ovrlpx, ovrlpy, ovrlpz
  
  ovrlpx = (i1>=lr%ns1 .and. i1<=lr%ns1+lr%d%n1)
  ovrlpy = (i2>=lr%ns2 .and. i2<=lr%ns2+lr%d%n2)
  ovrlpz = (i3>=lr%ns3 .and. i3<=lr%ns3+lr%d%n3)
  if(ovrlpx .and. ovrlpy .and. ovrlpz) then
      overlap_possible=.true.
  else
      overlap_possible=.true.
  end if

end subroutine check_grid_point_from_boxes


subroutine get_reverse_indices(n, indices, reverse_indices)
  use module_base
  implicit none
  
  ! Calling arguments
  integer,intent(in):: n
  integer,dimension(n),intent(in):: indices
  integer,dimension(n),intent(out):: reverse_indices

  ! Local variables
  integer:: i, j

  do i=1,n
      j=indices(i)
      reverse_indices(j)=i
  end do
  !!do i=1,n
  !!    do j=1,n
  !!        if(indices(j)==i) then
  !!            reverse_indices(i)=j
  !!        end if
  !!    end do
  !!end do

end subroutine get_reverse_indices
