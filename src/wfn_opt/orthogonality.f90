!! FUNCTION
!!    Orthogonality routine, for all the orbitals
!!    Uses wavefunctions in their transposed form
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine orthogonalize(iproc,nproc,orbs,comms,wfd,psi,input)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(input_variables), intent(in) :: input
  real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(inout) :: psi
  !local variables
  character(len=*), parameter :: subname='orthogonalize'
  integer :: i_stat,i_all,ierr,info
  integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor
  integer, dimension(:,:), allocatable :: ndimovrlp
  real(wp), dimension(:), allocatable :: ovrlp
  integer,dimension(:),allocatable:: norbArr
  character(len=20):: category


  ! Determine wheter we have close shell (nspin=1) or spin polarized (nspin=2)
  if (orbs%norbd>0) then 
     nspin=2 
  else 
     nspin=1 
  end if

  ! ndimovrlp describes the shape of the overlap matrix.
  allocate(ndimovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)
  
  ! Allocate norbArr which contains the number of up and down orbitals.
  allocate(norbArr(nspin), stat=i_stat)
  call memocc(i_stat,norbArr,'norbArr',subname)
  do ispin=1,nspin
     if(ispin==1) norbArr(ispin)=orbs%norbu
     if(ispin==2) norbArr(ispin)=orbs%norbd
  end do

  ! Choose which orthogonalization method shall be used:
  ! methOrtho==0: Cholesky orthonormalization (i.e. a pseudo Gram-Schmidt)
  ! methOrtho==1: hybrid Gram-Schmidt/Cholesky orthonormalization
  ! methOrtho==2: Loewdin orthonormalization
  if(input%methOrtho==0) then
     category='Chol'
     call timing(iproc, trim(category)//'_comput', 'ON')

     call dimension_ovrlp(nspin,orbs,ndimovrlp)

     ! Allocate the overlap matrix
     allocate(ovrlp(ndimovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)

     !print *,'there',iproc

     ! Make a loop over npsin; calculate the overlap matrix (for up/down, resp.) and orthogonalize (again for up/down, resp.).
     do ispin=1,nspin
        call getOverlap(iproc,nproc,nspin,norbArr(ispin),orbs,comms,&
             psi(1),ndimovrlp,ovrlp,norbArr,1,ispin,category)
        call cholesky(iproc,nproc,norbArr(ispin),psi(1),nspinor,nspin,orbs,comms,&
             ndimovrlp,ovrlp(1),norbArr,1,ispin)
     end do

     ! Deallocate the arrays.
     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)

  else if(input%methOrtho==1) then
       category='GS/Chol'
       call timing(iproc, trim(category)//'_comput', 'ON')
       
       ! Make a hybrid Gram-Schmidt/Cholesky orthonormalization.
       call gsChol(iproc,nproc,psi(1),input,nspinor,orbs,nspin,ndimovrlp,norbArr,comms)
  else if(input%methOrtho==2) then
     category='Loewdin'
     call timing(iproc,trim(category)//'_comput','ON')

     call dimension_ovrlp(nspin,orbs,ndimovrlp)
     
     ! Allocate the overlap matrix
     allocate(ovrlp(ndimovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
          
     ! Make a loop over npsin; calculate the overlap matrix (for up/down,resp.) and orthogonalize (again for up/down,resp.).
     do ispin=1,nspin
        call getOverlap(iproc,nproc,nspin,norbArr(ispin),orbs,comms,psi(1),ndimovrlp,ovrlp,norbArr,1,ispin,category)
        call loewdin(iproc,nproc,norbArr(ispin),orbs%nspinor,1,ispin,orbs,comms,nspin,psi,ovrlp,ndimovrlp,norbArr)
     end do
     
     ! Deallocate the arrays.
     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)
          
  else
     if(iproc==0) write(*,'(a)') 'ERROR: invalid choice for methOrtho.'
     if(iproc==0) write(*,'(a)') "Change it in 'input.perf' to 0, 1 or 2!"
     stop
  end if

  ! Deallocate the remaining arrays.
  i_all=-product(shape(norbArr))*kind(norbArr)
  deallocate(norbArr, stat=i_stat)
  call memocc(i_stat,i_all,'norbArr',subname)

  i_all=-product(shape(ndimovrlp))*kind(ndimovrlp)
  deallocate(ndimovrlp, stat=i_stat)
  call memocc(i_stat,i_all,'ndimovrlp',subname)


  call timing(iproc,trim(category)//'_comput','OF')
  
END SUBROUTINE orthogonalize
!!***


!!****f* BigDFT/orthoconstraint
!! FUNCTION
!!   Orthogonality routine, for all the orbitals
!!   Uses wavefunctions in their transposed form
!! SOURCE
!!
subroutine orthoconstraint(iproc,nproc,orbs,comms,wfd,psi,hpsi,scprsum)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(in) :: psi
  real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(out) :: hpsi
  real(dp), intent(out) :: scprsum
  !local variables
  character(len=*), parameter :: subname='orthoconstraint'
  integer :: i_stat,i_all,ierr,iorb,ise,jorb
  integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor
  real(dp) :: occ,tt
  integer, dimension(:,:), allocatable :: ndimovrlp
  real(wp), dimension(:), allocatable :: alag

  !separate the orthogonalisation procedure for up and down orbitals 
  !and for different k-points
  call timing(iproc,'LagrM_comput  ','ON')

  !number of components of the overlap matrix for parallel case
  !calculate the dimension of the overlap matrix for each k-point
  if (orbs%norbd > 0) then
     nspin=2
  else
     nspin=1
  end if

  !number of components for the overlap matrix in wp-kind real numbers

  allocate(ndimovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)

  call dimension_ovrlp(nspin,orbs,ndimovrlp)

  allocate(alag(ndimovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
  call memocc(i_stat,alag,'alag',subname)

  !put to zero all the k-points which are not needed
  call razero(ndimovrlp(nspin,orbs%nkpts),alag)

  !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
  ispsi=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

     do ispin=1,nspin

        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        if (nvctrp == 0) cycle

        if(nspinor==1) then
           call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
                max(1,nvctrp),hpsi(ispsi),max(1,nvctrp),0.0_wp,&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb)
        else
           !this part should be recheck in the case of nspinor == 2
           call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
                max(1,ncomp*nvctrp), &
                hpsi(ispsi),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb)
        end if
        ispsi=ispsi+nvctrp*norb*nspinor
     end do
  end do

  if (nproc > 1) then
     call timing(iproc,'LagrM_comput  ','OF')
     call timing(iproc,'LagrM_commun  ','ON')
     call mpiallred(alag(1),ndimovrlp(nspin,orbs%nkpts),MPI_SUM,MPI_COMM_WORLD,ierr)
     !call MPI_ALLREDUCE (alag(1,2),alag(1,1),ndimovrlp(nspin,orbs%nkpts),&
     !mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     call timing(iproc,'LagrM_commun  ','OF')
     call timing(iproc,'LagrM_comput  ','ON')
  end if

  !now each processors knows all the overlap matrices for each k-point
  !even if it does not handle it.
  !this is somehow redundant but it is one way of reducing the number of communications
  !without defining group of processors

  !calculate the sum of the diagonal of the overlap matrix, for each k-point
  scprsum=0.0_dp
  !for each k-point calculate the gradient
  ispsi=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

     do ispin=1,nspin
        if (ispin==1) ise=0
        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        if (nvctrp == 0) cycle

!!$        !correct the orthogonality constraint if there are some orbitals which have zero occupation number
!!$        do iorb=1,norb
!!$           do jorb=iorb+1,norb
!!$              if (orbs%occup((ikpt-1)*orbs%norb+iorb+ise) /= 0.0_gp .and. &
!!$                   orbs%occup((ikpt-1)*orbs%norb+jorb+ise) == 0.0_gp) then
!!$                 alag(ndimovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs) = 0.0_wp
!!$                 alag(ndimovrlp(ispin,ikpt-1)+jorb+(iorb-1)*norbs) = 0.0_wp
!!$                 !if (iproc ==0) print *,'i,j',iorb,jorb,alag(ndimovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs)
!!$              end if
!!$           end do
!!$        end do

        !calculate the scprsum if the k-point is associated to this processor
        if (orbs%ikptproc(ikpt) == iproc) then
           if(nspinor == 1) then
              do iorb=1,norb
                 occ=real(orbs%kwgts(ikpt)*orbs%occup((ikpt-1)*orbs%norb+iorb+ise),dp)
                 scprsum=scprsum+&
                      occ*real(alag(ndimovrlp(ispin,ikpt-1)+iorb+(iorb-1)*norbs),dp)
              enddo
           else if (nspinor == 4 .or. nspinor == 2) then
              !not sure about the imaginary part of the diagonal
              do iorb=1,norb
                 occ=real(orbs%kwgts(ikpt)*orbs%occup((ikpt-1)*orbs%norb+iorb+ise),dp)
                 scprsum=scprsum+&
                      occ*real(alag(ndimovrlp(ispin,ikpt-1)+2*iorb-1+(iorb-1)*norbs),dp)
                 scprsum=scprsum+&
                      occ*real(alag(ndimovrlp(ispin,ikpt-1)+2*iorb+(iorb-1)*norbs),dp)
              enddo
           end if
        end if
        ise=norb

        if(nspinor==1 .and. nvctrp /= 0) then
           call gemm('N','N',nvctrp,norb,norb,-1.0_wp,psi(ispsi),max(1,nvctrp),&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
                hpsi(ispsi),max(1,nvctrp))
        else if (nvctrp /= 0) then
           call c_gemm('N','N',ncomp*nvctrp,norb,norb,(-1.0_wp,0.0_wp),psi(ispsi),max(1,ncomp*nvctrp),&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb,(1.0_wp,0.0_wp),hpsi(ispsi),max(1,ncomp*nvctrp))
        end if
        ispsi=ispsi+nvctrp*norb*nspinor
     end do
  end do

  if (nproc > 1) then
     tt=scprsum
     call mpiallred(scprsum,1,MPI_SUM,MPI_COMM_WORLD,ierr)
     !call MPI_ALLREDUCE(tt,scprsum,1,mpidtypd,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  i_all=-product(shape(alag))*kind(alag)
  deallocate(alag,stat=i_stat)
  call memocc(i_stat,i_all,'alag',subname)

  i_all=-product(shape(ndimovrlp))*kind(ndimovrlp)
  deallocate(ndimovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ndimovrlp',subname)

  call timing(iproc,'LagrM_comput  ','OF')

END SUBROUTINE orthoconstraint
!!***


!!****f* BigDFT/subspace_diagonalisation
!! FUNCTION
!!   Found the linear combination of the wavefunctions which diagonalises
!!   the overlap matrix
!! SOURCE
!!
subroutine subspace_diagonalisation(iproc,nproc,orbs,comms,psi,hpsi,evsum)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(inout) :: orbs !eval is updated
  type(communications_arrays), intent(in) :: comms
  real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(in) :: hpsi
  real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(out) :: psi
  real(wp), intent(out) :: evsum
  !local variables
  character(len=*), parameter :: subname='subspace_diagonalisation'
  integer :: i_stat,i_all,ierr,info,iorb,n_lp,n_rp,npsiw,isorb,ise
  integer :: istart,ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor
  real(wp) :: occ
  integer, dimension(:,:), allocatable :: ndimovrlp
  real(wp), dimension(:), allocatable :: work_lp,work_rp,psiw
  real(wp), dimension(:,:), allocatable :: hamks

  !separate the diagonalisation procedure for up and down orbitals 
  !and for different k-points

  !number of components of the overlap matrix for parallel case
  istart=2
  if (nproc == 1) istart=1

  !calculate the dimension of the overlap matrix for each k-point
  if (orbs%norbd > 0) then
     nspin=2
  else
     nspin=1
  end if

  !number of components for the overlap matrix in wp-kind real numbers
  allocate(ndimovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)

  call dimension_ovrlp(nspin,orbs,ndimovrlp)

  allocate(hamks(ndimovrlp(nspin,orbs%nkpts),istart+ndebug),stat=i_stat)
  call memocc(i_stat,hamks,'hamks',subname)

  !put to zero all the k-points which are not needed
  call razero(ndimovrlp(nspin,orbs%nkpts)*istart,hamks)

  !dimension of the work arrays
  n_lp=0
  n_rp=0
  npsiw=0

  !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
  ispsi=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

     do ispin=1,nspin

        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        if (nvctrp == 0) cycle

        if(nspinor==1) then
           call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),max(1,nvctrp),hpsi(ispsi),&
                max(1,nvctrp),0.0_wp,&
                hamks(ndimovrlp(ispin,ikpt-1)+1,istart),norb)
        else
           !this part should be recheck in the case of nspinor == 2
           call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
                max(1,ncomp*nvctrp), &
                hpsi(ispsi),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
                hamks(ndimovrlp(ispin,ikpt-1)+1,istart),norb)
        end if
        ispsi=ispsi+nvctrp*norb*nspinor

        !dimensions of the work arrays
        n_lp=max(4*norbs,1000,n_lp)
        n_rp=max(3*norb+1,n_rp)
        npsiw=max(nvctrp*orbs%norb*nspinor,npsiw)

     end do
  end do

  if (nproc > 1) then
     call MPI_ALLREDUCE (hamks(1,2),hamks(1,1),ndimovrlp(nspin,orbs%nkpts),&
          mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  !now each processors knows all the overlap matrices for each k-point
  !even if it does not handle it.
  !this is somehow redundant but it is one way of reducing the number of communications
  !without defining group of processors

  allocate(work_lp(n_lp*2+ndebug),stat=i_stat)
  call memocc(i_stat,work_lp,'work_lp',subname)
  allocate(work_rp(n_rp+ndebug),stat=i_stat)
  call memocc(i_stat,work_rp,'work_rp',subname)

  allocate(psiw(npsiw+ndebug),stat=i_stat)
  call memocc(i_stat,psiw,'psiw',subname)


  !for each k-point now reorthogonalise wavefunctions
  ispsi=1
  evsum=0.0_wp
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
     isorb=1
     do ispin=1,nspin

        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        if (nvctrp == 0) cycle

        if(nspinor==1) then

           !shift to be add for eval
           call syev('V','U',norb,hamks(ndimovrlp(ispin,ikpt-1)+1,1),norb,&
                orbs%eval(isorb+(ikpt-1)*orbs%norb),work_lp(1),n_lp,info)
           if (info /= 0) write(*,*) 'SYEV ERROR',info

        else

           call  heev('V','U',norb,hamks(ndimovrlp(ispin,ikpt-1)+1,1),norb,&
                orbs%eval(isorb+(ikpt-1)*orbs%norb),work_lp(1),n_lp,work_rp(1),info)
           if (info /= 0) write(*,*) 'HEEV ERROR',info

        end if

        !here we have to add evsum and the KS orbitals written in terms of linear algebra
        !evsum should be corrected like the scprsum above

        !calculate the evsum if the k-point is associated to this processor
        if (orbs%ikptproc(ikpt) == iproc) then
           if (ispin==1) ise=0
           do iorb=1,norb
              occ=real(orbs%kwgts(ikpt)*orbs%occup((ikpt-1)*orbs%norb+iorb+ise),dp)
              evsum=evsum+orbs%eval(isorb+iorb-1+(ikpt-1)*orbs%norb)*occ
           enddo
           ise=norb
        end if

        !!do iorb=1,norb
        !!   occ=real(orbs%kwgts(ikpt)*orbs%occup((ikpt-1)*orbs%norb+iorb),wp)
        !!   evsum=evsum+orbs%eval(isorb+iorb-1+(ikpt-1)*orbs%norb)*occ
        !!   !if (iproc.eq.0) write(*,'(1x,a,i0,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
        !!enddo

!!        ! Transform to KS orbitals
!!        ! dgemm can be used instead of daxpy
!!        if(nspinor==1) then
!!           do iorb=1,norb
!!              call razero(nvctrp,psitt(1,iorb))
!!              do jorb=1,norb
!!                 alpha=hamks(jorb,iorb,1)
!!                 call axpy(nvctrp,alpha,psit(1,jorb),1,psitt(1,iorb),1)
!!              enddo
!!           enddo
!!        else
!!           do iorb=1,norb
!!              call razero(nvctrp*nspinor,psitt(1,iorb))
!!              do jorb=1,norb
!!                 call c_axpy(ncomp*nvctrp,hamks(2*jorb-1,iorb,1),psit(1,jorb),1,psitt(1,iorb),1)
!!              enddo
!!           enddo
!!        end if

        !sample of dgemm
        if (nspinor == 1) then
           call gemm('N','N',nvctrp,norb,norb,1.0_wp,psi(ispsi),max(1,nvctrp),&
                hamks(ndimovrlp(ispin,ikpt-1)+1,1),norb,0.0_wp,psiw(1),max(1,nvctrp))
        else
           call c_gemm('N','N',ncomp*nvctrp,norb,norb,(1.0_wp,0.0_wp),&
                psi(ispsi),max(1,ncomp*nvctrp),hamks(ndimovrlp(ispin,ikpt-1)+1,1),norb,&
                (0.0_wp,0.0_wp),psiw(1),max(1,ncomp*nvctrp))
        end if

        call DCOPY(nvctrp*norb*nspinor,psiw(1),1,psi(ispsi),1)

        ispsi=ispsi+nvctrp*norb*nspinor
        isorb=isorb+norb
     end do
  end do

  if (nproc > 1) then
     !evsumtmp=evsum
     call mpiallred(evsum,1,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  i_all=-product(shape(psiw))*kind(psiw)
  deallocate(psiw,stat=i_stat)
  call memocc(i_stat,i_all,'psiw',subname)

  i_all=-product(shape(work_lp))*kind(work_lp)
  deallocate(work_lp,stat=i_stat)
  call memocc(i_stat,i_all,'work_lp',subname)
  i_all=-product(shape(work_rp))*kind(work_rp)
  deallocate(work_rp,stat=i_stat)
  call memocc(i_stat,i_all,'work_rp',subname)

  i_all=-product(shape(hamks))*kind(hamks)
  deallocate(hamks,stat=i_stat)
  call memocc(i_stat,i_all,'hamks',subname)

  i_all=-product(shape(ndimovrlp))*kind(ndimovrlp)
  deallocate(ndimovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ndimovrlp',subname)

END SUBROUTINE subspace_diagonalisation
!!***


!!****f* BigDFT/orthon_virt_occup
!! DESCRIPTION
!!   Makes sure all psivirt/gradients are othogonal to the occupied states psi.
!!   This routine is almost the same as orthoconstraint_p. Difference:
!!   hpsi(:,norb) -->  psivirt(:,nvirte) , therefore rectangular alag.
!! 
!! WARNING
!!   Orthogonality to spin polarized channels is achieved in two calls,
!! SOURCE
!!
subroutine orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi_occ,psi_virt,msg)
  use module_base
  use module_types
  implicit none
  logical, intent(in) :: msg
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs,orbsv
  type(communications_arrays), intent(in) :: comms,commsv
  real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(in) :: psi_occ
  real(wp), dimension(sum(commsv%nvctr_par(iproc,1:orbsv%nkptsp))*orbsv%nspinor*orbsv%norb), intent(out) :: psi_virt
  !local variables
  character(len=*), parameter :: subname='orthon_virt_occup'
  integer :: i_stat,i_all,ierr,ispsiv,iorb,jorb,isorb
  integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor
  integer :: norbv,norbsv,ncompv,nvctrpv,nspinorv
  real(wp) :: scprsum,tt
  integer, dimension(:,:), allocatable :: ndimovrlp
  real(wp), dimension(:), allocatable :: alag

  !separate the orthogonalisation procedure for up and down orbitals 
  !and for different k-points
  call timing(iproc,'LagrM_comput  ','ON')

  !calculate the dimension of the overlap matrix for each k-point
  if (orbs%norbd > 0) then
     nspin=2
  else
     nspin=1
  end if

  !number of components for the overlap matrix in wp-kind real numbers

  allocate(ndimovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)

  call dimension_ovrlp_virt(nspin,orbs,orbsv,ndimovrlp)

  allocate(alag(ndimovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
  call memocc(i_stat,alag,'alag',subname)

  !put to zero all the k-points which are not needed
  call razero(ndimovrlp(nspin,orbs%nkpts),alag)

  !differentiate between real and complex wavefunctions
  !Lower triangle of overlap matrix using BLAS
  !     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; lower triangle


  !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
  ispsi=1
  ispsiv=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

     do ispin=1,nspin

        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        call orbitals_and_components(iproc,ikptp,ispin,orbsv,commsv,&
             nvctrpv,norbv,norbsv,ncompv,nspinorv)
        !there checks ensure that the component distribution scheme of virtual and occupied states is the same
        if (nvctrpv /= nvctrp) stop 'nvctrp'
        if (ncompv /= ncomp) stop 'ncomp'
        if (nspinorv /= nspinor) stop 'nspinor'

        if (nvctrp == 0) cycle
        
        norbv=orbsv%norbu
        if (ispin==2) norbv=orbsv%norbd

        !print *,'nvctrp',iproc,commsv%nvctr_par(iproc,ikptp),nvctrp,ikpt,orbs%nkpts,orbsv%nkpts,norbs,norbv,orbsv%nkptsp,orbs%nkptsp

        !print *,'iproc,nvctrp,nspin,norb,ispsi,ndimovrlp',iproc,nvctrp,nspin,norb,ispsi,ndimovrlp(ispin,ikpt-1)

        if(nspinor==1) then
           call gemm('T','N',norb,norbv,nvctrp,1.0_wp,psi_occ(ispsi),max(1,nvctrp),&
                psi_virt(ispsiv),max(1,nvctrp),0.0_wp,&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb)
        else
           !this part should be recheck in the case of nspinor == 2
           call c_gemm('C','N',norb,norbv,ncomp*nvctrp,(1.0_wp,0.0_wp),psi_occ(ispsi),&
                max(1,ncomp*nvctrp), &
                psi_virt(ispsiv),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb)
        end if

        ispsi=ispsi+nvctrp*norb*nspinor
        ispsiv=ispsiv+nvctrp*norbv*nspinor
     end do
  end do

  if (nproc > 1) then
     call timing(iproc,'LagrM_comput  ','OF')
     call timing(iproc,'LagrM_commun  ','ON')
     call mpiallred(alag(1),ndimovrlp(nspin,orbs%nkpts),MPI_SUM,MPI_COMM_WORLD,ierr)
     !call MPI_ALLREDUCE (alag(1,2),alag(1,1),ndimovrlp(nspin,orbs%nkpts),&
     !     mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     call timing(iproc,'LagrM_commun  ','OF')
     call timing(iproc,'LagrM_comput  ','ON')
  end if

  !now each processors knows all the overlap matrices for each k-point
  !even if it does not handle it.
  !this is somehow redundant but it is one way of reducing the number of communications
  !without defining group of processors

  !for each k-point now reorthogonalise wavefunctions
  ispsi=1
  ispsiv=1
  isorb=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

     do ispin=1,nspin

        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        if (nvctrp == 0) cycle

        norbv=orbsv%norbu
        if (ispin==2) norbv=orbsv%norbd


        if (msg) then
           write(*,'(1x,a)')'scalar products are'
           write(*,'(1x,a)')'iocc  ivirt       value'!               zero if<1d-12'

           scprsum=0.0_wp

           do iorb=1,norb
              do jorb=1,norbv
                 tt=alag(isorb+iorb-1+(jorb-1)*norbs)
                 write(*,'(1x,2i3,1pe21.14)')iorb,jorb,tt
                 scprsum=scprsum+tt**2
                 !if(abs(tt)<1d-12)alag(iorb,jorb,1)=0d0
                 !if(msg)write(*,'(2(i3),7x,2(1pe21.14))')iorb,jorb,tt,alag(iorb,jorb,1)
              end do
           enddo
           scprsum=sqrt(scprsum/real(norb,wp)/real(norbv,wp))
           write(*,'(1x,a,1pe21.14)')'sqrt sum squares is',scprsum
           write(*,'(1x)')
        end if

        if(nspinor==1 .and. nvctrp /= 0) then
           call gemm('N','N',nvctrp,norbv,norb,-1.0_wp,psi_occ(ispsi),max(1,nvctrp),&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
                psi_virt(ispsiv),max(1,nvctrp))
        else if (nvctrp /= 0) then
           call c_gemm('N','N',ncomp*nvctrp,norbv,norb,(-1.0_wp,0.0_wp),psi_occ(ispsi),&
                max(1,ncomp*nvctrp),&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb,(1.0_wp,0.0_wp),psi_virt(ispsi),&
                max(1,ncomp*nvctrp))
        end if
        ispsi=ispsi+nvctrp*norb*nspinor
        ispsiv=ispsiv+nvctrp*norbv*nspinor
        isorb=isorb+norbs*norbv
     end do
  end do


  i_all=-product(shape(alag))*kind(alag)
  deallocate(alag,stat=i_stat)
  call memocc(i_stat,i_all,'alag',subname)

  i_all=-product(shape(ndimovrlp))*kind(ndimovrlp)
  deallocate(ndimovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ndimovrlp',subname)

  call timing(iproc,'LagrM_comput  ','OF')

END SUBROUTINE orthon_virt_occup
!!***

subroutine complex_components(nspinor,norb,norbs,ncomp)
  implicit none
  integer, intent(in) :: nspinor,norb
  integer, intent(out) :: norbs,ncomp

  if(nspinor == 1) then
     norbs=norb
     ncomp=1 !useless
  else if (nspinor == 2) then
     norbs=2*norb
     ncomp=1
  else if (nspinor == 4) then
     norbs=2*norb
     ncomp=2
  end if
  
END SUBROUTINE complex_components

subroutine orbitals_and_components(iproc,ikptp,ispin,orbs,comms,nvctrp,norb,norbs,ncomp,nspinor)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,ikptp,ispin
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  integer, intent(out) :: nvctrp,norb,norbs,ncomp,nspinor

  nvctrp=comms%nvctr_par(iproc,ikptp)
  norb=orbs%norbu
  nspinor=orbs%nspinor
  if (ispin==2) norb=orbs%norbd

  call complex_components(nspinor,norb,norbs,ncomp)

END SUBROUTINE orbitals_and_components

subroutine dimension_ovrlp(nspin,orbs,ndimovrlp)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspin
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(nspin,0:orbs%nkpts), intent(out) :: ndimovrlp
  !local variables
  integer :: norb,norbs,ncomp,ikpt

  ndimovrlp(1,0)=0
  if (nspin == 2) then
     norb=orbs%norbu

     !this is first k-point
     call complex_components(orbs%nspinor,norb,norbs,ncomp)

     ndimovrlp(2,0)=norbs*norb
  end if

  do ikpt=1,orbs%nkpts
     !this part should be enhanced for real k-points
     norb=orbs%norbu
     if (nspin == 2) norb = orbs%norbd
     !this is ikpt k-point
     call complex_components(orbs%nspinor,norb,norbs,ncomp)

     ndimovrlp(1,ikpt)=ndimovrlp(nspin,ikpt-1)+norbs*norb
     if (orbs%norbd > 0) then
        norb=orbs%norbu
        !this is ikpt+1
        call complex_components(orbs%nspinor,norb,norbs,ncomp)
        if (ikpt == orbs%nkpts) then
           ndimovrlp(2,ikpt)=ndimovrlp(1,ikpt)
        else
           ndimovrlp(2,ikpt)=ndimovrlp(1,ikpt)+norbs*norb
        end if
     end if
  end do

END SUBROUTINE dimension_ovrlp

subroutine dimension_ovrlp_virt(nspin,orbs,orbsv,ndimovrlp)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspin
  type(orbitals_data), intent(in) :: orbs,orbsv
  integer, dimension(nspin,0:orbs%nkpts), intent(out) :: ndimovrlp
  !local variables
  integer :: norb,norbs,ncomp,ikpt,norbv

  ndimovrlp(1,0)=0
  if (nspin == 2) then
     norb=orbs%norbu
     norbv=orbsv%norbu
     call complex_components(orbs%nspinor,norb,norbs,ncomp)

     ndimovrlp(2,0)=norbs*norbv
  end if

  do ikpt=1,orbs%nkpts
     !this part should be enhanced for real k-points
     norb=orbs%norbu
     norbv=orbsv%norbu 
     if (nspin == 2) then
        norb=orbs%norbd
        norbv=orbsv%norbd
     end if


     call complex_components(orbs%nspinor,norb,norbs,ncomp)

     ndimovrlp(1,ikpt)=ndimovrlp(nspin,ikpt-1)+norbs*norbv
     if (orbs%norbd > 0) then

        norb=orbs%norbu
        norbv=orbsv%norbu
        
        call complex_components(orbs%nspinor,norb,norbs,ncomp)

        if (ikpt == orbs%nkpts) then
           ndimovrlp(2,ikpt)=ndimovrlp(1,ikpt)
        else
           ndimovrlp(2,ikpt)=ndimovrlp(1,ikpt)+norbs*norbv
        end if
     end if
  end do

END SUBROUTINE dimension_ovrlp_virt


subroutine orthoconstraint_p(iproc,nproc,norb,occup,nvctrp,psit,hpsit,scprsum,nspinor)
  !Effect of orthogonality constraints on gradient 
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norb,nvctrp,nspinor
  real(gp), dimension(norb), intent(in) :: occup
  real(wp), dimension(nspinor*nvctrp,norb), intent(in) :: psit 
  real(dp), intent(out) :: scprsum
  real(wp), dimension(nspinor*nvctrp,norb), intent(out) :: hpsit
  !local variables
  character(len=*), parameter :: subname='orthoconstraint_p'
  integer :: i_stat,i_all,istart,iorb,ierr,norbs,ncomp
  real(dp) :: occ
  real(wp), dimension(:,:,:), allocatable :: alag

  call timing(iproc,'LagrM_comput  ','ON')

  istart=2
  if (nproc == 1) istart=1

  if(nspinor == 1) then
     norbs=norb
  else if(nspinor == 2) then
     norbs=2*norb
     ncomp=1
  else if (nspinor == 4) then
     norbs=2*norb
     ncomp=2
  end if

  allocate(alag(norbs,norb,istart+ndebug),stat=i_stat)
  call memocc(i_stat,alag,'alag',subname)

  !initialise if nvctrp=0
  if (nvctrp == 0) then
     call razero(norbs*norb*istart,alag)
  end if

  !     alag(jorb,iorb,istart)=+psit(k,jorb)*hpsit(k,iorb)
  if(nspinor==1) then
     call GEMM('T','N',norb,norb,nvctrp,1.0_wp,psit(1,1),max(1,nvctrp),hpsit(1,1),max(1,nvctrp),0.0_wp,&
          alag(1,1,istart),norb)
  else
     !this part should be recheck in the case of nspinor == 2
     call C_GEMM('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psit(1,1),max(1,ncomp*nvctrp), &
          hpsit(1,1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),alag(1,1,istart),norb)
  end if

  if (nproc > 1) then
     call timing(iproc,'LagrM_comput  ','OF')
     call timing(iproc,'LagrM_commun  ','ON')
     call MPI_ALLREDUCE(alag(1,1,2),alag(1,1,1),norbs*norb,&
          mpidtypd,MPI_SUM,MPI_COMM_WORLD,ierr)
     call timing(iproc,'LagrM_commun  ','OF')
     call timing(iproc,'LagrM_comput  ','ON')
  end if
!          if (iproc.eq.0) then
!          write(*,*) 'ALAG',iproc,norb,norbs
!          do iorb=1,norb
!          write(*,'(10(1x,1pe10.3))') (alag(jorb,iorb,1),jorb=1,norbs)
!          enddo
!          endif



  scprsum=0.0_dp
  if(nspinor == 1) then
     do iorb=1,norb
        occ=real(occup(iorb),dp)
        scprsum=scprsum+occ*real(alag(iorb,iorb,1),dp)
     enddo
  else if (nspinor == 4 .or. nspinor == 2) then
     !not sure about the imaginary part of the diagonal
    do iorb=1,norb
       occ=real(occup(iorb),dp)
       scprsum=scprsum+occ*real(alag(2*iorb-1,iorb,1),dp)
       scprsum=scprsum+occ*real(alag(2*iorb,iorb,1),dp)
     enddo
  end if


!  if(iproc==0) print *,'ortho_p',scprsum

  ! hpsit(k,iorb)=-psit(k,jorb)*alag(jorb,iorb,1)
  if(nspinor==1) then
     call GEMM('N','N',nvctrp,norb,norb,-1.0_wp,psit(1,1),max(1,nvctrp),alag(1,1,1),norb,1.0_wp,&
          hpsit(1,1),max(1,nvctrp))
  else
     call C_GEMM('N','N',ncomp*nvctrp,norb,norb,(-1.0_wp,0.0_wp),psit(1,1),max(1,ncomp*nvctrp),&
          alag(1,1,1),norb,(1.0_wp,0.0_wp),hpsit(1,1),max(1,ncomp*nvctrp))
  end if


  i_all=-product(shape(alag))*kind(alag)
  deallocate(alag,stat=i_stat)
  call memocc(i_stat,i_all,'alag',subname)



  call timing(iproc,'LagrM_comput  ','OF')

END SUBROUTINE orthoconstraint_p

subroutine orthon_p(iproc,nproc,norb,nvctrp,nvctr_tot,psit,nspinor)
  ! Gram-Schmidt orthogonalisation
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norb,nvctrp,nvctr_tot,nspinor
  real(wp), dimension(nspinor*nvctrp,norb), intent(inout) :: psit
  !local variables
  character(len=*), parameter :: subname='orthon_p'
  integer :: info,i_all,i_stat,nvctr_eff,ierr,istart,norbs,ncomp
  real(wp) :: tt,ttLOC
  real(wp), dimension(:,:,:), allocatable :: ovrlp
  integer :: volta

  call timing(iproc,'GramS_comput  ','ON')

  do volta=1,2

  if (norb == 1) then 

     !for the inhomogeneous distribution this should  be changed
     nvctr_eff=nvctrp!min(nvctr_tot-iproc*nvctrp,nvctrp)

     if (nvctr_eff > 0) then
     !parallel treatment of a run with only one orbital
        if(nspinor==1) then
           tt=nrm2(nvctr_eff,psit(1,1),1)
        else
           !print *,'for one orbital the norm of the spinor must be calculated'
           !stop
           tt=nrm2(nvctr_eff*nspinor,psit(1,1),1) !NOT CORRECT
        end if
        ttLOC=tt**2
     else
        ttLOC=0.0_wp
     end if
     
     if (nproc > 1) then
        call MPI_ALLREDUCE(ttLOC,tt,1,mpidtypd,MPI_SUM,MPI_COMM_WORLD,ierr)
     else
        tt=ttLOC
     end if

     tt=1.0_wp/sqrt(tt)
     if(nspinor==1) then 
        !correct normalisation
        call vscal(nvctr_eff,tt,psit(1,1),1)
     else
        !not correct, to be adjusted
        call vscal(nvctr_eff*nspinor,tt,psit(1,1),1)
     end if

  else

     istart=2
     if (nproc == 1) istart=1

     if(nspinor==1) then
        norbs=norb
     else if (nspinor ==2) then
        norbs=2*norb
        ncomp=1
     else if (nspinor ==4) then
        norbs=2*norb
        ncomp=2
     end if

     allocate(ovrlp(norbs,norb,istart+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)

     call razero(norbs*norb*istart,ovrlp)

     ! Upper triangle of overlap matrix using BLAS
     !     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; upper triangle
     if(nspinor==1) then
        call syrk('L','T',norb,nvctrp,1.0_wp,psit(1,1),max(1,nvctrp),0.0_wp,ovrlp(1,1,istart),norb)
     else
!!        ovrlp=0.0d0
!!        do iorb=1,norb
!!           do jorb=1,norb
!!              ttr=ddot(nvctrp*nspinor,psit(1,iorb),1,psit(1,jorb),1)
!!              tti=ddot(nvctrp,psit(1,iorb),1,psit(nvctrp+1,jorb),1)
!!              tti=tti-ddot(nvctrp,psit(nvctrp+1,iorb),1,psit(1,jorb),1)
!!              tti=tti-ddot(nvctrp,psit(3*nvctrp+1,iorb),1,psit(2*nvctrp+1,jorb),1)
!!              tti=tti+ddot(nvctrp,psit(2*nvctrp+1,iorb),1,psit(3*nvctrp+1,jorb),1)
!!              ovrlp(2*iorb-1,jorb,1)=ttr
!!              ovrlp(2*iorb,jorb,1)=tti*0.0d0
!!              print *,iorb,norb,ttr,tti
!!           end do
!!        end do
!!        stop
        call herk('L','C',norb,ncomp*nvctrp,1.0_wp,psit(1,1),max(1,ncomp*nvctrp),&
             0.0_wp,ovrlp(1,1,istart),norb)
     end if

     if (nproc > 1) then
        call timing(iproc,'GramS_comput  ','OF')
        call timing(iproc,'GramS_commun  ','ON')
        call MPI_ALLREDUCE (ovrlp(1,1,2),ovrlp(1,1,1),norbs*norb,&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call timing(iproc,'GramS_commun  ','OF')
        call timing(iproc,'GramS_comput  ','ON')
     end if
     
!!     if (iproc==0) then
!!        write(*,*) 'parallel ovrlp'
!!        do i=1,norbs,min(2,nspinor)
!!           write(*,'(10(1x,1pe10.3))') (ovrlp(i,j,1),j=1,norb)
!!        enddo
!!     end if

     !to be excluded if nvctrp==0
     if(nspinor==1) then
        
        ! Cholesky factorization
        call potrf( 'L',norb,ovrlp(1,1,1),norb,info)
        if (info /= 0) then
           write(*,*) 'info Cholesky factorization',info
        end if
        
        ! calculate L^{-1}
        call trtri( 'L','N',norb,ovrlp(1,1,1),norb,info)
        if (info.ne.0) write(6,*) 'info L^-1',info
        
        ! new vectors   
        call trmm ('R','L','T','N',nvctrp,norb,1.0_wp,ovrlp(1,1,1),norb,psit(1,1),max(1,nvctrp))

     else

       ! Cholesky factorization
!!        do i=1,norb
!!           if(iproc==0) then
!!              write(*,*) 'parallel ovrlp',i
!!              write(*,'(10f10.3)') (ovrlp(j,i,1), j=1,norbs)
!!           end if
!!        end do
        call c_potrf( 'L',norb,ovrlp(1,1,1),norb,info )
        if (info /= 0) then
           write(*,*) 'info Cholesky factorization',info
        end if
        
        ! calculate L^{-1}
!!         do i=1,norb
!!           if(iproc==0) then
!!              write(*,*) 'parallel ovrlp2',i
!!              write(*,'(10f10.3)') (ovrlp(j,i,1), j=1,norbs)
!!           end if
!!        end do
       call c_trtri( 'L','N',norb,ovrlp(1,1,1),norb,info)
        if (info.ne.0) write(6,*) 'info L^-1',info
        
!!        do i=1,norb
!!           if(iproc==0) then
!!              write(*,'(10f10.3)') (ovrlp(j,i,1), j=1,norbs)
!!           end if
!!        end do
       ! new vectors   !!check if third argument should be transpose or conjugate
        call c_trmm ('R','L','C','N',ncomp*nvctrp,norb,(1.0_wp,0.0_wp),&
             ovrlp(1,1,1),norb,psit(1,1),max(1,ncomp*nvctrp))

        !if(nproc==1) call psitransspi(nvctrp,norb,psit,.true.)


     end if

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)

  end if

  enddo

  call timing(iproc,'GramS_comput  ','OF')

END SUBROUTINE orthon_p

!the loewe routines must be uniformised serial/parallel and nspinor should be added
subroutine loewe_p(iproc,nproc,norb,ndim,nvctrp,nvctr_tot,psit)
  ! loewdin orthogonalisation
  use module_base
  implicit real(kind=8) (a-h,o-z)
  logical, parameter :: parallel=.true.
  dimension psit(nvctrp,ndim)
  character(len=*), parameter :: subname='loewe_p'
  real(kind=8), allocatable :: ovrlp(:,:,:),evall(:),psitt(:,:)

  if (norb.eq.1) then

     nvctr_eff=min(nvctr_tot-iproc*nvctrp,nvctrp)

     if (nvctr_eff > 0) then
     !parallel treatment of a run with only one orbital
     tt=dnrm2(nvctr_eff,psit,1)     
     ttLOC=tt**2

     else
        ttLOC =0.d0
     end if
     
     call MPI_ALLREDUCE(ttLOC,tt,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

     tt=1.d0/sqrt(tt)
     call vscal(nvctr_eff,tt,psit(1,1),1)

     !stop 'more than one orbital needed for a parallel run'

  else

     allocate(ovrlp(norb,norb,3+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     allocate(evall(norb+ndebug),stat=i_stat)
     call memocc(i_stat,evall,'evall',subname)

     ! Upper triangle of overlap matrix using BLAS
     !     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; upper triangle
     call DSYRK('U','T',norb,nvctrp,1.d0,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

     ! Full overlap matrix using  BLAS
     !     ovrlap(jorb,iorb,2)=+psit(k,jorb)*psit(k,iorb)
     !      call DGEMM('T','N',norb,norb,nvctrp,1.d0,psit,&
     !                     nvctrp,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

     call MPI_ALLREDUCE(ovrlp(1,1,2),ovrlp(1,1,1),norb**2,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

     !       write(*,*) 'OVERLAP',iproc
     !       do i=1,norb
     !       write(*,'(10(x,e17.10))') (ovrlp(i,j,1),j=1,norb)
     !       enddo

     ! LAPACK
     call DSYEV('V','U',norb,ovrlp(1,1,1),norb,evall,ovrlp(1,1,3),norb**2,info)
     if (info.ne.0) write(6,*) 'info loewe', info
     !        if (iproc.eq.0) then 
     !          write(6,*) 'overlap eigenvalues'
     !77        format(8(1x,e10.3))
     !          if (norb.le.16) then
     !          write(6,77) evall
     !          else
     !          write(6,77) (evall(i),i=1,4), (evall(i),i=norb-3,norb)
     !          endif
     !        endif

     ! calculate S^{-1/2} ovrlp(*,*,3)
     do lorb=1,norb
        do jorb=1,norb
           ovrlp(jorb,lorb,2)=ovrlp(jorb,lorb,1)*sqrt(1.d0/evall(lorb))
        end do
     end do
     !        do 3985,j=1,norb
     !        do 3985,i=1,norb
     !        ovrlp(i,j,3)=0.d0
     !        do 3985,l=1,norb
     !3985    ovrlp(i,j,3)=ovrlp(i,j,3)+ovrlp(i,l,1)*ovrlp(j,l,2)
     ! BLAS:
     call DGEMM('N','T',norb,norb,norb,1.d0,ovrlp(1,1,1),norb,&
          ovrlp(1,1,2),norb,0.d0,ovrlp(1,1,3),norb)

     allocate(psitt(nvctrp,ndim+ndebug),stat=i_stat)
     call memocc(i_stat,psitt,'psitt',subname)
     ! new eigenvectors
     !   psitt(i,iorb)=psit(i,jorb)*ovrlp(jorb,iorb,3)
     call DGEMM('N','N',nvctrp,norb,norb,1.d0,psit,nvctrp,ovrlp(1,1,3),norb,0.d0,psitt,nvctrp)
     call DCOPY(nvctrp*ndim,psitt,1,psit,1)
     i_all=-product(shape(psitt))*kind(psitt)
     deallocate(psitt,stat=i_stat)
     call memocc(i_stat,i_all,'psitt',subname)

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)
     i_all=-product(shape(evall))*kind(evall)
     deallocate(evall,stat=i_stat)
     call memocc(i_stat,i_all,'evall',subname)

  end if

END SUBROUTINE loewe_p


subroutine loewe(norb,nvctrp,psi)
  ! loewdin orthogonalisation
  use module_base
  implicit real(kind=8) (a-h,o-z)
  dimension psi(nvctrp,norb)
  character(len=*), parameter :: subname='loewe'
  real(kind=8), allocatable :: ovrlp(:,:,:),evall(:),tpsi(:,:)

  if (norb.eq.1) then
     tt=0.d0
     do i=1,nvctrp
        tt=tt+psi(i,1)**2
     enddo
     tt=1.d0/sqrt(tt)
     do i=1,nvctrp
        psi(i,1)=psi(i,1)*tt
     enddo

  else

     allocate(ovrlp(norb,norb,3+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     allocate(evall(norb+ndebug),stat=i_stat)
     call memocc(i_stat,evall,'evall',subname)

     ! Overlap matrix using BLAS
     !     ovrlp(iorb,jorb)=psi(k,iorb)*psi(k,jorb) ; upper triangle
     call DSYRK('U','T',norb,nvctrp,1.d0,psi,nvctrp,0.d0,ovrlp(1,1,1),norb)

     !       write(*,*) 'OVERLAP'
     !       do i=1,norb
     !       write(*,'(10(1x,1pe17.10))') (ovrlp(i,j,1),j=1,norb)
     !       enddo


     ! LAPACK
     call DSYEV('V','U',norb,ovrlp(1,1,1),norb,evall,ovrlp(1,1,3),norb**2,info)
     if (info.ne.0) write(6,*) 'info loewe', info
     !          write(6,*) 'overlap eigenvalues'
     !77        format(8(1x,e10.3))
     !          if (norb.le.16) then
     !          write(6,77) evall
     !          else
     !          write(6,77) (evall(i),i=1,4), (evall(i),i=norb-3,norb)
     !          endif

     ! calculate S^{-1/2} ovrlp(*,*,3)
     do lorb=1,norb
        do jorb=1,norb
           ovrlp(jorb,lorb,2)=ovrlp(jorb,lorb,1)*sqrt(1.d0/evall(lorb))
        end do
     end do
     !        do 3985,j=1,norb
     !        do 3985,i=1,norb
     !        ovrlp(i,j,3)=0.d0
     !        do 3985,l=1,norb
     !3985    ovrlp(i,j,3)=ovrlp(i,j,3)+ovrlp(i,l,1)*ovrlp(j,l,2)
     ! BLAS:
     call DGEMM('N','T',norb,norb,norb,1.d0,ovrlp(1,1,1),norb,ovrlp(1,1,2),norb,0.d0,ovrlp(1,1,3),norb)

     ! new eigenvectors
     allocate(tpsi(nvctrp,norb+ndebug),stat=i_stat)
     call memocc(i_stat,tpsi,'tpsi',subname)
     !   tpsi(i,iorb)=psi(i,jorb)*ovrlp(jorb,iorb,3)
     call DGEMM('N','N',nvctrp,norb,norb,1.d0,psi(1,1),nvctrp,ovrlp(1,1,3),norb,0.d0,tpsi,nvctrp)
     call DCOPY(nvctrp*norb,tpsi,1,psi,1)
     i_all=-product(shape(tpsi))*kind(tpsi)
     deallocate(tpsi,stat=i_stat)
     call memocc(i_stat,i_all,'tpsi',subname)

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)
     i_all=-product(shape(evall))*kind(evall)
     deallocate(evall,stat=i_stat)
     call memocc(i_stat,i_all,'evall',subname)

  endif

END SUBROUTINE loewe


subroutine checkortho_p(iproc,nproc,norb,nvctrp,psit)
  use module_base
  implicit real(kind=8) (a-h,o-z)
  dimension psit(nvctrp,norb)
  character(len=*), parameter :: subname='checkortho_p'
  real(kind=8), allocatable :: ovrlp(:,:,:)

  allocate(ovrlp(norb,norb,2+ndebug),stat=i_stat)
  call memocc(i_stat,ovrlp,'ovrlp',subname)

  do iorb=1,norb
     do jorb=1,norb
        ovrlp(iorb,jorb,2)=ddot(nvctrp,psit(1,iorb),1,psit(1,jorb),1)
     end do
  end do

  call MPI_ALLREDUCE(ovrlp(1,1,2),ovrlp(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  toler=1.d-10
  dev=0.d0
  do iorb=1,norb
     do jorb=1,norb
        scpr=ovrlp(iorb,jorb,1)
        if (iorb.eq.jorb) then
           dev=dev+(scpr-1.d0)**2
        else
           dev=dev+scpr**2
        endif
        if (iproc == 0) then
           if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler) write(*,'(1x,a,2(1x,i0),1x,1pe12.6)')&
                'ERROR ORTHO',iorb,jorb,scpr
           if (iorb.ne.jorb .and. abs(scpr).gt.toler)      write(*,'(1x,a,2(1x,i0),1x,1pe12.6)')&
                'ERROR ORTHO',iorb,jorb,scpr
        end if
     end do
  end do

  if (dev.gt.toler) write(*,'(1x,a,i0,1pe13.5)') 'Deviation from orthogonality ',iproc,dev

  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp',subname)

END SUBROUTINE checkortho_p


subroutine checkortho(norb,nvctrp,psi)
  use module_base
  implicit real(kind=8) (a-h,o-z)
  dimension psi(nvctrp,norb)
  character(len=*), parameter :: subname='checkortho'
  real(kind=8), allocatable :: ovrlp(:,:,:)

  allocate(ovrlp(norb,norb,1+ndebug),stat=i_stat)
  call memocc(i_stat,ovrlp,'ovrlp',subname)

  do iorb=1,norb
     do jorb=1,norb
        ovrlp(iorb,jorb,1)=ddot(nvctrp,psi(1,iorb),1,psi(1,jorb),1)
     enddo
  enddo

  toler=1.d-10
  dev=0.d0
  do iorb=1,norb
     do jorb=1,norb
        scpr=ovrlp(iorb,jorb,1)
        if (iorb.eq.jorb) then
           dev=dev+(scpr-1.d0)**2
        else
           dev=dev+scpr**2
        endif
        if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler) write(*,'(1x,a,2(1x,i0),1x,1pe12.6)')&
             'ERROR ORTHO',iorb,jorb,scpr
        if (iorb.ne.jorb .and. abs(scpr).gt.toler)      write(*,'(1x,a,2(1x,i0),1x,1pe12.6)')&
             'ERROR ORTHO',iorb,jorb,scpr
     enddo
  enddo

  if (dev.gt.1.d-10) write(*,'(1x,a,i0,1pe13.5)') 'Deviation from orthogonality ',0,dev

  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp',subname)


END SUBROUTINE checkortho


subroutine KStrans_p(iproc,nproc,norb,nvctrp,occup,  & 
     hpsit,psit,evsum,eval,nspinor)
  ! at the start each processor has all the Psi's but only its part of the HPsi's
  ! at the end each processor has only its part of the Psi's
  !implicit real(kind=8) (a-h,o-z)
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norb,nvctrp,nspinor
  real(wp), intent(out) :: evsum
  real(gp), dimension(norb), intent(in) :: occup
  real(wp), dimension(nvctrp*nspinor,norb), intent(in) :: hpsit
  real(wp), dimension(norb), intent(out) :: eval
  real(wp), dimension(nvctrp*nspinor,norb), intent(out) :: psit
  !local variables
  character(len=*), parameter :: subname='KStrans_p'
  integer :: i_all,i_stat,ierr,iorb,jorb,n_lp,istart,info,norbs,ncomp
  real(wp) :: alpha
  ! arrays for KS orbitals
  real(wp), dimension(:), allocatable :: work_lp,work_rp
  real(wp), dimension(:,:), allocatable :: psitt
  real(wp), dimension(:,:,:), allocatable :: hamks

  if(nspinor==4) then
     norbs=2*norb
     ncomp=2
  else if (nspinor==2) then
     norbs=2*norb
     ncomp=1
  else
     norbs=norb
  end if
  ! set up Hamiltonian matrix
  allocate(hamks(norbs,norb,2+ndebug),stat=i_stat)
  call memocc(i_stat,hamks,'hamks',subname)

  do jorb=1,norb
     do iorb=1,norbs
        hamks(iorb,jorb,2)=0.0_wp
     enddo
  enddo
  if (nproc > 1) then
     istart=2
  else
     istart=1
  end if

  if(nspinor==1) then
!     do iorb=1,norb
!        do jorb=1,norb
!           scpr=ddot(nvctrp,psit(1,jorb),1,hpsit(1,iorb),1)
!           hamks(iorb,jorb,istart)=scpr
!        enddo
!     enddo
     call gemm('T','N',norb,norb,nvctrp,1.0_wp,psit(1,1),max(1,nvctrp),hpsit(1,1),max(1,nvctrp),0.0_wp,&
          hamks(1,1,istart),norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psit(1,1),max(1,ncomp*nvctrp), &
          hpsit(1,1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),hamks(1,1,istart),norb)
  end if

  if (nproc > 1) then
     call MPI_ALLREDUCE(hamks(1,1,2),hamks(1,1,1),norbs*norb,mpidtypw,&
          MPI_SUM,MPI_COMM_WORLD,ierr)
  end if
!  do iorb=1,norb
!     if(iproc==0) write(*,'(30f10.5)')(hamks(jorb,iorb,2),jorb=1,norbs)
!  end do
  !        write(*,*) 'KS Hamiltonian',iproc
  !        do iorb=1,norb
  !        write(*,'(10(1x,e10.3))') (hamks(iorb,jorb,1),jorb=1,norb)
  !        enddo

  n_lp=max(4*norbs,1000)
  allocate(work_lp(n_lp*2+ndebug),stat=i_stat)
  call memocc(i_stat,work_lp,'work_lp',subname)
  if(nspinor==1) then
     call  syev('V','U',norb,hamks(1,1,1),norb,eval(1),work_lp(1),n_lp,info)
  else
     allocate(work_rp(3*norb+1+ndebug),stat=i_stat)
     call memocc(i_stat,work_rp,'work_rp',subname)
     
     call  heev('V','U',norb,hamks(1,1,1),norb,eval(1),work_lp(1),n_lp,work_rp(1),info)
     
     i_all=-product(shape(work_rp))*kind(work_rp)
     deallocate(work_rp,stat=i_stat)
     call memocc(i_stat,i_all,'work_rp',subname)
  end if

  evsum=0.0_wp
  do iorb=1,norb
     evsum=evsum+eval(iorb)*real(occup(iorb),wp)
     !if (iproc.eq.0) write(*,'(1x,a,i0,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
  enddo
  i_all=-product(shape(work_lp))*kind(work_lp)
  deallocate(work_lp,stat=i_stat)
  call memocc(i_stat,i_all,'work_lp',subname)
  if (info.ne.0) write(*,*) 'DSYEV ERROR',info

  allocate(psitt(nvctrp*nspinor,norb+ndebug),stat=i_stat)
  call memocc(i_stat,psitt,'psitt',subname)
  ! Transform to KS orbitals
  ! dgemm can be used instead of daxpy
  if(nspinor==1) then
     do iorb=1,norb
        call razero(nvctrp,psitt(1,iorb))
        do jorb=1,norb
           alpha=hamks(jorb,iorb,1)
           call axpy(nvctrp,alpha,psit(1,jorb),1,psitt(1,iorb),1)
        enddo
     enddo
  else
     do iorb=1,norb
        call razero(nvctrp*nspinor,psitt(1,iorb))
        do jorb=1,norb
           call c_axpy(ncomp*nvctrp,hamks(2*jorb-1,iorb,1),psit(1,jorb),1,psitt(1,iorb),1)
        enddo
     enddo
  end if
  i_all=-product(shape(hamks))*kind(hamks)
  deallocate(hamks,stat=i_stat)
  call memocc(i_stat,i_all,'hamks',subname)

  call DCOPY(nvctrp*norb*nspinor,psitt,1,psit,1)
  i_all=-product(shape(psitt))*kind(psitt)
  deallocate(psitt,stat=i_stat)
  call memocc(i_stat,i_all,'psitt',subname)

END SUBROUTINE KStrans_p


! ********************************************************************************************************


subroutine gsChol(iproc, nproc, psi, input, nspinor, orbs, nspin,ndimovrlp,norbArr,comms)
!
! Purpose:
! =======
!  This subroutine orthonormalizes the orbitals psi in a parallel way. To do so, it first transposes the orbitals to all
!  processors using mpi_alltoallv. The orthonomalization is then done in this data layout using a combination of blockwise Gram-Schmidt
!  and Cholesky orthonomalization. At the end the vectors are again untransposed.
!
! Calling arguments:
! =================
!  Input arguments:
!    iproc     process ID
!    nproc     total number of processes
!    norb      total number of vectors that have to be orthonomalized, shared over all processes
!    input     data type containing many parameters
!  Input/Output arguments:
!    psi       on input: the vectors to be orthonormalized
!              on output: the orthonomalized vectors
!
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  !integer, intent(in) :: ikpt
  integer, intent(in) :: iproc, nproc, nspinor,nspin
  type(input_variables):: input
  type(orbitals_data):: orbs
  type(communications_arrays), intent(in) :: comms
  integer, dimension(nspin), intent(in) :: norbArr
  integer, dimension(nspin,0:orbs%nkpts), intent(in) :: ndimovrlp
  real(wp),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb),intent(inout):: psi
  
  ! Local variables
  integer:: iblock, jblock, ist, jst, iter, iter2, gcd, blocksize, blocksizeSmall, i_stat, i_all
  integer:: getBlocksize, ispin, ikptp, norbs, ncomp
  real(wp),dimension(:), allocatable :: ovrlp
  character(len=*), parameter:: subname='gsChol',category='GS/Chol'
   
  
  ! Make a loop over spin up/down.
  do ispin=1,nspin
     ! Get the blocksize.
     blocksize=getBlocksize(iproc, nproc, input, norbArr(ispin))
     
     ! There are two orthonormalization subroutines: gramschmidt orthogonalizes a given bunch of vectors to another bunch
     ! of already orthonormal vectors, and the subroutine cholesky orthonormalizes the given bunch.
     ! First determine how many bunches can be created for the given blocksize.
     iter=floor(real(norbArr(ispin))/real(blocksize))
     
     ! Get the dimensions of the overlap matrix for handling blocksize orbitals.
     call dimension_ovrlpFixedNorb(nspin,orbs,ndimovrlp,blocksize)
     allocate(ovrlp(ndimovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     
     ! Make a loop over all blocks.
     do iblock=1,iter
        ! ist is the starting orbitals of the current bunch of vectors.
        ist=(iblock-1)*blocksize+1
        ! Now orthogonalize this bunch to all previous ones.
        do jblock=1,iblock-1
           ! jst is the starting vector of the bunch to which the current bunch has to be orthogonalized.
           jst=blocksize*(jblock-1)+1
           call getOverlapDifferentPsi(iproc, nproc, nspin, blocksize,orbs, &
                comms, psi(1), ndimovrlp, ovrlp, norbArr, ist, jst, ispin, category)
           call gramschmidt(iproc, nproc, blocksize, psi(1), ndimovrlp, ovrlp, &
                orbs, nspin, nspinor, comms, norbArr, ist, jst, ispin)
        end do
    
        ! Orthonormalize the current bunch of vectors.
        call getOverlap(iproc, nproc, nspin, blocksize, orbs, comms, psi(1), &
             ndimovrlp, ovrlp, norbArr, ist, ispin, category)
        call cholesky(iproc, nproc, blocksize, psi(1), nspinor, nspin, orbs, &
             comms, ndimovrlp, ovrlp(1), norbArr, ist, ispin)
    
    end do

    i_all=-product(shape(ovrlp))*kind(ovrlp)
    deallocate(ovrlp, stat=i_stat)
    call memocc(i_stat,i_all,'ovrlp',subname)
    

    ! Orthonormalize the remaining vectors, if there are any.
    remainingIf: if(blocksize*iter/=norbArr(ispin)) then
        ! ist is the starting vector of the bunch that still havs to be orthonormalized.
        ist=blocksize*iter+1

        ! We have to find a new block size that matches both the remaining vectors and the already orthonomalized ones. This is done by determining
        ! the greatest common divisor of these two numbers.
        blocksizeSmall=gcd(blocksize*iter,norbArr(ispin)-ist+1)

        ! Get the dimensions of the overlap matrix for handling blocksize orbitals.
        call dimension_ovrlpFixedNorb(nspin,orbs,ndimovrlp,blocksizeSmall)
        allocate(ovrlp(ndimovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
        call memocc(i_stat,ovrlp,'ovrlp',subname)

        ! Determine how many blocks can be created with this new block size.
        iter2=(norbArr(ispin)-ist+1)/blocksizeSmall
        ! Now make a loop over all these blocks.
        do iblock=1,iter2
            ! ist is the starting vector of the current bunch.
            ist=iter*blocksize+blocksizeSmall*(iblock-1)+1
            ! Now orthogonalize this bunch to all previous ones.
            do jblock=1,(blocksize*iter)/blocksizeSmall+iblock-1
                ! jst is the starting vector of the bunch to which the current bunch has to be orthogonalized.
                jst=blocksizeSmall*(jblock-1)+1
                call getOverlapDifferentPsi(iproc, nproc, nspin, blocksizeSmall, &
                     orbs, comms, psi(1), ndimovrlp, ovrlp, norbArr, ist, jst, ispin, category)
                call gramschmidt(iproc, nproc, blocksizeSmall, psi(1), ndimovrlp, &
                     ovrlp, orbs, nspin, nspinor, comms, norbArr, ist, jst, ispin)
            end do
            ! Orthonormalize the current bunch of vectors.
            call getOverlap(iproc, nproc, nspin, blocksizeSmall, orbs, comms,&
                 psi(1), ndimovrlp, ovrlp, norbArr, ist, ispin, category)
            call cholesky(iproc, nproc, blocksizeSmall, psi(1), nspinor, nspin,&
                 orbs, comms, ndimovrlp, ovrlp(1), norbArr, ist, ispin)
        end do
        i_all=-product(shape(ovrlp))*kind(ovrlp)
        deallocate(ovrlp, stat=i_stat)
        call memocc(i_stat,i_all,'ovrlp',subname)
    end if remainingIf
    
end do

end subroutine gsChol


! ********************************************************************************************************


subroutine gramschmidt(iproc, nproc, norbIn, psit, ndimovrlp, ovrlp, orbs, nspin,&
     nspinor, comms, norbTot, block1, block2, ispinIn)
!
! Purpose:
! =======
!  This subroutine orthogonalizes a given bunch of vectors in psit to another bunch of equal size. These other vectors
!  are assumed to be orthonomal themselves. The starting indices of the two bunches are given by block1 and block2.
!  The orthonormalization is done in parallel, assuming that each process holds a small portion of each vector.
!
! Calling arguments:
! =================
!  Input arguments:
!    iproc       process ID
!    nproc       total number of processes
!    norbIn     number of orbitals to be orthonormalized
!    ndimovrlp      describes the shape of the overlap matrix
!    orbs       type that contains many parameters concerning the orbitals
!    nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!    nspinor    real wavefunction -> nspinor=1, complex wavefunction -> nspinor>1
!    comms      type containing parameters for communicating the wavefunstion between processors
!    norbTot    total number of orbitals (if nspin=2:
!                 norbTot(1)=total number of up orbitals
!                 norbTot(2)=total number of down orbitals)
!    block1     gives the starting orbital of the orbitals to be orthogonalized
!    block2     gives the starting orbital of the orbitals to which they shall be orthogonalized
!    ispinIn    indicates whether the up or down orbitals shall be handled
!  Input/Output arguments:
!    psit       the vectors that shall be orthonormalized
!    ovrlp      the overlap matrix which will be destroyed during this subroutine
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norbIn, nspin, nspinor, block1, block2, ispinIn
type(orbitals_data):: orbs
type(communications_arrays), intent(in) :: comms
real(wp),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb),intent(inout):: psit
integer,dimension(nspin,0:orbs%nkpts):: ndimovrlp
real(wp),dimension(ndimovrlp(nspin,orbs%nkpts)):: ovrlp
integer,dimension(nspin):: norbTot

! Local arguments
integer:: nvctrp, ist, ierr, i_stat, i_all, ncomp, ikptp, ikpt, ispin, norb, norbs, istThis, istOther
real(kind=8),dimension(:),allocatable:: A1D
character(len=*),parameter:: subname='gramschmidt'

! Initialize the starting indices. istThis is the starting index of the orbitals that shall be orthogonalized,
! istOther is the starting index of the orbitals to which they shall be orthogonalized.
istThis=1
istOther=1

! Make a loop over the number of k-points handled by the process.
do ikptp=1,orbs%nkptsp
    ! ikpt is the number of the k-point.
    ikpt=orbs%iskpts+ikptp
    ! Now make a loop over spin up and down.
    do ispin=1,nspin
        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for.
        ! In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
            nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we possibly treat only a part of the orbitals).
        norb=norbIn

        ! Allocate the matrix A which will hold some partial results.
        allocate(A1D(nvctrp*norb*nspinor), stat=i_stat)
        call memocc(i_stat, A1D, 'A1D', subname)

        ! Count up the starting indices.
        istThis=istThis+nvctrp*(block1-1)*nspinor
        istOther=istOther+nvctrp*(block2-1)*nspinor

        if(ispin==ispinIn) then
            ! Calculate matrix product psit*ovrlp=A. This will give the components that will be projected out of psit.
            ! We actually calculate -psit*ovrlp=-A, since this is better for further processing with daxpy.
            if(nspinor==1) then
                call dgemm('n', 'n', nvctrp, norb, norb, -1.d0, psit(istOther), nvctrp,&
                     ovrlp(ndimovrlp(ispin,ikpt-1)+1), norb, 0.d0, A1D(1), nvctrp)
            else
                call zgemm('n', 'n', nvctrp, norb, norb, (-1.d0,0.d0), psit(istOther), nvctrp, &
                     ovrlp(ndimovrlp(ispin,ikpt-1)+1), norb, (0.d0,0.d0), A1D(1), nvctrp)
            end if
            ! Now project out: psit=psit-A.
            ! Since we calculated -A, we have to put psit=psit+A and can use daxpy to perform psit=A+psit
            if(nspinor==1) then
                call daxpy(nvctrp*norb*nspinor,1.d0,A1D(1),1,psit(istThis),1)
            else
                call daxpy(nvctrp*norb*nspinor,1.d0,A1D(1),1,psit(istThis),1)
            end if
        end if

        ! Increase the starting indices. This will bring the starting index to the start of the the next spin case (up/down) and k-point.
        istThis=istThis+nvctrp*(norbTot(ispin)-block1+1)*nspinor
        istOther=istOther+nvctrp*(norbTot(ispin)-block2+1)*nspinor

        i_all=-product(shape(A1D))*kind(A1D)
        deallocate(A1D)
        call memocc(i_stat,i_all,'A1D',subname)
    end do
end do

end subroutine gramschmidt


! ********************************************************************************************************


subroutine cholesky(iproc, nproc, norbIn, psi, nspinor, nspin, orbs, comms, ndimL, Lc, norbTot, block1, ispinIn)
!
! Purpose:
! =======
!  This subroutine orthonormalizes a given bunch of vectors psi.
!  It first calculates the Cholesky composition S=L*L^T of the overlap matrix S.  This matrix L is then
!  inverted to get L^{-1} and the orthonormal vectors are finally given by psi=psi*L^{-1}.
!
! Calling arguments:
! =================
!  Input arguments:
!    iproc      process ID
!    nproc      total number of processes
!    norbIn     number of orbitals to be orthonormalized
!    nspinor    real wavefunction -> nspinor=1, complex wavefunction -> nspinor>1
!    nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!    orbs       type that contains many parameters concerning the orbitals
!    comms      type containing parameters for communicating the wavefunstion between processors
!    ndimL      describes the shape of the overlap matrix
!    norbTot    total number of orbitals (if nspin=2:
!                 norbTot(1)=total number of up orbitals
!                 norbTot(2)=total number of down orbitals)
!    block1     gives the starting orbital of the orbitals to be orthonormalized
!    ispinIn    indicates whether the up or down orbitals shall be handled
!  Input/Output arguments:
!    psi        the vectors that shall be orthonormalized
!    ovrlp      the overlap matrix which will be destroyed during this subroutine
!
use module_base
use module_types
implicit none

! Calling arguments
!integer:: iproc,nproc,nvctrp,norbIn, nspinor, nspin, norbTot, block1, ispinIn
integer:: iproc,nproc,nvctrp,norbIn, nspinor, nspin, block1, ispinIn
type(orbitals_data):: orbs
type(communications_arrays):: comms
real(kind=8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb),intent(in out):: psi
integer,dimension(nspin,0:orbs%nkpts):: ndimL
real(kind=8),dimension(ndimL(nspin,orbs%nkpts),1):: Lc
integer,dimension(nspin):: norbTot

! Local variables
integer:: ist, info, i_stat, i_all, ispin, ikptp, ikpt, ncomp, norbs, norb
character(len=*),parameter:: subname='cholesky'
  
 
! Set the starting index to 1.
ist=1
! Make a loop over the number of k-points handled by the process.
do ikptp=1,orbs%nkptsp
    ! ikpt is the number of the k-point.
    ikpt=orbs%iskpts+ikptp
    ! Now make a loop over spin up and down.
    do ispin=1,nspin
        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for.
        ! In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
            nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we possibly treat only a part of the orbitals).
        norb=norbIn
        ! Count up the starting index
        ist=ist+nvctrp*(block1-1)*nspinor
 
        ! The following part is only executed if ispin==ispinIn. Otherwise only the starting index ist
        ! is increased.
        if(ispin==ispinIn) then
            ! Make a Cholesky factorization of L.
            if(nspinor==1) then
                call dpotrf('l', norb, Lc(ndimL(ispin,ikpt-1)+1,1), norb, info)
            else
                call zpotrf('l', norb, Lc(ndimL(ispin,ikpt-1)+1,1), norb, info)
            end if
            
            ! Invert the Cholesky matrix: L^{-1}.
            if(nspinor==1) then
                call dtrtri('l', 'n', norb, Lc(ndimL(ispin,ikpt-1)+1,1), norb, info)
            else
                call ztrtri('l', 'n', norb, Lc(ndimL(ispin,ikpt-1)+1,1), norb, info)
            end if
 
            ! Calculate the matrix product psi*L^{-1}=psi. This will give the orthonormal orbitals.
            if(nspinor==1) then
                call dtrmm('r', 'l', 't', 'n', nvctrp, norb, 1.d0, &
                     Lc(ndimL(ispin,ikpt-1)+1,1), norb, psi(ist), nvctrp)
            else
                call ztrmm('r', 'l', 'c', 'n', ncomp*nvctrp, norb, (1.d0,0.d0),&
                     Lc(ndimL(ispin,ikpt-1)+1,1), norb, psi(ist), ncomp*nvctrp)
            end if
        end if
 
        ! Increase the starting index.
        ist=ist+nvctrp*(norbTot(ispin)-block1+1)*nspinor

    end do
end do         

end subroutine cholesky


! ********************************************************************************************************


subroutine loewdin(iproc,nproc, norbIn, nspinor, block1, ispinIn, orbs, comms, nspin, psit, ovrlp, ndimovrlp, norbTot)
!
! Purpose:
! =======
!   Orthonormalizes the vectors provided in psit by a loewdin orthonormalization.
!
! Calling arguments:
! =================
!  Input arguments:
!    iproc      process ID
!    nproc      number of processes
!    norbIn     number of orbitals to be orthonormalized
!    nspinor    real wavefunction -> nspinor=1, complex wavefunction -> nspinor>1
!    block1     gives the starting orbital of the orbitals to be orthonormalized
!    ispinIn    indicates whether the up or down orbitals shall be handled
!    orbs       type that contains many parameters concerning the orbitals
!    comms      type containing parameters for communicating the wavefunstion between processors
!    nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!    ndimovrlp  describes the shape of the overlap matrix
!    norbTot    total number of orbitals (if nspin=2:
!                 norbTot(1)=total number of up orbitals
!                 norbTot(2)=total number of down orbitals)
!  Input/output Arguments
!    psit       the orbitals to be orthonormalized
!    ovrlp      the overlap matrix which will be destroyed during this subroutine
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc,nproc,norbIn, nspinor, nspin, block1, ispinIn
type(orbitals_data),intent(in):: orbs
type(communications_arrays),intent(in):: comms
real(kind=8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb),intent(in out):: psit
integer,dimension(nspin,0:orbs%nkpts):: ndimovrlp
real(kind=8),dimension(ndimovrlp(nspin,orbs%nkpts)):: ovrlp
integer,dimension(nspin):: norbTot

! Local variables
integer:: jorb, lorb, i_stat, i_all, info, nvctrp, ispin, ist, ikptp, ikpt, ncomp, norbs, norb, lwork
real(kind=8),dimension(:),allocatable:: evall, psitt
real(kind=8),dimension(:,:),allocatable:: tempArr
character(len=*), parameter :: subname='loewdin'

! Allocate the work arrays.
lwork=nspinor*norbIn**2+10
allocate(tempArr(norbIn**2*nspinor,2), stat=i_stat)
call memocc(i_stat,tempArr,'tempArr',subname)

allocate(evall(norbIn), stat=i_stat)
call memocc(i_stat,evall,'evall',subname)

ist=1
! Make a loop over the number of k-points handled by the process.
do ikptp=1,orbs%nkptsp
    ! ikpt is the number of the k-point.
    ikpt=orbs%iskpts+ikptp
    ! Now make a loop over spin up and down.
    do ispin=1,nspin
        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for.
        ! In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
            nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we possibly treat only a part of the orbitals).
        norb=norbIn
        ! Count up the starting index
        ist=ist+nvctrp*(block1-1)*nspinor
 
        ! The following part is only executed if ispin==ispinIn. Otherwise only the starting index ist
        ! is increased.
        if(ispin==ispinIn) then

            ! Diagonalize the overlap matrix.
            if(nspinor==1) then
                call dsyev('v', 'l', norb, ovrlp(ndimovrlp(ispin,ikpt-1)+1), norb,&
                     evall, tempArr(1,1), lwork, info)
            else
                call zheev('v', 'l', norb,ovrlp(ndimovrlp(ispin,ikpt-1)+1), norb,&
                     evall, tempArr(1,1), lwork, tempArr(1,2), info)
            end if
            if (info/=0) then
                write(*,'(a,i0)') 'ERROR in dsyev (Loewdin); info=',info
                stop
            end if

            ! Calculate S^{-1/2}. 
            ! First calulate ovrlp*diag(evall) (ovrlp is the diagonalized overlap
            ! matrix and diag(evall) the diagonal matrix consisting of the eigenvalues...
            do lorb=1,norb
               do jorb=1,norb*nspinor
                  tempArr((lorb-1)*norb*nspinor+jorb,1)=&
                       ovrlp(ndimovrlp(ispin,ikpt-1)+(lorb-1)*norb*nspinor+jorb)*sqrt(1.d0/evall(lorb))
               end do
            end do

            ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
            ! This will give S^{-1/2}.
            if(nspinor==1) then
                call dgemm('n', 't', norb, norb, norb, 1.d0, ovrlp(ndimovrlp(ispin,ikpt-1)+1), norb,&
                     tempArr(1,1), norb, 0.d0, tempArr(1,2), norb)
            else
                call zgemm('n', 't', norb, norb, norb, (1.d0,0.d0), ovrlp(ndimovrlp(ispin,ikpt-1)+1), norb,&
                     tempArr(1,1), norb, (0.d0,0.d0), tempArr(1,2), norb)
            end if

            ! Now calculate the orthonormal orbitals by applying S^{-1/2} to the orbitals.
            ! This requires the use of a temporary variable psitt.
            allocate(psitt(nvctrp*norb*nspinor),stat=i_stat)
            call memocc(i_stat,psitt,'psitt',subname)
            if(nspinor==1) then
                call dgemm('n', 'n', nvctrp, norb, norb, 1.d0, psit(ist), &
                     nvctrp, tempArr(1,2), norb, 0.d0, psitt, nvctrp)
            else
                call zgemm('n', 'n', nvctrp, norb, norb, (1.d0,0.d0), &
                     psit(ist), nvctrp, tempArr(1,2), norb, (0.d0,0.d0), psitt, nvctrp)
            end if

            ! Now copy the orbitals from the temporary variable to psit.
            call dcopy(nvctrp*norb*nspinor, psitt(1), 1, psit(ist), 1)

            ! Deallocate the temporary variable psitt.
            i_all=-product(shape(psitt))*kind(psitt)
            deallocate(psitt,stat=i_stat)
            call memocc(i_stat,i_all,'psitt',subname)

        end if
        ! Increase the starting index.
        ist=ist+nvctrp*(norbTot(ispin)-block1+1)*nspinor

    end do
end do         


! Deallocate the remaining arrays.
i_all=-product(shape(tempArr))*kind(tempArr)
deallocate(tempArr,stat=i_stat)
call memocc(i_stat,i_all,'tempArr',subname)

i_all=-product(shape(evall))*kind(evall)
deallocate(evall,stat=i_stat)
call memocc(i_stat,i_all,'evall',subname)

end subroutine loewdin


! ********************************************************************************************************


subroutine getOverlap(iproc,nproc,nspin,norbIn,orbs,comms,&
     psi,ndimovrlp,ovrlp,norbTot,block1,ispinIn,category)
  !
  ! Purpose:
  ! =======
  !  This subroutine calculates the overlap matrix for a given bunch of orbitals. It also takes into 
  !  account k-points and spin.
  !
  ! Calling arguments:
  ! =================
  !  Input arguments:
  !    iproc      process ID
  !    nproc      total number of processes
  !    nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
  !    norbIn     number of orbitals to be orthonormalized
  !    orbs       type that contains many parameters concerning the orbitals
  !    comms      type containing parameters for communicating the wavefunstion between processors
  !    ndimovrlp      describes the shape of the overlap matrix
  !    norbTot    total number of orbitals (if nspin=2:
  !                 norbTot(1)=total number of up orbitals
  !                 norbTot(2)=total number of down orbitals)
  !    block1     gives the starting orbital of the orbitals to be orthonormalized
  !    ispinIn    indicates whether the up or down orbitals shall be handled
  !    catgeory   gives the category for the timing
  !  Output arguments:
  !    ovrlp      the overlap matrix of the orbitals given in psi
  !
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  character(len=*), intent(in) :: category
  integer,intent(in):: iproc,nproc,nspin,norbIn,block1,ispinIn
  type(orbitals_data),intent(in):: orbs
  type(communications_arrays),intent(in) :: comms
  real(wp),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb),intent(in) :: psi
  integer,dimension(nspin,0:orbs%nkpts),intent(in):: ndimovrlp
  real(wp),dimension(ndimovrlp(nspin,orbs%nkpts)),intent(out):: ovrlp
  integer,dimension(nspin),intent(in):: norbTot

  ! Local variables
  integer:: ispsi,ikptp,ikpt,ispin,nspinor,ncomp,norbs,ierr,nvctrp,norb



  ! Set the whole overlap matrix to zero. This is necessary since each process treats only a part
  ! of the matrix.
  call razero(ndimovrlp(nspin,orbs%nkpts),ovrlp)


  ispsi=1
  ! First make a loop over the k points handled by this process.
  do ikptp=1,orbs%nkptsp
     ! ikpt is the index of the k point.
     ikpt=orbs%iskpts+ikptp

     ! Now make also a loop over spin up/down.
     do ispin=1,nspin

        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for which the overlap
        ! matrix shall be calculated. In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we treat only a part of the orbitals).
        norb=norbIn
        ! Put the starting index to the right place. The current block of vector starts at the block1-th vector.
        ispsi=ispsi+nvctrp*(block1-1)*nspinor
        if(ispin==ispinIn) then
           if (nvctrp == 0) cycle

           ! Now calclulate one part of the overlap matrix. The starting index of this part is given by ndimovrlp(ispin,ikpt-1)+1.
           if(nspinor==1) then
              call syrk('L','T',norb,nvctrp,1.0_wp,psi(ispsi),max(1,nvctrp),&
                   0.0_wp,ovrlp(ndimovrlp(ispin,ikpt-1)+1),norb)
           else
              call herk('L','C',norb,ncomp*nvctrp,1.0_wp,psi(ispsi),&
                   max(1,ncomp*nvctrp),0.0_wp,ovrlp(ndimovrlp(ispin,ikpt-1)+1),norb)
           end if
        end if
        ! Move the starting indices to the end of the actual k point. This is necessary since nvctrp is
        ! different for the next k point and we cannot jump directly to the starting indices of our block for 
        ! the next k point.
        ispsi=ispsi+nvctrp*(norbTot(ispin)-block1+1)*nspinor
     end do
  end do

  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !print *,'here',iproc

  if (nproc > 1) then
     !call timing(iproc,'GramS_comput  ','OF')
     !call timing(iproc,'GramS_commun  ','ON')
     call timing(iproc, trim(category)//'_comput', 'OF')
     call timing(iproc, trim(category)//'_commun', 'ON')
     call mpiallred(ovrlp(1),ndimovrlp(nspin,orbs%nkpts),MPI_SUM,MPI_COMM_WORLD,ierr)
     !call MPI_ALLREDUCE (ovrlp(1,2),ovrlp(1,1),ndimovrlp(nspin,orbs%nkpts),mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     call timing(iproc, trim(category)//'_commun', 'OF')
     call timing(iproc, trim(category)//'_comput', 'ON')
     !call timing(iproc,'GramS_commun  ','OF')
     !call timing(iproc,'GramS_comput  ','ON')
  end if

  ! Now each processors knows all the overlap matrices for each k-point
  ! even if it does not handle it.
  ! This is somehow redundant but it is one way of reducing the number of communications
  ! without defining group of processors.

end subroutine getOverlap


! ********************************************************************************************************


subroutine getOverlapDifferentPsi(iproc, nproc, nspin, norbIn, orbs, comms,&
     psit, ndimovrlp, ovrlp, norbTot, block1, block2, ispinIn, category)
!
! Purpose:
! =======
!  This subroutine calculates the overlap matrix for a given bunch of orbitals. It also takes into 
!  account k-points and spin.
!
! Calling arguments:
! =================
!  Input arguments:
!    iproc      process ID
!    nproc      total number of processes
!    nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!    norbIn     number of orbitals to be orthonormalized
!    istart     second dimension of the overlpa matrix
!    orbs       type that contains many parameters concerning the orbitals
!    comms      type containing parameters for communicating the wavefunstion between processors
!    psit    the orbitals 
!    ndimovrlp  describes the shape of the overlap matrix
!    norbTot    total number of orbitals (if nspin=2:
!                 norbTot(1)=total number of up orbitals
!                 norbTot(2)=total number of down orbitals)
!    block1     gives the starting orbital of the orbitals to be orthogonalized
!    block2     gives the starting orbital of the orbitals to which the orbitals shall orthogonalized
!    ispinIn    indicates whether the up or down orbitals shall be handled
!    category   gives the category for the timing
!  Output arguments:
!    ovrlp      the overlap matrix of the orbitals given in psi
!
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  !integer,intent(in):: iproc, nproc, nspin, norbIn,  istart, norbTot, block1, block2
  character(len=*), intent(in) :: category
  integer,intent(in):: iproc, nproc, nspin, norbIn, block1, block2, ispinIn
  type(orbitals_data),intent(in):: orbs
  type(communications_arrays),intent(in) :: comms
  real(kind=8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb),intent(in) :: psit
  integer,dimension(nspin,0:orbs%nkpts),intent(in):: ndimovrlp
  real(kind=8),dimension(ndimovrlp(nspin,orbs%nkpts)):: ovrlp
  integer,dimension(nspin):: norbTot
  ! Local variables
  integer:: ikptp, ikpt, ispin, nspinor, ncomp, norbs, ierr, nvctrp, norb, ispsi1, ispsi2
  
  ! Set the whole overlap matrix to zero. This is necessary since each process treats only a part
  ! of the matrix.
  call razero(ndimovrlp(nspin,orbs%nkpts),ovrlp)

  ispsi1=1
  ispsi2=1
  ! First make a loop over the k points handled by this process.
  do ikptp=1,orbs%nkptsp
     ! ikpt is the index of the k point.
     ikpt=orbs%iskpts+ikptp
     
     ! Now make also a loop over spin up/down.
     do ispin=1,nspin
        
        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for which the overlap
        ! matrix shall be calculated. In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we treat only a part of the orbitals).
        norb=norbIn
        
        ! Put the starting index to the right place. The current block of vector starts at the block1-th and
        ! block2-th vector, respectively. 
        ispsi1=ispsi1+nvctrp*(block1-1)*nspinor
        ispsi2=ispsi2+nvctrp*(block2-1)*nspinor
        if(ispin==ispinIn) then
            if (nvctrp == 0) cycle
       
            ! Now calclulate one part of the overlap matrix. The starting index of this part is given by ndimovrlp(ispin,ikpt-1)+1.
            if(nspinor==1) then
               call gemm('t','n',norb,norb,ncomp*nvctrp,1.0_wp,psit(ispsi2),&
                    ncomp*nvctrp,psit(ispsi1),ncomp*nvctrp,0.d0,ovrlp(ndimovrlp(ispin,ikpt-1)+1),norb)
            else
               call c_gemm('c','n',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psit(ispsi2),&
                    ncomp*nvctrp,psit(ispsi1),ncomp*nvctrp,(0.d0,0.d0),ovrlp(ndimovrlp(ispin,ikpt-1)+1),norb)
            end if

        end if
        ! Move the starting indices to the end of the actual k point. This is necessary since nvctrp is
        ! different for the next k point and we cannot jump directly to the starting indices of our block for 
        ! the next k point.
        ispsi1=ispsi1+nvctrp*(norbTot(ispin)-block1+1)*nspinor
        ispsi2=ispsi2+nvctrp*(norbTot(ispin)-block2+1)*nspinor

     end do
  end do

  ! Sum up the overlap matrices from all processes.
  if (nproc > 1) then
     !call timing(iproc,'GramS_comput  ','OF')
     !call timing(iproc,'GramS_commun  ','ON')
     call timing(iproc,trim(category)//'_comput','OF')
     call timing(iproc,trim(category)//'_commun','ON')
     call mpiallred(ovrlp(1),ndimovrlp(nspin,orbs%nkpts),MPI_SUM,MPI_COMM_WORLD,ierr)
     !call mpi_allreduce(ovrlp(1,2),ovrlp(1,1),ndimovrlp(nspin,orbs%nkpts),mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
     call timing(iproc,trim(category)//'_commun','OF')
     call timing(iproc,trim(category)//'_comput','ON')
     !call timing(iproc,'GramS_commun  ','OF')
     !call timing(iproc,'GramS_comput  ','ON')
  end if
  
  ! Now each processors knows all the overlap matrices for each k-point even if it does not handle it.
  
end subroutine getOverlapDifferentPsi


! ********************************************************************************************************


subroutine dimension_ovrlpFixedNorb(nspin,orbs,ndimovrlp,norb)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspin,norb
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(nspin,0:orbs%nkpts), intent(out) :: ndimovrlp
  !local variables
  integer :: norbs,ncomp,ikpt

  ndimovrlp(1,0)=0
  if (nspin == 2) then
     !norb=orbs%norbu

     !this is first k-point
     call complex_components(orbs%nspinor,norb,norbs,ncomp)

     ndimovrlp(2,0)=norbs*norb
  end if

  do ikpt=1,orbs%nkpts
     !this part should be enhanced for real k-points
     !norb=orbs%norbu
     !if (nspin == 2) norb = orbs%norbd
     !this is ikpt k-point
     call complex_components(orbs%nspinor,norb,norbs,ncomp)

     ndimovrlp(1,ikpt)=ndimovrlp(nspin,ikpt-1)+norbs*norb
     if (orbs%norbd > 0) then
        !norb=orbs%norbu
        !this is ikpt+1
        call complex_components(orbs%nspinor,norb,norbs,ncomp)
        if (ikpt == orbs%nkpts) then
           ndimovrlp(2,ikpt)=ndimovrlp(1,ikpt)
        else
           ndimovrlp(2,ikpt)=ndimovrlp(1,ikpt)+norbs*norb
        end if
     end if
  end do

END SUBROUTINE dimension_ovrlpFixedNorb
