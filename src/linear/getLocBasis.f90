subroutine getLinearPsi(iproc,nproc,lzd,orbs,lorbs,llborbs,comsr,&
    mad,lbmad,op,lbop,comon,lbcomon,comgp,lbcomgp,at,rxyz,denspot,&
    GPU,updatePhi,&
    infoBasisFunctions,infoCoeff,itSCC,ebs,coeff,lphi,nlpspd,proj,communicate_lphi,coeff_proj,&
    ldiis,nit,nItInnerLoop,newgradient,orthpar,confdatarr,&
    methTransformOverlap,blocksize_pdgemm,convCrit,nItPrecond,&
    useDerivativeBasisFunctions,lphiRestart,comrp,blocksize_pdsyev,nproc_pdsyev,&
    hx,hy,hz,SIC,factor_enlarge)
!
! Purpose:
! ========
!   This subroutine creates the orbitals psi out of a linear combination of localized basis functions
!   phi. To do so, it proceeds as follows:
!    1. Create the basis functions (with subroutine 'getLocalizedBasis')
!    2. Write the Hamiltonian in this new basis.
!    3. Diagonalize this Hamiltonian matrix.
!    4. Build the new linear combinations. 
!   The basis functions are localized by adding a confining quartic potential to the ordinary DFT 
!   Hamiltonian. There is no self consistency cycle for the potential, i.e. the basis functionsi
!   are optimized with a fixed potential.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc           process ID
!     nproc           total number of processes
!     Glr             type describing the localization region
!     orbs            type describing the physical orbitals psi
!     at              type containing the paraneters for the atoms
!     lin             type containing parameters for the linear version
!     rxyz            the atomic positions
!     rxyzParab       the center of the confinement potential (at the moment identical rxyz)
!     nscatterarr     ???
!     ngatherarr      ???
!     nlpsp           ???
!     rhopot          the charge density
!     GPU             parameters for GPUs
!     input           type containing some very general parameters
!     pkernelseq      ???
!     n3p             ???
!     itSCC           iteration in the self consistency cycle
!  Input/Output arguments
!  ---------------------
!     phi             the localized basis functions. It is assumed that they have been initialized
!                     somewhere else
!   Output arguments
!   ----------------
!     psi             the physical orbitals, which will be a linear combinations of the localized
!                     basis functions phi
!     psit            psi transposed
!     infoBasisFunctions  indicated wheter the basis functions converged to the specified limit (value is the
!                         number of iterations it took to converge) or whether the iteration stopped due to 
!                         the iteration limit (value is -1). This info is returned by 'getLocalizedBasis'
!     infoCoeff           the same as infoBasisFunctions, just for the coefficients. This value is returned
!                         by 'optimizeCoefficients'
!
use module_base
use module_types
use module_interfaces, exceptThisOne => getLinearPsi
use Poisson_Solver
!use deallocatePointers
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, itSCC, nit, nItInnerLoop
integer,intent(in):: methTransformOverlap, blocksize_pdgemm, nItPrecond
integer,intent(in):: blocksize_pdsyev, nproc_pdsyev
type(local_zone_descriptors),intent(inout):: lzd
type(orbitals_data),intent(in) :: orbs
type(orbitals_data),intent(inout):: lorbs, llborbs
!type(p2pCommsSumrho),intent(inout):: comsr
type(p2pComms),intent(inout):: comsr
type(matrixDescriptors),intent(inout):: mad, lbmad
type(overlapParameters),intent(inout):: op, lbop
type(p2pComms),intent(inout):: comon, lbcomon
!type(p2pCommsGatherPot):: comgp, lbcomgp
type(p2pComms):: comgp, lbcomgp
type(atoms_data),intent(in):: at
real(8),dimension(3,at%nat),intent(in):: rxyz
type(DFT_local_fields), intent(inout) :: denspot
!integer,dimension(0:nproc-1,4),intent(inout):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!integer,dimension(0:nproc-1,2),intent(inout):: ngatherarr
!real(dp),dimension(max(lzd%Glr%d%n1i*lzd%Glr%d%n2i*n3p,1)*nspin),intent(inout) :: rhopot
type(GPU_pointers),intent(inout):: GPU
!real(dp), dimension(size_pkernel),intent(in):: pkernel
logical,intent(in):: updatePhi, newgradient, useDerivativeBasisFunctions
!real(dp),dimension(:),pointer,intent(in):: pkernelseq
integer,intent(out):: infoBasisFunctions, infoCoeff
real(8),intent(out):: ebs
real(8),intent(in):: convCrit, hx, hy, hz, factor_enlarge
real(8),dimension(llborbs%norb,orbs%norb),intent(in out):: coeff
!real(8),dimension(max(llborbs%npsidim_orbs,llborbs%npsidim_comp)),intent(inout):: lphi
real(8),dimension(:),pointer,intent(inout):: lphi
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
logical,intent(in):: communicate_lphi
real(8),dimension(lorbs%norb,orbs%norb),intent(inout):: coeff_proj
type(localizedDIISParameters),intent(inout):: ldiis
type(orthon_data),intent(in):: orthpar
type(confpot_data),dimension(lorbs%norbp),intent(in) :: confdatarr
!real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout)::lphiRestart
real(8),dimension(:),pointer,intent(inout)::lphiRestart
type(p2pCommsRepartition),intent(inout):: comrp
type(SIC_data),intent(in):: SIC

! Local variables 
integer:: istat, iall, ilr, istr, iorb, jorb, korb, tag, norbu, norbd, nspin, npsidim, norb, nlr
real(8),dimension(:),allocatable:: eval, lhphi, locrad
real(8),dimension(:,:),allocatable:: HamSmall, ovrlp, overlapmatrix, locregCenter
real(8),dimension(:,:,:),allocatable:: matrixElements
real(8):: epot_sum, ekin_sum, eexctX, eproj_sum, trace, tt, ddot, tt2, dnrm2, t1, t2, time,eSIC_DC
character(len=*),parameter:: subname='getLinearPsi' 
logical:: withConfinement
type(workarr_sumrho):: w
integer:: ist, ierr, iiorb, info, lorb, lwork, norbtot, k, l, ncnt, inc, jjorb, ii
!real(8),dimension(:),pointer:: lpot
type(confpot_data),dimension(:),allocatable :: confdatarrtmp
integer,dimension(:),allocatable:: orbsperlocreg
type(locreg_descriptors):: glr_tmp
integer:: ilrold, iiprocold, iiproc, jjlr, jjproc, i, gdim, ldim, ind, klr, kkorb, jst, iseg, istart, iend, ncount, indovrlp 



  ! Allocate the local arrays.  
  allocate(matrixElements(llborbs%norb,llborbs%norb,2), stat=istat)
  call memocc(istat, matrixElements, 'matrixElements', subname)
  allocate(HamSmall(llborbs%norb,llborbs%norb), stat=istat)
  call memocc(istat, HamSmall, 'HamSmall', subname)
  allocate(eval(llborbs%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)
  allocate(ovrlp(llborbs%norb,llborbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)


  ! This is a flag whether the basis functions shall be updated.
  if(updatePhi) then
      ! If we use the derivative basis functions, the trace minimizing orbitals of the last iteration are
      ! stored in lin%lphiRestart.
      ! WARNING: Will probably not work if we use the random input guess
      if(useDerivativeBasisFunctions) then
          !call dcopy(max(lorbs%npsidim_orbs,lorbs%npsidim_comp),lphiRestart(1),1,lphi(1),1)
          call dcopy(lorbs%npsidim_orbs,lphiRestart(1),1,lphi(1),1)
      end if


      ! Improve the trace minimizing orbitals.
      call getLocalizedBasis(iproc,nproc,at,lzd,lorbs,orbs,comon,op,comgp,mad,rxyz,&
           denspot,GPU,lphi,trace,&
          infoBasisFunctions,ovrlp,nlpspd,proj,coeff_proj,ldiis,nit,nItInnerLoop,newgradient,&
          orthpar,confdatarr,methTransformOverlap,blocksize_pdgemm,convCrit,&
          hx,hy,hz,SIC,nItPrecond,factor_enlarge)


  end if


  if(updatePhi .or. itSCC==0) then
      if(newgradient) then

          ! Cancel the communication of the potential
          if(useDerivativeBasisFunctions) then
              call cancelCommunicationPotential(iproc, nproc, lbcomgp)
              call deallocateCommunicationsBuffersPotential(lbcomgp, subname)
          else
              call cancelCommunicationPotential(iproc, nproc, comgp)
              call deallocateCommunicationsBuffersPotential(comgp, subname)
          end if

          ! Reallocate lphiRestart, since its size might have changed
          iall=-product(shape(lphiRestart))*kind(lphiRestart)
          deallocate(lphiRestart, stat=istat)
          call memocc(istat, iall, 'lphiRestart', subname)
          allocate(lphiRestart(lorbs%npsidim_orbs), stat=istat)
          call memocc(istat, lphiRestart, 'lphiRestart', subname)

          ! Create new types for large basis...
          call deallocate_orbitals_data(llborbs, subname)
          call deallocate_overlapParameters(lbop, subname)
          call deallocate_p2pComms(lbcomon, subname)
          call deallocate_matrixDescriptors(lbmad, subname)
          call deallocate_p2pComms(lbcomgp, subname)


          call nullify_orbitals_data(llborbs)
          call nullify_overlapParameters(lbop)
          call nullify_p2pComms(lbcomon)
          call nullify_matrixDescriptors(lbmad)
          call nullify_p2pComms(lbcomgp)
          tag=1
          if(.not.useDerivativeBasisFunctions) then
              norbu=lorbs%norb
          else
              norbu=4*lorbs%norb
          end if
          norb=norbu
          norbd=0
          nspin=1
          call orbitals_descriptors_forLinear(iproc, nproc, norb, norbu, norbd, nspin, lorbs%nspinor,&
               lorbs%nkpts, lorbs%kpts, lorbs%kwgts, llborbs)
          call repartitionOrbitals(iproc, nproc, llborbs%norb, llborbs%norb_par,&
               llborbs%norbp, llborbs%isorb_par, llborbs%isorb, llborbs%onWhichMPI)

          allocate(orbsperlocreg(lzd%nlr), stat=istat)
          call memocc(istat, orbsperlocreg, 'orbsperlocreg', subname)
          !!orbsperlocreg=0
          !!do iorb=1,lorbs%norbp
          !!    iiorb=lorbs%isorb+iorb
          !!    ilr=iiorb
          !!    orbsperlocreg(ilr)=orbsperlocreg(ilr)+1
          !!end do
          !!call mpiallred(orbsperlocreg(1), lzdlarge%nlr, mpi_sum, mpi_comm_world, ierr)
          do iorb=1,lzd%nlr
              if(useDerivativeBasisFunctions) then
                  orbsperlocreg(iorb)=4
              else
                  orbsperlocreg(iorb)=1
              end if
          end do

          iall=-product(shape(llborbs%inWhichLocreg))*kind(llborbs%inWhichLocreg)
          deallocate(llborbs%inWhichLocreg, stat=istat)
          call memocc(istat, iall, 'llborbs%inWhichLocreg', subname)

          allocate(locregCenter(3,lzd%nlr), stat=istat)
          call memocc(istat, locregCenter, 'locregCenter', subname)
          !!locregCenter=0.d0
          !!do iorb=1,lorbs%norbp
          !!    iiorb=lorbs%isorb+iorb
          !!    ilr=lorbs%inwhichlocreg(iiorb)
          !!    !locregCenter(:,iiorb)=confdatarr(iorb)%rxyzConf
          !!    locregCenter(:,ilr)=confdatarr(iorb)%rxyzConf
          !!end do
          !!call mpiallred(locregCenter(1,1), 3*lorbs%norb, mpi_sum, mpi_comm_world, ierr)
          do ilr=1,lzd%nlr
              locregCenter(:,ilr)=lzd%llr(ilr)%locregCenter
          end do

          call assignToLocreg2(iproc, nproc, llborbs%norb, llborbs%norb_par, 0, lzd%nlr, &
               nspin, orbsperlocreg, locregCenter, llborbs%inwhichlocreg)

          ! Assign inwhichlocreg manually
          if(useDerivativeBasisFunctions) then
              norb=4
          else
              norb=1
          end if
          ii=0
          do iorb=1,lorbs%norb
              do i=1,norb
                  ii=ii+1
                  llborbs%inwhichlocreg(ii)=lorbs%inwhichlocreg(iorb)
              end do
          end do

          ! Recreate lzd, since it has to contain the bounds also for the derivatives
          ! First copy to some temporary structure
          allocate(locrad(lzd%nlr), stat=istat)
          call memocc(istat, locrad, 'locrad', subname)
          call nullify_locreg_descriptors(glr_tmp)
          call copy_locreg_descriptors(lzd%glr, glr_tmp, subname)
          nlr=lzd%nlr
          locrad=lzd%llr(:)%locrad
          call deallocate_local_zone_descriptors(lzd, subname)
          call nullify_local_zone_descriptors(lzd)
          call initLocregs(iproc, nproc, nlr, locregCenter, hx, hy, hz, lzd, lorbs, glr_tmp, locrad, 's', llborbs)
          call nullify_locreg_descriptors(lzd%glr)
          call copy_locreg_descriptors(glr_tmp, lzd%glr, subname)
          call deallocate_locreg_descriptors(glr_tmp, subname)
          iall=-product(shape(locrad))*kind(locrad)
          deallocate(locrad, stat=istat)
          call memocc(istat, iall, 'locrad', subname)

          iall=-product(shape(locregCenter))*kind(locregCenter)
          deallocate(locregCenter, stat=istat)
          call memocc(istat, iall, 'locregCenter', subname)
          iall=-product(shape(orbsperlocreg))*kind(orbsperlocreg)
          deallocate(orbsperlocreg, stat=istat)
          call memocc(istat, iall, 'orbsperlocreg', subname)

          npsidim = 0
          do iorb=1,llborbs%norbp
           ilr=llborbs%inwhichlocreg(iorb+llborbs%isorb)
           npsidim = npsidim + lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
          end do
          allocate(llborbs%eval(llborbs%norb), stat=istat)
          call memocc(istat, llborbs%eval, 'llborbs%eval', subname)
          llborbs%eval=-.5d0
          llborbs%npsidim_orbs=max(npsidim,1)
          call initCommsOrtho(iproc, nproc, nspin, hx, hy, hz, lzd, llborbs, llborbs%inWhichLocreg,&
               's', lbop, lbcomon, tag)
          call initMatrixCompression(iproc, nproc, lzd%nlr, llborbs, &
               lbop%noverlaps, lbop%overlaps, lbmad)
          call initCompressedMatmul3(llborbs%norb, lbmad)

          call initializeCommunicationPotential(iproc, nproc, denspot%dpcom%nscatterarr, llborbs, &
               lzd, lbcomgp, llborbs%inWhichLocreg, tag)

          ! Communicate the potential
          if(useDerivativeBasisFunctions) then
              call allocateCommunicationsBuffersPotential(lbcomgp, subname)
              call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lbcomgp)
          else
              call allocateCommunicationsBuffersPotential(comgp, subname)
              call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, comgp)
          end if


          tag=1
          call deallocateCommunicationbufferSumrho(comsr, subname)
          call deallocate_p2pComms(comsr, subname)
          call nullify_p2pComms(comsr)
          call initializeCommsSumrho(iproc, nproc, denspot%dpcom%nscatterarr, lzd, llborbs, tag, comsr)
          call allocateCommunicationbufferSumrho(iproc, .false., comsr, subname)



      end if

  end if

  ! Calculate the derivative basis functions. Copy the trace minimizing orbitals to lin%lphiRestart.
  ! Keep the value of lphi for the next iteration
  if(updatePhi .or. itSCC==0) call dcopy(lorbs%npsidim_orbs, lphi(1), 1, lphiRestart(1), 1)

  if(updatePhi .and. newgradient) then

      ! Reallocate lphi, since it is now allocated without the derivatives
      iall=-product(shape(lphi))*kind(lphi)
      deallocate(lphi, stat=istat)
      call memocc(istat, iall, 'lphi', subname)

      allocate(lphi(llborbs%npsidim_orbs), stat=istat)
      call memocc(istat, lphi, 'lphi', subname)

      if(.not.useDerivativeBasisFunctions) call dcopy(lorbs%npsidim_orbs, lphiRestart(1), 1, lphi(1), 1)
  end if
      
  !!call dcopy(lorbs%npsidim_orbs, lphiRestart(1), 1, lphi(1), 1)


  if(useDerivativeBasisFunctions .and. (updatePhi .or. itSCC==0)) then
      !!call dcopy(max(lorbs%npsidim_orbs,lorbs%npsidim_comp),lphi(1),1,lin%lphiRestart(1),1)


      call deallocate_p2pCommsRepartition(comrp, subname)
      call nullify_p2pCommsRepartition(comrp)
      call initializeRepartitionOrbitals(iproc, nproc, tag, lorbs, llborbs, lzd, comrp)

      if(iproc==0) write(*,'(1x,a)',advance='no') 'calculating derivative basis functions...'
      call getDerivativeBasisFunctions(iproc,nproc,hx,lzd,lorbs,llborbs,comrp,&
           max(lorbs%npsidim_orbs,lorbs%npsidim_comp),lphiRestart,lphi)
      if(iproc==0) write(*,'(a)') 'done.'
  end if



  ! Get the overlap matrix.
  if(.not.useDerivativeBasisFunctions) then
      call getOverlapMatrix2(iproc, nproc, lzd, lorbs, comon, op, lphi, mad, ovrlp)
  else
      call getOverlapMatrix2(iproc, nproc, lzd, llborbs, lbcomon, lbop, lphi, lbmad, ovrlp)
  end if


  if(communicate_lphi) then
      ! Allocate the communication buffers for the calculation of the charge density.
      !call allocateCommunicationbufferSumrho(iproc, comsr, subname)
      ! Transform all orbitals to real space.
      ist=1
      istr=1
      do iorb=1,llborbs%norbp
          !ilr=llborbs%inWhichLocregp(iorb)
          iiorb=llborbs%isorb+iorb
          ilr=llborbs%inWhichLocreg(iiorb)
          call initialize_work_arrays_sumrho(lzd%Llr(ilr), w)
          !!write(*,'(a,5i12)') 'iorb, ist, istr, size(lphi), size(comsr%sendBuf)', iorb, ist, istr, size(lphi), size(comsr%sendBuf)
          call daub_to_isf(lzd%Llr(ilr), w, lphi(ist), comsr%sendBuf(istr))
          call deallocate_work_arrays_sumrho(w)
          ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
          istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
      end do
      if(istr/=comsr%nsendBuf+1) then
          write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=comsr%nsendBuf+1'
          stop
      end if
      
      ! Post the MPI messages for the communication of sumrho. Since we use non blocking point
      ! to point communication, the program will continue immediately. The messages will be gathered
      ! in the subroutine sumrhoForLocalizedBasis2.
      call postCommunicationSumrho2(iproc, nproc, comsr, comsr%sendBuf, comsr%recvBuf)
  end if
  

  if(iproc==0) write(*,'(1x,a)') '----------------------------------- Determination of the orbitals in this new basis.'

  ! Gather the potential (it has been posted in the subroutine linearScaling) if the basis functions
  ! have not been updated (in that case it was gathered there). If newgradient is true, it has to be
  ! gathered as well since the locregs changed.
  if(.not.updatePhi .or. newgradient) then
      call gatherPotential(iproc, nproc, comgp)
  end if
  ! If we use the derivative basis functions the potential has to be gathered anyway.
  if(useDerivativeBasisFunctions) call gatherPotential(iproc, nproc, lbcomgp)


  if(.not.useDerivativeBasisFunctions) then


     call local_potential_dimensions(lzd,lorbs,denspot%dpcom%ngatherarr(0,1))

     call full_local_potential(iproc,nproc,lorbs,Lzd,2,denspot%dpcom,denspot%rhov,denspot%pot_full,comgp)
     !call full_local_potential(iproc,nproc,&
     !     lzd%glr%d%n1i*lzd%glr%d%n2i*nscatterarr(iproc,2),&
     !     lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i,nspin,&
     !     lzd%glr%d%n1i*lzd%glr%d%n2i*nscatterarr(iproc,1)*nspin,0,&
     !     lorbs,lzd,2,ngatherarr,rhopot,lpot,comgp)
  else

     call local_potential_dimensions(lzd,llborbs,denspot%dpcom%ngatherarr(0,1))
      
     call full_local_potential(iproc,nproc,llborbs,Lzd,2,denspot%dpcom,denspot%rhov,denspot%pot_full,lbcomgp)
     !call full_local_potential(iproc,nproc,&
     !     lzd%glr%d%n1i*lzd%glr%d%n2i*nscatterarr(iproc,2),&
     !     lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i,nspin,&
     !     lzd%glr%d%n1i*lzd%glr%d%n2i*nscatterarr(iproc,1)*nspin,0,&
     !     llborbs,lzd,2,ngatherarr,rhopot,lpot,lbcomgp)
     
  end if

  ! Apply the Hamitonian to the orbitals. The flag withConfinement=.false. indicates that there is no
  ! confining potential added to the Hamiltonian.
  allocate(lhphi(max(llborbs%npsidim_orbs,llborbs%npsidim_comp)), stat=istat)
  call memocc(istat, lhphi, 'lhphi', subname)
  withConfinement=.false.
  !if(iproc==0) write(*,'(1x,a)',advance='no') 'Hamiltonian application...'
  allocate(lzd%doHamAppl(lzd%nlr), stat=istat)
  call memocc(istat, lzd%doHamAppl, 'lzd%doHamAppl', subname)
  lzd%doHamAppl=.true.
  if(.not.useDerivativeBasisFunctions) then
     allocate(confdatarrtmp(lorbs%norbp))
     call default_confinement_data(confdatarrtmp,lorbs%norbp)
     call FullHamiltonianApplication(iproc,nproc,at,lorbs,&
          hx,hy,hz,rxyz,&
          proj,lzd,nlpspd,confdatarrtmp,denspot%dpcom%ngatherarr,denspot%pot_full,lphi,lhphi,&
          ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,SIC,GPU,&
          pkernel=denspot%pkernelseq)
     deallocate(confdatarrtmp)

  else

     allocate(confdatarrtmp(llborbs%norbp))
     call default_confinement_data(confdatarrtmp,llborbs%norbp)
     call FullHamiltonianApplication(iproc,nproc,at,llborbs,&
          hx,hy,hz,rxyz,&
          proj,lzd,nlpspd,confdatarrtmp,denspot%dpcom%ngatherarr,denspot%pot_full,lphi,lhphi,&
          ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,SIC,GPU,&
          pkernel=denspot%pkernelseq)
     deallocate(confdatarrtmp)
  end if

  iall=-product(shape(lzd%doHamAppl))*kind(lzd%doHamAppl)
  deallocate(lzd%doHamAppl, stat=istat)
  call memocc(istat, iall, 'lzd%doHamAppl', subname)



  iall=-product(shape(denspot%pot_full))*kind(denspot%pot_full)
  deallocate(denspot%pot_full, stat=istat)
  call memocc(istat, iall, 'denspot%pot_full', subname)

  if(iproc==0) write(*,'(1x,a)') 'done.'

  ! Deallocate the buffers needed for the communication of the potential.
  call deallocateCommunicationsBuffersPotential(comgp, subname)
  if(useDerivativeBasisFunctions) call deallocateCommunicationsBuffersPotential(lbcomgp, subname)



  ! Calculate the matrix elements <phi|H|phi>.
  call allocateCommuncationBuffersOrtho(lbcomon, subname)
  if(.not. useDerivativeBasisFunctions) then
      call getMatrixElements2(iproc, nproc, lzd, llborbs, lbop, lbcomon, lphi, lhphi, mad, matrixElements)
      !call getMatrixElements2(iproc, nproc, lzd, lorbs, op, comon, lphi, lhphi, mad, matrixElements)
  else
      call getMatrixElements2(iproc, nproc, lzd, llborbs, lbop, lbcomon, lphi, lhphi, lbmad, matrixElements)
  end if
  call deallocateCommuncationBuffersOrtho(lbcomon, subname)


  !!!! TEST ########################################################
  !!call transpose_linear(iproc, 0, nproc-1, llborbs, lin%lb%collComms, lphi, mpi_comm_world, phiWork)
  !!call transpose_linear(iproc, 0, nproc-1, llborbs, lin%lb%collComms, lhphi, mpi_comm_world, phiWork)
  !!if(iproc==0) then
  !!    do ierr=1,orbs%npsidim
  !!        write(53,*) ierr, lin%lb%collComms%indexarray(ierr)
  !!    end do
  !!end if
  !!call calculate_overlap_matrix(iproc, llborbs, lin%lb%collComms, lphi, lhphi, matrixElements(1,1,2))
  !!call untranspose_linear(iproc, 0, nproc-1,  llborbs, lin%lb%collComms, lhphi, mpi_comm_world, phiWork)
  !!call untranspose_linear(iproc, 0, nproc-1,  llborbs, lin%lb%collComms, lphi, mpi_comm_world, phiWork)
  !!!! END TEST ####################################################
  !!do iorb=1,llborbs%norb
  !!    do jorb=1,llborbs%norb
  !!        write(400+iproc,*) iorb, jorb, matrixElements(jorb,iorb,1)
  !!        write(410+iproc,*) iorb, jorb, matrixElements(jorb,iorb,2)
  !!    end do
  !!end do

  !!tt=0.d0
  !!do iall=1,llborbs%norb
  !!     tt=tt+matrixElements(iall,iall,1)
  !!end do
  !!if(iproc==0) write(*,*) 'trace of H without confinement:',tt


  allocate(overlapmatrix(llborbs%norb,llborbs%norb), stat=istat)
  call memocc(istat, overlapmatrix, 'overlapmatrix', subname)
  overlapmatrix=ovrlp
  

  ! Diagonalize the Hamiltonian, either iteratively or with lapack.
  !!call mpi_barrier(mpi_comm_world, ierr) !To measure the time correctly.
  !!t1=mpi_wtime()
  ! Make a copy of the matrix elements since dsyev overwrites the matrix and the matrix elements
  ! are still needed later.
  call dcopy(llborbs%norb**2, matrixElements(1,1,1), 1, matrixElements(1,1,2), 1)
  if(blocksize_pdsyev<0) then
      if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, sequential version... '
      call diagonalizeHamiltonian2(iproc, nproc, llborbs, lbop%nsubmax, matrixElements(1,1,2), ovrlp, eval)
  else
      if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, parallel version... '
      call dsygv_parallel(iproc, nproc, blocksize_pdsyev, nproc_pdsyev, mpi_comm_world, 1, 'v', 'l', llborbs%norb,&
           matrixElements(1,1,2), llborbs%norb, ovrlp, llborbs%norb, eval, info)
  end if
  if(iproc==0) write(*,'(a)') 'done.'
  call dcopy(llborbs%norb*orbs%norb, matrixElements(1,1,2), 1, coeff(1,1), 1)
  !if(.not.updatePhi) call dcopy(llborbs%norb*(orbs%norb+lin%norbvirt), matrixElements(1,1,2), 1, lin%coeffall(1,1), 1)
  infoCoeff=0

  ! Write some eigenvalues. Don't write all, but only a few around the last occupied orbital.
  if(iproc==0) then
      write(*,'(1x,a)') '-------------------------------------------------'
      write(*,'(1x,a)') 'some selected eigenvalues:'
      do iorb=max(orbs%norb-8,1),min(orbs%norb+8,lorbs%norb)
          if(iorb==orbs%norb) then
              write(*,'(3x,a,i0,a,es12.5,a)') 'eval(',iorb,')=',eval(iorb),'  <-- last occupied orbital'
          else if(iorb==orbs%norb+1) then
              write(*,'(3x,a,i0,a,es12.5,a)') 'eval(',iorb,')=',eval(iorb),'  <-- first virtual orbital'
          else
              write(*,'(3x,a,i0,a,es12.5)') 'eval(',iorb,')=',eval(iorb)
          end if
      end do
      write(*,'(1x,a)') '-------------------------------------------------'
  end if
  !!t2=mpi_wtime()
  !!time=t2-t1
  !!if(iproc==0) write(*,'(1x,a,es10.3)') 'time for diagonalizing the Hamiltonian:',time


  ! Calculate the band structure energy with matrixElements.
  ebs=0.d0
  do iorb=1,orbs%norb
      do jorb=1,llborbs%norb
          do korb=1,llborbs%norb
              ebs = ebs + coeff(jorb,iorb)*coeff(korb,iorb)*matrixElements(korb,jorb,1)
          end do
      end do
  end do
  ! If closed shell multiply by two.
  if(orbs%nspin==1) ebs=2.d0*ebs



  ! Project the lb coefficients on the smaller subset
  if(useDerivativeBasisFunctions) then
      inc=4
  else
      inc=1
  end if
  do iorb=1,orbs%norb
      jjorb=1
      do jorb=1,llborbs%norb,inc
          tt=0.d0
          do korb=1,llborbs%norb
              tt = tt + coeff(korb,iorb)*overlapmatrix(korb,jorb)
          end do
          coeff_proj(jjorb,iorb)=tt
          !if(iproc==0) write(99,'(2i7,2es16.8)') iorb, jjorb,  coeff_proj(jjorb,iorb), coeff(jorb,iorb)
          jjorb=jjorb+1
      end do
  end do

  

  !!! Copy the basis functions for the next iterations
  !!call dcopy(max(lorbs%npsidim_orbs,lorbs%npsidim_comp), lphi(1), 1, lin%lphiold(1), 1)

  !!! Copy the Hamiltonian matrix for the next iteration
  !!call dcopy(lorbs%norb**2, matrixElements(1,1,1), 1, lin%hamold(1,1), 1)

  ! Deallocate all local arrays.
  iall=-product(shape(HamSmall))*kind(HamSmall)
  deallocate(HamSmall, stat=istat)
  call memocc(istat, iall, 'HamSmall', subname)

  iall=-product(shape(lhphi))*kind(lhphi)
  deallocate(lhphi, stat=istat)
  call memocc(istat, iall, 'lhphi', subname)

  iall=-product(shape(matrixElements))*kind(matrixElements)
  deallocate(matrixElements, stat=istat)
  call memocc(istat, iall, 'matrixElements', subname)
  
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)

  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)

  iall=-product(shape(overlapmatrix))*kind(overlapmatrix)
  deallocate(overlapmatrix, stat=istat)
  call memocc(istat, iall, 'overlapmatrix', subname)

end subroutine getLinearPsi


subroutine getLocalizedBasis(iproc,nproc,at,lzd,lorbs,orbs,comon,op,comgp,mad,rxyz,&
    denspot,GPU,lphi,trH,&
    infoBasisFunctions,ovrlp,nlpspd,proj,coeff,ldiis,nit,nItInnerLoop,newgradient,orthpar,&
    confdatarr,methTransformOverlap,blocksize_pdgemm,convCrit,hx,hy,hz,SIC,nItPrecond,factor_enlarge)
!
! Purpose:
! ========
!   Calculates the localized basis functions phi. These basis functions are obtained by adding a
!   quartic potential centered on the atoms to the ordinary Hamiltonian. The eigenfunctions are then
!   determined by minimizing the trace until the gradient norm is below the convergence criterion.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc           process ID
!     nproc           total number of processes
!     at              type containing the paraneters for the atoms
!     orbs            type describing the physical orbitals psi
!     Glr             type describing the localization region
!     input           type containing some very general parameters
!     lin             type containing parameters for the linear version
!     rxyz            the atomic positions
!     nspin           npsin==1 -> closed shell; npsin==2 -> spin polarized
!     nlpsp           ???
!     nscatterarr     ???
!     ngatherarr      ???
!     rhopot          the charge density
!     GPU             parameters for GPUs
!     pkernelseq      ???
!     rxyzParab       the center of the confinement potential (at the moment identical rxyz)
!     n3p             ???
!     itSCC           iteration in the self consistency cycle
!  Input/Output arguments
!  ---------------------
!     phi             the localized basis functions. It is assumed that they have been initialized
!                     somewhere else
!   Output arguments
!   ----------------
!     hphi            the modified Hamiltonian applied to phi
!     trH             the trace of the Hamiltonian
!     infoBasisFunctions  indicates wheter the basis functions converged to the specified limit (value is 0)
!                         or whether the iteration stopped due to the iteration limit (value is -1). This info
!                         is returned by 'getLocalizedBasis'
!
! Calling arguments:
!   Input arguments
!   Output arguments
!    phi   the localized basis functions
!
use module_base
use module_types
use module_interfaces, except_this_one => getLocalizedBasis, except_this_one_A => writeonewave
!  use Poisson_Solver
!use allocModule
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nit, nItInnerLoop, methTransformOverlap, blocksize_pdgemm
integer,intent(in):: nItPrecond
integer,intent(out):: infoBasisFunctions
type(atoms_data), intent(in) :: at
type(local_zone_descriptors),intent(inout):: lzd
type(orbitals_data):: lorbs, orbs
type(p2pComms):: comon
type(overlapParameters):: op
!type(p2pCommsGatherPot):: comgp
type(p2pComms):: comgp
type(matrixDescriptors),intent(inout):: mad
real(8),dimension(3,at%nat):: rxyz
type(DFT_local_fields), intent(inout) :: denspot
!integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
!real(dp), dimension(*), intent(inout) :: rhopot
type(GPU_pointers), intent(inout) :: GPU
!real(dp), dimension(:), pointer :: pkernelseq
!real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)):: lphi
real(8),dimension(:),pointer,intent(inout):: lphi
real(8),intent(out):: trH
real(8),intent(in):: convCrit, hx, hy, hz, factor_enlarge
real(8),dimension(lorbs%norb,lorbs%norb),intent(out):: ovrlp
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
real(8),dimension(lorbs%norb,orbs%norb),intent(in):: coeff
type(localizedDIISParameters),intent(inout):: ldiis
logical,intent(in):: newgradient
type(orthon_data),intent(in):: orthpar
type(confpot_data), dimension(lorbs%norbp),intent(inout) :: confdatarr
type(SIC_data) :: SIC !<parameters for the SIC methods

! Local variables
real(8) ::epot_sum,ekin_sum,eexctX,eproj_sum,eval_zero,t1tot,eSIC_DC
real(8) :: t2tot,timetot,tt1,tt2,tt3,tt4,tt5
real(8):: tt,ddot,fnrm,fnrmMax,meanAlpha,gnrm,gnrm_zero,gnrmMax,t1,t2, dnrm2
real(8) :: timecommunp2p, timeextract, timecommuncoll, timeoverlap, timecompress, energyconf_0, energyconf_trial
real(8):: trHold, factor, factor2
integer:: iorb, icountSDSatur, icountSwitch, idsx, icountDIISFailureTot, consecutive_rejections
integer :: icountDIISFailureCons,itBest, ncnt, lorb, ilrlarge2
integer:: istat,istart,ierr,ii,it,iall,ind1,ind2,jorb,ist,iiorb
integer:: gdim,ilr,ncount,offset,istsource,istdest,korb
integer,dimension(:),allocatable:: norbsPerAtom, inwhichlocreg_reference
real(8),dimension(:),allocatable:: alpha,fnrmOldArr,alphaDIIS
real(8),dimension(:,:),allocatable:: fnrmArr, fnrmOvrlpArr, lagmat, Umat
real(8),dimension(:,:),allocatable:: kernel, kernelold, locregCenter
logical:: withConfinement, resetDIIS, immediateSwitchToSD, variable_locregs
character(len=*),parameter:: subname='getLocalizedBasis'
real(8),dimension(5):: time
!real(8),dimension(:),pointer:: lpot
!real(8),external :: mpi_wtime1
character(len=3):: orbname, comment
real(8),dimension(:),allocatable:: lvphiovrlp, locrad, locrad_tmp
real(8),dimension(:),pointer:: phiWork
real(8),dimension(:),pointer:: lphilarge, lhphilarge, lhphilargeold, lphilargeold, lhphi, lhphiold, lphiold
real(8),dimension(:),pointer:: lphilarge2, lhphilarge2, lhphilarge2old, lphilarge2old
integer:: jst, istl, istg, nvctrp, ldim, nspin, norbu, norbd, tag, npsidim, ilrlarge, icenter, nl1, nl2, nl3
type(local_zone_descriptors):: lzdlarge, lzdlarge2
type(orbitals_data):: orbslarge, orbslarge2
type(overlapParameters):: oplarge, oplarge2
type(p2pComms):: comonlarge, comonlarge2
type(p2pComms):: comgplarge, comgplarge2
type(matrixDescriptors):: madlarge, madlarge2
type(localizedDIISParameters):: ldiis2
integer :: unitwf !DEBUG
character(len=10):: c1, c2, c3 !DEBUG
character(len=50):: filename !DEBUG
logical,parameter:: secondLocreg=.false.

! automatic array, for debugging
real(8),dimension(3,lzd%nlr):: locregCenterTemp


  ! Allocate all local arrays.
  call allocateLocalArrays()

  allocate(kernel(lorbs%norb,lorbs%norb), stat=istat)
  call memocc(istat, kernel, 'kernel', subname)

  call dgemm('n', 't', lorbs%norb, lorbs%norb, orbs%norb, 1.d0, coeff(1,1), lorbs%norb, &
       coeff(1,1), lorbs%norb, 0.d0, kernel(1,1), lorbs%norb)

  !!if(iproc==0) write(*,*) 'WARNING: MODIFYING KERNEL'
  !!do iorb=1,lorbs%norb
  !!    do jorb=1,lorbs%norb
  !!        if(jorb==iorb) then
  !!            kernel(jorb,iorb)=1.d0
  !!        else
  !!            kernel(jorb,iorb)=0.d0
  !!        end if
  !!    end do
  !!end do

  
  if(iproc==0) write(*,'(1x,a)') '======================== Creation of the basis functions... ========================'

  if(nItInnerLoop==-1 .and. factor_enlarge==1.d0) then
      variable_locregs=.false.
  else
      variable_locregs=.true.
  end if

  ! Initialize the arrays and variable needed for DIIS.
  if(newgradient .and. ldiis%isx>0) then
      if(variable_locregs) then
          if(iproc==0) write(*,'(1x,a)') 'ERROR: if newgradient is true, only steepest descent is &
                                          &allowed since the locreg shapes may change!'
          call mpi_barrier(mpi_comm_world, ierr)
          stop
      end if
  end if
  icountSDSatur=0
  icountSwitch=0
  icountDIISFailureTot=0
  icountDIISFailureCons=0
  ldiis%is=0
  ldiis%switchSD=.false.
  ldiis%trmin=1.d100
  ldiis%trold=1.d100
  alpha=ldiis%alphaSD
  alphaDIIS=ldiis%alphaDIIS

  ! Copy parameters to ldiis2 (needed for debugging)
  if(secondLocreg) then
      ldiis2%is=ldiis%is
      ldiis2%isx=ldiis%isx
      ldiis2%mis=ldiis%mis
      ldiis2%DIISHistMax=ldiis%DIISHistMax
      ldiis2%DIISHistMin=ldiis%DIISHistMin
      ldiis2%trmin=ldiis%trmin
      ldiis2%trold=ldiis%trold
      ldiis2%alphaSD=ldiis%alphaSD
      ldiis2%alphaDIIS=ldiis%alphaDIIS
      ldiis2%switchSD=ldiis%switchSD
      allocate(ldiis2%phiHist(1), stat=istat)
      call memocc(istat, ldiis2%phiHist, 'ldiis2%phiHist', subname)
      allocate(ldiis2%hphiHist(1), stat=istat)
      call memocc(istat, ldiis2%hphiHist, 'ldiis2%hphiHist', subname)
  end if



  if(.not.newgradient) then
      ! Gather the potential that each process needs for the Hamiltonian application for all its orbitals.
      ! The messages for this point ', to point communication have been posted in the subroutine linearScaling.
      call gatherPotential(iproc, nproc, comgp)

      ! Build the required potential
      call local_potential_dimensions(lzd,lorbs,denspot%dpcom%ngatherarr(0,1))

      call full_local_potential(iproc,nproc,lorbs,Lzd,2,denspot%dpcom,denspot%rhov,denspot%pot_full,comgp)
  end if


  allocate(lphiold(size(lphi)), stat=istat)
  call memocc(istat, lphiold, 'lphiold', subname)
  allocate(lagmat(lorbs%norb,lorbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat', subname)
  allocate(Umat(lorbs%norb,lorbs%norb), stat=istat)
  call memocc(istat, Umat, 'Umat', subname)
  allocate(kernelold(lorbs%norb,lorbs%norb), stat=istat)
  call memocc(istat, kernelold, 'kernelold', subname)
  allocate(locregCenter(3,lzd%nlr), stat=istat)
  call memocc(istat, locregCenter, 'locregCenter', subname)
  allocate(locrad(lzd%nlr), stat=istat)
  call memocc(istat, locrad, 'locrad', subname)
  allocate(locrad_tmp(lzd%nlr), stat=istat)
  call memocc(istat, locrad_tmp, 'locrad_tmp', subname)
  allocate(inwhichlocreg_reference(lorbs%norb), stat=istat)
  call memocc(istat, inwhichlocreg_reference, 'inwhichlocreg_reference', subname)

  time=0.d0
  resetDIIS=.false.
  immediateSwitchToSD=.false.
  t1tot=mpi_wtime()
  consecutive_rejections=0
  trHold=1.d100
 
  do ilr=1,lzd%nlr
      locrad(ilr)=lzd%llr(ilr)%locrad
  end do



  ! ratio of large locreg and standard locreg
  !factor=1.5d0
  factor=factor_enlarge
  factor2=200.0d0

  ! always use the same inwhichlocreg
  inwhichlocreg_reference = lorbs%inwhichlocreg


  ! Initialize largestructures if required
  if(newgradient) then
      do iorb=1,lorbs%norb
          ilr=lorbs%inwhichlocreg(iorb)
          locregCenter(:,ilr)=lzd%llr(ilr)%locregCenter
          !!if(iproc==0) write(*,'(a,i6,3f11.4)') 'ilr, locregCenter(:,ilr)', ilr, locregCenter(:,ilr)
      end do
      locregCenterTemp=locregCenter

      locrad_tmp=factor*locrad
      call create_new_locregs(iproc, nproc, lzd%nlr, hx, hy, hz, lorbs, lzd%glr, locregCenter, &
           locrad_tmp, denspot%dpcom%nscatterarr, .false., inwhichlocreg_reference, ldiis, &
           lzdlarge, orbslarge, oplarge, comonlarge, madlarge, comgplarge, &
           lphilarge, lhphilarge, lhphilargeold, lphilargeold)

      if(secondLocreg) then
          locrad_tmp=factor2*locrad
          call create_new_locregs(iproc, nproc, lzd%nlr, hx, hy, hz, lorbs, lzd%glr, locregCenter, &
               locrad_tmp, denspot%dpcom%nscatterarr, .false., inwhichlocreg_reference, ldiis2, &
               lzdlarge2, orbslarge2, oplarge2, comonlarge2, madlarge2, comgplarge2, &
               lphilarge2, lhphilarge2, lhphilarge2old, lphilarge2old)
      end if
  end if



  iterLoop: do it=1,nit






      fnrmMax=0.d0
      fnrm=0.d0
  
      if (iproc==0) then
          write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
      endif


      ! Orthonormalize the orbitals. If the localization regions are smaller that the global box (which
      ! will be the usual case), the orthogonalization can not be done exactly, but only approximately.
      if(iproc==0) then
          write(*,'(1x,a)',advance='no') 'Orthonormalization...'
      end if
      t1=mpi_wtime()
      if(.not.ldiis%switchSD) then
          !if(iproc==0) write(*,*) 'WARNING: COMMENTED'
          if(.not.newgradient) then
             call orthonormalizeLocalized(iproc, nproc, &
             orthpar%methTransformOverlap, orthpar%nItOrtho, &
             orthpar%blocksize_pdsyev, orthpar%blocksize_pdgemm, lorbs, op, comon, lzd, &
             mad, lphi, ovrlp)

         else

            lphilarge=0.d0
            istl=1
            istg=1
            do iorb=1,orbslarge%norbp
                ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
                ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
                ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                gdim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                nspin=1
                call Lpsi_to_global2(iproc, nproc, ldim, gdim, lorbs%norb, lorbs%nspinor, nspin, lzdlarge%llr(ilrlarge), &
                     lzd%llr(ilr), lphi(istl), lphilarge(istg))
                istl=istl+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                istg=istg+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
            end do
            call orthonormalizeLocalized(iproc, nproc, &
                 orthpar%methTransformOverlap, orthpar%nItOrtho, &
                 orthpar%blocksize_pdsyev, orthpar%blocksize_pdgemm, orbslarge, oplarge, comonlarge, lzdlarge, &
                 madlarge, lphilarge, ovrlp)

            !!!! Plot basis functions
            !!call plotOrbitals(iproc, lorbs, lzd%Glr, lphilarge, lzd%nlr, locregCenter, lorbs%inwhichlocreg, .5d0*hx, &
            !!    .5d0*hy, .5d0*hz, it)
            !!! plot the orbitals -- EXPERIMENTAL ##################################################
            !!allocate(lvphiovrlp(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f))
            !!ist=1
            !!write(comment,'(i3.3)') it
            !!do iorb=1,orbslarge%norbp
            !!    iiorb=iorb+orbslarge%isorb
            !!    ilr=orbslarge%inwhichlocreg(iiorb)
            !!    write(orbname,'(i3.3)') iiorb
            !!    write(*,'(a,i0)') 'plotting orbital ',iiorb
            !!    lvphiovrlp=0.d0
            !!    call Lpsi_to_global2(iproc, nproc, lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f, &
            !!         lzdlarge%glr%wfd%nvctr_c+7*lzdlarge%glr%wfd%nvctr_f, orbslarge%norb, orbslarge%nspinor, nspin, &
            !!         lzdlarge%Glr, lzdlarge%Llr(ilr), lphilarge(ist), lvphiovrlp(1))
            !!    call plot_wf(orbname//'_'//comment, 2, at, 1.d0, lzdlarge%glr, hx, hx, hx, rxyz, lvphiovrlp(1))
            !!    ist=ist+lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
            !!end do
            !!deallocate(lvphiovrlp)
            !!! ####################################################################################

            if(secondLocreg) then
                ! Go to even larger region
                lphilarge2=0.d0
                istl=1
                istg=1
                do iorb=1,orbslarge%norbp
                    ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
                    ilrlarge2 = orbslarge2%inWhichLocreg(orbslarge2%isorb+iorb)
                    ldim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                    gdim=lzdlarge2%llr(ilrlarge2)%wfd%nvctr_c+7*lzdlarge2%llr(ilrlarge2)%wfd%nvctr_f
                    nspin=1
                    call Lpsi_to_global2(iproc, nproc, ldim, gdim, lorbs%norb, lorbs%nspinor, nspin, lzdlarge2%llr(ilrlarge2), &
                         lzdlarge%llr(ilrlarge), lphilarge(istl), lphilarge2(istg))
                    istl=istl+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                    istg=istg+lzdlarge2%llr(ilrlarge2)%wfd%nvctr_c+7*lzdlarge2%llr(ilrlarge2)%wfd%nvctr_f
                end do
                if(istl/=orbslarge%npsidim_orbs+1) stop 'istl/=orbslarge%npsidim_orbs+1'
                if(istg/=orbslarge2%npsidim_orbs+1) stop 'istg/=orbslarge2%npsidim_orbs+1'
            end if

            if(.not.secondLocreg) then
                ! Update confdatarr...
                do iorb=1,orbslarge%norbp
                   iiorb=orbslarge%isorb+iorb
                   ilr=orbslarge%inWhichlocreg(iiorb)
                   icenter=orbslarge%inwhichlocreg(iiorb)
                   !confdatarr(iorb)%potorder=lin%confpotorder
                   !confdatarr(iorb)%prefac=lin%potentialprefac(at%iatype(icenter))
                   !confdatarr(iorb)%hh(1)=.5_gp*hx
                   !confdatarr(iorb)%hh(2)=.5_gp*hy
                   !confdatarr(iorb)%hh(3)=.5_gp*hz
                   confdatarr(iorb)%rxyzConf(1:3)=locregCenterTemp(1:3,icenter)
                   call my_geocode_buffers(lzdlarge%Llr(ilr)%geocode,nl1,nl2,nl3)
                   confdatarr(iorb)%ioffset(1)=lzdlarge%llr(ilr)%nsi1-nl1-1
                   confdatarr(iorb)%ioffset(2)=lzdlarge%llr(ilr)%nsi2-nl2-1
                   confdatarr(iorb)%ioffset(3)=lzdlarge%llr(ilr)%nsi3-nl3-1
                end do
                call MLWFnew(iproc, nproc, lzdlarge, orbslarge, at, oplarge, &
                     comonlarge, madlarge, rxyz, nItInnerLoop, kernel, &
                     newgradient, confdatarr, hx, locregCenterTemp, 3.d0, lphilarge, Umat, locregCenter)
            end if


            if(secondLocreg) then
                ! Update confdatarr...
                do iorb=1,orbslarge%norbp
                   iiorb=orbslarge%isorb+iorb
                   ilr=orbslarge%inWhichlocreg(iiorb)
                   icenter=orbslarge%inwhichlocreg(iiorb)
                   !confdatarr(iorb)%potorder=lin%confpotorder
                   !confdatarr(iorb)%prefac=lin%potentialprefac(at%iatype(icenter))
                   !confdatarr(iorb)%hh(1)=.5_gp*hx
                   !confdatarr(iorb)%hh(2)=.5_gp*hy
                   !confdatarr(iorb)%hh(3)=.5_gp*hz
                   confdatarr(iorb)%rxyzConf(1:3)=locregCenterTemp(1:3,icenter)
                   call my_geocode_buffers(lzdlarge%Llr(ilr)%geocode,nl1,nl2,nl3)
                   confdatarr(iorb)%ioffset(1)=lzdlarge%llr(ilr)%nsi1-nl1-1
                   confdatarr(iorb)%ioffset(2)=lzdlarge%llr(ilr)%nsi2-nl2-1
                   confdatarr(iorb)%ioffset(3)=lzdlarge%llr(ilr)%nsi3-nl3-1
                end do
                !!!!call random_number(lphilarge)

                !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
                !!     lzdlarge%llr(ilrlarge)%locrad/factor, confdatarr, lphilarge)
                call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
                     lzdlarge%llr(1)%locrad/factor2, confdatarr, lphilarge)
                !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
                !!     6.d0, confdatarr, lphilarge)
                call MLWFnew(iproc, nproc, lzdlarge, orbslarge, at, oplarge, &
                     comonlarge, madlarge, rxyz, nItInnerLoop, kernel, &
                     newgradient, confdatarr, hx, locregCenterTemp, 3.d0, lphilarge, Umat, locregCenter)
                call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
                     lzdlarge%llr(1)%locrad/factor2, confdatarr, lphilarge)
                !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
                !!     6.d0, confdatarr, lphilarge)
            end if

            !!     !!write(*,*) "ATTENTION HERE!!!"
            !!call plotOrbitals(iproc, lorbs, lzd%Glr, lphilarge, lzd%nlr, locregCenter, lorbs%inwhichlocreg, .5d0*hx, &
            !!    .5d0*hy, .5d0*hz, 100+it)
            !!! plot the orbitals -- EXPERIMENTAL ##################################################
            !!allocate(lvphiovrlp(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f))
            !!ist=1
            !!write(comment,'(i3.3)') 100+it
            !!do iorb=1,orbslarge%norbp
            !!    iiorb=iorb+orbslarge%isorb
            !!    ilr=orbslarge%inwhichlocreg(iiorb)
            !!    write(orbname,'(i3.3)') iiorb
            !!    write(*,'(a,i0)') 'plotting orbital ',iiorb
            !!    lvphiovrlp=0.d0
            !!    call Lpsi_to_global2(iproc, nproc, lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f, &
            !!         lzdlarge%glr%wfd%nvctr_c+7*lzdlarge%glr%wfd%nvctr_f, orbslarge%norb, orbslarge%nspinor, nspin, &
            !!         lzdlarge%Glr, lzdlarge%Llr(ilr), lphilarge(ist), lvphiovrlp(1))
            !!    call plot_wf(orbname//'_'//comment, 2, at, 1.d0, lzdlarge%glr, hx, hx, hx, rxyz, lvphiovrlp(1))
            !!    ist=ist+lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
            !!end do
            !!deallocate(lvphiovrlp)
            !!! ####################################################################################

                 !!locregCenter=locregCenterTemp
            do ilr=1,lzd%nlr
                tt = (locregCenter(1,ilr)-locregCenterTemp(1,ilr))**2 &
                    +(locregCenter(2,ilr)-locregCenterTemp(2,ilr))**2 &
                    +(locregCenter(3,ilr)-locregCenterTemp(3,ilr))**2 
                tt=sqrt(tt)
                !!if(iproc==0) write(*,'(a,i5,3es13.4)') 'ilr, new centers:', &
                !!     ilr, locregCenter(1,ilr), locregCenter(2,ilr), locregCenter(3,ilr)
                !!if(iproc==0) write(*,'(a,i5,es16.4)') '2: ilr, tt', ilr, tt
                if( floor(locregCenter(1,ilr)/hx) < 0 .or. ceiling(locregCenter(1,ilr)/hx) > lzd%glr%d%n1 ) then
                    if(iproc==0) then
                        write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
                            ' is outside of box in x direction! Box limits=',0,lzd%glr%d%n1,&
                            ', center=',floor(locregCenter(1,ilr)/hx),ceiling(locregCenter(1,ilr)/hx)
                    end if
                    call mpi_barrier(mpi_comm_world, ierr)
                    stop
                end if
                if( floor(locregCenter(2,ilr)/hx) < 0 .or. ceiling(locregCenter(2,ilr)/hy) > lzd%glr%d%n2 ) then
                    if(iproc==0) then
                        write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
                            'is outside of box in y direction! Box limits=',0,lzd%glr%d%n2,&
                            ', center=',floor(locregCenter(2,ilr)/hy),ceiling(locregCenter(2,ilr)/hy)
                    end if
                    call mpi_barrier(mpi_comm_world, ierr)
                    stop
                end if
                if( floor(locregCenter(3,ilr)/hz) < 0 .or. ceiling(locregCenter(3,ilr)/hz) > lzd%glr%d%n3 ) then
                    if(iproc==0) then
                        write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
                            'is outside of box in x direction! Box limits=',0,lzd%glr%d%n3,&
                            ', center=',floor(locregCenter(3,ilr)/hz),ceiling(locregCenter(3,ilr)/hz)
                    end if
                    call mpi_barrier(mpi_comm_world, ierr)
                    stop
                end if
            end do
            do ilr=1,lzd%nlr
                tt = (locregCenter(1,ilr)-locregCenterTemp(1,ilr))**2 &
                    +(locregCenter(1,ilr)-locregCenterTemp(1,ilr))**2 &
                    +(locregCenter(1,ilr)-locregCenterTemp(1,ilr))**2 
                tt=sqrt(tt)
                !!if(iproc==0) write(*,'(a,i5,es16.4)') '1: ilr, tt', ilr, tt
            end do

            !write(*,*) 'WARNING: CHECK HERE THE INDICES OF UMAT!!'
            if(nItInnerLoop>0) then                          
                 kernelold=kernel
                 do iorb=1,lorbs%norb
                     do jorb=1,lorbs%norb
                         tt=0.d0
                         do korb=1,lorbs%norb
                             do lorb=1,lorbs%norb
                                 tt=tt+kernelold(korb,lorb)*Umat(korb,iorb)*Umat(lorb,jorb)
                                 !tt=tt+kernelold(korb,lorb)*Umat(iorb,korb)*Umat(jorb,lorb)
                             end do
                         end do
                         kernel(jorb,iorb)=tt
                     end do
                 end do
             end if

            if(variable_locregs) then
                !!call deallocateCommunicationsBuffersPotential(comgp, subname)
                call destroy_new_locregs(lzd, lorbs, op, comon, mad, comgp, &
                     lphi, lhphi, lhphiold, lphiold)
                !!call mpi_barrier(mpi_comm_world, ierr)
                !!stop
                 !!do iorb=1,lorbs%norb
                 !!    write(*,'(a,2i8,3es12.4,4x,3es12.4)') 'iproc, iorb, locregCenterTemp(:,iorb), locregCenter(:,iorb)', &
                 !!        iproc, iorb, locregCenterTemp(:,iorb), locregCenter(:,iorb)
                 !!end do

                call create_new_locregs(iproc, nproc, lzdlarge%nlr, hx, hy, hz, orbslarge, lzdlarge%glr, locregCenter, &
                     locrad, denspot%dpcom%nscatterarr, .false., inwhichlocreg_reference, ldiis, &
                     lzd, lorbs, op, comon, mad, comgp, &
                     lphi, lhphi, lhphiold, lphiold)
                call allocateCommunicationsBuffersPotential(comgp, subname)

                !!     write(*,'(a,2i9)') 'sub 1: lorbs%npsidim_orbs, size(lphi)', lorbs%npsidim_orbs, size(lphi)
                !!do iorb=1,lorbs%norb
                !!    ilr=lorbs%inwhichlocreg(iorb)
                !!    write(*,'(a,2i6,3x,3es16.4)') '1: iorb, ilr, center', iorb, ilr, lzd%llr(ilr)%locregCenter
                !!end do
                !!call mpi_barrier(mpi_comm_world, ierr)
                !!stop
            end if


            !call deallocateCommunicationsBuffersPotential(comgp, subname)
            !call allocateCommunicationsBuffersPotential(comgp, subname)
            call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, comgp)

            if(secondLocreg) then
                ! Transform back to small locreg
                ind1=1
                ind2=1
                lphi=0.d0
                do iorb=1,lorbs%norbp
                    !ilr = lin%lig%orbsig%inWhichLocregp(iorb)
                    ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
                    ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
                    ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                    gdim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                    !write(*,'(a,7i12)') 'iproc, ind1, ind1+gdim-1, ind2, ind2+ldim-1, size(phi), size(lphilarge)', &
                    !      iproc, ind1, ind1+gdim-1, ind2, ind2+ldim-1, size(phi), size(lphilarge)
                    !call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdlarge%llr(ilr), lzd%glr, phi(ind1:ind1+gdim-1), lphilarge(ind2:ind2+ldim-1))
                    !!write(*,'(a,3i5,2(4x,3es16.4))') '1: iproc, iorb, ilr, lzd%llr(ilr)%locregCenter, locregCenter(:ilrlarge)', &
                    !!    iproc, iorb, ilr, lzd%llr(ilr)%locregCenter, locregCenter(:,ilrlarge)
                    call psi_to_locreg2(iproc, nproc, ldim, gdim, lzd%llr(ilr), lzdlarge%llr(ilrlarge), &
                         lphilarge(ind1:ind1+gdim-1), lphi(ind2:ind2+ldim-1))
                    !!do istat=ind2,ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1
                    !!    write(510+iproc,*) istat, lphi(istat)
                    !!end do
                    ind1=ind1+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                    ind2=ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                end do
            else
                ! Transform back to small locreg
                ind1=1
                ind2=1
                lphi=0.d0
                do iorb=1,lorbs%norbp
                    !ilr = lin%lig%orbsig%inWhichLocregp(iorb)
                    ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
                    ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
                    ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                    gdim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                    !write(*,'(a,7i12)') 'iproc, ind1, ind1+gdim-1, ind2, ind2+ldim-1, size(phi), size(lphilarge)', &
                    !      iproc, ind1, ind1+gdim-1, ind2, ind2+ldim-1, size(phi), size(lphilarge)
                    !call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdlarge%llr(ilr), lzd%glr, phi(ind1:ind1+gdim-1), lphilarge(ind2:ind2+ldim-1))
                    !write(*,'(a,3i5,2(4x,3es16.4))') '1: iproc, iorb, ilr, lzd%llr(ilr)%locregCenter, locregCenter(:ilrlarge)', &
                    !    iproc, iorb, ilr, lzd%llr(ilr)%locregCenter, locregCenter(:,ilrlarge)
                    call psi_to_locreg2(iproc, nproc, ldim, gdim, lzd%llr(ilr), lzdlarge%llr(ilrlarge), &
                         lphilarge(ind1:ind1+gdim-1), lphi(ind2:ind2+ldim-1))
                    !!do istat=ind2,ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1
                    !!    write(510+iproc,*) istat, lphi(istat)
                    !!end do
                    ind1=ind1+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                    ind2=ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                end do
            end if

            ! Update confdatarr...
            do iorb=1,orbslarge%norbp
               iiorb=orbslarge%isorb+iorb
               ilr=orbslarge%inWhichlocreg(iiorb)
               icenter=orbslarge%inwhichlocreg(iiorb)
               !confdatarr(iorb)%potorder=lin%confpotorder
               !confdatarr(iorb)%prefac=lin%potentialprefac(at%iatype(icenter))
               !confdatarr(iorb)%hh(1)=.5_gp*hx
               !confdatarr(iorb)%hh(2)=.5_gp*hy
               !confdatarr(iorb)%hh(3)=.5_gp*hz
               confdatarr(iorb)%rxyzConf(1:3)=locregCenter(1:3,icenter)
               call my_geocode_buffers(lzdlarge%Llr(ilr)%geocode,nl1,nl2,nl3)
               confdatarr(iorb)%ioffset(1)=lzdlarge%llr(ilr)%nsi1-nl1-1
               confdatarr(iorb)%ioffset(2)=lzdlarge%llr(ilr)%nsi2-nl2-1
               confdatarr(iorb)%ioffset(3)=lzdlarge%llr(ilr)%nsi3-nl3-1
            end do
            !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
            !!     lzdlarge%llr(ilrlarge)%locrad/factor, confdatarr, lphilarge)
            !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
            !!     6.d0, confdatarr, lphilarge)
 
            if(variable_locregs) then
                call destroy_new_locregs(lzdlarge, orbslarge, oplarge, comonlarge, madlarge, comgplarge, &
                     lphilarge, lhphilarge, lhphilargeold, lphilargeold)
                locrad_tmp=factor*locrad
                call create_new_locregs(iproc, nproc, lzd%nlr, hx, hy, hz, lorbs, lzd%glr, locregCenter, &
                     locrad_tmp, denspot%dpcom%nscatterarr, .false., inwhichlocreg_reference, ldiis, &
                     lzdlarge, orbslarge, oplarge, comonlarge, madlarge, comgplarge, &
                     lphilarge, lhphilarge, lhphilargeold, lphilargeold)
                locregCenterTemp=locregCenter
            end if

!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! do here another MLWF on the level of lzdlarge
!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


            if(secondLocreg) then
                call destroy_new_locregs(lzdlarge, orbslarge, oplarge, comonlarge, madlarge, comgplarge, &
                     lphilarge, lhphilarge, lhphilargeold, lphilargeold)
                locrad_tmp=factor2*locrad
                call create_new_locregs(iproc, nproc, lzd%nlr, hx, hy, hz, lorbs, lzd%glr, locregCenter, &
                     locrad_tmp, denspot%dpcom%nscatterarr, .false., inwhichlocreg_reference, ldiis2, &
                     lzdlarge, orbslarge, oplarge, comonlarge, madlarge, comgplarge, &
                     lphilarge, lhphilarge, lhphilargeold, lphilargeold)
            end if


        end if

      end if



      t2=mpi_wtime()
      time(1)=time(1)+t2-t1


      t1=mpi_wtime()
      !if(.not.ldiis%switchSD .and. .not.newgradient) then
      if(.not.ldiis%switchSD .and. newgradient) then
          !!!confdatarr%prefac=2.0d-2
          !!!call unitary_optimization(iproc, nproc, lzd, lorbs, at, op, &
          !!!                          comon, mad, rxyz, nItInnerLoop, kernel, &
          !!!                          newgradient, confdatarr, hx, lphi)
          !!!confdatarr%prefac=0.d0
      end if
      t2=mpi_wtime()
      time(5)=time(5)+t2-t1

      !!! Post the sends again to calculate the overlap matrix (will be needed for the orthoconstraint).
      !!call allocateSendBufferOrtho(comon, subname)
      !!call allocateRecvBufferOrtho(comon, subname)
      !!! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
      !!call extractOrbital3(iproc, nproc, lorbs, lorbs%npsidim, lorbs%inWhichLocreg, &
      !!     lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      !!! Post the send messages.
      !!call postCommsOverlapNew(iproc, nproc, lorbs, op, lzd, lphi, comon, timecommunp2p, timeextract)

  
      ! Calculate the unconstrained gradient by applying the Hamiltonian.
      !if(iproc==0) then
      !    write(*,'(1x,a)', advance='no') 'Hamiltonian application... '
      !end if
      t1=mpi_wtime()
      !withConfinement=.false.
      withConfinement=.true.
      allocate(lzd%doHamAppl(lorbs%norb), stat=istat)
      call memocc(istat, lzd%doHamAppl, 'lzd%doHamAppl', subname)
      lzd%doHamAppl=.true.



      if(newgradient) then
          ! Gather the potential that each process needs for the Hamiltonian application for all its orbitals.
          ! The messages for this point to point communication have been posted in the subroutine linearScaling.
          call gatherPotential(iproc, nproc, comgp)

          ! Build the required potential
          call local_potential_dimensions(lzd,lorbs,denspot%dpcom%ngatherarr(0,1))

          call full_local_potential(iproc,nproc,lorbs,Lzd,2,denspot%dpcom,denspot%rhov,denspot%pot_full,comgp)
      end if




      call FullHamiltonianApplication(iproc,nproc,at,lorbs,&
           hx,hy,hz,rxyz,&
           proj,lzd,nlpspd,confdatarr,denspot%dpcom%ngatherarr,denspot%pot_full,lphi,lhphi,&
           ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,SIC,GPU,&
           pkernel=denspot%pkernelseq)

      iall=-product(shape(lzd%doHamAppl))*kind(lzd%doHamAppl)
      deallocate(lzd%doHamAppl,stat=istat)
      call memocc(istat,iall,'lzd%doHamAppl',subname)
      ind2=1
      do iorb=1,lorbs%norbp
          ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
          !do istat=ind2,ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1
              !if(newgradient) write(530+iproc,*) istat, lhphi(istat)
          !end do
          ind2=ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      end do

   
      if(newgradient) then
          ! Deallocate potential
          iall=-product(shape(denspot%pot_full))*kind(denspot%pot_full)
          deallocate(denspot%pot_full, stat=istat)
          call memocc(istat, iall, 'denspot%pot_full', subname)

          !call allocateCommunicationsBuffersPotential(comgp, subname)
      end if


      !!!if(newgradient) then
      !!!    ! NEW: modify gradient
      !!!    call allocateSendBufferOrtho(comon, subname)
      !!!    call allocateRecvBufferOrtho(comon, subname)
      !!!    ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
      !!!    call extractOrbital3(iproc, nproc, lorbs, max(lorbs%npsidim_orbs,lorbs%npsidim_comp), &
      !!!         lorbs%inWhichLocreg, lzd, op, &
      !!!         lhphi, comon%nsendBuf, comon%sendBuf)
      !!!    call postCommsOverlapNew(iproc, nproc, lorbs, op, lzd, lhphi, comon, tt1, tt2)
      !!!    call collectnew(iproc, nproc, comon, mad, op, lorbs, lzd, comon%nsendbuf, &
      !!!         comon%sendbuf, comon%nrecvbuf, comon%recvbuf, tt3, tt4, tt5)
      !!!    call build_new_linear_combinations(iproc, nproc, lzd, lorbs, op, comon%nrecvbuf, &
      !!!         comon%recvbuf, kernel, .true., lhphi)
      !!!    call deallocateRecvBufferOrtho(comon, subname)
      !!!    call deallocateSendBufferOrtho(comon, subname)
      !!!    ind2=1
      !!!    do iorb=1,lorbs%norbp
      !!!        ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
      !!!        ind2=ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      !!!    end do
      !!!end if


      t2=mpi_wtime()

      ! Post the sends again to calculate the overlap matrix (will be needed for the orthoconstraint).
      call allocateSendBufferOrtho(comon, subname)
      call allocateRecvBufferOrtho(comon, subname)
      ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
      call extractOrbital3(iproc, nproc, lorbs, max(lorbs%npsidim_orbs,lorbs%npsidim_comp), lorbs%inWhichLocreg, &
           lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      ! Post the send messages.
      call postCommsOverlapNew(iproc, nproc, lorbs, op, lzd, lphi, comon, timecommunp2p, timeextract)


      time(2)=time(2)+t2-t1
      !if(iproc==0) then
      !    write(*,'(a)', advance='no') 'done. '
      !end if


  
      ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
      if(iproc==0) then
          write(*,'(a)', advance='no') ' Orthoconstraint... '
      end if
      !!!!!!!!!!!!!!!call orthoconstraintLocalized(iproc, nproc, lin, input, lphi, lhphi, trH)

      ! Gather the messages and calculate the overlap matrix.
      call collectnew(iproc, nproc, comon, mad, op, lorbs, lzd, comon%nsendbuf, &
           comon%sendbuf, comon%nrecvbuf, comon%recvbuf, timecommunp2p, timecommuncoll, timecompress)
      call calculateOverlapMatrix3(iproc, nproc, lorbs, op, lorbs%inWhichLocreg, comon%nsendBuf, &
           comon%sendBuf, comon%nrecvBuf, comon%recvBuf, mad, ovrlp)
      call deallocateRecvBufferOrtho(comon, subname)
      call deallocateSendBufferOrtho(comon, subname)

      t1=mpi_wtime()
      !!call flatten_at_edges(iproc, nproc, lin, at, input, hx, hy, hz, lorbs, lzd, rxyz, lhphi)
      !!if(.not.newgradient .or. it<=10) call orthoconstraintNonorthogonal(iproc, nproc, lzd, lorbs, op, comon, mad, ovrlp, &
      !!     methTransformOverlap, blocksize_pdgemm, lphi, lhphi, lagmat)

      !!if(iproc==0) write(*,*) 'WARNING: COMMENTED!'
      if(.not.newgradient) then
          call orthoconstraintNonorthogonal(iproc, nproc, lzd, lorbs, op, comon, mad, ovrlp, &
               methTransformOverlap, blocksize_pdgemm, lphi, lhphi, lagmat)
      else
          lphilarge=0.d0
          lhphilarge=0.d0
          istl=1
          istg=1
          do iorb=1,orbslarge%norbp
              ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
              ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
              ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
              gdim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
              !write(*,'(a,7i12)') 'iproc, ind1, ind1+gdim-1, ind2, ind2+ldim-1, size(phi), size(lphilarge)', &
              !      iproc, ind1, ind1+gdim-1, ind2, ind2+ldim-1, size(phi), size(lphilarge)
              !call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdlarge%llr(ilr), lzd%glr, phi(ind1:ind1+gdim-1), lphilarge(ind2:ind2+ldim-1))
              !call psi_to_locreg2(iproc, nproc, ldim, gdim, lzd%llr(ilr), lzdlarge%llr(ilr), lphi(ind1), lphilarge(ind2))
              nspin=1
              call Lpsi_to_global2(iproc, nproc, ldim, gdim, lorbs%norb, lorbs%nspinor, nspin, lzdlarge%llr(ilrlarge), &
                   lzd%llr(ilr), lphi(istl), lphilarge(istg))
              call Lpsi_to_global2(iproc, nproc, ldim, gdim, lorbs%norb, lorbs%nspinor, nspin, lzdlarge%llr(ilrlarge), &
                   lzd%llr(ilr), lhphi(istl), lhphilarge(istg))
              istl=istl+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
              istg=istg+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
          end do


          if(newgradient) then
              ! NEW: modify gradient
              call allocateSendBufferOrtho(comonlarge, subname)
              call allocateRecvBufferOrtho(comonlarge, subname)
              ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
              call extractOrbital3(iproc, nproc, orbslarge, max(orbslarge%npsidim_orbs,orbslarge%npsidim_comp), &
                   orbslarge%inWhichLocreg, lzdlarge, oplarge, &
                   lhphilarge, comonlarge%nsendBuf, comonlarge%sendBuf)
              call postCommsOverlapNew(iproc, nproc, orbslarge, oplarge, lzdlarge, lhphilarge, comonlarge, tt1, tt2)
              call collectnew(iproc, nproc, comonlarge, madlarge, oplarge, orbslarge, lzdlarge, comonlarge%nsendbuf, &
                   comonlarge%sendbuf, comonlarge%nrecvbuf, comonlarge%recvbuf, tt3, tt4, tt5)
              call build_new_linear_combinations(iproc, nproc, lzdlarge, orbslarge, oplarge, comonlarge%nrecvbuf, &
                   comonlarge%recvbuf, kernel, .true., lhphilarge)
              call deallocateRecvBufferOrtho(comonlarge, subname)
              call deallocateSendBufferOrtho(comonlarge, subname)
              ind2=1
              do iorb=1,orbslarge%norbp
                  ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
                  ind2=ind2+lzd%llr(ilrlarge)%wfd%nvctr_c+7*lzd%llr(ilrlarge)%wfd%nvctr_f
              end do
          end if


          call orthoconstraintNonorthogonal(iproc, nproc, lzdlarge, orbslarge, &
               oplarge, comonlarge, madlarge, ovrlp, &
               orthpar%methTransformOverlap, orthpar%blocksize_pdgemm, lphilarge, lhphilarge, lagmat)
      end if



      ! Calculate modified trace
      !!tt=0.d0
      !!do iorb=1,orbs%norb
      !!    do jorb=1,lorbs%norb
      !!        do korb=1,lorbs%norb
      !!            tt = tt + coeff(jorb,iorb)*coeff(korb,iorb)*W(korb,jorb)
      !!        end do
      !!    end do
      !!end do
      if(newgradient) then
          trH=0.d0
          do jorb=1,lorbs%norb
              do korb=1,lorbs%norb
                  trH = trH + kernel(korb,jorb)*lagmat(korb,jorb)
              end do
          end do
      else
          trH=0.d0
          do jorb=1,lorbs%norb
              trH = trH + lagmat(jorb,jorb)
          end do
      end if


      ! Flatten at the edges -  EXPERIMENTAL
      !!call flatten_at_edges(iproc, nproc, lin, at, input, hx, hy, hz, lorbs, lzd, rxyz, lhphi)

      ! Cycle if the trace increased (steepest descent only)
      !if(iproc==0) write(*,*) 'ldiis%switchSD, ldiis%isx', ldiis%switchSD, ldiis%isx
      if(.not. ldiis%switchSD .and. ldiis%isx==0) then
         !if(iproc==0) write(*,*) 'trH, trHold', trH, trHold
           !if(trH>trHold) then
           if(trH > trHold + 1.d-8*abs(trHold)) then
               consecutive_rejections=consecutive_rejections+1
               !!if(consecutive_rejections==300000) then
               !!    if(fnrmMax<convCrit .or. it>=nit) then
               !!        if(it>=nit) then
               !!            if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
               !!                ' iterations! Exiting loop due to consective failures of SD.'
               !!            if(iproc==0) write(*,'(1x,a,2es15.7,f12.7)') 'CHECK THIS Final values for fnrm, fnrmMax, trace: ', &
               !!                         fnrm, fnrmMax, trHold
               !!            infoBasisFunctions=-1
               !!        else
               !!            if(iproc==0) then
               !!                write(*,'(1x,a,i0,a,2es15.7,f12.7)') 'converged in ', it, ' iterations.'
               !!                write (*,'(1x,a,2es15.7,f12.7)') 'CHECK THIS Final values for fnrm, fnrmMax, trace: ', &
               !!                fnrm, fnrmMax, trHold
               !!            end if
               !!            infoBasisFunctions=it
               !!        end if
               !!        if(iproc==0) write(*,'(1x,a)') '============================= Basis functions created. &
               !!                     &============================='
               !!        !!if(lin%plotBasisFunctions) then
               !!        !!    call plotOrbitals(iproc, lorbs, Glr, phi, at%nat, rxyz, lin%onWhichAtom, .5d0*input%hx, &
               !!        !!        .5d0*input%hy, .5d0*input%hz, 1)
               !!        !!end if
               !!        exit iterLoop
               !!    end if
               !!end if
               if(iproc==0) write(*,'(a,2es20.12)') 'trH, trHold', trH, trHold
               if(consecutive_rejections<=3) then
                   alpha=alpha*.6d0
                   if(.not. newgradient) then
                       call dcopy(size(lphi), lphiold, 1, lphi, 1)
                   else
                       ! this makes no sense
                       !call dcopy(size(lphilarge), lphilargeold, 1, lphilarge, 1)
                   end if
                   !if(iproc==0) write(*,'(x,a)') 'trace increased; reject orbitals and cycle'
                   !cycle iterLoop
                   !if(iproc==0) write(*,'(x,a)') 'trace increased; reject orbitals and exit'
                   !exit iterLoop
               else
                   consecutive_rejections=0
               end if
           else
                   consecutive_rejections=0
           end if
      end if
      !consecutive_rejections=0


      t2=mpi_wtime()
      time(3)=time(3)+t2-t1
      !if(iproc==0) then
      !    write(*,'(a)', advance='no') 'done. '
      !end if
      !!if(iproc==0) write(*,*) 'trH',trH


  
      ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
      ! of the previous iteration (fnrmOvrlpArr).
      istart=1
      do iorb=1,lorbs%norbp
          !ilr=lorbs%inWhichLocregp(iorb)
          if(.not.newgradient) then
              iiorb=lorbs%isorb+iorb
              ilr=lorbs%inWhichLocreg(iiorb)
              ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
              if(it>1) fnrmOvrlpArr(iorb,1)=ddot(ncount, lhphi(istart), 1, lhphiold(istart), 1)
              fnrmArr(iorb,1)=ddot(ncount, lhphi(istart), 1, lhphi(istart), 1)
          else
              iiorb=orbslarge%isorb+iorb
              ilr=orbslarge%inWhichLocreg(iiorb)
              ncount=lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
              !if(it>1) fnrmOvrlpArr(iorb,1)=ddot(ncount, lhphilarge(istart), 1, lhphilargeold(istart), 1)
              ! The following line is nonsense, but it is not possible otherwise (since the shape of the locregs may change)
              if(it>1) fnrmOvrlpArr(iorb,1)=ddot(ncount, lhphilarge(istart), 1, lhphilarge(istart), 1)
              fnrmArr(iorb,1)=ddot(ncount, lhphilarge(istart), 1, lhphilarge(istart), 1)
          end if
          istart=istart+ncount
      end do

      ! Keep the gradient for the next iteration.
      if(it>1) then
          call dcopy(lorbs%norbp, fnrmArr(1,1), 1, fnrmOldArr(1), 1)
      end if
  
      ! Determine the gradient norm and its maximal component. In addition, adapt the
      ! step size for the steepest descent minimization (depending on the angle 
      ! between the current gradient and the one from the previous iteration).
      ! This is of course only necessary if we are using steepest descent and not DIIS.
      do iorb=1,lorbs%norbp
          fnrm=fnrm+fnrmArr(iorb,1)
          if(fnrmArr(iorb,1)>fnrmMax) fnrmMax=fnrmArr(iorb,1)
          if(it>1 .and. ldiis%isx==0 .and. .not.ldiis%switchSD) then
          ! Adapt step size for the steepest descent minimization.
              tt=fnrmOvrlpArr(iorb,1)/sqrt(fnrmArr(iorb,1)*fnrmOldArr(iorb))
              !if(tt>.7d0) then
              if(tt>.9d0 .and. trH<trHold) then
                  alpha(iorb)=alpha(iorb)*1.1d0
              else
                  alpha(iorb)=alpha(iorb)*.6d0
              end if
          end if
      end do
      call mpiallred(fnrm, 1, mpi_sum, mpi_comm_world, ierr)
      call mpiallred(fnrmMax, 1, mpi_max, mpi_comm_world, ierr)
      !fnrm=sqrt(fnrm)
      fnrm=sqrt(fnrm/dble(lorbs%norb))
      fnrmMax=sqrt(fnrmMax)
      ! Copy the gradient (will be used in the next iteration to adapt the step size).
      !call dcopy(max(lorbs%npsidim_orbs,lorbs%npsidim_comp), lhphi, 1, lhphiold, 1)
      call dcopy(lorbs%npsidim_orbs, lhphi, 1, lhphiold, 1)
      if(newgradient) call dcopy(max(orbslarge%npsidim_orbs,orbslarge%npsidim_comp), lhphilarge, 1, lhphilargeold, 1)
      trHold=trH
  
      ! Precondition the gradient.
      if(iproc==0) then
          !write(*,'(a)', advance='no') 'Preconditioning... '
          write(*,'(a)') 'Preconditioning.'
      end if
      gnrm=1.d3 ; gnrm_zero=1.d3
      t1=mpi_wtime()

      !! THIS WAS THE ORIGINAL
      !!ind2=1
      !!do iorb=1,lorbs%norbp
      !!    iiorb=lorbs%isorb+iorb
      !!    ilr = lorbs%inWhichLocreg(iiorb)
      !!    call choosePreconditioner2(iproc, nproc, lorbs, lzd%Llr(ilr), hx, hy, hz, &
      !!        nItPrecond, lhphi(ind2), at%nat, rxyz, at, confdatarr(iorb)%potorder, confdatarr(iorb)%prefac, it, iorb, eval_zero)
      !!    ind2=ind2+lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
      !!end do
      ind2=1
      do iorb=1,lorbs%norbp
          if(.not.newgradient) then
             iiorb=lorbs%isorb+iorb
             ilr = lorbs%inWhichLocreg(iiorb)
             call choosePreconditioner2(iproc, nproc, lorbs, lzd%llr(ilr), hx, hy, hz, &
                 nItPrecond, lhphi(ind2), at%nat, rxyz, at, confdatarr(iorb)%potorder, confdatarr(iorb)%prefac, it, iorb, eval_zero)
             ind2=ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
         else
             iiorb=orbslarge%isorb+iorb
             ilr = orbslarge%inWhichLocreg(iiorb)
             call choosePreconditioner2(iproc, nproc, orbslarge, lzdlarge%llr(ilr), hx, hy, hz, &
                 nItPrecond, lhphilarge(ind2), lzdlarge%nlr, rxyz, at, confdatarr(iorb)%potorder, &
                 confdatarr(iorb)%prefac, it, iorb, eval_zero)
             !!write(*,'(a,i7,es16.6)') 'iiorb, ddot1', iiorb, ddot(lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f, lphilarge(ind2), 1, lphilarge(ind2), 1)
             !!write(*,'(a,i7,es16.6)') 'iiorb, ddot2', iiorb, ddot(lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f, lhphilarge(ind2), 1, lhphilarge(ind2), 1)
             ind2=ind2+lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
         end if
      end do




      !!!end if
      !call preconditionall(iproc, nproc, lorbs, lzd%glr, input%hx, input%hy, input%hz, lin%nItPrecond, lhphi, tt, tt2)
      t2=mpi_wtime()
      time(4)=time(4)+t2-t1
      !if(iproc==0) then
      !    write(*,'(a)') 'done. '
      !end if


      !!!! plot the orbitals -- EXPERIMENTAL ##################################################
      !!!allocate(lvphiovrlp(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f))
      !!!ist=1
      !!!write(comment,'(i3.3)') it
      !!!do iorb=1,lorbs%norbp
      !!!    iiorb=iorb+lorbs%isorb
      !!!    ilr=lorbs%inwhichlocreg(iiorb)
      !!!    write(orbname,'(i3.3)') iiorb
      !!!    write(*,'(a,i0)') 'plotting orbital ',iiorb
      !!!    lvphiovrlp=0.d0
      !!!    call Lpsi_to_global2(iproc, nproc, lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f, &
      !!!         lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f, lorbs%norb, lorbs%nspinor, input%nspin, &
      !!!         lzd%Glr, lzd%Llr(ilr), lphi(ist), lvphiovrlp(1))
      !!!    call plot_wf(orbname//'_'//comment, 2, at, 1.d0, lzd%glr, input%hx, input%hy, input%hz, rxyz, lvphiovrlp(1))
      !!!    ist=ist+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      !!!end do
      !!!deallocate(lvphiovrlp)
      !!!! ####################################################################################


      ! Determine the mean step size for steepest descent iterations.
      tt=sum(alpha)
      meanAlpha=tt/dble(lorbs%norb)
  
      ! Write some informations to the screen.
      if(iproc==0) write(*,'(1x,a,i6,2es15.7,f17.10)') 'iter, fnrm, fnrmMax, trace', it, fnrm, fnrmMax, trH
      if(fnrmMax<convCrit .or. it>=nit) then
          if(it>=nit) then
              if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0) write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              infoBasisFunctions=-1
          else
              if(iproc==0) then
                  write(*,'(1x,a,i0,a,2es15.7,f12.7)') 'converged in ', it, ' iterations.'
                  write (*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              end if
              infoBasisFunctions=it
          end if
          if(iproc==0) write(*,'(1x,a)') '============================= Basis functions created. ============================='
          !!if(lin%plotBasisFunctions) then
          !!    call plotOrbitals(iproc, lorbs, Glr, phi, at%nat, rxyz, lin%onWhichAtom, .5d0*input%hx, &
          !!        .5d0*input%hy, .5d0*input%hz, 1)
          !!end if


          exit iterLoop
      end if
  
  
      ! Determine whether the basis functions shall be further optimized using DIIS or steepest descent.
      call DIISorSD()
      if(iproc==0) then
          if(ldiis%isx>0) then
              write(*,'(1x,3(a,i0))') 'DIIS informations: history length=',ldiis%isx, ', consecutive failures=', &
                  icountDIISFailureCons, ', total failures=', icountDIISFailureTot
          else
              !if(allowDIIS) then
              !    message='y'
              !else
              !    message='n'
              !end if
              write(*,'(1x,a,es9.3,a,i0,a)') 'steepest descent informations: mean alpha=', meanAlpha, &
              ', consecutive successes=', icountSDSatur, ', DIIS=y'
          end if
      end if


      ! Improve the orbitals, depending on the choice made above.
      if(.not.ldiis%switchSD) then
          call improveOrbitals()
      else
          if(iproc==0) write(*,'(x,a)') 'no improvement of the orbitals, recalculate gradient'
      end if

      
      if(newgradient) then
           !!confdatarr%prefac=1.0d-4
           !!call unitary_optimization(iproc, nproc, lzdlarge, orbslarge, at, oplarge, &
           !!                          comonlarge, madlarge, rxyz, nItInnerLoop, kernel, &
           !!                          newgradient, confdatarr, hx, lphilarge)
           !!confdatarr%prefac=0.0d0

           ! Update confdatarr...
           do iorb=1,orbslarge%norbp
              iiorb=orbslarge%isorb+iorb
              ilr=orbslarge%inWhichlocreg(iiorb)
              icenter=orbslarge%inwhichlocreg(iiorb)
              !confdatarr(iorb)%potorder=lin%confpotorder
              !confdatarr(iorb)%prefac=lin%potentialprefac(at%iatype(icenter))
              !confdatarr(iorb)%hh(1)=.5_gp*hx
              !confdatarr(iorb)%hh(2)=.5_gp*hy
              !confdatarr(iorb)%hh(3)=.5_gp*hz
              confdatarr(iorb)%rxyzConf(1:3)=locregCenterTemp(1:3,icenter)
              call my_geocode_buffers(lzdlarge%Llr(ilr)%geocode,nl1,nl2,nl3)
              confdatarr(iorb)%ioffset(1)=lzdlarge%llr(ilr)%nsi1-nl1-1
              confdatarr(iorb)%ioffset(2)=lzdlarge%llr(ilr)%nsi2-nl2-1
              confdatarr(iorb)%ioffset(3)=lzdlarge%llr(ilr)%nsi3-nl3-1
           end do

           !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
           !!     lzdlarge%llr(ilrlarge)%locrad/factor, confdatarr, lphilarge)
           !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
           !!     6.d0, confdatarr, lphilarge)

           !! EXPERIMENTAL: normalize lphilarge
           if(variable_locregs) then
               ist=1
               do iorb=1,orbslarge%norbp
                   iiorb=orbslarge%isorb+iorb
                   ilrlarge=orbslarge%inwhichlocreg(iiorb)
                   ncnt=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                   tt=dnrm2(ncnt, lphilarge(ist), 1)
                   call dscal(ncnt, 1/tt, lphilarge(ist), 1)
                   ist=ist+ncnt
               end do
           end if


           if(secondLocreg) then
                ! Go to even larger region
                lphilarge=0.d0
                istl=1
                istg=1
                do iorb=1,orbslarge%norbp
                    ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
                    ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
                    ldim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                    gdim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                    nspin=1
                    call Lpsi_to_global2(iproc, nproc, ldim, gdim, lorbs%norb, lorbs%nspinor, nspin, lzdlarge%llr(ilrlarge), &
                         lzdlarge%llr(ilrlarge), lphilarge(istl), lphilarge(istg))
                    istl=istl+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                    istg=istg+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                end do
                if(istl/=orbslarge%npsidim_orbs+1) stop 'istl/=orbslarge%npsidim_orbs+1'
                if(istg/=orbslarge%npsidim_orbs+1) stop 'istg/=orbslarge%npsidim_orbs+1'
            end if

           if(.not.secondLocreg) then
                ! Update confdatarr...
                do iorb=1,orbslarge%norbp
                   iiorb=orbslarge%isorb+iorb
                   ilr=orbslarge%inWhichlocreg(iiorb)
                   icenter=orbslarge%inwhichlocreg(iiorb)
                   !confdatarr(iorb)%potorder=lin%confpotorder
                   !confdatarr(iorb)%prefac=lin%potentialprefac(at%iatype(icenter))
                   !confdatarr(iorb)%hh(1)=.5_gp*hx
                   !confdatarr(iorb)%hh(2)=.5_gp*hy
                   !confdatarr(iorb)%hh(3)=.5_gp*hz
                   confdatarr(iorb)%rxyzConf(1:3)=locregCenterTemp(1:3,icenter)
                   call my_geocode_buffers(lzdlarge%Llr(ilr)%geocode,nl1,nl2,nl3)
                   confdatarr(iorb)%ioffset(1)=lzdlarge%llr(ilr)%nsi1-nl1-1
                   confdatarr(iorb)%ioffset(2)=lzdlarge%llr(ilr)%nsi2-nl2-1
                   confdatarr(iorb)%ioffset(3)=lzdlarge%llr(ilr)%nsi3-nl3-1
                end do
               call MLWFnew(iproc, nproc, lzdlarge, orbslarge, at, oplarge, &
                    comonlarge, madlarge, rxyz, nItInnerLoop, kernel, &
                    newgradient, confdatarr, hx, locregCenterTemp, 3.d0, lphilarge, Umat, locregCenter)
            end if

            if(secondLocreg) then
                ! Update confdatarr...
                do iorb=1,orbslarge2%norbp
                   iiorb=orbslarge2%isorb+iorb
                   ilr=orbslarge2%inWhichlocreg(iiorb)
                   icenter=orbslarge2%inwhichlocreg(iiorb)
                   !confdatarr(iorb)%potorder=lin%confpotorder
                   !confdatarr(iorb)%prefac=lin%potentialprefac(at%iatype(icenter))
                   !confdatarr(iorb)%hh(1)=.5_gp*hx
                   !confdatarr(iorb)%hh(2)=.5_gp*hy
                   !confdatarr(iorb)%hh(3)=.5_gp*hz
                   confdatarr(iorb)%rxyzConf(1:3)=locregCenterTemp(1:3,icenter)
                   call my_geocode_buffers(lzdlarge2%Llr(ilr)%geocode,nl1,nl2,nl3)
                   confdatarr(iorb)%ioffset(1)=lzdlarge2%llr(ilr)%nsi1-nl1-1
                   confdatarr(iorb)%ioffset(2)=lzdlarge2%llr(ilr)%nsi2-nl2-1
                   confdatarr(iorb)%ioffset(3)=lzdlarge2%llr(ilr)%nsi3-nl3-1
                end do
               !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
               !!     lzdlarge%llr(ilrlarge)%locrad/factor, confdatarr, lphilarge)
                call check_cutoff(iproc, nproc, orbslarge2, lzdlarge2, hx, hy, hz, &
                     lzdlarge2%llr(1)%locrad/factor2, confdatarr, lphilarge2)
               !!call check_cutoff(iproc, nproc, orbslarge2, lzdlarge2, hx, hy, hz, &
               !!     6.d0, confdatarr, lphilarge2)
               call MLWFnew(iproc, nproc, lzdlarge2, orbslarge2, at, oplarge2, &
                    comonlarge2, madlarge2, rxyz, nItInnerLoop, kernel, &
                    newgradient, confdatarr, hx, locregCenterTemp, 3.d0, lphilarge2, Umat, locregCenter)
               !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
               !!     lzdlarge%llr(ilrlarge)%locrad/factor, confdatarr, lphilarge)
                call check_cutoff(iproc, nproc, orbslarge2, lzdlarge2, hx, hy, hz, &
                     lzdlarge2%llr(1)%locrad/factor2, confdatarr, lphilarge2)
               !!call check_cutoff(iproc, nproc, orbslarge2, lzdlarge2, hx, hy, hz, &
               !!     6.d0, confdatarr, lphilarge2)
                     !!write(*,*) "ATTENTION HERE!!!"
                     !!locregCenter=locregCenterTemp
           end if
           do ilr=1,lzd%nlr
               tt = (locregCenter(1,ilr)-locregCenterTemp(1,ilr))**2 &
                   +(locregCenter(2,ilr)-locregCenterTemp(2,ilr))**2 &
                   +(locregCenter(3,ilr)-locregCenterTemp(3,ilr))**2 
               tt=sqrt(tt)
               !!if(iproc==0) write(*,'(a,i5,3es13.4)') 'ilr, new centers:', &
               !!     ilr, locregCenter(1,ilr), locregCenter(2,ilr), locregCenter(3,ilr)
               !!if(iproc==0) write(*,'(a,i5,es16.4)') '2: ilr, tt', ilr, tt
               if( floor(locregCenter(1,ilr)/hx) < 0 .or. ceiling(locregCenter(1,ilr)/hx) > lzd%glr%d%n1 ) then
                   if(iproc==0) then
                       write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
                           ' is outside of box in x direction! Box limits=',0,lzd%glr%d%n1,&
                           ', center=',floor(locregCenter(1,ilr)/hx),ceiling(locregCenter(1,ilr)/hx)
                   end if
                   call mpi_barrier(mpi_comm_world, ierr)
                   stop
               end if
               if( floor(locregCenter(2,ilr)/hx) < 0 .or. ceiling(locregCenter(2,ilr)/hy) > lzd%glr%d%n2 ) then
                   if(iproc==0) then
                       write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
                           'is outside of box in y direction! Box limits=',0,lzd%glr%d%n2,&
                           ', center=',floor(locregCenter(2,ilr)/hy),ceiling(locregCenter(2,ilr)/hy)
                   end if
                   call mpi_barrier(mpi_comm_world, ierr)
                   stop
               end if
               if( floor(locregCenter(3,ilr)/hz) < 0 .or. ceiling(locregCenter(3,ilr)/hz) > lzd%glr%d%n3 ) then
                   if(iproc==0) then
                       write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
                           'is outside of box in x direction! Box limits=',0,lzd%glr%d%n3,&
                           ', center=',floor(locregCenter(3,ilr)/hz),ceiling(locregCenter(3,ilr)/hz)
                   end if
                   call mpi_barrier(mpi_comm_world, ierr)
                   stop
               end if
           end do
           !write(*,*) 'WARNING CHECK INDICES OF UMAT!!'
           if(nItInnerLoop>0) then
               kernelold=kernel
               do iorb=1,lorbs%norb
                   do jorb=1,lorbs%norb
                       tt=0.d0
                       do korb=1,lorbs%norb
                           do lorb=1,lorbs%norb
                               !tt=tt+kernelold(korb,lorb)*Umat(jorb,korb)*Umat(iorb,lorb)
                               !tt=tt+kernelold(korb,lorb)*Umat(iorb,korb)*Umat(jorb,lorb)
                               tt=tt+kernelold(korb,lorb)*Umat(korb,iorb)*Umat(lorb,jorb)
                               !tt=tt+kernelold(korb,lorb)*Umat(iorb,korb)*Umat(jorb,lorb)
                           end do
                       end do
                       kernel(jorb,iorb)=tt
                   end do
               end do
               !!kernelold=kernel
               !!do iorb=1,lorbs%norb
               !!    do jorb=1,lorbs%norb
               !!        tt=0.d0
               !!        !!tt2=0.d0
               !!        do korb=1,lorbs%norb
               !!            do lorb=1,lorbs%norb
               !!                !tt=tt+kernelold(korb,lorb)*Umat(jorb,korb)*Umat(iorb,lorb)
               !!                tt=tt+kernelold(korb,jorb)*Umat(korb,lorb)*Umat(iorb,lorb)
               !!            end do
               !!        end do
               !!        kernel(jorb,iorb)=tt
               !!    end do
               !!end do
          end if


        if(variable_locregs) then
            call destroy_new_locregs(lzd, lorbs, op, comon, mad, comgp, &
                 lphi, lhphi, lhphiold, lphiold)
            call create_new_locregs(iproc, nproc, lzdlarge%nlr, hx, hy, hz, orbslarge, lzdlarge%glr, locregCenter, &
                 locrad, denspot%dpcom%nscatterarr, .false., inwhichlocreg_reference, ldiis, &
                 lzd, lorbs, op, comon, mad, comgp, &
                 lphi, lhphi, lhphiold, lphiold)
        end if
          !!       write(*,'(a,2i9)') 'sub 2: lorbs%npsidim_orbs, size(lphi)', lorbs%npsidim_orbs, size(lphi)
          !!if(iproc==0) then
          !!    do iorb=1,lorbs%norb
          !!        ilr=lorbs%inwhichlocreg(iorb)
          !!        write(*,'(a,2i8)') 'lorbs%inwhichlocreg(iorb), orbslarge%inwhichlocreg(iorb)', &
          !!            lorbs%inwhichlocreg(iorb), orbslarge%inwhichlocreg(iorb)
          !!        write(*,'(a,2i6,3x,3es16.4)') '2: iorb, ilr, center', iorb, ilr, lzd%llr(ilr)%locregCenter
          !!    end do
          !!end if

          if(secondLocreg) then
              ! Transform back to localization region
              ! Transform back to small locreg
              ind1=1
              ind2=1
              lphi=0.d0
              do iorb=1,lorbs%norbp
                  !ilr = lin%lig%orbsig%inWhichLocregp(iorb)
                  ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
                  ilrlarge2 = orbslarge2%inWhichLocreg(orbslarge2%isorb+iorb)
                  !write(*,'(a,4i8)') 'iproc, iorb, ilr, ilrlarge2', iproc, iorb, ilr, ilrlarge2
                  ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                  gdim=lzdlarge2%llr(ilrlarge2)%wfd%nvctr_c+7*lzdlarge2%llr(ilrlarge2)%wfd%nvctr_f
                  !write(*,'(a,7i12)') 'iproc, ind1, ind1+gdim-1, ind2, ind2+ldim-1, size(phi), size(lphilarge2)', &
                  !      iproc, ind1, ind1+gdim-1, ind2, ind2+ldim-1, size(phi), size(lphilarge2)
                  !call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdlarge2%llr(ilr), lzd%glr, phi(ind1:ind1+gdim-1), lphilarge2(ind2:ind2+ldim-1))
                   !!write(*,'(a,3i5,2(4x,3es16.4))') '2: iproc, iorb, ilr, lzd%llr(ilr)%locregCenter, locregCenter(:,ilrlarge2)', &
                   !!    iproc, iorb, ilr, lzd%llr(ilr)%locregCenter, locregCenter(:,ilrlarge2)
                  call psi_to_locreg2(iproc, nproc, ldim, gdim, lzd%llr(ilr), lzdlarge2%llr(ilrlarge2), &
                       lphilarge2(ind1:ind1+gdim-1), lphi(ind2:ind2+ldim-1))
                  ind1=ind1+lzdlarge2%llr(ilrlarge2)%wfd%nvctr_c+7*lzdlarge2%llr(ilrlarge2)%wfd%nvctr_f
                  ind2=ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
              end do
          else
              ind1=1
              ind2=1
              lphi=0.d0
              do iorb=1,lorbs%norbp
                  !ilr = lin%lig%orbsig%inWhichLocregp(iorb)
                  ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
                  ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
                  !write(*,'(a,4i8)') 'iproc, iorb, ilr, ilrlarge', iproc, iorb, ilr, ilrlarge
                  ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                  gdim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                  !write(*,'(a,7i12)') 'iproc, ind1, ind1+gdim-1, ind2, ind2+ldim-1, size(phi), size(lphilarge)', &
                  !      iproc, ind1, ind1+gdim-1, ind2, ind2+ldim-1, size(phi), size(lphilarge)
                  !call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdlarge%llr(ilr), lzd%glr, phi(ind1:ind1+gdim-1), lphilarge(ind2:ind2+ldim-1))
                   !!write(*,'(a,3i5,2(4x,3es16.4))') '2: iproc, iorb, ilr, lzd%llr(ilr)%locregCenter, locregCenter(:,ilrlarge)', &
                   !!    iproc, iorb, ilr, lzd%llr(ilr)%locregCenter, locregCenter(:,ilrlarge)
                  call psi_to_locreg2(iproc, nproc, ldim, gdim, lzd%llr(ilr), lzdlarge%llr(ilrlarge), &
                       lphilarge(ind1:ind1+gdim-1), lphi(ind2:ind2+ldim-1))
                  ind1=ind1+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
                  ind2=ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
              end do
          end if

           ! Update confdatarr...
           do iorb=1,orbslarge%norbp
              iiorb=orbslarge%isorb+iorb
              ilr=orbslarge%inWhichlocreg(iiorb)
              icenter=orbslarge%inwhichlocreg(iiorb)
              !confdatarr(iorb)%potorder=lin%confpotorder
              !confdatarr(iorb)%prefac=lin%potentialprefac(at%iatype(icenter))
              !confdatarr(iorb)%hh(1)=.5_gp*hx
              !confdatarr(iorb)%hh(2)=.5_gp*hy
              !confdatarr(iorb)%hh(3)=.5_gp*hz
              confdatarr(iorb)%rxyzConf(1:3)=locregCenter(1:3,icenter)
              call my_geocode_buffers(lzdlarge%Llr(ilr)%geocode,nl1,nl2,nl3)
              confdatarr(iorb)%ioffset(1)=lzdlarge%llr(ilr)%nsi1-nl1-1
              confdatarr(iorb)%ioffset(2)=lzdlarge%llr(ilr)%nsi2-nl2-1
              confdatarr(iorb)%ioffset(3)=lzdlarge%llr(ilr)%nsi3-nl3-1
           end do
           !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
           !!     lzdlarge%llr(ilrlarge)%locrad/factor, confdatarr, lphilarge)
           !!call check_cutoff(iproc, nproc, orbslarge, lzdlarge, hx, hy, hz, &
           !!     6.d0, confdatarr, lphilarge)

           if(variable_locregs) then
               call destroy_new_locregs(lzdlarge, orbslarge, oplarge, comonlarge, madlarge, comgplarge, &
                    lphilarge, lhphilarge, lhphilargeold, lphilargeold)
               locrad_tmp=factor*locrad
               call create_new_locregs(iproc, nproc, lzd%nlr, hx, hy, hz, lorbs, lzd%glr, locregCenter, &
                    locrad_tmp, denspot%dpcom%nscatterarr, .false., inwhichlocreg_reference, ldiis, &
                    lzdlarge, orbslarge, oplarge, comonlarge, madlarge, comgplarge, &
                    lphilarge, lhphilarge, lhphilargeold, lphilargeold)
               locregCenterTemp=locregCenter
           end if

!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! do here another MLWF on the level of lzdlarge
!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

           if(secondLocreg) then
               call destroy_new_locregs(lzdlarge2, orbslarge2, oplarge2, comonlarge2, madlarge2, comgplarge2, &
                    lphilarge2, lhphilarge2, lhphilarge2old, lphilarge2old)
               locrad_tmp=factor2*locrad
               call create_new_locregs(iproc, nproc, lzd%nlr, hx, hy, hz, lorbs, lzd%glr, locregCenter, &
                    locrad_tmp, denspot%dpcom%nscatterarr, .false., inwhichlocreg_reference, ldiis2, &
                    lzdlarge2, orbslarge2, oplarge2, comonlarge2, madlarge2, comgplarge2, &
                    lphilarge2, lhphilarge2, lhphilarge2old, lphilarge2old)
!          write(c3,'(i0)') it
!              write(c2,'(i0)') iorb
!              filename='TMB_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)
!              open(unitwf,file=trim(filename),status='unknown')
!              call writeonewave(unitwf,.true.,1,lzd%Glr%d%n1,lzd%Glr%d%n2,lzd%Glr%d%n3,hx,hy,hz,at%nat,rxyz,  &
!                     lzd%Glr%wfd%nseg_c,lzd%Glr%wfd%nvctr_c,lzd%Glr%wfd%keygloc(1,1),lzd%Glr%wfd%keyvloc(1),  &
!                     lzd%Glr%wfd%nseg_f,lzd%Glr%wfd%nvctr_f,lzd%Glr%wfd%keygloc(1,1+lzd%Glr%wfd%nseg_c),&
!                     lzd%Glr%wfd%keyvloc(1+lzd%Glr%wfd%nseg_c), &
!                     lvphiovrlp(1),lvphiovrlp(1+lzd%glr%wfd%nvctr_c),lorbs%eval)
!              close(unitwf)
           end if
      end if





      !!! plot the orbitals -- EXPERIMENTAL ##################################################
      !! Modify this (rxyz is not correct...)
      !!if(newgradient) then
      !!    allocate(lvphiovrlp(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f))
      !!    ist=1
      !!    write(comment,'(i3.3)') it
      !!    do iorb=1,lorbs%norbp
      !!        iiorb=iorb+lorbs%isorb
      !!        ilr=lorbs%inwhichlocreg(iiorb)
      !!        lvphiovrlp=0.d0
      !!        call Lpsi_to_global2(iproc, nproc, lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f, &
      !!             lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f, lorbs%norb, lorbs%nspinor, 1, &
      !!             lzd%Glr, lzd%Llr(ilr), lphi(ist), lvphiovrlp(1))
      !!        call plotOrbitals(iproc, orbs, lzd%Glr, lvphiovrlp, at%nat, rxyz, lorbs%inwhichlocreg, &
      !!             .5d0*hx, .5d0*hy, .5d0*hz, it)
      !!        ist=ist+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      !!    end do
      !!    deallocate(lvphiovrlp)
      !!end if
      !!! ####################################################################################




     ! Flush the standard output
     !flush(unit=6) 


  end do iterLoop




  if(newgradient) then
      !!! Create new logrecs taking incto account the derivatives
      !!ii=lorbs%npsidim_orbs
      !!call dcopy(ii, lphi(1), 1, lphilarge(1), 1)

      !!call destroy_new_locregs(lzd, lorbs, op, comon, mad, comgp, &
      !!     lphi, lhphi, lhphiold, lphiold)
      !!call create_new_locregs(iproc, nproc, lzdlarge%nlr, hx, hy, hz, orbslarge, lzdlarge%glr, locregCenterTemp, &
      !!     lzdlarge%llr(:)%locrad, denspot%dpcom%nscatterarr, .false., ldiis, &
      !!     lzd, lorbs, op, comon, mad, comgp, &
      !!     lphi, lhphi, lhphiold, lphiold)

      !!call dcopy(ii, lphilarge(1), 1, lphi(1), 1)

      call destroy_new_locregs(lzdlarge, orbslarge, oplarge, comonlarge, madlarge, comgplarge, &
           lphilarge, lhphilarge, lhphilargeold, lphilargeold)
      if(secondLocreg) then
          call destroy_new_locregs(lzdlarge2, orbslarge2, oplarge2, comonlarge2, madlarge2, comgplarge2, &
               lphilarge2, lhphilarge2, lhphilarge2old, lphilarge2old)
          call deallocateDIIS(ldiis2)
      end if

      ! Write the locreg centers
      if(iproc==0) then
          write(*,'(1x,a)') 'the new centers of the localization regions:'
          do ilr=1,lzd%nlr
              write(*,'(3x,i7,3es24.12)') ilr, lzd%llr(ilr)%locregCenter(1:3)
          end do
      end if
  end if

  ! thi svalue cannot be determined


  iall=-product(shape(lphiold))*kind(lphiold)
  deallocate(lphiold, stat=istat)
  call memocc(istat, iall, 'lphiold', subname)
  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)
  iall=-product(shape(Umat))*kind(Umat)
  deallocate(Umat, stat=istat)
  call memocc(istat, iall, 'Umat', subname)
  iall=-product(shape(kernelold))*kind(kernelold)
  deallocate(kernelold, stat=istat)
  call memocc(istat, iall, 'kernelold', subname)
  iall=-product(shape(locregCenter))*kind(locregCenter)
  deallocate(locregCenter, stat=istat)
  call memocc(istat, iall, 'locregCenter', subname)

  iall=-product(shape(locrad))*kind(locrad)
  deallocate(locrad, stat=istat)
  call memocc(istat, iall, 'locrad', subname)
  iall=-product(shape(locrad_tmp))*kind(locrad_tmp)
  deallocate(locrad_tmp, stat=istat)
  call memocc(istat, iall, 'locrad_tmp', subname)
  iall=-product(shape(inwhichlocreg_reference))*kind(inwhichlocreg_reference)
  deallocate(inwhichlocreg_reference, stat=istat)
  call memocc(istat, iall, 'inwhichlocreg_reference', subname)


!!$  iall=-product(shape(lorbs%ispot))*kind(lorbs%ispot)
!!$  deallocate(lorbs%ispot, stat=istat)
!!$  call memocc(istat, iall, 'lorbs%ispot', subname)

  !!! Deallocate potential
  !!iall=-product(shape(denspot%pot_full))*kind(denspot%pot_full)
  !!deallocate(denspot%pot_full, stat=istat)
  !!call memocc(istat, iall, 'denspot%pot_full', subname)
  if(.not. newgradient) then
      ! Deallocate potential
      iall=-product(shape(denspot%pot_full))*kind(denspot%pot_full)
      deallocate(denspot%pot_full, stat=istat)
      call memocc(istat, iall, 'denspot%pot_full', subname)
  end if


  ! Deallocate PSP stuff
  !call free_lnlpspd(lorbs, lzd)

  t2tot=mpi_wtime()
  timetot=t2tot-t1tot

  ! Sum up the timings.
  call mpiallred(time(1), 5, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(timetot, 1, mpi_sum, mpi_comm_world, ierr)
  time=time/dble(nproc)
  timetot=timetot/dble(nproc)
  if(iproc==0) then
      write(*,'(1x,a)') 'timings:'
      write(*,'(3x,a,es10.3)') '-total time:', timetot
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- orthonormalization:', time(1), '=', 100.d0*time(1)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- Hamiltonian application:', time(2),  '=', 100.d0*time(2)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- orthoconstraint:', time(3),  '=', 100.d0*time(3)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- preconditioning:', time(4),  '=', 100.d0*time(4)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- unitary optimization:', time(5),  '=', 100.d0*time(5)/timetot, '%'
      tt=time(1)+time(2)+time(3)+time(4)+time(5)
      tt=timetot-tt
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- other:', tt,  '=', 100.d0*tt/timetot, '%'
  end if


  ! Deallocate all quantities related to DIIS,
  !!!!if(ldiis%isx>0) call deallocateDIIS(ldiis)

  ! Deallocate all local arrays.
  call deallocateLocalArrays()

  iall=-product(shape(kernel))*kind(kernel)
  deallocate(kernel, stat=istat)
  call memocc(istat, iall, 'kernel', subname)


contains



    subroutine my_geocode_buffers(geocode,nl1,nl2,nl3)
      implicit none
      integer, intent(in) :: nl1,nl2,nl3
      character(len=1), intent(in) :: geocode
      !local variables
      logical :: perx,pery,perz
      integer :: nr1,nr2,nr3

      !conditions for periodicity in the three directions
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')

      call ext_buffers(perx,nl1,nr1)
      call ext_buffers(pery,nl2,nr2)
      call ext_buffers(perz,nl3,nr3)

    end subroutine my_geocode_buffers






    subroutine DIISorSD()
    !
    ! Purpose:
    ! ========
    !   This subroutine decides whether one should use DIIS or variable step size
    !   steepest descent to improve the orbitals. In the beginning we start with DIIS
    !   with history length lin%DIISHistMax. If DIIS becomes unstable, we switch to
    !   steepest descent. If the steepest descent iterations are successful, we switch
    !   back to DIIS, but decrease the DIIS history length by one. However the DIIS
    !   history length is limited to be larger or equal than lin%DIISHistMin.
    !


      ! First there are some checks whether the force is small enough to allow DIIS.

      ! Decide whether the force is small eneough to allow DIIS
      !if(fnrmMax<lin%startDIIS .and. .not.allowDIIS) then
      !    allowDIIS=.true.
      !    if(iproc==0) write(*,'(1x,a)') 'The force is small enough to allow DIIS.'
      !    ! This is to get the correct DIIS history 
      !    ! (it is chosen as max(lin%DIISHistMin,lin%DIISHistMax-icountSwitch).
      !    icountSwitch=icountSwitch-1
      !else if(fnrmMax>lin%startDIIS .and. allowDIIS) then
      !    allowDIIS=.false.
      !    if(iproc==0) write(*,'(1x,a)') 'The force is too large to allow DIIS.'
      !end if    

      !! Switch to SD if the flag indicating that we should start with SD is true.
      !! If this is the case, this flag is set to false, since this flag concerns only the beginning.
      !if(startWithSD .and. ldiis%isx>0) then
      !    call deallocateDIIS(ldiis)
      !    ldiis%isx=0
      !    ldiis%switchSD=.false.
      !    startWithSD=.false.
      !end if

      ! Decide whether we should switch from DIIS to SD in case we are using DIIS and it 
      ! is not allowed.
      !if(.not.startWithSD .and. .not.allowDIIS .and. ldiis%isx>0) then
      !if(.not.startWithSD .and. ldiis%isx>0) then
      !if(ldiis%isx>0) then
      !    if(iproc==0) write(*,'(1x,a,es10.3)') 'The force is too large, switch to SD with stepsize', alpha(1)
      !    call deallocateDIIS(ldiis)
      !    ldiis%isx=0
      !    ldiis%switchSD=.true.
      !end if

      ! If we swicthed to SD in the previous iteration, reset this flag.
      if(ldiis%switchSD) ldiis%switchSD=.false.

      ! Now come some checks whether the trace is descreasing or not. This further decides
      ! whether we should use DIIS or SD.

      ! Determine wheter the trace is decreasing (as it should) or increasing.
      ! This is done by comparing the current value with diisLIN%energy_min, which is
      ! the minimal value of the trace so far.
      if(trH<=ldiis%trmin .and. .not.resetDIIS) then
          ! Everything ok
          ldiis%trmin=trH
          ldiis%switchSD=.false.
          itBest=it
          icountSDSatur=icountSDSatur+1
          icountDIISFailureCons=0

          ! If we are using SD (i.e. diisLIN%idsx==0) and the trace has been decreasing
          ! for at least 10 iterations, switch to DIIS. However the history length is decreased.
          !if(icountSDSatur>=10 .and. ldiis%isx==0 .and. allowDIIS .or. immediateSwitchToSD) then
          if(icountSDSatur>=10 .and. ldiis%isx==0 .or. immediateSwitchToSD) then
              icountSwitch=icountSwitch+1
              idsx=max(ldiis%DIISHistMin,ldiis%DIISHistMax-icountSwitch)
              if(idsx>0) then
                  if(iproc==0) write(*,'(1x,a,i0)') 'switch to DIIS with new history length ', idsx
                  icountSDSatur=0
                  icountSwitch=0
                  icountDIISFailureTot=0
                  icountDIISFailureCons=0
                  ldiis%is=0
                  ldiis%switchSD=.false.
                  ldiis%trmin=1.d100
                  ldiis%trold=1.d100
                  alpha=ldiis%alphaSD
                  alphaDIIS=ldiis%alphaDIIS
                  !!!call initializeDIIS(lin%DIISHistMax, lzd, lorbs, lorbs%norb, ldiis)
                  icountDIISFailureTot=0
                  icountDIISFailureCons=0
                  immediateSwitchToSD=.false.
              end if
          end if
      else
          ! The trace is growing.
          ! Count how many times this occurs and (if we are using DIIS) switch to SD after 3 
          ! total failures or after 2 consecutive failures.
          icountDIISFailureCons=icountDIISFailureCons+1
          icountDIISFailureTot=icountDIISFailureTot+1
          icountSDSatur=0
          if((icountDIISFailureCons>=2 .or. icountDIISFailureTot>=3 .or. resetDIIS) .and. ldiis%isx>0) then
          !if((icountDIISFailureCons>=200 .or. icountDIISFailureTot>=300 .or. resetDIIS) .and. ldiis%isx>0) then
          !if((icountDIISFailureCons>=400 .or. icountDIISFailureTot>=600 .or. resetDIIS) .and. ldiis%isx>0) then
              ! Switch back to SD.
              alpha=ldiis%alphaSD
              if(iproc==0) then
                  if(icountDIISFailureCons>=2) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                      icountDIISFailureCons, ' times consecutively. Switch to SD with stepsize', alpha(1)
                  if(icountDIISFailureTot>=3) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                      icountDIISFailureTot, ' times in total. Switch to SD with stepsize', alpha(1)
                  if(resetDIIS) write(*,'(1x,a)') 'reset DIIS due to flag'
              end if
              if(resetDIIS) then
                  resetDIIS=.false.
                  immediateSwitchToSD=.true.
                  ldiis%trmin=1.d100
              end if
              ! Try to get back the orbitals of the best iteration. This is possible if
              ! these orbitals are still present in the DIIS history.
              if(it-itBest<ldiis%isx) then
                 if(iproc==0) then
                     if(iproc==0) write(*,'(1x,a,i0,a)')  'Recover the orbitals from iteration ', &
                         itBest, ' which are the best so far.'
                 end if
                 ii=modulo(ldiis%mis-(it-itBest),ldiis%mis)
                 offset=0
                 istdest=1
                 do iorb=1,lorbs%norbp
                     !ilr=lorbs%inWhichLocregp(iorb)
                     if(.not.newgradient) then
                         iiorb=lorbs%isorb+iorb
                         ilr=lorbs%inWhichLocreg(iiorb)
                         ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                     else
                         iiorb=orbslarge%isorb+iorb
                         ilr=orbslarge%inWhichLocreg(iiorb)
                         ncount=lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
                     end if
                     istsource=offset+ii*ncount+1
                     !write(*,'(a,4i9)') 'iproc, ncount, istsource, istdest', iproc, ncount, istsource, istdest
                     if(.not.newgradient) then
                         call dcopy(ncount, ldiis%phiHist(istsource), 1, lphi(istdest), 1)
                         call dcopy(ncount, ldiis%phiHist(istsource), 1, lphiold(istdest), 1)
                     else
                         call dcopy(ncount, ldiis%phiHist(istsource), 1, lphilarge(istdest), 1)
                         call dcopy(ncount, ldiis%phiHist(istsource), 1, lphilargeold(istdest), 1)
                     end if
                     offset=offset+ldiis%isx*ncount
                     istdest=istdest+ncount
                 end do
             else
                 ! else copy the orbitals of the last iteration to lphiold
                 if(.not.newgradient) then
                     call dcopy(size(lphi), lphi(1), 1, lphiold(1), 1)
                 else
                     call dcopy(size(lphilarge), lphilarge(1), 1, lphilargeold(1), 1)
                 end if
              end if
              !!!call deallocateDIIS(ldiis)
              ldiis%isx=0
              ldiis%switchSD=.true.
          end if
      end if

    end subroutine DIISorSD


    subroutine improveOrbitals()
    !
    ! Purpose:
    ! ========
    !   This subroutine improves the basis functions by following the gradient 
    ! For DIIS 
    !!if (diisLIN%idsx > 0) then
    !!   diisLIN%mids=mod(diisLIN%ids,diisLIN%idsx)+1
    !!   diisLIN%ids=diisLIN%ids+1
    !!end if
    if (ldiis%isx > 0) then
        ldiis%mis=mod(ldiis%is,ldiis%isx)+1
        ldiis%is=ldiis%is+1
    end if

    ! Follow the gradient using steepest descent.
    ! The same, but transposed
    
    ! steepest descent
    if(ldiis%isx==0) then
        call timing(iproc,'optimize_SD   ','ON')
        istart=1
        do iorb=1,lorbs%norbp
            !ilr=lorbs%inWhichLocregp(iorb)
            if(.not.newgradient) then
                iiorb=lorbs%isorb+iorb
                ilr=lorbs%inWhichLocreg(iiorb)
                ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                call daxpy(ncount, -alpha(iorb), lhphi(istart), 1, lphi(istart), 1)
            else
                iiorb=orbslarge%isorb+iorb
                ilr=orbslarge%inWhichLocreg(iiorb)
                ncount=lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
                call daxpy(ncount, -alpha(iorb), lhphilarge(istart), 1, lphilarge(istart), 1)
            end if
            istart=istart+ncount
        end do
        call timing(iproc,'optimize_SD   ','OF')
    else
        ! DIIS
        if(ldiis%alphaDIIS/=1.d0) then
            if(.not.newgradient) then
                call dscal(max(lorbs%npsidim_orbs,lorbs%npsidim_comp), ldiis%alphaDIIS, lhphi, 1)
            else
                call dscal(max(orbslarge%npsidim_orbs,orbslarge%npsidim_comp), ldiis%alphaDIIS, lhphilarge, 1)
            end if
        end if
        if(.not.newgradient) then
            call optimizeDIIS(iproc, nproc, lorbs, lorbs, lzd, lhphi, lphi, ldiis, it)
        else
            call optimizeDIIS(iproc, nproc, orbslarge, orbslarge, lzdlarge, lhphilarge, lphilarge, ldiis, it)
        end if
    end if
    end subroutine improveOrbitals



    subroutine allocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine allocates all local arrays.
    !
      allocate(alpha(lorbs%norb), stat=istat)
      call memocc(istat, alpha, 'alpha', subname)

      allocate(alphaDIIS(lorbs%norb), stat=istat)
      call memocc(istat, alphaDIIS, 'alphaDIIS', subname)

      allocate(fnrmArr(lorbs%norb,2), stat=istat)
      call memocc(istat, fnrmArr, 'fnrmArr', subname)

      allocate(fnrmOldArr(lorbs%norb), stat=istat)
      call memocc(istat, fnrmOldArr, 'fnrmOldArr', subname)

      allocate(fnrmOvrlpArr(lorbs%norb,2), stat=istat)
      call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)

      !allocate(lhphi(lin%Lorbs%npsidim), stat=istat)
      allocate(lhphi(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)), stat=istat)
      call memocc(istat, lhphi, 'lhphi', subname)
    
      !allocate(lhphiold(lin%Lorbs%npsidim), stat=istat)
      allocate(lhphiold(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)), stat=istat)
      call memocc(istat, lhphiold, 'lhphiold', subname)

    end subroutine allocateLocalArrays


    subroutine deallocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine deallocates all local arrays.
    !
      iall=-product(shape(alpha))*kind(alpha)
      deallocate(alpha, stat=istat)
      call memocc(istat, iall, 'alpha', subname)

      iall=-product(shape(alphaDIIS))*kind(alphaDIIS)
      deallocate(alphaDIIS, stat=istat)
      call memocc(istat, iall, 'alphaDIIS', subname)

      iall=-product(shape(fnrmArr))*kind(fnrmArr)
      deallocate(fnrmArr, stat=istat)
      call memocc(istat, iall, 'fnrmArr', subname)

      iall=-product(shape(fnrmOldArr))*kind(fnrmOldArr)
      deallocate(fnrmOldArr, stat=istat)
      call memocc(istat, iall, 'fnrmOldArr', subname)

      iall=-product(shape(fnrmOvrlpArr))*kind(fnrmOvrlpArr)
      deallocate(fnrmOvrlpArr, stat=istat)
      call memocc(istat, iall, 'fnrmOvrlpArr', subname)

      iall=-product(shape(lhphi))*kind(lhphi)
      deallocate(lhphi, stat=istat)
      call memocc(istat, iall, 'lhphi', subname)

      iall=-product(shape(lhphiold))*kind(lhphiold)
      deallocate(lhphiold, stat=istat)
      call memocc(istat, iall, 'lhphiold', subname)

      ! if diisLIN%idsx==0, these arrays have already been deallocated
      !if(diisLIN%idsx>0 .and. lin%DIISHistMax>0) call deallocate_diis_objects(diisLIN,subname)

    end subroutine deallocateLocalArrays


end subroutine getLocalizedBasis



subroutine create_new_locregs(iproc, nproc, nlr, hx, hy, hz, lorbs, glr, locregCenter, locrad, nscatterarr, withder, &
           inwhichlocreg_reference, ldiis, lzdlarge, orbslarge, oplarge, comonlarge, madlarge, comgplarge, &
           lphilarge, lhphilarge, lhphilargeold, lphilargeold)
use module_base
use module_types
use module_interfaces, except_this_one => create_new_locregs
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nlr
real(8),intent(in):: hx, hy, hz
type(orbitals_data),intent(in):: lorbs
type(locreg_descriptors),intent(in):: glr
real(8),dimension(3,nlr),intent(in):: locregCenter
real(8),dimension(nlr):: locrad
integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
logical,intent(in):: withder
integer,dimension(lorbs%norb),intent(in):: inwhichlocreg_reference
type(localizedDIISParameters),intent(inout):: ldiis
type(local_zone_descriptors),intent(out):: lzdlarge
type(orbitals_data),intent(out):: orbslarge
type(overlapParameters),intent(out):: oplarge
type(p2pComms),intent(out):: comonlarge
type(matrixDescriptors),intent(out):: madlarge
type(p2pComms),intent(out):: comgplarge
real(8),dimension(:),pointer,intent(out):: lphilarge, lhphilarge, lhphilargeold, lphilargeold

! Local variables
integer:: tag, norbu, norbd, nspin, iorb, iiorb, ilr, npsidim, ii, istat, iall, ierr, norb
integer,dimension(:),allocatable:: orbsperlocreg
character(len=*),parameter:: subname='create_new_locregs'


   if(iproc==0) write(*,'(x,a)') 'creating new locregs...'
   call nullify_local_zone_descriptors(lzdlarge)
   call nullify_orbitals_data(orbslarge)
   call nullify_overlapParameters(oplarge)
   call nullify_p2pComms(comonlarge)
   call nullify_matrixDescriptors(madlarge)
   call nullify_p2pComms(comgplarge)

   tag=1
   lzdlarge%nlr=nlr
   if(.not.withder) then
       norbu=lorbs%norb
   else
       norbu=4*lorbs%norb
   end if
   norb=norbu
   norbd=0
   nspin=1
   call orbitals_descriptors_forLinear(iproc, nproc, norb, norbu, norbd, nspin, lorbs%nspinor,&
        lorbs%nkpts, lorbs%kpts, lorbs%kwgts, orbslarge)
   call repartitionOrbitals(iproc, nproc, orbslarge%norb, orbslarge%norb_par,&
        orbslarge%norbp, orbslarge%isorb_par, orbslarge%isorb, orbslarge%onWhichMPI)

   !!!!allocate(orbsperlocreg(lzdlarge%nlr), stat=istat)
   !!!!call memocc(istat, orbsperlocreg, 'orbsperlocreg', subname)
   !!!!!!orbsperlocreg=0
   !!!!!!do iorb=1,lorbs%norbp
   !!!!!!    iiorb=lorbs%isorb+iorb
   !!!!!!    ilr=iiorb
   !!!!!!    orbsperlocreg(ilr)=orbsperlocreg(ilr)+1
   !!!!!!end do
   !!!!!!call mpiallred(orbsperlocreg(1), lzdlarge%nlr, mpi_sum, mpi_comm_world, ierr)
   !!!!do iorb=1,lzdlarge%nlr
   !!!!    orbsperlocreg(iorb)=1
   !!!!end do
   !!!!iall=-product(shape(orbslarge%inWhichLocreg))*kind(orbslarge%inWhichLocreg)
   !!!!deallocate(orbslarge%inWhichLocreg, stat=istat)
   !!!!call memocc(istat, iall, 'orbslarge%inWhichLocreg', subname)
   !!!!call assignToLocreg2(iproc, nproc, orbslarge%norb, orbslarge%norb_par, 0, lzdlarge%nlr, &
   !!!!     nspin, orbsperlocreg, locregCenter, orbslarge%inwhichlocreg)
   !!!iall=-product(shape(orbsperlocreg))*kind(orbsperlocreg)
   !!!deallocate(orbsperlocreg, stat=istat)
   !!!call memocc(istat, iall, 'orbsperlocreg', subname)
   orbslarge%inwhichlocreg = inwhichlocreg_reference

   call initLocregs(iproc, nproc, lzdlarge%nlr, locregCenter, hx, hy, hz, lzdlarge, orbslarge, Glr, locrad, 's')
   call nullify_locreg_descriptors(lzdlarge%Glr)
   call copy_locreg_descriptors(Glr, lzdlarge%Glr, subname)
   npsidim = 0
   do iorb=1,orbslarge%norbp
    ilr=orbslarge%inwhichlocreg(iorb+orbslarge%isorb)
    npsidim = npsidim + lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
   end do
   allocate(orbslarge%eval(orbslarge%norb), stat=istat)
   call memocc(istat, orbslarge%eval, 'orbslarge%eval', subname)
   orbslarge%eval=-.5d0
   orbslarge%npsidim_orbs=max(npsidim,1)
   call initCommsOrtho(iproc, nproc, nspin, hx, hy, hz, lzdlarge, orbslarge, orbslarge%inWhichLocreg,&
        's', oplarge, comonlarge, tag)
   call initMatrixCompression(iproc, nproc, lzdlarge%nlr, orbslarge, &
        oplarge%noverlaps, oplarge%overlaps, madlarge)
   call initCompressedMatmul3(orbslarge%norb, madlarge)

   call initializeCommunicationPotential(iproc, nproc, nscatterarr, orbslarge, lzdlarge, comgplarge, orbslarge%inWhichLocreg, tag)

   iall=-product(shape(ldiis%phiHist))*kind(ldiis%phiHist)
   deallocate(ldiis%phiHist, stat=istat)
   call memocc(istat, iall, 'ldiis%phiHist', subname)
   iall=-product(shape(ldiis%hphiHist))*kind(ldiis%hphiHist)
   deallocate(ldiis%hphiHist, stat=istat)
   call memocc(istat, iall, 'ldiis%hphiHist', subname)
   ii=0
   do iorb=1,orbslarge%norbp
       ilr=orbslarge%inwhichlocreg(orbslarge%isorb+iorb)
       ii=ii+ldiis%isx*(lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f)
   end do
   allocate(ldiis%phiHist(ii), stat=istat)
   call memocc(istat, ldiis%phiHist, 'ldiis%phiHist', subname)
   allocate(ldiis%hphiHist(ii), stat=istat)
   call memocc(istat, ldiis%hphiHist, 'ldiis%hphiHist', subname)

   allocate(lphilarge(orbslarge%npsidim_orbs), stat=istat)
   call memocc(istat, lphilarge, 'lphilarge', subname)
   allocate(lhphilarge(orbslarge%npsidim_orbs), stat=istat)
   call memocc(istat, lhphilarge, 'lhphilarge', subname)
   allocate(lhphilargeold(orbslarge%npsidim_orbs), stat=istat)
   call memocc(istat, lhphilargeold, 'lhphilargeold', subname)
   allocate(lphilargeold(orbslarge%npsidim_orbs), stat=istat)
   call memocc(istat, lphilargeold, 'lphilargeold', subname)

   lphilarge=0.d0
   lhphilarge=0.d0
   lhphilargeold=0.d0
   lphilargeold=0.d0

end subroutine create_new_locregs



subroutine destroy_new_locregs(lzdlarge, orbslarge, oplarge, comonlarge, madlarge, comgplarge, &
           lphilarge, lhphilarge, lhphilargeold, lphilargeold)
  use module_base
  use module_types
  use module_interfaces, except_this_one => destroy_new_locregs
  implicit none

  ! Calling arguments
  type(local_zone_descriptors),intent(inout):: lzdlarge
  type(orbitals_data),intent(inout):: orbslarge
  type(overlapParameters),intent(inout):: oplarge
  type(p2pComms),intent(inout):: comonlarge
  type(matrixDescriptors),intent(inout):: madlarge
  type(p2pComms),intent(inout):: comgplarge
  real(8),dimension(:),pointer,intent(inout):: lphilarge, lhphilarge, lhphilargeold, lphilargeold

  ! Local variables
  integer:: istat, iall
  character(len=*),parameter:: subname='destroy_new_locregs'

  call deallocate_local_zone_descriptors(lzdlarge, subname)
  call deallocate_orbitals_data(orbslarge, subname)
  call deallocate_overlapParameters(oplarge, subname)
  call deallocate_p2pComms(comonlarge, subname)
  call deallocate_matrixDescriptors(madlarge, subname)
  call deallocate_p2pComms(comgplarge, subname)

  iall=-product(shape(lphilarge))*kind(lphilarge)
  deallocate(lphilarge, stat=istat)
  call memocc(istat, iall, 'lphilarge', subname)
  iall=-product(shape(lhphilarge))*kind(lhphilarge)
  deallocate(lhphilarge, stat=istat)
  call memocc(istat, iall, 'lhphilarge', subname)
  iall=-product(shape(lhphilargeold))*kind(lhphilargeold)
  deallocate(lhphilargeold, stat=istat)
  call memocc(istat, iall, 'lhphilargeold', subname)
  iall=-product(shape(lphilargeold))*kind(lphilargeold)
  deallocate(lphilargeold, stat=istat)
  call memocc(istat, iall, 'lphilargeold', subname)

end subroutine destroy_new_locregs






subroutine transformHam(iproc, nproc, orbs, comms, phi, hphi, HamSmall)
!
! Purpose:
! =======
!   Builds the Hamiltonian in the basis of the localized basis functions phi. To do so, it gets all basis
!   functions |phi_i> and H|phi_i> and then calculates H_{ij}=<phi_i|H|phi_j>. The basis functions phi are
!   provided in the transposed form.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc      process ID
!     nproc      total number of processes
!     orbs       type describing the basis functions psi
!     comms      type containing the communication parameters for the physical orbitals phi
!     phi        basis functions 
!     hphi       the Hamiltonian applied to the basis functions 
!   Output arguments:
!   -----------------
!     HamSmall   Hamiltonian in small basis
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data), intent(in) :: orbs
type(communications_arrays), intent(in) :: comms
real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb), intent(in) :: phi, hphi
real(8),dimension(orbs%norb,orbs%norb),intent(out):: HamSmall

! Local variables
integer:: istat, ierr, nvctrp, iall
real(8),dimension(:,:),allocatable:: HamTemp
character(len=*),parameter:: subname='transformHam'



  ! Allocate a temporary array if there are several MPI processes
  if(nproc>1) then
      allocate(HamTemp(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, HamTemp, 'HamTemp', subname)
  end if
  
  ! nvctrp is the amount of each phi hold by the current process
  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  
  ! Build the Hamiltonian. In the parallel case, each process writes its Hamiltonian in HamTemp
  ! and a mpi_allreduce sums up the contribution from all processes.
  if(nproc==1) then
      call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi(1,1), nvctrp, &
                 hphi(1,1), nvctrp, 0.d0, HamSmall(1,1), orbs%norb)
  else
      call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi(1,1), nvctrp, &
                 hphi(1,1), nvctrp, 0.d0, HamTemp(1,1), orbs%norb)
  end if
  if(nproc>1) then
      call mpi_allreduce(HamTemp(1,1), HamSmall(1,1), orbs%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  end if
  
  if(nproc>1) then
     iall=-product(shape(HamTemp))*kind(HamTemp)
     deallocate(HamTemp,stat=istat)
     call memocc(istat, iall, 'HamTemp', subname)
  end if

end subroutine transformHam




subroutine diagonalizeHamiltonian(iproc, nproc, orbs, HamSmall, eval)
!
! Purpose:
! ========
!   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
!   the same result. This is done by requiring that the first entry of each vector
!   is positive.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc     process ID
!     nproc     number of MPI processes
!     orbs      type describing the physical orbitals psi
!   Input / Putput arguments
!     HamSmall  on input: the Hamiltonian
!               on exit: the eigenvectors
!   Output arguments
!     eval      the associated eigenvalues 
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc
type(orbitals_data), intent(inout) :: orbs
real(8),dimension(orbs%norb, orbs%norb):: HamSmall
real(8),dimension(orbs%norb):: eval

! Local variables
integer:: lwork, info, istat, iall, i, iorb, jorb
real(8),dimension(:),allocatable:: work
character(len=*),parameter:: subname='diagonalizeHamiltonian'

  ! Get the optimal work array size
  lwork=-1 
  allocate(work(1), stat=istat)
  call memocc(istat, work, 'work', subname)
  call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  lwork=work(1) 

  ! Deallocate the work array ane reallocate it with the optimal size
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
  call memocc(istat, work, 'work', subname)

  ! Diagonalize the Hamiltonian
  call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 

  ! Deallocate the work array.
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  
  ! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
  ! the first entry of each vector is positive.
  do iorb=1,orbs%norb
      if(HamSmall(1,iorb)<0.d0) then
          do jorb=1,orbs%norb
              HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
          end do
      end if
  end do


end subroutine diagonalizeHamiltonian



subroutine diagonalizeHamiltonian2(iproc, nproc, orbs, nsubmax, HamSmall, ovrlp, eval)
!
! Purpose:
! ========
!   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
!   the same result. This is done by requiring that the first entry of each vector
!   is positive.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc     process ID
!     nproc     number of MPI processes
!     orbs      type describing the physical orbitals psi
!   Input / Putput arguments
!     HamSmall  on input: the Hamiltonian
!               on exit: the eigenvectors
!   Output arguments
!     eval      the associated eigenvalues 
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc, nsubmax
type(orbitals_data), intent(inout) :: orbs
real(8),dimension(orbs%norb, orbs%norb):: HamSmall, ovrlp
real(8),dimension(orbs%norb):: eval

! Local variables
integer:: lwork, info, istat, iall, i, iorb, jorb, nsub
real(8),dimension(:),allocatable:: work
real(8),dimension(:,:),allocatable:: ham_band, ovrlp_band
character(len=*),parameter:: subname='diagonalizeHamiltonian'

  call timing(iproc,'diagonal_seq  ','ON')

  !! OLD VERSION #####################################################################################################
  ! Get the optimal work array size
  lwork=-1 
  allocate(work(1), stat=istat)
  call memocc(istat, work, 'work', subname)
  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  lwork=work(1) 

  ! Deallocate the work array ane reallocate it with the optimal size
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
  call memocc(istat, work, 'work', subname)

  ! Diagonalize the Hamiltonian
  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 

  ! Deallocate the work array.
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  
  ! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
  ! the first entry of each vector is positive.
  do iorb=1,orbs%norb
      if(HamSmall(1,iorb)<0.d0) then
          do jorb=1,orbs%norb
              HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
          end do
      end if
  end do
  !! #################################################################################################################

  !!!! NEW VERSION #####################################################################################################
  !!! Determine the maximal number of non-zero subdiagonals
  !!!!nsubmax=0
  !!!!do iorb=1,orbs%norb
  !!!!    nsub=0
  !!!!    do jorb=orbs%norb,iorb+1,-1
  !!!!        if(Hamsmall(jorb,iorb)/=0.d0) then
  !!!!            nsub=jorb-iorb
  !!!!            exit
  !!!!        end if
  !!!!    end do
  !!!!    if(iproc==0) write(*,*) 'iorb,nsub',iorb,nsub
  !!!!    nsubmax=max(nsub,nsubmax)
  !!!!end do
  !!!!if(iproc==0) write(*,*) 'nsubmax',nsubmax
  !!!!if(iproc==0) then
  !!!!      do iorb=1,orbs%norb
  !!!!           write(*,'(14es10.3)') (hamsmall(iorb,jorb), jorb=1,orbs%norb)
  !!!!      end do
  !!!!end if

  !!! Copy to banded format
  !!allocate(ham_band(nsubmax+1,orbs%norb), stat=istat)
  !!call memocc(istat, ham_band, 'ham_band', subname)
  !!allocate(ovrlp_band(nsubmax+1,orbs%norb), stat=istat)
  !!call memocc(istat, ovrlp_band, 'ovrlp_band', subname)
  !!do iorb=1,orbs%norb
  !!    do jorb=iorb,min(iorb+nsubmax,orbs%norb)
  !!        ham_band(1+jorb-iorb,iorb)=HamSmall(jorb,iorb)
  !!        ovrlp_band(1+jorb-iorb,iorb)=ovrlp(jorb,iorb)
  !!    end do
  !!end do
  !!!!if(iproc==0) then
  !!!!      write(*,*) '+++++++++++++++++++++++++++++'
  !!!!      do iorb=1,nsubmax+1
  !!!!           write(*,'(14es10.3)') (ham_band(iorb,jorb), jorb=1,orbs%norb)
  !!!!      end do
  !!!!end if



  !!!!! Get the optimal work array size
  !!!!lwork=-1 
  !!!!allocate(work(1), stat=istat)
  !!!!call memocc(istat, work, 'work', subname)
  !!!!call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  !!!!lwork=work(1) 

  !!!!! Deallocate the work array ane reallocate it with the optimal size
  !!!!iall=-product(shape(work))*kind(work)
  !!!!deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  !!!!call memocc(istat, iall, 'work', subname)
  !!allocate(work(3*orbs%norb), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
  !!call memocc(istat, work, 'work', subname)

  !!! Diagonalize the Hamiltonian
  !!!call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  !!call dsbgv('v', 'l', orbs%norb, nsubmax, nsubmax, ham_band(1,1), nsubmax+1, ovrlp_band(1,1), nsubmax+1, &
  !!     eval(1), HamSmall(1,1), orbs%norb, work, info)

  !!! Deallocate the work array.
  !!iall=-product(shape(work))*kind(work)
  !!deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  !!call memocc(istat, iall, 'work', subname)
  !!
  !!! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
  !!! the first entry of each vector is positive.
  !!do iorb=1,orbs%norb
  !!    if(HamSmall(1,iorb)<0.d0) then
  !!        do jorb=1,orbs%norb
  !!            HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
  !!        end do
  !!    end if
  !!end do


  !!iall=-product(shape(ham_band))*kind(ham_band)
  !!deallocate(ham_band, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating ham_band' 
  !!call memocc(istat, iall, 'ham_band', subname)

  !!iall=-product(shape(ovrlp_band))*kind(ovrlp_band)
  !!deallocate(ovrlp_band, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating ovrlp_band' 
  !!call memocc(istat, iall, 'ovrlp_band', subname)

  call timing(iproc,'diagonal_seq  ','OF')

end subroutine diagonalizeHamiltonian2




subroutine buildWavefunction(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, HamSmall)
!
! Purpose:
! =======
!   Builds the physical orbitals psi as a linear combination of the basis functions phi. The coefficients
!   for this linear combination are obtained by diagonalizing the Hamiltonian matrix HamSmall.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc      process ID
!     nproc      total number of processes
!     orbs       type describing the physical orbitals psi
!     orbsLIN    type describing the basis functions phi
!     comms      type containing the communication parameters for the physical orbitals psi
!     commsLIN   type containing the communication parameters for the basis functions phi
!     phi        the basis functions 
!     HamSmall   the  Hamiltonian matrix
!   Output arguments:
!   -----------------
!     psi        the physical orbitals 
!

use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc
type(orbitals_data), intent(in) :: orbs
type(orbitals_data), intent(in) :: orbsLIN
type(communications_arrays), intent(in) :: comms
type(communications_arrays), intent(in) :: commsLIN
real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbsLIN%norb) :: phi
real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb) :: psi
real(8),dimension(orbsLIN%norb,orbsLIN%norb):: HamSmall

! Local variables
integer:: nvctrp


  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  call dgemm('n', 'n', nvctrp, orbs%norb, orbsLIN%norb, 1.d0, phi(1,1), nvctrp, HamSmall(1,1), &
             orbsLIN%norb, 0.d0, psi(1,1), nvctrp)
  

end subroutine buildWavefunction





subroutine buildWavefunctionModified(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, coeff)

!
! Purpose:
! =======
!   Builds the physical orbitals psi as a linear combination of the basis functions phi. The coefficients
!   for this linear combination are obtained by diagonalizing the Hamiltonian matrix HamSmall.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc      process ID
!     nproc      total number of processes
!     orbs       type describing the physical orbitals psi
!     orbsLIN    type describing the basis functions phi
!     comms      type containing the communication parameters for the physical orbitals psi
!     commsLIN   type containing the communication parameters for the basis functions phi
!     phi        the basis functions 
!     coeff      the coefficients for the linear combination
!   Output arguments:
!   -----------------
!     psi        the physical orbitals 
!

use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc
type(orbitals_data), intent(in) :: orbs
type(orbitals_data), intent(in) :: orbsLIN
type(communications_arrays), intent(in) :: comms
type(communications_arrays), intent(in) :: commsLIN
real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbsLIN%norb) :: phi
real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb) :: psi
real(8),dimension(orbsLIN%norb,orbs%norb):: coeff

! Local variables
integer:: nvctrp


  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  call dgemm('n', 'n', nvctrp, orbs%norb, orbsLIN%norb, 1.d0, phi(1,1), nvctrp, coeff(1,1), &
             orbsLIN%norb, 0.d0, psi(1,1), nvctrp)
  

end subroutine buildWavefunctionModified


!!$subroutine getMatrixElements(iproc, nproc, Glr, orbs, comms, phi, hphi, matrixElements)
!!$!
!!$! Purpose:
!!$! ========
!!$!
!!$! Calling arguments:
!!$! ==================
!!$!
!!$use module_base
!!$use module_types
!!$use module_interfaces, exceptThisOne => getMatrixElements
!!$implicit none
!!$
!!$! Calling arguments
!!$integer,intent(in):: iproc, nproc
!!$type(locreg_descriptors),intent(in):: Glr
!!$type(orbitals_data),intent(in):: orbs
!!$type(communications_arrays),intent(in):: comms
!!$real(8),dimension(orbs%npsidim),intent(inout):: phi, hphi
!!$real(8),dimension(orbs%norb,orbs%norb,2),intent(out):: matrixElements
!!$
!!$! Local variables
!!$integer:: istart, jstart, nvctrp, iorb, jorb, istat, iall, ierr
!!$real(8):: ddot
!!$real(8),dimension(:),pointer:: phiWork
!!$character(len=*),parameter:: subname='getMatrixELements'
!!$real(8),dimension(:),allocatable:: phi2, hphi2
!!$
!!$
!!$
!!$  allocate(phi2(orbs%npsidim), stat=istat)
!!$  call memocc(istat, phi2, 'phi2', subname)
!!$  allocate(hphi2(orbs%npsidim), stat=istat)
!!$  call memocc(istat, hphi2, 'hphi2', subname)
!!$  allocate(phiWork(orbs%npsidim), stat=istat)
!!$  call memocc(istat, phiWork, 'phiWork', subname)
!!$  phiWork=0.d0
!!$  call vcopy(orbs%npsidim, phi(1), 1, phi2(1), 1)
!!$  call vcopy(orbs%npsidim, hphi(1), 1, hphi2(1), 1)
!!$
!!$  call transpose_v(iproc, nproc, orbs, Glr%wfd, comms, phi2, work=phiWork)
!!$  call transpose_v(iproc, nproc, orbs, Glr%wfd, comms, hphi2, work=phiWork)
!!$
!!$  matrixElements=0.d0
!!$
!!$  ! Calculate <phi_i|H_j|phi_j>
!!$  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!$  jstart=1
!!$  do jorb=1,orbs%norb
!!$      istart=1
!!$      do iorb=1,orbs%norb
!!$          matrixElements(iorb,jorb,2)=ddot(nvctrp, phi2(istart), 1, hphi2(jstart), 1)
!!$          !write(*,'(a,3i7,3es16.7)') 'iorb, jorb, iproc, ddot', iorb, jorb, iproc,  matrixElements(iorb,jorb,2), phi2(istart), hphi2(jstart)
!!$          istart=istart+nvctrp
!!$      end do
!!$      jstart=jstart+nvctrp
!!$  end do
!!$  call mpi_allreduce(matrixElements(1,1,2), matrixElements(1,1,1), orbs%norb**2, &
!!$      mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!$
!!$  call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, phi2, work=phiWork)
!!$  call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, hphi2, work=phiWork)
!!$
!!$
!!$
!!$  iall=-product(shape(phiWork))*kind(phiWork)
!!$  deallocate(phiWork, stat=istat)
!!$  call memocc(istat, iall, 'phiWork', subname)
!!$
!!$  iall=-product(shape(phi2))*kind(phi2)
!!$  deallocate(phi2, stat=istat)
!!$  call memocc(istat, iall, 'phi2', subname)
!!$
!!$  iall=-product(shape(hphi2))*kind(hphi2)
!!$  deallocate(hphi2, stat=istat)
!!$  call memocc(istat, iall, 'hphi2', subname)
!!$
!!$
!!$end subroutine getMatrixElements





subroutine modifiedBSEnergy(nspin, orbs, lin, HamSmall, matrixElements, ebsMod)
!
! Purpose:
! ========
!
! Calling arguments:
! ==================
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: nspin
type(orbitals_data),intent(in) :: orbs
type(linearParameters),intent(in):: lin
real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(in):: HamSmall, matrixElements
real(8),intent(out):: ebsMod

! Local variables
integer:: iorb, jorb, korb
real(8):: tt

  ! Calculate the modified band structure energy
  tt=0.d0
  do iorb=1,orbs%norb
      do jorb=1,lin%orbs%norb
          do korb=1,lin%orbs%norb
              tt=tt+HamSmall(korb,iorb)*HamSmall(jorb,iorb)*matrixElements(korb,jorb)
          end do
      end do
  end do
  if(nspin==1) then
      ebsMod=2.d0*tt ! 2 for closed shell
  else
      ebsMod=tt
  end if



end subroutine modifiedBSEnergy





subroutine modifiedBSEnergyModified(nspin, orbs, lin, coeff, matrixElements, ebsMod)
!
! Purpose:
! ========
!
! Calling arguments:
! ==================
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: nspin
type(orbitals_data),intent(in) :: orbs
type(linearParameters),intent(in):: lin
real(8),dimension(lin%orbs%norb,orbs%norb),intent(in):: coeff
real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(in):: matrixElements
real(8),intent(out):: ebsMod

! Local variables
integer:: iorb, jorb, korb
real(8):: tt

  ! Calculate the modified band structure energy
  tt=0.d0
  do iorb=1,orbs%norb
      do jorb=1,lin%orbs%norb
          do korb=1,lin%orbs%norb
              tt=tt+coeff(korb,iorb)*coeff(jorb,iorb)*matrixElements(korb,jorb)
          end do
      end do
  end do
  if(nspin==1) then
      ebsMod=2.d0*tt ! 2 for closed shell
  else
      ebsMod=tt
  end if



end subroutine modifiedBSEnergyModified












!!!subroutine optimizeCoefficients(iproc, orbs, lin, nspin, matrixElements, coeff, infoCoeff)
!!!!
!!!! Purpose:
!!!! ========
!!!!   Determines the optimal coefficients which minimize the modified band structure energy, i.e.
!!!!   E = sum_{i}sum_{k,l}c_{ik}c_{il}<phi_k|H_l|phi_l>.
!!!!   This is done by a steepest descen minimization using the gradient of the above expression with
!!!!   respect to the coefficients c_{ik}.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   ----------------
!!!!     iproc            process ID
!!!!     orbs             type describing the physical orbitals psi
!!!!     lin              type containing parameters for the linear version
!!!!     nspin            nspin==1 -> closed shell, npsin==2 -> open shell
!!!!     matrixElements   contains the matrix elements <phi_k|H_l|phi_l>
!!!!   Output arguments:
!!!!   -----------------
!!!!     coeff            the optimized coefficients 
!!!!     infoCoeff        if infoCoeff=0, the optimization converged
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nspin
!!!type(orbitals_data),intent(in):: orbs
!!!type(linearParameters),intent(in):: lin
!!!real(8),dimension(lin%lb%orbs%norb,lin%lb%orbs%norb),intent(in):: matrixElements
!!!real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(inout):: coeff
!!!integer,intent(out):: infoCoeff
!!!
!!!! Local variables
!!!integer:: it, iorb, jorb, k, l, istat, iall, korb, ierr
!!!real(8):: tt, fnrm, ddot, dnrm2, meanAlpha, cosangle, ebsMod, ebsModOld
!!!real(8),dimension(:,:),allocatable:: grad, gradOld, lagMat
!!!real(8),dimension(:),allocatable:: alpha
!!!character(len=*),parameter:: subname='optimizeCoefficients'
!!!logical:: converged
!!!
!!!
!!!! Allocate all local arrays.
!!!allocate(grad(lin%lb%orbs%norb,orbs%norb), stat=istat)
!!!call memocc(istat, grad, 'grad', subname)
!!!allocate(gradOld(lin%lb%orbs%norb,orbs%norb), stat=istat)
!!!call memocc(istat, gradOld, 'gradOld', subname)
!!!allocate(lagMat(orbs%norb,orbs%norb), stat=istat)
!!!call memocc(istat, lagMat, 'lagMat', subname)
!!!allocate(alpha(orbs%norb), stat=istat)
!!!call memocc(istat, alpha, 'alpha', subname)
!!!
!!!! trace of matrixElements
!!!if(iproc==0) then
!!!    tt=0.d0
!!!    do iorb=1,lin%lb%orbs%norb
!!!        do jorb=1,lin%lb%orbs%norb
!!!            if(iorb==jorb) tt=tt+matrixElements(iorb,jorb)
!!!            !write(777,*) iorb,jorb,matrixElements(jorb,iorb)
!!!        end do
!!!    end do
!!!    !write(*,*) 'trace',tt
!!!end if
!!!
!!!! Do everything only on the root process and then broadcast to all processes.
!!!! Maybe this part can be parallelized later, but at the moment it is not necessary since
!!!! it is fast enough.
!!!processIf: if(iproc==0) then
!!!    
!!!
!!!    !!! Orthonormalize the coefficient vectors (Gram-Schmidt).
!!!    !!do iorb=1,orbs%norb
!!!    !!    do jorb=1,iorb-1
!!!    !!        tt=ddot(lin%lb%orbs%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
!!!    !!        call daxpy(lin%lb%orbs%norb, -tt, coeff(1,jorb), 1, coeff(1,iorb), 1)
!!!    !!    end do
!!!    !!    tt=dnrm2(lin%lb%orbs%norb, coeff(1,iorb), 1)
!!!    !!    call dscal(lin%lb%orbs%norb, 1/tt, coeff(1,iorb), 1)
!!!    !!end do
!!!    
!!!    ! Initial step size for the optimization
!!!    alpha=5.d-3
!!!
!!!    ! Flag which checks convergence.
!!!    converged=.false.
!!!
!!!    if(iproc==0) write(*,'(1x,a)') '============================== optmizing coefficients =============================='
!!!
!!!    ! The optimization loop.
!!!    iterLoop: do it=1,lin%nItCoeff
!!!
!!!        if (iproc==0) then
!!!            write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
!!!        endif
!!!
!!!
!!!        ! Orthonormalize the coefficient vectors (Gram-Schmidt).
!!!        do iorb=1,orbs%norb
!!!            do jorb=1,iorb-1
!!!                tt=ddot(lin%lb%orbs%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
!!!                call daxpy(lin%lb%orbs%norb, -tt, coeff(1,jorb), 1, coeff(1,iorb), 1)
!!!            end do
!!!            tt=dnrm2(lin%lb%orbs%norb, coeff(1,iorb), 1)
!!!            call dscal(lin%lb%orbs%norb, 1/tt, coeff(1,iorb), 1)
!!!        end do
!!!
!!!
!!!        ! Calculate the gradient grad. At the same time we determine whether the step size shall be increased
!!!        ! or decreased (depending on gradient feedback).
!!!        meanAlpha=0.d0
!!!        grad=0.d0
!!!        do iorb=1,orbs%norb
!!!            do l=1,lin%lb%orbs%norb
!!!                do k=1,lin%lb%orbs%norb
!!!                    grad(l,iorb)=grad(l,iorb)+coeff(k,iorb)*(matrixElements(k,l)+matrixElements(l,k))
!!!                end do
!!!            end do
!!!            if(it>1) then
!!!                cosangle=ddot(lin%lb%orbs%norb, grad(1,iorb), 1, gradOld(1,iorb), 1)
!!!                cosangle=cosangle/dnrm2(lin%lb%orbs%norb, grad(1,iorb), 1)
!!!                cosangle=cosangle/dnrm2(lin%lb%orbs%norb, gradOld(1,iorb), 1)
!!!                !write(*,*) 'cosangle, ebsMod, ebsModOld', cosangle, ebsMod, ebsModOld
!!!                if(cosangle>.8d0 .and. ebsMod<ebsModOld+1.d-6*abs(ebsModOld)) then
!!!                    alpha(iorb)=max(alpha(iorb)*1.05d0,1.d-6)
!!!                else
!!!                    alpha(iorb)=max(alpha(iorb)*.5d0,1.d-6)
!!!                end if
!!!            end if
!!!            call dcopy(lin%lb%orbs%norb, grad(1,iorb), 1, gradOld(1,iorb), 1)
!!!            meanAlpha=meanAlpha+alpha(iorb)
!!!        end do
!!!        meanAlpha=meanAlpha/orbs%norb
!!!    
!!!    
!!!        ! Apply the orthoconstraint to the gradient. To do so first calculate the Lagrange
!!!        ! multiplier matrix.
!!!        lagMat=0.d0
!!!        do iorb=1,orbs%norb
!!!            do jorb=1,orbs%norb
!!!                do k=1,lin%lb%orbs%norb
!!!                    lagMat(iorb,jorb)=lagMat(iorb,jorb)+coeff(k,iorb)*grad(k,jorb)
!!!                end do
!!!            end do
!!!        end do
!!!
!!!        ! Now apply the orthoconstraint.
!!!        do iorb=1,orbs%norb
!!!            do k=1,lin%lb%orbs%norb
!!!                do jorb=1,orbs%norb
!!!                    grad(k,iorb)=grad(k,iorb)-.5d0*(lagMat(iorb,jorb)*coeff(k,jorb)+lagMat(jorb,iorb)*coeff(k,jorb))
!!!                end do
!!!            end do
!!!        end do
!!!    
!!!        
!!!        ! Calculate the modified band structure energy and the gradient norm.
!!!        if(it>1) then
!!!            ebsModOld=ebsMod
!!!        else
!!!            ebsModOld=1.d10
!!!        end if
!!!        ebsMod=0.d0
!!!        fnrm=0.d0
!!!        do iorb=1,orbs%norb
!!!            fnrm=fnrm+dnrm2(lin%lb%orbs%norb, grad(1,iorb), 1)
!!!            do jorb=1,lin%lb%orbs%norb
!!!                do korb=1,lin%lb%orbs%norb
!!!                    ebsMod=ebsMod+coeff(korb,iorb)*coeff(jorb,iorb)*matrixElements(korb,jorb)
!!!                end do
!!!            end do
!!!        end do
!!!    
!!!        ! Multiply the energy with a factor of 2 if we have a closed-shell system.
!!!        if(nspin==1) then
!!!            ebsMod=2.d0*ebsMod
!!!        end if
!!!
!!!        !if(iproc==0) write(*,'(1x,a,4x,i0,es12.4,3x,es10.3, es19.9)') 'iter, fnrm, meanAlpha, Energy', &
!!!        if(iproc==0) write(*,'(1x,a,es11.2,es22.13,es10.2)') 'fnrm, band structure energy, mean alpha', &
!!!            fnrm, ebsMod, meanAlpha
!!!        
!!!        ! Check for convergence.
!!!        if(fnrm<lin%convCritCoeff) then
!!!            if(iproc==0) write(*,'(1x,a,i0,a)') 'converged in ', it, ' iterations.'
!!!            if(iproc==0) write(*,'(3x,a,2es14.5)') 'Final values for fnrm, Energy:', fnrm, ebsMod
!!!            converged=.true.
!!!            infoCoeff=it
!!!            exit
!!!        end if
!!!  
!!!        if(it==lin%nItCoeff) then
!!!            if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
!!!                ' iterations! Exiting loop due to limitations of iterations.'
!!!            if(iproc==0) write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, Energy: ', fnrm, ebsMod
!!!            infoCoeff=-1
!!!            exit
!!!        end if
!!!
!!!        ! Improve the coefficients (by steepet descent).
!!!        do iorb=1,orbs%norb
!!!            do l=1,lin%lb%orbs%norb
!!!                coeff(l,iorb)=coeff(l,iorb)-alpha(iorb)*grad(l,iorb)
!!!            end do
!!!        end do
!!!    
!!!
!!!    end do iterLoop
!!!
!!!    !!if(.not.converged) then
!!!    !!    if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
!!!    !!        ' iterations! Exiting loop due to limitations of iterations.'
!!!    !!    if(iproc==0) write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, Energy: ', fnrm, ebsMod
!!!    !!    infoCoeff=-1
!!!    !!    ! Orthonormalize the coefficient vectors (Gram-Schmidt).
!!!    !!    do iorb=1,orbs%norb
!!!    !!        do jorb=1,iorb-1
!!!    !!            tt=ddot(lin%lb%orbs%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
!!!    !!            call daxpy(lin%lb%orbs%norb, -tt, coeff(1,jorb), 1, coeff(1,iorb), 1)
!!!    !!        end do
!!!    !!        tt=dnrm2(lin%lb%orbs%norb, coeff(1,iorb), 1)
!!!    !!        call dscal(lin%lb%orbs%norb, 1/tt, coeff(1,iorb), 1)
!!!    !!    end do
!!!    !!end if
!!!
!!!    if(iproc==0) write(*,'(1x,a)') '===================================================================================='
!!!
!!!end if processIf
!!!
!!!
!!!! Now broadcast the result to all processes
!!!call mpi_bcast(coeff(1,1), lin%lb%orbs%norb*orbs%norb, mpi_double_precision, 0, mpi_comm_world, ierr)
!!!call mpi_bcast(infoCoeff, 1, mpi_integer, 0, mpi_comm_world, ierr)
!!!
!!!
!!!! Deallocate all local arrays.
!!!iall=-product(shape(grad))*kind(grad)
!!!deallocate(grad, stat=istat)
!!!call memocc(istat, iall, 'grad', subname)
!!!
!!!iall=-product(shape(gradOld))*kind(gradOld)
!!!deallocate(gradOld, stat=istat)
!!!call memocc(istat, iall, 'gradOld', subname)
!!!
!!!iall=-product(shape(lagMat))*kind(lagMat)
!!!deallocate(lagMat, stat=istat)
!!!call memocc(istat, iall, 'lagMat', subname)
!!!
!!!iall=-product(shape(alpha))*kind(alpha)
!!!deallocate(alpha, stat=istat)
!!!call memocc(istat, iall, 'alpha', subname)
!!!
!!!end subroutine optimizeCoefficients







subroutine postCommunicationSumrho2(iproc, nproc, comsr, sendBuf, recvBuf)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
!type(p2pCommsSumrho),intent(inout):: comsr
type(p2pComms),intent(inout):: comsr
real(8),dimension(comsr%nsendBuf),intent(inout):: sendBuf
real(8),dimension(comsr%nrecvBuf),intent(out):: recvBuf

! Local variables
integer:: jproc, nreceives, nsends, iorb, mpisource, istsource, ncount, lrsource, mpidest, istdest, tag, ierr
integer:: ist, istr, ilr


! Communicate the orbitals for the calculation of the charge density.
! Since we use non blocking point to point communication, only post the message
! and continues with other calculations.
! Be aware that you must not modify the send buffer without checking whether
! the communications has completed.
if(iproc==0) write(*,'(1x,a)', advance='no') 'Posting sends / receives for the calculation of the charge density... '
nreceives=0
nsends=0
comsr%communComplete=.false.
procLoop1: do jproc=0,nproc-1
    orbsLoop1: do iorb=1,comsr%noverlaps(jproc)
        mpisource=comsr%comarr(1,iorb,jproc)
        istsource=comsr%comarr(2,iorb,jproc)
        ncount=comsr%comarr(3,iorb,jproc)
        lrsource=comsr%comarr(4,iorb,jproc)
        mpidest=comsr%comarr(5,iorb,jproc)
        istdest=comsr%comarr(6,iorb,jproc)
        tag=comsr%comarr(7,iorb,jproc)
        if(ncount==0) then
            ! No communication is needed. This should be improved in the initialization, i.e. this communication
            ! with 0 elements should be removed from comgp%noverlaps etc.
            comsr%comarr(8,iorb,jproc)=mpi_request_null
            comsr%comarr(9,iorb,jproc)=mpi_request_null
            comsr%communComplete(iorb,jproc)=.true.
            if(iproc==mpidest) then
                ! This is just to make the check at the end happy.
                nreceives=nreceives+1
            end if
        else
            if(mpisource/=mpidest) then
                ! The orbitals are on different processes, so we need a point to point communication.
                if(iproc==mpisource) then
                    !write(*,'(6(a,i0))') 'sumrho: process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
                    !call mpi_isend(lphi(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, comsr%comarr(8,iorb,jproc), ierr)
                    call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
                         comsr%comarr(8,iorb,jproc), ierr)
                    comsr%comarr(9,iorb,jproc)=mpi_request_null !is this correct?
                    nsends=nsends+1
                else if(iproc==mpidest) then
                   !write(*,'(6(a,i0))') 'sumrho: process ', mpidest, ' receives ', ncount, &
                   !     ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
                    call mpi_irecv(recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
                         comsr%comarr(9,iorb,jproc), ierr)
                    comsr%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
                    nreceives=nreceives+1
                else
                    comsr%comarr(8,iorb,jproc)=mpi_request_null
                    comsr%comarr(9,iorb,jproc)=mpi_request_null
                end if
            else
                ! The orbitals are on the same process, so simply copy them.
                if(iproc==mpisource) then
                    !write(*,'(6(a,i0))') 'sumrho: process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
                    call dcopy(ncount, sendBuf(istsource), 1, recvBuf(istdest), 1)
                    comsr%comarr(8,iorb,jproc)=mpi_request_null
                    comsr%comarr(9,iorb,jproc)=mpi_request_null
                    nsends=nsends+1
                    nreceives=nreceives+1
                    comsr%communComplete(iorb,mpisource)=.true.
                else
                    comsr%comarr(8,iorb,jproc)=mpi_request_null
                    comsr%comarr(9,iorb,jproc)=mpi_request_null
                    comsr%communComplete(iorb,mpisource)=.true.
                end if
            end if
        end if
    end do orbsLoop1
end do procLoop1
if(iproc==0) write(*,'(a)') 'done.'

if(nreceives/=comsr%noverlaps(iproc)) then
    write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comsr%noverlaps(iproc)', nreceives,&
         comsr%noverlaps(iproc)
    stop
end if
call mpi_barrier(mpi_comm_world, ierr)

end subroutine postCommunicationSumrho2





!subroutine initializeCommunicationPotential(iproc, nproc, nscatterarr, lin)
subroutine initializeCommunicationPotential(iproc, nproc, nscatterarr, orbs, lzd, comgp, onWhichAtomAll, tag)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!type(linearParameters),intent(inout):: lin
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
!type(p2pCommsGatherPot),intent(out):: comgp
type(p2pComms),intent(out):: comgp
integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
integer,intent(inout):: tag

! Local variables
integer:: is1, ie1, is2, ie2, is3, ie3, ilr, ii, iorb, iiorb, jproc, kproc, istat, iall
integer:: ioverlap, is3j, ie3j, is3k, ie3k, mpidest, istdest, ioffset, is3min, ie3max
integer,dimension(:,:),allocatable:: iStartEnd
character(len=*),parameter:: subname='setCommunicationPotential'


! Determine the bounds of the potential that we need for
! the orbitals on this process.
allocate(iStartEnd(6,0:nproc-1), stat=istat)
call memocc(istat, iStartEnd, 'iStartEnd', subname)
is1=0
ie1=0
is2=0
ie2=0
is3=0
ie3=0
iiorb=0
do jproc=0,nproc-1
    do iorb=1,orbs%norb_par(jproc,0)
        
        iiorb=iiorb+1 
        ilr=onWhichAtomAll(iiorb)
    
        ii=lzd%Llr(ilr)%nsi1
        if(ii < is1 .or. iorb==1) then
            is1=ii
        end if
        ii=lzd%Llr(ilr)%nsi1+lzd%Llr(ilr)%d%n1i
        if(ii > ie1 .or. iorb==1) then
            ie1=ii
        end if
    
        ii=lzd%Llr(ilr)%nsi2
        if(ii < is2 .or. iorb==1) then
            is2=ii
        end if
        ii=lzd%Llr(ilr)%nsi2+lzd%Llr(ilr)%d%n2i
        if(ii > ie2 .or. iorb==1) then
            ie2=ii
        end if
    
        ii=lzd%Llr(ilr)%nsi3
        if(ii < is3 .or. iorb==1) then
            is3=ii
        end if
        ii=lzd%Llr(ilr)%nsi3+lzd%Llr(ilr)%d%n3i
        if(ii > ie3 .or. iorb==1) then
            ie3=ii
        end if
    
    end do
    iStartEnd(1,jproc)=is1
    iStartEnd(2,jproc)=ie1
    iStartEnd(3,jproc)=is2
    iStartEnd(4,jproc)=ie2
    iStartEnd(5,jproc)=is3
    iStartEnd(6,jproc)=ie3
end do

! Determine how many slices each process receives.
allocate(comgp%noverlaps(0:nproc-1), stat=istat)
call memocc(istat, comgp%noverlaps, 'comgp%noverlaps', subname)
do jproc=0,nproc-1
    is3j=istartEnd(5,jproc)
    ie3j=istartEnd(6,jproc)
    mpidest=jproc
    ioverlap=0
    do kproc=0,nproc-1
        is3k=nscatterarr(kproc,3)+1
        ie3k=is3k+nscatterarr(kproc,2)-1
        if(is3j<=ie3k .and. ie3j>=is3k) then
            ioverlap=ioverlap+1
            !if(iproc==0) write(*,'(2(a,i0),a)') 'process ',jproc,' gets potential from process ',kproc,'.' 
        end if
    end do
    comgp%noverlaps(jproc)=ioverlap
    !if(iproc==0) write(*,'(2(a,i0),a)') 'Process ',jproc,' gets ',ioverlap,' potential slices.'
end do

! Determine the parameters for the communications.
allocate(comgp%overlaps(comgp%noverlaps(iproc)), stat=istat)
call memocc(istat, comgp%overlaps, 'comgp%overlaps', subname)
allocate(comgp%comarr(8,maxval(comgp%noverlaps),0:nproc-1))
call memocc(istat, comgp%comarr, 'comgp%comarr', subname)
allocate(comgp%ise3(2,0:nproc-1), stat=istat)
call memocc(istat, comgp%ise3, 'comgp%ise3', subname)
comgp%nrecvBuf = 0
is3min=0
ie3max=0
do jproc=0,nproc-1
    is3j=istartEnd(5,jproc)
    ie3j=istartEnd(6,jproc)
    mpidest=jproc
    ioverlap=0
    istdest=1
    do kproc=0,nproc-1
        is3k=nscatterarr(kproc,3)+1
        ie3k=is3k+nscatterarr(kproc,2)-1
        if(is3j<=ie3k .and. ie3j>=is3k) then
            is3=max(is3j,is3k) ! starting index in z dimension for data to be sent
            ie3=min(ie3j,ie3k) ! ending index in z dimension for data to be sent
            ioffset=is3-is3k ! starting index (in z direction) of data to be sent (actually it is the index -1)
            ioverlap=ioverlap+1
            tag=tag+1
            if(is3<is3min .or. ioverlap==1) then
                is3min=is3
            end if
            if(ie3>ie3max .or. ioverlap==1) then
                ie3max=ie3
            end if
            !write(*,'(a,8i8)') 'jproc, kproc, is3j, ie3j, is3k, ie3k, is3, ie3', jproc, kproc, is3j, ie3j, is3k, ie3k, is3, ie3
            call setCommunicationPotential(kproc, is3, ie3, ioffset, lzd%Glr%d%n1i, lzd%Glr%d%n2i, jproc,&
                 istdest, tag, comgp%comarr(1,ioverlap,jproc))
            !if(iproc==0) write(*,'(6(a,i0))') 'process ',comgp%comarr(1,ioverlap,jproc),' sends ',comgp%comarr(3,ioverlap,jproc),' elements from position ',&
            !                        comgp%comarr(2,ioverlap,jproc),' to position ',comgp%comarr(5,ioverlap,jproc),' on process ',&
            !                        comgp%comarr(4,ioverlap,jproc),'; tag=',comgp%comarr(6,ioverlap,jproc)
            istdest = istdest + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
            !write(*,'(a,4i8)') 'jproc, kproc, (ie3-is3+1),lzd%Glr%d%n1i*lzd%Glr%d%n2i', jproc, kproc, (ie3-is3+1),lzd%Glr%d%n1i*lzd%Glr%d%n2i
            if(iproc==jproc) then
                comgp%nrecvBuf = comgp%nrecvBuf + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
            end if
        end if
    end do
    comgp%ise3(1,jproc)=is3min
    comgp%ise3(2,jproc)=ie3max
    !if(iproc==0) write(*,'(a,3i8)') 'jproc, comgp%ise3(1,jproc), comgp%ise3(2,jproc)', jproc, comgp%ise3(1,jproc), comgp%ise3(2,jproc)
end do

!write(*,'(a,i4,i12)') 'iproc, comgp%nrecvBuf', iproc, comgp%nrecvBuf

allocate(comgp%communComplete(maxval(comgp%noverlaps),0:nproc-1), stat=istat)
call memocc(istat, comgp%communComplete, 'comgp%communComplete', subname)

iall=-product(shape(iStartEnd))*kind(iStartEnd)
deallocate(iStartEnd, stat=istat)
call memocc(istat, iall, 'iStartEnd', subname)

end subroutine initializeCommunicationPotential


subroutine setCommunicationPotential(mpisource, is3, ie3, ioffset, n1i, n2i, mpidest, istdest, tag, comarr)
use module_base
implicit none
! Calling arguments
integer,intent(in):: mpisource, is3, ie3, ioffset, n1i, n2i, mpidest, istdest, tag
integer,dimension(8),intent(out):: comarr

! Local variables
integer:: istsource, ncount

! From which MPI process shall the slice be sent
comarr(1)=mpisource

! Starting index on the sending process
istsource=ioffset*n1i*n2i+1
comarr(2)=istsource

! Amount of data to be sent
ncount=(ie3-is3+1)*n1i*n2i
comarr(3)=ncount

! To which MPI process shall the slice be sent
comarr(4)=mpidest

! Starting index on the receiving index
comarr(5)=istdest

! Tag for the communication
comarr(6)=tag

! comarr(7): this entry is used as request for the mpi_isend.

! comarr(8): this entry is used as request for the mpi_irecv.


end subroutine setCommunicationPotential




subroutine postCommunicationsPotential(iproc, nproc, ndimpot, pot, comgp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, ndimpot
real(8),dimension(ndimpot),intent(in):: pot
!type(p2pCommsGatherPot),intent(inout):: comgp
type(p2pComms),intent(inout):: comgp

! Local variables
integer:: jproc, kproc, nsends, nreceives, istat, mpisource, istsource, ncount, mpidest, istdest, tag, ierr


! Post the messages
if(iproc==0) write(*,'(1x,a)', advance='no') 'Posting sends / receives for communicating the potential... '
nreceives=0
nsends=0
comgp%communComplete=.false.
destLoop: do jproc=0,nproc-1
    sourceLoop: do kproc=1,comgp%noverlaps(jproc)
        mpisource=comgp%comarr(1,kproc,jproc)
        istsource=comgp%comarr(2,kproc,jproc)
        ncount=comgp%comarr(3,kproc,jproc)
        mpidest=comgp%comarr(4,kproc,jproc)
        istdest=comgp%comarr(5,kproc,jproc)
        tag=comgp%comarr(6,kproc,jproc)
        if(ncount==0) then
            ! No communication is needed. This should be improved in the initialization, i.e. this communication
            ! with 0 elements should be removed from comgp%noverlaps etc.
            comgp%comarr(7,kproc,jproc)=mpi_request_null
            comgp%comarr(8,kproc,jproc)=mpi_request_null
            comgp%communComplete(kproc,jproc)=.true.
            if(iproc==mpidest) then
                ! This is just to make the check at the end happy.
                nreceives=nreceives+1
            end if
        else
            if(mpisource/=mpidest) then
                if(iproc==mpisource) then
                    !write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
                    call mpi_isend(pot(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
                         comgp%comarr(7,kproc,jproc), ierr)
                    comgp%comarr(8,kproc,jproc)=mpi_request_null !is this correct?
                    nsends=nsends+1
                else if(iproc==mpidest) then
                   !write(*,'(6(a,i0))') 'process ', mpidest, ' receives ', ncount, &
                   !    ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
                    call mpi_irecv(comgp%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
                        comgp%comarr(8,kproc,jproc), ierr)
                    comgp%comarr(7,kproc,jproc)=mpi_request_null !is this correct?
                    nreceives=nreceives+1
                else
                    comgp%comarr(7,kproc,jproc)=mpi_request_null
                    comgp%comarr(8,kproc,jproc)=mpi_request_null
                end if
            else
                ! The orbitals are on the same process, so simply copy them.
                if(iproc==mpisource) then
                    !write(*,'(6(a,i0))') 'process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
                    !if(iproc==0) write(*,'(a,es22.12)') 'pot(istsource)', pot(istsource)
                    call dcopy(ncount, pot(istsource), 1, comgp%recvBuf(istdest), 1)
                    comgp%comarr(7,kproc,jproc)=mpi_request_null
                    comgp%comarr(8,kproc,jproc)=mpi_request_null
                    nsends=nsends+1
                    nreceives=nreceives+1
                    comgp%communComplete(kproc,iproc)=.true.
                else
                    comgp%comarr(7,kproc,jproc)=mpi_request_null
                    comgp%comarr(8,kproc,jproc)=mpi_request_null
                    comgp%communComplete(kproc,jproc)=.true.
                end if
            end if
        end if
    end do sourceLoop
end do destLoop
if(iproc==0) write(*,'(a)') 'done.'

if(nreceives/=comgp%noverlaps(iproc)) then
    write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comgp%noverlaps(iproc)',&
         nreceives, comgp%noverlaps(iproc)
    stop
end if


end subroutine postCommunicationsPotential





subroutine gatherPotential(iproc, nproc, comgp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
!type(p2pCommsGatherPot),intent(inout):: comgp
type(p2pComms),intent(inout):: comgp

! Local variables
integer:: kproc, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete

if(iproc==0) write(*,'(1x,a)',advance='no') 'Gathering the potential... '
! Check whether the communications have completed.
nfast=0
nsameproc=0
testLoop: do 
    do jproc=0,nproc-1
        do kproc=1,comgp%noverlaps(jproc)
           if(comgp%communComplete(kproc,jproc)) cycle
            call mpi_test(comgp%comarr(7,kproc,jproc), sendComplete, stat, ierr)      
            call mpi_test(comgp%comarr(8,kproc,jproc), receiveComplete, stat, ierr)   
            ! Attention: mpi_test is a local function.
            if(sendComplete .and. receiveComplete) comgp%communComplete(kproc,jproc)=.true.
            !!if(comgp%communComplete(kproc,jproc)) then
            !!    !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', korb
            !!    mpisource=comgp%comarr(1,kproc,jproc)
            !!    mpidest=comgp%comarr(4,kproc,jproc)
            !!    if(mpisource/=mpidest) then
            !!        nfast=nfast+1
            !!    else
            !!        nsameproc=nsameproc+1
            !!    end if
            !!end if
        end do
    end do
    ! If we made it until here, either all all the communication is
    ! complete or we better wait for each single orbital.
    exit testLoop
end do testLoop


! Since mpi_test is a local function, check whether the communication has completed on all processes.
call mpiallred(comgp%communComplete(1,0), nproc*maxval(comgp%noverlaps), mpi_land, mpi_comm_world, ierr)

! Wait for the communications that have not completed yet
nslow=0
do jproc=0,nproc-1
    do kproc=1,comgp%noverlaps(jproc)
        if(comgp%communComplete(kproc,jproc)) then
            mpisource=comgp%comarr(1,kproc,jproc)
            mpidest=comgp%comarr(4,kproc,jproc)
            if(mpisource==mpidest) then
                nsameproc=nsameproc+1
            else
                nfast=nfast+1
            end if
            cycle
        end if
        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
        nslow=nslow+1
        call mpi_wait(comgp%comarr(7,kproc,jproc), stat, ierr)  
        call mpi_wait(comgp%comarr(8,kproc,jproc), stat, ierr) 
        comgp%communComplete(kproc,jproc)=.true.
    end do
end do

call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
if (verbose > 3) then
   if(iproc==0) write(*,'(a,f5.1,a)') 'done. Communication overlap ratio:',100.d0*dble(nfast)/(dble(nfast+nslow)),'%'
else
   if(iproc==0) write(*,'(a,f5.1,a)') 'done.'
end if
!if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!                       nfast, ' could be overlapped with computation.'
!if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


end subroutine gatherPotential




!!!subroutine diagonalizeHamiltonianParallel(iproc, nproc, norb, ham, ovrlp, eval)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc, norb
!!!real(8),dimension(norb,norb),intent(inout):: ham, ovrlp
!!!real(8),dimension(norb),intent(out):: eval
!!!
!!!! Local variables
!!!integer:: ierr, mbrow, mbcol, i, j, istat, lwork, info, ii1, ii2, nproc_scalapack, iall
!!!integer:: nprocrow, nproccol, context, irow, icol, lnrow, lncol, numroc, jproc, liwork, neval_found, neval_computed
!!!real(8):: tt1, tt2
!!!real(8),dimension(:,:),allocatable:: lmat, loverlap, levec
!!!real(8),dimension(:),allocatable:: work, gap
!!!integer,dimension(9):: desc_levec, desc_lmat, desc_loverlap
!!!integer,dimension(:),allocatable:: iwork, ifail, icluster
!!!character(len=*),parameter:: subname='diagonalizeHamiltonianParallel'
!!!
!!!
!!!
!!!
!!!! Block size for scalapack
!!!mbrow=64
!!!mbcol=64
!!!
!!!! Number of processes that will be involved in the calculation
!!!tt1=dble(norb)/dble(mbrow)
!!!tt2=dble(norb)/dble(mbcol)
!!!ii1=ceiling(tt1)
!!!ii2=ceiling(tt2)
!!!nproc_scalapack = min(ii1*ii2,nproc)
!!!!nproc_scalapack = nproc
!!!if(iproc==0) write(*,'(a,i0,a)') 'scalapack will use ',nproc_scalapack,' processes.'
!!!
!!!! process grid: number of processes per row and column
!!!tt1=sqrt(dble(nproc_scalapack))
!!!ii1=ceiling(tt1)
!!!do i=ii1,nproc_scalapack
!!!    if(mod(nproc_scalapack,i)==0) then
!!!        nprocrow=i
!!!        exit
!!!    end if
!!!end do
!!!nproccol=nproc_scalapack/nprocrow
!!!if(iproc==0) write(*,'(a,i0,a,i0,a)') 'calculation is done on process grid with dimension ',nprocrow,' x ',nproccol,'.'
!!!
!!!
!!!! Initialize blacs context
!!!call blacs_get(-1, 0, context)
!!!call blacs_gridinit(context, 'r', nprocrow, nproccol )
!!!call blacs_gridinfo(context,nprocrow, nproccol, irow, icol)
!!!!write(*,*) 'iproc, irow, icol', iproc, irow, icol
!!!
!!!! Initialize the matrix mat to zero for processes that don't do the calculation.
!!!! For processes participating in the diagonalization, 
!!!! it will be partially (only at the position that process was working on) overwritten with the result. 
!!!! At the end we can the make an allreduce to get the correct result on all processes.
!!!if(irow==-1) ham=0.d0
!!!
!!!! Everything that follows is only done if the current process is part of the grid.
!!!processIf: if(irow/=-1) then
!!!    ! Determine the size of the matrix (lnrow x lncol):
!!!    lnrow = numroc(norb, mbrow, irow, 0, nprocrow)
!!!    lncol = numroc(norb, mbcol, icol, 0, nproccol)
!!!    write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local matrix of size ',lnrow,' x ',lncol
!!!
!!!    ! Initialize descriptor arrays.
!!!    call descinit(desc_lmat, norb, norb, mbrow, mbcol, 0, 0, context, lnrow, info)
!!!    call descinit(desc_loverlap, norb, norb, mbrow, mbcol, 0, 0, context, lnrow, info)
!!!    call descinit(desc_levec, norb, norb, mbrow, mbcol, 0, 0, context, lnrow, info)
!!!
!!!    ! Allocate the local array lmat
!!!    allocate(lmat(lnrow,lncol), stat=istat)
!!!    call memocc(istat, lmat, 'lmat', subname)
!!!    allocate(loverlap(lnrow,lncol), stat=istat)
!!!    call memocc(istat, loverlap, 'loverlap', subname)
!!!
!!!    ! Copy the global array mat to the local array lmat.
!!!    ! The same for loverlap and overlap, respectively.
!!!    !call dcopy(norb**2, ham(1,1), 1, mat(1,1), 1)
!!!    !call dcopy(norb**2, ovrlp(1,1), 1, overlap(1,1), 1)
!!!    do i=1,norb
!!!        do j=1,norb
!!!            call pdelset(lmat(1,1), j, i, desc_lmat, ham(j,i))
!!!            call pdelset(loverlap(1,1), j, i, desc_loverlap, ovrlp(j,i))
!!!        end do
!!!    end do
!!!
!!!
!!!    ! Solve the generalized eigenvalue problem.
!!!    allocate(levec(lnrow,lncol), stat=istat)
!!!    call memocc(istat, levec, 'levec', subname)
!!!    allocate(ifail(norb), stat=istat)
!!!    call memocc(istat, ifail, 'ifail', subname)
!!!    allocate(icluster(2*nprocrow*nproccol), stat=istat)
!!!    call memocc(istat, icluster, 'icluster', subname)
!!!    allocate(gap(nprocrow*nproccol), stat=istat)
!!!    call memocc(istat, gap, 'gap', subname)
!!!
!!!    ! workspace query
!!!    lwork=-1
!!!    liwork=-1
!!!    allocate(work(1), stat=istat)
!!!    call memocc(istat, work, 'work', subname)
!!!    allocate(iwork(1), stat=istat) ; if(istat/=0) stop 'ERROR in allocating'
!!!    call memocc(istat, iwork, 'iwork', subname)
!!!    call pdsygvx(1, 'v', 'a', 'l', norb, lmat(1,1), 1, 1, desc_lmat, loverlap(1,1), 1, 1, &
!!!                 desc_loverlap, 0.d0, 1.d0, 0, 1, -1.d0, neval_found, neval_computed, eval(1), &
!!!                 -1.d0, levec(1,1), 1, 1, desc_levec, work, lwork, iwork, liwork, &
!!!                 ifail, icluster, gap, info)
!!!    lwork=ceiling(work(1))
!!!    liwork=iwork(1)
!!!    !write(*,*) 'iproc, lwork, liwork', iproc, lwork, liwork
!!!    iall=-product(shape(work))*kind(work)
!!!    deallocate(work, stat=istat)
!!!    call memocc(istat, iall, 'work', subname)
!!!    iall=-product(shape(iwork))*kind(iwork)
!!!    deallocate(iwork, stat=istat)
!!!    call memocc(istat, iall, 'iwork', subname)
!!!
!!!    allocate(work(lwork), stat=istat)
!!!    call memocc(istat, work, 'work', subname)
!!!    allocate(iwork(liwork), stat=istat)
!!!    call memocc(istat, iwork, 'iwork', subname)
!!!
!!!    call pdsygvx(1, 'v', 'a', 'l', norb, lmat(1,1), 1, 1, desc_lmat, loverlap(1,1), 1, 1, &
!!!                 desc_loverlap, 0.d0, 1.d0, 0, 1, -1.d0, neval_found, neval_computed, eval(1), &
!!!                 -1.d0, levec(1,1), 1, 1, desc_levec, work, lwork, iwork, liwork, &
!!!                 ifail, icluster, gap, info)
!!!
!!!    ! Gather together the eigenvectors from all processes and store them in mat.
!!!    do i=1,norb
!!!        do j=1,norb
!!!            call pdelset2(ham(j,i), levec(1,1), j, i, desc_lmat, 0.d0)
!!!        end do
!!!    end do
!!!
!!!
!!!    iall=-product(shape(lmat))*kind(lmat)
!!!    deallocate(lmat, stat=istat)
!!!    call memocc(istat, iall, 'lmat', subname)
!!!
!!!    iall=-product(shape(levec))*kind(levec)
!!!    deallocate(levec, stat=istat)
!!!    call memocc(istat, iall, 'levec', subname)
!!!
!!!    iall=-product(shape(loverlap))*kind(loverlap)
!!!    deallocate(loverlap, stat=istat)
!!!    call memocc(istat, iall, 'loverlap', subname)
!!!
!!!    iall=-product(shape(work))*kind(work)
!!!    deallocate(work, stat=istat)
!!!    call memocc(istat, iall, 'work', subname)
!!!
!!!    iall=-product(shape(iwork))*kind(iwork)
!!!    deallocate(iwork, stat=istat)
!!!    call memocc(istat, iall, 'iwork', subname)
!!!
!!!    iall=-product(shape(ifail))*kind(ifail)
!!!    deallocate(ifail, stat=istat)
!!!    call memocc(istat, iall, 'ifail', subname)
!!!
!!!    iall=-product(shape(icluster))*kind(icluster)
!!!    deallocate(icluster, stat=istat)
!!!    call memocc(istat, iall, 'icluster', subname)
!!!
!!!    iall=-product(shape(gap))*kind(gap)
!!!    deallocate(gap, stat=istat)
!!!    call memocc(istat, iall, 'gap', subname)
!!!
!!!end if processIF
!!!
!!!! Gather the eigenvectors on all processes.
!!!call mpiallred(ham(1,1), norb**2, mpi_sum, mpi_comm_world, ierr)
!!!
!!!! Broadcast the eigenvalues if required. If nproc_scalapack==nproc, then all processes
!!!! diagonalized the matrix and therefore have the eigenvalues.
!!!if(nproc_scalapack/=nproc) then
!!!    call mpi_bcast(eval(1), norb, mpi_double_precision, 0, mpi_comm_world, ierr)
!!!end if
!!!
!!!
!!!
!!!
!!!end subroutine diagonalizeHamiltonianParallel





!!$subroutine prepare_lnlpspd(iproc, at, input, orbs, rxyz, radii_cf, locregShape, lzd)
!!$  use module_base
!!$  use module_types
!!$  use module_interfaces, exceptThisOne => prepare_lnlpspd
!!$  implicit none
!!$  
!!$  ! Calling arguments
!!$  integer,intent(in):: iproc
!!$  type(atoms_data),intent(in):: at
!!$  type(input_variables),intent(in):: input
!!$  type(orbitals_data),intent(in):: orbs
!!$  real(8),dimension(3,at%nat),intent(in):: rxyz
!!$  real(8),dimension(at%ntypes,3),intent(in):: radii_cf
!!$  character(len=1),intent(in):: locregShape
!!$  type(local_zone_descriptors),intent(inout):: lzd
!!$  
!!$  ! Local variables
!!$  integer:: ilr, istat, iorb
!!$  logical:: calc
!!$  character(len=*),parameter:: subname='prepare_lnlpspd'
!!$
!!$
!!$  allocate(Lzd%Lnlpspd(Lzd%nlr), stat=istat)
!!$  do ilr=1,Lzd%nlr
!!$      call nullify_nonlocal_psp_descriptors(Lzd%Lnlpspd(ilr))
!!$  end do
!!$
!!$  do ilr=1,Lzd%nlr
!!$
!!$      nullify(Lzd%Llr(ilr)%projflg) !to avoid problems when deallocating
!!$      calc=.false.
!!$      do iorb=1,orbs%norbp
!!$          if(ilr == orbs%inwhichLocreg(iorb+orbs%isorb)) calc=.true.
!!$      end do
!!$      if (.not. calc) cycle !calculate only for the locreg on this processor, without repeating for same locreg.
!!$      ! allocate projflg
!!$      allocate(Lzd%Llr(ilr)%projflg(at%nat),stat=istat)
!!$      call memocc(istat,Lzd%Llr(ilr)%projflg,'Lzd%Llr(ilr)%projflg',subname)
!!$
!!$      call nlpspd_to_locreg(input,iproc,Lzd%Glr,Lzd%Llr(ilr),rxyz,at,orbs,&
!!$           radii_cf,input%frmult,input%frmult,&
!!$           input%hx,input%hy,input%hz,locregShape,lzd%Gnlpspd,&
!!$           Lzd%Lnlpspd(ilr),Lzd%Llr(ilr)%projflg)
!!$  end do
!!$
!!$end subroutine prepare_lnlpspd


!!$subroutine free_lnlpspd(orbs, lzd)
!!$  use module_base
!!$  use module_types
!!$  !use deallocatePointers
!!$  use module_interfaces, exceptThisOne => free_lnlpspd
!!$  implicit none
!!$  
!!$  ! Calling arguments
!!$  type(orbitals_data),intent(in):: orbs
!!$  type(local_zone_descriptors),intent(inout):: lzd
!!$
!!$  ! Local variables
!!$  integer:: ilr, iorb, istat, iall
!!$  logical:: go
!!$  character(len=*),parameter:: subname='free_lnlpspd'
!!$
!!$  do ilr=1,lzd%nlr
!!$
!!$      go=.false.
!!$      do iorb=1,orbs%norbp
!!$         if(ilr == orbs%inwhichLocreg(iorb+orbs%isorb)) go=.true.
!!$      end do
!!$      if (.not. go) cycle !deallocate only for the locreg on this processor, without repeating for same locreg.
!!$
!!$      ! Deallocate projflg.
!!$      !call checkAndDeallocatePointer(lzd%llr(ilr)%projflg, 'lzd%llr(ilr)%projflg', subname)
!!$      iall=-product(shape(lzd%llr(ilr)%projflg))*kind(lzd%llr(ilr)%projflg)
!!$      deallocate(lzd%llr(ilr)%projflg, stat=istat)
!!$      call memocc(istat, iall, 'lzd%llr(ilr)%projflg', subname)
!!$
!!$      call deallocate_nonlocal_psp_descriptors(lzd%lnlpspd(ilr), subname)
!!$  end do
!!$
!!$!!$  deallocate(lzd%lnlpspd)
!!$!!$  nullify(lzd%lnlpspd)
!!$
!!$end subroutine free_lnlpspd




subroutine apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, confdatarr, hx, &
           psi, centralLocreg, vpsi)
use module_base
use module_types
use module_interfaces, except_this_one => apply_orbitaldependent_potential
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, centralLocreg
type(atoms_data),intent(in):: at
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(3,at%nat),intent(in):: rxyz
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),intent(in):: hx
!real(8),dimension(lzd%lpsidimtot),intent(in):: psi  !!!! ATENTION, intent should be in !
!real(8),dimension(lzd%lpsidimtot),intent(inout):: psi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: vpsi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
real(8),dimension(:,:),allocatable:: psir, vpsir
type(workarr_precond) :: work
type(workarrays_quartic_convolutions):: work_conv
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:,:,:),allocatable:: ypsitemp_c
real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
character(len=*),parameter:: subname='apply_orbitaldependent_potential'



  vpsi=0.d0
  ist_c=1
  ist_f=1
  do iorb=1,orbs%norbp
      iiorb=iorb+orbs%isorb
      ilr = orbs%inwhichlocreg(iiorb)
      if(centralLocreg<0) then
          !icenter=lin%orbs%inWhichLocregp(iorb)
          icenter=orbs%inWhichLocreg(iiorb)
      else
          icenter=centralLocreg
      end if
      !components of the potential
      npot=orbs%nspinor
      if (orbs%nspinor == 2) npot=1
      ist_f=ist_f+lzd%llr(ilr)%wfd%nvctr_c
      call allocate_workarrays_quartic_convolutions(lzd%llr(ilr), subname, work_conv)

      call uncompress_for_quartic_convolutions(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
           lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, &
           lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  & 
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  & 
           scal, psi(ist_c), psi(ist_f), &
           work_conv)

      if(confdatarr(iorb)%potorder==4) then
          call ConvolQuartic4(iproc,nproc,lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
               lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
               lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, & 
               hx, lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
               lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
               lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
               rxyz(1,ilr), confdatarr(iorb)%prefac, .false., 0.d0, &
               work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
               work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
               work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
               work_conv%y_c, work_conv%y_f)
      else if(confdatarr(iorb)%potorder==6) then
          call ConvolSextic(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
               lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
               lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, & 
               hx, lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
               lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
               lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
               rxyz(1,ilr), confdatarr(iorb)%prefac, .false., 0.d0, &
               work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
               work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
               work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
               work_conv%y_c, work_conv%y_f)
      else
          stop 'wronf conf pot'
      end if

      call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
           lzd%llr(Ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, &
           lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  & 
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  & 
           scal, work_conv%y_c, work_conv%y_f, vpsi(ist_c), vpsi(ist_f))

      call deallocate_workarrays_quartic_convolutions(lzd%llr(ilr), subname, work_conv)

      ist_f = ist_f + 7*lzd%llr(ilr)%wfd%nvctr_f
      ist_c = ist_c + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f

  end do



end subroutine apply_orbitaldependent_potential






!!!subroutine apply_orbitaldependent_potential_foronelocreg(iproc, nproc, lr, lin, at, input, orbs, rxyz, ndimpsi, &
!!!           psi, centralLocreg, vpsi)
!!!use module_base
!!!use module_types
!!!use module_interfaces
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc, centralLocreg, ndimpsi
!!!type(locreg_descriptors),intent(in):: lr
!!!type(linearParameters),intent(in):: lin
!!!type(atoms_data),intent(in):: at
!!!type(input_variables),intent(in):: input
!!!type(orbitals_data),intent(in):: orbs
!!!real(8),dimension(3,at%nat),intent(in):: rxyz
!!!real(8),dimension(ndimpsi),intent(inout):: psi
!!!real(8),dimension(ndimpsi),intent(out):: vpsi
!!!
!!!! Local variables
!!!integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
!!!real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
!!!type(workarr_sumrho):: work_sr
!!!real(8),dimension(:,:),allocatable:: psir, vpsir
!!!type(workarr_precond) :: work
!!!type(workarrays_quartic_convolutions):: work_conv
!!!character(len=*),parameter:: subname='apply_orbitaldependent_potential'
!!!real(8),dimension(0:3),parameter:: scal=1.d0
!!!real(8),dimension(:,:,:),allocatable:: ypsitemp_c
!!!real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
!!!
!!!
!!!
!!!  vpsi=0.d0
!!!  ist_c=1
!!!  ist_f=1
!!!     iiorb=iorb+orbs%isorb
!!!     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
!!!  
!!!     if(centralLocreg<0) then
!!!         !icenter=lin%orbs%inWhichLocregp(iorb)
!!!         icenter=lin%orbs%inWhichLocreg(iiorb)
!!!     else
!!!         icenter=centralLocreg
!!!     end if
!!!     ist_f=ist_f+lr%wfd%nvctr_c
!!!     call allocate_workarrays_quartic_convolutions(lr, subname, work_conv)
!!!     call uncompress_for_quartic_convolutions(lr%d%n1, lr%d%n2, lr%d%n3, &
!!!          lr%d%nfl1, lr%d%nfu1, &
!!!          lr%d%nfl2, lr%d%nfu2, &
!!!          lr%d%nfl3, lr%d%nfu3, &
!!!          lr%wfd%nseg_c, lr%wfd%nvctr_c, &
!!!          lr%wfd%keygloc, lr%wfd%keyv,  & 
!!!          lr%wfd%nseg_f, lr%wfd%nvctr_f, &
!!!          lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
!!!          lr%wfd%keyv(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)),  & 
!!!          scal, psi(ist_c), psi(ist_f), &
!!!          work_conv)
!!!
!!!     if(lin%confpotorder==4) then
!!!      call ConvolQuartic4(lr%d%n1, lr%d%n2, lr%d%n3, &
!!!           lr%d%nfl1, lr%d%nfu1, &
!!!           lr%d%nfl2, lr%d%nfu2, &
!!!           lr%d%nfl3, lr%d%nfu3, & 
!!!           input%hx, lr%ns1, lr%ns2, lr%ns3, &
!!!           lr%bounds%kb%ibyz_c, lr%bounds%kb%ibxz_c, lr%bounds%kb%ibxy_c, &
!!!           lr%bounds%kb%ibyz_f, lr%bounds%kb%ibxz_f, lr%bounds%kb%ibxy_f, &
!!!           rxyz(1,ilr), lin%potentialprefac(at%iatype(icenter)), .false., 0.d0, &
!!!           work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
!!!           work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
!!!           work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
!!!           work_conv%y_c, work_conv%y_f)
!!!      else if(lin%confpotorder==6) then
!!!
!!!      call ConvolSextic(lr%d%n1, lr%d%n2, lr%d%n3, &
!!!           lr%d%nfl1, lr%d%nfu1, &
!!!           lr%d%nfl2, lr%d%nfu2, &
!!!           lr%d%nfl3, lr%d%nfu3, & 
!!!           input%hx, lr%ns1, lr%ns2, lr%ns3, &
!!!           lr%bounds%kb%ibyz_c, lr%bounds%kb%ibxz_c, lr%bounds%kb%ibxy_c, &
!!!           lr%bounds%kb%ibyz_f, lr%bounds%kb%ibxz_f, lr%bounds%kb%ibxy_f, &
!!!           rxyz(1,ilr), lin%potentialprefac(at%iatype(icenter)), .false., 0.d0, &
!!!           work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
!!!           work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
!!!           work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
!!!           work_conv%y_c, work_conv%y_f)
!!!
!!!       else
!!!           stop 'wronf conf pot'
!!!
!!!       end if
!!!
!!!     call compress_forstandard(lr%d%n1, lr%d%n2, lr%d%n3, &
!!!          lr%d%nfl1, lr%d%nfu1, &
!!!          lr%d%nfl2, lr%d%nfu2, &
!!!          lr%d%nfl3, lr%d%nfu3, &
!!!          lr%wfd%nseg_c, lr%wfd%nvctr_c, &
!!!          lr%wfd%keygloc, lr%wfd%keyv,  & 
!!!          lr%wfd%nseg_f, lr%wfd%nvctr_f, &
!!!          lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
!!!          lr%wfd%keyv(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)),  & 
!!!          scal, work_conv%y_c, work_conv%y_f, vpsi(ist_c), vpsi(ist_f))
!!!
!!!     call deallocate_workarrays_quartic_convolutions(lr, subname, work_conv)
!!!
!!!
!!!end subroutine apply_orbitaldependent_potential_foronelocreg












!> Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
subroutine uncompress_for_quartic_convolutions(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  & 
     mseg_c, mvctr_c, keyg_c, keyv_c,  & 
     mseg_f, mvctr_f, keyg_f, keyv_f,  & 
     scal, psi_c, psi_f, &
     work)
  use module_base
  use module_types
  implicit none
  integer,intent(in):: n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, mseg_c, mvctr_c, mseg_f, mvctr_f
  integer,dimension(mseg_c),intent(in):: keyv_c
  integer,dimension(mseg_f),intent(in):: keyv_f
  integer,dimension(2,mseg_c),intent(in):: keyg_c
  integer,dimension(2,mseg_f),intent(in):: keyg_f
  real(wp),dimension(0:3),intent(in):: scal
  real(wp),dimension(mvctr_c),intent(in):: psi_c
  real(wp),dimension(7,mvctr_f),intent(in):: psi_f
  type(workarrays_quartic_convolutions),intent(out):: work
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !!!$omp parallel default(private) &
  !!!$omp shared(scal,psig_c,psig_f,x_f1,x_f2,x_f3) &
  !!!$omp shared(psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,mseg_c,mseg_f)
  !!! coarse part
  !!!$omp do
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        work%xx_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
        work%xy_c(i2,i,i3)=psi_c(i-i0+jj)*scal(0)
        work%xz_c(i3,i,i2)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo
  !!!$omp enddo
  !!! fine part
  !!!$omp do
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        work%xx_f1(i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xx_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xy_f(1,i2,i,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xz_f(1,i3,i,i2)=psi_f(1,i-i0+jj)*scal(1)

        work%xy_f2(i2,i,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xx_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xy_f(2,i2,i,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xz_f(2,i3,i,i2)=psi_f(2,i-i0+jj)*scal(1)

        work%xx_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        work%xy_f(3,i2,i,i3)=psi_f(3,i-i0+jj)*scal(2)
        work%xz_f(3,i3,i,i2)=psi_f(3,i-i0+jj)*scal(2)

        work%xz_f4(i3,i,i2)=psi_f(4,i-i0+jj)*scal(1)
        work%xx_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        work%xy_f(4,i2,i,i3)=psi_f(4,i-i0+jj)*scal(1)
        work%xz_f(4,i3,i,i2)=psi_f(4,i-i0+jj)*scal(1)

        work%xx_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        work%xy_f(5,i2,i,i3)=psi_f(5,i-i0+jj)*scal(2)
        work%xz_f(5,i3,i,i2)=psi_f(5,i-i0+jj)*scal(2)

        work%xx_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        work%xy_f(6,i2,i,i3)=psi_f(6,i-i0+jj)*scal(2)
        work%xz_f(6,i3,i,i2)=psi_f(6,i-i0+jj)*scal(2)

        work%xx_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
        work%xy_f(7,i2,i,i3)=psi_f(7,i-i0+jj)*scal(3)
        work%xz_f(7,i3,i,i2)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo
 !!!$omp enddo
 !!!$omp end parallel

END SUBROUTINE uncompress_for_quartic_convolutions





function dfactorial(n)
implicit none

! Calling arguments
integer,intent(in):: n
real(8):: dfactorial

! Local variables
integer:: i

  dfactorial=1.d0
  do i=1,n
      dfactorial=dfactorial*dble(i)
  end do

end function dfactorial

!!!!subroutine minimize_in_subspace(iproc, nproc, lin, at, input, lpot, GPU, ngatherarr, proj, rxyz, pkernelseq, nlpspd, lphi)
!!!!  use module_base
!!!!  use module_types
!!!!  use module_interfaces, exceptThisOne => minimize_in_subspace
!!!!  implicit none
!!!!  ! Calling arguments
!!!!  integer,intent(in):: iproc, nproc
!!!!  type(linearParameters),intent(inout):: lin
!!!!  type(atoms_data),intent(in):: at
!!!!  type(input_variables),intent(in):: input
!!!!  real(8),dimension(lin%lzd%ndimpotisf),intent(in):: lpot
!!!!  type(GPU_pointers),intent(inout):: GPU
!!!!  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
!!!!  type(nonlocal_psp_descriptors),intent(in):: nlpspd
!!!!  real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
!!!!  real(8),dimension(3,at%nat),intent(in):: rxyz
!!!!  real(dp), dimension(:), pointer :: pkernelseq
!!!!  real(8),dimension(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp)),intent(inout):: lphi
!!!!
!!!!  ! Local variables
!!!!  integer:: ndim_lhchi,jorb,jjorb,jlr,ii,istat,iall,jproc,iiorb,kk,iorb,norbTarget,nprocTemp
!!!!  integer:: is1,ie1,is2,ie2,is3,ie3,js1,je1,js2,je2,js3,je3,iat,ilr,ierr,tag,jlrold,nlocregPerMPI
!!!!  integer,dimension(:),allocatable:: onWhichAtomTemp,norb_parTemp,onWhichMPITemp
!!!!  logical,dimension(:),allocatable:: skip,doNotCalculate
!!!!  logical:: ovrlpx,ovrlpy,ovrlpz,check_whether_locregs_overlap,withConfinement
!!!!  real(8),dimension(:),allocatable:: lchi
!!!!  real(8),dimension(:,:),allocatable:: lhchi
!!!!  real(8):: ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,t1,t2,time,tt,tt1,tt2,tt3
!!!!  real(8),dimension(:,:,:),allocatable:: ham3
!!!!  character(len=*),parameter:: subname='minimize_in_subspace'
!!!!  type(confpot_data),dimension(:),allocatable :: confdatarr
!!!!
!!!!  !resetDIIS=.true.
!!!!  ! Apply the Hamiltonian for each atom.
!!!!  ! onWhichAtomTemp indicates that all orbitals feel the confining potential
!!!!  ! centered on atom iat.
!!!!  allocate(onWhichAtomTemp(lin%orbs%norb),stat=istat)
!!!!  call memocc(istat,onWhichAtomTemp,'onWhichAtomTemp',subname)
!!!!  allocate(doNotCalculate(lin%lzd%nlr),stat=istat)
!!!!  call memocc(istat,doNotCalculate,'doNotCalculate',subname)
!!!!  allocate(skip(lin%lzd%nlr), stat=istat)
!!!!  call memocc(istat, skip, 'skip', subname)
!!!!  !allocate(skipGlobal(lin%lig%lzdig%nlr,0:nproc-1), stat=istat)
!!!!  !call memocc(istat, skipGlobal, 'skipGlobal', subname)
!!!!
!!!!
!!!!  ! Determine for how many localization regions we need a Hamiltonian application.
!!!!  ndim_lhchi=0
!!!!  do iat=1,at%nat
!!!!     call getIndices(lin%lzd%llr(iat), is1, ie1, is2, ie2, is3, ie3)
!!!!     skip(iat)=.true.
!!!!     do jorb=1,lin%orbs%norbp
!!!!        jjorb=jorb+lin%orbs%isorb
!!!!        !onWhichAtomTemp(jorb)=iat
!!!!        onWhichAtomTemp(jjorb)=iat
!!!!        jlr=lin%orbs%inWhichLocreg(jjorb)
!!!!        if(lin%orbs%inWhichlocreg(jorb+lin%orbs%isorb)/=jlr) stop 'this should not happen'
!!!!        call getIndices(lin%lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!!!!        ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!!!        ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!!!        ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!!!        if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!!!           if(check_whether_locregs_overlap(lin%lzd%llr(iat), lin%lzd%llr(jlr), lin%lzd%glr)) then
!!!!              skip(iat)=.false.
!!!!           end if
!!!!        end if
!!!!     end do
!!!!     if(.not.skip(iat)) then
!!!!        ndim_lhchi=ndim_lhchi+1
!!!!     end if
!!!!  end do
!!!!
!!!!  allocate(lhchi(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp),ndim_lhchi),stat=istat)
!!!!  call memocc(istat, lhchi, 'lhchi', subname)
!!!!  lhchi=0.d0
!!!!
!!!!  if(iproc==0) write(*,'(1x,a)') 'Hamiltonian application for all atoms. This may take some time.'
!!!!  call mpi_barrier(mpi_comm_world, ierr)
!!!!  call cpu_time(t1)
!!!!
!!!!  allocate(lin%lzd%doHamAppl(lin%lzd%nlr), stat=istat)
!!!!  call memocc(istat, lin%lzd%doHamAppl, 'lin%lzd%doHamAppl', subname)
!!!!  withConfinement=.true.
!!!!  ii=0
!!!!  allocate(lchi(size(lphi)), stat=istat)
!!!!  call memocc(istat, lchi, 'lchi', subname)
!!!!  lchi=lphi
!!!!  do iat=1,at%nat
!!!!     doNotCalculate=.true.
!!!!     lin%lzd%doHamAppl=.false.
!!!!     !!call mpi_barrier(mpi_comm_world, ierr)
!!!!     call getIndices(lin%lzd%llr(iat), is1, ie1, is2, ie2, is3, ie3)
!!!!     skip(iat)=.true.
!!!!     do jorb=1,lin%orbs%norbp
!!!!        !onWhichAtomTemp(jorb)=iat
!!!!        onWhichAtomTemp(lin%orbs%isorb+jorb)=iat
!!!!        !jlr=onWhichAtomp(jorb)
!!!!        !jlr=lin%orbs%inWhichLocregp(jorb)
!!!!        jjorb=lin%orbs%isorb+jorb
!!!!        jlr=lin%orbs%inWhichLocreg(jjorb)
!!!!        call getIndices(lin%lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!!!!        ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!!!        ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!!!        ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!!!        if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!!!           doNotCalculate(jlr)=.false.
!!!!           lin%lzd%doHamAppl(jlr)=.true.
!!!!           skip(iat)=.false.
!!!!        else
!!!!           doNotCalculate(jlr)=.true.
!!!!           lin%lzd%doHamAppl(jlr)=.false.
!!!!        end if
!!!!     end do
!!!!     !write(*,'(a,2i4,4x,100l4)') 'iat, iproc, doNotCalculate', iat, iproc, doNotCalculate
!!!!     if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'Hamiltonian application for atom ', iat, '... '
!!!!     if(.not.skip(iat)) then
!!!!        ii=ii+1
!!!!        if(lin%nItInguess>0) then
!!!!!!$           call HamiltonianApplication3(iproc, nproc, at, lin%orbs, input%hx, input%hy, input%hz, rxyz, &
!!!!!!$                proj, lin%lzd, ngatherarr, lpot, lchi, lhchi(1,ii), &
!!!!!!$                ekin_sum, epot_sum, eexctX, eproj_sum, input%nspin, GPU, withConfinement, .false., &
!!!!!!$                pkernel=pkernelseq, lin=lin, confinementCenter=onWhichAtomTemp)
!!!!        
!!!!           !confdatarr to be initialized
!!!!           allocate(confdatarr(lin%orbs%norbp))
!!!!           call define_confinement_data(confdatarr,lin%orbs,rxyz,at,&
!!!!                input%hx,input%hy,input%hz,lin,lin%Lzd,onWhichAtomTemp)
!!!!
!!!!           call LocalHamiltonianApplication(iproc,nproc,at,lin%orbs,&
!!!!                input%hx,input%hy,input%hz,&
!!!!                lin%lzd,confdatarr,ngatherarr,Lpot,lchi,lhchi(1,ii),&
!!!!                ekin_sum,epot_sum,eexctX,eSIC_DC,input%SIC,GPU,&
!!!!                pkernel=pkernelseq)
!!!!
!!!!           call NonLocalHamiltonianApplication(iproc,at,lin%orbs,&
!!!!                input%hx,input%hy,input%hz,rxyz,&
!!!!                proj,lin%lzd,nlpspd,lchi,lhchi(1,ii),eproj_sum)
!!!!           deallocate(confdatarr)
!!!!           print *,'iproc,energies',ekin_sum,epot_sum,eproj_sum 
!!!!        end if
!!!!
!!!!     else
!!!!     end if
!!!!
!!!!
!!!!     if(iproc==0) write(*,'(a)') 'done.'
!!!!  end do
!!!!
!!!!
!!!!  !!iall=-product(shape(lpot))*kind(lpot)
!!!!  !!deallocate(lpot, stat=istat)
!!!!  !!call memocc(istat, iall, 'lpot', subname)
!!!!  if(ii/=ndim_lhchi) then
!!!!     write(*,'(a,i0,a,2(a2,i0))') 'ERROR on process ',iproc,': ii/=ndim_lhchi',ii,ndim_lhchi
!!!!     stop
!!!!  end if
!!!!  call mpi_barrier(mpi_comm_world, ierr)
!!!!  call cpu_time(t2)
!!!!  time=t2-t1
!!!!  if(iproc==0) write(*,'(1x,a,es10.3)') 'time for applying potential:', time
!!!!
!!!!
!!!!
!!!!  ! The input guess is possibly performed only with a subset of all processes.
!!!!  if(lin%norbsPerProcIG>lin%orbs%norb) then
!!!!     norbTarget=lin%orbs%norb
!!!!  else
!!!!     norbTarget=lin%norbsperProcIG
!!!!  end if
!!!!  nprocTemp=ceiling(dble(lin%orbs%norb)/dble(norbTarget))
!!!!  nprocTemp=min(nprocTemp,nproc)
!!!!  if(iproc==0) write(*,'(a,i0,a)') 'The minimization is performed using ', nprocTemp, ' processes.'
!!!!
!!!!  ! Create temporary norb_parTemp, onWhichMPITemp
!!!!  allocate(norb_parTemp(0:nprocTemp-1), stat=istat)
!!!!  call memocc(istat, norb_parTemp, 'norb_parTemp', subname)
!!!!  norb_parTemp=0
!!!!  tt=dble(lin%orbs%norb)/dble(nprocTemp)
!!!!  ii=floor(tt)
!!!!  ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
!!!!  norb_parTemp(0:nprocTemp-1)=ii
!!!!  kk=lin%orbs%norb-nprocTemp*ii
!!!!  norb_parTemp(0:kk-1)=ii+1
!!!!
!!!!  allocate(onWhichMPITemp(lin%orbs%norb), stat=istat)
!!!!  call memocc(istat, onWhichMPITemp, 'onWhichMPITemp', subname)
!!!!  iiorb=0
!!!!  do jproc=0,nprocTemp-1
!!!!     do iorb=1,norb_parTemp(jproc)
!!!!        iiorb=iiorb+1
!!!!        onWhichMPITemp(iiorb)=jproc
!!!!     end do
!!!!  end do
!!!!
!!!!  ! Calculate the number of different matrices that have to be stored on a given MPI process.
!!!!  jlrold=0
!!!!  nlocregPerMPI=0
!!!!  do jorb=1,lin%orbs%norb
!!!!     jlr=lin%orbs%inWhichLocreg(jorb)
!!!!     !jproc=lin%orbs%onWhichMPI(jorb)
!!!!     jproc=onWhichMPITemp(jorb)
!!!!     !if(iproc==0) write(*,'(a,5i7)') 'jorb, jlr, jlrold, jproc, nlocregPerMPI', jorb, jlr, jlrold, jproc, nlocregPerMPI
!!!!     if(iproc==jproc) then
!!!!        if(jlr/=jlrold) then
!!!!           nlocregPerMPI=nlocregPerMPI+1
!!!!           jlrold=jlr
!!!!        end if
!!!!     end if
!!!!  end do
!!!!
!!!!
!!!!
!!!!  ! Calculate the Hamiltonian matrix.
!!!!  call cpu_time(t1)
!!!!  allocate(ham3(lin%orbs%norb,lin%orbs%norb,nlocregPerMPI), stat=istat)
!!!!  call memocc(istat,ham3,'ham3',subname)
!!!!  if(lin%nItInguess>0) then
!!!!     if(iproc==0) write(*,*) 'calling getHamiltonianMatrix6'
!!!!     call getHamiltonianMatrix6(iproc, nproc, nprocTemp, lin%lzd, lin%orbs, lin%orbs, &
!!!!          onWhichMPITemp, input, lin%orbs%inWhichLocreg, ndim_lhchi, &
!!!!          nlocregPerMPI, lchi, lhchi, skip, lin%mad, &
!!!!          lin%memoryForCommunOverlapIG, lin%locregShape, tag, ham3)
!!!!  end if
!!!!
!!!!  iall=-product(shape(lhchi))*kind(lhchi)
!!!!  deallocate(lhchi, stat=istat)
!!!!  call memocc(istat, iall, 'lhchi',subname)
!!!!
!!!!
!!!!  ! Build the orbitals phi as linear combinations of the atomic orbitals.
!!!!  if(iproc==0) write(*,*) 'calling buildLinearCombinationsLocalized3'
!!!!  call buildLinearCombinationsLocalized3(iproc, nproc, lin%orbs, lin%orbs, lin%comms,&
!!!!       at, lin%lzd%Glr, input, lin%norbsPerType, &
!!!!       lin%orbs%inWhichLocreg, lchi, lphi, rxyz, lin%orbs%inWhichLocreg, &
!!!!       lin, lin%lzd, nlocregPerMPI, tag, ham3)
!!!!
!!!!  iall=-product(shape(lchi))*kind(lchi)
!!!!  deallocate(lchi, stat=istat)
!!!!  call memocc(istat, iall, 'lchi',subname)
!!!!
!!!!  iall=-product(shape(lin%lzd%doHamAppl))*kind(lin%lzd%doHamAppl)
!!!!  deallocate(lin%lzd%doHamAppl, stat=istat)
!!!!  call memocc(istat, iall, 'lin%lzd%doHamAppl',subname)
!!!!
!!!!  iall=-product(shape(norb_parTemp))*kind(norb_parTemp)
!!!!  deallocate(norb_parTemp, stat=istat)
!!!!  call memocc(istat, iall, 'norb_parTemp',subname)
!!!!
!!!!  iall=-product(shape(ham3))*kind(ham3)
!!!!  deallocate(ham3, stat=istat)
!!!!  call memocc(istat, iall, 'ham3',subname)
!!!!
!!!!  ! Deallocate all remaining local arrays.
!!!!  iall=-product(shape(onWhichAtomTemp))*kind(onWhichAtomTemp)
!!!!  deallocate(onWhichAtomTemp, stat=istat)
!!!!  call memocc(istat, iall, 'onWhichAtomTemp',subname)
!!!!
!!!!  iall=-product(shape(doNotCalculate))*kind(doNotCalculate)
!!!!  deallocate(doNotCalculate, stat=istat)
!!!!  call memocc(istat, iall, 'doNotCalculate',subname)
!!!!
!!!!  iall=-product(shape(skip))*kind(skip)
!!!!  deallocate(skip, stat=istat)
!!!!  call memocc(istat, iall, 'skip',subname)
!!!!
!!!!  iall=-product(shape(onWhichMPITemp))*kind(onWhichMPITemp)
!!!!  deallocate(onWhichMPITemp, stat=istat)
!!!!  call memocc(istat, iall, 'onWhichMPITemp',subname)
!!!!
!!!!
!!!!end subroutine minimize_in_subspace








!> Some description of the routine goes here
subroutine unitary_optimization(iproc, nproc, lzd, orbs, at, op, comon, mad, rxyz, nit, kernel, &
           newgradient, confdatarr, hx, lphi)
use module_base
use module_types
use module_interfaces, exceptThisOne => unitary_optimization
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nit
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(atoms_data),intent(in):: at
type(overlapParameters),intent(inout):: op
type(p2pComms),intent(inout):: comon
type(matrixDescriptors),intent(in):: mad
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(orbs%norb,orbs%norb),intent(in):: kernel
logical,intent(in):: newgradient
real(8),intent(in):: hx
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: lphi

! Local variables
integer:: it, info, lwork, k, istat, iorb, jorb, iall, ierr, ist, jst, ilrold, ncount, jjorb, iiorb, ilr, lorb, jlr
integer:: nlocregOnMPI, jlrold, jj, ii
real(8):: trace, lstep, dfactorial, energyconf_trial, energyconf_0, energyconf_der0, lstep_optimal, ddot
real(8):: tt1, tt2, tt3, tt4, tt5, tt
real(8):: t1, t2, t1_tot, t2_tot
real(8):: time_convol, time_commun, time_lincomb, time_linalg, time_matrixmodification, time_exponential, time_tot
real(8):: time_matrixelements
complex(8):: ttc
real(8),dimension(:,:),allocatable:: gmat, hamtrans, ttmat, Kmat
real(8),dimension(:,:,:),allocatable:: potmat, potmatsmall
complex(8),dimension(:,:),allocatable:: gmatc, omatc
complex(8),dimension(:,:,:),allocatable:: tempmatc
complex(8),dimension(:),allocatable:: work, expD_cmplx
real(8),dimension(:),allocatable:: rwork, eval, lphiovrlp, lvphi, recvbuf
real(8),dimension(:,:,:),allocatable:: tempmat3
character(len=*),parameter:: subname='unitary_optimization'
type(p2pComms):: comon_local

! Quick return if possible
if(nit==0) return

allocate(gmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, gmat, 'gmat', subname)
allocate(hamtrans(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, hamtrans, 'hamtrans', subname)
allocate(gmatc(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, gmatc, 'gmatc', subname)
allocate(omatc(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, omatc, 'omatc', subname)
allocate(tempmat3(orbs%norb,orbs%norb,3), stat=istat)
call memocc(istat, tempmat3, 'tempmat3', subname)
allocate(eval(orbs%norb), stat=istat)
call memocc(istat, eval, 'eval', subname)
allocate(expD_cmplx(orbs%norb), stat=istat)
call memocc(istat, expD_cmplx, 'expD_cmplx', subname)
allocate(tempmatc(orbs%norb,orbs%norb,2), stat=istat)
!call memocc(istat, tempmatc, 'tempmatc', subname)
allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
allocate(lvphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, lvphi, 'lvphi', subname)
allocate(Kmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Kmat, 'Kmat', subname)
allocate(recvbuf(comon%nrecvbuf), stat=istat)
call memocc(istat, recvbuf, 'recvbuf', subname)

allocate(potmat(orbs%norb,orbs%norb,at%nat), stat=istat)
call memocc(istat, potmat, 'potmat', subname)


! Count how many locregs each process handles
ilrold=-1
nlocregOnMPI=0
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inwhichlocreg(iiorb)
    !if(ilr>ilrold) then
    if(ilr/=ilrold) then
        nlocregOnMPI=nlocregOnMPI+1
    end if
    ilrold=ilr
end do
allocate(potmatsmall(orbs%norb,orbs%norb,nlocregOnMPI), stat=istat)
call memocc(istat, potmatsmall, 'potmatsmall', subname)



  call allocateSendBufferOrtho(comon, subname)
  call allocateRecvBufferOrtho(comon, subname)
  ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
  !if(nit>0) then
  !    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
  !    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
  !end if

  energyconf_trial=0.d0 !just to initialize this variable and make the compiler happy
  lstep=0.d0 !just to initialize this variable and make the compiler happy
  lstep_optimal=0.d0 !just to initialize this variable and make the compiler happy
  energyconf_der0=0.d0 !just to initialize this variable and make the compiler happy

  time_convol=0.d0
  time_lincomb=0.d0
  time_commun=0.d0
  time_linalg=0.d0
  time_exponential=0.d0
  time_matrixmodification=0.d0
  time_matrixElements=0.d0
  t1_tot=mpi_wtime()
  innerLoop: do it=1,nit

  !write(*,*) '1: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)

      t1=mpi_wtime()
      call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
           confdatarr, hx, lphi, -1, lvphi)
      t2=mpi_wtime()
      time_convol=time_convol+t2-t1

      t1=mpi_wtime()
      !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
      !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
      !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
      !deallocate(ttmat)
      !write(*,*) '2: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      t2=mpi_wtime()
      time_commun=time_commun+t2-t1

      t1=mpi_wtime()
      !call getMatrixElements2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%op, lin%lb%comon, lphi, lvphi, lin%mad, Kmat)
      !call deallocateRecvBufferOrtho(comon, subname)
      !call deallocateSendBufferOrtho(comon, subname)
      call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, Kmat)
      !do iorb=1,orbs%norb
      !    do jorb=1,orbs%norb
      !        if(iproc==0) write(66,*) iorb,jorb,Kmat(jorb,iorb)
      !    end do
      !end do
      !call allocateSendBufferOrtho(comon, subname)
      !call allocateRecvBufferOrtho(comon, subname)
      !write(*,*) '3: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      t2=mpi_wtime()
      time_matrixElements=time_matrixElements+t2-t1


      !!if(newgradient) then
      !!    call get_potential_matrices(iproc, nproc, at, orbs, lzd, op, comon, mad, rxyz, &
      !!         confdatarr, hx, lphi, potmat)
      !!    !call get_potential_matrices_new(iproc, nproc, lin, at, input, orbs, lzd, op, comon, rxyz, lphi, &
      !!    !     nlocregOnMPI, potmatsmall)
      !!end if


      !!if(.not.newgradient) then
          !energyconf_0=ddot(orbs%npsidim, lphi(1), 1, lvphi(1), 1)
          !call mpiallred(energyconf_0, 1, mpi_sum, mpi_comm_world, ierr)
          energyconf_0=0.d0
          do iorb=1,orbs%norb
              energyconf_0 = energyconf_0 + Kmat(iorb,iorb)
          end do
      !!else
      !!    energyconf_0=0.d0
      !!    do iorb=1,orbs%norb
      !!        do jorb=1,orbs%norb
      !!            energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(jorb,iorb)
      !!            !energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(iorb,jorb)
      !!        end do
      !!    end do
      !!end if
      if(iproc==0) write(*,'(a,i6,3es20.10,2es17.7)') &
                   'it, energyconf_0, energyvonf_trial, energyconf_der0, lstep, lstep_optimal', &
                   it, energyconf_0, energyconf_trial, energyconf_der0, lstep, lstep_optimal

      t1=mpi_wtime()
      !!if(.not.newgradient) then
          ! Construct antisymmtric matrix Gmat
          do iorb=1,orbs%norb
              do jorb=1,orbs%norb
                  gmat(jorb,iorb)=2.d0*(Kmat(jorb,iorb)-Kmat(iorb,jorb))
              end do
          end do 
      !!else
      !!    !!! THIS IS THE OLD VERSION #############################################################################
      !!    do iorb=1,orbs%norb
      !!        ilr=orbs%inwhichlocreg(iorb)
      !!        do jorb=1,orbs%norb
      !!            jlr=orbs%inwhichlocreg(jorb)
      !!            tt=0.d0
      !!            do lorb=1,orbs%norb
      !!                tt = tt + kernel(jorb,lorb)*Kmat(lorb,iorb) - kernel(iorb,lorb)*Kmat(lorb,jorb) + &
      !!                          kernel(jorb,lorb)*potmat(lorb,iorb,jlr) - kernel(iorb,lorb)*potmat(lorb,jorb,ilr)
      !!            end do
      !!            gmat(jorb,iorb)=-tt
      !!            !if(iproc==0) then
      !!            !    write(77,*) iorb, jorb, gmat(jorb,iorb)
      !!            !end if
      !!        end do
      !!    end do 
      !!    ! ########################################################################################################
      !!    !!! THIS IS THE NEW VERSION
      !!    !!gmat=0.d0
      !!    !!ii=0
      !!    !!ilrold=-1
      !!    !!do iorb=1,orbs%norbp
      !!    !!    iiorb=orbs%isorb+iorb
      !!    !!    ilr=orbs%inwhichlocreg(iiorb)
      !!    !!    if(ilr>ilrold) then
      !!    !!        ii=ii+1
      !!    !!    end if
      !!    !!    do jorb=1,orbs%norb
      !!    !!        jlr=orbs%inwhichlocreg(jorb)
      !!    !!        tt=0.d0
      !!    !!        do lorb=1,orbs%norb
      !!    !!            !tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
      !!    !!            !          - kernel(iiorb,lorb)*potmat(lorb,jorb,ilr)
      !!    !!            tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
      !!    !!                      - kernel(iiorb,lorb)*potmatsmall(lorb,jorb,ii)
      !!    !!        end do
      !!    !!        gmat(jorb,iiorb)=-tt
      !!    !!    end do
      !!    !!    ilrold=ilr
      !!    !!end do 
      !!    !!do iorb=1,orbs%norb
      !!    !!    ilr=orbs%inwhichlocreg(iorb)
      !!    !!    jlrold=-1
      !!    !!    jj=0
      !!    !!    do jorb=1,orbs%norbp
      !!    !!        jjorb=orbs%isorb+jorb
      !!    !!        jlr=orbs%inwhichlocreg(jjorb)
      !!    !!        if(jlr>jlrold) then
      !!    !!            jj=jj+1
      !!    !!        end if
      !!    !!        tt=0.d0
      !!    !!        do lorb=1,orbs%norb
      !!    !!            !tt = tt + kernel(jjorb,lorb)*potmat(lorb,iorb,jlr)
      !!    !!            tt = tt + kernel(jjorb,lorb)*potmatsmall(lorb,iorb,jj)
      !!    !!        end do
      !!    !!        gmat(jjorb,iorb)=gmat(jjorb,iorb)-tt
      !!    !!        jlrold=jlr
      !!    !!    end do
      !!    !!end do 
      !!    !!call mpiallred(gmat(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
      !!    !!do iorb=1,orbs%norb
      !!    !!    do jorb=1,orbs%norb
      !!    !!        if(iproc==0) then
      !!    !!            write(77,*) iorb, jorb, gmat(jorb,iorb)
      !!    !!        end if
      !!    !!    end do
      !!    !!end do
      !!end if
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1


      t1=mpi_wtime()
      !Build the complex matrix -iGmat
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              gmatc(jorb,iorb)=cmplx(0.d0,-gmat(jorb,iorb),kind=8)
          end do
      end do 
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1



      ! Diagonalize Gmatc
      t1=mpi_wtime()
      lwork=10*orbs%norb
      allocate(work(lwork), stat=istat) ! factor of 2 since it is assumed to be complex
      allocate(rwork(lwork), stat=istat)
      call zheev('v', 'l', orbs%norb, gmatc(1,1), orbs%norb, eval(1), work, lwork, rwork, info)
      if(info/=0) stop 'ERROR in zheev 1'
      deallocate(work)
      deallocate(rwork)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1


      ! Calculate step size
      if(it==1) then
          !!if(.not.newgradient) then
              lstep=5.d-2/(maxval(eval))
          !!else
          !!    lstep=1.d-4/(maxval(eval))
          !!end if
      else
          lstep=2.d0*lstep_optimal
          !lstep=1.d-3/(maxval(eval))
      end if

      t1=mpi_wtime()
      ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
      ! This is also a diagonal matrix, so only calculate the diagonal part.
      do iorb=1,orbs%norb
         ttc=cmplx(0.d0,-lstep*eval(iorb),kind=8)
         expD_cmplx(iorb)=(0.d0,0.d0)
          do k=0,50
              expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
          end do
      end do
      t2=mpi_wtime()
      time_exponential=time_exponential+t2-t1

      t1=mpi_wtime()
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              if(iorb==jorb) then
                  tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
              else
                  tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
              end if
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
           gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
      call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
           tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1

      t1=mpi_wtime()
      ! Build new lphi
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      !write(*,*) '5: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
           comon%recvbuf, tempmat3(1,1,1), .true., lphi)
      t2=mpi_wtime()
      time_lincomb=time_lincomb+t2-t1

      t1=mpi_wtime()
      call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
           confdatarr, hx, lphi, -1, lvphi)
      t2=mpi_wtime()
      time_convol=time_convol+t2-t1

      !!if(.not.newgradient) then
          energyconf_trial=ddot(max(orbs%npsidim_orbs,orbs%npsidim_comp), lphi(1), 1, lvphi(1), 1)
          call mpiallred(energyconf_trial, 1, mpi_sum, mpi_comm_world, ierr)
      !!else

      !!    call dcopy(comon%nrecvbuf, comon%recvbuf, 1, recvbuf, 1)
      !!    !call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      !!    !call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
      !!    !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
      !!    !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
      !!    !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
      !!    !deallocate(ttmat)
      !!    call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, Kmat)
      !!    call dcopy(comon%nrecvbuf, recvbuf, 1, comon%recvbuf, 1)

      !!    energyconf_trial=0.d0
      !!    do iorb=1,orbs%norb
      !!        do jorb=1,orbs%norb
      !!            energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(jorb,iorb)
      !!            !energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(iorb,jorb)
      !!        end do
      !!    end do
      !!end if

      ! Calculate the gradient of the confinement
      energyconf_der0=0.d0
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              energyconf_der0=energyconf_der0+gmat(jorb,iorb)**2
          end do
      end do
      energyconf_der0=-.5d0*energyconf_der0

      ! Calculate optimal step size
      lstep_optimal = -energyconf_der0*lstep**2/(2.d0*(energyconf_trial-energyconf_0-lstep*energyconf_der0))
      !!if(.not.newgradient) then
          lstep_optimal=min(lstep_optimal,lstep)
      !!else
      !!    if(lstep_optimal<0) then
      !!        lstep_optimal=lstep
      !!    else
      !!        lstep_optimal=min(lstep_optimal,lstep)
      !!    end if
      !!end if

      t1=mpi_wtime()
      ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
      ! This is also a diagonal matrix, so only calculate the diagonal part.
      do iorb=1,orbs%norb
         ttc=cmplx(0.d0,-lstep_optimal*eval(iorb),kind=8)
         expD_cmplx(iorb)=(0.d0,0.d0)
          do k=0,50
              expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
          end do
      end do
      t2=mpi_wtime()
      time_exponential=time_exponential+t2-t1

      t1=mpi_wtime()
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              if(iorb==jorb) then
                  tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
              else
                  tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
              end if
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
           gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
      call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
           tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1


      ! Build new lphi
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
          end do
      end do
      t1=mpi_wtime()
      call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
           comon%recvbuf, tempmat3(1,1,1), .true., lphi)
      t2=mpi_wtime()
      time_lincomb=time_lincomb+t2-t1


      !if(it<nit) then
      !    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      !    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
      !end if


  end do innerLoop

  t2_tot=mpi_wtime()
  time_tot=t2_tot-t1_tot
  call mpiallred(time_convol, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_commun, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_lincomb, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_linalg, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_matrixmodification, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_exponential, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_matrixelements, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_tot, 1, mpi_max, mpi_comm_world, ierr)
  !time_convol=time_convol/dble(nproc)
  !time_commun=time_commun/dble(nproc)
  !time_lincomb=time_lincomb/dble(nproc)
  !time_linalg=time_linalg/dble(nproc)
  !time_matrixmodification=time_matrixmodification/dble(nproc)
  !time_exponential_=time_exponential/dble(nproc)
  !time_tot=time_tot/dble(nproc)
  !!if(iproc==0) then
  !!    write(*,'(a,es16.6)') 'total time: ',time_tot
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'convolutions: ',time_convol,' (',time_convol/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear combinations: ',time_lincomb,' (',time_lincomb/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'communication: ',time_commun,' (',time_commun/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear algebra: ',time_linalg,' (',time_linalg/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix modification: ',time_matrixmodification, &
  !!                                   ' (',time_matrixmodification/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'building exponential: ',time_exponential,' (',time_exponential/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix elements ',time_matrixelements,' (',time_matrixelements/time_tot*100.d0,'%)'
  !!end if


  iall=-product(shape(gmat))*kind(gmat)
  deallocate(gmat, stat=istat)
  call memocc(istat, iall, 'gmat', subname)
  iall=-product(shape(gmatc))*kind(gmatc)
  deallocate(gmatc, stat=istat)
  !call memocc(istat, iall, 'gmatc', subname)
  iall=-product(shape(omatc))*kind(omatc)
  deallocate(omatc, stat=istat)
  !call memocc(istat, iall, 'omatc', subname)
  iall=-product(shape(tempmat3))*kind(tempmat3)
  deallocate(tempmat3, stat=istat)
  call memocc(istat, iall, 'tempmat3', subname)
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)
  iall=-product(shape(expD_cmplx))*kind(expD_cmplx)
  deallocate(expD_cmplx, stat=istat)
  call memocc(istat, iall, 'expD_cmplx', subname)
  iall=-product(shape(tempmatc))*kind(tempmatc)
  deallocate(tempmatc, stat=istat)
  !call memocc(istat, iall, 'tempmatc', subname)
  iall=-product(shape(hamtrans))*kind(hamtrans)
  deallocate(hamtrans, stat=istat)
  call memocc(istat, iall, 'hamtrans', subname)
  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  deallocate(lphiovrlp, stat=istat)
  call memocc(istat, iall, 'lphiovrlp', subname)
  iall=-product(shape(lvphi))*kind(lvphi)
  deallocate(lvphi, stat=istat)
  call memocc(istat, iall, 'lvphi', subname)
  iall=-product(shape(Kmat))*kind(Kmat)
  deallocate(Kmat, stat=istat)
  call memocc(istat, iall, 'Kmat', subname)
  iall=-product(shape(recvbuf))*kind(recvbuf)
  deallocate(recvbuf, stat=istat)
  call memocc(istat, iall, 'recvbuf', subname)

  iall=-product(shape(potmat))*kind(potmat)
  deallocate(potmat, stat=istat)
  call memocc(istat, iall, 'potmat', subname)
  iall=-product(shape(potmatsmall))*kind(potmatsmall)
  deallocate(potmatsmall, stat=istat)
  call memocc(istat, iall, 'potmatsmall', subname)

  call deallocateRecvBufferOrtho(comon, subname)
  call deallocateSendBufferOrtho(comon, subname)

end subroutine unitary_optimization



!!!!subroutine MLWF(iproc, nproc, lzd, orbs, at, op, comon, mad, rxyz, nit, kernel, &
!!!!           newgradient, confdatarr, hx, lphi, Umat)
!!!!use module_base
!!!!use module_types
!!!!use module_interfaces, exceptThisOne => MLWF
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc, nit
!!!!type(local_zone_descriptors),intent(in):: lzd
!!!!type(orbitals_data),intent(in):: orbs
!!!!type(atoms_data),intent(in):: at
!!!!type(overlapParameters),intent(inout):: op
!!!!type(p2pComms),intent(inout):: comon
!!!!type(matrixDescriptors),intent(in):: mad
!!!!real(8),dimension(3,at%nat),intent(in):: rxyz
!!!!real(8),dimension(orbs%norb,orbs%norb),intent(in):: kernel
!!!!logical,intent(in):: newgradient
!!!!real(8),intent(in):: hx
!!!!type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
!!!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: lphi
!!!!real(8),dimension(orbs%norb,orbs%norb),intent(out):: Umat
!!!!
!!!!! Local variables
!!!!integer:: it, info, lwork, k, istat, iorb, jorb, iall, ierr, ist, jst, ilrold, ncount, jjorb, iiorb, ilr, lorb, jlr
!!!!integer:: nlocregOnMPI, jlrold, jj, ii
!!!!real(8):: trace, lstep, dfactorial, energyconf_trial, energyconf_0, energyconf_der0, lstep_optimal, ddot
!!!!real(8):: tt1, tt2, tt3, tt4, tt5, tt
!!!!real(8):: t1, t2, t1_tot, t2_tot, omega
!!!!real(8):: time_convol, time_commun, time_lincomb, time_linalg, time_matrixmodification, time_exponential, time_tot
!!!!real(8):: time_matrixelements
!!!!complex(8):: ttc
!!!!real(8),dimension(:,:),allocatable:: gmat, hamtrans, ttmat, Kmat, X, Y, Z, Xd, Yd, Zd, commutX, commutY, commutZ
!!!!real(8),dimension(:,:),allocatable:: Xprime, Yprime, Zprime, Xprimesquare, Yprimesquare, Zprimesquare
!!!!real(8),dimension(:,:,:),allocatable:: potmat, potmatsmall
!!!!complex(8),dimension(:,:),allocatable:: gmatc, omatc
!!!!complex(8),dimension(:,:,:),allocatable:: tempmatc
!!!!complex(8),dimension(:),allocatable:: work, expD_cmplx
!!!!real(8),dimension(:),allocatable:: rwork, eval, lphiovrlp, lvphi, recvbuf, lxphi, lyphi, lzphi
!!!!real(8),dimension(:,:,:),allocatable:: tempmat3
!!!!character(len=*),parameter:: subname='unitary_optimization'
!!!!type(p2pComms):: comon_local
!!!!
!!!!! Quick return if possible
!!!!if(nit==0) return
!!!!
!!!!allocate(gmat(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, gmat, 'gmat', subname)
!!!!allocate(hamtrans(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, hamtrans, 'hamtrans', subname)
!!!!allocate(gmatc(orbs%norb,orbs%norb), stat=istat)
!!!!!call memocc(istat, gmatc, 'gmatc', subname)
!!!!allocate(omatc(orbs%norb,orbs%norb), stat=istat)
!!!!!call memocc(istat, omatc, 'omatc', subname)
!!!!allocate(tempmat3(orbs%norb,orbs%norb,3), stat=istat)
!!!!call memocc(istat, tempmat3, 'tempmat3', subname)
!!!!allocate(eval(orbs%norb), stat=istat)
!!!!call memocc(istat, eval, 'eval', subname)
!!!!allocate(expD_cmplx(orbs%norb), stat=istat)
!!!!call memocc(istat, expD_cmplx, 'expD_cmplx', subname)
!!!!allocate(tempmatc(orbs%norb,orbs%norb,2), stat=istat)
!!!!!call memocc(istat, tempmatc, 'tempmatc', subname)
!!!!allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
!!!!call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
!!!!allocate(lvphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!!call memocc(istat, lvphi, 'lvphi', subname)
!!!!allocate(Kmat(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Kmat, 'Kmat', subname)
!!!!allocate(recvbuf(comon%nrecvbuf), stat=istat)
!!!!call memocc(istat, recvbuf, 'recvbuf', subname)
!!!!
!!!!allocate(potmat(orbs%norb,orbs%norb,at%nat), stat=istat)
!!!!call memocc(istat, potmat, 'potmat', subname)
!!!!
!!!!allocate(lxphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!!call memocc(istat, lxphi, 'lxphi', subname)
!!!!allocate(lyphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!!call memocc(istat, lyphi, 'lyphi', subname)
!!!!allocate(lzphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!!call memocc(istat, lzphi, 'lzphi', subname)
!!!!
!!!!allocate(X(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, X, 'X', subname)
!!!!allocate(Y(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Y, 'Y', subname)
!!!!allocate(Z(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Z, 'Z', subname)
!!!!allocate(Xd(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Xd, 'Xd', subname)
!!!!allocate(Yd(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Yd, 'Yd', subname)
!!!!allocate(Zd(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Zd, 'Zd', subname)
!!!!allocate(Xprime(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Xprime, 'Xprime', subname)
!!!!allocate(Yprime(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Yprime, 'Yprime', subname)
!!!!allocate(Zprime(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Zprime, 'Zprime', subname)
!!!!allocate(Xprimesquare(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Xprimesquare, 'Xprimesquare', subname)
!!!!allocate(Yprimesquare(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Yprimesquare, 'Yprimesquare', subname)
!!!!allocate(Zprimesquare(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Zprimesquare, 'Zprimesquare', subname)
!!!!allocate(commutX(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, commutX, 'commutX', subname)
!!!!allocate(commutY(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, commutY, 'commutY', subname)
!!!!allocate(commutZ(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, commutZ, 'commutZ', subname)
!!!!
!!!!
!!!!! Count how many locregs each process handles
!!!!ilrold=-1
!!!!nlocregOnMPI=0
!!!!do iorb=1,orbs%norbp
!!!!    iiorb=orbs%isorb+iorb
!!!!    ilr=orbs%inwhichlocreg(iiorb)
!!!!    !if(ilr>ilrold) then
!!!!    if(ilr/=ilrold) then
!!!!        nlocregOnMPI=nlocregOnMPI+1
!!!!    end if
!!!!    ilrold=ilr
!!!!end do
!!!!allocate(potmatsmall(orbs%norb,orbs%norb,nlocregOnMPI), stat=istat)
!!!!call memocc(istat, potmatsmall, 'potmatsmall', subname)
!!!!
!!!!
!!!!
!!!!  call allocateSendBufferOrtho(comon, subname)
!!!!  call allocateRecvBufferOrtho(comon, subname)
!!!!  ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
!!!!  !if(nit>0) then
!!!!  !    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
!!!!  !    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
!!!!  !end if
!!!!
!!!!  energyconf_trial=0.d0 !just to initialize this variable and make the compiler happy
!!!!  lstep=0.d0 !just to initialize this variable and make the compiler happy
!!!!  lstep_optimal=0.d0 !just to initialize this variable and make the compiler happy
!!!!  energyconf_der0=0.d0 !just to initialize this variable and make the compiler happy
!!!!
!!!!  time_convol=0.d0
!!!!  time_lincomb=0.d0
!!!!  time_commun=0.d0
!!!!  time_linalg=0.d0
!!!!  time_exponential=0.d0
!!!!  time_matrixmodification=0.d0
!!!!  time_matrixElements=0.d0
!!!!  t1_tot=mpi_wtime()
!!!!
!!!!  ! Initialize Umat
!!!!  do iorb=1,orbs%norb
!!!!      do jorb=1,orbs%norb
!!!!          if(jorb==iorb) then
!!!!              Umat(jorb,iorb)=1.d0
!!!!          else
!!!!              Umat(jorb,iorb)=0.d0
!!!!          end if
!!!!      end do
!!!!  end do
!!!!
!!!!  innerLoop: do it=1,nit
!!!!
!!!!  !write(*,*) '1: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
!!!!
!!!!      t1=mpi_wtime()
!!!!      !!call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
!!!!      !!     confdatarr, hx, lphi, -1, lvphi)
!!!!      call apply_position_operators(iproc, nproc, orbs, lzd, hx, hx, hx, lphi, lxphi, lyphi, lzphi)
!!!!
!!!!
!!!!      t2=mpi_wtime()
!!!!      time_convol=time_convol+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
!!!!      !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
!!!!      !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
!!!!      !deallocate(ttmat)
!!!!      !write(*,*) '2: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
!!!!      t2=mpi_wtime()
!!!!      time_commun=time_commun+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      !call getMatrixElements2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%op, lin%lb%comon, lphi, lvphi, lin%mad, Kmat)
!!!!      !call deallocateRecvBufferOrtho(comon, subname)
!!!!      !call deallocateSendBufferOrtho(comon, subname)
!!!!      !!call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, Kmat)
!!!!
!!!!      ! Build the matrices X, Y, Z
!!!!      call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lxphi, mad, X)
!!!!      call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lyphi, mad, Y)
!!!!      call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lzphi, mad, Z)
!!!!
!!!!      ! Build the matrices Xd, Yd, Zd
!!!!      do iorb=1,orbs%norb
!!!!          do jorb=1,orbs%norb
!!!!              if(jorb==iorb) then
!!!!                  Xd(jorb,iorb)=X(jorb,iorb)
!!!!                  Yd(jorb,iorb)=Y(jorb,iorb)
!!!!                  Zd(jorb,iorb)=Z(jorb,iorb)
!!!!              else
!!!!                  Xd(jorb,iorb)=0.d0
!!!!                  Yd(jorb,iorb)=0.d0
!!!!                  Zd(jorb,iorb)=0.d0
!!!!              end if
!!!!          end do
!!!!      end do
!!!!
!!!!      ! Build the matrices Xprime, Yprime, Zprime
!!!!      Xprime=X-Xd
!!!!      Yprime=Y-Yd
!!!!      Zprime=Z-Zd
!!!!
!!!!      ! Calculate value of Omega
!!!!      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Xprime(1,1), orbs%norb, &
!!!!           Xprime(1,1), orbs%norb, 0.d0, Xprimesquare, orbs%norb)
!!!!      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Yprime(1,1), orbs%norb, &
!!!!           Yprime(1,1), orbs%norb, 0.d0, Yprimesquare, orbs%norb)
!!!!      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Zprime(1,1), orbs%norb, &
!!!!           Zprime(1,1), orbs%norb, 0.d0, Zprimesquare, orbs%norb)
!!!!      omega=0.d0
!!!!      do iorb=1,orbs%norb
!!!!          omega=omega+Xprimesquare(iorb,iorb)+Yprimesquare(iorb,iorb)+Zprimesquare(iorb,iorb)
!!!!      end do
!!!!      if(iproc==0) write(*,'(a,i7,2es16.7)') 'it, omega, lstep', it, omega, lstep
!!!!
!!!!      call commutator(orbs%norb, Xprime, Xd, commutX)
!!!!      call commutator(orbs%norb, Yprime, Yd, commutY)
!!!!      call commutator(orbs%norb, Zprime, Zd, commutZ)
!!!!
!!!!      gmat=-2.d0*(commutX+commutY+commutZ)
!!!!
!!!!
!!!!
!!!!      !do iorb=1,orbs%norb
!!!!      !    do jorb=1,orbs%norb
!!!!      !        if(iproc==0) write(66,*) iorb,jorb,Kmat(jorb,iorb)
!!!!      !    end do
!!!!      !end do
!!!!      !call allocateSendBufferOrtho(comon, subname)
!!!!      !call allocateRecvBufferOrtho(comon, subname)
!!!!      !write(*,*) '3: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
!!!!      t2=mpi_wtime()
!!!!      time_matrixElements=time_matrixElements+t2-t1
!!!!
!!!!
!!!!      !!if(newgradient) then
!!!!      !!    call get_potential_matrices(iproc, nproc, at, orbs, lzd, op, comon, mad, rxyz, &
!!!!      !!         confdatarr, hx, lphi, potmat)
!!!!      !!    !call get_potential_matrices_new(iproc, nproc, lin, at, input, orbs, lzd, op, comon, rxyz, lphi, &
!!!!      !!    !     nlocregOnMPI, potmatsmall)
!!!!      !!end if
!!!!
!!!!
!!!!      !!if(.not.newgradient) then
!!!!          !energyconf_0=ddot(orbs%npsidim, lphi(1), 1, lvphi(1), 1)
!!!!          !call mpiallred(energyconf_0, 1, mpi_sum, mpi_comm_world, ierr)
!!!!          !!$$energyconf_0=0.d0
!!!!          !!$$do iorb=1,orbs%norb
!!!!          !!$$    energyconf_0 = energyconf_0 + Kmat(iorb,iorb)
!!!!          !!$$end do
!!!!      !!else
!!!!      !!    energyconf_0=0.d0
!!!!      !!    do iorb=1,orbs%norb
!!!!      !!        do jorb=1,orbs%norb
!!!!      !!            energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(jorb,iorb)
!!!!      !!            !energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(iorb,jorb)
!!!!      !!        end do
!!!!      !!    end do
!!!!      !!end if
!!!!      !!$$if(iproc==0) write(*,'(a,i6,3es20.10,2es17.7)') &
!!!!      !!$$             'it, energyconf_0, energyvonf_trial, energyconf_der0, lstep, lstep_optimal', &
!!!!      !!$$             it, energyconf_0, energyconf_trial, energyconf_der0, lstep, lstep_optimal
!!!!
!!!!      t1=mpi_wtime()
!!!!      !!if(.not.newgradient) then
!!!!          ! Construct antisymmtric matrix Gmat
!!!!          !!$$do iorb=1,orbs%norb
!!!!          !!$$    do jorb=1,orbs%norb
!!!!          !!$$        gmat(jorb,iorb)=2.d0*(Kmat(jorb,iorb)-Kmat(iorb,jorb))
!!!!          !!$$    end do
!!!!          !!$$end do 
!!!!      !!else
!!!!      !!    !!! THIS IS THE OLD VERSION #############################################################################
!!!!      !!    do iorb=1,orbs%norb
!!!!      !!        ilr=orbs%inwhichlocreg(iorb)
!!!!      !!        do jorb=1,orbs%norb
!!!!      !!            jlr=orbs%inwhichlocreg(jorb)
!!!!      !!            tt=0.d0
!!!!      !!            do lorb=1,orbs%norb
!!!!      !!                tt = tt + kernel(jorb,lorb)*Kmat(lorb,iorb) - kernel(iorb,lorb)*Kmat(lorb,jorb) + &
!!!!      !!                          kernel(jorb,lorb)*potmat(lorb,iorb,jlr) - kernel(iorb,lorb)*potmat(lorb,jorb,ilr)
!!!!      !!            end do
!!!!      !!            gmat(jorb,iorb)=-tt
!!!!      !!            !if(iproc==0) then
!!!!      !!            !    write(77,*) iorb, jorb, gmat(jorb,iorb)
!!!!      !!            !end if
!!!!      !!        end do
!!!!      !!    end do 
!!!!      !!    ! ########################################################################################################
!!!!      !!    !!! THIS IS THE NEW VERSION
!!!!      !!    !!gmat=0.d0
!!!!      !!    !!ii=0
!!!!      !!    !!ilrold=-1
!!!!      !!    !!do iorb=1,orbs%norbp
!!!!      !!    !!    iiorb=orbs%isorb+iorb
!!!!      !!    !!    ilr=orbs%inwhichlocreg(iiorb)
!!!!      !!    !!    if(ilr>ilrold) then
!!!!      !!    !!        ii=ii+1
!!!!      !!    !!    end if
!!!!      !!    !!    do jorb=1,orbs%norb
!!!!      !!    !!        jlr=orbs%inwhichlocreg(jorb)
!!!!      !!    !!        tt=0.d0
!!!!      !!    !!        do lorb=1,orbs%norb
!!!!      !!    !!            !tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
!!!!      !!    !!            !          - kernel(iiorb,lorb)*potmat(lorb,jorb,ilr)
!!!!      !!    !!            tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
!!!!      !!    !!                      - kernel(iiorb,lorb)*potmatsmall(lorb,jorb,ii)
!!!!      !!    !!        end do
!!!!      !!    !!        gmat(jorb,iiorb)=-tt
!!!!      !!    !!    end do
!!!!      !!    !!    ilrold=ilr
!!!!      !!    !!end do 
!!!!      !!    !!do iorb=1,orbs%norb
!!!!      !!    !!    ilr=orbs%inwhichlocreg(iorb)
!!!!      !!    !!    jlrold=-1
!!!!      !!    !!    jj=0
!!!!      !!    !!    do jorb=1,orbs%norbp
!!!!      !!    !!        jjorb=orbs%isorb+jorb
!!!!      !!    !!        jlr=orbs%inwhichlocreg(jjorb)
!!!!      !!    !!        if(jlr>jlrold) then
!!!!      !!    !!            jj=jj+1
!!!!      !!    !!        end if
!!!!      !!    !!        tt=0.d0
!!!!      !!    !!        do lorb=1,orbs%norb
!!!!      !!    !!            !tt = tt + kernel(jjorb,lorb)*potmat(lorb,iorb,jlr)
!!!!      !!    !!            tt = tt + kernel(jjorb,lorb)*potmatsmall(lorb,iorb,jj)
!!!!      !!    !!        end do
!!!!      !!    !!        gmat(jjorb,iorb)=gmat(jjorb,iorb)-tt
!!!!      !!    !!        jlrold=jlr
!!!!      !!    !!    end do
!!!!      !!    !!end do 
!!!!      !!    !!call mpiallred(gmat(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
!!!!      !!    !!do iorb=1,orbs%norb
!!!!      !!    !!    do jorb=1,orbs%norb
!!!!      !!    !!        if(iproc==0) then
!!!!      !!    !!            write(77,*) iorb, jorb, gmat(jorb,iorb)
!!!!      !!    !!        end if
!!!!      !!    !!    end do
!!!!      !!    !!end do
!!!!      !!end if
!!!!      t2=mpi_wtime()
!!!!      time_matrixmodification=time_matrixmodification+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      !Build the complex matrix -iGmat
!!!!      do iorb=1,orbs%norb
!!!!          do jorb=1,orbs%norb
!!!!              gmatc(jorb,iorb)=cmplx(0.d0,-gmat(jorb,iorb),kind=8)
!!!!              if(iproc==0) write(999,'(a,2i8,2es16.7)') 'iorb, jorb, gmatc(jorb,iorb)', iorb, jorb, gmatc(jorb,iorb)
!!!!          end do
!!!!      end do 
!!!!      t2=mpi_wtime()
!!!!      time_matrixmodification=time_matrixmodification+t2-t1
!!!!
!!!!
!!!!
!!!!      ! Diagonalize Gmatc
!!!!      t1=mpi_wtime()
!!!!      lwork=10*orbs%norb
!!!!      allocate(work(lwork), stat=istat) ! factor of 2 since it is assumed to be complex
!!!!      allocate(rwork(lwork), stat=istat)
!!!!      call zheev('v', 'l', orbs%norb, gmatc(1,1), orbs%norb, eval(1), work, lwork, rwork, info)
!!!!      if(info/=0) stop 'ERROR in zheev'
!!!!      deallocate(work)
!!!!      deallocate(rwork)
!!!!      t2=mpi_wtime()
!!!!      time_linalg=time_linalg+t2-t1
!!!!
!!!!
!!!!      !!$$! Calculate step size
!!!!      !!$$if(it==1) then
!!!!      !!$$    !!if(.not.newgradient) then
!!!!      !!$$        lstep=5.d-2/(maxval(eval))
!!!!      !!$$    !!else
!!!!      !!$$    !!    lstep=1.d-4/(maxval(eval))
!!!!      !!$$    !!end if
!!!!      !!$$else
!!!!      !!$$    lstep=2.d0*lstep_optimal
!!!!      !!$$    !lstep=1.d-3/(maxval(eval))
!!!!      !!$$end if
!!!!      lstep=1.d-2/(maxval(eval))
!!!!
!!!!      t1=mpi_wtime()
!!!!      ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
!!!!      ! This is also a diagonal matrix, so only calculate the diagonal part.
!!!!      do iorb=1,orbs%norb
!!!!         ttc=cmplx(0.d0,-lstep*eval(iorb),kind=8)
!!!!         expD_cmplx(iorb)=(0.d0,0.d0)
!!!!          do k=0,50
!!!!              expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
!!!!          end do
!!!!      end do
!!!!      t2=mpi_wtime()
!!!!      time_exponential=time_exponential+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      do iorb=1,orbs%norb
!!!!          do jorb=1,orbs%norb
!!!!              if(iorb==jorb) then
!!!!                  tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
!!!!              else
!!!!                  tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
!!!!              end if
!!!!          end do
!!!!      end do
!!!!      t2=mpi_wtime()
!!!!      time_matrixmodification=time_matrixmodification+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
!!!!           gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
!!!!      call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
!!!!           tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
!!!!      t2=mpi_wtime()
!!!!      time_linalg=time_linalg+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      ! Build new lphi
!!!!      do iorb=1,orbs%norb
!!!!          do jorb=1,orbs%norb
!!!!              tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
!!!!          end do
!!!!      end do
!!!!      t2=mpi_wtime()
!!!!      time_matrixmodification=time_matrixmodification+t2-t1
!!!!
!!!!      ! Update Umat
!!!!      call dcopy(orbs%norb**2, Umat(1,1), 1, tempmat3(1,1,2), 1)
!!!!      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
!!!!           tempmat3(1,1,1), orbs%norb, 0.d0, Umat(1,1), orbs%norb)
!!!!
!!!!      t1=mpi_wtime()
!!!!      !write(*,*) '5: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
!!!!      call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
!!!!           comon%recvbuf, tempmat3(1,1,1), .true., lphi)
!!!!      t2=mpi_wtime()
!!!!      time_lincomb=time_lincomb+t2-t1
!!!!
!!!!      !!$$t1=mpi_wtime()
!!!!      !!$$call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
!!!!      !!$$     confdatarr, hx, lphi, -1, lvphi)
!!!!      !!$$t2=mpi_wtime()
!!!!      !!$$time_convol=time_convol+t2-t1
!!!!
!!!!      !!$$!!if(.not.newgradient) then
!!!!      !!$$    energyconf_trial=ddot(max(orbs%npsidim_orbs,orbs%npsidim_comp), lphi(1), 1, lvphi(1), 1)
!!!!      !!$$    call mpiallred(energyconf_trial, 1, mpi_sum, mpi_comm_world, ierr)
!!!!      !!$$!!else
!!!!
!!!!      !!$$!!    call dcopy(comon%nrecvbuf, comon%recvbuf, 1, recvbuf, 1)
!!!!      !!$$!!    !call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
!!!!      !!$$!!    !call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
!!!!      !!$$!!    !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
!!!!      !!$$!!    !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
!!!!      !!$$!!    !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
!!!!      !!$$!!    !deallocate(ttmat)
!!!!      !!$$!!    call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, Kmat)
!!!!      !!$$!!    call dcopy(comon%nrecvbuf, recvbuf, 1, comon%recvbuf, 1)
!!!!
!!!!      !!$$!!    energyconf_trial=0.d0
!!!!      !!$$!!    do iorb=1,orbs%norb
!!!!      !!$$!!        do jorb=1,orbs%norb
!!!!      !!$$!!            energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(jorb,iorb)
!!!!      !!$$!!            !energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(iorb,jorb)
!!!!      !!$$!!        end do
!!!!      !!$$!!    end do
!!!!      !!$$!!end if
!!!!
!!!!      !!$$! Calculate the gradient of the confinement
!!!!      !!$$energyconf_der0=0.d0
!!!!      !!$$do iorb=1,orbs%norb
!!!!      !!$$    do jorb=1,orbs%norb
!!!!      !!$$        energyconf_der0=energyconf_der0+gmat(jorb,iorb)**2
!!!!      !!$$    end do
!!!!      !!$$end do
!!!!      !!$$energyconf_der0=-.5d0*energyconf_der0
!!!!
!!!!      !!$$! Calculate optimal step size
!!!!      !!$$lstep_optimal = -energyconf_der0*lstep**2/(2.d0*(energyconf_trial-energyconf_0-lstep*energyconf_der0))
!!!!      !!$$!!if(.not.newgradient) then
!!!!      !!$$    lstep_optimal=min(lstep_optimal,lstep)
!!!!      !!$$!!else
!!!!      !!$$!!    if(lstep_optimal<0) then
!!!!      !!$$!!        lstep_optimal=lstep
!!!!      !!$$!!    else
!!!!      !!$$!!        lstep_optimal=min(lstep_optimal,lstep)
!!!!      !!$$!!    end if
!!!!      !!$$!!end if
!!!!
!!!!      !!$$t1=mpi_wtime()
!!!!      !!$$! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
!!!!      !!$$! This is also a diagonal matrix, so only calculate the diagonal part.
!!!!      !!$$do iorb=1,orbs%norb
!!!!      !!$$   ttc=cmplx(0.d0,-lstep_optimal*eval(iorb),kind=8)
!!!!      !!$$   expD_cmplx(iorb)=(0.d0,0.d0)
!!!!      !!$$    do k=0,50
!!!!      !!$$        expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
!!!!      !!$$    end do
!!!!      !!$$end do
!!!!      !!$$t2=mpi_wtime()
!!!!      !!$$time_exponential=time_exponential+t2-t1
!!!!
!!!!      !!$$t1=mpi_wtime()
!!!!      !!$$do iorb=1,orbs%norb
!!!!      !!$$    do jorb=1,orbs%norb
!!!!      !!$$        if(iorb==jorb) then
!!!!      !!$$            tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
!!!!      !!$$        else
!!!!      !!$$            tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
!!!!      !!$$        end if
!!!!      !!$$    end do
!!!!      !!$$end do
!!!!      !!$$t2=mpi_wtime()
!!!!      !!$$time_matrixmodification=time_matrixmodification+t2-t1
!!!!
!!!!      !!$$t1=mpi_wtime()
!!!!      !!$$call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
!!!!      !!$$     gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
!!!!      !!$$call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
!!!!      !!$$     tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
!!!!      !!$$t2=mpi_wtime()
!!!!      !!$$time_linalg=time_linalg+t2-t1
!!!!
!!!!
!!!!      !!$$! Build new lphi
!!!!      !!$$do iorb=1,orbs%norb
!!!!      !!$$    do jorb=1,orbs%norb
!!!!      !!$$        tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
!!!!      !!$$    end do
!!!!      !!$$end do
!!!!      !!$$t1=mpi_wtime()
!!!!      !!$$call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
!!!!      !!$$     comon%recvbuf, tempmat3(1,1,1), .true., lphi)
!!!!      !!$$t2=mpi_wtime()
!!!!      !!$$time_lincomb=time_lincomb+t2-t1
!!!!
!!!!
!!!!      !!$$!if(it<nit) then
!!!!      !!$$!    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
!!!!      !!$$!    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
!!!!      !!$$!end if
!!!!
!!!!
!!!!  end do innerLoop
!!!!
!!!!  t2_tot=mpi_wtime()
!!!!  time_tot=t2_tot-t1_tot
!!!!  call mpiallred(time_convol, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_commun, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_lincomb, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_linalg, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_matrixmodification, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_exponential, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_matrixelements, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_tot, 1, mpi_max, mpi_comm_world, ierr)
!!!!  !time_convol=time_convol/dble(nproc)
!!!!  !time_commun=time_commun/dble(nproc)
!!!!  !time_lincomb=time_lincomb/dble(nproc)
!!!!  !time_linalg=time_linalg/dble(nproc)
!!!!  !time_matrixmodification=time_matrixmodification/dble(nproc)
!!!!  !time_exponential_=time_exponential/dble(nproc)
!!!!  !time_tot=time_tot/dble(nproc)
!!!!  !!if(iproc==0) then
!!!!  !!    write(*,'(a,es16.6)') 'total time: ',time_tot
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'convolutions: ',time_convol,' (',time_convol/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear combinations: ',time_lincomb,' (',time_lincomb/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'communication: ',time_commun,' (',time_commun/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear algebra: ',time_linalg,' (',time_linalg/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix modification: ',time_matrixmodification, &
!!!!  !!                                   ' (',time_matrixmodification/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'building exponential: ',time_exponential,' (',time_exponential/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix elements ',time_matrixelements,' (',time_matrixelements/time_tot*100.d0,'%)'
!!!!  !!end if
!!!!
!!!!
!!!!  iall=-product(shape(gmat))*kind(gmat)
!!!!  deallocate(gmat, stat=istat)
!!!!  call memocc(istat, iall, 'gmat', subname)
!!!!  iall=-product(shape(gmatc))*kind(gmatc)
!!!!  deallocate(gmatc, stat=istat)
!!!!  !call memocc(istat, iall, 'gmatc', subname)
!!!!  iall=-product(shape(omatc))*kind(omatc)
!!!!  deallocate(omatc, stat=istat)
!!!!  !call memocc(istat, iall, 'omatc', subname)
!!!!  iall=-product(shape(tempmat3))*kind(tempmat3)
!!!!  deallocate(tempmat3, stat=istat)
!!!!  call memocc(istat, iall, 'tempmat3', subname)
!!!!  iall=-product(shape(eval))*kind(eval)
!!!!  deallocate(eval, stat=istat)
!!!!  call memocc(istat, iall, 'eval', subname)
!!!!  iall=-product(shape(expD_cmplx))*kind(expD_cmplx)
!!!!  deallocate(expD_cmplx, stat=istat)
!!!!  call memocc(istat, iall, 'expD_cmplx', subname)
!!!!  iall=-product(shape(tempmatc))*kind(tempmatc)
!!!!  deallocate(tempmatc, stat=istat)
!!!!  !call memocc(istat, iall, 'tempmatc', subname)
!!!!  iall=-product(shape(hamtrans))*kind(hamtrans)
!!!!  deallocate(hamtrans, stat=istat)
!!!!  call memocc(istat, iall, 'hamtrans', subname)
!!!!  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
!!!!  deallocate(lphiovrlp, stat=istat)
!!!!  call memocc(istat, iall, 'lphiovrlp', subname)
!!!!  iall=-product(shape(lvphi))*kind(lvphi)
!!!!  deallocate(lvphi, stat=istat)
!!!!  call memocc(istat, iall, 'lvphi', subname)
!!!!  iall=-product(shape(Kmat))*kind(Kmat)
!!!!  deallocate(Kmat, stat=istat)
!!!!  call memocc(istat, iall, 'Kmat', subname)
!!!!  iall=-product(shape(recvbuf))*kind(recvbuf)
!!!!  deallocate(recvbuf, stat=istat)
!!!!  call memocc(istat, iall, 'recvbuf', subname)
!!!!
!!!!  iall=-product(shape(potmat))*kind(potmat)
!!!!  deallocate(potmat, stat=istat)
!!!!  call memocc(istat, iall, 'potmat', subname)
!!!!  iall=-product(shape(potmatsmall))*kind(potmatsmall)
!!!!  deallocate(potmatsmall, stat=istat)
!!!!  call memocc(istat, iall, 'potmatsmall', subname)
!!!!
!!!!
!!!!  iall=-product(shape(lxphi))*kind(lxphi)
!!!!  deallocate(lxphi, stat=istat)
!!!!  call memocc(istat, iall, 'lxphi', subname)
!!!!  iall=-product(shape(lyphi))*kind(lyphi)
!!!!  deallocate(lyphi, stat=istat)
!!!!  call memocc(istat, iall, 'lyphi', subname)
!!!!  iall=-product(shape(lzphi))*kind(lzphi)
!!!!  deallocate(lzphi, stat=istat)
!!!!  call memocc(istat, iall, 'lzphi', subname)
!!!!
!!!!
!!!!
!!!!  iall=-product(shape(X))*kind(X)
!!!!  deallocate(X, stat=istat)
!!!!  call memocc(istat, iall, 'X', subname)
!!!!  iall=-product(shape(Y))*kind(Y)
!!!!  deallocate(Y, stat=istat)
!!!!  call memocc(istat, iall, 'Y', subname)
!!!!  iall=-product(shape(Z))*kind(Z)
!!!!  deallocate(Z, stat=istat)
!!!!  call memocc(istat, iall, 'Z', subname)
!!!!  iall=-product(shape(Xd))*kind(Xd)
!!!!  deallocate(Xd, stat=istat)
!!!!  call memocc(istat, iall, 'Xd', subname)
!!!!  iall=-product(shape(Yd))*kind(Yd)
!!!!  deallocate(Yd, stat=istat)
!!!!  call memocc(istat, iall, 'Yd', subname)
!!!!  iall=-product(shape(Zd))*kind(Zd)
!!!!  deallocate(Zd, stat=istat)
!!!!  call memocc(istat, iall, 'Zd', subname)
!!!!  iall=-product(shape(Xprime))*kind(Xprime)
!!!!  deallocate(Xprime, stat=istat)
!!!!  call memocc(istat, iall, 'Xprime', subname)
!!!!  iall=-product(shape(Yprime))*kind(Yprime)
!!!!  deallocate(Yprime, stat=istat)
!!!!  call memocc(istat, iall, 'Yprime', subname)
!!!!  iall=-product(shape(Zprime))*kind(Zprime)
!!!!  deallocate(Zprime, stat=istat)
!!!!  call memocc(istat, iall, 'Zprime', subname)
!!!!  iall=-product(shape(Xprimesquare))*kind(Xprimesquare)
!!!!  deallocate(Xprimesquare, stat=istat)
!!!!  call memocc(istat, iall, 'Xprimesquare', subname)
!!!!  iall=-product(shape(Yprimesquare))*kind(Yprimesquare)
!!!!  deallocate(Yprimesquare, stat=istat)
!!!!  call memocc(istat, iall, 'Yprimesquare', subname)
!!!!  iall=-product(shape(Zprimesquare))*kind(Zprimesquare)
!!!!  deallocate(Zprimesquare, stat=istat)
!!!!  call memocc(istat, iall, 'Zprimesquare', subname)
!!!!  iall=-product(shape(commutX))*kind(commutX)
!!!!  deallocate(commutX, stat=istat)
!!!!  call memocc(istat, iall, 'commutX', subname)
!!!!  iall=-product(shape(commutY))*kind(commutY)
!!!!  deallocate(commutY, stat=istat)
!!!!  call memocc(istat, iall, 'commutY', subname)
!!!!  iall=-product(shape(commutZ))*kind(commutZ)
!!!!  deallocate(commutZ, stat=istat)
!!!!  call memocc(istat, iall, 'commutZ', subname)
!!!!
!!!!
!!!!
!!!!  call deallocateRecvBufferOrtho(comon, subname)
!!!!  call deallocateSendBufferOrtho(comon, subname)
!!!!
!!!!end subroutine MLWF





subroutine MLWFnew(iproc, nproc, lzd, orbs, at, op, comon, mad, rxyz, nit, kernel, &
           newgradient, confdatarr, hx, locregCenters, maxDispl, lphi, Umat, centers)
use module_base
use module_types
use module_interfaces, exceptThisOne => MLWFnew
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nit
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(atoms_data),intent(in):: at
type(overlapParameters),intent(inout):: op
type(p2pComms),intent(inout):: comon
type(matrixDescriptors),intent(in):: mad
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(orbs%norb,orbs%norb),intent(in):: kernel
logical,intent(in):: newgradient
real(8),intent(in):: hx, maxDispl
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(3,lzd%nlr),intent(in):: locregCenters
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: lphi
real(8),dimension(orbs%norb,orbs%norb),intent(out):: Umat
real(8),dimension(3,lzd%nlr),intent(out):: centers

! Local variables
integer:: it, info, lwork, k, istat, iorb, jorb, iall, ierr, ist, jst, ilrold, ncount, jjorb, iiorb, ilr, lorb, jlr
integer:: nlocregOnMPI, jlrold, jj, ii, ncnt, order
real(8):: trace, lstep, dfactorial, energyconf_trial, energyconf_0, energyconf_der0, lstep_optimal, ddot
real(8):: tt1, tt2, tt3, tt4, tt5, tt, var
real(8):: t1, t2, t1_tot, t2_tot, omega
real(8):: time_convol, time_commun, time_lincomb, time_linalg, time_matrixmodification, time_exponential, time_tot
real(8):: time_matrixelements, rspread, locdiff
complex(8):: ttc
real(8),dimension(:,:),allocatable:: gmat, hamtrans, ttmat, Kmat, X, Y, Z, Xd, Yd, Zd, X2, Y2, Z2, R, R2
real(8),dimension(:,:),allocatable:: commutX1, commutY1, commutZ1, commutX2, commutY2, commutZ2, commutX3, commutY3, commutZ3
real(8),dimension(:,:),allocatable:: Xprime, Yprime, Zprime, Xprimesquare, Yprimesquare, Zprimesquare
real(8),dimension(:,:),allocatable:: X0, H, gHmat, Y0, M, gMmat, Z0, N, gNmat, centers_start, centers_end
real(8),dimension(:,:),allocatable:: HXd, HXdsquare, HXdX0, HX0, HX0square, Xdsquare, XdX0, X0square
real(8),dimension(:,:),allocatable:: MYd, MYdsquare, MYdY0, MY0, MY0square, Ydsquare, YdY0, Y0square
real(8),dimension(:,:),allocatable:: NZd, NZdsquare, NZdZ0, NZ0, NZ0square, Zdsquare, ZdZ0, Z0square
real(8),dimension(:,:,:),allocatable:: potmat, potmatsmall
complex(8),dimension(:,:),allocatable:: gmatc, omatc
complex(8),dimension(:,:,:),allocatable:: tempmatc
complex(8),dimension(:),allocatable:: work, expD_cmplx
real(8),dimension(:),allocatable:: rwork, eval, lphiovrlp, lvphi, lxphi, lyphi, lzphi, normarr
real(8),dimension(:,:,:),allocatable:: tempmat3
character(len=*),parameter:: subname='MLWFnew'
type(p2pComms):: comon_local

! Quick return if possible. In this way the localization regions will remain unchanged.
if(nit==-1) return

allocate(gmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, gmat, 'gmat', subname)
allocate(hamtrans(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, hamtrans, 'hamtrans', subname)
allocate(gmatc(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, gmatc, 'gmatc', subname)
allocate(omatc(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, omatc, 'omatc', subname)
allocate(tempmat3(orbs%norb,orbs%norb,3), stat=istat)
call memocc(istat, tempmat3, 'tempmat3', subname)
allocate(eval(orbs%norb), stat=istat)
call memocc(istat, eval, 'eval', subname)
allocate(expD_cmplx(orbs%norb), stat=istat)
call memocc(istat, expD_cmplx, 'expD_cmplx', subname)
allocate(tempmatc(orbs%norb,orbs%norb,2), stat=istat)
!call memocc(istat, tempmatc, 'tempmatc', subname)
allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
allocate(lvphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, lvphi, 'lvphi', subname)
allocate(Kmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Kmat, 'Kmat', subname)
!!allocate(recvbuf(comon%nrecvbuf), stat=istat)
!!call memocc(istat, recvbuf, 'recvbuf', subname)

allocate(potmat(orbs%norb,orbs%norb,at%nat), stat=istat)
call memocc(istat, potmat, 'potmat', subname)

allocate(lxphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, lxphi, 'lxphi', subname)
allocate(lyphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, lyphi, 'lyphi', subname)
allocate(lzphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, lzphi, 'lzphi', subname)

allocate(X(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, X, 'X', subname)
allocate(Y(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Y, 'Y', subname)
allocate(Z(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Z, 'Z', subname)
allocate(X2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, X2, 'X2', subname)
allocate(Y2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Y2, 'Y2', subname)
allocate(Z2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Z2, 'Z2', subname)
allocate(R(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, R, 'R', subname)
allocate(R2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, R2, 'R2', subname)
allocate(Xd(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Xd, 'Xd', subname)
allocate(Yd(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Yd, 'Yd', subname)
allocate(Zd(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Zd, 'Zd', subname)
allocate(Xprime(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Xprime, 'Xprime', subname)
allocate(Yprime(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Yprime, 'Yprime', subname)
allocate(Zprime(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Zprime, 'Zprime', subname)
allocate(Xprimesquare(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Xprimesquare, 'Xprimesquare', subname)
allocate(Yprimesquare(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Yprimesquare, 'Yprimesquare', subname)
allocate(Zprimesquare(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Zprimesquare, 'Zprimesquare', subname)

allocate(commutX1(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutX1, 'commutX1', subname)
allocate(commutY1(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutY1, 'commutY1', subname)
allocate(commutZ1(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutZ1, 'commutZ1', subname)
allocate(commutX2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutX2, 'commutX2', subname)
allocate(commutY2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutY2, 'commutY2', subname)
allocate(commutZ2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutZ2, 'commutZ2', subname)
allocate(commutX3(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutX3, 'commutX3', subname)
allocate(commutY3(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutY3, 'commutY3', subname)
allocate(commutZ3(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutZ3, 'commutZ3', subname)

!!allocate(X0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, X0, 'X0', subname)
!!allocate(H(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, H, 'H', subname)
!!allocate(gHmat(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, gHmat, 'gHmat', subname)
!!
!!allocate(Y0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Y0, 'Y0', subname)
!!allocate(M(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, M, 'M', subname)
!!allocate(gMmat(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, gMmat, 'gMmat', subname)
!!
!!allocate(Z0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Z0, 'Z0', subname)
!!allocate(N(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, N, 'N', subname)
!!allocate(gNmat(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, gNmat, 'gNmat', subname)
!!
!!allocate(HXd(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, HXd, 'HXd', subname)
!!allocate(HXdsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, HXdsquare, 'HXdsquare', subname)
!!allocate(HXdX0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, HXdX0, 'HXdX0', subname)
!!allocate(HX0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, HX0, 'HX0', subname)
!!allocate(HX0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, HX0square, 'HX0square', subname)
!!allocate(Xdsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Xdsquare, 'Xdsquare', subname)
!!allocate(XdX0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, XdX0, 'XdX0', subname)
!!allocate(X0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, X0square, 'X0square', subname)
!!
!!allocate(MYd(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, MYd, 'MYd', subname)
!!allocate(MYdsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, MYdsquare, 'MYdsquare', subname)
!!allocate(MYdY0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, MYdY0, 'MYdY0', subname)
!!allocate(MY0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, MY0, 'MY0', subname)
!!allocate(MY0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, MY0square, 'MY0square', subname)
!!allocate(Ydsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Ydsquare, 'Ydsquare', subname)
!!allocate(YdY0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, YdY0, 'YdY0', subname)
!!allocate(Y0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Y0square, 'Y0square', subname)
!!
!!allocate(NZd(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, NZd, 'NZd', subname)
!!allocate(NZdsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, NZdsquare, 'NZdsquare', subname)
!!allocate(NZdZ0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, NZdZ0, 'NZdZ0', subname)
!!allocate(NZ0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, NZ0, 'NZ0', subname)
!!allocate(NZ0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, NZ0square, 'NZ0square', subname)
!!allocate(Zdsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Zdsquare, 'Zdsquare', subname)
!!allocate(ZdZ0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, ZdZ0, 'ZdZ0', subname)
!!allocate(Z0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Z0square, 'Z0square', subname)

allocate(centers_start(3,lzd%nlr), stat=istat)
call memocc(istat, centers_start, 'centers_start', subname)
allocate(centers_end(3,lzd%nlr), stat=istat)
call memocc(istat, centers_end, 'centers_end', subname)
allocate(normarr(orbs%norb), stat=istat)
call memocc(istat, normarr, 'normarr', subname)


! Count how many locregs each process handles
ilrold=-1
nlocregOnMPI=0
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inwhichlocreg(iiorb)
    !if(ilr>ilrold) then
    if(ilr/=ilrold) then
        nlocregOnMPI=nlocregOnMPI+1
    end if
    ilrold=ilr
end do
allocate(potmatsmall(orbs%norb,orbs%norb,nlocregOnMPI), stat=istat)
call memocc(istat, potmatsmall, 'potmatsmall', subname)



  call allocateSendBufferOrtho(comon, subname)
  call allocateRecvBufferOrtho(comon, subname)
  ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
  !if(nit>0) then
  !    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
  !    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
  !end if

  energyconf_trial=0.d0 !just to initialize this variable and make the compiler happy
  lstep=0.d0 !just to initialize this variable and make the compiler happy
  lstep_optimal=0.d0 !just to initialize this variable and make the compiler happy
  energyconf_der0=0.d0 !just to initialize this variable and make the compiler happy

  time_convol=0.d0
  time_lincomb=0.d0
  time_commun=0.d0
  time_linalg=0.d0
  time_exponential=0.d0
  time_matrixmodification=0.d0
  time_matrixElements=0.d0
  t1_tot=mpi_wtime()

  ! Initialize Umat
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          if(jorb==iorb) then
              Umat(jorb,iorb)=1.d0
          else
              Umat(jorb,iorb)=0.d0
          end if
      end do
  end do

  !!! Initialize Lagragme multipliers to zero
  !!H=0.d0
  !!M=0.d0
  !!N=0.d0

  !!! Initialize X0 etc.
  !!X0=0.d0
  !!Y0=0.d0
  !!Z0=0.d0
  !!do iorb=1,orbs%norbp
  !!    iiorb=orbs%isorb+iorb
  !!    X0(iiorb,iiorb)=confdatarr(iorb)%rxyzConf(1)
  !!    Y0(iiorb,iiorb)=confdatarr(iorb)%rxyzConf(2)
  !!    Z0(iiorb,iiorb)=confdatarr(iorb)%rxyzConf(3)
  !!end do
  !!call mpiallred(X0(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
  !!call mpiallred(Y0(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
  !!call mpiallred(Z0(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)


  ! Calculate the norm of all basis functions
  normarr=0.d0
  ist=1
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      ncnt=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      normarr(iiorb)=ddot(ncnt, lphi(ist), 1, lphi(ist), 1)
      ist=ist+ncnt
  end do
  call mpiallred(normarr(1), orbs%norb, mpi_sum, mpi_comm_world, ierr)
  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        write(*,'(a,i6,es16.6)') 'new norm before loop: iorb, normarr(iorb)', iorb, normarr(iorb)
  !!    end do
  !!end if


  !!! OTHER CHECK
  !!call apply_rminusmu_operator(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, locregCenters, lvphi)
  !!call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R2)
  !!! New calculation of the spread
  !!var=0.d0
  !!do iorb=1,orbs%norb
  !!    tt = R2(iorb,iorb)/normarr(iorb)
  !!    var=var+tt
  !!end do
  !!!!if(iproc==0) write(*,'(a,es18.8)') 'FIRST: total variance', var

  innerLoop: do it=1,nit

  !write(*,*) '1: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)

      t1=mpi_wtime()
      !!call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
      !!     confdatarr, hx, lphi, -1, lvphi)

      if(it==1) then
          order=1
          call apply_position_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lxphi, lyphi, lzphi)
          ! Build the matrices X, Y, Z
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lxphi, mad, X)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lyphi, mad, Y)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lzphi, mad, Z)

          order=2
          call apply_position_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lxphi, lyphi, lzphi)
          ! Build the matrices X, Y, Z
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lxphi, mad, X2)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lyphi, mad, Y2)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lzphi, mad, Z2)

          order=1
          call apply_r_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lvphi)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R)

          order=2
          call apply_r_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lvphi)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R2)


          do iorb=1,orbs%norb
              ilr=orbs%inwhichlocreg(iorb)
              !!if(iproc==0) write(*,'(a,2i5,3f10.4,4x,3f10.4)') 'START: iorb, ilr, centers: ', &
              !!  iorb, ilr, X(iorb,iorb)/normarr(iorb), Y(iorb,iorb)/normarr(iorb), Z(iorb,iorb)/normarr(iorb), &
              !!  locregCenters(1,ilr),  locregCenters(2,ilr), locregCenters(3,ilr)
              centers_start(1,iorb)=X(iorb,iorb)/normarr(iorb)
              centers_start(2,iorb)=Y(iorb,iorb)/normarr(iorb)
              centers_start(3,iorb)=Z(iorb,iorb)/normarr(iorb)
          end do
      else
        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             X(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, X(1,1), orbs%norb)

        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             Y(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, Y(1,1), orbs%norb)

        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             Z(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, Z(1,1), orbs%norb)


        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             X2(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, X2(1,1), orbs%norb)

        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             Y2(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, Y2(1,1), orbs%norb)

        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             Z2(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, Z2(1,1), orbs%norb)


        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             R(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, R(1,1), orbs%norb)

        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             R2(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, R2(1,1), orbs%norb)
      end if

      ! New calculation of the spread
      !tt=0.d0
      var=0.d0
      do iorb=1,orbs%norb
          !!if(iproc==0) write(*,'(a,i7,9es15.6)') 'iorb, diagonal part of X,Y,Z,X2,Y2,Z2,R,R2,normarr', &
          !!    iorb, X(iorb,iorb), Y(iorb,iorb), Z(iorb,iorb), X2(iorb,iorb), Y2(iorb,iorb), Z2(iorb,iorb), &
          !!    R(iorb,iorb), R2(iorb,iorb), normarr(iorb)
          !!tt = X2(iorb,iorb)/normarr(iorb) + Y2(iorb,iorb)/normarr(iorb) + Z2(iorb,iorb)/normarr(iorb) &
          !!    -X(iorb,iorb)**2/normarr(iorb) - Y(iorb,iorb)**2/normarr(iorb) - Z(iorb,iorb)**2/normarr(iorb)
          !tt = R2(iorb,iorb)/normarr(iorb)-R(iorb,iorb)**2/normarr(iorb)
          tt = R2(iorb,iorb)/normarr(iorb)-(R(iorb,iorb)/normarr(iorb))**2
          var=var+tt
          !if(iproc==0) write(*,'(a,i8,es15.6)') 'iorb, tt', iorb, tt
      end do
      !if(iproc==0) write(*,'(a,i6,es18.8)') 'it, total variance', it, var
      !if(iproc==0) write(*,'(a,es18.8)') 'NEW SPREAD: ',tt

      ! Check the deviation from the center of the original localization region
      !do ilr=1,lzd%nlr
      do iorb=1,orbs%norb
          ilr=orbs%inwhichlocreg(iorb)
          tt = (locregCenters(1,ilr)-X(iorb,iorb)/normarr(iorb))**2 &
              +(locregCenters(2,ilr)-Y(iorb,iorb)/normarr(iorb))**2 &
              +(locregCenters(3,ilr)-Z(iorb,iorb)/normarr(iorb))**2 
          tt=sqrt(tt)
          if(tt>maxDispl) then
              if(iproc==0) write(*,'(a,i0,2x,es11.2)') 'WARNING: too large displacement for locreg ',ilr,tt
          end if
      end do

      !!h=0.d0
      !!m=0.d0
      !!n=0.d0
      !!y=0.d0
      !!z=0.d0
      !!y0=0.d0
      !!z0=0.d0

      !!tt=0.d0
      !!do iorb=1,orbs%norb
      !!    tt = tt + (X(iorb,iorb)-X0(iorb,iorb))**2 + (Y(iorb,iorb)-Y0(iorb,iorb))**2 + (Z(iorb,iorb)-Z0(iorb,iorb))**2
      !!end do
      !!if(iproc==0) write(*,'(a,es16.7)') 'difference to locreg center', tt

      ! Build the matrices Xd, Yd, Zd
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              if(jorb==iorb) then
                  Xd(jorb,iorb)=X(jorb,iorb)
                  Yd(jorb,iorb)=Y(jorb,iorb)
                  Zd(jorb,iorb)=Z(jorb,iorb)
              else
                  Xd(jorb,iorb)=0.d0
                  Yd(jorb,iorb)=0.d0
                  Zd(jorb,iorb)=0.d0
              end if
          end do
      end do

      ! Build the matrices Xprime, Yprime, Zprime
      Xprime=X-Xd
      Yprime=Y-Yd
      Zprime=Z-Zd

      ! Calculate value of Omega
      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Xprime(1,1), orbs%norb, &
           Xprime(1,1), orbs%norb, 0.d0, Xprimesquare, orbs%norb)
      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Yprime(1,1), orbs%norb, &
           Yprime(1,1), orbs%norb, 0.d0, Yprimesquare, orbs%norb)
      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Zprime(1,1), orbs%norb, &
           Zprime(1,1), orbs%norb, 0.d0, Zprimesquare, orbs%norb)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, H(1,1), orbs%norb, &
      !!     Xd(1,1), orbs%norb, 0.d0, HXd, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, M(1,1), orbs%norb, &
      !!     Yd(1,1), orbs%norb, 0.d0, MYd, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, N(1,1), orbs%norb, &
      !!     Zd(1,1), orbs%norb, 0.d0, NZd, orbs%norb)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, HXd(1,1), orbs%norb, &
      !!     Xd(1,1), orbs%norb, 0.d0, HXdsquare, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, MYd(1,1), orbs%norb, &
      !!     Yd(1,1), orbs%norb, 0.d0, MYdsquare, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, NZd(1,1), orbs%norb, &
      !!     Zd(1,1), orbs%norb, 0.d0, NZdsquare, orbs%norb)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, HXd(1,1), orbs%norb, &
      !!     X0(1,1), orbs%norb, 0.d0, HXdX0, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, MYd(1,1), orbs%norb, &
      !!     Y0(1,1), orbs%norb, 0.d0, MYdY0, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, NZd(1,1), orbs%norb, &
      !!     Z0(1,1), orbs%norb, 0.d0, NZdZ0, orbs%norb)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, H(1,1), orbs%norb, &
      !!     X0(1,1), orbs%norb, 0.d0, HX0, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, M(1,1), orbs%norb, &
      !!     Y0(1,1), orbs%norb, 0.d0, MY0, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, N(1,1), orbs%norb, &
      !!     Z0(1,1), orbs%norb, 0.d0, NZ0, orbs%norb)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, HX0(1,1), orbs%norb, &
      !!     X0(1,1), orbs%norb, 0.d0, HX0square, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, MY0(1,1), orbs%norb, &
      !!     Y0(1,1), orbs%norb, 0.d0, MY0square, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, NZ0(1,1), orbs%norb, &
      !!     Z0(1,1), orbs%norb, 0.d0, NZ0square, orbs%norb)
      omega=0.d0
      rspread=0.d0
      do iorb=1,orbs%norb
          !!omega = omega + Xprimesquare(iorb,iorb)+Yprimesquare(iorb,iorb)+Zprimesquare(iorb,iorb) &
          !!              - (HXdsquare(iorb,iorb)-2.d0*HXdX0(iorb,iorb)+HX0square(iorb,iorb)) &
          !!              - (MYdsquare(iorb,iorb)-2.d0*MYdY0(iorb,iorb)+MY0square(iorb,iorb)) &
          !!              - (NZdsquare(iorb,iorb)-2.d0*NZdZ0(iorb,iorb)+NZ0square(iorb,iorb))
          omega = omega + Xprimesquare(iorb,iorb)+Yprimesquare(iorb,iorb)+Zprimesquare(iorb,iorb) 
          rspread = rspread + Xprimesquare(iorb,iorb)+Yprimesquare(iorb,iorb)+Zprimesquare(iorb,iorb)
      end do
      !if(iproc==0) write(*,'(a,i7,2es16.7)') 'it, omega, lstep', it, omega, lstep

      call commutator(orbs%norb, Xprime, Xd, commutX1)
      call commutator(orbs%norb, Yprime, Yd, commutY1)
      call commutator(orbs%norb, Zprime, Zd, commutZ1)

      !!call commutator(orbs%norb, HXd, X, commutX2)
      !!call commutator(orbs%norb, MYd, Y, commutY2)
      !!call commutator(orbs%norb, NZd, Z, commutZ2)

      !!call commutator(orbs%norb, HX0, X, commutX3)
      !!call commutator(orbs%norb, MY0, Y, commutY3)
      !!call commutator(orbs%norb, NZ0, Z, commutZ3)

      !!gmat=2.d0*( (commutX1+commutY1+commutZ1) &
      !!             -commutX2+commutX3 &
      !!             -commutY2+commutY3 &
      !!             -commutZ2+commutZ3 )
      gmat=2.d0*( commutX1+commutY1+commutZ1 )


      !!! Calculate gradient for Lagrange multipliers
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Xd(1,1), orbs%norb, &
      !!     Xd(1,1), orbs%norb, 0.d0, Xdsquare(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Xd(1,1), orbs%norb, &
      !!     X0(1,1), orbs%norb, 0.d0, XdX0(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, X0(1,1), orbs%norb, &
      !!     X0(1,1), orbs%norb, 0.d0, X0square(1,1), orbs%norb)
      !!gHmat=Xdsquare-2.d0*XdX0+X0square
      !!H=H-1.d-3*gHmat

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Yd(1,1), orbs%norb, &
      !!     Yd(1,1), orbs%norb, 0.d0, Ydsquare(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Yd(1,1), orbs%norb, &
      !!     Y0(1,1), orbs%norb, 0.d0, YdY0(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Y0(1,1), orbs%norb, &
      !!     Y0(1,1), orbs%norb, 0.d0, Y0square(1,1), orbs%norb)
      !!gMmat=Ydsquare-2.d0*YdY0+Y0square
      !!M=M-1.d-3*gMmat

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Zd(1,1), orbs%norb, &
      !!     Zd(1,1), orbs%norb, 0.d0, Zdsquare(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Zd(1,1), orbs%norb, &
      !!     Z0(1,1), orbs%norb, 0.d0, ZdZ0(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Z0(1,1), orbs%norb, &
      !!     Z0(1,1), orbs%norb, 0.d0, Z0square(1,1), orbs%norb)
      !!gNmat=Zdsquare-2.d0*ZdZ0+Z0square
      !!N=N-1.d-3*gNmat

      locdiff=0.d0
      !!do iorb=1,orbs%norb
      !!    locdiff = locdiff + gHmat(iorb,iorb) + gMmat(iorb,iorb) + gNmat(iorb,iorb)
      !!end do
      if(iproc==0 .and. verbose >2) write(*,'(a,i8,3es16.7,es12.4)')&
           'it, rspread, locdiff, omega, lstep', it, rspread, locdiff, omega, lstep



      !do iorb=1,orbs%norb
      !    do jorb=1,orbs%norb
      !        if(iproc==0) write(66,*) iorb,jorb,Kmat(jorb,iorb)
      !    end do
      !end do
      !call allocateSendBufferOrtho(comon, subname)
      !call allocateRecvBufferOrtho(comon, subname)
      !write(*,*) '3: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      t2=mpi_wtime()
      time_matrixElements=time_matrixElements+t2-t1


      !!if(newgradient) then
      !!    call get_potential_matrices(iproc, nproc, at, orbs, lzd, op, comon, mad, rxyz, &
      !!         confdatarr, hx, lphi, potmat)
      !!    !call get_potential_matrices_new(iproc, nproc, lin, at, input, orbs, lzd, op, comon, rxyz, lphi, &
      !!    !     nlocregOnMPI, potmatsmall)
      !!end if


      !!if(.not.newgradient) then
          !energyconf_0=ddot(orbs%npsidim, lphi(1), 1, lvphi(1), 1)
          !call mpiallred(energyconf_0, 1, mpi_sum, mpi_comm_world, ierr)
          !!$$energyconf_0=0.d0
          !!$$do iorb=1,orbs%norb
          !!$$    energyconf_0 = energyconf_0 + Kmat(iorb,iorb)
          !!$$end do
      !!else
      !!    energyconf_0=0.d0
      !!    do iorb=1,orbs%norb
      !!        do jorb=1,orbs%norb
      !!            energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(jorb,iorb)
      !!            !energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(iorb,jorb)
      !!        end do
      !!    end do
      !!end if
      !!$$if(iproc==0) write(*,'(a,i6,3es20.10,2es17.7)') &
      !!$$             'it, energyconf_0, energyvonf_trial, energyconf_der0, lstep, lstep_optimal', &
      !!$$             it, energyconf_0, energyconf_trial, energyconf_der0, lstep, lstep_optimal

      t1=mpi_wtime()
      !!if(.not.newgradient) then
          ! Construct antisymmtric matrix Gmat
          !!$$do iorb=1,orbs%norb
          !!$$    do jorb=1,orbs%norb
          !!$$        gmat(jorb,iorb)=2.d0*(Kmat(jorb,iorb)-Kmat(iorb,jorb))
          !!$$    end do
          !!$$end do 
      !!else
      !!    !!! THIS IS THE OLD VERSION #############################################################################
      !!    do iorb=1,orbs%norb
      !!        ilr=orbs%inwhichlocreg(iorb)
      !!        do jorb=1,orbs%norb
      !!            jlr=orbs%inwhichlocreg(jorb)
      !!            tt=0.d0
      !!            do lorb=1,orbs%norb
      !!                tt = tt + kernel(jorb,lorb)*Kmat(lorb,iorb) - kernel(iorb,lorb)*Kmat(lorb,jorb) + &
      !!                          kernel(jorb,lorb)*potmat(lorb,iorb,jlr) - kernel(iorb,lorb)*potmat(lorb,jorb,ilr)
      !!            end do
      !!            gmat(jorb,iorb)=-tt
      !!            !if(iproc==0) then
      !!            !    write(77,*) iorb, jorb, gmat(jorb,iorb)
      !!            !end if
      !!        end do
      !!    end do 
      !!    ! ########################################################################################################
      !!    !!! THIS IS THE NEW VERSION
      !!    !!gmat=0.d0
      !!    !!ii=0
      !!    !!ilrold=-1
      !!    !!do iorb=1,orbs%norbp
      !!    !!    iiorb=orbs%isorb+iorb
      !!    !!    ilr=orbs%inwhichlocreg(iiorb)
      !!    !!    if(ilr>ilrold) then
      !!    !!        ii=ii+1
      !!    !!    end if
      !!    !!    do jorb=1,orbs%norb
      !!    !!        jlr=orbs%inwhichlocreg(jorb)
      !!    !!        tt=0.d0
      !!    !!        do lorb=1,orbs%norb
      !!    !!            !tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
      !!    !!            !          - kernel(iiorb,lorb)*potmat(lorb,jorb,ilr)
      !!    !!            tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
      !!    !!                      - kernel(iiorb,lorb)*potmatsmall(lorb,jorb,ii)
      !!    !!        end do
      !!    !!        gmat(jorb,iiorb)=-tt
      !!    !!    end do
      !!    !!    ilrold=ilr
      !!    !!end do 
      !!    !!do iorb=1,orbs%norb
      !!    !!    ilr=orbs%inwhichlocreg(iorb)
      !!    !!    jlrold=-1
      !!    !!    jj=0
      !!    !!    do jorb=1,orbs%norbp
      !!    !!        jjorb=orbs%isorb+jorb
      !!    !!        jlr=orbs%inwhichlocreg(jjorb)
      !!    !!        if(jlr>jlrold) then
      !!    !!            jj=jj+1
      !!    !!        end if
      !!    !!        tt=0.d0
      !!    !!        do lorb=1,orbs%norb
      !!    !!            !tt = tt + kernel(jjorb,lorb)*potmat(lorb,iorb,jlr)
      !!    !!            tt = tt + kernel(jjorb,lorb)*potmatsmall(lorb,iorb,jj)
      !!    !!        end do
      !!    !!        gmat(jjorb,iorb)=gmat(jjorb,iorb)-tt
      !!    !!        jlrold=jlr
      !!    !!    end do
      !!    !!end do 
      !!    !!call mpiallred(gmat(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
      !!    !!do iorb=1,orbs%norb
      !!    !!    do jorb=1,orbs%norb
      !!    !!        if(iproc==0) then
      !!    !!            write(77,*) iorb, jorb, gmat(jorb,iorb)
      !!    !!        end if
      !!    !!    end do
      !!    !!end do
      !!end if
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      !Build the complex matrix -iGmat
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              gmatc(jorb,iorb)=cmplx(0.d0,-gmat(jorb,iorb),kind=8)
              !!if(iproc==0) write(999,'(a,2i8,2es16.7)') 'iorb, jorb, gmatc(jorb,iorb)', iorb, jorb, gmatc(jorb,iorb)
          end do
      end do 
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1



      ! Diagonalize Gmatc
      t1=mpi_wtime()
      lwork=10*orbs%norb
      allocate(work(lwork), stat=istat) ! factor of 2 since it is assumed to be complex
      allocate(rwork(lwork), stat=istat)
      call zheev('v', 'l', orbs%norb, gmatc(1,1), orbs%norb, eval(1), work, lwork, rwork, info)
      if(info/=0) stop 'ERROR in zheev 3'
      deallocate(work)
      deallocate(rwork)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1


      !!$$! Calculate step size
      !!$$if(it==1) then
      !!$$    !!if(.not.newgradient) then
      !!$$        lstep=5.d-2/(maxval(eval))
      !!$$    !!else
      !!$$    !!    lstep=1.d-4/(maxval(eval))
      !!$$    !!end if
      !!$$else
      !!$$    lstep=2.d0*lstep_optimal
      !!$$    !lstep=1.d-3/(maxval(eval))
      !!$$end if
      lstep=-2.d-3/(maxval(eval))

      t1=mpi_wtime()
      ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
      ! This is also a diagonal matrix, so only calculate the diagonal part.
      do iorb=1,orbs%norb
         ttc=cmplx(0.d0,-lstep*eval(iorb),kind=8)
         expD_cmplx(iorb)=(0.d0,0.d0)
          do k=0,50
              expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
          end do
      end do
      t2=mpi_wtime()
      time_exponential=time_exponential+t2-t1

      t1=mpi_wtime()
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              if(iorb==jorb) then
                  tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
              else
                  tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
              end if
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
           gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
      call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
           tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1

      t1=mpi_wtime()
      ! Build new lphi
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      ! Update Umat
      call dcopy(orbs%norb**2, Umat(1,1), 1, tempmat3(1,1,2), 1)
      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
           tempmat3(1,1,1), orbs%norb, 0.d0, Umat(1,1), orbs%norb)

      !!t1=mpi_wtime()
      !!!write(*,*) '5: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      !!call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
      !!     comon%recvbuf, tempmat3(1,1,1), .true., lphi)
      !!t2=mpi_wtime()
      !!time_lincomb=time_lincomb+t2-t1

      !!$$t1=mpi_wtime()
      !!$$call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
      !!$$     confdatarr, hx, lphi, -1, lvphi)
      !!$$t2=mpi_wtime()
      !!$$time_convol=time_convol+t2-t1

      !!$$!!if(.not.newgradient) then
      !!$$    energyconf_trial=ddot(max(orbs%npsidim_orbs,orbs%npsidim_comp), lphi(1), 1, lvphi(1), 1)
      !!$$    call mpiallred(energyconf_trial, 1, mpi_sum, mpi_comm_world, ierr)
      !!$$!!else

      !!$$!!    call dcopy(comon%nrecvbuf, comon%recvbuf, 1, recvbuf, 1)
      !!$$!!    !call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      !!$$!!    !call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
      !!$$!!    !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
      !!$$!!    !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
      !!$$!!    !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
      !!$$!!    !deallocate(ttmat)
      !!$$!!    call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, Kmat)
      !!$$!!    call dcopy(comon%nrecvbuf, recvbuf, 1, comon%recvbuf, 1)

      !!$$!!    energyconf_trial=0.d0
      !!$$!!    do iorb=1,orbs%norb
      !!$$!!        do jorb=1,orbs%norb
      !!$$!!            energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(jorb,iorb)
      !!$$!!            !energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(iorb,jorb)
      !!$$!!        end do
      !!$$!!    end do
      !!$$!!end if

      !!$$! Calculate the gradient of the confinement
      !!$$energyconf_der0=0.d0
      !!$$do iorb=1,orbs%norb
      !!$$    do jorb=1,orbs%norb
      !!$$        energyconf_der0=energyconf_der0+gmat(jorb,iorb)**2
      !!$$    end do
      !!$$end do
      !!$$energyconf_der0=-.5d0*energyconf_der0

      !!$$! Calculate optimal step size
      !!$$lstep_optimal = -energyconf_der0*lstep**2/(2.d0*(energyconf_trial-energyconf_0-lstep*energyconf_der0))
      !!$$!!if(.not.newgradient) then
      !!$$    lstep_optimal=min(lstep_optimal,lstep)
      !!$$!!else
      !!$$!!    if(lstep_optimal<0) then
      !!$$!!        lstep_optimal=lstep
      !!$$!!    else
      !!$$!!        lstep_optimal=min(lstep_optimal,lstep)
      !!$$!!    end if
      !!$$!!end if

      !!$$t1=mpi_wtime()
      !!$$! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
      !!$$! This is also a diagonal matrix, so only calculate the diagonal part.
      !!$$do iorb=1,orbs%norb
      !!$$   ttc=cmplx(0.d0,-lstep_optimal*eval(iorb),kind=8)
      !!$$   expD_cmplx(iorb)=(0.d0,0.d0)
      !!$$    do k=0,50
      !!$$        expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
      !!$$    end do
      !!$$end do
      !!$$t2=mpi_wtime()
      !!$$time_exponential=time_exponential+t2-t1

      !!$$t1=mpi_wtime()
      !!$$do iorb=1,orbs%norb
      !!$$    do jorb=1,orbs%norb
      !!$$        if(iorb==jorb) then
      !!$$            tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
      !!$$        else
      !!$$            tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
      !!$$        end if
      !!$$    end do
      !!$$end do
      !!$$t2=mpi_wtime()
      !!$$time_matrixmodification=time_matrixmodification+t2-t1

      !!$$t1=mpi_wtime()
      !!$$call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
      !!$$     gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
      !!$$call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
      !!$$     tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
      !!$$t2=mpi_wtime()
      !!$$time_linalg=time_linalg+t2-t1


      !!$$! Build new lphi
      !!$$do iorb=1,orbs%norb
      !!$$    do jorb=1,orbs%norb
      !!$$        tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
      !!$$    end do
      !!$$end do
      !!$$t1=mpi_wtime()
      !!$$call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
      !!$$     comon%recvbuf, tempmat3(1,1,1), .true., lphi)
      !!$$t2=mpi_wtime()
      !!$$time_lincomb=time_lincomb+t2-t1


      !!$$!if(it<nit) then
      !!$$!    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      !!$$!    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
      !!$$!end if


  end do innerLoop


  !!!!!! EXPERIMENTAL
  !!!!write(*,*) 'WARNING: TRANSPOSE UMAT!!!!'
  !!!!call dcopy(orbs%norb**2, Umat(1,1), 1, tempmat3(1,1,1), 1)
  !!!!do iorb=1,orbs%norb
  !!!!    do jorb=1,orbs%norb
  !!!!        Umat(jorb,iorb)=tempmat3(iorb,jorb,1)
  !!!!    end do
  !!!!end do


  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        do jorb=1,orbs%norb
  !!            write(999,*) iorb,jorb,Umat(jorb,iorb)
  !!        end do
  !!    end do
  !!end if


  t1=mpi_wtime()
  !write(*,*) '5: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
  call extractOrbital3(iproc, nproc, orbs, orbs%npsidim_orbs, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
  call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
  call collectnew(iproc, nproc, comon, mad, op, orbs, lzd, comon%nsendbuf, &
       comon%sendbuf, comon%nrecvbuf, comon%recvbuf, tt3, tt4, tt5)
  call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
       comon%recvbuf, Umat(1,1), .true., lphi)
  t2=mpi_wtime()
  time_lincomb=time_lincomb+t2-t1


  ! Recalculate the norm of all basis functions
  normarr=0.d0
  ist=1
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      ncnt=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      normarr(iiorb)=ddot(ncnt, lphi(ist), 1, lphi(ist), 1)
      ist=ist+ncnt
  end do
  call mpiallred(normarr(1), orbs%norb, mpi_sum, mpi_comm_world, ierr)
  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        write(*,'(a,i6,es16.6)') 'new norm after loop: iorb, normarr(iorb)', iorb, normarr(iorb)
  !!    end do
  !!end if

  ! Recalculate X,Y,Z
  order=1
  call apply_position_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lxphi, lyphi, lzphi)
  ! Build the matrices X, Y, Z
  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lxphi, mad, X)
  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lyphi, mad, Y)
  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lzphi, mad, Z)


!!!  !!! CHECK #####################################
!!!  order=1
!!!  call apply_r_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lvphi)
!!!  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R)
!!!
!!!  order=2
!!!  call apply_r_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lvphi)
!!!  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R2)
!!!
!!!  call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
!!!       R(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
!!!  call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
!!!       tempmat3(1,1,1), orbs%norb, 0.d0, R(1,1), orbs%norb)
!!!
!!!  call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
!!!       R2(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
!!!  call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
!!!       tempmat3(1,1,1), orbs%norb, 0.d0, R2(1,1), orbs%norb)
!!!
!!!  ! New calculation of the spread
!!!  var=0.d0
!!!  do iorb=1,orbs%norb
!!!      tt = R2(iorb,iorb)/normarr(iorb)-(R(iorb,iorb)/normarr(iorb))**2
!!!      var=var+tt
!!!  end do
!!!  if(iproc==0) write(*,'(a,es18.8)') 'FINAL: total variance', var



  do iorb=1,orbs%norb
      ilr=orbs%inwhichlocreg(iorb)
      !!if(iproc==0) write(*,'(a,2i5,3f10.4)') 'END: iorb, ilr, centers: ', &
      !!    iorb, ilr, X(iorb,iorb)/normarr(iorb), Y(iorb,iorb)/normarr(iorb), Z(iorb,iorb)/normarr(iorb)
      !!centers(1,iorb)=X(iorb,iorb)/normarr(iorb)
      !!centers(2,iorb)=Y(iorb,iorb)/normarr(iorb)
      !!centers(3,iorb)=Z(iorb,iorb)/normarr(iorb)
      centers(1,ilr)=X(iorb,iorb)/normarr(iorb)
      centers(2,ilr)=Y(iorb,iorb)/normarr(iorb)
      centers(3,ilr)=Z(iorb,iorb)/normarr(iorb)
      centers_end(1,iorb)=X(iorb,iorb)/normarr(iorb)
      centers_end(2,iorb)=Y(iorb,iorb)/normarr(iorb)
      centers_end(3,iorb)=Z(iorb,iorb)/normarr(iorb)
  end do

  !!write(*,*) 'ATTENTION DEBUG'
  !!centers=locregCenters

  do iorb=1,lzd%nlr
      tt = (centers_start(1,iorb)-centers_end(1,iorb))**2 &
          +(centers_start(2,iorb)-centers_end(2,iorb))**2 &
          +(centers_start(3,iorb)-centers_end(3,iorb))**2 
      tt=sqrt(tt)
      !!if(iproc==0) write(*,'(a,i5,es9.2)') 'iorb, shift of the center: ', iorb, tt
  end do


!!  ! OTHER CHECK
!!  call apply_rminusmu_operator(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, centers, lvphi)
!!  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R2)
!!  ! New calculation of the spread
!!  var=0.d0
!!  do iorb=1,orbs%norb
!!      tt = R2(iorb,iorb)/normarr(iorb)
!!      var=var+tt
!!  end do
!!  if(iproc==0) write(*,'(a,es18.8)') 'FINAL2: total variance', var
!!
!!  call apply_rminusmu_operator(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, locregCenters, lvphi)
!!  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R2)
!!  ! New calculation of the spread
!!  var=0.d0
!!  do iorb=1,orbs%norb
!!      tt = R2(iorb,iorb)/normarr(iorb)
!!      var=var+tt
!!  end do
!!  if(iproc==0) write(*,'(a,es18.8)') 'FINAL3: total variance', var



  t2_tot=mpi_wtime()
  time_tot=t2_tot-t1_tot
  call mpiallred(time_convol, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_commun, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_lincomb, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_linalg, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_matrixmodification, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_exponential, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_matrixelements, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_tot, 1, mpi_max, mpi_comm_world, ierr)
  !time_convol=time_convol/dble(nproc)
  !time_commun=time_commun/dble(nproc)
  !time_lincomb=time_lincomb/dble(nproc)
  !time_linalg=time_linalg/dble(nproc)
  !time_matrixmodification=time_matrixmodification/dble(nproc)
  !time_exponential_=time_exponential/dble(nproc)
  !time_tot=time_tot/dble(nproc)
  !!if(iproc==0) then
  !!    write(*,'(a,es16.6)') 'total time: ',time_tot
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'convolutions: ',time_convol,' (',time_convol/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear combinations: ',time_lincomb,' (',time_lincomb/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'communication: ',time_commun,' (',time_commun/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear algebra: ',time_linalg,' (',time_linalg/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix modification: ',time_matrixmodification, &
  !!                                   ' (',time_matrixmodification/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'building exponential: ',time_exponential,' (',time_exponential/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix elements ',time_matrixelements,' (',time_matrixelements/time_tot*100.d0,'%)'
  !!end if


  iall=-product(shape(gmat))*kind(gmat)
  deallocate(gmat, stat=istat)
  call memocc(istat, iall, 'gmat', subname)
  iall=-product(shape(gmatc))*kind(gmatc)
  deallocate(gmatc, stat=istat)
  !call memocc(istat, iall, 'gmatc', subname)
  iall=-product(shape(omatc))*kind(omatc)
  deallocate(omatc, stat=istat)
  !call memocc(istat, iall, 'omatc', subname)
  iall=-product(shape(tempmat3))*kind(tempmat3)
  deallocate(tempmat3, stat=istat)
  call memocc(istat, iall, 'tempmat3', subname)
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)
  iall=-product(shape(expD_cmplx))*kind(expD_cmplx)
  deallocate(expD_cmplx, stat=istat)
  call memocc(istat, iall, 'expD_cmplx', subname)
  iall=-product(shape(tempmatc))*kind(tempmatc)
  deallocate(tempmatc, stat=istat)
  !call memocc(istat, iall, 'tempmatc', subname)
  iall=-product(shape(hamtrans))*kind(hamtrans)
  deallocate(hamtrans, stat=istat)
  call memocc(istat, iall, 'hamtrans', subname)
  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  deallocate(lphiovrlp, stat=istat)
  call memocc(istat, iall, 'lphiovrlp', subname)
  iall=-product(shape(lvphi))*kind(lvphi)
  deallocate(lvphi, stat=istat)
  call memocc(istat, iall, 'lvphi', subname)
  iall=-product(shape(Kmat))*kind(Kmat)
  deallocate(Kmat, stat=istat)
  call memocc(istat, iall, 'Kmat', subname)
  !!iall=-product(shape(recvbuf))*kind(recvbuf)
  !!deallocate(recvbuf, stat=istat)
  !!call memocc(istat, iall, 'recvbuf', subname)

  iall=-product(shape(potmat))*kind(potmat)
  deallocate(potmat, stat=istat)
  call memocc(istat, iall, 'potmat', subname)
  iall=-product(shape(potmatsmall))*kind(potmatsmall)
  deallocate(potmatsmall, stat=istat)
  call memocc(istat, iall, 'potmatsmall', subname)


  iall=-product(shape(lxphi))*kind(lxphi)
  deallocate(lxphi, stat=istat)
  call memocc(istat, iall, 'lxphi', subname)
  iall=-product(shape(lyphi))*kind(lyphi)
  deallocate(lyphi, stat=istat)
  call memocc(istat, iall, 'lyphi', subname)
  iall=-product(shape(lzphi))*kind(lzphi)
  deallocate(lzphi, stat=istat)
  call memocc(istat, iall, 'lzphi', subname)



  iall=-product(shape(X))*kind(X)
  deallocate(X, stat=istat)
  call memocc(istat, iall, 'X', subname)
  iall=-product(shape(Y))*kind(Y)
  deallocate(Y, stat=istat)
  call memocc(istat, iall, 'Y', subname)
  iall=-product(shape(Z))*kind(Z)
  deallocate(Z, stat=istat)
  call memocc(istat, iall, 'Z', subname)
  iall=-product(shape(X2))*kind(X2)
  deallocate(X2, stat=istat)
  call memocc(istat, iall, 'X2', subname)
  iall=-product(shape(R))*kind(R)
  deallocate(R, stat=istat)
  call memocc(istat, iall, 'R', subname)
  iall=-product(shape(R2))*kind(R2)
  deallocate(R2, stat=istat)
  call memocc(istat, iall, 'R2', subname)
  iall=-product(shape(Y2))*kind(Y2)
  deallocate(Y2, stat=istat)
  call memocc(istat, iall, 'Y2', subname)
  iall=-product(shape(Z2))*kind(Z2)
  deallocate(Z2, stat=istat)
  call memocc(istat, iall, 'Z2', subname)
  iall=-product(shape(Xd))*kind(Xd)
  deallocate(Xd, stat=istat)
  call memocc(istat, iall, 'Xd', subname)
  iall=-product(shape(Yd))*kind(Yd)
  deallocate(Yd, stat=istat)
  call memocc(istat, iall, 'Yd', subname)
  iall=-product(shape(Zd))*kind(Zd)
  deallocate(Zd, stat=istat)
  call memocc(istat, iall, 'Zd', subname)
  iall=-product(shape(Xprime))*kind(Xprime)
  deallocate(Xprime, stat=istat)
  call memocc(istat, iall, 'Xprime', subname)
  iall=-product(shape(Yprime))*kind(Yprime)
  deallocate(Yprime, stat=istat)
  call memocc(istat, iall, 'Yprime', subname)
  iall=-product(shape(Zprime))*kind(Zprime)
  deallocate(Zprime, stat=istat)
  call memocc(istat, iall, 'Zprime', subname)
  iall=-product(shape(Xprimesquare))*kind(Xprimesquare)
  deallocate(Xprimesquare, stat=istat)
  call memocc(istat, iall, 'Xprimesquare', subname)
  iall=-product(shape(Yprimesquare))*kind(Yprimesquare)
  deallocate(Yprimesquare, stat=istat)
  call memocc(istat, iall, 'Yprimesquare', subname)
  iall=-product(shape(Zprimesquare))*kind(Zprimesquare)
  deallocate(Zprimesquare, stat=istat)
  call memocc(istat, iall, 'Zprimesquare', subname)

  iall=-product(shape(commutX1))*kind(commutX1)
  deallocate(commutX1, stat=istat)
  call memocc(istat, iall, 'commutX1', subname)
  iall=-product(shape(commutY1))*kind(commutY1)
  deallocate(commutY1, stat=istat)
  call memocc(istat, iall, 'commutY1', subname)
  iall=-product(shape(commutZ1))*kind(commutZ1)
  deallocate(commutZ1, stat=istat)
  call memocc(istat, iall, 'commutZ1', subname)
  iall=-product(shape(commutX2))*kind(commutX2)
  deallocate(commutX2, stat=istat)
  call memocc(istat, iall, 'commutX2', subname)
  iall=-product(shape(commutY2))*kind(commutY2)
  deallocate(commutY2, stat=istat)
  call memocc(istat, iall, 'commutY2', subname)
  iall=-product(shape(commutZ2))*kind(commutZ2)
  deallocate(commutZ2, stat=istat)
  call memocc(istat, iall, 'commutZ2', subname)
  iall=-product(shape(commutX3))*kind(commutX3)
  deallocate(commutX3, stat=istat)
  call memocc(istat, iall, 'commutX3', subname)
  iall=-product(shape(commutY3))*kind(commutY3)
  deallocate(commutY3, stat=istat)
  call memocc(istat, iall, 'commutY3', subname)
  iall=-product(shape(commutZ3))*kind(commutZ3)
  deallocate(commutZ3, stat=istat)
  call memocc(istat, iall, 'commutZ3', subname)

  !!iall=-product(shape(X0))*kind(X0)
  !!deallocate(X0, stat=istat)
  !!call memocc(istat, iall, 'X0', subname)
  !!iall=-product(shape(H))*kind(H)
  !!deallocate(H, stat=istat)
  !!call memocc(istat, iall, 'H', subname)
  !!iall=-product(shape(gHmat))*kind(gHmat)
  !!deallocate(gHmat, stat=istat)
  !!call memocc(istat, iall, 'gHmat', subname)

  !!iall=-product(shape(Y0))*kind(Y0)
  !!deallocate(Y0, stat=istat)
  !!call memocc(istat, iall, 'Y0', subname)
  !!iall=-product(shape(M))*kind(M)
  !!deallocate(M, stat=istat)
  !!call memocc(istat, iall, 'M', subname)
  !!iall=-product(shape(gMmat))*kind(gMmat)
  !!deallocate(gMmat, stat=istat)
  !!call memocc(istat, iall, 'gMmat', subname)

  !!iall=-product(shape(Z0))*kind(Z0)
  !!deallocate(Z0, stat=istat)
  !!call memocc(istat, iall, 'Z0', subname)
  !!iall=-product(shape(N))*kind(N)
  !!deallocate(N, stat=istat)
  !!call memocc(istat, iall, 'N', subname)
  !!iall=-product(shape(gNmat))*kind(gNmat)
  !!deallocate(gNmat, stat=istat)
  !!call memocc(istat, iall, 'gNmat', subname)


  !!iall=-product(shape(HXd))*kind(HXd)
  !!deallocate(HXd, stat=istat)
  !!call memocc(istat, iall, 'HXd', subname)
  !!iall=-product(shape(HXdsquare))*kind(HXdsquare)
  !!deallocate(HXdsquare, stat=istat)
  !!call memocc(istat, iall, 'HXdsquare', subname)
  !!iall=-product(shape(HXdX0))*kind(HXdX0)
  !!deallocate(HXdX0, stat=istat)
  !!call memocc(istat, iall, 'HXdX0', subname)
  !!iall=-product(shape(HX0))*kind(HX0)
  !!deallocate(HX0, stat=istat)
  !!call memocc(istat, iall, 'HX0', subname)
  !!iall=-product(shape(HX0square))*kind(HX0square)
  !!deallocate(HX0square, stat=istat)
  !!call memocc(istat, iall, 'HX0square', subname)
  !!iall=-product(shape(Xdsquare))*kind(Xdsquare)
  !!deallocate(Xdsquare, stat=istat)
  !!call memocc(istat, iall, 'Xdsquare', subname)
  !!iall=-product(shape(XdX0))*kind(XdX0)
  !!deallocate(XdX0, stat=istat)
  !!call memocc(istat, iall, 'XdX0', subname)
  !!iall=-product(shape(X0square))*kind(X0square)
  !!deallocate(X0square, stat=istat)
  !!call memocc(istat, iall, 'X0square', subname)

  !!iall=-product(shape(MYd))*kind(MYd)
  !!deallocate(MYd, stat=istat)
  !!call memocc(istat, iall, 'MYd', subname)
  !!iall=-product(shape(MYdsquare))*kind(MYdsquare)
  !!deallocate(MYdsquare, stat=istat)
  !!call memocc(istat, iall, 'MYdsquare', subname)
  !!iall=-product(shape(MYdY0))*kind(MYdY0)
  !!deallocate(MYdY0, stat=istat)
  !!call memocc(istat, iall, 'MYdY0', subname)
  !!iall=-product(shape(MY0))*kind(MY0)
  !!deallocate(MY0, stat=istat)
  !!call memocc(istat, iall, 'MY0', subname)
  !!iall=-product(shape(MY0square))*kind(MY0square)
  !!deallocate(MY0square, stat=istat)
  !!call memocc(istat, iall, 'MY0square', subname)
  !!iall=-product(shape(Ydsquare))*kind(Ydsquare)
  !!deallocate(Ydsquare, stat=istat)
  !!call memocc(istat, iall, 'Ydsquare', subname)
  !!iall=-product(shape(YdY0))*kind(YdY0)
  !!deallocate(YdY0, stat=istat)
  !!call memocc(istat, iall, 'YdY0', subname)
  !!iall=-product(shape(Y0square))*kind(Y0square)
  !!deallocate(Y0square, stat=istat)
  !!call memocc(istat, iall, 'Y0square', subname)

  !!iall=-product(shape(NZd))*kind(NZd)
  !!deallocate(NZd, stat=istat)
  !!call memocc(istat, iall, 'NZd', subname)
  !!iall=-product(shape(NZdsquare))*kind(NZdsquare)
  !!deallocate(NZdsquare, stat=istat)
  !!call memocc(istat, iall, 'NZdsquare', subname)
  !!iall=-product(shape(NZdZ0))*kind(NZdZ0)
  !!deallocate(NZdZ0, stat=istat)
  !!call memocc(istat, iall, 'NZdZ0', subname)
  !!iall=-product(shape(NZ0))*kind(NZ0)
  !!deallocate(NZ0, stat=istat)
  !!call memocc(istat, iall, 'NZ0', subname)
  !!iall=-product(shape(NZ0square))*kind(NZ0square)
  !!deallocate(NZ0square, stat=istat)
  !!call memocc(istat, iall, 'NZ0square', subname)
  !!iall=-product(shape(Zdsquare))*kind(Zdsquare)
  !!deallocate(Zdsquare, stat=istat)
  !!call memocc(istat, iall, 'Zdsquare', subname)
  !!iall=-product(shape(ZdZ0))*kind(ZdZ0)
  !!deallocate(ZdZ0, stat=istat)
  !!call memocc(istat, iall, 'ZdZ0', subname)
  !!iall=-product(shape(Z0square))*kind(Z0square)
  !!deallocate(Z0square, stat=istat)
  !!call memocc(istat, iall, 'Z0square', subname)

  iall=-product(shape(centers_start))*kind(centers_start)
  deallocate(centers_start, stat=istat)
  call memocc(istat, iall, 'centers_start', subname)
  iall=-product(shape(centers_end))*kind(centers_end)
  deallocate(centers_end, stat=istat)
  call memocc(istat, iall, 'centers_end', subname)
  iall=-product(shape(normarr))*kind(normarr)
  deallocate(normarr, stat=istat)
  call memocc(istat, iall, 'normarr', subname)


  call deallocateRecvBufferOrtho(comon, subname)
  call deallocateSendBufferOrtho(comon, subname)

end subroutine MLWFnew










subroutine build_new_linear_combinations(iproc, nproc, lzd, orbs, op, nrecvbuf, recvbuf, omat, reset, lphi)
use module_base
use module_types
implicit none

!Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(overlapParameters),intent(in):: op
integer,intent(in):: nrecvbuf
real(8),dimension(nrecvbuf),intent(in):: recvbuf
real(8),dimension(orbs%norb,orbs%norb),intent(in):: omat
logical,intent(in):: reset
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: lphi

! Local variables
integer:: ist, jst, ilrold, iorb, iiorb, ilr, ncount, jorb, jjorb, i, ldim, ind, indout, gdim, iorbref, m, ii
integer:: istart, iend, iseg, start, ifine, igrid, irecv, kseg, kold, kstart, kend 
real(8):: tt

   call timing(iproc,'build_lincomb ','ON')

      ! Build new lphi
      if(reset) then
          lphi=0.d0
      end if

      indout=1
      ilrold=-1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          !if(ilr>ilrold) then
          if(ilr/=ilrold) then
              iorbref=iorb
          end if
          gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
          do jorb=1,op%noverlaps(iiorb)
              jjorb=op%overlaps(jorb,iiorb)
              !jjorb=op%overlaps(jorb,ilr)
              jst=op%indexInRecvBuf(iorbref,jjorb)
              ldim=op%wfd_overlap(jorb,iorbref)%nvctr_c+7*op%wfd_overlap(jorb,iorbref)%nvctr_f
              tt=omat(jjorb,iiorb)
              !!tt=tt*lzd%cutoffweight(jjorb,iiorb)
              !Coarse part
              kold=1
              do iseg=1,op%wfd_overlap(jorb,iorbref)%nseg_c
                  istart=op%wfd_overlap(jorb,iorbref)%keyglob(1,iseg)
                  iend=op%wfd_overlap(jorb,iorbref)%keyglob(2,iseg)
                  ncount=iend-istart+1
                  inner_loop: do kseg=kold,lzd%llr(ilr)%wfd%nseg_c
                     kstart = lzd%llr(ilr)%wfd%keyglob(1,kseg)
                     kend   = lzd%llr(ilr)%wfd%keyglob(2,kseg)
                     if(kstart <= iend .and. kend >= istart) then 
                        kold = kseg + 1
                        start = lzd%llr(ilr)%wfd%keyvglob(kseg) + max(0,istart-kstart)
                        call daxpy(ncount, tt, recvBuf(jst), 1, lphi(indout+start-1), 1)
                        jst=jst+ncount
                        exit inner_loop
                     end if
                  end do inner_loop
              end do
              ! Fine part
              kold = 1
              jst=op%indexInRecvBuf(iorbref,jjorb)              
              do iseg=1,op%wfd_overlap(jorb,iorbref)%nseg_f
                 istart=op%wfd_overlap(jorb,iorbref)%keyglob(1,iseg+op%wfd_overlap(jorb,iorbref)%nseg_c)
                 iend=op%wfd_overlap(jorb,iorbref)%keyglob(2,iseg+op%wfd_overlap(jorb,iorbref)%nseg_c)
                 start = op%wfd_overlap(jorb,iorbref)%keyvglob(iseg+op%wfd_overlap(jorb,iorbref)%nseg_c)
                 ncount=7*(iend-istart+1)
                 inner_loop2: do kseg=kold,lzd%llr(ilr)%wfd%nseg_f
                    kstart = lzd%llr(ilr)%wfd%keyglob(1,kseg+lzd%llr(ilr)%wfd%nseg_c)
                    kend   = lzd%llr(ilr)%wfd%keyglob(2,kseg+lzd%llr(ilr)%wfd%nseg_c)
                    if(kstart <= iend .and. kend >= istart) then 
                       kold = kseg + 1
                       start = lzd%llr(ilr)%wfd%nvctr_c+(lzd%llr(ilr)%wfd%keyvglob(kseg+lzd%llr(ilr)%wfd%nseg_c) +&
                              max(0,istart-kstart)-1)*7
                       call daxpy(ncount, tt, recvBuf(jst+op%wfd_overlap(jorb,iorbref)%nvctr_c), 1,&
                               lphi(indout+start), 1)
                       jst=jst+ncount
                       exit inner_loop2
                    end if
                  end do inner_loop2
              end do
          end do
          indout=indout+gdim
          ilrold=ilr

      end do

   call timing(iproc,'build_lincomb ','OF')
          

end subroutine build_new_linear_combinations




!!subroutine flatten_at_edges(iproc, nproc, lin, at, input,hx, hy, hz, orbs, lzd, rxyz, psi)
!!use module_base
!!use module_types
!!use module_interfaces
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!real(gp),intent(in):: hx, hy, hz
!!type(linearParameters),intent(in):: lin
!!type(atoms_data),intent(in):: at
!!type(input_variables),intent(in):: input
!!type(orbitals_data),intent(in):: orbs
!!type(local_zone_descriptors),intent(in):: lzd
!!real(8),dimension(3,at%nat),intent(in):: rxyz
!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
!!
!!! Local variables
!!integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, iiorb
!!real(8):: hxh, hyh, hzh, alpha
!!type(workarr_sumrho):: work_sr
!!real(8),dimension(:,:),allocatable:: psir
!!character(len=*),parameter:: subname='flatten_at_edges'
!!
!!
!!
!!  oidx = 0
!!  do iorb=1,orbs%norbp
!!     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
!!  
!!     !initialise the work arrays
!!     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)
!!
!!     ! Wavefunction in real space
!!     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     call memocc(i_stat,psir,'psir',subname)
!!     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)
!!
!!     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
!!     !the psir wavefunction is given in the spinorial form
!!
!!     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
!!     !apply the potential to the psir wavefunction and calculate potential energy
!!     hxh=.5d0*hx
!!     hyh=.5d0*hy
!!     hzh=.5d0*hz
!!     !icenter=confinementCenter(iorb)
!!     !icenter=lin%orbs%inWhichLocregp(iorb)
!!     iiorb=orbs%isorb+iorb
!!     icenter=lin%orbs%inWhichLocreg(iiorb)
!!     !components of the potential
!!     npot=orbs%nspinor
!!     if (orbs%nspinor == 2) npot=1
!!
!!     alpha=-log(1.d-5)/(1.d0-.7d0*lin%locrad(icenter))**2
!!     !write(*,*) 'iproc, iorb, alpha', iproc, iorb, alpha
!!     call flatten(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
!!          rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
!!          lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3, .7d0*lin%locrad(icenter), alpha, &
!!          lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!
!!
!!
!!     call isf_to_daub(lzd%llr(ilr), work_sr, psir, psi(1+oidx))
!!
!!     i_all=-product(shape(psir))*kind(psir)
!!     deallocate(psir,stat=i_stat)
!!     call memocc(i_stat,i_all,'psir',subname)
!!
!!     call deallocate_work_arrays_sumrho(work_sr)
!!
!!     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
!!
!!  enddo
!!
!!
!!end subroutine flatten_at_edges


subroutine flatten(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
     rxyzConfinement, hxh, hyh, hzh, potentialPrefac, confPotOrder, offsetx, offsety, offsetz, cut, alpha, &
     ibyyzz_r) !optional
use module_base
implicit none
integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor, confPotOrder, offsetx, offsety, offsetz
real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
real(8),dimension(3),intent(in):: rxyzConfinement
real(8),intent(in):: hxh, hyh, hzh, potentialPrefac, cut, alpha
!local variables
integer :: i1,i2,i3,i1s,i1e,ispinor, order
real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
real(gp) :: epot_p, epot, values, valuesold


  !the Tail treatment is allowed only in the Free BC case
  if (nbuf /= 0 .and. nl1*nl2*nl3 == 0) stop 'NONSENSE: nbuf/=0 only for Free BC'

  ! The order of the cofinement potential (we take order divided by two, 
  ! since later we calculate (r**2)**order.
  if(confPotOrder==2) then
      ! parabolic potential
      order=1
  else if(confPotOrder==4) then
      ! quartic potential
      order=2
  else if(confPotOrder==6) then
      ! sextic potential
      order=3
  end if
  
  values=0.d0
  valuesold=0.d0

!!!$omp parallel default(private)&
!!!$omp shared(psir,n1,n2,n3,epot,ibyyzz_r,nl1,nl2,nl3,nbuf,nspinor)
  !case without bounds
  i1s=-14*nl1
  i1e=2*n1+1+15*nl1
  !epot_p=0._gp
!!!$omp do
  do i3=-14*nl3,2*n3+1+15*nl3
     if (i3 >= -14+2*nbuf .and. i3 <= 2*n3+16-2*nbuf) then !check for the nbuf case
        do i2=-14*nl2,2*n2+1+15*nl2
           if (i2 >= -14+2*nbuf .and. i2 <= 2*n2+16-2*nbuf) then !check for the nbuf case
              !this if statement is inserted here for avoiding code duplication
              !it is to be seen whether the code results to be too much unoptimised
              if (present(ibyyzz_r)) then
                 !in this case we are surely in Free BC
                 !the min is to avoid to calculate for no bounds
                 do i1=-14+2*nbuf,min(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14-1
                    psir(i1,i2,i3,:)=0.0_wp
                 enddo
                 i1s=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf)
                 i1e=min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)
              end if
              !write(*,'(a,5i8)') 'iproc, i1, i2, i1s, i1e', iproc, i1, i2, i1s, i1e
              
              !here we put the branchments wrt to the spin
              if (nspinor == 4) then
                 stop 'this part is not yet implemented'
              else
                 do ispinor=1,nspinor
                    do i1=i1s,i1e
                       tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
                          (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2
                       tt=sqrt(tt)
                       tt=tt-cut

                       if(tt>0.d0) then
                           tt=exp(-alpha*tt**2)
                       else
                           tt=1.d0
                       end if
                       values=values+tt*abs(psir(i1,i2,i3,ispinor))
                       valuesold=valuesold+abs(psir(i1,i2,i3,ispinor))
     
                       tt=tt*psir(i1,i2,i3,ispinor)
                       psir(i1,i2,i3,ispinor)=tt
                    end do
                 end do
              end if
              
              if (present(ibyyzz_r)) then
                 !the max is to avoid the calculation for no bounds
                 do i1=max(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14+1,2*n1+16-2*nbuf
                    psir(i1,i2,i3,:)=0.0_wp
                 enddo
              end if

           else
              do i1=-14,2*n1+16
                 psir(i1,i2,i3,:)=0.0_wp
              enddo
           endif
        enddo
     else
        do i2=-14,2*n2+16
           do i1=-14,2*n1+16
              psir(i1,i2,i3,:)=0.0_wp
           enddo
        enddo
     endif
  enddo
!!!$omp end do

!!!$omp end parallel

write(*,*) 'values',values/valuesold

END SUBROUTINE flatten
!!***





subroutine get_potential_matrices(iproc, nproc, at, orbs, lzd, op, comon, mad, rxyz, &
           confdatarr, hx, psi, potmat)
use module_base
use module_types
use module_interfaces, eccept_this_one => get_potential_matrices
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(atoms_data),intent(in):: at
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
type(p2pComms),intent(inout):: comon
type(matrixDescriptors),intent(in):: mad
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),intent(in):: hx
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
real(8),dimension(orbs%norb,orbs%norb,at%nat),intent(out):: potmat

! Local variables
integer:: iorb, ilr, ilrold, istat, iall
real(8),dimension(:,:),allocatable:: ttmat
real(8):: tt1, tt2, tt3, tt4, tt5
real(8),dimension(:),allocatable:: vpsi
character(len=*),parameter:: subname='get_potential_matrices'

allocate(vpsi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, vpsi, 'vpsi', subname)


ilrold=-1
do iorb=1,orbs%norb
    ilr=orbs%inwhichlocreg(iorb)
    if(ilr==ilrold) cycle
    call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
         confdatarr, hx, psi, ilr, vpsi)

    !call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, vpsi, comon%nsendBuf, comon%sendBuf)
    !call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, vpsi, comon, tt1, tt2)
    !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
    !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
    !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
    !deallocate(ttmat)
    call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, psi, vpsi, mad, potmat(1,1,ilr))
    ilrold=ilr
    
end do

iall=-product(shape(vpsi))*kind(vpsi)
deallocate(vpsi, stat=istat)
call memocc(istat, iall, 'vpsi', subname)



end subroutine get_potential_matrices




!!subroutine get_potential_matrices_new(iproc, nproc, lin, at, input, orbs, lzd, op, comon, rxyz, psi, nlocregOnMPI, potmat)
!!use module_base
!!use module_types
!!use module_interfaces
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, nlocregOnMPI
!!type(linearParameters),intent(inout):: lin
!!type(atoms_data),intent(in):: at
!!type(input_variables),intent(in):: input
!!type(orbitals_data),intent(in):: orbs
!!type(local_zone_descriptors),intent(in):: lzd
!!type(overlapParameters),intent(inout):: op
!!type(p2pCommsOrthonormality),intent(inout):: comon
!!real(8),dimension(3,at%nat),intent(in):: rxyz
!!real(8),dimension(lzd%lpsidimtot),intent(inout):: psi
!!real(8),dimension(orbs%norb,orbs%norb,nlocregOnMPI),intent(out):: potmat
!!
!!! Local variables
!!integer:: iorb, ilr, ilrold, istat, iall, ii, iiorb
!!real(8),dimension(:,:),allocatable:: ttmat
!!real(8):: tt1, tt2, tt3, tt4, tt5
!!real(8),dimension(:),allocatable:: vpsi
!!character(len=*),parameter:: subname='get_potential_matrices'
!!integer,dimensioN(:),allocatable:: sendcounts, displs
!!
!!
!!allocate(sendcounts(0:nproc-1), stat=istat)
!!call memocc(istat, sendcounts, 'sendcounts', subname)
!!allocate(displs(0:nproc-1), stat=istat)
!!call memocc(istat, displs, 'displs', subname)
!!
!!call getCommunArraysMatrixCompression(iproc, nproc, orbs, lin%mad, sendcounts, displs)
!!
!!
!!allocate(vpsi(lzd%lpsidimtot), stat=istat)
!!call memocc(istat, vpsi, 'vpsi', subname)
!!
!!call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, vpsi, comon%nsendBuf, comon%sendBuf)
!!call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, vpsi, comon, tt1, tt2)
!!allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
!!call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
!!     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
!!deallocate(ttmat)
!!
!!ilrold=-1
!!ii=0
!!do iorb=1,orbs%norbp
!!    iiorb=orbs%isorb+iorb
!!    ilr=orbs%inwhichlocreg(iiorb)
!!    if(ilr==ilrold) cycle
!!    ii=ii+1
!!    call apply_orbitaldependent_potential(iproc, nproc, lin, at, input, lin%orbs, lin%lzd, rxyz, psi, ilr, vpsi)
!!    !call calculateOverlapMatrix3(iproc, nproc, orbs, op, orbs%inWhichLocreg, comon%nsendBuf, &
!!    !                             comon%sendBuf, comon%nrecvBuf, comon%recvBuf, lin%mad, potmat(1,1,ii))
!!    !call getMatrixElements2(iproc, nproc, lin%lzd, lin%orbs, lin%op, comon, psi, vpsi, lin%mad, potmat(1,1,ilr))
!!
!!    call calculateOverlapMatrix3Partial(iproc, nproc, orbs, op, orbs%inwhichlocreg, comon%nsendBuf, comon%sendBuf, &
!!         comon%nrecvBuf, comon%recvBuf, lin%mad, tempmat)
!!
!!    ilrold=ilr
!!    
!!end do
!!
!!iall=-product(shape(vpsi))*kind(vpsi)
!!deallocate(vpsi, stat=istat)
!!call memocc(istat, iall, 'vpsi', subname)
!!
!!iall=-product(shape(sendcounts))*kind(sendcounts)
!!deallocate(sendcounts, stat=istat)
!!call memocc(istat, iall, 'sendcounts', subname)
!!
!!iall=-product(shape(displs))*kind(displs)
!!deallocate(displs, stat=istat)
!!call memocc(istat, iall, 'displs', subname)
!!
!!end subroutine get_potential_matrices_new



!!!subroutine get_potential_matrices_new2(iproc, nproc, lin, at, input, orbs, lzd, op, comon, mad, rxyz, &
!!!           psi, nlocregOnMPI, potmat)
!!!use module_base
!!!use module_types
!!!use module_interfaces
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc, nlocregOnMPI
!!!type(linearParameters),intent(inout):: lin
!!!type(atoms_data),intent(in):: at
!!!type(input_variables),intent(in):: input
!!!type(orbitals_data),intent(in):: orbs
!!!type(local_zone_descriptors),intent(in):: lzd
!!!type(overlapParameters),intent(inout):: op
!!!type(p2pCommsOrthonormality),intent(inout):: comon
!!!type(matrixDescriptors),intent(in):: mad
!!!real(8),dimension(3,at%nat),intent(in):: rxyz
!!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
!!!real(8),dimension(orbs%norb,orbs%norb,nlocregOnMPI),intent(out):: potmat
!!!
!!!! Local variables
!!!integer:: iorb, ilr, ilrold, istat, iall, ii, iiorb, ncount, iilr, jjorb, ist, jst, jorb, iorbout
!!!real(8),dimension(:,:),allocatable:: ttmat
!!!real(8):: tt1, tt2, tt3, tt4, tt5, ddot
!!!real(8),dimension(:),allocatable:: vpsi
!!!character(len=*),parameter:: subname='get_potential_matrices_new2'
!!!integer,dimensioN(:),allocatable:: sendcounts, displs
!!!
!!!
!!!
!!!
!!!allocate(vpsi(comon%nrecvbuf), stat=istat)
!!!call memocc(istat, vpsi, 'vpsi', subname)
!!!
!!!call extractOrbital3(iproc, nproc, orbs, max(orbs%npsidim_orbs,orbs%npsidim_comp), orbs%inWhichLocreg, &
!!!     lzd, op, psi, comon%nsendBuf, comon%sendBuf)
!!!call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, psi, comon, tt1, tt2)
!!!!!allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
!!!call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, lin%lzd, comon%nsendbuf, &
!!!     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, tt3, tt4, tt5)
!!!!!deallocate(ttmat)
!!!
!!!! Now all other psi are in the receive buffer. Apply the potential to them.
!!!iilr=0
!!!ilrold=-1
!!!do iorbout=1,orbs%norbp
!!!    ilr=orbs%inwhichlocreg(orbs%isorb+iorbout)
!!!    if(ilr==ilrold) cycle
!!!    iilr=iilr+1
!!!    do iorb=1,orbs%norbp
!!!        iiorb=orbs%isorb+iorb
!!!        do jorb=1,op%noverlaps(iiorb)
!!!            jjorb=op%overlaps(jorb,iiorb)
!!!            call getStartingIndices(iorb, jorb, op, orbs, ist, jst)
!!!            ncount=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
!!!    
!!!            call apply_orbitaldependent_potential_foronelocreg(iproc, nproc, op%olr(jorb,iorb), lin, &
!!!                 at, input, orbs, rxyz, ncount, &
!!!                 comon%recvbuf, iilr, vpsi(1))
!!!            potmat(jjorb,iiorb,iilr)=ddot(ncount, comon%sendBuf(ist), 1, vpsi(1), 1)
!!!        end do
!!!    end do
!!!    ilrold=ilr
!!!end do
!!!
!!!!!ilrold=-1
!!!!!ii=0
!!!!!do iorb=1,orbs%norbp
!!!!!    iiorb=orbs%isorb+iorb
!!!!!    ilr=orbs%inwhichlocreg(iiorb)
!!!!!    if(ilr==ilrold) cycle
!!!!!    ii=ii+1
!!!!!    call apply_orbitaldependent_potential(iproc, nproc, lin, at, input, lin%orbs, lin%lzd, rxyz, psi, ilr, vpsi)
!!!!!    !call calculateOverlapMatrix3(iproc, nproc, orbs, op, orbs%inWhichLocreg, comon%nsendBuf, &
!!!!!    !                             comon%sendBuf, comon%nrecvBuf, comon%recvBuf, lin%mad, potmat(1,1,ii))
!!!!!    !call getMatrixElements2(iproc, nproc, lin%lzd, lin%orbs, lin%op, comon, psi, vpsi, lin%mad, potmat(1,1,ilr))
!!!!!
!!!!!    call calculateOverlapMatrix3Partial(iproc, nproc, orbs, op, orbs%inwhichlocreg, comon%nsendBuf, comon%sendBuf, &
!!!!!         comon%nrecvBuf, comon%recvBuf, lin%mad, tempmat)
!!!!!
!!!!!    ilrold=ilr
!!!!!    
!!!!!end do
!!!
!!!iall=-product(shape(vpsi))*kind(vpsi)
!!!deallocate(vpsi, stat=istat)
!!!call memocc(istat, iall, 'vpsi', subname)
!!!
!!!iall=-product(shape(sendcounts))*kind(sendcounts)
!!!deallocate(sendcounts, stat=istat)
!!!call memocc(istat, iall, 'sendcounts', subname)
!!!
!!!iall=-product(shape(displs))*kind(displs)
!!!deallocate(displs, stat=istat)
!!!call memocc(istat, iall, 'displs', subname)
!!!
!!!end subroutine get_potential_matrices_new2




subroutine apply_position_operators(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi, order, xpsi, ypsi, zpsi)
use module_base
use module_types
use module_interfaces, except_this_one => apply_position_operators
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, order
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),intent(in):: hx, hy, hz
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: xpsi, ypsi, zpsi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
real(8),dimension(:,:),allocatable:: psir, psirx, psiry, psirz
type(workarr_sumrho):: work_sr
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:,:,:),allocatable:: ypsitemp_c
real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
character(len=*),parameter:: subname='apply_position_operators'
integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation
interface
subroutine position_operator(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
     hxh, hyh, hzh, dir, &
     ibyyzz_r) !optional
use module_base
implicit none
integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor
real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
real(8),intent(in):: hxh, hyh, hzh
character(len=1),intent(in):: dir
integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
end subroutine
end interface

  ishift=(/0,0,0/)

  xpsi=0.d0
  ypsi=0.d0
  zpsi=0.d0
  oidx = 0
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)

     iiorb=orbs%isorb+iorb
     !!write(*,'(a,4i8,4x,3i6)') 'iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)', &
     !!    iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)
     !!write(*,'(a,3i8,6i6)') 'iproc, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3', &
     !!    1, lzd%llr(ilr)%d%n1i, 1, lzd%llr(ilr)%d%n2i, 1, lzd%llr(ilr)%d%n3i
  
     !initialise the work arrays
     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)

     ! Wavefunction in real space
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

     allocate(psirx(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psirx,'psirx',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirx)

     allocate(psiry(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psiry,'psiry',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psiry)

     allocate(psirz(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psirz,'psirz',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirz)

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form

     !psi(1+oidx+lzd%llr(ilr)%wfd%nvctr_c:1+oidx+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0

     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
     !!do i_stat=1,Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
     !!    write(1000+iproc,'(i9,es18.7,i9)') i_stat, psir(i_stat,1), Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
     !!end do
     !apply the potential to the psir wavefunction and calculate potential energy
     hxh=.5d0*hx
     hyh=.5d0*hy
     hzh=.5d0*hz
     !icenter=confinementCenter(iorb)
     !components of the potential
     npot=orbs%nspinor
     if (orbs%nspinor == 2) npot=1

     !!call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
     !!     rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
     !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     call position_operators(lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
                             psir, order, psirx, psiry, psirz, &
                             confdatarr(iorb), lzd%llr(ilr)%bounds%ibyyzz_r) !optional

     call isf_to_daub(lzd%llr(ilr), work_sr, psirx, xpsi(1+oidx))
     call isf_to_daub(lzd%llr(ilr), work_sr, psiry, ypsi(1+oidx))
     call isf_to_daub(lzd%llr(ilr), work_sr, psirz, zpsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'x', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, xpsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'y', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, ypsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'z', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, zpsi(1+oidx))

     !!iall=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot x', iorb, ddot(iall, xpsi(1+oidx), 1, psi(1+oidx), 1)
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot y', iorb, ddot(iall, ypsi(1+oidx), 1, psi(1+oidx), 1)
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot z', iorb, ddot(iall, zpsi(1+oidx), 1, psi(1+oidx), 1)


     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     i_all=-product(shape(psirx))*kind(psirx)
     deallocate(psirx,stat=i_stat)
     call memocc(i_stat,i_all,'psirx',subname)

     i_all=-product(shape(psiry))*kind(psiry)
     deallocate(psiry,stat=i_stat)
     call memocc(i_stat,i_all,'psiry',subname)

     i_all=-product(shape(psirz))*kind(psirz)
     deallocate(psirz,stat=i_stat)
     call memocc(i_stat,i_all,'psirz',subname)

     call deallocate_work_arrays_sumrho(work_sr)

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

  enddo


end subroutine apply_position_operators



!!!!!subroutine position_operator(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
!!!!!     hxh, hyh, hzh, ioffset, dir, &
!!!!!     ibyyzz_r) !optional
!!!!!use module_base
!!!!!implicit none
!!!!!integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor
!!!!!integer,dimension(3),intent(in):: ioffset
!!!!!real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
!!!!!integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
!!!!!real(8),intent(in):: hxh, hyh, hzh
!!!!!character(len=1),intent(in):: dir
!!!!!!local variables
!!!!!integer :: i1,i2,i3,i1s,i1e,ispinor, order
!!!!!real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
!!!!!real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
!!!!!real(gp) :: epot_p, epot
!!!!!
!!!!!
!!!!!  !the Tail treatment is allowed only in the Free BC case
!!!!!  if (nbuf /= 0 .and. nl1*nl2*nl3 == 0) stop 'NONSENSE: nbuf/=0 only for Free BC'
!!!!!
!!!!!
!!!!!
!!!!!!!!$omp parallel default(private)&
!!!!!!!!$omp shared(psir,n1,n2,n3,epot,ibyyzz_r,nl1,nl2,nl3,nbuf,nspinor)
!!!!!  !case without bounds
!!!!!  i1s=-14*nl1
!!!!!  i1e=2*n1+1+15*nl1
!!!!!  !epot_p=0._gp
!!!!!!!!$omp do
!!!!!  do i3=-14*nl3,2*n3+1+15*nl3
!!!!!     if (i3 >= -14+2*nbuf .and. i3 <= 2*n3+16-2*nbuf) then !check for the nbuf case
!!!!!        do i2=-14*nl2,2*n2+1+15*nl2
!!!!!           if (i2 >= -14+2*nbuf .and. i2 <= 2*n2+16-2*nbuf) then !check for the nbuf case
!!!!!              !this if statement is inserted here for avoiding code duplication
!!!!!              !it is to be seen whether the code results to be too much unoptimised
!!!!!              if (present(ibyyzz_r)) then
!!!!!                 !in this case we are surely in Free BC
!!!!!                 !the min is to avoid to calculate for no bounds
!!!!!                 do i1=-14+2*nbuf,min(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14-1
!!!!!                    psir(i1,i2,i3,:)=0.0_wp
!!!!!                 enddo
!!!!!                 i1s=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf)
!!!!!                 i1e=min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)
!!!!!              end if
!!!!!              !write(*,'(a,5i8)') 'iproc, i1, i2, i1s, i1e', iproc, i1, i2, i1s, i1e
!!!!!
!!!!!              !here we put the branchments wrt to the spin
!!!!!              if (nspinor == 4) then
!!!!!                 stop 'this part is not yet implemented'
!!!!!              else
!!!!!                 do ispinor=1,nspinor
!!!!!                    do i1=i1s,i1e
!!!!!                       if(dir=='x') then
!!!!!                           tt=dble(i1+ioffset(1))*hxh
!!!!!                       else if(dir=='y') then
!!!!!                           tt=dble(i2+ioffset(2))*hyh
!!!!!                       else if(dir=='z') then
!!!!!                           tt=dble(i3+ioffset(3))*hzh
!!!!!                       else
!!!!!                           stop 'wrong direction!'
!!!!!                       end if
!!!!!                       tt=tt*psir(i1,i2,i3,ispinor)
!!!!!                       if(dir=='x') write(100,'(4i9,es11.2,es16.5)') ioffset(1), i1, i2, i3, tt, psir(i1,i2,i3,ispinor)
!!!!!                       psir(i1,i2,i3,ispinor)=tt
!!!!!                    end do
!!!!!                 end do
!!!!!              end if
!!!!!
!!!!!              if (present(ibyyzz_r)) then
!!!!!                 !the max is to avoid the calculation for no bounds
!!!!!                 do i1=max(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14+1,2*n1+16-2*nbuf
!!!!!                    psir(i1,i2,i3,:)=0.0_wp
!!!!!                 enddo
!!!!!              end if
!!!!!
!!!!!           else
!!!!!              do i1=-14,2*n1+16
!!!!!                 psir(i1,i2,i3,:)=0.0_wp
!!!!!              enddo
!!!!!           endif
!!!!!        enddo
!!!!!     else
!!!!!        do i2=-14,2*n2+16
!!!!!           do i1=-14,2*n1+16
!!!!!              psir(i1,i2,i3,:)=0.0_wp
!!!!!           enddo
!!!!!        enddo
!!!!!     endif
!!!!!  enddo
!!!!!!!!$omp end do
!!!!!
!!!!!!!!$omp end parallel
!!!!!
!!!!!END SUBROUTINE position_operator




subroutine position_operators(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,order,&
     psirx, psiry, psirz, &
     confdata,ibyyzz_r) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,order
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(in) :: psir !< real-space wfn in lr
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(out) :: psirx, psiry, psirz !< x,y,z operator applied to real-space wfn in lr
  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
  real(wp):: ttx, tty, ttz, potx, poty, potz


  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))


  !$omp parallel default(none)&
  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)&
  !$omp private(ispinor,i1,i2,i3,epot_p,i1st,i1et)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4)

!!$  !$omp parallel default(private)&
!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
  !case without bounds


  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psirx(i1,i2,i3,ispinor)=0.0_wp 
             psiry(i1,i2,i3,ispinor)=0.0_wp 
             psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psirx(i1,i2,i3,ispinor)=0.0_wp 
              psiry(i1,i2,i3,ispinor)=0.0_wp 
              psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psirx(i1,i2,i3,ispinor)=0.0_wp 
              psiry(i1,i2,i3,ispinor)=0.0_wp 
              psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psirx(i1,i2,i3,ispinor)=0.0_wp 
              psiry(i1,i2,i3,ispinor)=0.0_wp 
              psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psirx(i1,i2,i3,ispinor)=0.0_wp 
              psiry(i1,i2,i3,ispinor)=0.0_wp 
              psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do

     end do
     !$omp end do
  end do


  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psirx(i1,i2,i3,ispinor)=0.0_wp 
              psiry(i1,i2,i3,ispinor)=0.0_wp 
              psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do


  !important part of the array
  if (nspinor==4) then
      stop 'not yet implemented for nspinor==4!'
     !!!$omp do
     !!do i3=i3s,i3e
     !!   do i2=i2s,i2e
     !!      !thanks to the optional argument the conditional is done at compile time
     !!      if (present(ibyyzz_r)) then
     !!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
     !!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
     !!      else
     !!         i1st=i1s
     !!         i1et=i1e
     !!      end if
     !!      !no need of setting up to zero values outside wavefunction bounds
     !!      do i1=i1st,i1et
     !!         !wavefunctions
     !!         psir1=psir(i1,i2,i3,1)
     !!         psir2=psir(i1,i2,i3,2)
     !!         psir3=psir(i1,i2,i3,3)
     !!         psir4=psir(i1,i2,i3,4)
     !!         !potentials + confining term
     !!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
     !!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
     !!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
     !!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)

     !!         !diagonal terms
     !!         tt11=pot1*psir1 !p1
     !!         tt22=pot1*psir2 !p2
     !!         tt33=pot4*psir3 !p3
     !!         tt44=pot4*psir4 !p4
     !!         !Rab*Rb
     !!         tt13=pot2*psir3 !p1
     !!         !Iab*Ib
     !!         tt14=pot3*psir4 !p1
     !!         !Rab*Ib
     !!         tt23=pot2*psir4 !p2
     !!         !Iab*Rb
     !!         tt24=pot3*psir3 !p2
     !!         !Rab*Ra
     !!         tt31=pot2*psir1 !p3
     !!         !Iab*Ia
     !!         tt32=pot3*psir2 !p3
     !!         !Rab*Ia
     !!         tt41=pot2*psir2 !p4
     !!         !Iab*Ra
     !!         tt42=pot3*psir1 !p4

     !!         !value of the potential energy
     !!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
     !!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3

     !!         !wavefunction update
     !!         !p1=h1p1+h2p3-h3p4
     !!         !p2=h1p2+h2p4+h3p3
     !!         !p3=h2p1+h3p2+h4p3
     !!         !p4=h2p2-h3p1+h4p4
     !!         psir(i1,i2,i3,1)=tt11+tt13-tt14
     !!         psir(i1,i2,i3,2)=tt22+tt23+tt24
     !!         psir(i1,i2,i3,3)=tt33+tt31+tt32
     !!         psir(i1,i2,i3,4)=tt44+tt41-tt42
     !!      end do
     !!   end do
     !!end do
     !!!$omp end do

  else !case with nspinor /=4
     do ispinor=1,nspinor
        !$omp do
        do i3=i3s,i3e
           do i2=i2s,i2e
              !thanks to the optional argument the conditional is done at compile time
              if (present(ibyyzz_r)) then
                 i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
                 i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
              else
                 i1st=i1s
                 i1et=i1e
              end if
              !no need of setting up to zero values outside wavefunction bounds
              do i1=i1st,i1et
                 psir1=psir(i1,i2,i3,ispinor)
                 !the local potential is always real (npot=1) + confining term
                 !!pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
                 potx=(confdata%hh(1)*real(i1+confdata%ioffset(1),wp))**order
                 poty=(confdata%hh(2)*real(i2+confdata%ioffset(2),wp))**order
                 potz=(confdata%hh(3)*real(i3+confdata%ioffset(3),wp))**order

                 ttx=potx*psir1
                 tty=poty*psir1
                 ttz=potz*psir1

                 psirx(i1,i2,i3,ispinor)=ttx
                 psiry(i1,i2,i3,ispinor)=tty
                 psirz(i1,i2,i3,ispinor)=ttz
              end do
           end do
        end do
        !$omp end do
     end do
  end if
  
  
  !$omp end parallel


END SUBROUTINE position_operators







subroutine commutator(norb, A, B, res)
implicit none

! Calling arguments
integer,intent(in):: norb
real(8),dimension(norb,norb),intent(in):: A, B
real(8),dimension(norb,norb),intent(out):: res

! Local variables
real(8),dimension(norb,norb):: AB, BA

call dgemm('n', 'n', norb, norb, norb, 1.d0, A, norb, B, norb, 0.d0, AB, norb)
call dgemm('n', 'n', norb, norb, norb, 1.d0, B, norb, A, norb, 0.d0, BA, norb)
res=AB-BA

end subroutine commutator








subroutine check_cutoff(iproc, nproc, orbs, lzd, hx, hy, hz, locrad, confdatarr, psi)
use module_base
use module_types
use module_interfaces, except_this_one => apply_position_operators
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),intent(in):: hx, hy, hz, locrad
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time, weight_in, weight_out
real(8),dimension(:,:),allocatable:: psir, psirx, psiry, psirz
type(workarr_sumrho):: work_sr
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:,:,:),allocatable:: ypsitemp_c
real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
character(len=*),parameter:: subname='apply_orbitaldependent_potential'
integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation

  ishift=(/0,0,0/)

  oidx = 0
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)

  
     !initialise the work arrays
     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)

     ! Wavefunction in real space
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
     !apply the potential to the psir wavefunction and calculate potential energy
     !icenter=confinementCenter(iorb)
     !components of the potential
     npot=orbs%nspinor
     if (orbs%nspinor == 2) npot=1

     write(*,'(a,2i8,3es16.7)') 'iproc, iorb, confdatarr(iorb)%rxyzConf', iproc, iorb, confdatarr(iorb)%rxyzConf

     call get_cutoff_weight(lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
                             psir, locrad, weight_in, weight_out, &
                             confdatarr(iorb), lzd%llr(ilr)%bounds%ibyyzz_r) !optional

     write(*,'(a,2i8,3es16.6)') 'iproc, iorb, weight_in, weight_out, ratio', &
         iproc, iorb, weight_in, weight_out, weight_in/(weight_in+weight_out)

     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     call deallocate_work_arrays_sumrho(work_sr)

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

  enddo


end subroutine check_cutoff



subroutine get_cutoff_weight(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,&
     cutoff, weight_in, weight_out, &
     confdata,ibyyzz_r) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(in) :: psir !< real-space wfn in lr
  real(8),intent(in):: cutoff
  real(8),intent(out):: weight_in, weight_out
  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
  real(wp):: ttx, tty, ttz, cutoff2

  ! Square of the cutoff radius
  cutoff2=cutoff**2

  ! Initialize return values
  weight_in=0.d0
  weight_out=0.d0


  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))


  !$omp parallel default(none)&
  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)&
  !$omp private(ispinor,i1,i2,i3,epot_p,i1st,i1et)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4)

!!$  !$omp parallel default(private)&
!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
  !case without bounds


  !!!put to zero the external part of psir if the potential is more little than the wavefunction
  !!!first part of the array
  !!do ispinor=1,nspinor
  !!   !$omp do 
  !!   do i3=1,i3s-1
  !!      do i2=1,n2i
  !!         do i1=1,n1i
  !!           psirx(i1,i2,i3,ispinor)=0.0_wp 
  !!           psiry(i1,i2,i3,ispinor)=0.0_wp 
  !!           psirz(i1,i2,i3,ispinor)=0.0_wp 
  !!         end do
  !!      end do
  !!   end do
  !!   !$omp end do
  !!end do

  !!!central part of the array
  !!do ispinor=1,nspinor
  !!   !$omp do 
  !!   do i3=i3s,i3e

  !!      !first part
  !!      do i2=1,i2s-1
  !!         do i1=1,n1i
  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
  !!         end do
  !!      end do
  !!      !central part
  !!      do i2=i2s,i2e
  !!         do i1=1,i1s-1
  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
  !!         end do
  !!         do i1=i1e+1,n1i
  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
  !!         end do
  !!      end do
  !!      !last part
  !!      do i2=i2e+1,n2i
  !!         do i1=1,n1i
  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
  !!         end do
  !!      end do

  !!   end do
  !!   !$omp end do
  !!end do


  !!!last part of the array
  !!do ispinor=1,nspinor
  !!   !$omp do 
  !!   do i3=i3e+1,n3i
  !!      do i2=1,n2i
  !!         do i1=1,n1i
  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
  !!         end do
  !!      end do
  !!   end do
  !!   !$omp end do
  !!end do


  !important part of the array
  if (nspinor==4) then
      stop 'not yet implemented for nspinor==4!'
     !!!$omp do
     !!do i3=i3s,i3e
     !!   do i2=i2s,i2e
     !!      !thanks to the optional argument the conditional is done at compile time
     !!      if (present(ibyyzz_r)) then
     !!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
     !!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
     !!      else
     !!         i1st=i1s
     !!         i1et=i1e
     !!      end if
     !!      !no need of setting up to zero values outside wavefunction bounds
     !!      do i1=i1st,i1et
     !!         !wavefunctions
     !!         psir1=psir(i1,i2,i3,1)
     !!         psir2=psir(i1,i2,i3,2)
     !!         psir3=psir(i1,i2,i3,3)
     !!         psir4=psir(i1,i2,i3,4)
     !!         !potentials + confining term
     !!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
     !!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
     !!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
     !!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)

     !!         !diagonal terms
     !!         tt11=pot1*psir1 !p1
     !!         tt22=pot1*psir2 !p2
     !!         tt33=pot4*psir3 !p3
     !!         tt44=pot4*psir4 !p4
     !!         !Rab*Rb
     !!         tt13=pot2*psir3 !p1
     !!         !Iab*Ib
     !!         tt14=pot3*psir4 !p1
     !!         !Rab*Ib
     !!         tt23=pot2*psir4 !p2
     !!         !Iab*Rb
     !!         tt24=pot3*psir3 !p2
     !!         !Rab*Ra
     !!         tt31=pot2*psir1 !p3
     !!         !Iab*Ia
     !!         tt32=pot3*psir2 !p3
     !!         !Rab*Ia
     !!         tt41=pot2*psir2 !p4
     !!         !Iab*Ra
     !!         tt42=pot3*psir1 !p4

     !!         !value of the potential energy
     !!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
     !!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3

     !!         !wavefunction update
     !!         !p1=h1p1+h2p3-h3p4
     !!         !p2=h1p2+h2p4+h3p3
     !!         !p3=h2p1+h3p2+h4p3
     !!         !p4=h2p2-h3p1+h4p4
     !!         psir(i1,i2,i3,1)=tt11+tt13-tt14
     !!         psir(i1,i2,i3,2)=tt22+tt23+tt24
     !!         psir(i1,i2,i3,3)=tt33+tt31+tt32
     !!         psir(i1,i2,i3,4)=tt44+tt41-tt42
     !!      end do
     !!   end do
     !!end do
     !!!$omp end do

  else !case with nspinor /=4
     do ispinor=1,nspinor
        !$omp do
        do i3=i3s,i3e
           do i2=i2s,i2e
              !thanks to the optional argument the conditional is done at compile time
              if (present(ibyyzz_r)) then
                 i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
                 i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
              else
                 i1st=i1s
                 i1et=i1e
              end if
              !no need of setting up to zero values outside wavefunction bounds
              do i1=i1st,i1et
                 psir1=psir(i1,i2,i3,ispinor)
                 !the local potential is always real (npot=1) + confining term
                 !!pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
                 ttx=(confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-confdata%rxyzConf(1))**2
                 tty=(confdata%hh(2)*real(i2+confdata%ioffset(2),wp)-confdata%rxyzConf(2))**2
                 ttz=(confdata%hh(3)*real(i3+confdata%ioffset(3),wp)-confdata%rxyzConf(3))**2
                 tt=ttx+tty+ttz
                 !write(1000,*) tt, cutoff2, psir1**2
                 if(tt>cutoff2) then
                     weight_out=weight_out+psir1**2
                 else
                     weight_in=weight_in+psir1**2
                 end if
              end do
           end do
        end do
        !$omp end do
     end do
  end if
  
  
  !$omp end parallel


END SUBROUTINE get_cutoff_weight



subroutine apply_r_operators(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi, order, vpsi)
use module_base
use module_types
use module_interfaces, except_this_one => apply_r_operators
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, order
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),intent(in):: hx, hy, hz
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: vpsi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
real(8),dimension(:,:),allocatable:: psir, psirx, psiry, psirz
type(workarr_sumrho):: work_sr
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:,:,:),allocatable:: ypsitemp_c
real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
character(len=*),parameter:: subname='apply_orbitaldependent_potential'
integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation

  ishift=(/0,0,0/)

  !xpsi=0.d0
  !ypsi=0.d0
  vpsi=0.d0
  oidx = 0
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)

     iiorb=orbs%isorb+iorb
     !!write(*,'(a,4i8,4x,3i6)') 'iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)', &
     !!    iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)
     !!write(*,'(a,3i8,6i6)') 'iproc, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3', &
     !!    1, lzd%llr(ilr)%d%n1i, 1, lzd%llr(ilr)%d%n2i, 1, lzd%llr(ilr)%d%n3i
  
     !initialise the work arrays
     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)

     ! Wavefunction in real space
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

     !!allocate(psirx(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     !!call memocc(i_stat,psirx,'psirx',subname)
     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirx)

     !!allocate(psiry(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     !!call memocc(i_stat,psiry,'psiry',subname)
     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psiry)

     !!allocate(psirz(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     !!call memocc(i_stat,psirz,'psirz',subname)
     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirz)

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form

     !psi(1+oidx+lzd%llr(ilr)%wfd%nvctr_c:1+oidx+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0

     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
     !!!write(*,*) 'WARNING DEBUG in r_operator'
     !!!psir=1.d0/sqrt(dble(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i))
     !apply the potential to the psir wavefunction and calculate potential energy
     hxh=.5d0*hx
     hyh=.5d0*hy
     hzh=.5d0*hz
     !icenter=confinementCenter(iorb)
     !components of the potential
     npot=orbs%nspinor
     if (orbs%nspinor == 2) npot=1

     !!call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
     !!     rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
     !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     call r_operator(lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
                             psir, order, &
                             confdatarr(iorb), lzd%llr(ilr)%bounds%ibyyzz_r) !optional

     call isf_to_daub(lzd%llr(ilr), work_sr, psir, vpsi(1+oidx))
     !!call isf_to_daub(lzd%llr(ilr), work_sr, psiry, ypsi(1+oidx))
     !!call isf_to_daub(lzd%llr(ilr), work_sr, psirz, zpsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'x', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, xpsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'y', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, ypsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'z', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, zpsi(1+oidx))

     !!iall=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot x', iorb, ddot(iall, xpsi(1+oidx), 1, psi(1+oidx), 1)
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot y', iorb, ddot(iall, ypsi(1+oidx), 1, psi(1+oidx), 1)
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot z', iorb, ddot(iall, zpsi(1+oidx), 1, psi(1+oidx), 1)


     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     !!i_all=-product(shape(psirx))*kind(psirx)
     !!deallocate(psirx,stat=i_stat)
     !!call memocc(i_stat,i_all,'psirx',subname)

     !!i_all=-product(shape(psiry))*kind(psiry)
     !!deallocate(psiry,stat=i_stat)
     !!call memocc(i_stat,i_all,'psiry',subname)

     !!i_all=-product(shape(psirz))*kind(psirz)
     !!deallocate(psirz,stat=i_stat)
     !!call memocc(i_stat,i_all,'psirz',subname)

     call deallocate_work_arrays_sumrho(work_sr)

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

  enddo


end subroutine apply_r_operators








subroutine r_operator(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,order,&
     confdata,ibyyzz_r) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,order
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
  real(wp):: ttx, tty, ttz, potx, poty, potz

  if(order/=1 .and. order/=2) stop 'wrong order'

  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))


  !$omp parallel default(none)&
  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)&
  !$omp private(ispinor,i1,i2,i3,epot_p,i1st,i1et)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4)

!!$  !$omp parallel default(private)&
!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
  !case without bounds


  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psir(i1,i2,i3,ispinor)=0.0_wp 
             !!psiry(i1,i2,i3,ispinor)=0.0_wp 
             !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do

     end do
     !$omp end do
  end do


  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do


  !important part of the array
  if (nspinor==4) then
      stop 'not yet implemented for nspinor==4!'
     !!!$omp do
     !!do i3=i3s,i3e
     !!   do i2=i2s,i2e
     !!      !thanks to the optional argument the conditional is done at compile time
     !!      if (present(ibyyzz_r)) then
     !!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
     !!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
     !!      else
     !!         i1st=i1s
     !!         i1et=i1e
     !!      end if
     !!      !no need of setting up to zero values outside wavefunction bounds
     !!      do i1=i1st,i1et
     !!         !wavefunctions
     !!         psir1=psir(i1,i2,i3,1)
     !!         psir2=psir(i1,i2,i3,2)
     !!         psir3=psir(i1,i2,i3,3)
     !!         psir4=psir(i1,i2,i3,4)
     !!         !potentials + confining term
     !!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
     !!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
     !!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
     !!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)

     !!         !diagonal terms
     !!         tt11=pot1*psir1 !p1
     !!         tt22=pot1*psir2 !p2
     !!         tt33=pot4*psir3 !p3
     !!         tt44=pot4*psir4 !p4
     !!         !Rab*Rb
     !!         tt13=pot2*psir3 !p1
     !!         !Iab*Ib
     !!         tt14=pot3*psir4 !p1
     !!         !Rab*Ib
     !!         tt23=pot2*psir4 !p2
     !!         !Iab*Rb
     !!         tt24=pot3*psir3 !p2
     !!         !Rab*Ra
     !!         tt31=pot2*psir1 !p3
     !!         !Iab*Ia
     !!         tt32=pot3*psir2 !p3
     !!         !Rab*Ia
     !!         tt41=pot2*psir2 !p4
     !!         !Iab*Ra
     !!         tt42=pot3*psir1 !p4

     !!         !value of the potential energy
     !!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
     !!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3

     !!         !wavefunction update
     !!         !p1=h1p1+h2p3-h3p4
     !!         !p2=h1p2+h2p4+h3p3
     !!         !p3=h2p1+h3p2+h4p3
     !!         !p4=h2p2-h3p1+h4p4
     !!         psir(i1,i2,i3,1)=tt11+tt13-tt14
     !!         psir(i1,i2,i3,2)=tt22+tt23+tt24
     !!         psir(i1,i2,i3,3)=tt33+tt31+tt32
     !!         psir(i1,i2,i3,4)=tt44+tt41-tt42
     !!      end do
     !!   end do
     !!end do
     !!!$omp end do

  else !case with nspinor /=4
     do ispinor=1,nspinor
        !$omp do
        do i3=i3s,i3e
           do i2=i2s,i2e
              !thanks to the optional argument the conditional is done at compile time
              if (present(ibyyzz_r)) then
                 i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
                 i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
              else
                 i1st=i1s
                 i1et=i1e
              end if
              !no need of setting up to zero values outside wavefunction bounds
              do i1=i1st,i1et
                 psir1=psir(i1,i2,i3,ispinor)
                 !the local potential is always real (npot=1) + confining term
                 !!pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
                 ttx=(confdata%hh(1)*real(i1+confdata%ioffset(1),wp))**2
                 tty=(confdata%hh(2)*real(i2+confdata%ioffset(2),wp))**2
                 ttz=(confdata%hh(3)*real(i3+confdata%ioffset(3),wp))**2

                 tt = ttx+tty+ttz

                 if(order==1) then
                     tt=sqrt(tt)
                 end if

                 psir(i1,i2,i3,ispinor)=tt*psir1
              end do
           end do
        end do
        !$omp end do
     end do
  end if
  
  
  !$omp end parallel


END SUBROUTINE r_operator



subroutine apply_rminusmu_operator(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi, centers, vpsi)
use module_base
use module_types
use module_interfaces, except_this_one => apply_rminusmu_operator
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),intent(in):: hx, hy, hz
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
real(8),dimension(3,lzd%nlr),intent(in):: centers
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: vpsi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
real(8),dimension(:,:),allocatable:: psir, psirx, psiry, psirz
type(workarr_sumrho):: work_sr
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:,:,:),allocatable:: ypsitemp_c
real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
real(8),dimension(3):: mu
character(len=*),parameter:: subname='apply_orbitaldependent_potential'
integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation

  ishift=(/0,0,0/)

  !xpsi=0.d0
  !ypsi=0.d0
  vpsi=0.d0
  oidx = 0
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)

     iiorb=orbs%isorb+iorb
     !!write(*,'(a,4i8,4x,3i6)') 'iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)', &
     !!    iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)
     !!write(*,'(a,3i8,6i6)') 'iproc, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3', &
     !!    1, lzd%llr(ilr)%d%n1i, 1, lzd%llr(ilr)%d%n2i, 1, lzd%llr(ilr)%d%n3i
  
     !initialise the work arrays
     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)

     ! Wavefunction in real space
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

     !!allocate(psirx(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     !!call memocc(i_stat,psirx,'psirx',subname)
     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirx)

     !!allocate(psiry(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     !!call memocc(i_stat,psiry,'psiry',subname)
     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psiry)

     !!allocate(psirz(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     !!call memocc(i_stat,psirz,'psirz',subname)
     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirz)

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form

     !psi(1+oidx+lzd%llr(ilr)%wfd%nvctr_c:1+oidx+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0

     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
     !!!write(*,*) 'WARNING DEBUG in rminusmu_operator'
     !!!psir=1.d0/sqrt(dble(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i))
     !apply the potential to the psir wavefunction and calculate potential energy
     hxh=.5d0*hx
     hyh=.5d0*hy
     hzh=.5d0*hz
     !icenter=confinementCenter(iorb)
     !components of the potential
     npot=orbs%nspinor
     if (orbs%nspinor == 2) npot=1

     !!call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
     !!     rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
     !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     mu(1:3) = centers(1:3,ilr)
     write(*,'(a,2i8,3es18.7)') 'iproc, iorb, mu', iproc, iorb, mu
     call rminusmu_operator(lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
                             psir, mu, &
                             confdatarr(iorb), lzd%llr(ilr)%bounds%ibyyzz_r) !optional

     call isf_to_daub(lzd%llr(ilr), work_sr, psir, vpsi(1+oidx))
     !!call isf_to_daub(lzd%llr(ilr), work_sr, psiry, ypsi(1+oidx))
     !!call isf_to_daub(lzd%llr(ilr), work_sr, psirz, zpsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'x', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, xpsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'y', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, ypsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'z', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, zpsi(1+oidx))

     !!iall=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot x', iorb, ddot(iall, xpsi(1+oidx), 1, psi(1+oidx), 1)
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot y', iorb, ddot(iall, ypsi(1+oidx), 1, psi(1+oidx), 1)
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot z', iorb, ddot(iall, zpsi(1+oidx), 1, psi(1+oidx), 1)


     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     !!i_all=-product(shape(psirx))*kind(psirx)
     !!deallocate(psirx,stat=i_stat)
     !!call memocc(i_stat,i_all,'psirx',subname)

     !!i_all=-product(shape(psiry))*kind(psiry)
     !!deallocate(psiry,stat=i_stat)
     !!call memocc(i_stat,i_all,'psiry',subname)

     !!i_all=-product(shape(psirz))*kind(psirz)
     !!deallocate(psirz,stat=i_stat)
     !!call memocc(i_stat,i_all,'psirz',subname)

     call deallocate_work_arrays_sumrho(work_sr)

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

  enddo


end subroutine apply_rminusmu_operator



subroutine rminusmu_operator(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,mu,&
     confdata,ibyyzz_r) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
  real(8),dimension(3):: mu
  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
  real(wp):: ttx, tty, ttz, potx, poty, potz


  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))


  !$omp parallel default(none)&
  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)&
  !$omp private(ispinor,i1,i2,i3,epot_p,i1st,i1et)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4)

!!$  !$omp parallel default(private)&
!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
  !case without bounds


  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psir(i1,i2,i3,ispinor)=0.0_wp 
             !!psiry(i1,i2,i3,ispinor)=0.0_wp 
             !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do

     end do
     !$omp end do
  end do


  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do


  !important part of the array
  if (nspinor==4) then
      stop 'not yet implemented for nspinor==4!'
     !!!$omp do
     !!do i3=i3s,i3e
     !!   do i2=i2s,i2e
     !!      !thanks to the optional argument the conditional is done at compile time
     !!      if (present(ibyyzz_r)) then
     !!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
     !!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
     !!      else
     !!         i1st=i1s
     !!         i1et=i1e
     !!      end if
     !!      !no need of setting up to zero values outside wavefunction bounds
     !!      do i1=i1st,i1et
     !!         !wavefunctions
     !!         psir1=psir(i1,i2,i3,1)
     !!         psir2=psir(i1,i2,i3,2)
     !!         psir3=psir(i1,i2,i3,3)
     !!         psir4=psir(i1,i2,i3,4)
     !!         !potentials + confining term
     !!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
     !!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
     !!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
     !!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)

     !!         !diagonal terms
     !!         tt11=pot1*psir1 !p1
     !!         tt22=pot1*psir2 !p2
     !!         tt33=pot4*psir3 !p3
     !!         tt44=pot4*psir4 !p4
     !!         !Rab*Rb
     !!         tt13=pot2*psir3 !p1
     !!         !Iab*Ib
     !!         tt14=pot3*psir4 !p1
     !!         !Rab*Ib
     !!         tt23=pot2*psir4 !p2
     !!         !Iab*Rb
     !!         tt24=pot3*psir3 !p2
     !!         !Rab*Ra
     !!         tt31=pot2*psir1 !p3
     !!         !Iab*Ia
     !!         tt32=pot3*psir2 !p3
     !!         !Rab*Ia
     !!         tt41=pot2*psir2 !p4
     !!         !Iab*Ra
     !!         tt42=pot3*psir1 !p4

     !!         !value of the potential energy
     !!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
     !!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3

     !!         !wavefunction update
     !!         !p1=h1p1+h2p3-h3p4
     !!         !p2=h1p2+h2p4+h3p3
     !!         !p3=h2p1+h3p2+h4p3
     !!         !p4=h2p2-h3p1+h4p4
     !!         psir(i1,i2,i3,1)=tt11+tt13-tt14
     !!         psir(i1,i2,i3,2)=tt22+tt23+tt24
     !!         psir(i1,i2,i3,3)=tt33+tt31+tt32
     !!         psir(i1,i2,i3,4)=tt44+tt41-tt42
     !!      end do
     !!   end do
     !!end do
     !!!$omp end do

  else !case with nspinor /=4
     do ispinor=1,nspinor
        !$omp do
        do i3=i3s,i3e
           do i2=i2s,i2e
              !thanks to the optional argument the conditional is done at compile time
              if (present(ibyyzz_r)) then
                 i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
                 i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
              else
                 i1st=i1s
                 i1et=i1e
              end if
              !no need of setting up to zero values outside wavefunction bounds
              do i1=i1st,i1et
                 psir1=psir(i1,i2,i3,ispinor)
                 !the local potential is always real (npot=1) + confining term
                 !!pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
                 ttx=(confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-mu(1))**2
                 tty=(confdata%hh(2)*real(i2+confdata%ioffset(2),wp)-mu(2))**2
                 ttz=(confdata%hh(3)*real(i3+confdata%ioffset(3),wp)-mu(3))**2

                 tt = ttx+tty+ttz

                 psir(i1,i2,i3,ispinor)=tt*psir1
              end do
           end do
        end do
        !$omp end do
     end do
  end if
  
  
  !$omp end parallel


END SUBROUTINE rminusmu_operator
