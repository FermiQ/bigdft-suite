!> This subroutine initializes all parameters needed for the linear scaling version
!! and allocate all arrays.
subroutine allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, nlpspd, lin, phi, &
    input, rxyz, nscatterarr, tag, coeff, lphi)
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
!     input           type containing some very general parameters
!     rxyz            the atomic positions
!     occuprForINguess  delete maybe
!  Output arguments
!  ---------------------
!     phi             the localized basis functions. They are only initialized here, but
!                       not normalized.
!
use module_base
use module_types
use module_interfaces, exceptThisOne => allocateAndInitializeLinear
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(locreg_descriptors),intent(in):: Glr
type(orbitals_data),intent(in):: orbs
type(atoms_data),intent(inout):: at
type(nonlocal_psp_descriptors),intent(in):: nlpspd
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(in):: rxyz
integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer,intent(inout):: tag
real(8),dimension(:),pointer,intent(out):: phi
real(8),dimension(:,:),pointer,intent(out):: coeff
real(8),dimension(:),pointer,intent(out):: lphi

! Local variables
integer:: norb, norbu, norbd, istat, iat, ityp, iall, ilr, iorb, iiorb, ii, jproc
integer,dimension(:),allocatable:: norbsPerAtom
character(len=*),parameter:: subname='allocateAndInitializeLinear'
character(len=20),dimension(:),allocatable:: atomNames
real(8):: t1, t2
integer :: npsidim



! Nullify all pointers
call nullify_linearParameters(lin)

! Allocate all local arrays.
allocate(atomNames(at%ntypes), stat=istat)
call memocc(istat, atomNames, 'atomNames', subname)
allocate(norbsPerAtom(at%nat), stat=istat)
call memocc(istat, norbsPerAtom, 'norbsPerAtom', subname)

! Number of localization regions.
lin%nlr=at%nat
lin%lzd%nlr=at%nat
!lin%lb%lzd%nlr=at%nat

! Allocate the basic arrays that are needed for reading the input parameters.
call allocateBasicArrays(at, lin)

! Read in all parameters related to the linear scaling version.
call readLinearParameters(iproc, nproc, lin, at, atomNames)

! Count the number of basis functions.
norb=0
do iat=1,at%nat
    ityp=at%iatype(iat)
    norbsPerAtom(iat)=lin%norbsPerType(ityp)
    norb=norb+norbsPerAtom(iat)
end do
!!write(*,'(a,2i9)') 'iproc, norb',iproc, norb
!!! ATTENTION: DEBUG
!norb=31
!!!!!!!!!!!!!!!!!!!

! Distribute the basis functions among the processors.
norbu=norb
norbd=0
!call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor,&
!     input%nkpt, input%kpt, input%wkpt, lin%orbs)
call orbitals_descriptors_forLinear(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor,&
     input%nkpt, input%kpt, input%wkpt, lin%orbs)
call repartitionOrbitals(iproc, nproc, lin%orbs%norb, lin%orbs%norb_par, lin%orbs%norbp, lin%orbs%isorb_par, lin%orbs%isorb, lin%orbs%onWhichMPI)
!call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor,&
!     input%nkpt, input%kpt, input%wkpt, lin%gorbs)
call orbitals_descriptors_forLinear(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor,&
     input%nkpt, input%kpt, input%wkpt, lin%gorbs)
call repartitionOrbitals(iproc, nproc, lin%gorbs%norb, lin%gorbs%norb_par, lin%gorbs%norbp, lin%gorbs%isorb_par, lin%gorbs%isorb, lin%gorbs%onWhichMPI)
!!ii=0
!!do jproc=0,nproc-1
!!    ii=ii+lin%orbs%norb_par(jproc)
!!    if(iproc==0) write(*,'(a,2i9)') 'jproc, lin%orbs%norb_par(jproc)', jproc, lin%orbs%norb_par(jproc)
!!end do
!!if(ii/=lin%orbs%norb) then
!!    write(*,'(a,2(2x,i0))') 'ERROR: ii/=lin%orbs%norb', ii, lin%orbs%norb
!!end if



! Do the same again, but take into acount that we may also use the derivatives of the basis functions with
! respect to x,y,z. These informations will be stored in lin%lb%orbs. If we don't use the derivtaive, then
! lin%lb%orbs will be identical to lin%orbs.
if(.not. lin%useDerivativeBasisFunctions) then
    norb=lin%orbs%norb
    norbu=norb
    norbd=0
else
    norb=4*lin%orbs%norb
    norbu=norb
    norbd=0
end if
!call orbitals_descriptors(iproc,nproc,norb,norbu,norbd,input%nspin,orbs%nspinor,input%nkpt,input%kpt,input%wkpt,lin%lb%orbs)
call orbitals_descriptors_forLinear(iproc,nproc,norb,norbu,norbd,input%nspin,orbs%nspinor,input%nkpt,input%kpt,input%wkpt,lin%lb%orbs)
call repartitionOrbitals(iproc, nproc, lin%lb%orbs%norb, lin%lb%orbs%norb_par, lin%lb%orbs%norbp, lin%lb%orbs%isorb_par, lin%lb%orbs%isorb, lin%lb%orbs%onWhichMPI)
!call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor, input%nkpt, input%kpt, input%wkpt, &
!     lin%lb%gorbs)
call orbitals_descriptors_forLinear(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor, input%nkpt, input%kpt, input%wkpt, &
     lin%lb%gorbs)
call repartitionOrbitals(iproc, nproc, lin%lb%gorbs%norb, lin%lb%gorbs%norb_par, lin%lb%gorbs%norbp, lin%lb%gorbs%isorb_par, lin%lb%gorbs%isorb, lin%lb%gorbs%onWhichMPI)




! Assign the parameters needed for the communication to lin%comms. Again distinguish
! between the 'normal' basis and the 'large' basis inlcuding the derivtaives.
call orbitals_communicators(iproc,nproc,Glr,lin%orbs,lin%comms)
call orbitals_communicators(iproc,nproc,Glr,lin%lb%orbs,lin%lb%comms)
call orbitals_communicators(iproc,nproc,Glr,lin%gorbs,lin%gcomms)
call orbitals_communicators(iproc,nproc,Glr,lin%lb%gorbs,lin%lb%gcomms)


! Write all parameters related to the linear scaling version to the screen.
if(iproc==0) call writeLinearParameters(iproc, nproc, at, lin, atomNames, lin%norbsPerType)

! Do some checks on the input parameters.
call checkLinearParameters(iproc, nproc, lin)


! Allocate (almost) all remaining arrays.
call allocateLinArrays(lin)

! Decide which orbital is centered on which atom, again for the 'normal' and
! the 'large' basis.

! This is the same as above, but with orbs%inWhichLocreg instead of lin%onWhichAtom
! The array inWhichLocreg has already been allocated in orbitals_descriptors. Since it will again be allocated
! in assignToLocreg2, deallocate it first.
iall=-product(shape(lin%orbs%inWhichLocreg))*kind(lin%orbs%inWhichLocreg)
deallocate(lin%orbs%inWhichLocreg, stat=istat)
call memocc(istat, iall, 'lin%orbs%inWhichLocreg', subname)
iall=-product(shape(lin%lb%orbs%inWhichLocreg))*kind(lin%lb%orbs%inWhichLocreg)
deallocate(lin%lb%orbs%inWhichLocreg, stat=istat)
call memocc(istat, iall, 'lin%lb%orbs%inWhichLocreg', subname)

!call assignToLocreg2(iproc, at%nat, lin%lzd%nlr, input%nspin, norbsPerAtom, lin%orbs)
call assignToLocreg2(iproc, at%nat, lin%lzd%nlr, input%nspin, norbsPerAtom, rxyz, lin%orbs)
if(lin%useDerivativeBasisFunctions) norbsPerAtom=4*norbsPerAtom
!call assignToLocreg2(iproc, at%nat, lin%lb%lzd%nlr, input%nspin, norbsPerAtom, lin%lb%lzd%orbs)
call assignToLocreg2(iproc, at%nat, lin%lzd%nlr, input%nspin, norbsPerAtom, rxyz, lin%lb%orbs)
if(lin%useDerivativeBasisFunctions) norbsPerAtom=norbsPerAtom/4

! Initialize the localization regions.
if(iproc==0) write(*,'(x,a)',advance='no') 'Initializing localization regions... '
t1=mpi_wtime()
call initLocregs(iproc, at%nat, rxyz, lin, input, Glr)
allocate(phi(lin%lb%gorbs%npsidim), stat=istat)
call memocc(istat, phi, 'phi', subname)
allocate(lphi(lin%lb%orbs%npsidim), stat=istat)
call memocc(istat, lphi, 'lphi', subname)

t2=mpi_wtime()
if(iproc==0) write(*,'(a,es9.3,a)') 'done in ',t2-t1,'s.'
npsidim = 0
do iorb=1,lin%orbs%norbp
 ilr=lin%orbs%inwhichlocreg(iorb+lin%orbs%isorb)
 npsidim = npsidim + lin%Lzd%Llr(ilr)%wfd%nvctr_c+7*lin%Lzd%Llr(ilr)%wfd%nvctr_f
end do
lin%lzd%Lpsidimtot = max(npsidim,1)

! The same for the lb type, i.e. with the derivatives.
npsidim = 0
do iorb=1,lin%lb%orbs%norbp
 ilr=lin%lb%orbs%inwhichlocreg(iorb+lin%lb%orbs%isorb)
 npsidim = npsidim + lin%Lzd%Llr(ilr)%wfd%nvctr_c+7*lin%Lzd%Llr(ilr)%wfd%nvctr_f
end do
lin%lzd%Lpsidimtot_der = max(npsidim,1)


! Maybe this could be moved to another subroutine? Or be omitted at all?
allocate(lin%orbs%eval(lin%orbs%norb), stat=istat)
call memocc(istat, lin%orbs%eval, 'lin%orbs%eval', subname)
lin%orbs%eval=-.5d0
allocate(lin%lb%orbs%eval(lin%lb%orbs%norb), stat=istat)
call memocc(istat, lin%lb%orbs%eval, 'lin%lb%orbs%eval', subname)
lin%lb%orbs%eval=-.5d0

! Initialize the coefficients.
call initCoefficients(iproc, orbs, lin, coeff)

! Initialize the parameters for the point to point communication for the
! calculation of the charge density.
if(iproc==0) write(*,'(x,a)',advance='no') 'Initializing communications sumrho... '
t1=mpi_wtime()
call initializeCommsSumrho2(iproc, nproc, nscatterarr, lin, tag)
t2=mpi_wtime()
if(iproc==0) write(*,'(a,es9.3,a)') 'done in ',t2-t1,'s.'
!call allocateCommunicationbufferSumrho(lin%comsr, subname)

! Copy Glr to lin%lzd
call nullify_locreg_descriptors(lin%lzd%Glr)
call copy_locreg_descriptors(Glr, lin%lzd%Glr, subname)
!call nullify_locreg_descriptors(lin%lb%lzd%Glr)
!call copy_locreg_descriptors(Glr, lin%lb%lzd%Glr, subname)

! Copy nlpspd to lin%lzd
call nullify_nonlocal_psp_descriptors(lin%lzd%Gnlpspd)
call copy_nonlocal_psp_descriptors(nlpspd, lin%lzd%Gnlpspd, subname)
!call nullify_nonlocal_psp_descriptors(lin%lb%lzd%Gnlpspd)
!call copy_nonlocal_psp_descriptors(nlpspd, lin%lb%lzd%Gnlpspd, subname)

! Set localnorb
do ilr=1,lin%lzd%nlr
    lin%lzd%Llr(ilr)%localnorb=0
    do iorb=1,lin%orbs%norbp
        if(lin%orbs%inWhichLocregp(iorb)==ilr) then
            lin%lzd%Llr(ilr)%localnorb = lin%lzd%Llr(ilr)%localnorb+1
        end if
    end do
end do

!! The same for the derivatives
!do ilr=1,lin%lzd%nlr
!    lin%lb%lzd%Llr(ilr)%localnorb=0
!    do iorb=1,lin%lb%orbs%norbp
!        if(lin%lb%orbs%inWhichLocregp(iorb)==ilr) then
!            lin%lb%lzd%Llr(ilr)%localnorb = lin%lb%lzd%Llr(ilr)%localnorb+1
!        end if
!    end do
!end do

! Initialize the parameters for the communication for the
! potential.
if(iproc==0) write(*,'(x,a)',advance='no') 'Initializing communications potential... '
t1=mpi_wtime()
call initializeCommunicationPotential(iproc, nproc, nscatterarr, lin%orbs, lin%lzd, lin%comgp, lin%orbs%inWhichLocreg, tag)
!call initializeCommunicationPotential(iproc, nproc, nscatterarr, lin%lb%orbs, lin%lb%lzd, lin%lb%comgp, &
!     lin%lb%lzd%orbs%inWhichLocreg, tag)
call initializeCommunicationPotential(iproc, nproc, nscatterarr, lin%lb%orbs, lin%lzd, lin%lb%comgp, &
     lin%lb%orbs%inWhichLocreg, tag)
t2=mpi_wtime()
if(iproc==0) write(*,'(a,es9.3,a)') 'done in ',t2-t1,'s.'

! Initialize the parameters for the communication for the orthonormalization.
if(iproc==0) write(*,'(x,a)',advance='no') 'Initializing communications orthonormalization... '
t1=mpi_wtime()
call initCommsOrtho(iproc, nproc, lin%lzd, lin%orbs, lin%orbs%inWhichLocreg, input, lin%locregShape, lin%op, lin%comon, tag)
call initCommsOrtho(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%orbs%inWhichLocreg, input, lin%locregShape, lin%lb%op, lin%lb%comon, tag)
t2=mpi_wtime()
if(iproc==0) write(*,'(a,es9.3,a)') 'done in ',t2-t1,'s.'

! Initialize the parameters for the repartitioning of the orbitals.
if(lin%useDerivativeBasisFunctions) call initializeRepartitionOrbitals(iproc, nproc, tag, lin)

! Restart array for the basis functions (only needed if we use the derivative basis functions).
allocate(lin%lphiRestart(lin%orbs%npsidim), stat=istat)
call memocc(istat, lin%lphiRestart, 'lin%lphiRestart', subname)
allocate(lin%lphiold(lin%orbs%npsidim), stat=istat)
call memocc(istat, lin%lphiold, 'lin%lphiold', subname)

! hamold stores the Hamiltonian matrix of the previous iterations.
allocate(lin%hamold(lin%orbs%norb,lin%orbs%norb), stat=istat)
call memocc(istat, lin%hamold, 'lin%hamold', subname)


!! Stores the Hamiltonian in the basis of the localized orbitals
!allocate(lin%hamold(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
!call memocc(istat, lin%hamold, 'lin%hamold', subname)

! Deallocate all local arrays.
iall=-product(shape(atomNames))*kind(atomNames)
deallocate(atomNames, stat=istat)
call memocc(istat, iall, 'atomNames', subname)

iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
deallocate(norbsPerAtom, stat=istat)
call memocc(istat, iall, 'norbsPerAtom', subname)

if(iproc==0) write(*,'(x,a)',advance='no') 'Initializing input guess... '
t1=mpi_wtime()
call initInputguessConfinement(iproc, nproc, at, Glr, input, lin, rxyz, nscatterarr, tag)
t2=mpi_wtime()
if(iproc==0) write(*,'(a,es9.3,a)') 'done in ',t2-t1,'s.'


! Estimate the memory requirements.
call estimateMemory(iproc, nproc, at%nat, lin, nscatterarr)


if(iproc==0) write(*,'(x,a)',advance='no') 'Initializing matrix compression... '
t1=mpi_wtime()
call initMatrixCompression(iproc, nproc, lin%orbs, lin%op, lin%mad)
call initCompressedMatmul3(lin%orbs%norb, lin%mad)
call initMatrixCompression(iproc, nproc, lin%lb%orbs, lin%lb%op, lin%lb%mad)
call initCompressedMatmul3(lin%lb%orbs%norb, lin%lb%mad)
t2=mpi_wtime()
if(iproc==0) write(*,'(a,es9.3,a)') 'done in ',t2-t1,'s.'
!!if(iproc==0) then
!!    do iall=1,lin%mad%nsegmatmul
!!        write(*,'(a,4i8)') 'iall, lin%mad%keyvmatmul(iall), lin%mad%keygmatmul(1,iall), lin%mad%keygmatmul(2,iall)', iall, lin%mad%keyvmatmul(iall), lin%mad%keygmatmul(1,iall), lin%mad%keygmatmul(2,iall)
!!    end do
!!end if


!!! Plot basis functions grid
!!do iorb=1,lin%orbs%norbp
!!    iiorb=lin%orbs%isorb+iorb
!!    ilr=lin%orbs%inWhichLocreg(iiorb)
!!    call plotGrid(iproc, nproc, lin%lb%orbs%norb, lin%orbs%nspinor, input%nspin, iiorb, lin%lzd%llr(ilr), lin%lzd%glr, at, rxyz, input%hx, input%hy, input%hz)
!!end do


end subroutine allocateAndInitializeLinear




subroutine readLinearParameters(iproc, nproc, lin, at, atomNames)
  use module_base
  use module_types
  implicit none
  
  integer,intent(in):: iproc, nproc
  type(linearParameters):: lin
  type(atoms_data),intent(inout):: at
  character(len=20),dimension(at%ntypes):: atomNames
  !integer,dimension(at%ntypes):: norbsPerType
  
  ! Local variables
  integer:: istat, itype, jtype, ierr, iall, iat, npt, ios
  logical:: fileExists, found
  character(len=*),parameter:: subname='readLinearParameters'
  character(len=20):: atomname
  real(8):: pp, lt
  real(8),dimension(:),allocatable:: locradType
  logical,dimension(at%ntypes):: parametersSpecified
  
  allocate(locradType(at%ntypes), stat=istat)
  call memocc(istat, locradType, 'locradType', subname)

  ! Open the input file and read in the parameters.
  inquire(file='input.lin', exist=fileExists)
  if(.not. fileExists) then
      if(iproc==0) write(*,'(x,a)') "ERROR: the file 'input.lin' must be present for the linear &
          & scaling version!"
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if
  open(unit=99, file='input.lin')
  read(99,*) lin%nItBasisFirst, lin%nItBasis, lin%fixBasis
  read(99,*) lin%convCrit
  read(99,*) lin%DIISHistMin, lin%DIISHistMax, lin%alphaDIIS, lin%alphaSD
  read(99,*) lin%startWithSD, lin%startDIIS
  read(99,*) lin%nItPrecond
  read(99,*) lin%getCoeff, lin%locregShape
  read(99,*) lin%blocksize_pdsyev, lin%blocksize_pdgemm
  read(99,*) lin%nproc_pdsyev, lin%nproc_pdgemm
  read(99,*) lin%methTransformOverlap, lin%nItOrtho, lin%convCritOrtho
  read(99,*) lin%correctionOrthoconstraint
  read(99,*) lin%nItCoeff, lin%convCritCoeff
  read(99,*) lin%mixingMethod
  read(99,*) lin%mixHist, lin%nItSCC, lin%alphaMix, lin%convCritMix
  read(99,*) lin%nItOuterSCC, lin%factorFixBasis, lin%minimalFixBasis, lin%convCritMixOut
  read(99,*) lin%useDerivativeBasisFunctions, lin%ConfPotOrder
  read(99,*) lin%nItInguess, lin%memoryForCommunOverlapIG
  read(99,*) lin%plotBasisFunctions
  read(99,*) lin%transformToGlobal
  read(99,*) lin%norbsPerProcIG

  ! This is maybe already done
  allocate(at%rloc(at%ntypes,3), stat=istat)
  call memocc(istat, at%rloc, 'at%rloc', subname)

  ! Now read in the parameters specific for each atom type.
  parametersSpecified=.false.
  do itype=1,at%ntypes
      read(99,*,iostat=ios) atomname, npt, pp, lt
      if(ios/=0) then
          ! The parameters where not specified for all atom types.
          if(iproc==0) then
              write(*,'(x,a)',advance='no') "ERROR: the file 'input.lin' does not contain the parameters&
                       & for the following atom types:"
              do jtype=1,at%ntypes
                  if(.not.parametersSpecified(jtype)) write(*,'(x,a)',advance='no') trim(at%atomnames(jtype))
              end do
          end if
          call mpi_barrier(mpi_comm_world, ierr)
          stop
      end if
      ! The reading was succesful. Check whether this atom type is actually present.
      found=.false.
      do jtype=1,at%ntypes
          if(trim(atomname)==trim(at%atomnames(jtype))) then
              found=.true.
              parametersSpecified(jtype)=.true.
              atomNames(jtype)=atomname
              lin%norbsPerType(jtype)=npt
              lin%potentialPrefac(jtype)=pp
              locradType(jtype)=lt
              at%rloc(jtype,:)=locradType(jtype)
          end if
      end do
      if(.not.found) then
          if(iproc==0) write(*,'(x,3a)') "ERROR: you specified informations about the atomtype '",trim(atomname), &
                     "', which is not present in the file containing the atomic coordinates."
          call mpi_barrier(mpi_comm_world, ierr)
          stop
      end if
  end do
  close(unit=99)
  
  ! Assign the localization radius to each atom.
  do iat=1,at%nat
      itype=at%iatype(iat)
      lin%locrad(iat)=locradType(itype)
  end do
  
  iall=-product(shape(locradType))*kind(locradType)
  deallocate(locradType, stat=istat)
  call memocc(istat, iall, 'locradType', subname)

end subroutine readLinearParameters




subroutine writeLinearParameters(iproc, nproc, at, lin, atomNames, norbsPerType)
!
! Purpose:
! ========
!   Write all parameters concerning the linear scaling version to the screen.
!
! Calling arguments:
! ==================
use module_base
use module_types
implicit none

integer,intent(in):: iproc, nproc
type(atoms_data),intent(in):: at
type(linearParameters),intent(in):: lin
character(len=20),dimension(at%ntypes),intent(in):: atomNames
integer,dimension(at%ntypes),intent(in):: norbsPerType

! Local variables
integer:: itype, jproc, len1, len2, space1, space2
logical:: written
character(len=8):: message1
character(len=14):: message2
character(len=2):: hist


write(*,'(x,a)') '################################# Input parameters #################################'
write(*,'(x,a)') '>>>> General parameters.'
write(*,'(4x,a)') '|           |    number of    |     prefactor for     | localization |'
write(*,'(4x,a)') '| atom type | basis functions | confinement potential |    radius    |'
do itype=1,at%ntypes
    write(*,'(4x,a,4x,a,a,a,a,i0,7x,a,7x,es9.3,6x,a,3x,f8.4,3x,a)') '| ', trim(atomNames(itype)), &
        repeat(' ', 6-len_trim(atomNames(itype))), '|', repeat(' ', 10-ceiling(log10(dble(norbsPerType(itype)+1)+1.d-10))), &
         norbsPerType(itype), '|', lin%potentialPrefac(itype), ' |', lin%locrad(itype), '|'
end do
close(unit=99)
write(*,'(4x,a)') '----------------------------------------------------------------------'
write(*,'(4x,a)') '| mixing | mixing | iterations in | alpha mix | convergence crit. |  iterations   | &
& factor for  |  minimal  | pot diff to |'
write(*,'(4x,a)') '| scheme | method |  in SC cycle  |           |    for mixing     | in outer loop | &
& fixing basis | fix basis | exit outer SCC |'
if(lin%mixHist==0) then
    message1=' linear '
else
    write(hist,'(i2)') lin%mixHist
    message1=' DIIS'//hist//' '
end if
write(*,'(4x,a,2x,a,2x,a,a,a,a,i0,5x,a,x,es9.3,x,a,5x,es9.3,5x,a,a,i0,4x,a,3x,es8.2,3x,a,2x,es9.3,2x,a,3x,es9.3,3x,a)') '|', &
     lin%mixingMethod, '|', message1, '|', repeat(' ', 10-ceiling(log10(dble(lin%nItSCC+1)+1.d-10))), &
     lin%nItSCC, '|', lin%alphaMix, '|', lin%convCritMix, '|', repeat(' ', 10-ceiling(log10(dble(lin%nItOuterSCC+1)+1.d-10))), &
     lin%nItOuterSCC, '|', lin%factorFixBasis, '|', lin%minimalFixBasis, '|', lin%convCritMixOut, '|'
write(*,'(4x,a)') '-------------------------------------------------------------------'
write(*,'(4x,a)') '| use the derivative | order of conf. | iterations in | IG: orbitals | IG: correction  | transform |'
write(*,'(4x,a)') '|  basis functions   |   potential    |  input guess  | per process  | orthoconstraint | to global |'
if(lin%correctionOrthoconstraint==0) then
    message1='  yes   '
else if(lin%correctionOrthoconstraint==1) then
    message1='  no    '
end if
write(*,'(4x,a,8x,l,10x,a,7x,i1,8x,a,a,i0,5x,a,a,i0,6x,a,5x,a,4x,a,5x,l,4x,a)')  '|', lin%useDerivativeBasisFunctions, '|', &
     lin%confPotOrder, '|', repeat(' ', 10-ceiling(log10(dble(lin%nItInguess+1)+1.d-10))), &
     lin%nItInguess, '|', repeat(' ', 8-ceiling(log10(dble(lin%norbsPerProcIG+1)+1.d-10))), lin%norbsPerProcIG, '|', &
     message1, '|', lin%transformToGlobal, '|'
write(*,'(4x,a)') '----------------------------------------------------------------------'
write(*,'(x,a)') '>>>> Parameters for the optimization of the basis functions.'
write(*,'(4x,a)') '| maximal number | convergence | iterations in  | get coef- | plot  |     stop     |'
write(*,'(4x,a)') '|  of iterations |  criterion  | preconditioner | ficients  | basis | optimization |'
write(*,'(4x,a)') '|  first   else  |             |                |           |       |              |'
if(trim(lin%getCoeff)=='diag') then
    !if(trim(lin%diagMethod)=='seq') then
    !    message1='diag seq'
    !else if(trim(lin%diagMethod)=='par') then
    !    message1='diag par'
    !end if
    message1='  diag  '
else if(trim(lin%getCoeff)=='min') then
    message1='   min  '
end if
write(*,'(4x,a,a,i0,3x,a,i0,2x,a,x,es9.3,x,a,a,i0,a,a,a,l,a,2x,es10.3,2x,a)') '| ', &
    repeat(' ', 5-ceiling(log10(dble(lin%nItBasisFirst+1)+1.d-10))), lin%nItBasisFirst, &
    repeat(' ', 5-ceiling(log10(dble(lin%nItBasis+1)+1.d-10))), lin%nItBasis, &
      '| ', lin%convCrit, ' | ', &
      repeat(' ', 8-ceiling(log10(dble(lin%nItPrecond+1)+1.d-10))), lin%nItPrecond, '       | ' , &
      message1, '  |  ', &
      lin%plotBasisFunctions, '   |', lin%fixBasis, '|'
write(*,'(4x,a)') '---------------------------------------------------------------------'
write(*,'(4x,a)') '| DIIS history | alpha DIIS | alpha SD |  start  | allow DIIS | orthonormalization: | transformation |'
write(*,'(4x,a)') '|  min   max   |            |          | with SD |            | nit max   conv crit | of overlap mat |'
if(lin%methTransformOverlap==0) then
    message2='    exact     '
else if(lin%methTransformOverlap==1) then
    message2='taylor appr. 1'
else if(lin%methTransformOverlap==2) then
    message2='taylor appr. 2'
else if(lin%methTransformOverlap==3) then
    message2='taylor appr. 3'
end if
write(*,'(4x,a,a,i0,3x,a,i0,3x,a,2x,es8.2,2x,a,x,es8.2,x,a,l3,a,x,es10.3,a,a,i0,7x,es7.1,2x,a,x,a,x,a)') '|', &
    repeat(' ', 4-ceiling(log10(dble(lin%DIISHistMin+1)+1.d-10))), lin%DIISHistMin, &
    repeat(' ', 3-ceiling(log10(dble(lin%DIISHistMax+1)+1.d-10))), lin%DIISHistMax, ' |', &
    lin%alphaDIIS, '|', lin%alphaSD, '|   ', lin%startWithSD, '    |', lin%startDIIS, ' |', &
    repeat(' ', 5-ceiling(log10(dble(lin%nItOrtho+1)+1.d-10))), lin%nItOrtho, lin%convCritOrtho, '|', message2, '|'
write(*,'(4x,a)') '------------------------------------------------------------------------------------------------------'
write(*,'(x,a)') '>>>> Parameters for the optimization of the coefficients.'
write(*,'(4x,a)') '| maximal number | convergence |'
write(*,'(4x,a)') '|  of iterations |  criterion  |'
write(*,'(4x,a,a,i0,5x,a,x,es9.3,x,a)') '| ', &
    repeat(' ', 9-ceiling(log10(dble(lin%nItCoeff+1)+1.d-10))), lin%nItCoeff, ' | ', lin%convCritCoeff, ' | '
write(*,'(4x,a)') '--------------------------------'
write(*,'(x,a)') '>>>> Performance options'
write(*,'(4x,a)') '| blocksize | blocksize | max proc | max proc | memory for |'
write(*,'(4x,a)') '|  pdsyev   |  pdgemm   |  pdsyev  |  pdgemm  | overlap IG |'
write(*,'(4x,a,a,i0,4x,a,a,i0,4x,a,a,i0,3x,a,a,i0,3x,a,a,i0,4x,a)') '|',repeat(' ', &
    6-ceiling(log10(dble(abs(lin%blocksize_pdgemm)+1)+1.d-10))),&
    lin%blocksize_pdsyev,'|',repeat(' ', 6-ceiling(log10(dble(abs(lin%blocksize_pdgemm)+1)+1.d-10))),lin%blocksize_pdgemm,&
    '|',repeat(' ', 6-ceiling(log10(dble(abs(lin%nproc_pdgemm)+1)+1.d-10))),lin%nproc_pdgemm,'|',&
    repeat(' ', 6-ceiling(log10(dble(abs(lin%nproc_pdgemm)+1)+1.d-10))),lin%nproc_pdgemm, '|',&
    repeat(' ', 8-ceiling(log10(dble(abs(lin%memoryForCommunOverlapIG)+1)+1.d-10))),lin%memoryForCommunOverlapIG, '|'
write(*,'(x,a)') 'lin%locregShape:',lin%locregShape


written=.false.
write(*,'(x,a)') '>>>> Partition of the basis functions among the processes.'
do jproc=1,nproc-1
    if(lin%orbs%norb_par(jproc)<lin%orbs%norb_par(jproc-1)) then
        len1=1+ceiling(log10(dble(jproc-1)+1.d-5))+ceiling(log10(dble(lin%orbs%norb_par(jproc-1)+1.d-5)))
        len2=ceiling(log10(dble(jproc)+1.d-5))+ceiling(log10(dble(nproc-1)+1.d-5))+&
             ceiling(log10(dble(lin%orbs%norb_par(jproc)+1.d-5)))
        if(len1>=len2) then
            space1=1
            space2=1+len1-len2
        else
            space1=1+len2-len1
            space2=1
        end if
        write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
            lin%orbs%norb_par(jproc-1), ' orbitals,', repeat(' ', space1), '|'
        write(*,'(4x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
            lin%orbs%norb_par(jproc),' orbitals.', repeat(' ', space2), '|'
        written=.true.
        exit
    end if
end do
if(.not.written) then
    write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
        ' treat ',lin%orbs%norbp,' orbitals. |'!, &
end if
write(*,'(x,a)') '-----------------------------------------------'


written=.false.
write(*,'(x,a)') '>>>> Partition of the basis functions including the derivatives among the processes.'
do jproc=1,nproc-1
    if(lin%lb%orbs%norb_par(jproc)<lin%lb%orbs%norb_par(jproc-1)) then
        len1=1+ceiling(log10(dble(jproc-1)+1.d-5))+ceiling(log10(dble(lin%lb%orbs%norb_par(jproc-1)+1.d-5)))
        len2=ceiling(log10(dble(jproc)+1.d-5))+ceiling(log10(dble(nproc-1)+1.d-5))+&
             ceiling(log10(dble(lin%lb%orbs%norb_par(jproc)+1.d-5)))
        if(len1>=len2) then
            space1=1
            space2=1+len1-len2
        else
            space1=1+len2-len1
            space2=1
        end if
        write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
            lin%lb%orbs%norb_par(jproc-1), ' orbitals,', repeat(' ', space1), '|'
        write(*,'(4x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
            lin%lb%orbs%norb_par(jproc),' orbitals.', repeat(' ', space2), '|'
        written=.true.
        exit
    end if
end do
if(.not.written) then
    write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
        ' treat ',lin%lb%orbs%norbp,' orbitals. |'!, &
end if
write(*,'(x,a)') '####################################################################################'


end subroutine writeLinearParameters




subroutine assignOrbitalsToAtoms(iproc, orbs, nat, norbsPerAt, onWhichAtom, onWhichAtomAll)
!
! Purpose:
! ========
!   Assigns the orbitals to the atoms, using the array lin%onWhichAtom.
!   If orbital i is centered on atom j, we have lin%onWhichAtom(i)=j.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc         process ID
!     nat           number of atoms
!     norbsPerAt    indicates how many orbitals are centered on each atom.
!  Input / Output arguments
!  ---------------------
!     lin           type containing parameters for the linear version
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nat
type(orbitals_data):: orbs
integer,dimension(nat):: norbsPerAt
integer,dimension(orbs%norbp):: onWhichAtom
integer,dimension(orbs%norb):: onWhichAtomAll

! Local variables
integer:: jproc, iiOrb, iorb, jorb, jat


  ! There are four counters:
  !   jproc: indicates which MPI process is handling the basis function which is being treated
  !   jat: counts the atom numbers
  !   jorb: counts the orbitals handled by a given process
  !   iiOrb: counts the number of orbitals for a given atoms thas has already been assigned
  jproc=0
  jat=1
  jorb=0
  iiOrb=0
  
  do iorb=1,orbs%norb
  
      ! Switch to the next MPI process if the numbers of orbitals for a given
      ! MPI process is reached.
      if(jorb==orbs%norb_par(jproc)) then
          jproc=jproc+1
          jorb=0
      end if
      
      ! Switch to the next atom if the number of basis functions for this atom is reached.
      if(iiOrb==norbsPerAt(jat)) then
          jat=jat+1
          iiOrb=0
      end if
      jorb=jorb+1
      iiOrb=iiOrb+1
      if(iproc==jproc) onWhichAtom(jorb)=jat

      ! Global assignment, i.e. without taking in account
      ! the various MPI processes.
      onWhichAtomAll(iorb)=jat
  end do    

end subroutine assignOrbitalsToAtoms



!!subroutine assignOrbitalsToAtoms2(iproc, orbs, nat, norbsPerAt, onWhichAtom, onWhichAtomAll)
!!!
!!! Purpose:
!!! ========
!!!   Assigns the orbitals to the atoms, using the array lin%onWhichAtom.
!!!   If orbital i is centered on atom j, we have lin%onWhichAtom(i)=j.
!!!
!!! Calling arguments:
!!! ==================
!!!   Input arguments:
!!!   ----------------
!!!     iproc         process ID
!!!     nat           number of atoms
!!!     norbsPerAt    indicates how many orbitals are centered on each atom.
!!!  Input / Output arguments
!!!  ---------------------
!!!     lin           type containing parameters for the linear version
!!!
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nat
!!type(orbitals_data):: orbs
!!integer,dimension(nat):: norbsPerAt
!!integer,dimension(orbs%norbp):: onWhichAtom
!!integer,dimension(orbs%norb):: onWhichAtomAll
!!
!!! Local variables
!!integer:: jproc, iiOrb, iorb, jorb, jat
!!
!!
!!  ! Count the total 'weighted orbitals', i.e. the sum of 1/(number of orbitals for this atom).
!!  weightLocreg=0.d0
!!  do iat=1,nat
!!      weightLocreg=weightLocreg+1.d0/dble(norbsPerAt(iat))
!!  end do
!!
!!
!!  ! There are four counters:
!!  !   jproc: indicates which MPI process is handling the basis function which is being treated
!!  !   jat: counts the atom numbers
!!  !   jorb: counts the orbitals handled by a given process
!!  !   iiOrb: counts the number of orbitals for a given atoms thas has already been assigned
!!  jproc=0
!!  jat=1
!!  jorb=0
!!  iiOrb=0
!!  
!!  do iorb=1,orbs%norb
!!  
!!      ! Switch to the next MPI process if the numbers of orbitals for a given
!!      ! MPI process is reached.
!!      if(jorb==orbs%norb_par(jproc)) then
!!          jproc=jproc+1
!!          jorb=0
!!      end if
!!      
!!      ! Get the atom number. First calculate the 'weighted localization region, i.e. 1/(number of basis functions for this atom).
!!      ! Switch to the next atom if the number of basis functions for this atom is reached.
!!      if(iiOrb==norbsPerAt(jat)) then
!!          jat=jat+1
!!          iiOrb=0
!!      end if
!!      jorb=jorb+1
!!      iiOrb=iiOrb+1
!!      if(iproc==jproc) onWhichAtom(jorb)=jat
!!
!!      ! Global assignment, i.e. without taking in account
!!      ! the various MPI processes.
!!      onWhichAtomAll(iorb)=jat
!!  end do    
!!
!!end subroutine assignOrbitalsToAtoms2





subroutine checkLinearParameters(iproc, nproc, lin)
!
! Purpose:
! ========
!  Checks some values contained in the variable lin on errors.
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(in):: lin

! Local variables
integer:: norbTarget, nprocIG, ierr


  if(lin%DIISHistMin>lin%DIISHistMax) then
      if(iproc==0) write(*,'(x,a,i0,a,i0,a)') 'ERROR: DIISHistMin must not be larger than &
      & DIISHistMax, but you chose ', lin%DIISHistMin, ' and ', lin%DIISHistMax, '!'
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if

  if(trim(lin%getCoeff)/='min' .and. trim(lin%getCoeff)/='diag') then
      if(iproc==0) write(*,'(x,a,a,a)') "ERROR: lin%getCoeff can have the values 'diag' or 'min', &
          & but we found '", trim(lin%getCoeff), "'!"
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if

  if(lin%methTransformOverlap<0 .or. lin%methTransformOverlap>2) then
      if(iproc==0) write(*,'(x,a,i0,a)') 'ERROR: lin%methTransformOverlap must be 0,1 or 2, but you specified ', &
                               lin%methTransformOverlap,'.'
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if

  !!if(trim(lin%getCoeff)=='diag') then
  !!    if(trim(lin%diagMethod)/='seq' .and. trim(lin%diagMethod)/='par') then
  !!        if(iproc==0) write(*,'(x,a,a,a)') "ERROR: lin%diagMethod can have the values 'seq' or 'par', &
  !!            & but we found '", trim(lin%diagMethod), "'!"
  !!        call mpi_barrier(mpi_comm_world, ierr)
  !!        stop
  !!    end if
  !!end if

  if(lin%confPotOrder/=4 .and. lin%confPotOrder/=6) then
      if(iproc==0) write(*,'(x,a,i0,a)') 'ERROR: lin%confPotOrder can have the values 4 or 6, &
          & but we found ', lin%confPotOrder, '!'
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if


  ! Determine the number of processes we need for the minimization of the trace in the input guess.
  if(lin%norbsPerProcIG>lin%orbs%norb) then
      norbTarget=lin%orbs%norb
  else
      norbTarget=lin%norbsperProcIG
  end if
  nprocIG=ceiling(dble(lin%orbs%norb)/dble(norbTarget))
  nprocIG=min(nprocIG,nproc)

  if( nprocIG/=nproc .and. ((lin%methTransformOverlap==0 .and. (lin%blocksize_pdsyev>0 .or. lin%blocksize_pdgemm>0)) .or. &
      (lin%methTransformOverlap==1 .and. lin%blocksize_pdgemm>0)) ) then
      if(iproc==0) then
          write(*,'(x,a)') 'ERROR: You want to use some routines from scalapack. This is only possible if all processes are &
                     &involved in these calls, which is not the case here.'
          write(*,'(x,a)') 'To avoid this problem you have several possibilities:'
          write(*,'(3x,a,i0,a)') "-set 'lin%norbsperProcIG' to a value not greater than ",floor(dble(lin%orbs%norb)/dble(nproc)), &
              ' (recommended; probably only little influence on performance)'
          write(*,'(3x,a)') "-if you use 'lin%methTransformOverlap==1': set 'lin%blocksize_pdgemm' to a negative value &
              &(may heavily affect performance)"
          write(*,'(3x,a)') "-if you use 'lin%methTransformOverlap==0': set 'lin%blocksize_pdsyev' and 'lin%blocksize_pdsyev' &
              &to negative values (may very heavily affect performance)"
      end if
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if

  if(lin%nproc_pdsyev>nproc) then
      if(iproc==0) write(*,'(x,a)') 'ERROR: lin%nproc_pdsyev can not be larger than nproc'
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if

  if(lin%nproc_pdgemm>nproc) then
      if(iproc==0) write(*,'(x,a)') 'ERROR: lin%nproc_pdgemm can not be larger than nproc'
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if

  if(lin%locregShape/='c' .and. lin%locregShape/='s') then
      if(iproc==0) write(*,*) "ERROR: lin%locregShape must be 's' or 'c'!"
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if



end subroutine checkLinearParameters






subroutine deallocateLinear(iproc, lin, phi, lphi, coeff)
!
! Purpose:
! ========
!   Deallocates all array related to the linear scaling version which have not been 
!   deallocated so far.
!
! Calling arguments:
! ==================
!
use module_base
use module_types
use module_interfaces, exceptThisOne => deallocateLinear
implicit none

! Calling arguments
integer,intent(in):: iproc
type(linearParameters),intent(inout):: lin
real(8),dimension(:),pointer,intent(inout):: phi, lphi
real(8),dimension(:,:),pointer,intent(inout):: coeff

! Local variables
integer:: istat, iall, iorb
character(len=*),parameter:: subname='deallocateLinear'


  iall=-product(shape(phi))*kind(phi)
  deallocate(phi, stat=istat)
  call memocc(istat, iall, 'phi', subname)

  iall=-product(shape(lphi))*kind(lphi)
  deallocate(lphi, stat=istat)
  call memocc(istat, iall, 'lphi', subname)

  iall=-product(shape(coeff))*kind(coeff)
  deallocate(coeff, stat=istat)
  call memocc(istat, iall, 'coeff', subname)

  call deallocate_linearParameters(lin, subname)


end subroutine deallocateLinear





subroutine randomWithinCutoff(iproc, orbs, Glr, at, lin, input, rxyz, phi)
!
! Purpose:
! ========
!   Initializes the basis functions phi to random number within a cutoff range around
!   the atom at which they are centered.
!   The cutoff radius id given by 
!     cut=1.d0/lin%potentialPrefac(at%iatype(iiAt))
!     cut=cut**.25d0
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc
type(orbitals_data), intent(inout) :: orbs
type(locreg_descriptors), intent(in) :: Glr
type(atoms_data),intent(in):: at
type(linearParameters),intent(in):: lin
type(input_variables), intent(in):: input
real(8),dimension(3,at%nat):: rxyz
real(8),dimension((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp):: phi

integer:: ix, iy, iz, ix0, iy0, iz0, iiAt, jj, iorb, i1, i2, i3, istart, ii, istat, iall
real(8),dimension(:),allocatable:: phir
real(8):: hx, hy, hz, hxh, hyh, hzh, kx, ky, kz, tt, tt2, cut
real :: ttreal
type(workarr_sumrho) :: w
type(workarr_locham):: w_lh
character(len=*),parameter:: subname='randomWithinCutoff'


allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)
call memocc(istat, phir, 'phir', subname)
phi=0.d0

call initialize_work_arrays_sumrho(Glr,w)
call initialize_work_arrays_locham(Glr, orbs%nspinor,w_lh)

hx=input%hx
hy=input%hy
hz=input%hz
hxh=.5d0*hx
hyh=.5d0*hy
hzh=.5d0*hz

! Initialize phi to zero.
phi=0.d0

istart=0

!!call the random number as many times as the number of orbitals before
!!so that to associate unambiguously a random number to a  component-orbital pair
do ii=1,orbs%isorb*orbs%nspinor*(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i)
   call random_number(ttreal)
end do

    orbLoop: do iorb=1,orbs%norbp
        call daub_to_isf(Glr,w,phi(istart+1),phir(1))
        !iiAt=lin%onWhichAtom(iorb)
        iiAt=lin%orbs%inWhichLocregp(iorb)
        ix0=nint(rxyz(1,iiAt)/hxh)
        iy0=nint(rxyz(2,iiAt)/hyh)
        iz0=nint(rxyz(3,iiAt)/hzh)
        cut=1.d0/lin%potentialPrefac(at%iatype(iiAt))
        cut=cut**(1.d0/dble(lin%confPotOrder))
        !cut=cut**.166666d0
!cut=80000.d0

        jj=0
        do i3=-14,Glr%d%n3i-15
            do i2=-14,Glr%d%n2i-15
                do i1=-14,Glr%d%n1i-15
                  jj=jj+1

                   tt=hxh**2*(i1-ix0)**2 + hyh**2*(i2-iy0)**2 + hzh**2*(i3-iz0)**2
                   tt=sqrt(tt)
                   if(tt<cut) then
                      call random_number(ttreal)
                      phir(jj)=real(ttreal,kind=8)
                   else
                      !call random_number(ttreal)
                      phir(jj)=0.d0
                   end if
                end do
            end do
        end do

        kx=orbs%kpts(1,orbs%iokpt(iorb))
        ky=orbs%kpts(2,orbs%iokpt(iorb))
        kz=orbs%kpts(3,orbs%iokpt(iorb))
        call isf_to_daub(Glr, w, phir(1), phi(istart+1))

        istart=istart+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor


    end do orbLoop

! Deallocate everything.
call deallocate_work_arrays_sumrho(w)
call deallocate_work_arrays_locham(Glr, w_lh)
iall=-product(shape(phir))*kind(phir)
deallocate(phir, stat=istat)
call memocc(istat, iall, 'phir', subname)


end subroutine randomWithinCutoff











subroutine plotOrbitals(iproc, orbs, Glr, phi, nat, rxyz, onWhichAtom, hxh, hyh, hzh, it)
!
! Plots the orbitals
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc
type(orbitals_data), intent(inout) :: orbs
type(locreg_descriptors), intent(in) :: Glr
real(8),dimension((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp):: phi
integer:: nat
real(8),dimension(3,nat):: rxyz
integer,dimension(orbs%norbp):: onWhichAtom
real(8):: hxh, hyh, hzh
integer:: it

integer:: ix, iy, iz, ix0, iy0, iz0, iiAt, jj, iorb, i1, i2, i3, istart, ii, istat
integer:: unit1, unit2, unit3
real(8),dimension(:),allocatable:: phir
type(workarr_sumrho) :: w
character(len=10):: c1, c2, c3
character(len=50):: file1, file2, file3

allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)

call initialize_work_arrays_sumrho(Glr,w)

istart=0

unit1=10*iproc+7
unit2=10*iproc+8
unit3=10*iproc+9

!write(*,*) 'write, orbs%nbasisp', orbs%norbp
    orbLoop: do iorb=1,orbs%norbp
        phir=0.d0
        call daub_to_isf(Glr,w,phi(istart+1),phir(1))
        iiAt=onWhichAtom(iorb)
        ix0=nint(rxyz(1,iiAt)/hxh)
        iy0=nint(rxyz(2,iiAt)/hyh)
        iz0=nint(rxyz(3,iiAt)/hzh)

        jj=0
        write(c1,'(i0)') iproc
        write(c2,'(i0)') iorb
        write(c3,'(i0)') it
        file1='orbs_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_x'
        file2='orbs_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_y'
        file3='orbs_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_z'
        open(unit=unit1, file=trim(file1))
        open(unit=unit2, file=trim(file2))
        open(unit=unit3, file=trim(file3))
        do i3=1,Glr%d%n3i
            do i2=1,Glr%d%n2i
                do i1=1,Glr%d%n1i
                   jj=jj+1
                   ! z component of point jj
                   iz=jj/(Glr%d%n2i*Glr%d%n1i)
                   ! Subtract the 'lower' xy layers
                   ii=jj-iz*(Glr%d%n2i*Glr%d%n1i)
                   ! y component of point jj
                   iy=ii/Glr%d%n1i
                   ! Subtract the 'lower' y rows
                   ii=ii-iy*Glr%d%n1i
                   ! x component
                   ix=ii
!if(phir(jj)>1.d0) write(*,'(a,3i7,es15.6)') 'WARNING: ix, iy, iz, phir(jj)', ix, iy, iz, phir(jj)
                   if(iy==ix0 .and. iz==iz0) write(unit1,*) ix, phir(jj)
                   ! Write along y-axis
                   if(ix==ix0 .and. iz==iz0) write(unit2,*) iy, phir(jj)
                   ! Write along z-axis
                   if(ix==ix0 .and. iy==iy0) write(unit3,*) iz, phir(jj)


                end do
            end do
        end do
        close(unit=unit1)
        close(unit=unit2)
        close(unit=unit3)

        istart=istart+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor

    end do orbLoop

call deallocate_work_arrays_sumrho(w)
deallocate(phir, stat=istat)


end subroutine plotOrbitals




subroutine cutoffOutsideLocreg(iproc, nproc, Glr, at, input, lin, rxyz, phi)
! Cut off everything outside the localization region by setting it to zero.
! Then do a orthonormalization.
use module_base
use module_types
use module_interfaces
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(locreg_descriptors),intent(in) :: Glr
type(atoms_data),intent(in):: at
type(input_variables),intent(in):: input
type(linearParameters),intent(in):: lin
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(lin%gorbs%npsidim),intent(inout):: phi

! Local variables
integer:: iorb, ist, i1, i2, i3, jj, iiAt, istat, iall, ierr
real(8):: tt, cut, hxh, hyh, hzh, ttIn, ttOut, ttIntot, ttOuttot
type(workarr_sumrho) :: w
real(8),dimension(:),allocatable:: phir
real(8),dimension(:),pointer:: phiWork
character(len=*),parameter:: subname='cutoffOutsideLocreg'

!write(*,*) 'in cutoffOutsideLocreg'

call initialize_work_arrays_sumrho(Glr, w)
hxh=input%hx*.5d0
hyh=input%hy*.5d0
hzh=input%hz*.5d0


allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)
call memocc(istat, phir, 'phir', subname)

ist=1
ttIntot=0.d0
ttOuttot=0.d0
do iorb=1,lin%orbs%norbp
    ! Transform the orbitals to real space.
    phir=0.d0
    call daub_to_isf(Glr, w, phi(ist), phir(1))
    
    !iiAt=lin%onWhichAtom(iorb)
    iiAt=lin%orbs%inWhichLocregp(iorb)
    cut=lin%locrad(iiAt)
    
    jj=0
    ttIn=0.d0
    ttOut=0.d0
    do i3=-14,Glr%d%n3i-15
        do i2=-14,Glr%d%n2i-15
            do i1=-14,Glr%d%n1i-15
               jj=jj+1
               tt = (hxh*i1-rxyz(1,iiAt))**2 + (hyh*i2-rxyz(2,iiAt))**2 + (hzh*i3-rxyz(3,iiAt))**2
               tt=sqrt(tt)
               if(tt>cut) then
                  !write(*,'(a,4i7,3es20.12)') 'iorb, i1, i2, i3, tt, cut, phir(jj)', iorb, i1, i2, i3, tt, cut, phir(jj)
                  ttOut=ttOut+phir(jj)**2
                  phir(jj)=0.d0
               else
                  ttIn=ttIn+phir(jj)**2
               end if
            end do
        end do
    end do
    
    call isf_to_daub(Glr, w, phir(1), phi(ist))

    !write(*,'(a,i7,2es20.12)') 'before: iorb, ttIn, ttOut', iorb, ttIn, ttOut
    ist=ist+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)

    ttIntot = ttIntot + ttIn
    ttOuttot = ttOuttot + ttOut

end do

call mpiallred(ttIntot, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(ttOuttot, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(x,a)') 'cutting of outside localization region:'
if(iproc==0) write(*,'(3x,a,2es17.8)') 'before cut; average weights in / out:', ttIntot/dble(lin%orbs%norb),&
             ttOuttot/dble(lin%orbs%norb)


call mpi_barrier(mpi_comm_world, ierr)
allocate(phiWork(lin%gorbs%npsidim), stat=istat)
call memocc(istat, phiWork, 'phiWork', subname)
call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
call orthogonalize(iproc, nproc, lin%orbs, lin%comms, Glr%wfd, phi, input)
call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
iall=-product(shape(phiWork))*kind(phiWork)
deallocate(phiWork, stat=istat)
call memocc(istat, iall, 'phiWork', subname)

! Check
ist=1
ttIntot=0.d0
ttOuttot=0.d0
do iorb=1,lin%orbs%norbp
    ! Transform the orbitals to real space.
    phir=0.d0
    call daub_to_isf(Glr, w, phi(ist), phir(1))
    
    !iiAt=lin%onWhichAtom(iorb)
    iiAt=lin%orbs%inWhichLocregp(iorb)
    cut=lin%locrad(iiAt)
    !write(*,'(a,2i8,es10.3)') 'iorb, iiAt, cut', iorb, iiAt, cut
    
    jj=0
    ttIn=0.d0
    ttOut=0.d0
    do i3=-14,Glr%d%n3i-15
        do i2=-14,Glr%d%n2i-15
            do i1=-14,Glr%d%n1i-15
               jj=jj+1
               tt = (hxh*i1-rxyz(1,iiAt))**2 + (hyh*i2-rxyz(2,iiAt))**2 + (hzh*i3-rxyz(3,iiAt))**2
               tt=sqrt(tt)
               if(tt>cut) then
                  ttOut = ttOut + phir(jj)**2
               else
                  ttIn = ttIn + phir(jj)**2
               end if
            end do
        end do
    end do
    
    !write(*,'(a,i7,2es20.12)') 'after: iorb, ttIn, ttOut', iorb, ttIn, ttOut
    ist=ist+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)

    ttIntot = ttIntot + ttIn
    ttOuttot = ttOuttot + ttOut

end do

call mpiallred(ttIntot, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(ttOuttot, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(3x,a,2es17.8)') 'after cut; average weights in / out:', ttIntot/dble(lin%orbs%norb),&
            ttOuttot/dble(lin%orbs%norb)

iall=-product(shape(phir))*kind(phir)
deallocate(phir, stat=istat)
call memocc(istat, iall, 'phir', subname)

call deallocate_work_arrays_sumrho(w)

end subroutine cutoffOutsideLocreg






subroutine initializeCommsSumrho2(iproc, nproc, nscatterarr, lin, tag)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
type(linearParameters),intent(inout):: lin
integer,intent(inout):: tag

! Local variables
integer:: istat, jproc, is, ie, ioverlap, i3s, i3e, ilr, iorb, is3ovrlp, n3ovrlp
character(len=*),parameter:: subname='initializeCommsSumrho'


! First count the number of overlapping orbitals for each slice.
allocate(lin%comsr%noverlaps(0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%noverlaps, 'lin%comsr%noverlaps', subname)
do jproc=0,nproc-1
    is=nscatterarr(jproc,3)-14
    ie=is+nscatterarr(jproc,1)-1
    !if(iproc==0) write(*,'(a,3i8)') 'jproc, is, ie', jproc, is, ie
    ioverlap=0
    do iorb=1,lin%lb%orbs%norb
        ilr=lin%lb%orbs%inWhichLocreg(iorb)
        i3s=2*lin%lzd%Llr(ilr)%ns3-14
        i3e=i3s+lin%lzd%Llr(ilr)%d%n3i-1
        if(i3s<=ie .and. i3e>=is) then
            ioverlap=ioverlap+1
        end if
    end do
    lin%comsr%noverlaps(jproc)=ioverlap
    !if(iproc==0) write(*,'(a,2i8)') 'jproc, lin%comsr%noverlaps(jproc)', jproc, lin%comsr%noverlaps(jproc)
end do
! Do the initialization concerning the calculation of the charge density.
allocate(lin%comsr%istarr(0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%istarr, 'lin%comsr%istarr', subname)
!allocate(lin%comsr%istrarr(lin%comsr%noverlaps(iproc)), stat=istat)
allocate(lin%comsr%istrarr(0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%istrarr, 'lin%comsr%istrarr', subname)
allocate(lin%comsr%overlaps(lin%comsr%noverlaps(iproc)), stat=istat)
call memocc(istat, lin%comsr%overlaps, 'lin%comsr%overlaps', subname)

allocate(lin%comsr%comarr(9,maxval(lin%comsr%noverlaps),0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%comarr, 'lin%coms%commsSumrho', subname)


lin%comsr%istarr=1
lin%comsr%istrarr=1
lin%comsr%nrecvBuf=0
do jproc=0,nproc-1
    is=nscatterarr(jproc,3)-14
    ie=is+nscatterarr(jproc,1)-1
    ioverlap=0
    do iorb=1,lin%lb%orbs%norb
        ilr=lin%lb%orbs%inWhichLocreg(iorb)
        i3s=2*lin%lzd%Llr(ilr)%ns3-14
        i3e=i3s+lin%lzd%Llr(ilr)%d%n3i-1
        if(i3s<=ie .and. i3e>=is) then
            ioverlap=ioverlap+1
            tag=tag+1
            is3ovrlp=max(is,i3s) !start of overlapping zone in z direction
            n3ovrlp=min(ie,i3e)-max(is,i3s)+1  !extent of overlapping zone in z direction
            is3ovrlp=is3ovrlp-2*lin%lzd%Llr(ilr)%ns3+15
            !call setCommunicationInformation2(jproc, iorb, is3ovrlp, n3ovrlp, lin%comsr%istrarr(jproc), tag, lin, lin%comsr%comarr(1,ioverlap,jproc))
            call setCommunicationInformation2(jproc, iorb, is3ovrlp, n3ovrlp, lin%comsr%istrarr(jproc), &
                 tag, lin%nlr, lin%lzd%Llr,&
                 lin%lb%orbs%inWhichLocreg, lin%lb%orbs, lin%comsr%comarr(1,ioverlap,jproc))
            if(iproc==jproc) then
                !lin%comsr%sizePhibuffr = lin%comsr%sizePhibuffr + lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*n3ovrlp
                lin%comsr%nrecvBuf = lin%comsr%nrecvBuf + lin%lzd%Llr(ilr)%d%n1i*lin%lzd%Llr(ilr)%d%n2i*n3ovrlp
                lin%comsr%overlaps(ioverlap)=iorb
                                                        !lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*lin%Llr(ilr)%d%n3i
            end if
            lin%comsr%istrarr(jproc) = lin%comsr%istrarr(jproc) + lin%lzd%Llr(ilr)%d%n1i*lin%lzd%Llr(ilr)%d%n2i*n3ovrlp
        end if
    end do
end do

! To avoid allocations with size 0.
lin%comsr%nrecvbuf=max(lin%comsr%nrecvbuf,1)


allocate(lin%comsr%communComplete(maxval(lin%comsr%noverlaps(:)),0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%communComplete, 'lin%comsr%communComplete', subname)
allocate(lin%comsr%computComplete(maxval(lin%comsr%noverlaps(:)),0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%computComplete, 'lin%comsr%computComplete', subname)


! Calculate the dimension of the wave function for each process.
! Do it for both the compressed ('npsidim') and for the uncompressed real space
! ('npsidimr') case.
lin%comsr%nsendBuf=0
do iorb=1,lin%lb%orbs%norbp
    ilr=lin%lb%orbs%inWhichLocregp(iorb)
    lin%comsr%nsendBuf=lin%comsr%nsendBuf+lin%lzd%Llr(ilr)%d%n1i*lin%lzd%Llr(ilr)%d%n2i*lin%lzd%Llr(ilr)%d%n3i*lin%lb%orbs%nspinor
end do

!!allocate(lin%comsr%sendBuf(lin%comsr%nsendBuf), stat=istat)
!!call memocc(istat, lin%comsr%sendBuf, 'lin%comsr%sendBuf', subname)
!!call razero(lin%comsr%nSendBuf, lin%comsr%sendBuf)
!!
!!allocate(lin%comsr%recvBuf(lin%comsr%nrecvBuf), stat=istat)
!!call memocc(istat, lin%comsr%recvBuf, 'lin%comsr%recvBuf', subname)
!!call razero(lin%comsr%nrecvBuf, lin%comsr%recvBuf)

end subroutine initializeCommsSumrho2



subroutine allocateCommunicationbufferSumrho(comsr, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(p2pCommsSumrho),intent(inout):: comsr
character(len=*),intent(in):: subname

! Local variables
integer:: istat

allocate(comsr%sendBuf(comsr%nsendBuf), stat=istat)
call memocc(istat, comsr%sendBuf, 'comsr%sendBuf', subname)
call razero(comsr%nSendBuf, comsr%sendBuf)

allocate(comsr%recvBuf(comsr%nrecvBuf), stat=istat)
call memocc(istat, comsr%recvBuf, 'comsr%recvBuf', subname)
call razero(comsr%nrecvBuf, comsr%recvBuf)

end subroutine allocateCommunicationbufferSumrho


subroutine deallocateCommunicationbufferSumrho(comsr, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(p2pCommsSumrho),intent(inout):: comsr
character(len=*),intent(in):: subname

! Local variables
integer:: istat, iall

iall=-product(shape(comsr%sendBuf))*kind(comsr%sendBuf)
deallocate(comsr%sendBuf, stat=istat)
call memocc(istat, iall, 'comsr%sendBuf', subname)

iall=-product(shape(comsr%recvBuf))*kind(comsr%recvBuf)
deallocate(comsr%recvBuf, stat=istat)
call memocc(istat, iall, 'comsr%recvBuf', subname)


end subroutine deallocateCommunicationbufferSumrho




subroutine allocateLinArrays(lin)
use module_base
use module_types
implicit none

! Calling arguments
type(linearParameters),intent(inout):: lin

! Local variables
integer:: istat
character(len=*),parameter:: subname='allocateLinArrays'


!allocate(lin%onWhichAtom(lin%orbs%norbp), stat=istat)
!call memocc(istat, lin%onWhichAtom, 'lin%onWhichAtom', subname)

!allocate(lin%onWhichAtomAll(lin%orbs%norb), stat=istat)
!call memocc(istat, lin%onWhichAtom, 'lin%onWhichAtomAll', subname)

!allocate(lin%lb%onWhichAtom(lin%lb%orbs%norbp), stat=istat)
!call memocc(istat, lin%lb%onWhichAtom, 'lin%lb%onWhichAtom', subname)

!allocate(lin%lb%onWhichAtomAll(lin%lb%orbs%norb), stat=istat)
!call memocc(istat, lin%lb%onWhichAtom, 'lin%lb%onWhichAtomAll', subname)


end subroutine allocateLinArrays




subroutine allocateBasicArrays(at, lin)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(atoms_data),intent(in):: at
  type(linearParameters),intent(inout):: lin
  
  ! Local variables
  integer:: istat
  character(len=*),parameter:: subname='allocateBasicArrays'
  
  allocate(lin%norbsPerType(at%ntypes), stat=istat)
  call memocc(istat, lin%norbsPerType, 'lin%norbsPerType', subname)
  
  allocate(lin%potentialPrefac(at%ntypes), stat=istat)
  call memocc(istat, lin%potentialPrefac, 'lin%potentialPrefac', subname)

  allocate(lin%locrad(lin%nlr),stat=istat)
  call memocc(istat,lin%locrad,'lin%locrad',subname)

end subroutine allocateBasicArrays



subroutine initLocregs(iproc, nat, rxyz, lin, input, Glr)
use module_base
use module_types
use module_interfaces, exceptThisOne => initLocregs
implicit none

! Calling arguments
integer,intent(in):: iproc, nat
real(8),dimension(3,nat),intent(in):: rxyz
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
type(locreg_descriptors),intent(in):: Glr

! Local variables
integer:: istat, npsidim, npsidimr, iorb, ilr, jorb, jjorb, jlr, iall
character(len=*),parameter:: subname='initLocregs'
logical,dimension(:),allocatable:: calculateBounds

! Allocate the array of localisation regions
!allocate(lin%Llr(lin%nlr),stat=istat)
allocate(lin%lzd%Llr(lin%lzd%nlr),stat=istat)
do ilr=1,lin%lzd%nlr
    call nullify_locreg_descriptors(lin%lzd%llr(ilr))
end do

 allocate(calculateBounds(lin%lzd%nlr), stat=istat)
 call memocc(istat, calculateBounds, 'calculateBounds', subname)
 calculateBounds=.false.
 do ilr=1,lin%lzd%nlr
     do jorb=1,lin%orbs%norbp
         jjorb=lin%orbs%isorb+jorb
         jlr=lin%orbs%inWhichLocreg(jjorb)
         if(jlr==ilr) then
             calculateBounds(ilr)=.true.
             exit
         end if
     end do
     do jorb=1,lin%lb%orbs%norbp
         jjorb=lin%lb%orbs%isorb+jorb
         jlr=lin%lb%orbs%inWhichLocreg(jjorb)
         if(jlr==ilr) then
             calculateBounds(ilr)=.true.
             exit
         end if
     end do
     lin%lzd%llr(ilr)%locrad=lin%locrad(ilr)
     lin%lzd%llr(ilr)%locregCenter=rxyz(:,ilr)
 end do

 if(lin%locregShape=='c') then
     call determine_locreg_periodic(iproc, lin%lzd%nlr, rxyz, lin%locrad, input%hx, input%hy, input%hz, &
          Glr, lin%lzd%Llr, calculateBounds)
 else if(lin%locregShape=='s') then
     call determine_locregSphere(iproc, lin%lzd%nlr, rxyz, lin%locrad, input%hx, input%hy, input%hz, &
          Glr, lin%lzd%Llr, calculateBounds)
 end if
 !call determine_locreg_periodic(iproc, lin%lb%lzd%nlr, rxyz, lin%locrad, input%hx, input%hy, input%hz, Glr, lin%lb%lzd%Llr, calculateBounds)

 if(iproc==0) write(*,'(a,2i9)') 'glr%wfd%nvctr_c, glr%wfd%nvctr_f', glr%wfd%nvctr_c, glr%wfd%nvctr_f
 do ilr=1,lin%lzd%nlr
     if(iproc==0) write(*,'(a,i4,2i9)') 'ilr, nvctrc_c, nvctr_f', ilr, lin%lzd%Llr(ilr)%wfd%nvctr_c, lin%lzd%Llr(ilr)%wfd%nvctr_f
 end do

 iall=-product(shape(calculateBounds))*kind(calculateBounds)
 deallocate(calculateBounds, stat=istat)
 call memocc(istat, iall, 'calculateBounds', subname)


! Calculate the dimension of the wave function for each process.
npsidim=0
do iorb=1,lin%orbs%norbp
    !ilr=lin%onWhichAtom(iorb)
    ilr=lin%orbs%inWhichLocregp(iorb)
    npsidim = npsidim + (lin%lzd%Llr(ilr)%wfd%nvctr_c+7*lin%lzd%Llr(ilr)%wfd%nvctr_f)*lin%orbs%nspinor
end do
!lin%Lorbs%npsidim=npsidim
lin%orbs%npsidim=max(npsidim,1)
lin%orbs%npsidim=max(npsidim,1)

if(.not. lin%useDerivativeBasisFunctions) then
    !lin%lb%Lorbs%npsidim=npsidim
    lin%lb%orbs%npsidim=max(npsidim,1)
    lin%lb%orbs%npsidim=max(npsidim,1)
else
    npsidim=0
    do iorb=1,lin%lb%orbs%norbp
        ilr=lin%lb%orbs%inWhichLocregp(iorb)
        npsidim = npsidim + (lin%lzd%Llr(ilr)%wfd%nvctr_c+7*lin%lzd%Llr(ilr)%wfd%nvctr_f)*lin%lb%orbs%nspinor
        !npsidimr = npsidimr + lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*lin%Llr(ilr)%d%n3i*lin%lb%orbs%nspinor
    end do
    !lin%lb%Lorbs%npsidim=npsidim
    lin%lb%orbs%npsidim=max(npsidim,1)
    lin%lb%orbs%npsidim=max(npsidim,1)
end if


lin%lzd%linear=.true.


!!allocate(phi(lin%lb%gorbs%npsidim), stat=istat)
!!call memocc(istat, phi, 'phi', subname)
!!
!!allocate(lphi(lin%lb%orbs%npsidim), stat=istat)
!!call memocc(istat, lphi, 'lphi', subname)

!allocate(lin%lphiold(lin%lb%orbs%npsidim), stat=istat)
!call memocc(istat, lin%lphiold, 'lin%lphiold', subname)

!allocate(lin%lhphiold(lin%lb%orbs%npsidim), stat=istat)
!call memocc(istat, lin%lhphiold, 'lin%lhphiold', subname)

end subroutine initLocregs







!> Does the same as initLocregs, but has as argumenst lzd instead of lin, i.e. all quantities are
!! are assigned to lzd%Llr etc. instead of lin%Llr. Can probably completely replace initLocregs.
!subroutine initLocregs2(iproc, nat, rxyz, lzd, input, Glr, locrad, phi, lphi)
subroutine initLocregs2(iproc, nat, rxyz, lzd, orbs, input, Glr, locrad, locregShape)
use module_base
use module_types
use module_interfaces, exceptThisOne => initLocregs2
implicit none

! Calling arguments
integer,intent(in):: iproc, nat
real(8),dimension(3,nat),intent(in):: rxyz
type(local_zone_descriptors),intent(inout):: lzd
type(orbitals_data),intent(inout):: orbs
type(input_variables),intent(in):: input
type(locreg_descriptors),intent(in):: Glr
real(8),dimension(lzd%nlr),intent(in):: locrad
character(len=1),intent(in):: locregShape
!real(8),dimension(:),pointer:: phi, lphi

! Local variables
integer:: istat, npsidim, npsidimr, iorb, ilr, jorb, jjorb, jlr, iall
character(len=*),parameter:: subname='initLocregs2'
logical,dimension(:),allocatable:: calculateBounds

! Allocate the array of localisation regions
allocate(lzd%Llr(lzd%nlr),stat=istat)

do ilr=1,lzd%nlr
    call nullify_locreg_descriptors(lzd%Llr(ilr))
end do
!! ATTENTION: WHAT ABOUT OUTOFZONE??


 allocate(calculateBounds(lzd%nlr), stat=istat)
 call memocc(istat, calculateBounds, 'calculateBounds', subname)
 calculateBounds=.false.
 do ilr=1,lzd%nlr
     do jorb=1,orbs%norbp
         jjorb=orbs%isorb+jorb
         jlr=orbs%inWhichLocreg(jjorb)
         if(jlr==ilr) then
             calculateBounds(ilr)=.true.
             exit
         end if
     end do
 end do

 if(locregShape=='c') then
     call determine_locreg_periodic(iproc, lzd%nlr, rxyz, locrad, input%hx, input%hy, input%hz, Glr, lzd%Llr, calculateBounds)
 else if(locregShape=='s') then
     call determine_locregSphere(iproc, lzd%nlr, rxyz, locrad, input%hx, input%hy, input%hz, &
          Glr, lzd%Llr, calculateBounds)
 end if


 iall=-product(shape(calculateBounds))*kind(calculateBounds)
 deallocate(calculateBounds, stat=istat)
 call memocc(istat, iall, 'calculateBounds', subname)

!do ilr=1,lin%nlr
!    if(iproc==0) write(*,'(x,a,i0)') '>>>>>>> zone ', ilr
!    if(iproc==0) write(*,'(3x,a,4i10)') 'nseg_c, nseg_f, nvctr_c, nvctr_f', lin%Llr(ilr)%wfd%nseg_c, lin%Llr(ilr)%wfd%nseg_f, lin%Llr(ilr)%wfd%nvctr_c, lin%Llr(ilr)%wfd%nvctr_f
!    if(iproc==0) write(*,'(3x,a,3i8)') 'lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i', lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i
!    if(iproc==0) write(*,'(a,6i8)') 'lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3',&
!    lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3
!end do


lzd%linear=.true.

! Calculate the dimension of the wave function for each process.
! Do it for both the compressed ('npsidim') and for the uncompressed real space
! ('npsidimr') case.
npsidim=0
do iorb=1,orbs%norbp
    ilr=orbs%inWhichLocregp(iorb)
    npsidim = npsidim + (lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
end do
!! WARNING: CHECHK THIS
orbs%npsidim=max(npsidim,1)


end subroutine initLocregs2


!> Allocate the coefficients for the linear combinations of the  orbitals and initialize
!! them at random.
!! Do this only on the root, since the calculations to determine coeff are not yet parallelized.
subroutine initCoefficients(iproc, orbs, lin, coeff)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc
  type(orbitals_data),intent(in):: orbs
  type(linearParameters),intent(in):: lin
  real(8),dimension(:,:),pointer,intent(out):: coeff
  
  ! Local variables
  integer:: iorb, jorb, istat
  real:: ttreal
  character(len=*),parameter:: subname='initCoefficients'
  
  
  allocate(coeff(lin%lb%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff, 'coeff', subname)
  
  call initRandomSeed(0, 1)
  if(iproc==0) then
      do iorb=1,orbs%norb
         do jorb=1,lin%lb%orbs%norb
            call random_number(ttreal)
            coeff(jorb,iorb)=real(ttreal,kind=8)
         end do
      end do
  end if

end subroutine initCoefficients




!!!!!> Copies lrin to lrout.
!!!! Can be used to copy Glr to lin%lzr%Glr, can probabaly be deleted
!!!! as soon as this initialization is done somewhere else.
!!subroutine allocateAndCopyLocreg(lrin, lrout)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!type(locreg_descriptors),intent(in):: lrin
!!type(locreg_descriptors),intent(out):: lrout
!!
!!
!!lrout%geocode = lrin%geocode
!!
!!lrout%hybrid_on = lrin%hyrbid_on
!!
!!lrout%ns1 = lrin%ns1 ; lrout%ns2 = lrin%ns2 ; lrout%ns3 = lrin%ns3
!!
!!lrout%nsi1 = lrin%nsi1 ; lrout%nsi2 = lrin%nsi2 ; lrout%nsi3 = lrin%nsi3
!!
!!lrout%Localnorb = lrin%Localnorb
!!
!!call dcopy(3, lrin%outofzone(1), 1, lrout%outofzone(1), 1)
!!
!!ii=size(lrin%projflg)
!!allocate(lrout%projflg(ii))
!!call dcopy(ii, lrin%projflg(1), 1, lrout%projflg(1), 1)
!!
!!lrout%Glr = lrin%Glr
!!
!!lrout%Gnlpspd = lrin%nlpspd
!!
!!lrout%orbs = lrin%orbs
!!
!!lrout%comms = lrin%comms
!!
!!
!!
!!!> Contains the information needed for describing completely a
!!!! wavefunction localisation region
!!  type, public :: locreg_descriptors
!!     character(len=1) :: geocode
!!     logical :: hybrid_on               !< interesting for global, periodic, localisation regions
!!     integer :: ns1,ns2,ns3             !< starting point of the localisation region in global coordinates
!!     integer :: nsi1,nsi2,nsi3          !< starting point of locreg for interpolating grid
!!     integer :: Localnorb               !< number of orbitals contained in locreg
!!     integer,dimension(3) :: outofzone  !< vector of points outside of the zone outside Glr for periodic systems
!!     integer,dimension(:),pointer :: projflg    !< atoms contributing nlpsp projectors to locreg
!!     type(grid_dimensions) :: d
!!     type(wavefunctions_descriptors) :: wfd
!!     type(convolutions_bounds) :: bounds
!!  end type locreg_descriptors
!!
!!
!!
!!end subroutine allocateAndCopyLocreg



subroutine estimateMemory(iproc, nproc, nat, lin, nscatterarr)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nat
type(linearParameters),intent(in):: lin
integer,dimension(0:nproc-1,4),intent(in):: nscatterarr

! Local variables
integer,parameter:: nsection=5, narray=12
integer,dimension(narray):: mem
logical,dimension(nsection,narray):: loc
integer:: mempeak, peaksection, isection, iarray, megabytes, memtot, iorb, ilr, ii, iimax
character(len=100),dimension(nsection):: section

if(iproc==0) then
    write(*,'(x,a)') '################################# Memory estimator ##################################'
    write(*,'(x,a)') 'WARNING: The memory requirements are underestimated by about 20-30%!'
    write(*,'(x,a)') 'Memory requirements of the largest arrays:'

    ! For all large arrays determine the memory the occupy and in which code segment they are allcoated.
    ! There are .. segments:
    section(1)='Optimization of the basis functions'
    section(2)='Calculation of the charge density'
    section(3)='Calculate the derivative basis functions'
    section(4)='Calculation the Hamiltonian matrix'
    section(5)='Input guess'

    ! the trace minimizing orbitals:
    mem(1)=8*lin%orbs%npsidim
    loc(1,1)=.true.
    loc(2,1)=.true.
    loc(3,1)=.true.
    loc(4,1)=.true.
    loc(5,1)=.true.
    write(*,'(3x,a,i0,a)') 'trace minimizing orbitals phi: ',megabytes(mem(1)),'MB'

    ! DIIS history of the trace minimizing orbitals
    mem(2)=8*lin%orbs%npsidim*lin%DIISHistMax
    loc(1,2)=.true.
    loc(2,2)=.false.
    loc(3,2)=.false.
    loc(4,2)=.false.
    loc(5,2)=.false.
    write(*,'(3x,a,i0,a)') 'DIIS history of the trace minimizing orbitals phi: ',megabytes(mem(2)),'MB'


    ! The Hamiltonian applied to the orbital, i.e. hphi
    mem(3)=8*lin%orbs%npsidim*lin%DIISHistMax
    loc(1,3)=.false.
    loc(2,3)=.false.
    loc(3,3)=.false.
    loc(4,3)=.true.
    loc(5,3)=.false.
    write(*,'(3x,a,i0,a)') 'The Hamiltonian applied to the orbital, i.e. hphi: ',megabytes(mem(3)),'MB'

    ! charge density / potential (including rhopotold for the mixing or the partial density for the
    ! input guess, therefore times 2)
    mem(4)=8*2*lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
    loc(1,4)=.true.
    loc(2,4)=.true.
    loc(3,4)=.true.
    loc(4,4)=.true.
    loc(5,4)=.true.
    write(*,'(3x,a,i0,a)') 'charge density / potential: ',megabytes(mem(4)),'MB'

    ! send / receive buffers for the charge density
    mem(5)=8*(lin%comsr%nrecvBuf+lin%comsr%nsendBuf)
    loc(1,5)=.false.
    loc(2,5)=.true.
    loc(3,5)=.false.
    loc(4,5)=.true.
    loc(5,5)=.false.
    write(*,'(3x,a,i0,a)') 'communication buffers sumrho: ',megabytes(mem(5)),'MB'

    ! send / receive buffers for the potential (used for the Hamiltonian application)
    mem(6)=8*lin%comgp%nrecvBuf
    loc(1,6)=.true.
    loc(2,6)=.true.
    loc(3,6)=.false.
    loc(4,6)=.true.
    loc(5,6)=.false.
    write(*,'(3x,a,i0,a)') 'communication buffers for gathering the potential: ',megabytes(mem(6)),'MB'

    ! send / receive buffers for the orthonormalization
    mem(7)=8*(lin%comon%nrecvBuf+lin%comon%nsendBuf+lin%op%ndim_lphiovrlp)
    loc(1,7)=.true.
    loc(2,7)=.false.
    loc(3,7)=.false.
    loc(4,7)=.false.
    loc(5,7)=.false.
    write(*,'(3x,a,i0,a)') 'communication buffers / workk arrays for orthonormalization: ',megabytes(mem(7)),'MB'

    ! auxiliary arrays for the orthonormalization (integer arrays)
    mem(8)=4*(lin%comon%nrecvBuf+lin%comon%nsendBuf)
    loc(1,8)=.true.
    loc(2,8)=.false.
    loc(3,8)=.false.
    loc(4,8)=.false.
    loc(5,8)=.false.
    write(*,'(3x,a,i0,a)') 'auxilliary arrays for orthonormalization: ',megabytes(mem(8)),'MB'

    ! full potential need for one localization region during Hamiltonian application and
    ! one orbital in real space (same size)
    iimax=0
    do iorb=1,lin%orbs%norbp
        ilr=lin%orbs%inWhichLocregp(iorb)
        ii=lin%lzd%Llr(ilr)%d%n1i*lin%lzd%Llr(ilr)%d%n2i*lin%lzd%Llr(ilr)%d%n3i
        if(ii>iimax) iimax=ii
    end do
    mem(9)=8*2*iimax
    loc(1,9)=.true.
    loc(2,9)=.false.
    loc(3,9)=.false.
    loc(4,9)=.true.
    loc(5,9)=.false.
    write(*,'(3x,a,i0,a)') 'potential / orbital in real space (Hamiltonian application): ',megabytes(mem(9)),'MB'

    ! Input guess: atomic orbitals (larger cutoff), atomic orbitals (smaller cutoff), hphi for all atoms
    mem(10)=8*( lin%lig%orbsGauss%npsidim + lin%lig%orbsig%npsidim + lin%lig%orbsig%npsidim*nat)
    loc(1,10)=.false.
    loc(2,10)=.false.
    loc(3,10)=.false.
    loc(4,10)=.false.
    loc(5,10)=.true.
    write(*,'(3x,a,i0,a)') 'input guess, all orbitals: ',megabytes(mem(10)),'MB'

    ! Input guess: Buffers for the orthonormalization communication (8 for double precicion arrays
    ! and 4 for single precision arrays -> factor 12)
    mem(11)=12*( lin%lig%comon%nrecvBuf + lin%lig%comon%nsendBuf) + 8*lin%lig%op%ndim_lphiovrlp
    loc(1,11)=.false.
    loc(2,11)=.false.
    loc(3,11)=.false.
    loc(4,11)=.false.
    loc(5,11)=.true.
    write(*,'(3x,a,i0,a)') 'input guess, communication buffers and auxilliary arrays for orthonormalization : ', &
         megabytes(mem(11)),'MB'

    ! Input guess: Buffers for the communicatin the potential
    mem(12)=8*lin%lig%comgp%nrecvBuf
    loc(1,12)=.false.
    loc(2,12)=.false.
    loc(3,12)=.false.
    loc(4,12)=.false.
    loc(5,12)=.true.
    write(*,'(3x,a,i0,a)') 'input guess, communication buffers for gathering the potential: ',megabytes(mem(12)),'MB'


    ! Calculate the memory peak
    mempeak=0
    do isection=1,nsection
        memtot=0
        do iarray=1,narray
            if(loc(isection,iarray)) then
                memtot=memtot+mem(iarray)
            end if
        end do
        if(memtot>mempeak) then
            mempeak=memtot
            peaksection=isection
        end if
    end do
    write(*,'(x,a,i0,a)') '>>> estimated memory peak: ',megabytes(mempeak),'MB'
    write(*,'(x,a,a)') '>>> peak section: ',trim(section(peaksection))

    write(*,'(x,a)') '#####################################################################################'

end if


end subroutine estimateMemory


function megabytes(bytes)
  implicit none
  
  integer,intent(in):: bytes
  integer:: megabytes
  
  megabytes=nint(dble(bytes)/1048576.d0)
  
end function megabytes




subroutine initMatrixCompression(iproc, nproc, orbs, op, mad)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(overlapParameters),intent(in):: op
  type(matrixDescriptors),intent(out):: mad
  
  ! Local variables
  integer:: jproc, iorb, jorb, iiorb, jjorb, ijorb, jjorbold, istat, iseg, nseg, ii, irow, irowold, isegline, i
  character(len=*),parameter:: subname='initMatrixCompression'
  
  
  mad%nseg=0
  mad%nvctr=0
  jjorbold=-1
  irowold=0
  allocate(mad%nsegline(orbs%norb), stat=istat)
  call memocc(istat, mad%nsegline, 'mad%nsegline', subname)
  mad%nsegline=0
  do jproc=0,nproc-1
      do iorb=1,orbs%norb_par(jproc)
          iiorb=orbs%isorb_par(jproc)+iorb
          ijorb=(iiorb-1)*orbs%norb
          do jorb=1,op%noverlaps(iiorb)
              jjorb=op%overlaps(jorb,iiorb)+ijorb
              ! Entry (iiorb,jjorb) is not zero.
              !if(iproc==0) write(300,*) iiorb,jjorb
              if(jjorb==jjorbold+1) then
                  ! There was no zero element in between, i.e. we are in the same segment.
                  jjorbold=jjorb
                  mad%nvctr=mad%nvctr+1

                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  if(irow/=irowold) then
                      ! We are in a new line
                      mad%nsegline(irow)=mad%nsegline(irow)+1
                      irowold=irow
                  end if

              else
                  ! There was a zero segment in between, i.e. we are in a new segment
                  mad%nseg=mad%nseg+1
                  mad%nvctr=mad%nvctr+1
                  jjorbold=jjorb
                  
                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  mad%nsegline(irow)=mad%nsegline(irow)+1
                  irowold=irow
              end if
          end do
      end do
  end do

  !if(iproc==0) write(*,*) 'mad%nseg, mad%nvctr',mad%nseg, mad%nvctr
  mad%nseglinemax=0
  do iorb=1,orbs%norb
      if(mad%nsegline(iorb)>mad%nseglinemax) then
          mad%nseglinemax=mad%nsegline(iorb)
      end if
  end do

  allocate(mad%keyv(mad%nseg), stat=istat)
  call memocc(istat, mad%keyv, 'mad%keyv', subname)
  allocate(mad%keyg(2,mad%nseg), stat=istat)
  call memocc(istat, mad%keyg, 'mad%keyg', subname)
  allocate(mad%keygline(2,mad%nseglinemax,orbs%norb), stat=istat)
  call memocc(istat, mad%keygline, 'mad%keygline', subname)


  nseg=0
  mad%keyv=0
  jjorbold=-1
  irow=0
  isegline=0
  irowold=0
  mad%keygline=0
  mad%keyg=0
  do jproc=0,nproc-1
      do iorb=1,orbs%norb_par(jproc)
          iiorb=orbs%isorb_par(jproc)+iorb
          ijorb=(iiorb-1)*orbs%norb
          do jorb=1,op%noverlaps(iiorb)
              jjorb=op%overlaps(jorb,iiorb)+ijorb
              ! Entry (iiorb,jjorb) is not zero.
              !if(iproc==0) write(300,*) iiorb,jjorb
              if(jjorb==jjorbold+1) then
                  ! There was no zero element in between, i.e. we are in the same segment.
                  mad%keyv(nseg)=mad%keyv(nseg)+1

                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  if(irow/=irowold) then
                      ! We are in a new line, so close the last segment and start the new one
                      mad%keygline(2,isegline,irowold)=mod(jjorbold-1,orbs%norb)+1
                      isegline=1
                      mad%keygline(1,isegline,irow)=mod(jjorb-1,orbs%norb)+1
                      irowold=irow
                  end if
                  jjorbold=jjorb
              else
                  ! There was a zero segment in between, i.e. we are in a new segment.
                  ! First determine the end of the previous segment.
                  if(jjorbold>0) then
                      mad%keyg(2,nseg)=jjorbold
                      mad%keygline(2,isegline,irowold)=mod(jjorbold-1,orbs%norb)+1
                  end if
                  ! Now add the new segment.
                  nseg=nseg+1
                  mad%keyg(1,nseg)=jjorb
                  jjorbold=jjorb
                  mad%keyv(nseg)=mad%keyv(nseg)+1

                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  if(irow/=irowold) then
                      ! We are in a new line
                      isegline=1
                      mad%keygline(1,isegline,irow)=mod(jjorb-1,orbs%norb)+1
                      irowold=irow
                  else
                      ! We are in the same line
                      isegline=isegline+1
                      mad%keygline(1,isegline,irow)=mod(jjorb-1,orbs%norb)+1
                      irowold=irow
                  end if
              end if
          end do
      end do
  end do
  ! Close the last segment
  mad%keyg(2,nseg)=jjorb
  mad%keygline(2,isegline,orbs%norb)=mod(jjorb-1,orbs%norb)+1

  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        write(*,'(a,2x,i0,2x,i0,3x,100i4)') 'iorb, mad%nsegline(iorb), mad%keygline(1,:,iorb)', iorb, mad%nsegline(iorb), mad%keygline(1,:,iorb)
  !!        write(*,'(a,2x,i0,2x,i0,3x,100i4)') 'iorb, mad%nsegline(iorb), mad%keygline(2,:,iorb)', iorb, mad%nsegline(iorb), mad%keygline(2,:,iorb)
  !!    end do
  !!end if

  !!if(iproc==0) then
  !!    do iseg=1,mad%nseg
  !!        write(*,'(a,4i8)') 'iseg, mad%keyv(iseg), mad%keyg(1,iseg), mad%keyg(2,iseg)', iseg, mad%keyv(iseg), mad%keyg(1,iseg), mad%keyg(2,iseg)
  !!    end do
  !!end if

  ! Some checks
  ii=0
  do iseg=1,mad%nseg
      ii=ii+mad%keyv(iseg)
  end do
  if(ii/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR: ii/=mad%nvctr',ii,mad%nvctr
      stop
  end if



end subroutine initMatrixCompression





subroutine initMatrixCompressionForInguess(iproc, nproc, nlr, orbs, noverlaps, overlaps, mad)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, nlr
  type(orbitals_data),intent(in):: orbs
  integer,dimension(nlr),intent(in):: noverlaps
  integer,dimension(maxval(noverlaps(:)),nlr),intent(in):: overlaps
  type(matrixDescriptors),intent(out):: mad
  
  ! Local variables
  integer:: jproc, iorb, jorb, iiorb, jjorb, ijorb, jjorbold, istat, iseg, nseg, ii, irow, irowold, isegline, i, ilr
  character(len=*),parameter:: subname='initMatrixCompression'
  
  
  mad%nseg=0
  mad%nvctr=0
  jjorbold=-1
  irowold=0
  allocate(mad%nsegline(orbs%norb), stat=istat)
  call memocc(istat, mad%nsegline, 'mad%nsegline', subname)
  mad%nsegline=0
  do jproc=0,nproc-1
      do iorb=1,orbs%norb_par(jproc)
          iiorb=orbs%isorb_par(jproc)+iorb
          ilr=orbs%inWhichLocreg(iiorb)
          ijorb=(iiorb-1)*orbs%norb
          !do jorb=1,noverlaps(iiorb)
          do jorb=1,noverlaps(ilr)
              !jjorb=overlaps(jorb,iiorb)+ijorb
              jjorb=overlaps(jorb,ilr)+ijorb
              ! Entry (iiorb,jjorb) is not zero.
              !if(iproc==0) write(300,*) iiorb,jjorb
              if(jjorb==jjorbold+1) then
                  ! There was no zero element in between, i.e. we are in the same segment.
                  jjorbold=jjorb
                  mad%nvctr=mad%nvctr+1

                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  if(irow/=irowold) then
                      ! We are in a new line
                      mad%nsegline(irow)=mad%nsegline(irow)+1
                      irowold=irow
                  end if

              else
                  ! There was a zero segment in between, i.e. we are in a new segment
                  mad%nseg=mad%nseg+1
                  mad%nvctr=mad%nvctr+1
                  jjorbold=jjorb
                  
                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  mad%nsegline(irow)=mad%nsegline(irow)+1
                  irowold=irow
              end if
          end do
      end do
  end do

  !if(iproc==0) write(*,*) 'mad%nseg, mad%nvctr',mad%nseg, mad%nvctr
  mad%nseglinemax=0
  do iorb=1,orbs%norb
      if(mad%nsegline(iorb)>mad%nseglinemax) then
          mad%nseglinemax=mad%nsegline(iorb)
      end if
  end do

  allocate(mad%keyv(mad%nseg), stat=istat)
  call memocc(istat, mad%keyv, 'mad%keyv', subname)
  allocate(mad%keyg(2,mad%nseg), stat=istat)
  call memocc(istat, mad%keyg, 'mad%keyg', subname)
  allocate(mad%keygline(2,mad%nseglinemax,orbs%norb), stat=istat)
  call memocc(istat, mad%keygline, 'mad%keygline', subname)


  nseg=0
  mad%keyv=0
  jjorbold=-1
  irow=0
  isegline=0
  irowold=0
  mad%keygline=0
  mad%keyg=0
  do jproc=0,nproc-1
      do iorb=1,orbs%norb_par(jproc)
          iiorb=orbs%isorb_par(jproc)+iorb
          ilr=orbs%inWhichLocreg(iiorb)
          ijorb=(iiorb-1)*orbs%norb
          !do jorb=1,noverlaps(iiorb)
          do jorb=1,noverlaps(ilr)
              !jjorb=overlaps(jorb,iiorb)+ijorb
              jjorb=overlaps(jorb,ilr)+ijorb
              ! Entry (iiorb,jjorb) is not zero.
              !if(iproc==0) write(300,*) iiorb,jjorb
              if(jjorb==jjorbold+1) then
                  ! There was no zero element in between, i.e. we are in the same segment.
                  mad%keyv(nseg)=mad%keyv(nseg)+1

                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  if(irow/=irowold) then
                      ! We are in a new line, so close the last segment and start the new one
                      mad%keygline(2,isegline,irowold)=mod(jjorbold-1,orbs%norb)+1
                      isegline=1
                      mad%keygline(1,isegline,irow)=mod(jjorb-1,orbs%norb)+1
                      irowold=irow
                  end if
                  jjorbold=jjorb
              else
                  ! There was a zero segment in between, i.e. we are in a new segment.
                  ! First determine the end of the previous segment.
                  if(jjorbold>0) then
                      mad%keyg(2,nseg)=jjorbold
                      mad%keygline(2,isegline,irowold)=mod(jjorbold-1,orbs%norb)+1
                  end if
                  ! Now add the new segment.
                  nseg=nseg+1
                  mad%keyg(1,nseg)=jjorb
                  jjorbold=jjorb
                  mad%keyv(nseg)=mad%keyv(nseg)+1

                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  if(irow/=irowold) then
                      ! We are in a new line
                      isegline=1
                      mad%keygline(1,isegline,irow)=mod(jjorb-1,orbs%norb)+1
                      irowold=irow
                  else
                      ! We are in the same line
                      isegline=isegline+1
                      mad%keygline(1,isegline,irow)=mod(jjorb-1,orbs%norb)+1
                      irowold=irow
                  end if
              end if
          end do
      end do
  end do
  ! Close the last segment
  mad%keyg(2,nseg)=jjorb
  mad%keygline(2,isegline,orbs%norb)=mod(jjorb-1,orbs%norb)+1

  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        write(*,'(a,2x,i0,2x,i0,3x,100i4)') 'iorb, mad%nsegline(iorb), mad%keygline(1,:,iorb)', iorb, mad%nsegline(iorb), mad%keygline(1,:,iorb)
  !!        write(*,'(a,2x,i0,2x,i0,3x,100i4)') 'iorb, mad%nsegline(iorb), mad%keygline(2,:,iorb)', iorb, mad%nsegline(iorb), mad%keygline(2,:,iorb)
  !!    end do
  !!end if

  !!if(iproc==0) then
  !!    do iseg=1,mad%nseg
  !!        write(*,'(a,4i8)') 'iseg, mad%keyv(iseg), mad%keyg(1,iseg), mad%keyg(2,iseg)', iseg, mad%keyv(iseg), mad%keyg(1,iseg), mad%keyg(2,iseg)
  !!    end do
  !!end if

  ! Some checks
  ii=0
  do iseg=1,mad%nseg
      ii=ii+mad%keyv(iseg)
  end do
  if(ii/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR: ii/=mad%nvctr',ii,mad%nvctr
      stop
  end if



end subroutine initMatrixCompressionForInguess






subroutine compressMatrix(norb, mad, mat, lmat)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: norb
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(norb**2),intent(in):: mat
  real(8),dimension(mad%nvctr),intent(out):: lmat
  
  ! Local variables
  integer:: iseg, jj, jorb, iiorb, jjorb
  
  
  jj=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jj=jj+1
          lmat(jj)=mat(jorb)
      end do
  end do
  if(jj/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix: jj/=mad%nvctr',jj,mad%nvctr
      stop
  end if
  
end subroutine compressMatrix



subroutine compressMatrix2(iproc, nproc, orbs, mad, mat, lmat, sendcounts, displs)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb**2),intent(in):: mat
  real(8),dimension(mad%nvctr),intent(out):: lmat
  integer,dimension(0:nproc-1),intent(out):: sendcounts, displs
  
  ! Local variables
  integer:: iseg, jj, jorb, iiorb, jjorb, jjproc, jjprocold, ncount
  
  sendcounts=0
  displs=0
  
  jj=0
  ncount=0
  jjprocold=0
  displs(0)=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jj=jj+1
          lmat(jj)=mat(jorb)
          
          ncount=ncount+1
          jjorb=(jorb-1)/orbs%norb+1
          jjproc=orbs%onWhichMPI(jjorb)
          if(jjproc>jjprocold) then
              ! This part of the matrix is calculated by a new MPI process.
              sendcounts(jjproc-1)=ncount-1
              displs(jjproc)=displs(jjproc-1)+sendcounts(jjproc-1)
              ncount=1
              jjprocold=jjproc
          end if
      end do
  end do
  sendcounts(nproc-1)=ncount
  if(jj/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix: jj/=mad%nvctr',jj,mad%nvctr
      stop
  end if

  if(sum(sendcounts)/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix2: sum(sendcounts)/=mad%nvctr',sum(sendcounts),mad%nvctr
      stop
  end if

  !if(iproc==0) then
  !    do jjproc=0,nproc-1
  !        write(*,'(a,3i8)') 'jjproc, displs(jjproc), sendcounts(jjproc)', jjproc, displs(jjproc), sendcounts(jjproc)
  !    end do
  !end if
  
end subroutine compressMatrix2



subroutine compressMatrixPerProcess(iproc, nproc, orbs, mad, mat, size_lmat, lmat)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, size_lmat
  type(orbitals_data),intent(in):: orbs
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb**2),intent(in):: mat
  real(8),dimension(size_lmat),intent(out):: lmat
  
  ! Local variables
  integer:: iseg, jj, jorb, jjorb, jjproc
  
  
  jj=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jjorb=(jorb-1)/orbs%norb+1
          jjproc=orbs%onWhichMPI(jjorb)
          if(iproc==jjproc) then
              jj=jj+1
              lmat(jj)=mat(jorb)
          end if
      end do
  end do
  
end subroutine compressMatrixPerProcess




subroutine getCommunArraysMatrixCompression(iproc, nproc, orbs, mad, sendcounts, displs)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(matrixDescriptors),intent(in):: mad
  integer,dimension(0:nproc-1),intent(out):: sendcounts, displs
  
  ! Local variables
  integer:: iseg, jj, jorb, iiorb, jjorb, jjproc, jjprocold, ncount
  
  sendcounts=0
  displs=0
  
  jj=0
  ncount=0
  jjprocold=0
  displs(0)=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jj=jj+1
          ncount=ncount+1
          jjorb=(jorb-1)/orbs%norb+1
          jjproc=orbs%onWhichMPI(jjorb)
          if(jjproc>jjprocold) then
              ! This part of the matrix is calculated by a new MPI process.
              sendcounts(jjproc-1)=ncount-1
              displs(jjproc)=displs(jjproc-1)+sendcounts(jjproc-1)
              ncount=1
              jjprocold=jjproc
          end if
      end do
  end do
  sendcounts(nproc-1)=ncount
  if(jj/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix: jj/=mad%nvctr',jj,mad%nvctr
      stop
  end if

  if(sum(sendcounts)/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix2: sum(sendcounts)/=mad%nvctr',sum(sendcounts),mad%nvctr
      stop
  end if

  
end subroutine getCommunArraysMatrixCompression





subroutine initCommsCompression(iproc, nproc, orbs, mad, mat, lmat, sendcounts, displs)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb**2),intent(in):: mat
  real(8),dimension(mad%nvctr),intent(out):: lmat
  integer,dimension(0:nproc-1),intent(out):: sendcounts, displs
  
  ! Local variables
  integer:: iseg, jj, jorb, iiorb, jjorb, jjproc, jjprocold, ncount
  
  sendcounts=0
  displs=0
  
  jj=0
  ncount=0
  jjprocold=0
  displs(0)=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jj=jj+1
          lmat(jj)=mat(jorb)
          
          ncount=ncount+1
          jjorb=(jorb-1)/orbs%norb+1
          jjproc=orbs%onWhichMPI(jjorb)
          if(jjproc>jjprocold) then
              ! This part of the matrix is calculated by a new MPI process.
              sendcounts(jjproc-1)=ncount-1
              displs(jjproc)=displs(jjproc-1)+sendcounts(jjproc-1)
              ncount=1
              jjprocold=jjproc
          end if
      end do
  end do
  sendcounts(nproc-1)=ncount
  if(jj/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix: jj/=mad%nvctr',jj,mad%nvctr
      stop
  end if

  if(sum(sendcounts)/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix2: sum(sendcounts)/=mad%nvctr',sum(sendcounts),mad%nvctr
      stop
  end if

  !if(iproc==0) then
  !    do jjproc=0,nproc-1
  !        write(*,'(a,3i8)') 'jjproc, displs(jjproc), sendcounts(jjproc)', jjproc, displs(jjproc), sendcounts(jjproc)
  !    end do
  !end if
  
end subroutine initCommsCompression



subroutine uncompressMatrix(norb, mad, lmat, mat)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: norb
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(mad%nvctr),intent(in):: lmat
  real(8),dimension(norb**2),intent(out):: mat
  
  ! Local variables
  integer:: iseg, jj, jorb, iiorb, jjorb
  
  mat=0.d0
  
  jj=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jj=jj+1
          mat(jorb)=lmat(jj)
      end do
  end do
  if(jj/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in uncompressMatrix: jj/=mad%nvctr',jj,mad%nvctr
      stop
  end if
  
end subroutine uncompressMatrix




subroutine initCompressedMatmul(iproc, nproc, norb, mad)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb
  type(matrixDescriptors),intent(inout):: mad
  
  ! Local variables
  integer:: iorb, jorb, ii, j, istat, iall, ij, iseg
  logical:: segment
  integer,dimension(:),allocatable:: row, column
  character(len=*),parameter:: subname='initCompressedMatmul'
  
  
  allocate(row(norb), stat=istat)
  call memocc(istat, row, 'row', subname)
  allocate(column(norb), stat=istat)
  call memocc(istat, column, 'column', subname)
  
  
  segment=.false.
  mad%nsegmatmul=0
  mad%nvctrmatmul=0
  do iorb=1,norb
      do jorb=1,norb
          ! Get an array of this line and column indicating whether
          ! there are nonzero numbers at these positions. Since the localization
          ! within the matrix is symmetric, we can use both time the same subroutine.
          call getRow(norb, mad, iorb, row) 
          call getRow(norb, mad, jorb, column) 
          !!if(iproc==0) write(*,'(a,i4,4x,100i4)') 'iorb, row', iorb, row
          !!if(iproc==0) write(*,'(a,i4,4x,100i4)') 'jorb, row', jorb, column
          ii=0
          do j=1,norb
              ii=ii+row(j)*column(j)
          end do
          if(ii>0) then
              ! This entry of the matrix will be different from zero.
              mad%nvctrmatmul=mad%nvctrmatmul+1
              if(.not. segment) then
                  ! This is the start of a new segment
                  segment=.true.
                  mad%nsegmatmul=mad%nsegmatmul+1
              end if
          else
              if(segment) then
                  ! We reached the end of a segment
                  segment=.false.
              end if
          end if
      end do
  end do
  
  allocate(mad%keygmatmul(2,mad%nsegmatmul), stat=istat)
  allocate(mad%keyvmatmul(mad%nsegmatmul), stat=istat)
  
  ! Now fill the descriptors.
  segment=.false.
  ij=0
  iseg=0
  do iorb=1,norb
      do jorb=1,norb
          ij=ij+1
          ! Get an array of this line and column indicating whether
          ! there are nonzero numbers at these positions. Since the localization
          ! within the matrix is symmetric, we can use both time the same subroutine.
          call getRow(norb, mad, iorb, row) 
          call getRow(norb, mad, jorb, column) 
          ii=0
          do j=1,norb
              ii=ii+row(j)*column(j)
          end do
          if(ii>0) then
              ! This entry of the matrix will be different from zero.
              if(.not. segment) then
                  ! This is the start of a new segment
                  segment=.true.
                  iseg=iseg+1
                  mad%keygmatmul(1,iseg)=ij
              end if
              mad%keyvmatmul(iseg)=mad%keyvmatmul(iseg)+1
          else
              if(segment) then
                  ! We reached the end of a segment
                  segment=.false.
                  mad%keygmatmul(2,iseg)=ij-1
              end if
          end if
      end do
  end do

  ! Close the last segment if required.
  if(segment) then
      mad%keygmatmul(2,iseg)=ij
  end if
  
  
  iall=-product(shape(row))*kind(row)
  deallocate(row, stat=istat)
  call memocc(istat, iall, 'row', subname)
  iall=-product(shape(column))*kind(column)
  deallocate(column, stat=istat)
  call memocc(istat, iall, 'column', subname)

end subroutine initCompressedMatmul



subroutine getRow(norb, mad, rowX, row)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: norb, rowX
  type(matrixDescriptors),intent(in):: mad
  integer,dimension(norb),intent(out):: row
  
  ! Local variables
  integer:: iseg, i, irow, icolumn
  
  row=0
  
  do iseg=1,mad%nseg
      do i=mad%keyg(1,iseg),mad%keyg(2,iseg)
      ! Get the row index of this element. Since the localization is symmetric, we can
      ! assume row or column ordering with respect to the segments.
          irow=(i-1)/norb+1
          if(irow==rowX) then
              ! Get the column index of this element.
              icolumn=i-(irow-1)*norb
              row(icolumn)=1
          end if
      end do
  end do

end subroutine getRow






subroutine initCompressedMatmul2(norb, nseg, keyg, nsegmatmul, keygmatmul, keyvmatmul)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: norb, nseg
  integer,dimension(2,nseg),intent(in):: keyg
  integer,intent(out):: nsegmatmul
  integer,dimension(:,:),pointer,intent(out):: keygmatmul
  integer,dimension(:),pointer,intent(out):: keyvmatmul

  ! Local variables
  integer:: iorb, jorb, ii, j, istat, iall, ij, iseg, i
  logical:: segment
  character(len=*),parameter:: subname='initCompressedMatmul2'
  real(8),dimension(:),allocatable:: mat1, mat2, mat3



  allocate(mat1(norb**2), stat=istat)
  call memocc(istat, mat1, 'mat1', subname)
  allocate(mat2(norb**2), stat=istat)
  call memocc(istat, mat2, 'mat2', subname)
  allocate(mat3(norb**2), stat=istat)
  call memocc(istat, mat2, 'mat2', subname)

  mat1=0.d0
  mat2=0.d0
  do iseg=1,nseg
      do i=keyg(1,iseg),keyg(2,iseg)
          ! the localization region is "symmetric"
          mat1(i)=1.d0
          mat2(i)=1.d0
      end do
  end do

  call dgemm('n', 'n', norb, norb, norb, 1.d0, mat1, norb, mat2, norb, 0.d0, mat3, norb)

  segment=.false.
  nsegmatmul=0
  do iorb=1,norb**2
      if(mat3(iorb)>0.d0) then
          ! This entry of the matrix will be different from zero.
          if(.not. segment) then
              ! This is the start of a new segment
              segment=.true.
              nsegmatmul=nsegmatmul+1
          end if
      else
          if(segment) then
              ! We reached the end of a segment
              segment=.false.
          end if
      end if
  end do


  allocate(keygmatmul(2,nsegmatmul), stat=istat)
  call memocc(istat, keygmatmul, 'keygmatmul', subname)
  allocate(keyvmatmul(nsegmatmul), stat=istat)
  call memocc(istat, keyvmatmul, 'keyvmatmul', subname)
  keyvmatmul=0
  ! Now fill the descriptors.
  segment=.false.
  ij=0
  iseg=0
  do iorb=1,norb**2
      ij=iorb
      if(mat3(iorb)>0.d0) then
          ! This entry of the matrix will be different from zero.
          if(.not. segment) then
              ! This is the start of a new segment
              segment=.true.
              iseg=iseg+1
              keygmatmul(1,iseg)=ij
          end if
          keyvmatmul(iseg)=keyvmatmul(iseg)+1
      else
          if(segment) then
              ! We reached the end of a segment
              segment=.false.
              keygmatmul(2,iseg)=ij-1
          end if
      end if
  end do
  ! Close the last segment if required.
  if(segment) then
      keygmatmul(2,iseg)=ij
  end if


iall=-product(shape(mat1))*kind(mat1)
deallocate(mat1, stat=istat)
call memocc(istat, iall, 'mat1', subname)
iall=-product(shape(mat2))*kind(mat2)
deallocate(mat2, stat=istat)
call memocc(istat, iall, 'mat2', subname)
iall=-product(shape(mat3))*kind(mat3)
deallocate(mat3, stat=istat)
call memocc(istat, iall, 'mat3', subname)


end subroutine initCompressedMatmul2




subroutine initCompressedMatmul3(norb, mad)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: norb
  type(matrixDescriptors),intent(inout):: mad

  ! Local variables
  integer:: iorb, jorb, ii, j, istat, iall, ij, iseg, i
  logical:: segment
  character(len=*),parameter:: subname='initCompressedMatmul3'
  real(8),dimension(:),allocatable:: mat1, mat2, mat3



  allocate(mat1(norb**2), stat=istat)
  call memocc(istat, mat1, 'mat1', subname)
  allocate(mat2(norb**2), stat=istat)
  call memocc(istat, mat2, 'mat2', subname)
  allocate(mat3(norb**2), stat=istat)
  call memocc(istat, mat2, 'mat2', subname)

  mat1=0.d0
  mat2=0.d0
  do iseg=1,mad%nseg
      do i=mad%keyg(1,iseg),mad%keyg(2,iseg)
          ! the localization region is "symmetric"
          mat1(i)=1.d0
          mat2(i)=1.d0
      end do
  end do

  call dgemm('n', 'n', norb, norb, norb, 1.d0, mat1, norb, mat2, norb, 0.d0, mat3, norb)

  segment=.false.
  mad%nsegmatmul=0
  do iorb=1,norb**2
      if(mat3(iorb)>0.d0) then
          ! This entry of the matrix will be different from zero.
          if(.not. segment) then
              ! This is the start of a new segment
              segment=.true.
              mad%nsegmatmul=mad%nsegmatmul+1
          end if
      else
          if(segment) then
              ! We reached the end of a segment
              segment=.false.
          end if
      end if
  end do


  allocate(mad%keygmatmul(2,mad%nsegmatmul), stat=istat)
  call memocc(istat, mad%keygmatmul, 'mad%keygmatmul', subname)
  allocate(mad%keyvmatmul(mad%nsegmatmul), stat=istat)
  call memocc(istat, mad%keyvmatmul, 'mad%keyvmatmul', subname)
  mad%keyvmatmul=0
  ! Now fill the descriptors.
  segment=.false.
  ij=0
  iseg=0
  do iorb=1,norb**2
      ij=iorb
      if(mat3(iorb)>0.d0) then
          ! This entry of the matrix will be different from zero.
          if(.not. segment) then
              ! This is the start of a new segment
              segment=.true.
              iseg=iseg+1
              mad%keygmatmul(1,iseg)=ij
          end if
          mad%keyvmatmul(iseg)=mad%keyvmatmul(iseg)+1
      else
          if(segment) then
              ! We reached the end of a segment
              segment=.false.
              mad%keygmatmul(2,iseg)=ij-1
          end if
      end if
  end do
  ! Close the last segment if required.
  if(segment) then
      mad%keygmatmul(2,iseg)=ij
  end if


iall=-product(shape(mat1))*kind(mat1)
deallocate(mat1, stat=istat)
call memocc(istat, iall, 'mat1', subname)
iall=-product(shape(mat2))*kind(mat2)
deallocate(mat2, stat=istat)
call memocc(istat, iall, 'mat2', subname)
iall=-product(shape(mat3))*kind(mat3)
deallocate(mat3, stat=istat)
call memocc(istat, iall, 'mat3', subname)


end subroutine initCompressedMatmul3





subroutine dgemm_compressed2(iproc, nproc, norb, nsegline, nseglinemax, keygline, nsegmatmul, keygmatmul, a, b, c)
!! ATTENTION: A MUST BE SYMMETRIC
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norb, nseglinemax, nsegmatmul
integer,dimension(2,nsegmatmul),intent(in):: keygmatmul
integer,dimension(norb):: nsegline
!integer,dimension(2,maxval(nsegline),norb):: keygline
integer,dimension(2,nseglinemax,norb):: keygline
real(8),dimension(norb,norb),intent(in):: a, b
real(8),dimension(norb,norb),intent(out):: c

! Local variables
integer:: iseg, i, irow, icolumn, k, iorb, jorb, korb, jseg, j, jrow, jcolumn, ii
integer:: ierr, istart, iend, iiseg, jjseg, ncount
real(8):: tt, ddot
logical:: iistop, jjstop




c=0.d0
ii=0
do iseg=1,nsegmatmul
    do i=keygmatmul(1,iseg),keygmatmul(2,iseg)
        ii=ii+1
        ! Get the row and column index
        irow=(i-1)/norb+1
        icolumn=i-(irow-1)*norb
        !c(irow,icolumn)=ddot(norb, a(1,irow), 1, b(1,icolumn), 1)
        iiseg=1
        jjseg=1
        iistop=.false.
        jjstop=.false.
        do
            istart=max(keygline(1,iiseg,irow),keygline(1,jjseg,icolumn))
            iend=min(keygline(2,iiseg,irow),keygline(2,jjseg,icolumn))
            ncount=iend-istart+1

            if(ncount>0) then
                tt=ddot(ncount, a(istart,irow), 1, b(istart,icolumn), 1)
            else
                tt=0.d0
            end if
            c(irow,icolumn) = c(irow,icolumn) + tt
            if(iiseg==nsegline(irow)) iistop=.true.
            if(jjseg==nsegline(icolumn)) jjstop=.true.
            if(iistop .and. jjstop) exit
            if((keygline(1,iiseg,irow)<=keygline(1,jjseg,icolumn) .or. jjstop) .and. .not.iistop) then
                iiseg=iiseg+1
            else
                jjseg=jjseg+1
            end if
        end do
        !if(iproc==0) write(*,'(3(a,i0),a,es15.6)') 'process ',iproc,': c(',irow,',',icolumn,')=',c(irow,icolumn)
    end do
end do
!write(*,*) 'ii, norb**2', ii, norb**2
!!do icolumn=1,norb
!!    do irow=1,norb
!!        if(iproc==0) write(200,*) icolumn, irow, c(irow,icolumn)
!!    end do
!!end do



end subroutine dgemm_compressed2



subroutine dgemm_compressed_parallel(iproc, nproc, norb, nsegline, nseglinemax, keygline, &
           nsegmatmul, keygmatmul, norb_par, isorb_par, norbp, a, b, c)
!! ATTENTION: A MUST BE SYMMETRIC
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norb, norbp, nseglinemax, nsegmatmul
integer,dimension(2,nsegmatmul),intent(in):: keygmatmul
integer,dimension(norb):: nsegline
!integer,dimension(2,maxval(nsegline),norb):: keygline
integer,dimension(2,nseglinemax,norb):: keygline
integer,dimension(0:nproc-1),intent(in):: norb_par, isorb_par
real(8),dimension(norb,norb),intent(in):: a, b
real(8),dimension(norb,norb),intent(out):: c

! Local variables
integer:: iseg, i, irow, icolumn, k, iorb, jorb, korb, jseg, j, jrow, jcolumn, ii
integer:: ierr, istart, iend, iiseg, jjseg, ncount, jproc, istat, iall, iirow, iicolumn
real(8):: tt, ddot
logical:: iistop, jjstop
integer,dimension(:),allocatable:: sendcounts, displs
real(8),dimension(:,:),allocatable:: c_loc
character(len=*),parameter:: subname='dgemm_compressed_parallel'


allocate(c_loc(norb,norbp), stat=istat)
call memocc(istat, c_loc, 'c_loc', subname)

!c=0.d0
c_loc=0.d0
ii=0
do iseg=1,nsegmatmul
    do i=keygmatmul(1,iseg),keygmatmul(2,iseg)
        ii=ii+1
        ! Get the row and column index
        !irow=(i-1)/norb+1
        !icolumn=i-(irow-1)*norb
        icolumn=(i-1)/norb+1
        irow=i-(icolumn-1)*norb
        !if(irow>isorb_par(iproc) .and. irow<=isorb_par(min(iproc+1,nproc-1))) then
        if((icolumn>isorb_par(iproc) .and. icolumn<=isorb_par(min(iproc+1,nproc-1))) .or. (iproc==nproc-1 .and. icolumn>isorb_par(iproc))) then
            !iirow=irow-isorb_par(iproc)
            iicolumn=icolumn-isorb_par(iproc)
            ! This process handles this entry of the matrix
            !c(irow,icolumn)=ddot(norb, a(1,irow), 1, b(1,icolumn), 1)
            iiseg=1
            jjseg=1
            iistop=.false.
            jjstop=.false.
            !write(*,'(a,3(i0,a))') 'process ',iproc,' calculates entry (',irow,',',iicolumn,')'
            do
                istart=max(keygline(1,iiseg,irow),keygline(1,jjseg,icolumn))
                iend=min(keygline(2,iiseg,irow),keygline(2,jjseg,icolumn))
                ncount=iend-istart+1

                if(ncount>0) then
                    tt=ddot(ncount, a(istart,irow), 1, b(istart,icolumn), 1)
                    !tt=ddot(ncount, a(istart,icolumn), 1, b(istart,irow), 1)
                else
                    tt=0.d0
                end if
                !c(irow,icolumn) = c(irow,icolumn) + tt
                !c_loc(icolumn,iirow) = c_loc(icolumn,iirow) + tt
                c_loc(irow,iicolumn) = c_loc(irow,iicolumn) + tt
                if(iiseg==nsegline(irow)) iistop=.true.
                if(jjseg==nsegline(icolumn)) jjstop=.true.
                if(iistop .and. jjstop) exit
                if((keygline(1,iiseg,irow)<=keygline(1,jjseg,icolumn) .or. jjstop) .and. .not.iistop) then
                    iiseg=iiseg+1
                else
                    jjseg=jjseg+1
                end if
            end do
            !write(*,'(5(a,i0),a,es15.6)') 'process ',iproc,': c_loc(',irow,',',iicolumn,')=c(',irow,',',icolumn,')=',c_loc(irow,iicolumn)
        end if
    end do
end do
!write(*,*) 'ii, norb**2', ii, norb**2

! Communicate the matrix.
allocate(sendcounts(0:nproc-1), stat=istat)
call memocc(istat, sendcounts, 'sendcounts', subname)
allocate(displs(0:nproc-1), stat=istat)
call memocc(istat, displs, 'displs', subname)

displs(0)=0
do jproc=0,nproc-1
    sendcounts(jproc)=norb*norb_par(jproc)
    if(jproc>0) displs(jproc)=displs(jproc-1)+sendcounts(jproc-1)
end do
call mpi_allgatherv(c_loc(1,1), sendcounts(iproc), mpi_double_precision, c(1,1), sendcounts, displs, &
     mpi_double_precision, mpi_comm_world, ierr)

iall=-product(shape(sendcounts))*kind(sendcounts)
deallocate(sendcounts, stat=istat)
call memocc(istat, iall, 'sendcounts', subname)
iall=-product(shape(displs))*kind(displs)
deallocate(displs, stat=istat)
call memocc(istat, iall, 'displs', subname)
iall=-product(shape(c_loc))*kind(c_loc)
deallocate(c_loc, stat=istat)
call memocc(istat, iall, 'c_loc', subname)

!!do icolumn=1,norb
!!    do irow=1,norb
!!        if(iproc==0) write(201,*) icolumn, irow, c(irow,icolumn)
!!    end do
!!end do


end subroutine dgemm_compressed_parallel






subroutine plotGrid(iproc, nproc, norb, nspinor, nspin, orbitalNumber, llr, glr, atoms, rxyz, hx, hy, hz)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb, nspinor, nspin, orbitalNumber
  type(locreg_descriptors),intent(in):: llr, glr
  type(atoms_data),intent(in)::atoms
  real(8),dimension(3,atoms%nat),intent(in):: rxyz
  real(8),intent(in):: hx, hy, hz
  
  ! Local variables
  integer:: iseg, jj, j0, j1, ii, i3, i2, i0, i1, i, ishift, iat, ldim, gdim, jjj, istat
  character(len=10):: num
  character(len=20):: filename
  real(8),dimension(:),allocatable:: lphi, phi


    ldim=llr%wfd%nvctr_c+7*llr%wfd%nvctr_f
    gdim=glr%wfd%nvctr_c+7*glr%wfd%nvctr_f
    allocate(lphi(ldim), stat=istat)
    allocate(phi(gdim), stat=istat)
    lphi=1.d0
    phi=0.d0
    call Lpsi_to_global2(iproc, nproc, ldim, gdim, norb, nspinor, nspin, glr, llr, lphi, phi)
  
    write(num,'(i0)') orbitalNumber
    filename='orbital_'//trim(num)
  
    open(unit=2000+iproc,file=trim(filename)//'.xyz',status='unknown')
    !write(2000+iproc,*) llr%wfd%nvctr_c+llr%wfd%nvctr_f+atoms%nat,' atomic'
    write(2000+iproc,*) glr%wfd%nvctr_c+glr%wfd%nvctr_f+llr%wfd%nvctr_c+llr%wfd%nvctr_f+atoms%nat,' atomic'
    if (atoms%geocode=='F') then
       write(2000+iproc,*)'complete simulation grid with low and high resolution points'
    else if (atoms%geocode =='S') then
       write(2000+iproc,'(a,2x,3(1x,1pe24.17))')'surface',atoms%alat1,atoms%alat2,atoms%alat3
    else if (atoms%geocode =='P') then
       write(2000+iproc,'(a,2x,3(1x,1pe24.17))')'periodic',atoms%alat1,atoms%alat2,atoms%alat3
    end if

   do iat=1,atoms%nat
      write(2000+iproc,'(a6,2x,3(1x,e12.5),3x)') trim(atoms%atomnames(atoms%iatype(iat))),rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
   end do

  
    jjj=0
    do iseg=1,glr%wfd%nseg_c
       jj=glr%wfd%keyv(iseg)
       j0=glr%wfd%keyg(1,iseg)
       j1=glr%wfd%keyg(2,iseg)
       ii=j0-1
       i3=ii/((glr%d%n1+1)*(glr%d%n2+1))
       ii=ii-i3*(glr%d%n1+1)*(glr%d%n2+1)
       i2=ii/(glr%d%n1+1)
       i0=ii-i2*(glr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
           jjj=jjj+1
           if(phi(jjj)==1.d0) write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  lg ',real(i,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  g ',real(i,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
       enddo
    enddo
  
    ! fine part
    !$omp do
    ishift=glr%wfd%nseg_c
    do iseg=1,glr%wfd%nseg_f
       jj=glr%wfd%keyv(ishift+iseg)
       j0=glr%wfd%keyg(1,ishift+iseg)
       j1=glr%wfd%keyg(2,ishift+iseg)
       ii=j0-1
       i3=ii/((glr%d%n1+1)*(glr%d%n2+1))
       ii=ii-i3*(glr%d%n1+1)*(glr%d%n2+1)
       i2=ii/(glr%d%n1+1)
       i0=ii-i2*(glr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
          jjj=jjj+1
          if(phi(jjj)==1.d0) write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  lG ',real(i,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
          write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  G ',real(i,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
          jjj=jjj+6
       enddo
    enddo
  
    close(unit=2000+iproc)

end subroutine plotGrid




subroutine repartitionOrbitals(iproc, nproc, norb, norb_par, norbp, isorb_par, isorb, onWhichMPI)
  use module_base
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb
  integer,dimension(0:nproc-1),intent(out):: norb_par, isorb_par
  integer,dimension(norb),intent(out):: onWhichMPI
  integer,intent(out):: norbp, isorb

  ! Local variables
  integer:: ii, kk, iiorb, mpiflag, iorb, ierr, jproc
  real(8):: tt

  ! Determine norb_par
  norb_par=0
  tt=dble(norb)/dble(nproc)
  ii=floor(tt)
  ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
  norb_par(0:nproc-1)=ii
  kk=norb-nproc*ii
  norb_par(0:kk-1)=ii+1

  ! Determine norbp
  norbp=norb_par(iproc)

  ! Determine isorb
  isorb=0
  do jproc=0,iproc-1
      isorb=isorb+norb_par(jproc)
  end do

  ! Determine onWhichMPI and isorb_par
  iiorb=0
  isorb_par=0
  do jproc=0,nproc-1
      do iorb=1,norb_par(jproc)
          iiorb=iiorb+1
          onWhichMPI(iiorb)=jproc
      end do
      if(iproc==jproc) then
          isorb_par(jproc)=isorb
      end if
  end do
  call MPI_Initialized(mpiflag,ierr)
  if(mpiflag /= 0) call mpiallred(isorb_par(0), nproc, mpi_sum, mpi_comm_world, ierr)


end subroutine repartitionOrbitals




subroutine repartitionOrbitals2(iproc, nproc, norb, norb_par, norbp, isorb)
  use module_base
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb
  integer,dimension(0:nproc-1),intent(out):: norb_par
  integer,intent(out):: norbp, isorb

  ! Local variables
  integer:: ii, kk, iiorb, mpiflag, iorb, ierr, jproc
  real(8):: tt

  ! Determine norb_par
  norb_par=0
  tt=dble(norb)/dble(nproc)
  ii=floor(tt)
  ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
  norb_par(0:nproc-1)=ii
  kk=norb-nproc*ii
  norb_par(0:kk-1)=ii+1

  ! Determine norbp
  norbp=norb_par(iproc)

  ! Determine isorb
  isorb=0
  do jproc=0,iproc-1
      isorb=isorb+norb_par(jproc)
  end do


end subroutine repartitionOrbitals2
