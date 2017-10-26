!> @file
!!  Module to handle the fragments of a system
!! @author
!!    Copyright (C) 2013-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Module which defines some important structures and methods to manipulate embedding systems
module module_fragments
  use module_base, only: gp,wp,bigdft_mpi,dp
  use module_types
  use dynamic_memory
  use module_atoms
  use module_input_keys, only: input_variables
  use fragment_base, only: fragmentInputParameters
  implicit none

  private

  !> information about the basis set metadata shared between cubic and ASF approaches
  type, public :: minimal_orbitals_data
     integer :: norb          !< Total number of orbitals per k point
     integer :: norbp         !< Total number of orbitals for the given processors
     integer :: isorb         !< Total number of orbitals for the given processors
     integer, dimension(:), pointer :: inwhichlocreg,onwhichatom !< associate the basis centers
     integer, dimension(:), pointer :: isorb_par,ispot
     integer, dimension(:,:), pointer :: norb_par
  end type minimal_orbitals_data

  type, public :: phi_array
     real(wp), dimension(:,:,:,:,:,:), pointer :: psig
  end type phi_array

  type, public :: fragment_basis
     integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
     type(local_zone_descriptors) :: Lzd
     type(minimal_orbitals_data) :: forbs
     real(wp), dimension(:), pointer :: psi_full
     type(phi_array), dimension(:), pointer :: phi
     real(wp), dimension(:), pointer :: density
  end type fragment_basis

  type, public :: fragment_orbitals
     integer :: norb, nphi

  end type fragment_orbitals

  !> Defines the minimal information to identify a system building block
  type, public :: system_fragment
     !integer :: nat_env !< environment atoms which complete fragment specifications
     !real(gp), dimension(:,:), pointer :: rxyz_env !< position of atoms in environment (AU), external reference frame
     type(atomic_structure) :: astruct_frg !< Number of atoms, positions, atom type etc for fragment
     type(atomic_structure) :: astruct_env !< Number of atoms, positions, atom type etc for fragment environment (includes fragment)
     integer :: nbasis_env !< Number of tmbs for environment (includes fragment)
     !integer, dimension(:,:), pointer :: env_mapping !< contains mapping information between fragment environment and full system
     type(fragment_basis) :: fbasis !< fragment basis, associated only if coherent with positions, pointer - do we really want this to be a pointer?
     ! add coeffs and or kernel
     integer :: nelec
     real(gp), dimension(:,:), pointer :: coeff
     real(gp), dimension(:,:,:), pointer :: kernel
     real(gp), dimension(:,:,:), pointer :: kernel_env
     real(gp), dimension(:), pointer :: eval
  end type system_fragment

  !public operator(*)

  public :: fragment_null, fragment_free, init_fragments, minimal_orbitals_data_null, fragmentInputParameters
  public :: calculate_fragment_density, init_minimal_orbitals_data


contains


  !> Need to check somewhere if fragment linking is consistent - initialize linear from file?!
  !! initializes reference fragments (already nullified), if it isn't a fragment calculation sets to appropriate dummy values
  !! ignoring environment for now
  subroutine init_fragments(in,orbs,astruct,ref_frags)
    implicit none
    type(input_variables), intent(in) :: in
    type(orbitals_data), intent(in) :: orbs ! orbitals of full system, needed to set 'dummy' values
    type(atomic_structure), intent(in) :: astruct ! atomic structure of full system, needed to set 'dummy' values
    type(system_fragment), dimension(in%frag%nfrag_ref), intent(inout) :: ref_frags

    ! local variables
    integer :: ifrag

    call f_routine(id='init_fragments')

    if (in%lin%fragment_calculation) then
        ! read fragment posinps and initialize fragment, except for psi and lzds
        do ifrag=1,in%frag%nfrag_ref
           call init_fragment_from_file(ref_frags(ifrag),trim(in%dir_output)//trim(in%frag%label(ifrag)),&
                in,astruct)
        end do

        ! check that fragments are sensible, i.e. correct number of atoms, atom types etc.
        call check_fragments(in,ref_frags,astruct)

     ! set appropriate default values as this is not a fragment calculation
     else
        ! nullify fragment
        ref_frags(1)=fragment_null()

        ! want all components of ref_frags(1)%fbasis%forbs to point towards appropriate component of orbs of full system
        call orbs_to_min_orbs_point(orbs,ref_frags(1)%fbasis%forbs)

        ! environment and coeffs
        call fragment_allocate(ref_frags(1))

        ! rest of fbasis isn't needed (I think!) so can remain nullified
  
        ! astruct - fill in other bits later
        ref_frags(1)%astruct_frg%nat=astruct%nat
        ref_frags(1)%astruct_env%nat=0

     end if

    call f_release_routine()

  end subroutine init_fragments


  !> Initializes all of fragment except lzd using the fragment posinp and tmb files
  subroutine init_fragment_from_file(frag,frag_name,input,astruct) ! switch this to pure if possible
    !use module_interfaces
    implicit none
    type(system_fragment), intent(inout) :: frag
    character(len=*), intent(in) :: frag_name
    type(input_variables), intent(in) :: input
    type(atomic_structure), intent(in) :: astruct ! atomic structure of full system
    
    logical :: env_exists
    integer :: iat, ityp, jtyp

    call f_routine(id='init_fragment_from_file')

    ! nullify fragment
    frag=fragment_null()

    ! read fragment positions
    call set_astruct_from_file(frag_name(1:len(frag_name)),bigdft_mpi%iproc,frag%astruct_frg)

    ! first check if an environment file exists
    inquire(FILE = frag_name(1:len(frag_name))//'_env.xyz', EXIST = env_exists)

    if (env_exists) then

       call set_astruct_from_file(frag_name(1:len(frag_name))//'_env',bigdft_mpi%iproc,frag%astruct_env)
       ! check that this contains at least 1 environment atom
       if (frag%astruct_env%nat < frag%astruct_frg%nat+1) then
          stop 'Fragment environment file missing some atoms'
       end if

       ! calculate nbasis_env here:
       frag%nbasis_env=0
       !example with iterators
!!$       atit=atoms_iter(frag%astruct_env)
!!$       do while (atoms_iter_next(atit))
!!$          jtyp= .index. (astruct%atomnames == atit%name) !to be defined
!!$          if (jtyp > 0) then
!!$             frag%nbasis_env=&
!!$                  frag%nbasis_env+input%lin%norbsPerType(jtyp)
!!$          else
!!$             call f_err_throw('Error in fragment_init_orbitals, atom type '+&
!!$                  atit%name+' does not exist in full structure',&
!!$                  err_name='BIGDFT_INPUT_FILE_ERROR')
!!$
!!$          end if
!!$       end do

       do iat=1,frag%astruct_env%nat
          ityp=frag%astruct_env%iatype(iat)
          do jtyp=1,astruct%ntypes
             if (astruct%atomnames(jtyp)==frag%astruct_env%atomnames(ityp)) exit
          end do
          ! if it's a ghost atom 'X' then we don't expect it to be in full structure
          if (jtyp==astruct%ntypes+1 .and. (trim(frag%astruct_env%atomnames(ityp))/='X')) then
             print*, 'Error in fragment_init_orbitals, atom type ',frag%astruct_env%atomnames(ityp),&
                  ' does not exist in full structure'
             stop
          end if
          if (trim(frag%astruct_env%atomnames(ityp))/='X') frag%nbasis_env=frag%nbasis_env+input%lin%norbsPerType(jtyp)
       end do

    else
       frag%astruct_env%nat=0
    end if

    ! iproc, nproc, nspinor not needed yet, add in later
    call init_minimal_orbitals_data(bigdft_mpi%iproc, bigdft_mpi%nproc, 1, input, frag%astruct_frg, &
         frag%fbasis%forbs,astruct)
    !call init_minimal_orbitals_data(iproc, nproc, nspinor, input, frag%astruct_frg, frag%fbasis%forbs)

    ! environment and coeffs
    call fragment_allocate(frag)

    ! allocate/initialize fragment basis...
    !currently orbitals are initialized via initAndUtils/init_orbitals_data_for_linear
    !which calls wavefunctions/orbitals_descriptors to assign parallel bits and superfluous stuff
    !and locreg_orbitals/assignToLocreg2 to give inwhichlocreg and onwhichatom - but we should really take onwhichatom from file?!

    !type, public :: fragment_basis
    !   integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
    !   integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
    !   type(local_zone_descriptors) :: Lzd
    !   type(minimal_orbitals_data) :: forbs

    call f_release_routine()

  end subroutine init_fragment_from_file


  !> Sanity check on fragment definitions
  subroutine check_fragments(input,ref_frags,astruct)
    implicit none
    type(input_variables), intent(in) :: input
    type(atomic_structure), intent(in) :: astruct ! atomic structure of full system, needed to check fragments are sensible
    type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags

    integer :: ifrag, ifrag_ref, tot_frag_ats, isfat, iat
    logical :: fragments_ok

    fragments_ok=.true.

    tot_frag_ats=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref = input%frag%frag_index(ifrag)
       tot_frag_ats = tot_frag_ats + ref_frags(ifrag_ref)%astruct_frg%nat
    end do

    if (tot_frag_ats /= astruct%nat) then
       fragments_ok=.false.
       write(*,*) 'Sum of atoms in fragments not equal to total atoms in system',tot_frag_ats,astruct%nat
    end if

    isfat=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref=input%frag%frag_index(ifrag)

       do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
          if (trim(ref_frags(ifrag_ref)%astruct_frg%atomnames(ref_frags(ifrag_ref)%astruct_frg%iatype(iat))) &
               /= trim(astruct%atomnames(astruct%iatype(iat+isfat)))) then
             fragments_ok=.false.
             write(*,*) 'Atom type for fragment ',ifrag,', reference fragment ',ifrag_ref,' atom number ',iat,&
                  ' does not match ',&
                  trim(ref_frags(ifrag_ref)%astruct_frg%atomnames(ref_frags(ifrag_ref)%astruct_frg%iatype(iat))),&
                  trim(astruct%atomnames(astruct%iatype(iat+isfat)))
          end if
       end do

       isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat     
    end do

    if (.not. fragments_ok) stop 'Problem in fragment definitions'

  end subroutine check_fragments


  !type, public :: minimal_orbitals_data
  !   integer :: norb          !< Total number of orbitals per k point
  !   integer :: norbp         !< Total number of orbitals for the given processors
  !   integer :: isorb         !< Total number of orbitals for the given processors
  !   integer, dimension(:), pointer :: inwhichlocreg,onwhichatom !< associate the basis centers
  !   integer, dimension(:), pointer :: isorb_par,ispot
  !   integer, dimension(:,:), pointer :: norb_par
  !> Just initializing norb for now, come back and do the rest later
  subroutine init_minimal_orbitals_data(iproc, nproc, nspinor, input, astruct, forbs, astruct_full)
    use module_base
    implicit none
  
    ! Calling arguments
    integer,intent(in) :: iproc, nproc, nspinor
    type(input_variables),intent(in) :: input
    type(atomic_structure),intent(in) :: astruct
    type(minimal_orbitals_data),intent(out) :: forbs
    type(atomic_structure), intent(in) :: astruct_full ! atomic structure of full system
  
    ! Local variables
    integer :: norb, ityp, iat, jtyp
    !integer :: norbu, norbd, nlr, ilr, iall, iorb, istat
    !integer,dimension(:),allocatable :: norbsPerLocreg, norbsPerAtom
    !real(kind=8),dimension(:,:),allocatable :: locregCenter
    !character(len=*),parameter :: subname='init_minimal_orbitals_data'

    call timing(iproc,'init_orbs_lin ','ON')
 
    ! Count the number of basis functions.
    !allocate(norbsPerAtom(astruct%nat), stat=istat)
    !call memocc(istat, norbsPerAtom, 'norbsPerAtom', subname)
    norb=0
    !nlr=0
    do iat=1,astruct%nat
       ityp=astruct%iatype(iat)
       !norbsPerAtom(iat)=input%lin%norbsPerType(ityp)
       do jtyp=1,astruct_full%ntypes
          if (astruct_full%atomnames(jtyp)==astruct%atomnames(ityp)) exit
       end do
       if (jtyp==astruct_full%ntypes+1 .and. trim(astruct%atomnames(ityp))/='X') then
          print*, 'Error in fragment_init_orbitals, atom type ',astruct%atomnames(ityp),' does not exist in full structure'
          stop
       end if
       norb=norb+input%lin%norbsPerType(jtyp)
       !nlr=nlr+input%lin%norbsPerType(ityp)
    end do

    ! Distribute the basis functions among the processors.
    !norbu=norb
    !norbd=0
 
    forbs%norb=norb

    !call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, nspinor,&
    !     input%nkpt, input%kpt, input%wkpt, lorbs,.true.) !simple repartition
 

    !allocate(locregCenter(3,nlr), stat=istat)
    !call memocc(istat, locregCenter, 'locregCenter', subname)
  
    !ilr=0
    !do iat=1,astruct%nat
    !    ityp=astruct%iatype(iat)
    !    do iorb=1,input%lin%norbsPerType(ityp)
    !        ilr=ilr+1
    !        locregCenter(:,ilr)=astruct%rxyz(:,iat)
    !        ! DEBUGLR write(10,*) iorb,locregCenter(:,ilr)
    !    end do
    !end do
 
    !allocate(norbsPerLocreg(nlr), stat=istat)
    !call memocc(istat, norbsPerLocreg, 'norbsPerLocreg', subname)
    !norbsPerLocreg=1 !should be norbsPerLocreg
    
    !iall=-product(shape(forbs%inWhichLocreg))*kind(forbs%inWhichLocreg)
    !deallocate(forbs%inWhichLocreg, stat=istat)
    !call memocc(istat, iall, 'forbs%inWhichLocreg', subname)
    !call assignToLocreg2(iproc, nproc, forbs%norb, forbs%norb_par, astruct%nat, nlr, &
    !     input%nspin, norbsPerLocreg, locregCenter, forbs%inwhichlocreg)

    !iall=-product(shape(forbs%onwhichatom))*kind(forbs%onwhichatom)
    !deallocate(forbs%onwhichatom, stat=istat)
    !call memocc(istat, iall, 'forbs%onwhichatom', subname)
    !call assignToLocreg2(iproc, nproc, forbs%norb, forbs%norb_par, astruct%nat, astruct%nat, &
    !     input%nspin, norbsPerAtom, astruct%rxyz, forbs%onwhichatom)
  

    !iall=-product(shape(norbsPerLocreg))*kind(norbsPerLocreg)
    !deallocate(norbsPerLocreg, stat=istat)
    !call memocc(istat, iall, 'norbsPerLocreg', subname)
  
    !iall=-product(shape(locregCenter))*kind(locregCenter)
    !deallocate(locregCenter, stat=istat)
    !call memocc(istat, iall, 'locregCenter', subname)

    !iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
    !deallocate(norbsPerAtom, stat=istat)
    !call memocc(istat, iall, 'norbsPerAtom', subname)


    call timing(iproc,'init_orbs_lin ','OF')

  end subroutine init_minimal_orbitals_data


  !> Point minimal orbs structure to a given full orbs structure
  subroutine orbs_to_min_orbs_point(orbs,forbs)
    implicit none
    ! Calling arguments
    type(orbitals_data),intent(in):: orbs
    type(minimal_orbitals_data),intent(inout):: forbs

    forbs%norb = orbs%norb
    forbs%norbp = orbs%norbp
    forbs%isorb = orbs%isorb

    !forbs%norb_par => orbs%norb_par
    !forbs%inwhichlocreg => orbs%inwhichlocreg
    !forbs%onwhichatom => orbs%onwhichatom
    !forbs%isorb_par => orbs%isorb_par
    !forbs%ispot => orbs%ispot
    forbs%norb_par      = f_malloc_ptr(src_ptr=orbs%norb_par,id='forbs%norb_par')
    forbs%inwhichlocreg = f_malloc_ptr(src_ptr=orbs%inwhichlocreg,id='forbs%inwhichlocreg')
    forbs%onwhichatom   = f_malloc_ptr(src_ptr=orbs%onwhichatom,id='forbs%onwhichatom')
    forbs%isorb_par     = f_malloc_ptr(src_ptr=orbs%isorb_par,id='forbs%isorb_par')
    forbs%ispot         = f_malloc_ptr(src_ptr=orbs%ispot,id='forbs%ispot')

!!$    if(associated(orbs%norb_par)) then
!!$        forbs%norb_par = f_malloc_ptr(src=orbs%norb_par,lbounds=lbound(orbs%norb_par),id='forbs%norb_par')
!!$    end if
!!$    if(associated(orbs%inwhichlocreg)) then
!!$        forbs%inwhichlocreg = f_malloc_ptr(src=orbs%inwhichlocreg,lbounds=lbound(orbs%inwhichlocreg),id='forbs%inwhichlocreg')
!!$    end if
!!$    if(associated(orbs%onwhichatom)) then
!!$        forbs%onwhichatom = f_malloc_ptr(src=orbs%onwhichatom,lbounds=lbound(orbs%onwhichatom),id='forbs%onwhichatom')
!!$    end if
!!$    if(associated(orbs%isorb_par)) then
!!$        forbs%isorb_par = f_malloc_ptr(src=orbs%isorb_par,lbounds=lbound(orbs%isorb_par),id='forbs%isorb_par')
!!$    end if
!!$    if(associated(orbs%ispot)) then
!!$        forbs%ispot = f_malloc_ptr(src=orbs%ispot,lbounds=lbound(orbs%ispot),id='forbs%ispot')
!!$    end if

  end subroutine orbs_to_min_orbs_point


  subroutine calculate_fragment_density(frag,ndimrho,tmb,iorb_start,charge,atoms,rxyz,denspot)
    use locreg_operations, only: Lpsi_to_global2,workarr_sumrho,initialize_work_arrays_sumrho,deallocate_work_arrays_sumrho
    implicit none
    type(system_fragment), intent(inout) :: frag
    integer, intent(in) :: ndimrho ! add to fragment structure?
    integer, intent(in) :: iorb_start ! the first tmb for this fragment
    integer, intent(in) :: charge ! charge on this fragment for calculating correct kernel?!
    type(dft_wavefunction), intent(in) :: tmb
    type(atoms_data), intent(in) :: atoms ! just for plotting
    real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz ! just for plotting
    type(DFT_local_fields), intent(in) :: denspot ! just for plotting

    ! local variables
    integer :: iiorb, ilr, ind, indg, indr, ipt, jjorb, iorb, itmb
    real(kind=gp),dimension(:,:), allocatable :: kernel, fcoeff
    real(kind=wp),dimension(:), allocatable :: gpsi
    real(kind=wp),dimension(:), allocatable :: psir
    real(kind=gp) :: total_charge, tt, tt1, factor, hxh, hyh, hzh
    type(workarr_sumrho) :: w

    ! calculate fragment density and store in frag%fbasis%density
    ! allocate to cover whole simulation cell or just the fragment region?
    call f_routine(id='calculate_fragment_density')
    frag%fbasis%density=f_malloc_ptr(ndimrho,id='frag%fbasis%density')
    kernel=f_malloc((/frag%fbasis%forbs%norb,frag%fbasis%forbs%norb/),id='kernel')
    fcoeff=f_malloc((/frag%fbasis%forbs%norb,frag%fbasis%forbs%norb/),id='fcoeff')

    ! calculate density kernel for this fragment (ideally would use routine but don't have access to correct orbs and doing monoproc)
    ! add occup
    do iorb=1,frag%fbasis%forbs%norb
       do itmb=1,frag%fbasis%forbs%norb
          if (iorb<=floor((frag%nelec-charge)/2.0_gp)) then
             fcoeff(itmb,iorb) = 2.0_gp*frag%coeff(itmb,iorb)
          else if (iorb<=ceiling((frag%nelec-charge)/2.0_gp)) then
             fcoeff(itmb,iorb) = 1.0_gp*frag%coeff(itmb,iorb)
          else
             fcoeff(itmb,iorb) = 0.0_gp
          end if
       end do
    end do

    call dgemm('n', 't', frag%fbasis%forbs%norb, frag%fbasis%forbs%norb, &
         nint((frag%nelec-charge)/2.0_gp), 1.d0, &
         frag%coeff(1,1), frag%fbasis%forbs%norb, &
         fcoeff(1,1), frag%fbasis%forbs%norb, 0.d0, &
         kernel(1,1), frag%fbasis%forbs%norb)

    call f_free(fcoeff)

    ! expand all tmbs to global box to simplify density calculation and convert to real space
    ! as we don't have psi_full allocated, just using tmb%psi for now, but need to know where to start
    ind=1
    indg=1
    indr=1

    gpsi=f_malloc0(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f,id='gpsi')
    !call f_zero(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f,gpsi)
    call initialize_work_arrays_sumrho(tmb%lzd%glr,.true.,w)
    psir=f_malloc(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*tmb%lzd%glr%d%n3i*frag%fbasis%forbs%norb,id='psir')

    do iiorb=1,tmb%orbs%norb
       ilr = tmb%orbs%inwhichlocreg(iiorb)
       if (iiorb < iorb_start) then
          ind = ind + tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f
          cycle
       else if (iiorb >= iorb_start + frag%fbasis%forbs%norb) then
          exit
       end if

       call Lpsi_to_global2(bigdft_mpi%iproc, tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f, &
            tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f, &
            1, 1, tmb%Lzd%glr, tmb%Lzd%Llr(ilr), tmb%psi(ind:), gpsi(indg:))

       call daub_to_isf(tmb%lzd%glr, w, gpsi(indg), psir(indr))

       !indg = indg + tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f
       ind = ind + tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f
       indr = indr + tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*tmb%lzd%glr%d%n3i
    end do

    call deallocate_work_arrays_sumrho(w)
    call f_free(gpsi)


    ! calculate density using kernel and psi
    hxh=.5d0*tmb%lzd%hgrids(1)
    hyh=.5d0*tmb%lzd%hgrids(2)
    hzh=.5d0*tmb%lzd%hgrids(3)
    factor=1.d0/(hxh*hyh*hzh)
    total_charge=0.d0
    !call to_zero(ndimrho,frag%fbasis%density)
    do ipt=1,ndimrho
       tt=1.e-20_dp
       do iiorb=1,frag%fbasis%forbs%norb
          tt1=psir(ipt+(iiorb-1)*ndimrho)
          tt=tt+kernel(iiorb,iiorb)*tt1**2
          do jjorb=iiorb+1,frag%fbasis%forbs%norb
             tt=tt+2.0_dp*kernel(jjorb,iiorb)*tt1*psir(ipt+(jjorb-1)*ndimrho)
          end do
       end do  
        tt=factor*tt
        total_charge=total_charge+tt
        ! apply a cut-off
        if (tt >= 1.e-6) then
           frag%fbasis%density(ipt)=tt
        else
           frag%fbasis%density(ipt)=0.0_gp
        end if
    end do

    print*,'Fragment density total charge',total_charge*hxh*hyh*hzh

    !call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'fragment_density.cube', &
    !     atoms,rxyz,denspot%dpbox,1,frag%fbasis%density)

    ! Define some constant factors.
    !hxh=.5d0*hx
    !hyh=.5d0*hy
    !hzh=.5d0*hz
    !factor=1.d0/(hxh*hyh*hzh)

    !do ipt=1,collcom_sr%nptsp_c
    !    ii=collcom_sr%norb_per_gridpoint_c(ipt)
    !    i0=collcom_sr%isptsp_c(ipt)
    !    tt=1.e-20_dp
    !    do i=1,ii
    !        iiorb=collcom_sr%indexrecvorbital_c(i0+i)
    !        tt1=collcom_sr%psit_c(i0+i)
    !        ind=denskern%matrixindex_in_compressed(iiorb,iiorb)
    !        tt=tt+denskern%matrix_compr(ind)*tt1*tt1
    !        do j=i+1,ii
    !            jjorb=collcom_sr%indexrecvorbital_c(i0+j)
    !            ind=denskern%matrixindex_in_compressed(jjorb,iiorb)
    !            if (ind==0) cycle
    !            tt=tt+2.0_dp*denskern%matrix_compr(ind)*tt1*collcom_sr%psit_c(i0+j)
    !        end do
    !    end do
    !    tt=factor*tt
    !    total_charge=total_charge+tt
    !    rho_local(ipt)=tt
    !end do

    call f_free(psir)
    call f_free(kernel)
    call f_release_routine()


  end subroutine calculate_fragment_density


  pure function minimal_orbitals_data_null() result(forbs)
    implicit none
    type(minimal_orbitals_data) :: forbs
    call nullify_minimal_orbitals_data(forbs)
  end function minimal_orbitals_data_null

  pure subroutine nullify_minimal_orbitals_data(forbs)
    implicit none
    type(minimal_orbitals_data), intent(out) :: forbs

    forbs%norb=0
    forbs%norbp=0
    forbs%isorb=0
    nullify(forbs%inwhichlocreg)
    nullify(forbs%onwhichatom)
    nullify(forbs%isorb_par)
    nullify(forbs%ispot)
    nullify(forbs%norb_par)
  end subroutine nullify_minimal_orbitals_data


  pure function fragment_null() result(frag)
    implicit none
    type(system_fragment) :: frag

    !frag%nat_env=0
    !nullify(frag%rxyz_env)
    frag%nelec=0
    frag%nbasis_env=0
    nullify(frag%coeff)
    nullify(frag%kernel)
    nullify(frag%kernel_env)
    nullify(frag%eval)
    !nullify(frag%env_mapping)
    call nullify_atomic_structure(frag%astruct_frg)
    call nullify_atomic_structure(frag%astruct_env)
    ! nullify fragment basis
    call nullify_fragment_basis(frag%fbasis)

  end function fragment_null


  pure function fragment_basis_null() result(basis)
    implicit none
    type(fragment_basis) :: basis
    call nullify_fragment_basis(basis)
  end function fragment_basis_null

  pure subroutine nullify_fragment_basis(basis)
    implicit none
    type(fragment_basis), intent(out) :: basis

    basis%npsidim_orbs=0
    basis%npsidim_comp=0
    call nullify_local_zone_descriptors(basis%lzd)
    call nullify_minimal_orbitals_data(basis%forbs)
    !basis%forbs=minimal_orbitals_data_null()
    nullify(basis%psi_full)
    nullify(basis%phi)
    nullify(basis%density)
  end subroutine nullify_fragment_basis
  

  !> This routine is dangerous as it frees orbs, when forbs points to it
  !! either garbage collectors or other techniques should be considered
  subroutine minimal_orbitals_data_free(forbs)
    implicit none
    type(minimal_orbitals_data), intent(inout) :: forbs

    if (associated(forbs%inwhichlocreg)) call f_free_ptr(forbs%inwhichlocreg)
    if (associated(forbs%onwhichatom)) call f_free_ptr(forbs%onwhichatom)
    if (associated(forbs%isorb_par)) call f_free_ptr(forbs%isorb_par)
    if (associated(forbs%ispot)) call f_free_ptr(forbs%ispot)
    if (associated(forbs%norb_par)) call f_free_ptr(forbs%norb_par)
    forbs=minimal_orbitals_data_null()
  end subroutine minimal_orbitals_data_free


  subroutine fragment_basis_free(basis)
    implicit none
    type(fragment_basis), intent(inout) :: basis
    character(len=*), parameter :: subname='fragment_basis_free'
    call deallocate_local_zone_descriptors(basis%lzd)
    call minimal_orbitals_data_free(basis%forbs)
    if (associated(basis%psi_full)) call f_free_ptr(basis%psi_full)
    if (associated(basis%density)) call f_free_ptr(basis%density)
    basis=fragment_basis_null()
  end subroutine fragment_basis_free


  subroutine fragment_free(frag)
    implicit none
    type(system_fragment), intent(inout) :: frag

    call deallocate_atomic_structure(frag%astruct_frg)
    frag%astruct_frg=atomic_structure_null()
    call deallocate_atomic_structure(frag%astruct_env)
    frag%astruct_env=atomic_structure_null()
    call minimal_orbitals_data_free(frag%fbasis%forbs)
    frag%fbasis%forbs = minimal_orbitals_data_null()
    !call f_free_ptr(frag%rxyz_env)
    call f_free_ptr(frag%coeff)
    call f_free_ptr(frag%kernel)
    call f_free_ptr(frag%kernel_env)
    call f_free_ptr(frag%eval)
    !call f_free_ptr(frag%env_mapping)
    call fragment_basis_free(frag%fbasis)
    frag=fragment_null()

  end subroutine fragment_free


  subroutine fragment_allocate(frag)
    implicit none
    type(system_fragment), intent(inout) :: frag

    call f_routine(id='fragment_allocate')

    !frag%rxyz_env=f_malloc_ptr((/3,min(1,frag%nat_env)/),id='frag%rxyz_env')
    frag%coeff=f_malloc_ptr((/frag%fbasis%forbs%norb,frag%fbasis%forbs%norb/),id='frag%coeff')
    frag%kernel=f_malloc_ptr((/frag%fbasis%forbs%norb,frag%fbasis%forbs%norb,1/),id='frag%kernel') !NEED SPIN HERE
    frag%kernel_env=f_malloc_ptr((/frag%nbasis_env,frag%nbasis_env,1/),id='frag%kernel_env') !NEED SPIN HERE
    frag%eval=f_malloc_ptr(frag%fbasis%forbs%norb,id='frag%eval')
    !frag%env_mapping=f_malloc0_ptr((/frag%nbasis_env,3/),id='frag%env_mapping')

    call f_release_routine()

  end subroutine fragment_allocate

  !function transform_fragment(trans,frag) result(frag_new)
  !  implicit none
  !  type(fragment_transformation), intent(in) :: trans
  !  type(system_fragment), intent(in) :: frag
  !  type(system_fragment) :: frag_new
  !
  !  ! local variables
  !  integer :: iat
  !
  !  frag_new=fragment_null()
  !  frag_new%astruct_frg%nat=frag%astruct_frg%nat
  !  frag_new%nat_env=frag%nat_env
  !  frag_new%astruct_frg%ntypes=frag%astruct_frg%ntypes
  !
  !  ! allocate arrays here, leave fragment_basis nullified
  !  call fragment_allocate(frag_new)
  !
  !  ! do fragment first, then environment
  !  do iat=1,frag%astruct_frg%nat
  !     frag_new%astruct_frg%rxyz(:,iat)=rotate_vector(trans%rot_axis,&
  !          trans%theta,frag%astruct_frg%rxyz(:,iat)-trans%rot_center(:))
  !     frag_new%astruct_frg%rxyz(:,iat)=frag_new%astruct_frg%rxyz(:,iat)+trans%rot_center(:)+trans%dr(:)
  !  end do
  !
  !  do iat=1,frag%nat_env
  !     frag_new%rxyz_env(:,iat)=rotate_vector(trans%rot_axis,trans%theta,frag%rxyz_env(:,iat)-trans%rot_center(:))
  !     frag_new%rxyz_env(:,iat)=frag_new%rxyz_env(:,iat)+trans%rot_center(:)+trans%dr(:)
  !  end do
  !
  !  ! to complete should copy across iatype and atomnames but only using as a check to see how effective trans is
  !
  !end function transform_fragment

end module module_fragments
