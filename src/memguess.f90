!> @file
!!   Program to guess the used memory by BigDFT
!! @author
!!   Copyright (C) 2007-2011 BigDFT group (LG)
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Test the input files and estimates the memory occupation versus the number
!! of processors
program memguess

  use module_base
  use module_types
  use module_interfaces
  use module_xc
  use m_ab6_symmetry

  implicit none
  character(len=*), parameter :: subname='memguess'
  character(len=20) :: tatonam, radical
  character(len=40) :: comment
  character(len=128) :: fileFrom, fileTo,filename_wfn
  logical :: optimise,GPUtest,atwf,convert=.false.,exportwf=.false.
  integer :: nelec,ntimes,nproc,i_stat,i_all,output_grid, i_arg,istat
  integer :: norbe,norbsc,nspin,iorb,norbu,norbd,nspinor,norb
  integer :: norbgpu,nspin_ig,ng
  real(gp) :: peakmem,hx,hy,hz
  type(input_variables) :: in
  type(atoms_data) :: atoms
  type(orbitals_data) :: orbs,orbstst
  type(communications_arrays) :: comms
  type(locreg_descriptors) :: Glr
  type(nonlocal_psp_descriptors) :: nlpspd
  type(gaussian_basis) :: G !basis for davidson IG
  real(gp), dimension(3) :: shift
  logical, dimension(:,:,:), allocatable :: logrid
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(gp), dimension(:,:), pointer :: rxyz
  real(wp), dimension(:), allocatable :: rhoexpo,psi
  real(wp), dimension(:,:), pointer :: rhocoeff
  real(kind=8), dimension(:,:), allocatable :: radii_cf
  logical, dimension(:,:,:), allocatable :: scorb
  real(kind=8), dimension(:), allocatable :: locrad
  real(gp), dimension(:), pointer :: gbd_occ
  !! By Ali
  integer :: ierror

  ! Get arguments

  !call getarg(1,tatonam)
  call get_command_argument(1, value = tatonam, status = istat)

  write(radical, "(A)") ""
  optimise=.false.
  GPUtest=.false.
  atwf=.false.
  if(trim(tatonam)=='' .or. istat>0) then
     write(*,'(1x,a)')&
          'Usage: ./memguess <nproc> [option]'
     write(*,'(1x,a)')&
          'Indicate the number of processes after the executable'
     write(*,'(1x,a)')&
          '[option] can be the following: '
     write(*,'(1x,a)')&
          '"y": grid to be plotted with V_Sim'
     write(*,'(1x,a)')&
          '"o" rotate the molecule such that the volume of the simulation box is optimised'
     write(*,'(1x,a)')&
          '"GPUtest <nrep>" case of a CUDAGPU calculation, to test the speed of 3d operators'
     write(*,'(1x,a)')&
          '         <nrep> is the number of repeats'
     write(*,'(1x,a)')&
          '"upgrade" upgrades input files older than 1.2 into actual format'
     write(*,'(1x,a)')&
          '"convert" <from.[cube,etsf]> <to.[cube,etsf]>" converts "from" to file "to" using the given formats'
     write(*,'(1x,a)')&
          '"exportwf" <n>[u,d] <from.[bin,formatted,etsf]> "'//&
          ' converts n-th wavefunction of file "from" to cube using BigDFT uncompression'
     write(*,'(1x,a)')&
          '"atwf" <ng> calculates the atomic wavefunctions of the first atom in the gatom basis and write their expression '
     write(*,'(1x,a)')&
          '            in the "gatom-wfn.dat" file '
     write(*,'(1x,a)')&
          '           <ng> is the number of gaussians used for the gatom calculation'
     stop
  else
     read(unit=tatonam,fmt=*) nproc
     i_arg = 2
     output_grid=0
     loop_getargs: do
        call get_command_argument(i_arg, value = tatonam, status = istat)
        !call getarg(i_arg,tatonam)
        if(trim(tatonam)=='' .or. istat > 0) then
           exit loop_getargs
        else if (trim(tatonam)=='y') then
           output_grid=1
           write(*,'(1x,a)') 'The system grid will be displayed in the "grid.xyz" file'
           exit loop_getargs
        else if (trim(tatonam)=='o') then
           optimise=.true.
           output_grid=1
           write(*,'(1x,a)')&
                'The optimised system grid will be displayed in the "grid.xyz" file'
           exit loop_getargs
        else if (trim(tatonam)=='GPUtest') then
           GPUtest=.true.
           write(*,'(1x,a)')&
                'Perform the test with GPU, if present.'
           i_arg = i_arg + 1
           call get_command_argument(i_arg, value = tatonam, status = istat)
           !call getarg(i_arg,tatonam)
           ntimes=1
           norbgpu=0
           read(tatonam,*,iostat=ierror)ntimes
           if (ierror==0) then
              write(*,'(1x,a,i0,a)')&
                   'Repeat each calculation ',ntimes,' times.'
              i_arg = i_arg + 1
              call get_command_argument(i_arg, value = tatonam)
              !call getarg(i_arg,tatonam)
              read(tatonam,*,iostat=ierror)norbgpu
           end if
           exit loop_getargs
        else if (trim(tatonam)=='convert') then
           convert=.true.
           i_arg = i_arg + 1
           call get_command_argument(i_arg, value = fileFrom)
           !call getarg(i_arg,fileFrom)
           i_arg = i_arg + 1
           call get_command_argument(i_arg, value = fileTo)
           !call getarg(i_arg,fileTo)
           write(*,'(1x,5a)')&
                'convert "', trim(fileFrom),'" file to "', trim(fileTo),'"'
           exit loop_getargs
        else if (trim(tatonam)=='exportwf') then
           exportwf=.true.
           i_arg = i_arg + 1
           call get_command_argument(i_arg, value = filename_wfn)
           !call getarg(i_arg,filename_wfn)
           write(*,'(1x,3a)')&
                'export wavefunction file: "', trim(filename_wfn),'" in .cube format'
           exit loop_getargs
        else if (trim(tatonam)=='atwf') then
           atwf=.true.
           write(*,'(1x,a)')&
                'Perform the calculation of atomic wavefunction of the first atom'
           i_arg = i_arg + 1
           call get_command_argument(i_arg, value = tatonam)
           !call getarg(i_arg,tatonam)
           read(tatonam,*,iostat=ierror)ng
           write(*,'(1x,a,i0,a)')&
                'Use gaussian basis of',ng,' elements.'
           exit loop_getargs
        else 
           ! Use value as radical for input files.
           if (trim(radical) /= "") then
              write(*,'(1x,a)')&
                   'Usage: ./memguess <nproc> [y]'
              write(*,'(1x,a)')&
                   'Indicate the number of processes after the executable'
              write(*,'(1x,a)')&
                   'ERROR: The only second argument which is accepted is "y", "o","convert", "GPUtest" or "atwf" ' 
              write(*,'(1x,a)')&
                   '       (type "memguess" without arguments to have an help)'
              stop
           end if
           write(radical, "(A)") trim(tatonam)
        end if
        i_arg = i_arg + 1
     end do loop_getargs
  end if

!!!  open(unit=1,file='input.memguess',status='old')
!!!  
!!!  !line number, to control the input values
!!!  iline=0
!!!  
!!!  !number of MPI proccessors
!!!  read(1,*) nproc
!!!  write(*,*) 'Number of mpi processes is: ',nproc
!!!  
!!!  read(1,*) optimise
!!!  if (optimise) write(*,*) 'Molecule will be rotated to minimize simulation box size and workarrays in BigDFT'
!!!  
!!!  !    "T"  If the system grid is to be displayed in the "grid.xyz" file
!!!  read(1,*) output_grid
!!!  write(*,*)  'output_grid= ',output_grid
!!!  
!!!  !    "T"   'Perform the test with GPU, if present.'   
!!!  read(1,*) GPUtest
!!!  if (GPUtest) write(*,*) 'Perform the test with GPU'
!!!!!! END of By Ali



  !welcome screen
  !call print_logo()

  if (convert) then
     atoms%geocode = "P"
     write(*,*) "Read density file..."
     call read_density(trim(fileFrom), atoms%geocode, Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, &
          & nspin, hx, hy, hz, rhocoeff, atoms%nat, rxyz, atoms%iatype, atoms%nzatom)
     atoms%ntypes = size(atoms%nzatom) - ndebug
     write(*,*) "Write new density file..."
     call plot_density(trim(fileTo), 0, 1, Glr%d%n1i / 2 - 1, Glr%d%n2i / 2 - 1, &
          & Glr%d%n3i / 2 - 1, Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, Glr%d%n3i, nspin, hx, hy, hz, &
          & atoms, rxyz, norbsc_arr, rhocoeff)
     write(*,*) "Done"
     stop
  end if

  !standard names
  call standard_inputfile_names(in, radical)
  if (trim(radical) == "") then
     call read_input_variables(0, "posinp", in, atoms, rxyz)
  else
     call read_input_variables(0, trim(radical), in, atoms, rxyz)
  end if
  !initialize memory counting
  !call memocc(0,0,'count','start')

  if (in%ixc < 0) then
     call xc_init(in%ixc, XC_MIXED, nspin)
  else
     call xc_init(in%ixc, XC_ABINIT, nspin)
  end if

  call print_general_parameters(nproc,in,atoms)
  call print_dft_parameters(in,atoms)
  call xc_dump()

  write(*,'(1x,a)')&
       '------------------------------------------------------------------ System Properties'

  ! store PSP parameters
  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

  call system_properties(0,nproc,in,atoms,orbs,radii_cf,nelec)

  if (optimise) then
     if (atoms%geocode =='F') then
        call optimise_volume(atoms,in%crmult,in%frmult,in%hx,in%hy,in%hz,rxyz,radii_cf)
     else
        call shift_periodic_directions(atoms,rxyz,radii_cf)
     end if
     write(*,'(1x,a)')'Writing optimised positions in file posopt.[xyz,ascii]...'
     write(comment,'(a)')'POSITIONS IN OPTIMIZED CELL '
     call write_atomic_file('posopt',0.d0,rxyz,atoms,trim(comment))
     !call wtxyz('posopt',0.d0,rxyz,atoms,trim(comment))
  end if

  !in the case in which the number of orbitals is not "trivial" check whether they are too many
  if ( max(orbs%norbu,orbs%norbd) /= ceiling(real(nelec,kind=4)/2.0) .or. .true.) then
     ! Allocations for readAtomicOrbitals (check inguess.dat and psppar files + give norbe)
     allocate(scorb(4,2,atoms%natsc+ndebug),stat=i_stat)
     call memocc(i_stat,scorb,'scorb',subname)
     allocate(norbsc_arr(atoms%natsc+1,in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
     allocate(locrad(atoms%nat+ndebug),stat=i_stat)
     call memocc(i_stat,locrad,'locrad',subname)

     !calculate the inputguess orbitals
     !spin for inputguess orbitals
     if (in%nspin==4) then
        nspin_ig=1
     else
        nspin_ig=in%nspin
     end if

     ! Read the inguess.dat file or generate the input guess via the inguess_generator
     call readAtomicOrbitals(atoms,norbe,norbsc,nspin_ig,orbs%nspinor,&
          scorb,norbsc_arr,locrad)

     if (in%nspin==4) then
        !in that case the number of orbitals doubles
        norbe=2*norbe
     end if

     ! De-allocations
     i_all=-product(shape(locrad))*kind(locrad)
     deallocate(locrad,stat=i_stat)
     call memocc(i_stat,i_all,'locrad',subname)
     i_all=-product(shape(scorb))*kind(scorb)
     deallocate(scorb,stat=i_stat)
     call memocc(i_stat,i_all,'scorb',subname)
     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
     deallocate(norbsc_arr,stat=i_stat)
     call memocc(i_stat,i_all,'norbsc_arr',subname)

     ! Check the maximum number of orbitals
     if (in%nspin==1 .or. in%nspin==4) then
        if (orbs%norb>norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals (',orbs%norb,&
                ') must not be greater than the number of orbitals (',norbe,&
                ') generated from the input guess.'
           stop
        end if
     else if (in%nspin == 2) then
        if (orbs%norbu > norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals up (',orbs%norbu,&
                ') must not be greater than the number of orbitals (',norbe,&
                ') generated from the input guess.'
           stop
        end if
        if (orbs%norbd > norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals down (',orbs%norbd,&
                ') must not be greater than the number of orbitals (',norbe,&
                ') generated from the input guess.'
           stop
        end if
     end if
  end if

  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space
  hx=in%hx
  hy=in%hy
  hz=in%hz

  call system_size(0,atoms,rxyz,radii_cf,in%crmult,in%frmult,hx,hy,hz,Glr,shift)

  ! Build and print the communicator scheme.
  call createWavefunctionsDescriptors(0,hx,hy,hz,&
       atoms,rxyz,radii_cf,in%crmult,in%frmult,Glr, output_grid = (output_grid > 0))
  call orbitals_communicators(0,nproc,Glr,orbs,comms)  

  if (exportwf) then

       allocate(psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
       call memocc(i_stat,psi,'psi',subname)

       call take_psi_from_file(filename_wfn,in%hx,in%hy,in%hz,Glr,atoms,rxyz,psi)

       call plot_wf(filename_wfn,1,atoms,Glr,in%hx,in%hy,in%hz,rxyz,psi,' ')
  
       i_all=-product(shape(psi))*kind(psi)
       deallocate(psi,stat=i_stat)
       call memocc(i_stat,i_all,'psi',subname)

  end if

  if (GPUtest) then
     !test the hamiltonian in CPU or GPU
     !create the orbitals data structure for one orbital
     !test orbitals
     nspin=1
     if (norbgpu == 0) then
        norb=orbs%norb
     else
        norb=norbgpu
     end if
     norbu=norb
     norbd=0
     nspinor=1

     call orbitals_descriptors(0,nproc,norb,norbu,norbd,in%nspin,nspinor, &
          & in%nkpt,in%kpt,in%wkpt,orbstst)
     allocate(orbstst%eval(orbstst%norbp+ndebug),stat=i_stat)
     call memocc(i_stat,orbstst%eval,'orbstst%eval',subname)
     do iorb=1,orbstst%norbp
        orbstst%eval(iorb)=-0.5_gp
     end do

     do iorb=1,orbstst%norb
        orbstst%occup(iorb)=1.0_gp
        orbstst%spinsgn(iorb)=1.0_gp
     end do

     call compare_cpu_gpu_hamiltonian(0,1,in%iacceleration,atoms,orbstst,nspin,in%ncong,in%ixc,&
          Glr,hx,hy,hz,rxyz,ntimes)

     call deallocate_orbs(orbstst,subname)


     i_all=-product(shape(orbstst%eval))*kind(orbstst%eval)
     deallocate(orbstst%eval,stat=i_stat)
     call memocc(i_stat,i_all,'orbstst%eval',subname)

  end if

  call deallocate_comms(comms,subname)

  ! determine localization region for all projectors, but do not yet fill the descriptor arrays
  allocate(nlpspd%nseg_p(0:2*atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nseg_p,'nseg_p',subname)
  allocate(nlpspd%nvctr_p(0:2*atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nvctr_p,'nvctr_p',subname)
  allocate(nlpspd%nboxp_c(2,3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nboxp_c,'nboxp_c',subname)
  allocate(nlpspd%nboxp_f(2,3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nboxp_f,'nboxp_f',subname)

  allocate(logrid(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid,'logrid',subname)

  call localize_projectors(0,Glr%d%n1,Glr%d%n2,Glr%d%n3,hx,hy,hz,&
       in%frmult,in%frmult,rxyz,radii_cf,logrid,atoms,orbs,nlpspd)
  
  !allocations for arrays holding the data descriptors
  !just for modularity
  allocate(nlpspd%keyg_p(2,nlpspd%nseg_p(2*atoms%nat)+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%keyg_p,'nlpspd%keyg_p',subname)
  allocate(nlpspd%keyv_p(nlpspd%nseg_p(2*atoms%nat)+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%keyv_p,'nlpspd%keyv_p',subname)

  if (atwf) then
     !here the treatment of the AE Core charge density
     !number of gaussians defined in the input of memguess
     !ng=31
     !plot the wavefunctions for the pseudo atom
     nullify(G%rxyz)
     call gaussian_pswf_basis(ng,.false.,0,in%nspin,atoms,rxyz,G,gbd_occ)
     !for the moment multiply the number of coefficients for each channel
     allocate(rhocoeff((ng*(ng+1))/2,4+ndebug),stat=i_stat)
     call memocc(i_stat,rhocoeff,'rhocoeff',subname)
     allocate(rhoexpo((ng*(ng+1))/2+ndebug),stat=i_stat)
     call memocc(i_stat,rhoexpo,'rhoexpo',subname)

     call plot_gatom_basis('gatom',1,ng,G,gbd_occ,rhocoeff,rhoexpo)

     if (associated(gbd_occ)) then
        i_all=-product(shape(gbd_occ))*kind(gbd_occ)
        deallocate(gbd_occ,stat=i_stat)
        call memocc(i_stat,i_all,'gbd_occ',subname)
        nullify(gbd_occ)
     end if
     !deallocate the gaussian basis descriptors
     call deallocate_gwf(G,subname)

!!$  !plot the wavefunctions for the AE atom
!!$  !not possible, the code should recognize the AE eleconf
!!$  call razero(35,atoms%psppar(0,0,atoms%iatype(1)))
!!$  atoms%psppar(0,0,atoms%iatype(1))=0.01_gp
!!$  nullify(G%rxyz)
!!$  call gaussian_pswf_basis(ng,.false.,0,in%nspin,atoms,rxyz,G,gbd_occ)
!!$  !for the moment multiply the number of coefficients for each channel
!!$  allocate(rhocoeff((ng*(ng+1))/2,4+ndebug),stat=i_stat)
!!$  call memocc(i_stat,rhocoeff,'rhocoeff',subname)
!!$  allocate(rhoexpo((ng*(ng+1))/2+ndebug),stat=i_stat)
!!$  call memocc(i_stat,rhoexpo,'rhoexpo',subname)
!!$  
!!$  call plot_gatom_basis('all-elec',1,ng,G,gbd_occ,rhocoeff,rhoexpo)
!!$
!!$  if (associated(gbd_occ)) then
!!$     i_all=-product(shape(gbd_occ))*kind(gbd_occ)
!!$     deallocate(gbd_occ,stat=i_stat)
!!$     call memocc(i_stat,i_all,'gbd_occ',subname)
!!$     nullify(gbd_occ)
!!$  end if
!!$  !deallocate the gaussian basis descriptors
!!$  call deallocate_gwf(G,subname)

     i_all=-product(shape(rhoexpo))*kind(rhoexpo)
     deallocate(rhoexpo,stat=i_stat)
     call memocc(i_stat,i_all,'rhoexpo',subname)
     i_all=-product(shape(rhocoeff))*kind(rhocoeff)
     deallocate(rhocoeff,stat=i_stat)
     call memocc(i_stat,i_all,'rhocoeff',subname)

  end if

  i_all=-product(shape(logrid))*kind(logrid)
  deallocate(logrid,stat=i_stat)
  call memocc(i_stat,i_all,'logrid',subname)

  call deallocate_proj_descr(nlpspd,subname)
  call deallocate_atoms_scf(atoms,subname) 

  call MemoryEstimator(nproc,in%idsx,Glr,&
       atoms%nat,orbs%norb,orbs%nspinor,orbs%nkpts,nlpspd%nprojel,&
       in%nspin,in%itrpmax,in%iscf,peakmem)

  !add the comparison between cuda hamiltonian and normal one if it is the case

  call deallocate_atoms(atoms,subname)

  call deallocate_lr(Glr,subname)

  call xc_end()

  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf',subname)
  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz',subname)

  ! De-allocations
  call deallocate_orbs(orbs,subname)
  call free_input_variables(in)  

  !finalize memory counting
  call memocc(0,0,'count','stop')

end program memguess

  
!>  Rotate the molecule via an orthogonal matrix in order to minimise the
!!  volume of the cubic cell
subroutine optimise_volume(atoms,crmult,frmult,hx,hy,hz,rxyz,radii_cf)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  real(gp), intent(in) :: crmult,frmult,hx,hy,hz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
  !local variables
  character(len=*), parameter :: subname='optimise_volume'
  integer :: iat,i_all,i_stat,it,i
  real(gp) :: x,y,z,vol,tx,ty,tz,tvol,s,diag,dmax
  type(locreg_descriptors) :: Glr
  real(gp), dimension(3) :: shift
  real(gp), dimension(3,3) :: urot
  real(gp), dimension(:,:), allocatable :: txyz

  allocate(txyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,txyz,'txyz',subname)
  call system_size(1,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,Glr,shift)
  !call volume(nat,rxyz,vol)
  vol=atoms%alat1*atoms%alat2*atoms%alat3
  write(*,'(1x,a,1pe16.8)')'Initial volume (Bohr^3)',vol

  it=0
  diag=1.d-2 ! initial small diagonal element allows for search over all angles
  loop_rotations: do  ! loop over all trial rotations
     diag=diag*1.0001_gp ! increase diag to search over smaller angles
     it=it+1
     if (diag > 100._gp) exit loop_rotations ! smaller angle rotations do not make sense

     ! create a random orthogonal (rotation) matrix
     call random_number(urot)
     urot(:,:)=urot(:,:)-.5_gp
     do i=1,3
        urot(i,i)=urot(i,i)+diag
     enddo

     s=urot(1,1)**2+urot(2,1)**2+urot(3,1)**2
     s=1._gp/sqrt(s)
     urot(:,1)=s*urot(:,1) 

     s=urot(1,1)*urot(1,2)+urot(2,1)*urot(2,2)+urot(3,1)*urot(3,2)
     urot(:,2)=urot(:,2)-s*urot(:,1)
     s=urot(1,2)**2+urot(2,2)**2+urot(3,2)**2
     s=1._gp/sqrt(s)
     urot(:,2)=s*urot(:,2) 

     s=urot(1,1)*urot(1,3)+urot(2,1)*urot(2,3)+urot(3,1)*urot(3,3)
     urot(:,3)=urot(:,3)-s*urot(:,1)
     s=urot(1,2)*urot(1,3)+urot(2,2)*urot(2,3)+urot(3,2)*urot(3,3)
     urot(:,3)=urot(:,3)-s*urot(:,2)
     s=urot(1,3)**2+urot(2,3)**2+urot(3,3)**2
     s=1._gp/sqrt(s)
     urot(:,3)=s*urot(:,3) 

     ! eliminate reflections
     if (urot(1,1) <= 0._gp) urot(:,1)=-urot(:,1)
     if (urot(2,2) <= 0._gp) urot(:,2)=-urot(:,2)
     if (urot(3,3) <= 0._gp) urot(:,3)=-urot(:,3)

     ! apply the rotation to all atomic positions! 
     do iat=1,atoms%nat
        x=rxyz(1,iat) 
        y=rxyz(2,iat) 
        z=rxyz(3,iat)

        txyz(:,iat)=x*urot(:,1)+y*urot(:,2)+z*urot(:,3)
     enddo

     call system_size(1,atoms,txyz,radii_cf,crmult,frmult,hx,hy,hz,Glr,shift)
     tvol=atoms%alat1*atoms%alat2*atoms%alat3
     !call volume(nat,txyz,tvol)
     if (tvol < vol) then
        write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol,it,diag
        rxyz(:,:)=txyz(:,:)
        vol=tvol
        dmax=max(atoms%alat1,atoms%alat2,atoms%alat3)
        ! if box longest along x switch x and z
        if (atoms%alat1 == dmax)  then
           do  iat=1,atoms%nat
              tx=rxyz(1,iat)
              tz=rxyz(3,iat)

              rxyz(1,iat)=tz
              rxyz(3,iat)=tx
           enddo
           ! if box longest along y switch y and z
        else if (atoms%alat2 == dmax .and. atoms%alat1 /= dmax)  then
           do  iat=1,atoms%nat
              ty=rxyz(2,iat) 
              tz=rxyz(3,iat)

              rxyz(2,iat)=tz 
              rxyz(3,iat)=ty
           enddo
        endif
     endif
  end do loop_rotations

  i_all=-product(shape(txyz))*kind(txyz)
  deallocate(txyz,stat=i_stat)
  call memocc(i_stat,i_all,'txyz',subname)

END SUBROUTINE optimise_volume


!>  Add a shift in the periodic directions such that the system
!!  uses as less as possible the modulo operation
subroutine shift_periodic_directions(at,rxyz,radii_cf)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: at
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  real(gp), dimension(3,at%nat), intent(inout) :: rxyz
  !local variables
  character(len=*), parameter :: subname='shift_periodic_directions'
  integer :: iat,i_all,i_stat,i,ityp
  real(gp) :: vol,tvol,maxsh,shiftx,shifty,shiftz
  real(gp), dimension(:,:), allocatable :: txyz

  !calculate maximum shift between these values
  !this is taken as five times the coarse radius around atoms
  maxsh=0.0_gp
  do ityp=1,at%ntypes
     maxsh=max(maxsh,5_gp*radii_cf(ityp,1))
  end do
  
  allocate(txyz(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,txyz,'txyz',subname)

  call calc_vol(at%geocode,at%nat,rxyz,vol)

  if (at%geocode /= 'F') then
     loop_shiftx: do i=1,5000 ! loop over all trial rotations
        ! create a random orthogonal (rotation) matrix
        call random_number(shiftx)

        !apply the shift to all atomic positions taking into account the modulo operation
        do iat=1,at%nat
           txyz(1,iat)=modulo(rxyz(1,iat)+shiftx*maxsh,at%alat1)
           txyz(2,iat)=rxyz(2,iat)
           txyz(3,iat)=rxyz(3,iat)
        end do

        call calc_vol(at%geocode,at%nat,txyz,tvol)
        !print *,'vol',tvol

        if (tvol < vol) then
           write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol
           rxyz(:,:)=txyz(:,:)
           vol=tvol
        endif
     end do loop_shiftx
  end if

  if (at%geocode == 'P') then
     loop_shifty: do i=1,5000 ! loop over all trial rotations
        ! create a random orthogonal (rotation) matrix
        call random_number(shifty)

        !apply the shift to all atomic positions taking into account the modulo operation
        do iat=1,at%nat
           txyz(1,iat)=rxyz(1,iat)
           txyz(2,iat)=modulo(rxyz(2,iat)+shifty*maxsh,at%alat2)
           txyz(3,iat)=rxyz(3,iat)
        end do

        call calc_vol(at%geocode,at%nat,txyz,tvol)

        if (tvol < vol) then
           write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol
           rxyz(:,:)=txyz(:,:)
           vol=tvol
        endif
     end do loop_shifty
  end if

    if (at%geocode /= 'F') then
     loop_shiftz: do i=1,5000 ! loop over all trial rotations
        ! create a random orthogonal (rotation) matrix
        call random_number(shiftz)

        !apply the shift to all atomic positions taking into account the modulo operation
        do iat=1,at%nat
           txyz(1,iat)=rxyz(1,iat)
           txyz(2,iat)=rxyz(2,iat)
           txyz(3,iat)=modulo(rxyz(3,iat)+shiftz*maxsh,at%alat3)
        end do

        call calc_vol(at%geocode,at%nat,txyz,tvol)

        if (tvol < vol) then
           write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol
           rxyz(:,:)=txyz(:,:)
           vol=tvol
        endif
     end do loop_shiftz
  end if

  


  i_all=-product(shape(txyz))*kind(txyz)
  deallocate(txyz,stat=i_stat)
  call memocc(i_stat,i_all,'txyz',subname)

END SUBROUTINE shift_periodic_directions


subroutine calc_vol(geocode,nat,rxyz,vol)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: nat
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(gp), intent(out) :: vol
  !local variables
  integer :: iat
  real(gp) :: cxmin,cxmax,cymin,cymax,czmin,czmax

  !calculate the extremes of the boxes taking into account the spheres around the atoms
  cxmax=-1.e10_gp 
  cxmin=1.e10_gp

  cymax=-1.e10_gp 
  cymin=1.e10_gp

  czmax=-1.e10_gp 
  czmin=1.e10_gp

  do iat=1,nat
     cxmax=max(cxmax,rxyz(1,iat)) 
     cxmin=min(cxmin,rxyz(1,iat))

     cymax=max(cymax,rxyz(2,iat)) 
     cymin=min(cymin,rxyz(2,iat))
     
     czmax=max(czmax,rxyz(3,iat)) 
     czmin=min(czmin,rxyz(3,iat))
  enddo
  !print *,cxmax,cxmin,cymax,cymin,czmax,czmin
  !now calculate the volume for the periodic part
  if (geocode == 'P') then
     vol=(cxmax-cxmin)*(cymax-cymin)*(czmax-czmin)
  else if (geocode == 'S') then
     vol=(cxmax-cxmin)*(czmax-czmin)
  end if

END SUBROUTINE calc_vol


subroutine compare_cpu_gpu_hamiltonian(iproc,nproc,iacceleration,at,orbs,nspin,ixc,ncong,&
     lr,hx,hy,hz,rxyz,ntimes)
  use module_base
  use module_types
  use module_interfaces
  use Poisson_Solver
  use module_xc

  implicit none
  integer, intent(in) :: iproc,nproc,nspin,ncong,ixc,ntimes,iacceleration
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  !local variables
  character(len=*), parameter :: subname='compare_cpu_gpu_hamiltonian'
  logical :: rsflag
  integer :: icoeff,i_stat,i_all,i1,i2,i3,ispin,j
  integer :: iorb,n3d,n3p,n3pi,i3xcsh,i3s,jproc,nrhotot,nspinn,nvctrp
  integer(kind=8) :: itsc0,itsc1
  real(kind=4) :: tt,t0,t1
  real(gp) :: ttd,x,y,z,r2,arg,sigma2,ekin_sum,epot_sum,ekinGPU,epotGPU,gnrm,gnrm_zero,gnrmGPU
  real(gp) :: Rden,Rham,Rgemm,Rsyrk,Rprec,eSIC_DC
  real(kind=8) :: CPUtime,GPUtime
  type(gaussian_basis) :: G
  type(GPU_pointers) :: GPU
  integer, dimension(:,:), allocatable :: nscatterarr
  real(wp), dimension(:,:,:,:), allocatable :: pot,rho
  real(wp), dimension(:,:), allocatable :: gaucoeffs,psi,hpsi
  real(wp), dimension(:,:,:), allocatable :: overlap
  real(wp), dimension(:), pointer :: gbd_occ
  real(wp), dimension(:), pointer :: fake_pkernelSIC

  !nullify pkernelSIC pointer
  nullify(fake_pkernelSIC)

  !nullify the G%rxyz pointer
  nullify(G%rxyz)
  !extract the gaussian basis from the pseudowavefunctions
  call gaussian_pswf_basis(21,.false.,iproc,nspin,at,rxyz,G,gbd_occ)
  
  allocate(gaucoeffs(G%ncoeff,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,gaucoeffs,'gaucoeffs',subname)

  !fill randomly the gaussian coefficients for the orbitals considered
  do iorb=1,orbs%norbp*orbs%nspinor
     do icoeff=1,G%ncoeff
        call random_number(tt)
        gaucoeffs(icoeff,iorb)=real(tt,wp)
     end do
  end do

  !allocate the wavefunctions
  allocate(psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  allocate(hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)

  call razero(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f*orbs%nspinor*orbs%norbp,psi)

  !convert the gaussians in wavelets
  call gaussians_to_wavelets(iproc,nproc,at%geocode,orbs,lr%d,&
       hx,hy,hz,lr%wfd,G,gaucoeffs,psi)

  i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
  deallocate(gaucoeffs,stat=i_stat)
  call memocc(i_stat,i_all,'gaucoeffs',subname)

  i_all=-product(shape(gbd_occ))*kind(gbd_occ)
  deallocate(gbd_occ,stat=i_stat)
  call memocc(i_stat,i_all,'gbd_occ',subname)

  !deallocate the gaussian basis descriptors
  call deallocate_gwf(G,subname)

  !allocate and initialise the potential and the density
  allocate(pot(lr%d%n1i,lr%d%n2i,lr%d%n3i,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,pot,'pot',subname)
  allocate(rho(lr%d%n1i,lr%d%n2i,lr%d%n3i,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,rho,'rho',subname)

  !here the potential can be used for building the density
  allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,nscatterarr,'nscatterarr',subname)
  !normally nproc=1
  do jproc=0,nproc-1
     call PS_dim4allocation(at%geocode,'D',jproc,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,ixc,&
          n3d,n3p,n3pi,i3xcsh,i3s)
     nscatterarr(jproc,1)=n3d
     nscatterarr(jproc,2)=n3p
     nscatterarr(jproc,3)=i3s+i3xcsh-1
     nscatterarr(jproc,4)=i3xcsh
  end do

  !components of the charge density
  if (orbs%nspinor ==4) then
     nspinn=4
  else
     nspinn=nspin
  end if

  !flag for toggling the REDUCE_SCATTER stategy
  rsflag = .not.xc_isgga()

  !calculate dimensions of the complete array to be allocated before the reduction procedure
  if (rsflag) then
     nrhotot=0
     do jproc=0,nproc-1
        nrhotot=nrhotot+nscatterarr(jproc,1)
     end do
  else
     nrhotot=lr%d%n3i
  end if

  !allocate the necessary objects on the GPU
  !set initialisation of GPU part 
  !initialise the acceleration strategy if required
  call init_material_acceleration(iproc,iacceleration,GPU)
  
  if (GPUconv .eqv. OCLconv) stop 'ERROR: One (and only one) acceleration should be present with GPUtest'
  
  !allocate arrays for the GPU if a card is present
  if (GPUconv) then
     call prepare_gpu_for_locham(lr%d%n1,lr%d%n2,lr%d%n3,nspin,&
          hx,hy,hz,lr%wfd,orbs,GPU)
  else if (OCLconv) then
     !the same with OpenCL, but they cannot exist at same time
     call allocate_data_OCL(lr%d%n1,lr%d%n2,lr%d%n3,lr%geocode,&
          nspin,hx,hy,hz,lr%wfd,orbs,GPU)
  end if
  if (iproc == 0) write(*,*)&
       'GPU data allocated'

  write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Density calculation'

  !for each of the orbitals treated by the processor build the partial densities
  !call cpu_time(t0)
  call nanosec(itsc0)
  do j=1,ntimes
     call tenminustwenty(lr%d%n1i*lr%d%n2i*nrhotot*nspinn,pot,nproc)
     call local_partial_density(iproc,nproc,rsflag,nscatterarr,&
          nrhotot,lr,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,nspin,orbs,&
          psi,pot)
  end do
  call nanosec(itsc1)
  !call cpu_time(t1)
  !CPUtime=real(t1-t0,kind=8)
  CPUtime=real(itsc1-itsc0,kind=8)*1.d-9

  !now the GPU part
  !for each of the orbitals treated by the processor build the partial densities
  !call cpu_time(t0)
  call nanosec(itsc0)
  do j=1,ntimes
     !switch between GPU/CPU treatment of the density
     if (GPUconv) then
        call local_partial_density_GPU(iproc,nproc,orbs,nrhotot,lr,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,nspin,psi,rho,GPU)
     else if (OCLconv) then
        call local_partial_density_OCL(iproc,nproc,orbs,nrhotot,lr,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,nspin,psi,rho,GPU)
     end if
  end do
  call nanosec(itsc1)
  !call cpu_time(t1)
  !GPUtime=real(t1-t0,kind=8)
  GPUtime=real(itsc1-itsc0,kind=8)*1.d-9

  i_all=-product(shape(nscatterarr))*kind(nscatterarr)
  deallocate(nscatterarr,stat=i_stat)
  call memocc(i_stat,i_all,'nscatterarr',subname)

  !compare the results between the different actions of the hamiltonian
  !check the differences between the results
  call compare_data_and_gflops(CPUtime,GPUtime,&
       real(lr%d%n1i*lr%d%n2i*lr%d%n3i,kind=8)*192.d0,pot,rho,&
       lr%d%n1i*lr%d%n2i*lr%d%n3i,ntimes*orbs%norbp,.false.,Rden)

  i_all=-product(shape(rho))*kind(rho)
  deallocate(rho,stat=i_stat)
  call memocc(i_stat,i_all,'rho',subname)


  !here the grid spacings are the small ones
  sigma2=0.125_gp*((lr%d%n1i*hx)**2+(lr%d%n2i*hy)**2+(lr%d%n3i*hz)**2)
  do ispin=1,nspin
     do i3=1,lr%d%n3i
        z=hz*real(i3-lr%d%n3i/2-1,gp)
        do i2=1,lr%d%n2i
           y=hy*real(i2-lr%d%n2i/2-1,gp)
           do i1=1,lr%d%n1i
              x=hx*real(i1-lr%d%n1i/2-1,gp)
              !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
              r2=x**2+y**2+z**2
              arg=0.5d0*r2/sigma2
              ttd=dexp(-arg)

              pot(i1,i2,i3,ispin)=ttd
           end do
        end do
     end do
  end do

  write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Local Hamiltonian calculation'

  !warm-up
  !call local_hamiltonian(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum) 

  !apply the CPU hamiltonian
  !take timings
  call nanosec(itsc0)
  do j=1,ntimes
     call local_hamiltonian(iproc,orbs,lr,hx,hy,hz,0,pot,psi,hpsi,fake_pkernelSIC,0,0.0_gp,ekin_sum,epot_sum,eSIC_DC)
  end do
  call nanosec(itsc1)
  CPUtime=real(itsc1-itsc0,kind=8)*1.d-9

  print *,'ekin,epot=',ekin_sum,epot_sum

  !WARNING: local hamiltonian overwrites the psis
  !warm-up
  !call gpu_locham(lr%d%n1,lr%d%n2,lr%d%n3,hx,hy,hz,orbs,GPU,ekinGPU,epotGPU)

  !apply the GPU hamiltonian and put the results in the hpsi_GPU array
  allocate(GPU%hpsi_ASYNC((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),stat=i_stat)
  call memocc(i_stat,GPU%hpsi_ASYNC,'GPU%hpsi_ASYNC',subname)

  !take timings
  call nanosec(itsc0)
  do j=1,ntimes
     if (GPUconv) then
        call local_hamiltonian_GPU(iproc,orbs,lr,hx,hy,hz,orbs%nspin,pot,psi,GPU%hpsi_ASYNC,ekinGPU,epotGPU,GPU)
     else if (OCLconv) then
        call local_hamiltonian_OCL(iproc,orbs,lr,hx,hy,hz,orbs%nspin,pot,psi,GPU%hpsi_ASYNC,ekinGPU,epotGPU,GPU)
     end if
  end do
  if(ASYNCconv .and. OCLconv) call finish_hamiltonian_OCL(orbs,ekinGPU,epotGPU,GPU)
  call nanosec(itsc1)
  GPUtime=real(itsc1-itsc0,kind=8)*1.d-9

  print *,'ekinGPU,epotGPU',ekinGPU,epotGPU



  !compare the results between the different actions of the hamiltonian
  !check the differences between the results
  call compare_data_and_gflops(CPUtime,GPUtime,&
       real(lr%d%n1i*lr%d%n2i*lr%d%n3i,kind=8)*real(192+46*3+192+2,kind=8),hpsi,GPU%hpsi_ASYNC,&
       orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),ntimes*orbs%norbp,.false.,Rham)

  i_all=-product(shape(pot))*kind(pot)
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot',subname)

  write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Linear Algebra (Blas)'
 
  !perform the scalar product between the hpsi wavefunctions
  !actually this is <hpsi|hpsi> it has no meaning.
  !this works only if nspinor==1
  allocate(overlap(orbs%norbp,orbs%norbp,2+ndebug),stat=i_stat)
  call memocc(i_stat,overlap,'overlap',subname)

  call nanosec(itsc0)
  do j=1,ntimes
     call DGEMM('T','N',orbs%norbp,orbs%norbp,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),1.0_wp,&
          psi(1,1),(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),&
          hpsi(1,1),(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),0.0_wp,&
          overlap(1,1,1),orbs%norbp)
  end do
  call nanosec(itsc1)
  CPUtime=real(itsc1-itsc0,kind=8)*1.d-9


  call nanosec(itsc0)
  do j=1,ntimes
     call GEMMSY('T','N',orbs%norbp,orbs%norbp,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),1.0_wp,&
          psi(1,1),(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),&
          hpsi(1,1),(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),0.0_wp,&
          overlap(1,1,2),orbs%norbp)
  end do
  call nanosec(itsc1)
  GPUtime=real(itsc1-itsc0,kind=8)*1.d-9

  !comparison between the results
  call compare_data_and_gflops(CPUtime,GPUtime,&
       real(orbs%norbp**2*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*2,kind=8),overlap(1,1,1),overlap(1,1,2),&
       orbs%norbp**2,ntimes,.false.,Rgemm)


  nvctrp=lr%wfd%nvctr_c+7*lr%wfd%nvctr_f

  call nanosec(itsc0)
  do j=1,ntimes
     call dsyrk('L','T',orbs%norbp,nvctrp,1.0_wp,psi(1,1),nvctrp,0.0_wp,&
          overlap(1,1,1),orbs%norbp)
  end do
  call nanosec(itsc1)
  CPUtime=real(itsc1-itsc0,kind=8)*1.d-9


  call nanosec(itsc0)
  do j=1,ntimes
     call syrk('L','T',orbs%norbp,nvctrp,1.0_wp,psi(1,1),nvctrp,0.0_wp,&
          overlap(1,1,2),orbs%norbp)
  end do
  call nanosec(itsc1)
  GPUtime=real(itsc1-itsc0,kind=8)*1.d-9

  call compare_data_and_gflops(CPUtime,GPUtime,&
       real(orbs%norbp*(orbs%norbp+1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),kind=8),overlap(1,1,1),overlap(1,1,2),&
       orbs%norbp**2,ntimes,.false.,Rsyrk)

  i_all=-product(shape(overlap))*kind(overlap)
  deallocate(overlap,stat=i_stat)
  call memocc(i_stat,i_all,'overlap',subname)


  !-------------------now the same for preconditioning
  write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Preconditioner'

  !the input function is psi
  call nanosec(itsc0)
  do j=1,ntimes
     call preconditionall(iproc,nproc,orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero)
  end do
  call nanosec(itsc1)

  CPUtime=real(itsc1-itsc0,kind=8)*1.d-9
  print *,'gnrm',gnrm


  !GPU data are already on the card, must be only copied back
  !the input function is GPU%hpsi in that case
  call nanosec(itsc0)
  do j=1,ntimes
     !Preconditions all orbitals belonging to iproc
     !and calculate the partial norm of the residue
     !switch between CPU and GPU treatment
     if (GPUconv) then
        call preconditionall_GPU(iproc,nproc,orbs,lr,hx,hy,hz,ncong,&
             GPU%hpsi_ASYNC,gnrmGPU,gnrm_zero,GPU)
     else if (OCLconv) then
        call preconditionall_OCL(iproc,nproc,orbs,lr,hx,hy,hz,ncong,&
             GPU%hpsi_ASYNC,gnrmGPU,gnrm_zero,GPU)
     end if
  end do
  call nanosec(itsc1)
  
  GPUtime=real(itsc1-itsc0,kind=8)*1.d-9
  print *,'gnrmGPU',gnrmGPU

  call compare_data_and_gflops(CPUtime,GPUtime,&
       real(lr%d%n1i*lr%d%n2i*lr%d%n3i,kind=8)*real((192+46*3+192+2-1+12)*(ncong+1),kind=8),hpsi,GPU%hpsi_ASYNC,&
       orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),ntimes*orbs%norbp,.false.,Rprec)

  i_all=-product(shape(GPU%hpsi_ASYNC))*kind(GPU%hpsi_ASYNC)
  deallocate(GPU%hpsi_ASYNC,stat=i_stat)
  call memocc(i_stat,i_all,'GPU%hpsi_ASYNC',subname)
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)
  i_all=-product(shape(hpsi))*kind(hpsi)
  deallocate(hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi',subname)


  !free the card at the end
  if (GPUconv) then
     call free_gpu(GPU,orbs%norbp)
  else if (OCLconv) then
     call free_gpu_OCL(GPU,orbs,nspin)
  end if

  !finalise the material accelearion usage
  call release_material_acceleration(GPU)



  write(*,'(1x,a,5(1x,f7.3))')'Ratios:',Rden,Rham,Rgemm,Rsyrk,Rprec
  
END SUBROUTINE compare_cpu_gpu_hamiltonian


subroutine compare_data_and_gflops(CPUtime,GPUtime,GFlopsfactor,&
     CPUdata,GPUdata,n,ntimes,dowrite,ratio)
  use module_base
  implicit none
  logical, intent(in) :: dowrite
  integer, intent(in) :: n,ntimes
  real(gp), intent(in) :: CPUtime,GPUtime,GFlopsfactor
  real(gp), intent(out) :: ratio
  real(wp), dimension(n), intent(in) :: CPUdata,GPUdata
  !local variables
  integer :: i
  real(gp) :: CPUGflops,GPUGflops,maxdiff,comp,threshold

  threshold=1.d-12
  !un-initialize valies which might suffer from fpe
  GPUGflops=-1.0_gp
  CPUGflops=-1.0_gp
  ratio=-1.0_gp

  if (CPUtime > 0.0_gp) CPUGflops=GFlopsfactor*real(ntimes,gp)/(CPUtime*1.d9)
  if (GPUtime > 0.0_gp) GPUGflops=GFlopsfactor*real(ntimes,gp)/(GPUtime*1.d9)

  maxdiff=0.0_gp

  rewind(17)

  do i=1,n
     if (dowrite) write(17,'(i6,2(1pe24.17))')i,CPUdata(i),GPUdata(i)
     comp=abs(CPUdata(i)-GPUdata(i))
     maxdiff=max(maxdiff,comp)
  end do
  if (GPUtime > 0.0_gp) ratio=CPUtime/GPUtime
  write(*,'(1x,a)')'| CPU: ms  |  Gflops  || GPU:  ms |  GFlops  || Ratio  | No. Elements | Max. Diff. |'

  write(*,'(1x,2(2(a,f10.2),a),a,f8.3,a,i14,a,1pe12.4,a)',advance='no')&
                                !'nbelem,REF/TEST ratio,Time,Gflops: REF,TEST',&
       '|',CPUtime*1.d3/real(ntimes,kind=8),'|',& ! Time CPU (ms)
       CPUGflops,'|',& !Gflops CPU (ms)
       '|',GPUtime*1.d3/real(ntimes,kind=8),'|',& ! Time GPU (ms)
       GPUGflops,'|',&!Gflops GPU (ms)
       '|',ratio,'|',& ! ratio
       n,'|',& !No. elements
       maxdiff,'|' ! maxdiff
  if (maxdiff <= threshold) then
     write(*,'(a)')''
  else
     write(*,'(a)')'<<<< WARNING' 
  end if

END SUBROUTINE compare_data_and_gflops

!> Extract the compressed wavefunction from the given file 
subroutine take_psi_from_file(filename,hx,hy,hz,lr,at,rxyz,psi)
  use module_base
  use module_types
  implicit none
  real(gp), intent(in) :: hx,hy,hz
  character(len=*), intent(in) :: filename
  type(locreg_descriptors), intent(in) :: lr
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='take_psi_form_file'
  logical :: perx,pery,perz,exists
  integer :: nb1,nb2,nb3,i_stat,isuffix,iorb_out,i_all
  real(gp) :: eval_fake
  real(wp), dimension(:,:,:), allocatable :: psifscf
  real(gp), dimension(:,:), allocatable :: rxyz_file

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  !buffers realted to periodicity
  !WARNING: the boundary conditions are not assumed to change between new and old
  call ext_buffers_coarse(perx,nb1)
  call ext_buffers_coarse(pery,nb2)
  call ext_buffers_coarse(perz,nb3)

  allocate(psifscf(-nb1:2*lr%d%n1+1+nb1,-nb2:2*lr%d%n2+1+nb2,-nb3:2*lr%d%n3+1+nb3+ndebug),stat=i_stat)
  call memocc(i_stat,psifscf,'psifscf',subname)

  allocate(rxyz_file(at%nat,3+ndebug),stat=i_stat)
  call memocc(i_stat,rxyz_file,'rxyz_file',subname)
     
  isuffix = index(filename, ".bin", back = .true.)
  exists=(isuffix > 0) !the file is written in binary format
  if (exists) then
     write(*,*) "Reading wavefunctions in BigDFT binary file format."
     open(unit=99,file=trim(filename),status='unknown',form="unformatted")
  else
     write(*,*) "Reading wavefunctions in plain text file format."
     open(unit=99,file=trim(filename),status='unknown')
  end if

  !find the value of iorb_out
  read(filename(index(filename, ".", back = .true.)+1:len(filename)),*)iorb_out

  !@ todo geocode should be passed in the localisation regions descriptors
  call readonewave(99, .not. exists,iorb_out,0,lr%d%n1,lr%d%n2,lr%d%n3, &
       hx,hy,hz,at,lr%wfd,rxyz_file,rxyz,&
       psi,eval_fake,psifscf)

  i_all=-product(shape(psifscf))*kind(psifscf)
  deallocate(psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'psifscf',subname)

  i_all=-product(shape(rxyz_file))*kind(rxyz_file)
  deallocate(rxyz_file,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz_file',subname)

end subroutine take_psi_from_file
