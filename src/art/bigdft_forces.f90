!> @file
!!   Information to interface BigDFT with ART
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2009-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! Modified by:
!! -EM 2010, see ~/AUTHORS
!! -Laurent Karim Beland, UdeM, 2011. For working with QM/MM !!

!> ART Module bigdft_forces
!! Module which contains information for Bigdft run inside art
module bigdft_forces

   use module_base!, only : gp,wp,dp,bohr2ang
   use module_types
   use module_interfaces
   use defs, only : iproc
   implicit none

   private

   ! Storage of required variables for a SCF loop calculation.
   logical :: initialised = .false.
   logical :: first_time  = .True.
   integer :: nproc, me
   type(atoms_data) :: at, atoms_all !at are the quantum atoms
   type(input_variables) :: in
   type(restart_objects) :: rst
   real(gp), parameter :: ht2ev = 27.2113834_gp

   real(kind=8) :: gnrm_l
   real(kind=8) :: gnrm_h  ! For lanczos

   public :: bigdft_init
   public :: calcforce_bigdft
   public :: mingeo
   public :: bigdft_finalise
   public :: init_all_atoms
   public :: copy_atoms_object
   public :: prepare_quantum_atoms_Si

   ! development : At some key points we do not use the previously 
   ! calculated wave function
   logical, public :: new_wf
   public :: check_force_clean_wf
   integer, dimension(:), allocatable, public :: in_system   ! Constraint over atoms

   contains

   !> ART init_all_atoms
   !! Routine to initialize all the positions. uses BigDFT read files
   subroutine init_all_atoms( nat, typa, posa, const_, boxl, boxtype, nproc_, me_, file )

      implicit none

      !Arguments
      integer,      intent(out)               :: nat
      integer,      pointer                   :: typa(:)
      real(kind=8), pointer                   :: posa(:)
      integer,      pointer                   :: const_(:)
      real(kind=8), dimension(3), intent(out) :: boxl
      character(len=1), intent(out)           :: boxtype
      integer,      intent(in)                :: nproc_
      integer,      intent(in)                :: me_
      character(len=*), intent(in)            :: file

      !Local variables
      integer                                 :: i
      character(len=2)                        :: symbol
      character(len=10)                       :: name
      real(gp), dimension(:,:), pointer       :: rxyz
      !_______________________

      nproc = nproc_
      me = me_

      call read_atomic_file(file,me_,atoms_all,rxyz)
      nat = atoms_all%nat
      boxtype = atoms_all%geocode

      allocate(posa(3 * nat))
      allocate(typa(nat))
      allocate(const_(nat))

      do i = 1, nat, 1
         posa(i)           = rxyz(1, i) * bohr2ang
         posa(i + nat)     = rxyz(2, i) * bohr2ang
         posa(i + 2 * nat) = rxyz(3, i) * bohr2ang
         typa(i) = atoms_all%iatype(i)
      end do

      boxl(1) = atoms_all%alat1 * bohr2ang
      boxl(2) = atoms_all%alat2 * bohr2ang
      boxl(3) = atoms_all%alat3 * bohr2ang
      ! Blocked atoms 
      const_ = 0                          ! Initialization, everyone is free.
      const_(:) = atoms_all%ifrztyp(:)

      if ( iproc == 0 ) then
         do i=1,atoms_all%nat
            name=trim(atoms_all%atomnames(atoms_all%iatype(i)))
            if (name(3:3)=='_') then
               symbol=name(1:2)
            else if (name(2:2)=='_') then
               symbol=name(1:1)
            else
               symbol=name(1:2)
            end if
            !write(9,'(i3, a2,4x,3(1x,1pe24.17))')i, symbol,(rxyz(j,i),j=1,3)
         enddo
      end if

   END SUBROUTINE init_all_atoms


   !> ART bigdft_init
   !! Routine to initialize all BigDFT stuff
   subroutine bigdft_init( nat, me_, my_gnrm,passivate,total_nb_atoms )

      implicit none

      !Arguments
      integer,      intent(in) :: nat
      integer,      intent(in)  :: me_
      real(kind=8), intent(in)  :: my_gnrm
      logical,      intent(in)  :: passivate
      integer,      intent(in)  :: total_nb_atoms

      !Local variables
      character(len=*), parameter :: subname='bigdft_init'
      real(gp), dimension(:,:), pointer     :: rxyz
      real(gp),dimension(3*total_nb_atoms) :: posquant
      integer :: natoms_calcul
      !_______________________

      me = me_

      if (nat .eq. total_nb_atoms .and. .not. passivate) then 
         ! we just reread all atoms
         call read_atomic_file("posinp",me_,at,rxyz)
      else 
         !uses the big object to prepare. everything should
         ! be alright in the object exept the length
         call prepare_quantum_atoms_Si(atoms_all,posquant,natoms_calcul)
         !we just copy it in a smaller vect
         call copy_atoms_object(atoms_all,at,rxyz,natoms_calcul,total_nb_atoms,posquant)
         call initialize_atomic_file(me_,at,rxyz)
      endif
      !standard names
      call standard_inputfile_names(in,'input')
      ! Read inputs.
      call read_input_parameters(me_, in, at, rxyz)

      call init_atomic_values(me_, at, in%ixc)

      ! Transfer "at" data to ART variables.
      gnrm_l = in%gnrm_cv
      if ( my_gnrm == 1.0d0 ) then 
         gnrm_h = in%gnrm_cv 
      else
         gnrm_h = my_gnrm
      end if
      ! The BigDFT restart structure.
      call init_restart_objects(me, in%iacceleration, at, rst, subname)

   END SUBROUTINE bigdft_init

   !> ART calcforce_bigdft
   !! Calculation of forces
   subroutine calcforce_bigdft( posa, forca, boxl, energy, evalf_number, conv )

      implicit none

      !Arguments

      real(kind=8), intent(in),  dimension(3*at%nat), target :: posa
      real(kind=8), intent(out), dimension(3*at%nat), target :: forca
      real(kind=8), dimension(3), intent(inout)           :: boxl
      real(kind=8), intent(out)                           :: energy
      integer,      intent(inout)                         :: evalf_number
      logical,      intent(in)                            :: conv

      !Local variables
      integer  :: infocode, i, ierror 
      real(gp) :: fnoise
      real(gp), allocatable :: xcart(:,:), fcart(:,:)
      !_______________________

      if ( conv ) then                    ! Convergence criterion for the wavefunction optimization
         in%gnrm_cv = gnrm_h              ! in Lanczos procedure.              
      else 
         in%gnrm_cv = gnrm_l                                    
      end if
      ! We transfer acell into 'at'
      at%alat1 = boxl(1)/bohr2ang
      at%alat2 = boxl(2)/bohr2ang
      at%alat3 = boxl(3)/bohr2ang
      ! Need to transform posa into xcart
      ! 1D -> 2D array
      allocate(xcart(3, at%nat))
      do i = 1, at%nat, 1
         xcart(:, i) = (/ posa(i), posa(at%nat + i), posa(2 * at%nat + i) /) / bohr2ang
      end do

      allocate(fcart(3, at%nat))

      if ( first_time ) then              ! This is done by default at the beginning.


         in%inputPsiId = 0
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
         call call_bigdft( nproc, me, at, xcart, in, energy, fcart, fnoise, rst, infocode )
         evalf_number = evalf_number + 1

         in%inputPsiId = 1
         initialised   = .true.
         first_time    = .False.
         new_wf        = .False.

      else 

         if ( .not. initialised ) then
            write(0,*) "No previous call to bigdft_init(). On strike, refuse to work."
            write(*,*) "No previous call to bigdft_init(). On strike, refuse to work."
            stop
         end if

         if ( new_wf ) then               ! if true,  we do not use the previously 
            ! calculated wave function.
            in%inputPsiId = 0
         else 
            in%inputPsiId = 1
         end if

         ! Get into BigDFT
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
         call call_bigdft( nproc, me, at, xcart, in, energy, fcart, fnoise, rst, infocode )
         evalf_number = evalf_number + 1

      end if
      ! Energy in eV 
      energy = energy * ht2ev
      ! box in ang
      boxl(1) = at%alat1 * bohr2ang
      boxl(2) = at%alat2 * bohr2ang
      boxl(3) = at%alat3 * bohr2ang

      ! zero forces for blocked atoms:
      ! This was already done in clean_forces (forces.f90).
      ! But, up to now, ART only works with totally frozen atoms
      ! ( i.e "f" ). Therefore, this is a safe action.
      do i = 1, at%nat, 1
         if ( at%ifrztyp(i) /= 0  .or. in_system(i) /= 0 ) fcart(:,i) = 0.0d0 
      end do 

      call center_f( fcart, at%nat )         ! We remove the net force over our free atomos.

      do i = 1, at%nat, 1                    ! Forces into ev/ang and in 1D array.
         forca( i )              = fcart(1, i) * ht2ev / bohr2ang
         forca( at%nat + i )     = fcart(2, i) * ht2ev / bohr2ang
         forca( 2 * at%nat + i ) = fcart(3, i) * ht2ev / bohr2ang
      end do

      deallocate(xcart)
      deallocate(fcart)

   END SUBROUTINE calcforce_bigdft
   !!***


   !!****f* bigdft_forces/mingeo
   !! FUNCTION
   !!   Minimise geometry
   !! SOURCE
   !!
   subroutine mingeo( posa, forca, boxl, evalf_number, total_energy, success )

      implicit none

      !Arguments
      real(kind=8), intent(inout), dimension(3*atoms_all%nat) :: posa
      real(kind=8), intent(in),    dimension(3*atoms_all%nat), target :: forca
      real(kind=8), intent(inout), dimension(3)     :: boxl
      integer,      intent(inout)                   :: evalf_number
      real(kind=8), intent(out)                     :: total_energy
      logical,      intent(out)                     :: success

      !Local variables
      integer :: i, ierror, ncount_bigdft
      real(gp), allocatable :: xcart(:,:), fcart(:,:)

      if ( .not. initialised ) then
         write(0,*) "No previous call to bigdft_init(). On strike, refuse to work."
         write(*,*) "No previous call to bigdft_init(). On strike, refuse to work."
         stop
      end if

      success = .True.                    ! success will be .False. if:
      !ncount_bigdft > in%ncount_cluster_x-1

      in%gnrm_cv = gnrm_l                 ! For relaxation, we use always the default value in input.dft

      at%alat1 = boxl(1)/bohr2ang
      at%alat2 = boxl(2)/bohr2ang
      at%alat3 = boxl(3)/bohr2ang
      ! Need to transform posa into xcart
      ! 1D -> 2D array
      allocate(xcart(3, at%nat))
      do i = 1, at%nat, 1
         xcart(:, i) = (/ posa(i), posa(at%nat + i), posa(2 * at%nat + i) /) / bohr2ang
      end do

      allocate(fcart(3, at%nat))
      do i = 1, at%nat, 1
         fcart(:, i) = (/ forca(i), forca(at%nat + i), forca(2 * at%nat + i) /) * bohr2ang / ht2ev
      end do

      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      call geopt( nproc, me, xcart, at, fcart, total_energy, rst, in, ncount_bigdft )
      evalf_number = evalf_number + ncount_bigdft 
      if (ncount_bigdft > in%ncount_cluster_x-1) success = .False.

      total_energy = total_energy * ht2ev
      ! box in ang
      boxl(1) = at%alat1 * bohr2ang
      boxl(2) = at%alat2 * bohr2ang
      boxl(3) = at%alat3 * bohr2ang
      ! Positions into ang.
      do i = 1, at%nat, 1
         posa(i)              = xcart(1, i) * bohr2ang
         posa(at%nat + i)     = xcart(2, i) * bohr2ang
         posa(2 * at%nat + i) = xcart(3, i) * bohr2ang
      end do

      deallocate(xcart)
      deallocate(fcart)

   END SUBROUTINE mingeo
   !!***


   !!****f* bigdft_forces/bigdft_finalise
   !! FUNCTION
   !!   Routine to finalise all BigDFT stuff
   !! SOURCE
   !!
   subroutine bigdft_finalise ( )

      implicit none

      !Local variable
      character(len=*), parameter :: subname='bigdft_finalise'
      ! Warning: for what this ??
      call free_restart_objects ( rst, subname )
      ! Warning, there is this note of Damian :
      ! To be completed
      ! but what ??? 

      call memocc( 0, 0, 'count', 'stop' )  ! finalize memory counting.

   END SUBROUTINE bigdft_finalise
   !!***

   !!****f* bigdft_forces/center_f
   !! FUNCTION
   !!   Removes the net force taking into account the blocked atoms
   !!
   !! SOURCE
   !!
   subroutine center_f( vector, natoms )

      implicit none

      !Arguments
      integer, intent(in) :: natoms 
      real(kind=8), dimension(3,natoms), intent(inout), target :: vector

      !Local variables
      integer :: i
      integer :: natoms_f                              ! degrees of freedom
      real(kind=8) :: xtotal, ytotal, ztotal
      logical, dimension(natoms) :: mask

      ! degrees of freedom 
      mask = at%ifrztyp .eq. 0 .and. in_system .eq. 0
      natoms_f = count(mask)

      xtotal = 0.0d0
      ytotal = 0.0d0
      ztotal = 0.0d0

      ! Do over free atoms ( although, the frozen ones add zero ) 
      do i = 1, natoms
         if ( mask(i) ) then
            xtotal = xtotal + vector(1,i)
            ytotal = ytotal + vector(2,i)
            ztotal = ztotal + vector(3,i)
         end if 
      enddo 

      if (iproc==0) then 
         write(*,'(a,1x,1p,e24.17,0p)') 'CENTER: 1) net force over free region ', sqrt(xtotal**2 + ytotal**2 + ztotal**2)
      end if 

      ! The average is only over the degrees of freedom in each direction 
      xtotal = xtotal / natoms_f 
      ytotal = ytotal / natoms_f
      ztotal = ztotal / natoms_f

      if (iproc==0) then 
         write(*,'(a,1x,1p,e24.17,0p)') 'CENTER: 1) residual force along x(pa)=', xtotal  
         write(*,'(a,1x,1p,e24.17,0p)') 'CENTER: 1) residual force along y(pa)=', ytotal  
         write(*,'(a,1x,1p,e24.17,0p)') 'CENTER: 1) residual force along z(pa)=', ztotal  
      end if

      ! Do over free atoms 
      do i = 1, natoms, 1
         if ( mask(i) ) then
            vector(1,i) = vector(1,i) - xtotal
            vector(2,i) = vector(2,i) - ytotal
            vector(3,i) = vector(3,i) - ztotal
         end if 
      end do 

   END SUBROUTINE center_f
   !!***


   subroutine copy_atoms_object(atoms1,atoms2,rxyz,nat,total_nb_atoms,posquant)
      use module_types
      use module_base
      implicit none


      type(atoms_data),intent(in) :: atoms1 
      type(atoms_data),intent(out) :: atoms2
      real(gp), dimension(:,:), pointer     :: rxyz
      integer, intent(in)         :: nat
      integer,intent(in) :: total_nb_atoms
      real(8), dimension(3*total_nb_atoms),intent(in) :: posquant
      integer :: i_stat

      character(len=*), parameter :: subname='copy_atoms_object'



      atoms2%units   = atoms1%units
      atoms2%nat  =    nat
      atoms2%alat1  = atoms1%alat1*bohr2ang
      atoms2%alat2  = atoms1%alat2*bohr2ang
      atoms2%alat3  = atoms1%alat3*bohr2ang

      atoms2%geocode = atoms1%geocode


      allocate(atoms2%natpol(atoms2%nat+ndebug),stat=i_stat)
      call memocc(i_stat,atoms2%natpol,'atoms%natpol',subname)

      !also the spin polarisation and the charge are is fixed to zero by default
      !this corresponds to the value of 100
      !RULE natpol=charge*1000 + 100 + spinpol

      atoms2%natpol(:)=100
      atoms2%natpol(:) = atoms1%natpol(1:nat)

      allocate(atoms2%ifrztyp(atoms2%nat+ndebug),stat=i_stat)
      call memocc(i_stat,atoms2%ifrztyp,'atoms%ifrztyp',subname)


      !this array is useful for frozen atoms
      !no atom is frozen by default
      atoms2%ifrztyp(:)=0
      atoms2%ifrztyp(:)= atoms1%ifrztyp(1:nat)

      allocate(rxyz(3,atoms2%nat+ndebug),stat=i_stat)
      call memocc(i_stat,rxyz,'rxyz',subname)


      rxyz(1,:) = posquant(1:nat)
      rxyz(2,:) = posquant(1+total_nb_atoms:nat+total_nb_atoms)
      rxyz(3,:) = posquant(1+total_nb_atoms+total_nb_atoms:nat+total_nb_atoms+total_nb_atoms)

      atoms2%ntypes  = atoms1%ntypes




      allocate(atoms2%iatype(atoms2%nat+ndebug),stat=i_stat)
      call memocc(i_stat,atoms2%iatype,'atoms%iatype',subname)
      atoms2%iatype(:) = atoms1%iatype(1:nat)

      allocate(atoms2%atomnames(atoms2%ntypes+ndebug),stat=i_stat)
      call memocc(i_stat,atoms2%atomnames,'atoms%atomnames',subname)
      atoms2%atomnames(1:atoms2%ntypes)=atoms1%atomnames(1:atoms2%ntypes)
      atoms2%format = atoms1%format

   END SUBROUTINE copy_atoms_object




   !Right now, this routine only prepares a region made mostly of Si
   !It only passivates the bottom
   !Modifications could be made to passivate the sides as well
   subroutine prepare_quantum_atoms_Si(atoms,posquant,nat)
      use module_types
      use module_base
      use defs
      implicit none

      type(atoms_data),intent(inout) :: atoms
      real(8), dimension(3*natoms),intent(out) :: posquant
      integer, intent(out)           :: nat

      real(8) :: x_min,x_max,y_min,y_max,z_min,z_max
      integer :: i,j,k
      logical, dimension(natoms) :: is_at_quantum
      real(8) :: xij,yij,zij,rij2
      logical :: have_hydro
      character(len=20), dimension(100) :: atomnames
      integer :: i_stat
      integer :: hydro_atom_type

      integer, dimension(natoms) :: numnei
      integer, dimension(natoms,maxnei) :: nei
      real(8), dimension(3) :: invbox

      character(len=*), parameter :: subname='prepare_quantum_atoms_Si'

      invbox = 1.0d0/box

      is_at_quantum = .false. !vectorial operation
      nat = 0
      posquant = 0.0d0

      !First : define a box which encompasses all of the atoms of the quantum region
      x_max = -10000000000000000.0d0
      y_max = -10000000000000000.0d0
      z_max = -10000000000000000.0d0

      x_min = 10000000000000000.0d0
      y_min = 10000000000000000.0d0
      z_min = 10000000000000000.0d0

      have_hydro = .false.
      if (passivate) then
         do i = 1,5   !!!type_name only goes to 5 and is hard coded!!!
            if (type_name(i) == "H") then
               have_hydro = .true.
               hydro_atom_type = i
            endif      
         end do
         if (.not. have_hydro) then
            atomnames(1:atoms%ntypes) = atoms%atomnames(1:atoms%ntypes)
            atoms%ntypes = atoms%ntypes +1
            deallocate(atoms%atomnames,stat = i_stat)
            allocate(atoms%atomnames(atoms%ntypes),stat = i_stat)
            call memocc(i_stat,atoms%atomnames,'atoms%atomnames',subname)
            atomnames(atoms%ntypes) = "H"
            atoms%atomnames(1:atoms%ntypes) = atomnames(1:atoms%ntypes)
            hydro_atom_type = atoms%ntypes
         endif
      endif


      call neighbours(natoms,pos,box,boundary,maxnei,numnei, nei)

      do i = 1,nbr_quantum
         is_at_quantum(i) = .true.
         nat = nat + 1
         if ( pos(i) < x_min) x_min = pos(i) -0.02d0
         if ( pos(i+natoms) < y_min) y_min = pos(i+natoms) -0.02d0
         if ( pos(i+natoms+natoms) < z_min) z_min = pos(i+natoms+natoms) -0.02d0

         if ( pos(i) > x_max) x_max = pos(i) +0.02d0
         if ( pos(i+natoms) > y_max) y_max = pos(i+natoms)  +0.02d0
         if ( pos(i+natoms+natoms) > z_max) z_max = pos(i+natoms+natoms) +0.02d0

         posquant(i) = pos(i)
         posquant(i+natoms) = pos(i+natoms)
         posquant(i+natoms+natoms) = pos(i+natoms+natoms)
      enddo

      do i = 1,nbr_quantum
         if (passivate) then 
            do j = 1,numnei(i)
               k = nei(i,j)
               if ( .not. is_at_quantum(k)) then 
                  xij = pos(k)-pos(i) - box(1) * nint((pos(k)-pos(i))*invbox(1))
                  yij = pos(k+natoms)-pos(i+natoms)
                  zij = pos(k+2*natoms)-pos(i+2*natoms) - box(3) * nint((pos(k+2*natoms)-pos(i+2*natoms))*invbox(3))
                  rij2 = xij*xij + yij*yij + zij*zij
                  if (rij2 .lt. 2.5d0*2.5d0) then
                     nat = nat + 1
                     posquant(nat) = pos(i) + 0.5d0*xij
                     posquant(nat+natoms) = pos(i+natoms) + 0.5d0*yij
                     posquant(nat+natoms+natoms) = pos(i+2*natoms) + 0.5d0*zij !we passivate with hydrogene at this distance
                     atoms%iatype(nat) = hydro_atom_type
                     atoms%ifrztyp(i) = 1  !this atom is frozen
                     atoms%ifrztyp(nat) = 1  !this one as well
                  endif
               endif   
            enddo
         endif
      enddo

   END SUBROUTINE prepare_quantum_atoms_Si


   !!****f* bigdft_forces/check_force_clean_wf
   !! FUNCTION
   !! SOURCE
   !!
   subroutine check_force_clean_wf( posa, boxl, evalf_number, total_energy, success )

      implicit none

      real(kind=8), intent(inout), dimension(3*atoms_all%nat) :: posa
      real(kind=8), intent(inout), dimension(3)     :: boxl
      integer,      intent(inout)                   :: evalf_number
      real(kind=8), intent(out)                     :: total_energy
      logical,      intent(out)                     :: success

      !Local variables
      integer      :: infocode, i, ierror, ncount_bigdft 
      real(kind=8) :: energy
      real(gp)     :: fnoise
      real(gp)     ::  fmax, fnrm
      real(gp), allocatable :: xcart(:,:), fcart(:,:)
      !_______________________

      in%inputPsiId = 0 
      in%gnrm_cv = gnrm_l 
      ! We transfer acell into 'at'
      at%alat1 = boxl(1)/bohr2ang
      at%alat2 = boxl(2)/bohr2ang
      at%alat3 = boxl(3)/bohr2ang

      allocate(xcart(3, at%nat))
      do i = 1, at%nat, 1
         xcart(:, i) = (/ posa(i), posa(at%nat + i), posa(2 * at%nat + i) /) / bohr2ang
      end do

      allocate(fcart(3, at%nat))

      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      call call_bigdft( nproc, me, at, xcart, in, energy, fcart, fnoise, rst, infocode )
      evalf_number = evalf_number + 1
      in%inputPsiId = 1

      call fnrmandforcemax(fcart,fnrm,fmax, at%nat)

      if ( fmax > in%forcemax ) then

         if ( iproc == 0 ) then
            write(*,*) 'BART:check_force_clean_wf'
            write(*,*) 'BART: fmax =', fmax,'H/Bohr. We relax again!' 
         end if

         call MPI_Barrier(MPI_COMM_WORLD,ierror)
         call geopt( nproc, me, xcart, at, fcart, total_energy, rst, in, ncount_bigdft )
         evalf_number = evalf_number + ncount_bigdft 
         if (ncount_bigdft > in%ncount_cluster_x-1) success = .False.

         ! and we clean again here
         in%inputPsiId = 0 
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
         call call_bigdft( nproc, me, at, xcart, in, energy, fcart, fnoise, rst, infocode )
         evalf_number = evalf_number + 1
         in%inputPsiId = 1

         total_energy = total_energy * ht2ev
         ! box in ang
         boxl(1) = at%alat1 * bohr2ang
         boxl(2) = at%alat2 * bohr2ang
         boxl(3) = at%alat3 * bohr2ang
         ! Positions into ang.
         do i = 1, at%nat, 1
            posa(i)              = xcart(1, i) * bohr2ang
            posa(at%nat + i)     = xcart(2, i) * bohr2ang
            posa(2 * at%nat + i) = xcart(3, i) * bohr2ang
         end do

      end if 

      deallocate(xcart)
      deallocate(fcart)

   END SUBROUTINE check_force_clean_wf


END MODULE bigdft_forces
!!***
