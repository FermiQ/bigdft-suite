!> @file
!! LOCalized REGion initialization for support functions
!! @author
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Initialization of the localized regions for the support functions
module locregs_init
  use locregs
  implicit none

  private

  !> Public routines
  public :: initLocregs
  public :: check_linear_inputguess
  public :: small_to_large_locreg
  public :: assign_to_atoms_and_locregs
  public :: lr_set


  contains

    !> form a locreg given atoms. Should be used to create the Glr localisation region
    subroutine lr_set(lr,iproc,OCLconv,dump,crmult,frmult,hgrids,rxyz,atoms,calculate_bounds,output_grid)
      use module_interfaces
      use locregs
      use module_atoms
      use module_defs, only: gp
      implicit none
      type(locreg_descriptors), intent(out) :: lr
      integer, intent(in) :: iproc
      logical, intent(in) :: OCLconv,dump
      type(atoms_data), intent(inout) :: atoms
      real(gp), intent(in) :: crmult,frmult
      real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
      real(gp), dimension(3), intent(inout) :: hgrids
      logical, intent(in) :: calculate_bounds,output_grid

      lr=locreg_null()
      call system_size(atoms,rxyz,crmult,frmult,&
           hgrids(1),hgrids(2),hgrids(3),OCLconv,lr)
      if (iproc == 0 .and. dump) call print_atoms_and_grid(lr, atoms, rxyz, &
           hgrids(1),hgrids(2),hgrids(3))
      call createWavefunctionsDescriptors(iproc,hgrids(1),hgrids(2),hgrids(3),atoms,&
           rxyz,crmult,frmult,calculate_bounds,lr, output_grid)
      if (iproc == 0 .and. dump) call print_wfd(lr%wfd)
    end subroutine lr_set


    ! lzd%llr already allocated, locregcenter and locrad already filled - could tidy this!
    subroutine initLocregs(iproc, nproc, lzd, hx, hy, hz, rxyz,locrad, orbs, Glr, locregShape, lborbs)
      use module_base
      use module_types
      implicit none

      ! Calling arguments
      integer, intent(in) :: iproc, nproc
      type(local_zone_descriptors), intent(inout) :: lzd
      real(kind=8), intent(in) :: hx, hy, hz
      real(gp), dimension(lzd%nlr), intent(in) :: locrad
      real(gp), dimension(3,lzd%nlr), intent(in) :: rxyz
      type(orbitals_data), intent(in) :: orbs
      type(locreg_descriptors), intent(in) :: Glr
      character(len=1), intent(in) :: locregShape
      type(orbitals_data),optional,intent(in) :: lborbs

      ! Local variables
      integer :: jorb, jjorb, jlr, ilr
      character(len=*), parameter :: subname='initLocregs'
      logical,dimension(:), allocatable :: calculateBounds
      real(gp), dimension(3) :: hgrids
      real(gp) :: maxlr, lrtol

      call f_routine(id=subname)


      calculateBounds = f_malloc(lzd%nlr,id='calculateBounds')
      calculateBounds=.false.

      do jorb=1,orbs%norbp
         jjorb=orbs%isorb+jorb
         jlr=orbs%inWhichLocreg(jjorb)
         calculateBounds(jlr)=.true.
      end do

      if(present(lborbs)) then
         do jorb=1,lborbs%norbp
            jjorb=lborbs%isorb+jorb
            jlr=lborbs%inWhichLocreg(jjorb)
            calculateBounds(jlr)=.true.
         end do
      end if

      ! make sure we have one locreg which is defined on all MPI so that we can use it for onsite overlap
      ! need to make sure that there are no other bigger locrads
      ! it would be helpful to make a variable indicating which locreg we have, but don't want to edit structures unnecessarily...
      ! just in case there is some noise in the locrads
      lrtol=1.0e-3
      jlr=1
      maxlr=Lzd%Llr(jlr)%locrad
      do ilr=2,Lzd%nlr
         if (Lzd%Llr(ilr)%locrad > maxlr + lrtol) then
            jlr=ilr
            maxlr=Lzd%Llr(jlr)%locrad
         end if
      end do
      Lzd%llr_on_all_mpi=jlr

      hgrids=[hx,hy,hz]
      if(locregShape=='c') then
         calculateBounds=.true.
         call determine_locreg_parallel(iproc,nproc,lzd%nlr,rxyz,locrad,&
              hx,hy,hz,Glr,lzd%Llr,orbs,calculateBounds)
         !stop 'locregShape c is deprecated'
      else if(locregShape=='s') then
         call determine_locregSphere_parallel(iproc, nproc, lzd%nlr, hx, hy, hz, &
              orbs, Glr, lzd%Llr, calculateBounds, Lzd%llr_on_all_mpi)
      end if
      call f_free(calculateBounds)

      !DEBUG
      !do ilr=1,lin%nlr
      !    if(iproc==0) write(*,'(1x,a,i0)') '>>>>>>> zone ', ilr
      !    if(iproc==0) write(*,'(3x,a,4i10)') 'nseg_c, nseg_f, nvctr_c, nvctr_f', lin%Llr(ilr)%wfd%nseg_c, lin%Llr(ilr)%wfd%nseg_f, lin%Llr(ilr)%wfd%nvctr_c, lin%Llr(ilr)%wfd%nvctr_f
      !    if(iproc==0) write(*,'(3x,a,3i8)') 'lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i', lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i
      !    if(iproc==0) write(*,'(a,6i8)') 'lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3',&
      !    lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3
      !end do
      !END DEBUG

      lzd%linear=.true.

      call f_release_routine()

    end subroutine initLocregs

    subroutine small_to_large_locreg(iproc, norb, norbp, isorb, in_which_locreg, &
               npsidim_orbs_small, npsidim_orbs_large, lzdsmall, lzdlarge, &
               phismall, philarge, to_global)
      use module_base
      use module_types, only: orbitals_data, local_zone_descriptors
      use locreg_operations, only: lpsi_to_global2
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, norb, norbp, isorb, npsidim_orbs_small, npsidim_orbs_large
      integer,dimension(norb),intent(in) :: in_which_locreg
      type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
      real(kind=8),dimension(npsidim_orbs_small),intent(in) :: phismall
      real(kind=8),dimension(npsidim_orbs_large),intent(out) :: philarge
      logical,intent(in),optional :: to_global

      ! Local variables
      integer :: ists, istl, iorb, ilr, sdim, ldim, nspin
      logical :: global

      call f_routine(id='small_to_large_locreg')

      if (present(to_global)) then
         global=to_global
      else
         global=.false.
      end if

      call timing(iproc,'small2large','ON') ! lr408t 
      ! No need to put arrays to zero, Lpsi_to_global2 will handle this.
      call f_zero(philarge)
      ists=1
      istl=1
      do iorb=1,norbp
         ilr = in_which_locreg(isorb+iorb)
         sdim=lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
         if (global) then
            ldim=lzdsmall%glr%wfd%nvctr_c+7*lzdsmall%glr%wfd%nvctr_f
         else
            ldim=lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
         end if
         nspin=1 !this must be modified later
         if (global) then
            call Lpsi_to_global2(iproc, sdim, ldim, norb, nspin, lzdsmall%glr, &
                 lzdsmall%llr(ilr), phismall(ists), philarge(istl))
         else
            call Lpsi_to_global2(iproc, sdim, ldim, norb, nspin, lzdlarge%llr(ilr), &
                 lzdsmall%llr(ilr), phismall(ists), philarge(istl))
         end if
         ists=ists+sdim
         istl=istl+ldim
      end do
      if(norbp>0 .and. ists/=npsidim_orbs_small+1) then
         write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',ists,'=ists /= npsidim_orbs_small+1=',npsidim_orbs_small+1
         stop
      end if
      if(norbp>0 .and. istl/=npsidim_orbs_large+1) then
         write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',istl,'=istl /= npsidim_orbs_large+1=',npsidim_orbs_large+1
         stop
      end if
      call timing(iproc,'small2large','OF') ! lr408t 
      call f_release_routine()
    end subroutine small_to_large_locreg


    subroutine determine_locregSphere_parallel(iproc,nproc,nlr,hx,hy,hz,orbs,Glr,Llr,calculateBounds,llr_on_all_mpi)!,outofzone)
    
      use module_base
      use module_types
      !use module_interfaces, except_this_one => determine_locregSphere_parallel
      use communications, only: communicate_locreg_descriptors_basics, communicate_locreg_descriptors_keys
      use bounds, only: locreg_bounds , ext_buffers
      implicit none
      integer, intent(in) :: iproc,nproc
      integer, intent(in) :: nlr
      real(gp), intent(in) :: hx,hy,hz
      !type(atomic_structure),intent(in) :: astruct
      type(orbitals_data),intent(in) :: orbs
      type(locreg_descriptors), intent(in) :: Glr
      type(locreg_descriptors), dimension(nlr), intent(inout) :: Llr
      logical,dimension(nlr),intent(in) :: calculateBounds
      integer, intent(in) :: llr_on_all_mpi
    !  integer, dimension(3,nlr),intent(out) :: outofzone
      !local variables
      character(len=*), parameter :: subname='determine_locreg'
      logical :: Gperx,Gpery,Gperz!>,Lperx,Lpery,Lperz
      logical :: warningx,warningy,warningz
!!$      integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
!!$      integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
      integer :: ilr!,isx,isy,isz,iex,iey,iez
!!$      integer :: ln1,ln2,ln3
      integer :: ii, iorb, jproc, iiorb
      ! integer :: iat, norb, norbu, norbd, nspin
!!$      integer, dimension(3) :: outofzone
      !! start and end points for each direction
      integer, dimension(2,3) :: nbox 
      integer, dimension(:),allocatable :: rootarr, norbsperatom, norbsperlocreg, onwhichmpi
      !real(8),dimension(:,:),allocatable :: locregCenter
      type(orbitals_data) :: orbsder
      logical :: perx, pery, perz
      real(gp), dimension(3) :: hgrids

      call f_routine(id='determine_locregSphere_parallel')
    

      rootarr = f_malloc(nlr,id='rootarr')
    
      ! Determine how many locregs one process handles at most
      ii=ceiling(dble(nlr)/dble(nproc))
      !determine the limits of the different localisation regions
      rootarr=1000000000
    
      onwhichmpi=f_malloc(nlr,id='onwhichmpi')
      iiorb=0
      do jproc=0,nproc-1
          do iorb=1,orbs%norb_par(jproc,0)
              iiorb=iiorb+1
              onWhichMPI(iiorb)=jproc
          end do
      end do
    
      ! Periodicity in the three directions
      Gperx=(Glr%geocode /= 'F')
      Gpery=(Glr%geocode == 'P')
      Gperz=(Glr%geocode /= 'F')
    
      call timing(iproc,'wfd_creation  ','ON')  
      do ilr=1,nlr
         !initialize out of zone and logicals
         !outofzone (:) = 0     
         warningx = .false.
         warningy = .false.
         warningz = .false. 
!!$         xperiodic = .false.
!!$         yperiodic = .false.
!!$         zperiodic = .false. 
    
    
         if(calculateBounds(ilr) .or. ilr==llr_on_all_mpi) then 
             ! This makes sure that each locreg is only handled once by one specific processor.
        
             ! Determine the extrema of this localization regions (using only the coarse part, since this is always larger or equal than the fine part).
             call determine_boxbounds_sphere(gperx, gpery, gperz, glr%d%n1, glr%d%n2, glr%d%n3, glr%ns1, glr%ns2, glr%ns3, &
                  hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
                  glr%wfd%nseg_c, glr%wfd%keygloc, glr%wfd%keyvloc, &
                  nbox(1,1),nbox(1,2),nbox(1,3),nbox(2,1),nbox(2,2),nbox(2,3))
             !!!>isx, isy, isz, iex, iey, iez)
             !write(*,'(a,3i7)') 'ilr, isx, iex', ilr, isx, iex
             hgrids = [hx,hy,hz]
             call lr_box(llr(ilr),Glr,hgrids,nbox,.false.)
            ! construct the wavefunction descriptors (wfd)
            if (calculateBounds(ilr)) rootarr(ilr)=iproc
            call determine_wfdSphere(ilr,nlr,Glr,hx,hy,hz,Llr)
         end if
      end do !on ilr
      call timing(iproc,'wfd_creation  ','OF') 
    
      ! Communicate the locregs
      ! This communication is uneffective. Instead of using bcast we should be using mpialltoallv.
      call timing(iproc,'comm_llr      ','ON')
      if (nproc > 1) then
         call fmpi_allreduce(rootarr(1), nlr, FMPI_MIN, comm=bigdft_mpi%mpi_comm)
         
         ! Communicate those parts of the locregs that all processes need.
         call communicate_locreg_descriptors_basics(iproc, nproc, nlr, rootarr, orbs, llr)
    
         ! SM: Seems not to be necessary to be called again
         !!do ilr=1,nlr
         !!   if(.not. calculateBounds(ilr)) call lr_box(llr(ilr),Glr,[hx,hy,hz])
         !!   !    write(*,*) 'iproc, nseg_c', iproc, llr(ilr)%wfd%nseg_c
         !!end do
    
         ! Now communicate those parts of the locreg that only some processes need (the keys).
         ! For this we first need to create orbsder that describes the derivatives.
         !call create_orbsder()
    
         !iiorb=0
         !onwhichmpider=f_malloc(orbsder%norb,id='onwhichmpider')
         !do jproc=0,nproc-1
         !   do iorb=1,orbsder%norb_par(jproc,0)
         !     iiorb=iiorb+1
         !     onWhichMPIder(iiorb)=jproc
         !   end do
         !end do
    
         ! Now communicate the keys
         call communicate_locreg_descriptors_keys(iproc, nproc, nlr, glr, llr, orbs, rootarr, onwhichmpi, llr_on_all_mpi)
    
         !call deallocate_orbitals_data(orbsder)
         !call f_free(onwhichmpider)
      end if
      call timing(iproc,'comm_llr      ','OF')
    
      !create the bound arrays for the locregs we need on the MPI tasks
      call timing(iproc,'calc_bounds   ','ON') 
      do ilr=1,nlr
             if (Llr(ilr)%geocode=='F' .and. (calculateBounds(ilr) .or. ilr==llr_on_all_mpi) ) then
                !write(*,*) 'calling locreg_bounds, ilr', ilr
                call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
                     Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
                     Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
             end if
      end do
    
      call timing(iproc,'calc_bounds   ','OF') 
    
      call f_free(rootarr)
    
      call f_free(onwhichmpi)
      call f_release_routine()
    END SUBROUTINE determine_locregSphere_parallel

    !> Determine a set of localisation regions from the centers and the radii.
    !! cut in cubes the global reference system
    subroutine determine_locreg_parallel(iproc,nproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,orbs,calculateBounds)!,outofzone)
      use module_base
      use module_types
      use bounds, only: locreg_bounds, ext_buffers
      use box
      implicit none
      integer, intent(in) :: iproc,nproc
      integer, intent(in) :: nlr
      real(gp), intent(in) :: hx,hy,hz
      type(locreg_descriptors), intent(in) :: Glr
      real(gp), dimension(nlr), intent(in) :: locrad
      real(gp), dimension(3,nlr), intent(in) :: cxyz
      type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
      type(orbitals_data),intent(in) :: orbs
      logical,dimension(nlr),intent(in) :: calculateBounds
    !  integer, dimension(3,nlr),intent(out) :: outofzone
      !local variables
      logical :: Gperx,Gpery,Gperz,Lperx,Lpery,Lperz
      logical :: warningx,warningy,warningz,calc
      integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
      integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
      integer :: ilr,isx,isy,isz,iex,iey,iez,iorb
      integer :: ln1,ln2,ln3
      integer,dimension(3) :: outofzone
      integer, dimension(2,3) :: nbox 
      real(gp) :: rx,ry,rz,cutoff
    !!  integer :: iilr,ierr
    !!  integer,dimension(0:nproc-1) :: nlr_par,islr_par
    
      !!if (iproc == 0) then
      !!   write(*,*)'Inside determine_locreg_periodic:'
      !!end if
    
    !  call parallel_repartition_locreg(iproc,nproc,nlr,nlr_par,islr_par)
    
      !initialize out of zone and logicals
!!$      outofzone (:) = 0
      warningx = .false.
      warningy = .false.
      warningz = .false.
      !determine the limits of the different localisation regions
      do ilr=1,nlr
         call nullify_locreg_descriptors(Llr(ilr))
         calc=.false.
         do iorb=1,orbs%norbp
            if(ilr == orbs%inwhichLocreg(iorb+orbs%isorb)) calc=.true.
         end do
         if (.not. calc) cycle         !calculate only for the locreg on this processor, without repeating for same locreg
    
         rx=cxyz(1,ilr)
         ry=cxyz(2,ilr)
         rz=cxyz(3,ilr)
    
         cutoff=locrad(ilr)

         nbox=box_nbox_from_cutoff(Glr%mesh_coarse,cxyz(1,ilr),locrad(ilr),inner=.false.)

!!$         nbox(1,1)=floor((rx-cutoff)/hx)
!!$         nbox(1,2)=floor((ry-cutoff)/hy)
!!$         nbox(1,3)=floor((rz-cutoff)/hz)
!!$    
!!$         nbox(2,1)=ceiling((rx+cutoff)/hx)
!!$         nbox(2,2)=ceiling((ry+cutoff)/hy)
!!$         nbox(2,3)=ceiling((rz+cutoff)/hz)

         call lr_box(Llr(ilr),Glr,[hx,hy,hz],nbox,.true.)
        
        ! construct the wavefunction descriptors (wfd)
         call determine_wfd_periodicity(ilr,nlr,Glr,Llr)
    
         ! Sould check if nfu works properly... also relative to locreg!!
         !if the localisation region is isolated build also the bounds
         if (Llr(ilr)%geocode=='F') then
            ! Check whether the bounds shall be calculated. Do this only if the currect process handles
            ! orbitals in the current localization region.
            if(calculateBounds(ilr)) then
               call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
                   Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
                   Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
            end if
         end if
      end do !on iilr
    
    !  call make_LLr_MpiType(Llr,nlr,mpiLlr)
    
    !  call MPI_ALLREDUCE(Llr(1),Llr(1),nlr,mpidtypg,FMPI_SUM,bigdft_mpi%mpi_comm,ierr)
      !after all localisation regions are determined draw them
      !call draw_locregs(nlr,hx,hy,hz,Llr)
    
    END SUBROUTINE determine_locreg_parallel

    subroutine assign_to_atoms_and_locregs(iproc, nproc, norb, nat, nspin, norb_per_atom, rxyz, &
         on_which_atom, in_which_locreg)
      use module_base
      use sort, only: QsortC
      implicit none

      integer,intent(in):: nat,iproc,nproc,nspin,norb
      integer,dimension(nat),intent(in):: norb_per_atom
      real(8),dimension(3,nat),intent(in):: rxyz
      integer,dimension(:),pointer,intent(out):: on_which_atom, in_which_locreg

      ! Local variables
      integer:: iat, jproc, iiOrb, iorb, jorb, jat, iiat, i_stat, i_all, ispin, iispin, istart, iend, ilr, jjat
      integer :: jjorb, norb_check, norb_nospin, isat, ii, iisorb, natp
      integer,dimension(:),allocatable :: irxyz_ordered, nat_par
      real(kind=8), parameter :: tol=1.0d-6 
      real(8):: tt, dmin, minvalue, xmin, xmax, ymin, ymax, zmin, zmax
      integer:: iatxmin, iatxmax, iatymin, iatymax, iatzmin, iatzmax, idir, nthread, ithread
      real(8),dimension(3):: diff
      real(kind=8),dimension(:),allocatable :: rxyz_dir, minvalue_thread
      integer,dimension(:),allocatable :: iiat_thread
      integer :: nbin, ibin, iatx, iiat_tot, i2dir, imin, jmin
      real(kind=8) :: dist, d, rmin, rmax, tmin
      integer,dimension(:),allocatable :: nat_per_bin, iat_per_bin
      integer,dimension(:,:),allocatable :: iat_bin
      logical,dimension(:),allocatable :: covered
      real(kind=8),parameter :: binwidth = 6.d0 ! should probably be more or less the same as the locrads... make it variable?
      character(len=3),parameter :: old='old', new='new'
      character(len=3),parameter :: mode = new !old

      call f_routine(id='assign_to_atoms_and_locregs')

      ! Determine in which direction the system has its largest extent
      xmin=1.d100
      ymin=1.d100
      zmin=1.d100
      xmax=-1.d100
      ymax=-1.d100
      zmax=-1.d100
      do iat=1,nat
         if(rxyz(1,iat)<xmin) then
            xmin=rxyz(1,iat)
            iatxmin=iat
         end if
         if(rxyz(1,iat)>xmax) then
            xmax=rxyz(1,iat)
            iatxmax=iat
         end if
         if(rxyz(2,iat)<ymin) then
            ymin=rxyz(2,iat)
            iatymin=iat
         end if
         if(rxyz(2,iat)>ymax) then
            ymax=rxyz(2,iat)
            iatymax=iat
         end if
         if(rxyz(3,iat)<zmin) then
            zmin=rxyz(3,iat)
            iatzmin=iat
         end if
         if(rxyz(3,iat)>zmax) then
            zmax=rxyz(3,iat)
            iatzmax=iat
         end if
      end do

      diff(1)=xmax-xmin
      diff(2)=ymax-ymin
      diff(3)=zmax-zmin
      !First 4 ifs control if directions the same length to disambiguate (was random before)
      !else, just choose the biggest
      if(abs(diff(1)-diff(2)) < tol .and. diff(1) > diff(3)) then
         idir=1
         iiat=iatxmin
         rmin = xmin
         rmax = xmax
         i2dir = 2
      else if(abs(diff(1)-diff(3)) < tol .and. diff(1) > diff(2)) then
         idir=1
         iiat=iatxmin
         rmin = xmin
         rmax = xmax
         i2dir = 3
      else if(abs(diff(2)-diff(3)) < tol .and. diff(2) > diff(1)) then
         idir=2
         iiat=iatymin
         rmin = ymin
         rmax = ymax
         i2dir = 3
      else if(abs(diff(1)-diff(3)) < tol .and. abs(diff(2)-diff(3)) < tol) then
         idir=1
         iiat=iatxmin
         rmin = xmin
         rmax = xmax
         i2dir = 3
      else
         if(maxloc(diff,1)==1) then
            idir=1
            iiat=iatxmin
            rmin = xmin
            rmax = xmax
         else if(maxloc(diff,1)==2) then
            idir=2
            iiat=iatymin
            rmin = ymin
            rmax = ymax
         else if(maxloc(diff,1)==3) then
            idir=3
            iiat=iatzmin
            rmin = zmin
            rmax = zmax
         else
            call f_err_throw('ERROR: not possible to determine the maximal extent')
         end if
         diff(idir) = -1.d100
         i2dir = maxloc(diff,1)
      end if


      irxyz_ordered = f_malloc(nat,id='irxyz_ordered')
      if (mode==old) then
         ! Lookup array which orders the atoms with respect to the direction idir
         rxyz_dir = f_malloc(nat,id='rxyz_dir')
         do ilr=1,nat
            rxyz_dir(ilr) = rxyz(idir,ilr)
            irxyz_ordered(ilr) = ilr
         end do
         !!do ilr=1,nat
         !!    if (iproc==0) write(*,*) 'B: rxyz_dir(ilr), irxyz_ordered(ilr)', rxyz_dir(ilr), irxyz_ordered(ilr)
         !!end do
         call QsortC(rxyz_dir, irxyz_ordered)
         !!do ilr=1,nat
         !!    if (iproc==0) write(*,*) 'A: rxyz_dir(ilr), irxyz_ordered(ilr)', rxyz_dir(ilr), irxyz_ordered(ilr)
         !!end do
         call f_free(rxyz_dir)
      end if


      if (mode==new) then
         ! # NEW ##################################################
         ! Calculate the number of bins (in the largest direction)
         nbin = max(1,ceiling((rmax-rmin)/binwidth))

         ! Assign the atoms to the bins. First count how many atoms per bin we have...
         nat_per_bin = f_malloc0(nbin,id='nat_per_bin')
         do iat=1,nat
            dist = rxyz(idir,iat)-rmin
            ibin = ceiling(dist/binwidth)
            ibin = max(ibin,1) ! to avoid bin 0 for the "minimal" atoms, i.e. the one which defines the minium value
            if (ibin<1 .or. ibin>nbin) then
               call f_err_throw('wrong bin (ibin='//trim(yaml_toa(ibin))//', nbin='//trim(yaml_toa(nbin))//')',&
                    err_name='BIGDFT_RUNTIME_ERROR')
            end if
            nat_per_bin(ibin) = nat_per_bin(ibin) + 1
         end do
         ! ... and now assign the atoms to the bins.
         iat_bin = f_malloc0((/maxval(nat_per_bin),nbin/),id='iat_bin')
         iat_per_bin = f_malloc0(nbin,id='iat_per_bin')
         do iat=1,nat
            dist = rxyz(idir,iat)-rmin
            ibin = ceiling(dist/binwidth)
            ibin = max(ibin,1) ! to avoid bin 0 for the "minimal" atoms, i.e. the one which defines the minium value
            if (ibin<1 .or. ibin>nbin) then
               call f_err_throw('wrong bin (ibin='//trim(yaml_toa(ibin))//', nbin='//trim(yaml_toa(nbin))//')',&
                    err_name='BIGDFT_RUNTIME_ERROR')
            end if
            iat_per_bin(ibin) = iat_per_bin(ibin) + 1
            iat_bin(iat_per_bin(ibin),ibin) = iat
         end do
         ! Check
         do ibin=1,nbin
            if (iat_per_bin(ibin)/=nat_per_bin(ibin)) then
               call f_err_throw('iat_per_bin(ibin)/=nat_per_bin(ibin)',err_name='BIGDFT_RUNTIME_ERROR')
            end if
         end do
         call f_free(iat_per_bin)

         ! Order the atoms by always choosing the closest one within a bin
         iiat_tot = 0
         iatx = iiat !starting atoms
         do ibin=1,nbin
            covered = f_malloc(nat_per_bin(ibin),id='covered')
            covered(:) = .false.
            ! Determine the minimum value in the i2dir direction
            tmin = huge(1.d0)
            do iat=1,nat_per_bin(ibin)
               jjat = iat_bin(iat,ibin)
               if (rxyz(i2dir,jjat)<tmin) then
                  tmin = rxyz(i2dir,jjat)
                  imin = jjat
               end if
            end do
            iatx = imin!jjat
            do iat=1,nat_per_bin(ibin)
               iiat_tot = iiat_tot + 1
               dmin = huge(1.d0)
               do jat=1,nat_per_bin(ibin)
                  jjat = iat_bin(jat,ibin)
                  if (covered(jat)) cycle
                  d = (rxyz(1,jjat)-rxyz(1,iatx))**2 + (rxyz(2,jjat)-rxyz(2,iatx))**2 + (rxyz(3,jjat)-rxyz(3,iatx))**2
                  if (d<dmin) then
                     dmin = d
                     jmin = jat
                  end if
               end do
               covered(jmin) = .true.
               !iatx = jjat !iat_bin(jmin,ibin)
               irxyz_ordered(iiat_tot) = iat_bin(jmin,ibin) !iatx
               !write(*,*) 'iiat_tot, irxyz_ordered(iiat_tot)', iiat_tot, irxyz_ordered(iiat_tot)
            end do
            call f_free(covered)
         end do

         call f_free(nat_per_bin)
         call f_free(iat_bin)
         ! # NEW ##################################################
      end if

      on_which_atom = f_malloc0_ptr(norb,id='on_which_atom')
      in_which_locreg = f_malloc0_ptr(norb,id='inWhichLocreg')

      ! Total number of orbitals without spin
      norb_nospin = 0
      do jat=1,nat
         norb_nospin = norb_nospin + norb_per_atom(jat)
      end do

      ! Parallelization over the atoms
      nat_par = f_malloc(0.to.nproc-1,id='nat_par')
      natp = nat/nproc
      nat_par(:) = natp
      ii = nat-natp*nproc
      nat_par(0:ii-1) = nat_par(0:ii-1) + 1
      isat = sum(nat_par(0:iproc-1))
      ! check
      ii = sum(nat_par)
      if (ii/=nat) call f_err_throw('ii/=nat',err_name='BIGDFT_RUNTIME_ERROR')
      iisorb = 0
      iiat = 0
      do jproc=0,iproc-1
         do iat=1,nat_par(jproc)
            iiat = iiat + 1
            jjat = irxyz_ordered(iiat)
            iisorb = iisorb + norb_per_atom(jjat)
         end do
      end do

      ! Calculate on_which_atom and in_which_locreg...
      norb_check = 0
      do ispin=1,nspin
         iiorb = iisorb + (ispin-1)*norb_nospin
         do iat=1,nat_par(iproc)
            iiat = isat + iat
            jjat = irxyz_ordered(iiat)
            jjorb = (ispin-1)*norb_nospin
            do jat=1,jjat-1
               jjorb = jjorb + norb_per_atom(jat)
            end do
            do iorb=1,norb_per_atom(jjat)
               iiorb = iiorb + 1
               jjorb = jjorb + 1
               on_which_atom(iiorb) = jjat
               in_which_locreg(iiorb) = jjorb
               norb_check = norb_check + 1
            end do
         end do
      end do
      call fmpi_allreduce(norb_check, 1, FMPI_SUM, comm=bigdft_mpi%mpi_comm)
      if (norb_check/=norb) call f_err_throw('norb_check (' // &
           & trim(yaml_toa(norb_check)) // ') /= norb (' // &
           & trim(yaml_toa(norb)) // ')',err_name='BIGDFT_RUNTIME_ERROR')

      call fmpi_allreduce(on_which_atom, FMPI_SUM, comm=bigdft_mpi%mpi_comm)
      call fmpi_allreduce(in_which_locreg, FMPI_SUM, comm=bigdft_mpi%mpi_comm)

      call f_free(nat_par)
      call f_free(irxyz_ordered)


      call f_release_routine()

    end subroutine assign_to_atoms_and_locregs

    ! SM: Don't really know what this is for
    !> Determine a set of localisation regions from the centers and the radii.
    !! cut in cubes the global reference system
    subroutine check_linear_inputguess(iproc,nlr,cxyz,locrad,hx,hy,hz,Glr,linear)
      use module_base
      use module_types
      use box
      implicit none
      integer, intent(in) :: iproc
      integer, intent(in) :: nlr
      logical,intent(out) :: linear
      real(gp), intent(in) :: hx,hy,hz
      type(locreg_descriptors), intent(in) :: Glr
      real(gp), dimension(nlr), intent(in) :: locrad
      real(gp), dimension(3,nlr), intent(in) :: cxyz
      !local variables
      integer, parameter :: START_=1,END_=2
      logical :: warningx,warningy,warningz
      integer :: ilr,isx,isy,isz,iex,iey,iez
      integer :: ln1,ln2,ln3
      integer, dimension(2,3) :: nbox
      real(gp) :: rx,ry,rz,cutoff

      !to check if the floor and the ceiling are meaningful in this context

      linear = .true.

      !determine the limits of the different localisation regions
      do ilr=1,nlr

         !initialize logicals
         warningx = .false.
         warningy = .false.
         warningz = .false.

         rx=cxyz(1,ilr)
         ry=cxyz(2,ilr)
         rz=cxyz(3,ilr)

         cutoff=locrad(ilr)

         nbox=box_nbox_from_cutoff(Glr%mesh_coarse,cxyz(1,ilr),locrad(ilr),inner=.false.)

         isx=nbox(START_,1)
         isy=nbox(START_,2)
         isz=nbox(START_,3)

         iex=nbox(END_,1)
         iey=nbox(END_,2)
         iez=nbox(END_,3)

!!$         isx=floor((rx-cutoff)/hx)
!!$         isy=floor((ry-cutoff)/hy)
!!$         isz=floor((rz-cutoff)/hz)
!!$
!!$         iex=ceiling((rx+cutoff)/hx)
!!$         iey=ceiling((ry+cutoff)/hy)
!!$         iez=ceiling((rz+cutoff)/hz)

         ln1 = iex-isx
         ln2 = iey-isy
         ln3 = iez-isz

         ! First check if localization region fits inside box
         if (iex - isx >= Glr%d%n1 - 14) then
            warningx = .true.
         end if
         if (iey - isy >= Glr%d%n2 - 14) then
            warningy = .true.
         end if
         if (iez - isz >= Glr%d%n3 - 14) then
            warningz = .true.
         end if

         !If not, then don't use linear input guess (set linear to false)
         if(warningx .and. warningy .and. warningz .and. (Glr%geocode .ne. 'F')) then
            linear = .false.
            if(iproc == 0) then
               write(*,*)'Not using the linear scaling input guess, because localization'
               write(*,*)'region greater or equal to simulation box.'
            end if
            exit 
         end if
      end do

    end subroutine check_linear_inputguess
   
end module locregs_init
