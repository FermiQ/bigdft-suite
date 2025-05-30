!> @file
!!  Files containing the self-consistent loop routines used
!!  for molecular dynamics (see libABINIT/72_geoptim)
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module determining the Self-Consistent Loop API
module scfloop_API

  use bigdft_run, only: run_objects
!!$  use module_base
!!$  use module_types
  implicit none

  private

  ! Storage of required variables for a SCF loop calculation.
  logical :: scfloop_initialised = .false.
  integer :: scfloop_nproc, itime_shift_for_restart
  character(len=*), parameter :: ab6file='ab6_moldyn'
  type(run_objects), pointer, save :: scfloop_obj

  public :: scfloop_init
  public :: scfloop_finalise
  public :: scfloop_initialised, scfloop_obj,itime_shift_for_restart
contains

  subroutine scfloop_init(nproc, obj)
    use f_utils, only: f_open_file
    use abi_defs_basis, only: abi_io_redirect
    implicit none
    integer, intent(in) :: nproc
    type(run_objects), intent(in), target :: obj
    !local variables
    integer :: unt,unt2
    character(len=10) :: suffix

    !this routine set up the file for dump of the moldyn
    unt=7
    call f_open_file(unt,file=trim(ab6log()))
    unt2=unt+1
    call f_open_file(unt2,file=trim(ab6out()))
    !reaffect the value of stdout in abinit
    call abi_io_redirect(new_ab_out=unt,new_std_out=unt2)

    scfloop_nproc = nproc
    scfloop_obj => obj

    scfloop_initialised = .true.
  END SUBROUTINE scfloop_init

  function ab6log() result(suffix)
    use module_base, only: bigdft_mpi
    use yaml_strings, only: yaml_toa,f_strcpy
    implicit none
    character(len=32) :: suffix
    if (bigdft_mpi%ngroup>1) then
       call f_strcpy(src=ab6file//'-'//trim(adjustl(yaml_toa(bigdft_mpi%iproc)))//'-'//&
            trim(adjustl(yaml_toa(bigdft_mpi%igroup)))//'.log',dest=suffix)
    else
       call f_strcpy(src=ab6file//'-'//trim(adjustl(yaml_toa(bigdft_mpi%iproc)))//'.log',dest=suffix)
    end if
  end function ab6log

  function ab6out() result(suffix)
    use module_base, only: bigdft_mpi
    use yaml_strings, only: yaml_toa,f_strcpy
    implicit none
    character(len=32) :: suffix
    if (bigdft_mpi%ngroup>1) then
       call f_strcpy(src=ab6file//'-'//trim(adjustl(yaml_toa(bigdft_mpi%iproc)))//'-'//&
            trim(adjustl(yaml_toa(bigdft_mpi%igroup)))//'.out',dest=suffix)
    else
       call f_strcpy(src=ab6file//'-'//trim(adjustl(yaml_toa(bigdft_mpi%iproc)))//'.out',dest=suffix)
    end if
  end function ab6out


  subroutine scfloop_finalise()
    use f_utils, only: f_close,f_file_unit
    implicit none
    !local variables
    integer :: unit

    call f_file_unit(trim(ab6out()),unit)
    if (unit >0 .and. unit /=6 ) call f_close(unit)

    call f_file_unit(trim(ab6log()),unit)
    if (unit >0 .and. unit /=6 ) call f_close(unit)


  END SUBROUTINE scfloop_finalise
end module scfloop_API


subroutine scfloop_main(acell, epot, fcart, grad, itime, me, natom, rprimd, xred)
  use scfloop_API
  use module_base
  use bigdft_run
  use yaml_output
  use module_input_keys, only: inputpsiid_set_policy
  use public_enums, only: ENUM_MEMORY
  use at_domain, only: domain_periodic_dims
  implicit none

  integer, intent(in) :: natom, itime, me
  real(dp), intent(out) :: epot
  real(dp), intent(in) :: acell(3)
  real(dp), intent(in) :: rprimd(3,3), xred(3,natom)
  real(dp), intent(out) :: grad(3, natom)
  real(dp), intent(out), target :: fcart(3, natom)

  character(len=*), parameter :: subname='scfloop_main'
  logical, dimension(3) :: peri
  integer :: infocode, i, j
  real(dp) :: favg(3)
  type(state_properties) :: outs

  if (.not. scfloop_initialised) then
     write(0,*) "No previous call to scfloop_init(). On strike, refuse to work."
     stop
  end if

  if (me == 0) then
     !write( *,'(1x,a,1x,i0)') &
     !     & 'SCFloop API, call force calculation step=', itime
     call yaml_map('SCFloop API, call force calculation step',itime)
  end if
!!$  if (scfloop_obj%atoms%astruct%geocode=='W') then
!!$     call f_err_throw("Wires bc has to be implemented here", &
!!$          err_name='BIGDFT_RUNTIME_ERROR')
!!$  end if

  ! We transfer acell into at
  scfloop_obj%atoms%astruct%cell_dim(1) = acell(1)
  scfloop_obj%atoms%astruct%cell_dim(2)= rprimd(2,2)
  scfloop_obj%atoms%astruct%cell_dim(3)= acell(3)
  !scfloop_obj%inputs%inputPsiId=1
  call inputpsiid_set_policy(ENUM_MEMORY,scfloop_obj%inputs%inputPsiId)
  ! need to transform xred into xcart
  !peri=bc_periodic_dims(geocode_to_bc(scfloop_obj%atoms%astruct%geocode))
  peri=domain_periodic_dims(scfloop_obj%atoms%astruct%dom)
  do i = 1, scfloop_obj%atoms%astruct%nat, 1
     do j=1,3
        if (peri(j)) then
           scfloop_obj%atoms%astruct%rxyz(j,i)=modulo(xred(j,i),1._gp)*acell(j)
        else
           scfloop_obj%atoms%astruct%rxyz(j,i)=xred(j,i)*acell(j)
        end if
!!$        if (scfloop_obj%atoms%astruct%geocode=='F') then
!!$           scfloop_obj%atoms%astruct%rxyz(j,i)=xred(j,i)*acell(j)
!!$        else
!!$           scfloop_obj%atoms%astruct%rxyz(j,i)=modulo(xred(j,i),1._gp)*acell(j)
!!$        end if
     end do
  end do
!!$  open(100+me)
!!$  write(100+me,*)xcart
!!$  close(100+me)
  call init_state_properties(outs, 1)
  call f_free_ptr(outs%fxyz)
  outs%fxyz => fcart
  !scfloop_obj%inputs%inputPsiId = 1
  call inputpsiid_set_policy(ENUM_MEMORY,scfloop_obj%inputs%inputPsiId)
  call bigdft_state(scfloop_obj,outs,infocode)
  epot = outs%energy
  nullify(outs%fxyz)
  
  call deallocate_state_properties(outs)
  ! need to transform the forces into reduced ones.
  favg(:) = real(0, dp)
  do i = 1, scfloop_obj%atoms%astruct%nat, 1
     favg(:) = favg(:) + fcart(:, i) / real(natom, dp)
     grad(:, i) = -fcart(:, i) / acell(:)
  end do
  do i = 1, scfloop_obj%atoms%astruct%nat, 1
     fcart(:, i) = fcart(:, i) - favg(:)
  end do
END SUBROUTINE scfloop_main


subroutine scfloop_output(acell, epot, ekin, fred, itime, me, natom, rprimd, vel, xred)
  use scfloop_API
  use module_base
  !use module_types
  use module_atoms, only: astruct_dump_to_file

  implicit none

  !Arguments
  integer, intent(in) :: natom, itime, me
  real(dp), intent(in) :: epot, ekin
  real(dp), dimension(3), intent(in) :: acell
  real(dp), dimension(3,natom), intent(in) :: xred
  real(dp), dimension(3,natom), intent(in) :: fred, vel, rprimd
  !Local variables
  character(len=*), parameter :: subname='scfloop_output'
  character(len = 5) :: fn5
  character(len = 40) :: comment
  real :: fnrm
  real(dp), dimension(:,:), allocatable :: xcart,fcart
  integer :: i

  if (me /= 0) return

  fnrm = real(0, dp)
  ! need to transform xred into xcart
  xcart = f_malloc((/ 3, natom /),id='xcart')
  fcart = f_malloc((/ 3, natom /),id='fcart')

  do i = 1, natom
     xcart(:, i) = xred(:, i) * acell(:)
     fnrm = fnrm + real(fred(1, i) * acell(1) * fred(1, i) * acell(1) + &
          & fred(2, i) * acell(2) * fred(2, i) * acell(2) + &
          & fred(3, i) * acell(3) * fred(3, i) * acell(3))
     fcart(:, i) = fred(:, i) * acell(:)
  end do

  write(fn5,'(i5.5)') itime+itime_shift_for_restart
  write(comment,'(a,1pe10.3)')'AB6MD:fnrm= ', sqrt(fnrm)
  call astruct_dump_to_file(scfloop_obj%atoms%astruct,&
       trim(scfloop_obj%inputs%dir_output)//'posmd_'//fn5, &
       trim(comment),&
       epot + ekin,rxyz=xcart,forces=fcart)
!!$  call write_atomic_file(trim(scfloop_obj%inputs%dir_output)//'posmd_'//fn5, &
!!$       & epot + ekin, xcart, scfloop_obj%atoms%astruct%ixyz_int, scfloop_obj%atoms, trim(comment),forces=fcart)

  !write velocities
  write(comment,'(a,i6.6)')'Timestep= ',itime+itime_shift_for_restart
  call wtvel('velocities.xyz',vel,scfloop_obj%atoms,comment)

  call f_free(xcart)
  call f_free(fcart)

  
  !To avoid warning from compiler
  fnrm=real(rprimd(1,1),kind=4)

END SUBROUTINE scfloop_output

!>    Read atomic positions
subroutine read_velocities(iproc,filename,atoms,vxyz)
  use scfloop_API
  use module_base
  use module_types
  use at_domain, only: domain_geocode
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(out) :: vxyz
  !local variables
  !n(c) character(len=*), parameter :: subname='read_velocities'
  character(len=2) :: symbol
  character(len=20) :: tatonam,units
  character(len=50) :: extra
  character(len=150) :: line
  logical :: lpsdbl,exists
  integer :: iat,i,ierrsfx,nat,unt
! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: vx,vy,vz,alat1,alat2,alat3
! case for which the atomic positions are given whithin general precision
  real(gp) :: vxd0,vyd0,vzd0,alat1d0,alat2d0,alat3d0

  !inquire whether the input file is present, otherwise put velocities to zero
  call f_file_exists(filename,exists)
  !inquire(file=filename,exist=exists)
  if (.not. exists) then  
     call f_zero(vxyz)
     return
  end if

  !controls if the positions are provided with machine precision
  if (atoms%astruct%units== 'atomicd0' .or. atoms%astruct%units== 'bohrd0') then
     lpsdbl=.true.
  else
     lpsdbl=.false.
  end if

  unt=99
  call f_open_file(unt,file=trim(filename),status='old',action='read')
  !open(unit=99,file=trim(filename),status='old',action='read')

  read(unt,*) nat,units,extra,itime_shift_for_restart
 
  !check whether the number of atoms is different 
  if (nat /= atoms%astruct%nat) then
     if (iproc ==0) write(*,*)' ERROR: the number of atoms in the velocities is different'
     stop
  end if

  !read from positions of .xyz format, but accepts also the old .ascii format
  read(unt,'(a150)')line

  if (lpsdbl) then
     read(line,*,iostat=ierrsfx) tatonam,alat1d0,alat2d0,alat3d0
  else
     read(line,*,iostat=ierrsfx) tatonam,alat1,alat2,alat3
  end if

  !convert the values of the cell sizes in bohr
  if (units=='angstroem' .or. units=='angstroemd0') then
  else if  (units=='atomic' .or. units=='bohr'  .or.&
       units== 'atomicd0' .or. units== 'bohrd0') then
  else if (units == 'reduced') then
     !assume that for reduced coordinates cell size is in bohr
  else
     write(*,*) 'length units in input file unrecognized'
     write(*,*) 'recognized units are angstroem or atomic = bohr'
     stop 
  endif
  do iat=1,atoms%astruct%nat
     !xyz input file, allow extra information
     read(unt,'(a150)')line 
     if (lpsdbl) then
        read(line,*,iostat=ierrsfx)symbol,vxd0,vyd0,vzd0
     else
        read(line,*,iostat=ierrsfx)symbol,vx,vy,vz
     end if
     tatonam=trim(symbol)
     if (lpsdbl) then
        vxyz(1,iat)=vxd0
        vxyz(2,iat)=vyd0
        vxyz(3,iat)=vzd0
     else
        vxyz(1,iat)=real(vx,gp)
        vxyz(2,iat)=real(vy,gp)
        vxyz(3,iat)=real(vz,gp)
     end if
 
     if (units=='angstroem' .or. units=='angstroemd0') then
        ! if Angstroem convert to Bohr
        do i=1,3 
           vxyz(i,iat)=vxyz(i,iat)/Bohr_Ang
        enddo
     else if (units == 'reduced') then 
        vxyz(1,iat)=vxyz(1,iat)*atoms%astruct%cell_dim(1)
        if (domain_geocode(atoms%astruct%dom) == 'P') vxyz(2,iat)=vxyz(2,iat)*atoms%astruct%cell_dim(2)
        vxyz(3,iat)=vxyz(3,iat)*atoms%astruct%cell_dim(3)
     endif
  enddo

  call f_close(unit=99)
END SUBROUTINE read_velocities

subroutine wtvel(filename,vxyz,atoms,comment)
  use module_base
  use module_types
  use module_atoms, only: write_xyz_bc
  use at_domain, only: domain_geocode
  implicit none
  character(len=*), intent(in) :: filename,comment
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: vxyz
  !Local variables
  integer, parameter :: iunit = 9
  character(len=2) :: symbol
  character(len=10) :: name,bcname
  character(len=11) :: units
  integer :: iat,j
  real(gp) :: factor

  open(unit=iunit,file=trim(filename),status='unknown',action='write')
  if (trim(atoms%astruct%units) == 'angstroem' .or. trim(atoms%astruct%units) == 'angstroemd0') then
     factor=Bohr_Ang
     units='angstroemd0'
  else
     factor=1.0_gp
     units='atomicd0'
  end if

  write(iunit,'(i6,2x,a,2x,a)') atoms%astruct%nat,trim(units),comment  

  call write_xyz_bc(iunit,domain_geocode(atoms%astruct%dom),factor,atoms%astruct%cell_dim)

!!$  if (atoms%astruct%geocode == 'P') then
!!$     write(iunit,'(a,3(1x,1pe24.17))') 'periodic',&
!!$       &  atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,atoms%astruct%cell_dim(3)*factor
!!$  else if (atoms%astruct%geocode == 'S') then
!!$     write(iunit,'(a,3(1x,1pe24.17))') 'surface',&
!!$       &  atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,atoms%astruct%cell_dim(3)*factor
!!$  else if (atoms%astruct%geocode == 'W') then
!!$     call f_err_throw("Wires bc has to be implemented here", &
!!$          err_name='BIGDFT_RUNTIME_ERROR')
!!$  else
!!$     write(9,*)'free'
!!$  end if
  do iat=1,atoms%astruct%nat
     name=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:2)
     end if

     write(iunit,'(a2,4x,3(1x,1pe24.17))') symbol,(vxyz(j,iat)*factor,j=1,3)

  enddo

  close(unit=iunit)
END SUBROUTINE wtvel
