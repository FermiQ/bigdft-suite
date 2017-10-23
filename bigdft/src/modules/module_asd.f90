!> @file
!! Handling of the constrained magnetic field of the system
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module module_asd
   use module_base
   !use abSpinlib
   implicit none

   private

   !> loop counter
   real(dp), parameter :: ha2tesla = 470103.892795094_dp
   !
   !
   !
   !>associated to an instance of an ASD simulation
   type, public :: asd_data
      !
      !> control flag for ASD evolution (predictor)
      logical :: do_pred
      !> control flag for ASD evolution (corrector)
      logical :: do_corr
      !> control flag for ASD evolution (corrector)
      logical :: in_progress
      !> Current timestep of the simulation
      integer :: mstep=0
      !> time in steps of the simulation (simulation time = nstep*delta_t)
      integer :: nstep=0
      !> temperature for the spin system
      real(dp) :: spin_temp
      !> ASD time step
      real(dp) :: delta_t
      !> Gilbert damping for ASD
      real(dp) :: damping
      !> Fudge factor for converting magnetic field to correct unit
      real(dp) :: bfield_scale = 1.0_dp
      !
      !!> value of the magnetic field close to each of the centers
      !real(gp), dimension(:,:), pointer :: B_at => null()
      !!> local magnetization of each of the centers
      !real(gp), dimension(:,:), pointer :: m_at => null()
      !
   end type asd_data

   type(asd_data) :: asd
   public :: asd_allocate,asd_free, asd
   public :: asd_read_external, asd_wrapper

contains

   pure function asd_data_null() result(asd)
      implicit none
      type(asd_data) :: asd
      call nullify_asd_data(asd)
   end function asd_data_null

   pure subroutine nullify_asd_data(asd)
      implicit none
      type(asd_data), intent(out) :: asd

      asd%mstep=0
      asd%nstep=0
      asd%delta_t=0.0_dp
      asd%spin_temp=0.0_dp
      asd%damping=0.0_dp
      ! bfield_scale put to zero below. Potential risk for bugs if kept like that
      asd%bfield_scale=1.0_dp
      !!!$ nullify(asd%B_at)
      !!!$ nullify(asd%m_at)

   end subroutine nullify_asd_data

   subroutine asd_free(asd)
      implicit none
      type(asd_data), intent(inout) :: asd

      !!!$ call f_free_ptr(asd%B_at)
      !!!$ call f_free_ptr(asd%m_at)
      !!!$ call f_free_ptr(asd%rho_at)
      !
      call nullify_asd_data(asd)

   end subroutine asd_free

   subroutine asd_allocate(asd,mstep, spin_temp,delta_t,damping,bfield_scale)
      implicit none
      type(asd_data), intent(inout) :: asd
      integer, intent(in),optional :: mstep
      real(dp), intent(in), optional :: spin_temp
      real(dp), intent(in), optional :: delta_t
      real(dp), intent(in), optional :: damping
      real(dp), intent(in), optional :: bfield_scale

      call asd_free(asd) !we can as the initial status of the data is defined
      asd%in_progress=.true.

      if(present(mstep)) asd%mstep=mstep
      if(present(spin_temp)) asd%spin_temp=spin_temp
      if(present(delta_t)) asd%delta_t=delta_t
      if(present(damping)) asd%damping=damping
      if(present(bfield_scale)) then
         asd%bfield_scale=bfield_scale
      else
         asd%bfield_scale=1.0d0
      end if
      ! Might modify scaling later on
      !asd%bfield_scale=asd%bfield_scale*ha2tesla

      !!!$ asd%B_at=f_malloc_ptr([3,nat],id='asd%B_at')
      !!!$ asd%m_at=f_malloc_ptr([3,nat],id='asd%m_at')
      !

   end subroutine asd_allocate

   !!!$ subroutine asd_dump_info(asd)
   !!!$    use yaml_output
   !!!$    implicit none
   !!!$    type(asd_data), intent(in) :: asd
   !!!$    !local variables
   !!!$    integer :: iat

   !!!$    call yaml_newline()
   !!!$    call yaml_sequence_open('Local information on the magnetic centers')
   !!!$    do iat=1,asd%nat
   !!!$       call yaml_newline()
   !!!$       call yaml_sequence(advance='no')
   !!!$       call yaml_mapping_open(flow=.true.)
   !!!$       call yaml_map('R',asd%rxyz(:,iat)) !position
   !!!$       call yaml_map('D',asd%radii(iat)) !radius
   !!!$       call yaml_newline()
   !!!$       call yaml_map('M',asd%m_at(:,iat),fmt='(1pe12.5)') !mag mom
   !!!$       call yaml_map('C',asd%rho_at(iat)) !charge
   !!!$       call yaml_mapping_close()
   !!!$    end do
   !!!$    call yaml_sequence_close()
   !!!$    call yaml_newline()

   !!!$ end subroutine asd_dump_info


   subroutine asd_read_external(asd,fname)
      !
      !
      implicit none
      !
      type(asd_data), intent(inout) :: asd  !< the currently used asd object
      character*30 :: fname !< name of external file containing constraining moment directions
      !
      integer :: iat
      !
      open(unit=22,file=fname,form='formatted',action='read',status='old')
      !
      read(22,*,END=100) asd%nstep,asd%spin_temp,asd%delta_t,asd%damping,asd%bfield_scale !asd_temp,asd_dt,asd_damp, asd_alpha!, xc_int
      print '(a,i6,4g14.6)','asd.in contents: ',asd%nstep,asd%spin_temp,asd%delta_t,asd%damping,asd%bfield_scale
      100   continue
      !
      close(22)
      !
      return
      !
   end subroutine asd_read_external

   subroutine asd_wrapper(asd,moments,bfield,nat)
      !
      !
      implicit none
      !
      type(asd_data), intent(inout) :: asd !< the currently used asd object
      integer, intent(in) :: nat
      real(dp), dimension(3,nat), intent(inout) :: moments
      real(dp), dimension(3,nat), intent(in) :: bfield
      !
      character*30 :: fname='asd.in'
      !
      ! Check if in beginning or end of simulation
      print *,'In asd_wrapper ',asd%mstep, asd%nstep, asd%in_progress
      if(asd%mstep==0 .and. .not. asd%in_progress) then
         asd%in_progress=.true.
         asd%do_pred=.true.
         asd%do_corr=.false.
         call allocate_asd(nat,1)
         call asd_read_external(asd,fname)
      else if(asd%mstep>asd%nstep) then
         asd%in_progress=.false.
         asd%do_pred=.false.
         asd%do_corr=.false.
         call allocate_asd(nat,-1)
         return
      end if

      !!!$   print *,'ASD scaling:',asd%bfield_scale,asd%bfield_scale*ha2tesla
      !!!$   print *,'asd_moments:'
      !!!$   print '(3f12.6)' , moments
      !!!$   print *,'asd_fields:'
      !!!$   print '(3g14.6)' , bfield
      !!!$   print *,'asd_fields_scaled 1:'
      !!!$   print '(3g14.6)' , bfield*asd%bfield_scale*ha2tesla
      !!!$   print *,'asd_fields_scaled 1:'
      !!!$   print '(3g14.6)' , bfield*ha2tesla
      !!!$   print *,'asd_do_pred:',asd%do_pred
      !!!$   print *,'asd_do_corr:',asd%do_corr

      if(asd%do_pred) then
         call asd_pred(moments,bfield*asd%bfield_scale*ha2tesla,asd%spin_temp,asd%delta_t,nat)
         print '(3f12.6)', moments
         asd%do_pred=.false.
         asd%do_corr=.true.
      else if(asd%do_corr) then
         call asd_corr(moments,bfield*asd%bfield_scale*ha2tesla,asd%spin_temp,asd%delta_t,nat)
         print '(3f12.6)', moments
         asd%do_pred=.true.
         asd%do_corr=.false.
         asd%mstep=asd%mstep+1
      end if
      !
      return
   end subroutine asd_wrapper

   !
end module module_asd
