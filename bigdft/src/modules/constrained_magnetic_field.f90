!> @file
!! Handling of the constrained magnetic field of the system
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module module_cfd
   use module_base
   implicit none

   private

   !> temporary output file number
   integer :: stdout = 666
   !> type of constraining algorithm (2=regular Lagrange,3=orthogonal Lagrange,4=PID,5=Ma-Dudarev)
   integer :: i_cons = 0
   !> dimensionality of magnetism (3=xyz)
   integer, parameter :: ncomp=3
   !> conversion of units from Ry to Tesla. Used for EOM solver (need to adjust for current energy unit)
   real(gp), parameter :: b2t = 235298.924212429_gp
   !> conversion of units from Ry to Tesla. Used for EOM solver (need to adjust for current energy unit)
   real(gp), parameter :: cfd_thresh = 5.0e-7_gp
   !
   !> Lagrange penalty factor (to be modified and moved)
   real(gp) :: lambda = 10.0_gp
   !> temporary Lagrange penalty factor (to be modified and moved)
   real(gp) :: lambda_t = 1.0_gp
   !> starting Lagrange penalty factor (to be modified and moved)
   real(gp) :: lambda_0 = 1.0_gp
   !> Threshold for determining induced moments (currently not constraining induced moments)
   real(gp) :: induced_mom_thresh = 0.5_gp
   !> prefactor for future use 
   real(gp) :: cfd_prefac = 1.0_gp
   !> mixing factor for constraining b-field
   real(gp) :: B_at_beta = 0.9_gp
   !> loop counter
   integer :: b_constr_iter
   !
   !


   !>associated to an instance of the cfd calculation
   type, public :: cfd_data
      !> number of magnetic centers
      integer :: nat=0
      !> position of the centers in the simulation domain
      real(gp), dimension(:,:), pointer :: rxyz => null()
      !> radius of each of the magnetic atoms
      real(gp), dimension(:), pointer :: radii => null()
      !> value of the magnetic field close to each of the centers
      real(gp), dimension(:,:), pointer :: B_at => null()
      !> local magnetization of each of the centers
      real(gp), dimension(:,:), pointer :: m_at => null()
      !> electronic charge inside each of the center
      real(gp), dimension(:), pointer :: rho_at => null()
      !
      !> reference magnetization of each of the centers
      real(gp), dimension(:,:), pointer :: m_at_ref => null()
      !> arrays for PID
      real(gp), dimension(:,:), pointer :: d_delta => null()
      real(gp), dimension(:,:), pointer :: s_delta => null()
      real(gp), dimension(:,:), pointer :: dd_delta => null()
      !> error in the constrained moments
      real(gp) :: constrained_mom_err
      !
   end type cfd_data

   public :: cfd_allocate,cfd_free,cfd_set_radius,cfd_dump_info
   public :: cfd_set_centers, cfd_read_external, cfd_field
   public :: cfd_is_converged

contains

   pure function cfd_data_null() result(cfd)
      implicit none
      type(cfd_data) :: cfd
      call nullify_cfd_data(cfd)
   end function cfd_data_null

   pure subroutine nullify_cfd_data(cfd)
      implicit none
      type(cfd_data), intent(out) :: cfd
      cfd%nat=0
      nullify(cfd%rxyz)
      nullify(cfd%radii)
      nullify(cfd%B_at)
      nullify(cfd%m_at)
      nullify(cfd%rho_at)
      !
      nullify(cfd%m_at_ref)
      nullify(cfd%d_delta)
      nullify(cfd%s_delta)
      nullify(cfd%dd_delta)
   end subroutine nullify_cfd_data

   subroutine cfd_free(cfd)
      implicit none
      type(cfd_data), intent(inout) :: cfd

      call f_free_ptr(cfd%rxyz)
      call f_free_ptr(cfd%radii)
      call f_free_ptr(cfd%B_at)
      call f_free_ptr(cfd%m_at)
      call f_free_ptr(cfd%rho_at)
      !
      call f_free_ptr(cfd%m_at_ref)
      call f_free_ptr(cfd%d_delta)
      call f_free_ptr(cfd%s_delta)
      call f_free_ptr(cfd%dd_delta)
      !
      call nullify_cfd_data(cfd)
      close(unit=stdout)
   end subroutine cfd_free

   subroutine cfd_allocate(cfd,nat)
      implicit none
      integer, intent(in) :: nat
      type(cfd_data), intent(inout) :: cfd

      call cfd_free(cfd) !we can as the initial status of the data is defined

      cfd%nat=nat
      cfd%rxyz=f_malloc_ptr([3,nat],id='cfd%rxyz')
      cfd%radii=f_malloc_ptr(nat,id='cfd%radii')
      cfd%B_at=f_malloc_ptr([3,nat],id='cfd%B_at')
      cfd%m_at=f_malloc_ptr([3,nat],id='cfd%m_at')
      cfd%rho_at=f_malloc_ptr(nat,id='cfd%rho_at')
      !
      cfd%m_at_ref=f_malloc_ptr([3,nat],id='cfd%m_at_ref')
      cfd%d_delta=f_malloc_ptr([3,nat],id='cfd%d_delta')
      cfd%s_delta=f_malloc_ptr([3,nat],id='cfd%s_delta')
      cfd%dd_delta=f_malloc_ptr([3,nat],id='cfd%dd_delta')

      open(file='cfd.log',unit=stdout)
   end subroutine cfd_allocate

   !!$  function cfd_get_centers_ptr(cfd) result(ptr)
   !!$    implicit none
   !!$    type(cfd_data), intent(inout) :: cfd
   !!$    real(gp), dimension(:,:), pointer :: ptr
   !!$
   !!$    ptr => cfd%rxyz
   !!$
   !!$  end function cfd_get_centers_ptr

   subroutine cfd_set_centers(cfd,rxyz)
      implicit none
      type(cfd_data), intent(inout) :: cfd
      real(gp), dimension(3,cfd%nat), intent(in) :: rxyz

      call f_memcpy(src=rxyz,dest=cfd%rxyz)

   end subroutine cfd_set_centers


   pure subroutine cfd_set_radius(cfd,iat,radius)
      implicit none
      integer, intent(in) :: iat
      real(gp), intent(in) :: radius
      type(cfd_data), intent(inout) :: cfd

      cfd%radii(iat)=radius
   end subroutine cfd_set_radius

   subroutine cfd_dump_info(cfd)
      use yaml_output
      implicit none
      type(cfd_data), intent(in) :: cfd
      !local variables
      integer :: iat

      call yaml_newline()
      call yaml_sequence_open('Local information on the magnetic centers')
      do iat=1,cfd%nat
         call yaml_newline()
         call yaml_sequence(advance='no')
         call yaml_mapping_open(flow=.true.)
         call yaml_map('R',cfd%rxyz(:,iat)) !position
         call yaml_map('D',cfd%radii(iat)) !radius
         call yaml_newline()
         call yaml_map('M',cfd%m_at(:,iat),fmt='(1pe12.5)') !mag mom
         call yaml_map('C',cfd%rho_at(iat)) !charge
         call yaml_mapping_close()
      end do
      call yaml_sequence_close()
      call yaml_newline()

   end subroutine cfd_dump_info



   !
   subroutine cfd_field(cfd,iproc)
      !
      implicit none
      !
      type(cfd_data), intent(inout) :: cfd
      integer, intent(in) :: iproc !< Label of the process,from 0 to nproc-1
      !
      ! arguments
      !integer, intent(in) :: ndim
      !real(gp), intent(in) :: cfd%m_at(ncomp,ndim)
      !real(gp), intent(in) :: cfd%m_at_ref(ncomp,ndim)
      !real(gp), intent(inout) :: cfd%B_at(ncomp,ndim)
      !
      !
      integer :: iidim,na
      real(gp), dimension(:,:), allocatable :: mom_tmp
      real(gp), dimension(3) :: e_i, e_out, c_in
      real(gp), dimension(3) :: m_delta, B_at_new, B_diff
      real(gp) :: etcon, ma, mnorm
      ! PID factors (e_k in ref = m_delta here)
      real(gp) :: K=-0.3d0, T_s= 0.3d0, T_i=3.0d0, T_d=-0.1d0
      !real(gp) :: K=-0.3d0, T_s= 0.3d0, T_i=1.5d0, T_d=-0.1d0 !best so far
      ! u_k = K*(e_k+I_k-1+T_s/T_i*e_k + T_d/T_s*(e_k-e_k_1))

      !!!   cfd_prefac=b2t*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
      if(i_cons==0) return

      if(cfd_is_converged(cfd)) then
         if(iproc==0) print *,"CFD converged, error:",cfd%constrained_mom_err
         return
      else
         if(iproc==0) print *,"CFD not converged, error:",cfd%constrained_mom_err
      end if

      allocate ( mom_tmp(ncomp,cfd%nat))


      cfd%constrained_mom_err=0.0_gp
      do na=1,cfd%nat
         if (i_cons==2) then
            ! Lagrange multiplier without orthogonalization 
            !
            mom_tmp(:,na) = cfd%m_at(:,na)/norm2(cfd%m_at(:,na)) - cfd%m_at_ref(:,na)
            ! Gramm-Schmidt step
            etcon = etcon + lambda_t * sum(mom_tmp(:,na)**2)
         else if (i_cons==3) then
            ! Lagrange multiplier with orthogonalization (b _|_ m)
            !
            mom_tmp(:,na) = cfd%m_at(:,na)/norm2(cfd%m_at(:,na)) - cfd%m_at_ref(:,na)
            ! Gramm-Schmidt step
            mom_tmp(:,na) = mom_tmp(:,na) - sum(mom_tmp(:,na)*cfd%m_at_ref(:,na)) * cfd%m_at_ref(:,na)
            etcon = etcon + lambda_t * sum(mom_tmp(:,na)**2)
         else if (i_cons==4) then
            ! i_cons = 4 means that we try to use a PID regulator
            !
            !
            !
            if(iproc==0) write (stdout,'(4x,a)') ' | AMN-PID noncolinear constraints '
            !
            ! Check moment magnitude
            mnorm = sqrt(cfd%m_at(1,na)**2+cfd%m_at(2,na)**2+cfd%m_at(3,na)**2)
            if(iproc==0) write (stdout,'(4x,a,i4)' ) " | - atom: ", na
            if (mnorm.lt.induced_mom_thresh) then
               if(iproc==0) write (stdout,'(2x,a,i4,a,f10.4)') ' | Local magnetization for atom ', na , ' is less than threshold',mnorm
               m_delta=0.0_gp
            else
               c_in=cfd%B_at(:,na)
               ! Direction only
               e_out=cfd%m_at(:,na)/mnorm
               e_i =cfd%m_at_ref(:,na)/norm2(cfd%m_at_ref(:,na))
               !! Direction and magnitude
               !e_out=cfd%m_at(:,na)
               !e_i =cfd%m_at_ref(:,na)
               ! P I D
               ! Full direction (untested)
               !m_delta=(e_i-e_out)
               ! Perp direction (works for bcc fe)
               !m_delta=-(cfd%m_at(:,na)-sum(cfd%m_at(:,na)*e_i)*e_i)
               m_delta=(cfd%m_at(:,na)-sum(cfd%m_at(:,na)*e_i)*e_i)
               !m_delta=-(e_out-norm2(e_out*e_i)*e_i)
            end if
            !
            ! Reducing the effect for first iteration (ie when cfd%d_delta=0)
            if(norm2(cfd%d_delta(:,na))<1e-15) m_delta=0.1_gp*m_delta
            ! e) m_delta=-lambda_t*(cfd%m_at(:,na)-sum(cfd%m_at(:,na)*e_i)*e_i)*10.0_gp
            ! others:lambda_t=0.1
            !gs
            !m_delta=-(e_out-norm2(e_out*e_i)*e_i)
            !
            ! Check to don't mix first iteration (dd_delta=derivative)
            if(norm2(cfd%d_delta(:,na))>1e-15) cfd%dd_delta(:,na)=m_delta-cfd%d_delta(:,na)
            !
            if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | Output moments     for atom ",na,cfd%m_at(:,na)
            if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | Target direction   for atom ",na,cfd%m_at_ref(:,na)/norm2(cfd%m_at_ref(:,na))
            if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | Outut direction    for atom ",na,e_out
            if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | Input field        for atom ",na,cfd%B_at(:,na)
            !
            ! Check to don't mix first iteration (s_delta=integral=I_k)
            if(norm2(cfd%d_delta(:,na))>1e-15) cfd%s_delta(:,na)=cfd%s_delta(:,na)+T_s/T_i*m_delta!*0.5_gp

            !cfd%B_at(:,na)=lambda_t*(1.20_gp*m_delta+0.35_gp*cfd%s_delta(:,na)+0.10_gp*cfd%dd_delta(:,na))
            !cfd%B_at(:,na)=lambda_t*(1.30_gp*m_delta+0.35_gp*cfd%s_delta(:,na)-0.10_gp*cfd%dd_delta(:,na))   !<-- use this for atoms
            !cfd%B_at(:,na)=lambda_t*(1.30_gp*m_delta+0.35_gp*cfd%s_delta(:,na)-0.10_gp*cfd%dd_delta(:,na))   !<-- use this for !atoms bigdft best
            !B_diff=lambda_t*(2.00_gp*m_delta+0.20_gp*cfd%s_delta(:,na)-0.20_gp*cfd%dd_delta(:,na)) 
            B_diff=K*(m_delta+cfd%s_delta(:,na)+T_d/T_s*cfd%dd_delta(:,na))
            cfd%constrained_mom_err=cfd%constrained_mom_err+sum((cfd%B_at(:,na)-B_diff)**2)
            cfd%B_at(:,na)=B_diff

            ! Calculate Zeeman-like constraining energy cost
            etcon = etcon + sum(cfd%B_at(:,na)*cfd%m_at(:,na))

            ! d_delta=error=e_k
            cfd%d_delta(:,na)=m_delta!*0.5_gp
            !
            !if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | P  contribution    for atom ",na,cfd%d_delta(:,na)
            if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | P  contribution    for atom ",na,m_delta*K
            if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | I  contribution    for atom ",na,cfd%s_delta(:,na)
            if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | D  contribution    for atom ",na,cfd%dd_delta(:,na)*T_d/T_s
            if(iproc==0) write (stdout,'(4x,a,i4,4f15.8)' ) " | Constraining field for atom ",na,cfd%B_at(:,na)
            if(iproc==0) write (stdout,'(4x,a,i4,3f15.4)' ) " | Constraining field for atom (t)",na,cfd_prefac*cfd%B_at(:,na)
            if(iproc==0 .and. na==cfd%nat) write (stdout,'(4x,a,2f15.4,g14.6)' ) " | Lambda prefactors ",lambda_t,lambda

         else if (i_cons==5) then
            ! i_cons = 5 means that we try the Ma-Dudarev approach
            ! which is very analogous to the normal Lagrange approach
            !
            !
            if(iproc==0) write (stdout,'(2x,a)') ' Ma-Dudarev constraints '
            !
            ! Check moment magnitude
            ma = dsqrt(cfd%m_at(1,na)**2+cfd%m_at(2,na)**2+cfd%m_at(3,na)**2)
            if(iproc==0) write (stdout,'(4x,a,i4)' ) " | - Atom: ", na
            !
            if (ma.lt.induced_mom_thresh) then
               if(iproc==0) write (stdout,'(2x,a,i4,a,f10.4)') ' | Local magnetization for atom ', na , ' is less than threshold',ma
               cfd%B_at(:,na)=0.0_gp
            else
               if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | Output moments     for atom ",na,cfd%m_at(:,na)
               if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | Input direction    for atom ",na,cfd%m_at_ref(:,na)
               if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | Input field        for atom ",na,cfd%B_at(:,na)
               !
               e_out=cfd%m_at(:,na)
               e_i =cfd%m_at_ref(:,na)/norm2(cfd%m_at_ref(:,na))*ma
               !e_dot=e_out(1)*e_i (1)+e_out(2)*e_i (2)+e_out(3)*e_i (3)
               ! new constraining field
               !cfd%B_at_new=-10.0_gp*lambda_t*(e_out-e_i)
               B_at_new=lambda_t*(e_out-e_i)
               ! gram-schmidt orthogonalization ( b _|_ m)
               !B_at_new=B_at_new-sum(e_out*B_at_new)*e_out/ma**2
               B_at_new=B_at_new-sum(e_i*B_at_new)*e_i/ma**2
               ! mixing field
               cfd%B_at(:,na)=(1.0_gp-B_at_beta)*cfd%B_at(:,na)-B_at_beta*B_at_new
               !
               ! calculate zeeman-like constraining energy cost
               etcon = etcon + lambda_t*(sqrt(sum(cfd%m_at(:,na)*cfd%m_at(:,na)))-sum(cfd%m_at(:,na)*e_i/ma))
               cfd%constrained_mom_err = cfd%constrained_mom_err + sum(e_out-e_i)**2
               ! if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | new field          for atom ",na,B_at_new
               if(iproc==0) write (stdout,'(4x,a,i4,4f15.8)' ) " | Output field       for atom ",na,cfd%B_at(:,na),&
                  sum(cfd%B_at(:,na)*e_i)
               if(iproc==0) write (stdout,'(4x,a,i4,4f15.8)' ) " | Output field (t)   for atom ",na, &
                  cfd%B_at(:,na)*cfd_prefac, &
                  sum(cfd%B_at(:,na)*e_i)
               if(iproc==0) write (stdout,'(4x,a,i4,3f15.8)' ) " | Output direction    for atom ",na,e_out/ma
            end if
         end if
      end do ! na
      !
      cfd%constrained_mom_err=sqrt(cfd%constrained_mom_err)/cfd%nat

      b_constr_iter=b_constr_iter+1
      !if(i_cons==5) lambda_t=min(lambda_t+4_gp,100.0_gp)
      !if(i_cons==5) lambda_t=min(lambda_t+1.0_gp-lambda_t/lambda,lambda)
      !if(i_cons==3) lambda_t=min(lambda_t+1.0_gp-lambda_t/lambda,lambda)
      ! works for moderate lambdas
      !if(i_cons==3) then
      !   if (b_constr_iter<=30.0_gp) then
      !      lambda_t=min(lambda_t*(2.0_gp-b_constr_iter/30.0_gp),lambda)
      !   else
      !      lambda_t=lambda
      !   end if
      !end if
      !
      !
      ! scale up lambda in case of Lagrangian formulation
      if(i_cons==5.and.etcon<1.0d-2) lambda_t=min(1.2_gp*lambda_t,1.0e4_gp)
      if(i_cons==3) lambda_t=min(lambda_t*(2.0_gp-lambda_t/lambda),lambda)
      !if(i_cons==3) lambda_t=min(lambda_t*(2.0_gp-lambda_t/lambda),lambda)
      !if(i_cons==5) lambda_t=min(lambda_t*(2.0_gp-lambda_t/lambda),lambda)
      !if(i_cons==5) lambda_t=min(lambda_t+2.0_gp,lambda)
     ! if(i_cons==5) lambda_t=min(lambda_t+1.0_gp,lambda)
      ! if(i_cons==4) lambda_t=min(lambda_t+1.0_gp,25.0_gp)
      !if(i_cons==4) lambda_t=lambda_t+1.0_gp
      !if(i_cons==5) lambda_t=lambda_t+lambda_t*min(0.1_gp,etcon**2)
      !if(i_cons==5) lambda_t=lambda_t*(1.0_gp+0.5_gp*(etcon))
      !if(i_cons==5) lambda_t=lambda_t*(1.0_gp+2.0_gp*min(cfd%constrained_mom_err,0.1_gp))
      if(i_cons==5.or.i_cons==4) then
         if(iproc==0) write (stdout,'(4x,a,f12.4a,g10.2)' ) " | New lambda_t: ", lambda_t, "     error: ", cfd%constrained_mom_err
         if(iproc==0) write (stdout,'(4x,a,f12.4a,g10.2)' ) " | New lambda  : ", lambda  
      end if
      if(iproc==0) write (stdout,'(4x,a)' ) " | -  "
      deallocate(mom_tmp)
      return
   end subroutine cfd_field

   subroutine cfd_read_external(cfd,fname)
      !
      !
      implicit none
      !
      type(cfd_data), intent(inout) :: cfd  !< the currently used cfd object
      character*30 :: fname !< name of external file containing constraining moment directions
      !
      integer :: iat
      !
      open(unit=22,file=fname,form='formatted',action='read',status='old')
      !
      !First read the reference directions 
      do iat=1,cfd%nat
         read(22,*) cfd%m_at_ref(1:3,iat)
         !print '(a,3f12.6)',' CFD ', cfd%m_at_ref(1:3,iat)
      end do
      !
      ! Then read controlling flags
      read(22,*,end=100) i_cons  ! which constraining scheme to choose 
      read(22,*,end=100) lambda_0,lambda  ! Lagrange factor for certains schemes
      read(22,*,end=100) cfd_prefac ! Fudge factor for testing the magnitude of the field
      lambda_t=lambda_0
      100 continue
      !
      close(22)
      !
      return
      !
   end subroutine cfd_read_external
   !
   function cfd_is_converged(cfd)
      !
      !
      implicit none
      !
      type(cfd_data), intent(inout) :: cfd  !< the currently used cfd object
      !
      logical :: cfd_is_converged
      real(gp), dimension(3) :: e_in, e_out
      integer :: iat
      !
      !cfd%constrained_mom_err=sqrt(sum((cfd%m_at(1:3,1:cfd%nat)-cfd%m_at_ref(1:3,1:cfd%nat))**2))
      !
      !!! cfd%constrained_mom_err=0.0_gp
      !!! do iat=1,cfd%nat
      !!!    e_in=cfd%m_at_ref(:,iat)
      !!!    e_in=e_in/norm2(e_in+1.0e-12_gp)
      !!!    e_out=cfd%m_at(:,iat)
      !!!    e_out=e_out/norm2(e_out+1.0e-12_gp)
      !!!    cfd%constrained_mom_err=cfd%constrained_mom_err+sum((e_in-e_out)**2)
      !!!    !print '(7f12.6)',e_in,e_out,cfd%constrained_mom_err
      !!! end do
      !!! cfd%constrained_mom_err=sqrt(cfd%constrained_mom_err)
      !
      !cfd_is_converged=cfd%constrained_mom_err<5.0d-4
      !cfd_is_converged=cfd%constrained_mom_err<1.0d-5
      cfd_is_converged=cfd%constrained_mom_err<cfd_thresh.and.cfd%constrained_mom_err>0.0_gp
      if(cfd_is_converged) write(999,'(3g14.6)') cfd%B_at
      if(cfd_is_converged) lambda_t=lambda_0
      !
      if(cfd_is_converged) cfd%constrained_mom_err=0.0_gp
      !print '(3f12.6)', cfd%m_at(:,:)
      !print *,'ref.  moments',cfd%constrained_mom_err
      !print '(3f12.6)', cfd%m_at_ref(:,:)
      !
      return
  end function cfd_is_converged
end module module_cfd
