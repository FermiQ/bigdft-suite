!> @file
!! Datatypes and associated methods relativ s to the nonlocal projectors
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining datatypes of the projectors as well as constructors and destructors
module psp_projectors
  use module_base
  use gaussians
  use locregs
  use psp_projectors_base
  implicit none

  private

  public :: projector_has_overlap,get_proj_locreg
  public :: DFT_PSP_projector_iter, DFT_PSP_projectors_iter_new
  public :: DFT_PSP_projectors_iter_next, DFT_PSP_projectors_iter_ensure
  public :: DFT_PSP_projectors_iter_apply
  public :: bounds_to_plr_limits,pregion_size
  public :: update_nlpsp

  !routines which are typical of the projector application or creation follow
contains

  !> converts the bound of the lr descriptors in local bounds of the plr locregs
  pure subroutine bounds_to_plr_limits(thatway,icoarse,plr,nl1,nl2,nl3,nu1,nu2,nu3)
    implicit none
    logical, intent(in) :: thatway !< if .true., the plr descriptors has to be filled
    !! if .false., the nl bounds are filled from the plr
    integer, intent(in) :: icoarse !<controls whether to assign coarse or fine
    !!limits. Can be 1 or 2.
    !!The 2 case cannot be doe before
    !!the case with 1 has been filled
    type(locreg_descriptors), intent(inout) :: plr !<projectors locreg
    integer, intent(inout) :: nl1,nl2,nl3,nu1,nu2,nu3 !<lower and upper bounds of locregs

    if (thatway) then
       if (icoarse==1) then !coarse limits (to be done first)
          plr%ns1=nl1
          plr%ns2=nl2
          plr%ns3=nl3

          plr%d%n1=nu1-nl1
          plr%d%n2=nu2-nl2
          plr%d%n3=nu3-nl3
       else if (icoarse == 2) then
          plr%d%nfl1=nl1-plr%ns1
          plr%d%nfl2=nl2-plr%ns2
          plr%d%nfl3=nl3-plr%ns3

          plr%d%nfu1=nu1-plr%ns1
          plr%d%nfu2=nu2-plr%ns2
          plr%d%nfu3=nu3-plr%ns3
!!$         else
!!$            stop 'WRONG icoarse'
       end if
    else
       if (icoarse==1) then !coarse limits
          nl1=plr%ns1
          nl2=plr%ns2
          nl3=plr%ns3

          nu1=plr%d%n1+plr%ns1
          nu2=plr%d%n2+plr%ns2
          nu3=plr%d%n3+plr%ns3
       else if (icoarse == 2) then
          nl1=plr%d%nfl1+plr%ns1
          nl2=plr%d%nfl2+plr%ns2
          nl3=plr%d%nfl3+plr%ns3

          nu1=plr%d%nfu1+plr%ns1
          nu2=plr%d%nfu2+plr%ns2
          nu3=plr%d%nfu3+plr%ns3
!!$         else
!!$            stop 'WRONG icoarse, false case'
       end if
    end if

  end subroutine bounds_to_plr_limits


  !> Finds the size of the smallest subbox that contains a localization region made
  !! out of atom centered spheres
  subroutine pregion_size(geocode,rxyz,radius,rmult,hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
    !use module_base, only: gp
    implicit none
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, intent(in) :: n1,n2,n3
    real(gp), intent(in) :: hx,hy,hz,rmult,radius
    real(gp), dimension(3), intent(in) :: rxyz
    integer, intent(out) :: nl1,nu1,nl2,nu2,nl3,nu3
    !Local variables
    double precision, parameter :: eps_mach=1.d-12
    !n(c) real(kind=8) :: onem
    real(gp) :: cxmax,cymax,czmax,cxmin,cymin,czmin,rad

    rad=radius*rmult
    cxmax=rxyz(1)+rad ; cxmin=rxyz(1)-rad
    cymax=rxyz(2)+rad ; cymin=rxyz(2)-rad
    czmax=rxyz(3)+rad ; czmin=rxyz(3)-rad
    !n(c) onem=1.d0-eps_mach

    nl1=ceiling(real(cxmin/hx,kind=8) - eps_mach)
    nl2=ceiling(real(cymin/hy,kind=8) - eps_mach)
    nl3=ceiling(real(czmin/hz,kind=8) - eps_mach)
    nu1=floor(real(cxmax/hx,kind=8) + eps_mach)
    nu2=floor(real(cymax/hy,kind=8) + eps_mach)
    nu3=floor(real(czmax/hz,kind=8) + eps_mach)

    !for non-free BC the projectors are not allowed to be also outside the box
    if (geocode == 'F') then
       if (nl1 < 0)   stop 'nl1: projector region outside cell'
       if (nl2 < 0)   stop 'nl2: projector region outside cell'
       if (nl3 < 0)   stop 'nl3: projector region outside cell'
       if (nu1 > n1)   stop 'nu1: projector region outside cell'
       if (nu2 > n2)   stop 'nu2: projector region outside cell'
       if (nu3 > n3)   stop 'nu3: projector region outside cell'
    else if (geocode == 'S') then
       !correct the extremes if they run outside the box
       if (nl1 < 0 .or. nu1 > n1) then
          nl1=0
          nu1=n1
       end if
       if (nl2 < 0)   stop 'nl2: projector region outside cell'
       if (nu2 > n2)   stop 'nu2: projector region outside cell'
       if (nl3 < 0 .or. nu3 > n3) then
          nl3=0
          nu3=n3
       end if
    else if (geocode == 'P') then
       !correct the extremes if they run outside the box
       if (nl1 < 0 .or. nu1 > n1) then
          nl1=0
          nu1=n1
       end if
       if (nl2 < 0 .or. nu2 > n2) then
          nl2=0
          nu2=n2
       end if
       if (nl3 < 0 .or. nu3 > n3) then
          nl3=0
          nu3=n3
       end if
    end if

  END SUBROUTINE pregion_size

  !> routine to update the PSP descriptors as soon as the localization regions
  ! are modified
  subroutine update_nlpsp(nl,nlr,lrs,Glr,lr_mask)
    use locreg_operations
    use compression
    implicit none
    integer, intent(in) :: nlr
    type(locreg_descriptors), intent(in) :: Glr
    !>logical array of the localization regions active on site
    !it is true for all the elements corresponding to localisation
    !! regions whose descriptors are calculated
    logical, dimension(nlr), intent(in) :: lr_mask
    type(locreg_descriptors), dimension(nlr), intent(in) :: lrs
    type(DFT_PSP_projectors), intent(inout) :: nl
    !local variables
    type(DFT_PSP_projector_iter) :: iter
    integer :: nbseg_dim,nkeyg_dim,ilr
    integer, dimension(:), allocatable :: nbsegs_cf,keyg_lin

    call f_routine(id='update_nlpsp')

    !find allocating dimensions for work arrays
    nbseg_dim=0
    call DFT_PSP_projectors_iter_new(iter, nl)
    do while(DFT_PSP_projectors_iter_next(iter))
       nbseg_dim=max(nbseg_dim,&
            iter%pspd%plr%wfd%nseg_c+iter%pspd%plr%wfd%nseg_f)
    end do
    nkeyg_dim=0
    do ilr=1,nlr
       nkeyg_dim=max(nkeyg_dim,lrs(ilr)%wfd%nseg_c+lrs(ilr)%wfd%nseg_f)
    end do

    !allocate the work arrays for building tolr array of structures
    nbsegs_cf=f_malloc(nbseg_dim,id='nbsegs_cf')
    keyg_lin=f_malloc(nkeyg_dim,id='keyg_lin')
    !reconstruct the projectors for any of the atoms
    call DFT_PSP_projectors_iter_new(iter, nl)
    do while(DFT_PSP_projectors_iter_next(iter))
       call free_tolr_ptr(iter%pspd%tolr)
       call f_free_ptr(iter%pspd%lut_tolr)
       if (iter%mproj > 0) then
          !then fill it again, if the locreg is demanded
          iter%pspd%nlr=nlr
          call set_wfd_to_wfd(Glr,iter%pspd%plr,&
               keyg_lin,nbsegs_cf,iter%pspd%noverlap,iter%pspd%lut_tolr,iter%pspd%tolr,lrs,lr_mask)
       end if
    end do

    call f_free(keyg_lin)
    call f_free(nbsegs_cf)

    call f_release_routine()

  end subroutine update_nlpsp

!!$  !> Calculate the scalar product with the projectors of a given set of 
!!$  !! orbitals (or support functions) given in the same localization region
!!$  subroutine calculate_cproj(ncplx_p,n_p,wfd_p,proj,&
!!$       ncplx_w,n_w,wfd_w,tolr,psi_pack,scpr,psi,pdpsi,hpdpsi)
!!$    use pseudopotentials, only: apply_hij_coeff
!!$    implicit none
!!$    integer, intent(in) :: ncplx_p !< number of complex components of the projector
!!$    integer, intent(in) :: n_p !< number of elements of the projector
!!$    integer, intent(in) :: ncplx_w !< number of complex components of the wavefunction
!!$    integer, intent(in) :: n_w !< number of complex components of the wavefunction
!!$    type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
!!$    type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
!!$    !> interaction between the wavefuntion and the psp projector
!!$    type(nlpsp_to_wfd), intent(in) :: tolr
!!$    !> components of the projectors, real and imaginary parts
!!$    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,ncplx_p,n_p), intent(in) :: proj
!!$    !> components of wavefunctions, real and imaginary parts
!!$    real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,ncplx_w,n_w), intent(in) :: psi
!!$    !> workspaces for the packing array
!!$    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w*ncplx_w), intent(inout) :: psi_pack
!!$    !> array of the scalar product between the projectors and the wavefunctions
!!$    real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(inout) :: scpr
!!$    !> array of the coefficients of the hgh projectors
!!$    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(inout) :: pdpsi
!!$    !> array of the coefficients of the hgh projectors multiplied by HGH matrix
!!$    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(inout) :: hpdpsi
!!$
!!$    !put to zero the array
!!$
!!$    call f_zero(psi_pack)
!!$
!!$    !here also the PSP application strategy can be considered
!!$    call proj_dot_psi(n_p*ncplx_p,wfd_p,proj,n_w*ncplx_w,wfd_w,psi,&
!!$         tolr%nmseg_c,tolr%nmseg_f,tolr%mask,psi_pack,scpr)
!!$
!!$    !first create the coefficients for the application of the matrix
!!$    !pdpsi = < p_i | psi >
!!$    call full_coefficients('C',ncplx_p,n_p,'N',ncplx_w,n_w,scpr,'N',pdpsi)
!!$
!!$    !then create the coefficients for the evaluation of the projector energy
!!$    !pdpsi= < psi | p_i> = conj(< p_i | psi >)
!!$    call full_coefficients('N',ncplx_p,n_p,'C',ncplx_w,n_w,scpr,'C',pdpsi)
!!$
!!$  end subroutine calculate_cproj


!!$    !call to_zero(max(ncplx_w,ncplx_p)*n_w*n_p,pdpsi(1,1,1))
!!$    pdpsi=0.0_wp

  !> find the locreg that is associated to the given projector of atom iat
  !! for a locreg of label ilr. Shoudl the locreg not be found, the result is zero.
  function get_proj_locreg(pspd,ilr) result(iilr)
    implicit none
    integer, intent(in) :: ilr
    type(nonlocal_psp_descriptors),intent(in) :: pspd
    integer :: iilr
    !local variables
    integer :: jlr

    iilr=0
    do jlr=1, pspd%noverlap
       if (pspd%lut_tolr(jlr) == ilr) then
          iilr=jlr
          exit
       end if
    end do

  end function get_proj_locreg

  function projector_has_overlap(ilr, llr, glr, pspd) result(overlap)
    use compression, only: wfd_to_wfd_skip
    implicit none
    ! Calling arguments
    integer,intent(in) :: ilr
    type(locreg_descriptors),intent(in) :: llr, glr
    type(nonlocal_psp_descriptors),intent(in) :: pspd
    logical :: overlap
    ! Local variables
    integer :: iilr

    overlap = .false.

    !no projector on this atom
    !if(pspd%mproj == 0) return

    ! Check whether the projectors of this atom have an overlap with locreg ilr
    iilr = get_proj_locreg(pspd, ilr)
    if (iilr == 0) return
    if (wfd_to_wfd_skip(pspd%tolr(iilr))) return

    call check_overlap(llr, pspd%plr, glr, overlap)

  end function projector_has_overlap

  recursive function DFT_PSP_projectors_iter_next(iter, ilr, lr, glr) result(ok)
    use compression, only: wfd_to_wfd_skip
    implicit none
    type(DFT_PSP_projector_iter), intent(inout) :: iter
    integer, intent(in), optional :: ilr
    type(locreg_descriptors), intent(in), optional :: lr, glr
    logical :: ok

    logical :: overlap
    integer :: iilr

    ok = .false.
    nullify(iter%coeff)
    iter%iat = iter%iat + 1
    if (iter%iat > size(iter%parent%projs)) return
    
    ok = .true.
    iter%current => iter%parent%projs(iter%iat)
    iter%pspd => iter%current%region
    iter%mproj = iter%current%mproj
    if (iter%mproj == 0) ok = DFT_PSP_projectors_iter_next(iter, ilr, lr, glr)
    if (present(ilr) .and. present(lr) .and. present(glr)) then
       iilr = get_proj_locreg(iter%pspd, ilr)
       if (iilr > 0) iter%tolr => iter%pspd%tolr(iilr)
       overlap = (iilr > 0)
       if (overlap) overlap = .not. wfd_to_wfd_skip(iter%pspd%tolr(iilr))
       if (overlap) call check_overlap(lr, iter%pspd%plr, glr, overlap)
       if (.not. overlap) ok = DFT_PSP_projectors_iter_next(iter, ilr, lr, glr)
    end if
  end function DFT_PSP_projectors_iter_next

  subroutine DFT_PSP_projectors_iter_ensure(iter, kpt, idir, nwarnings, glr)
    implicit none
    type(DFT_PSP_projector_iter), intent(inout) :: iter
    real(gp), dimension(3) :: kpt
    integer, intent(in) :: idir
    integer, intent(out) :: nwarnings
    type(locreg_descriptors), intent(in), optional :: glr
    
    type(projector_coefficients), pointer :: proj
    type(atomic_projector_iter) :: a_it

    if (all(kpt == 0.0_gp)) then
       iter%ncplx = 1
    else
       iter%ncplx = 2
    end if
    
    ! Ensure that projector exists for this kpoint, or build it if not.
    nullify(proj)
    if (.not. iter%parent%on_the_fly) then
       proj => iter%current%projs
       do while (associated(proj))
          if (proj%kpt(1) == kpt(1) .and. proj%kpt(2) == kpt(2) .and. proj%kpt(3) == kpt(3)) exit
          proj => proj%next
       end do
       if (.not. associated(proj)) then
          allocate(proj)
          proj%kpt = kpt
          proj%idir = idir
          nullify(proj%coeff)
          proj%next => iter%current%projs
          iter%current%projs => proj
       end if
       if (associated(proj%coeff) .and. proj%idir == idir) then
          iter%coeff => proj%coeff
          return
       end if
    end if

    ! Rebuild fallback.
    call atomic_projector_iter_new(a_it, iter%parent%pbasis(iter%iat), &
         & iter%pspd%plr, kpt)
    if (iter%parent%on_the_fly) then
       iter%coeff => iter%parent%shared_proj
    else
       if (f_err_raise(.not. associated(proj), "Non existing projector.", &
            & err_name='BIGDFT_RUNTIME_ERROR')) return
       proj%idir = idir
       if (.not. associated(proj%coeff)) proj%coeff = f_malloc_ptr(a_it%nproj * a_it%nc)
       iter%coeff => proj%coeff
    end if
    call atomic_projector_iter_set_destination(a_it, iter%coeff)
    if (PROJECTION_1D_SEPARABLE == iter%parent%method) then
       if (iter%parent%pbasis(iter%iat)%kind == PROJ_DESCRIPTION_GAUSSIAN .and. &
            & present(glr)) then
          call atomic_projector_iter_set_method(a_it, PROJECTION_1D_SEPARABLE, glr)
       else
          call atomic_projector_iter_set_method(a_it, PROJECTION_RS_COLLOCATION) 
       end if
    else
       call atomic_projector_iter_set_method(a_it, iter%parent%method)
    end if

    call atomic_projector_iter_start(a_it)
    ! Loop on shell.
    do while (atomic_projector_iter_next(a_it))
       call atomic_projector_iter_to_wavelets(a_it, idir, nwarnings)
    end do

    call atomic_projector_iter_free(a_it)
  end subroutine DFT_PSP_projectors_iter_ensure
  
  subroutine DFT_PSP_projectors_iter_apply(psp_it, psi_it, at, eproj, &
       & hcproj_in, hcproj_out, hpsi, paw)
    use module_atoms
    use module_types
    use orbitalbasis
    use pseudopotentials
    use compression
    use ao_inguess, only: lmax_ao
    implicit none
    type(DFT_PSP_projector_iter), intent(in) :: psp_it
    type(ket), intent(in) :: psi_it
    type(atoms_data), intent(in) :: at
    real(wp), intent(out) :: eproj
    real(wp), dimension(max(psp_it%ncplx, psp_it%ncplx), psi_it%n_ket, psp_it%mproj), intent(in), optional :: hcproj_in
    real(wp), dimension(max(psp_it%ncplx, psp_it%ncplx), psi_it%n_ket, psp_it%mproj), intent(out), optional :: hcproj_out
    !real(wp), dimension(psi_it%ob%orbs%npsidim_orbs), intent(inout), optional :: hpsi !the dimension of the hpsi should be specified differently
    real(wp), dimension(:), intent(inout), optional :: hpsi !the dimension of the hpsi should be specified differently
    type(paw_objects), intent(inout), optional :: paw

    logical :: usepaw
    integer :: ityp, nc, m, mm
    real(gp), dimension(3,3,4) :: hij
    type(atomic_proj_matrix) :: prj
    real(wp), dimension(:), pointer :: hpsi_ptr, spsi_ptr

    call pr_dot_psi(psp_it%ncplx, psp_it%mproj, psp_it%pspd%plr%wfd, &
         & psp_it%coeff, psi_it%ncplx, psi_it%n_ket, psi_it%lr%wfd, &
         & psi_it%phi_wvl, psp_it%tolr, psp_it%parent%wpack, &
         & psp_it%parent%scpr, psp_it%parent%cproj)

    !here the cproj can be extracted to update the density matrix for the atom iat 
    if (associated(psp_it%parent%iagamma) .and. .not. present(hcproj_in)) then
       call cproj_to_gamma(psp_it%parent%pbasis(psp_it%iat), &
            & psi_it%n_ket, psp_it%mproj, lmax_ao, max(psi_it%ncplx, psp_it%ncplx), &
            & psp_it%parent%cproj, psi_it%kwgt * psi_it%occup, &
            & psp_it%parent%iagamma(0, psp_it%iat), psp_it%parent%gamma_mmp(1,1,1,1,psi_it%ispin))
    end if
    usepaw = .false.
    if (present(paw)) usepaw = paw%usepaw
    if (usepaw) then
       ! Can be done in a better way I guess...
       mm = 1
       nc = max(psp_it%ncplx, psi_it%ncplx) * psi_it%n_ket
       do m = 1, psp_it%mproj
          paw%cprj(psp_it%iat, psi_it%iorb)%cp(1:nc, m) = &
               & psp_it%parent%cproj(mm:mm+nc-1)
          mm = mm + nc
       end do
    end if

    if (.not. present(hcproj_in)) then
       nc = max(psi_it%ncplx, psp_it%ncplx) * psi_it%n_ket
       ! Compute hcproj.
       if (.not. usepaw) then
          ityp = at%astruct%iatype(psp_it%iat)
          call hgh_hij_matrix(at%npspcode(ityp), at%psppar(0,0,ityp), hij)
          if (associated(at%gamma_targets) .and. psp_it%parent%apply_gamma_target) then
             call allocate_atomic_proj_matrix(hij, psp_it%iat, psi_it%ispin, prj, &
                  & at%gamma_targets)
          else
             call allocate_atomic_proj_matrix(hij, psp_it%iat, psi_it%ispin, prj)
          end if

          if (present(hcproj_out)) then
             call apply_hij_coeff(prj, nc, psp_it%mproj, psp_it%parent%cproj, &
                  & hcproj_out)
          else
             call apply_hij_coeff(prj, nc, psp_it%mproj, psp_it%parent%cproj, &
                  & psp_it%parent%hcproj)
          end if

          call free_atomic_proj_matrix(prj)
       else
          if (present(hcproj_out)) then
             call apply_paw_coeff(paw%paw_ij(psp_it%iat)%dij, &
                  & paw%paw_ij(psp_it%iat)%cplex_dij, nc, psp_it%mproj, &
                  & psp_it%parent%cproj, hcproj_out)
          else
             call apply_paw_coeff(paw%paw_ij(psp_it%iat)%dij, &
                  & paw%paw_ij(psp_it%iat)%cplex_dij, nc, psp_it%mproj, &
                  & psp_it%parent%cproj, psp_it%parent%hcproj)
          end if
       end if
    end if

    nullify(hpsi_ptr)
    if (present(hpsi)) hpsi_ptr => ob_ket_map(hpsi, psi_it)
    if (present(hcproj_in)) then
       call cproj_dot(psp_it%ncplx, psp_it%mproj, psi_it%ncplx, psi_it%n_ket, &
            & psp_it%parent%scpr, psp_it%parent%cproj, hcproj_in, eproj)
       if (associated(hpsi_ptr)) &
            & call cproj_pr_p_psi(hcproj_in, psp_it%ncplx, psp_it%mproj, &
            & psp_it%pspd%plr%wfd, psp_it%coeff, psi_it%ncplx, psi_it%n_ket, &
            & psi_it%lr%wfd, hpsi_ptr, psp_it%tolr, &
            & psp_it%parent%wpack, psp_it%parent%scpr)
    else if (present(hcproj_out)) then
       call cproj_dot(psp_it%ncplx, psp_it%mproj, psi_it%ncplx, psi_it%n_ket, &
            & psp_it%parent%scpr, psp_it%parent%cproj, hcproj_out, eproj)
       if (associated(hpsi_ptr)) &
            & call cproj_pr_p_psi(hcproj_out, psp_it%ncplx, psp_it%mproj, &
            & psp_it%pspd%plr%wfd, psp_it%coeff, psi_it%ncplx, psi_it%n_ket, &
            & psi_it%lr%wfd, hpsi_ptr, psp_it%tolr, &
            & psp_it%parent%wpack, psp_it%parent%scpr)
    else
       call cproj_dot(psp_it%ncplx, psp_it%mproj, psi_it%ncplx, psi_it%n_ket, &
            & psp_it%parent%scpr, psp_it%parent%cproj, psp_it%parent%hcproj, eproj)
       if (associated(hpsi_ptr)) &
            & call cproj_pr_p_psi(psp_it%parent%hcproj, psp_it%ncplx, psp_it%mproj, &
            & psp_it%pspd%plr%wfd, psp_it%coeff, psi_it%ncplx, psi_it%n_ket, &
            & psi_it%lr%wfd, hpsi_ptr, psp_it%tolr, &
            & psp_it%parent%wpack, psp_it%parent%scpr)
    end if

    if (usepaw .and. present(hpsi)) then
       ityp = at%astruct%iatype(psp_it%iat)
       nc = max(psi_it%ncplx, psp_it%ncplx) * psi_it%n_ket
       call apply_paw_coeff(at%pawtab(ityp)%sij, 1, nc, psp_it%mproj, &
            & psp_it%parent%cproj, psp_it%parent%hcproj)

       spsi_ptr => ob_ket_map(paw%spsi, psi_it)
       call cproj_pr_p_psi(psp_it%parent%hcproj, psp_it%ncplx, psp_it%mproj, &
            & psp_it%pspd%plr%wfd, psp_it%coeff, psi_it%ncplx, psi_it%n_ket, &
            & psi_it%lr%wfd, spsi_ptr, psp_it%tolr, &
            & psp_it%parent%wpack, psp_it%parent%scpr)
    end if
  end subroutine DFT_PSP_projectors_iter_apply

end module psp_projectors

!> calculate the density matrix of the system from the scalar product with the projectors
subroutine cproj_to_gamma(aproj,n_w,mproj,lmax,ncplx,cproj,factor,iagamma,gamma_mmp)
  use module_defs, only: wp
  use module_base, only: f_err_raise
  use gaussians
  use psp_projectors_base, only: atomic_projectors, PROJ_DESCRIPTION_GAUSSIAN
  implicit none
  integer, intent(in) :: mproj,ncplx,lmax,n_w
  real(wp), intent(in) :: factor
  type(atomic_projectors), intent(in) :: aproj
  integer, dimension(2*lmax+1), intent(in) :: iagamma
  real(wp), dimension(ncplx,n_w,mproj), intent(in) :: cproj
  real(wp), dimension(2*n_w,2*lmax+1,2*lmax+1,*), intent(inout) :: gamma_mmp 
  !local variables
  integer :: iproj
  type(gaussian_basis_iter) :: iter

  if (f_err_raise(aproj%kind /= PROJ_DESCRIPTION_GAUSSIAN, "Not implemented.", &
       & err_name = 'BIGDFT_RUNTIME_ERROR')) return
  
  call gaussian_iter_start(aproj%gbasis, 1, iter)

  ! Loop on shell.
  iproj=1
  do while (gaussian_iter_next_shell(aproj%gbasis, iter))
     if (iter%n ==1 .and. iagamma(iter%l)/=0) then
        call atomic_PSP_density_matrix_update('C',lmax,iter%l-1,ncplx,n_w,cproj(1,1,iproj),&
             factor,gamma_mmp(1,1,1,iagamma(iter%l)))
     end if
     iproj = iproj + (2*iter%l-1)!*ncplx
  end do

end subroutine cproj_to_gamma

!>calculate the density matrix for a atomic contribution
!!from the values of scalprod calculated in the code
subroutine atomic_PSP_density_matrix_update(transp,lmax,l,ncplx,n_w,sp,fac,gamma_mmp)
  use module_defs, only: wp,gp
  use yaml_strings
  implicit none
  !> scalprod coefficients are <p_i | psi> ('N') or <psi | p_i> ('C'),
  !! ignored if ncplx=1
  character(len=1), intent(in) :: transp
  integer, intent(in) :: ncplx,n_w
  integer, intent(in) :: lmax !< maximum value of the angular momentum considered
  integer, intent(in) :: l !<angular momentum of the density matrix, form 0 to l_max
  !> coefficients of the scalar products between projectos and orbitals
  real(gp), intent(in) :: fac !<rescaling factor
  real(wp), dimension(ncplx*n_w,2*l+1), intent(in) :: sp
  !>density matrix for this angular momenum and this spin
  real(wp), dimension(2*n_w,2*lmax+1,2*lmax+1), intent(inout) :: gamma_mmp
  !local variables
  integer :: m,mp,icplx
  real(wp) :: gamma_im
  real(wp), dimension(4) :: pauliv

  if (fac==0.0_gp) return

  if (n_w==2) then
     do m=1,2*l+1
        do mp=1,2*l+1
           !here we neglect the imaginary part of the 
           !results in the case m/=m'
           pauliv=pauli_representation(sp(1,mp),sp(1,m))
           gamma_mmp(1,m,mp)=gamma_mmp(1,m,mp)+pauliv(1)
           gamma_mmp(2,m,mp)=gamma_mmp(2,m,mp)+pauliv(2)
           gamma_mmp(3,m,mp)=gamma_mmp(3,m,mp)+pauliv(3)
           gamma_mmp(4,m,mp)=gamma_mmp(4,m,mp)+pauliv(4)
        end do
     end do
  else
     do m=1,2*l+1
        do mp=1,2*l+1
           do icplx=1,ncplx
              gamma_mmp(1,m,mp)=gamma_mmp(1,m,mp)+&
                   real(fac,wp)*sp(icplx,mp)*sp(icplx,m)
           end do
           if (ncplx==2) then
              gamma_im=real(fac,wp)*(sp(2,mp)*sp(1,m)-sp(1,mp)*sp(2,m))
              if (transp .eqv. 'N') then
                 gamma_mmp(2,m,mp)=gamma_mmp(2,m,mp)+gamma_im
              else if (transp .eqv. 'C') then
                 gamma_mmp(2,m,mp)=gamma_mmp(2,m,mp)-gamma_im
              end if
           end if
        end do
     end do
  end if

  contains

    !> for the density, where phi is conjugated
    pure function pauli_representation(psimp,psim) result(rho)
      implicit none
      real(wp), dimension(4), intent(in) :: psimp
      real(wp), dimension(4), intent(in) :: psim
      real(wp), dimension(4) :: rho
      !local variables
      real(wp), dimension(4) :: p,pp

      p=psim
      where( abs(p) < 1.e-10) p=0.0_wp

      pp=psimp
      where( abs(pp) < 1.e-10) pp=0.0_wp

      !density values
      rho(1)=pp(1)*p(1)+pp(2)*p(2)+pp(3)*p(3)+pp(4)*p(4)
      rho(2)=2.0_wp*(p(1)*pp(3)+p(2)*pp(4))
      rho(3)=2.0_wp*(p(1)*pp(4)-p(2)*pp(3)) !this seems with the opposite sign
      rho(4)=p(1)*pp(1)+p(2)*pp(2)-pp(3)*p(3)-pp(4)*p(4)

    end function pauli_representation

end subroutine atomic_PSP_density_matrix_update
