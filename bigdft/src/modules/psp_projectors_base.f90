module psp_projectors_base
  use module_base
  use gaussians
  use locregs
  use compression
  use locreg_operations
  use pspiof_m, only: pspiof_projector_t
  use m_pawrad, only: pawrad_type
  use m_pawtab, only: pawtab_type
  implicit none

  private

  !> Non local pseudopotential descriptors
  type, public :: nonlocal_psp_descriptors
     integer :: nlr !< total no. localization regions potentially interacting with the psp
     type(locreg_descriptors) :: plr !< localization region descriptor of a given projector (null if nlp=0)
     type(wfd_to_wfd), dimension(:), pointer :: tolr !<maskings for the locregs, dimension noverlap
     integer,dimension(:),pointer :: lut_tolr !< lookup table for tolr, dimension noverlap
     integer :: noverlap !< number of locregs which overlap with the projectors of the given atom
  end type nonlocal_psp_descriptors

  type, public :: projector_coefficients
     real(gp), dimension(3) :: kpt
     integer :: idir
     real(wp), dimension(:), pointer :: coeff
     type(projector_coefficients), pointer :: next
  end type projector_coefficients
  
  type, public :: daubechies_projectors
     type(nonlocal_psp_descriptors) :: region
     integer :: mproj !< number of projectors per k-point
     type(projector_coefficients), pointer :: projs
  end type daubechies_projectors

  integer, parameter :: PROJ_DESCRIPTION_GAUSSIAN = 1
  integer, parameter :: PROJ_DESCRIPTION_PSPIO = 2
  integer, parameter :: PROJ_DESCRIPTION_PAW = 3
  
  !> Description of the atomic functions
  type, public :: atomic_projectors
     integer :: iat !< Index of the atom this structure refers to.
     real(gp), dimension(3) :: rxyz !< Position of the center.
     real(gp), dimension(:), pointer :: normalized !< The normalisation value.
     real(gp) :: radius, fine_radius
     integer :: kind

     ! Gaussian specifics
     type(gaussian_basis_new) :: gbasis !< Gaussian description of the projectors.

     ! PSPIO specifics
     type(pspiof_projector_t), dimension(:), pointer :: rfuncs !< Radial projectors.

     ! libPAW specifics
     type(pawrad_type), pointer :: pawrad !< Radial mesh.
     type(pawtab_type), pointer :: pawtab !< Radial projector.
  end type atomic_projectors
  
  integer, parameter :: SEPARABLE_1D=0
  integer, parameter :: RS_COLLOCATION=1
  integer, parameter :: MP_COLLOCATION=2

  type(f_enumerator), public :: PROJECTION_1D_SEPARABLE=&
       f_enumerator('SEPARABLE_1D',SEPARABLE_1D,null())
  type(f_enumerator), public :: PROJECTION_RS_COLLOCATION=&
       f_enumerator('REAL_SPACE_COLLOCATION',RS_COLLOCATION,null())
  type(f_enumerator), public :: PROJECTION_MP_COLLOCATION=&
       f_enumerator('MULTIPOLE_PRESERVING_COLLOCATION',MP_COLLOCATION,null())

  !> describe the information associated to the non-local part of Pseudopotentials
  type, public :: DFT_PSP_projectors
     logical :: on_the_fly             !< strategy for projector creation
     logical :: apply_gamma_target     !< apply the target identified by the gamma_mmp value
     type(f_enumerator) :: method                 !< Prefered projection method
     integer :: nproj,nprojel,nregions   !< Number of projectors and number of elements
     real(gp) :: zerovol               !< Proportion of zero components.
     type(atomic_projectors), dimension(:), pointer :: pbasis !< Projectors in their own basis.
     type(daubechies_projectors), dimension(:), pointer :: projs !< Projectors in their region in daubechies.
     !> array to identify the order of the which are the atoms for which the density matrix is needed
     !! array of size natom,lmax
     integer, dimension(:,:), pointer :: iagamma
     !> density matrix for the required atoms, allocated from 1 to maxval(iagamma)
     real(wp), dimension(:,:,:,:,:), pointer :: gamma_mmp
     !>workspace for packing the wavefunctions in the case of multiple projectors
     real(wp), dimension(:), pointer :: wpack
     !> scalar product of the projectors and the wavefuntions, term by term (raw data)
     real(wp), dimension(:), pointer :: scpr
     !> full data of the scalar products
     real(wp), dimension(:), pointer :: cproj
     !> same quantity after application of the hamiltonian
     real(wp), dimension(:), pointer :: hcproj
     real(wp), dimension(:), pointer :: shared_proj
  end type DFT_PSP_projectors

  type, public :: atomic_projector_iter
     type(atomic_projectors), pointer :: parent
     real(gp), dimension(3) :: kpoint
     real(gp) :: normalisation
     
     integer :: cplx !< 1 for real coeff. 2 for complex ones.
     integer :: nc !< Number of components in one projector.
     integer :: nproj !< Total number of projectors.
     
     integer :: n, l !< Quantum number and orbital moment of current shell.
     integer :: istart_c !< Starting index in proj array of current shell.
     integer :: mproj !< Number of projectors for this shell.

     type(f_enumerator) :: method
     real(wp), dimension(:), pointer :: proj !< Subptr, pointing on the current
                                             !! mproj of this shell.
     real(wp), dimension(:), pointer :: proj_root
     real(wp), dimension(:), pointer :: proj_tmp

     ! Gaussian specific attributes.
     integer :: lmax
     type(gaussian_basis_iter) :: giter

     ! Radial functions specific attributes.
     integer :: riter

     ! Method specific attributes and work arrays.
     type(locreg_descriptors), pointer :: glr
     type(workarrays_projectors) :: wpr
     type(locreg_descriptors), pointer :: lr
     type(workarr_sumrho) :: wcol
  end type atomic_projector_iter

  type, public :: DFT_PSP_projector_iter
     type(DFT_PSP_projectors), pointer :: parent
     type(daubechies_projectors), pointer :: current

     type(nonlocal_psp_descriptors), pointer :: pspd
     type(wfd_to_wfd), pointer :: tolr
     integer :: iregion
     integer :: mproj
     integer :: ncplx
     real(wp), dimension(:), pointer :: coeff
  end type DFT_PSP_projector_iter

  public :: free_DFT_PSP_projectors
  public :: DFT_PSP_projectors_null
  public :: allocate_daubechies_projectors_ptr
  public :: nonlocal_psp_descriptors_null!,free_pspd_ptr
  public :: atomic_projectors_null
  public :: psp_update_positions

  public :: PROJ_DESCRIPTION_GAUSSIAN
  public :: allocate_atomic_projectors_ptr
  public :: rfunc_basis_from_pspio
  public :: rfunc_basis_from_paw

  public :: atomic_projector_iter_new, atomic_projector_iter_set_method
  public :: atomic_projector_iter_free, atomic_projector_iter_set_destination
  public :: atomic_projector_iter_start, atomic_projector_iter_next
  public :: atomic_projector_iter_to_wavelets, atomic_projector_iter_wnrm2

  public :: DFT_PSP_projectors_iter_new

contains

  pure function nonlocal_psp_descriptors_null() result(pspd)
    implicit none
    type(nonlocal_psp_descriptors) :: pspd
    call nullify_nonlocal_psp_descriptors(pspd)
  end function nonlocal_psp_descriptors_null

  pure subroutine nullify_nonlocal_psp_descriptors(pspd)
    use module_defs, only: UNINITIALIZED
    implicit none
    type(nonlocal_psp_descriptors), intent(out) :: pspd
    pspd%nlr=0
    call nullify_locreg_descriptors(pspd%plr)
    nullify(pspd%tolr)
    nullify(pspd%lut_tolr)
    pspd%noverlap=0
  end subroutine nullify_nonlocal_psp_descriptors

  pure function atomic_projectors_null() result(ap)
    implicit none
    type(atomic_projectors) :: ap
    call nullify_atomic_projectors(ap)
  end function atomic_projectors_null

  pure subroutine nullify_atomic_projectors(ap)
    use module_defs, only: UNINITIALIZED
    implicit none
    type(atomic_projectors), intent(out) :: ap
    ap%iat=0
    call nullify_gaussian_basis_new(ap%gbasis)
    nullify(ap%rfuncs)
    nullify(ap%pawrad)
    nullify(ap%pawtab)
    nullify(ap%normalized)
  end subroutine nullify_atomic_projectors

  subroutine allocate_atomic_projectors_ptr(aps, nat)
    implicit none
    type(atomic_projectors), dimension(:), pointer :: aps
    integer, intent(in) :: nat
    !local variables
    integer :: iat

    allocate(aps(nat))
    do iat = 1, nat
       call nullify_atomic_projectors(aps(iat))
    end do
  end subroutine allocate_atomic_projectors_ptr

  pure subroutine nullify_daubechies_projectors(proj)
    implicit none
    type(daubechies_projectors), intent(out) :: proj

    call nullify_nonlocal_psp_descriptors(proj%region)
    proj%mproj = 0
    nullify(proj%projs)
  end subroutine nullify_daubechies_projectors

  subroutine allocate_daubechies_projectors_ptr(projs, nat)
    implicit none
    type(daubechies_projectors), dimension(:), pointer :: projs
    integer, intent(in) :: nat
    !local variables
    integer :: iat

    allocate(projs(nat))
    do iat = 1, nat
       call nullify_daubechies_projectors(projs(iat))
    end do
  end subroutine allocate_daubechies_projectors_ptr

  pure function DFT_PSP_projectors_null() result(nl)
    implicit none
    type(DFT_PSP_projectors) :: nl
    call nullify_DFT_PSP_projectors(nl)
  end function DFT_PSP_projectors_null

  pure subroutine nullify_DFT_PSP_projectors(nl)
    implicit none
    type(DFT_PSP_projectors), intent(out) :: nl
    nl%on_the_fly=.true.
    nl%apply_gamma_target=.false.
!!$    nl%method = f_enumerator_null()
    nl%nproj=0
    nl%nprojel=0
    nl%nregions=0
    nl%zerovol=100.0_gp
    nullify(nl%pbasis)
    nullify(nl%iagamma)
    nullify(nl%gamma_mmp)
    nullify(nl%shared_proj)
    nullify(nl%projs)
    nullify(nl%wpack)
    nullify(nl%scpr)
    nullify(nl%cproj)
    nullify(nl%hcproj)
  end subroutine nullify_DFT_PSP_projectors

  !destructors
  subroutine deallocate_nonlocal_psp_descriptors(pspd)
    implicit none
    type(nonlocal_psp_descriptors), intent(inout) :: pspd
    !local variables
    call free_tolr_ptr(pspd%tolr)
    call f_free_ptr(pspd%lut_tolr)
    call deallocate_locreg_descriptors(pspd%plr)
  end subroutine deallocate_nonlocal_psp_descriptors

  subroutine free_pspd_ptr(pspd)
    implicit none
    type(nonlocal_psp_descriptors), dimension(:), pointer :: pspd
    !local variables
    integer :: iat

    if (.not. associated(pspd)) return
    do iat=lbound(pspd,1),ubound(pspd,1)
       call deallocate_nonlocal_psp_descriptors(pspd(iat))
    end do
    deallocate(pspd)
    nullify(pspd)

  end subroutine free_pspd_ptr

  subroutine deallocate_atomic_projectors(ap)
    use pspiof_m, only: pspiof_projector_free
    implicit none
    type(atomic_projectors), intent(inout) :: ap
    integer :: i
    call gaussian_basis_free(ap%gbasis)
    if (associated(ap%rfuncs)) then
       do i = lbound(ap%rfuncs, 1), ubound(ap%rfuncs, 1)
          call pspiof_projector_free(ap%rfuncs(i))
       end do
       deallocate(ap%rfuncs)
    end if
    call f_free_ptr(ap%normalized)
  end subroutine deallocate_atomic_projectors

  subroutine free_atomic_projectors_ptr(aps)
    implicit none
    type(atomic_projectors), dimension(:), pointer :: aps
    !local variables
    integer :: iat

    if (.not. associated(aps)) return
    do iat = lbound(aps, 1), ubound(aps, 1)
       call deallocate_atomic_projectors(aps(iat))
    end do
    deallocate(aps)
    nullify(aps)
  end subroutine free_atomic_projectors_ptr

  subroutine deallocate_daubechies_projectors(projs)
    implicit none
    type(daubechies_projectors), intent(inout) :: projs

    type(projector_coefficients), pointer :: proj, doomed
    
    call deallocate_nonlocal_psp_descriptors(projs%region)
    proj => projs%projs
    do while (associated(proj))
       !print *,'shape',shape(proj%coeff)
       call f_free_ptr(proj%coeff)
       doomed => proj
       proj => proj%next
       deallocate(doomed)
       nullify(doomed)
    end do
  end subroutine deallocate_daubechies_projectors

  subroutine free_daubechies_projectors_ptr(projs)
    implicit none
    type(daubechies_projectors), dimension(:), pointer :: projs
    !local variables
    integer :: iat

    if (.not. associated(projs)) return
    do iat = lbound(projs, 1), ubound(projs, 1)
       call deallocate_daubechies_projectors(projs(iat))
    end do
    deallocate(projs)
    nullify(projs)
  end subroutine free_daubechies_projectors_ptr

  subroutine deallocate_DFT_PSP_projectors(nl)
    implicit none
    type(DFT_PSP_projectors), intent(inout) :: nl

    call free_daubechies_projectors_ptr(nl%projs)
    call free_atomic_projectors_ptr(nl%pbasis)
    call f_free_ptr(nl%iagamma)
    call f_free_ptr(nl%gamma_mmp)
    call f_free_ptr(nl%shared_proj)
    call f_free_ptr(nl%wpack)
    call f_free_ptr(nl%scpr)
    call f_free_ptr(nl%cproj)
    call f_free_ptr(nl%hcproj)
  END SUBROUTINE deallocate_DFT_PSP_projectors

  subroutine free_DFT_PSP_projectors(nl)
    implicit none
    type(DFT_PSP_projectors), intent(inout) :: nl
    call deallocate_DFT_PSP_projectors(nl)
    call nullify_DFT_PSP_projectors(nl)
  end subroutine free_DFT_PSP_projectors

  subroutine atomic_projector_iter_new(iter, aproj, lr, kpoint)
    implicit none
    type(atomic_projector_iter), intent(out) :: iter
    type(atomic_projectors), intent(in), target :: aproj
    type(locreg_descriptors), intent(in), target, optional :: lr
    real(gp), dimension(3), intent(in), optional :: kpoint

    integer :: riter

    iter%parent => aproj

    iter%kpoint = 0._gp
    iter%cplx = 1
    if (present(kpoint)) then
       iter%kpoint = kpoint
       if (kpoint(1)**2 + kpoint(2)**2 + kpoint(3)**2 /= 0.0_gp) iter%cplx = 2
    end if
    iter%nc = 0
    nullify(iter%lr)
    if (present(lr)) then
       iter%lr => lr
       iter%nc = (lr%wfd%nvctr_c + 7 * lr%wfd%nvctr_f) * iter%cplx
    end if

    nullify(iter%proj)
    nullify(iter%proj_root)
    nullify(iter%proj_tmp)
    iter%method = PROJECTION_1D_SEPARABLE

    iter%lmax = 0
    iter%nproj = 0
    call atomic_projector_iter_start(iter)
    do while (atomic_projector_iter_next(iter))
       iter%nproj = iter%nproj + iter%mproj
       iter%lmax = max(iter%lmax, iter%l)
    end do
    ! Restart the iterator after use.
    call atomic_projector_iter_start(iter)

    call nullify_workarrays_projectors(iter%wpr)
    call nullify_work_arrays_sumrho(iter%wcol)
  end subroutine atomic_projector_iter_new

  subroutine atomic_projector_iter_set_destination(iter, proj)
    implicit none
    type(atomic_projector_iter), intent(inout) :: iter
    real(wp), dimension(:), intent(in), target :: proj

    iter%proj_root => proj
  end subroutine atomic_projector_iter_set_destination

  subroutine atomic_projector_iter_set_method(iter, method, glr)
    implicit none
    type(atomic_projector_iter), intent(inout) :: iter
    type(f_enumerator), intent(in) :: method
    type(locreg_descriptors), intent(in), target, optional :: glr

    iter%method = method
    nullify(iter%proj_tmp)
    nullify(iter%glr)
    call nullify_workarrays_projectors(iter%wpr)
    if (iter%method == PROJECTION_1D_SEPARABLE) then
       if (iter%lmax > 0 .and. iter%nc > 0) then
          iter%proj_tmp = f_malloc_ptr(iter%nc * (2*iter%lmax-1), id = 'proj_tmp')
       end if
       if (present(glr)) iter%glr => glr
       if (.not. associated(iter%glr)) &
            & call f_err_throw("Missing global region for 1D seprable method.", &
            & err_name = 'BIGDFT_RUNTIME_ERROR')
       ! Workarrays for the projector creation
       call allocate_workarrays_projectors(glr%d%n1, glr%d%n2, glr%d%n3, iter%wpr)
    else if (iter%method == PROJECTION_RS_COLLOCATION .or. &
         & iter%method == PROJECTION_MP_COLLOCATION) then
       iter%proj_tmp = f_malloc_ptr(iter%lr%mesh%ndim, id = 'proj_tmp')
       call initialize_work_arrays_sumrho(iter%lr, .true., iter%wcol)
    else
       call f_err_throw("Unknown projection method.", &
            & err_name = 'BIGDFT_RUNTIME_ERROR')
    end if
  end subroutine atomic_projector_iter_set_method

  subroutine atomic_projector_iter_free(iter)
    implicit none
    type(atomic_projector_iter), intent(inout) :: iter

    if (associated(iter%proj_tmp)) call f_free_ptr(iter%proj_tmp)
    call deallocate_workarrays_projectors(iter%wpr)
    call deallocate_work_arrays_sumrho(iter%wcol)
  end subroutine atomic_projector_iter_free

  subroutine atomic_projector_iter_start(iter)
    implicit none
    type(atomic_projector_iter), intent(inout) :: iter

    iter%n = -1
    iter%l = -1
    iter%mproj = 0
    iter%istart_c = 1
    if (iter%parent%kind == PROJ_DESCRIPTION_GAUSSIAN) then
       call gaussian_iter_start(iter%parent%gbasis, 1, iter%giter)
    else if (iter%parent%kind == PROJ_DESCRIPTION_PSPIO .or. &
         & iter%parent%kind == PROJ_DESCRIPTION_PAW) then
       iter%riter = 0
    else
       call f_err_throw("Unknown atomic projector kind.", &
            & err_name = 'BIGDFT_RUNTIME_ERROR')
    end if
  end subroutine atomic_projector_iter_start

  function atomic_projector_iter_next(iter) result(next)
    use pspiof_m, only: pspiof_qn_t, pspiof_qn_get_l, pspiof_qn_get_n, &
         & pspiof_projector_get_qn
    implicit none
    type(atomic_projector_iter), intent(inout) :: iter
    type(pspiof_qn_t) :: qn
    logical :: next

    iter%istart_c = iter%istart_c + iter%nc * iter%mproj ! add the previous shift.
    iter%n = -1
    iter%l = -1
    iter%mproj = 0
    if (iter%parent%kind == PROJ_DESCRIPTION_GAUSSIAN) then
       next = gaussian_iter_next_shell(iter%parent%gbasis, iter%giter)
       if (.not. next) return
       iter%n = iter%giter%n
       iter%l = iter%giter%l
       iter%normalisation = 1._gp
       if (associated(iter%parent%normalized)) iter%normalisation = &
            & iter%parent%normalized(iter%giter%ishell)
    else if (iter%parent%kind == PROJ_DESCRIPTION_PSPIO) then
       next = (iter%riter < size(iter%parent%rfuncs))
       if (.not. next) return
       iter%riter = iter%riter + 1
       qn = pspiof_projector_get_qn(iter%parent%rfuncs(iter%riter))
       iter%l = pspiof_qn_get_l(qn) + 1
       iter%n = pspiof_qn_get_n(qn)
       if (iter%n == 0) then
          iter%n = 1
       end if
       iter%normalisation = iter%parent%normalized(iter%riter)
    else if (iter%parent%kind == PROJ_DESCRIPTION_PAW) then
       next = (iter%riter < iter%parent%pawtab%basis_size)
       if (.not. next) return
       iter%riter = iter%riter + 1
       iter%l = iter%parent%pawtab%orbitals(iter%riter) + 1
       iter%n = 1
       iter%normalisation = iter%parent%normalized(iter%riter)
    else
       call f_err_throw("Unknown atomic projector kind.", &
            & err_name = 'BIGDFT_RUNTIME_ERROR')
    end if
    iter%mproj = 2 * iter%l - 1    
    next = .true.
       
    if (associated(iter%proj_root)) then
       if (iter%istart_c + iter%nc * iter%mproj > size(iter%proj_root) + 1) &
            & call f_err_throw('istart_c > nprojel+1', err_name = 'BIGDFT_RUNTIME_ERROR')
       iter%proj => f_subptr(iter%proj_root, &
            & from = iter%istart_c, size = iter%mproj * iter%nc)
    end if
  end function atomic_projector_iter_next

  function atomic_projector_iter_wnrm2(iter, m) result(nrm2)
    implicit none
    type(atomic_projector_iter), intent(in) :: iter
    integer, intent(in) :: m
    real(wp) :: nrm2

    if (f_err_raise(m <= 0 .or. m > iter%mproj, &
         & 'm > mproj', err_name = 'BIGDFT_RUNTIME_ERROR')) return

    call wnrm_wrap(iter%cplx, iter%lr%wfd%nvctr_c, iter%lr%wfd%nvctr_f, &
         & iter%proj(1 + (m - 1) * iter%nc), nrm2)
  end function atomic_projector_iter_wnrm2

  subroutine atomic_projector_iter_to_wavelets(iter, ider, nwarnings)
    use yaml_output, only: yaml_warning
    use yaml_strings, only: yaml_toa
    use compression, only: wnrm2
    implicit none
    type(atomic_projector_iter), intent(inout) :: iter
    integer, intent(in) :: ider
    integer, intent(inout), optional :: nwarnings

    integer :: np
    real(gp) :: scpr, gau_cut

    if (iter%parent%kind == PROJ_DESCRIPTION_GAUSSIAN) then
       if (iter%method == PROJECTION_1D_SEPARABLE) then
          gau_cut = UNINITIALIZED(gau_cut)
          if (associated(iter%parent%pawtab)) gau_cut = iter%parent%pawtab%rpaw
          call gaussian_iter_to_wavelets_separable(iter%parent%gbasis, iter%giter, &
               & ider, iter%glr%mesh_coarse, iter%lr, gau_cut, &
               & iter%parent%rxyz, iter%kpoint, iter%cplx, iter%wpr, &
               & iter%proj, iter%proj_tmp)
       else if (iter%method == PROJECTION_RS_COLLOCATION) then
          call gaussian_iter_to_wavelets_collocation(iter%parent%gbasis, iter%giter, &
               & ider, iter%lr, iter%parent%rxyz, iter%cplx, iter%proj, &
               & iter%proj_tmp, iter%wcol)
       else if (iter%method == PROJECTION_MP_COLLOCATION) then
          call gaussian_iter_to_wavelets_collocation(iter%parent%gbasis, iter%giter, &
               & ider, iter%lr, iter%parent%rxyz, iter%cplx, iter%proj, &
               & iter%proj_tmp, iter%wcol, 16)
       end if
    else if (iter%parent%kind == PROJ_DESCRIPTION_PSPIO) then
       if (iter%method == PROJECTION_1D_SEPARABLE) then
          call f_err_throw("1D separable projection is not possible for PSPIO.", &
               & err_name = 'BIGDFT_RUNTIME_ERROR')
       else if (iter%method == PROJECTION_RS_COLLOCATION) then
          call rfuncs_to_wavelets_collocation(iter%parent%rfuncs(iter%riter), &
               & ider, iter%lr, iter%parent%rxyz, iter%l, iter%n, iter%cplx, &
               & iter%proj, iter%proj_tmp, iter%wcol)
       else if (iter%method == PROJECTION_MP_COLLOCATION) then
          call f_err_throw("Multipole preserving projection is not implemented for PSPIO.", &
               & err_name = 'BIGDFT_RUNTIME_ERROR')
       end if
    else if (iter%parent%kind == PROJ_DESCRIPTION_PAW) then
       if (iter%method == PROJECTION_1D_SEPARABLE) then
          call f_err_throw("1D separable projection is not possible for PAW.", &
               & err_name = 'BIGDFT_RUNTIME_ERROR')
       else if (iter%method == PROJECTION_RS_COLLOCATION) then
          call paw_to_wavelets_collocation(iter%parent%pawrad, &
               & iter%parent%pawtab%tproj(1, iter%riter), &
               & ider, iter%lr, iter%parent%rxyz, iter%parent%radius, &
               & iter%l, iter%cplx, iter%proj, iter%proj_tmp, iter%wcol)
       else if (iter%method == PROJECTION_MP_COLLOCATION) then
          call f_err_throw("Multipole preserving projection is not implemented for PAW.", &
               & err_name = 'BIGDFT_RUNTIME_ERROR')
       end if
    else
       call f_err_throw("Unknown atomic projector kind.", &
            & err_name = 'BIGDFT_RUNTIME_ERROR')
    end if
    
    ! Check norm for each proj.
    if (ider == 0) then
       do np = 1, iter%mproj
          !here the norm should be done with the complex components
          scpr = atomic_projector_iter_wnrm2(iter, np)
          if (abs(iter%normalisation-scpr) > 1.d-2) then
             if (abs(iter%normalisation-scpr) > 1.d-1) then
                if (bigdft_mpi%iproc == 0) call yaml_warning( &
                     'Norm of the nonlocal PSP atom ' // trim(yaml_toa(iter%parent%iat)) // &
                     ' l=' // trim(yaml_toa(iter%l)) // &
                     ' m=' // trim(yaml_toa(np)) // ' is ' // trim(yaml_toa(scpr)) // &
                     ' while it is supposed to be about ' // &
                     & trim(yaml_toa(iter%normalisation)) //'.')
                !stop commented for the moment
                !restore the norm of the projector
                !call wscal_wrap(mbvctr_c,mbvctr_f,1.0_gp/sqrt(scpr),proj(istart_c))
             else if (present(nwarnings)) then
                nwarnings = nwarnings + 1
             end if
          end if
       end do
    end if
  end subroutine atomic_projector_iter_to_wavelets

  subroutine atomic_projectors_set_position(aproj, rxyz)
    implicit none
    type(atomic_projectors), intent(inout) :: aproj
    real(gp), dimension(3, 1), intent(in), target :: rxyz

    aproj%rxyz = rxyz(:, 1)
    if (aproj%kind == PROJ_DESCRIPTION_GAUSSIAN) then
       aproj%gbasis%rxyz => rxyz
    end if
  end subroutine atomic_projectors_set_position
  
  subroutine psp_update_positions(nlpsp, lr, lr0, rxyz)
    implicit none
    type(DFT_PSP_projectors), intent(inout) :: nlpsp
    type(locreg_descriptors), intent(in) :: lr, lr0
    real(gp), dimension(3, nlpsp%nregions), intent(in) :: rxyz

    integer :: ireg, iseg, j0, j1, ii, i0, i1, i2, i3, n1, n2, nb1, nb2, nbuf, nseg
    integer, dimension(:), allocatable :: nbsegs_cf,keyg_lin

    nb1 = lr%d%n1
    nb2 = lr%d%n2
    n1 = lr0%d%n1
    n2 = lr0%d%n2
    nbuf = (nb1 - n1) / 2
    nseg = 0
    do ireg = 1, nlpsp%nregions
       nseg = max(nseg, nlpsp%projs(ireg)%region%plr%wfd%nseg_c+nlpsp%projs(ireg)%region%plr%wfd%nseg_f)
    end do
    nbsegs_cf = f_malloc(nseg, id = 'nbsegs_cf')
    keyg_lin = f_malloc(lr%wfd%nseg_c + lr%wfd%nseg_f, id = 'keyg_lin')
    do ireg = 1, nlpsp%nregions
       call atomic_projectors_set_position(nlpsp%pbasis(ireg), rxyz(:, ireg))
       nseg = nlpsp%projs(ireg)%region%plr%wfd%nseg_c+nlpsp%projs(ireg)%region%plr%wfd%nseg_f
       do iseg = 1, nseg
          j0=nlpsp%projs(ireg)%region%plr%wfd%keyglob(1,iseg)
          j1=nlpsp%projs(ireg)%region%plr%wfd%keyglob(2,iseg)
          ii=j0-1
          i3=ii/((n1+1)*(n2+1))
          ii=ii-i3*(n1+1)*(n2+1)
          i2=ii/(n1+1)
          i0=ii-i2*(n1+1)
          i1=i0+j1-j0
          i3=i3+nbuf
          i2=i2+nbuf
          i1=i1+nbuf
          i0=i0+nbuf
          j0=i3*((nb1+1)*(nb2+1)) + i2*(nb1+1) + i0+1
          j1=i3*((nb1+1)*(nb2+1)) + i2*(nb1+1) + i1+1
          nlpsp%projs(ireg)%region%plr%wfd%keyglob(1,iseg)=j0
          nlpsp%projs(ireg)%region%plr%wfd%keyglob(2,iseg)=j1
       end do
       nlpsp%projs(ireg)%region%plr%mesh%ndims = lr%mesh%ndims
       call f_free_ptr(nlpsp%projs(ireg)%region%lut_tolr)
       call free_tolr_ptr(nlpsp%projs(ireg)%region%tolr)
       if (nseg > 0) then
          call set_wfd_to_wfd(lr, nlpsp%projs(ireg)%region%plr, &
               & keyg_lin, nbsegs_cf, nlpsp%projs(ireg)%region%noverlap, &
               & nlpsp%projs(ireg)%region%lut_tolr, nlpsp%projs(ireg)%region%tolr)
       end if
    end do
    call f_free(keyg_lin)
    call f_free(nbsegs_cf)
  end subroutine psp_update_positions

  subroutine rfunc_basis_from_pspio(aproj, pspio)
    use pspiof_m, only: pspiof_pspdata_t, pspiof_pspdata_get_n_projectors, &
         & pspiof_pspdata_get_projector, pspiof_projector_copy, PSPIO_SUCCESS, &
         & pspiof_projector_eval
    implicit none
    type(atomic_projectors), intent(inout) :: aproj
    type(pspiof_pspdata_t), intent(in) :: pspio

    integer :: i, ii, n
    real(gp) :: r
    real(gp), parameter :: eps = 1d-4

    aproj%kind = PROJ_DESCRIPTION_PSPIO
    n = pspiof_pspdata_get_n_projectors(pspio)
    allocate(aproj%rfuncs(n))
    aproj%normalized = f_malloc_ptr(n, id = "normalized")
    do i = 1, n
       if (f_err_raise(pspiof_projector_copy(pspiof_pspdata_get_projector(pspio, i), &
            & aproj%rfuncs(i)) /= PSPIO_SUCCESS, "Cannot copy projector " // &
            & trim(yaml_toa(i)), err_name = 'BIGDFT_RUNTIME_ERROR')) return
       ! Compute normalisation.
       aproj%normalized(i) = 0._gp
       do ii = 1, int(10.d0 / eps)
          r = eps * ii
          aproj%normalized(i) = aproj%normalized(i) + &
               & (pspiof_projector_eval(aproj%rfuncs(i), r) ** 2)  * r * r * eps
       end do

    end do
  end subroutine rfunc_basis_from_pspio

  subroutine rfuncs_to_wavelets_collocation(rfunc, ider, lr, rxyz, l, n, ncplx_p, psi, &
       & projector_real, w, mp_order)
    use box
    use pspiof_m, only: pspiof_projector_eval
    use f_utils, only: f_zero
    use gaussians, only: ylm_coefficients, ylm_coefficients_new
    implicit none
    type(pspiof_projector_t), intent(in) :: rfunc
    integer, intent(in) :: ider !<direction in which to perform the derivative (0 if any)
    type(locreg_descriptors), intent(in) :: lr !<projector descriptors for wavelets representation
    real(gp), dimension(3), intent(in) :: rxyz !<center of the Gaussian
    integer, intent(in) :: l !< angular momentum
    integer, intent(in) :: n !< quantum number
    integer, intent(in) :: ncplx_p !< 2 if the projector is supposed to be complex, 1 otherwise
    real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx_p, 2 * l - 1), intent(out) :: psi
    real(f_double), dimension(lr%mesh%ndim), intent(inout) :: projector_real
    type(workarr_sumrho), intent(inout) :: w
    integer, intent(in), optional :: mp_order

    !local variables
    real(gp) :: r, v, centre(3), factor
    type(box_iterator) :: boxit
    integer :: ithread, nthread
    type(ylm_coefficients) :: ylm
    integer, dimension(2,3) :: nbox
    double precision, parameter :: radius=10.d0
    double precision, parameter :: eps_mach=1.d-12
    !$ integer, external :: omp_get_thread_num, omp_get_num_threads

    call f_zero(psi)
    call ylm_coefficients_new(ylm, 1, l - 1)

    centre = rxyz - [cell_r(lr%mesh_coarse, lr%ns1 + 1, 1), &
         & cell_r(lr%mesh_coarse, lr%ns2 + 1, 2), &
         & cell_r(lr%mesh_coarse, lr%ns3 + 1, 3)]
    v = sqrt(lr%mesh%volume_element)
    do while(ylm_coefficients_next_m(ylm))
       call f_zero(projector_real)
       if (all(lr%mesh%bc == 0)) then
          boxit = lr%bit
       else
          nbox = box_nbox_from_cutoff(lr%mesh, rxyz, radius + &
               & maxval(lr%mesh%hgrids) * eps_mach)
          boxit = box_iter(lr%mesh, nbox)
       end if
       ithread=0
       nthread=1
       !$omp parallel default(shared)&
       !$omp private(ithread, r, factor) &
       !$omp firstprivate(boxit) 
       !$ ithread=omp_get_thread_num()
       !$ nthread=omp_get_num_threads()
       call box_iter_split(boxit,nthread,ithread)
       do while(box_next_point(boxit))
          factor = ylm_coefficients_at(ylm, boxit, centre, r) * v
          projector_real(boxit%ind) = factor * pspiof_projector_eval(rfunc, r)
       end do
       call box_iter_merge(boxit)
       !$omp end parallel
       call isf_to_daub(lr, w, projector_real, psi(1,1, ylm%m))
    end do
  end subroutine rfuncs_to_wavelets_collocation

  subroutine rfunc_basis_from_paw(aproj, pawrad, pawtab)
    use m_pawrad
    use m_pawtab
    use m_paw_numeric
    implicit none
    type(atomic_projectors), intent(inout) :: aproj
    type(pawrad_type), intent(in), target :: pawrad
    type(pawtab_type), intent(in), target :: pawtab

    real(gp) :: eps, r
    integer :: i, iproj, ierr, n
    real(dp), dimension(1) :: raux
    real(gp), dimension(:), allocatable :: d2
    integer, parameter :: nsteps = 100

    aproj%pawrad => pawrad
    aproj%pawtab => pawtab
    if (pawtab%has_wvl == 0) then
       aproj%kind = PROJ_DESCRIPTION_PAW
    else
       aproj%kind = PROJ_DESCRIPTION_GAUSSIAN
       call gaussian_basis_from_paw(1, [1], aproj%rxyz, [pawtab], 1, aproj%gbasis)
    end if
    aproj%normalized = f_malloc_ptr(size(pawtab%tproj, 2), id = "normalized")
    d2 = f_malloc(pawrad%mesh_size, id = "d2")
    eps = 1.05_gp * pawtab%rpaw / real(nsteps, gp)
    do iproj = 1, size(pawtab%tproj, 2)
       n = min(size(pawrad%rad), size(pawtab%tproj, 1))
       call paw_spline(pawrad%rad, pawtab%tproj(1, iproj), n, &
            & (pawtab%tproj(2, iproj) - pawtab%tproj(1, iproj)) / (pawrad%rad(2) - pawrad%rad(1)), 0._dp, d2)
       aproj%normalized(iproj) = 0._gp
       do i = 1, nsteps
          r = (i - 1) * eps
          call paw_splint(n, pawrad%rad, pawtab%tproj(1, iproj), d2, &
               & 1, [r], raux, ierr)
          aproj%normalized(iproj) = aproj%normalized(iproj) + raux(1) * raux(1) * eps
       end do
    end do
    call f_free(d2)
  end subroutine rfunc_basis_from_paw

  subroutine paw_to_wavelets_collocation(pawrad, tproj, &
       & ider, lr, rxyz, radius, l, ncplx_p, psi, projector_real, w, mp_order)
    use box
    use m_paw_numeric
    use f_utils, only: f_zero
    use gaussians, only: ylm_coefficients, ylm_coefficients_new
    implicit none
    type(pawrad_type), intent(in) :: pawrad
    real(dp), dimension(pawrad%mesh_size), intent(in) :: tproj
    integer, intent(in) :: ider !<direction in which to perform the derivative (0 if any)
    type(locreg_descriptors), intent(in) :: lr !<projector descriptors for wavelets representation
    real(gp), dimension(3), intent(in) :: rxyz !<center of the projector
    real(gp), intent(in) :: radius !<size of the projector
    integer, intent(in) :: l !< angular momentum
    integer, intent(in) :: ncplx_p !< 2 if the projector is supposed to be complex, 1 otherwise
    real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx_p, 2 * l - 1), intent(out) :: psi
    real(f_double), dimension(lr%mesh%ndim), intent(inout) :: projector_real
    type(workarr_sumrho), intent(inout) :: w
    integer, intent(in), optional :: mp_order

    !local variables
    real(gp) :: r, v, centre(3), factor
    type(box_iterator) :: boxit
    integer :: ithread, nthread, ierr
    type(ylm_coefficients) :: ylm
    real(dp), dimension(1) :: raux
    real(gp), dimension(:), allocatable :: d2
    integer, dimension(2,3) :: nbox
    double precision, parameter :: eps_mach=1.d-12
    !$ integer, external :: omp_get_thread_num, omp_get_num_threads

    call f_zero(psi)
    call ylm_coefficients_new(ylm, 1, l - 1)
    d2 = f_malloc(pawrad%mesh_size, id = "d2")
    call paw_spline(pawrad%rad, tproj, pawrad%mesh_size, &
         & (tproj(2) - tproj(1)) / (pawrad%rad(2) - pawrad%rad(1)), 0._dp, d2)
      
    centre = rxyz - [cell_r(lr%mesh_coarse, lr%ns1 + 1, 1), &
         & cell_r(lr%mesh_coarse, lr%ns2 + 1, 2), &
         & cell_r(lr%mesh_coarse, lr%ns3 + 1, 3)]
    v = sqrt(lr%mesh%volume_element)

    do while(ylm_coefficients_next_m(ylm))
       call f_zero(projector_real)
       if (all(lr%mesh%bc == 0)) then
          boxit = lr%bit
       else
          nbox = box_nbox_from_cutoff(lr%mesh, rxyz, radius + &
               & maxval(lr%mesh%hgrids) * eps_mach)
          boxit = box_iter(lr%mesh, nbox)
       end if
       ithread=0
       nthread=1
       !$omp parallel default(shared)&
       !$omp private(ithread, r, raux, ierr, factor) &
       !$omp firstprivate(boxit) 
       !$ ithread=omp_get_thread_num()
       !$ nthread=omp_get_num_threads()
       call box_iter_split(boxit,nthread,ithread)
       do while(box_next_point(boxit))
          factor = ylm_coefficients_at(ylm, boxit, centre, r)
          r = max(r, 1e-8_gp)
          call paw_splint(pawrad%mesh_size, pawrad%rad, tproj, d2, &
               & 1, [r], raux, ierr)
          projector_real(boxit%ind) = factor * raux(1) * v / r
       end do
       call box_iter_merge(boxit)
       !$omp end parallel
       call isf_to_daub(lr, w, projector_real, psi(1,1, ylm%m))
    end do
    call f_free(d2)
  end subroutine paw_to_wavelets_collocation

  subroutine DFT_PSP_projectors_iter_new(iter, nlpsp)
    implicit none
    type(DFT_PSP_projector_iter), intent(out) :: iter
    type(DFT_PSP_projectors), intent(in), target :: nlpsp

    iter%parent => nlpsp
    nullify(iter%current)
    iter%iregion = 0
    nullify(iter%pspd)
    nullify(iter%tolr)
    iter%mproj = 0
    nullify(iter%coeff)
    iter%ncplx = 0
  end subroutine DFT_PSP_projectors_iter_new

end module psp_projectors_base
