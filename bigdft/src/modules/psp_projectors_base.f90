module psp_projectors_base
  use module_base
  use gaussians
  use locregs
  use compression
  use locreg_operations
  use pspiof_m, only: pspiof_projector_t
  implicit none

  private

  !> Non local pseudopotential descriptors
  type, public :: nonlocal_psp_descriptors
     integer :: mproj !< number of projectors for this descriptor
     integer :: nlr !< total no. localization regions potentially interacting with the psp
     type(locreg_descriptors) :: plr !< localization region descriptor of a given projector (null if nlp=0)
     type(wfd_to_wfd), dimension(:), pointer :: tolr !<maskings for the locregs, dimension noverlap
     integer,dimension(:),pointer :: lut_tolr !< lookup table for tolr, dimension noverlap
     integer :: noverlap !< number of locregs which overlap with the projectors of the given atom
     logical :: shared !< True when the proj is shared
     real(wp), dimension(:), pointer :: proj !< A subptr of %proj when shared
  end type nonlocal_psp_descriptors

  integer, parameter :: PROJ_DESCRIPTION_GAUSSIAN = 1
  integer, parameter :: PROJ_DESCRIPTION_PSPIO = 2
  
  !> Description of the atomic functions
  type, public :: atomic_projectors
     integer :: iat !< Index of the atom this structure refers to.
     real(gp), dimension(3) :: rxyz !< Position of the center.
     logical :: normalized !< .true. if projectors are normalized to one.
     integer :: kind

     ! Gaussian specifics
     real(gp) :: gau_cut !< cutting radius for the gaussian description of projectors.
     type(gaussian_basis_new) :: gbasis !< Gaussian description of the projectors.

     ! PSPIO specifics
     type(pspiof_projector_t), dimension(:), pointer :: rfuncs !< Radial projectors.
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
     integer :: nproj,nprojel,natoms   !< Number of projectors and number of elements
     real(gp) :: zerovol               !< Proportion of zero components.
     type(atomic_projectors), dimension(:), pointer :: pbasis !< Projectors in their own basis.
     type(nonlocal_psp_descriptors), dimension(:), pointer :: pspd !<descriptor per projector, of size natom
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

  public :: free_DFT_PSP_projectors
  public :: DFT_PSP_projectors_null
  public :: nonlocal_psp_descriptors_null!,free_pspd_ptr
  public :: atomic_projectors_null
  public :: psp_update_positions

  public :: PROJ_DESCRIPTION_GAUSSIAN, PROJ_DESCRIPTION_PSPIO
  public :: allocate_atomic_projectors_ptr
  public :: rfunc_basis_from_pspio

  public :: atomic_projector_iter_new, atomic_projector_iter_set_method
  public :: atomic_projector_iter_free, atomic_projector_iter_set_destination
  public :: atomic_projector_iter_start, atomic_projector_iter_next
  public :: atomic_projector_iter_to_wavelets, atomic_projector_iter_wnrm2

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
    pspd%mproj=0
    pspd%nlr=0
    call nullify_locreg_descriptors(pspd%plr)
    nullify(pspd%tolr)
    nullify(pspd%lut_tolr)
    pspd%noverlap=0
    pspd%shared = .true.
    nullify(pspd%proj)
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
    ap%gau_cut = UNINITIALIZED(ap%gau_cut)
    call nullify_gaussian_basis_new(ap%gbasis)
    nullify(ap%rfuncs)
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
    nl%method = f_enumerator_null()
    nl%nproj=0
    nl%nprojel=0
    nl%natoms=0
    nl%zerovol=100.0_gp
    nullify(nl%pbasis)
    nullify(nl%iagamma)
    nullify(nl%gamma_mmp)
    nullify(nl%shared_proj)
    nullify(nl%pspd)
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
    if (.not. pspd%shared .and. associated(pspd%proj)) call f_free_ptr(pspd%proj)
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

  subroutine deallocate_DFT_PSP_projectors(nl)
    implicit none
    type(DFT_PSP_projectors), intent(inout) :: nl

    call free_pspd_ptr(nl%pspd)
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
    else if (iter%parent%kind == PROJ_DESCRIPTION_PSPIO) then
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
    real(gp) :: scpr

    if (iter%parent%kind == PROJ_DESCRIPTION_GAUSSIAN) then
       if (iter%method == PROJECTION_1D_SEPARABLE) then
          call gaussian_iter_to_wavelets_separable(iter%parent%gbasis, iter%giter, &
               & ider, iter%glr%mesh_coarse, iter%lr, iter%parent%gau_cut, &
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
    else
       call f_err_throw("Unknown atomic projector kind.", &
            & err_name = 'BIGDFT_RUNTIME_ERROR')
    end if
    
    ! Check norm for each proj.
    if (ider == 0 .and. iter%parent%normalized) then
       do np = 1, iter%mproj
          !here the norm should be done with the complex components
          scpr = atomic_projector_iter_wnrm2(iter, np)
          !print '(a,3(i6),1pe14.7,2(i6))','iat,l,m,scpr',iter%parent%iat,iter%l,np,scpr,ider,iter%istart_c
          if (abs(1.d0-scpr) > 1.d-2) then
             if (abs(1.d0-scpr) > 1.d-1) then
                if (bigdft_mpi%iproc == 0) call yaml_warning( &
                     'Norm of the nonlocal PSP atom ' // trim(yaml_toa(iter%parent%iat)) // &
                     ' l=' // trim(yaml_toa(iter%l)) // &
                     ' m=' // trim(yaml_toa(iter%n)) // ' is ' // trim(yaml_toa(scpr)) // &
                     ' while it is supposed to be about 1.0.')
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
  
  subroutine psp_update_positions(nlpsp, rxyz)
    implicit none
    type(DFT_PSP_projectors), intent(inout) :: nlpsp
    real(gp), dimension(3, nlpsp%natoms), intent(in) :: rxyz

    integer :: iat
    
    do iat = 1, nlpsp%natoms
       call atomic_projectors_set_position(nlpsp%pbasis(iat), rxyz(:, iat))
    end do
  end subroutine psp_update_positions

  subroutine rfunc_basis_from_pspio(pspio, rfuncs)
    use pspiof_m, only: pspiof_pspdata_t, pspiof_pspdata_get_n_projectors, &
         & pspiof_pspdata_get_projector, pspiof_projector_copy, PSPIO_SUCCESS
    implicit none
    type(pspiof_pspdata_t), intent(in) :: pspio
    type(pspiof_projector_t), dimension(:), pointer :: rfuncs

    integer :: i, n

    n = pspiof_pspdata_get_n_projectors(pspio)
    allocate(rfuncs(n))
    do i = 1, n
       if (f_err_raise(pspiof_projector_copy(pspiof_pspdata_get_projector(pspio, i), &
            & rfuncs(i)) /= PSPIO_SUCCESS, "Cannot copy projector " // &
            & trim(yaml_toa(i)), err_name = 'BIGDFT_RUNTIME_ERROR')) return
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
    real(gp) :: r, v, centre(3)
    type(box_iterator) :: boxit
    integer :: ithread, nthread
    type(ylm_coefficients) :: ylm

    call f_zero(psi)
    call ylm_coefficients_new(ylm, 1, l - 1)

    centre = rxyz - [cell_r(lr%mesh_coarse, lr%ns1 + 1, 1), &
         & cell_r(lr%mesh_coarse, lr%ns2 + 1, 2), &
         & cell_r(lr%mesh_coarse, lr%ns3 + 1, 3)]
    v = sqrt(lr%mesh%volume_element)
    do while(ylm_coefficients_next_m(ylm))
       call f_zero(projector_real)
       boxit = lr%bit
       ithread=0
       nthread=1
       !$omp parallel default(shared)&
       !$omp private(ithread, r) &
       !$omp firstprivate(boxit) 
       !$ ithread=omp_get_thread_num()
       !$ nthread=omp_get_num_threads()
       call box_iter_split(boxit,nthread,ithread)
       do while(box_next_point(boxit))
          projector_real(boxit%ind) = &
               & ylm_coefficients_at(ylm, boxit, centre, r) * &
               & pspiof_projector_eval(rfunc, r) * v
       end do
       call box_iter_merge(boxit)
       !$omp end parallel
       call isf_to_daub(lr, w, projector_real, psi(1,1, ylm%m))
    end do
  end subroutine rfuncs_to_wavelets_collocation
  
end module psp_projectors_base
