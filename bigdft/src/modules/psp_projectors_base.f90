module psp_projectors_base
  use module_base
  use gaussians
  use locregs
  use compression
  use locreg_operations
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
  end type nonlocal_psp_descriptors

  integer, parameter :: PROJ_DESCRIPTION_GAUSSIAN = 1
  
  !> Description of the atomic functions
  type, public :: atomic_projectors
     integer :: iat !< Index of the atom this structure refers to.
     real(gp), dimension(3) :: rxyz !< Position of the center.
     logical :: normalized !< .true. if projectors are normalized to one.
     integer :: kind

     ! Gaussian specifics
     real(gp) :: gau_cut !< cutting radius for the gaussian description of projectors.
     type(gaussian_basis_new) :: gbasis !< Gaussian description of the projectors.
  end type atomic_projectors
  
  !> describe the information associated to the non-local part of Pseudopotentials
  type, public :: DFT_PSP_projectors
     logical :: on_the_fly             !< strategy for projector creation
     logical :: apply_gamma_target     !< apply the target identified by the gamma_mmp value
     integer :: nproj,nprojel,natoms   !< Number of projectors and number of elements
     real(gp) :: zerovol               !< Proportion of zero components.
     type(atomic_projectors), dimension(:), pointer :: pbasis !< Projectors in their own basis.
     real(wp), dimension(:), pointer :: proj !<storage space of the projectors in wavelet basis
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
     integer :: nprojel !< Total size of proj array.
     integer :: istart !< Start shift in proj for this atom.
     real(wp), dimension(:), pointer :: proj
     real(wp), dimension(:), pointer :: proj_tmp

     ! Gaussian specific attributes.
     integer :: lmax
     type(gaussian_basis_iter) :: giter

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

  public :: PROJ_DESCRIPTION_GAUSSIAN

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
  end subroutine nullify_atomic_projectors

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
    nl%nproj=0
    nl%nprojel=0
    nl%natoms=0
    nl%zerovol=100.0_gp
    nullify(nl%pbasis)
    nullify(nl%iagamma)
    nullify(nl%gamma_mmp)
    nullify(nl%proj)
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
    implicit none
    type(atomic_projectors), intent(inout) :: ap
    call gaussian_basis_free(ap%gbasis)
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
    call f_free_ptr(nl%proj)
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

    type(gaussian_basis_iter) :: iterM
    integer :: mbseg_c, mbseg_f, mbvctr_c, mbvctr_f

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
       call plr_segs_and_vctrs(iter%lr, mbseg_c, mbseg_f, mbvctr_c, mbvctr_f)
       iter%nc = (mbvctr_c + 7 * mbvctr_f) * iter%cplx
    end if

    nullify(iter%proj)
    iter%istart = 1
    iter%method = PROJECTION_1D_SEPARABLE

    iter%nproj = 0
    if (iter%parent%kind == PROJ_DESCRIPTION_GAUSSIAN) then
       ! Specific treatment for gaussian basis set.
       call atomic_projector_iter_start(iter)

       ! Maximum number of terms for every projector.
       iter%lmax = 0
       iterM = iter%giter
       do
          if (.not. gaussian_iter_next_shell(iter%parent%gbasis, iterM)) exit
          if (iterM%ndoc > 1) iter%lmax = max(iter%lmax, iterM%l)
          iter%nproj = iter%nproj + 2 * iterM%l - 1
       end do
    else
       call f_err_throw("Unknown atomic projector kind.", &
            & err_name = 'BIGDFT_RUNTIME_ERROR')
    end if
  end subroutine atomic_projector_iter_new

  subroutine atomic_projector_iter_set_destination(iter, proj, nprojel, istart_c)
    implicit none
    type(atomic_projector_iter), intent(inout) :: iter
    integer, intent(in) :: istart_c, nprojel
    real(wp), dimension(nprojel), intent(in), target :: proj

    iter%istart  =  istart_c
    iter%nprojel =  nprojel
    iter%proj    => proj
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
    end if
  end subroutine atomic_projector_iter_set_method

  subroutine atomic_projector_iter_free(iter)
    implicit none
    type(atomic_projector_iter), intent(inout) :: iter

    if (iter%method == PROJECTION_1D_SEPARABLE) then
       if (associated(iter%proj_tmp)) call f_free_ptr(iter%proj_tmp)
       call deallocate_workarrays_projectors(iter%wpr)
    else if (iter%method == PROJECTION_RS_COLLOCATION .or. &
         & iter%method == PROJECTION_MP_COLLOCATION) then
       call deallocate_work_arrays_sumrho(iter%wcol)
       call f_free_ptr(iter%proj_tmp)
    end if
  end subroutine atomic_projector_iter_free

  subroutine atomic_projector_iter_start(iter)
    implicit none
    type(atomic_projector_iter), intent(inout) :: iter

    iter%n = -1
    iter%l = -1
    iter%mproj = 0
    iter%istart_c = iter%istart
    if (iter%parent%kind == PROJ_DESCRIPTION_GAUSSIAN) then
       call gaussian_iter_start(iter%parent%gbasis, 1, iter%giter)
    else
       call f_err_throw("Unknown atomic projector kind.", &
            & err_name = 'BIGDFT_RUNTIME_ERROR')
    end if
  end subroutine atomic_projector_iter_start

  function atomic_projector_iter_next(iter) result(next)
    implicit none
    type(atomic_projector_iter), intent(inout) :: iter
    logical :: next

    iter%istart_c = iter%istart_c + iter%nc * iter%mproj ! add the previous shift.
    iter%n = -1
    iter%l = -1
    iter%mproj = 0
    next = gaussian_iter_next_shell(iter%parent%gbasis, iter%giter)
    if (.not. next) return
    
    iter%n = iter%giter%n
    iter%l = iter%giter%l
    iter%mproj = 2 * iter%l - 1
    next = .true.
    if (iter%istart_c + iter%nc * iter%mproj > iter%nprojel+1) &
         & call f_err_throw('istart_c > nprojel+1', err_name = 'BIGDFT_RUNTIME_ERROR')
  end function atomic_projector_iter_next

  function atomic_projector_iter_wnrm2(iter, m) result(nrm2)
    implicit none
    type(atomic_projector_iter), intent(in) :: iter
    integer, intent(in) :: m
    real(wp) :: nrm2

    if (f_err_raise(m <= 0 .or. m > iter%mproj, &
         & 'm > mproj', err_name = 'BIGDFT_RUNTIME_ERROR')) return

    call wnrm_wrap(iter%cplx, iter%lr%wfd%nvctr_c, iter%lr%wfd%nvctr_f, &
         & iter%proj(iter%istart_c + (m - 1) * iter%nc), nrm2)
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
               & iter%proj(iter%istart_c:), iter%proj_tmp)
       else if (iter%method == PROJECTION_RS_COLLOCATION) then
          call gaussian_iter_to_wavelets_collocation(iter%parent%gbasis, iter%giter, &
               & ider, iter%lr, iter%parent%rxyz, iter%cplx, iter%proj(iter%istart_c:), &
               & iter%proj_tmp, iter%wcol)
       else if (iter%method == PROJECTION_MP_COLLOCATION) then
          call gaussian_iter_to_wavelets_collocation(iter%parent%gbasis, iter%giter, &
               & ider, iter%lr, iter%parent%rxyz, iter%cplx, iter%proj(iter%istart_c:), &
               & iter%proj_tmp, iter%wcol, 16)
       end if
    end if
    
    ! Check norm for each proj.
    if (ider == 0 .and. iter%parent%normalized) then
       do np = 1, iter%mproj
          !here the norm should be done with the complex components
          scpr = atomic_projector_iter_wnrm2(iter, np)
          !print '(a,3(i6),1pe14.7,2(i6))','iat,l,m,scpr',iat,l,m,scpr,idir,istart_c
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

    if (aproj%kind == PROJ_DESCRIPTION_GAUSSIAN) then
       aproj%gbasis%rxyz => rxyz
    else
       call f_err_throw("Unknown atomic projector kind.", &
            & err_name = 'BIGDFT_RUNTIME_ERROR')
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
end module psp_projectors_base
