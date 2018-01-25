!> @file
!!  Define operations over gaussian functions
!! @author
!!    Copyright (C) 2008-2013 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Modules which contains the handling of operations in Gaussian basis
!! Spherical harmonics are used in the cartesian form
module gaussians
  use module_base
  implicit none
  private

  integer, parameter :: NSD_=2,EXPO_=1,COEFF_=2    !< Positions of exponents and coefficients in the storage space
  integer, parameter :: NSHID_=3,DOC_=1,L_=2,N_=3  !< Positions of shell identification numbers in shell id space
  integer, parameter :: NTERM_MAX_OVERLAP=62       !< Maximum number of terms for the considered shells
  integer, parameter :: NTERM_MAX_KINETIC=190      !< Maximum number of terms for the considered shells in the case of laplacian
  integer, parameter :: L_MAX=3                    !< Maximum number of angular momentum considered
  integer, parameter :: NTERM_MAX_RS=NTERM_MAX_KINETIC !< Max terms for the real space expression of the Gaussians

  integer, parameter :: SEPARABLE_1D=0
  integer, parameter :: SEPARABLE_COLLOCATION=1
  integer, parameter :: RADIAL_COLLOCATION=2
  integer, parameter :: MULTIPOLE_PRESERVING_COLLOCATION=3

  type(f_enumerator), public :: PROJECTION_1D_SEPARABLE=f_enumerator('SEPARABLE_1D',SEPARABLE_1D,null())
  type(f_enumerator), public :: PROJECTION_RS_COLLOCATION=f_enumerator('REAL_SPACE_COLLOCATION',SEPARABLE_COLLOCATION,null())
  type(f_enumerator), public :: PROJECTION_MP_COLLOCATION=&
       f_enumerator('MULTIPOLE_PRESERVING_COLLOCATION',MULTIPOLE_PRESERVING_COLLOCATION,null())


  !> Structures of basis of gaussian functions
  type, public :: gaussian_basis
     integer :: nat     !< Number of centers
     integer :: ncoeff  !< Number of total basis elements
     integer :: nshltot !< Total number of shells (m quantum number ignored) 
     integer :: nexpo   !< Number of exponents (sum of the contractions)
     integer :: ncplx   !< Number of complex comp. (real or complex gaussians)
     !storage units
     integer, dimension(:), pointer :: nshell      !< Number of shells for any of the centers
     integer, dimension(:), pointer :: ndoc,nam    !< Degree of contraction, angular momentum of any shell
     real(gp), dimension(:,:), pointer :: xp,psiat !< Factors and values of the exponents (complex numbers are allowed)
     real(gp), dimension(:,:), pointer :: rxyz     !< Positions of the centers
  end type gaussian_basis

  !> Structures of basis of gaussian functions
  type, public :: gaussian_basis_new
     integer :: nat     !< Number of centers
     integer :: ncoeff  !< Number of total basis elements
     integer :: nshltot !< Total number of shells (m quantum number ignored) 
     integer :: nexpo   !< Number of exponents (sum of the contractions)
     integer :: ncplx   !< =1 if traditional, =2 if complex gaussians
     !storage units
     integer, dimension(:), pointer :: nshell  !< Number of shells for any of the centers
     integer, dimension(:,:), pointer :: shid  !< Degree of contraction, angular momentum and principal quantum number
     real(gp), dimension(:,:), pointer :: sd   !< Sigma and contraction coefficients the exponents (complex numbers are allowed)
     real(gp), dimension(:,:), pointer :: rxyz !< Positions of the centers
  end type gaussian_basis_new

  !>single-center gaussians for conversion in real space
  type, public :: gaussian_real_space
     integer :: discretization_method=RADIAL_COLLOCATION !<should be either RADIAL_COLLOCATION or MULTIPOLE_PRESERVING_COLLOCATION
     integer :: nterms=0
     integer :: mp_isf=0 !< family of the interpolating scaling fucntions used for the multipole preserving approach
     real(gp), dimension(NTERM_MAX_RS) :: factors=0.0_gp
     integer, dimension(3,NTERM_MAX_RS) :: lxyz=0
     integer, dimension(NTERM_MAX_RS) :: pows=0 !<powers of r
     real(gp), dimension(3) :: rxyz=0.0_gp
     real(gp) :: exponent=0.0_gp
     real(gp) :: sigma=0.0_gp
     real(gp) :: cutoff=0.0_gp
     real(dp), dimension(:), pointer :: mpx=>null(),mpy=>null(),mpz=>null()
  end type gaussian_real_space

  public :: gaudim_check,normalize_shell,gaussian_overlap,kinetic_overlap,gauint0,overlap,kinetic

  public :: nullify_gaussian_basis, deallocate_gwf, gaussian_basis_null, gaussian_basis_free

  public :: gaussian_basis_from_psp, gaussian_basis_from_paw,nullify_gaussian_basis_new,overlap_gain,kinetic_gain

  type, public :: gaussian_basis_iter
     integer :: nshell = 0 !< Number of shells to iter on, read only.
     integer :: ishell = 0 !< Current shell id, read only.
     integer :: ndoc       !< Number of gaussians for current shell, read only.
     integer :: idoc       !< Current gaussian id, read only.
     integer :: l, n       !< Quantum number of current shell, read only.

     integer :: ishell_s   !< Internal, may change.
     integer :: iexpo      !< Internal, may change.
  end type gaussian_basis_iter
  public :: gaussian_iter_start, gaussian_iter_next_shell, gaussian_iter_next_gaussian,three_dimensional_density
  public :: gaussian_real_space_set,gaussian_radial_value,gaussian_to_wavelets_locreg,set_box_around_gaussian

contains

  pure subroutine nullify_gaussian_real_space(g)
    implicit none
    type(gaussian_real_space), intent(out) :: g
    g%discretization_method=RADIAL_COLLOCATION
    g%mp_isf=0
    g%nterms=0
    g%factors=0.0_gp
    g%lxyz=0
    g%pows=0
    g%rxyz=0.0_gp
    g%exponent=0.0_gp
    g%sigma=0.0_gp
    g%cutoff=0.0_gp
    nullify(g%mpx,g%mpy,g%mpz)
  end subroutine nullify_gaussian_real_space
  

  !here the different treatment of the gaussian for multipole preserving can be triggered
  !pure 
  function gaussian_radial_value(g,rxyz,bit) result(f)
    use numerics
    use box
    use multipole_preserving
    implicit none
    real(gp), dimension(3), intent(in) :: rxyz
    type(gaussian_real_space), intent(in) :: g
    type(box_iterator) :: bit
    real(gp) :: f
    !local variables
    logical :: domp
    integer :: i
    real(gp) :: tt,r2,r,val
    real(gp), dimension(3) :: vect, rclosest
    integer, dimension(3) :: itmp

    select case(g%discretization_method)
    case(RADIAL_COLLOCATION)
       !r = distance(bit%mesh,bit%rxyz,rxyz)
       !r2=g%exponent*r**2
       bit%tmp=bit%mesh%hgrids*(bit%inext-2)-rxyz-bit%oxyz
       r2=square_gd(bit%mesh,bit%tmp)*g%exponent
       !bit%tmp=closest_r(bit%mesh,bit%rxyz,rxyz)
       !r2=square_gd(bit%mesh,bit%tmp)*g%exponent
       f=safe_exp(-r2,underflow=1.e-120_f_double)
       tt=0.0_gp
       do i=1,g%nterms
          !this should be in absolute coordinates
          rclosest = closest_r(bit%mesh,bit%rxyz,rxyz)
          vect=rxyz_ortho(bit%mesh,rclosest)
          val=product(vect**g%lxyz(:,i))
          tt=tt+g%factors(i)*val
       end do
       f=f*tt
    case(MULTIPOLE_PRESERVING_COLLOCATION)
!!$       f=0.0_gp
!!$       domp=g%mp_isf>0
!!$       !this should be in absolute cooordinates
!!$       itmp=nint((closest_r(bit%mesh,bit%rxyz,rxyz)+rxyz)/bit%mesh%hgrids)
!!$       !x=hgrid*j-x0
!!$       r = distance(bit%mesh,bit%rxyz,rxyz)
!!$       r2= sum((bit%mesh%hgrids*itmp-rxyz)**2)
!!$       if (abs(r2- r**2) > 1.e-10_gp) then
!!$          print *,'error',r,r**2,r2,rxyz
!!$          print *,'inext',bit%inext,'position',bit%mesh%hgrids*([bit%i,bit%j,bit%k]-1)
!!$          print *,'other',bit%rxyz+bit%oxyz,bit%i,bit%j,bit%k
!!$          print *,'other',sum(closest_r(bit%mesh,bit%rxyz,rxyz)**2)
!!$          stop
!!$       end if

       f=get_mp_exps_product(g%nterms,bit%inext-1,g%factors)
!!$       do i=1,g%nterms
          !here the distance is calculated in a way that might be better handled       
!!$          bit%tmp(3)=mp_exp(bit%mesh%hgrids(3),rxyz(3),g%exponent,itmp(3),g%lxyz(3,i),domp)
!!$          bit%tmp(2)=mp_exp(bit%mesh%hgrids(2),rxyz(2),g%exponent,itmp(2),g%lxyz(2,i),domp)
!!$          bit%tmp(1)=mp_exp(bit%mesh%hgrids(1),rxyz(1),g%exponent,itmp(1),g%lxyz(1,i),domp)
!!$          
!!$          f=f+product(bit%tmp)*g%factors(i)
!!$       end do
    end select

  end function gaussian_radial_value

  !>expand a power of the r2 into separable components of x^lx y^ly z^lz
  !!only works for even powers (pow%2 =0)
  subroutine expand_r2_power(pow,factors,lxyz,nterms)
    implicit none
    integer, intent(in) :: pow
    integer, intent(out) :: nterms
    real(gp), dimension(NTERM_MAX_RS), intent(out) :: factors
    integer, dimension(3,NTERM_MAX_RS), intent(out) :: lxyz
    !local variables
    integer :: n,p,q
    real(gp) :: facn,facq,facpq

    n=pow/2
    nterms=0
    do p=0,n
       facn=rfac(n-p+1,n)
       do q=0,p
          facq=rfac(1,q)
          facpq=rfac(1,p-q)
          nterms=nterms+1
          factors(nterms)=facn/(facpq*facq)          
          lxyz(1,nterms)=2*(n-p)
          lxyz(2,nterms)=2*q
          lxyz(3,nterms)=2*(p-q)
       end do
    end do
  end subroutine expand_r2_power

  !>multiply together separable terms and obtain a ploynomial which is the sum of the two
  subroutine multiply_separable_terms(nA,facA,lxyzA,nB,facB,lxyzB,nC,facC,lxyzC)
    implicit none
    integer, intent(in) :: nA,nB
    real(gp), dimension(nA), intent(in) :: facA
    real(gp), dimension(nB), intent(in) :: facB
    integer, dimension(3,nA), intent(in) :: lxyzA
    integer, dimension(3,nB), intent(in) :: lxyzB
    integer, intent(out) :: nC
    real(gp), dimension(NTERM_MAX_RS), intent(out) :: facC
    integer, dimension(3,NTERM_MAX_RS), intent(out) :: lxyzC
    !local variables
    integer :: iA,iB

    nc=0
    do iA=1,nA
       do iB=1,nB
          nC=nC+1
          facC(nC)=facA(iA)*facB(iB)
          lxyzC(:,nC)=lxyzA(:,iA)+lxyzB(:,iB)
       end do
    end do
  end subroutine multiply_separable_terms

  subroutine gaussian_real_space_set(g,sigma,nterms,factors,lxyz,pows,mp_isf_order)
    implicit none
    integer, intent(in) :: nterms
    integer, intent(in) :: mp_isf_order
    real(gp), intent(in) :: sigma
    real(gp), dimension(nterms), intent(in) :: factors
    type(gaussian_real_space), intent(out) :: g    
    integer, dimension(nterms), intent(in) :: pows
    integer, dimension(3,nterms), intent(in) :: lxyz
    !local variables
    integer :: nterms_tmp,nC,iterm
    real(gp), dimension(NTERM_MAX_RS) :: factors_tmp,facC
    integer, dimension(3,NTERM_MAX_RS) :: lxyz_tmp,lxyzC
    
    call nullify_gaussian_real_space(g)
    g%sigma=sigma!sqrt(onehalf/expo)
    g%exponent=onehalf/sigma**2
    g%cutoff=10.0_gp*sigma
    g%mp_isf=mp_isf_order
    if (mp_isf_order >0) g%discretization_method=MULTIPOLE_PRESERVING_COLLOCATION
    !g%discretization_method=MULTIPOLE_PRESERVING_COLLOCATION !to do the tests

    select case(g%discretization_method)
    case(RADIAL_COLLOCATION)
       !leave untouched the specified powers
       g%nterms=nterms
       g%factors(1:nterms)=factors
       g%lxyz(:,1:nterms)=lxyz
       g%pows(1:nterms)=pows
    case(MULTIPOLE_PRESERVING_COLLOCATION)
       !here the pow is not allowed (yet) therefore combine the powers 
       !into separable terms. Assume that the powers here are all even
       g%nterms=0
       do iterm=1,nterms
          call expand_r2_power(pows(iterm),factors_tmp,lxyz_tmp,nterms_tmp)
          call multiply_separable_terms(nterms_tmp,factors_tmp,lxyz_tmp,&
               nterms,factors,lxyz,nC,facC,lxyzC)
          g%factors(g%nterms+1:g%nterms+nC)=facC(1:nC)
          g%lxyz(:,g%nterms+1:g%nterms+nC)=lxyzC(:,1:nC)
          g%nterms=g%nterms+nC
       end do
    end select
    
  end subroutine gaussian_real_space_set

  !>prepare for the discretization of the gaussian in real space
  subroutine set_box_around_gaussian(bit,g,rxyz)
    use box
    use multipole_preserving
    implicit none
    type(box_iterator) :: bit
    type(gaussian_real_space), intent(in) :: g
    real(gp), dimension(3), intent(in) :: rxyz
    !local variables
    integer, dimension(2,3) :: nbox

    bit%tmp=rxyz+bit%oxyz
    nbox=gaussian_nbox(bit%tmp,bit%mesh,g)
    !we might increase the box in each direction for multipole preserving and pass such box to the set_nbox routine
    call box_iter_set_nbox(bit,nbox=nbox)
    if (g%discretization_method==MULTIPOLE_PRESERVING_COLLOCATION)&
         call mp_gaussian_workarrays(g%nterms,nbox,g%exponent,g%lxyz,rxyz+bit%oxyz,bit%mesh%hgrids)
  end subroutine set_box_around_gaussian

  !> set the nbox necessary for the Gaussian.
  !! based on the cutoff and on the possible discretization treatment
  pure function gaussian_nbox(oxyz,mesh,gaussian) result(nbox)
    use box
    use multipole_preserving
    implicit none
    real(gp), dimension(3), intent(in) :: oxyz
    type(gaussian_real_space), intent(in) :: gaussian
    type(cell), intent(in) :: mesh
    integer, dimension(2,3) :: nbox
    !first calculate the box that the gaussian would have without MP treatment
    nbox=box_nbox_from_cutoff(mesh,oxyz,gaussian%cutoff)

    !then further enlarge such box in each dimension in the case of separable MP treatment
    if (mp_initialized()) then
       nbox(1,:)=nbox(1,:)-gaussian%mp_isf
       nbox(2,:)=nbox(2,:)+gaussian%mp_isf
    end if
  end function gaussian_nbox
  
  !> accumulate the density array
  !!warning@ array should be initialized to zero at the first call
  subroutine three_dimensional_density(boxit,gaussian,prefactor,oxyz,density)
    use numerics
    use box
    implicit none
    type(box_iterator) :: boxit
    type(gaussian_real_space), intent(in) :: gaussian
    real(gp), intent(in) :: prefactor
    real(gp), dimension(3), intent(in) :: oxyz
    real(gp), dimension(*), intent(inout) :: density
    !local variables
    integer :: ithread,nthread
    real(gp) :: r
    !$ integer, external :: omp_get_thread_num,omp_get_num_threads

    !calculate nbox
    ithread=0
    nthread=1
    call set_box_around_gaussian(boxit,gaussian,oxyz)
    !we might increase the box in each direction for multipole preserving and pass such box to the set_nbox routine
    !call box_iter_set_nbox(boxit,nbox=gaussian_nbox(oxyz,boxit%mesh,gaussian))
    !select case(gaussian%discretization_method)
    !case(RADIAL_COLLOCATION)
       !radial case
       !$omp parallel default(shared)&
       !$omp private(ithread,r) &
       !$omp firstprivate(boxit) 
       !$ ithread=omp_get_thread_num()
       !$ nthread=omp_get_num_threads()
       call box_iter_split(boxit,nthread,ithread)
       do while(box_next_point(boxit))
          !r = distance(boxit%mesh,boxit%rxyz,oxyz)
          !f(bit%i,bit%j,bit%k) =factor*eval(func,r)
          density(boxit%ind) = density(boxit%ind)+ prefactor*gaussian_radial_value(gaussian,oxyz,boxit)
       end do
       call box_iter_merge(boxit)
       !$omp end parallel
!!$      case(MULTIPOLE_PRESERVING_COLLOCATION) !to be rewritten
!!$         !multipole preserving case
!!$         isx=boxit%nbox(1,1)
!!$         isy=boxit%nbox(1,2)
!!$         isz=boxit%nbox(1,3)
!!$         iex=boxit%nbox(2,1)
!!$         iey=boxit%nbox(2,2)
!!$         iez=boxit%nbox(2,3)
!!$         ! Check whether the temporary arrays are large enough
!!$         if (iex-isx>nmpx) then
!!$            call f_err_throw('Temporary array in x direction too small')
!!$         end if
!!$         if (iey-isy>nmpy) then
!!$            call f_err_throw('Temporary array in y direction too small')
!!$         end if
!!$         if (iez-isz>nmpz) then
!!$            call f_err_throw('Temporary array in z direction too small')
!!$         end if
!!$
!!$         !$omp parallel
!!$         !$omp do
!!$         do i1=isx,iex
!!$            mpx(i1-isx) = mp_exp(hxh,rx,expo,i1,0,.true.)
!!$         end do
!!$         !$omp end do nowait
!!$         !$omp do
!!$         do i2=isy,iey
!!$            mpy(i2-isy) = mp_exp(hyh,ry,expo,i2,0,.true.)
!!$         end do
!!$         !$omp end do nowait
!!$         !$omp do
!!$         do i3=isz,iez
!!$            mpz(i3-isz) = mp_exp(hzh,rz,expo,i3,0,.true.)
!!$         end do
!!$         !$omp end do
!!$         !$omp end parallel
!!$
!!$         !$omp parallel default(shared)&
!!$         !$omp private(ithread,fz,fy,fx,xp) &
!!$         !$omp firstprivate(boxit) 
!!$
!!$         !$ ithread=omp_get_thread_num()
!!$         !$ nthread=omp_get_num_threads()
!!$         call box_iter_split(boxit,nthread,ithread)  !check the +1 in the array indices
!!$         do while(box_next_z(boxit))
!!$            fz=mpz(boxit%k-boxit%nbox(1,3)+1)
!!$            do while(box_next_y(boxit))
!!$               fy=mpy(boxit%j-boxit%nbox(1,2)+1)
!!$               do while(box_next_x(boxit))
!!$                  fx=mpx(boxit%i-boxit%nbox(1,1)+1)
!!$                  xp=fx*fy*fz
!!$                  density(boxit%ind) = density(boxit%ind) - xp*charge
!!$               end do
!!$            end do
!!$         end do
!!$         call box_iter_merge(boxit)
!!$         !$omp end parallel
    !end select

    call box_iter_expand_nbox(boxit)

  end subroutine three_dimensional_density

  !>accumulate the coefficients of the expression of a given gaussian in wavelets on the array of
  !!compressed data
  subroutine gaussian_to_wavelets_locreg(mesh,ider,&
       ncplx_g,coeff,expo,distance_cutoff,n,l,rxyz,kpoint,&
       ncplx_p,lr,wpr,psi,method)
    use box
    use compression
    use f_functions
    use locregs
    use locreg_operations
    use f_utils, only: f_zero
    use f_enums, only: toi
    implicit none  
    integer, intent(in) :: ider !<direction in which to perform the derivative (0 if any)
    integer, intent(in) :: n !<principal quantum number
    integer, intent(in) :: l !<angular momentum of the shell
    integer, intent(in) :: ncplx_g !< 1 or 2 if the gaussian factor is real or complex respectively
    integer, intent(in) :: ncplx_p !< 2 if the projector is supposed to be complex, 1 otherwise
    real(gp), intent(in) :: distance_cutoff !< 1d-distance starting from which the gaussian is assumed to be zero
    type(cell), intent(in) :: mesh !<cell structure *of the wavelet box* (coarse grid)
    type(locreg_descriptors), intent(in) :: lr !<projector descriptors for wavelets representation
    real(gp), dimension(ncplx_g), intent(in) :: coeff !<prefactor of the gaussian
    real(gp), dimension(ncplx_g), intent(in) :: expo !<exponents (1/2sigma^2 for the first element) of the gaussian (real and imaginary part)
    real(gp), dimension(3), intent(in) :: rxyz !<center of the Gaussian
    real(gp), dimension(3), intent(in) :: kpoint !<coordinate of the kpoint in reciprocal space
    type(workarrays_projectors),intent(inout) :: wpr
    !> wavelet expression, @todo: create a routine that accumulates instead of overwriting
    real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx_p,2*l-1), intent(out) :: psi
    type(f_enumerator), intent(in), optional :: method
    !local variables
    integer, parameter :: nterm_max=20 !if GTH nterm_max=4 (this value should go in a module)
    integer :: i,m,meth,iterm
    type(box_iterator) :: bit
    type(workarr_sumrho) :: w
    type(gaussian_real_space) :: g
    real(gp), dimension(3) :: oxyz
    real(gp), dimension(ncplx_g) :: sigma_and_expo
    integer, dimension(2*l-1) :: nterms
    integer, dimension(nterm_max,3,2*l-1) :: lxyz
    integer, dimension(3,nterm_max) :: lxyz_gau
    real(gp), dimension(ncplx_g,nterm_max,2*l-1) :: factors
    type(f_function), dimension(3) :: funcs
    real(f_double), dimension(:), allocatable :: projector_real
    !local variables
    integer, dimension(1), parameter :: zero_v1=[0]

    meth=SEPARABLE_1D
    if (present(method)) meth=toi(method)

    select case(meth)
    case(SEPARABLE_1D)
       !here iat is useless
       call projector(cell_geocode(mesh), -1, ider,l,n, coeff, expo, &
            distance_cutoff, rxyz,&
            0,0,0,mesh%ndims(1)-1,mesh%ndims(2)-1,mesh%ndims(3)-1, &
            mesh%hgrids(1),mesh%hgrids(2),mesh%hgrids(3),&
            kpoint(1),kpoint(2),kpoint(3), ncplx_p,ncplx_g, &
            lr%wfd%nvctr_c,lr%wfd%nvctr_f,lr%wfd%nseg_c,lr%wfd%nseg_f,&
            lr%wfd%keyvglob,lr%wfd%keyglob, &
            wpr,psi) 
    case(SEPARABLE_COLLOCATION)
       call get_projector_coeffs(ncplx_g,l,n,ider,nterm_max,coeff,expo,&
            nterms,lxyz,sigma_and_expo,factors)

       projector_real=f_malloc(lr%mesh%ndim,id='projector_real')
       call f_zero(psi)
       call initialize_work_arrays_sumrho(lr,.true.,w)
       do m=1,2*l-1
          call f_zero(projector_real)
          !call gaussian_real_space_set(g,sqrt(onehalf/expo(1)),nterms(m),factors(1,1,m),lxyz(1,1,m))
          do iterm=1,nterms(m)
             do i=1,3
                lxyz_gau(i,iterm)=lxyz(iterm,i,m)
             end do
          end do
          call gaussian_real_space_set(g,sigma_and_expo(1),nterms(m),factors(1,1,m),lxyz_gau,zero_v1,0)
          oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
          bit=box_iter(lr%mesh,origin=oxyz) !use here the real space mesh of the projector locreg
          call three_dimensional_density(bit,g,sqrt(lr%mesh%volume_element),rxyz,projector_real)
          call isf_to_daub(lr,w,projector_real,psi(:,:,m))
       end do
       !print *,'testRS:',sum(projector_real**2)
       call deallocate_work_arrays_sumrho(w)
       call f_free(projector_real)
    case(MULTIPOLE_PRESERVING_COLLOCATION)
       call get_projector_coeffs(ncplx_g,l,n,ider,nterm_max,coeff,expo,&
            nterms,lxyz,sigma_and_expo,factors)

       projector_real=f_malloc(lr%mesh%ndim,id='projector_real')
       call f_zero(psi)
       call initialize_work_arrays_sumrho(lr,.true.,w)
       do m=1,2*l-1
          call f_zero(projector_real)
          !call gaussian_real_space_set(g,sqrt(onehalf/expo(1)),nterms(m),factors(1,1,m),lxyz(1,1,m))
          do iterm=1,nterms(m)
             do i=1,3
                lxyz_gau(i,iterm)=lxyz(iterm,i,m)
             end do
          end do
          call gaussian_real_space_set(g,sigma_and_expo(1),nterms(m),factors(1,1,m),lxyz_gau,zero_v1,16) !to be customized
          oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
          bit=box_iter(lr%mesh,origin=oxyz) !use here the real space mesh of the projector locreg
          call three_dimensional_density(bit,g,sqrt(lr%mesh%volume_element),rxyz,projector_real)
          !print *,'test:',sum(projector_real**2)
          call isf_to_daub(lr,w,projector_real,psi(:,:,m))
       end do
       call deallocate_work_arrays_sumrho(w)
       call f_free(projector_real)

       !new method, still separable

!!$         !for the moment only with s projectors (l=0,n=1)
!!$         oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
!!$         oxyz=rxyz-oxyz
!!$         bit=box_iter(lr%mesh,origin=oxyz) !use here the real space mesh of the projector locreg
!!$         do m=1,2*l-1
!!$            do i=1,3
!!$               funcs(i)=f_function_new(f_gaussian,exponent=expo(1))
!!$            end do
!!$            !here we do not consider the lxyz terms yet
!!$            !take the reference functions
!!$            !print *,size(projector_real),'real',lr%mesh%ndims,&
!!$            !     lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3],&
!!$            !     lr%mesh_coarse%hgrids*[lr%ns1,lr%ns2,lr%ns3],rxyz,oxyz
!!$            call separable_3d_function(bit,funcs,factors(1,1,m)*sqrt(lr%mesh%volume_element),projector_real)
!!$         end do !not correctly written, it should be used to define the functions
!!$         
!!$         call f_zero(psi)
!!$         call initialize_work_arrays_sumrho(lr,.true.,w)
!!$         !from real space to wavelet
!!$         call isf_to_daub(lr,w,projector_real,psi)
!!$         !free work arrays
!!$         call deallocate_work_arrays_sumrho(w)
!!$         call f_free(projector_real)
    end select

  end subroutine gaussian_to_wavelets_locreg


  !> Nullify the pointers of the structure gaussian_basis
  pure subroutine nullify_gaussian_basis(G)

    implicit none
    !Arguments
    type(gaussian_basis),intent(inout) :: G 

    G%nat=0     
    G%ncoeff=0  
    G%nshltot=0 
    G%nexpo=0   
    G%ncplx=1
    nullify(G%nshell)
    nullify(G%ndoc)
    nullify(G%nam)
    nullify(G%psiat)
    nullify(G%xp)
    nullify(G%rxyz)
  END SUBROUTINE nullify_gaussian_basis


!> De-Allocate gaussian_basis type
  subroutine deallocate_gwf(G)
    implicit none
    type(gaussian_basis), intent(inout) :: G

    !normally positions should be deallocated outside
    call f_free_ptr(G%ndoc)
    call f_free_ptr(G%nam)
    call f_free_ptr(G%nshell)
    call f_free_ptr(G%psiat)
    call f_free_ptr(G%xp)

  END SUBROUTINE deallocate_gwf


  !> Nullify the pointers of the structure gaussian_basis_new
  pure function gaussian_basis_null() result(G)
    implicit none
    type(gaussian_basis_new) :: G
    call nullify_gaussian_basis_new(G)
  end function gaussian_basis_null
  pure subroutine nullify_gaussian_basis_new(G)
    implicit none
    type(gaussian_basis_new), intent(out) :: G
    G%nat=0
    G%ncoeff=0
    G%nshltot=0
    G%nexpo=0
    G%ncplx=1
    nullify(G%nshell)
    nullify(G%shid)
    nullify(G%sd)
    nullify(G%rxyz)
  end subroutine nullify_gaussian_basis_new


  subroutine init_gaussian_basis(nat,nshell,rxyz,G)
    implicit none
    integer, intent(in) :: nat
    integer, dimension(nat), intent(in) :: nshell
    real(gp), dimension(3,nat), intent(in), target :: rxyz
    type(gaussian_basis_new),intent(out) :: G
    !local variables
    character(len=*), parameter :: subname='init_gaussian_basis'
    integer :: iat

    G=gaussian_basis_null()

    G%nat=nat
    G%rxyz => rxyz

    !number of shells per atoms
    G%nshell = f_malloc_ptr(G%nat,id='G%nshell')

    G%nshltot=0
    do iat=1,nat
       G%nshell(iat)=nshell(iat)
       G%nshltot=G%nshltot+nshell(iat)
    end do

    G%shid = f_malloc_ptr((/ NSHID_, G%nshltot /),id='G%shid')
  end subroutine init_gaussian_basis

  subroutine gaussian_basis_from_psp(nat,iatyp,rxyz,psppar,ntyp,G)
    implicit none
    integer, intent(in) :: nat, ntyp
    integer, dimension(nat) :: iatyp
    real(gp), dimension(0:4,0:6,ntyp), intent(in) :: psppar
    real(gp), dimension(3,nat), intent(in), target :: rxyz
    type(gaussian_basis_new),intent(out) :: G
    !local variables
    integer, dimension(nat) :: nshell
    integer :: iat, l, i, ishell, iexpo

    call f_routine(id='gaussian_basis_from_psp')
    ! Build nshell from psppar.
    do iat = 1, nat
       nshell(iat) = 0
       do l=1,4 !generic case, also for HGHs (for GTH it will stop at l=2)
          do i=1,3 !generic case, also for HGHs (for GTH it will stop at i=2)
             if (psppar(l,i,iatyp(iat)) /= 0.0_gp) then
                nshell(iat) = nshell(iat) + 1
             end if
          end do
       end do
    end do


    call init_gaussian_basis(nat,nshell,rxyz,G)

    ! Associate values in shid.
    ishell = 1
    do iat = 1, nat
       do l=1,4 !generic case, also for HGHs (for GTH it will stop at l=2)
          do i=1,3 !generic case, also for HGHs (for GTH it will stop at i=2)
             if (psppar(l,i,iatyp(iat)) /= 0.0_gp) then
                G%shid(DOC_, ishell) = 1
                G%shid(L_, ishell) = l-1
                G%shid(N_, ishell) = i
                G%nexpo=G%nexpo+1
                G%ncoeff=G%ncoeff+2*l-1
                ishell = ishell + 1
             end if
          end do
       end do
    end do

    !allocate storage space (real exponents and coeffs for the moment)
    G%sd = f_malloc_ptr((/ G%ncplx*NSD_, G%nexpo /),id='G%sd')

    ! Associate values in sd.
    iexpo = 1
    do iat = 1, nat
       do l=1,4 !generic case, also for HGHs (for GTH it will stop at l=2)
          do i=1,3 !generic case, also for HGHs (for GTH it will stop at i=2)
             if (psppar(l,i,iatyp(iat)) /= 0.0_gp) then
                G%sd(EXPO_,iexpo)  = 0.5_gp/(psppar(l,0,iatyp(iat))**2)
                G%sd(COEFF_,iexpo) = 1._gp
                iexpo = iexpo + 1
             end if
          end do
       end do
    end do

    call f_release_routine()

  end subroutine gaussian_basis_from_psp

  !> Initialise the gaussian basis from PAW datas.
  subroutine gaussian_basis_from_paw(nat, iatyp, rxyz, pawtab, ntyp, G)
    use m_pawtab, only : pawtab_type
    
    implicit none

    !Arguments ------------------------------------
    integer, intent(in) :: nat, ntyp
    integer, dimension(nat) :: iatyp
    type(pawtab_type), dimension(ntyp), intent(in) :: pawtab
    real(gp), dimension(3,nat), intent(in), target :: rxyz
    type(gaussian_basis_new),intent(out) :: G

    !Local variables-------------------------------
    integer, dimension(nat) :: nshell
    integer :: iat, i, ib, lprev, nprev, l, n, ishell, iexpo

    ! Build nshell from pawtab.
    do iat = 1, nat
       nshell(iat) = pawtab(iatyp(iat))%basis_size
    end do
    call init_gaussian_basis(nat,nshell,rxyz,G)

    ! Store stuff.
    G%ncplx=2 !Complex gaussians

    ! Associate values in shid.
    ishell = 1
    do iat = 1, nat
       ib = 1
       l = -1
       n = -1
       do i = 1, G%nshell(iat)
          ! Look for next (l,n) tuple.
          do
             if (ib > pawtab(iatyp(iat))%lmn_size) stop "indlmn impl."
             lprev = l
             nprev = n
             l = pawtab(iatyp(iat))%indlmn(1, ib)
             n = pawtab(iatyp(iat))%indlmn(3, ib)
             ib = ib + 1
             !    write(*,*)ll,pawtab(itypat)%indlmn(2,ib),nn
             if (l /= lprev .or. n /= nprev) exit
             !    write(*,*)jb,pawtab(itypat)%indlmn(1:3,ib)
          end do

          G%shid(DOC_, ishell) = pawtab(iatyp(iat))%wvl%pngau(i)
          G%shid(L_, ishell) = l! + 1 ! 1 is added due to BigDFT convention
          G%shid(N_, ishell) = 1
          G%nexpo  = G%nexpo  + G%shid(DOC_, ishell)
          G%ncoeff = G%ncoeff + 2*G%shid(L_, ishell)-1
          ishell = ishell + 1
       end do
    end do
    !allocate storage space (real exponents and coeffs for the moment)
    G%sd = f_malloc_ptr((/ G%ncplx*NSD_, G%nexpo /),id='G%sd')
    ! Copy coefficients and factors.
    iexpo = 0
    do iat = 1, nat
       do i = 1, pawtab(iatyp(iat))%wvl%ptotgau
          ! real part of exponent is positive for BigDFT convention.
          G%sd((EXPO_ - 1) * 2 + 1,iexpo + i)  = abs(pawtab(iatyp(iat))%wvl%parg(1,i))
          G%sd((EXPO_ - 1) * 2 + 2,iexpo + i)  = pawtab(iatyp(iat))%wvl%parg(2,i)
          G%sd((COEFF_ - 1) * 2 + 1,iexpo + i) = pawtab(iatyp(iat))%wvl%pfac(1,i)
          G%sd((COEFF_ - 1) * 2 + 2,iexpo + i) = pawtab(iatyp(iat))%wvl%pfac(2,i)
!!$          write(*,*) iat, iexpo, i, G%sd((EXPO_ - 1) * 2 + 1,iexpo + i), &
!!$               & G%sd((EXPO_ - 1) * 2 + 2,iexpo + i), G%sd((COEFF_ - 1) * 2 + 1,iexpo + i), &
!!$               & G%sd((COEFF_ - 1) * 2 + 2,iexpo + i)
       end do
       iexpo = iexpo + pawtab(iatyp(iat))%wvl%ptotgau
    end do
  end subroutine gaussian_basis_from_paw

  !> Start a new iterator on gaussian basis for the shells of a given atom.
  subroutine gaussian_iter_start(G, iat, iter)
    implicit none
    type(gaussian_basis_new), intent(in) :: G
    integer, intent(in) :: iat
    type(gaussian_basis_iter), intent(out) :: iter

    integer :: i, j

    ! Compute offset in shell array and exposant array.
    iter%ishell_s = 0
    iter%iexpo    = 0
    do i = 1, iat - 1
       do j = iter%ishell_s + 1, iter%ishell_s + G%nshell(i)
          iter%iexpo = iter%iexpo + G%shid(DOC_, j)
       end do
       iter%ishell_s = iter%ishell_s + G%nshell(i)
    end do

    iter%ishell = 0
    iter%nshell = G%nshell(iat)
  end subroutine gaussian_iter_start

  !> Go to the next shell of the current iterator.
  function gaussian_iter_next_shell(G, iter)
    implicit none
    type(gaussian_basis_new), intent(in) :: G
    type(gaussian_basis_iter), intent(inout) :: iter
    logical :: gaussian_iter_next_shell

    gaussian_iter_next_shell = .false.
    ! End of loop case
    if (iter%ishell >= iter%nshell) return

    iter%ishell = iter%ishell + 1
    iter%l      = G%shid(L_, iter%ishell_s + iter%ishell)+1
    iter%n      = G%shid(N_, iter%ishell_s + iter%ishell)
    iter%ndoc   = G%shid(DOC_, iter%ishell_s + iter%ishell)
    iter%idoc   = 0
    if (iter%ishell > 1) iter%iexpo = iter%iexpo + G%shid(DOC_, iter%ishell_s + iter%ishell - 1)
    gaussian_iter_next_shell = .true.
  end function gaussian_iter_next_shell

  !> Go to the next gaussian of the current shell.
  function gaussian_iter_next_gaussian(G, iter, coeff, expo)
    implicit none
    type(gaussian_basis_new), intent(in) :: G
    type(gaussian_basis_iter), intent(inout) :: iter
    real(gp), dimension(G%ncplx), intent(out), optional :: coeff, expo
    logical :: gaussian_iter_next_gaussian

    gaussian_iter_next_gaussian = .false.
    ! End of loop case
    if (iter%idoc >= iter%ndoc) return

    iter%idoc  = iter%idoc + 1
    coeff(1) = G%sd(G%ncplx * (COEFF_ - 1) + 1, iter%iexpo + iter%idoc)
    expo(1)  = G%sd(G%ncplx * (EXPO_ - 1) + 1, iter%iexpo + iter%idoc)
    coeff(G%ncplx) = G%sd(G%ncplx * (COEFF_ - 1) + G%ncplx, iter%iexpo + iter%idoc)
    expo(G%ncplx)  = G%sd(G%ncplx * (EXPO_ - 1) + G%ncplx, iter%iexpo + iter%idoc)
    gaussian_iter_next_gaussian = .true.
  end function gaussian_iter_next_gaussian


  subroutine gaussian_basis_convert(G,Gold)
    implicit none
    type(gaussian_basis_new), intent(out) :: G
    type(gaussian_basis), intent(in) :: Gold
    !local variables
    character(len=*), parameter :: subname='gaussian_basis_convert'
    integer :: ishell,iexpo

    !G=init_gaussian_basis(Gold%nat,Gold%nshell,Gold%rxyz)
    call init_gaussian_basis(Gold%nat,Gold%nshell,Gold%rxyz,G)
    G%ncplx=1
    G%nexpo=0
    do ishell=1,G%nshltot
       G%shid(DOC_,ishell)=Gold%ndoc(ishell)
       G%shid(L_,ishell)=Gold%nam(ishell)-1 !<traditional convention for l restored
       G%shid(N_,ishell)=1 !<old structure only has N=1
       G%nexpo=G%nexpo+Gold%ndoc(ishell)
       G%ncoeff=G%ncoeff+2*Gold%nam(ishell)-1
    end do
    !allocate storage space (real exponents and coeffs for the moment)
    G%sd = f_malloc_ptr((/ G%ncplx*NSD_, G%nexpo /),id='G%sd')

    do iexpo=1,G%nexpo
       G%sd(EXPO_,iexpo)=0.5_gp/Gold%xp(1,iexpo)**2
       G%sd(COEFF_,iexpo)=Gold%psiat(1,iexpo)
    end do

  end subroutine gaussian_basis_convert

  subroutine gaussian_basis_free(G)
    implicit none
    type(gaussian_basis_new), intent(inout) :: G

    !do not deallocate the atomic centers
    if (associated(G%rxyz)) nullify(G%rxyz)

    call f_free_ptr(G%sd)
    call f_free_ptr(G%shid)
    call f_free_ptr(G%nshell)

    G=gaussian_basis_null()

  end subroutine gaussian_basis_free

  !> Overlap matrix between two different basis structures
  subroutine gaussian_overlap(A,B,ovrlp)
    implicit none
    type(gaussian_basis), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp
    !only lower triangular part for A%ncoeff=B%ncoeff
    !local variables
    integer, parameter :: niw=18,nrw=6
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo
    integer :: ngA,ngB,lA,lB,mA,mB
    real(gp) :: dx,dy,dz
    integer, dimension(niw) :: iw
    real(gp), dimension(nrw) :: rw
    type(gaussian_basis_new) :: G,H

    call gaussian_basis_convert(G,A)
    call gaussian_basis_convert(H,B)

    call overlap(G,H,ovrlp)

    call gaussian_basis_free(G)
    call gaussian_basis_free(H)

    return

    !deprecated version, not used
    iovrlp=0
    ishell=0
    iexpo=1
    icoeff=1

    !loop on each shell (intensive calculation)
    do iat=1,A%nat
       do isat=1,A%nshell(iat)
          ishell=ishell+1
          ngA=A%ndoc(ishell)
          lA=A%nam(ishell)
          do mA=1,2*lA-1
             iovrlp=iovrlp+1

             jovrlp=0
             jshell=0
             jexpo=1
             jcoeff=1

             do jat=1,B%nat
                dx=B%rxyz(1,jat)-A%rxyz(1,iat)
                dy=B%rxyz(2,jat)-A%rxyz(2,iat)
                dz=B%rxyz(3,jat)-A%rxyz(3,iat)
                do jsat=1,B%nshell(jat)
                   jshell=jshell+1
                   ngB=B%ndoc(jshell)
                   lB=B%nam(jshell)
                   do mB=1,2*lB-1
                      jovrlp=jovrlp+1
                      !if ((jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) .or. &
                      !     A%ncoeff /= B%ncoeff ) then
                      call gbasovrlp(A%xp(1,iexpo:),A%psiat(1,iexpo:),&
                           B%xp(1,jexpo:),B%psiat(1,jexpo:),&
                           ngA,ngB,lA,mA,lB,mB,dx,dy,dz,&
                           niw,nrw,iw,rw,ovrlp(iovrlp,jovrlp))
                      !end if
                   end do
                   jexpo=jexpo+ngB
                   jcoeff=jcoeff+2*lB-1
                end do
             end do
          end do
          iexpo=iexpo+ngA
          icoeff=icoeff+2*lA-1
       end do
    end do

    call gaudim_check(iexpo,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
    call gaudim_check(jexpo,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)

  END SUBROUTINE gaussian_overlap


  !> Overlap kinetic matrix between two different basis structures
  !! the kinetic operator is applicated on the A basis structure
  subroutine kinetic_overlap(A,B,ovrlp)
    implicit none
    type(gaussian_basis), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp 
    !only lower triangular part for A%ncoeff=B%ncoeff
    !local variables
    integer, parameter :: niw=18,nrw=6
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo
    integer :: ngA,ngB,lA,lB,mA,mB
    real(gp) :: dx,dy,dz
    integer, dimension(niw) :: iw
    real(gp), dimension(nrw) :: rw
    type(gaussian_basis_new) :: G,H

    call gaussian_basis_convert(G,A)
    call gaussian_basis_convert(H,B)

    call kinetic(G,H,ovrlp)

    call gaussian_basis_free(G)
    call gaussian_basis_free(H)

    return

    iovrlp=0
    ishell=0
    iexpo=1
    icoeff=1

    !loop on each shell (intensive calculation)
    do iat=1,A%nat
       do isat=1,A%nshell(iat)
          ishell=ishell+1
          ngA=A%ndoc(ishell)
          lA=A%nam(ishell)
          do mA=1,2*lA-1
             iovrlp=iovrlp+1

             jovrlp=0
             jshell=0
             jexpo=1
             jcoeff=1

             do jat=1,B%nat
                dx=B%rxyz(1,jat)-A%rxyz(1,iat)
                dy=B%rxyz(2,jat)-A%rxyz(2,iat)
                dz=B%rxyz(3,jat)-A%rxyz(3,iat)
                do jsat=1,B%nshell(jat)
                   jshell=jshell+1
                   ngB=B%ndoc(jshell)
                   lB=B%nam(jshell)
                   do mB=1,2*lB-1
                      jovrlp=jovrlp+1
                      if (jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff .or. &
                           A%ncoeff /= B%ncoeff ) then
                         call kineticovrlp(A%xp(1,iexpo:),A%psiat(1,iexpo:),&
                              B%xp(1,jexpo:),B%psiat(1,jexpo:),&
                              ngA,ngB,lA,mA,lB,mB,dx,dy,dz,&
                              niw,nrw,iw,rw,ovrlp(iovrlp,jovrlp))
                      end if
                   end do
                   jexpo=jexpo+ngB
                   jcoeff=jcoeff+2*lB-1
                end do
             end do
          end do
          iexpo=iexpo+ngA
          icoeff=icoeff+2*lA-1
       end do
    end do

    call gaudim_check(iexpo,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
    call gaudim_check(jexpo,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)

  END SUBROUTINE kinetic_overlap


  !> Calculates the scalar product between two shells
  !! by considering only the nonzero coefficients
  !! actual building block for calculating overlap matrix
  !! inserted work arrays for calculation
  subroutine gbasovrlp(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,m1,l2,m2,dx,dy,dz,&
       niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: ng1,ng2,l1,m1,l2,m2,niw,nrw
    real(gp), intent(in) :: dx,dy,dz
    integer, dimension(niw) :: iw
    real(gp), dimension(nrw) :: rw
    real(gp), dimension(ng1), intent(in) :: expo1,coeff1
    real(gp), dimension(ng2), intent(in) :: expo2,coeff2
    real(gp), intent(out) :: ovrlp
    !local variables
    integer :: i1,i2
    real(gp) :: a1,a2,c1,c2,govrlpr

    ovrlp=0.0_gp
    do i1=1,ng1
       a1=expo1(i1)
       a1=0.5_gp/a1**2
       c1=coeff1(i1)
       do i2=1,ng2
          a2=expo2(i2)
          a2=0.5_gp/a2**2
          c2=coeff2(i2)
          call gprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,govrlpr)
          govrlpr=c1*govrlpr*c2
          !print *,c1,c2,govrlpr
          ovrlp=ovrlp+govrlpr
       end do
    end do

  END SUBROUTINE gbasovrlp


  !> Calculates the scalar product between two shells
  !! by considering only the nonzero coefficients
  !! actual building block for calculating overlap matrix
  !! inserted work arrays for calculation
  subroutine kineticovrlp(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,m1,l2,m2,dx,dy,dz,&
       niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: ng1,ng2,l1,m1,l2,m2,niw,nrw
    real(gp), intent(in) :: dx,dy,dz
    integer, dimension(niw) :: iw
    real(gp), dimension(nrw) :: rw
    real(gp), dimension(ng1), intent(in) :: expo1,coeff1
    real(gp), dimension(ng2), intent(in) :: expo2,coeff2
    real(gp), intent(out) :: ovrlp
    !local variables
    integer :: i1,i2
    real(gp) :: a1,a2,c1,c2,govrlpr

    ovrlp=0.0_gp
    do i1=1,ng1
       a1=expo1(i1)
       a1=0.5_gp/a1**2
       c1=coeff1(i1)
       do i2=1,ng2
          a2=expo2(i2)
          a2=0.5_gp/a2**2
          c2=coeff2(i2)
          call kinprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,govrlpr)
          govrlpr=c1*govrlpr*c2
          !print *,c1,c2,govrlpr
          ovrlp=ovrlp+govrlpr
       end do
    end do

  END SUBROUTINE kineticovrlp


  !> Calculates a dot product between two differents gaussians times spherical harmonics
  !! valid only for shell which belongs to different atoms, and with also dy/=0/=dx dz/=0
  !! to be rearranged when only some of them is zero
  subroutine gprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: l1,l2,m1,m2,niw,nrw 
    real(gp), intent(in) :: a1,a2,dx,dy,dz
    integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
    real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
    real(gp), intent(out) :: ovrlp
    !local variables
    integer, parameter :: nx=3
    integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
    real(gp) :: fx,fy,fz,fa,fb!,govrlp

    !calculates the number of different couples
    call calc_coeff_inguess(l1,m1,nx,n1,&
         iw(1),iw(nx+1),iw(2*nx+1),rw(1))
    call calc_coeff_inguess(l2,m2,nx,n2,&
         iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))
    ovrlp=0.0_gp
    do i2=1,n2
       qx=iw(3*nx+i2)
       qy=iw(4*nx+i2)
       qz=iw(5*nx+i2)
       fb=rw(n1+i2)
       do i1=1,n1
          px=iw(i1)
          py=iw(nx+i1)
          pz=iw(2*nx+i1)
          fa=rw(i1)

          fx=govrlp(a1,a2,dx,px,qx)
          fy=govrlp(a1,a2,dy,py,qy)
          fz=govrlp(a1,a2,dz,pz,qz)

          ovrlp=ovrlp+fa*fb*fx*fy*fz
          !print *,i1,i2,fx,fy,fz,fa,fb
       end do
    end do

  END SUBROUTINE gprod


  !> Evaluate the wavefunction for a given grid mesh
  !! @todo TO BE verified and optimized
  subroutine wavefunction(j1,j2,j3,G,h1,h2,h3,coeff,wvfnct)
    use multipole_preserving, only: scfdotf
    implicit none
    integer, intent(in) :: j1,j2,j3
    real(gp), intent(in) :: h1,h2,h3
    type(gaussian_basis_new), intent(in) :: G
    real(gp), dimension(G%ncoeff), intent(in) :: coeff
    real(gp), intent(out) :: wvfnct
    !local variables
    integer :: iovrlp,ishell,iexpo,icoeff,iat,isat,ng,l,n,m,i
    integer :: itpd,ntpdsh,ig
    real(gp) :: s,d,wfn,f,psi

    integer, dimension(2*L_MAX+1) :: ntpd
    integer, dimension(3,NTERM_MAX_OVERLAP) :: pow
    real(gp), dimension(NTERM_MAX_OVERLAP) :: ftpd
    real(gp), dimension(3) :: r

    iovrlp=0
    ishell=0
    iexpo=1
    icoeff=1

    !loop on each shell (intensive calculation)
    wvfnct=0.0_gp
    psi=0.0_gp
    do iat=1,G%nat
       r(1)=G%rxyz(1,iat)
       r(2)=G%rxyz(2,iat)
       r(3)=G%rxyz(3,iat)
       do isat=1,G%nshell(iat)
          ishell=ishell+1
          ng=G%shid(DOC_,ishell)
          l=G%shid(L_,ishell)
          n=G%shid(N_,ishell)
          call tensor_product_decomposition(n,l,ntpdsh,ntpd,pow,ftpd)
          itpd=1
          !evaluate the wavefunction shell by shell
          do m=1,2*l+1
             !here the entire array should be extracted
             iovrlp=iovrlp+1
             do ig=0,ng-1
                s=G%sd(EXPO_,ig+iexpo)
                d=G%sd(COEFF_,ig+iexpo)
                wfn=0.0_gp
                do i=0,ntpd(m)-1
                   f=  scfdotf(j1,h1,s,r(1),pow(1,i+itpd))
                   f=f*scfdotf(j2,h2,s,r(2),pow(2,i+itpd))
                   f=f*scfdotf(j3,h3,s,r(3),pow(3,i+itpd))
                   wfn=wfn+ftpd(i+itpd)*f
                end do
                psi=psi+d*wfn
             end do
             itpd=itpd+ntpd(m)
             wvfnct=wvfnct+coeff(iovrlp)*psi
          end do
          iexpo=iexpo+ng
          icoeff=icoeff+2*l+1
       end do
    end do
    call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)

  end subroutine wavefunction


  !> Overlap matrix between two different basis structures
  subroutine overlap(A,B,ovrlp)
    implicit none
    type(gaussian_basis_new), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff), intent(out) :: ovrlp 
    !local variables
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo,itpdA,itpdB,ig1,ig2
    integer :: ngA,ngB,lA,lB,mA,mB,nA,nB,ntpdshA,ntpdshB,i1,i2
    real(gp) :: overlp,integral,s1,s2,d1,d2,f
    integer, dimension(2*L_MAX+1) :: ntpdA,ntpdB
    integer, dimension(3,NTERM_MAX_OVERLAP) :: powA,powB
    real(gp), dimension(NTERM_MAX_OVERLAP) :: ftpdA,ftpdB
    real(gp), dimension(3) :: dr,rA

    iovrlp=0
    ishell=0
    iexpo=1
    icoeff=1

    !loop on each shell (intensive calculation)
    do iat=1,A%nat
       rA(1)=A%rxyz(1,iat)
       rA(2)=A%rxyz(2,iat)
       rA(3)=A%rxyz(3,iat)
       do isat=1,A%nshell(iat)
          ishell=ishell+1
          ngA=A%shid(DOC_,ishell)
          lA=A%shid(L_,ishell)
          nA=A%shid(N_,ishell)
          call tensor_product_decomposition(nA,lA,ntpdshA,ntpdA,powA,ftpdA)
          itpdA=1
          do mA=1,2*lA+1
             !here the entire array should be extracted
             iovrlp=iovrlp+1

             jovrlp=0
             jshell=0
             jexpo=1
             jcoeff=1

             do jat=1,B%nat
                !here boundary conditions should be considered
                dr(1)=B%rxyz(1,jat)-rA(1)
                dr(2)=B%rxyz(2,jat)-rA(2)
                dr(3)=B%rxyz(3,jat)-rA(3)
                do jsat=1,B%nshell(jat)
                   jshell=jshell+1
                   ngB=B%shid(DOC_,jshell)
                   lB=B%shid(L_,jshell)
                   nB=B%shid(N_,jshell)
                   call tensor_product_decomposition(nB,lB,ntpdshB,ntpdB,powB,ftpdB)
                   itpdB=1
                   do mB=1,2*lB+1
                      jovrlp=jovrlp+1
                      !if ((jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) .or. &
                      !     A%ncoeff /= B%ncoeff ) then
!!$                      ovrlp(iovrlp,jovrlp)=&
!!$                           gdot(ngA,A%sd(1,iexpo),ntpdA(mA),powA(1,itpdA),ftpdA(itpdA),&
!!$                           ngB,B%sd(1,jexpo),ntpdB(mB),powB(1,itpdB),ftpdB(itpdB),dr)
                      !inline gdot
                      overlp=0.0_gp
                      do ig2=0,ngB-1
                         s2=B%sd(EXPO_,ig2+jexpo)
                         d2=B%sd(COEFF_,ig2+jexpo)
                         do ig1=0,ngA-1
                            s1=A%sd(EXPO_,ig1+iexpo)
                            d1=A%sd(COEFF_,ig1+iexpo)
                            integral=0.0_gp
                            do i2=0,ntpdB(mB)-1
                               do i1=0,ntpdA(mA)-1
                                  f=  govrlp(s1,s2,dr(1),powA(1,i1+itpdA),powB(1,i2+itpdB))
                                  f=f*govrlp(s1,s2,dr(2),powA(2,i1+itpdA),powB(2,i2+itpdB))
                                  f=f*govrlp(s1,s2,dr(3),powA(3,i1+itpdA),powB(3,i2+itpdB))
                                  integral=integral+ftpdA(i1+itpdA)*ftpdB(i2+itpdB)*f
                               end do
                            end do
                            overlp=overlp+d1*d2*integral
                         end do
                      end do
                      ovrlp(iovrlp,jovrlp)=overlp
                      !end if
                      itpdB=itpdB+ntpdB(mB)
                   end do
                   jexpo=jexpo+ngB
                   jcoeff=jcoeff+2*lB+1
                end do
             end do
             itpdA=itpdA+ntpdA(mA)
          end do
          iexpo=iexpo+ngA
          icoeff=icoeff+2*lA+1
       end do
    end do

    call gaudim_check(iexpo,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
    call gaudim_check(jexpo,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)

  END SUBROUTINE overlap

  !> Overlap matrix between two different basis structures
  subroutine overlap_gain(A,B,ovrlp)
    implicit none
    type(gaussian_basis_new), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff), intent(out) :: ovrlp 
    !local variables
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo,itpdA,itpdB,ig1,ig2
    integer :: ngA,ngB,lA,lB,mA,mB,nA,nB,ntpdshA,ntpdshB,i1,i2
    real(gp) :: overlp,integral,s1,s2,d1,d2,f,ss
    integer, dimension(2*L_MAX+1) :: ntpdA,ntpdB
    integer, dimension(3,NTERM_MAX_OVERLAP) :: powA,powB
    real(gp), dimension(NTERM_MAX_OVERLAP) :: ftpdA,ftpdB
    real(gp), dimension(3) :: dr,rA
    real(f_double), external :: c_overlap_c

    iovrlp=0
    ishell=0
    iexpo=1
    icoeff=1

    !loop on each shell (intensive calculation)
    do iat=1,A%nat
       rA(1)=A%rxyz(1,iat)
       rA(2)=A%rxyz(2,iat)
       rA(3)=A%rxyz(3,iat)
       do isat=1,A%nshell(iat)
          ishell=ishell+1
          ngA=A%shid(DOC_,ishell)
          lA=A%shid(L_,ishell)
          nA=A%shid(N_,ishell)
          call tensor_product_decomposition(nA,lA,ntpdshA,ntpdA,powA,ftpdA)
          itpdA=1
          do mA=1,2*lA+1
             !here the entire array should be extracted
             iovrlp=iovrlp+1

             jovrlp=0
             jshell=0
             jexpo=1
             jcoeff=1

             do jat=1,B%nat
                !here boundary conditions should be considered
                dr(1)=B%rxyz(1,jat)-rA(1)
                dr(2)=B%rxyz(2,jat)-rA(2)
                dr(3)=B%rxyz(3,jat)-rA(3)
                do jsat=1,B%nshell(jat)
                   jshell=jshell+1
                   ngB=B%shid(DOC_,jshell)
                   lB=B%shid(L_,jshell)
                   nB=B%shid(N_,jshell)
                   call tensor_product_decomposition(nB,lB,ntpdshB,ntpdB,powB,ftpdB)
                   itpdB=1
                   do mB=1,2*lB+1
                      jovrlp=jovrlp+1
                      !if ((jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) .or. &
                      !     A%ncoeff /= B%ncoeff ) then
!!$                      ovrlp(iovrlp,jovrlp)=&
!!$                           gdot(ngA,A%sd(1,iexpo),ntpdA(mA),powA(1,itpdA),ftpdA(itpdA),&
!!$                           ngB,B%sd(1,jexpo),ntpdB(mB),powB(1,itpdB),ftpdB(itpdB),dr)
                      !inline gdot
                      overlp=0.0_gp
                      do ig2=0,ngB-1
                         s2=B%sd(EXPO_,ig2+jexpo)
                         d2=B%sd(COEFF_,ig2+jexpo)
                         do ig1=0,ngA-1
                            s1=A%sd(EXPO_,ig1+iexpo)
                            d1=A%sd(COEFF_,ig1+iexpo)
                            integral=0.0_gp
                            do i2=0,ntpdB(mB)-1
                               do i1=0,ntpdA(mA)-1
                                  ss=c_overlap_c(s1,rA,powA(1,i1+itpdA),powA(2,i1+itpdA),powA(3,i1+itpdA),&
                                       s2,B%rxyz(1,jat),powB(1,i2+itpdB),powB(2,i2+itpdB),powB(3,i2+itpdB))
                                  integral=integral+ftpdA(i1+itpdA)*ftpdB(i2+itpdB)*ss
                               end do
                            end do
                            overlp=overlp+d1*d2*integral
                         end do
                      end do
                      ovrlp(iovrlp,jovrlp)=overlp
                      !end if
                      itpdB=itpdB+ntpdB(mB)
                   end do
                   jexpo=jexpo+ngB
                   jcoeff=jcoeff+2*lB+1
                end do
             end do
             itpdA=itpdA+ntpdA(mA)
          end do
          iexpo=iexpo+ngA
          icoeff=icoeff+2*lA+1
       end do
    end do

    call gaudim_check(iexpo,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
    call gaudim_check(jexpo,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)

  END SUBROUTINE overlap_gain


  !> Overlap matrix between two different basis structures
  subroutine kinetic_gain(A,B,ovrlp)
    implicit none
    type(gaussian_basis_new), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff), intent(out) :: ovrlp 
    !local variables
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo,itpdA,itpdB,ig1,ig2
    integer :: ngA,ngB,lA,lB,mA,mB,nA,nB,ntpdshA,ntpdshB,i1,i2
    real(gp) :: overlp,integral,s1,s2,d1,d2,f,ss
    integer, dimension(2*L_MAX+1) :: ntpdA,ntpdB
    integer, dimension(3,NTERM_MAX_OVERLAP) :: powA,powB
    real(gp), dimension(NTERM_MAX_OVERLAP) :: ftpdA,ftpdB
    real(gp), dimension(3) :: dr,rA
    real(f_double), external :: c_laplacian_c

    iovrlp=0
    ishell=0
    iexpo=1
    icoeff=1

    !loop on each shell (intensive calculation)
    do iat=1,A%nat
       rA(1)=A%rxyz(1,iat)
       rA(2)=A%rxyz(2,iat)
       rA(3)=A%rxyz(3,iat)
       do isat=1,A%nshell(iat)
          ishell=ishell+1
          ngA=A%shid(DOC_,ishell)
          lA=A%shid(L_,ishell)
          nA=A%shid(N_,ishell)
          call tensor_product_decomposition(nA,lA,ntpdshA,ntpdA,powA,ftpdA)
          itpdA=1
          do mA=1,2*lA+1
             !here the entire array should be extracted
             iovrlp=iovrlp+1

             jovrlp=0
             jshell=0
             jexpo=1
             jcoeff=1

             do jat=1,B%nat
                !here boundary conditions should be considered
                dr(1)=B%rxyz(1,jat)-rA(1)
                dr(2)=B%rxyz(2,jat)-rA(2)
                dr(3)=B%rxyz(3,jat)-rA(3)
                do jsat=1,B%nshell(jat)
                   jshell=jshell+1
                   ngB=B%shid(DOC_,jshell)
                   lB=B%shid(L_,jshell)
                   nB=B%shid(N_,jshell)
                   call tensor_product_decomposition(nB,lB,ntpdshB,ntpdB,powB,ftpdB)
                   itpdB=1
                   do mB=1,2*lB+1
                      jovrlp=jovrlp+1
                      !if ((jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) .or. &
                      !     A%ncoeff /= B%ncoeff ) then
!!$                      ovrlp(iovrlp,jovrlp)=&
!!$                           gdot(ngA,A%sd(1,iexpo),ntpdA(mA),powA(1,itpdA),ftpdA(itpdA),&
!!$                           ngB,B%sd(1,jexpo),ntpdB(mB),powB(1,itpdB),ftpdB(itpdB),dr)
                      !inline gdot
                      overlp=0.0_gp
                      do ig2=0,ngB-1
                         s2=B%sd(EXPO_,ig2+jexpo)
                         d2=B%sd(COEFF_,ig2+jexpo)
                         do ig1=0,ngA-1
                            s1=A%sd(EXPO_,ig1+iexpo)
                            d1=A%sd(COEFF_,ig1+iexpo)
                            integral=0.0_gp
                            do i2=0,ntpdB(mB)-1
                               do i1=0,ntpdA(mA)-1
                                  ss=c_laplacian_c(s2,B%rxyz(1,jat),powB(1,i2+itpdB),powB(2,i2+itpdB),powB(3,i2+itpdB),&
                                       s1,rA,powA(1,i1+itpdA),powA(2,i1+itpdA),powA(3,i1+itpdA))
                                  integral=integral+ftpdA(i1+itpdA)*ftpdB(i2+itpdB)*ss
                               end do
                            end do
                            overlp=overlp+d1*d2*integral
                         end do
                      end do
                      ovrlp(iovrlp,jovrlp)=overlp
                      !end if
                      itpdB=itpdB+ntpdB(mB)
                   end do
                   jexpo=jexpo+ngB
                   jcoeff=jcoeff+2*lB+1
                end do
             end do
             itpdA=itpdA+ntpdA(mA)
          end do
          iexpo=iexpo+ngA
          icoeff=icoeff+2*lA+1
       end do
    end do

    call gaudim_check(iexpo,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
    call gaudim_check(jexpo,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)

  END SUBROUTINE kinetic_gain



  !> Calculates a dot product between two basis functions
  !! Basis function is identified by its tensor product decompositions and sigmas
  !! the contraction coefficients are also given
  pure function gdot(ng1,sd1,ntpd1,pws1,ftpd1,ng2,sd2,ntpd2,pws2,ftpd2,dr) result(ovrlp)
    implicit none
    integer, intent(in) :: ng1,ng2 !< degrees of contractions of the basis
    integer, intent(in) :: ntpd1,ntpd2 !<number of tensor product decopositions
    integer, dimension(3,ntpd1), intent(in) :: pws1 !<power coefficients for each term and direction
    integer, dimension(3,ntpd2), intent(in) :: pws2 !<power coefficients for each term and directio
    real(gp), dimension(ntpd1), intent(in) :: ftpd1 !<factors of tensor product decompositions
    real(gp), dimension(ntpd2), intent(in) :: ftpd2 !<factors of tensor product decompositions
    real(gp), dimension(3), intent(in) :: dr !<separations between basis functions
    real(gp), dimension(NSD_,*), intent(in) :: sd1 !<exponents and coefficient
    real(gp), dimension(NSD_,*), intent(in) :: sd2 !<exponents and coefficient
    real(gp):: ovrlp !<scalar product
    !local variables
    integer :: ig1,ig2,i2,i1
    real(gp) ::  f,s1,s2,d1,d2,integral

    ovrlp=0.0_gp
    do ig2=1,ng2
       s2=sd2(EXPO_,ig2)
       d2=sd2(COEFF_,ig2)
       do ig1=1,ng1
          s1=sd1(EXPO_,ig1)
          d1=sd1(COEFF_,ig1)
          integral=0.0_gp
          do i2=1,ntpd2
             do i1=1,ntpd1
                f=  govrlp(s1,s2,dr(1),pws1(1,i1),pws2(1,i2))
                f=f*govrlp(s1,s2,dr(2),pws1(2,i1),pws2(2,i2))
                f=f*govrlp(s1,s2,dr(3),pws1(3,i1),pws2(3,i2))
                integral=integral+ftpd1(i1)*ftpd2(i2)*f
             end do
          end do
          ovrlp=ovrlp+d1*d2*integral
       end do
    end do
  end function gdot


!!$  !>calculate the density kernel matrix between two shells for a set of spatial points
!!$  pure subroutine density_kernel_shell
!!$    integer, intent(in) :: l1,l2 !<angular momenta of the shell
!!$    integer, intent(in) :: ntpdsh1,ntpdsh2 !<total number of shell tpd
!!$    real(gp), dimension(NSD_), intent(in) :: sd1,sd2 !<exponents and coefficient
!!$    integer, dimension(2*l1+1), intent(in) :: ntpd1 !<number of terms
!!$    integer, dimension(2*l2+1), intent(in) :: ntpd2 !<number of terms
!!$    integer, dimension(3,ntpdsh1), intent(in) :: pws1 !<exponents
!!$    integer, dimension(3,ntpdsh2), intent(in) :: pws2 !<exponents
!!$    real(gp), dimension(ntpdsh1), intent(in) :: ftpd1 !<factors of tpd, first shell
!!$    real(gp), dimension(ntpdsh2), intent(in) :: ftpd2 !<factors of tpd, second shell
!!$
!!$  end subroutine density_kernel_shell

  !> Performs the gaussian product for all the terms in the shell
  pure subroutine gdot_shell(sd1,l1,ntpdsh1,ntpd1,pws1,ftpd1,&
       sd2,l2,ntpdsh2,ntpd2,pws2,ftpd2,dr,overlap)
    implicit none
    integer, intent(in) :: l1,l2 !<angular momenta of the shell
    integer, intent(in) :: ntpdsh1,ntpdsh2 !<total number of shell tpd
    real(gp), dimension(NSD_), intent(in) :: sd1,sd2 !<exponents and coefficient
    integer, dimension(2*l1+1), intent(in) :: ntpd1 !<number of terms
    integer, dimension(2*l2+1), intent(in) :: ntpd2 !<number of terms
    integer, dimension(3,ntpdsh1), intent(in) :: pws1 !<exponents
    integer, dimension(3,ntpdsh2), intent(in) :: pws2 !<exponents
    real(gp), dimension(ntpdsh1), intent(in) :: ftpd1 !<factors of tpd, first shell
    real(gp), dimension(ntpdsh2), intent(in) :: ftpd2 !<factors of tpd, second shell
    real(gp), dimension(3), intent(in) :: dr !<separations between basis functions
    real(gp), dimension(2*l1+1,2*l2+1), intent(inout) :: overlap !<overlap of the shell
    !local variables
    integer :: m1,m2,i1,i2,itpd1,itpd2
    real(gp) :: f,integral

    itpd2=0
    do m2=1,2*l2+1
       itpd1=0
       do m1=1,2*l1+1
          integral=0.0_gp
          do i2=1,ntpd2(m2)
             do i1=1,ntpd1(m1)
                f=  govrlp(sd1(EXPO_),sd2(EXPO_),&
                     dr(1),pws1(1,i1+itpd1),pws2(1,i2+itpd2))
                f=f*govrlp(sd1(EXPO_),sd2(EXPO_),&
                     dr(2),pws1(2,i1+itpd1),pws2(2,i2+itpd2))
                f=f*govrlp(sd1(EXPO_),sd2(EXPO_),&
                     dr(3),pws1(3,i1+itpd1),pws2(3,i2+itpd2))
                integral=integral+ftpd1(i1+itpd1)*ftpd2(i2+itpd2)*f
             end do
          end do
          overlap(m1,m2)=overlap(m1,m2)+sd1(COEFF_)*sd2(COEFF_)*integral
          itpd1=itpd1+ntpd1(m1)
       end do
       itpd2=itpd2+ntpd2(m2)
    end do

  end subroutine gdot_shell


  !> Overlap matrix between two different basis structures
  !! laplacian is applied to the first one
  subroutine kinetic(A,B,ovrlp)
    implicit none
    type(gaussian_basis_new), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff), intent(out) :: ovrlp
    !local variables
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo,igA,igB
    integer :: ngA,ngB,lA,lB,mA,mB,nA,nB,ntpdshA,ntpdshB
    integer, dimension(2*L_MAX+1) :: ntpdA,ntpdB
    integer, dimension(3,NTERM_MAX_KINETIC) :: powA
    real(gp), dimension(NTERM_MAX_KINETIC) :: ftpdA
    integer, dimension(3,NTERM_MAX_OVERLAP) :: powB
    real(gp), dimension(NTERM_MAX_OVERLAP) :: ftpdB
    real(gp), dimension((2*L_MAX+1)*(2*L_MAX+1)) :: shell_overlap
    real(gp), dimension(3) :: dr,rA

    iovrlp=0
    ishell=0
    iexpo=0
    icoeff=1
    !loop on each shell (intensive calculation)
    do iat=1,A%nat
       rA(1)=A%rxyz(1,iat)
       rA(2)=A%rxyz(2,iat)
       rA(3)=A%rxyz(3,iat)
       do isat=1,A%nshell(iat)
          ishell=ishell+1
          ngA=A%shid(DOC_,ishell)
          lA=A%shid(L_,ishell)
          nA=A%shid(N_,ishell)

          jovrlp=0
          jshell=0
          jexpo=0
          jcoeff=1
          do jat=1,B%nat
             !here boundary conditions should be considered
             dr(1)=B%rxyz(1,jat)-rA(1)
             dr(2)=B%rxyz(2,jat)-rA(2)
             dr(3)=B%rxyz(3,jat)-rA(3)
             do jsat=1,B%nshell(jat)
                jshell=jshell+1
                ngB=B%shid(DOC_,jshell)
                lB=B%shid(L_,jshell)
                nB=B%shid(N_,jshell)
                !calculation of shell decomposition is independent of exponent
                call tensor_product_decomposition(nB,lB,ntpdshB,ntpdB,powB,ftpdB)

                shell_overlap=0.0_gp !overlap of the different terms
                do igA=1,ngA
                   call tensor_product_decomposition_laplacian(A%sd(EXPO_,igA+iexpo),nA,lA,&
                        ntpdshA,ntpdA,powA,ftpdA)
!                   call tensor_product_decomposition(nA,lA,ntpdshA,ntpdA,powA,ftpdA)
                   do igB=1,ngB
                      call gdot_shell(A%sd(1:,igA+iexpo),lA,ntpdshA,ntpdA,powA,ftpdA,&
                           B%sd(1:,igB+jexpo),lB,ntpdshB,ntpdB,powB,ftpdB,dr,&
                           shell_overlap)
                   end do
                end do

                !here the entire array should be copied in the right place
                do mB=1,2*lB+1
                   do mA=1,2*lA+1
                      ovrlp(iovrlp+mA,jovrlp+mB)=-0.5_gp*&
                           shell_overlap(mA+(mB-1)*(2*lA+1))
                   end do
                end do
                jexpo=jexpo+ngB
                jcoeff=jcoeff+2*lB+1
                jovrlp=jovrlp+2*lB+1
             end do
          end do
          iexpo=iexpo+ngA
          icoeff=icoeff+2*lA+1
          iovrlp=iovrlp+2*lA+1
       end do
    end do

    call gaudim_check(iexpo+1,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
    call gaudim_check(jexpo+1,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)

  END SUBROUTINE kinetic
!!$
!!$  !> Overlap matrix between two different basis structures
!!$  !! laplacian is applied to the first one
!!$  subroutine kinetic(A,B,ovrlp)
!!$    implicit none
!!$    type(gaussian_basis_new), intent(in) :: A,B
!!$    real(gp), dimension(A%ncoeff,B%ncoeff), intent(out) :: ovrlp
!!$    !local variables
!!$    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
!!$    integer :: iovrlp,jovrlp,jcoeff,jexpo,igA,igB,itpd1,itpd2,i1,i2
!!$    integer :: ngA,ngB,lA,lB,mA,mB,nA,nB,ntpdshA,ntpdshB
!!$    real(gp) :: f,integral
!!$    integer, dimension(2*L_MAX+1) :: ntpdA,ntpdB
!!$    integer, dimension(3,NTERM_MAX_KINETIC) :: powA
!!$    real(gp), dimension(NTERM_MAX_KINETIC) :: ftpdA
!!$    integer, dimension(3,NTERM_MAX_OVERLAP) :: powB
!!$    real(gp), dimension(NTERM_MAX_OVERLAP) :: ftpdB
!!$    real(gp), dimension((2*L_MAX+1)*(2*L_MAX+1)) :: shell_overlap
!!$    real(gp), dimension(3) :: dr,rA
!!$
!!$    iovrlp=0
!!$    ishell=0
!!$    iexpo=0
!!$    icoeff=1
!!$    !loop on each shell (intensive calculation)
!!$    do iat=1,A%nat
!!$       rA(1)=A%rxyz(1,iat)
!!$       rA(2)=A%rxyz(2,iat)
!!$       rA(3)=A%rxyz(3,iat)
!!$       do isat=1,A%nshell(iat)
!!$          ishell=ishell+1
!!$          ngA=A%shid(DOC_,ishell)
!!$          lA=A%shid(L_,ishell)
!!$          nA=A%shid(N_,ishell)
!!$
!!$          jovrlp=0
!!$          jshell=0
!!$          jexpo=0
!!$          jcoeff=1
!!$          do jat=1,B%nat
!!$             !here boundary conditions should be considered
!!$             dr(1)=B%rxyz(1,jat)-rA(1)
!!$             dr(2)=B%rxyz(2,jat)-rA(2)
!!$             dr(3)=B%rxyz(3,jat)-rA(3)
!!$             do jsat=1,B%nshell(jat)
!!$                jshell=jshell+1
!!$                ngB=B%shid(DOC_,jshell)
!!$                lB=B%shid(L_,jshell)
!!$                nB=B%shid(N_,jshell)
!!$                !calculation of shell decomposition is independent of exponent
!!$                call tensor_product_decomposition(nB,lB,ntpdshB,ntpdB,powB,ftpdB)
!!$
!!$                do mB=1,2*lB+1
!!$                   do mA=1,2*lA+1
!!$                      ovrlp(iovrlp+mA,jovrlp+mB)=0.0_gp
!!$                   end do
!!$                end do
!!$                do igA=1,ngA
!!$                   call tensor_product_decomposition_laplacian(A%sd(EXPO_,igA+iexpo),nA,lA,&
!!$                        ntpdshA,ntpdA,powA,ftpdA)
!!$!                   call tensor_product_decomposition(nA,lA,ntpdshA,ntpdA,powA,ftpdA)
!!$                   do igB=1,ngB
!!$                      itpd2=0
!!$                      do mB=1,2*lB+1
!!$                         itpd1=0
!!$                         do mA=1,2*lA+1
!!$                            integral=0.0_gp
!!$                            do i2=1,ntpdB(mB)
!!$                               do i1=1,ntpdA(mA)
!!$                                  f=  govrlp(A%sd(EXPO_,igA+iexpo),B%sd(EXPO_,igB+jexpo),&
!!$                                       dr(1),powA(1,i1+itpd1),powB(1,i2+itpd2))
!!$                                  f=f*govrlp(A%sd(EXPO_,igA+iexpo),B%sd(EXPO_,igB+jexpo),&
!!$                                       dr(2),powA(2,i1+itpd1),powB(2,i2+itpd2))
!!$                                  f=f*govrlp(A%sd(EXPO_,igA+iexpo),B%sd(EXPO_,igB+jexpo),&
!!$                                       dr(3),powA(3,i1+itpd1),powB(3,i2+itpd2))
!!$                                  integral=integral+ftpdA(i1+itpd1)*ftpdB(i2+itpd2)*f
!!$                               end do
!!$                            end do
!!$                            ovrlp(iovrlp+mA,jovrlp+mB)=ovrlp(iovrlp+mA,jovrlp+mB)-0.5_gp*&
!!$                                 A%sd(COEFF_,igA+iexpo)*B%sd(COEFF_,igB+jexpo)*integral
!!$                            itpd1=itpd1+ntpdA(mA)
!!$                         end do
!!$                         itpd2=itpd2+ntpdB(mB)
!!$                      end do
!!$
!!$!!!$                      call gdot_shell(A%sd(1:,igA+iexpo),lA,ntpdshA,ntpdA,powA,ftpdA,&
!!$!!!$                           B%sd(1:,igB+jexpo),lB,ntpdshB,ntpdB,powB,ftpdB,dr,&
!!$!!!$                           shell_overlap)
!!$                   end do
!!$                end do
!!$
!!$!!!$                !here the entire array should be copied in the right place
!!$!!!$                do mB=1,2*lB+1
!!$!!!$                   do mA=1,2*lA+1
!!$!!!$                      ovrlp(iovrlp+mA,jovrlp+mB)=-0.5_gp*&
!!$!!!$                           shell_overlap(mA+(mB-1)*(2*lA+1))
!!$!!!$                   end do
!!$!!!$                end do
!!$                jexpo=jexpo+ngB
!!$                jcoeff=jcoeff+2*lB+1
!!$                jovrlp=jovrlp+2*lB+1
!!$             end do
!!$          end do
!!$          iexpo=iexpo+ngA
!!$          icoeff=icoeff+2*lA+1
!!$          iovrlp=iovrlp+2*lA+1
!!$       end do
!!$    end do
!!$
!!$    call gaudim_check(iexpo+1,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
!!$    call gaudim_check(jexpo+1,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)
!!$
!!$  END SUBROUTINE kinetic
  


  !> Calculates @f$\int e^{-a1*x^2} x^l1 \exp^{-a2*(x-d)^2} (x-d)^l2 dx@f$
  !! Uses the function gauint0 if d==0.0
  pure function govrlp(a1,a2,d,l1,l2)
    implicit none
    integer, intent(in) :: l1,l2
    real(gp), intent(in) :: a1,a2,d
    real(gp) :: govrlp
    !local variables
    integer :: p
    real(gp) :: prefac,stot,aeff,ceff,tt,fsum!,gauint,gauint0

    !quick check
    if (d == 0.0_gp) then
       govrlp=gauint0(a1+a2,l1+l2)
       return
    end if

    !build the prefactor
    prefac=a1+a2
    prefac=a2/prefac
    prefac=a1*prefac
    prefac=-d**2*prefac
    prefac=exp(prefac)

    !build the effective exponent and coefficients
    aeff=a1+a2
    ceff=a2*d
    ceff=ceff/aeff

    !build the first term in the sum
    stot=(-d)**l2
    stot=gauint(aeff,ceff,l1)*stot

    !perform the sum
    do p=1,l2/2
       tt=rfac(1,p)
       fsum=rfac(l2-p+1,l2)
       fsum=fsum/tt
       tt=(-d)**(l2-p)
       fsum=fsum*tt
       tt=gauint(aeff,ceff,l1+p)
       fsum=fsum*tt
       stot=stot+fsum
    end do
    do p=l2/2+1,l2
       tt=rfac(1,l2-p)
       fsum=rfac(p+1,l2)
       fsum=fsum/tt
       tt=(-d)**(l2-p)
       fsum=fsum*tt
       tt=gauint(aeff,ceff,l1+p)
       fsum=fsum*tt
       stot=stot+fsum
    end do

    !final result
    govrlp=prefac*stot
  END FUNCTION govrlp
  

  !> Kinetic overlap between gaussians, based on cartesian coordinates
  !! calculates a dot product between two differents gaussians times spherical harmonics
  !! valid only for shell which belongs to different atoms, and with also dy/=0/=dx dz/=0
  !! to be rearranged when only some of them is zero
  subroutine kinprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: l1,l2,m1,m2,niw,nrw 
    real(gp), intent(in) :: a1,a2,dx,dy,dz
    integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
    real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
    real(gp), intent(out) :: ovrlp
    !local variables
    integer, parameter :: nx=3
    integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
    real(gp) :: fx,fy,fz,fa,fb,d2fx,d2fy,d2fz!,govrlp,kinovrlp

    !calculates the number of different couples
    call calc_coeff_inguess(l1,m1,nx,n1,&
         iw(1),iw(nx+1),iw(2*nx+1),rw(1))
    call calc_coeff_inguess(l2,m2,nx,n2,&
         iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))
    ovrlp=0.0_gp
    do i2=1,n2
       qx=iw(3*nx+i2)
       qy=iw(4*nx+i2)
       qz=iw(5*nx+i2)
       fb=rw(n1+i2)
       do i1=1,n1
          px=iw(i1)
          py=iw(nx+i1)
          pz=iw(2*nx+i1)
          fa=rw(i1)

          fx=govrlp(a1,a2,dx,px,qx)
          fy=govrlp(a1,a2,dy,py,qy)
          fz=govrlp(a1,a2,dz,pz,qz)

          d2fx=kinovrlp(a1,a2,dx,px,qx)
          d2fy=kinovrlp(a1,a2,dy,py,qy)
          d2fz=kinovrlp(a1,a2,dz,pz,qz)

          ovrlp=ovrlp-0.5_gp*fa*fb*(d2fx*fy*fz+fx*d2fy*fz+fx*fy*d2fz)
          !print *,i1,i2,fx,fy,fz,fa,fb
       end do
    end do

  END SUBROUTINE kinprod


  !> Calculates @f$\int d^2/dx^2(\exp^{-a1*x^2} x^l1) \exp^{-a2*(x-d)^2} (x-d)^l2 dx@f$
  !! in terms of the govrlp function below
  pure function kinovrlp(a1,a2,d,l1,l2)
    implicit none
    integer, intent(in) :: l1,l2
    real(gp), intent(in) :: a1,a2,d
    real(gp) :: kinovrlp
    !local variables
    real(gp) :: fac,ovrlp!govrlp

    !case l1+2
    fac=4._gp*a1**2
    ovrlp=govrlp(a1,a2,d,l1+2,l2)
    kinovrlp=fac*ovrlp
    !case l1
    fac=2._gp*a1*real(2*l1+1,gp)
    ovrlp=govrlp(a1,a2,d,l1,l2)
    kinovrlp=kinovrlp-fac*ovrlp
    !case l1-2 (if applicable)
    if (l1 >=2) then
       fac=real(l1*(l1-1),gp)
       ovrlp=govrlp(a1,a2,d,l1-2,l2)
       kinovrlp=kinovrlp+fac*ovrlp
    end if
  END FUNCTION kinovrlp


  !> Calculates @f$ I(l) = \int \exp^{-a*x^2} x^l dx @f$ i.e. the moments of the gaussian
  !! @f$ I(0)  = \sqrt{\pi} @f$
  !! @f$ I(2p) = (p-1/2) I(2p-2) @f$
  !! this works for all l
  pure function gauint0(a,l)
    implicit none
    !Arguments
    integer, intent(in) :: l
    real(gp), intent(in) :: a
    real(gp) :: gauint0
    !local variables
    real(gp), parameter :: gammaonehalf=1.772453850905516027298_gp ! i.e. sqrt(pi)
    integer :: p
    real(gp) :: prefac,tt
    !build the prefactor
    prefac=1.0_gp/sqrt(a)
    prefac=gammaonehalf*prefac**(l+1)

    p=l/2
    if (2*p < l) then
       ! l is odd
       gauint0=0.0_gp
    else
       tt=xfac(1,p,-0.5_gp)
       !final result
       gauint0=prefac*tt
    end if

  END FUNCTION gauint0


  !> Calculates @f$\int \exp^{-a*(x-c)^2} x^l dx@f$
  !! this works ONLY when c /= 0.d0
  pure function gauint(a,c,l)
    implicit none
    integer, intent(in) :: l
    real(gp), intent(in) :: a,c
    real(gp) :: gauint
    !local variables
    real(gp), parameter :: gammaonehalf=1.772453850905516027298_gp
    integer :: p
    real(gp) :: prefac,stot,fsum,tt!,firstprod

    !quick check
    !if (c==0.0_gp) then
    !   stop 'gauint0 should be called'
    !end if

    !build the prefactor
    prefac=sqrt(a)
    prefac=1.0_gp/prefac
    prefac=gammaonehalf*prefac

    !the first term of the sum is one
    !but we have to multiply for the prefactor
    stot=c**l

    !if (c==0.0_gp .and. l==0) then
    !   do p=0,20
    !      print *,'stot,p',stot,a,p,gauint0(a,p)
    !   end do
    !end if

    !calculate the sum
    do p=1,l/4
       tt=rfac(p+1,2*p)
       fsum=rfac(l-2*p+1,l)
       fsum=fsum/tt
       tt=firstprod(p)
       fsum=fsum*tt
       tt=c**(l-2*p)
       tt=tt/a**p
       fsum=fsum*tt
       stot=stot+fsum
    end do
    do p=l/4+1,l/3
       tt=rfac(p+1,l-2*p)
       fsum=rfac(2*p+1,l)
       fsum=fsum/tt
       tt=firstprod(p)
       fsum=fsum*tt
       tt=c**(l-2*p)
       tt=tt/a**p
       fsum=fsum*tt
       stot=stot+fsum
    end do
    do p=l/3+1,l/2
       tt=rfac(l-2*p+1,p)
       fsum=rfac(2*p+1,l)
       fsum=fsum*tt
       tt=firstprod(p)
       fsum=fsum*tt
       tt=c**(l-2*p)
       tt=tt/a**p
       fsum=fsum*tt
       stot=stot+fsum
    end do

    !final result
    gauint=stot*prefac

  END FUNCTION gauint


  pure function firstprod(p)
    implicit none
    integer, intent(in) :: p
    real(gp) :: firstprod
    !local variables
    integer :: i
    real(gp) :: tt
    firstprod=1.0_gp
    do i=1,p
       tt=real(2*i,gp)
       tt=1.0_gp/tt
       tt=1.0_gp-tt
       firstprod=firstprod*tt
    end do
  END FUNCTION firstprod


  subroutine gaudim_check(iexpo,icoeff,ishell,nexpo,ncoeff,nshltot)
    implicit none
    integer, intent(in) :: iexpo,icoeff,ishell,nexpo,ncoeff,nshltot
    !check of the dimensions
    if (iexpo /= nexpo+1) then
       write(*,*)' ERROR: nexpo+1 <> iexpo',nexpo,iexpo
       stop
    else if (icoeff /= ncoeff+1) then
       write(*,*)' ERROR: ncoeff+1 <> icoeff',ncoeff,icoeff
       stop
    else if (ishell /= nshltot) then
       write(*,*)' ERROR: nshltot <> ishell',nshltot,ishell
       stop
    end if
  END SUBROUTINE gaudim_check


  !> Normalize a given atomic shell following the angular momentum
  pure subroutine normalize_shell(ng,l,expo,coeff)
    implicit none
    integer, intent(in) :: ng,l
    real(gp), dimension(ng), intent(in) :: expo
    real(gp), dimension(ng), intent(inout) :: coeff
    !local variables
    integer :: i,j
    real(gp) :: norm,tt,e1,ex,c1,c2!,gauint0

    norm=0.0_gp
    do i=1,ng
       e1=expo(i)
       c1=coeff(i)
       do j=1,ng
          ex=expo(j)+e1
          c2=coeff(j)
          tt=gauint0(ex,2*l+2)
          norm=norm+c1*tt*c2
       end do
    end do
    norm=sqrt(0.5_gp*norm)
    norm=1.0_gp/norm
    do i=1,ng
       coeff(i)=coeff(i)*norm
    end do

    !print *,'l=',l,'norm=',norm

  END SUBROUTINE normalize_shell


  !> Factorial (float)
  pure function rfac(is,ie)
    implicit none
    integer, intent(in) :: is,ie
    real(gp) :: rfac
    !local variables
    integer :: i
    real(gp) :: tt
    rfac=1.0_gp
    do i=is,ie
       tt=real(i,gp)
       rfac=rfac*tt
    end do
  END FUNCTION rfac


  !> Partial factorial, with real shift
  !! With this function n! = xfac(1,n,0.d0)
  !! @f$ \prod_1^n (n - s) @f$
  pure function xfac(is,ie,sh)
    implicit none
    !Arguments
    integer, intent(in) :: is,ie
    real(gp), intent(in) :: sh
    real(gp) :: xfac
    !local variables
    integer :: i
    real(gp) :: tt
    xfac=1.0_gp
    do i=is,ie
       tt=real(i,gp)+sh
       xfac=xfac*tt
    end do
  END FUNCTION xfac


  !> Routine to extract the coefficients from the quantum numbers and the operation
  pure subroutine tensor_product_decomposition(n,l,ntpd_shell,ntpd,pow,ftpd)
    implicit none
    !Arguments
    integer, intent(in) :: n,l
    integer, intent(out) :: ntpd_shell !< No. of terms for the whole shell
    integer, dimension(2*l+1), intent(out) :: ntpd !< number of terms per shell element
    integer, dimension(3,NTERM_MAX_OVERLAP), intent(out) :: pow !< tensor product decompositions
    real(gp), dimension(NTERM_MAX_OVERLAP), intent(out) :: ftpd !<factors
    !All data alltogether
    real(gp), parameter :: rsp=0.564189583547756286948079451560772585844050629328998856844086_gp !1/sqrt(pi)
    select case(n)
    case(1)
       select case(l)
       case(0)
          ntpd_shell=1 ! l=0,n=1
          ntpd(1:1)=(/1/) ! l=0,n=1
          pow(1:3,1:1)=reshape((/0,0,0/),(/3,1/)) ! l=0,n=1
          ftpd(1:1)=rsp*(/0.5_gp/) ! m=1, l=0, n=1
       case(1)
          ntpd_shell=3 ! l=1,n=1
          ntpd(1:3)=(/1,1,1/) ! l=1,n=1
          pow(1:3,1:3)=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/)) ! l=1,n=1
          ftpd(1:1)=rsp*(/sqrt(0.75_gp)/) ! m=1, l=1, n=1
          ftpd(2:2)=rsp*(/sqrt(0.75_gp)/) ! m=2, l=1, n=1
          ftpd(3:3)=rsp*(/sqrt(0.75_gp)/) ! m=3, l=1, n=1
       case(2)
          ntpd_shell=8 ! l=2,n=1
          ntpd(1:5)=(/1,1,1,2,3/) ! l=2,n=1
          pow(1:3,1:8)=reshape((/0,1,1,1,0,1,1,1,0,2,0,0,0,2,0,2,0,0,0,2,0,0,0,2/),(/3,8/)) ! l=2,n=1
          ftpd(1:1)=rsp*(/sqrt(3.75_gp)/) ! m=1, l=2, n=1
          ftpd(2:2)=rsp*(/sqrt(3.75_gp)/) ! m=2, l=2, n=1
          ftpd(3:3)=rsp*(/sqrt(3.75_gp)/) ! m=3, l=2, n=1
          ftpd(4:5)=rsp*(/sqrt(0.9375_gp),-sqrt(0.9375_gp)/) ! m=4, l=2, n=1
          ftpd(6:8)=rsp*(/-sqrt(0.3125_gp),-sqrt(0.3125_gp),sqrt(1.25_gp)/) ! m=5, l=2, n=1
       case(3)
          ntpd_shell=16 ! l=3,n=1
          ntpd(1:7)=(/3,3,3,2,2,2,1/) ! l=3,n=1
          pow(1:3,1:16)=reshape((/3,0,0,1,2,0,1,0,2,2,1,0,0,3,0,0,1,2,2,0,1,0,2,1,0,0,3,3,0,0,1,&
               2,0,2,1,0,0,3,0,2,0,1,0,2,1,1,1,1/),(/3,16/)) ! l=3,n=1
          ftpd(1:3)=rsp*(/sqrt(0.65625_gp),sqrt(0.65625_gp),-sqrt(10.5_gp)/) ! m=1, l=3, n=1
          ftpd(4:6)=rsp*(/sqrt(0.65625_gp),sqrt(0.65625_gp),-sqrt(10.5_gp)/) ! m=2, l=3, n=1
          ftpd(7:9)=rsp*(/sqrt(3.9375_gp),sqrt(3.9375_gp),-sqrt(1.75_gp)/) ! m=3, l=3, n=1
          ftpd(10:11)=rsp*(/sqrt(1.09375_gp),-sqrt(9.84375_gp)/) ! m=4, l=3, n=1
          ftpd(12:13)=rsp*(/-sqrt(9.84375_gp),sqrt(1.09375_gp)/) ! m=5, l=3, n=1
          ftpd(14:15)=rsp*(/sqrt(6.5625_gp),-sqrt(6.5625_gp)/) ! m=6, l=3, n=1
          ftpd(16:16)=rsp*(/sqrt(26.25_gp)/) ! m=7, l=3, n=1
       end select
    case(2)
       select case(l)
       case(0)
          ntpd_shell=3 ! l=0,n=2
          ntpd(1:1)=(/3/) ! l=0,n=2
          pow(1:3,1:3)=reshape((/2,0,0,0,2,0,0,0,2/),(/3,3/)) ! l=0,n=2
          ftpd(1:3)=rsp*(/0.5_gp,0.5_gp,0.5_gp/) ! m=1, l=0, n=2
       case(1)
          ntpd_shell=9 ! l=1,n=2
          ntpd(1:3)=(/3,3,3/) ! l=1,n=2
          pow(1:3,1:9)=reshape((/3,0,0,1,2,0,1,0,2,2,1,0,0,3,0,0,1,2,2,0,1,0,2,1,0,0,3/),(/3,9/)) ! l=1,n=2
          ftpd(1:3)=rsp*(/sqrt(0.75_gp),sqrt(0.75_gp),sqrt(0.75_gp)/) ! m=1, l=1, n=2
          ftpd(4:6)=rsp*(/sqrt(0.75_gp),sqrt(0.75_gp),sqrt(0.75_gp)/) ! m=2, l=1, n=2
          ftpd(7:9)=rsp*(/sqrt(0.75_gp),sqrt(0.75_gp),sqrt(0.75_gp)/) ! m=3, l=1, n=2
       case(2)
          ntpd_shell=19 ! l=2,n=2
          ntpd(1:5)=(/3,3,3,4,6/) ! l=2,n=2
          pow(1:3,1:19)=reshape((/2,1,1,0,3,1,0,1,3,3,0,1,1,2,1,1,0,3,3,1,0,1,3,0,1,1,2,4,0,0,0,4,&
               0,2,0,2,0,2,2,4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4/),(/3,19/)) ! l=2,n=2
          ftpd(1:3)=rsp*(/sqrt(3.75_gp),sqrt(3.75_gp),sqrt(3.75_gp)/) ! m=1, l=2, n=2
          ftpd(4:6)=rsp*(/sqrt(3.75_gp),sqrt(3.75_gp),sqrt(3.75_gp)/) ! m=2, l=2, n=2
          ftpd(7:9)=rsp*(/sqrt(3.75_gp),sqrt(3.75_gp),sqrt(3.75_gp)/) ! m=3, l=2, n=2
          ftpd(10:13)=rsp*(/sqrt(0.9375_gp),-sqrt(0.9375_gp),sqrt(0.9375_gp),-sqrt(0.9375_gp)/) ! m=4, l=2, n=2
          ftpd(14:19)=rsp*(/-sqrt(0.3125_gp),-sqrt(1.25_gp),-sqrt(0.3125_gp),sqrt(0.3125_gp),sqrt(0.3125_gp),&
               sqrt(1.25_gp)/) ! m=5, l=2, n=2
       case(3)
          ntpd_shell=35 ! l=3,n=2
          ntpd(1:7)=(/6,6,6,5,5,4,3/) ! l=3,n=2
          pow(1:3,1:35)=reshape((/5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,0,1,4,4,&
               0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5,5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,&
               4,0,1,0,4,1,2,0,3,0,2,3,3,1,1,1,3,1,1,1,3/),(/3,35/)) ! l=3,n=2
          ftpd(1:6)=rsp*(/sqrt(0.65625_gp),sqrt(2.625_gp),sqrt(0.65625_gp),-sqrt(5.90625_gp),-sqrt(5.90625_gp),&
               -sqrt(10.5_gp)/) ! m=1, l=3, n=2
          ftpd(7:12)=rsp*(/sqrt(0.65625_gp),sqrt(2.625_gp),sqrt(0.65625_gp),-sqrt(5.90625_gp),-sqrt(5.90625_gp),&
               -sqrt(10.5_gp)/) ! m=2, l=3, n=2
          ftpd(13:18)=rsp*(/sqrt(3.9375_gp),sqrt(15.75_gp),sqrt(3.9375_gp),sqrt(0.4375_gp),sqrt(0.4375_gp),&
               -sqrt(1.75_gp)/) ! m=3, l=3, n=2
          ftpd(19:23)=rsp*(/sqrt(1.09375_gp),-sqrt(4.375_gp),-sqrt(9.84375_gp),sqrt(1.09375_gp),-sqrt(9.84375_gp)/) ! m=4, l=3, n=2
          ftpd(24:28)=rsp*(/-sqrt(9.84375_gp),-sqrt(4.375_gp),sqrt(1.09375_gp),-sqrt(9.84375_gp),sqrt(1.09375_gp)/) ! m=5, l=3, n=2
          ftpd(29:32)=rsp*(/sqrt(6.5625_gp),-sqrt(6.5625_gp),sqrt(6.5625_gp),-sqrt(6.5625_gp)/) ! m=6, l=3, n=2
          ftpd(33:35)=rsp*(/sqrt(26.25_gp),sqrt(26.25_gp),sqrt(26.25_gp)/) ! m=7, l=3, n=2
       end select
    case(3)
       select case(l)
       case(0)
          ntpd_shell=6 ! l=0,n=3
          ntpd(1:1)=(/6/) ! l=0,n=3
          pow(1:3,1:6)=reshape((/4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4/),(/3,6/)) ! l=0,n=3
          ftpd(1:6)=rsp*(/0.5_gp,1._gp,0.5_gp,1._gp,1._gp,0.5_gp/) ! m=1, l=0, n=3
       case(1)
          ntpd_shell=18 ! l=1,n=3
          ntpd(1:3)=(/6,6,6/) ! l=1,n=3
          pow(1:3,1:18)=reshape((/5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,0,1,4,4,&
               0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5/),(/3,18/)) ! l=1,n=3
          ftpd(1:6)=rsp*(/sqrt(0.75_gp),sqrt(3._gp),sqrt(0.75_gp),sqrt(3._gp),sqrt(3._gp),sqrt(0.75_gp)/) ! m=1, l=1, n=3
          ftpd(7:12)=rsp*(/sqrt(0.75_gp),sqrt(3._gp),sqrt(0.75_gp),sqrt(3._gp),sqrt(3._gp),sqrt(0.75_gp)/) ! m=2, l=1, n=3
          ftpd(13:18)=rsp*(/sqrt(0.75_gp),sqrt(3._gp),sqrt(0.75_gp),sqrt(3._gp),sqrt(3._gp),sqrt(0.75_gp)/) ! m=3, l=1, n=3
       case(2)
          ntpd_shell=33 ! l=2,n=3
          ntpd(1:5)=(/6,6,6,8,7/) ! l=2,n=3
          pow(1:3,1:33)=reshape((/4,1,1,2,3,1,0,5,1,2,1,3,0,3,3,0,1,5,5,0,1,3,2,1,1,4,1,3,0,3,1,2,3,1,0,5,5,1,&
               0,3,3,0,1,5,0,3,1,2,1,3,2,1,1,4,6,0,0,4,2,0,2,4,0,0,6,0,4,0,2,0,4,2,2,0,4,0,2,4,6,0,0,4,2,0,2,4,&
               0,0,6,0,2,0,4,0,2,4,0,0,6/),(/3,33/)) ! l=2,n=3
          ftpd(1:6)=rsp*(/sqrt(3.75_gp),sqrt(15._gp),sqrt(3.75_gp),sqrt(15._gp),sqrt(15._gp),sqrt(3.75_gp)/) ! m=1, l=2, n=3
          ftpd(7:12)=rsp*(/sqrt(3.75_gp),sqrt(15._gp),sqrt(3.75_gp),sqrt(15._gp),sqrt(15._gp),sqrt(3.75_gp)/) ! m=2, l=2, n=3
          ftpd(13:18)=rsp*(/sqrt(3.75_gp),sqrt(15._gp),sqrt(3.75_gp),sqrt(15._gp),sqrt(15._gp),sqrt(3.75_gp)/) ! m=3, l=2, n=3
          ftpd(19:26)=rsp*(/sqrt(0.9375_gp),sqrt(0.9375_gp),-sqrt(0.9375_gp),-sqrt(0.9375_gp),sqrt(3.75_gp),&
               -sqrt(3.75_gp),sqrt(0.9375_gp),-sqrt(0.9375_gp)/) ! m=4, l=2, n=3
          ftpd(27:33)=rsp*(/-sqrt(0.3125_gp),-sqrt(2.8125_gp),-sqrt(2.8125_gp),-sqrt(0.3125_gp),sqrt(2.8125_gp),&
               sqrt(2.8125_gp),sqrt(1.25_gp)/) ! m=5, l=2, n=3
       case(3)
          ntpd_shell=62 ! l=3,n=3
          ntpd(1:7)=(/10,10,10,9,9,8,6/) ! l=3,n=3
          pow(1:3,1:62)=reshape((/7,0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,4,2,3,0,4,1,2,4,1,0,6,6,1,0,4,3,0,2,5,0,0,&
               7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,3,4,0,1,6,6,0,1,4,2,1,2,4,1,0,6,1,4,0,3,2,2,3,0,4,3,2,0,5,0,2,5,0,0,7,7,&
               0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,4,2,3,0,4,1,2,4,6,1,0,4,3,0,2,5,0,0,7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,&
               3,4,6,0,1,4,2,1,2,4,1,0,6,1,4,0,3,0,4,3,2,0,5,0,2,5,5,1,1,3,3,1,1,5,1,3,1,3,1,3,3,1,1,5/),(/3,62/)) ! l=3,n=3
          ftpd(1:10)=rsp*(/sqrt(0.65625_gp),sqrt(5.90625_gp),sqrt(5.90625_gp),sqrt(0.65625_gp),-sqrt(2.625_gp),&
               -sqrt(10.5_gp),-sqrt(2.625_gp),-sqrt(32.15625_gp),-sqrt(32.15625_gp),-sqrt(10.5_gp)/) ! m=1, l=3, n=3
          ftpd(11:20)=rsp*(/sqrt(0.65625_gp),sqrt(5.90625_gp),sqrt(5.90625_gp),sqrt(0.65625_gp),-sqrt(2.625_gp),&
               -sqrt(10.5_gp),-sqrt(2.625_gp),-sqrt(32.15625_gp),-sqrt(32.15625_gp),-sqrt(10.5_gp)/) ! m=2, l=3, n=3
          ftpd(21:30)=rsp*(/sqrt(3.9375_gp),sqrt(35.4375_gp),sqrt(35.4375_gp),sqrt(3.9375_gp),sqrt(7._gp),&
               sqrt(28._gp),sqrt(7._gp),-sqrt(0.4375_gp),-sqrt(0.4375_gp),-sqrt(1.75_gp)/) ! m=3, l=3, n=3
          ftpd(31:39)=rsp*(/sqrt(1.09375_gp),-sqrt(1.09375_gp),-sqrt(27.34375_gp),-sqrt(9.84375_gp),sqrt(4.375_gp),&
               -sqrt(17.5_gp),-sqrt(39.375_gp),sqrt(1.09375_gp),-sqrt(9.84375_gp)/) ! m=4, l=3, n=3
          ftpd(40:48)=rsp*(/-sqrt(9.84375_gp),-sqrt(27.34375_gp),-sqrt(1.09375_gp),sqrt(1.09375_gp),-sqrt(39.375_gp),&
               -sqrt(17.5_gp),sqrt(4.375_gp),-sqrt(9.84375_gp),sqrt(1.09375_gp)/) ! m=5, l=3, n=3
          ftpd(49:56)=rsp*(/sqrt(6.5625_gp),sqrt(6.5625_gp),-sqrt(6.5625_gp),-sqrt(6.5625_gp),sqrt(26.25_gp),&
               -sqrt(26.25_gp),sqrt(6.5625_gp),-sqrt(6.5625_gp)/) ! m=6, l=3, n=3
          ftpd(57:62)=rsp*(/sqrt(26.25_gp),sqrt(105._gp),sqrt(26.25_gp),sqrt(105._gp),sqrt(105._gp),&
               sqrt(26.25_gp)/) ! m=7, l=3, n=3
       end select
    end select
  end subroutine tensor_product_decomposition


  !> Routine to extract the coefficients from the quantum numbers and the operation
  !! it provides the tensor product decomposition of the laplacian of a given shell
  pure subroutine tensor_product_decomposition_laplacian(a,n,l,ntpd_shell,ntpd,pow,ftpd)
    implicit none
    integer, intent(in) :: n,l
    real(gp), intent(in) :: a !<exponent of the gaussian
    integer, intent(out) :: ntpd_shell !< No. of terms for the whole shell
    integer, dimension(2*l+1), intent(out) :: ntpd !< number of terms per shell element
    integer, dimension(3,NTERM_MAX_KINETIC), intent(out) :: pow !< tensor product decompositions
    real(gp), dimension(NTERM_MAX_KINETIC), intent(out) :: ftpd !<factors
    !All data alltogether
    real(gp), parameter :: rsp=0.564189583547756286948079451560772585844050629328998856844086_gp !1/sqrt(pi)
    select case(n)
    case(1)
       select case(l)
       case(0)
          ntpd_shell=4 ! l=0,n=1
          ntpd(1:1)=(/4/) ! l=0,n=1
          pow(1:3,1:4)=reshape((/2,0,0,0,2,0,0,0,2,0,0,0/),(/3,4/)) ! l=0,n=1
          ftpd(1:4)=rsp*(/2._gp*a**2,2._gp*a**2,2._gp*a**2,-3._gp*a/) ! m=1, l=0, n=1
       case(1)
          ntpd_shell=12 ! l=1,n=1
          ntpd(1:3)=(/4,4,4/) ! l=1,n=1
          pow(1:3,1:12)=reshape((/3,0,0,1,2,0,1,0,2,1,0,0,2,1,0,0,3,0,0,1,2,0,1,0,2,0,1,0,2,1,0,0,3,0,0,1/),(/3,12/)) ! l=1,n=1
          ftpd(1:4)=rsp*(/sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(75._gp)*a/) ! m=1, l=1, n=1
          ftpd(5:8)=rsp*(/sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(75._gp)*a/) ! m=2, l=1, n=1
          ftpd(9:12)=rsp*(/sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(75._gp)*a/) ! m=3, l=1, n=1
       case(2)
          ntpd_shell=27 ! l=2,n=1
          ntpd(1:5)=(/4,4,4,6,9/) ! l=2,n=1
          pow(1:3,1:27)=reshape((/2,1,1,0,3,1,0,1,3,0,1,1,3,0,1,1,2,1,1,0,3,1,0,1,3,1,0,1,3,0,1,1,2,1,1,0,&
               4,0,0,0,4,0,2,0,2,0,2,2,2,0,0,0,2,0,4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4,2,0,0,0,2,0,0,0,2/),(/3,27/)) ! l=2,n=1
          ftpd(1:4)=rsp*(/sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(735._gp)*a/) ! m=1, l=2, n=1
          ftpd(5:8)=rsp*(/sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(735._gp)*a/) ! m=2, l=2, n=1
          ftpd(9:12)=rsp*(/sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(735._gp)*a/) ! m=3, l=2, n=1
          ftpd(13:18)=rsp*(/sqrt(15._gp)*a**2,-sqrt(15._gp)*a**2,sqrt(15._gp)*a**2,-sqrt(15._gp)*a**2,&
               -sqrt(183.75_gp)*a,sqrt(183.75_gp)*a/) ! m=4, l=2, n=1
          ftpd(19:27)=rsp*(/-sqrt(5._gp)*a**2,-sqrt(20._gp)*a**2,-sqrt(5._gp)*a**2,sqrt(5._gp)*a**2,&
               sqrt(5._gp)*a**2,sqrt(20._gp)*a**2,sqrt(61.25_gp)*a,sqrt(61.25_gp)*a,-sqrt(245._gp)*a/) ! m=5, l=2, n=1
       case(3)
          ntpd_shell=51 ! l=3,n=1
          ntpd(1:7)=(/9,9,9,7,7,6,4/) ! l=3,n=1
          pow(1:3,1:51)=reshape((/5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,3,0,0,1,2,0,1,0,2,4,1,0,2,3,0,0,5,0,&
               2,1,2,0,3,2,0,1,4,2,1,0,0,3,0,0,1,2,4,0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5,2,0,1,0,2,1,0,0,3,5,&
               0,0,3,2,0,1,4,0,3,0,2,1,2,2,3,0,0,1,2,0,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,2,1,0,0,3,0,4,0,1,0,4,1,2,&
               0,3,0,2,3,2,0,1,0,2,1,3,1,1,1,3,1,1,1,3,1,1,1/),(/3,51/)) ! l=3,n=1
          ftpd(1:9)=rsp*(/sqrt(10.5_gp)*a**2,sqrt(42._gp)*a**2,sqrt(10.5_gp)*a**2,-sqrt(94.5_gp)*a**2,&
               -sqrt(94.5_gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(212.625_gp)*a,-sqrt(212.625_gp)*a,&
               sqrt(3402._gp)*a/) ! m=1, l=3, n=1
          ftpd(10:18)=rsp*(/sqrt(10.5_gp)*a**2,sqrt(42._gp)*a**2,sqrt(10.5_gp)*a**2,&
               -sqrt(94.5_gp)*a**2,-sqrt(94.5_gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(212.625_gp)*a,&
               -sqrt(212.625_gp)*a,sqrt(3402._gp)*a/) ! m=2, l=3, n=1
          ftpd(19:27)=rsp*(/sqrt(63._gp)*a**2,sqrt(252._gp)*a**2,sqrt(63._gp)*a**2,sqrt(7._gp)*a**2,&
               sqrt(7._gp)*a**2,-sqrt(28._gp)*a**2,-sqrt(1275.75_gp)*a,-sqrt(1275.75_gp)*a,sqrt(567._gp)*a/) ! m=3, l=3, n=1
          ftpd(28:34)=rsp*(/sqrt(17.5_gp)*a**2,-sqrt(70._gp)*a**2,-sqrt(157.5_gp)*a**2,sqrt(17.5_gp)*a**2,&
               -sqrt(157.5_gp)*a**2,-sqrt(354.375_gp)*a,sqrt(3189.375_gp)*a/) ! m=4, l=3, n=1
          ftpd(35:41)=rsp*(/-sqrt(157.5_gp)*a**2,-sqrt(70._gp)*a**2,sqrt(17.5_gp)*a**2,-sqrt(157.5_gp)*a**2,&
               sqrt(17.5_gp)*a**2,sqrt(3189.375_gp)*a,-sqrt(354.375_gp)*a/) ! m=5, l=3, n=1
          ftpd(42:47)=rsp*(/sqrt(105._gp)*a**2,-sqrt(105._gp)*a**2,sqrt(105._gp)*a**2,-sqrt(105._gp)*a**2,&
               -sqrt(2126.25_gp)*a,sqrt(2126.25_gp)*a/) ! m=6, l=3, n=1
          ftpd(48:51)=rsp*(/sqrt(420._gp)*a**2,sqrt(420._gp)*a**2,sqrt(420._gp)*a**2,-sqrt(8505._gp)*a/) ! m=7, l=3, n=1
       end select
    case(2)
       select case(l)
       case(0)
          ntpd_shell=10 ! l=0,n=2
          ntpd(1:1)=(/10/) ! l=0,n=2
          pow(1:3,1:10)=reshape((/0,0,0,4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4,2,0,0,0,2,0,0,0,2/),(/3,10/)) ! l=0,n=2
          ftpd(1:10)=rsp*(/3._gp,2._gp*a**2,4._gp*a**2,2._gp*a**2,4._gp*a**2,4._gp*a**2,2._gp*a**2,&
               -7._gp*a,-7._gp*a,-7._gp*a/) ! m=1, l=0, n=2
       case(1)
          ntpd_shell=30 ! l=1,n=2
          ntpd(1:3)=(/10,10,10/) ! l=1,n=2
          pow(1:3,1:30)=reshape((/1,0,0,5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,3,0,0,1,2,0,1,0,2,0,1,0,4,1,0,2,3,&
               0,0,5,0,2,1,2,0,3,2,0,1,4,2,1,0,0,3,0,0,1,2,0,0,1,4,0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5,2,0,1,0,2,&
               1,0,0,3/),(/3,30/)) ! l=1,n=2
          ftpd(1:10)=rsp*(/sqrt(75._gp),sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,&
               sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(243._gp)*a,-sqrt(243._gp)*a,-sqrt(243._gp)*a/) ! m=1, l=1, n=2
          ftpd(11:20)=rsp*(/sqrt(75._gp),sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,&
               sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(243._gp)*a,-sqrt(243._gp)*a,-sqrt(243._gp)*a/) ! m=2, l=1, n=2
          ftpd(21:30)=rsp*(/sqrt(75._gp),sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,&
               sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(243._gp)*a,-sqrt(243._gp)*a,-sqrt(243._gp)*a/) ! m=3, l=1, n=2
       case(2)
          ntpd_shell=60 ! l=2,n=2
          ntpd(1:5)=(/10,10,10,14,16/) ! l=2,n=2
          pow(1:3,1:60)=reshape((/0,1,1,4,1,1,2,3,1,0,5,1,2,1,3,0,3,3,0,1,5,2,1,1,0,3,1,0,1,3,1,0,1,5,0,1,3,2,1,1,4,1,3,&
               0,3,1,2,3,1,0,5,3,0,1,1,2,1,1,0,3,1,1,0,5,1,0,3,3,0,1,5,0,3,1,2,1,3,2,1,1,4,3,1,0,1,3,0,1,1,2,2,0,0,0,2,&
               0,6,0,0,4,2,0,2,4,0,0,6,0,4,0,2,0,4,2,2,0,4,0,2,4,4,0,0,0,4,0,2,0,2,0,2,2,2,0,0,0,2,0,0,0,2,6,0,0,4,2,0,2,&
               4,0,0,6,0,2,0,4,0,2,4,0,0,6,4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4/),(/3,60/)) ! l=2,n=2
          ftpd(1:10)=rsp*(/sqrt(735._gp),sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,&
               sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a/) ! m=1, l=2, n=2
          ftpd(11:20)=rsp*(/sqrt(735._gp),sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,&
               sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a/) ! m=2, l=2, n=2
          ftpd(21:30)=rsp*(/sqrt(735._gp),sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,&
               sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a/) ! m=3, l=2, n=2
          ftpd(31:44)=rsp*(/sqrt(183.75_gp),-sqrt(183.75_gp),sqrt(15._gp)*a**2,sqrt(15._gp)*a**2,-sqrt(15._gp)*a**2,&
               -sqrt(15._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(60._gp)*a**2,sqrt(15._gp)*a**2,-sqrt(15._gp)*a**2,&
               -sqrt(453.75_gp)*a,sqrt(453.75_gp)*a,-sqrt(453.75_gp)*a,sqrt(453.75_gp)*a/) ! m=4, l=2, n=2
          ftpd(45:60)=rsp*(/-sqrt(61.25_gp),-sqrt(61.25_gp),sqrt(245._gp),-sqrt(5._gp)*a**2,-sqrt(45._gp)*a**2,&
               -sqrt(45._gp)*a**2,-sqrt(5._gp)*a**2,sqrt(45._gp)*a**2,sqrt(45._gp)*a**2,sqrt(20._gp)*a**2,&
               sqrt(151.25_gp)*a,sqrt(605._gp)*a,sqrt(151.25_gp)*a,-sqrt(151.25_gp)*a,-sqrt(151.25_gp)*a,&
               -sqrt(605._gp)*a/) ! m=5, l=2, n=2
       case(3)
          ntpd_shell=113 ! l=3,n=2
          ntpd(1:7)=(/19,19,19,16,16,14,10/) ! l=3,n=2
          pow(1:3,1:113)=reshape((/3,0,0,1,2,0,1,0,2,7,0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,4,2,3,0,4,1,2,4,1,0,6,&
               5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,2,1,0,0,3,0,0,1,2,6,1,0,4,3,0,2,5,0,0,7,0,4,1,2,2,3,2,0,5,2,2,&
               1,4,0,3,4,0,1,6,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,0,1,4,2,0,1,0,2,1,0,0,3,6,0,1,4,2,1,2,4,1,0,6,1,4,0,&
               3,2,2,3,0,4,3,2,0,5,0,2,5,0,0,7,4,0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5,3,0,0,1,2,0,7,0,0,5,2,0,3,4,0,&
               1,6,0,5,0,2,3,2,2,1,4,2,3,0,4,1,2,4,5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,2,1,0,0,3,0,6,1,0,4,3,0,2,5,0,0,&
               7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,3,4,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,2,0,1,0,2,1,6,0,1,4,2,1,2,4,1,0,6,&
               1,4,0,3,0,4,3,2,0,5,0,2,5,4,0,1,0,4,1,2,0,3,0,2,3,1,1,1,5,1,1,3,3,1,1,5,1,3,1,3,1,3,3,1,1,5,3,1,1,&
               1,3,1,1,1,3/),(/3,113/)) ! l=3,n=2
          ftpd(1:19)=rsp*(/sqrt(212.625_gp),sqrt(212.625_gp),-sqrt(3402._gp),sqrt(10.5_gp)*a**2,sqrt(94.5_gp)*a**2,&
               sqrt(94.5_gp)*a**2,sqrt(10.5_gp)*a**2,-sqrt(42._gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(42._gp)*a**2,&
               -sqrt(514.5_gp)*a**2,-sqrt(514.5_gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(443.625_gp)*a,-sqrt(1774.5_gp)*a,&
               -sqrt(443.625_gp)*a,sqrt(3992.625_gp)*a,sqrt(3992.625_gp)*a,sqrt(7098._gp)*a/) ! m=1, l=3, n=2
          ftpd(20:38)=rsp*(/sqrt(212.625_gp),sqrt(212.625_gp),-sqrt(3402._gp),sqrt(10.5_gp)*a**2,sqrt(94.5_gp)*a**2,&
               sqrt(94.5_gp)*a**2,sqrt(10.5_gp)*a**2,-sqrt(42._gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(42._gp)*a**2,&
               -sqrt(514.5_gp)*a**2,-sqrt(514.5_gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(443.625_gp)*a,-sqrt(1774.5_gp)*a,&
               -sqrt(443.625_gp)*a,sqrt(3992.625_gp)*a,sqrt(3992.625_gp)*a,sqrt(7098._gp)*a/) ! m=2, l=3, n=2
          ftpd(39:57)=rsp*(/sqrt(1275.75_gp),sqrt(1275.75_gp),-sqrt(567._gp),sqrt(63._gp)*a**2,sqrt(567._gp)*a**2,&
               sqrt(567._gp)*a**2,sqrt(63._gp)*a**2,sqrt(112._gp)*a**2,sqrt(448._gp)*a**2,sqrt(112._gp)*a**2,&
               -sqrt(7._gp)*a**2,-sqrt(7._gp)*a**2,-sqrt(28._gp)*a**2,-sqrt(2661.75_gp)*a,-sqrt(10647._gp)*a,&
               -sqrt(2661.75_gp)*a,-sqrt(295.75_gp)*a,-sqrt(295.75_gp)*a,sqrt(1183._gp)*a/) ! m=3, l=3, n=2
          ftpd(58:73)=rsp*(/sqrt(354.375_gp),-sqrt(3189.375_gp),sqrt(17.5_gp)*a**2,-sqrt(17.5_gp)*a**2,&
               -sqrt(437.5_gp)*a**2,-sqrt(157.5_gp)*a**2,sqrt(70._gp)*a**2,-sqrt(280._gp)*a**2,-sqrt(630._gp)*a**2,&
               sqrt(17.5_gp)*a**2,-sqrt(157.5_gp)*a**2,-sqrt(739.375_gp)*a,sqrt(2957.5_gp)*a,sqrt(6654.375_gp)*a,&
               -sqrt(739.375_gp)*a,sqrt(6654.375_gp)*a/) ! m=4, l=3, n=2
          ftpd(74:89)=rsp*(/-sqrt(3189.375_gp),sqrt(354.375_gp),-sqrt(157.5_gp)*a**2,-sqrt(437.5_gp)*a**2,-sqrt(17.5_gp)*a**2,&
               sqrt(17.5_gp)*a**2,-sqrt(630._gp)*a**2,-sqrt(280._gp)*a**2,sqrt(70._gp)*a**2,-sqrt(157.5_gp)*a**2,&
               sqrt(17.5_gp)*a**2,sqrt(6654.375_gp)*a,sqrt(2957.5_gp)*a,-sqrt(739.375_gp)*a,sqrt(6654.375_gp)*a,&
               -sqrt(739.375_gp)*a/) ! m=5, l=3, n=2
          ftpd(90:103)=rsp*(/sqrt(2126.25_gp),-sqrt(2126.25_gp),sqrt(105._gp)*a**2,sqrt(105._gp)*a**2,-sqrt(105._gp)*a**2,&
               -sqrt(105._gp)*a**2,sqrt(420._gp)*a**2,-sqrt(420._gp)*a**2,sqrt(105._gp)*a**2,-sqrt(105._gp)*a**2,&
               -sqrt(4436.25_gp)*a,sqrt(4436.25_gp)*a,-sqrt(4436.25_gp)*a,sqrt(4436.25_gp)*a/) ! m=6, l=3, n=2
          ftpd(104:113)=rsp*(/sqrt(8505._gp),sqrt(420._gp)*a**2,sqrt(1680._gp)*a**2,sqrt(420._gp)*a**2,sqrt(1680._gp)*a**2,&
               sqrt(1680._gp)*a**2,sqrt(420._gp)*a**2,-sqrt(17745._gp)*a,-sqrt(17745._gp)*a,-sqrt(17745._gp)*a/) ! m=7, l=3, n=2
       end select
    case(3)
       select case(l)
       case(0)
          ntpd_shell=19 ! l=0,n=3
          ntpd(1:1)=(/19/) ! l=0,n=3
          pow(1:3,1:19)=reshape((/2,0,0,0,2,0,0,0,2,6,0,0,4,2,0,2,4,0,0,6,0,4,0,2,2,2,2,0,4,2,2,0,4,0,2,4,0,0,6,4,&
               0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4/),(/3,19/)) ! l=0,n=3
          ftpd(1:19)=rsp*(/10._gp,10._gp,10._gp,2._gp*a**2,6._gp*a**2,6._gp*a**2,2._gp*a**2,6._gp*a**2,12._gp*a**2,&
               6._gp*a**2,6._gp*a**2,6._gp*a**2,2._gp*a**2,-11._gp*a,-22._gp*a,-11._gp*a,-22._gp*a,-22._gp*a,&
               -11._gp*a/) ! m=1, l=0, n=3
       case(1)
          ntpd_shell=57 ! l=1,n=3
          ntpd(1:3)=(/19,19,19/) ! l=1,n=3
          pow(1:3,1:57)=reshape((/3,0,0,1,2,0,1,0,2,7,0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,4,2,3,0,4,1,2,4,1,0,6,5,0,0,&
               3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,2,1,0,0,3,0,0,1,2,6,1,0,4,3,0,2,5,0,0,7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,3,4,&
               0,1,6,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,0,1,4,2,0,1,0,2,1,0,0,3,6,0,1,4,2,1,2,4,1,0,6,1,4,0,3,2,2,3,0,4,3,2,&
               0,5,0,2,5,0,0,7,4,0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5/),(/3,57/)) ! l=1,n=3
          ftpd(1:19)=rsp*(/sqrt(588._gp),sqrt(588._gp),sqrt(588._gp),sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,sqrt(432._gp)*a**2,sqrt(108._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(108._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(507._gp)*a,-sqrt(2028._gp)*a,-sqrt(507._gp)*a,-sqrt(2028._gp)*a,&
               -sqrt(2028._gp)*a,-sqrt(507._gp)*a/) ! m=1, l=1, n=3
          ftpd(20:38)=rsp*(/sqrt(588._gp),sqrt(588._gp),sqrt(588._gp),sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,sqrt(432._gp)*a**2,sqrt(108._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(108._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(507._gp)*a,-sqrt(2028._gp)*a,-sqrt(507._gp)*a,-sqrt(2028._gp)*a,&
               -sqrt(2028._gp)*a,-sqrt(507._gp)*a/) ! m=2, l=1, n=3
          ftpd(39:57)=rsp*(/sqrt(588._gp),sqrt(588._gp),sqrt(588._gp),sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(108._gp)*a**2,sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,sqrt(432._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(108._gp)*a**2,sqrt(108._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(507._gp)*a,-sqrt(2028._gp)*a,-sqrt(507._gp)*a,&
               -sqrt(2028._gp)*a,-sqrt(2028._gp)*a,-sqrt(507._gp)*a/) ! m=3, l=1, n=3
       case(2)
          ntpd_shell=109 ! l=2,n=3
          ntpd(1:5)=(/19,19,19,24,28/) ! l=2,n=3
          pow(1:3,1:109)=reshape((/2,1,1,0,3,1,0,1,3,6,1,1,4,3,1,2,5,1,0,7,1,4,1,3,2,3,3,0,5,3,2,1,5,0,3,5,0,1,7,&
               4,1,1,2,3,1,0,5,1,2,1,3,0,3,3,0,1,5,3,0,1,1,2,1,1,0,3,7,0,1,5,2,1,3,4,1,1,6,1,5,0,3,3,2,3,1,4,3,3,&
               0,5,1,2,5,1,0,7,5,0,1,3,2,1,1,4,1,3,0,3,1,2,3,1,0,5,3,1,0,1,3,0,1,1,2,7,1,0,5,3,0,3,5,0,1,7,0,5,1,&
               2,3,3,2,1,5,2,3,1,4,1,3,4,1,1,6,5,1,0,3,3,0,1,5,0,3,1,2,1,3,2,1,1,4,4,0,0,0,4,0,2,0,2,0,2,2,8,0,0,&
               6,2,0,2,6,0,0,8,0,6,0,2,4,2,2,2,4,2,0,6,2,4,0,4,0,4,4,2,0,6,0,2,6,6,0,0,4,2,0,2,4,0,0,6,0,4,0,2,0,&
               4,2,2,0,4,0,2,4,4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4,8,0,0,6,2,0,4,4,0,2,6,0,0,8,0,6,0,2,4,2,2,2,4,&
               2,0,6,2,4,0,4,2,2,4,0,4,4,2,0,6,0,2,6,0,0,8,6,0,0,4,2,0,2,4,0,0,6,0,2,0,4,0,2,4,0,0,6/),(/3,109/)) ! l=2,n=3
          ftpd(1:19)=rsp*(/sqrt(4860._gp),sqrt(4860._gp),sqrt(4860._gp),sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(540._gp)*a**2,sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,sqrt(2160._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(540._gp)*a**2,sqrt(540._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(3375._gp)*a,-sqrt(13500._gp)*a,&
               -sqrt(3375._gp)*a,-sqrt(13500._gp)*a,-sqrt(13500._gp)*a,-sqrt(3375._gp)*a/) ! m=1, l=2, n=3
          ftpd(20:38)=rsp*(/sqrt(4860._gp),sqrt(4860._gp),sqrt(4860._gp),sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(540._gp)*a**2,sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,sqrt(2160._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(540._gp)*a**2,sqrt(540._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(3375._gp)*a,-sqrt(13500._gp)*a,&
               -sqrt(3375._gp)*a,-sqrt(13500._gp)*a,-sqrt(13500._gp)*a,-sqrt(3375._gp)*a/) ! m=2, l=2, n=3
          ftpd(39:57)=rsp*(/sqrt(4860._gp),sqrt(4860._gp),sqrt(4860._gp),sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,sqrt(2160._gp)*a**2,sqrt(540._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(540._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(3375._gp)*a,-sqrt(13500._gp)*a,-sqrt(3375._gp)*a,-sqrt(13500._gp)*a,&
               -sqrt(13500._gp)*a,-sqrt(3375._gp)*a/) ! m=3, l=2, n=3
          ftpd(58:81)=rsp*(/sqrt(1215._gp),-sqrt(1215._gp),sqrt(1215._gp),-sqrt(1215._gp),sqrt(15._gp)*a**2,sqrt(60._gp)*a**2,&
               -sqrt(60._gp)*a**2,-sqrt(15._gp)*a**2,sqrt(135._gp)*a**2,sqrt(135._gp)*a**2,-sqrt(135._gp)*a**2,&
               -sqrt(135._gp)*a**2,sqrt(135._gp)*a**2,-sqrt(135._gp)*a**2,sqrt(15._gp)*a**2,-sqrt(15._gp)*a**2,&
               -sqrt(843.75_gp)*a,-sqrt(843.75_gp)*a,sqrt(843.75_gp)*a,sqrt(843.75_gp)*a,-sqrt(3375._gp)*a,sqrt(3375._gp)*a,&
               -sqrt(843.75_gp)*a,sqrt(843.75_gp)*a/) ! m=4, l=2, n=3
          ftpd(82:109)=rsp*(/-sqrt(405._gp),-sqrt(1620._gp),-sqrt(405._gp),sqrt(405._gp),sqrt(405._gp),sqrt(1620._gp),&
               -sqrt(5._gp)*a**2,-sqrt(80._gp)*a**2,-sqrt(180._gp)*a**2,-sqrt(80._gp)*a**2,-sqrt(5._gp)*a**2,&
               -sqrt(5._gp)*a**2,-sqrt(45._gp)*a**2,-sqrt(45._gp)*a**2,-sqrt(5._gp)*a**2,sqrt(45._gp)*a**2,&
               sqrt(180._gp)*a**2,sqrt(45._gp)*a**2,sqrt(125._gp)*a**2,sqrt(125._gp)*a**2,sqrt(20._gp)*a**2,&
               sqrt(281.25_gp)*a,sqrt(2531.25_gp)*a,sqrt(2531.25_gp)*a,sqrt(281.25_gp)*a,-sqrt(2531.25_gp)*a,&
               -sqrt(2531.25_gp)*a,-sqrt(1125._gp)*a/) ! m=5, l=2, n=3
       case(3)
          ntpd_shell=190 ! l=3,n=3
          ntpd(1:7)=(/31,31,31,27,27,24,19/) ! l=3,n=3
          pow(1:3,1:190)=reshape((/5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,9,0,0,7,2,0,5,4,0,3,6,0,1,8,0,7,0,2,5,&
               2,2,3,4,2,1,6,2,5,0,4,3,2,4,1,4,4,3,0,6,1,2,6,1,0,8,7,0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,4,2,3,&
               0,4,1,2,4,1,0,6,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,0,1,4,8,1,0,6,3,0,4,5,0,2,7,0,0,9,0,6,1,2,4,3,2,2,5,2,&
               0,7,2,4,1,4,2,3,4,0,5,4,2,1,6,0,3,6,0,1,8,6,1,0,4,3,0,2,5,0,0,7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,3,4,0,1,6,4,&
               0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5,8,0,1,6,2,1,4,4,1,2,6,1,0,8,1,6,0,3,4,2,3,2,4,3,0,6,3,4,0,5,2,2,5,0,4,5,2,&
               0,7,0,2,7,0,0,9,6,0,1,4,2,1,2,4,1,0,6,1,4,0,3,2,2,3,0,4,3,2,0,5,0,2,5,0,0,7,5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,9,0,&
               0,5,4,0,3,6,0,1,8,0,7,0,2,5,2,2,3,4,2,1,6,2,5,0,4,3,2,4,1,4,4,3,0,6,1,2,6,7,0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,&
               4,2,3,0,4,1,2,4,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,8,1,0,6,3,0,4,5,0,0,9,0,6,1,2,4,3,2,2,5,2,0,7,2,4,1,4,2,3,4,0,5,4,2,&
               1,6,0,3,6,6,1,0,4,3,0,2,5,0,0,7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,3,4,4,0,1,0,4,1,2,0,3,0,2,3,8,0,1,6,2,1,2,6,1,0,8,1,6,0,&
               3,4,2,3,2,4,3,0,6,3,4,0,5,0,4,5,2,0,7,0,2,7,6,0,1,4,2,1,2,4,1,0,6,1,4,0,3,0,4,3,2,0,5,0,2,5,3,1,1,1,3,1,1,1,3,7,1,1,&
               5,3,1,3,5,1,1,7,1,5,1,3,3,3,3,1,5,3,3,1,5,1,3,5,1,1,7,5,1,1,3,3,1,1,5,1,3,1,3,1,3,3,1,1,5/),(/3,190/)) ! l=3,n=3
          ftpd(1:31)=rsp*(/sqrt(1270.5_gp),sqrt(5082._gp),sqrt(1270.5_gp),-sqrt(11434.5_gp),-sqrt(11434.5_gp),&
               -sqrt(20328._gp),sqrt(10.5_gp)*a**2,sqrt(168._gp)*a**2,sqrt(378._gp)*a**2,sqrt(168._gp)*a**2,&
               sqrt(10.5_gp)*a**2,-sqrt(10.5_gp)*a**2,-sqrt(94.5_gp)*a**2,-sqrt(94.5_gp)*a**2,-sqrt(10.5_gp)*a**2,&
               -sqrt(850.5_gp)*a**2,-sqrt(3402._gp)*a**2,-sqrt(850.5_gp)*a**2,-sqrt(1270.5_gp)*a**2,-sqrt(1270.5_gp)*a**2,&
               -sqrt(168._gp)*a**2,-sqrt(758.625_gp)*a,-sqrt(6827.625_gp)*a,-sqrt(6827.625_gp)*a,-sqrt(758.625_gp)*a,&
               sqrt(3034.5_gp)*a,sqrt(12138._gp)*a,sqrt(3034.5_gp)*a,sqrt(37172.625_gp)*a,sqrt(37172.625_gp)*a,&
               sqrt(12138._gp)*a/) ! m=1, l=3, n=3
          ftpd(32:62)=rsp*(/sqrt(1270.5_gp),sqrt(5082._gp),sqrt(1270.5_gp),-sqrt(11434.5_gp),-sqrt(11434.5_gp),&
               -sqrt(20328._gp),sqrt(10.5_gp)*a**2,sqrt(168._gp)*a**2,sqrt(378._gp)*a**2,sqrt(168._gp)*a**2,sqrt(10.5_gp)*a**2,&
               -sqrt(10.5_gp)*a**2,-sqrt(94.5_gp)*a**2,-sqrt(94.5_gp)*a**2,-sqrt(10.5_gp)*a**2,-sqrt(850.5_gp)*a**2,&
               -sqrt(3402._gp)*a**2,-sqrt(850.5_gp)*a**2,-sqrt(1270.5_gp)*a**2,-sqrt(1270.5_gp)*a**2,-sqrt(168._gp)*a**2,&
               -sqrt(758.625_gp)*a,-sqrt(6827.625_gp)*a,-sqrt(6827.625_gp)*a,-sqrt(758.625_gp)*a,sqrt(3034.5_gp)*a,&
               sqrt(12138._gp)*a,sqrt(3034.5_gp)*a,sqrt(37172.625_gp)*a,sqrt(37172.625_gp)*a,sqrt(12138._gp)*a/) ! m=2, l=3, n=3
          ftpd(63:93)=rsp*(/sqrt(7623._gp),sqrt(30492._gp),sqrt(7623._gp),sqrt(847._gp),sqrt(847._gp),-sqrt(3388._gp),&
               sqrt(63._gp)*a**2,sqrt(1008._gp)*a**2,sqrt(2268._gp)*a**2,sqrt(1008._gp)*a**2,sqrt(63._gp)*a**2,sqrt(343._gp)*a**2,&
               sqrt(3087._gp)*a**2,sqrt(3087._gp)*a**2,sqrt(343._gp)*a**2,sqrt(63._gp)*a**2,sqrt(252._gp)*a**2,sqrt(63._gp)*a**2,&
               -sqrt(63._gp)*a**2,-sqrt(63._gp)*a**2,-sqrt(28._gp)*a**2,-sqrt(4551.75_gp)*a,-sqrt(40965.75_gp)*a,&
               -sqrt(40965.75_gp)*a,-sqrt(4551.75_gp)*a,-sqrt(8092._gp)*a,-sqrt(32368._gp)*a,-sqrt(8092._gp)*a,&
               sqrt(505.75_gp)*a,sqrt(505.75_gp)*a,&
               sqrt(2023._gp)*a/) ! m=3, l=3, n=3
          ftpd(94:120)=rsp*(/sqrt(2117.5_gp),-sqrt(8470._gp),-sqrt(19057.5_gp),sqrt(2117.5_gp),-sqrt(19057.5_gp),&
               sqrt(17.5_gp)*a**2,-sqrt(630._gp)*a**2,-sqrt(1120._gp)*a**2,-sqrt(157.5_gp)*a**2,sqrt(157.5_gp)*a**2,&
               -sqrt(157.5_gp)*a**2,-sqrt(3937.5_gp)*a**2,-sqrt(1417.5_gp)*a**2,sqrt(157.5_gp)*a**2,&
               -sqrt(630._gp)*a**2,-sqrt(1417.5_gp)*a**2,sqrt(17.5_gp)*a**2,-sqrt(157.5_gp)*a**2,&
               -sqrt(1264.375_gp)*a,sqrt(1264.375_gp)*a,sqrt(31609.375_gp)*a,sqrt(11379.375_gp)*a,&
               -sqrt(5057.5_gp)*a,sqrt(20230._gp)*a,sqrt(45517.5_gp)*a,-sqrt(1264.375_gp)*a,sqrt(11379.375_gp)*a/) ! m=4, l=3, n=3
          ftpd(121:147)=rsp*(/-sqrt(19057.5_gp),-sqrt(8470._gp),sqrt(2117.5_gp),-sqrt(19057.5_gp),sqrt(2117.5_gp),&
               -sqrt(157.5_gp)*a**2,-sqrt(1120._gp)*a**2,-sqrt(630._gp)*a**2,sqrt(17.5_gp)*a**2,-sqrt(1417.5_gp)*a**2,&
               -sqrt(3937.5_gp)*a**2,-sqrt(157.5_gp)*a**2,sqrt(157.5_gp)*a**2,-sqrt(1417.5_gp)*a**2,-sqrt(630._gp)*a**2,&
               sqrt(157.5_gp)*a**2,-sqrt(157.5_gp)*a**2,sqrt(17.5_gp)*a**2,sqrt(11379.375_gp)*a,sqrt(31609.375_gp)*a,&
               sqrt(1264.375_gp)*a,-sqrt(1264.375_gp)*a,sqrt(45517.5_gp)*a,sqrt(20230._gp)*a,-sqrt(5057.5_gp)*a,&
               sqrt(11379.375_gp)*a,-sqrt(1264.375_gp)*a/) ! m=5, l=3, n=3
          ftpd(148:171)=rsp*(/sqrt(12705._gp),-sqrt(12705._gp),sqrt(12705._gp),-sqrt(12705._gp),sqrt(105._gp)*a**2,&
               sqrt(420._gp)*a**2,-sqrt(420._gp)*a**2,-sqrt(105._gp)*a**2,sqrt(945._gp)*a**2,sqrt(945._gp)*a**2,&
               -sqrt(945._gp)*a**2,-sqrt(945._gp)*a**2,sqrt(945._gp)*a**2,-sqrt(945._gp)*a**2,sqrt(105._gp)*a**2,&
               -sqrt(105._gp)*a**2,-sqrt(7586.25_gp)*a,-sqrt(7586.25_gp)*a,sqrt(7586.25_gp)*a,sqrt(7586.25_gp)*a,&
               -sqrt(30345._gp)*a,sqrt(30345._gp)*a,-sqrt(7586.25_gp)*a,sqrt(7586.25_gp)*a/) ! m=6, l=3, n=3
          ftpd(172:190)=rsp*(/sqrt(50820._gp),sqrt(50820._gp),sqrt(50820._gp),sqrt(420._gp)*a**2,sqrt(3780._gp)*a**2,&
               sqrt(3780._gp)*a**2,sqrt(420._gp)*a**2,sqrt(3780._gp)*a**2,sqrt(15120._gp)*a**2,sqrt(3780._gp)*a**2,&
               sqrt(3780._gp)*a**2,sqrt(3780._gp)*a**2,sqrt(420._gp)*a**2,-sqrt(30345._gp)*a,-sqrt(121380._gp)*a,&
               -sqrt(30345._gp)*a,-sqrt(121380._gp)*a,-sqrt(121380._gp)*a,-sqrt(30345._gp)*a/) ! m=7, l=3, n=3
       end select
    end select
  end subroutine tensor_product_decomposition_laplacian

end module gaussians
