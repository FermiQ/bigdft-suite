!> @file
!! Multipole preserving quadrature scheme
!! @author
!!    Copyright (C) 2016-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module multipole_preserving

  use dynamic_memory
  use f_precisions, only: gp => f_double
  use f_arrays
  use dictionaries, only: f_err_throw
  use yaml_output, only: yaml_mapping_open,yaml_mapping_close,yaml_map
  implicit none

  private

  integer :: itype_scf=0   !< Type of the interpolating SCF, 0= data unallocated
  integer :: n_scf=-1      !< Number of points of the allocated data
  integer :: nrange_scf=0  !< Range of the integration
  real(gp), dimension(:), allocatable :: scf_data !< Values for the interpolating scaling functions points
  type(f_matrix), dimension(3), save :: mp_exps !< Refcounted arrays to precalculate the coefficients

  !> log of the minimum value of the scf data
  !! to avoid floating point exceptions while multiplying with it
  real(gp) :: mn_scf = 0.0_gp
  real(gp), parameter :: mp_tiny = 1.e-30_gp !<put zero when the value is lower than this


  public :: initialize_real_space_conversion,finalize_real_space_conversion,scfdotf,mp_exp
  public :: mp_gaussian_workarrays,get_mp_exps_product,mp_initialized

  contains

    !> Prepare the array for the evaluation with the interpolating Scaling Functions
    !! one might add also the function to be converted and the
    !! prescription for integrating knowing the scaling relation of the function
    subroutine initialize_real_space_conversion(npoints,isf_m,rloc,nmoms,verbose)
      implicit none
      !Arguments
      integer, intent(in), optional :: npoints !< Number of points (only 2**x)
      !> rloc of a given set of gaussian functions in order to determine npoints
      real(kind=8), dimension(:), intent(in), optional :: rloc
      integer, intent(in), optional :: isf_m !< Type of interpolatig scaling function
      integer, intent(in), optional :: nmoms !< Number of preserved moments if /= 0
                                             !! (see ISF_family in scaling_function.f90)
      logical, intent(in), optional :: verbose !< If .true., display some information
      !Local variables
      character(len=*), parameter :: subname='initialize_real_space_conversion'
      integer :: n_range,i,nmm,np
      real(gp) :: tt,r
      real(gp), dimension(:), allocatable :: x_scf !< to be removed in a future implementation

      !Check if already allocated
      if (itype_scf /= 0) &
         call f_err_throw('initialize_real_space_conversion already called.')

      if (present(isf_m)) then
         itype_scf=isf_m
      else
         itype_scf=16
      end if

      if (present(nmoms)) then
         nmm=nmoms
      else
         nmm=0
      end if

      np = 6
      !Determine the length of the array
      if (present(npoints)) then
         np=ceiling(log(dble(npoints))/log(2.d0))
      else if (present(rloc)) then
         r = minval(rloc)
         if (r /= 0.d0) then
           np=ceiling(log(10.d0/r)/log(2.0))
         else
           np=10
         end if
      end if
      np=max(6,np)
      !2**6 is a min to have 2**12=2048 points
      !otherwise trouble in scaling_function.f90
      n_scf=2*(itype_scf+nmm)*(2**np)

      !allocations for scaling function data array
      x_scf = f_malloc(0.to.n_scf,id='x_scf')

      scf_data = f_malloc(0.to.n_scf,id='scf_data')

      !Build the scaling function external routine coming from Poisson Solver. To be customized accordingly
      !call scaling_function(itype_scf,n_scf,n_range,x_scf,scf_data)
      !call wavelet_function(itype_scf,n_scf,x_scf,scf_data)
      call ISF_family(itype_scf,nmm,n_scf,n_range,x_scf,scf_data)
      call f_free(x_scf)

      if (present(verbose)) then
        if (verbose) then
          call yaml_mapping_open('Multipole preserving (ionic potential)',flow=.true.)
          if (present(rloc)) call yaml_map('rloc',r)
          call yaml_map('itype_scf',itype_scf)
          call yaml_map('npoints',n_scf)
          call yaml_mapping_close()
        end if
      end if

      nrange_scf=n_range
      !Define the log of the smallest non zero value as the  cutoff for multiplying with it.
      !This means that the values which are
      !lower than scf_data squared will be considered as zero.
      mn_scf=epsilon(1.d0)**2 !just to put a "big" value
      do i=0,n_scf
         tt=scf_data(i)
         if (tt /= 0.0_gp .and. abs(tt) > sqrt(tiny(1.0))) then
            mn_scf=min(mn_scf,tt**2)
         end if
      end do

    end subroutine initialize_real_space_conversion


    !> Check if scf_data is initialized
    pure function mp_initialized()
      implicit none
      logical :: mp_initialized
      mp_initialized=allocated(scf_data)
    end function mp_initialized


    !> Build the work arrays
    subroutine mp_gaussian_workarrays(nterms,nbox,expo,lxyz,rxyz,hgrids)
      implicit none
      !Arguments
      integer, intent(in) :: nterms
      real(gp), intent(in) :: expo
      real(gp), dimension(3), intent(in) :: rxyz,hgrids
      integer, dimension(2,3), intent(in) :: nbox
      integer, dimension(nterms,3), intent(in) :: lxyz
      !Local variables
      integer :: i,iterm,it
      !matrix to be added to allocate the arrays
      do i=1,3
         mp_exps(i)=f_malloc_ptr([1.to.nterms,nbox(1,i).to.nbox(2,i)],id='mp_Exp_i')
         do it=nbox(1,i),nbox(2,i)
            do iterm=1,nterms
               mp_exps(i)%ptr(iterm,it)=&
                    mp_exp(hgrids(i),rxyz(i),expo,it-1,lxyz(iterm,i),.true.)
            end do
         end do
      end do

    end subroutine mp_gaussian_workarrays


    pure function get_mp_exps_product(nterms,ival,factors) result(f)
      implicit none
      integer, intent(in) :: nterms
      real(gp), dimension(nterms), intent(in) :: factors
      integer, dimension(3), intent(in) :: ival
      real(gp) :: f
      !local variables
      integer :: iterm
      f=0.0_gp
      do iterm=1,nterms
         f=f+mp_exps(1)%ptr(iterm,ival(1))*&
              mp_exps(2)%ptr(iterm,ival(2))*&
              mp_exps(3)%ptr(iterm,ival(3))*factors(iterm)
      end do
    end function get_mp_exps_product


    !> Deallocate scf_data
    subroutine finalize_real_space_conversion()
      implicit none

      itype_scf=0 !< To indicate need initialization
      n_scf=-1
      mn_scf=0.0_gp
      nrange_scf=0
      call f_free(scf_data)
      call f_array_free(mp_exps)

    end subroutine finalize_real_space_conversion


    !> multipole-preserving gaussian function
    !! chooses between traditional exponential and scfdotf
    !! according to the value of the exponent in units of the grid spacing
    !! the function is supposed to be x**pow*exp(-expo*x**2)
    !! where x=hgrid*j-x0
    !! @warning
    !! this function is also elemental to ease its evaluation, though
    !! the usage for vector argument is discouraged: dedicated routines has to be
    !! written to meet performance
    !! @todo  Optimize it!
    elemental pure function mp_exp(hgrid,x0,expo,j,pow,modified)
      use numerics, only: safe_exp
      implicit none
      real(gp), intent(in) :: hgrid   !< Hgrid
      real(gp), intent(in) :: x0      !< X value
      real(gp), intent(in) :: expo    !< Exponent of the gaussian
      logical, intent(in) :: modified !< Switch to scfdotf if true
      integer, intent(in) :: j        !< Location of the scf from x0
      integer, intent(in) :: pow      !< Exp(-expo*x**2)*(x**pow)
      real(gp) :: mp_exp
      !local variables
      real(gp) :: x

      !added failsafe to avoid segfaults
      if (modified .and. allocated(scf_data)) then
         mp_exp=scfdotf(j,hgrid,expo,x0,pow)
      else
         x=hgrid*j-x0
         mp_exp=safe_exp(-expo*x**2,underflow=mp_tiny)
         if (pow /= 0) mp_exp=mp_exp*(x**pow)
      end if
    end function mp_exp


    !> This function calculates the scalar product between a ISF and a
    !! input function, which is a gaussian times a power centered
    !! @f$g(x) = (x-x_0)^{pow} e^{-pgauss (x-x_0)}@f$
    !! here pure specifier is redundant
    !! we should add here the threshold from which the
    !! normal function can be evaluated
    elemental pure function scfdotf(j,hgrid,pgauss,x0,pow) result(gint)
      use numerics, only: safe_exp
      implicit none
      !Arguments
      integer, intent(in) :: j !<value of the input result in the hgrid reference
      integer, intent(in) :: pow
      real(gp), intent(in) :: hgrid,pgauss,x0
      real(gp) :: gint
      !local variables
      integer :: i
      real(gp) :: x,absci,fabsci,dx
      gint=0.0_gp

      !Step grid for the integration
      !dx = real(2*itype_scf,gp)/real(n_scf,gp)
      dx = real(nrange_scf,gp)/real(n_scf,gp)
      !starting point for the x coordinate for integration
      !x  = real(j-itype_scf+1,gp)-dx
      x  = real(j-nrange_scf/2+1,gp)-dx

      !the loop can be unrolled to maximize performances
      if (pow /= 0) then
         do i=0,n_scf
            x=x+dx
            absci = x*hgrid - x0
            !here evaluate the function
            fabsci = absci**pow
            absci = -pgauss*absci*absci
            fabsci = fabsci*safe_exp(absci)!,underflow=mn_scf)
            !calculate the integral
            gint = gint + scf_data(i)*fabsci
            !       print *,'test',i,scf_data(i),fabsci,pgauss,pow,absci
         end do
      else
         do i=0,n_scf
            x=x+dx
            absci = x*hgrid - x0
            !          !here evaluate the function
            absci = -pgauss*absci*absci
            fabsci = safe_exp(absci)!,underflow=mn_scf)
            !calculate the integral
            !          fabsci= safe_gaussian(x0,x*hgrid,pgauss)
            !print *,'test',i,scf_data(i),fabsci,pgauss,absci,log(tiny(1.d0)),tiny(1.d0)
            gint = gint + scf_data(i)*fabsci

         end do
      end if
      gint = gint*dx
      if (abs(gint) < mp_tiny) gint=0.0_gp
    end function scfdotf

end module multipole_preserving


!> Creates the charge density of a Gaussian function, to be used for the local part
!! of the pseudopotentials (gives the error function term when later processed by the Poisson solver).
subroutine mp_gaussian_density(rloc, zion, multipole_preservingl, use_iterator,boxit,&
     mp_isf,nmpx, nmpy, nmpz, mpx, mpy, mpz, nrho, density)
  !use gaussians, only: mp_exp
  use box
  use f_precisions, only: dp=>f_double
  use multipole_preserving
  use dynamic_memory
  use numerics, only: pi
  use dictionaries, only: f_err_throw
  implicit none
  ! Calling arguments
  logical,intent(in) :: multipole_preservingl, use_iterator
  integer,intent(in) :: nrho
  real(dp),intent(in) :: rloc
  integer,intent(in) :: zion !< ionic charge (integer!)
  integer,intent(in) :: mp_isf !< interpolating scaling function order for the multipole preserving
  integer,intent(in) :: nmpx, nmpy, nmpz !< sizes of the temporary arrays; if too small the code stops
  type(box_iterator), intent(inout) :: boxit
  real(kind=8),dimension(0:nmpx),intent(inout) :: mpx !< temporary array for the exponetials in x direction
  real(kind=8),dimension(0:nmpy),intent(inout) :: mpy !< temporary array for the exponetials in y direction
  real(kind=8),dimension(0:nmpz),intent(inout) :: mpz !< temporary array for the exponetials in z direction
  real(dp),dimension(nrho),intent(inout) :: density
  ! Local variables
  real(dp),parameter :: mp_tiny = 1.e-30_dp
  logical :: perx, pery, perz,gox, goy, goz
  logical, dimension(3) :: peri
  integer :: i3s, n3pi,n1i,n2i,n3i
  real(dp) :: rlocinv2sq, charge, cutoff, xp, yp, zp, rx, ry, rz, hxh, hyh, hzh,fx,fy,fz
  integer :: i1, i2, i3, isx, iex, isy, iey, isz, iez, j1, j2, j3, ind

  call f_routine(id='mp_gaussian_density')

  rx=boxit%oxyz(1)
  ry=boxit%oxyz(2)
  rz=boxit%oxyz(3)

  i3s=boxit%i3s
  n3pi=boxit%i3e-i3s+1

  hxh=boxit%mesh%hgrids(1)
  hyh=boxit%mesh%hgrids(2)
  hzh=boxit%mesh%hgrids(3)

  peri=cell_periodic_dims(boxit%mesh)
  perx=peri(1)
  pery=peri(2)
  perz=peri(3)

  n1i=boxit%mesh%ndims(1)
  n2i=boxit%mesh%ndims(2)
  n3i=boxit%mesh%ndims(3)

  rlocinv2sq=0.5_dp/rloc**2
  charge=real(zion,dp)/(2.0_dp*pi*sqrt(2.0_dp*pi)*rloc**3)

  !cutoff of the range
  cutoff=10.0_dp*rloc
  if (multipole_preservingl) then
     !We want to have a good accuracy of the last point rloc*10
     !cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=dp)
     cutoff=cutoff+max(hxh,hyh,hzh)*real(mp_isf,kind=dp)
  end if

  isx=boxit%nbox(1,1)
  isy=boxit%nbox(1,2)
  isz=boxit%nbox(1,3)
  iex=boxit%nbox(2,1)
  iey=boxit%nbox(2,2)
  iez=boxit%nbox(2,3)

  ! Check whether the temporary arrays are large enough
  if (iex-isx>nmpx) then
     call f_err_throw('Temporary array in x direction too small')
  end if
  if (iey-isy>nmpy) then
     call f_err_throw('Temporary array in y direction too small')
  end if
  if (iez-isz>nmpz) then
     call f_err_throw('Temporary array in z direction too small')
  end if

  !$omp parallel
  !$omp do
  do i1=isx,iex
     mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,multipole_preservingl)
  end do
  !$omp end do nowait
  !$omp do
  do i2=isy,iey
     mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,multipole_preservingl)
  end do
  !$omp end do nowait
  !$omp do
  do i3=isz,iez
     mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,multipole_preservingl)
  end do
  !$omp end do
  !$omp end parallel


  if (use_iterator) then
     do while(box_next_z(boxit))
        !fz=mpz(boxit%inext(3)-boxit%nbox(1,3))
        fz=mpz(boxit%k-boxit%nbox(1,3)+1)
        do while(box_next_y(boxit))
           fy=mpy(boxit%j-boxit%nbox(1,2)+1)
           do while(box_next_x(boxit))
              fx=mpx(boxit%i-boxit%nbox(1,1)+1)
              xp=fx*fy*fz
              density(boxit%ind) = density(boxit%ind) - xp*charge
           end do
        end do
     end do
  else
     !$omp parallel do default(shared) &
     !$omp private(i3,i2,i1,zp,yp,xp,ind,j1,j2,j3,gox,goz,goy)
     do i3=isz,iez
        zp = mpz(i3-isz)
        if (abs(zp) < mp_tiny) cycle
        call ind_positions_new(perz,i3,n3i,j3,goz)
        do i2=isy,iey
           yp = zp*mpy(i2-isy)
           if (abs(yp) < mp_tiny) cycle
           call ind_positions_new(pery,i2,n2i,j2,goy)
           do i1=isx,iex
              xp = yp*mpx(i1-isx)
              if (abs(xp) < mp_tiny) cycle
              call ind_positions_new(perx,i1,n1i,j1,gox)
              if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                 ind=j1+(j2-1)*n1i+(j3-i3s)*n1i*n2i
                 density(ind)=density(ind)-xp*charge
              endif
           enddo
        enddo
     enddo
     !$omp end parallel do
  end if
  call f_release_routine()

contains

  !> Determine the index in which the potential must be inserted, following the BC
  !! Determine also whether the index is inside or outside the box for free BC
  pure subroutine ind_positions_new(periodic,i,ni,j,go)
    implicit none
    logical, intent(in) :: periodic
    integer, intent(in) :: i,ni
    logical, intent(out) :: go
    integer, intent(out) :: j

    if (periodic) then
       go=.true.
       j=modulo(i,ni)
    else
       j=i
       if (i >= -14 .and. i <= ni-15) then
          go=.true.
       else
          go=.false.
       end if
    end if

  END SUBROUTINE ind_positions_new

end subroutine mp_gaussian_density
