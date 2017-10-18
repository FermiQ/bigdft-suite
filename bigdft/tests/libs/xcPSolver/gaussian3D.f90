!> @file
!!  Use integral form for Poisson solver
!! @author
!!    Copyright (c) 2013-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Program testing new ideas like momentum-preserving gaussian integrals 3D
program gaussian3D
  use module_base, only: gp
  implicit none
  real(gp), parameter :: hgrid = 0.4_gp
  real(gp) :: rloc,shift
  integer :: i

  call f_lib_initialize()

  !call MP_gaussian()
  !call test_gain()

  rloc = 1.0_gp
  shift = 0.0_gp
  do
    call MP_rloc(hgrid,shift,rloc,nmoms=3)
    rloc = rloc * 0.1_gp
    shift = shift + 0.05_gp
    if (rloc < 1.0e-10_gp) exit
  end do

shift = 0.0_gp
  do i=0,3
    !shift = real(i,gp)*0.1_gp
    call MP_rloc_zero(hgrid,shift,5 )
  end do

  call f_lib_finalize()

contains

  !> Check multipole preserving for small rloc up to rloc==zero
  subroutine MP_rloc(hgrid,shift,rloc,nmoms)
    use module_base
    use gaussians, only: gauint0
    use multipole_preserving
    use yaml_output, only: yaml_map,yaml_mapping_open,yaml_mapping_close,yaml_comment
    implicit none
    !Arguments
    real(gp), intent(in) :: hgrid !< step size of the 3D mesh
    real(gp), intent(in) :: shift !< shift of the gaussian in the mesh
    real(gp), intent(in) :: rloc !< rloc of erf function
    integer, intent(in) :: nmoms !< Number of calculated moments
    !Local variables
    logical, parameter :: multipole_preserving = .true.
    integer, parameter :: mp_isf=16
    !integer, parameter :: nmoms = 7
    real(gp), parameter :: gammaonehalf=1.772453850905516027298_gp ! i.e. sqrt(pi)
    real(gp), parameter :: e_diff = 1.0e-13_gp
    real(gp), dimension(:,:,:), allocatable :: moments
    real(gp), dimension(:), allocatable :: mp,rx
    real(gp), dimension(:,:,:), allocatable :: array
    real(gp) :: x,y,z,reference,diff,factor,rlocinv2sq,cutoff,vol
    integer :: is,ie,n,i,ix,iy,iz,mx,my,mz

    call f_routine(id='MP_rloc')

    rlocinv2sq=0.5_gp/rloc**2
    cutoff=10.0_gp*rloc+hgrid*real(mp_isf,kind=gp)

    !Normalization factor
    !factor = 1.0_gp/gauint0(rlocinv2sq,0)
    !factor = sqrt(rlocinv2sq)/gammaonehalf
    factor = sqrt(0.5_gp)/(rloc*gammaonehalf)
    call yaml_mapping_open('MP_rloc',flow=.true.)
    call yaml_map('Normalized factor',factor*factor*factor)
    call yaml_map('rloc',rloc)
    call yaml_map('hgrid',hgrid)
    call yaml_map('shift', shift)
    call yaml_map('range',cutoff)

    moments = f_malloc( (/ 0.to.nmoms, 0.to.nmoms, 0.to.nmoms /),id="moments")

    !Separable function: do 1-D integrals before and store it.
    is=floor((-cutoff)/hgrid)
    ie=ceiling((cutoff)/hgrid)
    n = (ie - is + 1)

    call yaml_map('is',is)
    call yaml_map('ie',ie)

    !Separable function: do 1-D integrals before and store it.
    rx = f_malloc( is .to. ie,id='mpx')
    mp = f_malloc( is .to. ie,id='mpx')

    call yaml_map('mp size',size(mp))

    !Initialize the work arrays needed to integrate with is
    call initialize_real_space_conversion(isf_m=mp_isf,rlocs=(/rloc/),verbose=.true.)
    call yaml_mapping_close()

    do i=is,ie
       rx(i) =  real(i,gp)*hgrid
       mp(i) = factor*mp_exp(hgrid,shift,rlocinv2sq,i,0,multipole_preserving)
       print *,i,mp(i)
    end do

    !Calculate the corresponding 3D array
    array = f_malloc( (/ is.to.ie, is.to.ie, is.to.ie /),id='array')

    do iz=is,ie
      do iy=is,ie
        do ix=is,ie
          array(ix,iy,iz) = mp(ix)*mp(iy)*mp(iz)
        end do
      end do
    end do

    !Calculate the moments
    vol=hgrid*hgrid*hgrid
    do mz=0,nmoms
      do my=0,nmoms
        do mx=0,nmoms
          moments(mx,my,mz)=0.0_gp
          do iz=is,ie
            z=rx(iz)
            do iy=is,ie
              y=rx(iy)
              do ix=is,ie
                x=rx(ix)
                !Center to the shift to avoid complicated reference
                moments(mx,my,mz)=moments(mx,my,mz)&
                   +( (x-shift)**mx*(y-shift)**my*(z-shift)**mz)*array(ix,iy,iz)
              end do
            end do
          end do
          moments(mx,my,mz)=moments(mx,my,mz)*vol
        end do
      end do
    end do

    !print moments value
    call yaml_comment("Null integral with diff <"//trim(yaml_toa(e_diff))//" is not displayed.")
    call yaml_comment("Diff, moment, reference")
    call yaml_mapping_open('Moments')
    do mx=0,nmoms
      do my=0,nmoms
        do mz=0,nmoms
          reference=factor**3*gauint0(rlocinv2sq,mx)*gauint0(rlocinv2sq,my)*gauint0(rlocinv2sq,mz)
          diff = abs(reference-moments(mx,my,mz))
          if (diff < e_diff) then
          else
            call yaml_map(&
                trim(yaml_toa(mx))+"-"+trim(yaml_toa(my))+"+"+trim(yaml_toa(mz)),&
                (/ diff, moments(mx,my,mz), reference /),fmt="(1pe10.3)" )
          end if
        end do
      end do
    end do
    call yaml_mapping_close()

    !De-allocate
    call f_free(rx,mp)
    call f_free(array)
    call f_free(moments)

    call finalize_real_space_conversion()

    call f_release_routine()

  end subroutine MP_rloc


  !> Check multipole preserving for rloc == zero using a ISF directly
  subroutine MP_rloc_zero(hgrid,shift,nmoms)
    use module_base
    use gaussians, only: gauint0
    use multipole_preserving
    use yaml_output, only: yaml_map,yaml_mapping_open,yaml_mapping_close,yaml_comment
    implicit none
    !Arguments
    real(gp), intent(in) :: hgrid !< step size of the 3D mesh
    real(gp), intent(in) :: shift !< shift of the grid to delta or delat to grid
    integer, intent(in) :: nmoms !< Number of calculated moments
    !Local variables
    logical, parameter :: multipole_preserving = .true.
    integer, parameter :: mp_isf=16
    !integer, parameter :: nmoms = 7
    real(gp), parameter :: gammaonehalf=1.772453850905516027298_gp ! i.e. sqrt(pi)
    real(gp), parameter :: e_diff = 1.0e-12_gp
    real(gp), dimension(:,:,:), allocatable :: moments
    real(gp), dimension(:), allocatable :: mp,rx,x_scf,scf_data
    real(gp), dimension(:,:,:), allocatable :: array
    real(gp) :: x,y,z,reference,diff,cutoff
    real(gp) :: a1,a2,aa,dx,v,r0,rr,vol,x1,x2
    integer :: i,is,ie,n,ix,iy,iz,mx,my,mz
    integer :: np,n_scf,n_range,i1,ir

    call f_routine(id='MP_rloc')

    cutoff=hgrid*real(mp_isf,kind=gp)
    vol = hgrid*hgrid*hgrid

    call yaml_mapping_open("MP_rloc_zero",flow=.true.)
    call yaml_map('rloc',0.0)
    call yaml_map('hgrid',hgrid)
    call yaml_map('shift',shift)
    call yaml_map('range',cutoff)

    moments = f_malloc( (/ 0.to.nmoms, 0.to.nmoms, 0.to.nmoms /),id="moments")

    !Separable function: do 1-D integrals before and store it.
    is=floor((-cutoff)/hgrid)
    ie=ceiling((cutoff)/hgrid)
    n = (ie - is + 1)

    call yaml_map('is',is)
    call yaml_map('ie',ie)

    !Separable function: do 1-D integrals before and store it.
    rx = f_malloc( is .to. ie,id='rx')
    mp = f_malloc( is .to. ie,id='mpx')

    do i=is,ie
       rx(i) =  real(i,gp)*hgrid
       !mp(i) = factor*mp_exp(hgrid,0.0_gp,rlocinv2sq,i,0,multipole_preserving)
    end do

    call yaml_map('mp size',size(mp))

    np = 11
    !2**6 is a min to have 2**12=2048 points
    !otherwise trouble in scaling_function.f90
    n_scf=2*mp_isf*(2**np)

    !allocations for scaling function data array
    x_scf = f_malloc(0.to.n_scf,id='x_scf')
    scf_data = f_malloc(0.to.n_scf,id='scf_data')
    !Build the scaling function external routine coming from Poisson Solver. To be customized accordingly
    !call scaling_function(itype_scf,n_scf,n_range,x_scf,scf_data)
    !call wavelet_function(itype_scf,n_scf,x_scf,scf_data)
    call ISF_family(mp_isf,0,n_scf,n_range,x_scf,scf_data)

    call yaml_map('n_scf',n_scf)
    call yaml_map('n_range',n_range)
    call yaml_mapping_close()

    !Build mp
    ir = n_scf/n_range  !Step for index
    aa = real(ir,gp) ! Inverse of the step
    dx = 1.0_gp/aa ! grid step
    r0 = real(-n_range/2 + 1,gp) !Starting point

    !The first index is 'is' out.
    do i=is,ie
      !We have a ISF centered on rx(i) with a step of hgrid and a delta on 0
      !so need value for x_scf(i)*hgrid == -rx(i) (symmetric so rx(i))
      rr = (rx(i)-shift)/hgrid
      i1 = int((rr-r0)*ir)
      if (i1 < 0) then
        print *,i,i1,rr
        mp(i) = 0.0_gp
        cycle
      end if
      x1 = real(i1,gp)*dx+r0
      x2 = x1+dx
      !x1 = x_scf(i1)
      !x2 = x_scf(i1+1)
      a1 = aa*(rr-x1)
      a2 = aa*(x2-rr)
      v = a2*scf_data(i1)+a1*scf_data(i1+1)
      print *,i,i1,rr,x1,x2,v
      mp(i) = v
    end do

    !Calculate the corresponding 3D array
    array = f_malloc( (/ is.to.ie, is.to.ie, is.to.ie /),id='array')

    do iz=is,ie
      do iy=is,ie
        do ix=is,ie
          array(ix,iy,iz) = mp(ix)*mp(iy)*mp(iz)
        end do
      end do
    end do

    !Calculate the moments
    do mz=0,nmoms
      do my=0,nmoms
        do mx=0,nmoms
          moments(mx,my,mz)=0.0_gp
          do iz=is,ie
            z=rx(iz)
            do iy=is,ie
              y=rx(iy)
              do ix=is,ie
                x=rx(ix)
                moments(mx,my,mz)=moments(mx,my,mz)&
                 +(x**mx*y**my*z**mz)*array(ix,iy,iz)
              end do
            end do
          end do
          !moments(mx,my,mz)=moments(mx,my,mz)*vol
        end do
      end do
    end do

    !print moments value
    call yaml_comment("Null integral with diff <"//trim(yaml_toa(e_diff))//" is not displayed.")
    call yaml_comment("Diff, moment, reference")
    call yaml_mapping_open('Moments')
    do mx=0,nmoms
      do my=0,nmoms
        do mz=0,nmoms
          !if (mx==0.and.my==0.and.mz==0) then
          !  reference = 1.0_gp
          !else
          !  reference = 0.0_gp
          !end if
          reference = shift**(mx+my+mz)
          diff = abs(reference-moments(mx,my,mz))
          if (diff < e_diff) then
          else
            call yaml_map(&
                trim(yaml_toa(mx))+"-"+trim(yaml_toa(my))+"-"+trim(yaml_toa(mz)),&
                (/ diff, moments(mx,my,mz), reference /),fmt="(1pe10.3)" )
          end if
        end do
      end do
    end do
    call yaml_mapping_close()

    !De-allocate
    call f_free(rx,mp)
    call f_free(x_scf,scf_data)
    call f_free(array)
    call f_free(moments)

    call f_release_routine()

  end subroutine MP_rloc_zero


!> Program testing new ideas like momentum-preserving gaussian integrals 3D
subroutine test_gain()
  use module_base
  use gaussians
  use multipole_preserving
  use yaml_output
  use yaml_parse
  use pseudopotentials
  implicit none
  integer, parameter :: nmoms=1         !< Number of calculated moments
  integer, parameter :: nstep=3         !< Number of resolution to calculate the moments
  integer, parameter :: nsigma=1        !< Number of different gaussian functions
  integer, parameter :: npts=50         !< Arrays from -npts to npts
  real(gp), parameter :: hgrid = 0.8_gp !< Grid step
  real(gp), parameter :: sigma = 0.2_gp !< Sigma gaussian
  integer :: nat,ntyp,iat,i
  integer(f_long) :: t0,t1
  real(gp) :: diff
  type(gaussian_basis_new) :: G
  type(dictionary), pointer :: dict,types
  integer, dimension(:), allocatable :: iatype
  real(gp), dimension(:,:), allocatable :: rxyz,Sab,S2ab,Tab,T2ab
  real(gp), dimension(:,:,:), allocatable :: psppar

  types=>yaml_load('[ Zn, Er]')
  !extract a set of gaussian basis of two PSP
  nat=dict_len(types)
  ntyp=dict_len(types)
  iatype=f_malloc(nat,id='iatype') !all atoms are different
  call f_memcpy(src=[(iat,iat=1,nat)],dest=iatype)
  psppar=f_malloc0([0.to.4,0.to.6,1.to.ntyp],id='psppar')
  !put two atoms far apart
  rxyz=f_malloc0([3,nat],id='rxyz')
  rxyz(3,2)=3.d0

  !retrieve the parameters of the PSP by default
  call dict_init(dict)
  call psp_dict_fill_all(dict,dict_value( types // 0), 11, 15.d0, 5.d0, 8.d0)!a PSP-rich atom
  call psp_dict_fill_all(dict,dict_value( types // 1), 1, 15.d0, 5.d0, 8.d0)!a PSP-richer atom

!  call update_psp_dict(dict,'C')
!  call update_psp_dict(dict,'N')

  !print the dictionary
  call yaml_map('PSP dictionary',dict)

  !then retrieve the psppar components
  call psp_set_from_dict(dict //("psppar."+dict_value(types//0)), psppar=psppar(1:,:,1))
  call psp_set_from_dict(dict //("psppar."+dict_value(types//1)), psppar=psppar(1:,:,2))

  call yaml_map('Psppar for '+dict_value(types//0),psppar(:,:,1))
  call yaml_map('Psppar for '+dict_value(types//1),psppar(:,:,2))

  call gaussian_basis_from_psp(nat,iatype,rxyz,psppar,ntyp,G)

  Sab=f_malloc([G%ncoeff,G%ncoeff],id='Sab')
  S2ab=f_malloc([G%ncoeff,G%ncoeff],id='S2ab')
  Tab=f_malloc([G%ncoeff,G%ncoeff],id='Tab')
  T2ab=f_malloc([G%ncoeff,G%ncoeff],id='T2ab')

  call yaml_mapping_open('Basis set generated, calculating overlap')
  call yaml_map('Number of basis elements',G%ncoeff)
  !calculate the overlap matrix of the basis
  t0=f_time()
  do i=1,1000
     call overlap(G,G,Sab)
  end do
  t1=f_time()
  call yaml_map('Overlap matrix',Sab,fmt='(1pg12.3)')
  call yaml_map('Elapsed time',real(t1-t0,f_double)*1.e-9)

  t0=f_time()
  do i=1,1000
     call overlap_gain(G,G,S2ab)
  end do
  t1=f_time()
  call yaml_map('Overlap matrix with GaIn library',S2ab,fmt='(1pg12.3)')
  call yaml_map('Elapsed time',real(t1-t0,f_double)*1.e-9)

  call f_diff(int(G%ncoeff**2,f_long),Sab,S2ab,diff)
  call yaml_map('Maxdiff of both objects',diff)

  call yaml_mapping_close()

  call yaml_mapping_open('Basis set generated, calculating kinetic term')
  call yaml_map('Number of basis elements',G%ncoeff)
  !calculate the overlap matrix of the basis
  t0=f_time()
  do i=1,1000
     call kinetic(G,G,Tab)
  end do
  t1=f_time()
  call yaml_map('Laplacian matrix',Tab,fmt='(1pg12.3)')
  call yaml_map('Elapsed time',real(t1-t0,f_double)*1.e-9)

  t0=f_time()
  do i=1,1000
     call kinetic_gain(G,G,T2ab)
  end do
  t1=f_time()
  call yaml_map('Laplacian matrix with GaIn library',T2ab,fmt='(1pg12.3)')
  call yaml_map('Elapsed time',real(t1-t0,f_double)*1.e-9)

  call f_diff(int(G%ncoeff**2,f_long),Tab,-0.5d0*T2ab,diff)
  call yaml_map('Maxdiff of both objects',diff)

  call yaml_mapping_close()

  call f_free(Sab)
  call f_free(S2ab)
  call f_free(Tab)
  call f_free(T2ab)
  call f_free(iatype)
  call f_free(psppar)
  call dict_free(dict,types)

  !as the basis set is now generated we can use it to play with the Gaussian operations

  call f_free(rxyz)
  call gaussian_basis_free(G)


end subroutine test_gain


subroutine MP_gaussian()
  use module_base
  use gaussians
  use multipole_preserving
  use yaml_output
  use yaml_parse
  implicit none
  integer, parameter :: iunit=16        !< File unit for the plot
  integer, parameter :: nmoms=1         !< Number of calculated moments
  integer, parameter :: nstep=3         !< Number of resolution to calculate the moments
  integer, parameter :: nsigma=1        !< Number of different gaussian functions
  integer, parameter :: npts=50         !< Arrays from -npts to npts
  real(gp), parameter :: hgrid = 0.8_gp !< Grid step
  real(gp), parameter :: sigma = 0.2_gp !< Sigma gaussian
  integer :: j
  integer :: imoms,pow,istep,isigma
  real(gp) :: pgauss,x0,y0,z0,reference,max_fj
  real(gp), dimension(0:nmoms,2) :: moments
  real(gp), dimension(3,2,0:nmoms) :: avgmaxmin
  real(gp), dimension(:,:,:), allocatable :: fj_phi,fj_coll


  pow=0

  !pgauss=0.5_gp/((0.1_gp*hgrid)**2)!8.0e-3_dp*1.25_dp**(6*(8-1))
  !array where we have to write the value of the discretization
  fj_phi=f_malloc( (/ -npts .to. npts, -npts .to. npts, -npts .to. npts/), id='fj_phi')
  fj_coll=f_malloc((/ -npts .to. npts, -npts .to. npts, -npts .to. npts/), id='fj_coll')
  call initialize_real_space_conversion() !initialize the work arrays needed to integrate with isf

  ! Calculate for different nsigma sigma
  do isigma=1,nsigma
     pgauss=0.5_gp/((sigma+0.01_gp*(isigma-1)*hgrid)**2)
     call yaml_map('sigma/h',sqrt(0.5_gp/pgauss)/hgrid)
     !plot raw function (fort.iunit)
     do j=-npts,npts
        if (pow /= 0) then
           write(iunit,*) j,exp(-pgauss*(j*hgrid)**2)*((j*hgrid)**pow)
        else
           write(iunit,*) j,exp(-pgauss*(j*hgrid)**2)
        end if
     end do

     avgmaxmin=0.0_gp
     avgmaxmin(3,:,:)=1.d100
     max_fj=0.0_gp
     do istep=1,nstep
        x0=(-0.5_gp+real(istep-1,gp)/real(nstep,gp))*hgrid
        y0=(-0.5_gp+real(istep-1,gp)/real(nstep,gp))*hgrid
        z0=(-0.5_gp+real(istep-1,gp)/real(nstep,gp))*hgrid
        call yaml_map('x0',x0,advance='no')
        call yaml_comment('Step No.'//trim(yaml_toa(istep)),tabbing=70)
        call evaluate_moments3D(nmoms,npts,hgrid,pgauss,pow,x0,y0,z0,fj_phi,fj_coll,moments)
        max_fj=max(max_fj,maxval(abs(fj_coll-fj_phi)))
!!$  !print moments value
!!$  do imoms=0,nmoms
!!$     reference=gauint0(pgauss,imoms+pow)
!!$     if (reference /=0.0_gp) then
!!$        call yaml_map('Mom No.'//trim(yaml_toa(imoms)),&
!!$             (moments(imoms,:)-reference)/reference,fmt='(1pe22.14)',advance='no')
!!$     else
!!$        call yaml_map('Mom No.'//trim(yaml_toa(imoms)),&
!!$             moments(imoms,:),fmt='(1pe22.14)',advance='no')
!!$     end if
!!$     call yaml_comment('Ref: '//trim(yaml_toa(reference,fmt='(1pe22.14)')))
!!$  end do

        !calculate the average, maximum and minimum of each moment in function of the reference
        !j=1 use the elemental property of the mp_exp function with fj_phi
        !j=2 collocation array with fj_coll
        do j=1,2
           do imoms=0,nmoms
              reference=gauint0(pgauss,imoms+pow)**3
              print *,j,imoms,reference,moments(imoms,j)
              if (reference /= 0.0_gp) then
                 !x^even
                 moments(imoms,j) = abs((moments(imoms,j)-reference))!/reference)
              else
                 !x^odd
                 moments(imoms,j) = abs(moments(imoms,j))
              end if
              avgmaxmin(1,j,imoms) = avgmaxmin(1,j,imoms)+moments(imoms,j)/real(nstep,gp)
              avgmaxmin(2,j,imoms) = max(moments(imoms,j),avgmaxmin(2,j,imoms))
              avgmaxmin(3,j,imoms) = min(moments(imoms,j),avgmaxmin(3,j,imoms))
           end do
        end do
     end do

     !Plot fort.(iunit+1)
     write(iunit+1,'(104(1pe14.5))') sqrt(0.5_gp/pgauss)/hgrid,avgmaxmin
     call yaml_map('maxdiff' // trim(yaml_toa(isigma)), (/ sqrt(0.5_gp/pgauss)/hgrid, max_fj /) )
     !print *,'maxdiff',sqrt(0.5_gp/pgauss)/hgrid,max_fj
  end do

  call yaml_map('Results',reshape(avgmaxmin,(/6,nmoms+1/)),fmt='(1pe14.5)')

  call finalize_real_space_conversion()

  call f_free(fj_phi)
  call f_free(fj_coll)


end subroutine MP_gaussian


!> Classify the quality of a multipole extraction in both cases
subroutine evaluate_moments3D(nmoms,npts,hgrid,pgauss,pow,x0,y0,z0,fj_phi,fj_coll,moments)
  use module_base, only: gp
  use multipole_preserving, only: mp_exp
  implicit none
  !Arguments
  integer, intent(in) :: npts,pow,nmoms
  real(gp), intent(in) :: hgrid,pgauss,x0,y0,z0
  real(gp), dimension(0:nmoms,2), intent(out) :: moments
  real(gp), dimension(-npts:npts,-npts:npts,-npts:npts), intent(out) :: fj_phi,fj_coll
  !local variables
  integer :: j,jy,jz

  !use the elemental property of the mp_exp function
  do jz=-npts,npts
     do jy=-npts,npts
        fj_phi(:,jy,jz)=mp_exp(hgrid,x0,pgauss,(/(j,j=-npts,npts)/),pow,.true.) &
        & *mp_exp(hgrid,y0,pgauss,jy,pow,.true.)*mp_exp(hgrid,z0,pgauss,jz,pow,.true.)
     end do
  end do
  !scfdotf((/(j,j=-npts,npts)/),hgrid,pgauss,x0,pow)
  call moments_3d(2*npts+1,2*npts+1,2*npts+1,fj_phi, &
  & x0+hgrid*(npts+1),y0+hgrid*(npts+1),z0+hgrid*(npts+1),hgrid,nmoms,moments(0,1))

  !collocation array
  do jz=-npts,npts
     do jy=-npts,npts
        fj_coll(:,jy,jz)=mp_exp(hgrid,x0,pgauss,(/(j,j=-npts,npts)/),pow,.false.) &
        & *mp_exp(hgrid,y0,pgauss,jy,pow,.false.)*mp_exp(hgrid,z0,pgauss,jz,pow,.false.)
     end do
  end do
  !if (pow /=0) then
  !   fj_coll=(/(exp(-pgauss*(j*hgrid-x0)**2)*(j*hgrid-x0)**pow,j=-npts,npts)/)
  !else
  !   fj_coll=(/(exp(-pgauss*(j*hgrid-x0)**2),j=-npts,npts)/)
  !end if
  call moments_3d(2*npts+1,2*npts+1,2*npts+1,fj_coll,&
       & x0+hgrid*(npts+1),y0+hgrid*(npts+1),z0+hgrid*(npts+1),hgrid,nmoms,moments(0,2))

end subroutine evaluate_moments3D


!> Calculate the moments of an array with respect to a reference point
subroutine moments_3d(nx,ny,nz,array,x0,y0,z0,h,nmoms,moments)
  use module_base, only: gp
  implicit none
  !Arguments
  integer, intent(in) :: nmoms,nx,ny,nz
  real(gp), intent(in) :: x0,y0,z0,h !< grid spacing
  real(gp), dimension(nx,ny,nz), intent(in) :: array
  real(gp), dimension(0:nmoms), intent(out) :: moments
  !local variables
  integer :: j,kx,ky,kz
  real(gp) :: x,y,z

  do j=0,nmoms
     moments(j)=0.0_gp
     do kz =1,nz
        z=real(kz,gp)*h-z0
        do ky =1,ny
           y=real(ky,gp)*h-y0
           do kx=1,nx
              x=real(kx,gp)*h-x0
              moments(j)=moments(j)+(x**j*y**j*z**j)*array(kx,ky,kz)
           end do
        end do
     end do
     moments(j)=moments(j)*h*h*h
  end do

end subroutine moments_3d

end program gaussian3D
