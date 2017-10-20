!> test program for the function projection in wavelets
program kinetic_operator
  use module_defs, only: UNINITIALIZED
  use futile
  use box
  use f_functions
  use locregs
  use gaussians
  use locreg_operations
  use f_trees
  use BigDFT_API, only: bigdft_init_errors,bigdft_init_timing_categories
  use numerics
  use compression, only: wnrm2
  use multipole_preserving
  use numerics, only: pi
  use f_blas
  implicit none
  real(f_double) :: crmult,frmult,maxdiff,sigma,maxdiff1
  type(locreg_descriptors) :: lr
  real(f_double), dimension(3) :: kpoint,oxyz,hgrids,hgrids_init,acell_cur
  type(f_tree) :: dict_posinp
  type(workarrays_projectors) :: wp
  type(workarr_sumrho) :: w
  type(workarr_locham) :: wl
  real(f_double), dimension(:), allocatable :: psi,tpsi
  real(f_double), dimension(:), allocatable :: hpsi,psir
  type(dictionary), pointer :: options
  real(f_double), dimension(:), allocatable :: projector_real,gaussian
  real(f_double), dimension(3) :: rxyz,angdeg,angrad,k
  integer, parameter :: n=1 !<principal quantum number
  integer, parameter :: l=1 !<angular momentum of the shell
  integer, parameter :: ider=0 !<direction in which to perform the derivative (0 if any)
  integer, parameter :: nterm_max=20
  integer, parameter :: ncplx_g=1
  real(f_double), dimension(ncplx_g) :: expo 
  real(f_double), dimension(ncplx_g) :: coeff !<prefactor of the gaussian
  real(f_double) :: ekin,alpha,beta,gamma
  real(f_double), dimension(6) :: k_strten
  integer :: nlr,nspinor,nstress,is,ii,ii_fin,i,iproc,unit3,nn,ni
  real(f_double) :: hx,hy,hz,kx,ky,kz,da,tr,hgrid,scal,enea,enea2
  real(f_double), dimension(3,3) :: stress_3_3
  logical :: volstress=.false.
  real(f_double), dimension(:), allocatable :: ene_acell,dene,volele,detgd,acell_var
  real(f_double), dimension(:,:), allocatable :: stress_ana,val
  real(f_double), dimension(:,:,:), allocatable :: stress_kin
  integer, parameter :: nord = 16
  type(cell) :: mesh
  real(f_double), parameter :: acell=20.0_f_double
  real(f_double), parameter :: acelli=20.0_f_double
  logical :: wrtfiles=.true.
  logical :: wrtfunc=.false.

  iproc=0
  call f_lib_initialize()
 
  call bigdft_init_errors()
  call bigdft_init_timing_categories()

  call yaml_argparse(options,&
       '- {name: hgrid, shortname: g, default: 0.4, help_string: hgrid}'//f_cr//&
       '- {name: nstress, shortname: n, default: 1, help_string: nstress}'//f_cr//&
       '- {name: sigma, shortname: s, default: 1.4  , help_string: sigma}')

  !hgrids=0.5_f_double
  hgrids=options//'hgrid'
  sigma=options//'sigma'
  nstress=options//'nstress'

!  hgrids=[0.2_f_double,0.2_f_double,0.2_f_double]

  hgrids_init=hgrids
  angdeg(1)=90.0_f_double
  angdeg(2)=90.0_f_double
  angdeg(3)=90.0_f_double
  alpha = angdeg(1)/180.0_f_double*pi!2.0_dp*datan(1.0_dp) !to be modified
  beta  = angdeg(2)/180.0_f_double*pi!2.0_dp*datan(1.0_dp)
  gamma = angdeg(3)/180.0_f_double*pi!2.0_dp*datan(1.0_dp)
  angrad(1) = angdeg(1)/180.0_f_double*pi
  angrad(2) = angdeg(2)/180.0_f_double*pi
  angrad(3) = angdeg(3)/180.0_f_double*pi

  ene_acell = f_malloc((/ nstress /),id='ene_acell')
  acell_var = f_malloc((/ nstress /),id='acell_var')
  dene = f_malloc((/ nstress /),id='dene')
  val = f_malloc((/ nstress, 3 /),id='val')
  stress_ana = f_malloc((/ nstress, 3 /),id='stress_ana')
  stress_kin = f_malloc((/ nstress, 6, 3 /),id='stress_ps')
  volele = f_malloc((/ nstress /),id='volele')
  detgd = f_malloc((/ nstress /),id='detgd')

  unit3=205
  if (wrtfiles) call f_open_file(unit=unit3,file='func_ene_acell.dat')

  da=1.d0
  if (nstress.gt.1)  da=acelli/real(nstress-1,kind=8)

! Start of the stress code
  if (volstress) then
   ii_fin=1
  else
   ii_fin=3
  end if

  !do ii=1,ii_fin ! loop on the three x, y, z components.
  do ii=1,1

   do is=1,nstress

      if (iproc==0) then
        call yaml_comment('Stress itetation',hfill='-')
        call yaml_mapping_open('Kinetic stress input')
        call yaml_map('Kinetic stress iteration', is)
        call yaml_map('direction', ii)
      end if
      hx=hgrids_init(1)
      hy=hgrids_init(2)
      hz=hgrids_init(3)
      k=1.0_f_double
   
   !    da=acell/real(nstress-1,kind=8)
   !    acell_var=acell+real((is-1),kind=8)*da
   !    acell_var=acell*(1.0_f_double+real((is-1),kind=8)/real(nstress-1,kind=8))
   !    acell_var=acell*(1.0_f_double+real((is-1),kind=8)/real(nstress-1,kind=8))
   
      acell_cur=acell
      acell_var(is)=acelli
      if (nstress.gt.1) then
       scal=1.0_f_double+real((is-1),kind=8)/real(nstress-1,kind=8)
       acell_var(is)=acelli*scal
       if (volstress) then
   !     hx=acell_var/real(n01,kind=8)
   !     hy=acell_var/real(n02,kind=8)
   !     hz=acell_var/real(n03,kind=8)
       else
        k(ii)=1.0_f_double/scal
        acell_cur(ii)=acell_var(is)
        if (ii.eq.1) then
          !hx=hgrids_init(1) + real((is-1),kind=8)*da
          hx=hgrids_init(1)*scal
          hy=hgrids_init(2)
          hz=hgrids_init(3)
        else if (ii.eq.2) then
          hx=hgrids_init(1)
          hy=hgrids_init(2)*scal
          hz=hgrids_init(3)
        else if (ii.eq.3) then
          hx=hgrids_init(1)
          hy=hgrids_init(2)
          hz=hgrids_init(3)*scal
        end if
       end if
      end if
   
      hgrids=(/hx,hy,hz/)
      call yaml_map('hgrids',hgrids)
      call yaml_map('k',k)
      call yaml_map('acell_var',acell_var(is))
      call yaml_map('acell',acell_cur)
      call yaml_map('da',da)
      !grid for the free BC case
      hgrid=max(hx,hy,hz)
      !mesh=cell_new(geocode,ndims,hgrids,angrad) 
   
      !dict_posinp=f_tree_load('{positions: [{ C: [0.0, 0.0, 0.0]}], cell: [10,10,10]}')
      dict_posinp=f_tree_load('{positions: [{ C: [0.0, 0.0, 0.0]}], cell:'//trim(yaml_toa(acell_cur))//' }')
    ! dict_posinp=f_tree_load('{positions: [{ C: [0.0, 0.0, 0.0]}], cell:'+yaml_toa(acell)+'}')
    
    !  sigma=acelli*0.07d0
      crmult=1000.0_f_double
      frmult=1000.0_f_double
      angrad=onehalf*pi
      oxyz=5.0_f_double
      kpoint=0.0_f_double
       call yaml_map('sigma',sigma)
    
      coeff=[1.0_f_double]
      expo=[0.5_f_double/sigma**2]
      rxyz=acelli/2.0_f_double
      call yaml_map('rxyz',rxyz)
      rxyz=rxyz*(1.0_f_double/k)
      call yaml_map('rxyz',rxyz)
      call dict_free(options)
    !---------------------------------------------------------------------------------
    ! Check the projectors in Daubechies 
      
      call define_lr(lr,dict_posinp,crmult,frmult,hgrids)
    
      call f_tree_free(dict_posinp)
      call allocate_workarrays_projectors(lr%d%n1, lr%d%n2, lr%d%n3, wp)
      psi=f_malloc0(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,id='psi')
      tpsi=f_malloc0(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,id='tpsi')
    
    !!$  call project(psi,PROJECTION_1D_SEPARABLE)
    !!$
    !!$  call yaml_mapping_open('Traditional separable Projection')
    !!$  call yaml_map('Norm of the calculated projector',wnrm2(1,lr%wfd,psi))
    !!$  call yaml_mapping_close()
    !!$
      call initialize_real_space_conversion(isf_m=16)
    !!$
    !!$  call project(tpsi,PROJECTION_RS_COLLOCATION)
    !!$
    !!$  call yaml_mapping_open('Collocation-based separable Projection')
    !!$  call yaml_map('Norm of the calculated projector',wnrm2(1,lr%wfd,tpsi))
    !!$  call yaml_mapping_close()
    !!$
    !!$  call f_zero(tpsi) !reset
    !!$  call project(tpsi,PROJECTION_MP_COLLOCATION)
    !!$
    !!$  call yaml_mapping_open('Multipole-preserving-based separable Projection')
    !!$  call yaml_map('Norm of the calculated projector',wnrm2(1,lr%wfd,tpsi))
    !!$  call yaml_mapping_close()
    !!$
    !!$  !calculate the difference of the two arrays
    !!$  call f_diff(f_size(psi),psi,tpsi,maxdiff)
    !!$  call yaml_map('Maximum difference of the two arrays',maxdiff)
    
    !---------------------------------------------------------------------------------
    ! Check input/output gaussian in ISF real space
    
      ! compare input analytical gaussian and the roundtrip one to the daubechies
      projector_real=f_malloc(lr%mesh%ndim,id='projector_real') 
      call initialize_work_arrays_sumrho(lr,.true.,w)
    !  call daub_to_isf(lr,w,tpsi,projector_real)
    
      !build up the input gaussian as done in gaussian_to_wavelets_locreg
      gaussian=f_malloc(lr%mesh%ndim,id='gaussian')
      oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
      call yaml_map('oxyz',oxyz)
      call real_space_gaussian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,k,rxyz,oxyz,lr,gaussian)
    
      !calculate the difference of the two arrays
    !  call f_diff(f_size(gaussian),gaussian,projector_real,maxdiff)
    
    !  call yaml_map('Maximum difference of the in/out gaussian',maxdiff)
    
    
    !---------------------------------------------------------------------------------
    ! Check of the kinetic operator Daubechies
    
    ! Input analytical phi
    ! hpsi -> phi in Daubechies in isf_to_daub_kinetic -> real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nspinor), intent(inout) :: hpsi
    ! psir -> isf potential     in isf_to_daub_kinetic -> real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(in) :: psir
    
      nlr=1
      nspinor=1
      hpsi=f_malloc0(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,id='hpsi')
      psir=f_malloc(lr%mesh%ndim,id='psir') 
      hx=lr%mesh%hgrids(1)
      hy=lr%mesh%hgrids(2)
      hz=lr%mesh%hgrids(3)
      call yaml_map('lr%mesh%hgrids',lr%mesh%hgrids)
      kx=0.0_f_double
      ky=0.0_f_double
      kz=0.0_f_double
    
      call initialize_work_arrays_locham(lr,nspinor,.true.,wl)
      call f_zero(hpsi)
    !  call yaml_mapping_open('Checking of input func in isf')
    !  call yaml_map('size func',lr%mesh%ndim)
    !  call yaml_map('func maxval',maxval(gaussian))
    !  call yaml_map('func at 1',gaussian(1))
    !  call yaml_map('func at /4',gaussian(lr%mesh%ndim/4))
    !  call yaml_map('func at /3',gaussian(lr%mesh%ndim/3))
    !  call yaml_map('func at /2',gaussian(lr%mesh%ndim/2))
    !  call yaml_map('func at end',gaussian(lr%mesh%ndim))
    !  call yaml_mapping_close()
    if (wrtfunc) then 
      call yaml_map('total n',lr%mesh%ndim)
      ni=lr%mesh%ndim**(1.0_f_double/3.0_f_double)
      call yaml_map('ni',ni)
      call print_vect(ni,ni,ni,13,gaussian)
    end if 
      call isf_to_daub(lr,w,gaussian,hpsi)
      !hpsi=tpsi
    !!  nn=lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
    !  call yaml_mapping_open('Checking of input func in Daubechies')
    !  call yaml_map('size func',nn)
    !  call yaml_map('func maxval',maxval(hpsi))
    !  call yaml_map('func at 1',hpsi(1))
    !  call yaml_map('func at /4',hpsi(nn/4))
    !  call yaml_map('func at /3',hpsi(nn/3))
    !  call yaml_map('func at /2',hpsi(nn/2))
    !  call yaml_map('func at end',hpsi(nn))
    !  call yaml_mapping_close()
      
      call daub_to_isf_locham(nspinor,lr,wl,hpsi,psir)
    
      call f_zero(psir)
      call f_zero(hpsi)
    
      call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,nspinor,lr,wl,psir,hpsi,ekin,k_strten)
    
    !  call yaml_mapping_open('Checking of output of isf_to_daub_kinetic in Daubechies')
    !  call yaml_map('size func',nn)
    !  call yaml_map('func maxval',maxval(hpsi))
    !  call yaml_map('func at 1',hpsi(1))
    !  call yaml_map('func at /4',hpsi(nn/4))
    !  call yaml_map('func at /3',hpsi(nn/3))
    !  call yaml_map('func at /2',hpsi(nn/2))
    !  call yaml_map('func at end',hpsi(nn))
    !  call yaml_mapping_close()
    
      call f_zero(projector_real)
      call daub_to_isf(lr,w,hpsi,projector_real)
    !  call yaml_mapping_open('Checking of isf_to_daub_kinetic in isf')
    !  call yaml_map('size func',lr%mesh%ndim)
    !  call yaml_map('func maxval',maxval(projector_real))
    !  call yaml_map('func at 1',projector_real(1))
    !  call yaml_map('func at /4',projector_real(lr%mesh%ndim/4))
    !  call yaml_map('func at /3',projector_real(lr%mesh%ndim/3))
    !  call yaml_map('func at /2',projector_real(lr%mesh%ndim/2))
    !  call yaml_map('func at end',projector_real(lr%mesh%ndim))
    !  call yaml_mapping_close()
    
    !  call yaml_map('max diff of pr',projector_real(2084897))
    
      !build up the analytical kinetic operator applyed to a gaussian
      call f_zero(gaussian)
      oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
      call real_space_laplacian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,k,rxyz,oxyz,lr,gaussian)
    !  call yaml_mapping_open('Checking of analytical laplacian in isf')
    !  call yaml_map('size func',lr%mesh%ndim)
    !  call yaml_map('laplacian maxval',maxval(gaussian))
    !  call yaml_map('laplacian func at 1',gaussian(1))
    !  call yaml_map('laplacian func at /4',gaussian(lr%mesh%ndim/4))
    !  call yaml_map('laplacian func at /3',gaussian(lr%mesh%ndim/3))
    !  call yaml_map('laplacian func at /2',gaussian(lr%mesh%ndim/2))
    !  call yaml_map('laplacian func at end',gaussian(lr%mesh%ndim))
    !  call yaml_mapping_close()
    !  call yaml_map('max diff of ga',gaussian(2084897))
    
      call f_diff(f_size(gaussian),gaussian,projector_real,maxdiff1)
    
    !!  diffc=0.0d0
    !!  icur=1
    !!  do i=1,lr%mesh%ndim
    !!   diffcur=abs(gaussian(i)-projector_real(i))
    !!   if (diffcur.gt.diffc) then
    !!    diffc=diffcur
    !!    icur=i
    !!   end if
    !!  end do
    !!  call yaml_map('projector_real',projector_real(icur))
    !!  call yaml_map('gaussian',gaussian(icur))
    !!  ni=lr%mesh%ndim**(1.0_f_double/3.0_f_double)
    !!  call yaml_map('total n',lr%mesh%ndim)
    !!  call yaml_map('ni',ni)
    if (wrtfunc) then 
      call print_vect(ni,ni,ni,11,projector_real)
      call print_vect(ni,ni,ni,12,gaussian)
      call print_vect(ni,ni,ni,14,gaussian-projector_real)
    end if 
      call f_zero(psi)
      call isf_to_daub(lr,w,gaussian,psi)
    !  nn=lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
    !  call yaml_mapping_open('Checking of analytical laplacian from isf_to_daub in Daubechies')
    !  call yaml_map('size func',nn)
    !  call yaml_map('laplacian maxval',maxval(psi))
    !  call yaml_map('laplacian func at 1',psi(1))
    !  call yaml_map('laplacian func at /4',psi(nn/4))
    !  call yaml_map('laplacian func at /3',psi(nn/3))
    !  call yaml_map('laplacian func at /2',psi(nn/2))
    !  call yaml_map('laplacian func at end',psi(nn))
    !  call yaml_mapping_close()
    
      !calculate the difference of the two arrays
    
      call f_diff(f_size(hpsi),hpsi,psi,maxdiff)
      call yaml_map('Maximum difference of the kinetic operator in Daubechies',maxdiff)
      call yaml_map('Maximum difference of the kinetic operator in isf',maxdiff1)
    !  call yaml_map('Maximum difference of the kinetic operator in isf',diffc)
    !  call yaml_map('Maximum difference of the kinetic operator in isf',icur)
    
      call f_zero(gaussian)
      oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
      call real_space_gaussian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,k,rxyz,oxyz,lr,gaussian)
      call f_zero(projector_real)
      call real_space_laplacian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,k,rxyz,oxyz,lr,projector_real)
      enea=0.0_f_double
      do i=1,lr%mesh%ndim
       enea=enea+gaussian(i)*projector_real(i)
      end do
      enea=enea*lr%mesh%volume_element
      enea2=f_dot(gaussian,projector_real)
      enea2=enea2*lr%mesh%volume_element
      call yaml_map('Analytical kinetic energy',enea)
      call yaml_map('Analytical kinetic energy f_dot',enea2)


      mesh=lr%mesh
    
      call f_free(psi)
      call f_free(tpsi)
      call f_free(hpsi)
      call f_free(psir)
    
      call deallocate_work_arrays_locham(wl)
    
      call finalize_real_space_conversion()
      call deallocate_workarrays_projectors(wp)
      call deallocate_work_arrays_sumrho(w)
      call deallocate_locreg_descriptors(lr)
      call f_free(projector_real)
      call f_free(gaussian)
   
      do i=1,6
       stress_kin(is,i,ii)=k_strten(i)
      end do
   
      stress_3_3(1,1)=k_strten(1)
      stress_3_3(2,2)=k_strten(2)
      stress_3_3(3,3)=k_strten(3)
      stress_3_3(2,3)=k_strten(4)
      stress_3_3(1,3)=k_strten(5)
      stress_3_3(1,2)=k_strten(6)
      stress_3_3(3,2)=stress_3_3(2,3)
      stress_3_3(3,1)=stress_3_3(1,3)
      stress_3_3(2,1)=stress_3_3(1,2)
     
      val(is,ii)=0.d0
      do i=1,3
       val(is,ii)=val(is,ii)+mesh%gu(i,ii)*stress_3_3(i,ii)
      end do
      detgd(is)=mesh%detgd
      volele(is)=mesh%volume_element
      if (iproc==0) then
       call yaml_map('Controvariant Metric',mesh%gu)
       call yaml_map('Stress 3x3',stress_3_3)
       call yaml_map('sum_i mesh%gu(i,ii)*stress_3_3(i,ii)',val(is,ii))
       call yaml_map('mesh%detgd',mesh%detgd)
       call yaml_map('mesh%volume_element',mesh%volume_element)
      end if
   
      if (iproc == 0) then
         call yaml_mapping_open('Parallel calculation')
         call yaml_map('Ekinetic',ekin)
         call yaml_map('Stress tensor',k_strten)
         call yaml_map('Max diff',maxdiff)
         call yaml_mapping_close()
      end if
   
      ene_acell(is)=ekin
   
      if (iproc==0) call yaml_mapping_close()
   
   end do ! End of the stress loop

!post-processing of stress calculation

   call FD_first_der('F',nstress,da,ene_acell,dene,nord)
   !call fssnord1DmatNabla('F',nstress,da,ene_acell,dene,nord)

   do is=1,nstress
    if (volstress) then
     stress_ana(is,ii)=-dene(is)/(acell_var(is)*acell_var(is))
    else
     !stress_ana(is,ii)=-dene(is)/(acell*acell)/mesh%detgd
     stress_ana(is,ii)=-dene(is)/(acelli*acelli)/sqrt(mesh%detgd)
    end if
    if (wrtfiles) write(unit3,'(1(1x,i8),8(1x,1pe26.14e3))')is,acell_var(is),&
                  ene_acell(is),dene(is),stress_ana(is,ii),val(is,ii),stress_kin(is,1,1),&
                  stress_kin(is,2,1),stress_kin(is,3,1) !,detgd(is),volele(is),&
!                  stress_ps(is,1,ii),stress_ps(is,2,ii),stress_ps(is,3,ii),&
!                  stress_ps(is,4,ii),stress_ps(is,5,ii),stress_ps(is,6,ii)
   end do


  end do ! loop external ii to the stress one, for the three x,y,z directions.

  if (nstress.gt.1) then
   if (iproc==0) then
    call yaml_map('Total stress iterations', nstress)
   end if
   if (volstress) then
    tr=stress_kin(nstress/2,1,1)+stress_kin(nstress/2,2,1)+stress_kin(nstress/2,3,1)
    if (iproc == 0) then
     call yaml_comment('Stress post-processing',hfill='-')
     call yaml_mapping_open('Comparison between analytical vs psolver varing x,y,z concurrently')
     call yaml_map('Comparison at nstress/2',nstress/2)
     call yaml_map('stress analytical ',stress_ana(nstress/2,1))
     call yaml_map('stress psolver trace',tr)
     call yaml_mapping_close()
    end if
   else  
    if (iproc==0) then
     call yaml_map('Angles',[alpha,beta,gamma]*180.0_f_double*oneopi)
     call yaml_map('Contravariant Metric',mesh%gu)
     call yaml_map('Covariant Metric',mesh%gd)
     call yaml_map('Product of the two',matmul(mesh%gu,mesh%gd))
     call yaml_map('Covariant determinant',mesh%detgd)
    end if
 
    if (iproc == 0) then
     call yaml_comment('Stress post-processing',hfill='-')
     call yaml_mapping_open('Comparison between analytical vs psolver varing x,y,z individully')
     call yaml_map('Comparison at nstress/2',nstress/2)
     call yaml_map('stress analytical x',stress_ana(nstress/2,1))
     call yaml_map('stress psolver x',val(nstress/2,1))
     call yaml_map('stress analytical y',stress_ana(nstress/2,2))
     call yaml_map('stress psolver y',val(nstress/2,2))
     call yaml_map('stress analytical z',stress_ana(nstress/2,3))
     call yaml_map('stress psolver z',val(nstress/2,3))
     call yaml_mapping_close()
    end if
   end if
  end if

  if (wrtfiles) call f_close(unit3)

  call f_free(val)
  call f_free(dene)
  call f_free(ene_acell)
  call f_free(stress_ana)
  call f_free(stress_kin)
  call f_free(volele)
  call f_free(detgd)
  call f_free(acell_var)

  call f_lib_finalize()

  contains

    subroutine project(psi,method)
      use f_enums
      implicit none
      type(f_enumerator), intent(in) :: method
      real(f_double), dimension(*) :: psi
      call gaussian_to_wavelets_locreg(lr%mesh_coarse,0,&
           1,coeff,expo,UNINITIALIZED(1.0_f_double),&
           1,1,rxyz,kpoint,&
           1,lr,wp,psi,method=method)
    end subroutine project    

end program kinetic_operator

    subroutine real_space_gaussian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,k,rxyz,oxyz,lr,gaussian)
      use gaussians
      use futile
      use locregs
      use box
      use f_functions
      use numerics, only: pi
      implicit none
      integer, intent(in) :: ncplx_g !< 1 or 2 if the gaussian factor is real or complex respectively
      integer, intent(in) :: n !<principal quantum number
      integer, intent(in) :: l !<angular momentum of the shell
      integer, intent(in) :: ider !<direction in which to perform the derivative (0 if any)
      integer, intent(in) :: nterm_max !if GTH nterm_max=4 (this value should go in a module)
      real(f_double), dimension(ncplx_g), intent(in) :: coeff !<prefactor of the gaussian
      real(f_double), dimension(ncplx_g), intent(in) :: expo 
      real(f_double), intent(in) :: acelli
      real(f_double), dimension(3), intent(in) :: k,rxyz,oxyz
      type(locreg_descriptors), intent(in) :: lr
      real(f_double), dimension(lr%mesh%ndim), intent(out) :: gaussian
      ! Local variables
      integer, dimension(2*l-1) :: nterms
      integer, dimension(nterm_max,3,2*l-1) :: lxyz
      real(f_double), dimension(ncplx_g) :: sigma_and_expo
      real(f_double), dimension(ncplx_g,nterm_max,2*l-1) :: factors
      type(box_iterator) :: bit
      type(gaussian_real_space) :: g
      type(f_function), dimension(3) :: funcs
      real(f_double), dimension(3) :: noxyz,ex
      real(f_double) :: fact,sumv
      integer :: m,i
      real(f_double), dimension(3) :: coeffs

      call get_projector_coeffs(ncplx_g,l,n,ider,nterm_max,coeff,expo,&
           nterms,lxyz,sigma_and_expo,factors)

      !call gaussian_real_space_set(g,sqrt(onehalf/expo(1)),1,factors,lxyz)
      call gaussian_real_space_set(g,sigma_and_expo(1),1,factors,lxyz(1,:,1),[0],16)
      call f_zero(gaussian)
      bit=box_iter(lr%mesh)
      call three_dimensional_density(bit,g,sqrt(lr%mesh%volume_element),rxyz,gaussian)

      !for the moment only with s projectors (l=0,n=1)
      noxyz=rxyz-oxyz

      bit=box_iter(lr%mesh,origin=noxyz) !use here the real space mesh of the projector locreg

     coeffs=[1.0_f_double,0.0_f_double,0.0_f_double] 
!     call yaml_mapping_open('Data analytical gaussian')
!     call yaml_map('rxyz',rxyz)
!     call yaml_map('oxyz',oxyz)
!     call yaml_map('expo',expo)
!     call yaml_map('factor',factors)
!     call yaml_map('coeffs',coeffs)

      do m=1,2*l-1
         do i=1,3
            ex(i)=expo(1)*k(i)*k(i)
!            call yaml_map('expo',i)
!            call yaml_map('expo',ex)
            !funcs(i)=f_function_new(f_cosine,length=acelli,frequency=2.0_f_double*k(i))
            funcs(i)=f_function_new(f_gaussian,exponent=ex(i))
            !funcs(i)=f_function_new(f_polynomial,coefficients=coeffs)
         end do
!     call yaml_mapping_close()
!      fact=1.0_f_double
      fact=(sqrt(ex(1)*ex(2)*ex(3)))/(pi**(3.0_f_double/2.0_f_double))
         !here we do not consider the lxyz terms yet
         !take the reference functions
         !print *,size(gaussian),'real',lr%mesh%ndims,&
         !     lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3],&
         !     lr%mesh_coarse%hgrids*[lr%ns1,lr%ns2,lr%ns3],rxyz,noxyz
         !call separable_3d_function(bit,funcs,factors(1,1,m)*sqrt(lr%mesh%volume_element),gaussian)
         call separable_3d_function(bit,funcs,1.0_f_double*fact,gaussian)
         !call separable_3d_laplacian(bit,funcs,factors(1,1,m)*sqrt(lr%mesh%volume_element),gaussian)
      end do !not correctly written, it should be used to define the functions

      sumv=0.0_f_double
      do i=1,lr%mesh%ndim
       sumv=sumv+gaussian(i)
      end do
     sumv=sumv*lr%mesh%volume_element
     call yaml_map('lr%mesh%ndim',lr%mesh%ndim)
     call yaml_map('lr%mesh%volume_element',lr%mesh%volume_element)
     call yaml_map('gaussian fact',fact)
     call yaml_map('Integral of gaussian',sumv)

    end subroutine real_space_gaussian


    subroutine real_space_laplacian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,k,rxyz,oxyz,lr,gaussian)
      use gaussians
      use futile
      use locregs
      use box
      use f_functions
      use numerics, only: pi
      implicit none
      integer, intent(in) :: ncplx_g !< 1 or 2 if the gaussian factor is real or complex respectively
      integer, intent(in) :: n !<principal quantum number
      integer, intent(in) :: l !<angular momentum of the shell
      integer, intent(in) :: ider !<direction in which to perform the derivative (0 if any)
      integer, intent(in) :: nterm_max !if GTH nterm_max=4 (this value should go in a module)
      real(f_double), dimension(ncplx_g), intent(in) :: coeff !<prefactor of the gaussian
      real(f_double), dimension(ncplx_g), intent(in) :: expo 
      real(f_double), intent(in) :: acelli
      real(f_double), dimension(3), intent(in) :: k,rxyz,oxyz 
      type(locreg_descriptors), intent(in) :: lr
      real(f_double), dimension(lr%mesh%ndim), intent(out) :: gaussian
      ! Local variables
      integer, dimension(2*l-1) :: nterms
      integer, dimension(nterm_max,3,2*l-1) :: lxyz
      real(f_double), dimension(ncplx_g) :: sigma_and_expo
      real(f_double), dimension(ncplx_g,nterm_max,2*l-1) :: factors
      integer :: m,i
      type(box_iterator) :: bit
      type(f_function), dimension(3) :: funcs
      real(f_double), dimension(3) :: noxyz,ex
      real(f_double) :: fact
      real(f_double), dimension(3) :: coeffs

      call get_projector_coeffs(ncplx_g,l,n,ider,nterm_max,coeff,expo,&
           nterms,lxyz,sigma_and_expo,factors)

     coeffs=[1.0_f_double,0.0_f_double,0.0_f_double] 
      !call gaussian_real_space_set(g,sqrt(onehalf/expo(1)),1,factors,lxyz)
!      call gaussian_real_space_set(g,sigma_and_expo(1),1,factors,lxyz(1,:,1),[0],16)
!      call f_zero(gaussian)
!      bit=box_iter(lr%mesh)
!      call three_dimensional_density(bit,g,sqrt(lr%mesh%volume_element),rxyz,gaussian)
!     call yaml_mapping_open('Data analytical laplacian')
!     call yaml_map('rxyz',rxyz)
!     call yaml_map('oxyz',oxyz)
!     call yaml_map('expo',expo)
!     call yaml_map('factor',factors)
!     call yaml_map('coeffs',coeffs)

      !for the moment only with s projectors (l=0,n=1)
      noxyz=rxyz-oxyz

      bit=box_iter(lr%mesh,origin=noxyz) !use here the real space mesh of the projector locreg
      do m=1,2*l-1
         do i=1,3
            ex(i)=expo(1)*k(i)*k(i)
!            call yaml_map('expo',i)
!            call yaml_map('expo',ex)
            !funcs(i)=f_function_new(f_cosine,length=acelli,frequency=2.0_f_double*k(i))
            funcs(i)=f_function_new(f_gaussian,exponent=ex(i))
            !funcs(i)=f_function_new(f_polynomial,coefficients=coeffs)
         end do
      !fact=1.0_f_double
      fact=(sqrt(ex(1)*ex(2)*ex(3)))/(pi**(3.0_f_double/2.0_f_double))
!     call yaml_mapping_close()

         !here we do not consider the lxyz terms yet
         !take the reference functions
         !print *,size(gaussian),'real',lr%mesh%ndims,&
         !     lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3],&
         !     lr%mesh_coarse%hgrids*[lr%ns1,lr%ns2,lr%ns3],rxyz,noxyz
         !call separable_3d_function(bit,funcs,factors(1,1,m)*sqrt(lr%mesh%volume_element),gaussian)
         !call separable_3d_laplacian(bit,funcs,-2.0_f_double*factors(1,1,m)*sqrt(lr%mesh%volume_element),gaussian)
         call separable_3d_laplacian(bit,funcs,-2.0_f_double*fact,gaussian)
      end do !not correctly written, it should be used to define the functions
    end subroutine real_space_laplacian

!subroutine print_vect(nn,psi)
!  implicit none
!  integer, intent(in) :: nn
!  real(f_double), dimension(nn), intent(in) :: psi
!
!  call yaml_map('size func',nn)
!  call yaml_map('laplacian maxval',maxval(psi))
!  call yaml_map('laplacian func at 1',psi(1))
!  call yaml_map('laplacian func at /4',psi(nn/4))
!  call yaml_map('laplacian func at /3',psi(nn/3))
!  call yaml_map('laplacian func at /2',psi(nn/2))
!  call yaml_map('laplacian func at end',psi(nn))
!
!end subroutine get_diff

subroutine print_vect(n01,n02,n03,uni,psi)
  use futile
  implicit none
  integer, intent(in) :: n01,n02,n03,uni
  real(f_double), dimension(n01,n02,n03), intent(in) :: psi
  integer :: i1,i2,i3

  i3=n03/2
  i2=n02/2
  do i1=1,n01
  write(uni,'(1(1x,I4),1x,e14.7)')i1,psi(i1,i2,i3)
  end do

end subroutine print_vect

subroutine fssnord1DmatNabla(geocode,n01,hx,u,du,nord)
      implicit none
!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      character(len=1), intent(in) :: geocode
      integer, intent(in) :: n01,nord
      real(kind=8), intent(in) :: hx
      real(kind=8), dimension(n01) :: u
      real(kind=8), dimension(n01) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,i1,ii
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D,c1DF
      logical :: perx

      n = nord+1
      m = nord/2
      n_cell = n01

      !buffers associated to the geocode
      !conditions for periodicity
      perx=(geocode /= 'F')

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
        c1DF(i,j)=0.d0
       end do
      end do

      include 'FiniteDiffCorff.inc'
      
      do i1=1,n01
   
       du(i1) = 0.0d0
   
       if (i1.le.m) then
        if (perx) then
         do j=-m,m
          ii=modulo(i1 + j + n01 - 1, n01 ) + 1
          du(i1) = du(i1) + c1D(j,0)*u(ii)
         end do
        else
         do j=-m,m
          du(i1) = du(i1) + c1D(j,i1-m-1)*u(j+m+1)
         end do
        end if
       else if (i1.gt.n01-m) then
        if (perx) then
         do j=-m,m
          ii=modulo(i1 + j - 1, n01 ) + 1
          du(i1) = du(i1) + c1D(j,0)*u(ii)
         end do
        else
         do j=-m,m
          du(i1) = du(i1) + c1D(j,i1-n01+m)*u(n01 + j - m)
         end do
        end if
       else
        do j=-m,m
         du(i1) = du(i1) + c1D(j,0)*u(i1 + j)
        end do
       end if
       du(i1)=du(i1)/hx
   
      end do

end subroutine fssnord1DmatNabla
