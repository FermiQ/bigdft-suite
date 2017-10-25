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
  integer :: nlr,nspinor,nstress,is,ii,ii_fin,i,iproc,unit3,ni
  real(f_double) :: hx,hy,hz,kx,ky,kz,da,tr,hgrid,scal,enea2,enea
  real(f_double), dimension(3,3) :: stress_3_3
  logical :: volstress
  real(f_double), dimension(:), allocatable :: ene_acell,dene,volele,detgd,acell_var
  real(f_double), dimension(:,:), allocatable :: stress_ana,val
  real(f_double), dimension(:,:,:), allocatable :: stress_kin
  integer, parameter :: nord = 16
  real(f_double), parameter :: acell=20.0_f_double
  real(f_double), parameter :: acelli=20.0_f_double
  logical :: wrtfiles=.true.
  logical :: wrtfunc=.false.

  ! To run the test:
  ! ./kinetic_operator                -> test only the kinetic operator without its stress
  ! ./kinetic_operator -n 50          -> directions varied one by one
  ! ./kinetic_operator -n 50 -v yes   -> all directions varied concurrently

  iproc=0
  call f_lib_initialize()
 
  call bigdft_init_errors()
  call bigdft_init_timing_categories()

  call yaml_argparse(options,&
       '- {name: hgrid, shortname: g, default: 0.8, help_string: hgrid}'//f_cr//&
       '- {name: nstress, shortname: n, default: 1, help_string: nstress}'//f_cr//&
       '- {name: volstress, shortname: v, default: none, help_string: volstress}'//f_cr//&
       '- {name: sigma, shortname: s, default: 1.4  , help_string: sigma}')

  hgrids=options//'hgrid'
  sigma=options//'sigma'
  nstress=options//'nstress'
  volstress=options//'volstress'

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

  do ii=1,ii_fin ! loop on the three x, y, z components.

   do is=1,nstress

      if (iproc==0) then
        call yaml_comment('Stress itetation',hfill='-')
        call yaml_mapping_open('Kinetic stress input')
        call yaml_map('Kinetic stress iteration', is)
        call yaml_map('Check all directions concurrently', volstress)
        if (.not. volstress) call yaml_map('direction', ii)
        if (volstress) call yaml_map('direction', 'all')
      end if
      hx=hgrids_init(1)
      hy=hgrids_init(2)
      hz=hgrids_init(3)
      k=1.0_f_double
   
      acell_cur=acell
      acell_var(is)=acelli
      if (nstress.gt.1) then
       scal=1.0_f_double+real((is-1),kind=8)/real(nstress-1,kind=8)
       acell_var(is)=acelli*scal
       if (volstress) then
        k(1:3)=1.0_f_double/scal
        acell_cur(1:3)=acell_var(is)
          hx=hgrids_init(1)*scal
          hy=hgrids_init(2)*scal
          hz=hgrids_init(3)*scal
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
      !grid for the free BC case
      hgrid=max(hx,hy,hz)
   
      !dict_posinp=f_tree_load('{positions: [{ C: [0.0, 0.0, 0.0]}], cell: [10,10,10]}')
      dict_posinp=f_tree_load('{positions: [{ C: [0.0, 0.0, 0.0]}], cell:'//trim(yaml_toa(acell_cur))//' }')
    
      crmult=1000.0_f_double
      frmult=1000.0_f_double
      angrad=onehalf*pi
      oxyz=5.0_f_double
      kpoint=0.0_f_double
      if (iproc==0) call yaml_map('sigma',sigma)
    
      coeff=[1.0_f_double]
      expo=[0.5_f_double/sigma**2]
      rxyz=acelli/2.0_f_double
      rxyz=rxyz*(1.0_f_double/k)
      if (iproc==0) call yaml_map('Gaussian center',rxyz)
      call dict_free(options)
      
      call define_lr(lr,dict_posinp,crmult,frmult,hgrids)
    
      call f_tree_free(dict_posinp)
      call allocate_workarrays_projectors(lr%d%n1, lr%d%n2, lr%d%n3, wp)
      psi=f_malloc0(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,id='psi')
      tpsi=f_malloc0(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,id='tpsi')
    
      call initialize_real_space_conversion(isf_m=16)
    
      projector_real=f_malloc(lr%mesh%ndim,id='projector_real') 
      call initialize_work_arrays_sumrho(lr,.true.,w)

      if (iproc==0) then
         call yaml_mapping_open('Some checks')
         call yaml_map('hgrids',hgrids)
         call yaml_map('k',k)
         call yaml_map('acell',acell_cur)
         call yaml_map('da',da)
         call yaml_map('Controvariant Metric',lr%mesh%gu)
         call yaml_map('mesh%detgd',lr%mesh%detgd)
         call yaml_map('lr%mesh%ndim',lr%mesh%ndim)
         call yaml_map('lr%mesh%volume_element',lr%mesh%volume_element)
         call yaml_map('lr%mesh_coarse%volume_element',lr%mesh%volume_element)
         call yaml_mapping_close()
      end if
    
      !build up the input gaussian as done in gaussian_to_wavelets_locreg
      gaussian=f_malloc(lr%mesh%ndim,id='gaussian')
      oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
      call real_space_gaussian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,acell_cur,k,rxyz,oxyz,lr,gaussian,1,iproc)
    
      nlr=1
      nspinor=1
      hpsi=f_malloc0(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,id='hpsi')
      psir=f_malloc(lr%mesh%ndim,id='psir') 
      hx=lr%mesh%hgrids(1)
      hy=lr%mesh%hgrids(2)
      hz=lr%mesh%hgrids(3)
      kx=0.0_f_double
      ky=0.0_f_double
      kz=0.0_f_double
    
      call initialize_work_arrays_locham(lr,nspinor,.true.,wl)
      call f_zero(hpsi)

      if (wrtfunc) then 
        if (iproc==0) call yaml_map('total n',lr%mesh%ndim)
        ni=lr%mesh%ndim**(1.0_f_double/3.0_f_double)
        if (iproc==0) call yaml_map('ni',ni)
        call print_vect(ni,ni,ni,13,gaussian)
      end if 

      ! move the input gaussian in Daubechies space
      call isf_to_daub(lr,w,gaussian,hpsi)
      call f_zero(tpsi)      
      tpsi=hpsi

      ! to be done for the proper use of isf_to_daub_kinetic
      call daub_to_isf_locham(nspinor,lr,wl,hpsi,psir)

      ! import to avoid spurius cumulations
      call f_zero(psir)
      call f_zero(hpsi)

      ! calculate the kinetic operator on the input gaussian in Daubechies    
      call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,nspinor,lr,wl,psir,hpsi,ekin,k_strten)

      ekin=ekin*0.25_f_double*lr%mesh%volume_element 
      ene_acell(is)=ekin

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
       val(is,ii)=val(is,ii)+lr%mesh%gu(i,ii)*stress_3_3(i,ii)
      end do
      detgd(is)=lr%mesh%detgd
      volele(is)=lr%mesh%volume_element

      ! calculate the kinetic energy in Daubechies (to be compared to the bigdft ekin)
      enea=f_dot(hpsi,tpsi)
      enea=0.25_f_double*enea*lr%mesh%volume_element
    
      ! move the laplacian in isf space (important -> nullify projector_real)
      call f_zero(projector_real)
      call daub_to_isf(lr,w,hpsi,projector_real)
    
      !build up the analytical kinetic operator applyed to a gaussian
      call f_zero(gaussian)
      oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
      call real_space_laplacian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,k,rxyz,oxyz,lr,gaussian)
    
      call f_diff(f_size(gaussian),gaussian,projector_real,maxdiff1)
    
      if (wrtfunc) then 
       call print_vect(ni,ni,ni,11,projector_real)
       call print_vect(ni,ni,ni,12,gaussian)
       call print_vect(ni,ni,ni,14,gaussian-projector_real)
      end if 
      call f_zero(psi)
      ! move the analytical laplacian in the Dauchies space for comparison with
      ! the bigdft output
      call isf_to_daub(lr,w,gaussian,psi)
    
      !calculate the difference of the two arrays
    
      call f_diff(f_size(hpsi),hpsi,psi,maxdiff)

      ! compute the kinetic energy in isf space    
      call f_zero(gaussian)
      oxyz=lr%mesh%hgrids*[lr%nsi1,lr%nsi2,lr%nsi3]
      call real_space_gaussian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,acell_cur,k,rxyz,oxyz,lr,gaussian,2,iproc)
      call f_zero(projector_real)
      call real_space_laplacian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,k,rxyz,oxyz,lr,projector_real)
      enea2=f_dot(gaussian,projector_real)
      enea2=0.25_f_double*enea2*lr%mesh%volume_element
      if (iproc==0) then
         call yaml_mapping_open('Parallel calculation')
         call yaml_map('Numerical kinetic energy f_dot in Daubechies',enea)
         call yaml_map('Numerical kinetic energy f_dot in isf',enea2)
         call yaml_map('Bigdft Ekinetic',ekin)
         call yaml_map('Stress tensor',k_strten)
         call yaml_map('Stress 3x3',stress_3_3)
         call yaml_map('sum_i mesh%gu(i,ii)*stress_3_3(i,ii)',val(is,ii))
         call yaml_map('Maximum difference of the kinetic operator in Daubechies',maxdiff)
         call yaml_map('Maximum difference of the kinetic operator in isf',maxdiff1)
         call yaml_mapping_close()
      end if
    
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
   
   
      if (iproc==0) call yaml_mapping_close()
   
   end do ! End of the stress loop

!post-processing of stress calculation

   call FD_first_der('F',nstress,da,ene_acell,dene,nord)
   !call fssnord1DmatNabla('F',nstress,da,ene_acell,dene,nord)

   do is=1,nstress
    if (volstress) then
     stress_ana(is,ii)=-2.d0*real(lr%mesh%ndim,kind=8)*dene(is)/(acell_var(is)*acell_var(is))
    else
     !stress_ana(is,ii)=-dene(is)/(acell*acell)/mesh%detgd
     stress_ana(is,ii)=-2.d0*real(lr%mesh%ndim,kind=8)*dene(is)/(acelli*acelli)/sqrt(lr%mesh%detgd)
    end if
    if (wrtfiles) write(unit3,'(1(1x,i8),8(1x,1pe26.14e3))')is,acell_var(is),&
                  ene_acell(is),dene(is),stress_ana(is,ii),val(is,ii),stress_kin(is,1,1),&
                  stress_kin(is,2,1),stress_kin(is,3,1)
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
     call yaml_map('Contravariant Metric',lr%mesh%gu)
     call yaml_map('Covariant Metric',lr%mesh%gd)
     call yaml_map('Product of the two',matmul(lr%mesh%gu,lr%mesh%gd))
     call yaml_map('Covariant determinant',lr%mesh%detgd)
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

    subroutine real_space_gaussian(ncplx_g,n,l,ider,nterm_max,coeff,expo,acelli,acell,k,rxyz,oxyz,lr,gaussian,pt,iproc)
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
      real(f_double), dimension(3), intent(in) :: acell
      real(f_double), dimension(3), intent(in) :: k,rxyz,oxyz
      type(locreg_descriptors), intent(in) :: lr
      real(f_double), dimension(lr%mesh%ndim), intent(out) :: gaussian
      integer, intent(in) :: pt,iproc
      ! Local variables
      integer, dimension(2*l-1) :: nterms
      integer, dimension(nterm_max,3,2*l-1) :: lxyz
      real(f_double), dimension(ncplx_g) :: sigma_and_expo
      real(f_double), dimension(ncplx_g,nterm_max,2*l-1) :: factors
      type(box_iterator) :: bit
      type(gaussian_real_space) :: g
      type(f_function), dimension(3) :: funcs
      real(f_double), dimension(3) :: noxyz,ex
      real(f_double) :: fact,sumv,length,Itot
      integer :: m,i,ii,j,kk
      real(f_double), dimension(3) :: coeffs,Inte,acell_cur

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
      length=acelli
      do m=1,2*l-1
         do i=1,3
            ex(i)=expo(1)*k(i)*k(i)
            !funcs(i)=f_function_new(f_cosine,length=acelli,frequency=2.0_f_double*k(i))
            funcs(i)=f_function_new(f_gaussian,exponent=ex(i))
            !funcs(i)=f_function_new(f_polynomial,coefficients=coeffs)
         end do
      !fact=sqrt(ex(1)*ex(2)*ex(3)/(pi**3.0_f_double))
      fact=(8.0_f_double*ex(1)*ex(2)*ex(3)/(pi**3.0_f_double))**0.25_f_double

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
       sumv=sumv+gaussian(i)*gaussian(i)
      end do
     sumv=sumv*lr%mesh%volume_element
     if (pt==1 .and. iproc==0) then
        call yaml_mapping_open('Input gaussian data')
        call yaml_map('gaussian fact',fact)
        call yaml_map('Integral of gaussian**2',sumv)
     end if
     ! Compute analytical total energy as Int g \Delta d dxdydz
     acell_cur=0.5_f_double*acell

     Itot=0.0_f_double
     do ii=1,3
        i=mod(ii-1,3)+1
        j=mod(ii,3)+1
        kk=mod(ii+1,3)+1
       
        Inte(i)=-2.0_f_double*acell_cur(i)*ex(i)*(fact**2.0_f_double)*exp(-2.0_f_double*(acell_cur(i)**2.0_f_double)*ex(i))&
             - sqrt(pi/2.0_f_double)*sqrt(ex(i))*(fact**2.0_f_double)*erf(sqrt(2.0_f_double)*acell_cur(i)*sqrt(ex(i)))
        Inte(j)=sqrt(pi/2.0_f_double)*erf(sqrt(2.0_f_double)*acell_cur(j)*sqrt(ex(j)))/sqrt(ex(j))
        Inte(kk)=sqrt(pi/2.0_f_double)*erf(sqrt(2.0_f_double)*acell_cur(kk)*sqrt(ex(kk)))/sqrt(ex(kk))
 
        Itot = Itot + Inte(i) * Inte(j) * Inte(kk)
     end do
     Itot=-0.5_f_double*Itot
     if (pt==1 .and. iproc==0) call yaml_map('Analytical kinetic energy [finite box] ',Itot)

     Itot=0.0_f_double
     do ii=1,3
        i=mod(ii-1,3)+1
        j=mod(ii,3)+1
        kk=mod(ii+1,3)+1
       
        Inte(i)= - sqrt(pi/2.0_f_double)*sqrt(ex(i))*(fact**2.0_f_double)
        Inte(j)=sqrt(pi/2.0_f_double)/sqrt(ex(j))
        Inte(kk)=sqrt(pi/2.0_f_double)/sqrt(ex(kk))

        Itot = Itot + Inte(i) * Inte(j) * Inte(kk)
     end do
     Itot=-0.5_f_double*Itot
     if (pt==1 .and. iproc==0) then
      call yaml_map('Analytical kinetic energy [all R^3,Infinity box]',Itot)
      call yaml_mapping_close()
     end if

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
      real(f_double) :: fact,length
      real(f_double), dimension(3) :: coeffs

      call get_projector_coeffs(ncplx_g,l,n,ider,nterm_max,coeff,expo,&
           nterms,lxyz,sigma_and_expo,factors)

     coeffs=[1.0_f_double,0.0_f_double,0.0_f_double] 

      length=acelli
      !for the moment only with s projectors (l=0,n=1)
      noxyz=rxyz-oxyz

      bit=box_iter(lr%mesh,origin=noxyz) !use here the real space mesh of the projector locreg
      do m=1,2*l-1
         do i=1,3
            ex(i)=expo(1)*k(i)*k(i)
            !funcs(i)=f_function_new(f_cosine,length=acelli,frequency=2.0_f_double*k(i))
            funcs(i)=f_function_new(f_gaussian,exponent=ex(i))
            !funcs(i)=f_function_new(f_polynomial,coefficients=coeffs)
         end do
      !fact=sqrt(ex(1)*ex(2)*ex(3)/(pi**3.0_f_double))
      fact=(8.0_f_double*ex(1)*ex(2)*ex(3)/(pi**3.0_f_double))**0.25_f_double

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
