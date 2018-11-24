!> @file
!!  Routines to reformat wavefunctions
!! @author
!!    Copyright (C) 2010-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Reformat one wavefunction
subroutine reformatonewave(displ,wfd,at,hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,& !n(c) iproc (arg:1)
     rxyz_old,psigold,hx,hy,hz,n1,n2,n3,rxyz,psifscf,psi)
  use module_base
  use module_types
  use box
  use at_domain
  use bounds, only: ext_buffers_coarse
  use compression
  use numerics, only: onehalf,pi
  implicit none
  integer, intent(in) :: n1_old,n2_old,n3_old,n1,n2,n3  !n(c) iproc
  real(gp), intent(in) :: hx,hy,hz,displ,hx_old,hy_old,hz_old
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz_old,rxyz
  real(wp), dimension(0:n1_old,2,0:n2_old,2,0:n3_old,2), intent(in) :: psigold
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
  real(wp), dimension(*), intent(out) :: psifscf !this supports different BC
  !local variables
  character(len=*), parameter :: subname='reformatonewave'
  logical :: cif1,cif2,cif3,perx,pery,perz
  integer :: i1,i2,i3,j1,j2,j3,l1,l2,iat,nb1,nb2,nb3,ind,jj1,jj2,jj3a,jj3b,jj3c
  integer :: ind2,ind3,n1nb1,n2nb2,n1o7,n2o7,n3o7,n1nb1o,n2nb2o,n3nb3o
  real(gp) :: hxh,hyh,hzh,hxh_old,hyh_old,hzh_old,x,y,z,dx,dy,dz,xold,yold,zold!,mindist
  real(wp) :: zr,yr,xr,ym1,y00,yp1
  real(wp), dimension(-1:1,-1:1) :: xya
  real(wp), dimension(-1:1) :: xa
  real(wp), dimension(:), allocatable :: ww,wwold
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psig
  real(wp), dimension(:,:,:), allocatable :: psifscfold
  real(gp), dimension(3) :: rd
  type(cell) :: mesh
  type(domain) :: dom
  logical, dimension(3) :: peri
  !real(kind=4) :: t0, t1
  !real(kind=8) :: time

  call f_routine(id='reformatonewave')


  !mesh=cell_new(at%astruct%geocode,[n1,n2,n3],[hx,hy,hz])
  dom=domain_new(units=ATOMIC_UNITS,bc=geocode_to_bc_enum(at%astruct%geocode),&
            alpha_bc=onehalf*pi,beta_ac=onehalf*pi,gamma_ab=onehalf*pi,acell=[n1,n2,n3]*[hx,hy,hz])
  mesh=cell_new(dom,[n1,n2,n3],[hx,hy,hz])

  !conditions for periodicity in the three directions
!!$  perx=(at%astruct%geocode /= 'F')
!!$  pery=(at%astruct%geocode == 'P')
!!$  perz=(at%astruct%geocode /= 'F')
  peri=domain_periodic_dims(mesh%dom)
  perx=peri(1)
  pery=peri(2)
  perz=peri(3)

  !buffers related to periodicity
  !WARNING: the boundary conditions are not assumed to change between new and old
  call ext_buffers_coarse(perx,nb1)
  call ext_buffers_coarse(pery,nb2)
  call ext_buffers_coarse(perz,nb3)

  psifscfold = f_malloc((/ -nb1.to.2*n1_old+1+nb1, -nb2.to.2*n2_old+1+nb2, -nb3.to.2*n3_old+1+nb3 /),id='psifscfold')
  wwold = f_malloc((2*n1_old+2+2*nb1)*(2*n2_old+2+2*nb2)*(2*n3_old+2+2*nb3),id='wwold')

  if (domain_geocode(at%astruct%dom)=='F') then
     call synthese_grow(n1_old,n2_old,n3_old,wwold,psigold,psifscfold)
  else if (domain_geocode(at%astruct%dom) =='S') then
     call synthese_slab(n1_old,n2_old,n3_old,wwold,psigold,psifscfold)
  else if (domain_geocode(at%astruct%dom) =='P') then
     call synthese_per(n1_old,n2_old,n3_old,wwold,psigold,psifscfold)
  else if (domain_geocode(at%astruct%dom) =='W') then
     call synthese_wire(n1_old,n2_old,n3_old,wwold,psigold,psifscfold)
  end if

  call f_free(wwold)

  !write(*,*) iproc,' displ ',displ
  if (hx == hx_old .and. hy == hy_old .and. hz == hz_old .and. &
       n1_old==n1 .and. n2_old==n2 .and. n3_old==n3 .and. &
       displ<= 1.d-2) then
     !if (iproc==0) write(*,*) iproc,' orbital just copied'
     call vcopy((2*n1+2+2*nb1)*(2*n2+2+2*nb2)*(2*n3+2+2*nb3),psifscfold(-nb1,-nb2,-nb3),1,&
          psifscf(1),1)
!!$     do i3=-nb3,2*n3+1+nb3
!!$        do i2=-nb2,2*n2+1+nb2
!!$           do i1=-nb1,2*n1+1+nb1
!!$              ind=i1+nb1+1+(2*n1+2+2*nb1)*(i2+nb2)+&
!!$                   (2*n1+2+2*nb1)*(2*n2+2+2*nb2)*(i3+nb3)
!!$              psifscf(ind)=psifscfold(i1,i2,i3)
!!$           enddo
!!$        enddo
!!$     enddo

  else

    rd=0.0_gp
     do iat=1,at%astruct%nat
       rd=rd+closest_r(mesh%dom,rxyz(:,iat),center=rxyz_old(:,iat))
     end do
     rd=rd/real(at%astruct%nat,gp)
     dx=rd(1)
     dy=rd(2)
     dz=rd(3)
!!$     dx=0.0_gp
!!$     dy=0.0_gp
!!$     dz=0.0_gp
!!$     !Calculate average shift
!!$     !Take into account the modulo operation which should be done for non-isolated BC
!!$     do iat=1,at%astruct%nat
!!$        dx=dx+mindist(perx,at%astruct%cell_dim(1),rxyz(1,iat),rxyz_old(1,iat))
!!$        dy=dy+mindist(pery,at%astruct%cell_dim(2),rxyz(2,iat),rxyz_old(2,iat))
!!$        dz=dz+mindist(perz,at%astruct%cell_dim(3),rxyz(3,iat),rxyz_old(3,iat))
!!$     enddo
!!$     dx=dx/real(at%astruct%nat,gp)
!!$     dy=dy/real(at%astruct%nat,gp)
!!$     dz=dz/real(at%astruct%nat,gp)

     ! transform to new structure
     !if (iproc==0) write(*,*) iproc,' orbital fully transformed'
     hxh=.5_gp*hx
     hxh_old=.5_gp*hx_old
     hyh=.5_gp*hy
     hyh_old=.5_gp*hy_old
     hzh=.5_gp*hz
     hzh_old=.5_gp*hz_old

     call f_zero((2*n1+2+2*nb1)*(2*n2+2+2*nb2)*(2*n3+2+2*nb3),psifscf(1))
     !call cpu_time(t0)

     n1nb1=(2*n1+2+2*nb1)
     n2nb2=(2*n2+2+2*nb2)
     n1nb1o=2*n1_old+1+2*nb1+1
     n2nb2o=2*n2_old+1+2*nb2+1
     n3nb3o=2*n3_old+1+2*nb3+1
     n1o7=2*n1_old+7
     n2o7=2*n2_old+7
     n3o7=2*n3_old+7

     !$omp parallel do default(none) &
     !$omp shared(nb1,nb2,nb3,n1,n2,n3,hxh,hyh,hzh,hxh_old,hyh_old,hzh_old,dx,dy,dz,perx,pery,perz) &
     !$omp shared(n1o7,n2o7,n3o7,n1nb1o,n2nb2o,n3nb3o,n1nb1,n2nb2,psifscf,psifscfold) &
     !$omp private(i1,i2,i3,j1,j2,j3,x,y,z,xold,yold,zold,cif1,cif2,cif3,ind,ind2,ind3,xr,yr,zr) &
     !$omp private(l1,l2,jj1,jj2,jj3a,jj3b,jj3c,ym1,y00,yp1,xya,xa)
     do i3=-nb3,2*n3+1+nb3
        z=real(i3,gp)*hzh
        zold=z-dz
        j3=nint(zold/hzh_old)
        cif3=(j3 >= -6 .and. j3 <= n3o7) .or. perz
        if (.not.cif3) cycle

        zr =real((zold-real(j3,gp)*hzh_old)/hzh_old,wp)
        jj3a=modulo(j3-1+nb3,n3nb3o)-nb3
        jj3b=modulo(j3  +nb3,n3nb3o)-nb3
        jj3c=modulo(j3+1+nb3,n3nb3o)-nb3
        ind3=n1nb1*n2nb2*(i3+nb3)+nb1+1
        do i2=-nb2,2*n2+1+nb2
           y=real(i2,gp)*hyh
           yold=y-dy

           j2=nint(yold/hyh_old)
           cif2=(j2 >= -6 .and. j2 <= n2o7) .or. pery
           if (.not. cif2) cycle

           yr = real((yold-real(j2,gp)*hyh_old)/hyh_old,wp)
           ind2=n1nb1*(i2+nb2)
           do i1=-nb1,2*n1+1+nb1
              x=real(i1,gp)*hxh
              xold=x-dx

              j1=nint(xold/hxh_old)
              cif1=(j1 >= -6 .and. j1 <= n1o7) .or. perx

              if (cif1) then
                 ind=i1+ind2+ind3
                 do l2=-1,1
                    jj2=modulo(j2+l2+nb2,n2nb2o)-nb2
                    do l1=-1,1
                       !the modulo has no effect on free BC thanks to the
                       !if statement above
                       jj1=modulo(j1+l1+nb1,n1nb1o)-nb1

                       ym1=psifscfold(jj1,jj2,jj3a)
                       y00=psifscfold(jj1,jj2,jj3b)
                       yp1=psifscfold(jj1,jj2,jj3c)

                       xya(l1,l2)=ym1 + &
                            (1.0_wp + zr)*(y00 - ym1 + zr*(.5_wp*ym1 - y00  + .5_wp*yp1))
                    enddo
                 enddo

                 do l1=-1,1
                    ym1=xya(l1,-1)
                    y00=xya(l1,0)
                    yp1=xya(l1,1)
                    xa(l1)=ym1 + &
                         (1.0_wp + yr)*(y00 - ym1 + yr*(.5_wp*ym1 - y00  + .5_wp*yp1))
                 enddo

                 xr = real((xold-real(j1,gp)*hxh_old)/hxh_old,wp)
                 ym1=xa(-1)
                 y00=xa(0)
                 yp1=xa(1)
                 psifscf(ind)=ym1 + &
                      (1.0_wp + xr)*(y00 - ym1 + xr*(.5_wp*ym1 - y00  + .5_wp*yp1))

              endif

           enddo
        enddo
     enddo
     !$omp end parallel do

     !call cpu_time(t1)
     !time=real(t1-t0,kind=8)
     !if (bigdft_mpi%iproc==0) print*,'time',time

  endif

  !write(100+iproc,*) 'norm of psifscf ',dnrm2((2*n1+16)*(2*n2+16)*(2*n3+16),psifscf,1)

  call f_free(psifscfold)

  psig = f_malloc((/ 0.to.n1, 1.to.2, 0.to.n2, 1.to.2, 0.to.n3, 1.to.2 /),id='psig')
  ww = f_malloc((2*n1+2+2*nb1)*(2*n2+2+2*nb2)*(2*n3+2+2*nb3),id='ww')

  if (domain_geocode(at%astruct%dom) =='F') then
     call analyse_shrink(n1,n2,n3,ww,psifscf,psig)
  else if (domain_geocode(at%astruct%dom) == 'S') then
     call analyse_slab(n1,n2,n3,ww,psifscf,psig)
  else if (domain_geocode(at%astruct%dom) == 'P') then
     call analyse_per(n1,n2,n3,ww,psifscf,psig)
  else if (domain_geocode(at%astruct%dom) =='W') then
     call analyse_wire(n1,n2,n3,ww,psifscf,psig)
  end if

  !write(100+iproc,*) 'norm new psig ',dnrm2(8*(n1+1)*(n2+1)*(n3+1),psig,1)
  call compress_plain(n1,n2,0,n1,0,n2,0,n3,  &
       wfd%nseg_c,wfd%nvctr_c,wfd%keygloc(1,1),wfd%keyvloc(1),   &
       wfd%nseg_f,wfd%nvctr_f,&
       wfd%keygloc(1,wfd%nseg_c+min(1,wfd%nseg_f)),&
       wfd%keyvloc(wfd%nseg_c+min(1,wfd%nseg_f)),   &
       psig,psi(1),psi(wfd%nvctr_c+min(1,wfd%nvctr_f)))

  !write(100+iproc,*) 'norm of reformatted psi ',dnrm2(nvctr_c+7*nvctr_f,psi,1)

  call f_free(psig)
  call f_free(ww)

  call f_release_routine()

END SUBROUTINE reformatonewave


!> Calculates the minimum difference between two coordinates
!! knowing that there could have been a modulo operation
function mindist(periodic,alat,r,r_old)
  use module_base
  implicit none
  logical, intent(in) :: periodic
  real(gp), intent(in) :: r,r_old,alat
  real(gp) :: mindist

  !for periodic BC calculate mindist only if the center of mass can be defined without the modulo
  if (periodic) then
     if (r_old > 0.5_gp*alat) then
        if (r < 0.5_gp*alat) then
           !mindist=r+alat-r_old
           mindist=0.0_gp
        else
           mindist=r-r_old
        end if
     else
        if (r > 0.5_gp*alat) then
           !mindist=r-alat-r_old
           mindist=0.0_gp
        else
           mindist=r-r_old
        end if
     end if
  else
     mindist=r-r_old
  end if

end function mindist

!> Read one wavefunction
subroutine readonewave(unitwf,useFormattedInput,iorb,iproc,n1,n2,n3,&
     & hx,hy,hz,at,wfd,rxyz_old,rxyz,psi,eval,psifscf)
  use module_base
  use module_types
  use io, only: io_read_descr, io_error, read_psi_compress
  use yaml_output
  use box
  use at_domain
  use compression
  use numerics, only: onehalf,pi
  implicit none
  logical, intent(in) :: useFormattedInput
  integer, intent(in) :: unitwf,iorb,iproc,n1,n2,n3
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(atoms_data), intent(in) :: at
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(wp), intent(out) :: eval
  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
  real(wp), dimension(*), intent(out) :: psifscf !this supports different BC
  !local variables
  character(len=*), parameter :: subname='readonewave'
  character(len = 256) :: error
  logical :: lstat
  integer :: iorb_old,n1_old,n2_old,n3_old,iat,iel,nvctr_c_old,nvctr_f_old,i1,i2,i3
  real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
  real(gp) :: tx,ty,tz,displ,hx_old,hy_old,hz_old!,mindist
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psigold
  type(cell) :: mesh
  type(domain) :: dom

  call f_routine(id=subname)

  !mesh=cell_new(at%astruct%geocode,[n1,n2,n3],[hx,hy,hz])
  dom=domain_new(units=ATOMIC_UNITS,bc=geocode_to_bc_enum(at%astruct%geocode),&
            alpha_bc=onehalf*pi,beta_ac=onehalf*pi,gamma_ab=onehalf*pi,acell=[n1,n2,n3]*[hx,hy,hz])
  mesh=cell_new(dom,[n1,n2,n3],[hx,hy,hz])

  !write(*,*) 'INSIDE readonewave'
  call io_read_descr(unitwf, useFormattedInput, iorb_old, eval, n1_old, n2_old, n3_old, &
       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, at%astruct%nat)
  if (.not. lstat) call io_error(trim(error))
  if (iorb_old /= iorb) stop 'readonewave'


  
  displ=0.0_gp
  do iat=1,at%astruct%nat
    displ=displ+distance(mesh%dom,rxyz(:,iat),rxyz_old(:,iat))**2
  enddo
  displ=sqrt(displ)

  if (hx_old == hx .and. hy_old == hy .and. hz_old == hz .and.&
       nvctr_c_old == wfd%nvctr_c .and. nvctr_f_old == wfd%nvctr_f .and. &
       n1_old == n1  .and. n2_old == n2 .and. n3_old == n3 .and. displ <= 1.d-3) then

     if (iproc == 0) call yaml_map('Need to reformat wavefunctions',.false.)
     call read_psi_compress(unitwf, useFormattedInput, wfd%nvctr_c, wfd%nvctr_f, psi, lstat, error)
     if (.not. lstat) call io_error(trim(error))

  else

     if (iproc == 0 .and. iorb == 1) then
        call yaml_map('Need to reformat wavefunctions',.true.)
        if (hx_old /= hx .or. hy_old /= hy .or. hz_old /= hz) &
           & call yaml_comment('because hgrid_old /= hgrid' // &
             & trim(yaml_toa((/ hx_old,hy_old,hz_old,hx,hy,hz /), fmt='(f14.10)')))
        if (nvctr_c_old /= wfd%nvctr_c) &
           & call yaml_comment('because nvctr_c_old /= nvctr_c' // trim(yaml_toa((/ nvctr_c_old,wfd%nvctr_c /))))
        if (nvctr_f_old /= wfd%nvctr_f) &
           & call yaml_comment('because nvctr_f_old /= nvctr_f' // trim(yaml_toa((/ nvctr_f_old,wfd%nvctr_f /))))
        if (n1_old /= n1  .or. n2_old /= n2 .or. n3_old /= n3 ) &
           call yaml_comment('because cell size has changed' // trim(yaml_toa((/ n1_old,n1,n2_old,n2,n3_old,n3 /))))
        if (displ > 1.d-3 ) call yaml_comment('large displacement of molecule' // trim(yaml_toa(displ)))
     end if

     psigold = f_malloc0((/ 0.to.n1_old, 1.to.2, 0.to.n2_old, 1.to.2, 0.to.n3_old, 1.to.2 /),id='psigold')

     do iel=1,nvctr_c_old
        if (useFormattedInput) then
           read(unitwf,*) i1,i2,i3,tt
        else
           read(unitwf) i1,i2,i3,tt
        end if
        psigold(i1,1,i2,1,i3,1)=tt
     enddo
     do iel=1,nvctr_f_old
        if (useFormattedInput) then
           read(unitwf,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
        else
           read(unitwf) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
        end if
        psigold(i1,2,i2,1,i3,1)=t1
        psigold(i1,1,i2,2,i3,1)=t2
        psigold(i1,2,i2,2,i3,1)=t3
        psigold(i1,1,i2,1,i3,2)=t4
        psigold(i1,2,i2,1,i3,2)=t5
        psigold(i1,1,i2,2,i3,2)=t6
        psigold(i1,2,i2,2,i3,2)=t7
     enddo

     ! I put nat = 1 here, since only one position is saved in wavefunction files.
     call reformatonewave(displ,wfd,at,hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,&
          rxyz_old,psigold,hx,hy,hz,n1,n2,n3,rxyz,psifscf,psi)

     call f_free(psigold)

  endif

  call f_release_routine()

END SUBROUTINE readonewave

subroutine readwavedescr(lstat, filename, iorb, ispin, ikpt, ispinor, nspinor, fileRI)
  use module_base
  use module_types

  implicit none

  character(len = *), intent(in) :: filename
  integer, intent(out) :: iorb, ispin, ikpt, nspinor, ispinor
  logical, intent(out) :: lstat
  character(len = 1024), intent(out) :: fileRI

  logical :: exists
  integer :: i, i_stat
  character(len = 1) :: code

  lstat = .false.

  !find the value of iorb
  read(filename(index(filename, ".", back = .true.)+2:len(filename)),*,iostat = i_stat) iorb
  if (i_stat /= 0) return
  i = index(filename, "-k", back = .true.)+2
  read(filename(i:i+2),*,iostat = i_stat) ikpt
  if (i_stat /= 0) return
  i = index(filename, "-", back = .true.)+1
  read(filename(i:i),*,iostat = i_stat) code
  if (i_stat /= 0) return
  if (code == "U" .or. code == "N") ispin = 1
  if (code == "D") ispin = 2
  ! Test file for other spinor part.
  nspinor = 1
  ispinor = 1
  read(filename(i+1:i+1),*,iostat = i_stat) code
  if (i_stat /= 0) return
  if (code == "R") then
     write(fileRI, "(A)") filename
     fileRI(i+1:i+1) = "I"
     inquire(file=trim(fileRI), exist=exists)
     if (exists) then
        ispinor = 1
        nspinor = 2
     end if
  end if
  if (code == "I") then
     write(fileRI, "(A)") filename
     fileRI(i+1:i+1) = "R"
     inquire(file=trim(fileRI), exist=exists)
     if (exists) then
        ispinor = 2
        nspinor = 2
     end if
  end if

  lstat = .true.
END SUBROUTINE readwavedescr


!> Given a translation vector, find the inverse one
subroutine find_inverse(nin,iout,t0_field,t0_l,k1)
  use module_base
  use yaml_output
  implicit none
  integer, intent(in) :: iout                      !< Point of the new grid from which the inverse has to be found
  integer, intent(in) :: nin                       !< Number of points of the input grid
  real(gp), dimension(nin), intent(in) :: t0_field !< Array of displacements of the input grid
  integer, intent(out) :: k1                       !< Starting point of the input grid from which the interplation should be calculated
  real(gp), intent(out) :: t0_l                    !< Resulting shift from the starting point, from which the filter has to be calculated
  !local variables
  real(gp), parameter  :: tol=1.e-14_gp
  integer :: l,k2
  real(gp) :: ksh1,ksh2,k,kold,alpha

  kold=-1000.0_gp
  find_trans: do l=1,nin
     k=real(l,gp)+t0_field(l)
     if (k-real(iout,gp) > tol) exit find_trans
     kold=k
  end do find_trans
  ! want to use either l or l-1 to give us point i - pick closest
  if (k-real(iout,gp) < -kold+real(iout,gp)) then
     ksh1=k-real(iout,gp)
     ksh2=-kold+real(iout,gp)
     k1=l
     k2=l-1
     if (k2==0) then
        k2=1
        ksh2=ksh1
     end if
     if (k1==nin+1) then
        k1=nin
        ksh1=ksh2
     end if
  else
     ksh1=-kold+real(iout,gp)
     ksh2=k-real(iout,gp)
     k1=l-1
     k2=l
     if (k1==0) then
        k1=1
        ksh1=ksh2
     end if
     if (k2==nin+1) then
        k2=nin
        ksh2=ksh1
     end if
  end if

  if (ksh1==0.0_gp .or. k1==k2) then !otherwise already have exactly on point
     ksh2=1.0_gp
     ksh1=0.0_gp
  end if

  alpha=ksh2/(ksh1+ksh2)

  t0_l=alpha*t0_field(k1)+(1.0_gp-alpha)*t0_field(k2)

end subroutine find_inverse


 !> routine which directly applies the 3D transformation of the rototranslation
subroutine field_rototranslation3D_interpolation(da,newz,centre_old,centre_new,&
     sint,cost,onemc,hgrids_old,ndims_old,f_old,&
     hgrids_new,ndims_new,f_new)
  use module_base
  use yaml_output
  implicit none
  real(gp), intent(in) :: sint,cost,onemc !< rotation wrt newzeta vector
  real(gp), dimension(3), intent(in) :: da !<coordinates of rigid shift vector
  real(gp), dimension(3), intent(in) :: newz !<coordinates of new z vector (should be of norm one)
  real(gp), dimension(3), intent(in) :: centre_old,centre_new !<centre of rotation
  real(gp), dimension(3), intent(in) :: hgrids_old,hgrids_new !<dimension of old and new box
  integer, dimension(3), intent(in) :: ndims_old,ndims_new !<dimension of old and new box
  real(gp), dimension(ndims_old(1),ndims_old(2),ndims_old(3)), intent(in) :: f_old
  real(gp), dimension(ndims_new(1),ndims_new(2),ndims_new(3)), intent(out) :: f_new
  !local variables
  integer :: i,j,k,it
  real(gp), dimension(3) :: dt
  real(gp), dimension(27) :: coeffs

  call f_routine(id='field_rototranslation3D_interpolation')

  !loop on the coordinates of the new domain
  do k=1,ndims_new(3)
     do j=1,ndims_new(2)
        do i=1,ndims_new(1)
           do it=1,3
              call shift_and_start(i,j,k,coeffs,dt)
           end do
!           print *,'i,j,k',i,j,k,coeffs
!           print *,'dt',dt
           f_new(i,j,k)=interpolate(dt,coeffs)
!           print *,'interpolate',f_new(i,j,k)
        end do
     end do
  end do

  call f_release_routine()
!stop
contains

  function interpolate(dt,aijk)
    implicit none
    real(gp), dimension(3), intent(in) :: dt
    real(gp), dimension(0:2,0:2,0:2), intent(inout) :: aijk
    real(gp) :: interpolate
    !local variables
    integer :: px,py,pz,ix,iy,iz,info
    real(gp) :: x,y,z
    integer, dimension(27) :: ipiv
    real(gp), dimension(-1:1,-1:1,-1:1,0:2,0:2,0:2) :: bijk

    if (maxval(abs(aijk)) == 0.0_gp) then
       interpolate=0.0_gp
       return
    end if

    do iz=-1,1
       z=dt(3)+real(iz,gp)
       z=hgrids_old(3)*z
       do iy=-1,1
          y=dt(2)+real(iy,gp)
          y=hgrids_old(2)*y
          do ix=-1,1
             x=dt(1)+real(ix,gp)
             x=hgrids_old(1)*x
             do pz=0,2
                do py=0,2
                   do px=0,2
                      bijk(ix,iy,iz,px,py,pz)=(x**px)*(y**py)*(z**pz)
                   end do
                end do
             end do
          end do
       end do
    end do

    !here the linear system has to be solved to find the coefficients aijk
    !some pragma has to be passed to MKL to ensure a monothread execution
    call dgesv(27,1,bijk,27,ipiv,aijk,27,info)
    if (info /=0) then
       print *,'error', info, dt
       call f_err_severe()
    end if
    interpolate=aijk(0,0,0)

  end function interpolate

  pure subroutine shift_and_start(j1,j2,j3,fijk,dt)
    implicit none
    integer, intent(in) :: j1,j2,j3
    real(gp), dimension(-1:1,-1:1,-1:1), intent(out) :: fijk
    real(gp), dimension(3), intent(out) :: dt
    !local variables
    integer :: ix,iy,iz
    integer, dimension(3) :: istart,istart_shift
    real(gp), dimension(3) :: t0_l

    !define the coordinates in the reference frame
    !which depends of the transformed variables
    dt(1)=-centre_new(1)+real(j1-1,gp)*hgrids_new(1)
    dt(2)=-centre_new(2)+real(j2-1,gp)*hgrids_new(2)
    dt(3)=-centre_new(3)+real(j3-1,gp)*hgrids_new(3)

    !define the value of the shift of the variable we are going to transform
    t0_l=coord(newz,cost,sint,onemc,dt(1),dt(2),dt(3))-da
    istart=nint((t0_l+centre_old+hgrids_old)/hgrids_old)

    !!doubts about that
    !t0_l=(dt-t0_l)/hgrids_old
    !!identify shift
    !dt(1)=(real(istart(1),gp)+t0_l(1))-real(j1,gp)*hgrids_new(1)/hgrids_old(1)
    !dt(2)=(real(istart(2),gp)+t0_l(2))-real(j2,gp)*hgrids_new(2)/hgrids_old(2)
    !dt(3)=(real(istart(3),gp)+t0_l(3))-real(j3,gp)*hgrids_new(3)/hgrids_old(3)
    !!end of doubts


    !this shift brings the old point in the new reference frame
    dt=real(istart,gp)-(t0_l+centre_new+hgrids_new)/hgrids_old

    !purify the shift to be a inferior than multiple of the grid spacing
    istart_shift=nint(dt)
    dt=dt-real(istart_shift,gp)
    istart=istart-istart_shift

    !fill array if it is inside the old box
    fijk=0.0_gp
    do iz=-1,1
       if (istart(3)+iz >= 1 .and. istart(3)+iz <= ndims_old(3)) then
          do iy=-1,1
             if (istart(2)+iy >= 1 .and. istart(2)+iy <= ndims_old(2)) then
             do ix=-1,1
                if (istart(1)+ix >= 1 .and. istart(1)+ix <= ndims_old(1)) then
                   fijk(ix,iy,iz)=&
                        f_old(istart(1)+ix,istart(2)+iy,istart(3)+iz)
                end if
             end do
          end if
          end do
       end if
    end do

!    if (maxval(abs(fijk)) /= 0.0_gp) then
!       write(17,*)j1,j2,j3,dt,istart,fijk
!    end if


  end subroutine shift_and_start

  pure function coord(u,C,S,onemc,x,y,z)
    use module_base, only: gp
    implicit none
    real(gp), intent(in) :: C,S,onemc !<trigonometric functions of the theta angle
    real(gp), intent(in) :: x,y,z !<coordinates to be used for the mapping
    real(gp), dimension(3), intent(in) :: u !<axis of rotation
    real(gp), dimension(3) :: coord

    coord(1)=u(1)**2*x + u(1)*u(2)*y + S*u(3)*y - S*u(2)*z + u(1)*u(3)*z - C*((-1 + u(1)**2)*x + u(1)*(u(2)*y + u(3)*z))
    coord(2)=-(S*u(3)*x) + (C + onemc*u(2)**2)*y + onemc*u(2)*u(3)*z + u(1)*(u(2)*onemc*x + S*z)
    coord(3)=S*(u(2)*x - u(1)*y) + C*z + u(3)*(onemc*u(1)*x + onemc*u(2)*y + u(3)*z - C*u(3)*z)

  end function coord

end subroutine field_rototranslation3D_interpolation
