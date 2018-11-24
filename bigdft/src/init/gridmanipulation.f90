!> @file
!!  Routines to manipulate the grid
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculates the overall size of the simulation cell 
!! and shifts the atoms such that their position is the most symmetric possible wrt the simulation cell.
!! Assign these values to the global localisation region descriptor.
subroutine system_size(atoms,rxyz,crmult,frmult,hx,hy,hz,OCLconv,Glr)
   use module_base
   use module_types
   use yaml_strings, only: yaml_toa
   use locregs
   use at_domain
   use box, only: cell,box_nbox_from_cutoff,cell_new
   use numerics, only: onehalf,pi
   implicit none
   type(atoms_data), intent(inout) :: atoms
   real(gp), intent(in) :: crmult,frmult
   real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
   !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
   real(gp), intent(inout) :: hx,hy,hz
   logical, intent(in) :: OCLconv
   type(locreg_descriptors), intent(out) :: Glr
   !Local variables
   !character(len=*), parameter :: subname='system_size'
   type(cell) :: mesh_coarse
   integer, parameter :: lupfil=14,START_=1,END_=2
   real(gp), parameter :: eps_mach=1.e-12_gp
   integer :: iat,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,i!,n1i,n2i,n3i
   real(gp) :: ri,rad,cxmin,cxmax,cymin,cymax,czmin,czmax!,alatrue1,alatrue2,alatrue3
   real(gp), dimension(3) :: hgridsh,alat
   logical, dimension(3) :: peri
   integer, dimension(3) :: ndims,ndims1
   integer, dimension(2,3) :: nbox,nbox_tmp
   type(domain) :: dom

   !check the geometry code with the grid spacings
   if (domain_geocode(atoms%astruct%dom) == 'F' .and. (hx/=hy .or. hx/=hz .or. hy/=hz)) then
      call f_err_throw('Grid spacings must be equal' // &
           ' in the Free BC case, while hgrids = '//&
           trim(yaml_toa((/ hx, hy, hz/),fmt='(f7.4)')),&
           err_name='BIGDFT_INPUT_VARIABLES_ERROR')
      !write(*,'(1x,a,3(1x,F7.4))') 'ERROR: The values of the grid spacings must be equal' // &
      !     & ' in the Free BC case, while hgrids = ', hx, hy, hz
      return
   end if

   !Special case if no atoms (no posinp by error or electron gas)
   if (atoms%astruct%nat == 0) then
      ri = 0.0_gp
   else
      ri = 1.e10_gp
   end if

   !calculate the extremes of the boxes taking into account the spheres around the atoms
   cxmax = -ri 
   cxmin =  ri

   cymax = -ri 
   cymin =  ri

   czmax = -ri 
   czmin =  ri

   !this treatment is a pre-evaluation of the box dimensions
   !as if all the directions were Free BC.
   !therefore there is no need to intriduce the metric here
   do iat=1,atoms%astruct%nat

      rad=atoms%radii_cf(atoms%astruct%iatype(iat),1)*crmult

      cxmax=max(cxmax,rxyz(1,iat)+rad) 
      cxmin=min(cxmin,rxyz(1,iat)-rad)

      cymax=max(cymax,rxyz(2,iat)+rad) 
      cymin=min(cymin,rxyz(2,iat)-rad)

      czmax=max(czmax,rxyz(3,iat)+rad) 
      czmin=min(czmin,rxyz(3,iat)-rad)
   enddo

   !eliminate epsilon form the grid size calculation
   !!  cxmax=cxmax+eps_mach 
   !!  cymax=cymax+eps_mach  
   !!  czmax=czmax+eps_mach  
   !!
   !!  cxmin=cxmin-eps_mach
   !!  cymin=cymin-eps_mach
   !!  czmin=czmin-eps_mach

   !peri=bc_periodic_dims(geocode_to_bc(atoms%astruct%geocode))
   peri=domain_periodic_dims(atoms%astruct%dom)
   alat(1)=(cxmax-cxmin)
   alat(2)=(cymax-cymin)
   alat(3)=(czmax-czmin)
   hgridsh(1)=hx
   hgridsh(2)=hy
   hgridsh(3)=hz
   !assign the shift to the atomic positions
   atoms%astruct%shift(1)=cxmin
   atoms%astruct%shift(2)=cymin
   atoms%astruct%shift(3)=czmin
   do i=1,3
      if (peri(i)) then
         call correct_grid(atoms%astruct%cell_dim(i),hgridsh(i),ndims(i))
         atoms%astruct%shift(i)=0.0_gp
      else
         atoms%astruct%cell_dim(i)=alat(i)
         ndims(i)=int(atoms%astruct%cell_dim(i)/hgridsh(i))
         alat(i)=ndims(i)*hgridsh(i)
         atoms%astruct%shift(i)=atoms%astruct%shift(i)+&
              0.5_gp*(atoms%astruct%cell_dim(i)-alat(i))
         atoms%astruct%cell_dim(i)=alat(i)
      end if
   end do
!!$   alatrue1=alat(1)
!!$   alatrue2=alat(2)
!!$   alatrue3=alat(3)
   hx=hgridsh(1)
   hy=hgridsh(2)
   hz=hgridsh(3)

   n1=ndims(1)
   n2=ndims(2)
   n3=ndims(3)

!!$   !correct the box sizes for the isolated case
!!$   select case(atoms%astruct%geocode)
!!$   case('F')
!!$      atoms%astruct%cell_dim(1)=alatrue1
!!$      atoms%astruct%cell_dim(2)=alatrue2
!!$      atoms%astruct%cell_dim(3)=alatrue3
!!$   case('S')
!!$      cxmin=0.0_gp
!!$      atoms%astruct%cell_dim(2)=alatrue2
!!$      czmin=0.0_gp
!!$   case('P')
!!$      !for the moment we do not put the shift, at the end it will be tested
!!$      !here we should put the center of mass
!!$      cxmin=0.0_gp
!!$      cymin=0.0_gp
!!$      czmin=0.0_gp
!!$   end select
!!$
!!$   !assign the shift to the atomic positions
!!$   atoms%astruct%shift(1)=cxmin
!!$   atoms%astruct%shift(2)=cymin
!!$   atoms%astruct%shift(3)=czmin

   !here we can put a modulo operation for periodic directions
   do iat=1,atoms%astruct%nat
      rxyz(1,iat)=rxyz(1,iat)-atoms%astruct%shift(1)
      rxyz(2,iat)=rxyz(2,iat)-atoms%astruct%shift(2)
      rxyz(3,iat)=rxyz(3,iat)-atoms%astruct%shift(3)
   enddo

   !put the maximum and the miniumum inverted in the beginning
   nbox(START_,:)=ndims
   nbox(END_,:)=0

   ! fine grid size (needed for creation of input wavefunction, preconditioning)
   if (atoms%astruct%nat == 0) then
      !For homogeneous gaz, we fill the box with the fine grid
      nbox(START_,:)=0
      nbox(END_,:)=ndims
!!$      nfl1=0 
!!$      nfl2=0 
!!$      nfl3=0
!!$
!!$      nfu1=n1
!!$      nfu2=n2
!!$      nfu3=n3
!!$   else
!!$      !we start with nfl max to find th emin and nfu min to find the max
!!$      nfl1=n1 
!!$      nfl2=n2 
!!$      nfl3=n3
!!$
!!$      nfu1=0 
!!$      nfu2=0 
!!$      nfu3=0
   end if

   !here we will put the nformations about the cell angles
   !we will also have to decide where the mesh should be in astruct or not
   !mesh_coarse=cell_new(atoms%astruct%geocode,ndims+1,hgridsh)
   ndims1=ndims+1
   dom=domain_new(units=ATOMIC_UNITS,bc=geocode_to_bc_enum(domain_geocode(atoms%astruct%dom)),&
            alpha_bc=onehalf*pi,beta_ac=onehalf*pi,gamma_ab=onehalf*pi,acell=ndims1*hgridsh)
   mesh_coarse=cell_new(dom,ndims1,hgridsh)
  
   do iat=1,atoms%astruct%nat
      rad=atoms%radii_cf(atoms%astruct%iatype(iat),2)*frmult
      if (rad <= 0.0_gp) cycle
      nbox_tmp=box_nbox_from_cutoff(mesh_coarse,rxyz(:,iat),rad+eps_mach*maxval(hgridsh))
      do i=1,3
         nbox(START_,i)=min(nbox(START_,i),nbox_tmp(START_,i))
         nbox(END_,i)=max(nbox(END_,i),nbox_tmp(END_,i))
      end do
!!$      if (rad > 0.0_gp) then
!!$         nfl1=min(nfl1,ceiling((rxyz(1,iat)-rad)/hx - eps_mach))
!!$         nfu1=max(nfu1,floor((rxyz(1,iat)+rad)/hx + eps_mach))
!!$
!!$         nfl2=min(nfl2,ceiling((rxyz(2,iat)-rad)/hy - eps_mach))
!!$         nfu2=max(nfu2,floor((rxyz(2,iat)+rad)/hy + eps_mach))
!!$
!!$         nfl3=min(nfl3,ceiling((rxyz(3,iat)-rad)/hz - eps_mach)) 
!!$         nfu3=max(nfu3,floor((rxyz(3,iat)+rad)/hz + eps_mach))
!!$      end if
   enddo

   do i=1,3
      if (nbox(START_,i) < 0 .or. nbox(END_,i) > ndims(i)) then
         nbox(START_,i)=0
         nbox(END_,i)=ndims(i)
      end if
      if (nbox(START_,i) == ndims(i) .and. nbox(END_,i) == 0) nbox(:,i)=ndims(i)/2
   end do

!!$   !correct the values of the delimiter if they go outside the box
!!$   if (nfl1 < 0 .or. nfu1 > n1) then
!!$      nfl1=0
!!$      nfu1=n1
!!$   end if
!!$   if (nfl2 < 0 .or. nfu2 > n2) then
!!$      nfl2=0
!!$      nfu2=n2
!!$   end if
!!$   if (nfl3 < 0 .or. nfu3 > n3) then
!!$      nfl3=0
!!$      nfu3=n3
!!$   end if
!!$
!!$   !correct the values of the delimiter if there are no wavelets
!!$   if (nfl1 == n1 .and. nfu1 == 0) then
!!$      nfl1=n1/2
!!$      nfu1=n1/2
!!$   end if
!!$   if (nfl2 == n2 .and. nfu2 == 0) then
!!$      nfl2=n2/2
!!$      nfu2=n2/2
!!$   end if
!!$   if (nfl3 == n3 .and. nfu3 == 0) then
!!$      nfl3=n3/2
!!$      nfu3=n3/2
!!$   end if

   hgridsh(1)=0.5_gp*hx
   hgridsh(2)=0.5_gp*hy
   hgridsh(3)=0.5_gp*hz

   nfl1=nbox(START_,1)
   nfu1=nbox(END_,1)   
   nfl2=nbox(START_,2) 
   nfu2=nbox(END_,2)   
   nfl3=nbox(START_,3) 
   nfu3=nbox(END_,3) 

   !assign the values
!   call init_lr(Glr,mesh_coarse,nbox_fine,&
   call init_lr(Glr,domain_geocode(atoms%astruct%dom),hgridsh,n1,n2,n3,&
        nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,&
        hybrid_flag=.not. OCLconv)

END SUBROUTINE system_size


!> Here the dimensions should be corrected in order to 
!! allow the fft for the preconditioner and for Poisson Solver
subroutine correct_grid(a,h,n)
   use module_base
   use Poisson_Solver, except_dp => dp, except_gp => gp
   implicit none
   real(gp), intent(in) :: a
   integer, intent(out) :: n
   real(gp), intent(inout) :: h
   !local variables
   integer :: m,m2,nt

   n=ceiling(a/h)-1
   nt=n+1
   do
      !correct the direct dimension
      call fourier_dim(nt,m)

      !control if the double of this dimension is compatible with the FFT
      call fourier_dim(2*m,m2)
      !if this check is passed both the preconditioner and the PSolver works
      if (m2==2*m .and. mod(m,2) ==0) exit !only even dimensions are considered so far

      nt=m+1
   end do
   n=m-1

   !!!  !here the dimensions should be corrected in order to 
   !!!  !allow the fft for the preconditioner
   !!!  m=2*n+2
   !!!  do 
   !!!     call fourier_dim(m,m)
   !!!     if ((m/2)*2==m) then
   !!!        n=(m-2)/2
   !!!        exit
   !!!     else
   !!!        m=m+1
   !!!     end if
   !!!  end do

   h=a/real(n+1,gp)

END SUBROUTINE correct_grid

subroutine export_grids(fname, atoms, rxyz, hx, hy, hz, n1, n2, n3, logrid_c, logrid_f)
  use module_defs, only: gp
  use module_types
  use module_base, only: f_err_throw
  use module_atoms, only: write_xyz_bc
  use at_domain, only: domain_geocode
  implicit none
  character(len = *), intent(in) :: fname
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
  real(gp), intent(in) :: hx, hy, hz
  integer, intent(in) :: n1, n2, n3
  logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid_c
  logical, dimension(0:n1,0:n2,0:n3), intent(in), optional :: logrid_f

  integer :: nvctr, iat, i3, i2, i1

  nvctr = 0
  do i3=0,n3  
     do i2=0,n2  
        do i1=0,n1
           if (logrid_c(i1,i2,i3)) nvctr = nvctr + 1
        enddo
     enddo
  end do
  if (present(logrid_f)) then
     do i3=0,n3  
        do i2=0,n2  
           do i1=0,n1
              if (logrid_f(i1,i2,i3)) nvctr = nvctr + 1
           enddo
        enddo
     end do
  end if

  ! Create the file grid.xyz to visualize the grid of functions
  open(unit=22,file=fname,status='unknown')
  write(22,*) nvctr+atoms%astruct%nat,' atomic'

  call write_xyz_bc(22,domain_geocode(atoms%astruct%dom),1.0_gp,atoms%astruct%cell_dim)

!!$  if (atoms%astruct%geocode=='F') then
!!$     write(22,*)'complete simulation grid with low and high resolution points'
!!$  else if (atoms%astruct%geocode =='S') then
!!$     write(22,'(a,2x,3(1x,1pe24.17))')'surface',atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),atoms%astruct%cell_dim(3)
!!$  else if (atoms%astruct%geocode =='P') then
!!$     write(22,'(a,2x,3(1x,1pe24.17))')'periodic',atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
!!$          atoms%astruct%cell_dim(3)
!!$  else if (atoms%astruct%geocode =='W') then
!!$     call f_err_throw("Wires bc has to be implemented here", &
!!$          err_name='BIGDFT_RUNTIME_ERROR')
!!$  end if
!!$
  do iat=1,atoms%astruct%nat
     write(22,'(a6,2x,3(1x,e12.5),3x)') &
          &   trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
  enddo
  do i3=0,n3  
     do i2=0,n2  
        do i1=0,n1
           if (logrid_c(i1,i2,i3))&
                &   write(22,'(a4,2x,3(1x,e10.3))') &
                &   '  g ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
        enddo
     enddo
  end do
  if (present(logrid_f)) then
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              if (logrid_f(i1,i2,i3))&
                   &   write(22,'(a4,2x,3(1x,e10.3))') &
                   &   '  G ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           enddo
        enddo
     enddo
  end if
  close(22)
END SUBROUTINE export_grids
