!> @file
!!  Routines to do tail calculation (correct effect of finite size)
!! @author
!!    Copyright (C) 2007-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculate the finite size corrections over wavefunctions
!! Conceived only for isolated Boundary Conditions, no SIC correction
subroutine CalculateTailCorrection(iproc,nproc,at,rbuf,orbs,&
     Glr,nlpsp,ncongt,pot,hgrid,rxyz,crmult,frmult,nspin,&
     psi,output_denspot,ekin_sum,epot_sum,eproj_sum,paw)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, only: orbitals_descriptors
  use psp_projectors_base
  use gaussians, only: gaussian_basis
  use public_enums
  use bounds, only: make_bounds, make_all_ib
  use locregs
  use compression
  use orbitalbasis
  use psp_projectors
  implicit none
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: Glr
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  integer, intent(in) :: iproc,nproc,ncongt,nspin
  logical, intent(in) :: output_denspot
  real(kind=8), dimension(3), intent(in) :: hgrid
  real(kind=8), intent(in) :: crmult,frmult,rbuf
  !real(kind=8), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
  real(kind=8), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(kind=8), dimension(Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nspin), intent(in) :: pot
  real(kind=8), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbs%norbp), intent(in) :: psi
  real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
  type(paw_objects),intent(inout)::paw
  !local variables
  type(locreg_descriptors), target :: lr
  character(len=*), parameter :: subname='CalculateTailCorrection'
  integer :: iseg,i0,j0,i1,j1,i2,i3,ii,iat,iorb,npt,ipt,i,ierr,nbuf,ispin
  integer :: nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3
  integer :: n1,n2,n3,nsegb_c,nsegb_f,nvctrb_c,nvctrb_f
  real(kind=8) :: alatb1,alatb2,alatb3,ekin,epot,eproj,tt,cprecr,sum_tail !n(c) eproj1 epot1,ekin1
  type(orbitals_data), target :: orbsb
  type(orbital_basis), target :: ob
  type(ket) :: psi_it
  type(DFT_PSP_projector_iter) :: psp_it
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f
  integer, dimension(:,:,:), allocatable :: ibbyz_c,ibbyz_f,ibbxz_c,ibbxz_f,ibbxy_c,ibbxy_f
  real(kind=8), dimension(:,:), allocatable :: wrkallred
  real(kind=8), dimension(:), allocatable, target :: psib,hpsib,psir
  real(kind=8), dimension(:,:), pointer :: txyz

  !for shrink:
  integer, allocatable, dimension(:,:,:) :: ibbzzx_c,ibbyyzz_c
  integer, allocatable, dimension(:,:,:) :: ibbxy_ff,ibbzzx_f,ibbyyzz_f

  !for grow:
  integer, allocatable, dimension(:,:,:) :: ibbzxx_c,ibbxxyy_c
  integer, allocatable, dimension(:,:,:) :: ibbyz_ff,ibbzxx_f,ibbxxyy_f

  !real space border:
  integer, allocatable, dimension(:,:,:) :: ibbyyzz_r 

  integer nw1,nw2, nwarnings

  real(kind=8), dimension(:,:,:), allocatable::x_c!input 
  real(kind=8), dimension(:,:,:,:), allocatable :: x_f ! input
  real(kind=8), dimension(:,:,:), allocatable :: x_f1,x_f2,x_f3 ! internal
  real(kind=8), dimension(:), allocatable :: w1,w2
  real(kind=8), dimension(:,:,:), allocatable::y_c!output 
  real(kind=8), dimension(:,:,:,:), allocatable :: y_f! output

  !SM This routine had as argument hgrid(1), but must now have all three dur to
  !the modified calling sequence of convolutions routine. Still only use the
  !first entry throughout the routine.

  call f_routine(id='CalculateTailCorrection')

  n1=Glr%d%n1
  n2=Glr%d%n2
  n3=Glr%d%n3

  nbuf=nint(rbuf/hgrid(1))
  !    --- new grid sizes n1,n2,n3
  nb1=n1+2*nbuf
  nb2=n2+2*nbuf
  nb3=n3+2*nbuf

  ! Create new structure with modified grid sizes
  call init_lr(lr,Glr%geocode,0.5*hgrid,nb1,nb2,nb3,&
       Glr%d%nfl1,Glr%d%nfl2,Glr%d%nfl3,Glr%d%nfu1,Glr%d%nfu2,Glr%d%nfu3,&
       .true.,bnds=Glr%bounds)
 
  alatb1=real(nb1,kind=8)*hgrid(1) 
  alatb2=real(nb2,kind=8)*hgrid(1) 
  alatb3=real(nb3,kind=8)*hgrid(1)

  if (iproc == 0) then
     call yaml_comment('Finite-Size correction',hfill='-')
     call yaml_mapping_open('Estimation of Finite-Size Corrections')
     call yaml_map('Effective AU space more around each external atom',rbuf,fmt='(f6.3)')
     call yaml_map('Adding grid points around cell',nbuf)
     call yaml_map('Effective box size (AU)', (/ alatb1,alatb2,alatb3 /), fmt='(1pe12.5)')
     call yaml_map('Grid spacing units', (/ nb1,nb2,nb3 /))
     !write(*,'(1x,a)')&
     !     '---------------------------------------------- Estimation of Finite-Size Corrections'
     !write(*,'(1x,a,f6.3,a)') &
     !     'F-S Correction for an effective space of ',rbuf,' AU more around each external atom'
     !write(*,'(1x,a,i6,a)') &
     !     '                  requires the adding of ',nbuf,' additional grid points around cell'
     !write(*,'(1x,a,19x,a)') &
     !     '   Effective box size,   Atomic Units:','grid spacing units:'
     !write(*,'(1x,a,3(1x,1pe12.5),3x,3(1x,i9))')&
     !     '            ',alatb1,alatb2,alatb3,nb1,nb2,nb3
  end if

!end do

  !---reformat wavefunctions

  ! fine grid size (needed for creation of input wavefunction, preconditioning)
  nbfl1=Glr%d%nfl1+nbuf ; nbfl2=Glr%d%nfl2+nbuf ; nbfl3=Glr%d%nfl3+nbuf
  nbfu1=Glr%d%nfu1+nbuf ; nbfu2=Glr%d%nfu2+nbuf ; nbfu3=Glr%d%nfu3+nbuf
  if (iproc == 0) then
     call yaml_map('Extremes for the new high resolution grid points', &
          & (/ nbfl1, nbfu1, nbfl2,nbfu2,nbfl3,nbfu3 /))

     !write(*,'(1x,a,3x,3(3x,i4,a1,i0))')&
     !     '  Extremes for the new high resolution grid points:',&
     !     nbfl1,'<',nbfu1,nbfl2,'<',nbfu2,nbfl3,'<',nbfu3
  endif

  ! change atom coordinates according to the enlarged box
  txyz = f_malloc_ptr((/ 3, at%astruct%nat /),id='txyz')
  do iat=1,at%astruct%nat
     txyz(1,iat)=rxyz(1,iat)+real(nbuf,kind=8)*hgrid(1)
     txyz(2,iat)=rxyz(2,iat)+real(nbuf,kind=8)*hgrid(1)
     txyz(3,iat)=rxyz(3,iat)+real(nbuf,kind=8)*hgrid(1)
  enddo

  ! determine localization region for all orbitals, but do not yet fill the descriptor arrays
  logrid_c = f_malloc((/ 0.to.nb1, 0.to.nb2, 0.to.nb3 /),id='logrid_c')
  logrid_f = f_malloc((/ 0.to.nb1, 0.to.nb2, 0.to.nb3 /),id='logrid_f')
  ibbyz_c = f_malloc((/ 1.to.2, 0.to.nb2, 0.to.nb3 /),id='ibbyz_c')
  ibbxz_c = f_malloc((/ 1.to.2, 0.to.nb1, 0.to.nb3 /),id='ibbxz_c')
  ibbxy_c = f_malloc((/ 1.to.2, 0.to.nb1, 0.to.nb2 /),id='ibbxy_c')
  ibbyz_f = f_malloc((/ 1.to.2, 0.to.nb2, 0.to.nb3 /),id='ibbyz_f')
  ibbxz_f = f_malloc((/ 1.to.2, 0.to.nb1, 0.to.nb3 /),id='ibbxz_f')
  ibbxy_f = f_malloc((/ 1.to.2, 0.to.nb1, 0.to.nb2 /),id='ibbxy_f')

  !   allocate for grow
  ibbzxx_c = f_malloc((/ 1.to.2, 0.to.nb3, -14.to.2*nb1+16 /),id='ibbzxx_c')
  ibbxxyy_c = f_malloc((/ 1.to.2, -14.to.2*nb1+16, -14.to.2*nb2+16 /),id='ibbxxyy_c')
  ibbyz_ff = f_malloc((/ 1.to.2, 0.to.nb2, 0.to.nb3 /),id='ibbyz_ff')
  ibbzxx_f = f_malloc((/ 1.to.2, 0.to.nb3, -14.to.2*nb1+16 /),id='ibbzxx_f')
  ibbxxyy_f = f_malloc((/ 1.to.2, -14.to.2*nb1+16, -14.to.2*nb2+16 /),id='ibbxxyy_f')

  !allocate for shrink
  ibbzzx_c = f_malloc((/ 1.to.2, -14.to.2*nb3+16, 0.to.nb1 /),id='ibbzzx_c')
  ibbyyzz_c = f_malloc((/ 1.to.2, -14.to.2*nb2+16, -14.to.2*nb3+16 /),id='ibbyyzz_c')
  ibbxy_ff = f_malloc((/ 1.to.2, 0.to.nb1, 0.to.nb2 /),id='ibbxy_ff')
  ibbzzx_f = f_malloc((/ 1.to.2, -14.to.2*nb3+16, 0.to.nb1 /),id='ibbzzx_f')
  ibbyyzz_f = f_malloc((/ 1.to.2, -14.to.2*nb2+16, -14.to.2*nb3+16 /),id='ibbyyzz_f')

  !allocate for real space
  ibbyyzz_r = f_malloc((/ 1.to.2, -14.to.2*nb2+16, -14.to.2*nb3+16 /),id='ibbyyzz_r')

  ! coarse grid quantities
  call fill_logrid('F',nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,nbuf,at%astruct%nat,at%astruct%ntypes,at%astruct%iatype,txyz, & 
       at%radii_cf(1,1),crmult,hgrid(1),hgrid(1),hgrid(1),logrid_c)
  call num_segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_c,nsegb_c,nvctrb_c)
  lr%wfd%nseg_c = nsegb_c
  lr%wfd%nvctr_c = nvctrb_c

  if (iproc == 0) then
     call yaml_mapping_open('Coarse resolution grid',flow=.true.)
     call yaml_map('Segments',nsegb_c)
     call yaml_map('Points',nvctrb_c)
     call yaml_mapping_close()
     !write(*,'(2(1x,a,i10))') &
     !     'Coarse resolution grid: Number of segments= ',nsegb_c,'points=',nvctrb_c
     !write(*,'(1x,a,2(1x,i10))') 'BIG: orbitals have coarse segment, elements',nsegb_c,nvctrb_c
  end if
  call make_bounds(nb1,nb2,nb3,logrid_c,ibbyz_c,ibbxz_c,ibbxy_c)

  ! fine grid quantities
  call fill_logrid('F',nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,0,at%astruct%nat,at%astruct%ntypes,at%astruct%iatype,txyz, & 
       at%radii_cf(1,2),frmult,hgrid(1),hgrid(1),hgrid(1),logrid_f)
  call num_segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_f,nsegb_f,nvctrb_f)
  lr%wfd%nseg_f = nsegb_f
  lr%wfd%nvctr_f = nvctrb_f
  if (iproc == 0) then
     !Bug in yaml_output solved
     call yaml_mapping_open('Fine resolution grid',flow=.true.)
     call yaml_map('Segments',nsegb_f)
     call yaml_map('Points',nvctrb_f)
     call yaml_mapping_close()
     !write(*,'(2(1x,a,i10))') &
     !     '  Fine resolution grid: Number of segments= ',nsegb_f,'points=',nvctrb_f
     !write(*,'(1x,a,2(1x,i10))') 'BIG: orbitals have fine   segment, elements',nsegb_f,7*nvctrb_f
  end if
  call make_bounds(nb1,nb2,nb3,logrid_f,ibbyz_f,ibbxz_f,ibbxy_f)

! Create the file grid.xyz to visualize the grid of functions
  if (iproc ==0 .and. output_denspot) then
     call yaml_comment('Writing the file describing the new atomic positions of the effective system')
     !write(*,'(1x,a)')&
     !     'Writing the file describing the new atomic positions of the effective system'
     open(unit=22,file='grid_tail.xyz',status='unknown') !here the output directory can be passed
     write(22,*) nvctrb_c+nvctrb_f,' atomic' 
     write(22,*)'complete simulation grid for the tail correction'
     do iat=1,at%astruct%nat
        write(22,'(a6,2x,3(1x,e12.5),3x)') &
             trim(at%astruct%atomnames(at%astruct%iatype(iat))),txyz(1,iat),txyz(2,iat),txyz(3,iat)
     enddo
     do i3=0,nb3  
        do i2=0,nb2  
           do i1=0,nb1
              if (logrid_c(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   '  g ',real(i1,kind=8)*hgrid(1),real(i2,kind=8)*hgrid(1),real(i3,kind=8)*hgrid(1)
           enddo
        enddo
     end do
     do i3=0,nb3 
        do i2=0,nb2 
           do i1=0,nb1
              if (logrid_f(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   '  G ',real(i1,kind=8)*hgrid(1),real(i2,kind=8)*hgrid(1),real(i3,kind=8)*hgrid(1)
           enddo
        enddo
     enddo
     close(22)
  endif

  call make_all_ib(nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,&
       ibbxy_c,ibbzzx_c,ibbyyzz_c,ibbxy_f,ibbxy_ff,ibbzzx_f,ibbyyzz_f,&
       ibbyz_c,ibbzxx_c,ibbxxyy_c,ibbyz_f,ibbyz_ff,ibbzxx_f,ibbxxyy_f,ibbyyzz_r)

  ! now fill the wavefunction descriptor arrays
  call allocate_wfd(lr%wfd)
  ! coarse grid quantities
  call segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_c,nsegb_c,&
       & lr%wfd%keyglob(1,1),lr%wfd%keyvglob(1))

  ! fine grid quantities
  call segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_f,nsegb_f,&
       & lr%wfd%keyglob(1,lr%wfd%nseg_c+1),lr%wfd%keyvglob(lr%wfd%nseg_c+1))

  call f_free(logrid_c)
  call f_free(logrid_f)


  ! allocations for arrays holding the wavefunction
  !if (iproc == 0) &
  !     write(*,'(1x,a,i0)') 'Allocate words for psib and hpsib ',2*(nvctrb_c+7*nvctrb_f)
  psib = f_malloc(nvctrb_c+7*nvctrb_f,id='psib')
  hpsib = f_malloc(nvctrb_c+7*nvctrb_f,id='hpsib')
  !if (iproc == 0) write(*,*) 'Allocation done'

  ! work arrays applylocpotkin
  psir = f_malloc0((2*nb1+31)*(2*nb2+31)*(2*nb3+31),id='psir')

  if (iproc == 0) then
     call yaml_map('Wavefunction memory occupation in the extended grid (Bytes):',(nvctrb_c+7*nvctrb_f)*8)
     call yaml_comment('Calculating tail corrections, orbitals are processed separately')
     !write(*,'(1x,a,i0)') &
     !     'Wavefunction memory occupation in the extended grid (Bytes): ',&
     !     (nvctrb_c+7*nvctrb_f)*8
     !write(*,'(1x,a,i0,a)') &
     !     'Calculating tail corrections on ',orbs%norbp,' orbitals per processor.'
     !write(*,'(1x,a)',advance='no') &
     !     '     orbitals are processed separately'
  end if

  nw1=max(2*(nb3+1)*(2*nb1+31)*(2*nb2+31),&   ! shrink convention: nw1>nw2
       2*(nb1+1)*(2*nb2+31)*(2*nb3+31))
  nw2=max(4*(nb2+1)*(nb3+1)*(2*nb1+31),&
       4*(nb1+1)*(nb2+1)*(2*nb3+31))

  x_c = f_malloc0((/ 0.to.nb1, 0.to.nb2, 0.to.nb3 /),id='x_c')
  y_c = f_malloc0((/ 0.to.nb1, 0.to.nb2, 0.to.nb3 /),id='y_c')
  x_f = f_malloc0((/ 1.to.7, nbfl1.to.nbfu1, nbfl2.to.nbfu2, nbfl3.to.nbfu3 /),id='x_f')
  y_f = f_malloc0((/ 1.to.7, nbfl1.to.nbfu1, nbfl2.to.nbfu2, nbfl3.to.nbfu3 /),id='y_f')
  w1 = f_malloc(nw1,id='w1')
  w2 = f_malloc(nw2,id='w2')
  x_f1 = f_malloc0((/ nbfl1.to.nbfu1, nbfl2.to.nbfu2, nbfl3.to.nbfu3 /),id='x_f1')
  x_f2 = f_malloc0((/ nbfl1.to.nbfu1, nbfl2.to.nbfu2, nbfl3.to.nbfu3 /),id='x_f2')
  x_f3 = f_malloc0((/ nbfl1.to.nbfu1, nbfl2.to.nbfu2, nbfl3.to.nbfu3 /),id='x_f3')
  !put to zero the arrays for the hamiltonian procedure
  !call to_zero((nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),x_f1)
  !call to_zero((nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),x_f2)
  !call to_zero((nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),x_f3)
  !call to_zero((nb1+1)*(nb2+1)*(nb3+1),x_c)
  !call to_zero(7*(nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),x_f)
  !call to_zero((nb1+1)*(nb2+1)*(nb3+1),y_c)
  !call to_zero(7*(nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),y_f)
  !call to_zero((2*nb1+31)*(2*nb2+31)*(2*nb3+31),psir)
  ekin_sum=0.d0
  epot_sum=0.d0
  eproj_sum=0.d0

  !allocate the fake orbital structure for the application of projectors
  call orbitals_descriptors(0,1,1,1,0,1,1,1, &
       reshape((/0._gp,0._gp,0._gp/),(/3,1/)),(/1._gp /),orbsb,LINEAR_PARTITION_NONE)

  !change positions in gaussian projectors
  call psp_update_positions(nlpsp, lr, Glr, txyz)

  do iorb=1,orbs%norbp

     !build the compressed wavefunction in the enlarged box
     call transform_fortail(n1,n2,n3,nb1,nb2,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,&
        & Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc,Glr%wfd%keyvloc,&
        & Glr%wfd%nseg_f,Glr%wfd%nvctr_f,Glr%wfd%keygloc(1,Glr%wfd%nseg_c+1),Glr%wfd%keyvloc(Glr%wfd%nseg_c+1),  &
        & lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keyglob,lr%wfd%keyvglob,&
        & lr%wfd%nseg_f,lr%wfd%nvctr_f,lr%wfd%keyglob(1,nsegb_c+1),lr%wfd%keyvglob(nsegb_c+1),&
        & nbuf,psi(1,iorb),psi(Glr%wfd%nvctr_c+1,iorb),  & 
        & x_c,x_f,psib(1),psib(nvctrb_c+1))

     !write(*,*) 'transform_fortail finished',iproc,iorb

     if(orbs%spinsgn(iorb+orbs%isorb)>0.0d0) then
        ispin=1
     else
        ispin=2
     end if
     psi_it%ncplx = 1
     psi_it%n_ket = 1
     psi_it%lr => lr
     psi_it%phi_wvl => psib
     psi_it%ispin = ispin
     psi_it%ispsi = 1
     psi_it%nphidim = nvctrb_c+7*nvctrb_f
     psi_it%ob => ob
     ob%orbs => orbsb
     orbsb%npsidim_orbs = nvctrb_c+7*nvctrb_f
       
     npt=2
     tail_adding: do ipt=1,npt

        !for the tail application leave the old-fashioned hamiltonian
        !since it only deals with Free BC and thus no k-points or so

        !calculate gradient
        call applylocpotkinone(nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,nbuf, &
             hgrid,lr%wfd%nseg_c,lr%wfd%nseg_f,lr%wfd%nvctr_c,lr%wfd%nvctr_f, &
             lr%wfd%keyglob,lr%wfd%keyvglob, &
             ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f,y_c,y_f,psir, &
             psib,pot(1,1,1,ispin),hpsib,epot,ekin, &
             x_c,x_f1,x_f2,x_f3,x_f,w1,w2,&
             ibbzzx_c,ibbyyzz_c,ibbxy_ff,ibbzzx_f,ibbyyzz_f,&
             ibbzxx_c,ibbxxyy_c,ibbyz_ff,ibbzxx_f,ibbxxyy_f,nw1,nw2,ibbyyzz_r,1,1)
        !write(*,'(a,3i3,2f12.8)') 'applylocpotkinone finished',iproc,iorb,ipt,epot,ekin

        call DFT_PSP_projectors_iter_new(psp_it, nlpsp)
        loop_proj: do while (DFT_PSP_projectors_iter_next(psp_it, ilr = 1, lr = lr, glr = lr))
           call DFT_PSP_projectors_iter_ensure(psp_it, [0._gp, 0._gp, 0._gp], 0, nwarnings, lr)
           call DFT_PSP_projectors_iter_apply(psp_it, psi_it, at, eproj, hpsi = hpsib, paw = paw)
        end do loop_proj

        !calculate residue for the single orbital
        tt=0.d0
        do i=1,nvctrb_c+7*nvctrb_f
           hpsib(i)=hpsib(i)-orbs%eval(iorb+orbs%isorb)*psib(i)
           tt=tt+hpsib(i)**2
        enddo
        tt=sqrt(tt)

        if (ipt == npt) exit tail_adding

        !calculate tail using the preconditioner as solver for the green function application
        cprecr=-orbs%eval(iorb+orbs%isorb)
        call precong(nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3, &
             lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keyglob,lr%wfd%keyvglob, &
             ncongt,cprecr,hgrid(1),ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f,hpsib)
        !call plot_wf(10,nb1,nb2,nb3,hgrid,nsegb_c,nvctrb_c,keyg,keyv,nsegb_f,nvctrb_f,  & 
        !      txyz(1,1),txyz(2,1),txyz(3,1),psib)

        ! add tail to the bulk wavefunction
        sum_tail=0.d0
        do i=1,nvctrb_c+7*nvctrb_f
           psib(i)=psib(i)-hpsib(i)
           sum_tail=sum_tail+psib(i)**2
        enddo
        sum_tail=sqrt(sum_tail)
        !write(*,'(1x,a,i3,3(1x,1pe13.6),1x,1pe9.2)') &
        !     'BIG: iorb,ekin,epot,eproj,gnrm',iorb,ekin,epot,eproj,tt

        !values of the energies before tail application
        !n(c) ekin1=ekin
        !n(c) epot1=epot
        !n(c) eproj1=eproj
        !write(*,'(1x,a,1x,i0,f18.14)') 'norm orbital + tail',iorb,sum_tail
        !call plot_wf(20,nb1,nb2,nb3,hgrid,nsegb_c,nvctrb_c,keyg,keyv,nsegb_f,nvctrb_f,  & 
        !      txyz(1,1),txyz(2,1),txyz(3,1),psib)

        sum_tail=1.d0/sum_tail
        do i=1,nvctrb_c+7*nvctrb_f
           psib(i)=psib(i)*sum_tail
        enddo

     end do tail_adding

!!!     write(*,'(1x,a,i3,3(1x,1pe13.6),2(1x,1pe9.2))') &
!!!          'BIG: iorb,denergies,gnrm,dnorm',&
!!!          iorb,ekin-ekin1,epot-epot1,eproj-eproj1,tt,sum_tail-1.d0

     if (iproc == 0) then
        !write(*,'(a)',advance='no') &
        !     repeat('.',((iorb+orbs%isorb)*40)/orbs%norbp-((iorb-1)*40)/orbs%norbp)
     end if
     ekin_sum=ekin_sum+ekin*orbs%occup(iorb+orbs%isorb)
     epot_sum=epot_sum+epot*orbs%occup(iorb+orbs%isorb)
     eproj_sum=eproj_sum+eproj*orbs%occup(iorb+orbs%isorb)

  end do

  if (iproc == 0) then
     !write(*,'(1x,a)')'done.'
  end if
  call deallocate_orbs(orbsb)

  call f_free_ptr(txyz)
  call f_free(psir)
  call f_free(psib)
  call f_free(hpsib)

  call deallocate_wfd(lr%wfd)

  call f_free(ibbyz_c)
  call f_free(ibbxz_c)
  call f_free(ibbxy_c)
  call f_free(ibbyz_f)
  call f_free(ibbxz_f)
  call f_free(ibbxy_f)
  call f_free(x_c)
  call f_free(y_c)
  call f_free(w1)
  call f_free(w2)
  call f_free(x_f1)
  call f_free(x_f2)
  call f_free(x_f3)
  call f_free(x_f)
  call f_free(y_f)
  call f_free(ibbzzx_c)
  call f_free(ibbyyzz_c)
  call f_free(ibbxy_ff)
  call f_free(ibbzzx_f)
  call f_free(ibbyyzz_f)
  call f_free(ibbzxx_c)
  call f_free(ibbxxyy_c)
  call f_free(ibbyz_ff)
  call f_free(ibbzxx_f)
  call f_free(ibbxxyy_f)
  call f_free(ibbyyzz_r)

  if (nproc > 1) then
     !if (iproc == 0) then
     !   write(*,'(1x,a,f27.14)')'Tail calculation ended'
     !endif
     wrkallred = f_malloc((/ 3, 2 /),id='wrkallred')
     wrkallred(1,2)=ekin_sum
     wrkallred(2,2)=epot_sum 
     wrkallred(3,2)=eproj_sum 
     call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
          MPI_DOUBLE_PRECISION,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
     ekin_sum=wrkallred(1,1) 
     epot_sum=wrkallred(2,1) 
     eproj_sum=wrkallred(3,1)
     call f_free(wrkallred)
  endif

  call f_release_routine()

END SUBROUTINE CalculateTailCorrection


subroutine transform_fortail(n1,n2,n3,nb1,nb2,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,& 
     mseg_c,mvctr_c,keyg_c,keyv_c,mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     msegb_c,mvctrb_c,keybg_c,keybv_c,msegb_f,mvctrb_f,keybg_f,keybv_f,  & 
     nbuf,psi_c,psi_f,psig_c,psig_f,psib_c,psib_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nb1,nb2,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3
  integer, intent(in) :: mseg_c,mvctr_c,mseg_f,mvctr_f,msegb_c,mvctrb_c,msegb_f,mvctrb_f,nbuf
  integer, intent(in) :: keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  integer, intent(in) :: keybg_c(2,msegb_c),keybv_c(msegb_c),keybg_f(2,msegb_f),keybv_f(msegb_f)
  real(kind=8) :: psi_c(mvctr_c),psi_f(7,mvctr_f)
  real(kind=8) :: psib_c(mvctrb_c),psib_f(7,mvctrb_f)
  real(kind=8) :: psig_c(0:n1+2*nbuf,0:n2+2*nbuf,0:n3+2*nbuf)
  real(kind=8) :: psig_f(7,nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3)
  !Local variables
  integer :: iseg,jj,j0,j1,i0,i1,i2,i3,ii,i

  call f_zero(psig_c)
  call f_zero(psig_f)

  ! coarse part
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_c(i+nbuf,i2+nbuf,i3+nbuf)=psi_c(i-i0+jj)
     enddo
  enddo

  do iseg=1,msegb_c
     jj=keybv_c(iseg)
     j0=keybg_c(1,iseg)
     j1=keybg_c(2,iseg)
     ii=j0-1
     i3=ii/((nb1+1)*(nb2+1))
     ii=ii-i3*(nb1+1)*(nb2+1)
     i2=ii/(nb1+1)
     i0=ii-i2*(nb1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psib_c(i-i0+jj)=psig_c(i,i2,i3)
     enddo
  enddo


  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_f(1,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(1,i-i0+jj)
        psig_f(2,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(2,i-i0+jj)
        psig_f(3,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(3,i-i0+jj)
        psig_f(4,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(4,i-i0+jj)
        psig_f(5,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(5,i-i0+jj)
        psig_f(6,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(6,i-i0+jj)
        psig_f(7,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(7,i-i0+jj)
     enddo
  enddo

  do iseg=1,msegb_f
     jj=keybv_f(iseg)
     j0=keybg_f(1,iseg)
     j1=keybg_f(2,iseg)
     ii=j0-1
     i3=ii/((nb1+1)*(nb2+1))
     ii=ii-i3*(nb1+1)*(nb2+1)
     i2=ii/(nb1+1)
     i0=ii-i2*(nb1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psib_f(1,i-i0+jj)=psig_f(1,i,i2,i3)
        psib_f(2,i-i0+jj)=psig_f(2,i,i2,i3)
        psib_f(3,i-i0+jj)=psig_f(3,i,i2,i3)
        psib_f(4,i-i0+jj)=psig_f(4,i,i2,i3)
        psib_f(5,i-i0+jj)=psig_f(5,i,i2,i3)
        psib_f(6,i-i0+jj)=psig_f(6,i,i2,i3)
        psib_f(7,i-i0+jj)=psig_f(7,i,i2,i3)
     enddo
  enddo

END SUBROUTINE transform_fortail


subroutine transform_fortail_prev(n1,n2,n3,nb1,nb2,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,& 
     mseg_c,mvctr_c,keyg_c,keyv_c,mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     msegb_c,mvctrb_c,keybg_c,keybv_c,msegb_f,mvctrb_f,keybg_f,keybv_f,  & 
     nbuf,psi_c,psi_f,psig_c,psig_fc,psig_f,psib_c,psib_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nb1,nb2,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3
  integer, intent(in) :: mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, intent(in) :: msegb_c,mvctrb_c,msegb_f,mvctrb_f,nbuf
  integer, intent(in) :: keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  integer, intent(in) :: keybg_c(2,msegb_c),keybv_c(msegb_c),keybg_f(2,msegb_f),keybv_f(msegb_f)
  real(kind=8) :: psi_c(mvctr_c),psi_f(7,mvctr_f)
  real(kind=8) :: psib_c(mvctrb_c),psib_f(7,mvctrb_f)
  real(kind=8) :: psig_c(0:n1+2*nbuf,0:n2+2*nbuf,0:n3+2*nbuf)
  real(kind=8) :: psig_fc(0:n1+2*nbuf,0:n2+2*nbuf,0:n3+2*nbuf,3),  & 
        & psig_f(7,nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3)
  !Local variables
  integer :: iseg,j0,jj,j1,i0,i1,i2,i3,ii,i

  call f_zero(psig_c)
  call f_zero(psig_fc)
  call f_zero(psig_f)

  ! coarse part
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_c(i+nbuf,i2+nbuf,i3+nbuf)=psi_c(i-i0+jj)
     enddo
  enddo

  do iseg=1,msegb_c
     jj=keybv_c(iseg)
     j0=keybg_c(1,iseg)
     j1=keybg_c(2,iseg)
     ii=j0-1
     i3=ii/((nb1+1)*(nb2+1))
     ii=ii-i3*(nb1+1)*(nb2+1)
     i2=ii/(nb1+1)
     i0=ii-i2*(nb1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psib_c(i-i0+jj)=psig_c(i,i2,i3)
     enddo
  enddo


  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_f(1,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(1,i-i0+jj)
        psig_fc(i+nbuf,i2+nbuf,i3+nbuf,1)=psig_f(1,i+nbuf,i2+nbuf,i3+nbuf)
        psig_f(2,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(2,i-i0+jj)
        psig_fc(i+nbuf,i2+nbuf,i3+nbuf,2)=psig_f(2,i+nbuf,i2+nbuf,i3+nbuf)
        psig_f(3,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(3,i-i0+jj)
        psig_f(4,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(4,i-i0+jj)
        psig_fc(i+nbuf,i2+nbuf,i3+nbuf,3)=psig_f(4,i+nbuf,i2+nbuf,i3+nbuf)
        psig_f(5,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(5,i-i0+jj)
        psig_f(6,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(6,i-i0+jj)
        psig_f(7,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(7,i-i0+jj)
     enddo
  enddo

  do iseg=1,msegb_f
     jj=keybv_f(iseg)
     j0=keybg_f(1,iseg)
     j1=keybg_f(2,iseg)
     ii=j0-1
     i3=ii/((nb1+1)*(nb2+1))
     ii=ii-i3*(nb1+1)*(nb2+1)
     i2=ii/(nb1+1)
     i0=ii-i2*(nb1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psib_f(1,i-i0+jj)=psig_f(1,i,i2,i3)
        psib_f(2,i-i0+jj)=psig_f(2,i,i2,i3)
        psib_f(3,i-i0+jj)=psig_f(3,i,i2,i3)
        psib_f(4,i-i0+jj)=psig_f(4,i,i2,i3)
        psib_f(5,i-i0+jj)=psig_f(5,i,i2,i3)
        psib_f(6,i-i0+jj)=psig_f(6,i,i2,i3)
        psib_f(7,i-i0+jj)=psig_f(7,i,i2,i3)
     enddo
  enddo

END SUBROUTINE transform_fortail_prev


subroutine applylocpotkinone(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
     hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
     y_c,y_f,psir,  &
     psi,pot,hpsi,epot,ekin,x_c,x_f1,x_f2,x_f3,x_f,w1,w2,&
     ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,nw1,nw2,ibyyzz_r,nspinor,npot)!
  !  Applies the local potential and kinetic energy operator to one wavefunction 
  ! Input: pot,psi
  ! Output: hpsi,epot,ekin
  use module_base
  use module_interfaces, only: apply_potential
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf,nw1,nw2
  integer, intent(in) :: nseg_c,nseg_f,nvctr_c,nvctr_f,nspinor,npot
  real(gp), dimension(3), intent(in) :: hgrid
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  integer, dimension(2,-14:2*n3+16,0:n1), intent(in) :: ibzzx_c
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in) :: ibyyzz_c
  integer, dimension(2,nfl1:nfu1,nfl2:nfu2), intent(in) :: ibxy_ff
  integer, dimension(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(in) :: ibzzx_f
  integer, dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16), intent(in) :: ibyyzz_f
  integer, dimension(2,0:n3,-14:2*n1+16), intent(in) :: ibzxx_c
  integer, dimension(2,-14:2*n1+16,-14:2*n2+16), intent(in) :: ibxxyy_c
  integer, dimension(2,nfl2:nfu2,nfl3:nfu3), intent(in) :: ibyz_ff
  integer, dimension(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16), intent(in) :: ibzxx_f
  integer, dimension(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16), intent(in) :: ibxxyy_f
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in) :: ibyyzz_r
  real(wp), dimension(nvctr_c+7*nvctr_f,nspinor), intent(in) :: psi
  real(wp), dimension((2*n1+31)*(2*n2+31)*(2*n3+31),npot), intent(in) :: pot
  real(wp), dimension(nw1), intent(inout) :: w1
  real(wp), dimension(nw2), intent(inout) :: w2
  real(wp), dimension(0:n1,0:n2,0:n3,nspinor), intent(inout) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,nspinor), intent(inout) :: y_f
  real(wp), dimension((2*n1+31)*(2*n2+31)*(2*n3+31),nspinor), intent(inout) :: psir
  real(wp), dimension(0:n1,0:n2,0:n3,nspinor), intent(inout) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,nspinor), intent(inout) :: x_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,NSPINOR),intent(inout) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3,NSPINOR),intent(inout) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2,NSPINOR),intent(inout) :: x_f3
  real(gp), intent(out) :: epot,ekin
  real(wp), dimension(nvctr_c+7*nvctr_f,nspinor), intent(out) :: hpsi
  !local variables
  integer :: i,idx
  real(gp) :: ekino
  real(wp), dimension(0:3) :: scal

  do i=0,3
     scal(i)=1.0_wp
  enddo

  !call to_zero((2*n1+31)*(2*n2+31)*(2*n3+31)*nspinor,psir)

  do idx=1,nspinor  
     call uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
          nseg_c,nvctr_c,keyg(1,1),keyv(1),  & 
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
          scal,psi(1,IDX),psi(nvctr_c+1,IDX),  &
          x_c(0,0,0,idx),x_f(1,nfl1,nfl2,nfl3,idx),&
          x_f1(nfl1,nfl2,nfl3,idx),x_f2(nfl2,nfl1,nfl3,idx),x_f3(nfl3,nfl1,nfl2,idx))
     
     call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          w1,w2,x_c(0,0,0,idx),x_f(1,nfl1,nfl2,nfl3,idx), & 
          psir(1,IDX),ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
     
  end do

  call apply_potential(n1,n2,n3,1,1,1,nbuf,nspinor,npot,psir,pot,epot,&
       ibyyzz_r) !optional

!!  epot=0.0_gp
!!  if (nspinor==1 .or. nspinor == 2) then
!!     do ispinor=1,nspinor
!!        if (nbuf == 0) then
!!           call realspace(ibyyzz_r,pot,psir(1,ispinor),epots,n1,n2,n3)
!!        else
!!           !this is for the tails. In principle it should work only for 
!!           call realspace_nbuf(ibyyzz_r,pot,psir(1,ispinor),epot,n1,n2,n3,nbuf)
!!        endif
!!           epot=epot+epots
!!        end do
!!  else
!!     call realspaceINPLACE(ibyyzz_r,pot,psir,epot,n1,n2,n3)
!!  end if
  
  ekin=0.0_gp
  do idx=1,nspinor
     call comb_shrink(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          w1,w2,psir(1,IDX),&
          ibxy_c,ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
          y_c(0,0,0,IDX),y_f(1,nfl1,nfl2,nfl3,IDX))!,ibyz_c,ibyz_f)
     
     call ConvolkineticT(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
          hgrid(1),hgrid(2),hgrid(3),ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
          x_c(0,0,0,IDX),x_f(1,nfl1,nfl2,nfl3,IDX),&
          y_c(0,0,0,IDX),y_f(1,nfl1,nfl2,nfl3,IDX),EKINO, &
          x_f1(nfl1,nfl2,nfl3,IDX),x_f2(nfl2,nfl1,nfl3,IDX),x_f3(nfl3,nfl1,nfl2,IDX),111)
     ekin=ekin+ekino
     
     call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
          nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
          scal,y_c(0,0,0,IDX),y_f(1,nfl1,nfl2,nfl3,IDX),hpsi(1,IDX),hpsi(nvctr_c+1,IDX))
  end do
  
END SUBROUTINE applylocpotkinone
