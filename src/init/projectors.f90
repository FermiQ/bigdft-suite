!> @file
!!  Routines to handel projectors
!! @author
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


subroutine localize_projectors(iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,rxyz,radii_cf,&
     logrid,at,orbs,nlpspd)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  logical, dimension(0:n1,0:n2,0:n3), intent(inout) :: logrid
  !Local variables
  !n(c) logical :: cmplxprojs
  integer :: istart,ityp,natyp,iat,mproj,nl1,nu1,nl2,nu2,nl3,nu3,mvctr,mseg,nprojelat,i,l
  integer :: ikpt,nkptsproj,ikptp
  real(gp) :: maxfullvol,totfullvol,totzerovol,zerovol,fullvol,maxrad,maxzerovol,rad
  
  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------------ PSP Projectors Creation'
     write(*,'(1x,a4,4x,a4,2(1x,a))')&
          'Type','Name','Number of atoms','Number of projectors'
  end if
  
  nlpspd%nseg_p(0)=0 
  nlpspd%nvctr_p(0)=0 

  istart=1
  nlpspd%nproj=0
  nlpspd%nprojel=0

  if (iproc ==0) then
     !print the number of projectors to be created
     do ityp=1,at%ntypes
        call numb_proj(ityp,at%ntypes,at%psppar,at%npspcode,mproj)
        natyp=0
        do iat=1,at%nat
           if (at%iatype(iat) == ityp) natyp=natyp+1
        end do
        write(*,'(1x,i4,2x,a6,1x,i15,i21)')&
             ityp,trim(at%atomnames(ityp)),natyp,mproj
     end do
  end if

  do iat=1,at%nat

     call numb_proj(at%iatype(iat),at%ntypes,at%psppar,at%npspcode,mproj)
     if (mproj /= 0) then 

        !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))')&
        !     'projector descriptors for atom with mproj ',iat,mproj
        nlpspd%nproj=nlpspd%nproj+mproj

        ! coarse grid quantities
        call pregion_size(at%geocode,rxyz(1,iat),radii_cf(at%iatype(iat),3),cpmult, &
             hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        nlpspd%nboxp_c(1,1,iat)=nl1
        nlpspd%nboxp_c(1,2,iat)=nl2       
        nlpspd%nboxp_c(1,3,iat)=nl3       

        nlpspd%nboxp_c(2,1,iat)=nu1
        nlpspd%nboxp_c(2,2,iat)=nu2
        nlpspd%nboxp_c(2,3,iat)=nu3

        call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%ntypes,at%iatype(iat),rxyz(1,iat),radii_cf(1,3),cpmult,hx,hy,hz,logrid)
        call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)

        nlpspd%nseg_p(2*iat-1)=nlpspd%nseg_p(2*iat-2) + mseg
        nlpspd%nvctr_p(2*iat-1)=nlpspd%nvctr_p(2*iat-2) + mvctr
        istart=istart+mvctr*mproj

        nprojelat=mvctr*mproj

        !print *,'iat,mvctr',iat,mvctr,mseg,mproj

        ! fine grid quantities
        call pregion_size(at%geocode,rxyz(1,iat),radii_cf(at%iatype(iat),2),fpmult,&
             hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        nlpspd%nboxp_f(1,1,iat)=nl1
        nlpspd%nboxp_f(1,2,iat)=nl2
        nlpspd%nboxp_f(1,3,iat)=nl3

        nlpspd%nboxp_f(2,1,iat)=nu1
        nlpspd%nboxp_f(2,2,iat)=nu2
        nlpspd%nboxp_f(2,3,iat)=nu3

        call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%ntypes,at%iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hx,hy,hz,logrid)
        call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
        !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))') 'mseg,mvctr, fine  projectors ',mseg,mvctr
        nlpspd%nseg_p(2*iat)=nlpspd%nseg_p(2*iat-1) + mseg
        nlpspd%nvctr_p(2*iat)=nlpspd%nvctr_p(2*iat-1) + mvctr

        istart=istart+7*mvctr*mproj
        nprojelat=nprojelat+7*mvctr*mproj
        nlpspd%nprojel=max(nlpspd%nprojel,nprojelat)

        !print *,'iat,nprojelat',iat,nprojelat,mvctr,mseg

     else  !(atom has no nonlocal PSP, e.g. H)
        nlpspd%nseg_p(2*iat-1)=nlpspd%nseg_p(2*iat-2) 
        nlpspd%nvctr_p(2*iat-1)=nlpspd%nvctr_p(2*iat-2) 
        nlpspd%nseg_p(2*iat)=nlpspd%nseg_p(2*iat-1) 
        nlpspd%nvctr_p(2*iat)=nlpspd%nvctr_p(2*iat-1) 

        !! the following is necessary to the creation of preconditioning projectors
        !! coarse grid quantities ( when used preconditioners are applied to all atoms
        !! even H if present )
        call pregion_size(at%geocode,rxyz(1,iat),radii_cf(at%iatype(iat),3),cpmult, &
             hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        nlpspd%nboxp_c(1,1,iat)=nl1
        nlpspd%nboxp_c(1,2,iat)=nl2       
        nlpspd%nboxp_c(1,3,iat)=nl3       

        nlpspd%nboxp_c(2,1,iat)=nu1
        nlpspd%nboxp_c(2,2,iat)=nu2
        nlpspd%nboxp_c(2,3,iat)=nu3

        ! fine grid quantities
        call pregion_size(at%geocode,rxyz(1,iat),radii_cf(at%iatype(iat),2),fpmult,&
             hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        nlpspd%nboxp_f(1,1,iat)=nl1
        nlpspd%nboxp_f(1,2,iat)=nl2
        nlpspd%nboxp_f(1,3,iat)=nl3

        nlpspd%nboxp_f(2,1,iat)=nu1
        nlpspd%nboxp_f(2,2,iat)=nu2
        nlpspd%nboxp_f(2,3,iat)=nu3


     endif
  enddo

  !control the strategy to be applied following the memory limit
  !if the projectors takes too much memory allocate only one atom at the same time
  !control the memory of the projectors expressed in GB
  if (memorylimit /= 0.e0 .and. .not. DistProjApply .and. &
       real(istart-1,kind=4) > memorylimit*134217728.0e0) then
     if (iproc == 0) then
        write(*,'(44x,a)') '------ On-the-fly projectors application'
     end if
     DistProjApply =.true.
  end if

  !calculate the fraction of the projector array used for allocate zero values
  !control the hardest and the softest gaussian
  totzerovol=0.0_gp
  maxfullvol=0.0_gp
  totfullvol=0.0_gp
  do iat=1,at%nat
     ityp=at%iatype(iat)
     maxrad=min(maxval(at%psppar(1:4,0,ityp)),cpmult/15.0_gp*radii_cf(ityp,3))
     zerovol=0.0_gp
     fullvol=0.0_gp
     do l=1,4
        do i=1,3
           if (at%psppar(l,i,ityp) /= 0.0_gp) then
              rad=min(at%psppar(l,0,ityp),cpmult/15.0_gp*radii_cf(ityp,3))
              zerovol=zerovol+(maxrad**3-rad**3)
              fullvol=fullvol+maxrad**3
           end if
        end do
     end do
     if (fullvol >= maxfullvol .and. fullvol > 0.0_gp) then
        maxzerovol=zerovol/fullvol
        maxfullvol=fullvol
     end if
     totzerovol=totzerovol+zerovol
     totfullvol=totfullvol+fullvol
  end do

  !assign the total quantity per atom
  zerovol=0.d0
  if (totfullvol /= 0.0_gp) then
     if (DistProjApply) then
        zerovol=maxzerovol
     else
        zerovol=totzerovol/totfullvol
     end if
  end if

  !here is the point in which the projector strategy should be decided
  !DistProjApply shoud never change after this point

  !number of elements of the projectors
  if (.not. DistProjApply) nlpspd%nprojel=istart-1

  !Compute the multiplying coefficient for nprojel in case of imaginary k points.
  !activate the complex projector if there are kpoints
  !TO BE COMMENTED OUT
  !n(c) cmplxprojs= (orbs%kpts(1,1)**2+orbs%kpts(2,1)**2+orbs%kpts(3,1)**2 >0) .or. orbs%nkpts>1
  nkptsproj=1
  if ((.not.DistProjApply) .and. orbs%norbp > 0) then
     nkptsproj = 0
     do ikptp=1,orbs%nkptsp
        ikpt=orbs%iskpts+ikptp
        if (orbs%kpts(1,ikpt)**2+orbs%kpts(2,ikpt)**2+orbs%kpts(3,ikpt)**2 >0 .and. &
             &  orbs%nspinor > 1) then
           nkptsproj = nkptsproj + 2
        else
           nkptsproj = nkptsproj + 1
        end if
     end do
  else if (DistProjApply) then
     do ikptp=1,orbs%nkptsp
        ikpt=orbs%iskpts+ikptp
        if (orbs%kpts(1,ikpt)**2+orbs%kpts(2,ikpt)**2+orbs%kpts(3,ikpt)**2 >0 .and. &
             &  orbs%nspinor > 1) then
           nkptsproj = max(nkptsproj, 2)
        end if
     end do
  end if
  nlpspd%nprojel=nkptsproj*nlpspd%nprojel

  !print *,'iproc,nkptsproj',iproc,nkptsproj,nlpspd%nprojel,orbs%iskpts,orbs%iskpts+orbs%nkptsp

  if (iproc == 0) then
     if (DistProjApply) then
        write(*,'(44x,a)') '------  On-the-fly projectors application'
     else
        write(*,'(44x,a)') '------'
     end if
     write(*,'(1x,a,i21)') 'Total number of projectors =',nlpspd%nproj
     write(*,'(1x,a,i21)') 'Total number of components =',nlpspd%nprojel
     write(*,'(1x,a,i21)') 'Percent of zero components =',nint(100.0_gp*zerovol)
  end if

END SUBROUTINE localize_projectors


!> Fill the proj array with the PSP projectors or their derivatives, following idir value
subroutine fill_projectors(iproc,lr,hx,hy,hz,at,orbs,rxyz,nlpspd,proj,idir)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,idir
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors),intent(in) :: lr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(out) :: proj
  !Local variables
  !n(c) integer, parameter :: nterm_max=20 !if GTH nterm_max=4
  integer :: istart_c,iat,iproj,nwarnings,ikpt,iskpt,iekpt

  if (iproc.eq.0 .and. nlpspd%nproj /=0 .and. idir==0) write(*,'(1x,a)',advance='no') &
       'Calculating wavelets expansion of projectors...'
  !warnings related to the projectors norm
  nwarnings=0
  !allocate these vectors up to the maximum size we can get
  istart_c=1

  !create projectors for any of the k-point hosted by the processor
  !starting kpoint
  if (orbs%norbp > 0) then
     iskpt=orbs%iokpt(1)
     iekpt=orbs%iokpt(orbs%norbp)
  else
     iskpt=1
     iekpt=1
  end if

  do ikpt=iskpt,iekpt
     iproj=0
     do iat=1,at%nat
        !this routine is defined to uniformise the call for on-the-fly application
        call atom_projector(ikpt,iat,idir,istart_c,iproj,&
             lr,hx,hy,hz,rxyz,at,orbs,nlpspd,proj,nwarnings)
     enddo
     if (iproj /= nlpspd%nproj) stop 'incorrect number of projectors created'
     ! projector part finished
  end do
  if (istart_c-1 /= nlpspd%nprojel) stop 'incorrect once-and-for-all psp generation'
 
  if (iproc == 0 .and. nlpspd%nproj /=0 .and. idir == 0) then
     if (nwarnings == 0) then
        write(*,'(1x,a)')'done.'
     else
        write(*,'(1x,a,i0,a)')'found ',nwarnings,' warnings.'
        write(*,'(1x,a)')'Some projectors may be too rough.'
        write(*,'(1x,a,f6.3)')&
             'Consider the possibility of modifying hgrid and/or the localisation radii.'
     end if
  end if

END SUBROUTINE fill_projectors


subroutine atom_projector(ikpt,iat,idir,istart_c,iproj,&
     lr,hx,hy,hz,rxyz,at,orbs,nlpspd,proj,nwarnings)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iat,idir,ikpt
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors),intent(in) :: lr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  integer, intent(inout) :: istart_c,iproj,nwarnings
  real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
  !Local variables
  integer :: ityp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,jseg_c,l,i,ncplx
  real(gp) :: kx,ky,kz

  !features of the k-point ikpt
  kx=orbs%kpts(1,ikpt)
  ky=orbs%kpts(2,ikpt)
  kz=orbs%kpts(3,ikpt)

  !evaluate the complexity of the k-point
  if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
     ncplx=1
  else
     ncplx=2
  end if

  ityp=at%iatype(iat)
  mbvctr_c=nlpspd%nvctr_p(2*iat-1)-nlpspd%nvctr_p(2*iat-2)
  mbvctr_f=nlpspd%nvctr_p(2*iat  )-nlpspd%nvctr_p(2*iat-1)

  mbseg_c=nlpspd%nseg_p(2*iat-1)-nlpspd%nseg_p(2*iat-2)
  mbseg_f=nlpspd%nseg_p(2*iat  )-nlpspd%nseg_p(2*iat-1)
  jseg_c=nlpspd%nseg_p(2*iat-2)+1
  !decide the loop bounds
  do l=1,4 !generic case, also for HGHs (for GTH it will stop at l=2)
     do i=1,3 !generic case, also for HGHs (for GTH it will stop at i=2)
        if (at%psppar(l,i,ityp) /= 0.0_gp) then
           call projector(at%geocode,at%atomnames(ityp),iat,idir,l,i,&
                at%psppar(l,0,ityp),rxyz(1,iat),lr,&
                hx,hy,hz,kx,ky,kz,ncplx,&
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                nlpspd%keyv_p(jseg_c),nlpspd%keyg_p(1,jseg_c),proj(istart_c),nwarnings)
           iproj=iproj+2*l-1
           istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx
           !print *,'iproc,istart_c,nlpspd%nprojel',istart_c,nlpspd%nprojel,ncplx, kx, ky, kz, ikpt
           if (istart_c > nlpspd%nprojel+1) stop 'istart_c > nprojel+1'
        endif
     enddo
  enddo
END SUBROUTINE atom_projector


subroutine deallocate_proj_descr(nlpspd,subname)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: subname
  type(nonlocal_psp_descriptors), intent(inout) :: nlpspd

  integer :: i_all, i_stat
  i_all=-product(shape(nlpspd%nboxp_c))*kind(nlpspd%nboxp_c)
  deallocate(nlpspd%nboxp_c,stat=i_stat)
  call memocc(i_stat,i_all,'nboxp_c',subname)
  i_all=-product(shape(nlpspd%nboxp_f))*kind(nlpspd%nboxp_f)
  deallocate(nlpspd%nboxp_f,stat=i_stat)
  call memocc(i_stat,i_all,'nboxp_f',subname)
  i_all=-product(shape(nlpspd%keyg_p))*kind(nlpspd%keyg_p)
  deallocate(nlpspd%keyg_p,stat=i_stat)
  call memocc(i_stat,i_all,'keyg_p',subname)
  i_all=-product(shape(nlpspd%keyv_p))*kind(nlpspd%keyv_p)
  deallocate(nlpspd%keyv_p,stat=i_stat)
  call memocc(i_stat,i_all,'keyv_p',subname)
  i_all=-product(shape(nlpspd%nvctr_p))*kind(nlpspd%nvctr_p)
  deallocate(nlpspd%nvctr_p,stat=i_stat)
  call memocc(i_stat,i_all,'nvctr_p',subname)
  i_all=-product(shape(nlpspd%nseg_p))*kind(nlpspd%nseg_p)
  deallocate(nlpspd%nseg_p,stat=i_stat)
  call memocc(i_stat,i_all,'nseg_p',subname)
END SUBROUTINE deallocate_proj_descr


subroutine projector(geocode,atomname,iat,idir,l,i,gau_a,rxyz,lr,&
     hx,hy,hz,kx,ky,kz,ncplx,&
     mbvctr_c,mbvctr_f,mseg_c,mseg_f,keyv_p,keyg_p,proj,nwarnings)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=20), intent(in) :: atomname
  integer, intent(in) :: iat,idir,l,i,mbvctr_c,mbvctr_f,mseg_c,mseg_f,ncplx
  type(locreg_descriptors),intent(in) :: lr
  real(gp), intent(in) :: hx,hy,hz,gau_a,kx,ky,kz
  !integer, dimension(2,3), intent(in) :: nboxp_c,nboxp_f
  integer, dimension(mseg_c+mseg_f), intent(in) :: keyv_p
  integer, dimension(2,mseg_c+mseg_f), intent(in) :: keyg_p

  real(gp), dimension(3), intent(in) :: rxyz
  integer, intent(inout) :: nwarnings
  real(wp), dimension((mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx), intent(out) :: proj
  !Local variables
  integer, parameter :: nterm_max=20 !if GTH nterm_max=4
  integer :: m,iterm
  !integer :: nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f
  integer :: istart_c,nterm
  real(gp) :: fpi,factor,rx,ry,rz
  real(dp) :: scpr
  integer, dimension(3) :: nterm_arr
  integer, dimension(nterm_max) :: lx,ly,lz
  integer, dimension(3,nterm_max,3) :: lxyz_arr
  real(gp), dimension(nterm_max) :: factors
  real(gp), dimension(nterm_max,3) :: fac_arr

  !this value can also be inserted as a parameter
  fpi=(4.0_gp*atan(1.0_gp))**(-.75_gp)

  rx=rxyz(1) 
  ry=rxyz(2) 
  rz=rxyz(3)

  istart_c=1
  !start of the projectors expansion routine
  factor=sqrt(2.0_gp)*fpi/(sqrt(gau_a)**(2*(l-1)+4*i-1))
  do m=1,2*l-1
    
     if (idir==0) then !normal projector calculation case
        call calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,factors)
        
        factors(1:nterm)=factor*factors(1:nterm)
     else !calculation of projector derivative
        call calc_coeff_derproj(l,i,m,nterm_max,gau_a,nterm_arr,lxyz_arr,fac_arr)
        
        nterm=nterm_arr(idir)
        do iterm=1,nterm
           factors(iterm)=factor*fac_arr(iterm,idir)
           lx(iterm)=lxyz_arr(1,iterm,idir)
           ly(iterm)=lxyz_arr(2,iterm,idir)
           lz(iterm)=lxyz_arr(3,iterm,idir)
        end do
     end if
     
     call crtproj(geocode,nterm,lr,hx,hy,hz,kx,ky,kz,ncplx,&
          gau_a,factors,rx,ry,rz,lx,ly,lz,&
          mbvctr_c,mbvctr_f,mseg_c,mseg_f,keyv_p,keyg_p,proj(istart_c))

     ! testing
     if (idir == 0) then
        !here the norm should be done with the complex components
        call wnrm_wrap(ncplx,mbvctr_c,mbvctr_f,proj(istart_c),scpr)
        if (abs(1.d0-scpr) > 1.d-2) then
           if (abs(1.d0-scpr) > 1.d-1) then
              !if (iproc == 0) then
                 write(*,'(1x,a)')'error found!'
                 write(*,'(1x,a,i4,a,a6,a,i1,a,i1,a,f6.3)')&
                      'The norm of the nonlocal PSP for atom n=',iat,&
                      ' (',trim(atomname),&
                      ') labeled by l=',l,' m=',m,' is ',scpr
                 write(*,'(1x,a)')&
                      'while it is supposed to be about 1.0. Control PSP data or reduce grid spacing.'
              !end if
              stop
           else
!!!              write(*,'(1x,a,i4,a,a6,a,i1,a,i1,a,f4.3)')&
!!!                   'The norm of the nonlocal PSP for atom n=',iat,&
!!!                   ' (',trim(atomname),&
!!!                   ') labeled by l=',l,' m=',m,' is ',scpr
              nwarnings=nwarnings+1
           end if
        end if
        !do iterm=1,nterm
        !   if (iproc.eq.0) write(*,'(1x,a,i0,1x,a,1pe10.3,3(1x,i0))') &
        !        'projector: iat,atomname,gau_a,lx,ly,lz ', & 
        !        iat,trim(at%atomnames(at%iatype(iat))),gau_a,lx(iterm),ly(iterm),lz(iterm)
        !enddo
     end if
     !end testing
     istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
  enddo
END SUBROUTINE projector


!>   Determines the number of projectors (valid for GTH and HGH pseudopotentials)
subroutine numb_proj(ityp,ntypes,psppar,npspcode,mproj)
  use module_base
  implicit none
  integer, intent(in) :: ityp,ntypes
  integer, dimension(ntypes), intent(in) :: npspcode
  real(gp), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  integer, intent(out) :: mproj
  !Local variables
  integer :: l,i

  mproj=0
  if (npspcode(ityp) == 2) then !GTH
     do l=1,2 
        do i=1,2 
           if (psppar(l,i,ityp) /= 0.0_gp) mproj=mproj+2*l-1
        enddo
     enddo
  else if (npspcode(ityp) == 3 .or. npspcode(ityp) == 10) then !HGH and HGH-K
     do l=1,4 
        do i=1,3 
           if (psppar(l,i,ityp) /= 0.0_gp) mproj=mproj+2*l-1
        enddo
     enddo
  end if

END SUBROUTINE numb_proj


!>   Returns the compressed form of a Gaussian projector 
!!   @f$ x^lx * y^ly * z^lz * exp (-1/(2*gau_a^2) *((x-rx)^2 + (y-ry)^2 + (z-rz)^2 )) @f$
!!   in the arrays proj_c, proj_f
subroutine crtproj(geocode,nterm,lr, & 
     hx,hy,hz,kx,ky,kz,ncplx,gau_a,fac_arr,rx,ry,rz,lx,ly,lz, & 
     mvctr_c,mvctr_f,mseg_c,mseg_f,keyv_p,keyg_p,proj)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: nterm,mvctr_c,mvctr_f,mseg_c,mseg_f,ncplx
  real(gp), intent(in) :: hx,hy,hz,gau_a,rx,ry,rz,kx,ky,kz
  integer, dimension(nterm), intent(in) :: lx,ly,lz
  real(gp), dimension(nterm), intent(in) :: fac_arr
  integer, dimension(mseg_c+mseg_f), intent(in) :: keyv_p
  integer, dimension(2,mseg_c+mseg_f), intent(in) :: keyg_p
  real(wp), dimension((mvctr_c+7*mvctr_f)*ncplx), intent(out) :: proj
  type(locreg_descriptors), intent(in) :: lr
  !Local variables
  character(len=*), parameter :: subname='crtproj'
  integer, parameter :: nw=65536
  logical :: perx,pery,perz !variables controlling the periodicity in x,y,z
  integer :: iterm,n_gau,ml1,ml2,ml3,mu1,mu2,mu3,i1,i2,i3
  integer :: ns1,ns2,ns3,n1,n2,n3
  integer :: mvctr,i_all,i_stat,j1,i0,j0,jj,ii,i,iseg,ind_f,ind_c
  integer :: counter !test
  real(wp) :: re_cmplx_prod,im_cmplx_prod
  real(gp) :: factor !n(c) err_norm
  real(wp), allocatable, dimension(:,:,:) :: work
  real(wp), allocatable, dimension(:,:,:,:) :: wprojx,wprojy,wprojz
!$  integer :: ithread,nthread,omp_get_thread_num,omp_get_num_threads

  ! rename region boundaries
  ns1 = lr%ns1
  ns2 = lr%ns2
  ns3 = lr%ns3
  n1  = lr%d%n1
  n2  = lr%d%n2
  n3  = lr%d%n3

  allocate(work(0:nw,2,2+ndebug),stat=i_stat)  !always use complex value
  call memocc(i_stat,work,'work',subname)

  allocate(wprojx(ncplx,0:n1,2,nterm+ndebug),stat=i_stat)
  call memocc(i_stat,wprojx,'wprojx',subname)
  allocate(wprojy(ncplx,0:n2,2,nterm+ndebug),stat=i_stat)
  call memocc(i_stat,wprojy,'wprojy',subname)
  allocate(wprojz(ncplx,0:n3,2,nterm+ndebug),stat=i_stat)
  call memocc(i_stat,wprojz,'wprojz',subname)

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  ! make sure that the coefficients returned by CALL GAUSS_TO_DAUB are zero outside [ml:mr] 
  !n(c) err_norm=0.0_gp 
  do iterm=1,nterm

     factor=fac_arr(iterm)
     n_gau=lx(iterm) 
     call gauss_to_daub_k(hx,kx*hx,ncplx,factor,rx,gau_a,n_gau,ns1,n1,ml1,mu1,&
          wprojx(1,0,1,iterm),work,nw,perx) 
     
     n_gau=ly(iterm) 
     call gauss_to_daub_k(hy,ky*hy,ncplx,1.d0,ry,gau_a,n_gau,ns2,n2,ml2,mu2,&
          wprojy(1,0,1,iterm),work,nw,pery) 

     n_gau=lz(iterm) 
     call gauss_to_daub_k(hz,kz*hz,ncplx,1.d0,rz,gau_a,n_gau,ns3,n3,ml3,mu3,&
          wprojz(1,0,1,iterm),work,nw,perz)

  end do
  !the filling of the projector should be different if ncplx==1 or 2
  !split such as to avoid intensive call to if statements
  if (ncplx == 1) then

     !$omp parallel default(private) shared(mseg_c,keyv_p,keyg_p,n3,n2) &
     !$omp shared(n1,proj,wprojx,wprojy,wprojz,mvctr_c) &
     !$omp shared(nterm,mvctr_f,mseg_f)
     !$  ithread=omp_get_thread_num()
     !$  nthread=omp_get_num_threads()

     !$  if(ithread .eq. 0) then
     mvctr=0
     ! coarse part
     do iseg=1,mseg_c
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/((n1+1)*(n2+1))
        ii=ii-i3*(n1+1)*(n2+1)
        i2=ii/(n1+1)
        i0=ii-i2*(n1+1)
        i1=i0+j1-j0
        do i=i0,i1
           ind_c=i-i0+jj
           mvctr=mvctr+1
           proj(ind_c)=&
                wprojx(1,i,1,1)*wprojy(1,i2,1,1)*wprojz(1,i3,1,1)
        enddo
     enddo

     if (mvctr /=  mvctr_c) then
        !$  write(*,'(1x,a,i0,1x,i0)')' ithread,nthread: ',ithread,nthread
        write(*,'(1x,a,i0,1x,i0)')' ERROR (crtproj 1): mvctr /= mvctr_c ',mvctr,mvctr_c
        stop
     end if
     !$  end if

     !$  if(ithread .eq. 1 .or. nthread .eq. 1) then
     ! First term: fine projector components
     ! fine part
     mvctr=0
     do iseg=mseg_c+1,mseg_c+mseg_f
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/((n1+1)*(n2+1))
        ii=ii-i3*(n1+1)*(n2+1)
        i2=ii/(n1+1)
        i0=ii-i2*(n1+1)
        i1=i0+j1-j0
        do i=i0,i1
           ind_f=mvctr_c+7*(i-i0+jj-1)
           mvctr=mvctr+1
           proj(ind_f+1)=wprojx(1,i,2,1)*wprojy(1,i2,1,1)*wprojz(1,i3,1,1)
           proj(ind_f+2)=wprojx(1,i,1,1)*wprojy(1,i2,2,1)*wprojz(1,i3,1,1)
           proj(ind_f+3)=wprojx(1,i,2,1)*wprojy(1,i2,2,1)*wprojz(1,i3,1,1)
           proj(ind_f+4)=wprojx(1,i,1,1)*wprojy(1,i2,1,1)*wprojz(1,i3,2,1)
           proj(ind_f+5)=wprojx(1,i,2,1)*wprojy(1,i2,1,1)*wprojz(1,i3,2,1)
           proj(ind_f+6)=wprojx(1,i,1,1)*wprojy(1,i2,2,1)*wprojz(1,i3,2,1)
           proj(ind_f+7)=wprojx(1,i,2,1)*wprojy(1,i2,2,1)*wprojz(1,i3,2,1)
        enddo
     enddo

     if (mvctr /= mvctr_f) then
        !$ write(*,'(1x,a,i0,1x,i0)')' ithread,nthread: ',ithread,nthread
        write(*,'(1x,a,i0,1x,i0)')' ERROR (crtproj 1): mvctr /= mvctr_f ',mvctr,mvctr_f
        stop 
     end if
     !$  end if

     if (nterm >= 2) then

        !$  if(ithread .eq. 0) then
        ! Other terms: coarse projector components
        ! coarse part
        do iseg=1,mseg_c
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/((n1+1)*(n2+1))
           ii=ii-i3*(n1+1)*(n2+1)
           i2=ii/(n1+1)
           i0=ii-i2*(n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              ind_c=i-i0+jj
              do iterm=2,nterm
                 proj(ind_c)=proj(ind_c)+&
                      wprojx(1,i,1,iterm)*wprojy(1,i2,1,iterm)*wprojz(1,i3,1,iterm)
              end do
           end do
        end do

        !$  end if

        !$  if(ithread .eq. 1 .or. nthread .eq. 1) then
        ! Other terms: fine projector components
        do iseg=mseg_c+1,mseg_c+mseg_f
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/((n1+1)*(n2+1))
           ii=ii-i3*(n1+1)*(n2+1)
           i2=ii/(n1+1)
           i0=ii-i2*(n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              ind_f=mvctr_c+7*(i-i0+jj-1)
              do iterm=2,nterm
                 proj(ind_f+1)=proj(ind_f+1)+&
                      wprojx(1,i,2,iterm)*wprojy(1,i2,1,iterm)*wprojz(1,i3,1,iterm)
                 proj(ind_f+2)=proj(ind_f+2)+&
                      wprojx(1,i,1,iterm)*wprojy(1,i2,2,iterm)*wprojz(1,i3,1,iterm)
                 proj(ind_f+3)=proj(ind_f+3)+&
                      wprojx(1,i,2,iterm)*wprojy(1,i2,2,iterm)*wprojz(1,i3,1,iterm)
                 proj(ind_f+4)=proj(ind_f+4)+&
                      wprojx(1,i,1,iterm)*wprojy(1,i2,1,iterm)*wprojz(1,i3,2,iterm)
                 proj(ind_f+5)=proj(ind_f+5)+&
                      wprojx(1,i,2,iterm)*wprojy(1,i2,1,iterm)*wprojz(1,i3,2,iterm)
                 proj(ind_f+6)=proj(ind_f+6)+&
                      wprojx(1,i,1,iterm)*wprojy(1,i2,2,iterm)*wprojz(1,i3,2,iterm)
                 proj(ind_f+7)=proj(ind_f+7)+&
                      wprojx(1,i,2,iterm)*wprojy(1,i2,2,iterm)*wprojz(1,i3,2,iterm)
                 !! proj_f(1,i-i0+jj)=proj_f(1,i-i0+jj)+&
                 !!      wprojx(i,2,iterm)*wprojy(i2,1,iterm)*wprojz(i3,1,iterm)
                 !! proj_f(2,i-i0+jj)=proj_f(2,i-i0+jj)+&
                 !!      wprojx(i,1,iterm)*wprojy(i2,2,iterm)*wprojz(i3,1,iterm)
                 !! proj_f(3,i-i0+jj)=proj_f(3,i-i0+jj)+&
                 !!      wprojx(i,2,iterm)*wprojy(i2,2,iterm)*wprojz(i3,1,iterm)
                 !! proj_f(4,i-i0+jj)=proj_f(4,i-i0+jj)+&
                 !!      wprojx(i,1,iterm)*wprojy(i2,1,iterm)*wprojz(i3,2,iterm)
                 !! proj_f(5,i-i0+jj)=proj_f(5,i-i0+jj)+&
                 !!      wprojx(i,2,iterm)*wprojy(i2,1,iterm)*wprojz(i3,2,iterm)
                 !! proj_f(6,i-i0+jj)=proj_f(6,i-i0+jj)+&
                 !!      wprojx(i,1,iterm)*wprojy(i2,2,iterm)*wprojz(i3,2,iterm)
                 !! proj_f(7,i-i0+jj)=proj_f(7,i-i0+jj)+&
                 !!      wprojx(i,2,iterm)*wprojy(i2,2,iterm)*wprojz(i3,2,iterm)
              end do
           end do
        end do

        !$  end if
     end if
     !$omp end parallel

  else if (ncplx == 2) then

     !$omp parallel default(private) shared(mseg_c,keyv_p,keyg_p,n3,n2) &
     !$omp shared(n1,proj,wprojx,wprojy,wprojz,mvctr_c) &
     !$omp shared(nterm,mvctr_f,mseg_f)
     !$  ithread=omp_get_thread_num()
     !$  nthread=omp_get_num_threads()

     !part with real and imaginary part
     !modify the openMP statements such as to benefit from parallelisation

     !$  if(ithread .eq. 0) then 
     mvctr=0
     ! coarse part
     do iseg=1,mseg_c
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/((n1+1)*(n2+1))
        ii=ii-i3*(n1+1)*(n2+1)
        i2=ii/(n1+1)
        i0=ii-i2*(n1+1)
        i1=i0+j1-j0
        do i=i0,i1
           ind_c=i-i0+jj
           mvctr=mvctr+1
           proj(ind_c)=&
                re_cmplx_prod(wprojx(1,i,1,1),wprojy(1,i2,1,1),wprojz(1,i3,1,1))
        enddo
     enddo
     if (mvctr /=  mvctr_c) then
        !$ write(*,'(1x,a,i0,1x,i0)')' ithread,nthread: ',ithread,nthread
        write(*,'(1x,a,i0,1x,i0)')' ERROR (crtproj 2): mvctr /= mvctr_c ',mvctr,mvctr_c
        stop
     end if
     !$  end if

     !$  if(ithread .eq. 1 .or. nthread .eq. 1) then
     ! First term: fine projector components
     ! fine part
     mvctr=0
     do iseg=mseg_c+1,mseg_c+mseg_f
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/((n1+1)*(n2+1))
        ii=ii-i3*(n1+1)*(n2+1)
        i2=ii/(n1+1)
        i0=ii-i2*(n1+1)
        i1=i0+j1-j0
        do i=i0,i1
           ind_f=mvctr_c+7*(i-i0+jj-1)
           mvctr=mvctr+1
           proj(ind_f+1)=re_cmplx_prod(wprojx(1,i,2,1),wprojy(1,i2,1,1),wprojz(1,i3,1,1))
           proj(ind_f+2)=re_cmplx_prod(wprojx(1,i,1,1),wprojy(1,i2,2,1),wprojz(1,i3,1,1))
           proj(ind_f+3)=re_cmplx_prod(wprojx(1,i,2,1),wprojy(1,i2,2,1),wprojz(1,i3,1,1))
           proj(ind_f+4)=re_cmplx_prod(wprojx(1,i,1,1),wprojy(1,i2,1,1),wprojz(1,i3,2,1))
           proj(ind_f+5)=re_cmplx_prod(wprojx(1,i,2,1),wprojy(1,i2,1,1),wprojz(1,i3,2,1))
           proj(ind_f+6)=re_cmplx_prod(wprojx(1,i,1,1),wprojy(1,i2,2,1),wprojz(1,i3,2,1))
           proj(ind_f+7)=re_cmplx_prod(wprojx(1,i,2,1),wprojy(1,i2,2,1),wprojz(1,i3,2,1))
        enddo
     enddo
     if (mvctr /= mvctr_f) then
        !$ write(*,'(1x,a,i0,1x,i0)')' ithread,nthread: ',ithread,nthread
        write(*,'(1x,a,i0,1x,i0)')' ERROR (crtproj 2): mvctr /= mvctr_f ',mvctr,mvctr_f
        stop 
     end if
     !$  end if

     if (nterm >= 2) then

        !$  if(ithread .eq. 0) then
        ! Other terms: coarse projector components
        ! coarse part
        do iseg=1,mseg_c
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/((n1+1)*(n2+1))
           ii=ii-i3*(n1+1)*(n2+1)
           i2=ii/(n1+1)
           i0=ii-i2*(n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              ind_c=i-i0+jj
              do iterm=2,nterm
                 proj(ind_c)=proj(ind_c)+re_cmplx_prod(&
                      wprojx(1,i,1,iterm),wprojy(1,i2,1,iterm),wprojz(1,i3,1,iterm))
              end do
           end do
        end do

        !$  end if

        !$  if(ithread .eq. 1 .or. nthread .eq. 1) then
        ! Other terms: fine projector components
        do iseg=mseg_c+1,mseg_c+mseg_f
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/((n1+1)*(n2+1))
           ii=ii-i3*(n1+1)*(n2+1)
           i2=ii/(n1+1)
           i0=ii-i2*(n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              ind_f=mvctr_c+7*(i-i0+jj-1)
              do iterm=2,nterm
                 proj(ind_f+1)=proj(ind_f+1)+re_cmplx_prod(&
                      wprojx(1,i,2,iterm),wprojy(1,i2,1,iterm),wprojz(1,i3,1,iterm))
                 proj(ind_f+2)=proj(ind_f+2)+re_cmplx_prod(&
                      wprojx(1,i,1,iterm),wprojy(1,i2,2,iterm),wprojz(1,i3,1,iterm))
                 proj(ind_f+3)=proj(ind_f+3)+re_cmplx_prod(&
                      wprojx(1,i,2,iterm),wprojy(1,i2,2,iterm),wprojz(1,i3,1,iterm))
                 proj(ind_f+4)=proj(ind_f+4)+re_cmplx_prod(&
                      wprojx(1,i,1,iterm),wprojy(1,i2,1,iterm),wprojz(1,i3,2,iterm))
                 proj(ind_f+5)=proj(ind_f+5)+re_cmplx_prod(&
                      wprojx(1,i,2,iterm),wprojy(1,i2,1,iterm),wprojz(1,i3,2,iterm))
                 proj(ind_f+6)=proj(ind_f+6)+re_cmplx_prod(&
                      wprojx(1,i,1,iterm),wprojy(1,i2,2,iterm),wprojz(1,i3,2,iterm))
                 proj(ind_f+7)=proj(ind_f+7)+re_cmplx_prod(&
                      wprojx(1,i,2,iterm),wprojy(1,i2,2,iterm),wprojz(1,i3,2,iterm))
              end do
           end do
        end do

        !$  end if
     end if

     !now the imaginary part
     !$  if((ithread == 0 .and. nthread <= 2) .or. ithread == 2) then 
     ! coarse part
     mvctr=0
     do iseg=1,mseg_c
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/((n1+1)*(n2+1))
        ii=ii-i3*(n1+1)*(n2+1)
        i2=ii/(n1+1)
        i0=ii-i2*(n1+1)
        i1=i0+j1-j0
        do i=i0,i1
           ind_c=mvctr_c+7*mvctr_f+i-i0+jj
           !write(17,*)ind_c,(mvctr_c+7*mvctr_f)*ncplx
           mvctr=mvctr+1
           proj(ind_c)=&
                im_cmplx_prod(wprojx(1,i,1,1),wprojy(1,i2,1,1),wprojz(1,i3,1,1))
        enddo
     enddo
     if (mvctr /=  mvctr_c) then
        !$ write(*,'(1x,a,i0,1x,i0)')' ithread,nthread: ',ithread,nthread
        write(*,'(1x,a,i0,1x,i0)')' ERROR (crtproj 3): mvctr /= mvctr_c ',mvctr,mvctr_c
        stop
     end if
     !$  end if


     !$  if((ithread .eq. 1 .and. nthread <=3) .or. nthread .eq. 1 .or. ithread == 3) then
     ! First term: fine projector components
     ! fine part
     mvctr=0
     do iseg=mseg_c+1,mseg_c+mseg_f
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/((n1+1)*(n2+1))
        ii=ii-i3*(n1+1)*(n2+1)
        i2=ii/(n1+1)
        i0=ii-i2*(n1+1)
        i1=i0+j1-j0
        do i=i0,i1
           ind_f=mvctr_c+7*mvctr_f+mvctr_c+7*(i-i0+jj-1)
           mvctr=mvctr+1
           proj(ind_f+1)=im_cmplx_prod(wprojx(1,i,2,1),wprojy(1,i2,1,1),wprojz(1,i3,1,1))
           proj(ind_f+2)=im_cmplx_prod(wprojx(1,i,1,1),wprojy(1,i2,2,1),wprojz(1,i3,1,1))
           proj(ind_f+3)=im_cmplx_prod(wprojx(1,i,2,1),wprojy(1,i2,2,1),wprojz(1,i3,1,1))
           proj(ind_f+4)=im_cmplx_prod(wprojx(1,i,1,1),wprojy(1,i2,1,1),wprojz(1,i3,2,1))
           proj(ind_f+5)=im_cmplx_prod(wprojx(1,i,2,1),wprojy(1,i2,1,1),wprojz(1,i3,2,1))
           proj(ind_f+6)=im_cmplx_prod(wprojx(1,i,1,1),wprojy(1,i2,2,1),wprojz(1,i3,2,1))
           proj(ind_f+7)=im_cmplx_prod(wprojx(1,i,2,1),wprojy(1,i2,2,1),wprojz(1,i3,2,1))
        enddo
     enddo
     if (mvctr /= mvctr_f) then
        !$ write(*,'(1x,a,i0,1x,i0)')' ithread,nthread: ',ithread,nthread
        write(*,'(1x,a,i0,1x,i0)')' ERROR (crtproj 3): mvctr /= mvctr_f ',mvctr,mvctr_f
        stop 
     end if
     !$  end if

     if (nterm >= 2) then

        !$  if((ithread == 0 .and. nthread <= 2) .or. ithread == 2) then 
        ! Other terms: coarse projector components
        ! coarse part
        do iseg=1,mseg_c
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/((n1+1)*(n2+1))
           ii=ii-i3*(n1+1)*(n2+1)
           i2=ii/(n1+1)
           i0=ii-i2*(n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              ind_c=mvctr_c+7*mvctr_f+i-i0+jj
              do iterm=2,nterm
                 proj(ind_c)=proj(ind_c)+im_cmplx_prod(&
                      wprojx(1,i,1,iterm),wprojy(1,i2,1,iterm),wprojz(1,i3,1,iterm))
              end do
           end do
        end do
        !$  end if

        !$  if((ithread == 1 .and. nthread <=3) .or. nthread == 1 .or. ithread == 3) then
        ! Other terms: fine projector components
        do iseg=mseg_c+1,mseg_c+mseg_f
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/((n1+1)*(n2+1))
           ii=ii-i3*(n1+1)*(n2+1)
           i2=ii/(n1+1)
           i0=ii-i2*(n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              ind_f=mvctr_c+7*mvctr_f+mvctr_c+7*(i-i0+jj-1)
              do iterm=2,nterm
                 proj(ind_f+1)=proj(ind_f+1)+im_cmplx_prod(&
                      wprojx(1,i,2,iterm),wprojy(1,i2,1,iterm),wprojz(1,i3,1,iterm))
                 proj(ind_f+2)=proj(ind_f+2)+im_cmplx_prod(&
                      wprojx(1,i,1,iterm),wprojy(1,i2,2,iterm),wprojz(1,i3,1,iterm))
                 proj(ind_f+3)=proj(ind_f+3)+im_cmplx_prod(&
                      wprojx(1,i,2,iterm),wprojy(1,i2,2,iterm),wprojz(1,i3,1,iterm))
                 proj(ind_f+4)=proj(ind_f+4)+im_cmplx_prod(&
                      wprojx(1,i,1,iterm),wprojy(1,i2,1,iterm),wprojz(1,i3,2,iterm))
                 proj(ind_f+5)=proj(ind_f+5)+im_cmplx_prod(&
                      wprojx(1,i,2,iterm),wprojy(1,i2,1,iterm),wprojz(1,i3,2,iterm))
                 proj(ind_f+6)=proj(ind_f+6)+im_cmplx_prod(&
                      wprojx(1,i,1,iterm),wprojy(1,i2,2,iterm),wprojz(1,i3,2,iterm))
                 proj(ind_f+7)=proj(ind_f+7)+im_cmplx_prod(&
                      wprojx(1,i,2,iterm),wprojy(1,i2,2,iterm),wprojz(1,i3,2,iterm))
              end do
           end do
        end do
        !$  end if

     end if

     !$omp end parallel
     
  end if

  i_all=-product(shape(wprojx))*kind(wprojx)
  deallocate(wprojx,stat=i_stat)
  call memocc(i_stat,i_all,'wprojx',subname)
  i_all=-product(shape(wprojy))*kind(wprojy)
  deallocate(wprojy,stat=i_stat)
  call memocc(i_stat,i_all,'wprojy',subname)
  i_all=-product(shape(wprojz))*kind(wprojz)
  deallocate(wprojz,stat=i_stat)
  call memocc(i_stat,i_all,'wprojz',subname)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)

END SUBROUTINE crtproj


!> Real part of the complex product
function re_cmplx_prod(a,b,c)
  use module_base
  implicit none
  real(wp), dimension(2), intent(in) :: a,b,c
  real(wp) :: re_cmplx_prod
  
  re_cmplx_prod=a(1)*b(1)*c(1) &
       -a(1)*b(2)*c(2) &
       -a(2)*b(1)*c(2) &
       -a(2)*b(2)*c(1)
  
END FUNCTION re_cmplx_prod


!>   Imaginary part of the complex product
function im_cmplx_prod(a,b,c)
  use module_base
  implicit none
  real(wp), dimension(2), intent(in) :: a,b,c
  real(wp) :: im_cmplx_prod
  
  im_cmplx_prod=-a(2)*b(2)*c(2) &
       +a(2)*b(1)*c(1) &
       +a(1)*b(2)*c(1) &
       +a(1)*b(1)*c(2)
  
END FUNCTION im_cmplx_prod


!> Finds the size of the smallest subbox that contains a localization region made 
!! out of atom centered spheres
subroutine pregion_size(geocode,rxyz,radius,rmult,hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz,rmult,radius
  real(gp), dimension(3), intent(in) :: rxyz
  integer, intent(out) :: nl1,nu1,nl2,nu2,nl3,nu3
  !Local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  !n(c) real(kind=8) :: onem
  real(gp) :: cxmax,cymax,czmax,cxmin,cymin,czmin,rad

  rad=radius*rmult
  cxmax=rxyz(1)+rad ; cxmin=rxyz(1)-rad
  cymax=rxyz(2)+rad ; cymin=rxyz(2)-rad
  czmax=rxyz(3)+rad ; czmin=rxyz(3)-rad
  !n(c) onem=1.d0-eps_mach

  nl1=ceiling(real(cxmin/hx,kind=8) - eps_mach)   
  nl2=ceiling(real(cymin/hy,kind=8) - eps_mach)   
  nl3=ceiling(real(czmin/hz,kind=8) - eps_mach)   
  nu1=floor(real(cxmax/hx,kind=8) + eps_mach)  
  nu2=floor(real(cymax/hy,kind=8) + eps_mach)  
  nu3=floor(real(czmax/hz,kind=8) + eps_mach)  

  !for non-free BC the projectors are not allowed to be also outside the box
  if (geocode == 'F') then
     if (nl1 < 0)   stop 'nl1: projector region outside cell'
     if (nl2 < 0)   stop 'nl2: projector region outside cell'
     if (nl3 < 0)   stop 'nl3: projector region outside cell'
     if (nu1 > n1)   stop 'nu1: projector region outside cell'
     if (nu2 > n2)   stop 'nu2: projector region outside cell'
     if (nu3 > n3)   stop 'nu3: projector region outside cell'
  else if (geocode == 'S') then
     !correct the extremes if they run outside the box
     if (nl1 < 0 .or. nu1 > n1) then
        nl1=0
        nu1=n1
     end if
     if (nl2 < 0)   stop 'nl2: projector region outside cell'
     if (nu2 > n2)   stop 'nu2: projector region outside cell'
     if (nl3 < 0 .or. nu3 > n3) then
        nl3=0
        nu3=n3
     end if
  else if (geocode == 'P') then
     !correct the extremes if they run outside the box
     if (nl1 < 0 .or. nu1 > n1) then
        nl1=0
        nu1=n1
     end if
     if (nl2 < 0 .or. nu2 > n2) then
        nl2=0
        nu2=n2
     end if
     if (nl3 < 0 .or. nu3 > n3) then
        nl3=0
        nu3=n3
     end if
  end if

END SUBROUTINE pregion_size


subroutine calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,fac_arr)
  use module_base
  implicit none
  integer, intent(in) :: l,i,m,nterm_max
  integer, intent(out) :: nterm
  integer, dimension(nterm_max), intent(out) :: lx,ly,lz
  real(gp), dimension(nterm_max), intent(out) :: fac_arr

  if (l.eq.1 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=0.7071067811865475244008444_gp
  else if (l.eq.1 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3651483716701107423046465_gp
     fac_arr(2)=0.3651483716701107423046465_gp
     fac_arr(3)=0.3651483716701107423046465_gp
  else if (l.eq.1 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=2 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=4 ; lz(3)=0
     lx(4)=2 ; ly(4)=0 ; lz(4)=2
     lx(5)=0 ; ly(5)=2 ; lz(5)=2
     lx(6)=0 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.09200874124564722903948358_gp
     fac_arr(2)=0.1840174824912944580789672_gp
     fac_arr(3)=0.09200874124564722903948358_gp
     fac_arr(4)=0.1840174824912944580789672_gp
     fac_arr(5)=0.1840174824912944580789672_gp
     fac_arr(6)=0.09200874124564722903948358_gp
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=1.000000000000000000000000_gp
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.2) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.000000000000000000000000_gp
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.3) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.000000000000000000000000_gp
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3380617018914066310038473_gp
     fac_arr(2)=0.3380617018914066310038473_gp
     fac_arr(3)=0.3380617018914066310038473_gp
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.2) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3380617018914066310038473_gp
     fac_arr(2)=0.3380617018914066310038473_gp
     fac_arr(3)=0.3380617018914066310038473_gp
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.3) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.3380617018914066310038473_gp
     fac_arr(2)=0.3380617018914066310038473_gp
     fac_arr(3)=0.3380617018914066310038473_gp
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     lx(6)=1 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.06795295885835007261827187_gp
     fac_arr(2)=0.1359059177167001452365437_gp
     fac_arr(3)=0.06795295885835007261827187_gp
     fac_arr(4)=0.1359059177167001452365437_gp
     fac_arr(5)=0.1359059177167001452365437_gp
     fac_arr(6)=0.06795295885835007261827187_gp
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.2) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     lx(6)=0 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.06795295885835007261827187_gp
     fac_arr(2)=0.1359059177167001452365437_gp
     fac_arr(3)=0.06795295885835007261827187_gp
     fac_arr(4)=0.1359059177167001452365437_gp
     fac_arr(5)=0.1359059177167001452365437_gp
     fac_arr(6)=0.06795295885835007261827187_gp
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.3) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=2 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=4 ; lz(3)=1
     lx(4)=2 ; ly(4)=0 ; lz(4)=3
     lx(5)=0 ; ly(5)=2 ; lz(5)=3
     lx(6)=0 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.06795295885835007261827187_gp
     fac_arr(2)=0.1359059177167001452365437_gp
     fac_arr(3)=0.06795295885835007261827187_gp
     fac_arr(4)=0.1359059177167001452365437_gp
     fac_arr(5)=0.1359059177167001452365437_gp
     fac_arr(6)=0.06795295885835007261827187_gp
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=1.414213562373095048801689_gp
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.2) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.414213562373095048801689_gp
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.3) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.414213562373095048801689_gp
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.4) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.7071067811865475244008444_gp
     fac_arr(2)=-0.7071067811865475244008444_gp
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.5) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=-0.4082482904638630163662140_gp
     fac_arr(2)=-0.4082482904638630163662140_gp
     fac_arr(3)=0.8164965809277260327324280_gp
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=1
     lx(2)=0 ; ly(2)=3 ; lz(2)=1
     lx(3)=0 ; ly(3)=1 ; lz(3)=3
     fac_arr(1)=0.3563483225498991795794046_gp
     fac_arr(2)=0.3563483225498991795794046_gp
     fac_arr(3)=0.3563483225498991795794046_gp
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.2) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=1
     lx(2)=1 ; ly(2)=2 ; lz(2)=1
     lx(3)=1 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.3563483225498991795794046_gp
     fac_arr(2)=0.3563483225498991795794046_gp
     fac_arr(3)=0.3563483225498991795794046_gp
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.3) then
     nterm=3
     lx(1)=3 ; ly(1)=1 ; lz(1)=0
     lx(2)=1 ; ly(2)=3 ; lz(2)=0
     lx(3)=1 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3563483225498991795794046_gp
     fac_arr(2)=0.3563483225498991795794046_gp
     fac_arr(3)=0.3563483225498991795794046_gp
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.4) then
     nterm=4
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=4 ; lz(2)=0
     lx(3)=2 ; ly(3)=0 ; lz(3)=2
     lx(4)=0 ; ly(4)=2 ; lz(4)=2
     fac_arr(1)=0.1781741612749495897897023_gp
     fac_arr(2)=-0.1781741612749495897897023_gp
     fac_arr(3)=0.1781741612749495897897023_gp
     fac_arr(4)=-0.1781741612749495897897023_gp
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.5) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=2 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=4 ; lz(3)=0
     lx(4)=2 ; ly(4)=0 ; lz(4)=2
     lx(5)=0 ; ly(5)=2 ; lz(5)=2
     lx(6)=0 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=-0.1028688999747279401740630_gp
     fac_arr(2)=-0.2057377999494558803481260_gp
     fac_arr(3)=-0.1028688999747279401740630_gp
     fac_arr(4)=0.1028688999747279401740630_gp
     fac_arr(5)=0.1028688999747279401740630_gp
     fac_arr(6)=0.2057377999494558803481260_gp
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=1
     lx(2)=2 ; ly(2)=3 ; lz(2)=1
     lx(3)=0 ; ly(3)=5 ; lz(3)=1
     lx(4)=2 ; ly(4)=1 ; lz(4)=3
     lx(5)=0 ; ly(5)=3 ; lz(5)=3
     lx(6)=0 ; ly(6)=1 ; lz(6)=5
     fac_arr(1)=0.05959868750235655989526993_gp
     fac_arr(2)=0.1191973750047131197905399_gp
     fac_arr(3)=0.05959868750235655989526993_gp
     fac_arr(4)=0.1191973750047131197905399_gp
     fac_arr(5)=0.1191973750047131197905399_gp
     fac_arr(6)=0.05959868750235655989526993_gp
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.2) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=1
     lx(2)=3 ; ly(2)=2 ; lz(2)=1
     lx(3)=1 ; ly(3)=4 ; lz(3)=1
     lx(4)=3 ; ly(4)=0 ; lz(4)=3
     lx(5)=1 ; ly(5)=2 ; lz(5)=3
     lx(6)=1 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.05959868750235655989526993_gp
     fac_arr(2)=0.1191973750047131197905399_gp
     fac_arr(3)=0.05959868750235655989526993_gp
     fac_arr(4)=0.1191973750047131197905399_gp
     fac_arr(5)=0.1191973750047131197905399_gp
     fac_arr(6)=0.05959868750235655989526993_gp
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.3) then
     nterm=6
     lx(1)=5 ; ly(1)=1 ; lz(1)=0
     lx(2)=3 ; ly(2)=3 ; lz(2)=0
     lx(3)=1 ; ly(3)=5 ; lz(3)=0
     lx(4)=3 ; ly(4)=1 ; lz(4)=2
     lx(5)=1 ; ly(5)=3 ; lz(5)=2
     lx(6)=1 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.05959868750235655989526993_gp
     fac_arr(2)=0.1191973750047131197905399_gp
     fac_arr(3)=0.05959868750235655989526993_gp
     fac_arr(4)=0.1191973750047131197905399_gp
     fac_arr(5)=0.1191973750047131197905399_gp
     fac_arr(6)=0.05959868750235655989526993_gp
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.4) then
     nterm=8
     lx(1)=6 ; ly(1)=0 ; lz(1)=0
     lx(2)=4 ; ly(2)=2 ; lz(2)=0
     lx(3)=2 ; ly(3)=4 ; lz(3)=0
     lx(4)=0 ; ly(4)=6 ; lz(4)=0
     lx(5)=4 ; ly(5)=0 ; lz(5)=2
     lx(6)=0 ; ly(6)=4 ; lz(6)=2
     lx(7)=2 ; ly(7)=0 ; lz(7)=4
     lx(8)=0 ; ly(8)=2 ; lz(8)=4
     fac_arr(1)=0.02979934375117827994763496_gp
     fac_arr(2)=0.02979934375117827994763496_gp
     fac_arr(3)=-0.02979934375117827994763496_gp
     fac_arr(4)=-0.02979934375117827994763496_gp
     fac_arr(5)=0.05959868750235655989526993_gp
     fac_arr(6)=-0.05959868750235655989526993_gp
     fac_arr(7)=0.02979934375117827994763496_gp
     fac_arr(8)=-0.02979934375117827994763496_gp
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.5) then
     nterm=7
     lx(1)=6 ; ly(1)=0 ; lz(1)=0
     lx(2)=4 ; ly(2)=2 ; lz(2)=0
     lx(3)=2 ; ly(3)=4 ; lz(3)=0
     lx(4)=0 ; ly(4)=6 ; lz(4)=0
     lx(5)=2 ; ly(5)=0 ; lz(5)=4
     lx(6)=0 ; ly(6)=2 ; lz(6)=4
     lx(7)=0 ; ly(7)=0 ; lz(7)=6
     fac_arr(1)=-0.01720465913641697233541246_gp
     fac_arr(2)=-0.05161397740925091700623738_gp
     fac_arr(3)=-0.05161397740925091700623738_gp
     fac_arr(4)=-0.01720465913641697233541246_gp
     fac_arr(5)=0.05161397740925091700623738_gp
     fac_arr(6)=0.05161397740925091700623738_gp
     fac_arr(7)=0.03440931827283394467082492_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.1) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3162277660168379331998894_gp
     fac_arr(2)=0.3162277660168379331998894_gp
     fac_arr(3)=-1.264911064067351732799557_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.2) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3162277660168379331998894_gp
     fac_arr(2)=0.3162277660168379331998894_gp
     fac_arr(3)=-1.264911064067351732799557_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.3) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.7745966692414833770358531_gp
     fac_arr(2)=0.7745966692414833770358531_gp
     fac_arr(3)=-0.5163977794943222513572354_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.4) then
     nterm=2
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.4082482904638630163662140_gp
     fac_arr(2)=-1.224744871391589049098642_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.5) then
     nterm=2
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     fac_arr(1)=-1.224744871391589049098642_gp
     fac_arr(2)=0.4082482904638630163662140_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.6) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     fac_arr(1)=1.000000000000000000000000_gp
     fac_arr(2)=-1.000000000000000000000000_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.7) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=2.000000000000000000000000_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.1) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     lx(6)=1 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.06356417261637282102978506_gp
     fac_arr(2)=0.1271283452327456420595701_gp
     fac_arr(3)=0.06356417261637282102978506_gp
     fac_arr(4)=-0.1906925178491184630893552_gp
     fac_arr(5)=-0.1906925178491184630893552_gp
     fac_arr(6)=-0.2542566904654912841191402_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.2) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     lx(6)=0 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.06356417261637282102978506_gp
     fac_arr(2)=0.1271283452327456420595701_gp
     fac_arr(3)=0.06356417261637282102978506_gp
     fac_arr(4)=-0.1906925178491184630893552_gp
     fac_arr(5)=-0.1906925178491184630893552_gp
     fac_arr(6)=-0.2542566904654912841191402_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.3) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=2 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=4 ; lz(3)=1
     lx(4)=2 ; ly(4)=0 ; lz(4)=3
     lx(5)=0 ; ly(5)=2 ; lz(5)=3
     lx(6)=0 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.1556997888323045941832351_gp
     fac_arr(2)=0.3113995776646091883664703_gp
     fac_arr(3)=0.1556997888323045941832351_gp
     fac_arr(4)=0.05189992961076819806107838_gp
     fac_arr(5)=0.05189992961076819806107838_gp
     fac_arr(6)=-0.1037998592215363961221568_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.4) then
     nterm=5
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     fac_arr(1)=0.08206099398622182182282711_gp
     fac_arr(2)=-0.1641219879724436436456542_gp
     fac_arr(3)=-0.2461829819586654654684813_gp
     fac_arr(4)=0.08206099398622182182282711_gp
     fac_arr(5)=-0.2461829819586654654684813_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.5) then
     nterm=5
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     fac_arr(1)=-0.2461829819586654654684813_gp
     fac_arr(2)=-0.1641219879724436436456542_gp
     fac_arr(3)=0.08206099398622182182282711_gp
     fac_arr(4)=-0.2461829819586654654684813_gp
     fac_arr(5)=0.08206099398622182182282711_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.6) then
     nterm=4
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=4 ; lz(2)=1
     lx(3)=2 ; ly(3)=0 ; lz(3)=3
     lx(4)=0 ; ly(4)=2 ; lz(4)=3
     fac_arr(1)=0.2010075630518424150978747_gp
     fac_arr(2)=-0.2010075630518424150978747_gp
     fac_arr(3)=0.2010075630518424150978747_gp
     fac_arr(4)=-0.2010075630518424150978747_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.7) then
     nterm=3
     lx(1)=3 ; ly(1)=1 ; lz(1)=1
     lx(2)=1 ; ly(2)=3 ; lz(2)=1
     lx(3)=1 ; ly(3)=1 ; lz(3)=3
     fac_arr(1)=0.4020151261036848301957494_gp
     fac_arr(2)=0.4020151261036848301957494_gp
     fac_arr(3)=0.4020151261036848301957494_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.1) then
     nterm=10
     lx(1)=7 ; ly(1)=0 ; lz(1)=0
     lx(2)=5 ; ly(2)=2 ; lz(2)=0
     lx(3)=3 ; ly(3)=4 ; lz(3)=0
     lx(4)=1 ; ly(4)=6 ; lz(4)=0
     lx(5)=5 ; ly(5)=0 ; lz(5)=2
     lx(6)=3 ; ly(6)=2 ; lz(6)=2
     lx(7)=1 ; ly(7)=4 ; lz(7)=2
     lx(8)=3 ; ly(8)=0 ; lz(8)=4
     lx(9)=1 ; ly(9)=2 ; lz(9)=4
     lx(10)=1 ; ly(10)=0 ; lz(10)=6
     fac_arr(1)=0.009103849893318918298413687_gp
     fac_arr(2)=0.02731154967995675489524106_gp
     fac_arr(3)=0.02731154967995675489524106_gp
     fac_arr(4)=0.009103849893318918298413687_gp
     fac_arr(5)=-0.01820769978663783659682737_gp
     fac_arr(6)=-0.03641539957327567319365475_gp
     fac_arr(7)=-0.01820769978663783659682737_gp
     fac_arr(8)=-0.06372694925323242808889581_gp
     fac_arr(9)=-0.06372694925323242808889581_gp
     fac_arr(10)=-0.03641539957327567319365475_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.2) then
     nterm=10
     lx(1)=6 ; ly(1)=1 ; lz(1)=0
     lx(2)=4 ; ly(2)=3 ; lz(2)=0
     lx(3)=2 ; ly(3)=5 ; lz(3)=0
     lx(4)=0 ; ly(4)=7 ; lz(4)=0
     lx(5)=4 ; ly(5)=1 ; lz(5)=2
     lx(6)=2 ; ly(6)=3 ; lz(6)=2
     lx(7)=0 ; ly(7)=5 ; lz(7)=2
     lx(8)=2 ; ly(8)=1 ; lz(8)=4
     lx(9)=0 ; ly(9)=3 ; lz(9)=4
     lx(10)=0 ; ly(10)=1 ; lz(10)=6
     fac_arr(1)=0.009103849893318918298413687_gp
     fac_arr(2)=0.02731154967995675489524106_gp
     fac_arr(3)=0.02731154967995675489524106_gp
     fac_arr(4)=0.009103849893318918298413687_gp
     fac_arr(5)=-0.01820769978663783659682737_gp
     fac_arr(6)=-0.03641539957327567319365475_gp
     fac_arr(7)=-0.01820769978663783659682737_gp
     fac_arr(8)=-0.06372694925323242808889581_gp
     fac_arr(9)=-0.06372694925323242808889581_gp
     fac_arr(10)=-0.03641539957327567319365475_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.3) then
     nterm=10
     lx(1)=6 ; ly(1)=0 ; lz(1)=1
     lx(2)=4 ; ly(2)=2 ; lz(2)=1
     lx(3)=2 ; ly(3)=4 ; lz(3)=1
     lx(4)=0 ; ly(4)=6 ; lz(4)=1
     lx(5)=4 ; ly(5)=0 ; lz(5)=3
     lx(6)=2 ; ly(6)=2 ; lz(6)=3
     lx(7)=0 ; ly(7)=4 ; lz(7)=3
     lx(8)=2 ; ly(8)=0 ; lz(8)=5
     lx(9)=0 ; ly(9)=2 ; lz(9)=5
     lx(10)=0 ; ly(10)=0 ; lz(10)=7
     fac_arr(1)=0.02229978693352242055222348_gp
     fac_arr(2)=0.06689936080056726165667044_gp
     fac_arr(3)=0.06689936080056726165667044_gp
     fac_arr(4)=0.02229978693352242055222348_gp
     fac_arr(5)=0.02973304924469656073629797_gp
     fac_arr(6)=0.05946609848939312147259594_gp
     fac_arr(7)=0.02973304924469656073629797_gp
     fac_arr(8)=-0.007433262311174140184074493_gp
     fac_arr(9)=-0.007433262311174140184074493_gp
     fac_arr(10)=-0.01486652462234828036814899_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.4) then
     nterm=9
     lx(1)=7 ; ly(1)=0 ; lz(1)=0
     lx(2)=5 ; ly(2)=2 ; lz(2)=0
     lx(3)=3 ; ly(3)=4 ; lz(3)=0
     lx(4)=1 ; ly(4)=6 ; lz(4)=0
     lx(5)=5 ; ly(5)=0 ; lz(5)=2
     lx(6)=3 ; ly(6)=2 ; lz(6)=2
     lx(7)=1 ; ly(7)=4 ; lz(7)=2
     lx(8)=3 ; ly(8)=0 ; lz(8)=4
     lx(9)=1 ; ly(9)=2 ; lz(9)=4
     fac_arr(1)=0.01175301967439877980816756_gp
     fac_arr(2)=-0.01175301967439877980816756_gp
     fac_arr(3)=-0.05876509837199389904083778_gp
     fac_arr(4)=-0.03525905902319633942450267_gp
     fac_arr(5)=0.02350603934879755961633511_gp
     fac_arr(6)=-0.04701207869759511923267022_gp
     fac_arr(7)=-0.07051811804639267884900533_gp
     fac_arr(8)=0.01175301967439877980816756_gp
     fac_arr(9)=-0.03525905902319633942450267_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.5) then
     nterm=9
     lx(1)=6 ; ly(1)=1 ; lz(1)=0
     lx(2)=4 ; ly(2)=3 ; lz(2)=0
     lx(3)=2 ; ly(3)=5 ; lz(3)=0
     lx(4)=0 ; ly(4)=7 ; lz(4)=0
     lx(5)=4 ; ly(5)=1 ; lz(5)=2
     lx(6)=2 ; ly(6)=3 ; lz(6)=2
     lx(7)=0 ; ly(7)=5 ; lz(7)=2
     lx(8)=2 ; ly(8)=1 ; lz(8)=4
     lx(9)=0 ; ly(9)=3 ; lz(9)=4
     fac_arr(1)=-0.03525905902319633942450267_gp
     fac_arr(2)=-0.05876509837199389904083778_gp
     fac_arr(3)=-0.01175301967439877980816756_gp
     fac_arr(4)=0.01175301967439877980816756_gp
     fac_arr(5)=-0.07051811804639267884900533_gp
     fac_arr(6)=-0.04701207869759511923267022_gp
     fac_arr(7)=0.02350603934879755961633511_gp
     fac_arr(8)=-0.03525905902319633942450267_gp
     fac_arr(9)=0.01175301967439877980816756_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.6) then
     nterm=8
     lx(1)=6 ; ly(1)=0 ; lz(1)=1
     lx(2)=4 ; ly(2)=2 ; lz(2)=1
     lx(3)=2 ; ly(3)=4 ; lz(3)=1
     lx(4)=0 ; ly(4)=6 ; lz(4)=1
     lx(5)=4 ; ly(5)=0 ; lz(5)=3
     lx(6)=0 ; ly(6)=4 ; lz(6)=3
     lx(7)=2 ; ly(7)=0 ; lz(7)=5
     lx(8)=0 ; ly(8)=2 ; lz(8)=5
     fac_arr(1)=0.02878890113916869875409405_gp
     fac_arr(2)=0.02878890113916869875409405_gp
     fac_arr(3)=-0.02878890113916869875409405_gp
     fac_arr(4)=-0.02878890113916869875409405_gp
     fac_arr(5)=0.05757780227833739750818811_gp
     fac_arr(6)=-0.05757780227833739750818811_gp
     fac_arr(7)=0.02878890113916869875409405_gp
     fac_arr(8)=-0.02878890113916869875409405_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.7) then
     nterm=6
     lx(1)=5 ; ly(1)=1 ; lz(1)=1
     lx(2)=3 ; ly(2)=3 ; lz(2)=1
     lx(3)=1 ; ly(3)=5 ; lz(3)=1
     lx(4)=3 ; ly(4)=1 ; lz(4)=3
     lx(5)=1 ; ly(5)=3 ; lz(5)=3
     lx(6)=1 ; ly(6)=1 ; lz(6)=5
     fac_arr(1)=0.05757780227833739750818811_gp
     fac_arr(2)=0.1151556045566747950163762_gp
     fac_arr(3)=0.05757780227833739750818811_gp
     fac_arr(4)=0.1151556045566747950163762_gp
     fac_arr(5)=0.1151556045566747950163762_gp
     fac_arr(6)=0.05757780227833739750818811_gp

  else
     stop 'PSP format error'
  endif
  
END SUBROUTINE calc_coeff_proj


subroutine localize_projectors_paw(iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,rxyz,radii_cf,&
     logrid,at,orbs,PAWD)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs

  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  logical, dimension(0:n1,0:n2,0:n3), intent(inout) :: logrid
  type(PAWproj_data_type) ::PAWD

  !Local variables
  integer :: istart,ityp,natyp,iat,mproj,nl1,nu1,nl2,nu2,nl3,nu3,mvctr,mseg,nprojelat,i,l
  integer :: ikpt,nkptsproj,ikptp
  real(gp) :: maxfullvol,totfullvol,totzerovol,zerovol,fullvol,maxrad,maxzerovol,rad
  integer :: natpaw

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------------ PSP Projectors Creation'
     write(*,'(1x,a4,4x,a4,2(1x,a))')&
          'Type','Name','Number of atoms','Number of paw projectors per atom'
  end if
  
  PAWD%paw_nlpspd%nseg_p(0)=0 
  PAWD%paw_nlpspd%nvctr_p(0)=0 

  istart=1
  PAWD%paw_nlpspd%nproj=0
  PAWD%paw_nlpspd%nprojel=0

  if (iproc ==0) then
     !print the number of projectors to be created
     do ityp=1,at%ntypes
        natyp=0
        mproj=0
        if(  at%paw_NofL(ityp).gt.0) then
           do iat=1,at%nat
              if (at%iatype(iat) == ityp) then
                 if(natyp.eq.0) then
                    call numb_proj_paw(ityp,mproj)                    
                 endif
                 natyp=natyp+1
              endif
           end do
           write(*,'(1x,i4,2x,a6,1x,i15,i21)')&
                ityp,trim(at%atomnames(ityp)),natyp,mproj
        end if
     end do
  end if

  natpaw=0

  do iat=1,at%nat

     if(  at%paw_NofL(at%iatype(iat)).gt.0) then

        call numb_proj_paw(at%iatype(iat),mproj)

        if (mproj /= 0) then 
           natpaw=natpaw+1
           PAWD%paw_nlpspd%nproj=PAWD%paw_nlpspd%nproj+mproj



           ! coarse grid quantities
           call pregion_size(at%geocode,rxyz(1,iat),radii_cf(at%iatype(iat),3),cpmult, &
                hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

           PAWD%paw_nlpspd%nboxp_c(1,1,natpaw)=nl1
           PAWD%paw_nlpspd%nboxp_c(1,2,natpaw)=nl2       
           PAWD%paw_nlpspd%nboxp_c(1,3,natpaw)=nl3       

           PAWD%paw_nlpspd%nboxp_c(2,1,natpaw)=nu1
           PAWD%paw_nlpspd%nboxp_c(2,2,natpaw)=nu2
           PAWD%paw_nlpspd%nboxp_c(2,3,natpaw)=nu3

           call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
                at%ntypes,at%iatype(iat),rxyz(1,iat),radii_cf(1,3),cpmult,hx,hy,hz,logrid)
           call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)

           PAWD%paw_nlpspd%nseg_p (2*natpaw-1)=PAWD%paw_nlpspd%nseg_p (2*natpaw-2) + mseg
           PAWD%paw_nlpspd%nvctr_p(2*natpaw-1)=PAWD%paw_nlpspd%nvctr_p(2*natpaw-2) + mvctr
           istart=istart+mvctr*mproj
 


           nprojelat=mvctr*mproj

           !print *,'iat,mvctr',iat,mvctr,mseg,mproj

           ! fine grid quantities


           call pregion_size(at%geocode,rxyz(1,iat),radii_cf(at%iatype(iat),2),fpmult,&
                hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

           PAWD%paw_nlpspd%nboxp_f(1,1,natpaw)=nl1
           PAWD%paw_nlpspd%nboxp_f(1,2,natpaw)=nl2
           PAWD%paw_nlpspd%nboxp_f(1,3,natpaw)=nl3

           PAWD%paw_nlpspd%nboxp_f(2,1,natpaw)=nu1
           PAWD%paw_nlpspd%nboxp_f(2,2,natpaw)=nu2
           PAWD%paw_nlpspd%nboxp_f(2,3,natpaw)=nu3

           call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
                at%ntypes,at%iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hx,hy,hz,logrid)
           call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
           !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))') 'mseg,mvctr, fine  projectors ',mseg,mvctr
           PAWD%paw_nlpspd%nseg_p (2*natpaw)=PAWD%paw_nlpspd%nseg_p (2*natpaw-1) + mseg
           PAWD%paw_nlpspd%nvctr_p(2*natpaw)=PAWD%paw_nlpspd%nvctr_p(2*natpaw-1) + mvctr



           istart=istart+7*mvctr*mproj
           nprojelat=nprojelat+7*mvctr*mproj



           PAWD%paw_nlpspd%nprojel=max(PAWD%paw_nlpspd%nprojel,nprojelat)

           !print *,'iat,nprojelat',iat,nprojelat,mvctr,mseg

        else  !(atom has no nonlocal PSP, e.g. H)
           PAWD%paw_nlpspd%nseg_p(2*natpaw-1)=PAWD%paw_nlpspd%nseg_p(2*natpaw-2) 
           PAWD%paw_nlpspd%nvctr_p(2*natpaw-1)=PAWD%paw_nlpspd%nvctr_p(2*natpaw-2) 
           PAWD%paw_nlpspd%nseg_p(2*natpaw)=PAWD%paw_nlpspd%nseg_p(2*natpaw-1) 
           PAWD%paw_nlpspd%nvctr_p(2*natpaw)=PAWD%paw_nlpspd%nvctr_p(2*natpaw-1) 

           !! the following is necessary to the creati of preconditioning projectors

           ! coarse grid quantities
           call pregion_size(at%geocode,rxyz(1,iat),radii_cf(at%iatype(iat),3),cpmult, &
                hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

           PAWD%paw_nlpspd%nboxp_c(1,1,natpaw)=nl1
           PAWD%paw_nlpspd%nboxp_c(1,2,natpaw)=nl2       
           PAWD%paw_nlpspd%nboxp_c(1,3,natpaw)=nl3       

           PAWD%paw_nlpspd%nboxp_c(2,1,natpaw)=nu1
           PAWD%paw_nlpspd%nboxp_c(2,2,natpaw)=nu2
           PAWD%paw_nlpspd%nboxp_c(2,3,natpaw)=nu3

           ! fine grid quantities
           call pregion_size(at%geocode,rxyz(1,iat),radii_cf(at%iatype(iat),2),fpmult,&
                hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

           PAWD%paw_nlpspd%nboxp_f(1,1,natpaw)=nl1
           PAWD%paw_nlpspd%nboxp_f(1,2,natpaw)=nl2
           PAWD%paw_nlpspd%nboxp_f(1,3,natpaw)=nl3

           PAWD%paw_nlpspd%nboxp_f(2,1,natpaw)=nu1
           PAWD%paw_nlpspd%nboxp_f(2,2,natpaw)=nu2
           PAWD%paw_nlpspd%nboxp_f(2,3,natpaw)=nu3

        endif
     endif
  enddo
  
  
  !   if (memorylimit /= 0.e0 .and. .not. DistProjApply .and. &
  !        real(istart-1,kind=4) > memorylimit*134217728.0e0) then
  !      if (iproc == 0) then
  !         write(*,'(44x,a)') '------ On-the-fly paw projectors application'
  !      end if
  !      DistProjApply =.true.
  !   end if
  
  !calculate the fraction of the projector array used for allocate zero values
  !control the hardest and the softest gaussian
  totzerovol=0.0_gp
  maxfullvol=0.0_gp
  totfullvol=0.0_gp
  do iat=1,at%nat
     if(  at%paw_NofL(at%iatype(iat)).gt.0) then
        ityp=at%iatype(iat)
        maxrad=min(maxval(at%psppar(1:4,0,ityp)),cpmult/15.0_gp*radii_cf(ityp,3))
        zerovol=0.0_gp
        fullvol=0.0_gp
        do l=1,4
           do i=1,3
              if (at%psppar(l,i,ityp) /= 0.0_gp) then
                 rad=min(at%psppar(l,0,ityp),cpmult/15.0_gp*radii_cf(ityp,3))
                 zerovol=zerovol+(maxrad**3-rad**3)
                 fullvol=fullvol+maxrad**3
              end if
           end do
        end do
        if (fullvol >= maxfullvol .and. fullvol > 0.0_gp) then
           maxzerovol=zerovol/fullvol
           maxfullvol=fullvol
        end if
        totzerovol=totzerovol+zerovol
        totfullvol=totfullvol+fullvol
     endif
  end do

  !assign the total quantity per atom
  zerovol=0.d0
  if (totfullvol /= 0.0_gp) then
     if (PAWD%DistProjApply) then
        zerovol=maxzerovol
     else
        zerovol=totzerovol/totfullvol
     end if
  end if

  !here is the point in which the projector strategy should be decided
  !DistProjApply shoud never change after this point

  !number of elements of the projectors
  if (.not. PAWD%DistProjApply) PAWD%paw_nlpspd%nprojel=istart-1

  nkptsproj=1
  if ((.not.PAWD%DistProjApply) .and. orbs%norbp > 0) then
     nkptsproj = 0
     !the new solution did not work when there is no orbital on the processor
     do ikptp=1,orbs%nkptsp! orbs%iokpt(1), orbs%iokpt(orbs%norbp)
        ikpt=orbs%iskpts+ikptp
!!$         print *, " k points ", orbs%kpts

        if (orbs%kpts(1,ikpt)**2+orbs%kpts(2,ikpt)**2+orbs%kpts(3,ikpt)**2 >0 .and. &
             &  orbs%nspinor > 1) then
           nkptsproj = nkptsproj + 2
        else
           nkptsproj = nkptsproj + 1
        end if
     end do
  else if (PAWD%DistProjApply) then
     !the new solution did not work when there is no orbital on the processor
     do ikptp=1,orbs%nkptsp! orbs%iokpt(1), orbs%iokpt(orbs%norbp)
        ikpt=orbs%iskpts+ikptp
        if (orbs%kpts(1,ikpt)**2+orbs%kpts(2,ikpt)**2+orbs%kpts(3,ikpt)**2 >0 .and. &
             &  orbs%nspinor > 1) then
           nkptsproj = max(nkptsproj, 2)
        end if
     end do
  end if
  !   print *, " nkptsproj EST    ", nkptsproj
  !   print *, " PAWD%paw_nlpspd%nprojel EST  ", PAWD%paw_nlpspd%nprojel

  PAWD%paw_nlpspd%nprojel=nkptsproj*PAWD%paw_nlpspd%nprojel
  if (iproc == 0) then
     if (PAWD%DistProjApply) then
        write(*,'(44x,a)') '------  PAWD: On-the-fly projectors application'
     else
        write(*,'(44x,a)') '------'
     end if
     write(*,'(1x,a,i21)') 'Total number of projectors =',PAWD%paw_nlpspd%nproj
     write(*,'(1x,a,i21)') 'Total number of components =',PAWD%paw_nlpspd%nprojel
     write(*,'(1x,a,i21)') 'Percent of zero components =',nint(100.0_gp*zerovol)
  end if
contains
  
subroutine numb_proj_paw(ityp,mproj)
  integer , intent(in):: ityp
  integer, intent(out):: mproj
  

  integer :: il,jtyp

  mproj=0
  il=0
  do jtyp=1,ityp-1
     il=il+at%paw_NofL(jtyp)
  enddo
  do i =1, at%paw_NofL(ityp)
     il=il+1
     if( at%paw_l(il).ge.0) then
        mproj=mproj+at%paw_nofchannels(il)*(2*at%paw_l(il) +1)
     else
        mproj=mproj+at%paw_nofchannels(il)*(-2*at%paw_l(il) -1)        
     endif
  enddo
end subroutine numb_proj_paw

END subroutine localize_projectors_paw


