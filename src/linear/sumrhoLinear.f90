
!!>   Here starts the routine for building partial density inside the localisation region
!!!   This routine should be treated as a building-block for the linear scaling code
subroutine local_partial_densityLinear(iproc,nproc,rsflag,nscatterarr,&
     nrhotot,Lzd,hxh,hyh,hzh,nspin,orbs,mapping,psi,rho)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => local_partial_densityLinear
  use module_xc
  use Poisson_Solver
  implicit none
  logical, intent(in) :: rsflag
  integer, intent(in) :: iproc,nproc
  integer,intent(inout):: nrhotot
  integer, intent(in) :: nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(local_zone_descriptors), intent(in) :: Lzd
  type(orbitals_data),intent(in) :: orbs
  integer,dimension(orbs%norb),intent(in):: mapping
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  real(dp),dimension(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot,1),max(nspin,orbs%nspinor)),intent(out):: rho
  !local variables
  character(len=*), parameter :: subname='local_partial_densityLinear'
  integer :: iorb,i_stat,i_all,ii, ind, indSmall, indLarge, orbtot
  integer :: oidx,sidx,nspinn,npsir,ncomplex, i1, i2, i3, ilr, ispin
  integer :: nspincomp,i3s,i3e,ii1,ii2,ii3
  real(gp) :: hfac,spinval
  type(workarr_sumrho) :: w
  real(wp), dimension(:,:), allocatable :: psir
  real(dp), dimension(:),allocatable :: rho_p
  real(8):: dnrm2
  integer, dimension(:,:), allocatable :: Lnscatterarr
  integer :: n3d,n3p,n3pi,i3xcsh,i3tmp,jproc 
  character(len=8) :: filename
  character(len=3) :: numb
 !components of wavefunction in real space which must be considered simultaneously
  !and components of the charge density
  if (orbs%nspinor ==4) then
     npsir=4
     nspinn=4
     ncomplex=0
  else
     npsir=1
     nspinn=nspin
     ncomplex=orbs%nspinor-1
  end if
  nspincomp = 1
  if (nspin > 1) nspincomp = 2

 !allocate and define Lnscatterarr which is just a fake
  allocate(Lnscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,Lnscatterarr,'Lnscatterarr',subname)
  Lnscatterarr(:,3) = 0
  Lnscatterarr(:,4) = 0

  !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
  !otherwise use libXC routine
  call xc_init_rho(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot,1)*max(nspin,orbs%nspinor),rho,nproc)

  ind=1
  orbitalsLoop: do ii=1,orbs%norbp

     iorb = ii + orbs%isorb
     ilr = orbs%inwhichLocreg(iorb)

     Lnscatterarr(:,1) = Lzd%Llr(ilr)%d%n3i 
     Lnscatterarr(:,2) = Lzd%Llr(ilr)%d%n3i 


     call initialize_work_arrays_sumrho(Lzd%Llr(ilr),w)
     allocate(rho_p(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspinn), stat=i_stat) !must redefine the size of rho_p?
     call memocc(i_stat,rho_p,'rho_p',subname)
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,npsir+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
  
     if (Lzd%Llr(ilr)%geocode == 'F') then
        call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*npsir,psir)
     end if
 
     !Need to zero rho_p
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspinn, rho_p)

     !print *,'norbp',orbs%norbp,orbs%norb,orbs%nkpts,orbs%kwgts,orbs%iokpt,orbs%occup
     !hfac=orbs%kwgts(orbs%iokpt(ii))*(orbs%occup(iorb)/(hxh*hyh*hzh))
     hfac=orbs%kwgts(orbs%iokpt(ii))*(orbs%occup(mapping(iorb))/(hxh*hyh*hzh))
     spinval=orbs%spinsgn(iorb)



     if (hfac /= 0.d0) then

        !sum for complex function case, npsir=1 in that case
        do oidx=0,ncomplex

           do sidx=1,npsir
              call daub_to_isf(Lzd%Llr(ilr),w,psi(ind),psir(1,sidx))
              ind=ind+Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f
           end do
           

           select case(Lzd%Llr(ilr)%geocode)
           case('F')
              !write(*,*) 'WARNING: MODIFIED CALLING SEQUENCE OF partial_density_free!!!!'
              call partial_density_free((rsflag .and. .not. Lzd%linear),nproc,Lzd%Llr(ilr)%d%n1i,&
                   Lzd%Llr(ilr)%d%n2i,Lzd%Llr(ilr)%d%n3i,npsir,nspinn,Lzd%Llr(ilr)%d%n3i,&!nrhotot,&
                   hfac,Lnscatterarr,spinval,psir,rho_p,Lzd%Llr(ilr)%bounds%ibyyzz_r)
           case('P')

              call partial_density(rsflag,nproc,Lzd%Llr(ilr)%d%n1i,Lzd%Llr(ilr)%d%n2i,Lzd%Llr(ilr)%d%n3i,&
                   npsir,nspinn,Lzd%Llr(ilr)%d%n3i,&!nrhotot,&
                   hfac,nscatterarr,spinval,psir,rho_p)

           case('S')

              call partial_density(rsflag,nproc,Lzd%Llr(ilr)%d%n1i,Lzd%Llr(ilr)%d%n2i,Lzd%Llr(ilr)%d%n3i,&
                   npsir,nspinn,Lzd%Llr(ilr)%d%n3i,&!nrhotot,&
                   hfac,nscatterarr,spinval,psir,rho_p)

           end select

           ! Copy rho_p to the correct place in rho
           indSmall=0
           do ispin=1,nspinn
               do i3=1,Lzd%Llr(ilr)%d%n3i !min(Lzd%Llr(ilr)%d%n3i,nscatterarr(iproc,1)) 
                   ii3 = i3 + Lzd%Llr(ilr)%nsi3 - 1
                   if(ii3 < 0 .and. Lzd%Glr%geocode /='F') ii3=ii3+Lzd%Glr%d%n3i
                   if(ii3+1 > Lzd%Glr%d%n3i .and. Lzd%Glr%geocode /='F') ii3 = modulo(ii3+1,Lzd%Glr%d%n3i+1)
                   do i2=1,Lzd%Llr(ilr)%d%n2i
                       ii2 = i2 + Lzd%Llr(ilr)%nsi2 - 1
                       if(ii2 < 0 .and. Lzd%Glr%geocode =='P') ii2=ii2+Lzd%Glr%d%n2i
                       if(ii2+1 > Lzd%Glr%d%n2i .and. Lzd%Glr%geocode =='P') ii2 = modulo(ii2+1,Lzd%Glr%d%n2i+1)
                       do i1=1,Lzd%Llr(ilr)%d%n1i
                           ii1=i1 + Lzd%Llr(ilr)%nsi1-1
                           if(ii1<0 .and. Lzd%Glr%geocode /= 'F') ii1=ii1+Lzd%Glr%d%n1i
                           if(ii1+1 > Lzd%Glr%d%n1i.and.Lzd%Glr%geocode/='F') ii1 = modulo(ii1+1,Lzd%Glr%d%n1i+1)
                           ! indSmall is the index in the currect localization region
                           indSmall=indSmall+1
                           ! indLarge is the index in the whole box. 
                           indLarge=ii3*Lzd%Glr%d%n2i*Lzd%Glr%d%n1i +&
                               ii2*Lzd%Glr%d%n1i + ii1 + 1
                           rho(indLarge,ispin)=rho(indLarge,ispin)+rho_p(indSmall)
                       end do
                   end do
               end do
           end do
        end do
     else
        ind=ind+(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*max(ncomplex,1)*npsir
     end if

     i_all=-product(shape(rho_p))*kind(rho_p)
     deallocate(rho_p,stat=i_stat)
     call memocc(i_stat,i_all,'rho_p',subname)
     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)
     call deallocate_work_arrays_sumrho(w)
  end do orbitalsLoop
 
  i_all=-product(shape(Lnscatterarr))*kind(Lnscatterarr)
  deallocate(Lnscatterarr,stat=i_stat)
  call memocc(i_stat,i_all,'Lnscatterarr',subname)
 

END SUBROUTINE local_partial_densityLinear
!
!!!
!!!
subroutine partial_density_linear(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
     hfac,nscatterarr,spinsgn,psir,rho_p,ibyyzz_r) 
  use module_base
  use module_types
  implicit none
  logical, intent(in) :: rsflag
  integer, intent(in) :: nproc,n1i,n2i,n3i,nrhotot,nspinn,npsir
  real(gp), intent(in) :: hfac,spinsgn
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
  real(wp), dimension(n1i,n2i,n3i,npsir), intent(in) :: psir
  real(dp), dimension(n1i,n2i,nrhotot,nspinn), intent(inout) :: rho_p
  integer, dimension(:,:,:),pointer :: ibyyzz_r 
  !local variables
  integer :: i3s,jproc,i3off,n3d,isjmp,i1,i2,i3,i1s,i1e,j3,i3sg
  real(gp) :: hfac2
  real(dp) :: psisq,p1,p2,p3,p4,r1,r2,r3,r4
!  integer :: ncount0,ncount1,ncount_rate,ncount_max
!!!  integer :: ithread,nthread,omp_get_thread_num,omp_get_num_threads
  !sum different slices by taking into account the overlap
  i3sg=0
!$omp parallel default(private) shared(n1i,nproc,rsflag,nspinn,nscatterarr,spinsgn) &
!$omp shared(n2i,npsir,hfac,psir,rho_p,n3i,i3sg,ibyyzz_r)
  i3s=0
!!!   ithread=omp_get_thread_num()
!!!   nthread=omp_get_num_threads()
  hfac2=2.0_gp*hfac

!  call system_clock(ncount0,ncount_rate,ncount_max)

  !case without bounds
  i1s=1
  i1e=n1i
  loop_xc_overlap: do jproc=0,nproc-1
     !case for REDUCE_SCATTER approach, not used for GGA since it enlarges the 
     !communication buffer
     if (rsflag) then
        i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
        n3d=nscatterarr(jproc,1)
        if (n3d==0) exit loop_xc_overlap
     else
        i3off=0
        n3d=n3i
     end if
     !here the condition for the MPI_ALLREDUCE should be entered
     if(spinsgn > 0.0_gp) then
        isjmp=1
     else
        isjmp=2
     end if
     do i3=i3off+1,i3off+n3d
        !this allows the presence of GGA with non-isolated BC. If i3 is between 1 and n3i
        !j3=i3. This is useful only when dealing with rsflags and GGA, so we can comment it out
        !j3=modulo(i3-1,n3i)+1 
        j3=i3
        i3s=i3s+1
!!!    if(mod(i3s,nthread) .eq. ithread) then
     !$omp do
        do i2=1,n2i
              i1s=ibyyzz_r(1,i2-15,j3-15)+1
              i1e=ibyyzz_r(2,i2-15,j3-15)+1
           if (npsir == 1) then
              do i1=i1s,i1e
                 !conversion between the different types
                 psisq=real(psir(i1,i2,j3,1),dp)
                 psisq=psisq*psisq
                 rho_p(i1,i2,i3s,isjmp)=rho_p(i1,i2,i3s,isjmp)+real(hfac,dp)*psisq
              end do
           else !similar loop for npsir=4
              do i1=i1s,i1e
                 !conversion between the different types
                 p1=real(psir(i1,i2,j3,1),dp)
                 p2=real(psir(i1,i2,j3,2),dp)
                 p3=real(psir(i1,i2,j3,3),dp)
                 p4=real(psir(i1,i2,j3,4),dp)

                 !density values
                 r1=p1*p1+p2*p2+p3*p3+p4*p4
                 r2=p1*p3+p2*p4
                 r3=p1*p4-p2*p3
                 r4=p1*p1+p2*p2-p3*p3-p4*p4

                 rho_p(i1,i2,i3s,1)=rho_p(i1,i2,i3s,1)+real(hfac,dp)*r1
                 rho_p(i1,i2,i3s,2)=rho_p(i1,i2,i3s,2)+real(hfac2,dp)*r2
                 rho_p(i1,i2,i3s,3)=rho_p(i1,i2,i3s,3)+real(hfac2,dp)*r3
                 rho_p(i1,i2,i3s,4)=rho_p(i1,i2,i3s,4)+real(hfac,dp)*r4
              end do
           end if
        end do
     !$omp enddo
!!!    end if

!$omp critical
        i3sg=max(i3sg,i3s)
!$omp end critical

     end do
     if (.not. rsflag) exit loop_xc_overlap !the whole range is already done
  end do loop_xc_overlap
!$omp end parallel

  if (i3sg /= nrhotot) then
     write(*,'(1x,a,i0,1x,i0)')'ERROR: problem with rho_p: i3s,nrhotot,',i3sg,nrhotot
     stop
  end if

!  call system_clock(ncount1,ncount_rate,ncount_max)
!  write(*,*) 'TIMING:PDF',real(ncount1-ncount0)/real(ncount_rate)
END SUBROUTINE partial_density_linear


subroutine sumrhoForLocalizedBasis2(iproc,nproc,norb,lzd,input,hx,hy,hz,orbs,&
     comsr,ld_coeff,coeff,nrho,rho,at,nscatterarr)
!
use module_base
use module_types
use libxc_functionals
use module_interfaces, exceptThisOne => sumrhoForLocalizedBasis2
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nrho, norb, ld_coeff
real(gp),intent(in):: hx, hy, hz
type(local_zone_descriptors),intent(in):: lzd
type(input_variables),intent(in):: input
type(orbitals_data),intent(in):: orbs
!type(p2pCommsSumrho),intent(inout):: comsr
type(p2pComms),intent(inout):: comsr
!real(8),dimension(orbs%norb,norb),intent(in):: coeff
real(8),dimension(ld_coeff,norb),intent(in):: coeff
real(8),dimension(nrho),intent(out),target:: rho
type(atoms_data),intent(in):: at
integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh

! Local variables
integer:: iorb, jorb, korb, istat, indLarge, i1, i2, i3, ilr, jlr
integer:: i1s, i1e, i2s, i2e, i3s, i3e, i1d, j1d, i2d, j2d, i3d, j3d, indri, indrj, ldim, iall, istr, istri, istrj
integer:: indi2, indi3, indj2, indj3, indl2, indl3, mpisource, mpidest, iiorb, jjorb
integer:: ierr, jproc, is, ie, nreceives
integer:: nfast, nslow, nsameproc, m, i1d0, j1d0, indri0, indrj0, indLarge0
integer:: azones,bzones,ii,izones,jzones,x,y,z,ishift1,ishift2,ishift3,jshift1,jshift2,jshift3
integer,allocatable :: astart(:,:), aend(:,:), bstart(:,:),bend(:,:)
real(8):: tt, hxh, hyh, hzh, factor, totalCharge, tt0, tt1, tt2, tt3, factorTimesDensKern, t1, t2, time
real(8),dimension(:,:),allocatable:: densKern
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete
character(len=*),parameter:: subname='sumrhoForLocalizedBasis2'


if(iproc==0) write(*,'(1x,a)') 'Calculating charge density...'

!lin%comsr%communComplete=.false.
!lin%comsr%computComplete=.false.


! Allocate the density kernel.
allocate(densKern(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, densKern, 'densKern', subname)

!call mpi_barrier(mpi_comm_world, ierr)
!call cpu_time(t1)
! Calculate the density kernel.
if(iproc==0) write(*,'(3x,a)',advance='no') 'calculating the density kernel... '
call timing(iproc,'sumrho_TMB    ','ON')
!!call dgemm('n', 't', orbs%norb, orbs%norb, norb, 1.d0, coeff(1,1), orbs%norb, &
!!     coeff(1,1), orbs%norb, 0.d0, densKern(1,1), orbs%norb)
call dgemm('n', 't', orbs%norb, orbs%norb, norb, 1.d0, coeff(1,1), ld_coeff, &
     coeff(1,1), ld_coeff, 0.d0, densKern(1,1), orbs%norb)
call timing(iproc,'sumrho_TMB    ','OF')
if(iproc==0) write(*,'(a)') 'done.'
!call mpi_barrier(mpi_comm_world, ierr)
!call cpu_time(t2)
!time=t2-t1
!if(iproc==0) write(*,'(a,es12.4)') 'time for kernel:',time


! Define some constant factors.
hxh=.5d0*hx
hyh=.5d0*hy
hzh=.5d0*hz
if(input%nspin==1) then
    factor=2.d0/(hxh*hyh*hzh)
else
    factor=1.d0/(hxh*hyh*hzh)
end if

! Initialize rho.
if (libxc_functionals_isgga()) then
    call razero(nrho, rho)
else
    ! There is no mpi_allreduce, therefore directly initialize to
    ! 10^-20 and not 10^-20/nproc.
    rho=1.d-20
    !call tenminustwenty(nrho, rho, nproc)
end if
call timing(iproc,'p2pSumrho_wait','ON')


call wait_p2p_communication(iproc, nproc, comsr)



! Now calculate the charge density. Each process calculates only one slice of the total charge density.
! Such a slice has the full extent in the x and y direction, but is limited in the z direction.
! The bounds of the slice are given by nscatterarr. To do so, each process has received all orbitals that
! extend into this slice. The number of these orbitals is given by lin%comsr%noverlaps(iproc).
!call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)

call timing(iproc,'p2pSumrho_wait','OF')

call timing(iproc,'sumrho_TMB    ','ON')

! Bounds of the slice in global coordinates.
is=nscatterarr(iproc,3) 
ie=is+nscatterarr(iproc,1)-1

totalCharge=0.d0
do iorb=1,comsr%noverlaps(iproc)
    iiorb=comsr%overlaps(iorb) !global index of orbital iorb
    ilr=orbs%inwhichlocreg(iiorb) !localization region of orbital iorb
    istri=comsr%comarr(5,iorb,iproc)-1 !starting index of orbital iorb in the receive buffer
    do jorb=1,comsr%noverlaps(iproc)
        jjorb=comsr%overlaps(jorb) !global indes of orbital jorb
        jlr=orbs%inwhichlocreg(jjorb) !localization region of orbital jorb
        istrj=comsr%comarr(5,jorb,iproc)-1 !starting index of orbital jorb in the receive buffer

        azones = 1
        bzones = 1
        !Calculate the number of regions to cut alr and blr
        do ii=1,2
           if(lzd%llr(ilr)%outofzone(ii) > 0) azones = azones * 2
           if(lzd%llr(jlr)%outofzone(ii) > 0) bzones = bzones * 2
        end do
      
        !allocate astart and aend
        allocate(astart(3,azones),stat=istat)
        call memocc(istat,astart,'astart',subname)
        allocate(aend(3,azones),stat=istat)
        call memocc(istat,aend,'aend',subname)
       
        !FRACTURE THE FIRST LOCALIZATION REGION
        call fracture_periodic_zone_ISF(azones,lzd%Glr,lzd%Llr(ilr),lzd%Llr(ilr)%outofzone(:),astart,aend)
       
        !allocate bstart and bend
        allocate(bstart(3,bzones),stat=istat)
        call memocc(istat,bstart,'bstart',subname)
        allocate(bend(3,bzones),stat=istat)
        call memocc(istat,bend,'bend',subname)
       
        !FRACTURE SECOND LOCREG
        call fracture_periodic_zone_ISF(bzones,lzd%Glr,lzd%Llr(jlr),lzd%Llr(jlr)%outofzone(:),bstart,bend)

        do izones=1,azones
           do jzones=1,bzones
              ! Bounds of the overlap of orbital iorb and jorb in global coordinates.
              i1s=max(astart(1,izones),bstart(1,jzones))
              i1e=min(aend(1,izones)-1,bend(1,jzones)-1)
              i2s=max(astart(2,izones),bstart(2,jzones))
              i2e=min(aend(2,izones)-1,bend(2,jzones)-1)
              i3s=max(comsr%startingindex(iorb,1),comsr%startingindex(jorb,1))
              i3e=min(comsr%startingindex(iorb,2),comsr%startingindex(jorb,2))
              call transform_ISFcoordinates(1,i1s,i2s,i3s,lzd%Glr,lzd%Llr(ilr),x,y,z,ishift1, ishift2, ishift3)
              call transform_ISFcoordinates(1,i1s,i2s,i3s,lzd%Glr,lzd%Llr(jlr),x,y,z,jshift1, jshift2, jshift3)
              factorTimesDensKern = factor*densKern(iiorb,jjorb)
              ! Now loop over all points in the box in which the orbitals overlap.
              do i3=i3s,i3e !bounds in z direction
                  i3d=i3 -max(is,-ishift3) !z coordinate of orbital iorb with respect to the overlap box
                  j3d=i3 -max(is,-jshift3) !z coordinate of orbital jorb with respect to the overlap box
                  indi3=i3d*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i !z-part of the index of orbital iorb in the 1-dim receive buffer
                  indj3=j3d*lzd%llr(jlr)%d%n2i*lzd%llr(jlr)%d%n1i !z-part of the index of orbital jorb in the 1-dim receive buffer
                  !indl3=(i3-is)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
                  indl3=(modulo(i3-1,Lzd%Glr%d%n3i)-is+1)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
                  do i2=i2s,i2e !bounds in y direction
                      i2d=i2 + ishift2 !y coordinate of orbital iorb with respect to the overlap box
                      j2d=i2 + jshift2 !y coordinate of orbital jorb with respect to the overlap box
                      indi2=i2d*lzd%llr(ilr)%d%n1i !y-part of the index of orbital iorb in the 1-dim receive buffer
                      indj2=j2d*lzd%llr(jlr)%d%n1i !y-part of the index of orbital jorb in the 1-dim receive buffer
                      !indl2=i2*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
                      indl2=(modulo(i2-1,Lzd%Glr%d%n2i)+1)*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
                      m=mod(i1e-i1s+1,4)
                      if(m/=0) then
                          ! The following five variables hold some intermediate results to speed up the code.
                          i1d0= ishift1 
                          j1d0= jshift1
                          indri0 = indi3 + indi2 + istri + 1
                          indrj0 = indj3 + indj2 + istrj + 1
                          indLarge0 = indl3 + indl2 + 1 
                          do i1=i1s,i1s+m-1
                              i1d=i1d0+i1 !x coordinate of orbital iorb with respect to the overlap box
                              j1d=j1d0+i1 !x coordinate of orbital jorb with respect to the overlap box
                              indri = indri0 + i1d !index of orbital iorb in the 1-dim receive buffer
                              indrj = indrj0 + j1d !index of orbital jorb in the 1-dim receive buffer
                              !indLarge = indLarge0 + i1 !index for which the charge density is beeing calculated
                              indLarge = indLarge0 + modulo(i1-1,Lzd%Glr%d%n1i)+1 !index for which the charge density is beeing calculated
                              tt = factorTimesDensKern*comsr%recvBuf(indri)*comsr%recvBuf(indrj)
                              rho(indLarge) = rho(indLarge) + tt !update the charge density at point indLarge
                              totalCharge = totalCharge + tt !add the contribution to the total charge
                          end do
                      end if
                      ! This is the same again, this time with unrolled loops.
                      if(i1e-i1s+1>4) then
                          i1d0= ishift1 
                          j1d0= jshift1
                          indri0 = indi3 + indi2 + istri + 1
                          indrj0 = indj3 + indj2 + istrj + 1
                          indLarge0 = indl3 + indl2 + 1
                          do i1=i1s+m,i1e,4
                              i1d=i1d0+i1
                              j1d=j1d0+i1
                              indri = indri0 + i1d
                              indrj = indrj0 + j1d
                              !indLarge = indLarge0 + i1
                              tt0 = factorTimesDensKern*comsr%recvBuf(indri  )*comsr%recvBuf(indrj  )
                              tt1 = factorTimesDensKern*comsr%recvBuf(indri+1)*comsr%recvBuf(indrj+1)
                              tt2 = factorTimesDensKern*comsr%recvBuf(indri+2)*comsr%recvBuf(indrj+2)
                              tt3 = factorTimesDensKern*comsr%recvBuf(indri+3)*comsr%recvBuf(indrj+3)
                              indLarge = indLarge0 + modulo(i1-1,Lzd%Glr%d%n1i)+1
                              rho(indLarge  ) = rho(indLarge  ) + tt0
                              indLarge = indLarge0 +modulo(i1,Lzd%Glr%d%n1i)+1
                              rho(indLarge) = rho(indLarge) + tt1
                              !rho(indLarge+1) = rho(indLarge+1) + tt1
                              indLarge = indLarge0 + modulo(i1+1,Lzd%Glr%d%n1i)+1
                              rho(indLarge) = rho(indLarge) + tt2
                              !rho(indLarge+2) = rho(indLarge+2) + tt1
                              indLarge = indLarge0 + modulo(i1+2,Lzd%Glr%d%n1i)+1
                              rho(indLarge) = rho(indLarge) + tt3
                              !rho(indLarge+3) = rho(indLarge+3) + tt1
                              totalCharge = totalCharge + tt0 + tt1 + tt2 + tt3
                          end do
                      end if
                  end do
              end do
          end do !jzones
       end do !izones
       iall=-product(shape(astart))*kind(astart)
       deallocate(astart, stat=istat)
       call memocc(istat, iall, 'astart', subname)
       iall=-product(shape(bstart))*kind(bstart)
       deallocate(bstart, stat=istat)
       call memocc(istat, iall, 'bstart', subname)
       iall=-product(shape(aend))*kind(aend)
       deallocate(aend, stat=istat)
       call memocc(istat, iall, 'aend', subname)
       iall=-product(shape(bend))*kind(bend)
       deallocate(bend, stat=istat)
       call memocc(istat, iall, 'bend', subname)
    end do
end do
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1

call timing(iproc,'sumrho_TMB    ','OF')

call mpiallred(totalCharge, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(3x,a,es20.12)') 'Calculation finished. TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh


iall=-product(shape(densKern))*kind(densKern)
deallocate(densKern, stat=istat)
call memocc(istat, iall, 'densKern', subname)

end subroutine sumrhoForLocalizedBasis2


!> Initializes the parameters needed for the communication of the orbitals
!! when calculating the charge density.
!!
!! input arguments
!!  @param jproc        process to which the orbital shall be sent
!!  @param iorb         orbital that is to be sent
!!  @param istDest      the position on the MPI process to which it should be sent
!!  @param tag          communication tag
!!  @param lin          type containing the parameters for the linear scaling version
!! output arguments
!!  @param commsSumrho  contains the parameters
subroutine setCommunicationInformation2(jproc, iorb, is3ovrlp, n3ovrlp, istDest, tag, nlr, Llr, &
           onWhichAtomAll, orbs, commsSumrho)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: jproc, iorb, is3ovrlp, n3ovrlp, istDest, tag, nlr
type(locreg_descriptors),dimension(nlr),intent(in):: Llr
type(orbitals_data):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
integer,dimension(6),intent(out):: commsSumrho

! Local variables
integer:: mpisource, ist, jorb, jlr

! on which MPI process is the orbital that has to be sent to jproc
mpisource=orbs%onWhichMPI(iorb)
commsSumrho(1)=mpisource

! starting index of the orbital on that MPI process
ist=1
do jorb=orbs%isorb_par(mpisource)+1,iorb-1
    jlr=onWhichAtomAll(jorb)
    !ist=ist+lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f
    ist = ist + Llr(jlr)%d%n1i*Llr(jlr)%d%n2i*Llr(jlr)%d%n3i
end do
jlr=onWhichAtomAll(iorb)
ist = ist + Llr(jlr)%d%n1i*Llr(jlr)%d%n2i*(is3ovrlp-1)
commsSumrho(2)=ist

! amount of data to be sent
jlr=onWhichAtomAll(iorb)
!commsSumrho(3)=lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f
commsSumrho(3)=Llr(jlr)%d%n1i*Llr(jlr)%d%n2i*n3ovrlp

!!! localization region to which this orbital belongs to
!!commsSumrho(4)=onWhichAtomAll(iorb)

! to which MPI process should this orbital be sent
!commsSumrho(5)=jproc
commsSumrho(4)=jproc

! the position on the MPI process to which it should be sent
!commsSumrho(6)=istDest
commsSumrho(5)=istDest

! the tag for this communication
!commsSumrho(7)=tag
commsSumrho(6)=tag

! commsSumrho(8): this entry is used as request for the mpi_isend.

! commsSumrho(9): this entry is used as request for the mpi_irecv.


end subroutine setCommunicationInformation2


