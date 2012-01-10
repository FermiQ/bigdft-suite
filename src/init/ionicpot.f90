!> @file
!!  Routines for the ionic energy contribution
!! @author
!!    Copyright (C) 2007-2011 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>    Calculte the ionic contribution to the energy and the forces
subroutine IonicEnergyandForces(iproc,nproc,at,hxh,hyh,hzh,elecfield,&
     rxyz,eion,fion,psoffset,nvacancy,n1,n2,n3,n1i,n2i,n3i,i3s,n3pi,pot_ion,pkernel)
  use module_base
  use module_types
  use Poisson_Solver
  implicit none
  type(atoms_data), intent(in) :: at
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n1i,n2i,n3i,i3s,n3pi,nvacancy
  real(gp), intent(in) :: hxh,hyh,hzh,elecfield(3)
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(*), intent(in) :: pkernel
  real(gp), intent(out) :: eion,psoffset
  real(gp), dimension(3,at%nat), intent(out) :: fion
  real(dp), dimension(*), intent(out) :: pot_ion
  !local variables
  character(len=*), parameter :: subname='IonicEnergyandForces'
  logical :: slowion=.false.
  logical :: perx,pery,perz,gox,goy,goz
  integer :: iat,ii,i_all,i_stat,ityp,jat,jtyp,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
  integer :: isx,iex,isy,iey,isz,iez,i1,i2,i3,j1,j2,j3,ind,ierr
  real(gp) :: ucvol,rloc,twopitothreehalf,pi,atint,shortlength,charge,eself,rx,ry,rz
  real(gp) :: fxion,fyion,fzion,dist,fxerf,fyerf,fzerf,cutoff
  real(gp) :: hxx,hxy,hxz,hyy,hyz,hzz,chgprod,evacancy
  real(gp) :: x,y,z,xp,Vel,prefactor,r2,arg,ehart
  !real(gp) :: Mz,cmassy
  real(gp), dimension(3,3) :: gmet,rmet,rprimd,gprimd
  !other arrays for the ewald treatment
  real(gp), dimension(:,:), allocatable :: fewald,xred

  pi=4.d0*datan(1.d0)
  psoffset=0.0_gp
  if (at%geocode == 'P') then
     !here we insert the calculation of the ewald forces
     allocate(fewald(3,at%nat+ndebug),stat=i_stat)
     call memocc(i_stat,fewald,'fewald',subname)
     allocate(xred(3,at%nat+ndebug),stat=i_stat)
     call memocc(i_stat,xred,'xred',subname)

     !calculate rprimd
     rprimd(:,:)=0.0_gp

     rprimd(1,1)=at%alat1
     rprimd(2,2)=at%alat2
     rprimd(3,3)=at%alat3

     !calculate the metrics and the volume
     call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

     !calculate reduced coordinates
     do iat=1,at%nat
        do ii=1,3
           xred(ii,iat)= gprimd(1,ii)*rxyz(1,iat)+gprimd(2,ii)*rxyz(2,iat)+&
                gprimd(3,ii)*rxyz(3,iat)
        end do

        !in the case of a vacancy calculate the hartree energy by *not* extending 
        !results outside supercell, to compensate the electronic part
        if (iat == nvacancy) then
           rx=rxyz(1,iat) 
           ry=rxyz(2,iat)
           rz=rxyz(3,iat)
           ityp=at%iatype(iat)
           evacancy=0.0_gp
           do jat=1,iat-1
              dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
              jtyp=at%iatype(jat)
              chgprod=real(at%nelpsp(jtyp),gp)*real(at%nelpsp(ityp),gp)
              evacancy=evacancy+chgprod/dist
           enddo
           do jat=iat+1,at%nat
              dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
              jtyp=at%iatype(jat)
              chgprod=real(at%nelpsp(jtyp),gp)*real(at%nelpsp(ityp),gp)
              evacancy=evacancy+chgprod/dist
           end do
           if (iproc == 0) write(*,*)'Ionic energy of the vacancy, to be subtracted:',evacancy
        end if
     end do

     !calculate ewald energy and forces
     call ewald(eion,gmet,fewald,at%nat,at%ntypes,rmet,at%iatype,ucvol,&
          xred,real(at%nelpsp,kind=8))

     !make forces dimensional
     do iat=1,at%nat
        do ii=1,3
           fion(ii,iat)= - (gprimd(ii,1)*fewald(1,iat)+&
                gprimd(ii,2)*fewald(2,iat)+&
                gprimd(ii,3)*fewald(3,iat))
        end do
        !if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
     end do

     i_all=-product(shape(xred))*kind(xred)
     deallocate(xred,stat=i_stat)
     call memocc(i_stat,i_all,'xred',subname)
     i_all=-product(shape(fewald))*kind(fewald)
     deallocate(fewald,stat=i_stat)
     call memocc(i_stat,i_all,'fewald',subname)

     !now calculate the integral of the local psp
     !this is the offset to be applied in the Poisson Solver to have a neutralizing background
     psoffset=0.0_gp
     shortlength=0.0_gp
     charge=0.0_gp
     twopitothreehalf=2.0_gp*pi*sqrt(2.0_gp*pi)
     do iat=1,at%nat
        ityp=at%iatype(iat)
        rloc=at%psppar(0,0,ityp)
        atint=at%psppar(0,1,ityp)+3.0_gp*at%psppar(0,2,ityp)+&
             15.0_gp*at%psppar(0,3,ityp)+105.0_gp*at%psppar(0,4,ityp)
        psoffset=psoffset+rloc**3*atint
        shortlength=shortlength+real(at%nelpsp(ityp),gp)*rloc**2
        charge=charge+real(at%nelpsp(ityp),gp)
     end do
     psoffset=twopitothreehalf*psoffset
     shortlength=shortlength*2.d0*pi

     !print *,'psoffset',psoffset,'pspcore',(psoffset+shortlength)*charge/(at%alat1*at%alat2*at%alat3)
     !if (iproc ==0) print *,'eion',eion,charge/ucvol*(psoffset+shortlength)
     !correct ionic energy taking into account the PSP core correction
     eion=eion+charge/ucvol*(psoffset+shortlength)

!!!     !in the surfaces case, correct the energy term following (J.Chem.Phys. 111(7)-3155, 1999)
!!!     if (at%geocode == 'S') then
!!!        !calculate the Mz dipole component (which in our case corresponds to y direction)
!!!        !first calculate the center of mass
!!!        cmassy=0.0_gp
!!!        do iat=1,at%nat
!!!           cmassy=cmassy+rxyz(2,iat)
!!!        end do
!!!        
!!!        Mz=0.0_gp
!!!        do iat=1,at%nat
!!!           ityp=at%iatype(iat)
!!!           Mz=Mz+real(at%nelpsp(ityp),gp)*(rxyz(2,iat)-cmassy)
!!!        end do
!!!        
!!!        !correct energy and forces in the y direction
!!!        eion=eion+0.5_gp/ucvol*Mz**2
!!!        do iat=1,at%nat
!!!           ityp=at%iatype(iat)
!!!           fion(2,iat)=fion(2,iat)-real(at%nelpsp(ityp),gp)/ucvol*Mz
!!!           if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
!!!        end do
!!!
!!!     end if

  else if (at%geocode == 'F') then

     eion=0.0_gp
     eself=0.0_gp
     do iat=1,at%nat
        ityp=at%iatype(iat)
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)
        !inizialization of the forces
        fxion=0.0_gp
        fyion=0.0_gp
        fzion=0.0_gp
        !initialisation of the hessian
        hxx=0.0_gp
        hxy=0.0_gp
        hxz=0.0_gp
        hyy=0.0_gp
        hyz=0.0_gp
        hzz=0.0_gp

        !    ion-ion interaction
        do jat=1,iat-1
           dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
           jtyp=at%iatype(jat)
           chgprod=real(at%nelpsp(jtyp),gp)*real(at%nelpsp(ityp),gp)
           eion=eion+chgprod/dist
           !forces
           fxion=fxion+chgprod/(dist**3)*(rx-rxyz(1,jat))
           fyion=fyion+chgprod/(dist**3)*(ry-rxyz(2,jat))
           fzion=fzion+chgprod/(dist**3)*(rz-rxyz(3,jat))
           !hessian matrix
           hxx=hxx+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))**2-chgprod/(dist**3)
           hxy=hxy+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(ry-rxyz(2,jat))
           hxz=hxz+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(rz-rxyz(3,jat))
           hyy=hyy+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))**2-chgprod/(dist**3)
           hyz=hyz+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))*(rz-rxyz(3,jat))
           hzz=hzz+3.0_gp*chgprod/(dist**5)*(rz-rxyz(3,jat))**2-chgprod/(dist**3)
        enddo
        do jat=iat+1,at%nat
           dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
           jtyp=at%iatype(jat)
           chgprod=real(at%nelpsp(jtyp),gp)*real(at%nelpsp(ityp),gp)

           !forces
           fxion=fxion+chgprod/(dist**3)*(rx-rxyz(1,jat))
           fyion=fyion+chgprod/(dist**3)*(ry-rxyz(2,jat))
           fzion=fzion+chgprod/(dist**3)*(rz-rxyz(3,jat))
           !hessian matrix
           hxx=hxx+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))**2-chgprod/(dist**3)
           hxy=hxy+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(ry-rxyz(2,jat))
           hxz=hxz+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(rz-rxyz(3,jat))
           hyy=hyy+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))**2-chgprod/(dist**3)
           hyz=hyz+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))*(rz-rxyz(3,jat))
           hzz=hzz+3.0_gp*chgprod/(dist**5)*(rz-rxyz(3,jat))**2-chgprod/(dist**3)
        end do

        fion(1,iat)=fxion
        fion(2,iat)=fyion
        fion(3,iat)=fzion

        !if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
        !energy which comes from the self-interaction of the spread charge
        eself=eself+real(at%nelpsp(ityp)**2,gp)*0.5_gp*sqrt(1.d0/pi)/at%psppar(0,0,ityp)
     end do

     !if (nproc==1 .and. slowion) print *,'eself',eself

  end if

  !for the surfaces BC,
  !activate for the moment only the slow calculation of the ionic energy and forces
  !if (at%geocode == 'S' .or. at%geocode == 'P') slowion=.true.
  if (at%geocode == 'S') slowion=.true.
  !slowion=.true.

  if (slowion) then

     !case of slow ionic calculation
     !conditions for periodicity in the three directions
     perx=(at%geocode /= 'F')
     pery=(at%geocode == 'P')
     perz=(at%geocode /= 'F')

     call ext_buffers(perx,nbl1,nbr1)
     call ext_buffers(pery,nbl2,nbr2)
     call ext_buffers(perz,nbl3,nbr3)

     !the ions corresponds to gaussian charges disposed in the same way as the pseudopotentials

     !first calculate the self-energy and the forces
     !(the latter are zero for a symmetric grid distribution)

     !self energy initialisation
     eself=0.0_gp
     do iat=1,at%nat

        fion(1,iat)=0.0_gp
        fion(2,iat)=0.0_gp
        fion(3,iat)=0.0_gp

        ityp=at%iatype(iat)
        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
        prefactor=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**5)
        cutoff=10.0_gp*rloc

        !calculate the self energy of the isolated bc
        eself=eself+real(at%nelpsp(ityp),gp)**2/rloc

        !in the case of a vacancy calculate the hartree energy by *not* extending 
        !results outside supercell, to compensate the electronic part
        if (iat == nvacancy) then
           evacancy=0.0_gp
           do jat=1,iat-1
              dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
              jtyp=at%iatype(jat)
              chgprod=real(at%nelpsp(jtyp),gp)*real(at%nelpsp(ityp),gp)
              evacancy=evacancy+chgprod/dist
           enddo
           do jat=iat+1,at%nat
              dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
              jtyp=at%iatype(jat)
              chgprod=real(at%nelpsp(jtyp),gp)*real(at%nelpsp(ityp),gp)
              evacancy=evacancy+chgprod/dist
           end do
           print *,'Ionic energy of the vacancy, to be subtracted:',0.5_gp*evacancy
        end if

     enddo

     eself=0.5_gp/sqrt(pi)*eself

     !if (nproc==1) 
     !print *,'iproc,eself',iproc,eself
     call to_zero(n1i*n2i*n3pi,pot_ion(1))

     if (n3pi >0 ) then
        !then calculate the hartree energy and forces of the charge distributions
        !(and save the values for the ionic potential)

        do iat=1,at%nat
           ityp=at%iatype(iat)
           rx=rxyz(1,iat) 
           ry=rxyz(2,iat)
           rz=rxyz(3,iat)

           rloc=at%psppar(0,0,ityp)
           charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
           cutoff=10.0_gp*rloc

           isx=floor((rx-cutoff)/hxh)
           isy=floor((ry-cutoff)/hyh)
           isz=floor((rz-cutoff)/hzh)

           iex=ceiling((rx+cutoff)/hxh)
           iey=ceiling((ry+cutoff)/hyh)
           iez=ceiling((rz+cutoff)/hzh)

           !these nested loops will be used also for the actual ionic forces, to be recalculated
           do i3=isz,iez
              z=real(i3,gp)*hzh-rz
              call ind_positions(perz,i3,n3,j3,goz) 
              j3=j3+nbl3+1
              do i2=isy,iey
                 y=real(i2,gp)*hyh-ry
                 call ind_positions(pery,i2,n2,j2,goy)
                 do i1=isx,iex
                    x=real(i1,gp)*hxh-rx
                    call ind_positions(perx,i1,n1,j1,gox)
                    r2=x**2+y**2+z**2
                    arg=r2/rloc**2
                    xp=exp(-.5_gp*arg)
                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                       pot_ion(ind)=pot_ion(ind)-xp*charge
                    endif
                 enddo
              enddo
           enddo

        enddo

     end if

     !now call the Poisson Solver for the global energy forces
     call H_potential(at%geocode,'D',iproc,nproc,&
          n1i,n2i,n3i,hxh,hyh,hzh,&
          pot_ion,pkernel,pot_ion,ehart,-2.0_gp*psoffset,.false.)

!!$     call PSolver(at%geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
!!$          pot_ion,pkernel,pot_ion,ehart,zero,zero,-2.0_gp*psoffset,.false.,1)

     eion=ehart-eself

     !print *,'ehart,eself',iproc,ehart,eself

     !if (nproc==1) 
     !print *,'iproc,eion',iproc,eion

     do iat=1,at%nat
        ityp=at%iatype(iat)
        !coordinates of the center
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat) 
        rz=rxyz(3,iat)
        !inizialization of the forces

        fxerf=0.0_gp
        fyerf=0.0_gp
        fzerf=0.0_gp

        !local part
        rloc=at%psppar(0,0,ityp)
        prefactor=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**5)
        !maximum extension of the gaussian
        cutoff=10.0_gp*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        !calculate the forces near the atom due to the error function part of the potential
        !calculate forces for all atoms only in the distributed part of the simulation box
        if (n3pi >0 ) then
           do i3=isz,iez
              z=real(i3,gp)*hzh-rz
              call ind_positions(perz,i3,n3,j3,goz) 
              j3=j3+nbl3+1
              do i2=isy,iey
                 y=real(i2,gp)*hyh-ry
                 call ind_positions(pery,i2,n2,j2,goy)
                 do i1=isx,iex
                    x=real(i1,gp)*hxh-rx
                    call ind_positions(perx,i1,n1,j1,gox)
                    r2=x**2+y**2+z**2
                    arg=r2/rloc**2
                    xp=exp(-.5_gp*arg)
                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                       !error function part
                       Vel=pot_ion(ind)
                       fxerf=fxerf+xp*Vel*x
                       fyerf=fyerf+xp*Vel*y
                       fzerf=fzerf+xp*Vel*z
                    endif
                 end do
              end do
           end do
        end if
        !final result of the forces

        fion(1,iat)=fion(1,iat)+(hxh*hyh*hzh*prefactor)*fxerf
        fion(2,iat)=fion(2,iat)+(hxh*hyh*hzh*prefactor)*fyerf
        fion(3,iat)=fion(3,iat)+(hxh*hyh*hzh*prefactor)*fzerf

        !if (nproc==1) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)

!!!        write(10+iat,'(1x,f8.3,i5,(1x,3(1x,1pe12.5)))',advance='no') &
!!!             hxh,iat,(fion(j1,iat),j1=1,3)


     end do

     if (nproc > 1) then
        call mpiallred(fion(1,1),3*at%nat,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     !if (iproc ==0) print *,'eion',eion,psoffset,shortlength

  end if

  ! Add contribution from constant electric field to the forces
  do iat=1,at%nat
     ityp=at%iatype(iat)
     charge=real(at%nelpsp(ityp),gp)
     fion(1:3,iat)=fion(1:3,iat)+(charge*elecfield(1:3))
     !ry=rxyz(2,iat) 
     !eion=eion-(charge*elecfield)*ry
     eion=eion-charge*sum(elecfield(1:3)*rxyz(1:3,iat))
  enddo

  if (iproc == 0) then
     if(all(elecfield(1:3)==0._gp)) then 
           write(*,'(1x,a,1pe22.14)') 'ion-ion interaction energy',eion
     else 
           write(*,'(1x,a,1pe22.14)') 'ion-ion and ion-electric field interaction energy',eion
     endif
     if (nvacancy /= 0) then
        open(unit=22,file='eion_corr.tmp',status='unknown')
        write(22,*)eion-0.5_gp*evacancy,psoffset
        close(unit=22)
     end if
  end if
END SUBROUTINE IonicEnergyandForces


subroutine createIonicPotential(geocode,iproc,nproc,at,rxyz,&
     hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,psoffset,nvacancy,&
     correct_offset)
  use module_base
  use module_types
!  use module_interfaces, except_this_one => createIonicPotential
  use Poisson_Solver
  implicit none
  character(len=1), intent(in) :: geocode
  logical,intent(in) :: correct_offset
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,nvacancy
  real(gp), intent(in) :: hxh,hyh,hzh,psoffset
  type(atoms_data), intent(in) :: at
  real(gp), intent(in) :: elecfield(3)
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(*), intent(in) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  !local variables
  character(len=*), parameter :: subname='createIonicPotential'
  logical :: perx,pery,perz,gox,goy,goz,htoobig=.false.,efwrite,check_potion=.false.
  integer :: iat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp !n(c) nspin
  integer :: ind,i_all,i_stat,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,nloc,iloc
  real(kind=8) :: pi,rholeaked,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot,potxyz,offset
  real(wp) :: maxdiff
  real(gp) :: ehart
  real(dp), dimension(2) :: charges_mpi
  integer, dimension(:,:), allocatable :: ngatherarr
  real(dp), dimension(:), allocatable :: potion_corr
  real(dp), dimension(:), pointer :: pkernel_ref

  call timing(iproc,'CrtLocPot     ','ON')

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '----------------------------------------------------------- Ionic Potential Creation'
  end if

  pi=4.d0*atan(1.d0)
  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  !Creates charge density arising from the ionic PSP cores
  call to_zero(n1i*n2i*n3pi,pot_ion(1))

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  if (n3pi >0 .and. .not. htoobig) then

     do iat=1,at%nat
        ityp=at%iatype(iat)
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)

        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
        cutoff=10.d0*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        do i3=isz,iez
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=isy,iey
              y=real(i2,kind=8)*hyh-ry
              call ind_positions(pery,i2,n2,j2,goy)
              do i1=isx,iex
                 x=real(i1,kind=8)*hxh-rx
                 call ind_positions(perx,i1,n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5d0*arg)
                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                    pot_ion(ind)=pot_ion(ind)-xp*charge
                 else if (.not. goz ) then
                    rholeaked=rholeaked+xp*charge
                 endif
              enddo
           enddo
        enddo

     enddo

  end if

  ! Check
  tt=0.d0
  do j3=1,n3pi
     do i2= -nbl2,2*n2+1+nbr2
        do i1= -nbl1,2*n1+1+nbr1
           ind=i1+1+nbl1+(i2+nbl2)*n1i+(j3-1)*n1i*n2i
           tt=tt+pot_ion(ind)
        enddo
     enddo
  enddo

  tt=tt*hxh*hyh*hzh
  rholeaked=rholeaked*hxh*hyh*hzh

  !print *,'test case input_rho_ion',iproc,i3start,i3end,n3pi,2*n3+16,tt

  if (nproc > 1) then
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked

     call mpiallred(charges_mpi(1),2,MPI_SUM,MPI_COMM_WORLD,ierr)

     tt_tot=charges_mpi(1)
     rholeaked_tot=charges_mpi(2)
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
  end if

  if (iproc == 0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
       'total ionic charge, leaked charge ',tt_tot,rholeaked_tot

  if (.not. htoobig) then
     call timing(iproc,'CrtLocPot     ','OF')
     !here the value of the datacode must be kept fixed
     !n(c) nspin=1

     call H_potential(geocode,'D',iproc,nproc,&
          n1i,n2i,n3i,hxh,hyh,hzh,&
          pot_ion,pkernel,pot_ion,ehart,-psoffset,.false.)

     call timing(iproc,'CrtLocPot     ','ON')
     
     if (check_potion) then
        if (iproc == 0) write(*,'(1x,a)',advance='no') &
             'Check the ionic potential...'
          
        allocate(potion_corr(n1i*n2i*n3pi+ndebug),stat=i_stat)
        call memocc(i_stat,potion_corr,'potion_corr',subname)

        call razero(n1i*n2i*n3pi,potion_corr)

        !calculate pot_ion with an explicit error function to correct in the case of big grid spacings
        !for the moment works only in the isolated BC case
        do i3=1,n3pi
           z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,n2i
              y=real(i2-nbl2-1,gp)*hyh
              do i1=1,n1i
                 x=real(i1-nbl1-1,gp)*hxh
                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                 !if (i1==49 .and. i2==46 .and. i3==44) then
                    call sum_erfcr(at%nat,at%ntypes,x,y,z,at%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
                 !   stop
                 !end if
                 potion_corr(ind)=potion_corr(ind)+potxyz
                 !write(18,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind)
              end do
           end do
        end do

        !then calculate the maximum difference in the sup norm
        maxdiff=0.0_wp
        do i3=1,n3pi
           do i2=1,n2i
              do i1=1,n1i
                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                 maxdiff=max(maxdiff,abs(potion_corr(ind)-pot_ion(ind)))
                 !write(17,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind),maxdiff
              end do
           end do
        end do

        call mpiallred(maxdiff,1,MPI_MAX,MPI_COMM_WORLD,ierr)

        if (iproc == 0) write(*,'(1x,a,1pe24.17)')'...done. MaxDiff=',maxdiff

        stop

        i_all=-product(shape(potion_corr))*kind(potion_corr)
        deallocate(potion_corr,stat=i_stat)
        call memocc(i_stat,i_all,'potion_corr',subname)

     end if

  end if


!!!  !calculate the value of the offset to be put
!!!  tt_tot=0.d0
!!!  do ind=1,n1i*n2i*n3i
!!!     tt_tot=tt_tot+pot_ion(ind)
!!!  end do
!!!  print *,'previous offset',tt_tot*hxh*hyh*hzh

  if (n3pi > 0) then
     do iat=1,at%nat
        ityp=at%iatype(iat)

        rx=rxyz(1,iat)
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)

        ! determine number of local terms
        nloc=0
        do iloc=1,4
           if (at%psppar(0,iloc,ityp) /= 0.d0) nloc=iloc
        enddo
        rloc=at%psppar(0,0,ityp)
        cutoff=10.d0*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)
        
        !do not add the local part for the vacancy
        if (nloc /= 0) then

           do i3=isz,iez
              z=real(i3,kind=8)*hzh-rz
              call ind_positions(perz,i3,n3,j3,goz) 
              j3=j3+nbl3+1
              if (goz .and. j3 >= i3s .and. j3 <=  i3s+n3pi-1) then
                 do i2=isy,iey
                    y=real(i2,kind=8)*hyh-ry
                    call ind_positions(pery,i2,n2,j2,goy)
                    if (goy) then
                       do i1=isx,iex
                          x=real(i1,kind=8)*hxh-rx
                          call ind_positions(perx,i1,n1,j1,gox)
                          if (gox) then
                             r2=x**2+y**2+z**2
                             arg=r2/rloc**2
                             xp=exp(-.5d0*arg)
                             tt=at%psppar(0,nloc,ityp)
                             do iloc=nloc-1,1,-1
                                tt=arg*tt+at%psppar(0,iloc,ityp)
                             enddo
                             ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                             pot_ion(ind)=pot_ion(ind)+xp*tt
                          end if
                       enddo
                    end if
                 enddo
              end if
           end do

        end if

     enddo

     if (htoobig) then
        !add to pot_ion an explicit error function to correct in the case of big grid spacing
        !for the moment works only in the isolated BC case
        do i3=1,n3pi
           z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,n2i
              y=real(i2-nbl2-1,gp)*hyh
              do i1=1,n1i
                 x=real(i1-nbl1-1,gp)*hxh
                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                 call sum_erfcr(at%nat,at%ntypes,x,y,z,at%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
                 pot_ion(ind)=pot_ion(ind)+potxyz
              end do
           end do
        end do
     end if
     
  end if
  
  if (nvacancy /= 0) then
     !for a vacancy reference calculation, save the ionic potential to be used
     !in the following run

     !first calculate the kernel in isolated BC
     call timing(iproc,'CrtLocPot     ','OF')
     call createKernel(iproc,nproc,'F',n1i,n2i,n3i,hxh,hyh,hzh,16,pkernel_ref)
     call timing(iproc,'CrtLocPot     ','ON')


     !calculate the ionic potential correction in the global data distribution
     allocate(potion_corr(n1i*n2i*n3i+ndebug),stat=i_stat)
     call memocc(i_stat,potion_corr,'potion_corr',subname)

     call razero(n1i*n2i*n3i,potion_corr)

     iat=nvacancy
     ityp=at%iatype(iat)
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)

     rloc=at%psppar(0,0,ityp)
     charge=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
     cutoff=10.d0*rloc

     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)

     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

     do i3=isz,iez
        z=real(i3,kind=8)*hzh-rz
        call ind_positions(perz,i3,n3,j3,goz) 
        j3=j3+nbl3+1
        do i2=isy,iey
           y=real(i2,kind=8)*hyh-ry
           call ind_positions(pery,i2,n2,j2,goy)
           do i1=isx,iex
              x=real(i1,kind=8)*hxh-rx
              call ind_positions(perx,i1,n1,j1,gox)
              r2=x**2+y**2+z**2
              arg=r2/rloc**2
              xp=exp(-.5d0*arg)
              if (goz  .and. goy  .and. gox ) then
                 ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-1)*n1i*n2i
                 potion_corr(ind)=xp*charge !the sign is inverted here
              endif
           enddo
        enddo
     enddo

!!$     !plot the ionic potential in a .pot file
!!$     !allocate the arrays for plotting
!!$     !its values are ignored in the datacode='G' case
!!$     allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
!!$     call memocc(i_stat,nscatterarr,'nscatterarr',subname)
!!$     allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
!!$     call memocc(i_stat,ngatherarr,'ngatherarr',subname)
!!$     !create the descriptors for the density and the potential
!!$     !these descriptors should take into account the localisation regions
!!$     call createDensPotDescriptors(iproc,nproc,at%geocode,'D',n1i,n2i,n3i,0,&
!!$          n3d_fake,n3p_fake,n3pi_fake,i3xcsh_fake,i3s_fake,nscatterarr,ngatherarr)
!!$
!!$     i_all=-product(shape(nscatterarr))*kind(nscatterarr)
!!$     deallocate(nscatterarr,stat=i_stat)
!!$     call memocc(i_stat,i_all,'nscatterarr',subname)
!!$
!!$
!!$     call plot_density(at%geocode,'gaupotion.pot',iproc,1,n1,n2,n3,n1i,n2i,n3i,n3i,&
!!$          at%alat1,at%alat2,at%alat3,ngatherarr,potion_corr)



     call timing(iproc,'CrtLocPot     ','OF')
     !here the value of the datacode must be kept fixed
     call H_potential('F','G',iproc,nproc,&
          n1i,n2i,n3i,hxh,hyh,hzh,&
          potion_corr,pkernel_ref,potion_corr,ehart,0.0_gp,.false.)

!!$     call PSolver('F','G',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
!!$          potion_corr,pkernel_ref,potion_corr,ehart,eexcu,vexcu,0.0_gp,.false.,1)
     call timing(iproc,'CrtLocPot     ','ON')


     i_all=-product(shape(pkernel_ref))*kind(pkernel_ref)
     deallocate(pkernel_ref,stat=i_stat)
     call memocc(i_stat,i_all,'pkernel_ref',subname)


!!!     call plot_density(at%geocode,'deltapotion.pot',iproc,1,n1,n2,n3,n1i,n2i,n3i,n3i,&
!!!          at%alat1,at%alat2,at%alat3,ngatherarr,potion_corr)


     iat=nvacancy
     ityp=at%iatype(iat)

     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)

     ! determine number of local terms
     nloc=0
     do iloc=1,4
        if (at%psppar(0,iloc,ityp) /= 0.d0) nloc=iloc
     enddo
     rloc=at%psppar(0,0,ityp)
     cutoff=10.d0*rloc

     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)

     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

     !do not add the local part for the vacancy
     if (nloc /= 0) then

        do i3=isz,iez
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,n3,j3,goz) 
           j3=j3+nbl3+1
           if (goz) then
              do i2=isy,iey
                 y=real(i2,kind=8)*hyh-ry
                 call ind_positions(pery,i2,n2,j2,goy)
                 if (goy) then
                    do i1=isx,iex
                       x=real(i1,kind=8)*hxh-rx
                       call ind_positions(perx,i1,n1,j1,gox)
                       if (gox) then
                          r2=x**2+y**2+z**2
                          arg=r2/rloc**2
                          xp=exp(-.5d0*arg)
                          tt=at%psppar(0,nloc,ityp)
                          do iloc=nloc-1,1,-1
                             tt=arg*tt+at%psppar(0,iloc,ityp)
                          enddo
                          ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-1)*n1i*n2i
                          potion_corr(ind)=&
                               potion_corr(ind)-xp*tt ! the sign has changed here
                       end if
                    enddo
                 end if
              enddo
           end if
        end do

     end if


     !call plot_density(at%geocode,'deltapotion_final.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3pi,&
     !     at%alat1,at%alat2,at%alat3,ngatherarr,potion_corr)

     !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !add the periodic pot_ion
     ind=1+(i3s-1)*n1i*n2i
     call axpy(n1i*n2i*n3pi,1.0_dp,pot_ion(1),1,potion_corr(ind),1)

     if (correct_offset) then
        !calculate the offset
        tt=0.d0
        do i1=1,n3pi*n2i*n1i
           tt=tt+potion_corr(ind-1+i1)
        enddo
        !tt=tt*hxh*hyh*hzh
        if (nproc > 1) then
           call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
                MPI_SUM,MPI_COMM_WORLD,ierr)
        else
           offset=tt
        end if

        if (iproc==0) print *,'offset for potion',offset

        !now potion_corr has zero integral
        potion_corr=potion_corr-offset/real(n1i*n2i*n3i,dp)

        !calculate the offset
        tt=0.d0
        do i1=1,n3pi*n2i*n1i
           tt=tt+potion_corr(ind-1+i1)
        enddo
        !tt=tt*hxh*hyh*hzh
        if (nproc > 1) then
           call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
                MPI_SUM,MPI_COMM_WORLD,ierr)
        else
           offset=tt
        end if

        if (iproc==0) print *,'offset recheck',offset
     end if

     !here put nproc=1 
     call plot_density('potion_corr.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3pi,&
          at%alat1,at%alat2,at%alat3,ngatherarr,potion_corr(ind))



!!!     !reread file from disk
!!!     !overwrite pot_ion with the potential previously created
!!!     call read_potfile(at%geocode,'potion_corr.pot',n1,n2,n3,n1i,n2i,n3i,n3pi,i3s,1,potion_corr)
!!!
!!!     !calculate the offset
!!!     tt=0.d0
!!!     do ind=1,n3pi*n2i*n1i
!!!        tt=tt+potion_corr(ind)
!!!     enddo
!!!     !tt=tt*hxh*hyh*hzh
!!!     if (nproc > 1) then
!!!        call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
!!!             MPI_SUM,MPI_COMM_WORLD,ierr)
!!!     else
!!!        offset=tt
!!!     end if
!!!
!!!     if (iproc==0) print *,'offset reread',offset
!!!
!!!     call plot_density(at%geocode,'potion_corr_2.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,&
!!!          n3pi,at%alat1,at%alat2,at%alat3,ngatherarr,potion_corr)
          

     i_all=-product(shape(ngatherarr))*kind(ngatherarr)
     deallocate(ngatherarr,stat=i_stat)
     call memocc(i_stat,i_all,'ngatherarr',subname)
     i_all=-product(shape(potion_corr))*kind(potion_corr)
     deallocate(potion_corr,stat=i_stat)
     call memocc(i_stat,i_all,'potion_corr',subname)


  end if

!!!  !calculate the value of the offset to be put
!!!  tt_tot=0.d0
!!!  do ind=1,n1i*n2i*n3i
!!!     tt_tot=tt_tot+pot_ion(ind)
!!!  end do
!!!  print *,'actual offset',tt_tot*hxh*hyh*hzh

  !use rhopot to calculate the potential from a constant electric field along y direction
  if (.not. all(elecfield(1:3) == 0.0_gp)) then
     !constant electric field allowed only for surface and free BC
     if (geocode == 'P') then
     !if (iproc == 0) 
           write(*,'(1x,a)') &
          'The constant electric field is not allowed for Fully Periodic BC.'
          !'The constant electric field is allowed only for Free and Surfaces BC'
     stop
      !constant electric field allowed for surface BC only normal to the surface
     elseif (geocode == 'S' .and. (elecfield(1) /= 0.0_gp .or. elecfield(3) /= 0.0_gp) ) then
     !if (iproc == 0) 
           write(*,'(1x,a)') &
          'Only normal constant electric field (Ex=Ez=0) is allowed for Surface BC.'
     stop
     end if
     if (iproc == 0) write(*,'(1x,a,"(",es10.2,", ",es10.2,", ",es10.2,") ", a)') &
          'Constant electric field ',elecfield(1:3),' Ha/Bohr'
!or         'Parabolic confining potential: rprb=',elecfield,&
!           ';  v_conf(r)= 1/(2*rprb**4) * r**2'

     !write or not electric field in a separate file
     efwrite=.true.

     if (n3pi > 0) then
        do i3=1,n3pi
           z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,n2i
              y=real(i2-nbl2-1,gp)*hyh
                 do i1=1,n1i
                    x=real(i1-nbl1-1,gp)*hxh
                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    pot_ion(ind)=pot_ion(ind)+elecfield(1)*x+elecfield(2)*y+elecfield(3)*z
!                    parabola: these two lines replace the above line 
!                              comment out the if case and calculate x, z
!                    r2=(x-rx)**2+(y-ry)**2+(z-rz)**2
!                    pot_ion(ind)=pot_ion(ind)+0.5_gp/(elecfield**4)*r2
                 end do
           end do
        end do

        if (efwrite .and. iproc == 0) then
           open(unit=17,file='elecpotential_x',status='unknown')
           write(17,*) "# x , external electric potential(x,y=0,z=0)"
           do i1=nbl1+1,n1i-nbr1-1
              x=real(i1-nbl1-1,gp)*hxh
              write(17,*)x,-elecfield(1)*x
           end do
           close(17)
           open(unit=17,file='elecpotential_y',status='unknown')
           write(17,*) "# y , external electric potential(x=0,y,z=0)"
           do i2=nbl2+1,n2i-nbr2-1
              y=real(i2-nbl2-1,gp)*hyh
              write(17,*)y,-elecfield(2)*y
           end do
           close(17)
           open(unit=17,file='elecpotential_z',status='unknown')
           write(17,*) "# z , external electric potential(x=0,y=0,z)"
           do i3=1,n3pi
              z=real(i3+i3s-1-nbl3-1,gp)*hzh
              write(17,*)z,-elecfield(3)*z
           end do
           close(17)
        end if

     end if
  end if

  call timing(iproc,'CrtLocPot     ','OF')

END SUBROUTINE createIonicPotential


!>   Determine the index in which the potential must be inserted, following the BC
!!   Determine also whether the index is inside or outside the box for free BC
subroutine ind_positions(periodic,i,n,j,go)
  implicit none
  logical, intent(in) :: periodic
  integer, intent(in) :: i,n
  logical, intent(out) :: go
  integer, intent(out) :: j

  if (periodic) then
     go=.true.
     j=modulo(i,2*n+2)
  else
     j=i
     if (i >= -14 .and. i <= 2*n+16) then
        go=.true.
     else
        go=.false.
     end if
  end if

END SUBROUTINE ind_positions


subroutine sum_erfcr(nat,ntypes,x,y,z,iatype,nelpsp,psppar,rxyz,potxyz)
  use module_base
  implicit none
  integer, intent(in) :: nat,ntypes
  real(gp) :: x,y,z
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(gp), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(wp), intent(out) :: potxyz
  !local variables
  integer :: iat,ityp
  real(wp) :: pi,charge
  real(gp) :: r,sq2rl,rx,ry,rz,derf_val
  
  pi=4.0_wp*atan(1.0_wp)

  potxyz =0.0_wp

  do iat=1,nat

     ityp=iatype(iat)
     sq2rl=sqrt(2.0_gp)*psppar(0,0,ityp)
     charge=real(nelpsp(ityp),wp)

     rx=rxyz(1,iat)-x 
     ry=rxyz(2,iat)-y
     rz=rxyz(3,iat)-z

     r=sqrt(rx**2+ry**2+rz**2)

     if (r == 0.0_gp) then
        potxyz = potxyz - charge*2.0_wp/(sqrt(pi)*real(sq2rl,wp))
     else
        call derf_ab(derf_val,r/sq2rl)
        potxyz = potxyz - charge*real(derf_val/r,wp)
     end if

     !write(*,'(a,1x,i0,3(1pe24.17))')'iat,r,derf_val,derf_val/r',iat,r/sq2rl,derf_val,derf_val/r

  end do

END SUBROUTINE sum_erfcr


subroutine ext_buffers(periodic,nl,nr)
  implicit none
  logical, intent(in) :: periodic
  integer, intent(out) :: nl,nr

  if (periodic) then
     nl=0
     nr=0
  else
     nl=14
     nr=15
  end if
END SUBROUTINE ext_buffers


subroutine CounterIonPotential(geocode,iproc,nproc,in,shift,&
     hxh,hyh,hzh,grid,n3pi,i3s,pkernel,pot_ion)
  use module_base
  use module_types
  use module_interfaces, except_this_one => CounterIonPotential
  use Poisson_Solver
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,n3pi,i3s
  real(gp), intent(in) :: hxh,hyh,hzh
  real(gp), dimension(3), intent(in) :: shift
  type(input_variables), intent(in) :: in
  type(grid_dimensions), intent(in) :: grid
  real(dp), dimension(*), intent(in) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  !local variables
  character(len=*), parameter :: subname='CounterIonPotential'
  logical :: htoobig=.false.,check_potion=.false.
  logical :: perx,pery,perz,gox,goy,goz
  integer :: iat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp,nspin
  integer :: ind,i_all,i_stat,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
  real(kind=8) :: pi,rholeaked,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot,potxyz
  real(wp) :: maxdiff
  real(gp) :: ehart
  type(atoms_data) :: at
  real(dp), dimension(2) :: charges_mpi
  real(dp), dimension(:), allocatable :: potion_corr
  real(gp), dimension(:,:), allocatable :: radii_cf
  real(gp), dimension(:,:), pointer :: rxyz

  call timing(iproc,'CrtLocPot     ','ON')

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '--------------------------------------------------- Counter Ionic Potential Creation'
  end if

  !read the positions of the counter ions from file
  call read_atomic_file('posinp_ci',iproc,at,rxyz)
  ! Read associated pseudo files.
  call init_atomic_values((iproc == 0), at, in%ixc)

  allocate(radii_cf(at%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

  !read the specifications of the counter ions from pseudopotentials
  call read_atomic_variables('input.occup',iproc,in,at,radii_cf)

  pi=4.d0*atan(1.d0)
  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  !Creates charge density arising from the ionic PSP cores
  call razero(grid%n1i*grid%n2i*n3pi,pot_ion)


  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  if (n3pi >0 .and. .not. htoobig) then

     do iat=1,at%nat
        ityp=at%iatype(iat)
        !shift the positions of the counter_ion wrt the box
        rx=rxyz(1,iat)-shift(1)
        ry=rxyz(2,iat)-shift(2)
        rz=rxyz(3,iat)-shift(3)

        if (iproc == 0) then
           write(*,'(1x,a,i6,3(1pe14.7))')'counter ion No. ',iat,rx,ry,rz
        end if

        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
        cutoff=10.d0*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        !print *,'rloc,iat,nelpsp',isx,iex,isy,iey,isz,iez,shift(:),iproc

        do i3=isz,iez
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,grid%n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=isy,iey
              y=real(i2,kind=8)*hyh-ry
              call ind_positions(pery,i2,grid%n2,j2,goy)
              do i1=isx,iex
                 x=real(i1,kind=8)*hxh-rx
                 call ind_positions(perx,i1,grid%n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5d0*arg)
                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*grid%n1i+(j3-i3s+1-1)*grid%n1i*grid%n2i
                    pot_ion(ind)=pot_ion(ind)-xp*charge
                 else if (.not. goz ) then
                    rholeaked=rholeaked+xp*charge
                 endif
              enddo
           enddo
        enddo

     enddo

  end if

  ! Check
  tt=0.d0
  do j3=1,n3pi
     do i2= -nbl2,2*grid%n2+1+nbr2
        do i1= -nbl1,2*grid%n1+1+nbr1
           ind=i1+1+nbl1+(i2+nbl2)*grid%n1i+(j3-1)*grid%n1i*grid%n2i
           tt=tt+pot_ion(ind)
        enddo
     enddo
  enddo

  tt=tt*hxh*hyh*hzh
  rholeaked=rholeaked*hxh*hyh*hzh

  if (nproc > 1) then
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked

     call mpiallred(charges_mpi(1),2,MPI_SUM,MPI_COMM_WORLD,ierr)

     tt_tot=charges_mpi(1)
     rholeaked_tot=charges_mpi(2)
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
  end if

  if (iproc == 0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
       'total ionic charge, leaked charge ',tt_tot,rholeaked_tot

  if (.not. htoobig) then
     call timing(iproc,'CrtLocPot     ','OF')
     !here the value of the datacode must be kept fixed
     nspin=1

     call H_potential(geocode,'D',iproc,nproc,&
          grid%n1i,grid%n2i,grid%n3i,hxh,hyh,hzh,&
          pot_ion,pkernel,pot_ion,ehart,0.0_gp,.false.)

     call timing(iproc,'CrtLocPot     ','ON')
     
     if (check_potion) then
        if (iproc == 0) write(*,'(1x,a)',advance='no') &
             'Check the ionic potential...'
          
        allocate(potion_corr(grid%n1i*grid%n2i*n3pi+ndebug),stat=i_stat)
        call memocc(i_stat,potion_corr,'potion_corr',subname)

        call razero(grid%n1i*grid%n2i*n3pi,potion_corr)

        !calculate pot_ion with an explicit error function to correct in the case of big grid spacings
        !for the moment works only in the isolated BC case
        do i3=1,n3pi
           z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,grid%n2i
              y=real(i2-nbl2-1,gp)*hyh
              do i1=1,grid%n1i
                 x=real(i1-nbl1-1,gp)*hxh
                 ind=i1+(i2-1)*grid%n1i+(i3-1)*grid%n1i*grid%n2i
                 !if (i1==49 .and. i2==46 .and. i3==44) then
                    call sum_erfcr(at%nat,at%ntypes,x,y,z,at%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
                 !   stop
                 !end if
                 potion_corr(ind)=potion_corr(ind)+potxyz
                 !write(18,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind)
              end do
           end do
        end do

        !then calculate the maximum difference in the sup norm
        maxdiff=0.0_wp
        do i3=1,n3pi
           do i2=1,grid%n2i
              do i1=1,grid%n1i
                 ind=i1+(i2-1)*grid%n1i+(i3-1)*grid%n1i*grid%n2i
                 maxdiff=max(maxdiff,abs(potion_corr(ind)-pot_ion(ind)))
                 !write(17,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind),maxdiff
              end do
           end do
        end do

        call mpiallred(maxdiff,1,MPI_MAX,MPI_COMM_WORLD,ierr)

        if (iproc == 0) write(*,'(1x,a,1pe24.17)')'...done. MaxDiff=',maxdiff

        stop

        i_all=-product(shape(potion_corr))*kind(potion_corr)
        deallocate(potion_corr,stat=i_stat)
        call memocc(i_stat,i_all,'potion_corr',subname)

     end if

  end if

  if (n3pi > 0 .and. htoobig) then
     !add to pot_ion an explicit error function to correct in the case of big grid spacing
     !for the moment works only in the isolated BC case
     do i3=1,n3pi
        z=real(i3+i3s-1-nbl3-1,gp)*hzh
        do i2=1,grid%n2i
           y=real(i2-nbl2-1,gp)*hyh
           do i1=1,grid%n1i
              x=real(i1-nbl1-1,gp)*hxh
              ind=i1+(i2-1)*grid%n1i+(i3-1)*grid%n1i*grid%n2i
              call sum_erfcr(at%nat,at%ntypes,x,y,z,at%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
              pot_ion(ind)=pot_ion(ind)+potxyz
           end do
        end do
     end do
  end if

  !deallocations
  call deallocate_atoms(at,subname) 

  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf',subname)

  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz',subname)


  call timing(iproc,'CrtLocPot     ','OF')

END SUBROUTINE CounterIonPotential
