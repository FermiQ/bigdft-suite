!!****p* PSolver/PS_Check
!! FUNCTION
!!    Performs a check of the Poisson Solver suite by running with different regimes
!!    and for different choices of the XC functionals
!! COPYRIGHT
!!    Copyright (C) 2002-2009 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2007
!! SOURCE
!!
program PS_Check

  use module_base
  use Poisson_Solver

  implicit none
  !include 'mpif.h'
  !Length of the box
  character(len=*), parameter :: subname='PS_Check'
  real(kind=8), parameter :: a_gauss = 1.0d0,a2 = a_gauss**2
  real(kind=8), parameter :: acell = 10.d0
  character(len=50) :: chain
  character(len=1) :: geocode
  real(kind=8), dimension(:), allocatable :: density,rhopot,potential,pot_ion,xc_pot
  real(kind=8), pointer :: pkernel(:)
  real(kind=8) :: hx,hy,hz,offset
  real(kind=8) :: ehartree,eexcu,vexcu
  integer :: n01,n02,n03,itype_scf,i_all,i_stat
  integer :: iproc,nproc,ierr,ispden
  integer :: n_cell,ixc
  integer, dimension(4) :: nxyz

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  !initialize memory counting and timings
  !call memocc(0,iproc,'count','start')
  call timing(iproc,'parallel      ','IN')

  !the first proc read the data and then send them to the others
  if (iproc==0) then
     !Use arguments
     call getarg(1,chain)
     read(unit=chain,fmt=*) nxyz(1)
     call getarg(2,chain)
     read(unit=chain,fmt=*) nxyz(2)
     call getarg(3,chain)
     read(unit=chain,fmt=*) nxyz(3)
     call getarg(4,chain)
     read(unit=chain,fmt=*) nxyz(4)
     call getarg(5,chain)
     read(unit=chain,fmt=*) geocode
  end if

  call MPI_BCAST(nxyz,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(geocode,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

  n01=nxyz(1)
  n02=nxyz(2)
  n03=nxyz(3)
  ixc=nxyz(4)

  !print *,iproc,n01,n02,n03

  !Step size
  n_cell = max(n01,n02,n03)
  hx=acell/real(n01,kind=8)
  hy=acell/real(n02,kind=8)
  hz=acell/real(n03,kind=8)

  !order of the scaling functions chosen
  itype_scf=16

  !calculate the kernel in parallel for each processor
  call createKernel(iproc,nproc,geocode,n01,n02,n03,hx,hy,hz,itype_scf,pkernel,quiet='yes')

  !Allocations, considering also spin density
  !Density
  allocate(density(n01*n02*n03*2+ndebug),stat=i_stat)
  call memocc(i_stat,density,'density',subname)
  !Density then potential
  allocate(potential(n01*n02*n03+ndebug),stat=i_stat)
  call memocc(i_stat,potential,'potential',subname)
  !ionic potential
  allocate(pot_ion(n01*n02*n03+ndebug),stat=i_stat)
  call memocc(i_stat,pot_ion,'pot_ion',subname)
  !XC potential
  allocate(xc_pot(n01*n02*n03*2+ndebug),stat=i_stat)
  call memocc(i_stat,xc_pot,'xc_pot',subname)

  do ispden=1,2
     if (iproc == 0) write(unit=*,fmt="(1x,a,i0)")  '===================== nspden:  ',ispden
     !then assign the value of the analytic density and the potential
     !allocate the rhopot also for complex routines
     allocate(rhopot(n01*n02*n03*2+ndebug),stat=i_stat)
     call memocc(i_stat,rhopot,'rhopot',subname)
  
     call test_functions(geocode,ixc,n01,n02,n03,ispden,acell,a_gauss,hx,hy,hz,&
          density,potential,rhopot,pot_ion,offset)
     !calculate the Poisson potential in parallel
     !with the global data distribution (also for xc potential)

     call XC_potential(geocode,'G',iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
          rhopot,eexcu,vexcu,ispden,xc_pot)
     call H_potential(geocode,'G',iproc,nproc,n01,n02,n03,hx,hy,hz,&
          rhopot,pkernel,xc_pot,ehartree,offset,.false.) !optional argument
!!$     call PSolver(geocode,'G',iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
!!$          rhopot,pkernel,xc_pot,ehartree,eexcu,vexcu,offset,.false.,ispden)
     
     if (iproc == 0) write(unit=*,fmt="(1x,a,3(1pe20.12))") 'Energies:',ehartree,eexcu,vexcu
     !stop
     if (iproc == 0) then
        !compare the values of the analytic results (nproc == -1 indicates that it is serial)
        call compare(0,-1,n01,n02,n03,1,potential,rhopot,'ANALYTIC')
     end if
     !if the latter test pass, we have a reference for all the other calculations
     !build the reference quantities (based on the numerical result, not the analytic)
     potential(:)=rhopot(1:n01*n02*n03)

     !now the parallel calculation part
     i_all=-product(shape(rhopot))*kind(rhopot)
     deallocate(rhopot,stat=i_stat)
     call memocc(i_stat,i_all,'rhopot',subname)

     call compare_with_reference(iproc,nproc,geocode,'G',n01,n02,n03,ixc,ispden,hx,hy,hz,&
          offset,ehartree,eexcu,vexcu,&
          density,potential,pot_ion,xc_pot,pkernel)
     
     call compare_with_reference(iproc,nproc,geocode,'D',n01,n02,n03,ixc,ispden,hx,hy,hz,&
          offset,ehartree,eexcu,vexcu,&
          density,potential,pot_ion,xc_pot,pkernel)

     !test for the serial solver (always done to have a simpler comparison)
     !if (iproc == 0 .and. nproc > 1 ) then
     if (iproc == 0) then
        i_all=-product(shape(pkernel))*kind(pkernel)
        deallocate(pkernel,stat=i_stat)
        call memocc(i_stat,i_all,'pkernel',subname)

        !calculate the kernel 
        call createKernel(0,1,geocode,n01,n02,n03,hx,hy,hz,itype_scf,pkernel,quiet='yes')

        call compare_with_reference(0,1,geocode,'G',n01,n02,n03,ixc,ispden,hx,hy,hz,&
             offset,ehartree,eexcu,vexcu,&
             density,potential,pot_ion,xc_pot,pkernel)

        call compare_with_reference(0,1,geocode,'D',n01,n02,n03,ixc,ispden,hx,hy,hz,&
             offset,ehartree,eexcu,vexcu,&
             density,potential,pot_ion,xc_pot,pkernel)
     end if

     if (ixc == 0) exit
  end do

  if (ixc == 0) then
     !compare the calculations in complex
     call compare_cplx_calculations(iproc,nproc,geocode,'G',n01,n02,n03,hx,hy,hz,ehartree,&
          density,potential,pkernel)
     
     call compare_cplx_calculations(iproc,nproc,geocode,'D',n01,n02,n03,hx,hy,hz,ehartree,&
          density,potential,pkernel)
  end if

  i_all=-product(shape(pkernel))*kind(pkernel)
  deallocate(pkernel,stat=i_stat)
  call memocc(i_stat,i_all,'pkernel',subname)

  i_all=-product(shape(density))*kind(density)
  deallocate(density,stat=i_stat)
  call memocc(i_stat,i_all,'density',subname)
  i_all=-product(shape(potential))*kind(potential)
  deallocate(potential,stat=i_stat)
  call memocc(i_stat,i_all,'potential',subname)
  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',subname)
  i_all=-product(shape(xc_pot))*kind(xc_pot)
  deallocate(xc_pot,stat=i_stat)
  call memocc(i_stat,i_all,'xc_pot',subname)

  call timing(iproc,'              ','RE')
  !finalize memory counting
  call memocc(0,0,'count','stop')

  call MPI_FINALIZE(ierr)  

!end program PS_Check


contains

  subroutine compare_cplx_calculations(iproc,nproc,geocode,distcode,n01,n02,n03,hx,hy,hz,ehref,&
       density,potential,pkernel)
    use module_base
    use Poisson_Solver
    implicit none
    character(len=1), intent(in) :: geocode,distcode
    integer, intent(in) :: iproc,nproc,n01,n02,n03
    real(kind=8), intent(in) :: hx,hy,hz,ehref
    real(kind=8), dimension(n01*n02*n03), intent(in) :: potential
    real(kind=8), dimension(n01*n02*n03*2), intent(in) :: density
    real(kind=8), dimension(:), pointer :: pkernel
    !local varaibles
    character(len=*), parameter :: subname='compare_cplx_calculations'
    character(len=20) :: message
    integer :: n3d,n3p,n3pi,i3xcsh,i3s,i3sd,i3,i2,i1,istden,istpot,isp,i
    real(kind=8) :: eexcu,vexcu,ehartree,offset
    real(kind=8), dimension(:,:,:,:), allocatable :: rhopot

    offset=0.d0

    !this is performed always without XC since a complex
    !charge density makes no sense
    write(message,'(1x,a,1x,a)') geocode,distcode

    call PS_dim4allocation(geocode,distcode,iproc,nproc,n01,n02,n03,0,&
         n3d,n3p,n3pi,i3xcsh,i3s)

    !starting point of the three-dimensional arrays
    if (distcode == 'D') then
       istden=n01*n02*(i3s-1)+1
       istpot=n01*n02*(i3s+i3xcsh-1)+1
       i3sd=i3s
    else if (distcode == 'G') then
       istden=1
       istpot=1
       i3sd=1
    end if

    !input poisson solver, complex distribution
    allocate(rhopot(n01,n02,n3d,2+ndebug),stat=i_stat)
    call memocc(i_stat,rhopot,'rhopot',subname)

    !do isp=1,2
       isp=1
       do i3=1,n3d
          do i2=1,n02
             do i1=1,n01
                i=i1+(i2-1)*n01+(modulo(i3sd+i3-2,n03))*n01*n02!+(isp-1)*n01*n02*n03
                rhopot(i1,i2,i3,isp)=density(i)
             end do
          end do
       end do
    !end do
   
    !perform the calculation in complex, with distributed and gathered distribution
     call H_potential(geocode,distcode,iproc,nproc,n01,n02,n03,hx,hy,hz,&
          rhopot,pkernel,rhopot,ehartree,offset,.false.,quiet='YES')
!!$    call PSolver(geocode,distcode,iproc,nproc,n01,n02,n03,0,hx,hy,hz,&
!!$         rhopot,pkernel,rhopot,ehartree,eexcu,vexcu,offset,.false.,1,quiet='YES')

    call compare(iproc,-1,n01,n02,n3d,1,potential(istpot),rhopot,'CPLXREAL'//trim(message))

       isp=2
       do i3=1,n3d
          do i2=1,n02
             do i1=1,n01
                i=i1+(i2-1)*n01+(modulo(i3sd+i3-2,n03))*n01*n02!+(isp-1)*n01*n02*n03
                rhopot(i1,i2,i3,isp)=density(i)
             end do
          end do
       end do

       !perform the calculation in complex, with distributed and gathered distribution
       call H_potential(geocode,distcode,iproc,nproc,n01,n02,n03,hx,hy,hz,&
            rhopot(1,1,1,2),pkernel,rhopot,ehartree,offset,.false.,quiet='YES')
!!$       call PSolver(geocode,distcode,iproc,nproc,n01,n02,n03,0,hx,hy,hz,&
!!$            rhopot(1,1,1,2),pkernel,rhopot,ehartree,eexcu,vexcu,offset,.false.,1,quiet='YES')
    
    call compare(iproc,-1,n01,n02,n3d,1,potential(istpot),rhopot(1,1,1,2),'CPLXIMAG'//trim(message))


    if (iproc==0) write(unit=*,fmt="(1x,a,1pe20.12)")'Energy diff:',ehref-ehartree

    
    i_all=-product(shape(rhopot))*kind(rhopot)
    deallocate(rhopot,stat=i_stat)
    call memocc(i_stat,i_all,'rhopot',subname)
    
  end subroutine compare_cplx_calculations



  subroutine compare_with_reference(iproc,nproc,geocode,distcode,n01,n02,n03,&
       ixc,nspden,hx,hy,hz,offset,ehref,excref,vxcref,&
       density,potential,pot_ion,xc_pot,pkernel)
    use module_base
    use Poisson_Solver
    implicit none
    character(len=1), intent(in) :: geocode,distcode
    integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc,nspden
    real(kind=8), intent(in) :: hx,hy,hz,offset,ehref,excref,vxcref
    real(kind=8), dimension(n01*n02*n03), intent(in) :: potential
    real(kind=8), dimension(n01*n02*n03*nspden), intent(in) :: density
    real(kind=8), dimension(n01*n02*n03), intent(inout) :: pot_ion
    real(kind=8), dimension(n01*n02*n03*nspden), target, intent(inout) :: xc_pot
    real(kind=8), dimension(:), pointer :: pkernel
    !local variables
    character(len=*), parameter :: subname='compare_with_reference'
    character(len=100) :: message
    integer :: n3d,n3p,n3pi,i3xcsh,i3s,istden,istpot,i_all,i_stat,istpoti,i
    integer :: istxc,i1,i2,i3,isp,i3sd
    real(kind=8) :: eexcu,vexcu,ehartree
    real(kind=8), dimension(:), allocatable :: test,test_xc
    real(kind=8), dimension(:,:,:,:), allocatable :: rhopot
    real(kind=8), dimension(:), pointer :: xc_temp

    call PS_dim4allocation(geocode,distcode,iproc,nproc,n01,n02,n03,ixc,&
         n3d,n3p,n3pi,i3xcsh,i3s)

    !starting point of the three-dimensional arrays
    if (distcode == 'D') then
       istden=n01*n02*(i3s-1)+1
       istpot=n01*n02*(i3s+i3xcsh-1)+1
       i3sd=i3s
    else if (distcode == 'G') then
       istden=1
       istpot=1
       i3sd=1
    end if
    istpoti=n01*n02*(i3s+i3xcsh-1)+1

    !test arrays for comparison
    allocate(test(n01*n02*n03*nspden+ndebug),stat=i_stat)
    call memocc(i_stat,test,'test',subname)
    !XC potential
    allocate(test_xc(n01*n02*n03*nspden+ndebug),stat=i_stat)
    call memocc(i_stat,test_xc,'test_xc',subname)
    !input poisson solver
    allocate(rhopot(n01,n02,n3d,nspden+ndebug),stat=i_stat)
    call memocc(i_stat,rhopot,'rhopot',subname)

    if (iproc == 0) write(message,'(1x,a,1x,i0,1x,a,1x,i0)') geocode,ixc,distcode,nspden

    if (ixc /= 0) then
       if (nspden == 1) then
          test(1:n01*n02*n03)=&
               potential(1:n01*n02*n03)+pot_ion(1:n01*n02*n03)+xc_pot(1:n01*n02*n03)
       else
          if (distcode == 'G') then
             do i=1,n01*n02*n03
                test(i)=potential(i)+pot_ion(i)+xc_pot(i)
                test(i+n01*n02*n03)=potential(i)+pot_ion(i)+&
                     xc_pot(i+n01*n02*n03)
             end do
          else
             do i=1,n01*n02*n3p
                test(i+istpot-1)=potential(i+istpot-1)+pot_ion(i+istpot-1)+xc_pot(i+istpot-1)
                test(i+istpot-1+n01*n02*n3p)=potential(i+istpot-1)+pot_ion(i+istpot-1)+&
                     xc_pot(i+istpot-1+n01*n02*n03)
             end do
          end if
       end if
    else
       test(1:n01*n02*n03)=potential(1:n01*n02*n03)!+pot_ion
    end if

    do isp=1,nspden
       !add the initialisation of the density for the periodic GGA case
       do i3=1,n3d
          do i2=1,n02
             do i1=1,n01
                i=i1+(i2-1)*n01+(modulo(i3sd+i3-2,n03))*n01*n02+(isp-1)*n01*n02*n03
                rhopot(i1,i2,i3,isp)=density(i)
             end do
          end do
       end do
    end do

    if (nspden == 2 .and. distcode == 'D') then
       allocate(xc_temp(n01*n02*n3p*nspden+ndebug),stat=i_stat)
       call memocc(i_stat,xc_temp,'xc_temp',subname)
       !toggle the components of xc_pot in the distributed case
       do i=1,n01*n02*n3p
          xc_temp(i)=xc_pot(i+istpot-1)
          xc_temp(i+n01*n02*n3p)=xc_pot(i+istpot-1+n01*n02*n03)
       end do
       istxc=1
    else
       xc_temp => xc_pot
       istxc=istpot
    end if

     call XC_potential(geocode,distcode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
          rhopot(1,1,1,1),eexcu,vexcu,nspden,test_xc)
     call H_potential(geocode,distcode,iproc,nproc,n01,n02,n03,hx,hy,hz,&
          rhopot(1,1,1,1),pkernel,rhopot,ehartree,offset,.false.,quiet='yes') !optional argument
     !compare the values of the analytic results (no dependence on spin)
     call compare(iproc,nproc,n01,n02,n3p,1,potential(istpot),rhopot(1,1,1,1),&
          'ANACOMPLET '//message)

!!$    call PSolver(geocode,distcode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
!!$         rhopot(1,1,1,1),pkernel,test_xc,ehartree,eexcu,vexcu,offset,.false.,nspden,quiet='yes')
!!$    
!!$    !compare the values of the analytic results (no dependence on spin)
!!$    call compare(iproc,nproc,n01,n02,n3p,1,potential(istpot),rhopot(1,1,i3xcsh+1,1),&
!!$         'ANACOMPLET '//message)

    !compare also the xc_potential
    if (ixc/=0) call compare(iproc,nproc,n01,n02,nspden*n3p,1,xc_temp(istxc:),&
         test_xc(1),&
         'XCCOMPLETE '//message)
    if (iproc==0) write(unit=*,fmt="(1x,a,3(1pe20.12))") &
         'Energies diff:',ehref-ehartree,excref-eexcu,vxcref-vexcu

    do isp=1,nspden
       do i3=1,n3d
          do i2=1,n02
             do i1=1,n01
                i=i1+(i2-1)*n01+(modulo(i3sd+i3-2,n03))*n01*n02+(isp-1)*n01*n02*n03
                rhopot(i1,i2,i3,isp)=density(i)
             end do
          end do
       end do
    end do

    if (nspden == 2 .and. distcode == 'D') then
       i_all=-product(shape(xc_temp))*kind(xc_temp)
       deallocate(xc_temp,stat=i_stat)
       call memocc(i_stat,i_all,'xc_temp',subname)
    end if

    call XC_potential(geocode,distcode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
            rhopot(1,1,1,1),eexcu,vexcu,nspden,test_xc)

    call H_potential(geocode,distcode,iproc,nproc,n01,n02,n03,hx,hy,hz,&
          rhopot(1,1,1,1),pkernel,pot_ion(istpoti),ehartree,offset,ixc /= 0,quiet='yes') !optional argument
    !fill the other part, for spin, polarised
    if (nspden == 2) then
       !the starting point is not so simple
       if (n3d > n3p) then
          call dcopy(n01*n02*n3p,rhopot(1,1,1,1),1,rhopot(1,1,n3p+1,1),1)
       else
          call dcopy(n01*n02*n3p,rhopot(1,1,1,1),1,rhopot(1,1,1,nspden),1)
       end if
    end if
    !spin up and down together with the XC part
    call axpy(n01*n02*n3p*nspden,1.0_dp,test_xc(1),1,rhopot(1,1,1,1),1)
    !then compare again, but the complete result
    call compare(iproc,-1,n01,n02,nspden*n3p,1,test(istpot),&
         rhopot(1,1,1,1),'COMPLETE   '//message)

!!$    call PSolver(geocode,distcode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
!!$         rhopot(1,1,1,1),pkernel,pot_ion(istpoti),ehartree,eexcu,vexcu,offset,.true.,nspden,quiet='yes')
!!$    !then compare again, but the complete result
!!$    call compare(iproc,nproc,n01,n02,nspden*n3p,1,test(istpot),&
!!$         rhopot(1,1,i3xcsh+1,1),'COMPLETE   '//message)
    if (iproc==0) write(unit=*,fmt="(1x,a,3(1pe20.12))") &
         'Energies diff:',ehref-ehartree,excref-eexcu,vxcref-vexcu

    i_all=-product(shape(test))*kind(test)
    deallocate(test,stat=i_stat)
    call memocc(i_stat,i_all,'test',subname)
    i_all=-product(shape(test_xc))*kind(test_xc)
    deallocate(test_xc,stat=i_stat)
    call memocc(i_stat,i_all,'test_xc',subname)

    i_all=-product(shape(rhopot))*kind(rhopot)
    deallocate(rhopot,stat=i_stat)
    call memocc(i_stat,i_all,'rhopot',subname)

  end subroutine compare_with_reference


! Compare arrays potential and density
! if nproc == -1: serial version i.e. special comparison
  subroutine compare(iproc,nproc,n01,n02,n03,nspden,potential,density,description)
    implicit none
    include 'mpif.h'
    character(len=*), intent(in) :: description
    integer, intent(in) :: iproc,nproc,n01,n02,n03,nspden
    real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,density
    !local variables
    integer :: i1,i2,i3,ierr,i1_max,i2_max,i3_max
    real(kind=8) :: factor,diff_par,max_diff
    max_diff = 0.d0
    i1_max = 1
    i2_max = 1
    i3_max = 1
    do i3=1,n03
       do i2=1,n02 
          do i1=1,n01
             factor=abs(real(nspden,kind=8)*potential(i1,i2,i3)-density(i1,i2,i3))
             if (max_diff < factor) then
                max_diff = factor
                i1_max = i1
                i2_max = i2
                i3_max = i3
             end if
          end do
       end do
    end do

!!!  print *,'iproc,i3xcsh,i3s,max_diff',iproc,i3xcsh,i3s,max_diff

    if (nproc > 1) then
       !extract the max
       call MPI_ALLREDUCE(max_diff,diff_par,1,MPI_double_precision,  &
            MPI_MAX,MPI_COMM_WORLD,ierr)
    else
       diff_par=max_diff
    end if

    if (iproc == 0) then
       if (nproc == -1) then
          if (diff_par > 1.e-10) then
             write(unit=*,fmt="(1x,a,1pe20.12,a)") &
                  trim(description) // '    Max diff:',diff_par,'   <<<< WARNING'
             write(unit=*,fmt="(1x,a,1pe20.12)") &
                  '      result:',density(i1_max,i2_max,i3_max),&
                  '    original:',potential(i1_max,i2_max,i3_max)
             write(*,'(a,3(i0,1x))') '  Max diff at: ',i1_max,i2_max,i3_max
!!!           i3=i3_max
!!!           i1=i1_max
!!!           do i2=1,n02
!!!              !do i1=1,n01
!!!                 write(20,*)i1,i2,potential(i1,i2,i3),density(i1,i2,i3)
!!!              !end do
!!!           end do
!!!           stop
          else
             write(unit=*,fmt="(1x,a,1pe20.12)") &
                  trim(description) // '    Max diff:',diff_par
          end if
       else
          if (diff_par > 1.e-10) then
             write(unit=*,fmt="(1x,a,1pe20.12,a)") &
                  trim(description) // '    Max diff:',diff_par,'   <<<< WARNING'
          else
             write(unit=*,fmt="(1x,a,1pe20.12)") &
                  trim(description) //'    Max diff:',diff_par
          end if
       end if
    end if

    max_diff=diff_par

  end subroutine compare
!!***

!!****f* PS_Check/test_functions
!! FUNCTION
!!    This subroutine builds some analytic functions that can be used for 
!!    testing the poisson solver.
!!    The default choice is already well-tuned for comparison.
!!    WARNING: not all the test functions can be used for all the boundary conditions of
!!    the poisson solver, in order to have a reliable analytic comparison.
!!    The parameters of the functions must be adjusted in order to have a sufficiently localized
!!    function in the isolated direction and an explicitly periodic function in the periodic ones.
!!    Beware of the high-frequency components that may false the results when hgrid is too high.
!! SOURCE
!!
  subroutine test_functions(geocode,ixc,n01,n02,n03,nspden,acell,a_gauss,hx,hy,hz,&
       density,potential,rhopot,pot_ion,offset)
    implicit none
    character(len=1), intent(in) :: geocode
    integer, intent(in) :: n01,n02,n03,ixc,nspden
    real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz
    real(kind=8), intent(out) :: offset
    real(kind=8), dimension(n01,n02,n03), intent(out) :: pot_ion,potential
    real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: density,rhopot
    !local variables
    integer :: i1,i2,i3,ifx,ify,ifz,i
    real(kind=8) :: x1,x2,x3,length,denval,pi,a2,derf,factor,r,r2
    real(kind=8) :: fx,fx2,fy,fy2,fz,fz2,a,ax,ay,az,bx,by,bz,tt

    if (trim(geocode) == 'P') then
       !parameters for the test functions
       length=acell
       a=0.5d0/a_gauss**2
       !test functions in the three directions
       ifx=5
       ify=5
       ifz=5
       !parameters of the test functions
       ax=length
       ay=length
       az=length
       bx=2.d0!real(nu,kind=8)
       by=2.d0!real(nu,kind=8)
       bz=2.d0

!!!     !plot of the functions used
!!!     do i1=1,n03
!!!        x = hx*real(i1,kind=8)!valid if hy=hz
!!!        y = hz*real(i1,kind=8) 
!!!        call functions(x,ax,bx,fx,fx2,ifx)
!!!        call functions(y,az,bz,fz,fz2,ifz)
!!!        write(20,*)i1,fx,fx2,fz,fz2
!!!     end do

       !Initialization of density and potential
       denval=0.d0 !value for keeping the density positive
       do i3=1,n03
          x3 = hz*real(i3-n03/2-1,kind=8)
          call functions(x3,az,bz,fz,fz2,ifz)
          do i2=1,n02
             x2 = hy*real(i2-n02/2-1,kind=8)
             call functions(x2,ay,by,fy,fy2,ify)
             do i1=1,n01
                x1 = hx*real(i1-n01/2-1,kind=8)
                call functions(x1,ax,bx,fx,fx2,ifx)
                do i=1,nspden
                   density(i1,i2,i3,i) = 1.d0/real(nspden,kind=8)*(fx2*fy*fz+fx*fy2*fz+fx*fy*fz2)
                end do
                potential(i1,i2,i3) = -16.d0*datan(1.d0)*fx*fy*fz
                denval=max(denval,-density(i1,i2,i3,1))
             end do
          end do
       end do

       if (ixc==0) denval=0.d0

    else if (trim(geocode) == 'S') then
       !parameters for the test functions
       length=acell
       a=0.5d0/a_gauss**2
       !test functions in the three directions
       ifx=5
       ifz=5
       !non-periodic dimension
       ify=6
       !parameters of the test functions
       ax=length
       az=length
       bx=2.d0!real(nu,kind=8)
       bz=2.d0!real(nu,kind=8)
       !non-periodic dimension
       ay=length
       by=a

       !Initialisation of density and potential
       denval=0.d0 !value for keeping the density positive
       do i3=1,n03
          x3 = hz*real(i3-n03/2-1,kind=8)
          call functions(x3,az,bz,fz,fz2,ifz)
          do i2=1,n02
             x2 = hy*real(i2-n02/2-1,kind=8)
             call functions(x2,ay,by,fy,fy2,ify)
             do i1=1,n01
                x1 = hx*real(i1-n02/2-1,kind=8)
                call functions(x1,ax,bx,fx,fx2,ifx)
                do i=1,nspden
                   density(i1,i2,i3,i) = &
                        1.d0/real(nspden,kind=8)/(16.d0*datan(1.d0))*(fx2*fy*fz+fx*fy2*fz+fx*fy*fz2)
                end do
                potential(i1,i2,i3) = -fx*fy*fz
                denval=max(denval,-density(i1,i2,i3,1))
             end do
          end do
       end do

       if (ixc==0) denval=0.d0

    else if (trim(geocode) == 'F') then

       !grid for the free BC case
       !hgrid=max(hx,hy,hz)

       pi = 4.d0*atan(1.d0)
       a2 = a_gauss**2

       !Normalization
       factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
       !gaussian function
       do i3=1,n03
          x3 = hz*real(i3-n03/2,kind=8)
          do i2=1,n02
             x2 = hy*real(i2-n02/2,kind=8)
             do i1=1,n01
                x1 = hx*real(i1-n01/2,kind=8)
                r2 = x1*x1+x2*x2+x3*x3
                do i=1,nspden
                   density(i1,i2,i3,i) = 1.d0/real(nspden,kind=8)*max(factor*exp(-r2/a2),1d-24)
                end do
                r = sqrt(r2)
                !Potential from a gaussian
                if (r == 0.d0) then
                   potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
                else
                   potential(i1,i2,i3) = derf(r/a_gauss)/r
                end if
             end do
          end do
       end do
       
       denval=0.d0

    else

       print *,'geometry code not admitted',geocode
       stop

    end if

  ! For ixc/=0 the XC potential is added to the solution, and an analytic comparison is no more
  ! possible. In that case the only possible comparison is between the serial and the parallel case
  ! To ease the comparison between the serial and the parallel case we add a random pot_ion
  ! to the potential.

    if (denval /= 0.d0) then
       rhopot(:,:,:,:) = density(:,:,:,:) + denval +1.d-14
    else
       rhopot(:,:,:,:) = density(:,:,:,:) 
    end if

    offset=0.d0
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01
             tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
             pot_ion(i1,i2,i3)=tt
             offset=offset+potential(i1,i2,i3)
             !add the case for offset in the surfaces case 
             !(for periodic case it is absorbed in offset)
             if (geocode == 'S' .and. denval /= 0.d0) then
                x2 = hy*real(i2-1,kind=8)-0.5d0*acell+0.5d0*hy
                potential(i1,i2,i3)=potential(i1,i2,i3)&
                     -8.d0*datan(1.d0)*denval*real(nspden,kind=8)*(x2**2+0.25d0*acell**2)
                !this stands for
                !denval*2pi*Lx*Lz/Ly^2(y^2-Ly^2/4), less accurate in hgrid
             end if

!!!           if (rhopot(i1,i2,i3,1) <= 0.d0) then
!!!              print *,i1,i2,i3,rhopot(i1,i2,i3,1),denval
!!!           end if
          end do
       end do
    end do
    if (denval /= 0.d0) density=rhopot
    offset=offset*hx*hy*hz

    !print *,'offset',offset

  end subroutine test_functions


  subroutine functions(x,a,b,f,f2,whichone)
    implicit none
    integer, intent(in) :: whichone
    real(kind=8), intent(in) :: x,a,b
    real(kind=8), intent(out) :: f,f2
    !local variables
    real(kind=8) :: r,r2,y,yp,ys,factor,pi,g,h,g1,g2,h1,h2
    real(kind=8) :: length,frequency,nu,sigma,agauss

    pi = 4.d0*datan(1.d0)
    select case(whichone)
    case(1)
       !constant
       f=1.d0
       f2=0.d0
    case(2)
       !gaussian of sigma s.t. a=1/(2*sigma^2)
       r2=a*x**2
       f=dexp(-r2)
       f2=(-2.d0*a+4.d0*a*r2)*dexp(-r2)
    case(3)
       !gaussian "shrinked" with a=length of the system
       length=a
       r=pi*x/length
       y=dtan(r)
       yp=pi/length*1.d0/(dcos(r))**2
       ys=2.d0*pi/length*y*yp
       factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
       f2=factor*dexp(-y**2)
       f=dexp(-y**2)
    case(4)
       !cosine with a=length, b=frequency
       length=a
       frequency=b
       r=frequency*pi*x/length
       f=dcos(r)
       f2=-(frequency*pi/length)**2*dcos(r)
    case(5)
       !exp of a cosine, a=length
       nu=2.d0
       r=pi*nu/a*x
       y=dcos(r)
       yp=dsin(r)
       f=dexp(y)
       factor=(pi*nu/a)**2*(-y+yp**2)
       f2=factor*f
    case(6)
       !gaussian times "shrinked" gaussian, sigma=length/10
       length=a
       r=pi*x/length
       y=dtan(r)
       yp=pi/length*1.d0/(dcos(r))**2
       ys=2.d0*pi/length*y*yp
       factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
       g=dexp(-y**2)
       g1=-2.d0*y*yp*g
       g2=factor*dexp(-y**2)
       
       sigma=length/10
       agauss=0.5d0/sigma**2
       r2=agauss*x**2
       h=dexp(-r2)
       h1=-2.d0*agauss*x*h
       h2=(-2.d0*agauss+4.d0*agauss*r2)*dexp(-r2)
       f=g*h
       f2=g2*h+g*h2+2.d0*g1*h1
    case(7)
       !sine with a=length, b=frequency
       length=a
       frequency=b
       r=frequency*pi*x/length
       f=dsin(r)
       f2=-(frequency*pi/length)**2*dsin(r)
    end select

  end subroutine functions


end program PS_Check
!!***


