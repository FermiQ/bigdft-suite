!> @file
!!  Routines to build the kernel used by the Poisson solver
!! @author
!!    Copyright (C) 2006-2017 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Build the kernel of the Poisson operator with
!! surfaces Boundary conditions
!! in an interpolating scaling functions basis.
!! @warning
!!  Beware of the fact that the nonperiodic direction is y!
subroutine Periodic_Kernel(n1,n2,n3,nker1,nker2,nker3,h1,h2,h3,itype_scf,karray,iproc,nproc,&
     !,mu0_screening,alpha,beta,gamma,
     n3pr2,n3pr1)
  use Poisson_Solver, only: dp
  use memory_profiling
  use dynamic_memory
  implicit none
  !Arguments
  integer, intent(in) :: n1,n2,n3          !< Dimensions for the FFT
  integer, intent(in) :: nker1,nker2,nker3 !< Dimensions of the kernel (nker3=n3/2+1) nker(1,2)=n(1,2)/2+1
  integer, intent(in) :: itype_scf         !< Order of the scaling function
  integer, intent(in) :: iproc             !< Process Id
  integer, intent(in) :: nproc             !< Number of processes
  integer, intent(in) :: n3pr1,n3pr2
  real(dp), intent(in) :: h1,h2,h3         !< Mesh steps in the three dimensions
  real(dp), dimension(nker1,nker2,nker3/nproc), intent(out) :: karray !< output array
  !real(dp), intent(in) :: mu0_screening,alpha,beta,gamma
  !Local variables
  character(len=*), parameter :: subname='Periodic_Kernel'
  real(dp), parameter :: pi=3.14159265358979323846_dp
  integer :: i1,i2,i3,j3,iproc1
  real(dp) :: p1,p2,mu3,h
  real(dp), dimension(:), allocatable :: fourISFx,fourISFy,fourISFz
  !metric for triclinic lattices
  !! ABINIT stuff /// to be fixed
  !scalars
  integer :: iout,id1,id2,id3
  !arrays
!!$   integer(kind=8) :: id(3),ii,ing,ig
!  real(kind=8),allocatable :: gq(:,:)
  !! end of ABINIT stuff
  h=h1*h2*h3

  !!! PSolver n1-n2 plane mpi partitioning !!!

  if (n3pr1 >1) then
   iproc1=mod(iproc,n3pr2)
  else
   iproc1=iproc
  endif

  !first control that the domain is not shorter than the scaling function
  !add also a temporary flag for the allowed ISF types for the kernel
  if (itype_scf > min(n1,n2,n3) .or. itype_scf /= 16) then
     print *,'ERROR: dimension of the box are too small for the ISF basis chosen',&
          itype_scf,n1,n2,n3
     stop
  end if

!!$  detg = 1.0_dp - dcos(alpha)**2 - dcos(beta)**2 - dcos(gamma)**2 + 2.0_dp*dcos(alpha)*dcos(beta)*dcos(gamma)
!!$
  iout=1 !turn the switch on

!!$  rprimd(1,1) = h1*n1
!!$  rprimd(2,1) = 0.0d0
!!$  rprimd(3,1) = 0.0d0
!!$
!!$  rprimd(1,2) = h2*n3*dcos(alpha)
!!$  rprimd(2,2) = h2*n3*dsin(alpha)
!!$  rprimd(3,2) = 0.0d0
!!$
!!$  rprimd(1,3) = n2*h3*dcos(beta)
!!$  rprimd(2,3) = n2*h3*(-dcos(alpha)*dcos(beta)+dcos(gamma))/dsin(alpha)
!!$  rprimd(3,3) = n2*h3*sqrt(detg)/dsin(alpha)
!!$

  !acerioni
  id1=int(n1/2)+2
  id2=int(n2/2)+2
  id3=int(n3/2)+2
  !acerioni

  !calculate the FFT of the ISF for the three dimensions
  fourISFx = f_malloc(0.to.nker1-1,id='fourISFx')
  fourISFy = f_malloc(0.to.nker2-1,id='fourISFy')
  fourISFz = f_malloc(0.to.nker3-1,id='fourISFz')

  call fourtrans_isf(n1/2,fourISFx)
  call fourtrans_isf(n2/2,fourISFy)
  call fourtrans_isf(n3/2,fourISFz)

!!  fourISFx=0._dp
!!  fourISFy=0._dp
!!  fourISFz=0._dp

  !write(*,*) 'screening = ', mu0_screening

  !calculate directly the reciprocal space components of the kernel function
  do i3=1,nker3/nproc
     j3=iproc1*(nker3/nproc)+i3
     if (j3 <= n3/2+1) then
!!$        p3=real(j3-1,dp)/real(n3,dp)
        mu3=real(j3-1,dp)/real(n3,dp)
        mu3=(mu3/h2)**2 !beware of the exchanged dimension

!!$        !acerioni
!!$        ig3=i3-int(i3/id3)*n3-1
!!$        gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
!!$        gqgm23=gq(3,i3)*gmet(2,3)*2
!!$        gqgm13=gq(3,i3)*gmet(1,3)*2
!!$        !acerioni
        do i2=1,nker2
           p2=real(i2-1,dp)/real(n2,dp)
!!$           !acerioni
!!$           ig2=i2-int(i2/id2)*n2-1
!!$           gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
!!$           gqgm12=gq(2,i2)*gmet(1,2)*2
!!$           gqg2p3=gqgm13+gqgm12
!!$           !acerioni
           do i1=1,nker1
              p1=real(i1-1,dp)/real(n1,dp)
              !beware of the exchanged dimension
!              ker = pi*((p1/h1)**2+(p2/h3)**2+mu3)+mu0_screening**2/16.0_dp/datan(1.0_dp)

!!$              !triclinic cell
!!$              !acerioni
!!$              ig1=i1-int(i1/id1)*n1-1
!!$              gs=gs2 + gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
!!$              gs = gs * 4.0d0*datan(1.0d0)
!!$              !gs=gs/16.0d0*datan(1.0d0)
!!$              !write(16,*) i1,i2,i3,gs
!!$
!!$
!!$                b11=gprimd(1,1)*real(ig1,kind=8)
!!$                b21=gprimd(2,1)*real(ig1,kind=8)
!!$                b31=gprimd(3,1)*real(ig1,kind=8)
!!$                b12=gprimd(1,2)*real(ig2,kind=8)
!!$                b22=gprimd(2,2)*real(ig2,kind=8)
!!$                b32=gprimd(3,2)*real(ig2,kind=8)
!!$                b13=gprimd(1,3)*real(ig3,kind=8)
!!$                b23=gprimd(2,3)*real(ig3,kind=8)
!!$                b33=gprimd(3,3)*real(ig3,kind=8)
!!$
!!$                !g2cart(ifft)=( &
!!$                !     &     (b11+b12+b13)**2&
!!$                !     &     +(b21+b22+b23)**2&
!!$                !     &     +(b31+b32+b33)**2&
!!$                !     &     )
!!$
!!$                ker = mu0_screening**2/16.0_dp/datan(1.0_dp)
!!$                ker = ker + pi*((b11+b12+b13)**2 &
!!$                    &     +(b21+b22+b23)**2 &
!!$                    &     +(b31+b32+b33)**2 &
!!$                    &     )
!!$
!!$                ker = ker + pi*(gd(1,1)*(p1/h1)**2+gd(3,3)*(p2/h3)**2+gd(2,2)*(p3/h2)**2)
!!$                !ker = ker + 2.0_dp*pi*(gd(1,3)*(p1/h1)*(p2/h3)+gd(2,3)*(p2/h3)*(p3/h2)+gd(1,2)*(p1/h1)*(p3/h2))
!!$
!!$              if (i3 == nker3/2) write(16,*) i1,i2,gs,ker


!              if (ker/=0._dp) then
                 !karray(i1,i2,i3)=1._dp/ker*fourISFx(i1-1)*fourISFy(i2-1)*fourISFz(j3-1)
                 karray(i1,i2,i3)=fourISFx(i1-1)*fourISFy(i2-1)*fourISFz(j3-1)
!              else
!                 karray(i1,i2,i3)=0._dp
!              end if
              !write(1717,*) i1,i2,i3,karray(i1,i2,i3)
           end do
        end do
     else
        do i2=1,nker2
           do i1=1,nker1
              karray(i1,i2,i3)=0._dp
           end do
        end do
     end if
  end do

  ! do i3=1,nker3
  !    do i2=1,nker2
  !       do i1=1,nker1
  !          write(1717,*) i1,i2,i3,1/karray(i1,i2,i3)
  !       end do
  !    end do
  ! end do

  call f_free(fourISFx)
  call f_free(fourISFy)
  call f_free(fourISFz)

END SUBROUTINE Periodic_Kernel


!> Calculate the Fourier transform
!! Suppose the output symmetric and real
subroutine fourtrans_isf(n,ftisf)
  use Poisson_Solver, only: dp
  implicit none
  integer, intent(in) :: n
  real(dp), dimension(0:n), intent(out) :: ftisf
  !local variables
  real(dp), parameter :: twopi=6.28318530717958647688_dp
  integer :: i,j
  real(dp) :: p,pointval,hval,q,htp

  !zero fourier component
  ftisf(0)=1._dp
  !non-zero components, use the previous calculated values for powers of two
  loop_points: do i=1,n
     !do nothing if the point can be divided by two
     if (2*(i/2) == i) then
        cycle loop_points
     end if
     p=real(i,dp)*twopi/real(2*n,dp)
     !calculate the values of the given point
     pointval=1._dp
     q=p
     loop_calc: do
        q=0.5_dp*q
        call fourtrans(q,htp)
        hval=htp
        if (abs(hval - 1._dp) <= 1.d-16) then
           exit loop_calc
        end if
        pointval=pointval*hval
     end do loop_calc
     ftisf(i)=pointval
     !calculate the other points on a dyadic grid until needed
     j=i
     q=p
     loop_dyadic: do
        j=2*j
        if (j > n) then
           exit loop_dyadic
        end if
        call fourtrans(q,htp)
        ftisf(j)=htp*ftisf(j/2)
        q=2._dp*q
     end do loop_dyadic
  end do loop_points

END SUBROUTINE fourtrans_isf


!> Transform the wavelet filters
subroutine fourtrans(p,htp)
  use Poisson_Solver, only: dp
  implicit none
  real(dp), intent(in) :: p
  real(dp), intent(out) :: htp
  !local variables
  integer :: i,j
  real(dp) :: cp,x
  !include the filters for a given scaling function
  include 'lazy_16.inc'

  htp=0._dp
  do j=m-3,1,-2
     x=real(j,dp)
     cp=cos(p*x)
     htp=htp+ch(j)*cp
  end do
  !this is the value divided by two
  htp=0.5_dp+htp

END SUBROUTINE fourtrans


!> Build the kernel of the Poisson operator with surfaces Boundary conditions
!! in an interpolating scaling functions basis.
!! @warning Beware of the fact that the nonperiodic direction is y!
!! SYNOPSIS
subroutine Surfaces_Kernel(iproc,nproc,mpi_comm,inplane_comm,n1,n2,n3,m3,nker1,nker2,nker3,&
     mesh,itype_scf,karray,mu0_screening)!,alpha)!,beta,gamma)!,n3pr2,n3pr1)
  use Poisson_Solver, only: dp
  use wrapper_mpi
  use dynamic_memory
  use f_utils, only: f_zero
  use module_fft_sg, only: p_index
  use box
  use numerics, only: twopi
  use at_domain, only: square_gu
  implicit none
  include 'perfdata.inc'

  !Arguments
  integer, intent(in) :: n1,n2,n3          !< Dimensions for the FFT
  integer, intent(in) :: m3                !< Actual dimension in non-periodic direction
  integer, intent(in) :: nker1,nker2,nker3 !< Dimensions of the kernel (nker3=n3/2+1) nker(1,2)=n(1,2)/2+1
  integer, intent(in) :: itype_scf         !< Order of the scaling function
  integer, intent(in) :: iproc             !< Process Id
  integer, intent(in) :: nproc             !< Number of processes
  integer, intent(in) :: mpi_comm,inplane_comm!n3pr1,n3pr2
  type(cell), intent(in) :: mesh
  real(dp), dimension(nker1,nker2,nker3/nproc), intent(out) :: karray !< Output array
  real(dp), intent(in) :: mu0_screening!,alpha!,beta,gamma

  !Local variables
  character(len=*), parameter :: subname='Surfaces_Kernel'
  !Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
  integer, parameter :: n_points = 2**6
  !Maximum number of points for FFT (should be same number in fft3d routine)
  !n(c) integer, parameter :: nfft_max=24000

  real(dp), dimension(:), allocatable :: kernel_scf
  real(dp), dimension(:), allocatable :: x_scf ,y_scf
  !FFT arrays
  real(dp), dimension(:,:,:), allocatable :: halfft_cache,kernel
  real(dp), dimension(:,:,:,:), allocatable :: kernel_mpi
  real(dp), dimension(:,:), allocatable :: cossinarr,btrig
  integer, dimension(:), allocatable :: after,now,before
  real(dp) :: h1,h2,h3
  real(dp), dimension(3) :: p
  real(dp) :: pi,dx,mu1,ponx,pony
  real(dp) :: a,b,c,d,feR,foR,foI,fR,cp,sp,pion,x,value !n(c) ,diff, feI
  integer :: n_scf,ncache,imu,ierr,ntrig
  integer :: n_range,n_cell,num_of_mus,shift,istart,iend,ireim,jreim,j2st,j2nd,nact2
  integer :: i,i1,i2,i3
  integer :: j2,ind1,ind2,jnd1,ic,inzee,nfft,ipolyord,jp2

  !metric for monoclinic lattices
!  real(dp), dimension(2,2) :: gus
!  real(dp) :: oneodetg

  !coefficients for the polynomial interpolation
  real(dp), dimension(9,8) :: cpol

  !assign the values of the coefficients
  cpol(:,:)=1._dp

  cpol(1,2)=.25_dp

  cpol(1,3)=1._dp/3._dp

  cpol(1,4)=7._dp/12._dp
  cpol(2,4)=8._dp/3._dp

  cpol(1,5)=19._dp/50._dp
  cpol(2,5)=3._dp/2._dp

  cpol(1,6)=41._dp/272._dp
  cpol(2,6)=27._dp/34._dp
  cpol(3,6)=27._dp/272._dp

  cpol(1,7)=751._dp/2989._dp
  cpol(2,7)=73._dp/61._dp
  cpol(3,7)=27._dp/61._dp

  cpol(1,8)=-989._dp/4540._dp
  cpol(2,8)=-1472._dp/1135._dp
  cpol(3,8)=232._dp/1135._dp
  cpol(4,8)=-2624._dp/1135._dp

  !renormalize values
  cpol(1,1)=.5_dp*cpol(1,1)
  cpol(1:2,2)=2._dp/3._dp*cpol(1:2,2)
  cpol(1:2,3)=3._dp/8._dp*cpol(1:2,3)
  cpol(1:3,4)=2._dp/15._dp*cpol(1:3,4)
  cpol(1:3,5)=25._dp/144._dp*cpol(1:3,5)
  cpol(1:4,6)=34._dp/105._dp*cpol(1:4,6)
  cpol(1:4,7)=2989._dp/17280._dp*cpol(1:4,7)
  cpol(1:5,8)=-454._dp/2835._dp*cpol(1:5,8)

  !assign the complete values
  cpol(2,1)=cpol(1,1)

  cpol(3,2)=cpol(1,2)

  cpol(3,3)=cpol(2,3)
  cpol(4,3)=cpol(1,3)

  cpol(4,4)=cpol(2,4)
  cpol(5,4)=cpol(1,4)

  cpol(4,5)=cpol(3,5)
  cpol(5,5)=cpol(2,5)
  cpol(6,5)=cpol(1,5)

  cpol(5,6)=cpol(3,6)
  cpol(6,6)=cpol(2,6)
  cpol(7,6)=cpol(1,6)

  cpol(5,7)=cpol(4,7)
  cpol(6,7)=cpol(3,7)
  cpol(7,7)=cpol(2,7)
  cpol(8,7)=cpol(1,7)

  cpol(6,8)=cpol(4,8)
  cpol(7,8)=cpol(3,8)
  cpol(8,8)=cpol(2,8)
  cpol(9,8)=cpol(1,8)

  h1=mesh%hgrids(1)
  h2=mesh%hgrids(3)
  h3=mesh%hgrids(2)

  !write(*,*) ' '
  !write(*,*) 'mu0_screening = ', mu0_screening

  !Number of integration points : 2*itype_scf*n_points
  n_scf=2*itype_scf*n_points
  !Allocations
  x_scf = f_malloc(0.to.n_scf,id='x_scf')
  y_scf = f_malloc(0.to.n_scf,id='y_scf')
  !Build the scaling function
  call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
  !Grid step for the integration
  dx = real(n_range,dp)/real(n_scf,dp)
  !Extend the range (no more calculations because fill in by 0._dp)
  n_cell = m3
  n_range = max(n_cell,n_range)

  ntrig=n3/2

  !Allocations
  ncache=ncache_optimal

  !the HalFFT must be performed only in the third dimension,
  !and nker3=n3/2+1, hence
  if (ncache <= (nker3-1)*4) ncache=(nker3-1)*4

  !enlarge the second dimension of the kernel to be compatible with nproc
  nact2=nker2
  enlarge_ydim: do
     if (nproc*(nact2/nproc) /= nact2) then
        nact2=nact2+1
     else
        exit enlarge_ydim
     end if
  end do enlarge_ydim

  !array for the MPI procedure
  kernel = f_malloc((/ nker1, nact2/nproc, nker3 /),id='kernel')
  kernel_mpi = f_malloc((/ nker1, nact2/nproc, nker3/nproc, nproc /),id='kernel_mpi')
  kernel_scf = f_malloc(n_range,id='kernel_scf')
  halfft_cache = f_malloc((/ 2, ncache/4, 2 /),id='halfft_cache')
  cossinarr = f_malloc((/ 2, n3/2-1 /),id='cossinarr')
  btrig = f_malloc((/ 2, ntrig /),id='btrig')
  after = f_malloc(7,id='after')
  now = f_malloc(7,id='now')
  before = f_malloc(7,id='before')

  !constants
  pi=4._dp*datan(1._dp)

!!$  !monoclinic cell
!!$  oneodetg = 1.0_dp/dsin(alpha)**2
!!$
!!$  call f_zero(gus)
!!$  gus(1,1)=oneodetg
!!$  gus(2,1)=-cos(alpha)*oneodetg
!!$  gus(1,2)=-cos(alpha)*oneodetg
!!$  gus(2,2)=oneodetg


  !arrays for the halFFT
  call ctrig_sg(n3/2,ntrig,btrig,after,before,now,1,ic)

  !build the phases for the HalFFT reconstruction
  pion=2._dp*pi/real(n3,dp)
  do i3=2,n3/2
     x=real(i3-1,dp)*pion
     cossinarr(1,i3-1)= dcos(x)
     cossinarr(2,i3-1)=-dsin(x)
  end do
  !kernel=0._dp
  !kernel_mpi=0._dp

  !calculate the limits of the FFT calculations
  !that can be performed in a row remaining inside the cache
  num_of_mus=ncache/(2*n3)

  !n(c) diff=0._dp
  !order of the polynomial to be used for integration (must be a power of two)
  ipolyord=8 !this part should be incorporated inside the numerical integration
  !here we have to choice the piece of the x-y grid to cover

  !let us now calculate the fraction of mu that will be considered
  j2st=iproc*(nact2/nproc)
  j2nd=min((iproc+1)*(nact2/nproc),nker2)!n2/2+1) !here we might write nker2 instead of n2/2+1

  do ind2=(n1/2+1)*j2st+1,(n1/2+1)*j2nd,num_of_mus
     istart=ind2
     iend=min(ind2+(num_of_mus-1),(n1/2+1)*j2nd)
     nfft=iend-istart+1
     shift=0

     !initialization of the interesting part of the cache array
     halfft_cache(:,:,:)=0._dp

     if (istart == 1) then
        !i2=1
        shift=1

        if (mu0_screening == 0.0_dp) then
           call calculates_green_opt_muzero(n_range,n_scf,ipolyord,x_scf,y_scf,&
                cpol(1,ipolyord),dx,kernel_scf)

           !copy of the first zero value
           halfft_cache(1,1,1)=0._dp

           do i3=1,m3

              value=0.5_dp*h3*kernel_scf(i3)
              !index in where to copy the value of the kernel
              call indices(ireim,num_of_mus,n3/2+i3,1,ind1)
              !index in where to copy the symmetric value
              call indices(jreim,num_of_mus,n3/2+2-i3,1,jnd1)
              halfft_cache(ireim,ind1,1) = value
              halfft_cache(jreim,jnd1,1) = value

           end do

        else
           mu1=abs(mu0_screening)*h3
           call calculates_green_opt(n_range,n_scf,itype_scf,ipolyord,x_scf,y_scf,&
                cpol(1,ipolyord),mu1,dx,kernel_scf)

           !copy of the first zero value
           halfft_cache(1,1,1)=0._dp

           do i3=1,m3

              value=-0.5_dp*h3/mu1*kernel_scf(i3)
              !index in where to copy the value of the kernel
              call indices(ireim,num_of_mus,n3/2+i3,1,ind1)
              !index in where to copy the symmetric value
              call indices(jreim,num_of_mus,n3/2+2-i3,1,jnd1)
              halfft_cache(ireim,ind1,1) = value
              halfft_cache(jreim,jnd1,1) = value

           end do


        end if

        ! !copy of the first zero value
        ! halfft_cache(1,1,1)=0._dp

        ! do i3=1,m3

        !    value=0.5_dp*h3*kernel_scf(i3)
        !    !index in where to copy the value of the kernel
        !    call indices(ireim,num_of_mus,n3/2+i3,1,ind1)
        !    !index in where to copy the symmetric value
        !    call indices(jreim,num_of_mus,n3/2+2-i3,1,jnd1)
        !    halfft_cache(ireim,ind1,1) = value
        !    halfft_cache(jreim,jnd1,1) = value

        ! end do

     end if

     loopimpulses : do imu=istart+shift,iend

        !here there is the value of mu associated to hgrid
        !note that we have multiplicated mu for hgrid to be comparable
        !with mu0ref

        !calculate the proper value of mu taking into account the periodic dimensions
        !corresponding value of i1 and i2
        i1=mod(imu,n1/2+1)
        if (i1==0) i1=n1/2+1
        i2=(imu-i1)/(n1/2+1)+1
        !with these conventions imu=i1+(i2-1)*(n1/2+1), with i1=0,...,n1/2 and (hopefully) i2=0,..,n2/2,-n2/2,-n2/2+1,...,-1
        ponx=real(i1-1,dp)/real(n1,dp)
        pony=real(p_index(i2,n2),dp)/real(n2,dp)!real(i2-1,dp)/real(n2,dp) !here we should put the value of p

        !print *,'here',i1,i2,i1-1,p_index(i2,n2)

        !acerioni --- adding the mu0_screening
        !mu1=2._dp*pi*sqrt((ponx/h1)**2+(pony/h2)**2)*h3
        !old:
        !mu1=2._dp*pi*sqrt((ponx/h1)**2+(pony/h2)**2)
!!$        !new:
!!$        mu1 = 2._dp*pi*sqrt(gus(1,1)*(ponx/h1)**2+gus(2,2)*(pony/h2)**2+2.0_dp*gus(1,2)*(ponx/h1)*(pony/h2))
!!$        mu1 = h3*sqrt(mu1**2+mu0_screening**2)
        !acerioni
        !new3:
        p=[ponx/h1,mu0_screening/twopi,pony/h2]
        mu1=h3*twopi*sqrt(square_gu(mesh%dom,p))

        call calculates_green_opt(n_range,n_scf,itype_scf,ipolyord,x_scf,y_scf,&
             cpol(1,ipolyord),mu1,dx,kernel_scf)

        !readjust the coefficient and define the final kernel

        !copy of the first zero value
        halfft_cache(1,imu-istart+1,1) = 0._dp
        do i3=1,m3
           value=-0.5_dp*h3/mu1*kernel_scf(i3)
           !write(80,*)mu1,i3,kernel_scf(i03)
           !index in where to copy the value of the kernel
           call indices(ireim,num_of_mus,n3/2+i3,imu-istart+1,ind1)
           !index in where to copy the symmetric value
           call indices(jreim,num_of_mus,n3/2+2-i3,imu-istart+1,jnd1)
           halfft_cache(ireim,ind1,1)=value
           halfft_cache(jreim,jnd1,1)=value
        end do

     end do loopimpulses

     !now perform the FFT of the array in cache
     inzee=1
     do i=1,ic
        call fftstp_sg(num_of_mus,nfft,n3/2,num_of_mus,n3/2,&
             halfft_cache(1,1,inzee),halfft_cache(1,1,3-inzee),&
             ntrig,btrig,after(i),now(i),before(i),1)
        inzee=3-inzee
     enddo
     !assign the values of the FFT array
     !and compare with the good results
     do imu=istart,iend

        !corresponding value of i1 and i2
        i1=mod(imu,n1/2+1)
        if (i1==0) i1=n1/2+1
        i2=(imu-i1)/(n1/2+1)+1

        j2=i2-j2st

        a=halfft_cache(1,imu-istart+1,inzee)
        b=halfft_cache(2,imu-istart+1,inzee)
        kernel(i1,j2,1)=a+b
        kernel(i1,j2,n3/2+1)=a-b

        do i3=2,n3/2
           ind1=imu-istart+1+num_of_mus*(i3-1)
           jnd1=imu-istart+1+num_of_mus*(n3/2+2-i3-1)
           cp=cossinarr(1,i3-1)
           sp=cossinarr(2,i3-1)
           a=halfft_cache(1,ind1,inzee)
           b=halfft_cache(2,ind1,inzee)
           c=halfft_cache(1,jnd1,inzee)
           d=halfft_cache(2,jnd1,inzee)
           feR=.5_dp*(a+c)
           !n(c) feI=.5_dp*(b-d)
           foR=.5_dp*(a-c)
           foI=.5_dp*(b+d)
           fR=feR+cp*foI-sp*foR
           kernel(i1,j2,i3)=fR
        end do
     end do

  end do

  call f_free(cossinarr)


  !give to each processor a slice of the third dimension
  if (nproc > 1) then
     call MPI_ALLTOALL(kernel,nker1*(nact2/nproc)*(nker3/nproc), &
          MPI_double_precision, &
          kernel_mpi,nker1*(nact2/nproc)*(nker3/nproc), &
          MPI_double_precision,mpi_comm,ierr)

!!     !Maximum difference
!!     max_diff = 0._dp
!!     i1_max = 1
!!     i2_max = 1
!!     i3_max = 1
!!     do i3=1,nker3/nproc
!!        do i2=1,nact2/nproc
!!           do i1=1,nker1
!!              factor=abs(kernel(i1,i2,i3+iproc*(nker3/nproc))&
!!                   -kernel_mpi(i1,i2,i3,iproc+1))
!!              if (max_diff < factor) then
!!                 max_diff = factor
!!                 i1_max = i1
!!                 i2_max = i2
!!                 i3_max = i3
!!              end if
!!           end do
!!        end do
!!     end do
!!     write(*,*) '------------------'
!!     print *,'iproc=',iproc,'difference post-mpi, at',i1_max,i2_max,i3_max
!!     write(unit=*,fmt="(1x,a,1pe12.4)") 'Max diff: ',max_diff,&
!!          'calculated',kernel(i1_max,i2_max,i3_max+iproc*(nker3/nproc)),'post-mpi',kernel_mpi(i1_max,i2_max,i3_max,iproc+1)

     do jp2=1,nproc
        do i3=1,nker3/nproc
           do i2=1,nact2/nproc
              j2=i2+(jp2-1)*(nact2/nproc)
              if (j2 <= nker2) then
                 do i1=1,nker1
                    karray(i1,j2,i3)=&
                         kernel_mpi(i1,i2,i3,jp2)
                 end do
              end if
           end do
        end do
     end do

     call mpi_barrier(mpi_comm,ierr)

     !!! PSolver n1-n2 plane mpi partitioning !!!
     !! Broadcast the result to all the other processors of the internal group
     !! actually the result can be scattered since only a section of the
     !! kernel should be multiplied
     !if (iproc < n3pr1*n3pr2) then
     if (inplane_comm /= MPI_COMM_NULL) then
        call MPI_Bcast(karray,nker1*nker2*nker3/nproc,MPI_double_precision,&
             0,inplane_comm,ierr)
        !exception should be added
     endif

  else
     karray(1:nker1,1:nker2,1:nker3)=kernel(1:nker1,1:nker2,1:nker3)
  endif


  !De-allocations
  call f_free(kernel)
  call f_free(kernel_mpi)
  call f_free(btrig)
  call f_free(after)
  call f_free(now)
  call f_free(before)
  call f_free(halfft_cache)
  call f_free(kernel_scf)
  call f_free(x_scf)
  call f_free(y_scf)

END SUBROUTINE Surfaces_Kernel


subroutine calculates_green_opt(n,n_scf,itype_scf,intorder,xval,yval,c,mu,hres,g_mu)
  use Poisson_Solver, only: dp
  use memory_profiling
  use dynamic_memory
  implicit none
  real(dp), parameter :: mu_max=0.2_dp
  integer, intent(in) :: n,n_scf,intorder,itype_scf
  real(dp), intent(in) :: hres,mu
  real(dp), dimension(0:n_scf), intent(in) :: xval,yval
  real(dp), dimension(intorder+1), intent(in) :: c
  real(dp), dimension(n), intent(out) :: g_mu
  !local variables
  character(len=*), parameter :: subname='calculates_green_opt'
  integer :: izero,ivalue,i,iend,ikern,n_iter,nrec
  real(dp) :: f,x,filter,gleft,gright,gltmp,grtmp,fl,fr,ratio,mu0 !n(c) x0, x1
  real(dp), dimension(:), allocatable :: green,green1

  !We calculate the number of iterations to go from mu0 to mu0_ref
  if (mu <= mu_max) then
     n_iter = 0
     mu0 = mu
  else
     n_iter=1
     loop_iter: do
        ratio=real(2**n_iter,dp)
        mu0=mu/ratio
        if (mu0 <= mu_max) then
           exit loop_iter
        end if
        n_iter=n_iter+1
     end do loop_iter
  end if

  !dimension needed for the correct calculation of the recursion
  nrec=2**n_iter*n

  green = f_malloc(-nrec.to.nrec,id='green')


  !initialization of the branching value
  ikern=0
  izero=0
  initialization: do
     if (xval(izero)>=real(ikern,dp) .or. izero==n_scf) exit initialization
     izero=izero+1
  end do initialization
  green=0._dp
  !now perform the interpolation in right direction
  ivalue=izero
  gright=0._dp
  loop_right: do
     if(ivalue >= n_scf-intorder-1) exit loop_right
     do i=1,intorder+1
        x=xval(ivalue)-real(ikern,dp)
        f=yval(ivalue)*exp(-mu0*x)
        filter=real(intorder,dp)*c(i)
        gright=gright+filter*f
        ivalue=ivalue+1
     end do
     ivalue=ivalue-1
  end do loop_right
  iend=n_scf-ivalue
  do i=1,iend
     x=xval(ivalue)-real(ikern,dp)
     f=yval(ivalue)*exp(-mu0*x)
     filter=real(intorder,dp)*c(i)
     gright=gright+filter*f
     ivalue=ivalue+1
  end do
  gright=hres*gright

  !the scaling function is symmetric, so the same for the other direction
  gleft=gright

  green(ikern)=gleft+gright

  !now the loop until the last value
  do ikern=1,nrec
     gltmp=0._dp
     grtmp=0._dp
     ivalue=izero
     !n(c) x0=xval(izero)
     loop_integration: do
        if (izero==n_scf)  exit loop_integration
        do i=1,intorder+1
           x=xval(ivalue)
           fl=yval(ivalue)*exp(mu0*x)
           fr=yval(ivalue)*exp(-mu0*x)
           filter=real(intorder,dp)*c(i)
           gltmp=gltmp+filter*fl
           grtmp=grtmp+filter*fr
           ivalue=ivalue+1
           if (xval(izero)>=real(ikern,dp) .or. izero==n_scf) then
              !n(c) x1=xval(izero)
              exit loop_integration
           end if
           izero=izero+1
        end do
        ivalue=ivalue-1
        izero=izero-1
     end do loop_integration
     gleft=exp(-mu0)*(gleft+hres*exp(-mu0*real(ikern-1,dp))*gltmp)
     if (izero == n_scf) then
        gright=0._dp
     else
        gright=exp(mu0)*(gright-hres*exp(mu0*real(ikern-1,dp))*grtmp)
     end if
     green(ikern)=gleft+gright
     green(-ikern)=gleft+gright
     if (abs(green(ikern)) <= 1.d-20) then
        nrec=ikern
        exit
     end if
     !print *,ikern,izero,n_scf,gltmp,grtmp,gleft,gright,x0,x1,green(ikern)
  end do
  !now we must calculate the recursion
  green1 = f_malloc(-nrec.to.nrec,id='green1')
  !Start the iteration to go from mu0 to mu
  call scf_recursion(itype_scf,n_iter,nrec,green(-nrec),green1(-nrec))

  do i=1,min(n,nrec)
     g_mu(i)=green(i-1)
  end do
  do i=min(n,nrec)+1,n
     g_mu(i)=0._dp
  end do

  call f_free(green)
  call f_free(green1)

END SUBROUTINE calculates_green_opt


subroutine calculates_green_opt_muzero(n,n_scf,intorder,xval,yval,c,hres,green)
  use Poisson_Solver, only: dp
  implicit none
  integer, intent(in) :: n,n_scf,intorder
  real(dp), intent(in) :: hres
  real(dp), dimension(0:n_scf), intent(in) :: xval,yval
  real(dp), dimension(intorder+1), intent(in) :: c
  real(dp), dimension(n), intent(out) :: green
  !local variables
  !n(c) character(len=*), parameter :: subname='calculates_green_opt_muzero'
  integer :: izero,ivalue,i,iend,ikern
  real(dp) :: x,y,filter,gl0,gl1,gr0,gr1,c0,c1

  !initialization of the branching value
  ikern=0
  izero=0
  initialization: do
     if (xval(izero)>=real(ikern,dp) .or. izero==n_scf) exit initialization
     izero=izero+1
  end do initialization
  green=0._dp
  !first case, ikern=0
  !now perform the interpolation in right direction
  ivalue=izero
  gr1=0._dp
  loop_right: do
     if(ivalue >= n_scf-intorder-1) exit loop_right
     do i=1,intorder+1
        x=xval(ivalue)
        y=yval(ivalue)
        filter=real(intorder,dp)*c(i)
        gr1=gr1+filter*x*y
        ivalue=ivalue+1
     end do
     ivalue=ivalue-1
  end do loop_right
  iend=n_scf-ivalue
  do i=1,iend
     x=xval(ivalue)
     y=yval(ivalue)
     filter=real(intorder,dp)*c(i)
     gr1=gr1+filter*x*y
     ivalue=ivalue+1
  end do
  gr1=hres*gr1
  !the scaling function is symmetric
  gl1=-gr1
  gl0=0.5_dp
  gr0=0.5_dp

  green(1)=2._dp*gr1

  !now the loop until the last value
  do ikern=1,n-1
     c0=0._dp
     c1=0._dp
     ivalue=izero
     loop_integration: do
        if (izero==n_scf)  exit loop_integration
        do i=1,intorder+1
           x=xval(ivalue)
           y=yval(ivalue)
           filter=real(intorder,dp)*c(i)
           c0=c0+filter*y
           c1=c1+filter*y*x
           ivalue=ivalue+1
           if (xval(izero)>=real(ikern,dp) .or. izero==n_scf) then
              exit loop_integration
           end if
           izero=izero+1
        end do
        ivalue=ivalue-1
        izero=izero-1
     end do loop_integration
     c0=hres*c0
     c1=hres*c1

     gl0=gl0+c0
     gl1=gl1+c1
     gr0=gr0-c0
     gr1=gr1-c1
     !general case
     green(ikern+1)=real(ikern,dp)*(gl0-gr0)+gr1-gl1
     !print *,ikern,izero,n_scf,gltmp,grtmp,gleft,gright,x0,x1,green(ikern)
  end do

END SUBROUTINE calculates_green_opt_muzero


subroutine indices(nimag,nelem,intrn,extrn,nindex)
  implicit none
  !arguments
  integer, intent(in) :: intrn,extrn,nelem
  integer, intent(out) :: nimag,nindex
  !local
  integer :: i
  !real or imaginary part
  nimag=2-mod(intrn,2)
  !actual index over half the length
  i=(intrn+1)/2
  !check
  if (2*(i-1)+nimag /= intrn) then
     print *,'error, index=',intrn,'nimag=',nimag,'i=',i
  end if
  !complete index to be assigned
  nindex=extrn+nelem*(i-1)
END SUBROUTINE indices


!> Build the kernel (karray) of a gaussian function
!! for interpolating scaling functions
!! @f$ K(j) = \sum_k \omega_k \int \int \phi(x) g_k(x'-x) \delta(x'- j) dx dx' @f$
!!
!! Do the parallel HalFFT of the symmetrized function and stores into
!! memory only 1/8 of the grid divided by the number of processes nproc
!!
!! MODIFICATION
!!    Different calculation of the gaussian times ISF integral, LG, Dec 2009
subroutine Free_Kernel(n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,&
     hx,hy,hz,itype_scf,iproc,nproc,karray,mu0_screening,n3pr2,n3pr1)
  use Poisson_Solver, only: dp, gp
  use memory_profiling
  use dynamic_memory
  use f_utils, only: f_zero
 implicit none
 !Arguments
 integer, intent(in) :: n01, n02, n03       !< Mesh dimensions of the density
 integer, intent(in) :: nfft1, nfft2, nfft3 !< Dimensions of the FFT grid (HalFFT in the third direction)
 integer, intent(in) :: n1k, n2k, n3k       !< Dimensions of the kernel FFT
 integer, intent(in) :: itype_scf           !< Order of the scaling function
 integer, intent(in) :: iproc, nproc
 integer, intent(in) :: n3pr2,n3pr1
 real(dp), intent(in) :: hx,hy,hz           !< Mesh steps
 real(dp), dimension(n1k,n2k,n3k/nproc), intent(out) :: karray
 real(dp), intent(in) :: mu0_screening
 !Local variables
 character(len=*), parameter :: subname='Free_Kernel'
 !Do not touch !!!!
 integer :: n_gauss
 !Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
 integer, parameter :: n_points = 2**6
 !Better p_gauss for calculation
 !(the support of the exponential should be inside [-n_range/2,n_range/2])
 !real(dp), dimension(89) :: p_gauss,w_gauss! = 0.0_dp, w_gauss = 0.0_dp
 !real(dp), dimension(n_gauss_Yukawa) :: p_gauss_Yukawa, w_gauss_Yukawa
 real(dp), dimension(:), allocatable :: fwork, p_gauss, w_gauss
 real(dp), dimension(:,:), allocatable :: kernel_scf, fftwork
 real(dp) :: ur_gauss, dr_gauss, acc_gauss, pgauss, a_range
 real(dp) :: factor, factor2 !mu0_screening = 1.0_dp !n(c) ,dx
 real(dp) :: a1,a2,a3,wg,k2,k3
 integer :: n_scf, nker2, nker3 !n(c) nker1
 integer :: i_gauss, n_range, n_cell
 integer :: i1, i2, i3, iproc1
 integer :: i03, iMin, iMax
 logical, parameter :: high_accuracy=.false.

 !!! PSolver n1-n2 plane mpi partitioning !!!

 if (n3pr1 >1) then
   iproc1=mod(iproc,n3pr2)
 else
   iproc1=iproc
 endif

 !substitute 'iproc1' for all 'iproc' here after
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!print *,'arguments',n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,&
!     hx,hy,hz,itype_scf,iproc,nproc,mu0_screening
!print '(3(1pe25.17))',hx,hy,hz
 !Number of integration points : 2*itype_scf*n_points
 n_scf=2*itype_scf*n_points
 !Set karray

 !here we must set the dimensions for the fft part, starting from the nfft
 !remember that actually nfft2 is associated to n03 and viceversa

 !Auxiliary dimensions only for building the FFT part
 !n(c) nker1=nfft1
 nker2=nfft2
 nker3=nfft3/2+1

 !adjusting the last two dimensions to be multiples of nproc
 do
    if(modulo(nker2,nproc) == 0) exit
    nker2=nker2+1
 end do
 do
    if(modulo(nker3,nproc) == 0) exit
    nker3=nker3+1
 end do

!!$ !Build the scaling function
!!$ call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
 n_range=2*itype_scf
 !Grid step for the integration
 !n(c) dx = real(n_range,dp)/real(n_scf,dp)
 !Extend the range (no more calculations because fill in by 0._dp)
 n_cell = max(n01,n02,n03)
 n_range = max(n_cell,n_range)

 !Lengths of the box (use box dimension)
 a1 = hx * real(n01,dp)
 a2 = hy * real(n02,dp)
 a3 = hz * real(n03,dp)

 if (mu0_screening == 0.0_dp) then

 !We divide the p_gauss by a_range**2 and a_gauss by a_range
    a_range = sqrt(a1*a1+a2*a2+a3*a3)
    factor = 1._dp/a_range
    !factor2 = factor*factor
    factor2 = 1._dp/(a1*a1+a2*a2+a3*a3)

    if (.not.high_accuracy) then
     n_gauss = 89
     p_gauss = f_malloc(1.to.n_gauss,id='p_gauss')
     w_gauss = f_malloc(1.to.n_gauss,id='w_gauss')

     !Initialization of the gaussian (Beylkin)
     call gequad(p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
     !In order to have a range from a_range=sqrt(a1*a1+a2*a2+a3*a3)
     !(biggest length in the cube)
     !do i_gauss=1,n_gauss
     !   p_gauss(i_gauss) = factor2*p_gauss(i_gauss)
     !end do
     !do i_gauss=1,n_gauss
     !   w_gauss(i_gauss) = factor*w_gauss(i_gauss)
     !end do
   else if (high_accuracy) then
     n_gauss = 173
     p_gauss = f_malloc(1.to.n_gauss,id='p_gauss')
     w_gauss = f_malloc(1.to.n_gauss,id='w_gauss')

     !Initialization of the gaussian (Beylkin)
     call HighAcc_gequad(p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
   end if
 else

    n_gauss = 90
    p_gauss = f_malloc(1.to.n_gauss,id='p_gauss')
    w_gauss = f_malloc(1.to.n_gauss,id='w_gauss')

    !Initialization of the gaussian (Mirone)
    call Yukawa_gequad(p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
    !In order to have a range from a_range=sqrt(a1*a1+a2*a2+a3*a3)
    !(biggest length in the cube)
    !We divide the p_gauss by a_range**2 and a_gauss by a_range
    do i_gauss=1,n_gauss
       p_gauss(i_gauss) = mu0_screening**2*p_gauss(i_gauss)
    end do
    do i_gauss=1,n_gauss
       ! we do not put the 1/(4\pi) factor here because
       ! it is already accounted for in 'scal'
       w_gauss(i_gauss) = mu0_screening*w_gauss(i_gauss)
    end do
    factor2=1.0_gp
    factor=1.0_gp
 end if

 call f_zero(karray)
!!$  do i3=1,n3k/nproc
!!$    !$omp parallel do default(shared) private(i2, i1)
!!$    do i2=1,n2k
!!$      do i1=1,n1k
!!$        karray(i1,i2,i3) = 0.0_dp
!!$      end do
!!$    end do
!!$    !$omp end parallel do
!!$  end do
 
 
 fwork = f_malloc(0.to.n_range,id='fwork')
 fftwork = f_malloc((/ 2, max(nfft1, nfft2, nfft3)*2 /),id='fftwork')

 kernel_scf = f_malloc((/ max(n1k, n2k, n3k), 3 /),id='kernel_scf')
 iMin = iproc1 * (nker3/nproc) + 1
 iMax = min((iproc1+1)*(nker3/nproc),nfft3/2+1)

  do i_gauss=n_gauss,1,-1
    !Gaussian
    pgauss = p_gauss(i_gauss) * factor2

!!$    if (i_gauss == 71 .or. .true.) then
!!$       print *,'pgauss,wgauss',pgauss,w_gauss(i_gauss)
!!$       !take the timings
!!$       call cpu_time(t0)
!!$       do itimes=1,ntimes
!!$          !this routine can be substituted by the wofz calculation
!!$          call gauss_conv_scf(itype_scf,pgauss,hx,dx,n_range,n_scf,x_scf,y_scf,&
!!$               kernel_scf,kern_1_scf)
!!$       end do
!!$       call cpu_time(t1)
!!$       told=real(t1-t0,dp)/real(ntimes,dp)
!!$
!!$       !take the timings
!!$       call cpu_time(t0)
!!$       do itimes=1,ntimes
!!$          call analytic_integral(sqrt(pgauss)*hx,n_range,itype_scf,fwork)
!!$       end do
!!$       call cpu_time(t1)
!!$       tnew=real(t1-t0,dp)/real(ntimes,dp)
!!$
!!$       !calculate maxdiff
!!$       maxdiff=0.0_dp
!!$       do i=0,n_range
!!$          !write(17,*)i,kernel_scf(i,1),kern_1_scf(i)
!!$          maxdiff=max(maxdiff,abs(kernel_scf(i,1)-(fwork(i))))
!!$       end do
!!$
!!$       do i=0,n_range
!!$          write(18,'(i4,3(1pe25.17))')i,kernel_scf(i,1),fwork(i)
!!$       end do
!!$
!!$       write(*,'(1x,a,i3,2(1pe12.5),1pe24.17)')'time,i_gauss',i_gauss,told,tnew,maxdiff
!!$       !stop
!!$    end if
!!$    !STOP

!    fwork = 0.0_dp

    call gauconv_ffts(itype_scf,pgauss,hx,hy,hz,nfft1,nfft2,nfft3,n1k,n2k,n3k,n_range,&
         fwork,fftwork,kernel_scf)

    !Add to the kernel (only the local part)
    wg=w_gauss(i_gauss) * factor
!tt=0.d0
!print *,'aaa',iMin,iMax,n2k,n3k,sum(kernel_scf(:,:))
    do i03 = iMin, iMax
       i3=i03-iproc1*(nker3/nproc)
       k3=kernel_scf(i03,3) * wg

       !$omp parallel do default(shared) private(i2, k2, i1)
       do i2=1,n2k
          k2=kernel_scf(i2,2)*k3
          do i1=1,n1k
             karray(i1,i2,i3) = karray(i1,i2,i3) + kernel_scf(i1,1) * k2
 !            tt=tt+kernel_scf(i1,1) * k2
          end do
       end do
       !$omp end parallel do
    end do
!print *,'igauss,tt',i_gauss,tt
!!$    do i3=1,nker3/nproc
!!$       if (iproc*(nker3/nproc)+i3  <= nfft3/2+1) then
!!$          i03=iproc*(nker3/nproc)+i3
!!$          do i2=1,n2k
!!$             do i1=1,n1k
!!$                karray(i1,i2,i3) = karray(i1,i2,i3) + w_gauss(i_gauss)* &
!!$                     kernel_scf(i1,1)*kernel_scf(i2,2)*kernel_scf(i03,3)
!!$             end do
!!$          end do
!!$       end if
!!$    end do

    !!!!ALAM: here the print statement can be added print *,'igauss',i_gauss
!!$
 end do
!!$stop
!!$

!!$ !De-allocations
 call f_free(kernel_scf)
 call f_free(fwork)
 call f_free(fftwork)
 call f_free(p_gauss)
 call f_free(w_gauss)

END SUBROUTINE Free_Kernel


subroutine gauconv_ffts(itype_scf,pgauss,hx,hy,hz,n1,n2,n3,nk1,nk2,nk3,n_range,fwork,fftwork,kffts)
  use Poisson_Solver, only: dp
  use dynamic_memory
  !n(c) use module_fft_sg
  implicit none
  integer, intent(in) :: itype_scf,n1,n2,n3,nk1,nk2,nk3,n_range
  real(dp), intent(in) :: pgauss,hx,hy,hz
  real(dp), dimension(0:n_range), intent(inout) :: fwork
  real(dp), dimension(2,max(n1,n2,n3)*2), intent(inout) :: fftwork
  real(dp), dimension(max(nk1,nk2,nk3),3), intent(out) :: kffts
  !local variables
  integer :: j,inzee,idir,n,nk
  real(dp) :: h
  integer, dimension(3) :: ndims,ndimsk
  real(dp), dimension(3) :: hgrids
  !!acerioni
  real(dp), dimension(:), allocatable :: x_scf, y_scf
  real(dp), dimension(-n_range:n_range) :: fwork_tmp
  integer :: n_points, n_scf, nrange
  real(dp) :: dx
  real(dp), dimension(-n_range:n_range) :: work
  !real(dp), dimension(-n_range:n_range) :: kernel_scf
  !real(dp), dimension(:), allocatable :: fwork2
  logical :: an_int
  !!acerioni

  ndims(1)=n1
  ndims(2)=n2
  ndims(3)=n3

  ndimsk(1)=nk1
  ndimsk(2)=nk2
  ndimsk(3)=nk3

  !beware of the inversion of second and third dimensions
  hgrids(1)=hx
  hgrids(2)=hz
  hgrids(3)=hy

  an_int = .false.

  if (.not. an_int) then
     !!acerioni: setup for gauss_conv_scf
     n_points = 16
     fwork_tmp = 0.0_dp

     !Number of integration points : 2*itype_scf*n_points
     n_scf=2*itype_scf*n_points

     !Other allocations
     x_scf = f_malloc(0.to.n_scf,id='x_scf')
     y_scf = f_malloc(0.to.n_scf,id='y_scf')
     !allocate(gaussian(0:n_scf),stat=i_stat)

     !Build the scaling function
     call scaling_function(itype_scf,n_scf,nrange,x_scf,y_scf)

     !Grid step for the integration
     dx = real(nrange,dp)/real(n_scf,dp)
  else
  end if

  !write(*,*) 'n_range = ', n_range

  !Allocations
  !allocate(work(-n_range:n_range), stat=i_stat)
  !allocate(fwork2(-n_range:n_range), stat=i_stat)
  !!acerioni


  if (hx == hy .and. hy == hz) then

     if(an_int) then
        !call analytic_integral(sqrt(pgauss)*hx,n_range,itype_scf,fwork)
     else
        call gauss_conv_scf(itype_scf, pgauss, hx, dx, n_range, n_scf, x_scf, y_scf, fwork_tmp, work)
        fwork(0:n_range) = fwork_tmp(0:n_range)
     end if

     !!open(unit=52, file = 'integral_comparison.plot', status = 'replace', position = 'rewind')
     !!do j= 0,n_range
     !!   write (52,*) j,fwork(j)
     !! end do
     !!close(52)

     do idir=1,3
        n=ndims(idir)
        nk=ndimsk(idir)
        !copy the values on the real part of the fftwork array

        ! fftwork=0.0_dp
        ! do j=0,min(n_range,n/2)
        !    fftwork(1,n/2+1+j)=fwork(j)
        !    fftwork(1,n/2+1-j)=fftwork(1,n/2+1+j)
        ! end do

        fftwork=0.0_dp
        do j=0,min(n_range,n/2)-1
           !fftwork(1,n/2+1+j)=kernel_scf(j)
           !fftwork(1,n/2+1-j)=kernel_scf(j)
           fftwork(1,n/2+1+j)=fwork(j)
           fftwork(1,n/2+1-j)=fwork(j)
        end do
        !old version of the loop, after Cray compiler bug.
        !do j=0,min(n_range,n/2)
        !   fftwork(1,n/2+1+j)=fwork(j)
        !   fftwork(1,n/2+1-j)=fftwork(1,n/2+1+j)
        !end do
        !calculate the fft
        call fft_1d_ctoc(1,1,n,fftwork,inzee)
        !copy the real part on the kfft array
        call dcopy(nk,fftwork(1,n*(inzee-1)+1),2,kffts(1,idir),1)
        !write(17+idir-1,'(1pe24.17)')kffts(:,idir)
     end do
  else
     do idir=1,3
        h=hgrids(idir)
        n=ndims(idir)
        nk=ndimsk(idir)
        if(an_int) then
           !call analytic_integral(sqrt(pgauss)*h,n_range,itype_scf,fwork)
        else
           call gauss_conv_scf(itype_scf, pgauss, h, dx, n_range, n_scf, x_scf, y_scf, fwork_tmp, work)
           fwork(0:n_range) = fwork_tmp(0:n_range)
        end if
        !copy the values on the real part of the fftwork array
        fftwork=0.0_dp
        do j=0,min(n_range,n/2)-1
           fftwork(1,n/2+1+j)=fwork(j)
           fftwork(1,n/2+1-j)=fwork(j)
        end do

        ! fftwork=0.0_dp
        ! do j=0,min(n_range,n/2)
           !fftwork(1,n/2+1+j)=fwork(j)
           !fftwork(1,n/2+1-j)=fftwork(1,n/2+1+j)
        ! end do

        !calculate the fft
        call fft_1d_ctoc(1,1,n,fftwork,inzee)
        !copy the real part on the kfft array
        call dcopy(nk,fftwork(1,n*(inzee-1)+1),2,kffts(1,idir),1)
        !write(17+idir-1,'(1pe24.17)')kffts(:,idir)
     end do
  end if

  if(.not. an_int) then
     call f_free(x_scf)
     call f_free(y_scf)
  else
  end if

END SUBROUTINE gauconv_ffts


!!!!> Here alpha correspondes to sqrtalpha in mathematica
!!!!! the final result is fwork(j+m)-fwork(j-m)
!!!subroutine analytic_integral(alpha,ntot,m,fwork)
!!!  use module_base
!!!  implicit none
!!!  integer, intent(in) :: ntot,m
!!!  real(dp), intent(in) :: alpha
!!!  real(dp), dimension(0:ntot), intent(inout) :: fwork
!!!  !integer, optional, intent(in) :: argument_nf
!!!  !local variables
!!!  integer :: nf
!!!  real(dp), parameter :: pi=3.1415926535897932384_dp
!!!  logical :: flag,flag1,flag2
!!!  integer :: j,q,jz
!!!  real(dp) :: if,r1,r2,res,ypm,ymm,erfcpm,erfcmm,factor,re,ro,factorend
!!!  !fourier transform, from mathematica
!!!  real(dp), dimension(:), pointer :: fISF
!!!!commented since inc files do not fulfill fortran norm
!!!!commenting it out, these include files do not fulfill fortran norm!!!
!!!!!$  include 'lazy_ISF_8_2048.inc'
!!!!!$  include 'lazy_ISF_14_2048.inc'
!!!!!$  include 'lazy_ISF_16_2048.inc'
!!!!!$  include 'lazy_ISF_20_2048.inc'
!!!!!$  include 'lazy_ISF_24_2048.inc'
!!!!!$  include 'lazy_ISF_30_2048.inc'
!!!!!$  include 'lazy_ISF_40_2048.inc'
!!!!!$  include 'lazy_ISF_50_2048.inc'
!!!!!$  include 'lazy_ISF_60_2048.inc'
!!!!!$  include 'lazy_ISF_100_2048.inc'
!!!!!$
!!!!!$  ! real(dp), dimension(0:nf) :: fISF = (/&
!!!!!$  !      1._dp,&
!!!!!$  !      0.99999999999999999284235222189_dp,&
!!!!!$  !      0.99999999999956013105290196791_dp,&
!!!!!$  !      0.99999999974047362436549957134_dp,&
!!!!!$  !      0.999999977723831277163111033987_dp,&
!!!!!$  !      0.999999348120402080917750802356_dp,&
!!!!!$  !      0.99999049783895018128985387648_dp,&
!!!!!$  !      0.99991555832357702478243846513_dp,&
!!!!!$  !      0.999483965311871663439701778468_dp,&
!!!!!$  !      0.997656434461567612067080906608_dp,&
!!!!!$  !      0.99166060711872037183053190362_dp,&
!!!!!$  !      0.975852499662821752376101449171_dp,&
!!!!!$  !      0.941478752030026285801304437187_dp,&
!!!!!$  !      0.878678036869166599071794039064_dp,&
!!!!!$  !      0.780993377358505853552754551915_dp,&
!!!!!$  !      0.65045481899429616671929018797_dp,&
!!!!!$  !      0.499741982655935831719850889234_dp,&
!!!!!$  !      0.349006378709142477173505869256_dp,&
!!!!!$  !      0.218427566876086979858296964531_dp,&
!!!!!$  !      0.120744342131771503808889607743_dp,&
!!!!!$  !      0.0580243447291221387538310922609_dp,&
!!!!!$  !      0.0237939962805936437124240819547_dp,&
!!!!!$  !      0.00813738655512823360783743450662_dp,&
!!!!!$  !      0.00225348982730563717615393684127_dp,&
!!!!!$  !      0.000485814732465831887324674087566_dp,&
!!!!!$  !      0.0000771922718944619737021087074916_dp,&
!!!!!$  !      8.34911217930991992166687474335e-6_dp,&
!!!!!$  !      5.4386155057958302499753464285e-7_dp,&
!!!!!$  !      1.73971967107293590801788194389e-8_dp,&
!!!!!$  !      1.8666847551067074314481563463e-10_dp,&
!!!!!$  !      2.8611022063937216744058802299e-13_dp,&
!!!!!$  !      4.12604249873737563657665492649e-18_dp,&
!!!!!$  !      0._dp,&
!!!!!$  !      3.02775647387618521868673172298e-18_dp,&
!!!!!$  !      1.53514570268556433704096392259e-13_dp,&
!!!!!$  !      7.27079116834250613985176986757e-11_dp,&
!!!!!$  !      4.86563325394873292266790060414e-9_dp,&
!!!!!$  !      1.0762051057772619178393861559e-7_dp,&
!!!!!$  !      1.1473008487467664271213583755e-6_dp,&
!!!!!$  !      7.20032699735536213944939155489e-6_dp,&
!!!!!$  !      0.0000299412827430274803875806741358_dp,&
!!!!!$  !      0.00008893448268856997565968618255_dp,&
!!!!!$  !      0.000198412101719649392215603405001_dp,&
!!!!!$  !      0.000344251643242443482702827827289_dp,&
!!!!!$  !      0.000476137217995145280942423549555_dp,&
!!!!!$  !      0.000534463964629735808538669640174_dp,&
!!!!!$  !      0.000493380569658768192772778551283_dp,&
!!!!!$  !      0.000378335841074046383032625565151_dp,&
!!!!!$  !      0.000242907366232915943662337043783_dp,&
!!!!!$  !      0.000131442744929435769666863922254_dp,&
!!!!!$  !      0.0000602917442687924479053419936816_dp,&
!!!!!$  !      0.0000235549783294245003536257246338_dp,&
!!!!!$  !      7.86058640769338837604478652921e-6_dp,&
!!!!!$  !      2.23813489015410093932264202671e-6_dp,&
!!!!!$  !      5.39326427012172337715252302537e-7_dp,&
!!!!!$  !      1.07974438850159507475448146609e-7_dp,&
!!!!!$  !      1.73882195410933428198534789337e-8_dp,&
!!!!!$  !      2.14124508643951633300134563969e-9_dp,&
!!!!!$  !      1.86666701805198449182670610067e-10_dp,&
!!!!!$  !      1.02097507219447295331839243873e-11_dp,&
!!!!!$  !      2.86110214266058470148547621716e-13_dp,&
!!!!!$  !      2.81016133980507633418466387683e-15_dp,&
!!!!!$  !      4.12604249873556074813981250712e-18_dp,&
!!!!!$  !      5.97401377952175312560963543305e-23_dp,&
!!!!!$  !      0._dp&
!!!!!$  !    /)
!!!!!$
!!!!!$
!!!!!$  !Only itype=8,14,16,20,24,30,40,50,60,100
!!!!!$  select case(m)
!!!!!$  case(8)
!!!!!$     fISF => fISF8
!!!!!$  case(14)
!!!!!$     fISF => fISF14
!!!!!$  case(16)
!!!!!$     fISF => fISF16
!!!!!$  case(20)
!!!!!$     fISF => fISF20
!!!!!$  case(24)
!!!!!$     fISF => fISF24
!!!!!$  case(30)
!!!!!$     fISF => fISF30
!!!!!$  case(40)
!!!!!$     fISF => fISF40
!!!!!$  case(50)
!!!!!$     fISF => fISF50
!!!!!$  case(60)
!!!!!$     fISF => fISF60
!!!!!$  case(100)
!!!!!$     fISF => fISF100
!!!!!$  case default
!!!!!$     print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
!!!!!$     stop
!!!!!$  end select
!!!
!!!
!!!  ! if(present(argument_nf)) then
!!!  !    nf=argument_nf
!!!  ! else
!!!  !    nf = 64 ! "default value"
!!!  ! endif
!!!
!!!  nf = 256
!!!
!!!  flag=.false.
!!!  factor=pi/real(2*m,dp)/alpha
!!!  factorend=sqrt(pi)/alpha/real(4*m,dp)
!!!  !fill work array
!!!  !the calculation for j=0 can be separated from the rest
!!!  !since it only contains error functions
!!!  loop_nonzero: do j=0,ntot
!!!     ypm=alpha*real(j+m,dp)
!!!     ymm=alpha*real(j-m,dp)
!!!     call derfcf(erfcpm,ypm)
!!!     call derfcf(erfcmm,ymm)
!!!     !assume nf even
!!!     res=0._dp
!!!     !reso=0._dp
!!!     do q=nf,2,-2
!!!        !the sign of q only influences the imaginary part
!!!        !so we can multiply by a factor of two
!!!
!!!        !call wofz_mod(alpha,m,-q,-j,r1,if,flag1)
!!!        call wofz_mod(alpha,m,q,-j-m,r1,if,flag1)
!!!        call wofz_mod(alpha,m,q,-j+m,r2,if,flag2)
!!!!!$        call wofz_mod(-xe,-y,r1,if,flag1)
!!!!!$        call wofz_mod(xe,-y,r2,if,flag2)
!!!        re=r1-r2
!!!        flag=flag1 .or. flag2 .or. flag
!!!        !if (flag) then
!!!        !   print *,'here',xe,y,q,j
!!!        !   stop
!!!        !end if
!!!        !call wofz_mod(alpha,m,-q+1,-j,r1,if,flag1)
!!!        call wofz_mod(alpha,m,q-1,-j-m,r1,if,flag1)
!!!        call wofz_mod(alpha,m,q-1,-j+m,r2,if,flag2)
!!!!!$        call wofz_mod(-xo,-y,r1,if,flag1)
!!!!!$        call wofz_mod(xo,-y,r2,if,flag2)
!!!        ro=r1-r2
!!!        flag=flag1 .or. flag2 .or. flag
!!!        !if (flag) then
!!!        !   print *,'there',xo,y
!!!        !   stop
!!!        !end if
!!!        !write(16,'(2(i4),6(1pe15.7))')j,q,re,ro,erfcmm-erfcpm
!!!        re=re*fISF(q)
!!!        ro=ro*fISF(q-1)
!!!        res=res+re-ro
!!!     end do
!!!     !q=0
!!!     !fwork(j)=derf(y)+rese-reso
!!!     fwork(j)=erfcmm-erfcpm+2.0_dp*res!e-reso
!!!     fwork(j)=factorend*fwork(j)
!!!     !exit form the loop if it is close to zero
!!!     if (abs(fwork(j)) < 1.e-25_dp) exit loop_nonzero
!!!     !write(17,'(i4,8(1pe15.7))')j,derf(y),erfcsgn,rese,reso,derf(y)+rese-reso,&
!!!     !     -erfcsgn+rese-reso,erfcsgn+rese-reso
!!!  end do loop_nonzero
!!!
!!!  !check flag
!!!  if (flag) then
!!!     write (*,*)'value of alpha',alpha
!!!     stop 'problem occurred in wofz'
!!!  end if
!!!
!!!  !put to zero the rest
!!!  do jz=j+1,ntot
!!!     fwork(jz)=0.0_dp
!!!  end do
!!!
!!!END SUBROUTINE analytic_integral


subroutine gauss_conv_scf(itype_scf,pgauss,hgrid,dx,n_range,n_scf,x_scf,y_scf,kernel_scf,work)
  use Poisson_Solver, only: dp
  implicit none
  integer, intent(in) :: n_range,n_scf,itype_scf
  real(dp), intent(in) :: pgauss,hgrid,dx
  real(dp), dimension(0:n_scf), intent(in) :: x_scf
  real(dp), dimension(0:n_scf), intent(in) :: y_scf
  real(dp), dimension(-n_range:n_range), intent(inout) :: work
  real(dp), dimension(-n_range:n_range), intent(inout) :: kernel_scf
  !local variables
  real(dp), parameter :: p0_ref = 1.0_dp
  real(dp), parameter :: mx_expo=634.71369555645470d0! = -log(tiny(1.0_dp)*1.d32)
  integer :: n_iter,i_kern,i
  real(dp) :: p0_cell,p0gauss,absci,kern,x_limit

  !Grid step for the integration
  !dx = real(n_range,dp)/real(n_scf,dp)

  !To have a correct integration
  p0_cell = p0_ref/(hgrid*hgrid)

  !We calculate the number of iterations to go from pgauss to p0_ref
  n_iter = nint((log(pgauss) - log(p0_cell))/log(4.0_dp))
  if (n_iter <= 0)then
     n_iter = 0
     p0gauss = pgauss
  else
     p0gauss = pgauss/4._dp**n_iter
  end if
  x_limit= sqrt(mx_expo/p0gauss)/hgrid !to avoid illegal operations
  !Stupid integration
  !Do the integration with the exponential centered in i_kern
  kernel_scf(:) = 0.0_dp
  do i_kern=0,n_range
     kern = 0.0_dp
     do i=0,n_scf
        absci = x_scf(i) - real(i_kern,dp)
        if (abs(absci) < x_limit) then
           absci = absci*absci*hgrid**2
           !print *,'done',x_limit,abs(absci),exp(-p0gauss*absci),y_scf(i)
           kern = kern + y_scf(i)*exp(-p0gauss*absci)
        end if
     end do
     kernel_scf(i_kern) = kern*dx
     kernel_scf(-i_kern) = kern*dx
     if (abs(kern) < 1.d-18) then
        !Too small not useful to calculate
        exit
     end if
  end do

  !Start the iteration to go from p0gauss to pgauss
  call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,work)

END SUBROUTINE gauss_conv_scf


subroutine inserthalf(n1,n3,lot,nfft,i1,zf,zw)
  use Poisson_Solver, only: dp
  implicit none
  integer, intent(in) :: n1,n3,lot,nfft,i1
  real(dp), dimension(n1/2+1,n3/2+1), intent(in) :: zf
  real(dp), dimension(2,lot,n3/2), intent(out) :: zw
  !local variables
  integer :: l1,l3,i01,i03r,i03i,i3

  i3=0
  do l3=1,n3,2
     i3=i3+1
     i03r=abs(l3-n3/2-1)+1
     i03i=abs(l3-n3/2)+1
     do l1=1,nfft
        i01=abs(l1-1+i1-n1/2-1)+1
        zw(1,l1,i3)=zf(i01,i03r)
        zw(2,l1,i3)=zf(i01,i03i)
     end do
  end do

END SUBROUTINE inserthalf


!> The conversion from d0 to dp type should be finished
subroutine gequad(p,w,urange,drange,acc)

   use Poisson_Solver, only: dp
   implicit none

   !Arguments
   real(dp), intent(out) :: urange,drange,acc
   real(dp), intent(out) :: p(*),w(*)

   ! range [10^(-9),1] and accuracy ~10^(-8);
   p(1)=4.96142640560223544d19
   p(2)=1.37454269147978052d19
   p(3)=7.58610013441204679d18
   p(4)=4.42040691347806996d18
   p(5)=2.61986077948367892d18
   p(6)=1.56320138155496681d18
   p(7)=9.35645215863028402d17
   p(8)=5.60962910452691703d17
   p(9)=3.3666225119686761d17
   p(10)=2.0218253197947866d17
   p(11)=1.21477756091902017d17
   p(12)=7.3012982513608503d16
   p(13)=4.38951893556421099d16
   p(14)=2.63949482512262325d16
   p(15)=1.58742054072786174d16
   p(16)=9.54806587737665531d15
   p(17)=5.74353712364571709d15
   p(18)=3.455214877389445d15
   p(19)=2.07871658520326804d15
   p(20)=1.25064667315629928d15
   p(21)=7.52469429541933745d14
   p(22)=4.5274603337253175d14
   p(23)=2.72414006900059548d14
   p(24)=1.63912168349216752d14
   p(25)=9.86275802590865738d13
   p(26)=5.93457701624974985d13
   p(27)=3.5709554322296296d13
   p(28)=2.14872890367310454d13
   p(29)=1.29294719957726902d13
   p(30)=7.78003375426361016d12
   p(31)=4.68148199759876704d12
   p(32)=2.8169955024829868d12
   p(33)=1.69507790481958464d12
   p(34)=1.01998486064607581d12
   p(35)=6.13759486539856459d11
   p(36)=3.69320183828682544d11
   p(37)=2.22232783898905102d11
   p(38)=1.33725247623668682d11
   p(39)=8.0467192739036288d10
   p(40)=4.84199582415144143d10
   p(41)=2.91360091170559564d10
   p(42)=1.75321747475309216d10
   p(43)=1.0549735552210995d10
   p(44)=6.34815321079006586d9
   p(45)=3.81991113733594231d9
   p(46)=2.29857747533101109d9
   p(47)=1.38313653595483694d9
   p(48)=8.32282908580025358d8
   p(49)=5.00814519374587467d8
   p(50)=3.01358090773319025d8
   p(51)=1.81337994217503535d8
   p(52)=1.09117589961086823d8
   p(53)=6.56599771718640323d7
   p(54)=3.95099693638497164d7
   p(55)=2.37745694710665991d7
   p(56)=1.43060135285912813d7
   p(57)=8.60844290313506695d6
   p(58)=5.18000974075383424d6
   p(59)=3.116998193057466d6
   p(60)=1.87560993870024029d6
   p(61)=1.12862197183979562d6
   p(62)=679132.441326077231_dp
   p(63)=408658.421279877969_dp
   p(64)=245904.473450669789_dp
   p(65)=147969.568088321005_dp
   p(66)=89038.612357311147_dp
   p(67)=53577.7362552358895_dp
   p(68)=32239.6513926914668_dp
   p(69)=19399.7580852362791_dp
   p(70)=11673.5323603058634_dp
   p(71)=7024.38438577707758_dp
   p(72)=4226.82479307685999_dp
   p(73)=2543.43254175354295_dp
   p(74)=1530.47486269122675_dp
   p(75)=920.941785160749482_dp
   p(76)=554.163803906291646_dp
   p(77)=333.46029740785694_dp
   p(78)=200.6550575335041_dp
   p(79)=120.741366914147284_dp
   p(80)=72.6544243200329916_dp
   p(81)=43.7187810415471025_dp
   p(82)=26.3071631447061043_dp
   p(83)=15.8299486353816329_dp
   p(84)=9.52493152341244004_dp
   p(85)=5.72200417067776041_dp
   p(86)=3.36242234070940928_dp
   p(87)=1.75371394604499472_dp
   p(88)=0.64705932650658966_dp
   p(89)=0.072765905943708247_dp

   w(1)=47.67445484528304247d10
   w(2)=11.37485774750442175d9
   w(3)=78.64340976880190239d8
   w(4)=46.27335788759590498d8
   w(5)=24.7380464827152951d8
   w(6)=13.62904116438987719d8
   w(7)=92.79560029045882433d8
   w(8)=52.15931216254660251d8
   w(9)=31.67018011061666244d8
   w(10)=1.29291036801493046d8
   w(11)=1.00139319988015862d8
   w(12)=7.75892350510188341d7
   w(13)=6.01333567950731271d7
   w(14)=4.66141178654796875d7
   w(15)=3.61398903394911448d7
   w(16)=2.80225846672956389d7
   w(17)=2.1730509180930247d7
   w(18)=1.68524482625876965d7
   w(19)=1.30701489345870338d7
   w(20)=1.01371784832269282d7
   w(21)=7.86264116300379329d6
   w(22)=6.09861667912273717d6
   w(23)=4.73045784039455683d6
   w(24)=3.66928949951594161d6
   w(25)=2.8462050836230259d6
   w(26)=2.20777394798527011d6
   w(27)=1.71256191589205524d6
   w(28)=1.32843556197737076d6
   w(29)=1.0304731275955989d6
   w(30)=799345.206572271448_dp
   w(31)=620059.354143595343_dp
   w(32)=480986.704107449333_dp
   w(33)=373107.167700228515_dp
   w(34)=289424.08337412132_dp
   w(35)=224510.248231581788_dp
   w(36)=174155.825690028966_dp
   w(37)=135095.256919654065_dp
   w(38)=104795.442776800312_dp
   w(39)=81291.4458222430418_dp
   w(40)=63059.0493649328682_dp
   w(41)=48915.9040455329689_dp
   w(42)=37944.8484018048756_dp
   w(43)=29434.4290473253969_dp
   w(44)=22832.7622054490044_dp
   w(45)=17711.743950151233_dp
   w(46)=13739.287867104177_dp
   w(47)=10657.7895710752585_dp
   w(48)=8267.42141053961834_dp
   w(49)=6413.17397520136448_dp
   w(50)=4974.80402838654277_dp
   w(51)=3859.03698188553047_dp
   w(52)=2993.51824493299154_dp
   w(53)=2322.1211966811754_dp
   w(54)=1801.30750964719641_dp
   w(55)=1397.30379659817038_dp
   w(56)=1083.91149143250697_dp
   w(57)=840.807939169209188_dp
   w(58)=652.228524366749422_dp
   w(59)=505.944376983506128_dp
   w(60)=392.469362317941064_dp
   w(61)=304.444930257324312_dp
   w(62)=236.162932842453601_dp
   w(63)=183.195466078603525_dp
   w(64)=142.107732186551471_dp
   w(65)=110.23530215723992_dp
   w(66)=85.5113346705382257_dp
   w(67)=66.3325469806696621_dp
   w(68)=51.4552463353841373_dp
   w(69)=39.9146798429449273_dp
   w(70)=30.9624728409162095_dp
   w(71)=24.018098812215013_dp
   w(72)=18.6312338024296588_dp
   w(73)=14.4525541233150501_dp
   w(74)=11.2110836519105938_dp
   w(75)=8.69662175848497178_dp
   w(76)=6.74611236165731961_dp
   w(77)=5.23307018057529994_dp
   w(78)=4.05937850501539556_dp
   w(79)=3.14892659076635714_dp
   w(80)=2.44267408211071604_dp
   w(81)=1.89482240522855261_dp
   w(82)=1.46984505907050079_dp
   w(83)=1.14019261330527007_dp
   w(84)=0.884791217422925293_dp
   w(85)=0.692686387080616483_dp
   w(86)=0.585244576897023282_dp
   w(87)=0.576182522545327589_dp
   w(88)=0.596688817388997178_dp
   w(89)=0.607879901151108771_dp

   urange = 1._dp
   drange=1d-08
   acc   =1d-08

END SUBROUTINE gequad

!> The conversion from d0 to dp type should be finished
subroutine HighAcc_gequad(p,w,urange,drange,acc)

   use Poisson_Solver, only: dp
   implicit none

   !Arguments
   real(dp), intent(out) :: urange,drange,acc
   real(dp), intent(out) :: p(*),w(*)

   ! range [10^(-9),1] and accuracy ~10^(-14);
    p(1)=  0.95133594396161155d-4
    p(2)=  0.83736214493632534d-3
    p(3)=  0.22225444433217846d-2
    p(4)=  0.40609515900839489d-2
    p(5)=  0.61346196197642824d-2
    p(6)=  0.84564463671216172d-2
    p(7)=  0.11396643432652089d-1
    p(8)=  0.15332071501619455d-1
    p(9)=  0.20625729634019541d-1
    p(10)= 0.27747108469207471d-1
    p(11)= 0.37327262694132474d-1
    p(12)= 0.50215125723209117d-1
    p(13)= 0.67552739456408911d-1
    p(14)= 0.90876454899651224d-1
    p(15)= 0.12225307399202322_dp
    p(16)= 0.16446299668052608_dp
    p(17)= 0.22124660259179399_dp
    p(18)= 0.29763570010523072_dp
    p(19)= 0.40039941377347954_dp
    p(20)= 0.53864402184772930_dp
    p(21)= 0.72461989776150504_dp
    p(22)= 0.97480706168556297_dp
    p(23)= 0.13113755369505424d1
    p(24)= 0.17641499189991365d1
    p(25)= 0.23732522446939215d1
    p(26)= 0.31926573565471994d1
    p(27)= 0.42949758160367129d1
    p(28)= 0.57778881979022039d1
    p(29)= 0.77728009323841558d1
    p(30)= 0.10456490722061620d2
    p(31)= 0.14066769388756272d2
    p(32)= 0.18923557271366974d2
    p(33)= 0.25457232567482787d2
    p(34)= 0.34246768760304711d2
    p(35)= 0.46071039631381019d2
    p(36)= 0.61977838188827008d2
    p(37)= 0.83376725537225923d2
    p(38)= 0.11216393737597498d3
    p(39)= 0.15089041655953446d3
    p(40)= 0.20298786171522153d3
    p(41)= 0.27307282293478642d3
    p(42)= 0.36735579159993426d3
    p(43)= 0.49419153532626041d3
    p(44)= 0.66481944526980817d3
    p(45)= 0.89435950074916093d3
    p(46)= 0.12031521073450110d4
    p(47)= 0.16185605365584393d4
    p(48)= 0.21773956879694165d4
    p(49)= 0.29291780411679165d4
    p(50)= 0.39405258512574896d4
    p(51)= 0.53010584424009157d4
    p(52)= 0.71313377123977598d4
    p(53)= 0.95935515748118687d4
    p(54)= 0.12905886038543824d5
    p(55)= 0.17361859488742211d5
    p(56)= 0.23356332452229686d5
    p(57)= 0.31420497670356286d5
    p(58)= 0.42268951080915365d5
    p(59)= 0.56863014845445825d5
    p(60)= 0.76495923712980191d5
    p(61)= 0.10290742340354310d6
    p(62)= 0.13843793600414865d6
    p(63)= 0.18623595355151905d6
    p(64)= 0.25053703772501458d6
    p(65)= 0.33703914885940461d6
    p(66)= 0.45340756358966930d6
    p(67)= 0.60995412377475039d6
    p(68)= 0.82055100749600702d6
    p(69)= 0.11038599948073409d7
    p(70)= 0.14849861580872252d7
    p(71)= 0.19977025166995658d7
    p(72)= 0.26874427909605918d7
    p(73)= 0.36153274545693672d7
    p(74)= 0.48635798491140092d7
    p(75)= 0.65428123028831985d7
    p(76)= 0.88018279043078329d7
    p(77)= 0.11840806501954885d8
    p(78)= 0.15929043392011497d8
    p(79)= 0.21428812584910858d8
    p(80)= 0.28827469264697671d8
    p(81)= 0.38780636160504386d8
    p(82)= 0.52170300736567020d8
    p(83)= 0.70182971410764083d8
    p(84)= 0.94414818517454311d8
    p(85)= 0.12701311694985677d9
    p(86)= 0.17086652424519673d9
    p(87)= 0.22986105536768723d9
    p(88)= 0.30922443707537574d9
    p(89)= 0.41598935640329200d9
    p(90)= 0.55961665344915390d9
    p(91)= 0.75283368191280985d9
    p(92)= 0.10127621276622870d10
    p(93)= 0.13624352255612867d10
    p(94)= 0.18328388208343513d10
    p(95)= 0.24656571410751729d10
    p(96)= 0.33169665920581636d10
    p(97)= 0.44622048984605818d10
    p(98)= 0.60028559236984396d10
    p(99)= 0.80754425358440866d10
    p(100)=0.10863624411218838d11
    p(101)=0.14614472807427057d11
    p(102)=0.19660364474535069d11
    p(103)=0.26448434806018219d11
    p(104)=0.35580200183687508d11
    p(105)=0.47864860601247337d11
    p(106)=0.64391005911968369d11
    p(107)=0.86623079860117981d11
    p(108)=0.11653114993591792d12
    p(109)=0.15676548244781943d12
    p(110)=0.21089139256424255d12
    p(111)=0.28370518026816327d12
    p(112)=0.38165914849498779d12
    p(113)=0.51343336590554791d12
    p(114)=0.69070484033885205d12
    p(115)=0.92918226229045325d12
    p(116)=0.12499980109182966d13
    p(117)=0.16815807734514434d13
    p(118)=0.22621747178334521d13
    p(119)=0.30432284519412812d13
    p(120)=0.40939540777712070d13
    p(121)=0.55074603354892686d13
    p(122)=0.74090033182542949d13
    p(123)=0.99670858846100293d13
    p(124)=0.13408389329024725d14
    p(125)=0.18037860461933641d14
    p(126)=0.24265734090812223d14
    p(127)=0.32643885465719379d14
    p(128)=0.43914734015917680d14
    p(129)=0.59077031921155641d14
    p(130)=0.79474367289759625d14
    p(131)=0.10691422454224173d15
    p(132)=0.14382815238772038d15
    p(133)=0.19348723247851362d15
    p(134)=0.26029194222890300d15
    p(135)=0.35016209762998294d15
    p(136)=0.47106143035652962d15
    p(137)=0.63370328391172888d15
    p(138)=0.85249996319286688d15
    p(139)=0.11468398629050685d16
    p(140)=0.15428055459640645d16
    p(141)=0.20754850172613612d16
    p(142)=0.27920810034323050d16
    p(143)=0.37560937635744200d16
    p(144)=0.50529480854673600d16
    p(145)=0.67975630965426680d16
    p(146)=0.91445356789579920d16
    p(147)=0.12301839879392896d17
    p(148)=0.16549256269670064d17
    p(149)=0.22263164353000016d17
    p(150)=0.29949894963985972d17
    p(151)=0.40290598143695328d17
    p(152)=0.54201602400568144d17
    p(153)=0.72915614017732912d17
    p(154)=0.98090951782034928d17
    p(155)=0.13195849683395963d18
    p(156)=0.17751937941606634d18
    p(157)=0.23881092028441696d18
    p(158)=0.32126439284930816d18
    p(159)=0.43218630869104928d18
    p(160)=0.58140587496608563d18
    p(161)=0.78214599733355878d18
    p(162)=0.10521950112399391d19
    p(163)=0.14154829730671122d19
    p(164)=0.19042021922169725d19
    p(165)=0.25616599124375373d19
    p(166)=0.34461159291859302d19
    p(167)=0.46359452087019561d19
    p(168)=0.62365829878403523d19
    p(169)=0.83898677860153580d19
    p(170)=0.11286610248602075d20
    p(171)=0.15183501594171152d20
    p(172)=0.20425859986505535d20
    p(173)=0.27478230472770691d20

    w(1)=   0.21950860876721736d-1
    w(2)=   0.21219095663737503d-1
    w(3)=   0.19744429381095819d-1
    w(4)=   0.17584620938882533d-1
    w(5)=   0.15543228346269854d-1
    w(6)=   0.15699043444815337d-1
    w(7)=   0.17873729612986982d-1
    w(8)=   0.20719864890700004d-1
    w(9)=   0.24031984638887470d-1
    w(10)=  0.27873649959237631d-1
    w(11)=  0.32329429984634325d-1
    w(12)=  0.37497494754359741d-1
    w(13)=  0.43491707509891335d-1
    w(14)=  0.50444133255223095d-1
    w(15)=  0.58507948424236150d-1
    w(16)=  0.67860815676888606d-1
    w(17)=  0.78708798178009937d-1
    w(18)=  0.91290899598435685d-1
    w(19)=  0.10588432986415428_dp
    w(20)=  0.12281061266892306_dp
    w(21)=  0.14244266931156566_dp
    w(22)=  0.16521303492966022_dp
    w(23)=  0.19162338815039942_dp
    w(24)=  0.22225560411666404_dp
    w(25)=  0.25778457440952840_dp
    w(26)=  0.29899307631685168_dp
    w(27)=  0.34678901904888998_dp
    w(28)=  0.40222544687103023_dp
    w(29)=  0.46652373986441475_dp
    w(30)=  0.54110052347550464_dp
    w(31)=  0.62759887972808004_dp
    w(32)=  0.72792454774583493_dp
    w(33)=  0.84428791115842061_dp
    w(34)=  0.97925269746107102_dp
    w(35)=  0.11357924622763562d1
    w(36)=  0.13173561029839187d1
    w(37)=  0.15279438451200911d1
    w(38)=  0.17721953756864086d1
    w(39)=  0.20554920651272290d1
    w(40)=  0.23840755301398366d1
    w(41)=  0.27651851495032318d1
    w(42)=  0.32072175626854924d1
    w(43)=  0.37199116653170043d1
    w(44)=  0.43145631773650566d1
    w(45)=  0.50042735113948140d1
    w(46)=  0.58042384240023983d1
    w(47)=  0.67320828100131882d1
    w(48)=  0.78082490156604107d1
    w(49)=  0.90564472260915654d1
    w(50)=  0.10504177850178792d2
    w(51)=  0.12183337411861066d2
    w(52)=  0.14130921297065543d2
    w(53)=  0.16389838839188403d2
    w(54)=  0.19009858701169911d2
    w(55)=  0.22048705382898401d2
    w(56)=  0.25573331012288250d2
    w(57)=  0.29661390440244212d2
    w(58)=  0.34402952138924391d2
    w(59)=  0.39902482597958780d2
    w(60)=  0.46281147939014723d2
    w(61)=  0.53679483458068560d2
    w(62)=  0.62260489910967401d2
    w(63)=  0.72213224759915548d2
    w(64)=  0.83756967503520372d2
    w(65)=  0.97146050861305866d2
    w(66)=  0.11267546425377405d3
    w(67)=  0.13068735303434275d3
    w(68)=  0.15157855666479423d3
    w(69)=  0.17580935191597666d3
    w(70)=  0.20391359372466630d3
    w(71)=  0.23651047713081883d3
    w(72)=  0.27431817943523964d3
    w(73)=  0.31816968314280149d3
    w(74)=  0.36903112830364876d3
    w(75)=  0.42802309859277733d3
    w(76)=  0.49644531010462623d3
    w(77)=  0.57580524680831729d3
    w(78)=  0.66785137356240239d3
    w(79)=  0.77461165844094342d3
    w(80)=  0.89843825309819852d3
    w(81)=  0.10420593155733814d4
    w(82)=  0.12086391173001150d4
    w(83)=  0.14018477586030970d4
    w(84)=  0.16259420286597663d4
    w(85)=  0.18858591914407307d4
    w(86)=  0.21873257639284129d4
    w(87)=  0.25369836832249584d4
    w(88)=  0.29425366422741617d4
    w(89)=  0.34129198182779392d4
    w(90)=  0.39584967332783622d4
    w(91)=  0.45912875841548321d4
    w(92)=  0.53252340726214452d4
    w(93)=  0.61765065699817542d4
    w(94)=  0.71638603841216891d4
    w(95)=  0.83090489780439875d4
    w(96)=  0.96373032440103052d4
    w(97)=  0.11177887392700828d5
    w(98)=  0.12964743704786575d5
    w(99)=  0.15037240350136301d5
    w(100)= 0.17441038750676147d5
    w(101)= 0.20229099596711181d5
    w(102)= 0.23462849681346550d5
    w(103)= 0.27213535260806380d5
    w(104)= 0.31563791757995041d5
    w(105)= 0.36609464393144051d5
    w(106)= 0.42461719853838506d5
    w(107)= 0.49249495528908665d5
    w(108)= 0.57122340267917702d5
    w(109)= 0.66253709254107307d5
    w(110)= 0.76844785583709789d5
    w(111)= 0.89128912748989474d5
    w(112)= 0.10337673568186787d6
    w(113)= 0.11990216362602370d6
    w(114)= 0.13906928621197681d6
    w(115)= 0.16130039511072749d6
    w(116)= 0.18708528799967427d6
    w(117)= 0.21699205982660220d6
    w(118)= 0.25167961912575946d6
    w(119)= 0.29191220514660433d6
    w(120)= 0.33857622563778306d6
    w(121)= 0.39269978625783697d6
    w(122)= 0.45547534188632516d6
    w(123)= 0.52828596889089234d6
    w(124)= 0.61273583718311926d6
    w(125)= 0.71068555342615815d6
    w(126)= 0.82429315407857893d6
    w(127)= 0.95606165143667767d6
    w(128)= 0.11088941802138176d7
    w(129)= 0.12861579596507163d7
    w(130)= 0.14917584803753824d7
    w(131)= 0.17302255505040844d7
    w(132)= 0.20068130967580504d7
    w(133)= 0.23276149194226521d7
    w(134)= 0.26996989514724654d7
    w(135)= 0.31312629798701671d7
    w(136)= 0.36318152599042701d7
    w(137)= 0.42123840018765945d7
    w(138)= 0.48857603455672013d7
    w(139)= 0.56667801757110348d7
    w(140)= 0.65726509874694990d7
    w(141)= 0.76233310034236396d7
    w(142)= 0.88419688948269971d7
    w(143)= 0.10255413795095254d8
    w(144)= 0.11894806842191380d8
    w(145)= 0.13796267282818770d8
    w(146)= 0.16001688254729975d8
    w(147)= 0.18559659779891673d8
    w(148)= 0.21526539303970236d8
    w(149)= 0.24967693368357640d8
    w(150)= 0.28958937771356728d8
    w(151)= 0.33588207948280476d8
    w(152)= 0.38957496372425713d8
    w(153)= 0.45185099661896639d8
    w(154)= 0.52408224900732599d8
    w(155)= 0.60786012597021684d8
    w(156)= 0.70503042880080983d8
    w(157)= 0.81773402185524419d8
    w(158)= 0.94845400025772005d8
    w(159)= 0.11000703976630102d9
    w(160)= 0.12759236394022685d9
    w(161)= 0.14798881390172765d9
    w(162)= 0.17164576596684095d9
    w(163)= 0.19908443211058137d9
    w(164)= 0.23090934335338539d9
    w(165)= 0.26782166883985958d9
    w(166)= 0.31063466405685753d9
    w(167)= 0.36029158854735017d9
    w(168)= 0.41788648788472962d9
    w(169)= 0.48468829777769488d9
    w(170)= 0.56216880136942244d9
    w(171)= 0.65203505568868458d9
    w(172)= 0.75626700167512512d9
    w(173)= 0.87716108640598845d9

   urange = 1._dp
   drange=1d-09
   acc   =1d-14

END SUBROUTINE HighAcc_gequad


!> Build the kernel of the Poisson operator with wires Boundary conditions
!! in an interpolating scaling functions basis.
!! The periodic direction is z
subroutine Wires_Kernel(iproc,nproc,n01,n02,n03,n1,n2,n3,nker1,nker2,nker3,h1,h2,h3,itype_scf,karray, &
                        mu0_screening)
  use Poisson_Solver, only: dp
  use memory_profiling
  use dynamic_memory
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc              !< Number of processes
  integer, intent(in) :: n01,n02,n03
  integer, intent(in) :: n1,n2,n3           !< Dimensions for the FFT
  integer, intent(in) :: nker1,nker2,nker3  !< Dimensions of the kernel nker(1,2,3)=n(1,2,3)/2+1
  integer, intent(in) :: itype_scf          !< Order of the scaling function
  real(dp), intent(in) :: h1,h2,h3          !< Mesh steps in the three dimensions
  real(dp), dimension(nker1,nker2,nker3/nproc), intent(out) :: karray !< Output array
  real(dp), intent(in) :: mu0_screening
  !Local variables
  character(len=*), parameter :: subname='Wires_Kernel'
  real(dp), parameter :: pi=3.14159265358979323846_dp
  integer, parameter :: n_gauss = 144
  integer :: i1, i2, i3, n_range, n_cell, k,i3s,i3e,nn
  real(dp) :: mu, t0, t1
  !real(dp), dimension(:), allocatable :: fourISFx,fourISFy,fourISFz
  real(dp), dimension(:), allocatable :: fwork
  real(dp), dimension(:,:), allocatable :: kernel_scf,fftwork
  real(dp), dimension(:), pointer :: alpha, w

  !load alpha(:) and w(:) coefficients from an .inc file
  include 'gaussfit_wires.inc'

  !write(*,*) "Entering the Wires_Kernel subroutine..."
  call cpu_time(t0)

  nn=n03

  w2 = -w2 !because we actually need -K0(mu*r)

  n_range = 2*itype_scf

  n_cell = max(n01,n02)
  n_range = max(n_cell,n_range)

  kernel_scf = f_malloc((/ max(nker1, nker2, nker3), 3 /),id='kernel_scf')
  fwork = f_malloc(0.to.n_range,id='fwork')
  fftwork = f_malloc((/ 2, max(n1, n2, n3)*2 /),id='fftwork')

  ! initialization
  karray = 0.0_dp

  ! case i2 = 1 (namely mu = 0)

  !parallel region for the nker3
  i3s=iproc*(nker3/nproc)+1
  i3e=min((iproc+1)*(nker3/nproc),nker3)

  if (mu0_screening == 0.0_dp) then ! .and. i3s==1) then
     !loads the coefficients alpha(:) and w(:) of the Gaussian fit for log(x):
     alpha => p1
     w => w1

     ! we introduce kind of an 'effective' mu:
     mu = 1.0_dp/sqrt((h1*(n01+2*itype_scf))**2+(h2*(n02+2*itype_scf))**2)
     !mu = 1.0_dp/sqrt((h1*(n_range))**2+(h2*(n_range))**2)
     !mu = 2.0_dp*1.0_dp/(2.0d0*h1*(n01/2))
     !mu = 2.0d0/max(h1*n01,h2*n02)

     ! because of the scaling properties of the log function we have to add the following:
     if (i3s == 1) karray(1,1,1) = karray(1,1,1) - (n1*n3)*log(mu)
  else
     alpha => p2
     w => w2
     mu = mu0_screening
  end if

  do k = 1, n_gauss
     fwork = 0.0_dp
     call gauconv_ffts(itype_scf,alpha(k)*mu**2,h1,h2,h3,n1,n2,n3,nker1,nker2,nker3,n_range,fwork,fftwork,kernel_scf)
     do i3 = i3s,i3e !1, nker3
        do i1 = 1, nker1
           karray(i1,1,i3-i3s+1) = karray(i1,1,i3-i3s+1) + w(k)*kernel_scf(i1,1)*kernel_scf(i3,3)
        end do
     end do
  end do

  ! case i2 != 1 (namely mu != 0)
  ! loads the coefficients alpha(:) and w(:) of the Gaussian fit for -BesselK0
  alpha => p2
  w => w2

  do i2 = 2, nker2
     if (i2 <= n2/2+1) then
        mu = 2.0_dp*pi/real(n2,dp)*real(i2-1,dp)/h3
     else
        mu = 2.0_dp*pi/real(n2,dp)*real(-i2+n2+1,dp)/h3
     end if

     mu = sqrt(mu**2 + mu0_screening**2)

     do k = 1, n_gauss
        fwork = 0.0_dp
        call gauconv_ffts(itype_scf,alpha(k)*mu**2,h1,h2,h3,n1,n2,n3,nker1,nker2,nker3,n_range,fwork,fftwork,kernel_scf)
        do i3 = i3s,i3e !1, nker3
           do i1 = 1, nker1
             karray(i1,i2,i3-i3s+1) = karray(i1,i2,i3-i3s+1) + w(k)*kernel_scf(i1,1)*kernel_scf(i3,3)
           end do
        end do
     end do
  end do



  call cpu_time(t1)
  !write(*,*) "Exiting the Wires_Kernel subroutine..."
  !write(*,*) "Elapsed time = ", t1-t0


  call f_free(kernel_scf)
  call f_free(fwork)
  call f_free(fftwork)


END SUBROUTINE Wires_Kernel


!> The conversion from d0 to dp type should be finished
subroutine Yukawa_gequad(p,w,urange,drange,acc)

  use Poisson_Solver, only: dp
  implicit none

!Arguments
  real(dp), intent(out) :: urange,drange,acc
  real(dp), dimension(1:90), intent(out) :: p, w
!
! range [10^(-9), 15] and accuracy ~?;
!

  include 'gaussfit_Yukawa.inc'

  p = p1
  w = w1

  urange = 15._dp
  drange = 1.d-09
  acc = 1.d-06 !relative error
!
END SUBROUTINE Yukawa_gequad


subroutine test_g2cart(n1,n2,n3,gprimd,g2cart)
   implicit none
   integer, intent(in) :: n1,n2,n3
   real(kind=8), dimension(3,3), intent(in) :: gprimd
   real(kind=8), dimension(n1*n2*n3), intent(out) :: g2cart
   !local variables
   integer :: count, i1,i2,i3,id1,id2,id3,ifft,ig1,ig2,ig3,ii1
   real(kind=8) :: b11,b12,b13,b21,b22,b23,b31,b32,b33


   id1=int(n1/2)+2
   id2=int(n2/2)+2
   id3=int(n3/2)+2
   ifft=0
   count=0
   do i3=1,n3
      ifft=(i3-1)*n1*n2
      ig3=i3-int(i3/id3)*n3-1
      do i2=1,n2
         ig2=i2-int(i2/id2)*n2-1
         ii1=1
         do i1=ii1,n1
            ig1=i1-int(i1/id1)*n1-1
            ifft=ifft+1

            b11=gprimd(1,1)*real(ig1,kind=8)
            b21=gprimd(2,1)*real(ig1,kind=8)
            b31=gprimd(3,1)*real(ig1,kind=8)
            b12=gprimd(1,2)*real(ig2,kind=8)
            b22=gprimd(2,2)*real(ig2,kind=8)
            b32=gprimd(3,2)*real(ig2,kind=8)
            b13=gprimd(1,3)*real(ig3,kind=8)
            b23=gprimd(2,3)*real(ig3,kind=8)
            b33=gprimd(3,3)*real(ig3,kind=8)

            g2cart(ifft)=( &
                 &     (b11+b12+b13)**2&
                 &     +(b21+b22+b23)**2&
                 &     +(b31+b32+b33)**2&
                 &     )
!!$     do ifunc=1,nfunc
!!$!     compute the laplacien in fourrier space
!!$!     that is * (i x 2pi x G)**2
!!$      laplacerdfuncg(1,ifft,ifunc) = -rdfuncg(1,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
!!$      laplacerdfuncg(2,ifft,ifunc) = -rdfuncg(2,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
!!$     end do
         end do
      end do
   end do
end subroutine test_g2cart
