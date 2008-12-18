!!****f* BigDFT/system_size
!! FUNCTION
!!   Calculates the overall size of the simulation cell 
!!   and shifts the atoms such that their position is the most symmetric possible.
!!   Assign these values to the global localisation region descriptor.
!! SOURCE
!!
!!***
subroutine system_size(iproc,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,Glr)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: iproc
  real(gp), intent(in) :: crmult,frmult
  real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  real(gp), intent(inout) :: hx,hy,hz
  type(locreg_descriptors), intent(out) :: Glr
  !local variables
  integer, parameter :: lupfil=14
  real(gp), parameter ::eps_mach=1.e-12_gp,onem=1.0_gp-eps_mach
  integer :: iat,j,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
  real(gp) :: rad,cxmin,cxmax,cymin,cymax,czmin,czmax,alatrue1,alatrue2,alatrue3

  !check the geometry code with the grid spacings
  if (atoms%geocode == 'F' .and. (hx/=hy .or. hx/=hz .or. hy/=hz)) then
     write(*,'(1x,a)')'ERROR: The values of the grid spacings must be equal in the Free BC case'
     stop
  end if


  !calculate the extremes of the boxes taking into account the spheres around the atoms
  cxmax=-1.e10_gp 
  cxmin=1.e10_gp

  cymax=-1.e10_gp 
  cymin=1.e10_gp

  czmax=-1.e10_gp 
  czmin=1.e10_gp

  do iat=1,atoms%nat

     rad=radii_cf(atoms%iatype(iat),1)*crmult

     cxmax=max(cxmax,rxyz(1,iat)+rad) 
     cxmin=min(cxmin,rxyz(1,iat)-rad)

     cymax=max(cymax,rxyz(2,iat)+rad) 
     cymin=min(cymin,rxyz(2,iat)-rad)
     
     czmax=max(czmax,rxyz(3,iat)+rad) 
     czmin=min(czmin,rxyz(3,iat)-rad)
  enddo

  cxmax=cxmax+eps_mach 
  cymax=cymax+eps_mach  
  czmax=czmax+eps_mach  

  cxmin=cxmin-eps_mach
  cymin=cymin-eps_mach
  czmin=czmin-eps_mach


  !define the box sizes for free BC, and calculate dimensions for the fine grid with ISF
  if (atoms%geocode == 'F') then
     atoms%alat1=(cxmax-cxmin)
     atoms%alat2=(cymax-cymin)
     atoms%alat3=(czmax-czmin)

     ! grid sizes n1,n2,n3
     n1=int(atoms%alat1/hx)
     n2=int(atoms%alat2/hy)
     n3=int(atoms%alat3/hz)
     alatrue1=real(n1,gp)*hx
     alatrue2=real(n2,gp)*hy
     alatrue3=real(n3,gp)*hz

     n1i=2*n1+31
     n2i=2*n2+31
     n3i=2*n3+31

  else if (atoms%geocode == 'P') then 
     !define the grid spacings, controlling the FFT compatibility
     call correct_grid(atoms%alat1,hx,n1)
     call correct_grid(atoms%alat2,hy,n2)
     call correct_grid(atoms%alat3,hz,n3)
     alatrue1=(cxmax-cxmin)
     alatrue2=(cymax-cymin)
     alatrue3=(czmax-czmin)

     n1i=2*n1+2
     n2i=2*n2+2
     n3i=2*n3+2

  else if (atoms%geocode == 'S') then
     call correct_grid(atoms%alat1,hx,n1)
     atoms%alat2=(cymax-cymin)
     call correct_grid(atoms%alat3,hz,n3)

     alatrue1=(cxmax-cxmin)
     n2=int(atoms%alat2/hy)
     alatrue2=real(n2,gp)*hy
     alatrue3=(czmax-czmin)

     n1i=2*n1+2
     n2i=2*n2+31
     n3i=2*n3+2

  end if

  !balanced shift taking into account the missing space
  cxmin=cxmin+0.5_gp*(atoms%alat1-alatrue1)
  cymin=cymin+0.5_gp*(atoms%alat2-alatrue2)
  czmin=czmin+0.5_gp*(atoms%alat3-alatrue3)

  !correct the box sizes for the isolated case
  if (atoms%geocode == 'F') then
     atoms%alat1=alatrue1
     atoms%alat2=alatrue2
     atoms%alat3=alatrue3
  else if (atoms%geocode == 'S') then
     cxmin=0.0_gp
     atoms%alat2=alatrue2
     czmin=0.0_gp
  else if (atoms%geocode == 'P') then
     !for the moment we do not put the shift, at the end it will be tested
     !here we should put the center of mass
     cxmin=0.0_gp
     cymin=0.0_gp
     czmin=0.0_gp
  end if

  do iat=1,atoms%nat
     rxyz(1,iat)=rxyz(1,iat)-cxmin
     rxyz(2,iat)=rxyz(2,iat)-cymin
     rxyz(3,iat)=rxyz(3,iat)-czmin
  enddo

  ! fine grid size (needed for creation of input wavefunction, preconditioning)
  nfl1=n1 
  nfl2=n2 
  nfl3=n3

  nfu1=0 
  nfu2=0 
  nfu3=0

  do iat=1,atoms%nat
     rad=radii_cf(atoms%iatype(iat),2)*frmult
     nfl1=min(nfl1,ceiling((rxyz(1,iat)-rad)/hx - eps_mach))
     nfu1=max(nfu1,floor((rxyz(1,iat)+rad)/hx + eps_mach))

     nfl2=min(nfl2,ceiling((rxyz(2,iat)-rad)/hy - eps_mach))
     nfu2=max(nfu2,floor((rxyz(2,iat)+rad)/hy + eps_mach))

     nfl3=min(nfl3,ceiling((rxyz(3,iat)-rad)/hz - eps_mach)) 
     nfu3=max(nfu3,floor((rxyz(3,iat)+rad)/hz + eps_mach))
  enddo

  !correct the values of the delimiter if they go outside the box
  if (nfl1 < 0 .or. nfu1 > n1) then
     nfl1=0
     nfu1=n1
  end if
  if (nfl2 < 0 .or. nfu2 > n2) then
     nfl2=0
     nfu2=n2
  end if
  if (nfl3 < 0 .or. nfu3 > n3) then
     nfl3=0
     nfu3=n3
  end if

  if (iproc == 0) then
     write(*,'(1x,a,19x,a)') 'Shifted atomic positions, Atomic Units:','grid spacing units:'
     do iat=1,atoms%nat
        write(*,'(1x,i5,1x,a6,3(1x,1pe12.5),3x,3(1x,0pf9.3))') &
             iat,trim(atoms%atomnames(atoms%iatype(iat))),&
             (rxyz(j,iat),j=1,3),rxyz(1,iat)/hx,rxyz(2,iat)/hy,rxyz(3,iat)/hz
     enddo
     write(*,'(1x,a,3(1x,1pe12.5),a,3(1x,0pf5.2))') &
          '   Shift of=',-cxmin,-cymin,-czmin,' Grid Spacings=',hx,hy,hz
     write(*,'(1x,a,3(1x,1pe12.5),3x,3(1x,i9))')&
          '  Box Sizes=',atoms%alat1,atoms%alat2,atoms%alat3,n1,n2,n3
     write(*,'(1x,a,3x,3(3x,i4,a1,i0))')&
          '      Extremes for the high resolution grid points:',&
          nfl1,'<',nfu1,nfl2,'<',nfu2,nfl3,'<',nfu3
  endif

  !assign the values
  Glr%d%n1  =n1  
  Glr%d%n2  =n2  
  Glr%d%n3  =n3  
  Glr%d%n1i =n1i 
  Glr%d%n2i =n2i 
  Glr%d%n3i =n3i 
  Glr%d%nfl1=nfl1
  Glr%d%nfl2=nfl2
  Glr%d%nfl3=nfl3
  Glr%d%nfu1=nfu1
  Glr%d%nfu2=nfu2
  Glr%d%nfu3=nfu3

  Glr%ns1=0
  Glr%ns2=0
  Glr%ns3=0

  !evaluate if the conditiond for the hybrid evaluation if periodic BC hold
  Glr%hybrid_on=                   (nfu1-nfl1+lupfil < n1+1)
  Glr%hybrid_on=(Glr%hybrid_on.and.(nfu2-nfl2+lupfil < n2+1))
  Glr%hybrid_on=(Glr%hybrid_on.and.(nfu3-nfl3+lupfil < n3+1))

end subroutine system_size

!!****f* BigDFT/correct_grid
!! FUNCTION
!!   Here the dimensions should be corrected in order to 
!!   allow the fft for the preconditioner and for Poisson Solver
!! SOURCE
!!
!!***
subroutine correct_grid(a,h,n)
  use module_base
  use Poisson_Solver
  implicit none
  real(gp), intent(in) :: a
  integer, intent(inout) :: n
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
     if (m2==2*m) exit

     nt=m+1
  end do
  n=m-1

!!$  !here the dimensions should be corrected in order to 
!!$  !allow the fft for the preconditioner
!!$  m=2*n+2
!!$  do 
!!$     call fourier_dim(m,m)
!!$     if ((m/2)*2==m) then
!!$        n=(m-2)/2
!!$        exit
!!$     else
!!$        m=m+1
!!$     end if
!!$  end do

  h=a/real(n+1,gp)
  
end subroutine correct_grid

! Calculates the length of the keys describing a wavefunction data structure
subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
  implicit none
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3
  logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid 
  integer, intent(out) :: mseg,mvctr
  !local variables
  logical :: plogrid
  integer :: i1,i2,i3,nsrt,nend

  mvctr=0
  nsrt=0
  nend=0
  do i3=nl3,nu3 
     do i2=nl2,nu2
        plogrid=.false.
        do i1=nl1,nu1
           if (logrid(i1,i2,i3)) then
              mvctr=mvctr+1
              if (plogrid .eqv. .false.) then
                 nsrt=nsrt+1
              endif
           else
              if (plogrid .eqv. .true.) then
                 nend=nend+1
              endif
           endif
           plogrid=logrid(i1,i2,i3)
        enddo
        if (plogrid .eqv. .true.) then
           nend=nend+1
        endif
     enddo
  enddo
  if (nend.ne.nsrt) then 
     write(*,*)' ERROR: nend <> nsrt',nend,nsrt
     stop 
  endif
  mseg=nend
  
end subroutine num_segkeys

! Calculates the keys describing a wavefunction data structure
subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
  !implicit real(kind=8) (a-h,o-z)
  implicit none
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,mseg
  logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid  
  integer, dimension(mseg), intent(out) :: keyv
  integer, dimension(2,mseg), intent(out) :: keyg
  !local variables
  logical :: plogrid
  integer :: mvctr,nsrt,nend,i1,i2,i3,ngridp

  mvctr=0
  nsrt=0
  nend=0
  do i3=nl3,nu3 
     do i2=nl2,nu2
     plogrid=.false.
     do i1=nl1,nu1
        ngridp=i3*((n1+1)*(n2+1)) + i2*(n1+1) + i1+1
        if (logrid(i1,i2,i3)) then
           mvctr=mvctr+1
           if (plogrid .eqv. .false.) then
              nsrt=nsrt+1
              keyg(1,nsrt)=ngridp
              keyv(nsrt)=mvctr
           endif
        else
           if (plogrid .eqv. .true.) then
              nend=nend+1
              keyg(2,nend)=ngridp-1
           endif
        endif
        plogrid=logrid(i1,i2,i3)
     enddo
     if (plogrid .eqv. .true.) then
        nend=nend+1
        keyg(2,nend)=ngridp
     endif
  enddo; enddo
  if (nend /= nsrt) then 
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif
  !mseg=nend
end subroutine segkeys

subroutine fill_logrid(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
     ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)
  ! set up an array logrid(i1,i2,i3) that specifies whether the grid point
  ! i1,i2,i3 is the center of a scaling function/wavelet
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,ntypes
  real(gp), intent(in) :: rmult,hx,hy,hz
  integer, dimension(nat), intent(in) :: iatype
  real(gp), dimension(ntypes), intent(in) :: radii
  real(gp), dimension(3,nat), intent(in) :: rxyz
  logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12,onem=1.d0-eps_mach
  integer :: i1,i2,i3,iat,ml1,ml2,ml3,mu1,mu2,mu3,j1,j2,j3
  real(gp) :: dx,dy2,dz2,rad

  !some checks
  if (geocode /='F') then
     !the nbuf value makes sense only in the case of free BC
     if (nbuf /=0) then
        write(*,'(1x,a)')'ERROR: a nonzero value of nbuf is allowed only for Free BC (tails)'
        stop
     end if
     !the grid spacings must be the same
     if (hx/= hy .or. hy /=hz .or. hx/=hz) then
!        write(*,'(1x,a)')'ERROR: For Free BC the grid spacings must be the same'
     end if
  end if

  if (geocode == 'F') then
     do i3=nl3,nu3 
        do i2=nl2,nu2 
           do i1=nl1,nu1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
  else !
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
  end if

  do iat=1,nat
     rad=radii(iatype(iat))*rmult+real(nbuf,gp)*hx
     !        write(*,*) 'iat,nat,rad',iat,nat,rad
     ml1=ceiling((rxyz(1,iat)-rad)/hx - eps_mach)  
     ml2=ceiling((rxyz(2,iat)-rad)/hy - eps_mach)   
     ml3=ceiling((rxyz(3,iat)-rad)/hz - eps_mach)   
     mu1=floor((rxyz(1,iat)+rad)/hx + eps_mach)
     mu2=floor((rxyz(2,iat)+rad)/hy + eps_mach)
     mu3=floor((rxyz(3,iat)+rad)/hz + eps_mach)
     !for Free BC, there must be no incoherences with the previously calculated delimiters
     if (geocode == 'F') then
        if (ml1 < nl1) stop 'ml1 < nl1'
        if (ml2 < nl2) stop 'ml2 < nl2'
        if (ml3 < nl3) stop 'ml3 < nl3'

        if (mu1 > nu1) stop 'mu1 > nu1'
        if (mu2 > nu2) stop 'mu2 > nu2'
        if (mu3 > nu3) stop 'mu3 > nu3'
     end if
     !what follows works always provided the check before
     do i3=ml3,mu3
        dz2=(real(i3,gp)*hz-rxyz(3,iat))**2
        j3=modulo(i3,n3+1)
        do i2=ml2,mu2
           dy2=(real(i2,gp)*hy-rxyz(2,iat))**2
           j2=modulo(i2,n2+1)
           do i1=ml1,mu1
              j1=modulo(i1,n1+1)
              dx=real(i1,gp)*hx-rxyz(1,iat)
              if (dx**2+(dy2+dz2) <= rad**2) then 
                 logrid(j1,j2,j3)=.true.
              endif
           enddo
        enddo
     enddo
  enddo

END SUBROUTINE fill_logrid

subroutine make_bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
  implicit none
  integer, intent(in) :: n1,n2,n3
  logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid
  integer, dimension(2,0:n2,0:n3), intent(out) :: ibyz
  integer, dimension(2,0:n1,0:n3), intent(out) :: ibxz
  integer, dimension(2,0:n1,0:n2), intent(out) :: ibxy
  !local variables
  integer :: i1,i2,i3

  do i3=0,n3 
     do i2=0,n2 
        ibyz(1,i2,i3)= 1000
        ibyz(2,i2,i3)=-1000

        loop_i1s: do i1=0,n1
           if (logrid(i1,i2,i3)) then 
              ibyz(1,i2,i3)=i1
              exit loop_i1s
           endif
        enddo loop_i1s

        loop_i1e: do i1=n1,0,-1
           if (logrid(i1,i2,i3)) then 
              ibyz(2,i2,i3)=i1
              exit loop_i1e
           endif
        enddo loop_i1e
     end do
  end do


  do i3=0,n3 
     do i1=0,n1
        ibxz(1,i1,i3)= 1000
        ibxz(2,i1,i3)=-1000

        loop_i2s: do i2=0,n2 
           if (logrid(i1,i2,i3)) then 
              ibxz(1,i1,i3)=i2
              exit loop_i2s
           endif
        enddo loop_i2s

        loop_i2e: do i2=n2,0,-1
           if (logrid(i1,i2,i3)) then 
              ibxz(2,i1,i3)=i2
              exit loop_i2e
           endif
        enddo loop_i2e

     end do
  end do


  do i2=0,n2 
     do i1=0,n1 
        ibxy(1,i1,i2)= 1000
        ibxy(2,i1,i2)=-1000

        loop_i3s: do i3=0,n3
           if (logrid(i1,i2,i3)) then 
              ibxy(1,i1,i2)=i3
              exit loop_i3s
           endif
        enddo loop_i3s

        loop_i3e: do i3=n3,0,-1
           if (logrid(i1,i2,i3)) then 
              ibxy(2,i1,i2)=i3
              exit loop_i3e
           endif
        enddo loop_i3e
     end do
  end do

end subroutine make_bounds

