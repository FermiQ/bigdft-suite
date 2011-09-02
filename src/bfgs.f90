!> @file
!!  Routines to do BFGS geometry optimisation
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
subroutine bfgsdriver(nproc,iproc,rxyz,fxyz,epot,at,rst,in,ncount_bigdft)
    use module_base
    use module_types
    use minpar
    implicit none
    integer, intent(in) :: nproc,iproc
    integer, intent(inout) :: ncount_bigdft
    type(atoms_data), intent(inout) :: at
    type(input_variables), intent(inout) :: in
    type(restart_objects), intent(inout) :: rst
    real(gp), intent(inout) :: epot
    real(gp), dimension(3*at%nat), intent(inout) :: rxyz
    real(gp), dimension(3*at%nat), intent(inout) :: fxyz
    real(gp) :: fluct=0.0_gp,fnrm,fmax,fnoise
    integer :: infocode,i,ixyz,iat,istat,icall,icheck
    character(len=4) :: fn4
    character(len=40) :: comment
    logical :: move_this_coordinate
    integer:: nr
    integer:: nwork
    real(gp),allocatable:: x(:),f(:),work(:)
    !character(len=4) :: fn4
    !character(len=40) :: comment
    !real(gp), dimension(3*at%nat) :: rxyz0,rxyzwrite
    !character(len=*), parameter :: subname='bfgs'

    in%inputPsiId=1
    in%output_grid=0
    in%output_wf=.false.
    icheck=0
    !if(iproc==0) write(*,*) 'EPOT=',epot
    !return

    nr=0
    do i=1,3*at%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(at%ifrztyp(iat),ixyz)) nr=nr+1
    enddo
    parmin%iflag=0
    nwork=nr*nr+3*nr+3*nr*nr+3*nr
    allocate(work(nwork),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating work.'
    allocate(x(nr),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating x.'
    allocate(f(nr),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating f.'
    icall=0
    do 
        !call nebforce(n,np,x,f,fnrmtot,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
        !do ip=1,np-1
        !    call atomic_copymoving_forward(atoms,n,f(1,ip),nr,fa(1,ip))
        !enddo
        !if(icall/=0) then
            call call_bigdft(nproc,iproc,at,rxyz,in,epot,fxyz,fnoise,rst,infocode)
            ncount_bigdft=ncount_bigdft+1
        !endif
        call atomic_copymoving_forward(at,3*at%nat,fxyz,nr,f)
        call atomic_copymoving_forward(at,3*at%nat,rxyz,nr,x)
        call fnrmandforcemax(fxyz,fnrm,fmax,at%nat)
        if(fmax<3.d-1) call updatefluctsum(at%nat,fnoise,fluct)
        call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,icheck)
        if(iproc==0) write(*,*) 'ICHECK ',icheck
        if(icheck>5) parmin%converged=.true.
        !call calmaxforcecomponentanchors(atoms,np,f(1,1),fnrm,fspmax)
        !call checkconvergence(parmin,fspmax)
        !if(ncount_bigdft>in%ncount_cluster_x-1)
        !if(iproc==0) write(*,*) 'nr=',nr,f(1)
        if (iproc == 0) then
           write(fn4,'(i4.4)') ncount_bigdft
           write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
           call  write_atomic_file('posout_'//fn4,epot,rxyz,at,trim(comment))
        endif
        call bfgs_reza(iproc,nr,x,epot,f,nwork,work,in%betax,sqrt(fnrm),fmax, &
            ncount_bigdft,fluct*in%frac_fluct,fluct)
        !x(1:nr)=x(1:nr)+1.d-2*f(1:nr)
        call atomic_copymoving_backward(at,nr,x,3*at%nat,rxyz)
        if(parmin%converged) then
           if(iproc==0) write(16,'(a,i0,a)') "   BFGS converged in ",icall," iterations"
           if(iproc==0) then
              write(fn4,'(i4.4)') ncount_bigdft
              write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
              call  write_atomic_file('posout_'//fn4,epot,rxyz,at,trim(comment))
           endif
        endif
        !if(ncount_bigdft>in%ncount_cluster_x-1)
        !do ip=1,np-1
        !    call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
        !enddo
        if(parmin%converged) exit
        if(parmin%iflag<=0) exit
        icall=icall+1
        if(icall>in%ncount_cluster_x) exit
    enddo
    deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    deallocate(x,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating x.'
    deallocate(f,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating f.'
END SUBROUTINE bfgsdriver
!*****************************************************************************************
!subroutine inithess(iproc,nat,rat,at)
!    use module_types
!    implicit none
!    integer::iproc,nat,iat
!    real(8)::rat(3,nat)
!    type(atoms_data), intent(inout) :: at
!    integer, allocatable::ita(:)
!    allocate(ita(nat))
!    do iat=1,nat
!        if(trim(atoms%atomnames(atoms%iatype(iat)))=='H') then
!            ita(iat)=1
!        elseif(trim(atoms%atomnames(atoms%iatype(iat)))=='C') then
!            ita(iat)=2
!        elseif(trim(atoms%atomnames(atoms%iatype(iat)))=='N') then
!            ita(iat)=3
!        elseif(trim(atoms%atomnames(atoms%iatype(iat)))=='O') then
!            ita(iat)=4
!        else
!            write(*,*) 'This subroutine works only for systems which '
!            write(*,*) 'contain only ogganic elements, namely H,C,N,O.'
!            stop
!        endif
!    enddo
!    do iat=1,nat
!        do jat=iat+1,nat
!            dx=rat(1,jat)-rat(1,iat)
!            dy=rat(2,jat)-rat(2,iat)
!            dz=rat(3,jat)-rat(3,iat)
!            r=sqrt(dx**2+dy**2+dz**2)
!            if(r<1.35d0*r0(ita(iat),ita(iat))) then
!                nln(iat)=nln(iat)+1
!                nln(jat)=nln(jat)+1
!            endif
!        enddo
!    enddo
!    deallocate(ita)
!end subroutine inithess
!*****************************************************************************************
subroutine bfgs_reza(iproc,nr,x,epot,f,nwork,work,alphax,fnrm,fmax,ncount_bigdft,flt1,flt2)
    use minpar
    implicit none
    integer::iproc,nr,nwork,mf,my,ms,nrsqtwo,iw1,iw2,iw3,iw4,info,i,j,l,mx
    integer::ncount_bigdft
    real(8)::x(nr),f(nr),epot,work(nwork),alphax,flt1,flt2
    integer, allocatable::ipiv(:)
    !real(8), allocatable::eval(:),umat(:)
    !type(parameterminimization)::parmin
    real(8)::DDOT,DNRM2,tt1,tt2,de,fnrm,fmax,beta
    real(8)::tt3,tt4,tt5,tt6
    real(8), save::epotold,alpha,alphamax,zeta
    logical, save::reset
    integer, save::isatur
    if(nwork/=nr*nr+3*nr+3*nr*nr+3*nr) then
        stop 'ERROR: size of work array is insufficient.'
    endif
    nrsqtwo=nr*nr*2
    mf=nr*nr+1       !for force of previous iteration in wiki notation
    my=mf+nr         !for y_k in wiki notation
    ms=my+nr         !for s_k in wiki notation
    iw1=ms+nr        !work array to keep the hessian untouched
    iw2=iw1+nr*nr    !for work array of DSYTRF
    iw3=iw2+nrsqtwo  !for p_k in wiki notation
    mx =iw3+nr       !for position of previous iteration
    iw4=mx+nr        !for eigenvalues of inverse og hessian
    if(parmin%iflag==0) then
        parmin%iflag=1
        parmin%iter=0
        epotold=epot
        alpha=8.d-1
        reset=.false.
        alphamax=0.9d0
        zeta=1.d0
        isatur=0
        if(iproc==0) then
        open(unit=1390,file='bfgs_eigenvalues.dat',status='replace')
        close(1390)
        endif
    else
        parmin%iter=parmin%iter+1
    endif
    if(fnrm<min(6.d-2,max(1.d-2,2.d-3*sqrt(real(nr,8))))) then
        if(isatur<99) isatur=isatur+1
    else
        isatur=0
    endif
    de=epot-epotold
    !fnrm=calnorm(nr,f);fmax=calmaxforcecomponent(nr,f)
    if(iproc==0) then
    !write(*,'(a10,i5,es23.15,es11.3,2es12.5,2es12.4,i3)') &
    !    'GEOPT_BFGS',parmin%iter,epot,de,fnrm,fmax,zeta,alpha,isatur
    !       '(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)'
    write(*,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a7,i3)') &
        ncount_bigdft,parmin%iter,'GEOPT_BFGS',epot,de,fmax,fnrm,flt1,flt2,'isatur=',isatur
    write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a7,i3)') &
        ncount_bigdft,parmin%iter,'GEOPT_BFGS',epot,de,fmax,fnrm,flt1,flt2,'isatur=',isatur
    endif
    close(16)
    open(unit=16,file='geopt.mon',status='unknown',position='APPEND')
    !if(parmin%iter==602) then
    !    do i=1,nr/3
    !        write(31,*) x(i*3-2),x(i*3-1),x(i*3-0)
    !    enddo
    !    stop
    !endif
    !if(fmax<parmin%fmaxtol) then
    if(parmin%converged) then
        !parmin%converged=.true.
        parmin%iflag=0
        if(iproc==0) then
        write(*,'(a,i4,es23.15,2es12.5)') &
            'BFGS FINISHED: itfire,epot,fnrm,fmax ',parmin%iter,epot,fnrm,fmax
        endif
        return
    endif

    !if(de>0.d0 .and. zeta>1.d-1) then
    if(de>5.d-2) then
        epot=epotold
        x(1:nr)=work(mx:mx-1+nr)
        f(1:nr)=work(mf:mf-1+nr)
        reset=.true.
        !alpha=max(alpha*0.5d0/1.1d0,1.d-2)
        zeta=max(zeta*2.d-1,1.d-3)
        isatur=0
    else
        !zeta=1.d0
        !if(zeta>1.d-1) zeta=min(zeta*1.1d0,1.d0)
        zeta=min(zeta*1.1d0,1.d0)
        !isatur=isatur+1
    endif
    if(parmin%iter==0 .or. reset) then
        reset=.false.
        !if(isatur>=10) then
        !    reset=.false.
        !    !alpha=5.d-1
        !endif
        work(1:nr*nr)=0.d0
        do i=1,nr
            work(i+(i-1)*nr)=zeta*alphax
        enddo
        work(iw3:iw3-1+nr)=zeta*alphax*f(1:nr)
    else
        work(ms:ms-1+nr)=x(1:nr)-work(mx:mx-1+nr)
        work(my:my-1+nr)=work(mf:mf-1+nr)-f(1:nr)
        tt1=DDOT(nr,work(my),1,work(ms),1)
        do i=1,nr
            tt2=0.d0
            do j=1,nr
                tt2=tt2+work(i+(j-1)*nr)*work(my-1+j)
            enddo
            work(iw2-1+i)=tt2
        enddo
        tt2=DDOT(nr,work(my),1,work(iw2),1)
        !write(21,*) parmin%iter,tt1,tt2
        !tt1=max(tt1,1.d-2)
        do i=1,nr
            do j=i,nr
                l=i+(j-1)*nr
                work(l)=work(l)+(tt1+tt2)*work(ms-1+i)*work(ms-1+j)/tt1**2- &
                    (work(iw2-1+i)*work(ms-1+j)+work(iw2-1+j)*work(ms-1+i))/tt1
                work(j+(i-1)*nr)=work(l)
            enddo
        enddo
        !do i=1,nr
        !    tt2=0.d0
        !    do j=1,nr
        !        tt2=tt2+work(j+(i-1)*nr)*f(j)
        !    enddo
        !    work(iw3-1+i)=tt2
        !enddo
        !write(31,*) zeta
        work(iw1:iw1-1+nr*nr)=work(1:nr*nr)
        call DSYEV('V','L',nr,work(iw1),nr,work(iw4),work(iw2),nrsqtwo,info)
        if(info/=0) stop 'ERROR: DSYEV in bfgs_reza failed.'
        tt1=work(iw4+0)    ; tt2=work(iw4+1)    ; tt3=work(iw4+2)
        tt4=work(iw4+nr-3) ; tt5=work(iw4+nr-2) ; tt6=work(iw4+nr-1)
        if(iproc==0) then
        open(unit=1390,file='bfgs_eigenvalues.dat',status='old',position='append')
        write(1390,'(i5,6es15.5)') parmin%iter,tt1,tt2,tt3,tt4,tt5,tt6
        close(1390)
        endif
        work(iw3:iw3-1+nr)=0.d0
        if(isatur<3) then
            beta=1.d0/alphax
        elseif(isatur<6) then
            beta=1.d-1/alphax
        elseif(isatur<10) then
            beta=1.d-2/alphax
        else
            beta=1.d-3/alphax
        endif
        !do j=1,nr
        !    if(work(iw4-1+j)>0.d0) then
        !        tt3=work(iw4-1+j)
        !        exit
        !    enddo
        !enddo
        tt3=alphax*0.5d0
        do j=1,nr
            tt1=DDOT(nr,work(iw1+nr*(j-1)),1,f,1)
            if(work(iw4-1+j)<tt3) then
                tt4=tt3
            else
                tt4=work(iw4-1+j)
            endif
            tt2=1.d0/sqrt(1.d0/tt4**2+beta**2)
            do i=1,nr
                work(iw3-1+i)=work(iw3-1+i)+tt1*work(iw1-1+i+nr*(j-1))*tt2
            enddo
        enddo
    endif
    epotold=epot
    work(mf:mf-1+nr)=f(1:nr)
    work(mx:mx-1+nr)=x(1:nr)
    alpha=min(alphamax,alpha*1.1d0)
    x(1:nr)=x(1:nr)+alpha*work(iw3:iw3-1+nr)
end subroutine bfgs_reza
!*****************************************************************************************

!>   Driver for the LBFGS routine found on the Nocedal Homepage
!!   The subroutines have only been modified slightly, so a VIMDIFF will show all modifications!
!!   This is helpfull when we are looking for the source of problems during BFGS runs
subroutine lbfgsdriver(nproc,iproc,rxyz,fxyz,etot,at,rst,in,ncount_bigdft,fail) 
  use module_base
  use module_types
!  use par_driver
  use minpar
  implicit none
!  type(driverparameters)::par
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), intent(inout) :: etot
  real(gp), dimension(3*at%nat), intent(inout) :: rxyz
  logical, intent(out) :: fail
  real(gp), dimension(3*at%nat), intent(out) :: fxyz

  real(gp), dimension(3*at%nat):: txyz, sxyz
  real(gp) :: fluct,fnrm, fnoise
  real(gp) :: fmax
!  logical :: check
  integer :: check
  integer :: infocode,i,ixyz,iat,nitsd
  real(gp) :: fnormmax_sw,etotprev
  character(len=4) :: fn4
  character(len=40) :: comment
  logical :: move_this_coordinate

  integer:: n,nr,ndim
  integer:: NWORK
  real(gp),allocatable:: X(:),G(:),DIAG(:),W(:)
  real(gp):: F,EPS!,XTOL,GTOL,,STPMIN,STPMAX
  real(gp), dimension(3*at%nat) :: rxyz0,rxyzwrite
  INTEGER:: IPRINT(2),IFLAG,ICALL,M
  character(len=*), parameter :: subname='bfgs'
  integer :: i_stat,i_all

  check=0

!  call init_driver(par)     !Initialize the parameters
  parmin%finstep=0
  parmin%alpha=1.d0
  fail=.false.
  fnrm=1.d10
  nitsd=10!500                 !Maximum number of SD steps before entering BFGS
  fnormmax_sw=in%forcemax!1.e-2_gp      !SD till the max force comp is less than this value
  
  !Dummy variables
  txyz=0._gp
  sxyz=0._gp
  
  if (iproc==0)    write(*,*) 'Maximum number of SD steps used in the beginning: ',nitsd

  call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_bigdft,fnrm,fnoise,in,&
       fnormmax_sw,nitsd,fluct)
  etotprev=etot
  rxyz0=rxyz     !Save initial positions, since the unconstrained degrees of freedom will be updated upon them
  rxyzwrite=rxyz
  call fnrmandforcemax(fxyz,fnrm,fmax,at%nat)
  !call fnrmandforcemax(fxyz,fnrm,fmax,at)
  !check if the convergence is reached after SD
  call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)

  if (check.gt.5) then
     if (iproc.eq.0) write(*,*) 'Converged before entering BFGS'
     return
  endif


  !Make a list of all degrees of freedom that should be passed to bfgs
  n=3*at%nat
  nr=0
  do i=1,3*at%nat
     iat=(i-1)/3+1
     ixyz=mod(i-1,3)+1
     if(move_this_coordinate(at%ifrztyp(iat),ixyz)) nr=nr+1
  enddo
  if(iproc==0) write(*,*) 'DOF: n,nr ',n,nr
     NDIM=nr
     NWORK=NDIM*(2*parmin%MSAVE +1)+2*parmin%MSAVE
      
     allocate(X(NDIM),stat=i_stat)
     call memocc(i_stat,X,'X',subname)
     allocate(G(NDIM),stat=i_stat)
     call memocc(i_stat,G,'G',subname)
     allocate(DIAG(NDIM),stat=i_stat)
     call memocc(i_stat,DIAG,'DIAG',subname)
     allocate(W(NWORK),stat=i_stat)
     call memocc(i_stat,W,'W',subname)

     call atomic_copymoving_forward(at,n,rxyz,nr,X)

     N=nr
     M=parmin%MSAVE
     IPRINT(1)= 1
     IPRINT(2)= 0
     F=etot
!     We do not wish to provide the diagonal matrices Hk0, and 
!     therefore set DIAGCO to FALSE.

     EPS=0.0_gp
     ICALL=0
     IFLAG=0

 20   CONTINUE
        if (parmin%IWRITE) then
           if (iproc == 0) then
              write(fn4,'(i4.4)') ncount_bigdft
              write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
              call  write_atomic_file('posout_'//fn4,etot,rxyz,at,trim(comment))
           endif
           parmin%IWRITE=.false.
        endif
        rxyzwrite=rxyz

        if (fmax < 3.d-1) call updatefluctsum(at%nat,fnoise,fluct)
   if (iproc==0.and.ICALL.ne.0.and.parmin%verbosity > 0) & 
              &write(16,'(I5,1x,I5,2x,a11,1x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)')&
              &ncount_bigdft,ICALL,"GEOPT_LBFGS",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct&
              &,"BFGS-it=",parmin%finstep,"alpha=",parmin%alpha
   if (iproc==0.and.ICALL.ne.0.and.parmin%verbosity > 0) & 
              & write(* ,'(I5,1x,I5,2x,a11,1x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)')&
              &ncount_bigdft,ICALL,"GEOPT_LBFGS",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct&
              &,"BFGS-it=",parmin%finstep,"alpha=",parmin%alpha
              etotprev=etot
              if (iproc==0.and.ICALL.ne.0.and.parmin%verbosity > 0) write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
                           'FORCES norm(Ha/Bohr): maxval=',fmax,'fnrm2=',fnrm,'fluct=', fluct
              call convcheck(fnrm,fmax,fluct*in%frac_fluct, in%forcemax,check)
              if (ncount_bigdft >= in%ncount_cluster_x) goto 50
              close(16)
              open(unit=16,file='geopt.mon',status='unknown',position='APPEND')

      if(check.gt.5) then
         if(iproc==0)  write(16,'(a,i0,a)') "   BFGS converged in ",ICALL," iterations"
         if (iproc == 0) then
            write(fn4,'(i4.4)') ncount_bigdft
            write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
            call  write_atomic_file('posout_'//fn4,etot,rxyz,at,trim(comment))
         endif
         goto 100
      endif

      
      rxyz=rxyz0
      call atomic_copymoving_backward(at,nr,X,n,rxyz)
!      txyz=rxyz
!      alpha=0._gp
!      call atomic_axpy(at,txyz,alpha,sxyz,rxyz)
      in%inputPsiId=1
      in%output_grid=0
      in%output_wf=.false.
!      if(ICALL.ne.0) call call_bigdft(nproc,iproc,at,rxyz,in,F,fxyz,rst,infocode)
      if(ICALL.ne.0) call call_bigdft(nproc,iproc,at,rxyz,in,F,fxyz,fnoise,rst,infocode)
      if(ICALL.ne.0) ncount_bigdft=ncount_bigdft+1
      call atomic_copymoving_forward(at,n,fxyz,nr,G)
      etot=F
      G=-G
      call fnrmandforcemax(fxyz,fnrm,fmax,at%nat)
!      call fnrmandforcemax(fxyz,fnrm,fmax,at)

      CALL LBFGS(IPROC,IN,PARMIN,N,M,X,F,G,DIAG,IPRINT,EPS,W,IFLAG)
      IF(IFLAG.LE.0) GO TO 50
      ICALL=ICALL + 1
!     We allow at most the given number of evaluations of F and G
      if(ncount_bigdft>in%ncount_cluster_x-1)  then
        goto 100
      endif
      close(16)
      open(unit=16,file='geopt.mon',status='unknown',position='append')
      GO TO 20
  50  CONTINUE
        if (iproc==0) write(*,*) "# Error in BFGS, switching to SD and CG"
        if (iproc==0) write(16,*) "Error in BFGS, switching to SD and CG"
        rxyz(:)=rxyzwrite(:)
        fail=.true.
 100  CONTINUE
        
      i_all=-product(shape(X))*kind(X)
      deallocate(X,stat=i_stat)
      call memocc(i_stat,i_all,'X',subname)
      i_all=-product(shape(G))*kind(G)
      deallocate(G,stat=i_stat)
      call memocc(i_stat,i_all,'G',subname)
      i_all=-product(shape(DIAG))*kind(DIAG)
      deallocate(DIAG,stat=i_stat)
      call memocc(i_stat,i_all,'DIAG',subname)
      i_all=-product(shape(W))*kind(W)
      deallocate(W,stat=i_stat)
      call memocc(i_stat,i_all,'W',subname)

END SUBROUTINE lbfgsdriver

subroutine atomic_copymoving_forward(atoms,n,x,nr,xa)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer::n,nr,i,iat,ixyz,ir
    real(8)::x(n),xa(nr)
    logical::move_this_coordinate
    ir=0
    do i=1,3*atoms%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%ifrztyp(iat),ixyz)) then
            ir=ir+1
            xa(ir)=x(i)
        endif
    enddo
    if(ir/=nr) stop 'ERROR: inconsistent number of relaxing DOF'
END SUBROUTINE atomic_copymoving_forward


subroutine atomic_copymoving_backward(atoms,nr,xa,n,x)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer::n,nr,i,iat,ixyz,ir
    real(8)::x(n),xa(nr)
    logical::move_this_coordinate
    ir=0
    do i=1,3*atoms%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%ifrztyp(iat),ixyz)) then
            ir=ir+1
            x(i)=xa(ir)
        endif
    enddo
    if(ir/=nr) stop 'ERROR: inconsistent number of relaxing DOF'
END SUBROUTINE atomic_copymoving_backward
