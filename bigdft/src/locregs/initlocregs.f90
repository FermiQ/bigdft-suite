!> @file
!! Basic operations relative to the localization region iniitalization.
!! All these routines are appelled from the module locregs_init, they are not used from outside
!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine determine_boxbounds_sphere(gperx, gpery, gperz, n1glob, n2glob, n3glob, nl1glob, nl2glob, nl3glob, &
     hx, hy, hz, locrad, locregCenter, &
     nsegglob, keygglob, keyvglob, ixmin, iymin, izmin, ixmax, iymax, izmax)
  use dynamic_memory
  implicit none
  logical,intent(in) :: gperx, gpery, gperz
  integer, intent(in) :: n1glob, n2glob, n3glob, nl1glob, nl2glob, nl3glob, nsegglob
  real(kind=8),intent(in) :: hx, hy, hz, locrad
  real(kind=8),dimension(3),intent(in) :: locregCenter
  integer,dimension(2,nsegglob),intent(in) :: keygglob
  integer,dimension(nsegglob),intent(in) :: keyvglob
  integer,intent(out) :: ixmin, iymin, izmin, ixmax, iymax, izmax
  !local variables
  integer :: i, i1, i2, i3, iseg, jj, j0, j1, ii, i0, ii1, ii2, ii3, n1p1, np
  integer :: ij1, ij2 ,ij3, jj1, jj2, jj3
  integer :: ijs1, ije1, ijs2, ije2, ijs3, ije3
  real(kind=8) :: cut, dx,dy, dz
  !debug
  integer :: iiimin, isegmin

  call f_routine(id='determine_boxbounds_sphere')

  ! For perdiodic boundary conditions, one has to check also in the neighboring
  ! cells (see in the loop below)
  if (gperx) then
     ijs1 = -1
     ije1 = 1
  else
     ijs1 = 0
     ije1 = 0
  end if
  if (gpery) then
     ijs2 = -1
     ije2 = 1
  else
     ijs2 = 0
     ije2 = 0
  end if
  if (gperz) then
     ijs3 = -1
     ije3 = 1
  else
     ijs3 = 0
     ije3 = 0
  end if

  iiimin=0
  isegmin=0

  ! Initialize the return values
  ixmax=0
  iymax=0
  izmax=0
  ixmin=nl1glob+n1glob
  iymin=nl2glob+n2glob
  izmin=nl3glob+n3glob

  cut=locrad**2
  n1p1=n1glob+1
  np=n1p1*(n2glob+1)
  !$omp parallel default(none) &
  !$omp shared(nsegglob,keygglob,n1glob,n2glob,n3glob,nl1glob,nl2glob,nl3glob,locregCenter) &
  !$omp shared(ixmin,iymin,izmin,ixmax,iymax,izmax,hx,hy,hz,cut,n1p1,np,ijs1,ije1,ijs2,ije2,ijs3,ije3) &
  !$omp private(iseg,jj,j0,j1,ii,i3,i2,i0,i1,ii2,ii3,ii1,i,dx,dy,dz,iiimin,isegmin) &
  !$omp private(ij1, ij2, ij3, jj1, jj2, jj3)
  !$omp do reduction(max:ixmax,iymax,izmax) reduction(min:ixmin,iymin,izmin)
  do iseg=1,nsegglob
     j0=keygglob(1,iseg)
     j1=keygglob(2,iseg)
     ii=j0-1
     i3=ii/np
     ii=ii-i3*np
     i2=ii/n1p1
     i0=ii-i2*n1p1
     i1=i0+j1-j0

     ii2=i2+nl2glob
     ii3=i3+nl3glob

     !dz=((ii3*hz)-locregCenter(3))**2
     !dy=((ii2*hy)-locregCenter(2))**2
     do ij3=ijs3,ije3!-1,1
        jj3=ii3+ij3*(n3glob+1)
        dz=((jj3*hz)-locregCenter(3))**2
        if(dz<=cut) then
           ! From the viewpoint of the z coordinate we are inside the cutoff, check now also the y and x dimensions
           do ij2=ijs2,ije2!-1,1
              jj2=ii2+ij2*(n2glob+1)
              dy=((jj2*hy)-locregCenter(2))**2
              if(dy+dz<=cut) then
                 ! From the viewpoint of the y and z coordinate we are inside the cutoff, check now also the x dimension
                 do i=i0,i1
                    ii1=i+nl1glob
                    do ij1=ijs1,ije1!-1,1
                       jj1=ii1+ij1*(n1glob+1)
                       dx=((jj1*hx)-locregCenter(1))**2
                       if(dx+dy+dz<=cut) then
                          ixmax=max(jj1,ixmax)
                          iymax=max(jj2,iymax)
                          izmax=max(jj3,izmax)
                          ixmin=min(jj1,ixmin)
                          iymin=min(jj2,iymin)
                          izmin=min(jj3,izmin)
                       end if
                    end do
                 end do
              end if
           end do
        end if
     end do
     !dx=((ii1*hx)-locregCenter(1))**2
     !!dx=((ii1*hx)-locregCenter(1))**2
     !!if(dx+dy+dz<=cut) then
     !!    ixmax=max(ii1,ixmax)
     !!    iymax=max(ii2,iymax)
     !!    izmax=max(ii3,izmax)
     !!    ixmin=min(ii1,ixmin)
     !!    !if(ii1<ixmin) iiimin=j0-1 ; isegmin=iseg
     !!    iymin=min(ii2,iymin)
     !!    izmin=min(ii3,izmin)
     !!end if
  end do
  !$omp enddo
  !$omp end parallel

  call f_release_routine()

END SUBROUTINE determine_boxbounds_sphere


subroutine num_segkeys_sphere(perx, pery, perz, n1, n2, n3, nl1glob, nl2glob, nl3glob, hx, hy, hz, &
     locrad, locregCenter, &
     nsegglob, keygglob, keyvglob, nseg, nvctr)
  use module_defs
  use dynamic_memory
  implicit none
  logical,intent(in) :: perx, pery, perz 
  integer, intent(in) :: n1, n2, n3, nl1glob, nl2glob, nl3glob, nsegglob
  real(kind=8),intent(in) :: hx, hy, hz, locrad
  real(kind=8),dimension(3),intent(in) :: locregCenter
  integer,dimension(2,nsegglob),intent(in) :: keygglob
  integer,dimension(nsegglob),intent(in) :: keyvglob
  integer,intent(out) :: nseg, nvctr
  !local variables
  logical :: segment, inside
  integer :: i, i1, i2, i3, nstart, nend, iseg, jj, j0, j1, ii, i0, ii1, ii2, ii3, n1p1, np
  integer :: ij1, ij2, ij3, jj1, jj2, jj3, ijs1, ijs2, ijs3, ije1, ije2, ije3
  real(kind=8) :: cut, dx,dy, dz

  call f_routine(id='num_segkeys_sphere')

  nvctr=0
  nstart=0
  nend=0

  cut=locrad**2
  n1p1=n1+1
  np=n1p1*(n2+1)

  ! For perdiodic boundary conditions, one has to check also in the neighboring
  ! cells (see in the loop below)
  if (perx) then
     ijs1 = -1
     ije1 = 1
  else
     ijs1 = 0
     ije1 = 0
  end if
  if (pery) then
     ijs2 = -1
     ije2 = 1
  else
     ijs2 = 0
     ije2 = 0
  end if
  if (perz) then
     ijs3 = -1
     ije3 = 1
  else
     ijs3 = 0
     ije3 = 0
  end if

  !$omp parallel default(none) &
  !$omp shared(nsegglob,keygglob,nl1glob,nl2glob,nl3glob,locregCenter) &
  !$omp shared(hx,hy,hz,cut,n1p1,np,nstart,nvctr,nend, n1, n2, n3, ijs1, ijs2, ijs3, ije1, ije2, ije3) &
  !$omp private(iseg,jj,j0,j1,ii,i3,i2,i0,i1,ii2,ii3,ii1,i,dx,dy,dz,segment) &
  !$omp private(inside, ij1, ij2, ij3, jj1, jj2, jj3)
  segment=.false.
  !$omp do schedule(dynamic,50) reduction(+:nstart,nvctr,nend)
  do iseg=1,nsegglob
     j0=keygglob(1,iseg)
     j1=keygglob(2,iseg)
     ii=j0-1
     i3=ii/np
     ii=ii-i3*np
     i2=ii/n1p1
     i0=ii-i2*n1p1
     i1=i0+j1-j0

     ii2=i2+nl2glob
     ii3=i3+nl3glob

     ! First just check the z dimension. If inside is false, proceed directly,
     ! otherwise check also the other dimensions.
     inside=.false.
     do ij3=ijs3,ije3!-1,1
        jj3=ii3+ij3*(n3+1)
        dz=((jj3*hz)-locregCenter(3))**2
        if(dz<=cut) then
           inside=.true.
        end if
     end do
     check_z_if: if (inside) then
        ! May be inside the sphere, so check also the other dimensions.
        ! Since each line in y (and thus also each plane in the z dimensions) starts
        ! a new segment, the following does not have to be done.
        inside = .false.

        !dz=((ii3*hz)-locregCenter(3))**2
        !dy=((ii2*hy)-locregCenter(2))**2
        do i=i0,i1
           ii1=i+nl1glob
           !dx=((ii1*hx)-locregCenter(1))**2
           inside=.false.
           do ij3=ijs3,ije3!-1,1
              jj3=ii3+ij3*(n3+1)
              dz=((jj3*hz)-locregCenter(3))**2
              do ij2=ijs2,ije2!-1,1
                 jj2=ii2+ij2*(n2+1)
                 dy=((jj2*hy)-locregCenter(2))**2
                 do ij1=ijs1,ije1!-1,1
                    jj1=ii1+ij1*(n1+1)
                    dx=((jj1*hx)-locregCenter(1))**2
                    !write(*,'(a,6i7,4es12.4)') 'ii1, ii2, ii3, jj1, jj2, jj3, dx, dy, dz, cut', ii1, ii2, ii3, jj1, jj2, jj3, dx, dy, dz, cut
                    if(dx+dy+dz<=cut) then
                       !write(*,'(a,6i7,4es12.4)') 'ii1, ii2, ii3, jj1, jj2, jj3, dx, dy, dz, cut', ii1, ii2, ii3, jj1, jj2, jj3, dx, dy, dz, cut
                       inside=.true.
                    end if
                 end do
              end do
           end do
           if(inside) then
              nvctr=nvctr+1
              if(.not.segment) then
                 nstart=nstart+1
                 segment=.true.
              end if
           else
              if(segment) then
                 nend=nend+1
                 segment=.false.
              end if
           end if
        end do
        if(segment) then
           ! Always start a new segment if we come to a new line in y direction.
           nend=nend+1
           segment=.false.
        end if
     end if check_z_if
  end do
  !$omp enddo
  !$omp end parallel

  nseg=nstart

  !check
  if (nend /= nstart) then
     write(*,*) 'nend , nstart',nend,nstart
     stop 'nend <> nstart'
  endif

  call f_release_routine()

END SUBROUTINE num_segkeys_sphere


subroutine segkeys_Sphere(perx, pery, perz, n1, n2, n3, nl1glob, nl2glob, nl3glob, &
     nl1, nu1, nl2, nu2, nl3, nu3, nseg, hx, hy, hz, &
     locrad, locregCenter, &
     nsegglob, keygglob, keyvglob, nvctr_loc, keyg_loc, keyg_glob, keyv_loc, keyv_glob, keygloc)
  use dynamic_memory
  use dictionaries
  use sparsematrix_init, only: distribute_on_threads
  implicit none
  logical,intent(in) :: perx, pery, perz
  integer,intent(in) :: n1, n2, n3, nl1glob, nl2glob, nl3glob, nl1, nu1, nl2, nu2, nl3, nu3, nseg, nsegglob, nvctr_loc
  real(kind=8) :: hx, hy, hz, locrad
  real(kind=8),dimension(3) :: locregCenter
  integer,dimension(2,nsegglob),intent(in) :: keygglob
  integer,dimension(nsegglob),intent(in) :: keyvglob
  integer,dimension(2,nseg),intent(out) :: keyg_loc, keyg_glob
  integer,dimension(nseg),intent(out) :: keyv_loc, keyv_glob
  integer,dimension(2,nseg),intent(inout) :: keygloc !tmp
  !local variables
  character(len=*),parameter :: subname = 'segkeys_Sphere'
  integer :: i, i1, i2, i3, nstart, nend, nvctr, igridpoint, igridglob, iseg, jj, j0, j1, ii, i0, n1l, n2l, n3l
  integer :: i1l, i2l, i3l, ii1, ii2, ii3, loc, n1p1, np, n1lp1, nlp, igridgloba
  !integer :: igridpointa
  integer :: ij1, ij2, ij3, jj1, jj2, jj3, ii1mod, ii2mod, ii3mod, ivctr, jvctr, kvctr, ijs1, ijs2, ijs3, ije1, ije2, ije3
  integer :: nsegglob_start, nsegglob_end
  real(kind=8) :: cut, dx, dy, dz
  logical :: segment, inside
  integer,dimension(:,:),pointer :: ise
  integer :: ithread, jthread, nthread, ivctr_tot, jvctr_tot, nstart_tot, nend_tot, kthread, j, offset
  integer,dimension(:),allocatable :: nstartarr, keyv_last
  integer,dimension(:,:,:),allocatable :: keygloc_work, keyg_glob_work
  integer,dimension(:,:),allocatable :: keyv_glob_work
  !$ integer :: omp_get_thread_num
  !integer, allocatable :: keygloc(:,:)


  call f_routine('segkeys_Sphere')

  !keygloc = f_malloc((/ 2, nseg /),id='keygloc')

  !dimensions of the localisation region (O:nIl)
  ! must be smaller or equal to simulation box dimensions
  !n1l=i1ec-i1sc
  !n2l=i2ec-i2sc
  !n3l=i3ec-i3sc
  n1l=nu1-nl1
  n2l=nu2-nl2
  n3l=nu3-nl3


  ! For perdiodic boundary conditions, one has to check also in the neighboring
  ! cells (see in the loop below)
  if (perx) then
     ijs1 = -1
     ije1 = 1
  else
     ijs1 = 0
     ije1 = 0
  end if
  if (pery) then
     ijs2 = -1
     ije2 = 1
  else
     ijs2 = 0
     ije2 = 0
  end if
  if (perz) then
     ijs3 = -1
     ije3 = 1
  else
     ijs3 = 0
     ije3 = 0
  end if


!!! Count the number of segments which must be checked from the point of view of the z coordinate
  !!nseg_check = nseg_check + 1
  !!do iseg=1,nsegglob
  !!    j0=keygglob(1,iseg)
  !!    j1=keygglob(2,iseg)
  !!    ii=j0-1
  !!    i3=ii/np
  !!    ii3=i3+nl3glob

  !!    inside=.false.
  !!    do ij3=ijs3,ije3
  !!        jj3=ii3+ij3*(n3+1)
  !!        dz=((jj3*hz)-locregCenter(3))**2
  !!        if(dz<=cut) then
  !!            inside=.true.
  !!        end if
  !!    end do
  !!    if (inside) then
  !!        nseg_check = nseg_check + 1
  !!    end if
  !!end do
  !!iseg_lookup = f_malloc(nseg_check,id='iseg_lookup')
  !!nseg_check = nseg_check + 1
  !!do iseg=1,nsegglob
  !!    j0=keygglob(1,iseg)
  !!    j1=keygglob(2,iseg)
  !!    ii=j0-1
  !!    i3=ii/np
  !!    ii3=i3+nl3glob

  !!    inside=.false.
  !!    do ij3=ijs3,ije3
  !!        jj3=ii3+ij3*(n3+1)
  !!        dz=((jj3*hz)-locregCenter(3))**2
  !!        if(dz<=cut) then
  !!            inside=.true.
  !!        end if
  !!    end do
  !!    if (inside) then
  !!        iseg_check = iseg_check + 1
  !!        iseg_lookup(nseg_check) = iseg
  !!    end if
  !!end do



  ! Count how many global segments have an overlap from the viewpoint of the z coordinate
  nsegglob_start = huge(nsegglob_start)
  nsegglob_end = 0
  cut=locrad**2
  n1p1=n1+1
  np=n1p1*(n2+1)
  !$omp parallel if (nsegglob>1000) &
  !$omp default(none) &
  !$omp shared(nsegglob, keygglob, np, nl3glob, ijs3, ije3, n3) &
  !$omp shared(hz, locregCenter, cut, nsegglob_start, nsegglob_end) &
  !$omp private(iseg, j0, j1, ii, i3, ii3, inside, ij3, jj3, dz)
  !$omp do reduction(min: nsegglob_start) reduction(max: nsegglob_end)
  do iseg=1,nsegglob
     j0=keygglob(1,iseg)
     j1=keygglob(2,iseg)
     ii=j0-1
     i3=ii/np
     ii3=i3+nl3glob

     inside=.false.
     do ij3=ijs3,ije3
        jj3=ii3+ij3*(n3+1)
        !this is not trivial for non-orthorhombic case
        !should employ mesh-like quantities
        dz=((jj3*hz)-locregCenter(3))**2
        if(dz<=cut) then
           inside=.true.
        end if
     end do
     if (inside) then
        nsegglob_start = min(nsegglob_start,iseg)
        nsegglob_end = max(nsegglob_end,iseg)
     end if
  end do
  !$omp end do
  !$omp end parallel




  call distribute_on_threads(nsegglob_start, nsegglob_end, nthread, ise)

  keygloc_work = f_malloc((/1.to.2,1.to.nseg,0.to.nthread-1/),id='keygloc_work')
  keyg_glob_work = f_malloc((/1.to.2,1.to.nseg,0.to.nthread-1/),id='keyg_glob_work')
  keyv_glob_work = f_malloc((/1.to.nseg,0.to.nthread-1/),id='keyv_glob_work')
  keyv_last = f_malloc(0.to.nthread-1,id='keyv_last')


  nstartarr = f_malloc(0.to.nthread-1,id='nstartarr')

  !can add openmp here too as segment always ends at end of y direction? 
  !problem is need nend value - can do a pre-scan to find seg value only as with init_collcom.
  !for now just do omp section
  !!cut=locrad**2
  !!n1p1=n1+1
  !!np=n1p1*(n2+1)
  n1lp1=n1l+1
  nlp=n1lp1*(n2l+1)
  ivctr=0
  jvctr=0
  kvctr=0
  nvctr=0
  nstart=0
  nend=0
  ivctr_tot = 0
  jvctr_tot = 0
  nstart_tot = 0
  nend_tot = 0
  segment=.false.
  ithread = 0
  !$omp parallel &
  !$omp default(none) &
  !$omp shared(ise, hx, hy, hz, keygglob, np, n1p1, nl1glob, nl2glob, nl3glob, locregCenter, nsegglob) &
  !$omp shared(keygloc_work, keyg_glob_work, keyv_glob_work, nstartarr, nl1, nl2, nl3, nu1, nu2, nu3) &
  !$omp shared(ijs3, ije3, ijs2, ije2, ijs1, ije1, n1, n2, n3, cut, n1lp1, nlp, nthread) &
  !$omp shared(keygloc, keyg_glob, keyv_glob, ivctr_tot, jvctr_tot, nstart_tot, nend_tot, keyv_last) &
  !$omp firstprivate(ithread, ivctr, jvctr, kvctr, nvctr, nstart, nend, segment) &
  !$omp private(iseg, j0, j1, ii, i3, i2, i1, i0, ii2, ii3, dz, dy, igridgloba, jj1) &
  !$omp private(i, ii1, dx, i1l, igridglob, inside, ij3, jj3, ij2, jj2, ij1, i2l, i3l) &
  !$omp private(ii1mod, ii2mod, ii3mod, igridpoint, offset, j, kthread,jthread)
  !jj1, )
  !$ ithread = omp_get_thread_num()
  do iseg=ise(1,ithread),ise(2,ithread)
     !!omp do schedule(dynamic,50)
     !do iseg=1,nsegglob
     j0=keygglob(1,iseg)
     j1=keygglob(2,iseg)
     ii=j0-1
     i3=ii/np
     ii=ii-i3*np
     i2=ii/n1p1
     i0=ii-i2*n1p1
     i1=i0+j1-j0
     ii2=i2+nl2glob
     ii3=i3+nl3glob

!!! First just check the z dimension. If inside is false, proceed directly,
!!! otherwise check also the other dimensions.
     !!inside=.false.
     !!do ij3=ijs3,ije3!-1,1
     !!    jj3=ii3+ij3*(n3+1)
     !!    dz=((jj3*hz)-locregCenter(3))**2
     !!    if(dz<=cut) then
     !!        inside=.true.
     !!    end if
     !!end do
     !!check_z_if: if (inside) then
     ! May be inside the sphere, so check also the other dimensions.
     ! Since each line in y (and thus also each plane in the z dimensions) starts
     ! a new segment, the following does not have to be done.
     inside = .false.

     igridgloba=ii3*np+ii2*n1p1+1 
     do i=i0,i1
        ii1=i+nl1glob
        dx=((ii1*hx)-locregCenter(1))**2
        i1l=ii1-nl1
        !igridpoint=igridpointa+i1l
        igridglob=igridgloba+ii1 
        inside=.false.
        do ij3=ijs3,ije3!-1,1
           jj3=ii3+ij3*(n3+1)
           dz=((jj3*hz)-locregCenter(3))**2
           do ij2=ijs2,ije2!-1,1
              jj2=ii2+ij2*(n2+1)
              dy=((jj2*hy)-locregCenter(2))**2
              do ij1=ijs1,ije1!-1,1
                 jj1=ii1+ij1*(n1+1)
                 dx=((jj1*hx)-locregCenter(1))**2
                 if(dx+dy+dz<=cut) then
                    if (inside) call f_err_throw('twice inside',err_name='BIGDFT_RUNTIME_ERROR')
                    inside=.true.
                    ii1mod=jj1
                    ii2mod=jj2
                    ii3mod=jj3
                    i1l=jj1-nl1
                    i2l=jj2-nl2
                    i3l=jj3-nl3
                    igridpoint=i3l*nlp+i2l*n1lp1+i1l+1
                    !write(*,'(a,4i8)') 'i1l, i2l, i3l, igridpoint', i1l, i2l, i3l, igridpoint
                 end if
              end do
           end do
        end do
        !write(*,*) 'ii1, ii2, ii3, inside', ii1, ii2, ii3, inside
        if(inside) then
           ! Check that we are not outside of the locreg region
           ivctr=ivctr+1
           kvctr=kvctr+1
           !write(*,*) 'inside: kvctr, igridpoint', kvctr, igridpoint
           if(ii1mod<nl1) then
              write(*,'(a,i0,a,i0,a)') 'ERROR: ii1mod=',ii1mod,'<',nl1,'=nl1'
              stop
           end if
           if(ii2mod<nl2) then
              write(*,'(a,i0,a,i0,a)') 'ERROR: ii2mod=',ii2mod,'<',nl2,'=nl2'
              stop
           end if
           if(ii3mod<nl3) then
              write(*,'(a,i0,a,i0,a)') 'ERROR: ii3mod=',ii3mod,'<',nl3,'=nl3'
              stop
           end if
           if(ii1mod>nu1) then
              write(*,'(a,i0,a,i0,a)') 'ERROR: ii1mod=',ii1mod,'>',nu1,'=nu1'
              stop
           end if
           if(ii2mod>nu2) then
              write(*,'(a,i0,a,i0,a)') 'ERROR: ii2mod=',ii2mod,'>',nu2,'=nu2'
              stop
           end if
           if(ii3mod>nu3) then
              write(*,'(a,i0,a,i0,a)') 'ERROR: ii3mod=',ii3mod,'>',nu3,'=nu3'
              stop
           end if
           nvctr=nvctr+1
           if(.not.segment) then
              nstart=nstart+1
              keygloc_work(1,nstart,ithread)=igridpoint
              keyg_glob_work(1,nstart,ithread)=igridglob
              keyv_glob_work(nstart,ithread)=nvctr
              segment=.true.
           end if
        else
           if(segment) then
              nend=nend+1
              keygloc_work(2,nend,ithread)=igridpoint!-1
              keyg_glob_work(2,nend,ithread)=igridglob-1
              !write(*,'(a,4i7)') 'outside: kvctr, igridpoint, keygloc(1:2,nend)', kvctr, igridpoint, keygloc(1:2,nend)
              segment=.false.
              jvctr=jvctr+keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1
              if (kvctr/=keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1) then
                 write(*,*) 'kvctr, keygloc(2,nend)-keygloc(1,nend)+1', &
                      kvctr, keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1
                 stop 'kvctr/=keygloc(2,nend)-keygloc(1,nend)+1'
              end if
              kvctr=0
           end if
        end if
     end do
     if(segment) then
        ! Close the segment
        nend=nend+1
        keygloc_work(2,nend,ithread)=igridpoint
        keyg_glob_work(2,nend,ithread)=igridglob
        segment=.false.
        jvctr=jvctr+keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1
        if (kvctr/=keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1) then
           write(*,*) 'kvctr, keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1', &
                kvctr, keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1
           stop 'kvctr/=keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1'
        end if
        kvctr=0
     end if
     !!end if check_z_if
  end do
  ! Some checks
  if (nstart/=nend) call f_err_throw('nstart/=nend',err_name='BIGDFT_RUNTIME_ERROR')
  ! Number of segments calculated by ithread
  nstartarr(ithread) = nstart
  ! Number of elements calculated by ithread
  if (nstart>0) then
     keyv_last(ithread) = keyv_glob_work(nstart,ithread)+keyg_glob_work(2,nstart,ithread)-keyg_glob_work(1,nstart,ithread)
  else
     keyv_last(ithread) = 0
  end if
  !$omp barrier
  ii = 1
  do jthread=0,nthread-1
     if (ithread==jthread) then
        if (nstartarr(jthread)>0) then
           call f_memcpy(n=2*nstartarr(jthread), src=keygloc_work(1,1,ithread), dest=keygloc(1,ii))
           call f_memcpy(n=2*nstartarr(jthread), src=keyg_glob_work(1,1,ithread), dest=keyg_glob(1,ii))
           offset = 0
           do kthread=0,jthread-1
              offset = offset + keyv_last(kthread)
           end do
           do j=1,nstartarr(jthread)
              keyv_glob(ii+j-1) = keyv_glob_work(j,ithread) + offset
           end do
           !call f_memcpy(n=nstartarr(jthread), src=keyv_glob_work(1,ithread), dest=keyv_glob(ii))
        end if
     end if
     ii = ii + nstartarr(jthread)
  end do

  !$omp critical
  ivctr_tot = ivctr_tot + ivctr
  jvctr_tot = jvctr_tot + jvctr
  nstart_tot = nstart_tot + nstart
  nend_tot = nend_tot + nend
  !nseg_tot = nseg_tot + nseg
  !$omp end critical
  !$omp end parallel

  !write(*,*) 'nstartarr',nstartarr
  !do ii=1,nseg
  !    write(*,*) 'ii, keygloc(:,ii)', ii, keygloc(:,ii)
  !    write(*,*) 'ii, keyg_glob(:,ii)', ii, keyg_glob(:,ii)
  !    write(*,*) 'ii, keyv_glob(ii)', ii, keyv_glob(ii)
  !end do


  ! Some checks
  if (ivctr_tot/=nvctr_loc) then
     write(*,*) 'ivctr_tot, nvctr_loc', ivctr_tot, nvctr_loc
     stop 'ivctr_tot/=nvctr_loc'
  end if

  if (jvctr_tot/=nvctr_loc) then
     write(*,*) 'jvctr_tot, nvctr_loc', jvctr_tot, nvctr_loc
     stop 'jvctr_tot/=nvctr_loc'
  end if

  if (nend_tot /= nstart_tot) then
     write(*,*) 'nend_tot , nstart_tot',nend_tot,nstart_tot
     stop 'nend_tot <> nstart_tot'
  endif
  if (nseg /= nstart_tot) then
     write(*,*) 'nseg , nstart_tot',nseg,nstart_tot
     stop 'nseg <> nstart_tot'
  endif

  ! Now build the keyvloc where we replace the segments in order for the loc
  ivctr=0
  ii = maxval(keygloc)
  do iseg=1,nseg
     !sorting the keyg_loc
     loc = minloc(keygloc(1,:),1)
     keyg_loc(1,iseg) = keygloc(1,loc)
     keyg_loc(2,iseg) = keygloc(2,loc)
     !    print *,'iseg,keygloc,keyg_loc',iseg,keygloc(1,loc),keygloc(2,loc),keyg_loc(1,iseg),keyg_loc(2,iseg)
     keyv_loc(iseg) = keyv_glob(loc)
     !keygloc(1,loc) = maxval(keygloc) + 1
     keygloc(1,loc) = ii+iseg !just put to the maximal value
     !write(*,'(a,7i8)') 'iseg,keyglob,keyvglob,keygloc,keyvloc',iseg,keyg_glob(1:2,iseg),keyv_glob(iseg),keyg_loc(1:2,iseg),keyv_loc(iseg)
     ivctr=ivctr+keyg_loc(2,iseg)-keyg_loc(1,iseg)+1
  end do
  !call f_free(keygloc)
  if (ivctr/=nvctr_loc) then
     write(*,*) 'ivctr, nvctr_loc', ivctr, nvctr_loc
     stop 'rearrangement check: ivctr/=nvctr_loc'
  end if

  ! Some checks
  ivctr=0
  !write(*,*) 'nlp, n1lp1', nlp, n1lp1
  !$omp parallel &
  !$omp default(none) &
  !$omp shared(nseg, keyg_loc, nlp, n1lp1, n1l, n2l, n3l, ivctr) &
  !$omp private(iseg, j0, j1, ii, i3, i2, i1, i0, i)
  !$omp do reduction(+:ivctr)
  do iseg=1,nseg
     j0=keyg_loc(1,iseg)
     j1=keyg_loc(2,iseg)
     ii=j0-1
     i3=ii/nlp
     ii=ii-i3*nlp
     i2=ii/n1lp1
     i0=ii-i2*n1lp1
     i1=i0+j1-j0
     !if (i2<nl2) then
     !    write(*,'(a,2(i0,a))') 'ERROR: i2=',i2,'<',nl2,'=nl2' ; stop
     !end if
     if (i2>n2l) then
        write(*,'(a,2(i0,a))') 'ERROR: i2=',i2,'>',n2l,'=n2l' ; stop
     end if
     !if (i3<nl3) then
     !    write(*,'(a,2(i0,a))') 'ERROR: i3=',i3,'<',nl3,'=nl3' ; stop
     !end if
     if (i3>n3l) then
        write(*,'(a,2(i0,a))') 'ERROR: i3=',i3,'>',n3l,'=n3l' ; stop
     end if
     do i=i0,i1
        ivctr=ivctr+1
        !write(*,'(a,6i8)') 'j0, j1, ii, i, i2, i3', j0, j1, ii, i, i2, i3
        !if (i<nl1) then
        !    write(*,'(a,2(i0,a))') 'ERROR: i=',i,'<',nl1,'=nl1' ; stop
        !end if
        if (i>n1l) then
           write(*,'(a,2(i0,a))') 'ERROR: i=',i,'>',n1l,'=n1l' ; stop
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel

  if (ivctr/=nvctr_loc) then
     write(*,*) 'ivctr, nvctr_loc', ivctr, nvctr_loc
     stop 'second check: ivctr/=nvctr_loc'
  end if

  call f_free(keygloc_work)
  call f_free(keyg_glob_work)
  call f_free(keyv_glob_work)
  call f_free(keyv_last)
  call f_free(nstartarr)
  call f_free_ptr(ise)

  call f_release_routine()

END SUBROUTINE segkeys_Sphere


!> Determines the the wavefunction descriptors,wfd, and fine grid upper limit of locreg 
!! taking into account the pediodicity
!!          
!! @warning
!!    We assign Llr%nfl and llr%nfu with respect to the origin of the local zone, like in determine_locreg. 
subroutine determine_wfdSphere(ilr,nlr,Glr,hx,hy,hz,Llr)!,outofzone)

  use module_base
  use locregs, only: allocate_wfd,locreg_descriptors
  use at_domain, only: domain_periodic_dims
  implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in) :: ilr,nlr
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  real(kind=8),intent(in) :: hx, hy, hz
  type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr  ! Localization grid descriptors 

  !Subroutine Array Arguments
  !  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr

  !local variables
  integer :: ii
  integer,dimension(3) :: Gife,Gifs,iedir,isdir,Lifs,Life,period,Gics,Gice
  character(len=*), parameter :: subname='determine_wfdSphere'
  !!  integer :: nseg_c,nseg_f,nvctr_c,nvctr_f      ! total number of sgements and elements
  integer, allocatable :: keygloc_tmp(:,:)
  logical :: perx, pery, perz
  logical, dimension(3) :: peri

  call f_routine(id=subname)

  ! periodicity in the three directions
!!$  perx=(glr%geocode /= 'F')
!!$  pery=(glr%geocode == 'P')
!!$  perz=(glr%geocode /= 'F')

  peri=domain_periodic_dims(glr%mesh%dom)
  perx=peri(1)
  pery=peri(2)
  perz=peri(3)

  !starting point of locreg (can be outside the simulation box)
  isdir(1) = Llr(ilr)%ns1
  isdir(2) = Llr(ilr)%ns2
  isdir(3) = Llr(ilr)%ns3
  !ending point of locreg (can be outside the simulation box)
  iedir(1) = Llr(ilr)%ns1 + Llr(ilr)%d%n1
  iedir(2) = Llr(ilr)%ns2 + Llr(ilr)%d%n2
  iedir(3) = Llr(ilr)%ns3 + Llr(ilr)%d%n3
  ! starting and ending point of coarse grid in Global region
  Gics(1) = Glr%ns1
  Gics(2) = Glr%ns2
  Gics(3) = Glr%ns3
  Gice(1) = Glr%ns1 + Glr%d%n1
  Gice(2) = Glr%ns2 + Glr%d%n2
  Gice(3) = Glr%ns3 + Glr%d%n3
  ! starting and ending point of fine grid in Global region
  Gifs(1) = Glr%d%nfl1 + Glr%ns1
  Gifs(2) = Glr%d%nfl2 + Glr%ns2
  Gifs(3) = Glr%d%nfl3 + Glr%ns3
  Gife(1) = Glr%d%nfu1 + Glr%ns1
  Gife(2) = Glr%d%nfu2 + Glr%ns2
  Gife(3) = Glr%d%nfu3 + Glr%ns3
  ! periodicity
  period(1) = Glr%d%n1+1
  period(2) = Glr%d%n2+1
  period(3) = Glr%d%n3+1

!!! Determine starting point of the fine grid in locreg
  !!do ii=1,3
  !!   if (Llr(ilr)%outofzone(ii) > 0) then
  !!      ! When periodicity, we must check for 2 different situations:
  !!      ! (1) : starting of locreg before or in fine grid zone
  !!      if (isdir(ii) < Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
  !!      ! (2) : starting point after fine grid
  !!      if (isdir(ii) > Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii)+period(ii))-isdir(ii)
  !!   else
  !!       Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
  !!   end if 
  !!end do

!!! Determine ending point of the fine grid in locreg
  !!do ii=1,3
  !!   if(Llr(ilr)%outofzone(ii) > 0) then
  !!      !When periodicity, we must check for three different situations:
  !!      ! (1) : ending of locreg before fine grid zone
  !!      if(iedir(ii) < (Gifs(ii) + period(ii))) Life(ii) = Gife(ii)-isdir(ii)
  !!      ! (2) : ending of locreg in fine grid zone
  !!      if(iedir(ii) > (Gifs(ii) + period(ii)) .and. iedir(ii) < (Gife(ii) + period(ii))) then
  !!        Life(ii) = iedir(ii)-isdir(ii)
  !!      end if
  !!      ! (3) : ending of locreg after ending of fine grid zone
  !!      if(iedir(ii) > (Gife(ii)+period(ii))) Life(ii) = Gife(ii) + period(ii)-isdir(ii)
  !!   else
  !!      Life(ii) = min(iedir(ii),Gife(ii))-isdir(ii)
  !!   end if
  !!end do

  do ii=1,3
     ! Determine starting point of the fine grid in locreg. There are two possibilities:
     if (isdir(ii)<gics(ii)) then
        ! Start of the locreg locreg outside of the global box
        lifs(ii) = max(isdir(ii)+period(ii),gifs(ii)) - period(ii) - isdir(ii)
     else if(isdir(ii)>=gics(ii)) then
        ! Start of locreg inside of the global box
        lifs(ii) = max(isdir(ii),gifs(ii)) - isdir(ii)
     else
        stop 'cannot determine start of fine grid'
     end if

     ! Determine ending point of the fine grid in locreg. There are two possibilities:
     if (iedir(ii)>gice(ii)) then
        ! End of the locreg outside of the global box
        life(ii) = min(iedir(ii)-period(ii),gife(ii)) + period(ii) - isdir(ii)
     else if(iedir(ii)<=gice(ii)) then
        ! End of the locreg inside of the global box
        life(ii) = min(iedir(ii),gife(ii)) - isdir(ii)
     else
        stop 'cannot determine start of fine grid'
     end if
  end do


  ! Assign values to Llr
  Llr(ilr)%d%nfl1 = Lifs(1)
  Llr(ilr)%d%nfl2 = Lifs(2)
  Llr(ilr)%d%nfl3 = Lifs(3)
  Llr(ilr)%d%nfu1 = Life(1)
  Llr(ilr)%d%nfu2 = Life(2)
  Llr(ilr)%d%nfu3 = Life(3)

  ! define the wavefunction descriptors inside the localisation region
!!!coarse part
  !!call num_segkeys_sphere(perx, pery, perz, Glr%d%n1, Glr%d%n2, Glr%d%n3, &
  !!     glr%ns1, glr%ns2, glr%ns3, &
  !!     hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
  !!     Glr%wfd%nseg_c, Glr%wfd%keygloc, &
  !!     Glr%wfd%keyvloc, &
  !!     llr(ilr)%wfd%nseg_c, llr(ilr)%wfd%nvctr_c)

!!!fine part
  !!call num_segkeys_sphere(perx, pery, perz, Glr%d%n1, Glr%d%n2, Glr%d%n3, &
  !!     glr%ns1, glr%ns2, glr%ns3, &
  !!     hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
  !!     glr%wfd%nseg_f, Glr%wfd%keygloc(1:,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):), &
  !!     Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
  !!     llr(ilr)%wfd%nseg_f, llr(ilr)%wfd%nvctr_f)
  call get_num_segkeys(perx, pery, perz, glr%d%n1, glr%d%n2, glr%d%n3, &
       glr%ns1, glr%ns2, glr%ns3, &
       hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
       glr%wfd%nseg_c, glr%wfd%nseg_f, glr%wfd%keygloc,  Glr%wfd%keyvloc, &
       llr(ilr)%wfd%nseg_c, llr(ilr)%wfd%nvctr_c, llr(ilr)%wfd%nseg_f, llr(ilr)%wfd%nvctr_f)

  !write(*,'(a,2i8)') 'llr(ilr)%wfd%nvctr_c, llr(ilr)%wfd%nvctr_f', llr(ilr)%wfd%nvctr_c, llr(ilr)%wfd%nvctr_f

  !allocate the wavefunction descriptors following the needs
  call allocate_wfd(Llr(ilr)%wfd)

  !Now, fill the descriptors:

  !keygloc_tmp = f_malloc((/ 2, (llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f) /),id='keygloc_tmp')

!!$omp parallel default(private) &
!!$omp shared(Glr,llr,hx,hy,hz,ilr,keygloc_tmp,perx,pery,perz)  
!!$omp sections
!!$omp section

!!!!coarse part
!!!call segkeys_Sphere(perx, pery, perz, Glr%d%n1, Glr%d%n2, Glr%d%n3, &
!!!     glr%ns1, glr%ns2, glr%ns3, &
!!!     llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
!!!     llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
!!!     llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
!!!     llr(ilr)%wfd%nseg_c, hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
!!!     Glr%wfd%nseg_c, Glr%wfd%keygloc(1:,1:), &
!!!     Glr%wfd%keyvloc(1:), llr(ilr)%wfd%nvctr_c, &
!!!     llr(ilr)%wfd%keygloc(1:,1:),llr(ilr)%wfd%keyglob(1:,1:), &
!!!     llr(ilr)%wfd%keyvloc(1:), llr(ilr)%wfd%keyvglob(1:), &
!!!     keygloc_tmp(1:,1:))

!!!!!$omp section
!!!!fine part
!!!call segkeys_Sphere(perx, pery, perz, Glr%d%n1, Glr%d%n2, Glr%d%n3, &
!!!     glr%ns1, glr%ns2, glr%ns3, &
!!!     llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
!!!     llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
!!!     llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
!!!     llr(ilr)%wfd%nseg_f, hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
!!!     Glr%wfd%nseg_f, Glr%wfd%keygloc(1:,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):),&
!!!     Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):), llr(ilr)%wfd%nvctr_f, &
!!!     llr(ilr)%wfd%keygloc(1:,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f):), &
!!!     llr(ilr)%wfd%keyglob(1:,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f):), &
!!!     llr(ilr)%wfd%keyvloc(llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f):), &
!!!     llr(ilr)%wfd%keyvglob(llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f):), &
!!!     keygloc_tmp(1:,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f):))  
!!!!!$omp end sections
!!!!!$omp end parallel
  call get_segkeys(perx, pery, perz, glr%d%n1, glr%d%n2, glr%d%n3, &
       glr%ns1, glr%ns2, glr%ns3, &
       llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
       llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
       llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
       hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
       glr%wfd%nseg_c, glr%wfd%nseg_f, glr%wfd%keygloc, glr%wfd%keyvloc, &
       llr(ilr)%wfd%nseg_c, llr(ilr)%wfd%nseg_f, llr(ilr)%wfd%nvctr_c, llr(ilr)%wfd%nvctr_f, &
       llr(ilr)%wfd%keygloc, llr(ilr)%wfd%keyglob, llr(ilr)%wfd%keyvloc, llr(ilr)%wfd%keyvglob)

  !call f_free(keygloc_tmp)

  call f_release_routine()


END SUBROUTINE determine_wfdSphere


subroutine get_num_segkeys(perx, pery, perz, n1, n2, n3, ns1, ns2, ns3, hx, hy, hz, locrad, &
     locregCenter, nseg_c_glob, nseg_f_glob, keyg_glob, keyv_glob, &
     nseg_c, nvctr_c, nseg_f, nvctr_f)
  use module_base
  implicit none

  ! Calling arguments
  logical,intent(in) :: perx, pery, perz
  integer,intent(in) :: n1, n2, n3, ns1, ns2, ns3, nseg_c_glob, nseg_f_glob
  real(kind=8),intent(in) :: hx, hy, hz, locrad
  real(kind=8),dimension(3),intent(in) :: locregCenter
  integer,dimension(2,nseg_c_glob+nseg_f_glob),intent(in) :: keyg_glob
  integer,dimension(nseg_c_glob+nseg_f_glob),intent(in) :: keyv_glob
  integer,intent(out) :: nseg_c, nvctr_c, nseg_f, nvctr_f

  call f_routine(id='get_num_segkeys')

  !coarse part
  call num_segkeys_sphere(perx, pery, perz, n1, n2, n3, &
       ns1, ns2, ns3, &
       hx, hy, hz, locrad, locregCenter, &
       nseg_c_glob, keyg_glob, keyv_glob, &
       nseg_c, nvctr_c)

  !fine part
  call num_segkeys_sphere(perx, pery, perz, n1, n2, n3, &
       ns1, ns2, ns3, &
       hx, hy, hz, locrad, locregCenter, &
       nseg_f_glob, keyg_glob(1,nseg_c_glob+min(1,nseg_f_glob)), &
       keyv_glob(nseg_c_glob+min(1,nseg_f_glob)), &
       nseg_f, nvctr_f)

  call f_release_routine()

end subroutine get_num_segkeys


subroutine get_segkeys(perx, pery, perz, &
     n1_glob, n2_glob, n3_glob, nl1_glob, nl2_glob, nl3_glob, &
     nl1, nu1, nl2, nu2, nl3, nu3, hx, hy, hz, locrad, locregCenter, &
     nseg_c_glob, nseg_f_glob, keyg_glob, keyv_glob, &
     nseg_c, nseg_f, nvctr_c, nvctr_f, &
     keygloc, keygglob, keyvloc, keyvglob)
  use module_base
  implicit none

  ! Calling arguments
  logical,intent(in) :: perx, pery, perz
  integer,intent(in) :: n1_glob, n2_glob, n3_glob, nl1_glob, nl2_glob, nl3_glob
  integer,intent(in) :: nl1, nl2, nl3, nu1, nu2, nu3
  integer,intent(in) :: nseg_c_glob, nseg_c, nseg_f_glob, nseg_f
  integer,intent(in) :: nvctr_c, nvctr_f
  real(kind=8),intent(in) :: hx, hy, hz, locrad
  real(kind=8),dimension(3),intent(in) :: locregCenter
  integer,dimension(2,nseg_c_glob+nseg_f_glob),intent(in) :: keyg_glob
  integer,dimension(nseg_c_glob+nseg_f_glob),intent(in) :: keyv_glob
  integer,dimension(2,nseg_c+nseg_f),intent(out) :: keygloc, keygglob
  integer,dimension(nseg_c+nseg_f),intent(out) :: keyvloc, keyvglob

  integer, allocatable :: keygloc_tmp(:,:)

  call f_routine(id='get_segkeys')

  keygloc_tmp = f_malloc((/2,nseg_c+nseg_f/),id='keygloc_tmp')

  !coarse part
  call segkeys_Sphere(perx, pery, perz, n1_glob, n2_glob, n3_glob, &
       nl1_glob, nl2_glob, nl3_glob, &
       nl1, nu1, nl2, nu2, nl3, nu3, &
       nseg_c, hx, hy, hz, locrad, locregCenter, &
       nseg_c_glob, keyg_glob(1,1), &
       keyv_glob(1), nvctr_c, &
       keygloc(1,1),keygglob(1,1), &
       keyvloc(1), keyvglob(1), &
       keygloc_tmp(1,1))

  !fine part
  call segkeys_Sphere(perx, pery, perz, n1_glob, n2_glob, n3_glob, &
       nl1_glob, nl2_glob, nl3_glob, &
       nl1, nu1, nl2, nu2, nl3, nu3, &
       nseg_f, hx, hy, hz, locrad, locregCenter, &
       nseg_f_glob, keyg_glob(1,nseg_c_glob+min(1,nseg_f_glob)),&
       keyv_glob(nseg_c_glob+min(1,nseg_f_glob)), nvctr_f, &
       keygloc(1,nseg_c+min(1,nseg_f)), &
       keygglob(1,nseg_c+min(1,nseg_f)), &
       keyvloc(nseg_c+min(1,nseg_f)), &
       keyvglob(nseg_c+min(1,nseg_f)), &
       keygloc_tmp(1,nseg_c+min(1,nseg_f)))  

  call f_free(keygloc_tmp)

  call f_release_routine()

end subroutine get_segkeys

!> Determines the the wavefunction descriptors,wfd, and fine grid upper limit of locreg 
!! taking into account the pediodicity
!! @warning
!!    We assign Llr%nfl and llr%nfu with respect to the origin of the local zone, like in determine_locreg. 
subroutine determine_wfd_periodicity(ilr,nlr,Glr,Llr)!,outofzone)

  use module_base
  use locregs, only: allocate_wfd,locreg_descriptors
  implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in) :: ilr,nlr
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr  ! Localization grid descriptors 

  !Subroutine Array Arguments
  !  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr

  !local variables
  integer :: ii
  integer,dimension(3) :: Gife,Gifs,iedir,isdir,Lifs,Life,period
  integer :: nseg_c,nseg_f,nvctr_c,nvctr_f      ! total number of sgements and elements
  character(len=*), parameter :: subname='determine_wfd_periodicity'

  call f_routine(id='determine_wfd_periodicity')

  !starting point of locreg (always inside locreg)
  isdir(1) = Llr(ilr)%ns1
  isdir(2) = Llr(ilr)%ns2
  isdir(3) = Llr(ilr)%ns3
  !ending point of locreg (can be outside the simulation box)
  iedir(1) = Llr(ilr)%ns1 + Llr(ilr)%d%n1
  iedir(2) = Llr(ilr)%ns2 + Llr(ilr)%d%n2
  iedir(3) = Llr(ilr)%ns3 + Llr(ilr)%d%n3
  ! starting and ending point of fine grid in Global region
  Gifs(1) = Glr%d%nfl1 + Glr%ns1
  Gifs(2) = Glr%d%nfl2 + Glr%ns2
  Gifs(3) = Glr%d%nfl3 + Glr%ns3
  Gife(1) = Glr%d%nfu1 + Glr%ns1
  Gife(2) = Glr%d%nfu2 + Glr%ns2
  Gife(3) = Glr%d%nfu3 + Glr%ns3
  ! periodicity
  period(1) = Glr%d%n1
  period(2) = Glr%d%n2
  period(3) = Glr%d%n3

  ! Determine starting point of the fine grid in locreg
  do ii=1,3
     if (Llr(ilr)%outofzone(ii) > 0) then
        ! When periodicity, we must check for 2 different situations:
        ! (1) : starting of locreg before or in fine grid zone
        if (isdir(ii) < Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
        ! (2) : starting point after fine grid
        if (isdir(ii) > Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii)+period(ii))-isdir(ii)
     else
        Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
     end if
  end do

  ! Determine ending point of the fine grid in locreg
  do ii=1,3
     if(Llr(ilr)%outofzone(ii) > 0) then
        !When periodicity, we must check for three different situations:
        ! (1) : ending of locreg before fine grid zone
        if(iedir(ii) < (Gifs(ii) + period(ii))) Life(ii) = Gife(ii)-isdir(ii)
        ! (2) : ending of locreg in fine grid zone
        if(iedir(ii) > (Gifs(ii) + period(ii)) .and. iedir(ii) < (Gife(ii) + period(ii))) then
           Life(ii) = iedir(ii)-isdir(ii)
        end if
        ! (3) : ending of locreg after ending of fine grid zone
        if(iedir(ii) > (Gife(ii)+period(ii))) Life(ii) = Gife(ii) + period(ii)-isdir(ii)
     else
        Life(ii) = min(iedir(ii),Gife(ii))-isdir(ii)
     end if
  end do

  ! Assign values to Llr
  Llr(ilr)%d%nfl1 = Lifs(1)
  Llr(ilr)%d%nfl2 = Lifs(2)
  Llr(ilr)%d%nfl3 = Lifs(3)
  Llr(ilr)%d%nfu1 = Life(1)
  Llr(ilr)%d%nfu2 = Life(2)
  Llr(ilr)%d%nfu3 = Life(3)

  ! define the wavefunction descriptors inside the localisation region
  !coarse part
  call num_segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),isdir(2),&
       iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_c,Glr%wfd%nvctr_c,&
       Glr%wfd%keygloc(1:,1:),Glr%wfd%keyvloc(1:),nseg_c,nvctr_c,Llr(ilr)%outofzone(:))
  !fine part
  call num_segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),isdir(2),&
       iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
       Glr%wfd%keygloc(1:,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):),&
       Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):),nseg_f,nvctr_f,Llr(ilr)%outofzone(:))

  ! Assign the values to Llr
  Llr(ilr)%wfd%nseg_c = nseg_c
  Llr(ilr)%wfd%nseg_f = nseg_f
  Llr(ilr)%wfd%nvctr_c= nvctr_c
  Llr(ilr)%wfd%nvctr_f= nvctr_f

  !allocate the wavefunction descriptors following the needs
  call allocate_wfd(Llr(ilr)%wfd)

  !Now, fill the descriptors:
  !coarse part
  call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
       isdir(2),iedir(2),isdir(3),iedir(3),&
       Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1:,1:),Glr%wfd%keyvloc(1:),&
       Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nvctr_c,&
       Llr(ilr)%wfd%keygloc(1:,1:),Llr(ilr)%wfd%keyglob(1:,1:),Llr(ilr)%wfd%keyvloc(1:),&
       Llr(ilr)%wfd%keyvglob(1:),&
       Llr(ilr)%outofzone(:))

  !fine part
  call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
       isdir(2),iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
       Glr%wfd%keygloc(1:,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):),&
       Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):),&
       Llr(ilr)%wfd%nseg_f,Llr(ilr)%wfd%nvctr_f,&
       Llr(ilr)%wfd%keygloc(1:,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f):),&
       Llr(ilr)%wfd%keyglob(1:,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f):),&
       Llr(ilr)%wfd%keyvloc(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f):),&
       Llr(ilr)%wfd%keyvglob(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f):),&
       Llr(ilr)%outofzone(:))

  call f_release_routine()

END SUBROUTINE determine_wfd_periodicity


!> Calculates the number of segments and elements in localisation region
subroutine num_segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,outofzone)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, intent(out) :: nseg_loc,nvctr_loc
  integer, dimension(3),intent(in) :: outofzone
  !local variables
  logical :: lseg,go1,go2,go3
  integer :: iseg,j0,j1,ii,i1,i2,i3,i0,i,nsrt,nend,nvctr_check,n1p1,np

  nvctr_loc=0
  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0

  n1p1=n1+1
  np=n1p1*(n2+1)
  do iseg=1,nseg
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/np
     ii=ii-i3*np
     i2=ii/n1p1
     i0=ii-i2*n1p1
     i1=i0+j1-j0
     lseg=.false.
     ! overlap conditions if zone completely inside simulation box
     go2 = (i2sc <= i2 .and. i2 <= i2ec)
     go3 = (i3sc <= i3 .and. i3 <= i3ec)
     ! overlap conditions if zone as components in other periodic cells
     if(outofzone(2) > 0) go2 = (i2 <= outofzone(2) .or. i2 >= i2sc)
     if(outofzone(3) > 0) go3 = (i3 <= outofzone(3) .or. i3 >= i3sc)

     do i=i0,i1
        nvctr_check=nvctr_check+1
        go1 = (i1sc <= i .and. i <= i1ec)
        if(outofzone(1) > 0) go1 = (i <= outofzone(1) .or. i >= i1sc)

        if (go1 .and. go2 .and. go3 ) then
           nvctr_loc=nvctr_loc+1
           if (.not. lseg) then
              nsrt=nsrt+1
           end if
           lseg=.true.
        else
           if (lseg) then
              nend=nend+1
              lseg=.false.
           end if
        end if
     end do
     if (lseg) then
        nend=nend+1
     end if
  end do
  nseg_loc=nend

  !check
  if (nend /= nsrt) then
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif

  if (nvctr_check /= nvctr) then
     write(*,'(1x,a,2(i8))')&
          'ERROR: incorrect number of coarse points examined for reducing the localisation region',&
          nvctr_check,nvctr
     stop
  end if

END SUBROUTINE num_segkeys_periodic


subroutine segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,keygloc,keyglob,keyvloc,keyvglob,outofzone)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,nseg_loc,nvctr_loc
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, dimension(3), intent(in) :: outofzone
  integer, dimension(nseg_loc), intent(out) :: keyvglob
  integer, dimension(nseg_loc), intent(out) :: keyvloc
  integer, dimension(2,nseg_loc), intent(out) :: keygloc
  integer, dimension(2,nseg_loc), intent(out) :: keyglob
  !local variables
  character(len=*),parameter :: subname = 'segkeys_periodic'
  logical :: go1,go2,go3,lseg
  integer :: iseg,j0,j1,ii,i1,i2,i3,i0,i,nsrt,nend,nvctr_check,n1l,n2l,n3l,i1l,i2l,i3l,n1p1,np,n1lp1,nlp
  integer :: ngridp,ngridlob,loc
  integer, allocatable :: keyg_loc(:,:)

  call f_routine('segkeys_periodic')

  !should be initialized
  ngridp=-1000
  ngridlob=-1000

  !dimensions of the localisation region (O:nIl)
  ! must be smaller or equal to simulation box dimensions
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc


  keyg_loc = f_malloc((/ 2, nseg_loc /),id='keyg_loc')

  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0
  n1p1=n1+1
  np=n1p1*(n2+1)
  n1lp1=n1l+1
  nlp=n1lp1*(n2l+1)
  do iseg=1,nseg
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/np
     ii=ii-i3*np
     i2=ii/n1p1
     i0=ii-i2*n1p1
     i1=i0+j1-j0
     lseg=.false.

     ! intersection condition if zone inside simulation box
     go2 = (i2sc <= i2 .and. i2 <= i2ec)
     go3 = (i3sc <= i3 .and. i3 <= i3ec)
     ! intersection condition if zone has components outside simulation box (periodic)
     if(outofzone(2) > 0) go2 = (i2 <= outofzone(2) .or. i2 >= i2sc)
     if(outofzone(3) > 0) go3 = (i3 <= outofzone(3) .or. i3 >= i3sc)

     do i=i0,i1
        go1 = (i1sc <= i .and. i <= i1ec)
        if(outofzone(1) > 0) go1 = (i <= outofzone(1) .or. i >= i1sc)
        if (go1 .and. go2 .and. go3) then
           !index of the compressed function
           i1l=i-i1sc
           if(outofzone(1) > 0 .and. i <= outofzone(1))i1l = i - i1sc + n1 + 1
           i2l=i2-i2sc
           if(outofzone(2) > 0 .and. i2 <= outofzone(2))i2l = i2 - i2sc + n2 + 1
           i3l=i3-i3sc
           if(outofzone(3) > 0 .and. i3 <= outofzone(3))i3l = i3 - i3sc + n3 + 1
           ngridp=i3l*nlp + i2l*n1lp1 + i1l+1
           ngridlob = i3 * np + i2 * n1p1 + i + 1

           nvctr_check=nvctr_check+1
           if (.not. lseg) then
              !             print *,'         check:',i,i2,i3,i1l,i2l,i3l,ngridp
              nsrt=nsrt+1
              keyg_loc(1,nsrt)=ngridp
              keyglob(1,nsrt)=ngridlob
              keyvglob(nsrt)=nvctr_check
           end if
           lseg=.true.
        else
           if (lseg) then
              !              print *,'in        else:',i,i2,i3,i1l,i2l,i3l,ngridp
              nend=nend+1
              keyg_loc(2,nend)=ngridp
              keyglob(2,nend)=ngridlob
              lseg=.false.
           end if
        end if
     end do
     if (lseg) then
        !        print *,'in second else:',i,i2,i3,i1l,i2l,i3l,ngridp
        nend=nend+1
        keyg_loc(2,nend)=ngridp
        keyglob(2,nend)=ngridlob
     end if
  end do

  !check
  if (nvctr_check /= nvctr_loc .or. nend /= nsrt .or. nend /= nseg_loc) then
     print *,'global region statistics:',nseg,nvctr
     write(*,*)&
          'ERROR: problem in segkeys_periodic  ',&
          'nvctr_check:',nvctr_check,'nvctr_loc:',nvctr_loc,&
          'nend:',nend,'nsrt:',nsrt,'nseg_loc:',nseg_loc
     stop
  end if

  ! Now build the keyvloc where we replace the segments in order for the loc
  do iseg = 1, nseg_loc
     !sorting the keyg_loc
     loc = minloc(keyg_loc(1,:),1)
     keygloc(1,iseg) = keyg_loc(1,loc)
     keygloc(2,iseg) = keyg_loc(2,loc)
     !print *,'iseg,keygloc,keyg_loc',iseg,keygloc(1,loc),keygloc(2,loc),keyg_loc(1,iseg),keyg_loc(2,iseg)
     keyg_loc(1,loc) = maxval(keyg_loc) + 1
     keyvloc(iseg) = keyvglob(loc)
     !    print *,'iseg,keyglob,keyvglob,keygloc,keyvloc',iseg,keyglob(1,iseg),keyvglob(iseg),keygloc(1,iseg),keyvloc(iseg)
  end do

  call f_free(keyg_loc)

  call f_release_routine()

END SUBROUTINE segkeys_periodic

!> Divides the locreg into zones contained inside the simulation box, by applying the primitive vectors
!! It returns: astart(3,nzones) which is the starting points of the different zones (max. 8)
!!             aend(3,nzones) which is the ending points of the different zones (max. 8)
subroutine fracture_periodic_zone(nzones,Glr,Llr,outofzone,astart,aend)

  use module_base
  use locregs

  implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in) :: nzones
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Localization grid descriptors 

  !Subroutine Array Arguments
  integer,dimension(3),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  integer,dimension(3,nzones),intent(out) :: astart !
  integer,dimension(3,nzones),intent(out) :: aend !

  !local variables
  integer :: ii,index,jj
  integer,dimension(3) :: alrs,alre,Gend,Gstart,period

  ! Start and end of Global region
  Gstart(1) = Glr%ns1 
  Gstart(2) = Glr%ns2
  Gstart(3) = Glr%ns3  
  Gend(1) = Glr%ns1 + Glr%d%n1
  Gend(2) = Glr%ns2 + Glr%d%n2
  Gend(3) = Glr%ns3 + Glr%d%n3

  ! Periodicity of the system
  period(1) = Glr%d%n1 + 1
  period(2) = Glr%d%n2 + 1
  period(3) = Glr%d%n3 + 1

  ! Start and end of local region
  alrs(1) = Llr%ns1
  alrs(2) = Llr%ns2
  alrs(3) = Llr%ns3
  alre(1) = Llr%ns1 + Llr%d%n1
  alre(2) = Llr%ns2 + Llr%d%n2
  alre(3) = Llr%ns3 + Llr%d%n3

  !assign the first zone (necessarily without shift) and initialize the rest
  do ii=1,3
     astart(ii,:) = alrs(ii)
     aend(ii,:) = min(Gend(ii),alre(ii))
  end do

  !assign the other zones
  index = 2
  do ii=1,3
     if(outofzone(ii) > 0) then    !Translation: X,Y,Z
        astart(ii,index) = Gstart(ii)
        aend(ii,index) = modulo(alre(ii),period(ii))
        index = index + 1
     end if
     do jj=ii+1,3
        if(outofzone(ii) > 0 .and. outofzone(jj) > 0) then  !Translation: X+Y,X+Z,Y+Z
           astart(ii,index) = Gstart(ii)
           astart(jj,index) = Gstart(jj)
           aend(ii,index) = modulo(alre(ii),period(ii))
           aend(jj,index) = modulo(alre(jj),period(jj))
           index = index + 1
        end if
     end do
  end do

  if(outofzone(1) > 0 .and. outofzone(2) > 0 .and. outofzone(3) > 0 ) then ! Translation: X+Y+Z
     astart(1,index) = Gstart(1)
     astart(2,index) = Gstart(2)
     astart(3,index) = Gstart(3)
     aend(1,index) = modulo(alre(1),period(1))
     aend(2,index) = modulo(alre(2),period(2))
     aend(3,index) = modulo(alre(3),period(3))
  end if

END SUBROUTINE fracture_periodic_zone

!> routine moved as external to the module to avoid the compiler to create temporary arrays in the stack
subroutine transform_keyglob_to_keygloc(gmesh,Llr,nseg,keyglob,keygloc)
  use module_base
  use locregs, only: locreg_descriptors
  use box
  !use module_interfaces
  implicit none
  type(cell), intent(in) :: gmesh
  type(locreg_descriptors),intent(in) :: Llr
  integer, intent(in) :: nseg
  integer, dimension(2,nseg),intent(in) :: keyglob
  integer, dimension(2,nseg),intent(out) :: keygloc
  !local variables
  integer :: i, j, j0, ii, iz, iy, ix, n1p1, np, span

  call f_routine(id='transform_keyglob_to_keygloc')

  if (all(Llr%mesh%dom%bc == 0)) then
     n1p1=gmesh%ndims(1)
     np=n1p1*gmesh%ndims(2)
     do j = 1, nseg
        ! Writing keyglob in cartesian coordinates
        j0 = keyglob(1,j)
        span = keyglob(2, j) - keyglob(1, j)
        ii = j0-1
        iz = ii/np
        ii = ii-iz*np
        iy = ii/n1p1
        ix = ii-iy*n1p1

        ! Apply periodicity
        call withPer(ix, Llr%ns1, Llr%d%n1, gmesh%ndims(1), gmesh%dom%bc(1) == 1)
        call withPer(iy, Llr%ns2, Llr%d%n2, gmesh%ndims(2), gmesh%dom%bc(2) == 1)
        call withPer(iz, Llr%ns3, Llr%d%n3, gmesh%ndims(3), gmesh%dom%bc(3) == 1)
        ! Checking consistency
        if(iz < Llr%ns3 .or. iy < Llr%ns2 .or. ix < Llr%ns1) &
             & stop 'transform_keyglob_to_keygloc : minimum overflow'
        if(iz > Llr%ns3+Llr%d%n3 .or. iy > Llr%ns2+Llr%d%n2 .or. ix > Llr%ns1+Llr%d%n1 .or. ix + span > Llr%ns1+Llr%d%n1) &
             & stop 'transform_keyglob_to_keygloc : maximum overflow'

        keygloc(1,j) = (iz-Llr%ns3)*(Llr%d%n1+1)*(Llr%d%n2+1) + (iy-Llr%ns2)*(Llr%d%n1+1) + (ix-Llr%ns1) + 1
!!$        if (ix + span > Llr%ns1 + Llr%d%n1) span = Llr%ns1 + Llr%d%n1 - ix
        keygloc(2,j) = keygloc(1,j) + span
     end do
  else
     call f_memcpy(n = nseg * 2, src = keyglob, dest = keygloc(1,1))
  end if

  call f_release_routine()

contains
  subroutine withPer(i, is, in, N, per)
    implicit none
    integer, intent(inout) :: i
    integer, intent(in) :: is, in, N
    logical, intent(in) :: per

    if (.not. per) return
    do while (i > is + in)
       i = i - N
    end do
    do while (i < is)
       i = i + N
    end do
  end subroutine withPer
end subroutine transform_keyglob_to_keygloc


!from here the routines which are used in the cubic code
!> Calculates the length of the keys describing a wavefunction data structure
subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
  use dynamic_memory
  implicit none
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3
  logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid 
  integer, intent(out) :: mseg,mvctr
  !local variables
  logical :: plogrid
  integer :: i1,i2,i3,nsrt,nend,nsrti,nendi,mvctri

  call f_routine(id='num_segkeys')

  mvctr=0
  nsrt=0
  nend=0
  !$omp parallel default(private) shared(nl3,nu3,nl2,nu2,nl1,nu1,logrid,mvctr,nsrt,nend)
  mvctri=0
  nsrti=0
  nendi=0
  !$omp do  
  do i3=nl3,nu3 
     do i2=nl2,nu2
        plogrid=.false.
        do i1=nl1,nu1
           if (logrid(i1,i2,i3)) then
              mvctri=mvctri+1
              if (.not. plogrid) then
                 nsrti=nsrti+1
              endif
           else
              if (plogrid) then
                 nendi=nendi+1
              endif
           endif
           plogrid=logrid(i1,i2,i3)
        enddo
        if (plogrid) then
           nendi=nendi+1
        endif
     enddo
  enddo
  !$omp enddo
  !$omp critical
  mvctr=mvctr+mvctri
  nsrt=nsrt+nsrti
  nend=nend+nendi
  !$omp end critical
  !$omp end parallel
  if (nend /= nsrt) then 
     write(*,*)' ERROR: nend <> nsrt',nend,nsrt
     stop 
  endif
  mseg=nend

  call f_release_routine()

END SUBROUTINE num_segkeys


!> Calculates the keys describing a wavefunction data structure
subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
  use dynamic_memory
  implicit none
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,mseg
  logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid  
  integer, dimension(mseg), intent(out) :: keyv
  integer, dimension(2,mseg), intent(out) :: keyg
  !local variables
  logical :: plogrid
  integer :: mvctr,nsrt,nend,i1,i2,i3,ngridp,np,n1p1

  call f_routine(id='segkeys')

  mvctr=0
  nsrt=0
  nend=0
  n1p1=n1+1
  np=n1p1*(n2+1)
  do i3=nl3,nu3 
     do i2=nl2,nu2
        plogrid=.false.
        do i1=nl1,nu1
           ngridp=i3*np + i2*n1p1 + i1+1
           if (logrid(i1,i2,i3)) then
              mvctr=mvctr+1
              if (.not. plogrid) then
                 nsrt=nsrt+1
                 keyg(1,nsrt)=ngridp
                 keyv(nsrt)=mvctr
              endif
           else
              if (plogrid) then
                 nend=nend+1
                 keyg(2,nend)=ngridp-1
              endif
           endif
           plogrid=logrid(i1,i2,i3)
        enddo
        if (plogrid) then
           nend=nend+1
           keyg(2,nend)=ngridp
        endif
     enddo
  enddo
  if (nend /= nsrt) then 
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif
  !mseg=nend

  call f_release_routine()

END SUBROUTINE segkeys

!> Set up an array logrid(i1,i2,i3) that specifies whether the grid point
!! i1,i2,i3 is the center of a scaling function/wavelet
subroutine fill_logrid(mesh,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
     &   ntypes,iatype,rxyz,radii,rmult,logrid)
  use module_base
  use sparsematrix_init, only: distribute_on_tasks
  use box
  use at_domain, only: domain_geocode
  implicit none
  !Arguments
  type(cell), intent(in) :: mesh
  integer, intent(in) :: nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,ntypes
  real(gp), intent(in) :: rmult
  integer, dimension(nat), intent(in) :: iatype
  real(gp), dimension(ntypes), intent(in) :: radii
  real(gp), dimension(3,nat), intent(in) :: rxyz
  logical, dimension(0:mesh%ndims(1)-1,0:mesh%ndims(2)-1,0:mesh%ndims(3)-1), intent(out) :: logrid
  !local variables
  integer, parameter :: START_=1,END_=2
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: i1,i2,i3,iat,ml1,ml2,ml3,mu1,mu2,mu3,j1,j2,j3,i1s,i1e,i2s,i2e,i3s,i3e,i
  integer :: natp, isat, iiat
  !$ integer, external:: omp_get_num_threads,omp_get_thread_num
  !$ integer :: ithread,nthread
  real(gp) :: dx,dy2,dz2,rad,dy2pdz2,radsq
  logical :: parallel
  integer, dimension(2,3) :: nbox_limit,nbox,nbox_tmp
!  logical, dimension(0:n1,0:n2,0:n3) :: logrid_tmp
  type(box_iterator) :: bit
  !logical, dimension(0:n1,0:n2,0:n3) :: logrid_tmp

  call f_routine(id='fill_logrid')

  nbox_limit(START_,1)=nl1
  nbox_limit(END_,1)=nu1
  nbox_limit(START_,2)=nl2
  nbox_limit(END_,2)=nu2
  nbox_limit(START_,3)=nl3
  nbox_limit(END_,3)=nu3

  if (any(mesh%dom%bc /= 0) .and. nbuf /=0) then
        call f_err_throw('ERROR: a nonzero value of nbuf is allowed only for Free BC (tails)',&
             err_name='BIGDFT_RUNTIME_ERROR')
        return
  end if

  if (all(mesh%dom%bc == 0)) then
     !$omp parallel default(none) &
     !$omp shared(nl3, nu3, nl2, nu2, nl1, nu1, logrid) &
     !$omp private(i3, i2, i1)
     !$omp do schedule(static)
     do i3=nl3,nu3 
        do i2=nl2,nu2 
           do i1=nl1,nu1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
     !$omp end do
     !$omp end parallel
  else !
     !Special case if no atoms (homogeneous electron gas): all points are used (TD)
     if (nat == 0) then
        !$omp parallel default(shared) &
        !$omp private(i3, i2, i1)
        !$omp do schedule(static)
        do i3=0,mesh%ndims(3) - 1
           do i2=0,mesh%ndims(2) - 1
              do i1=0,mesh%ndims(1) - 1
                 logrid(i1,i2,i3)=.true.
              enddo
           enddo
        enddo
        !$omp end do
        !$omp end parallel
     else
        !$omp parallel default(shared) &
        !$omp private(i3, i2, i1)
        !$omp do schedule(static)
        do i3=0,mesh%ndims(3) - 1
           do i2=0,mesh%ndims(2) - 1
              do i1=0,mesh%ndims(1) - 1
                 logrid(i1,i2,i3)=.false.
                 !logrid_tmp(i1,i2,i3)=.false.
              enddo
           enddo
        enddo
        !$omp end do
        !$omp end parallel
     end if
  end if

  ! MPI parallelization over the atoms, ony if there are many atoms.
  ! Maybe 200 is too low, but in this way there is a test for this feature.
  if (nat>2000) then
     call distribute_on_tasks(nat, bigdft_mpi%iproc, bigdft_mpi%nproc, natp, isat)
     parallel = .true.
  else
     natp = nat
     isat = 0
     parallel = .false.
  end if

!  logrid_tmp=.false.

  do iat=1,natp
     iiat = iat + isat
     if (radii(iatype(iiat)) == 0.0_gp) cycle
     rad=radii(iatype(iiat))*rmult+(real(nbuf,gp)+eps_mach)*maxval(mesh%hgrids)
     nbox=box_nbox_from_cutoff(mesh,rxyz(:,iiat),rad)

     if (domain_geocode(mesh%dom) == 'F') then
        if (any( (nbox(START_,:) < nbox_limit(START_,:)) .or. (nbox(END_,:) > nbox_limit(END_,:)) )) then
           call f_err_throw('The box of the atom '//trim(yaml_toa(iat))//' is outside the limit; '//&
                trim(yaml_toa(reshape(nbox,[6])))//', '//trim(yaml_toa(reshape(nbox_limit,[6]))),&
                err_name='BIGDFT_RUNTIME_ERROR')
        end if
     end if

     !in the case of periodic bc there should be no need to wrap around multiple times here
     do i=1,3
        nbox(START_,i)=max(nbox(START_,i),-(mesh%ndims(i)-1)/2-1)
        nbox(END_,i)=min(nbox(END_,i),mesh%ndims(i)+(mesh%ndims(i)-1)/2)
     end do

!!$     ml1=ceiling((rxyz(1,iiat)-rad)/hx - eps_mach)  
!!$     ml2=ceiling((rxyz(2,iiat)-rad)/hy - eps_mach)   
!!$     ml3=ceiling((rxyz(3,iiat)-rad)/hz - eps_mach)   
!!$     mu1=floor((rxyz(1,iiat)+rad)/hx + eps_mach)
!!$     mu2=floor((rxyz(2,iiat)+rad)/hy + eps_mach)
!!$     mu3=floor((rxyz(3,iiat)+rad)/hz + eps_mach)
!!$     i3s=max(ml3,-n3/2-1)
!!$     i3e=min(mu3,n3+n3/2+1)
!!$     i2s=max(ml2,-n2/2-1)
!!$     i2e=min(mu2,n2+n2/2+1)
!!$     i1s=max(ml1,-n1/2-1)
!!$     i1e=min(mu1,n1+n1/2+1)
!!$
!!$
!!$
!!$     !print *,'limitold',ml3,mu3,i3s,i3e
!!$     !print *,'limitnew',nbox(:,3)
!!$

     bit=box_iter(mesh,nbox=nbox+1) !add here a plus one for the convention of ndims

     !split the iterator for openmp parallelisation
     !$omp parallel firstprivate(bit) private(ithread)
     !$ nthread=omp_get_num_threads()
     !$ ithread=omp_get_thread_num()
     !$ call box_iter_split(bit,nthread,ithread)
     call fill_logrid_with_spheres(bit,rxyz(1,iiat),rad,logrid)
     !$ call box_iter_merge(bit)
     !$omp end parallel  

!!$     call fill_logrid_with_spheres(bit,rxyz(1,iiat),rad,logrid_tmp)
!!$
!!$        ml1=ceiling((rxyz(1,iiat)-rad)/hx - eps_mach)  
!!$        ml2=ceiling((rxyz(2,iiat)-rad)/hy - eps_mach)   
!!$        ml3=ceiling((rxyz(3,iiat)-rad)/hz - eps_mach)   
!!$        mu1=floor((rxyz(1,iiat)+rad)/hx + eps_mach)
!!$        mu2=floor((rxyz(2,iiat)+rad)/hy + eps_mach)
!!$        mu3=floor((rxyz(3,iiat)+rad)/hz + eps_mach)
!!$
!!$        !for Free BC, there must be no incoherences with the previously calculated delimiters
!!$        if (geocode(1:1) == 'F') then
!!$           if (ml1 < nl1) then
!!$              write(*,'(a,i0,3x,i0)')  'ERROR: ml1 < nl1  ', ml1, nl1
!!$              stop
!!$           end if
!!$           if (ml2 < nl2) then
!!$              write(*,'(a,i0,3x,i0)')  'ERROR: ml2 < nl2  ', ml2, nl2
!!$              stop
!!$           end if
!!$           if (ml3 < nl3) then
!!$              write(*,'(a,i0,3x,i0)')  'ERROR: ml3 < nl3  ', ml3, nl3
!!$              stop
!!$           end if
!!$
!!$           if (mu1 > nu1) then
!!$              write(*,'(a,i0,3x,i0)')  'ERROR: mu1 > nu1  ', mu1, nu1
!!$              stop
!!$           end if
!!$           if (mu2 > nu2) then
!!$              write(*,'(a,i0,3x,i0)')  'ERROR: mu2 > nu2  ', mu2, nu2
!!$              stop
!!$           end if
!!$           if (mu3 > nu3) then
!!$              write(*,'(a,i0,3x,i0)')  'ERROR: mu3 > nu3  ', mu3, nu3
!!$              stop
!!$           end if
!!$        end if
!!$
!!$        i3s=max(ml3,-n3/2-1)
!!$        i3e=min(mu3,n3+n3/2+1)
!!$        i2s=max(ml2,-n2/2-1)
!!$        i2e=min(mu2,n2+n2/2+1)
!!$        i1s=max(ml1,-n1/2-1)
!!$        i1e=min(mu1,n1+n1/2+1)
!!$        radsq=rad**2
!!$
!!$        nbox_tmp(START_,3)=i3s
!!$        nbox_tmp(END_,3)=i3e
!!$        nbox_tmp(START_,2)=i2s
!!$        nbox_tmp(END_,2)=i2e
!!$        nbox_tmp(START_,1)=i1s
!!$        nbox_tmp(END_,1)=i1e
!!$
!!$
!!$        call f_assert(all(nbox_tmp == nbox),id='box different')
!!$
!!$        !what follows works always provided the check before
!!$        !$omp parallel default(shared) private(i3,dz2,j3,i2,dy2,j2,i1,j1,dx,dy2pdz2)
!!$        !$omp do schedule(static,1)
!!$        do i3=i3s,i3e
!!$           dz2=(real(i3,gp)*hz-rxyz(3,iiat))**2-eps_mach
!!$           if (dz2>radsq) cycle
!!$           j3=modulo(i3,n3+1)
!!$           do i2=i2s,i2e
!!$              dy2=(real(i2,gp)*hy-rxyz(2,iiat))**2
!!$              dy2pdz2=dy2+dz2
!!$              if (dy2pdz2>radsq) cycle
!!$              j2=modulo(i2,n2+1)
!!$              do i1=i1s,i1e
!!$                 j1=modulo(i1,n1+1)
!!$                 dx=real(i1,gp)*hx-rxyz(1,iiat)
!!$                 if (dx**2+dy2pdz2 <= radsq) then 
!!$                    logrid(j1,j2,j3)=.true.
!!$                 endif
!!$!!!                 if ((logrid(j1,j2,j3) .and. .not. logrid_tmp(j1,j2,j3)) .and. j3==0) then
!!$!!!                    print *,'j1,j2,j3',j1,j2,j3,radsq,nbox
!!$!!!                    print *,logrid(j1,j2,j3),logrid_tmp(j1,j2,j3)
!!$!!!                    print *,'BB',i1,i2,i3,dx**2+dy2pdz2,radsq
!!$!!!                    stop
!!$!!!                 end if
!!$              enddo
!!$           enddo
!!$        enddo
!!$        !$omp enddo
!!$        !$omp end parallel
  enddo

  if (parallel) then
     call fmpi_allreduce(logrid,FMPI_LOR, comm=bigdft_mpi%mpi_comm)
  end if

!!$  print *,'count',count(logrid .neqv. logrid_tmp)
!!$
!!$
!!$  call f_assert(all(logrid .eqv. logrid_tmp),'logrid')
 
  call f_release_routine()

END SUBROUTINE fill_logrid

!> Tick .true. inside the sphere of radius rad and center rxyz0.
!! Used to build up the localization regions around each atoms.
subroutine fill_logrid_with_spheres(bit,rxyz0,rad,logrid)
  use module_defs, only: gp
  use box
  use at_domain, only: square_gd
  implicit none
  type(box_iterator) :: bit
  real(gp), dimension(3), intent(in) :: rxyz0
  real(gp), intent(in) :: rad
  logical, dimension(bit%mesh%ndims(1),bit%mesh%ndims(2),bit%mesh%ndims(3)), intent(inout) :: logrid
  !local variables
  do while(box_next_point(bit))
     ! Tick .true. inside the sphere of radius rad and center rxyz0
     bit%tmp=bit%mesh%hgrids*(bit%inext-2)-rxyz0-bit%oxyz
!!$     if (bit%k==1 .or. bit%k==24) then
!!$        print *,'AA',bit%tmp,square_gd(bit%mesh,bit%tmp),rad**2
!!$        print *,'ii',bit%inext-2
!!$        print *,bit%i,bit%j,bit%k
!!$        print *,bit%nbox
!!$     end if
     if (square_gd(bit%mesh%dom,bit%tmp) <= rad**2) then
        logrid(bit%i,bit%j,bit%k)=.true.
     end if
  end do

end subroutine fill_logrid_with_spheres
