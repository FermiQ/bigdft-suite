!> @file
!!    Modulefile for handling fundamental data structed and methods of the simulation box
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2016-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module box

  use f_precisions, gp=>f_double
  use numerics, only: onehalf,pi
  use at_domain

  private

  !>parameter for the definition of the bc
!!$  integer, parameter :: NULL_BC=-100
  integer, parameter :: FREE=0
  integer, parameter :: PERIODIC=1

  !to ease readiness
  integer, parameter :: START_=1,END_=2
  integer, parameter :: X_=1,Y_=2,Z_=3

  !> data type which stores all informations of the simulation box. 
  !! It contains also the metric for nonorthorhombic cells.
  type, public :: cell
     type(domain) :: dom !< data type for the simulation domain
     !logical :: orthorhombic !<true if the cell is orthorhombic
     !integer, dimension(3) :: bc !< boundary conditions on each direction (FREE=0, PERIODIC=1)
     integer, dimension(3) :: ndims !< number of grid points on each direction
     real(gp), dimension(3) :: hgrids !< real space grid on each direction
     !real(gp), dimension(3) :: angrad !<angles between the dimensions in radiant (alpha_bc,beta_ac,gamma_bc)
     !derived data
     integer(f_long) :: ndim !< product of the dimension, long integer to avoid overflow
     real(gp) :: volume_element !< volume element of the primitive cell
     real(gp), dimension(3,3) :: habc !<primitive volume elements in the translation vectors direction
     !real(gp), dimension(3,3) :: uabc !<matrix of the normalized translation vectors direction
     !real(gp), dimension(3,3) :: gd !<covariant metric needed for non-orthorhombic operations
     !real(gp), dimension(3,3) :: gu !<controvariant metric needed for non-orthorhombic operations
     !real(gp) :: detgd !<determinant of the covariant matrix
  end type cell

  !> defines the object to iterate around the real-space grid points.
  !! given a cell type, it might iterate on a section of this gris provided by the extremes nbox
  !! it also provides a facility to parallelize over the
  type, public :: box_iterator
     integer :: i3s=-1 !<starting point in the dimension z
     integer :: i3e=-1 !<ending point in the dimension z
     integer :: i23=-1 !<collapsed index in 23 dimension (in relative conventions)
     integer :: ind=-1 !<one-dimensional index for arrays (in relative conventions)
     !> indices in absolute coordinates in the given box,
     !! from nbox(1,:) to nbox(2,:). To be intended as private
     integer, dimension(3)  :: inext=0
     !> actual index inside the box,from 1 to mesh%ndims(:)
     integer :: i,j,k !better as scalars
     !> Sub-box to iterate over the points (ex. around atoms)
     !! start and end points for each direction
     integer, dimension(2,3) :: nbox=-1
     real(gp), dimension(3) :: oxyz=-1.0_gp !<origin of the coordinate system
     real(gp), dimension(3) :: rxyz=-1.0_gp !<coordinates of the grid point which is internal at the box
     real(gp), dimension(3) :: rxyz_nbox=-1.0_gp !<coordinates of the grid point which is internal at the box
     real(gp), dimension(3) :: tmp=0.0_gp !< size 3 array buffer to avoid the creation of temporary arrays
     logical :: whole=.false. !<to assess if we run over the entire box or not (no check over the internal point)
     integer, dimension(2,3) :: subbox=-1 !<box of the local task
     !>reference mesh from which it starts
     type(cell), pointer :: mesh=>null()
  end type box_iterator

!!$  interface box_iter
!!$     module procedure box_iter_c,box_iter_base
!!$  end interface box_iter

  public :: cell_r,cell_new,box_iter,box_next_point
  public :: box_next_x,box_next_y,box_next_z,cell_null,nullify_box_iterator
  public :: box_iter_rewind,box_iter_split,box_iter_merge,box_iter_set_nbox,box_iter_expand_nbox,box_nbox_from_cutoff
  public :: box_iter_square_gd,box_iter_closest_r,box_iter_distance

contains

  !> Nullify the cell type
  pure function cell_null() result(me)
   implicit none
   type(cell) :: me
   me%dom=domain_null()
!!$   me%orthorhombic=.true.
!!$   me%bc=0
   me%ndims=0
   me%hgrids=0.0_gp
!!$   me%angrad=0.0_gp
   !derived data
   me%ndim=0
   me%volume_element=0.0_gp
   me%habc=0.0_gp
!!$   me%uabc=0.0_gp
!!$   me%gd=0.0_gp
!!$   me%gu=0.0_gp
!!$   me%detgd=0.0_gp
  end function cell_null

  !> Nullify the iterator dpbox type
  pure subroutine nullify_box_iterator(boxit)
    implicit none
    type(box_iterator), intent(out) :: boxit
    boxit%inext(X_)=-1
  end subroutine nullify_box_iterator

!!$  function box_iter_c(mesh,origin) result(boxit)
!!$    type(cell), intent(in), target :: mesh
!!$    !> starting point of the box in the z direction
!!$    integer, intent(in), optional :: i3s
!!$    !> number of planes of the box to be considered
!!$    integer, intent(in), optional :: n3p
!!$    !> Box of start and end points which have to be considered
!!$    integer, dimension(2,3), intent(in), optional :: nbox
!!$    !> real coordinates of the origin in the reference frame of the
!!$    !box (the first point has the 000 coordinate)
!!$    real(gp), dimension(3), intent(in), optional :: origin
!!$    type(box_iterator) :: boxit
!!$
!!$  end function box_iter_c

  !>define an iterator over the cell points
  function box_iter(mesh,nbox,origin,i3s,n3p,centered,cutoff) result(boxit)
    use f_utils, only: f_zero
    implicit none
    type(cell), intent(in), target :: mesh
    !>when true the origin is placed at the center of the box, origin is ignored
    logical, intent(in), optional :: centered
    !> starting point of the box in the z direction
    integer, intent(in), optional :: i3s
    !> number of planes of the box to be considered
    integer, intent(in), optional :: n3p
    !> Box of start and end points which have to be considered
    integer, dimension(2,3), intent(in), optional :: nbox
    real(gp), intent(in), optional :: cutoff !< determine the box around the origin
    !> real coordinates of the origin in the reference frame of the
    !! box (the first point has the 000 coordinate)
    real(gp), dimension(3), intent(in), optional :: origin

    type(box_iterator) :: boxit

    call nullify_box_iterator(boxit)

    !if the mesh is invalid (e.g. no dims, return)
    if (mesh%ndim==0) return
    !associate the mesh
    boxit%mesh => mesh

    call f_zero(boxit%oxyz)
    if (present(origin)) boxit%oxyz=origin
    if (present(centered)) then
       if (centered) boxit%oxyz=0.5_gp*real(boxit%mesh%ndims)*boxit%mesh%hgrids
    end if

    if (present(i3s)) then
       boxit%i3s=i3s
    else
       boxit%i3s=1
    end if

    if (present(n3p)) then
       boxit%i3e=boxit%i3s+n3p-1
    else
       boxit%i3e=boxit%i3s+mesh%ndims(Z_)-1
    end if

    call box_iter_set_nbox(boxit,nbox,boxit%oxyz,cutoff)

    call probe_iterator(boxit)

  end function box_iter

  pure subroutine box_iter_set_nbox(bit,nbox,oxyz,cutoff)
    implicit none
    type(box_iterator), intent(inout) :: bit
    real(gp), dimension(3), optional, intent(in) :: oxyz
    real(gp), intent(in), optional :: cutoff
    integer, dimension(2,3), intent(in), optional :: nbox

    if(present(nbox)) then
       bit%nbox=nbox
       bit%whole=.false.
       call set_subbox(bit%mesh%dom%bc,bit%mesh%ndims,bit%nbox,bit%subbox)
       call box_iter_rewind(bit)
    else if (present(cutoff)) then
!!$
!!$       bit%nbox(START_,:)=floor((oxyz-cutoff)/bit%mesh%hgrids)
!!$       bit%nbox(END_,:)=ceiling((oxyz+cutoff)/bit%mesh%hgrids)
       bit%nbox=box_nbox_from_cutoff(bit%mesh,oxyz,cutoff)
       bit%whole=.false.
       call set_subbox(bit%mesh%dom%bc,bit%mesh%ndims,bit%nbox,bit%subbox)
       call box_iter_rewind(bit)
    else
       call box_iter_expand_nbox(bit)
    end if

    call box_iter_rewind(bit)

  end subroutine box_iter_set_nbox

  !> this function has to be genralized for non-orthorhombic grids
  pure function box_nbox_from_cutoff(mesh,oxyz,cutoff,inner) result(nbox)
    implicit none
    type(cell), intent(in) :: mesh
    real(gp), dimension(3), intent(in) :: oxyz
    real(gp), intent(in) :: cutoff
    integer, dimension(2,3) :: nbox
    logical, intent(in), optional :: inner
    !local variables
    logical :: inner_
    real(gp), dimension(2,3) :: rbox
    !for non-orthorhombic cells the concept of distance has to be inserted here (the box should contain the sphere)
!!$    nbox(START_,:)=floor((oxyz-cutoff)/mesh%hgrids)
!!$    nbox(END_,:)=ceiling((oxyz+cutoff)/mesh%hgrids)

    inner_=.true.
    if (present(inner)) inner_=inner

    rbox=cell_cutoff_extrema(mesh,oxyz,cutoff)

    if (inner_) then
       nbox(START_,:)=ceiling(rbox(START_,:)/mesh%hgrids)
       nbox(END_,:)=floor(rbox(END_,:)/mesh%hgrids)
    else
       nbox(START_,:)=floor(rbox(START_,:)/mesh%hgrids)
       nbox(END_,:)=ceiling(rbox(END_,:)/mesh%hgrids)
    end if

  end function box_nbox_from_cutoff

  pure function cell_cutoff_extrema(mesh,oxyz,cutoff) result(rbox)
    !use yaml_strings
    implicit none
    type(cell), intent(in) :: mesh
    real(gp), dimension(3), intent(in) :: oxyz
    real(gp), intent(in) :: cutoff
    real(gp), dimension(2,3) :: rbox
    !local variables
    real(gp), dimension(3) :: fac
    real(gp) :: aa,b,c,p,area,h,l,p2,area2,hf
    integer :: i,i1,i2
    !for non-orthorhombic cells the concept of distance has to be inserted here (the box should contain the sphere)
    ! compute the inverse of mesh%uabc

    if (mesh%dom%orthorhombic) then
        rbox(START_,:)=oxyz-cutoff
        rbox(END_,:)=oxyz+cutoff
    else
        rbox(START_,:)=1.0d10
        rbox(END_,:) =-1.0d10
        fac(:)=1.0_gp
        do i=1,3
         i1=mod(i,3)+1
         i2=mod(i+1,3)+1
         if (abs(mesh%dom%angrad(i1) - onehalf*pi) .lt. 1.0d-15) then
          if (abs(mesh%dom%angrad(i2) - onehalf*pi) .lt. 1.0d-15) then
           hf=1.0_gp
          else
           c=cos(mesh%dom%angrad(i2))/sin(mesh%dom%angrad(i))
           hf=sqrt(1.0_gp-c**2)
          end if
         else if (abs(mesh%dom%angrad(i2) - onehalf*pi) .lt. 1.0d-15) then
          if (abs(mesh%dom%angrad(i1) - onehalf*pi) .lt. 1.0d-15) then
           hf=1.0_gp
          else
           c=cos(mesh%dom%angrad(i1))/sin(mesh%dom%angrad(i))
           hf=sqrt(1.0_gp-c**2)
          end if
         else
          aa=1.0_gp/cos(mesh%dom%angrad(i2))
          b=1.0_gp/cos(mesh%dom%angrad(i1))
          c=sqrt(aa**2+b**2-2.0_gp*aa*b*cos(mesh%dom%angrad(i)))
          p=(tan(mesh%dom%angrad(i2))+tan(mesh%dom%angrad(i1))+c)*0.5_gp
          area=sqrt(p*(p-tan(mesh%dom%angrad(i2)))*(p-tan(mesh%dom%angrad(i1)))*(p-c))
          h=2.0_gp*area/c
          l=sqrt(1.0_gp+h**2)
          p2=(1.0_gp+h+l)*0.5_gp
          area2=sqrt(p2*(p2-1.0_gp)*(p2-h)*(p2-l))
          hf=2.0_gp*area2/l
         end if
         fac(i)=1.0_gp/hf
        end do
        do i=1,3
         i1=mod(i,3)+1
         i2=mod(i+1,3)+1
         rbox(START_,i1)=min(rbox(START_,i1),oxyz(i1)-cutoff*fac(i1))
         rbox(END_,i1)  =max(rbox(END_,i1),oxyz(i1)+cutoff*fac(i1))
         rbox(START_,i2)=min(rbox(START_,i2),oxyz(i2)-cutoff*fac(i2))
         rbox(END_,i2)  =max(rbox(END_,i2),oxyz(i2)+cutoff*fac(i2))
        end do
        do i=1,3
         if (rbox(START_,i) .lt. 0.0_gp) rbox(START_,i) = 0.0_gp
         if (rbox(END_,i) .gt. mesh%hgrids(i)*(mesh%ndims(i)-1)) rbox(END_,i) = mesh%hgrids(i)*(mesh%ndims(i)-1)
        end do
    end if

  end function cell_cutoff_extrema

  pure subroutine box_iter_expand_nbox(bit)
    implicit none
    type(box_iterator), intent(inout) :: bit
    bit%whole=.true.
    bit%nbox(START_,:)=1
    bit%nbox(END_,:)=bit%mesh%ndims
    call set_subbox(bit%mesh%dom%bc,bit%mesh%ndims,bit%nbox,bit%subbox)
    call box_iter_rewind(bit)
  end subroutine box_iter_expand_nbox

  pure subroutine set_subbox(bc,ndims,nbox,subbox)
    implicit none
    integer, dimension(3), intent(in) :: bc,ndims
    integer, dimension(2,3), intent(in) :: nbox
    integer, dimension(2,3), intent(out) :: subbox
    !local variables
    integer :: i

    do i=1,3
       if (bc(i)==PERIODIC) then
          subbox(:,i)=nbox(:,i)
       else
          subbox(START_,i)=max(1,nbox(START_,i))
          subbox(END_,i)=min(ndims(i),nbox(END_,i))
       end if
    end do
  end subroutine set_subbox

  !>verify if the iterator can be used as expected
  subroutine probe_iterator(bit)
    use f_precisions
    use yaml_strings
    use dictionaries
    use dynamic_memory
    use f_arrays
    use f_utils, only: f_assert
    implicit none
    type(box_iterator), intent(inout) :: bit
    !local variables
    logical :: impossible_wrap_on_z
    integer :: iz,iy,ix,i,jx,jy,jz
    integer(f_long) :: icnt,itgt,itgt_inner
    logical(f_byte), dimension(:), allocatable :: lxyz
    integer, dimension(:), allocatable :: lx,ly,lz
    integer, dimension(3) :: subdims

    do i=1,3
       subdims(i)=bit%subbox(END_,i)-bit%subbox(START_,i)+1
    end do
    impossible_wrap_on_z=bit%mesh%dom%bc(Z_) /= PERIODIC.or. bit%subbox(END_,Z_)-bit%subbox(START_,Z_) < bit%mesh%ndims(Z_)
    if (impossible_wrap_on_z) subdims(3)=min(subdims(3),bit%i3e-bit%i3s+1)

    !first, count if the iterator covers all the required points
    itgt=product(int(subdims,f_long))
    itgt_inner=int(subdims(1),f_long)*int(subdims(2),f_long)*min(subdims(3),bit%i3e-bit%i3s+1)!int(bit%i3e-bit%i3s+1,f_long)

    !allocate array of values corresponding to the expected grid
!!$    lx=f_malloc(bit%mesh%ndims(1),id='lx')
!!$    ly=f_malloc(bit%mesh%ndims(2),id='ly')
!!$    lz=f_malloc(bit%i3e-bit%i3s+1,id='lz')
    lx=f_malloc(subdims(X_),id='lx')
    ly=f_malloc(subdims(Y_),id='ly')
    lz=f_malloc(subdims(Z_),id='lz')

    lxyz=f_malloc0(itgt,id='lxyz')

!!$    do iz=bit%i3s,bit%i3e
!!$       lz(iz-bit%i3s+1)=iz
!!$    end do
!!$    do iy=1,bit%mesh%ndims(2)
!!$       ly(iy)=iy
!!$    end do
!!$    do ix=1,bit%mesh%ndims(1)
!!$       lx(ix)=ix
!!$    end do

    do iz=1,subdims(Z_)
       jz=iz+bit%subbox(START_,Z_)-1
       if (bit%mesh%dom%bc(Z_) == PERIODIC) jz=modulo(jz-1,bit%mesh%ndims(Z_))+1
       lz(iz)=jz
    end do
    do iy=1,subdims(Y_)
       jy=iy+bit%subbox(START_,Y_)-1
       if (bit%mesh%dom%bc(Y_) == PERIODIC) jy=modulo(jy-1,bit%mesh%ndims(Y_))+1
       ly(iy)=jy
    end do
    do ix=1,subdims(X_)
       jx=ix+bit%subbox(START_,X_)-1
       if (bit%mesh%dom%bc(X_) == PERIODIC) jx=modulo(jx-1,bit%mesh%ndims(X_))+1
       lx(ix)=jx
    end do

    !separable mode
    iz=bit%subbox(START_,Z_)-1 !0
    icnt=0
    jz=0
    do while(box_next_z(bit))
       iz=iz+1
       jz=iz-bit%subbox(START_,Z_)+1
       !print *,'bit',bit%k,bit%inext(Z_)
       call f_assert(iz+bit%i3s-1==bit%inext(Z_)-1,'A')!,&
       !'Error iz='+iz+', inext(Z)='+bit%inext(Z_))
       iy=bit%subbox(START_,Y_)-1!0
       jy=0
       do while(box_next_y(bit))
          iy=iy+1
          jy=iy-bit%subbox(START_,Y_)+1
          call f_assert(iy==bit%inext(Y_)-1,'B')!,&
          !'Error iy='+iy+', inext(Y)='+bit%inext(Y_))
          ix=bit%subbox(START_,X_)-1!0
          jx=0
          do while(box_next_x(bit))
             ix=ix+1
             jx=ix-bit%subbox(START_,X_)+1
             call f_assert(ix==bit%inext(X_)-1,'C')!,&
             !'Error ix='+ix+', inext(X)='+bit%inext(X_))

             icnt=icnt+1
             call f_assert(lx(jx) == bit%i,'D')!,&
                 !'Error value, ix='+bit%i+', expected='+lx(jx))
             !convert the value of the logical array
             !if (lxyz(bit%ind)) &
             if (lxyz(icnt)) &
                  call f_err_throw('Error point ind='+bit%ind+&
               ', i,j,k='+yaml_toa([bit%i,bit%j,bit%k]))
             !lxyz(bit%ind)=f_T
             lxyz(icnt)=f_T
          end do
          call f_assert(jx == subdims(X_),'E')!,&
          !'Error boxit, ix='+ix+', itgtx='+subdims(X_))
          call f_assert(ly(jy) == bit%j,'F')!,&
          !'Error value, iy='+bit%j+', expected='+ly(iy))
       end do
       call f_assert(jy == subdims(Y_),'G')!,&
       !'Error boxit, iy='+iy+', itgty='+subdims(Y_))

       call f_assert(lz(jz)+bit%i3s-1 == bit%k,'H') !&
       !yaml_toa([lz(jz),bit%k,bit%i3s]))!,&
    end do
    if (impossible_wrap_on_z) then 
         call f_assert(jz == min(subdims(Z_),bit%i3e-bit%i3s+1),'I') !,&
    !'Error boxit, iz='+iz+', itgtz='+subdims(Z_))
         call f_assert(icnt == itgt_inner,'J')!,&
    !'Error sep boxit, icnt='+icnt+', itgt='+itgt)
      end if
    !complete mode
    if (all( bit%subbox(END_,:)-bit%subbox(START_,:) < bit%mesh%ndims) .or. all(bit%mesh%dom%bc == FREE)) then
       icnt=int(0,f_long)
       do while(box_next_point(bit))
          icnt=icnt+1
          !here we might see if there are points from which
          !we passed twice
          !print *,bit%i,bit%j,bit%k
          !if (.not. lxyz(bit%ind)) &
          if (.not. lxyz(icnt)) &
               call f_err_throw('Error point (2) ind='+bit%ind+&
               ', i,j,k='+yaml_toa([bit%i,bit%j,bit%k]))
          !lxyz(bit%ind)=f_F
          lxyz(icnt)=f_F
       end do
       call f_assert(icnt == itgt,'Error boxit, icnt='+icnt+&
            ', itgt='+itgt)
       !no points have to be left behind
       if (any(lxyz)) call f_err_throw('Error boxit, points not covered')
    end if



    call f_free(lxyz)
    call f_free(lx,ly,lz)

  end subroutine probe_iterator

  pure subroutine box_iter_rewind(bit)
    implicit none
    type(box_iterator), intent(inout) :: bit
    !local variables

    bit%inext=bit%subbox(START_,:)
!!$    if (bit%whole) then
!!$       bit%i=1
!!$       bit%j=1
!!$       bit%k=1
!!$    else
!!$       bit%i=-1
!!$       bit%j=-1
!!$       bit%k=-1
!!$    end if
    bit%k=bit%subbox(START_,Z_)-1
    bit%ind=0
    bit%i23=0

    if (bit%whole) bit%whole=bit%i3s == 1 .and. bit%i3e==bit%mesh%ndims(3)
  end subroutine box_iter_rewind

  !find the first z value which is available from the starting point
  function box_next_z(bit) result(ok)
    implicit none
    type(box_iterator), intent(inout) :: bit
    logical :: ok

!!    print *,'here',bit%inext(Z_),bit%k,bit%i3s,bit%i3e
!!$    call increment_dim(bit,3,bit%k,ok)
    ok = bit%i3e >= bit%i3s ! to be removed
    ok = bit%subbox(END_,Z_) >= bit%subbox(START_,Z_)
    if (.not. ok) return !there is nothing to explore
    ok= bit%inext(Z_) <= bit%subbox(END_,Z_)
!    print *,'ok1',ok,bit%inext(Z_),bit%subbox(END_,Z_),bit%k
    do while(ok)
       if (bit%whole) then
          bit%k=bit%inext(Z_)
       else
          call internal_point(bit%mesh%dom%bc(Z_),bit%inext(Z_),bit%mesh%ndims(Z_),&
               bit%k,bit%i3s,bit%i3e,ok)
!          print *,'ok1.2',ok,bit%inext(Z_),bit%k
          if (.not. ok) bit%inext(Z_)=bit%inext(Z_)+1
       end if
!       print *,'ok1.3',ok,bit%inext(Z_),bit%k
       if (ok) then
          bit%inext(Z_)=bit%inext(Z_)+1
          exit
       end if
       ok = bit%inext(Z_) <= bit%subbox(END_,Z_)
!       print *,'ok1.5',ok,bit%inext(Z_),bit%k
    end do
!    print *,'ok2',ok
    !reset x and y
    if (ok) then
       call update_boxit_z(bit)
       bit%inext(Y_)=bit%subbox(START_,Y_)
       bit%inext(X_)=bit%subbox(START_,X_)
    end if

    !in the case the z_direction is over, make the iterator ready for new use
!    print *,'o3',ok

    if (.not. ok) call box_iter_rewind(bit)

  end function box_next_z

  !find the first y value which is available from the starting point
  function box_next_y(bit) result(ok)
    implicit none
    type(box_iterator), intent(inout) :: bit
    logical :: ok

    call increment_dim(bit,2,bit%j,ok)
    !reset x
    if (ok) then
       call update_boxit_y(bit)
       bit%inext(X_)=bit%subbox(START_,X_)
    end if

  end function box_next_y

  !find the first x value which is available from the starting point
  function box_next_x(bit) result(ok)
    implicit none
    type(box_iterator), intent(inout) :: bit
    logical :: ok

    call increment_dim(bit,1,bit%i,ok)
    if (ok) call update_boxit_x(bit)
  end function box_next_x

  pure subroutine increment_dim(bit,idim,indi,ok)
    implicit none
    integer, intent(in) :: idim
    integer, intent(inout) :: indi
    type(box_iterator), intent(inout) :: bit
    logical, intent(out) :: ok

    ok= bit%inext(idim) <= bit%subbox(END_,idim)
    if (.not. ok) return
    if (bit%whole) then
       indi=bit%inext(idim)
    else
       if (bit%mesh%dom%bc(idim) == PERIODIC) then
          indi=modulo(bit%inext(idim)-1,bit%mesh%ndims(idim))+1
       else
          indi=bit%inext(idim)
       end if
    end if
    bit%inext(idim)=bit%inext(idim)+1

!!$    do !while(ok)
!!$       if (bit%whole) then
!!$          indi=bit%inext(idim)
!!$       else
!!$          if (bit%mesh%bc(idim) == PERIODIC) then
!!$             indi=modulo(bit%inext(idim)-1,bit%mesh%ndims(idim))+1
!!$          else
!!$             indi=bit%inext(idim)
!!$             !ok=indi >= 1 .and. indi <= bit%mesh%ndims(idim)
!!$             !if (ok) ok= indi <= bit%mesh%ndims(idim)
!!$          end if
!!$
!!$          !if (.not. ok) bit%inext(idim)=bit%inext(idim)+1
!!$       end if
!!$       if (ok) then
!!$          bit%inext(idim)=bit%inext(idim)+1
!!$          exit
!!$       end if
!!$       ok = bit%inext(idim) <= bit%subbox(END_,idim)
!!$       if (.not. ok) exit
!!$    end do
  end subroutine increment_dim

  pure subroutine update_boxit_x(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit

    !one dimensional index (to be corrected)
    boxit%ind = boxit%i+boxit%mesh%ndims(1)*boxit%i23

    !the position associated to the coordinates
    boxit%rxyz(X_)=cell_r(boxit%mesh,boxit%i,X_)-boxit%oxyz(X_)
    !and the position associated to the subbox
    boxit%rxyz_nbox(X_)=boxit%mesh%hgrids(X_)*&
         (boxit%inext(X_)-2)-boxit%oxyz(X_)

  end subroutine update_boxit_x

  pure subroutine update_boxit_y(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit
    !here we have the indices      boxit%inext as well as boxit%ixyz
    !we might then calculate the related quantities
    !two dimensional index, last two elements
    !to be corrected
    boxit%i23=(boxit%j-1)+&
         boxit%mesh%ndims(2)*(boxit%k-boxit%i3s)
    !the position associated to the coordinates
    boxit%rxyz(Y_)=cell_r(boxit%mesh,boxit%j,Y_)-boxit%oxyz(Y_)
    !and the position associated to the subbox
    boxit%rxyz_nbox(Y_)=boxit%mesh%hgrids(Y_)*&
         (boxit%inext(Y_)-2)-boxit%oxyz(Y_)

  end subroutine update_boxit_y

  pure subroutine update_boxit_z(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit

    !the position associated to the coordinates
    boxit%rxyz(Z_)=cell_r(boxit%mesh,boxit%k,Z_)-boxit%oxyz(Z_)

    !and the position associated to the subbox
    boxit%rxyz_nbox(Z_)=boxit%mesh%hgrids(Z_)*&
         (boxit%inext(Z_)-2)-boxit%oxyz(Z_)

  end subroutine update_boxit_z

!!!>  !this routine should not use inext as it is now prepared for the next step
!!!>  pure subroutine update_boxit(boxit)
!!!>    implicit none
!!!>    type(box_iterator), intent(inout) :: boxit
!!!>
!!!>    call update_boxit_x(boxit)
!!!>    call update_boxit_y(boxit)
!!!>    call update_boxit_z(boxit)
!!!>!!$
!!!>!!$    !one dimensional index (to be corrected)
!!!>!!$    boxit%ind = boxit%i+boxit%mesh%ndims(1)*boxit%i23
!!!>!!$    !here we have the indices      boxit%inext as well as boxit%ixyz
!!!>!!$    !we might then calculate the related quantities
!!!>!!$    !two dimensional index, last two elements
!!!>!!$    !to be corrected
!!!>!!$    boxit%i23=(boxit%j-1)+&
!!!>!!$         boxit%mesh%ndims(2)*(boxit%k-boxit%i3s)
!!!>!!$
!!!>!!$    !the position associated to the coordinates
!!!>!!$    boxit%rxyz(X_)=cell_r(boxit%mesh,boxit%i,X_)-boxit%oxyz(X_)
!!!>!!$    !the position associated to the coordinates
!!!>!!$    boxit%rxyz(Y_)=cell_r(boxit%mesh,boxit%j,Y_)-boxit%oxyz(Y_)
!!!>!!$    !the position associated to the coordinates
!!!>!!$    boxit%rxyz(Z_)=cell_r(boxit%mesh,boxit%k,Z_)-boxit%oxyz(Z_)
!!!>
!!!>  end subroutine update_boxit

  function box_next_point(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit
    logical :: box_next_point
    !local variables
    logical :: go

    box_next_point=associated(boxit%mesh)
    if (.not. box_next_point) return
    !this put the starting point
    !if (boxit%k==boxit%subbox(START_,Z_)-1) then
    if (boxit%inext(Z_)==boxit%subbox(START_,Z_)) then
       go=box_next_z(boxit)
       if (go) go=box_next_y(boxit)
       !if this one fails then there are no slices available
       if (.not. go) box_next_point=.false.
    end if
    !simulate loop
    flattened_loop: do
       if (box_next_x(boxit)) exit flattened_loop
       if (box_next_y(boxit)) cycle flattened_loop !and then redo the check for x
       box_next_point =box_next_z(boxit)
       if (box_next_point) box_next_point =box_next_y(boxit)
       if (.not. box_next_point) exit flattened_loop
    end do flattened_loop

  end function box_next_point

  !>split the box iterator in different tasks
  !!after the call to this routine the iterator will only run on the
  !!part corresponding to the given task.
  !!After the splitting finished the routine box_iter_merge have to be called.
  !!One cannot split an iterator more than once.
  pure subroutine box_iter_split(boxit,ntasks,itask) !other options may follow
    implicit none
    type(box_iterator), intent(inout) :: boxit
    integer, intent(in) :: ntasks,itask
    !local variables
    integer :: n,np,is

    if (ntasks==1) return
    !the present strategy splits the iterator in the direction Y
    !we might add in the following different approaches
    n=boxit%subbox(END_,Y_)-boxit%subbox(START_,Y_)+1
    call distribute_on_tasks(n,itask,ntasks,np,is)

    !then define the subbox on which the iteration has to be done
    boxit%subbox(START_,Y_)=boxit%subbox(START_,Y_)+is
    boxit%subbox(END_,Y_)=boxit%subbox(START_,Y_)+np-1

    call box_iter_rewind(boxit)
  end subroutine box_iter_split

  ! Parallelization a number n over nproc nasks
  ! this routine might go on a lower level module like f_utils
  pure subroutine distribute_on_tasks(n, iproc, nproc, np, is)
    implicit none
    ! Calling arguments
    integer,intent(in) :: n, iproc, nproc
    integer,intent(out) :: np, is

    ! Local variables
    integer :: ii

    ! First distribute evenly... (LG: if n is, say, 34 and nproc is 7 - thus 8 MPI processes)
    np = n/nproc                !(LG: we here have np=4)
    is = iproc*np               !(LG: is=iproc*4 : 0,4,8,12,16,20,24,28)
    ! ... and now distribute the remaining objects.
    ii = n-nproc*np             !(LG: ii=34-28=6)
    if (iproc<ii) np = np + 1   !(LG: the first 6 tasks (iproc<=5) will have np=5)
    is = is + min(iproc,ii)     !(LG: update is, so (iproc,np,is): (0,5,0),(1,5,5),(2,5,10),(3,5,15),(4,5,20),(5,5,25),(6,4,30),(7,4,34))

  end subroutine distribute_on_tasks

  !> Terminate the splitting section. As after the call to this routine
  !! the iterator will run on the entire nbox
  pure subroutine box_iter_merge(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit
    boxit%subbox=boxit%nbox
    call set_subbox(boxit%mesh%dom%bc,boxit%mesh%ndims,boxit%nbox,boxit%subbox)
    call box_iter_rewind(boxit)
  end subroutine box_iter_merge

  pure subroutine internal_point(bc,ipoint,npoint,jpoint,ilow,ihigh,go)
    implicit none
    integer, intent(in) :: bc
    integer, intent(in) :: npoint,ilow,ihigh,ipoint
    logical, intent(out) :: go
    integer, intent(out) :: jpoint

    if (bc == PERIODIC) then
       jpoint=modulo(ipoint-1,npoint)+1
    else
       jpoint=ipoint
    end if
    go=jpoint >= ilow
    if (go) go= jpoint <= ihigh
    if (.not. go .and. ihigh == npoint) go = ilow==1

  end subroutine internal_point

  function cell_new(dom,ndims,hgrids) result(mesh)
  !function cell_new(geocode,ndims,hgrids,alpha_bc,beta_ac,gamma_ab,abc) result(mesh)
    use numerics, only: onehalf,pi
    use wrapper_linalg, only: det_3x3
    use f_utils, only: f_assert
    use dictionaries, only: f_err_throw
    implicit none
    !character(len=1), intent(in) :: geocode
    type(domain), intent(in) :: dom !< data type for the simulation domain
    integer, dimension(3), intent(in) :: ndims
    real(gp), dimension(3), intent(in) :: hgrids
    !real(gp), dimension(3), intent(in), optional :: angrad
    !real(gp), intent(in), optional :: alpha_bc,beta_ac,gamma_ab
    !> arrays of the unit vectors of the cell. Normalized, in fortran order a_i=abc(i,1), b_i=abc(i,2)
    !real(gp), dimension(3,3), intent(in), optional :: abc
    type(cell) :: mesh
    !local variables
    !real(gp) :: aa,cc,a2,cosang
    integer :: i,j

!!$    mesh%bc=geocode_to_bc(geocode)

    mesh%dom=dom
    mesh%ndims=ndims
    mesh%hgrids=hgrids
    mesh%ndim=product(int(ndims,f_long))

!!$    !default orthorhombic
!!$    mesh%angrad=onehalf*pi
!!$
!!$    if (present(alpha_bc)) mesh%angrad(1)=alpha_bc
!!$    if (present(beta_ac)) mesh%angrad(2)=beta_ac
!!$    if (present(gamma_ab)) mesh%angrad(3)=gamma_ab
!!$
!!$    call f_assert(all(mesh%angrad > 0.0_gp),'Error, Cell new, some of the angles are not positive')
!!$
!!$    if (geocode == 'S') then
!!$       call f_assert(mesh%angrad(1)-onehalf*pi,id='Alpha angle invalid')
!!$       call f_assert(mesh%angrad(3)-onehalf*pi,id='Gamma angle invalid')
!!$    end if
!!$
!!$    mesh%orthorhombic=all(mesh%dom%angrad==onehalf*pi)
!!$
!!$    if ((geocode == 'F' .or. geocode== 'W') .and. (.not. mesh%orthorhombic)) &
!!$         call f_err_throw('For geocode="F","W" the cell must be orthorhombic')
!!$
    if (.not. dom%orthorhombic) then
!!$       !some consistency check on the angles should be performed
!!$       !1) sum(angrad) < twopi
!!$       if (all(mesh%dom%angrad==mesh%dom%angrad(1))) then
!!$          !Treat the case of equal angles (except all right angles) :
!!$          !generates trigonal symmetry wrt third axis
!!$          cosang=cos(mesh%dom%angrad(1))
!!$          a2=2.0_gp/3.0_gp*(1.0_gp-cosang)
!!$          aa=sqrt(a2)
!!$          cc=sqrt(1.0_gp-a2)
!!$          mesh%habc(1,1)=aa; mesh%habc(2,1)=0.0_gp; mesh%habc(3,1)=cc
!!$          mesh%habc(1,2)=-0.5_gp*aa ; mesh%habc(2,2)=sqrt(3.0_gp)*0.5_gp*aa ; mesh%habc(3,2)=cc
!!$          mesh%habc(1,3)=-0.5_gp*aa ; mesh%habc(2,3)=-sqrt(3.0_gp)*0.5_gp*aa ; mesh%habc(3,3)=cc
!!$          !Set the covariant metric
!!$          mesh%gd(1,1) = 1.0_gp
!!$          mesh%gd(1,2) = cos(mesh%dom%angrad(3)) !gamma_ab
!!$          mesh%gd(1,3) = cos(mesh%dom%angrad(2)) !beta_ac
!!$          mesh%gd(2,2) = 1.0_gp
!!$          mesh%gd(2,3) = cos(mesh%dom%angrad(1)) !alpha_bc
!!$          mesh%gd(3,3) = 1.0_gp
!!$          !Set the determinant of the covariant metric
!!$          mesh%detgd = 1.0_gp - cos(mesh%dom%angrad(1))**2 - cos(mesh%dom%angrad(2))**2 - cos(mesh%dom%angrad(3))**2 +&
!!$               2.0_gp*cos(mesh%dom%angrad(1))*cos(mesh%dom%angrad(2))*cos(mesh%dom%angrad(3))
!!$          !Set the contravariant metric
!!$          mesh%gu(1,1) = (sin(mesh%dom%angrad(1))**2)/mesh%detgd
!!$          mesh%gu(1,2) = (cos(mesh%dom%angrad(2))*cos(mesh%dom%angrad(1))-cos(mesh%dom%angrad(3)))/mesh%detgd
!!$          mesh%gu(1,3) = (cos(mesh%dom%angrad(3))*cos(mesh%dom%angrad(1))-cos(mesh%dom%angrad(2)))/mesh%detgd
!!$          mesh%gu(2,2) = (sin(mesh%dom%angrad(2))**2)/mesh%detgd
!!$          mesh%gu(2,3) = (cos(mesh%dom%angrad(3))*cos(mesh%dom%angrad(2))-cos(mesh%dom%angrad(1)))/mesh%detgd
!!$          mesh%gu(3,3) = (sin(mesh%dom%angrad(3))**2)/mesh%detgd
!!$       else if (geocode == 'P') then
!!$          mesh%habc=0.0_gp
!!$          mesh%habc(1,1)=1.0_gp
!!$          mesh%habc(1,2)=cos(mesh%dom%angrad(3))
!!$          mesh%habc(2,2)=sin(mesh%dom%angrad(3))
!!$          mesh%habc(1,3)=cos(mesh%dom%angrad(2))
!!$          mesh%habc(2,3)=(cos(mesh%dom%angrad(1))-mesh%habc(1,2)*mesh%habc(1,3))/mesh%habc(2,2)
!!$          mesh%habc(3,3)=sqrt(1.0_gp-mesh%habc(1,3)**2-mesh%habc(2,3)**2)
!!$          !Set the covariant metric
!!$          mesh%gd(1,1) = 1.0_gp
!!$          mesh%gd(1,2) = cos(mesh%dom%angrad(3)) !gamma_ab
!!$          mesh%gd(1,3) = cos(mesh%dom%angrad(2)) !beta_ac
!!$          mesh%gd(2,2) = 1.0_gp
!!$          mesh%gd(2,3) = cos(mesh%dom%angrad(1)) !alpha_bc
!!$          mesh%gd(3,3) = 1.0_gp
!!$          !Set the determinant of the covariant metric
!!$          mesh%detgd = 1.0_gp - cos(mesh%dom%angrad(1))**2 - cos(mesh%dom%angrad(2))**2 - cos(mesh%dom%angrad(3))**2 +&
!!$               2.0_gp*cos(mesh%dom%angrad(1))*cos(mesh%dom%angrad(2))*cos(mesh%dom%angrad(3))
!!$          !Set the contravariant metric
!!$          mesh%gu(1,1) = (sin(mesh%dom%angrad(1))**2)/mesh%detgd
!!$          mesh%gu(1,2) = (cos(mesh%dom%angrad(2))*cos(mesh%dom%angrad(1))-cos(mesh%dom%angrad(3)))/mesh%detgd
!!$          mesh%gu(1,3) = (cos(mesh%dom%angrad(3))*cos(mesh%dom%angrad(1))-cos(mesh%dom%angrad(2)))/mesh%detgd
!!$          mesh%gu(2,2) = (sin(mesh%dom%angrad(2))**2)/mesh%detgd
!!$          mesh%gu(2,3) = (cos(mesh%dom%angrad(3))*cos(mesh%dom%angrad(2))-cos(mesh%dom%angrad(1)))/mesh%detgd
!!$          mesh%gu(3,3) = (sin(mesh%dom%angrad(3))**2)/mesh%detgd
!!$       else !only Surfaces is possible here
!!$          mesh%habc=0.0_gp
!!$          mesh%habc(1,1)=1.0_gp
!!$          mesh%habc(2,2)=1.0_gp
!!$          mesh%habc(1,3)=cos(mesh%dom%angrad(2))
!!$          mesh%habc(3,3)=sin(mesh%dom%angrad(2))
!!$          !Set the covariant metric
!!$          mesh%gd=0.0_gp
!!$          mesh%gd(1,1) = 1.0_gp
!!$          mesh%gd(1,3) = cos(mesh%dom%angrad(2)) !beta_ac
!!$          mesh%gd(2,2) = 1.0_gp
!!$          mesh%gd(3,3) = 1.0_gp
!!$          !Set the determinant of the covariant metric
!!$          mesh%detgd = sin(mesh%dom%angrad(2))**2
!!$          !Set the contravariant metric
!!$          mesh%gu=0.0_gp
!!$          mesh%gu(1,1) = 1.0_gp/mesh%detgd
!!$          mesh%gu(1,3) = -cos(mesh%dom%angrad(2))/mesh%detgd
!!$          mesh%gu(2,2) = 1.0_gp!/mesh%detgd
!!$          mesh%gu(3,3) = 1.0_gp/mesh%detgd
!!$       end if
!!$       mesh%uabc=0.0_gp
!!$       mesh%uabc(1:3,1:3)=mesh%habc(1:3,1:3)
!!$
!!$       !Rescale habc using hgrid
!!$       mesh%habc(:,1)=hgrids*mesh%habc(:,1)
!!$       mesh%habc(:,2)=hgrids*mesh%habc(:,2)
!!$       mesh%habc(:,3)=hgrids*mesh%habc(:,3)
       ! here we assume the dom%uabc = mesh%uabc
       mesh%habc(:,1)=hgrids*dom%uabc(:,1)
       mesh%habc(:,2)=hgrids*dom%uabc(:,2)
       mesh%habc(:,3)=hgrids*dom%uabc(:,3)
       !the volume element
       !Compute unit cell volume
       mesh%volume_element=det_3x3(mesh%habc)
    else
       mesh%habc=0.0_gp
!!$       mesh%uabc=0.0_gp
       do i=1,3
          mesh%habc(i,i)=hgrids(i)
!!$          mesh%uabc(i,i)=1.0_gp
       end do
!!$       mesh%dom%angrad=onehalf*pi
       mesh%volume_element=product(mesh%hgrids)
!!$       mesh%gd(1,1) = 1.0_gp
!!$       mesh%gd(1,2) = 0.0_gp
!!$       mesh%gd(1,3) = 0.0_gp
!!$       mesh%gd(2,2) = 1.0_gp
!!$       mesh%gd(2,3) = 0.0_gp
!!$       mesh%gd(3,3) = 1.0_gp
!!$       mesh%detgd = 1.0_gp
!!$       !Set the contravariant metric
!!$       mesh%gu(1,1) = 1.0_gp
!!$       mesh%gu(1,2) = 0.0_gp
!!$       mesh%gu(1,3) = 0.0_gp
!!$       mesh%gu(2,2) = 1.0_gp
!!$       mesh%gu(2,3) = 0.0_gp
!!$       mesh%gu(3,3) = 1.0_gp
    end if
!!$    mesh%gd(2,1) = mesh%gd(1,2)
!!$    mesh%gd(3,1) = mesh%gd(1,3)
!!$    mesh%gd(3,2) = mesh%gd(2,3)
!!$
!!$    mesh%gu(2,1) = mesh%gu(1,2)
!!$    mesh%gu(3,1) = mesh%gu(1,3)
!!$    mesh%gu(3,2) = mesh%gu(2,3)
    do i=1,3
       do j=1,3
          if (abs(mesh%habc(i,j)).lt.1.0d-15) mesh%habc(i,j)=0.0_gp
!!$          if (abs(mesh%uabc(i,j)).lt.1.0d-15) mesh%uabc(i,j)=0.0_gp
!!$          if (abs(mesh%gd(i,j)).lt.1.0d-15) mesh%gd(i,j)=0.0_gp
!!$          if (abs(mesh%gu(i,j)).lt.1.0d-15) mesh%gu(i,j)=0.0_gp
       end do
    end do

    !here we should verify that the the inverse metric times the metric is the identity

  end function cell_new

  !>gives the value of the coordinate from the grid point
  elemental pure function cell_r(mesh,i,dim) result(t)
    implicit none
    integer, intent(in) :: i
    type(cell), intent(in) :: mesh
    integer, intent(in) :: dim
    real(gp) :: t

    t=mesh%hgrids(dim)*(i-1)
  end function cell_r

  pure function box_iter_square_gd(bit) 
    implicit none
    type(box_iterator), intent(in) :: bit
    real(gp) :: box_iter_square_gd

    box_iter_square_gd=square_gd(bit%mesh%dom,bit%rxyz)
    
  end function box_iter_square_gd

  pure function box_iter_closest_r(bit,rxyz,orthorhombic) result (r)
    !use f_utils, only: f_get_option
    implicit none
    type(box_iterator), intent(in) :: bit
    real(gp), dimension(3), intent(in) :: rxyz
    logical, intent(in), optional :: orthorhombic
    real(gp), dimension(3) :: r
    ! local variables
    logical :: orthorhombic_

    !orthorhombic_=f_get_option(orthorhombic,.false.)
    orthorhombic_=.false.
    if (present(orthorhombic)) orthorhombic_=orthorhombic

    r=closest_r(bit%mesh%dom,bit%rxyz,rxyz)

    if (orthorhombic_) r=rxyz_ortho(bit%mesh%dom,r)
    
  end function box_iter_closest_r

  pure function box_iter_distance(bit,rxyz0) 
    implicit none
    type(box_iterator), intent(in) :: bit
    real(gp), dimension(3), intent(in) :: rxyz0
    real(gp) :: box_iter_distance

    box_iter_distance=distance(bit%mesh%dom,bit%rxyz,rxyz0)
    
  end function box_iter_distance


end module box
