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

  private

  !>parameter for the definition of the bc
  integer, parameter :: FREE=0
  integer, parameter :: PERIODIC=1
  
  !to ease readiness
  integer, parameter :: START_=1,END_=2
  integer, parameter :: X_=1,Y_=2,Z_=3

  !> data type which stores all informations of the simulation box. 
  !! It contains also the metric for nonorthorhombic cells.
  type, public :: cell
     logical :: orthorhombic !<true if the cell is orthorhombic
     integer, dimension(3) :: bc !< boundary conditions on each direction (FREE=0, PERIODIC=1)
     integer, dimension(3) :: ndims !< number of grid points on each direction
     real(gp), dimension(3) :: hgrids !< real space grid on each direction
     real(gp), dimension(3) :: angrad !<angles between the dimensions in radiant (alpha_bc,beta_ac,gamma_bc)
     !derived data
     integer(f_long) :: ndim !< product of the dimension, long integer to avoid overflow
     real(gp) :: volume_element !< volume element of the primitive cell
     real(gp), dimension(3,3) :: habc !<primitive volume elements in the translation vectors direction
     real(gp), dimension(3,3) :: uabc !<matrix of the normalized translation vectors direction
     real(gp), dimension(3,3) :: gd !<covariant metric needed for non-orthorhombic operations
     real(gp), dimension(3,3) :: gu !<controvariant metric needed for non-orthorhombic operations
     real(gp) :: detgd !<determinant of the covariant matrix
  end type cell

  !> defines the object to iterate around the real-space grid points.
  !! given a cell type, it might iterate on a section of this gris provided by the extremes nbox
  !! it also provides a facility to parallelize over the 
  type, public :: box_iterator
     integer :: i3s !<starting point in the dimension z
     integer :: i3e !<ending point in the dimension z
     integer :: i23 !<collapsed index in 23 dimension (in relative conventions)
     integer :: ind !<one-dimensional index for arrays (in relative conventions)
     !> indices in absolute coordinates in the given box,
     !! from nbox(1,:) to nbox(2,:). To be intended as private 
     integer, dimension(3)  :: inext  
     !> actual index inside the box,from 1 to mesh%ndims(:)
     integer :: i,j,k !better as scalars
     !> Sub-box to iterate over the points (ex. around atoms)
     !! start and end points for each direction
     integer, dimension(2,3) :: nbox 
     real(gp), dimension(3) :: oxyz !<origin of the coordinate system
     real(gp), dimension(3) :: rxyz !<coordinates of the grid point
     real(gp), dimension(3) :: tmp !< size 3 array buffer to avoid the creation of temporary arrays
     logical :: whole !<to assess if we run over the entire box or not (no check over the internal point)
     integer, dimension(2,3) :: subbox !<box of the local task
     !>reference mesh from which it starts
     type(cell), pointer :: mesh
  end type box_iterator

!!$  interface box_iter
!!$     module procedure box_iter_c,box_iter_base
!!$  end interface box_iter

  interface dotp_gu
     module procedure dotp_gu,dotp_gu_add1,dotp_gu_add2
  end interface dotp_gu

  interface square_gu
     module procedure square,square_add
  end interface square_gu

  interface dotp_gd
     module procedure dotp_gd,dotp_gd_add1,dotp_gd_add2
  end interface dotp_gd


  interface square_gd
     module procedure square_gd,square_gd_add
  end interface square_gd

  public :: cell_r,cell_periodic_dims,rxyz_ortho,distance,closest_r,square_gu,square_gd,cell_new,box_iter,box_next_point
  public :: cell_geocode,box_next_x,box_next_y,box_next_z,dotp_gu,dotp_gd,cell_null,nullify_box_iterator
  public :: box_iter_rewind,box_iter_split,box_iter_merge,box_iter_set_nbox,box_iter_expand_nbox,box_nbox_from_cutoff


contains

  !> Nullify the cell type
  pure function cell_null() result(me)
   implicit none
   type(cell) :: me
   me%orthorhombic=.true.
   me%bc=0
   me%ndims=0
   me%hgrids=0.0_gp
   me%angrad=0.0_gp
   !derived data
   me%ndim=0
   me%volume_element=0.0_gp
   me%habc=0.0_gp
   me%uabc=0.0_gp
   me%gd=0.0_gp
   me%gu=0.0_gp
   me%detgd=0.0_gp
  end function cell_null

  !> Nullify the iterator dpbox type
  pure subroutine nullify_box_iterator(boxit)
    implicit none
    type(box_iterator), intent(out) :: boxit
    boxit%i3s =-1 
    boxit%i3e =-1
    boxit%i23 =-1
    boxit%ind =-1
    boxit%i=-1
    boxit%j=-1
    boxit%k=-1
    boxit%inext=0
    boxit%inext(X_)=-1
    boxit%nbox=-1 
    boxit%oxyz=-1.0_gp
    boxit%rxyz=-1.0_gp
    boxit%tmp=0.0_gp
    nullify(boxit%mesh)
    boxit%whole=.false.
    boxit%subbox=-1
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
       call set_subbox(bit%mesh%bc,bit%mesh%ndims,bit%nbox,bit%subbox)
       call box_iter_rewind(bit)
    else if (present(cutoff)) then
!!$       
!!$       bit%nbox(START_,:)=floor((oxyz-cutoff)/bit%mesh%hgrids)
!!$       bit%nbox(END_,:)=ceiling((oxyz+cutoff)/bit%mesh%hgrids)
       bit%nbox=box_nbox_from_cutoff(bit%mesh,oxyz,cutoff)
       bit%whole=.false.
       call set_subbox(bit%mesh%bc,bit%mesh%ndims,bit%nbox,bit%subbox)
       call box_iter_rewind(bit)
    else
       call box_iter_expand_nbox(bit)
    end if

    call box_iter_rewind(bit)

  end subroutine box_iter_set_nbox

  !> this function has to be genralized for non-orthorhombic grids
  pure function box_nbox_from_cutoff(mesh,oxyz,cutoff) result(nbox)
    implicit none
    type(cell), intent(in) :: mesh
    real(gp), dimension(3), intent(in) :: oxyz
    real(gp), intent(in) :: cutoff
    integer, dimension(2,3) :: nbox
    !for non-orthorhombic cells the concept of distance has to be inserted here (the box should contain the sphere)
    nbox(START_,:)=floor((oxyz-cutoff)/mesh%hgrids)
    nbox(END_,:)=ceiling((oxyz+cutoff)/mesh%hgrids)

  end function box_nbox_from_cutoff

  pure subroutine box_iter_expand_nbox(bit)
    implicit none
    type(box_iterator), intent(inout) :: bit
    bit%whole=.true.
    bit%nbox(START_,:)=1
    bit%nbox(END_,:)=bit%mesh%ndims
    call set_subbox(bit%mesh%bc,bit%mesh%ndims,bit%nbox,bit%subbox)
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
    integer :: iz,iy,ix,i,jx,jy,jz
    integer(f_long) :: icnt,itgt
    logical(f_byte), dimension(:), allocatable :: lxyz
    integer, dimension(:), allocatable :: lx,ly,lz
    integer, dimension(3) :: subdims

    do i=1,3
       subdims(i)=bit%subbox(END_,i)-bit%subbox(START_,i)+1
    end do
    subdims(3)=min(subdims(3),bit%i3e-bit%i3s+1)

    !first, count if the iterator covers all the required points
    itgt=product(int(subdims,f_long))
    !!!itgt=int(bit%mesh%ndims(1),f_long)*int(bit%mesh%ndims(2),f_long)*&
    !!!     int(bit%i3e-bit%i3s+1,f_long)

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
       lz(iz)=iz+bit%subbox(START_,Z_)-1
    end do
    do iy=1,subdims(Y_)
       ly(iy)=iy+bit%subbox(START_,Y_)-1
    end do
    do ix=1,subdims(X_)
       lx(ix)=ix+bit%subbox(START_,X_)-1
    end do

    !separable mode
    iz=bit%subbox(START_,Z_)-1 !0
    icnt=0
    do while(box_next_z(bit))
       iz=iz+1
       jz=iz-bit%subbox(START_,Z_)+1
       !print *,'bit',bit%k,bit%inext(Z_)
       call f_assert(iz+bit%i3s-1==bit%inext(Z_)-1,'A')!,&
       !'Error iz='+iz+', inext(Z)='+bit%inext(Z_))
       iy=bit%subbox(START_,Y_)-1!0
       do while(box_next_y(bit))
          iy=iy+1
          jy=iy-bit%subbox(START_,Y_)+1
          call f_assert(iy==bit%inext(Y_)-1,'B')!,&
          !'Error iy='+iy+', inext(Y)='+bit%inext(Y_))
          ix=bit%subbox(START_,X_)-1!0
          do while(box_next_x(bit))
             ix=ix+1
             jx=ix-bit%subbox(START_,X_)+1
             call f_assert(ix==bit%inext(X_)-1,'C')!,&
             !'Error ix='+ix+', inext(X)='+bit%inext(X_))

             icnt=icnt+1
             call f_assert(lx(jx) == bit%i,'D')!,&
             !'Error value, ix='+bit%i+', expected='+lx(ix))
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

       call f_assert(lz(jz)+bit%i3s-1 == bit%k,&
            yaml_toa([lz(jz),bit%k,bit%i3s]))!,&
    end do
    call f_assert(jz == subdims(Z_),'I')!,&
    !'Error boxit, iz='+iz+', itgtz='+subdims(Z_))
    call f_assert(icnt == itgt,'J')!,&
    !'Error sep boxit, icnt='+icnt+', itgt='+itgt)

    !complete mode
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

    if (any(lxyz)) call f_err_throw('Error boxit, points not covered')

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

!!$    call increment_dim(bit,3,bit%k,ok)
    ok = bit%i3e >= bit%i3s ! to be removed
    ok = bit%subbox(END_,Z_) >= bit%subbox(START_,Z_) 
    if (.not. ok) return !there is nothing to explore
    ok= bit%inext(Z_) <= bit%subbox(END_,Z_)
    do while(ok)
       if (bit%whole) then
          bit%k=bit%inext(Z_)
       else 
          call internal_point(bit%mesh%bc(Z_),bit%inext(Z_),bit%mesh%ndims(Z_),&
               bit%k,bit%i3s,bit%i3e,ok)
          if (.not. ok) bit%inext(Z_)=bit%inext(Z_)+1
       end if
       if (ok) then
          bit%inext(Z_)=bit%inext(Z_)+1
          exit
       end if
       ok = bit%inext(Z_) <= bit%subbox(END_,Z_)
    end do
    !reset x and y
    if (ok) then
       call update_boxit_z(bit)
       bit%inext(Y_)=bit%subbox(START_,Y_)
       bit%inext(X_)=bit%subbox(START_,X_)
    end if

    !in the case the z_direction is over, make the iterator ready for new use
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
       if (bit%mesh%bc(idim) == PERIODIC) then
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

  end subroutine update_boxit_y

  pure subroutine update_boxit_z(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit 

    !the position associated to the coordinates
    boxit%rxyz(Z_)=cell_r(boxit%mesh,boxit%k,Z_)-boxit%oxyz(Z_)

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
    if (boxit%k==boxit%subbox(START_,Z_)-1) then
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
    call set_subbox(boxit%mesh%bc,boxit%mesh%ndims,boxit%nbox,boxit%subbox)
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

  end subroutine internal_point

  function cell_new(geocode,ndims,hgrids,alpha_bc,beta_ac,gamma_ab,abc) result(mesh)
    use numerics, only: onehalf,pi
    use wrapper_linalg, only: det_3x3
    use f_utils, only: f_assert
    use dictionaries, only: f_err_throw
    implicit none
    character(len=1), intent(in) :: geocode
    integer, dimension(3), intent(in) :: ndims
    real(gp), dimension(3), intent(in) :: hgrids
    !real(gp), dimension(3), intent(in), optional :: angrad
    real(gp), intent(in), optional :: alpha_bc,beta_ac,gamma_ab
    !> arrays of the unit vectors of the cell. Normalized, in fortran order a_i=abc(i,1), b_i=abc(i,2)
    real(gp), dimension(3,3), intent(in), optional :: abc
    type(cell) :: mesh
    !local variables
    real(gp) :: aa,cc,a2,cosang
    integer :: i,j


    select case(geocode)
    case('P')
       mesh%bc=PERIODIC
    case('S')
       mesh%bc=PERIODIC
       mesh%bc(2)=FREE
    case('F')
       mesh%bc=FREE
    case('W')
       mesh%bc=FREE
       mesh%bc(3)=PERIODIC
    case default
       call f_err_throw('Invalid specification of the variable "geocode"')
    end select
    mesh%ndims=ndims
    mesh%hgrids=hgrids
    mesh%ndim=product(int(ndims,f_long))

    !default orthorhombic
    mesh%angrad=onehalf*pi

    if (present(alpha_bc)) mesh%angrad(1)=alpha_bc
    if (present(beta_ac)) mesh%angrad(2)=beta_ac
    if (present(gamma_ab)) mesh%angrad(3)=gamma_ab

    call f_assert(all(mesh%angrad > 0.0_gp),'Error, Cell new, some of the angles are not positive')

    if (geocode == 'S') then
       call f_assert(mesh%angrad(1)-onehalf*pi,id='Alpha angle invalid')
       call f_assert(mesh%angrad(3)-onehalf*pi,id='Gamma angle invalid')
    end if

    mesh%orthorhombic=all(mesh%angrad==onehalf*pi)

    if ((geocode == 'F' .or. geocode== 'W') .and. (.not. mesh%orthorhombic)) &
         call f_err_throw('For geocode="F","W" the cell must be orthorhombic')

    if (.not. mesh%orthorhombic) then
       !some consistency check on the angles should be performed
       !1) sum(angrad) < twopi
       if (all(mesh%angrad==mesh%angrad(1))) then
          !Treat the case of equal angles (except all right angles) :
          !generates trigonal symmetry wrt third axis
          cosang=cos(mesh%angrad(1))
          a2=2.0_gp/3.0_gp*(1.0_gp-cosang)
          aa=sqrt(a2)
          cc=sqrt(1.0_gp-a2)
          mesh%habc(1,1)=aa; mesh%habc(2,1)=0.0_gp; mesh%habc(3,1)=cc
          mesh%habc(1,2)=-0.5_gp*aa ; mesh%habc(2,2)=sqrt(3.0_gp)*0.5_gp*aa ; mesh%habc(3,2)=cc
          mesh%habc(1,3)=-0.5_gp*aa ; mesh%habc(2,3)=-sqrt(3.0_gp)*0.5_gp*aa ; mesh%habc(3,3)=cc
          !Set the covariant metric
          mesh%gd(1,1) = 1.0_gp
          mesh%gd(1,2) = cos(mesh%angrad(3)) !gamma_ab
          mesh%gd(1,3) = cos(mesh%angrad(2)) !beta_ac
          mesh%gd(2,2) = 1.0_gp
          mesh%gd(2,3) = cos(mesh%angrad(1)) !alpha_bc
          mesh%gd(3,3) = 1.0_gp 
          !Set the determinant of the covariant metric
          mesh%detgd = 1.0_gp - cos(mesh%angrad(1))**2 - cos(mesh%angrad(2))**2 - cos(mesh%angrad(3))**2 +&
               2.0_gp*cos(mesh%angrad(1))*cos(mesh%angrad(2))*cos(mesh%angrad(3))
          !Set the contravariant metric
          mesh%gu(1,1) = (sin(mesh%angrad(1))**2)/mesh%detgd
          mesh%gu(1,2) = (cos(mesh%angrad(2))*cos(mesh%angrad(1))-cos(mesh%angrad(3)))/mesh%detgd
          mesh%gu(1,3) = (cos(mesh%angrad(3))*cos(mesh%angrad(1))-cos(mesh%angrad(2)))/mesh%detgd
          mesh%gu(2,2) = (sin(mesh%angrad(2))**2)/mesh%detgd
          mesh%gu(2,3) = (cos(mesh%angrad(3))*cos(mesh%angrad(2))-cos(mesh%angrad(1)))/mesh%detgd
          mesh%gu(3,3) = (sin(mesh%angrad(3))**2)/mesh%detgd
       else if (geocode == 'P') then
          mesh%habc=0.0_gp
          mesh%habc(1,1)=1.0_gp
          mesh%habc(1,2)=cos(mesh%angrad(3))
          mesh%habc(2,2)=sin(mesh%angrad(3))
          mesh%habc(1,3)=cos(mesh%angrad(2))
          mesh%habc(2,3)=(cos(mesh%angrad(1))-mesh%habc(1,2)*mesh%habc(1,3))/mesh%habc(2,2)
          mesh%habc(3,3)=sqrt(1.0_gp-mesh%habc(1,3)**2-mesh%habc(2,3)**2)
          !Set the covariant metric
          mesh%gd(1,1) = 1.0_gp
          mesh%gd(1,2) = cos(mesh%angrad(3)) !gamma_ab
          mesh%gd(1,3) = cos(mesh%angrad(2)) !beta_ac
          mesh%gd(2,2) = 1.0_gp
          mesh%gd(2,3) = cos(mesh%angrad(1)) !alpha_bc
          mesh%gd(3,3) = 1.0_gp 
          !Set the determinant of the covariant metric
          mesh%detgd = 1.0_gp - cos(mesh%angrad(1))**2 - cos(mesh%angrad(2))**2 - cos(mesh%angrad(3))**2 +&
               2.0_gp*cos(mesh%angrad(1))*cos(mesh%angrad(2))*cos(mesh%angrad(3))
          !Set the contravariant metric
          mesh%gu(1,1) = (sin(mesh%angrad(1))**2)/mesh%detgd
          mesh%gu(1,2) = (cos(mesh%angrad(2))*cos(mesh%angrad(1))-cos(mesh%angrad(3)))/mesh%detgd
          mesh%gu(1,3) = (cos(mesh%angrad(3))*cos(mesh%angrad(1))-cos(mesh%angrad(2)))/mesh%detgd
          mesh%gu(2,2) = (sin(mesh%angrad(2))**2)/mesh%detgd
          mesh%gu(2,3) = (cos(mesh%angrad(3))*cos(mesh%angrad(2))-cos(mesh%angrad(1)))/mesh%detgd
          mesh%gu(3,3) = (sin(mesh%angrad(3))**2)/mesh%detgd
       else !only Surfaces is possible here
          mesh%habc=0.0_gp
          mesh%habc(1,1)=1.0_gp
          mesh%habc(2,2)=1.0_gp
          mesh%habc(1,3)=cos(mesh%angrad(2))
          mesh%habc(3,3)=sin(mesh%angrad(2))
          !Set the covariant metric
          mesh%gd=0.0_gp
          mesh%gd(1,1) = 1.0_gp
          mesh%gd(1,3) = cos(mesh%angrad(2)) !beta_ac
          mesh%gd(2,2) = 1.0_gp
          mesh%gd(3,3) = 1.0_gp 
          !Set the determinant of the covariant metric
          mesh%detgd = sin(mesh%angrad(2))**2 
          !Set the contravariant metric
          mesh%gu=0.0_gp
          mesh%gu(1,1) = 1.0_gp/mesh%detgd
          mesh%gu(1,3) = -cos(mesh%angrad(2))/mesh%detgd
          mesh%gu(2,2) = 1.0_gp!/mesh%detgd
          mesh%gu(3,3) = 1.0_gp/mesh%detgd
       end if
       mesh%uabc=0.0_gp
       mesh%uabc(1:3,1:3)=mesh%habc(1:3,1:3)

       !Rescale habc using hgrid
       mesh%habc(:,1)=hgrids*mesh%habc(:,1)
       mesh%habc(:,2)=hgrids*mesh%habc(:,2)
       mesh%habc(:,3)=hgrids*mesh%habc(:,3)
       !the volume element
       !Compute unit cell volume
       mesh%volume_element=det_3x3(mesh%habc)
    else
       mesh%habc=0.0_gp
       mesh%uabc=0.0_gp
       do i=1,3
          mesh%habc(i,i)=hgrids(i)
          mesh%uabc(i,i)=1.0_gp
       end do
       mesh%angrad=onehalf*pi
       mesh%volume_element=product(mesh%hgrids)
       mesh%gd(1,1) = 1.0_gp
       mesh%gd(1,2) = 0.0_gp
       mesh%gd(1,3) = 0.0_gp
       mesh%gd(2,2) = 1.0_gp
       mesh%gd(2,3) = 0.0_gp
       mesh%gd(3,3) = 1.0_gp 
       mesh%detgd = 1.0_gp 
       !Set the contravariant metric
       mesh%gu(1,1) = 1.0_gp
       mesh%gu(1,2) = 0.0_gp
       mesh%gu(1,3) = 0.0_gp
       mesh%gu(2,2) = 1.0_gp
       mesh%gu(2,3) = 0.0_gp 
       mesh%gu(3,3) = 1.0_gp
    end if
    mesh%gd(2,1) = mesh%gd(1,2)
    mesh%gd(3,1) = mesh%gd(1,3)
    mesh%gd(3,2) = mesh%gd(2,3)

    mesh%gu(2,1) = mesh%gu(1,2)
    mesh%gu(3,1) = mesh%gu(1,3)
    mesh%gu(3,2) = mesh%gu(2,3)
    do i=1,3
       do j=1,3
          if (abs(mesh%habc(i,j)).lt.1.0d-15) mesh%habc(i,j)=0.0_gp
          if (abs(mesh%uabc(i,j)).lt.1.0d-15) mesh%uabc(i,j)=0.0_gp
          if (abs(mesh%gd(i,j)).lt.1.0d-15) mesh%gd(i,j)=0.0_gp
          if (abs(mesh%gu(i,j)).lt.1.0d-15) mesh%gu(i,j)=0.0_gp
       end do
    end do

    !here we should verify that the the inverse metric times the metric is the identity



  end function cell_new

  !> returns a logical array of size 3 which is .true. for all the periodic dimensions
  pure function cell_periodic_dims(mesh) result(peri)
    implicit none
    type(cell), intent(in) :: mesh
    logical, dimension(3) :: peri
    !local variables

    peri= mesh%bc == PERIODIC

  end function cell_periodic_dims

  !>give the associated geocode, 'X' for unknown
  pure function cell_geocode(mesh)
    implicit none
    type(cell), intent(in) :: mesh
    character(len=1) :: cell_geocode
    !local variables
    logical, dimension(3) :: peri

    peri=cell_periodic_dims(mesh)
    if (all(peri)) then
       cell_geocode='P'
    else if (.not. any(peri)) then
       cell_geocode='F'
    else if (peri(1) .and. .not. peri(2) .and. peri(3)) then
       cell_geocode='S'
    else if (.not. peri(1) .and. .not. peri(2) .and. peri(3)) then
       cell_geocode='W'
    else
       cell_geocode='X'
    end if

  end function cell_geocode


  !>gives the value of the coordinate from the grid point
  elemental pure function cell_r(mesh,i,dim) result(t)
    implicit none
    integer, intent(in) :: i
    type(cell), intent(in) :: mesh
    integer, intent(in) :: dim
    real(gp) :: t

    t=mesh%hgrids(dim)*(i-1)
  end function cell_r

  !>gives the value of the coordinates for an orthorhombic reference system
  pure function rxyz_ortho(mesh,rxyz)
    implicit none
    type(cell), intent(in) :: mesh
    real(gp), dimension(3), intent(in) :: rxyz
    real(gp), dimension(3) :: rxyz_ortho
    ! local variables
    integer :: i,j

    if (mesh%orthorhombic) then
     rxyz_ortho(1:3)=rxyz(1:3)
    else
     do i=1,3
      rxyz_ortho(i)=0.0_gp
      do j=1,3
       rxyz_ortho(i)=rxyz_ortho(i)+mesh%uabc(i,j)*rxyz(j)
      end do
     end do
    end if

  end function rxyz_ortho

  pure function distance(mesh,r,c) result(d)
    use dictionaries, only: f_err_throw
    implicit none
    real(gp), dimension(3), intent(in) :: r,c
    type(cell), intent(in) :: mesh
    real(gp) :: d
    !local variables
    integer :: i !,j,k,ii
    real(gp) :: d2!,dold
    real(gp), dimension(3) :: rt!,ri,ci

    rt=closest_r(mesh,r,c)
    d2=square_gd(mesh,rt)
    d=sqrt(d2)

    d=0.0_gp
    if (mesh%orthorhombic) then
       d2=0.0_gp
       do i=1,3
          d2=d2+r_wrap(mesh%bc(i),mesh%hgrids(i)*mesh%ndims(i),&
               r(i),c(i))**2
       end do
       d=sqrt(d2)
    else
       rt=closest_r(mesh,r,c)
       d2=square_gd(mesh,rt)
       d=sqrt(d2)
!       dold=1.0d100 !huge_number
!       do ii=1,3
!        if (mesh%bc(ii)==PERIODIC) then
!          ri(ii)=mod(r(ii),mesh%ndims(ii)*mesh%hgrids(ii))
!          ci(ii)=mod(c(ii),mesh%ndims(ii)*mesh%hgrids(ii))
!        else
!          ri(ii)=r(ii)
!          ci(ii)=c(ii)
!        end if
!       end do
!       do i=-mesh%bc(1),mesh%bc(1)
!        do j=-mesh%bc(2),mesh%bc(2)
!         do k=-mesh%bc(3),mesh%bc(3)
!            rt(1)=ri(1)+real(i,kind=8)*mesh%ndims(1)*mesh%hgrids(1)
!            rt(2)=ri(2)+real(j,kind=8)*mesh%ndims(2)*mesh%hgrids(2)
!            rt(3)=ri(3)+real(k,kind=8)*mesh%ndims(3)*mesh%hgrids(3)
!            d2=square_gd(mesh,rt-ci)
!            d=sqrt(d2)
!            if (d.lt.dold) then
!               dold=d
!            end if
!         end do
!        end do
!       end do
!       d=dold
    end if

  end function distance

!!$  pure function min_dist(bc,alat,r,r_old)
!!$    implicit none
!!$    integer, intent(in) :: bc
!!$    real(gp), intent(in) :: r,r_old,alat
!!$    real(gp) :: min_dist
!!$
!!$    !for periodic BC calculate mindist only if the center of mass can be defined without the modulo
!!$    min_dist=abs(r-r_old)
!!$    if (bc==PERIODIC) then
!!$       if (min_dist > 0.5_gp*alat) then
!!$          if (r < 0.5_gp*alat) then
!!$             min_dist=abs(r+alat-r_old)
!!$          else
!!$             min_dist=abs(r-alat-r_old)
!!$          end if
!!$       end if
!!$    end if
!!$
!!$  end function min_dist

  !> Calculates the minimum difference between two coordinates
  !!@warning: this is only valid if the coordinates wrap once.
  pure function r_wrap(bc,alat,ri,ci)
    implicit none
    integer, intent(in) :: bc
    real(gp), intent(in) :: ri,ci,alat
    real(gp) :: r_wrap
    ! local variables
    real(gp) :: r,c

    !for periodic BC calculate mindist only if the center of mass can be defined without the modulo
    r=ri
    c=ci
    r_wrap=r-c
    if (bc==PERIODIC) then
      r=mod(ri,alat)
      c=mod(ci,alat)
      r_wrap=r-c
       if (abs(r_wrap) > 0.5_gp*alat) then
          if (r < 0.5_gp*alat) then
             r_wrap=r+alat-c
          else
             r_wrap=r-alat-c
          end if
       end if
    end if

  end function r_wrap

  !>find the closest center according to the periodiciy of the
  !! box and provide the vector
  pure function closest_r(mesh,v,center) result(r)
    implicit none
    real(gp), dimension(3), intent(in) :: v,center
    type(cell), intent(in) :: mesh
    real(gp), dimension(3) :: r
    !local variables
    integer :: i,j,k,ii,icurr,jcurr,kcurr
    real(gp) :: d,d2,dold
    real(gp), dimension(3) :: rt,ri,ci!,c_ortho,r_ortho

    if (mesh%orthorhombic) then
       do i=1,3
          r(i)=r_wrap(mesh%bc(i),mesh%hgrids(i)*mesh%ndims(i),&
               v(i),center(i))
       end do
    else
       dold=1.0d100 !huge_number
       icurr=0
       jcurr=0
       kcurr=0
       do ii=1,3
        if (mesh%bc(ii)==PERIODIC) then
          ri(ii)=mod(v(ii),mesh%ndims(ii)*mesh%hgrids(ii))
          ci(ii)=mod(center(ii),mesh%ndims(ii)*mesh%hgrids(ii))
        else
          ri(ii)=v(ii)
          ci(ii)=center(ii)
        end if
       end do
       do i=-mesh%bc(1),mesh%bc(1)
        do j=-mesh%bc(2),mesh%bc(2)
         do k=-mesh%bc(3),mesh%bc(3)
            rt(1)=ri(1)+real(i,gp)*mesh%ndims(1)*mesh%hgrids(1)
            rt(2)=ri(2)+real(j,gp)*mesh%ndims(2)*mesh%hgrids(2)
            rt(3)=ri(3)+real(k,gp)*mesh%ndims(3)*mesh%hgrids(3)
            d2=square_gd(mesh,rt-ci)
            d=sqrt(d2)
            if (d.lt.dold) then
               dold=d
               icurr=i
               jcurr=j
               kcurr=k
            end if
         end do
        end do
       end do
       d=dold
       r(1)=ri(1)+real(icurr,gp)*mesh%ndims(1)*mesh%hgrids(1) - ci(1)
       r(2)=ri(2)+real(jcurr,gp)*mesh%ndims(2)*mesh%hgrids(2) - ci(2)
       r(3)=ri(3)+real(kcurr,gp)*mesh%ndims(3)*mesh%hgrids(3) - ci(3)
    end if

  end function closest_r

  !> Calculates the square of the vector r in the cell defined by mesh
  !! Takes into account the non-orthorhombicity of the box
  !! with the controvariant metric (mesh%gu)
  pure function square(mesh,v)
    implicit none
    !> array of coordinate in the mesh reference frame
    real(gp), dimension(3), intent(in) :: v
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: square

    if (mesh%orthorhombic) then
       square=v(1)**2+v(2)**2+v(3)**2
    else
       square=dotp_gu(mesh,v,v)
    end if

  end function square

  function square_add(mesh,v_add) result(square)
    implicit none
    !> array of coordinate in the mesh reference frame
    real(gp) :: v_add
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: square

    if (mesh%orthorhombic) then
       call dotp_external_ortho(v_add,v_add,square)
    else
       call dotp_external_nonortho(mesh%gu,v_add,v_add,square)
    end if

  end function square_add

  !> Calculates the square of the vector r in the cell defined by mesh
  !! Takes into account the non-orthorhombicity of the box
  !! with the covariant metric (mesh%gd)
  pure function square_gd(mesh,v)
    implicit none
    !> array of coordinate in the mesh reference frame
    real(gp), dimension(3), intent(in) :: v
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: square_gd

    if (mesh%orthorhombic) then
       square_gd=v(1)**2+v(2)**2+v(3)**2
    else
       square_gd=dotp_gd(mesh,v,v)
    end if

  end function square_gd

  function square_gd_add(mesh,v_add) result(square)
    implicit none
    !> array of coordinate in the mesh reference frame
    real(gp) :: v_add
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: square

    if (mesh%orthorhombic) then
       call dotp_external_ortho(v_add,v_add,square)
    else
       call dotp_external_nonortho(mesh%gd,v_add,v_add,square)
    end if

  end function square_gd_add

  pure function dotp_gu(mesh,v1,v2)
    implicit none
    real(gp), dimension(3), intent(in) :: v1,v2
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: dotp_gu
    !local variables
    integer :: i,j

    if (mesh%orthorhombic) then
       dotp_gu=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
    else
       dotp_gu=0.0_gp
       do i=1,3
          do j=1,3
             dotp_gu=dotp_gu+mesh%gu(i,j)*v1(i)*v2(j)
          end do
       end do
    end if

  end function dotp_gu

  function dotp_gu_add2(mesh,v1,v2_add) result(dotp)
    implicit none
    real(gp), dimension(3), intent(in) :: v1
    real(gp) :: v2_add !<intent in, cannot be declared as such
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: dotp

    if (mesh%orthorhombic) then
       call dotp_external_ortho(v1,v2_add,dotp)
    else
       call dotp_external_nonortho(mesh%gu,v1,v2_add,dotp)
    end if

  end function dotp_gu_add2

  function dotp_gu_add1(mesh,v1_add,v2) result(dotp)
    implicit none
    real(gp), dimension(3), intent(in) :: v2
    real(gp) :: v1_add !<intent in, cannot be declared as such
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: dotp

    if (mesh%orthorhombic) then
       call dotp_external_ortho(v1_add,v2,dotp)
    else
       call dotp_external_nonortho(mesh%gu,v1_add,v2,dotp)
    end if

  end function dotp_gu_add1

  pure function dotp_gd(mesh,v1,v2)
    implicit none
    real(gp), dimension(3), intent(in) :: v1,v2
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: dotp_gd
    !local variables
    integer :: i,j

    if (mesh%orthorhombic) then
       dotp_gd=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
    else
       dotp_gd=0.0_gp
       do i=1,3
          do j=1,3
             dotp_gd=dotp_gd+mesh%gd(i,j)*v1(i)*v2(j)
          end do
       end do
    end if

  end function dotp_gd

  function dotp_gd_add2(mesh,v1,v2_add) result(dotp)
    implicit none
    real(gp), dimension(3), intent(in) :: v1
    real(gp) :: v2_add !<intent in, cannot be declared as such
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: dotp

    if (mesh%orthorhombic) then
       call dotp_external_ortho(v1,v2_add,dotp)
    else
       call dotp_external_nonortho(mesh%gd,v1,v2_add,dotp)
    end if

  end function dotp_gd_add2

  function dotp_gd_add1(mesh,v1_add,v2) result(dotp)
    implicit none
    real(gp), dimension(3), intent(in) :: v2
    real(gp) :: v1_add !<intent in, cannot be declared as such
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: dotp

    if (mesh%orthorhombic) then
       call dotp_external_ortho(v1_add,v2,dotp)
    else
       call dotp_external_nonortho(mesh%gd,v1_add,v2,dotp)
    end if

  end function dotp_gd_add1

end module box

subroutine dotp_external_ortho(v1,v2,dotp)  
  use f_precisions, only: gp=>f_double
  implicit none
  real(gp), dimension(3), intent(in) :: v1,v2
  real(gp), intent(out) :: dotp

  dotp=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
end subroutine dotp_external_ortho

subroutine dotp_external_nonortho(g,v1,v2,dotp)  
  use f_precisions, only: gp=>f_double
  implicit none
  real(gp), dimension(3,3), intent(in) :: g
  real(gp), dimension(3), intent(in) :: v1,v2
  real(gp), intent(out) :: dotp
  !local variables
  integer :: i,j

       dotp=0.0_gp
       do i=1,3
          do j=1,3
             dotp=dotp+g(i,j)*v1(i)*v2(j)
          end do
       end do
end subroutine dotp_external_nonortho
