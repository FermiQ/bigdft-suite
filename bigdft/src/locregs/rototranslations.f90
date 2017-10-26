!> @file
!!  Module to define the operator of rototranslation and related objects
!!  according to a transformation
!!
!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module rototranslations
  use module_defs, only: gp
  implicit none
  
  private

  !> Contains the rotation and translation (possibly deformation) which have to be applied to a given fragment
  type, public :: rototranslation
     !real(gp), dimension(3) ::  !< translation of fragment center
     real(gp), dimension(3) :: rot_center_new=0.0_gp !< positions of the centers
     real(gp), dimension(3) :: rot_center=0.0_gp !< center of rotation in original coordinates (might be fragment center or not)
     real(gp), dimension(3) :: rot_axis=[1.0_gp,0.0_gp,0.0_gp] !< unit rotation axis (should have modulus one)
     real(gp) :: theta=0.0_gp !< angle of rotation
     real(gp) :: Werror=0.0_gp !< Wahba error associated with transformation
     !> rotation matrix. Should be applied to the fragment
     !! functions to reformat
     real(gp), dimension(3,3) :: Rmat=reshape([1.0_gp,0.0_gp,0.0_gp,0.0_gp,1.0_gp,0.0_gp,&
          0.0_gp,0.0_gp,1.0_gp],[3,3])
  end type rototranslation

  public :: set_rototranslation
  
contains

  !> Defines an identity transformation
  pure function rototranslation_identity() result(ft)
    implicit none
    type(rototranslation) :: ft
    ft%rot_center    = 0.0_gp 
    ft%Werror        = 0.0_gp
    call set_no_translation(ft)
    call set_no_rotation(ft)
  end function rototranslation_identity

  pure subroutine set_translation(ft,src,dest)
    implicit none
    real(gp), dimension(3), intent(in) :: src,dest
    type(rototranslation), intent(inout) :: ft

    ft%rot_center=src
    ft%rot_center_new=dest
  end subroutine set_translation

  pure subroutine set_no_translation(ft)
    implicit none
    real(gp), dimension(3) :: src,dest
    type(rototranslation), intent(inout) :: ft

    ft%rot_center_new=ft%rot_center
  end subroutine set_no_translation

  pure subroutine set_rotation(ft,R)
    implicit none
    real(gp), dimension(3,3), intent(in) :: R
    type(rototranslation), intent(inout) :: ft

    ft%Rmat=R
    !find the angle from R matrix
    ft%theta=theta_from_r(transpose(R))
    !find rot_axis
    ft%rot_axis=axis_from_r(transpose(R))
  end subroutine set_rotation

  pure subroutine set_no_rotation(ft)
    implicit none
    type(rototranslation), intent(inout) :: ft
    !local variables
    real(gp), dimension(3,3) :: eye
    
    eye          = 0.0_gp
    eye(1,1)     = 1.0_gp
    eye(2,2)     = 1.0_gp
    eye(3,3)     = 1.0_gp

    call set_rotation(ft,eye)
    
  end subroutine set_no_rotation

  subroutine find_and_set_rotation(ft,nat,src,dest)
    implicit none
    integer, intent(in) :: nat !< fragment size
    real(gp), dimension(3,nat), intent(in) :: src,dest !<coordinates measured wrt rot_center
    type(rototranslation), intent(inout) :: ft
    !local variables
    real(gp) :: J
    real(gp), dimension(3,3) :: R
    
    call find_rotation(nat,src,dest,R,J)
    ft%Werror=J
    call set_rotation(ft,R)
  end subroutine find_and_set_rotation

  subroutine set_rototranslation(ft,nat,src,dest)
    use dynamic_memory
    implicit none
    integer, intent(in) :: nat !< fragment size
    real(gp), dimension(3,nat), intent(in) :: src,dest
    type(rototranslation), intent(out) :: ft
    !local variables
    integer :: iat
    real(gp), dimension(3) :: src_center,dest_center
    real(gp), dimension(:,:), allocatable :: src_0,dest_0

    ft=rototranslation_identity()
    
    src_center=frag_center(nat,src)
    dest_center=frag_center(nat,dest)

    call set_translation(ft,src=src_center,dest=dest_center)

    src_0=f_malloc(src=src,id='src_0')
    dest_0=f_malloc(src=dest,id='dest_0')
    do iat=1,nat
       src_0(:,iat)=src_0(:,iat)-src_center
       dest_0(:,iat)=dest_0(:,iat)-dest_center
    end do

    call find_and_set_rotation(ft,nat,src=src_0,dest=dest_0)
    
    call f_free(src_0)
    call f_free(dest_0)
  end subroutine set_rototranslation
  
  !> Express the coordinates of a vector into a rotated reference frame
  pure function rotate_vector(newz,theta,vec) result(vecn)
    !use module_base
    implicit none
    real(gp), intent(in) :: theta
    real(gp), dimension(3), intent(in) :: newz,vec
    real(gp), dimension(3) :: vecn
    !local variables
    real(gp) :: sint,cost,onemc,x,y,z

    !save recalculation
    sint=sin(theta)
    cost=cos(theta)
    onemc=1.0_gp-cost
    x=vec(1)
    y=vec(2)
    z=vec(3)

    vecn(1)=x*(cost + onemc*newz(1)**2) + y*(onemc*newz(1)*newz(2) - sint*newz(3)) &
         + z*(sint*newz(2) + onemc*newz(1)*newz(3))
    vecn(2)=y*(cost + onemc*newz(2)**2) + x*(onemc*newz(1)*newz(2) + sint*newz(3)) &
         + z*(-(sint*newz(1)) + onemc*newz(2)*newz(3))
    vecn(3)=z*(cost + onemc*newz(3)**2) + x*(onemc*newz(1)*newz(3) - sint*newz(2)) &
         + y*(sint*newz(1) + onemc*newz(2)*newz(3))

  end function rotate_vector


  !pure function transform_fragment_basis(trans,basis) result(basis_new)
  !  implicit none
  !  type(rototranslation), intent(in) :: trans
  !  type(fragment_basis), intent(in) :: basis
  !  type(fragment_basis) :: basis_new

  !  basis_new=fragment_basis_null()

  !  ! minimal_orbitals_data should remain the same
  !! want to have defined new lzd already?  in which case not a pure function...

  !!   integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
  !!   integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
  !!   type(local_zone_descriptors) :: Lzd
  !!   type(minimal_orbitals_data) :: forbs
  !!   real(wp), dimension(:), pointer :: psi

  !end function transform_fragment_basis

  subroutine find_rotation(nat,rxyz_ref,rxyz_new,R,J)
    use module_base
    use yaml_output
    implicit none
    integer, intent(in) :: nat !< fragment size
    real(gp), dimension(3,nat), intent(in) :: rxyz_ref,rxyz_new !<coordinates measured wrt rot_center
    real(gp), intent(out) :: J !< Wahba cost function, i.e. error in transformation - now included in structure
    real(gp), dimension(3,3), intent(out) :: R

    !local variables
    integer, parameter :: lwork=7*3
    integer :: info,iat!,i_stat,i
    real(gp) :: dets, J0
    real(gp), dimension(3) :: SM_arr !< array of SVD and M array
    real(gp), dimension(lwork) :: work !< array of SVD and M array
    real(gp), dimension(3,nat) :: J_arr !< matrix for calculating Wahba's cost function
    real(gp), dimension(3,3) :: B_mat,R_mat,U_mat,VT_mat !<matrices of Wahba's problem
    !character(len=100) :: subname

    !subname='find_frag_trans'

    !find the error with the identity transformation
    J0=0.0_gp
    do iat=1,nat
       J0=J0+(rxyz_new(1,iat)-rxyz_ref(1,iat))**2+&
            (rxyz_new(2,iat)-rxyz_ref(2,iat))**2+&
            (rxyz_new(3,iat)-rxyz_ref(3,iat))**2
    end do
    !if this is already little, go ahead
    if (J0 < 1.e-6_gp) then
       J=J0
       R=0.0_gp
       R(1,1)=1.0_gp
       R(2,2)=1.0_gp
       R(3,3)=1.0_gp
       return
    end if
    
    B_mat=0.0_gp
    R_mat=0.0_gp

    !all positions are of weight one for the moment
    call dgemm('N','T',3,3,nat,1.0_gp,rxyz_new,3,rxyz_ref,3,0.0_gp,B_mat,3)

    !find matrix of svd
    call dgesvd('A','A',3,3,B_mat,3,SM_arr,U_mat,3,VT_mat,3,work,lwork,info)
    if (f_err_raise(info/=0,'Problem in DGESVD')) return

    !multiply last line of VT_mat by det(U)*det(V)
    dets=det_3x3(U_mat)*det_3x3(VT_mat)
    VT_mat(3,:)=VT_mat(3,:)*dets

    !find rotation matrix
    call dgemm('N','N',3,3,3,1.0_gp,U_mat,3,VT_mat,3,0.0_gp,R_mat,3)

    R=transpose(R_mat)
    !call set_rotation(frag_trans,transpose(R_mat))
!!$    !store rotation matrix in columns
!!$    frag_trans%Rmat=transpose(R_mat)
!!$
!!$    !find the angle from R matrix
!!$    frag_trans%theta=theta_from_r(R_mat)
!!$    !find rot_axis
!!$    frag_trans%rot_axis=axis_from_r(R_mat)

!!$    call yaml_map('Rmat found',frag_trans%Rmat)

    !print*,'Rmat:',frag_trans%theta
    !do i=1,3
    !   write(*,'(3(F12.6,2x))') R_mat(i,:)
    !end do

    !print*,'Rcalc:',frag_trans%theta
    !write(*,'(3(F12.6,2x))') cos(frag_trans%theta)+frag_trans%rot_axis(1)**2*(1.0_gp-cos(frag_trans%theta)),&
    !     frag_trans%rot_axis(1)*frag_trans%rot_axis(2)*(1.0_gp-cos(frag_trans%theta))&
    !     -frag_trans%rot_axis(3)*sin(frag_trans%theta),&
    !     frag_trans%rot_axis(1)*frag_trans%rot_axis(3)*(1.0_gp-cos(frag_trans%theta))&
    !     +frag_trans%rot_axis(2)*sin(frag_trans%theta)
    !write(*,'(3(F12.6,2x))') frag_trans%rot_axis(2)*frag_trans%rot_axis(1)*(1.0_gp-cos(frag_trans%theta))&
    !     +frag_trans%rot_axis(3)*sin(frag_trans%theta),&
    !     cos(frag_trans%theta)+frag_trans%rot_axis(2)**2*(1.0_gp-cos(frag_trans%theta)),&
    !     frag_trans%rot_axis(2)*frag_trans%rot_axis(3)*(1.0_gp-cos(frag_trans%theta))&
    !     -frag_trans%rot_axis(1)*sin(frag_trans%theta)
    !write(*,'(3(F12.6,2x))') frag_trans%rot_axis(3)*frag_trans%rot_axis(1)*(1.0_gp-cos(frag_trans%theta))&
    !     -frag_trans%rot_axis(2)*sin(frag_trans%theta),&
    !     frag_trans%rot_axis(3)*frag_trans%rot_axis(2)*(1.0_gp-cos(frag_trans%theta))&
    !     +frag_trans%rot_axis(1)*sin(frag_trans%theta),&
    !     cos(frag_trans%theta)+frag_trans%rot_axis(3)**2*(1.0_gp-cos(frag_trans%theta))

    !to be verified that the cost function of Wahba's problem is little
    J_arr=rxyz_new
    call dgemm('N','N',3,nat,3,-1.0_gp,R_mat,3,rxyz_ref,3,1.0_gp,J_arr,3)
    J=0.0_gp!frag_trans%Werror=0.0_gp
    do iat=1,nat
       !frag_trans%Werror=frag_trans%Werror+J_arr(1,iat)**2+J_arr(2,iat)**2+J_arr(3,iat)**2
       J=J+J_arr(1,iat)**2+J_arr(2,iat)**2+J_arr(3,iat)**2
    end do

!!$    !here yaml output
!!$    !make it optional whether to print the warning from here or leave it to external function
!!$    if (J>1.0e-3) then
!!$       write(*,'(a,2es18.8)') "Error, Wahba's cost function is too big",J,frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp)
!!$    end if

    !check the pertinence of the suggested rotation
    !if (abs(frag_trans%theta) > 60.d0*(4.0_gp*atan(1.d0)/180.0_gp)) print*,'frag_trans%theta=',frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp)
    !if  (f_err_raise(abs(frag_trans%theta) > 60.d0*(4.0_gp*atan(1.d0)/180.0_gp),'Angle frag_trans%theta not optimal (frag_trans%theta= '//&
    !      yaml_toa(frag_trans%theta)//' )')) return

!!$    !check if we could achieve same error without a rotation
!!$    ! - want to avoid unecessary rotation by e.g. 180 degrees
!!$    if (frag_trans%theta/=0.0d0) then
!!$       J0=0.0_gp
!!$       do iat=1,nat
!!$          J0=J0+(rxyz_new(1,iat)-rxyz_ref(1,iat))**2+(rxyz_new(2,iat)-rxyz_ref(2,iat))**2+(rxyz_new(3,iat)-rxyz_ref(3,iat))**2
!!$       end do
!!$
!!$       !replace with no rotation
!!$       if (J0<J .or. J0-J<1e-6) then
!!$          !write(*,'(a,6(es12.4,2x))') 'replacing suggested transformation with zero transformation ',&
!!$          !     J0,J,frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp),frag_trans%rot_axis
!!$          !else 
!!$          !   write(*,'(a,6(es12.4,2x))') 'NOT replacing suggested transformation with zero transformation ',&
!!$          !        J0,J,frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp),frag_trans%rot_axis
!!$       end if
!!$    end if

  end subroutine find_rotation


  pure function theta_from_r(R_mat) result(theta)
    implicit none
    real(gp), dimension(3,3), intent(in) :: R_mat

    real(gp) :: theta

    !local variables
    real(gp) :: tracem1

    tracem1=R_mat(1,1)+R_mat(2,2)+R_mat(3,3)-1.0_gp

    if (abs(tracem1) - 2.0_gp > -1.e-14_gp) then
       if (tracem1 > 0.0_gp) theta = 0.0_gp
       if (tracem1 <= 0.0_gp) theta = 180.0_gp*(4.0_gp*atan(1.d0)/180.0_gp)
    else
       theta=acos(0.5_gp*tracem1)
    end if

  end function theta_from_r


  pure function axis_from_r(R_mat) result(rot_axis)
    implicit none
    real(gp), dimension(3,3), intent(in) :: R_mat
    real(gp), dimension(3) :: rot_axis

    !local variables
    real(gp) :: dnrm2, norm
    integer :: i

    rot_axis(1)=R_mat(3,2)-R_mat(2,3)    
    rot_axis(2)=R_mat(1,3)-R_mat(3,1)
    rot_axis(3)=R_mat(2,1)-R_mat(1,2)    

    !normalize it

    norm=rot_axis(1)**2+rot_axis(2)**2+rot_axis(3)**2
    norm=sqrt(norm)
    !print*,'axis_from_r',norm,rot_axis
!!$    if (norm>=1.e-5_gp) then
!!$       !call dscal(3,1.0_gp/norm,rot_axis,1)
!!$       !print*,'axis_from_r2',norm,rot_axis
!!$    else
    if (norm < 1.e-5_gp) then
       ! squares of rot_axis are diag 0.5(R+I), signs as before
       !this is not good as it gives a array of modulus bigger than one
       !rot_axis(1:2)=0.0_gp
       !rot_axis(3)=1.0_gp
       do i=1,3
          if (R_mat(i,i)<-1.0_gp.and.R_mat(i,i)>-1.0_gp-1.0e-5_gp) then
             rot_axis(i)=0.0_gp
          else if (R_mat(i,i)<=-1.0_gp-1.0e-5_gp) then
             !stop 'Problem in assigning axis from Rmat'
          else
             rot_axis(i)=sign(dsqrt(0.5_gp*(R_mat(i,i)+1.0_gp)),rot_axis(i))
          end if
       end do

       !print*,'axis_from_r3a',rot_axis,R_mat(1,1),R_mat(2,2),R_mat(3,3)
       !norm=dnrm2(3,rot_axis,1)
       norm=rot_axis(1)**2+rot_axis(2)**2+rot_axis(3)**2
       norm=sqrt(norm)

       !call dscal(3,1.0_gp/norm,rot_axis,1)
       !print*,'axis_from_r3',norm,rot_axis
    end if
    rot_axis=rot_axis/norm
  end function axis_from_r

  pure function frag_center(nat,rxyz) result(cen)
    implicit none
    integer, intent(in) :: nat 
    real(gp), dimension(3,nat), intent(in) :: rxyz
    real(gp), dimension(3) :: cen
    !local variables
    integer :: iat,i

    cen=0.0_gp
    if (nat > 0) then
       do iat=1,nat
          do i=1,3
             cen(i)=cen(i)+rxyz(i,iat)
          end do
       end do
       cen=cen/real(nat,gp)
    end if

  end function frag_center

  
end module rototranslations
