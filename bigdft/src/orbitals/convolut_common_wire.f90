!> @file
!!  Common convolutions
!! @author
!!    Copyright (C) 2018-2018 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!>  A analysis wavelet transformation where the size of the data is forced to shrink
!!  The input array y is overwritten
subroutine analyse_wire_self(n1,n2,n3,y,x)
  implicit none
  integer,intent(in)::n1,n2,n3
  real(kind=8),dimension((2*n1+16)*(2*n2+16)*(2*n3+2))::x,y
  integer nt

  ! i1,I2,i3 -> I2,i3,i1
  nt=(2*n2+16)*(2*n3+2)
  call  ana_rot_shrink(n1,nt,y,x)
  ! I2,i3,i1 -> i3,i1,i2
  nt=(2*n3+2)*(2*n1+16)
  call  ana_rot_shrink(n2,nt,x,y)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_per(n3,nt,y,x)

END SUBROUTINE analyse_wire_self

!>   An analysis wavelet transformation where the size of the data is forced to shrink
!!   The input array y is not overwritten
subroutine analyse_wire(n1,n2,n3,ww,y,x)
  implicit none
  !Arguments
  integer, intent(in) :: n1,n2,n3
  real(kind=8) :: x(0:n1,2,0:n2,2,0:n3,2)
  real(kind=8) :: y (-7:2*n1+8,-7:2*n2+8,0:2*n3+1)
  real(kind=8) :: ww(-7:2*n1+8,-7:2*n2+8,0:2*n3+1)
  !Local variables
  integer :: nt

  ! i1,I2,i3 -> I2,i3,i1
  nt=(2*n2+16)*(2*n3+2)
  call  ana_rot_shrink(n1,nt,y,ww)
  ! I2,i3,i1 -> i3,i1,i2
  nt=(2*n3+2)*(2*n1+16)
  call  ana_rot_shrink(n2,nt,ww,y)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_per(n3,nt,y,x)

END SUBROUTINE analyse_wire

!>   A synthesis wavelet transformation where the size of the data is allowed to grow
!!   The input array x is overwritten
subroutine synthese_wire_self(n1,n2,n3,x,y)
  implicit none
  integer,intent(in)::n1,n2,n3
  real(kind=8),dimension((2*n1+16)*(2*n2+16)*(2*n3+2))::x,y
  integer nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_grow(n1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,I2
  nt=(2*n3+2)*(2*n1+16)
  call  syn_rot_grow(n2,nt,y,x)
  ! i3,i1,I2  -> i1,I2,i3
  nt=(2*n1+16)*(2*n2+16)
  call  syn_rot_per(n3,nt,x,y)

END SUBROUTINE synthese_wire_self

!>   A synthesis wavelet transformation where the size of the data is allowed to grow
!!   The input array x is not overwritten
subroutine synthese_wire(n1,n2,n3,ww,x,y)
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: x(0:n1,2,0:n2,2,0:n3,2)
  real(kind=8) :: y (-7:2*n1+8,-7:2*n2+8,0:2*n3+1)
  real(kind=8) :: ww(-7:2*n1+8,-7:2*n2+8,0:2*n3+1)

  ! i1,i2,i3 -> i2,i3,i1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_grow(n1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,I2
  nt=(2*n3+2)*(2*n1+16)
  call  syn_rot_grow(n2,nt,y,ww)
  ! i3,i1,I2  -> i1,I2,i3
  nt=(2*n1+16)*(2*n2+16)
  call  syn_rot_per(n3,nt,ww,y)

END SUBROUTINE synthese_wire

!>   Applies the magic filter matrix in wire BC ( no transposition)
!!   The input array x is overwritten
!!   this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_wire_self(n1,n2,n3,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(-7:n1+8,-7:n2+8,0:n3), intent(inout) :: x
  real(wp), dimension(-7:n1+8,-7:n2+8,0:n3), intent(inout) :: y
  !local variables
  integer :: ndat

  !  (i1,i2*i3) -> (i2*i3,i1)
  ndat=(n2+1)*(n3+1)
  call convrot_grow(n1,ndat,x,y)
  !  (i2,i3*i1) -> (i3*i1,I2)
  ndat=(n3+1)*(n1+16)
  call convrot_grow(n2,ndat,y,x)
  !  (i3,i1*I2) -> (i1*I2,i3)
  ndat=(n1+16)*(n2+16)
  call convrot_n_per(n3,ndat,x,y)

END SUBROUTINE convolut_magic_n_wire_self


!>   Applies the magic filter matrix in periodic BC ( no transposition)
!!   The input array x is not overwritten
!!   this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_wire(n1,n2,n3,x,y,ww)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(-7:n1+8,-7:n2+8,0:n3), intent(inout) :: y
  !local variables
  integer :: ndat
  real(wp), dimension(-7:n1+8,-7:n2+8,0:n3):: ww ! work array

  !  (i1,i2*i3) -> (i2*i3,i1)
  ndat=(n2+1)*(n3+1)
  call convrot_grow(n1,ndat,x,y)
  !  (i2,i3*i1) -> (i3*i1,I2)
  ndat=(n3+1)*(n1+16)
  call convrot_grow(n2,ndat,y,ww)
  !  (i3,i1*I2) -> (i1*I2,i3)
  ndat=(n1+16)*(n2+16)
  call convrot_n_per(n3,ndat,ww,y)

END SUBROUTINE convolut_magic_n_wire


!>   Applies the magic filter matrix transposed in periodic BC 
!!   The input array x is overwritten
!!   this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_t_wire_self(n1,n2,n3,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(-7:n1+8,-7:n2+8,0:n3), intent(inout) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer :: ndat

  !  (i1,I2*i3) -> (I2*i3,i1)
  ndat=(n2+16)*(n3+1)
  call convrot_shrink(n1,ndat,x,y)
  !  (I2,i3*i1) -> (i3*i1,i2)
  ndat=(n3+1)*(n1+16)
  call convrot_shrink(n2,ndat,y,x)
  !  (i3,i1*i2) -> (i1*i2,i3)
  ndat=(n1+1)*(n2+1)
  call convrot_t_per(n3,ndat,x,y)


END SUBROUTINE convolut_magic_t_wire_self


!> Applies the modified kinetic energy operator onto x to get y. 
!! Computes the kinetic energy too.
!! Works for the slab BC.
!! Modified kinetic energy operator:
!! A=-1/2 exp(Ikr) Delta exp(-Ikr)
!! where k=(k1,k2,k3); r=(x,y,z)
subroutine convolut_kinetic_wire_T_k(n1,n2,n3,hgrid,x,y,ener,k1,k2,k3)

  use module_defs, only: wp,gp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::k1,k2,k3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(inout) :: y
  real(wp), intent(out) :: ener
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt1,tt2,c,enerp
  real(wp), dimension(3) :: scale,scale1
  real(wp), dimension(2,lowfil:lupfil,3) :: fil  

  scale (:)=real(-.5_gp/hgrid(:)**2,wp)
  scale1(1)=real(k1/hgrid(1),wp)
  scale1(2)=real(k2/hgrid(2),wp)
  scale1(3)=real(k3/hgrid(3),wp)
  c=.5_wp*(k1*k1+k2*k2+k3*k3)

  ! second derivative filters for Daubechies 16
  fil(1,0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(1,2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(1,3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(1,4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(1,5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(1,6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(1,7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(1,8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(1,9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(1,10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(1,11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(1,12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(1,13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(1,14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do i=1,14
     fil(1,-i,:)=fil(1,i,:)
  enddo

  fil(2,0,:)= 0._wp
  fil(2,1,:)= 0.8834460460908270942785856e0_wp*scale1(:)
  fil(2,2,:)= -0.3032593514765938346887962e0_wp*scale1(:)
  fil(2,3,:)= 0.1063640682894442760934532e0_wp*scale1(:)
  fil(2,4,:)= -0.03129014783948023634381564e0_wp*scale1(:)
  fil(2,5,:)= 0.006958379116450707495020408e0_wp*scale1(:)
  fil(2,6,:)= -0.001031530213375445369097965e0_wp*scale1(:)
  fil(2,7,:)= 0.00007667706908380351933901775e0_wp*scale1(:)
  fil(2,8,:)= 2.451992111053665419191564e-7_wp*scale1(:)
  fil(2,9,:)= 3.993810456408053712133667e-8_wp*scale1(:)
  fil(2,10,:)=-7.207948238588481597101904e-8_wp*scale1(:)
  fil(2,11,:)=-9.697184925637300947553069e-10_wp*scale1(:)
  fil(2,12,:)=-7.252206916665149851135592e-13_wp*scale1(:)
  fil(2,13,:)=1.240078536096648534547439e-14_wp*scale1(:)
  fil(2,14,:)=-1.585464751677102510097179e-19_wp*scale1(:)

  do i=1,14
     fil(2,-i,:)=-fil(2,i,:)
  enddo

  ener=0._wp
!$omp parallel default (private) shared(x,y,ener,fil,c,n1,n2,n3)
  enerp=0._wp

!$omp do 

  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,n1
           tt1=x(1,i1,i2,i3)*c
           tt2=x(2,i1,i2,i3)*c
           !do l=lowfil,lupfil
           !   j=modulo(i1+l,n1+1)
           do l=max(lowfil,-i1),min(lupfil,n1-i1)
           j=i1+l 
              tt1=tt1+x(1,j,i2,i3)*fil(1,l,1)-x(2,j,i2,i3)*fil(2,l,1)
              tt2=tt2+x(2,j,i2,i3)*fil(1,l,1)+x(1,j,i2,i3)*fil(2,l,1)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2
         enerp=enerp+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt1=0._wp
           tt2=0._wp
           !do l=lowfil,lupfil
           !   j=modulo(i2+l,n2+1)
           do l=max(lowfil,-i2),min(lupfil,n2-i2)
           j=i2+l 
              tt1=tt1+x(1,i1,j,i3)*fil(1,l,2)-x(2,i1,j,i3)*fil(2,l,2)
              tt2=tt2+x(2,i1,j,i3)*fil(1,l,2)+x(1,i1,j,i3)*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2
         enerp=enerp+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
        enddo
     enddo
     
  enddo
!$omp enddo

!$omp do 
  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt1=0._wp
           tt2=0._wp
           do l=lowfil,lupfil
              j=modulo(i3+l,n3+1)
              tt1=tt1+x(1,i1,i2,j)*fil(1,l,3)-x(2,i1,i2,j)*fil(2,l,3)
              tt2=tt2+x(2,i1,i2,j)*fil(1,l,3)+x(1,i1,i2,j)*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2
         enerp=enerp+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
        enddo
     enddo
  enddo
!$omp enddo
!$omp critical
ener=ener+enerp
!$omp end critical
!  ener=ener*.5_wp
!$omp end parallel  
END SUBROUTINE convolut_kinetic_wire_T_k


subroutine convolut_kinetic_wire_T(n1,n2,n3,hgrid,x,y,ekin)
!   applies the kinetic energy operator onto x to get y. Works for surface BC
!   y:=y-1/2Delta x
  use module_defs, only: wp,gp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
  real(wp),intent(out)::ekin
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  
  integer :: mod_arr1(lowfil:n1+lupfil)   
  integer :: mod_arr3(lowfil:n3+lupfil)   
  !integer :: ncount0,ncount1,ncount_max,ncount_rate

  call fill_mod_arr(mod_arr3,lowfil,n3+lupfil,n3+1)
   
  scale(:)=real(-.5_gp/hgrid(:)**2,wp)

  ! second derivative filters for Daubechies 16
  fil(0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do i=1,14
     fil(-i,:)=fil(i,:)
  enddo
  ekin=0.0_wp

  !call system_clock(ncount0,ncount_rate,ncount_max)

  !call conv_kin_x(x,y,(n2+1)*(n3+1))   
  call conv_kin_x_new_wire(n1,x,y,(n2+1)*(n3+1),lowfil,lupfil,fil,ekin)!,mod_arr1)

  !call conv_kin_y
  call conv_kin_y_new(n1,n2,n3,x,y,lowfil,lupfil,fil,ekin)

  !call conv_kin_z(x,y,(n1+1)*(n2+1))
  call conv_kin_z_new(n3,x,y,(n1+1)*(n2+1),lowfil,lupfil,fil,mod_arr3,ekin)
  
  !call system_clock(ncount1,ncount_rate,ncount_max)
  !write(*,*) 'TIMING:convolut_kinetic_slab_T',real(ncount1-ncount0)/real(ncount_rate)
  
END SUBROUTINE convolut_kinetic_wire_T


!>   Applies the kinetic energy operator onto x to get y. Works for wire BC
subroutine convolut_kinetic_wire_sdc(n1,n2,n3,x,y,cprecr,modul3,a,b,c,e)
  use module_defs, only: wp,gp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in) :: cprecr
  real(wp), dimension(8,0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(8,0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,l,j
  real(wp) :: tt111,tt112,tt121,tt122,tt211,tt212,tt221,tt222
  integer,intent(in) :: modul3(lowfil:n3+lupfil)
  real(gp),intent(in) :: a(lowfil:lupfil,3)
  real(gp),intent(in) :: b(lowfil:lupfil,3)
  real(gp),intent(in) :: c(lowfil:lupfil,3)
  real(gp),intent(in) :: e(lowfil:lupfil,3)
  
!dee
!$omp parallel default(private) &
!$omp shared (n1,n2,n3,x,y,cprecr,a,b,c,e,modul3) 


!$omp do  
  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,n1
           tt111=x(1,i1,i2,i3)*cprecr
           tt211=x(2,i1,i2,i3)*cprecr
           tt121=x(3,i1,i2,i3)*cprecr
           tt221=x(4,i1,i2,i3)*cprecr
           tt112=x(5,i1,i2,i3)*cprecr
           tt212=x(6,i1,i2,i3)*cprecr
           tt122=x(7,i1,i2,i3)*cprecr
           tt222=x(8,i1,i2,i3)*cprecr
           !do l=lowfil,lupfil
           !   j=modul1(i1+l)
           do l=max(lowfil,-i1),min(lupfil,n1-i1) ! so that 0 =< j =< n2
              j=i1+l
              tt111=tt111+x(1,j,i2,i3)*a(l,1)+x(2,j,i2,i3)*b(l,1)
              tt211=tt211+x(2,j,i2,i3)*e(l,1)+x(1,j,i2,i3)*c(l,1)
              tt121=tt121+x(3,j,i2,i3)*a(l,1)+x(4,j,i2,i3)*b(l,1)
              tt221=tt221+x(4,j,i2,i3)*e(l,1)+x(3,j,i2,i3)*c(l,1)
              tt112=tt112+x(5,j,i2,i3)*a(l,1)+x(6,j,i2,i3)*b(l,1)
              tt212=tt212+x(6,j,i2,i3)*e(l,1)+x(5,j,i2,i3)*c(l,1)
              tt122=tt122+x(7,j,i2,i3)*a(l,1)+x(8,j,i2,i3)*b(l,1)
              tt222=tt222+x(8,j,i2,i3)*e(l,1)+x(7,j,i2,i3)*c(l,1)
           enddo
           y(1,i1,i2,i3)=tt111
           y(2,i1,i2,i3)=tt211
           y(3,i1,i2,i3)=tt121
           y(4,i1,i2,i3)=tt221
           y(5,i1,i2,i3)=tt112
           y(6,i1,i2,i3)=tt212
           y(7,i1,i2,i3)=tt122
           y(8,i1,i2,i3)=tt222
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt111=0.e0_wp
           tt211=0.e0_wp
           tt121=0.e0_wp
           tt221=0.e0_wp
           tt112=0.e0_wp
           tt212=0.e0_wp
           tt122=0.e0_wp
           tt222=0.e0_wp
           do l=max(lowfil,-i2),min(lupfil,n2-i2) ! so that 0 =< j =< n2
              j=i2+l
              tt111=tt111+x(1,i1,j,i3)*a(l,2)+x(3,i1,j,i3)*b(l,2)
              tt211=tt211+x(2,i1,j,i3)*a(l,2)+x(4,i1,j,i3)*b(l,2)
              tt121=tt121+x(3,i1,j,i3)*e(l,2)+x(1,i1,j,i3)*c(l,2)
              tt221=tt221+x(4,i1,j,i3)*e(l,2)+x(2,i1,j,i3)*c(l,2)
              tt112=tt112+x(5,i1,j,i3)*a(l,2)+x(7,i1,j,i3)*b(l,2)
              tt212=tt212+x(6,i1,j,i3)*a(l,2)+x(8,i1,j,i3)*b(l,2)
              tt122=tt122+x(7,i1,j,i3)*e(l,2)+x(5,i1,j,i3)*c(l,2)
              tt222=tt222+x(8,i1,j,i3)*e(l,2)+x(6,i1,j,i3)*c(l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt111
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt211
           y(3,i1,i2,i3)=y(3,i1,i2,i3)+tt121
           y(4,i1,i2,i3)=y(4,i1,i2,i3)+tt221
           y(5,i1,i2,i3)=y(5,i1,i2,i3)+tt112
           y(6,i1,i2,i3)=y(6,i1,i2,i3)+tt212
           y(7,i1,i2,i3)=y(7,i1,i2,i3)+tt122
           y(8,i1,i2,i3)=y(8,i1,i2,i3)+tt222
        enddo
     enddo
     
  enddo
!$omp enddo

!$omp do  
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt111=0.e0_wp
           tt211=0.e0_wp
           tt121=0.e0_wp
           tt221=0.e0_wp
           tt112=0.e0_wp
           tt212=0.e0_wp
           tt122=0.e0_wp
           tt222=0.e0_wp
           do l=lowfil,lupfil
              j=modul3(i3+l)
              tt111=tt111+x(1,i1,i2,j)*a(l,3)+x(5,i1,i2,j)*b(l,3)
              tt211=tt211+x(2,i1,i2,j)*a(l,3)+x(6,i1,i2,j)*b(l,3)
              tt121=tt121+x(3,i1,i2,j)*a(l,3)+x(7,i1,i2,j)*b(l,3)
              tt221=tt221+x(4,i1,i2,j)*a(l,3)+x(8,i1,i2,j)*b(l,3)
              tt112=tt112+x(5,i1,i2,j)*e(l,3)+x(1,i1,i2,j)*c(l,3)
              tt212=tt212+x(6,i1,i2,j)*e(l,3)+x(2,i1,i2,j)*c(l,3)
              tt122=tt122+x(7,i1,i2,j)*e(l,3)+x(3,i1,i2,j)*c(l,3)
              tt222=tt222+x(8,i1,i2,j)*e(l,3)+x(4,i1,i2,j)*c(l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt111
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt211
           y(3,i1,i2,i3)=y(3,i1,i2,i3)+tt121
           y(4,i1,i2,i3)=y(4,i1,i2,i3)+tt221
           y(5,i1,i2,i3)=y(5,i1,i2,i3)+tt112
           y(6,i1,i2,i3)=y(6,i1,i2,i3)+tt212
           y(7,i1,i2,i3)=y(7,i1,i2,i3)+tt122
           y(8,i1,i2,i3)=y(8,i1,i2,i3)+tt222
        enddo
     enddo
  enddo
!$omp enddo 
!$omp end parallel

END SUBROUTINE convolut_kinetic_wire_sdc


!> Applies the modified kinetic energy operator onto x to get y. 
!! Works for the wire BC.
!! Modified kinetic energy operator:
!! A=-1/2 exp(Ikr) Delta exp(-Ikr)+C
!! where k=(k1,k2,k3); r=(x,y,z)
subroutine convolut_kinetic_wire_c_k(n1,n2,n3,hgrid,x,y,c_in,k1,k2,k3)

  use module_defs, only: wp,gp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::c_in,k1,k2,k3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt1,tt2,c
  real(wp), dimension(3) :: scale,scale1
  real(wp), dimension(2,lowfil:lupfil,3) :: fil  

  scale (:)=real(-.5_gp/hgrid(:)**2,wp)
  scale1(1)=real(k1/hgrid(1),wp)
  scale1(2)=real(k2/hgrid(2),wp)
  scale1(3)=real(k3/hgrid(3),wp)
  c=c_in+.5_wp*(k1*k1+k2*k2+k3*k3)

  ! second derivative filters for Daubechies 16
  fil(1,0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(1,2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(1,3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(1,4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(1,5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(1,6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(1,7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(1,8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(1,9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(1,10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(1,11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(1,12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(1,13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(1,14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do i=1,14
     fil(1,-i,:)=fil(1,i,:)
  enddo

  fil(2,0,:)= 0._wp
  fil(2,1,:)= 0.8834460460908270942785856e0_wp*scale1(:)
  fil(2,2,:)= -0.3032593514765938346887962e0_wp*scale1(:)
  fil(2,3,:)= 0.1063640682894442760934532e0_wp*scale1(:)
  fil(2,4,:)= -0.03129014783948023634381564e0_wp*scale1(:)
  fil(2,5,:)= 0.006958379116450707495020408e0_wp*scale1(:)
  fil(2,6,:)= -0.001031530213375445369097965e0_wp*scale1(:)
  fil(2,7,:)= 0.00007667706908380351933901775e0_wp*scale1(:)
  fil(2,8,:)= 2.451992111053665419191564e-7_wp*scale1(:)
  fil(2,9,:)= 3.993810456408053712133667e-8_wp*scale1(:)
  fil(2,10,:)=-7.207948238588481597101904e-8_wp*scale1(:)
  fil(2,11,:)=-9.697184925637300947553069e-10_wp*scale1(:)
  fil(2,12,:)=-7.252206916665149851135592e-13_wp*scale1(:)
  fil(2,13,:)=1.240078536096648534547439e-14_wp*scale1(:)
  fil(2,14,:)=-1.585464751677102510097179e-19_wp*scale1(:)

  do i=1,14
     fil(2,-i,:)=-fil(2,i,:)
  enddo
!$omp parallel default (private) shared(n1,n2,n3,x,y,c,fil)
!$omp do

  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,n1
           tt1=x(1,i1,i2,i3)*c
           tt2=x(2,i1,i2,i3)*c
!           do l=lowfil,lupfil
!              j=modulo(i1+l,n1+1)
           do l=max(lowfil,-i1),min(lupfil,n1-i1)
           j=i1+l 
              tt1=tt1+x(1,j,i2,i3)*fil(1,l,1)-x(2,j,i2,i3)*fil(2,l,1)
              tt2=tt2+x(2,j,i2,i3)*fil(1,l,1)+x(1,j,i2,i3)*fil(2,l,1)
           enddo
           y(1,i1,i2,i3)=tt1
           y(2,i1,i2,i3)=tt2
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt1=0._wp
           tt2=0._wp
!           do l=lowfil,lupfil
!              j=modulo(i2+l,n2+1)
           do l=max(lowfil,-i2),min(lupfil,n2-i2)
           j=i2+l 
              tt1=tt1+x(1,i1,j,i3)*fil(1,l,2)-x(2,i1,j,i3)*fil(2,l,2)
              tt2=tt2+x(2,i1,j,i3)*fil(1,l,2)+x(1,i1,j,i3)*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2
        enddo
     enddo
     
  enddo
!$omp enddo
!$omp do

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt1=0._wp
           tt2=0._wp
           do l=lowfil,lupfil
              j=modulo(i3+l,n3+1)
              tt1=tt1+x(1,i1,i2,j)*fil(1,l,3)-x(2,i1,i2,j)*fil(2,l,3)
              tt2=tt2+x(2,i1,i2,j)*fil(1,l,3)+x(1,i1,i2,j)*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2
        enddo
     enddo
  enddo
!$omp enddo
!$omp end parallel  
END SUBROUTINE convolut_kinetic_wire_c_k

subroutine prepare_sdc_wire(n3,modul3,a,b,c,e,hx,hy,hz)
  use module_defs, only: wp,gp
  implicit none
  integer,intent(in)::n3
  real(gp),intent(in)::hx,hy,hz
  
  integer, parameter :: lowfil=-14,lupfil=14
  
  integer,intent(out)::modul3(lowfil:n3+lupfil)
  real(gp),intent(out)::a(lowfil:lupfil,3)
  real(gp),intent(out)::b(lowfil:lupfil,3)
  real(gp),intent(out)::c(lowfil:lupfil,3)
  real(gp),intent(out)::e(lowfil:lupfil,3)
  
  real(gp)::hgrid(3)
  integer::i
  real(gp)::scale(3)

  call fill_mod_arr(modul3,lowfil,n3+lupfil,n3+1)
  
  hgrid(1)=hx
  hgrid(2)=hy
  hgrid(3)=hz

  scale(:)=real(-.5_gp/hgrid(:)**2,wp)
  
  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0,:)=   -3.5536922899131901941296809374_wp*scale(:)
  a(1,:)=    2.2191465938911163898794546405_wp*scale(:)
  a(2,:)=   -0.6156141465570069496314853949_wp*scale(:)
  a(3,:)=    0.2371780582153805636239247476_wp*scale(:)
  a(4,:)=   -0.0822663999742123340987663521_wp*scale(:)
  a(5,:)=    0.02207029188482255523789911295638968409_wp*scale(:)
  a(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  a(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  a(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  a(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  a(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  a(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  a(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  a(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  a(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)
  do i=1,14
     a(-i,:)=a(i,:)
  enddo
  !  <phi|D^2|psi_i>
  c(-14,:)=     -3.869102413147656535541850057188e-18_wp*scale(:)
  c(-13,:)=      1.5130616560866154733900029272077362e-13_wp*scale(:)
  c(-12,:)=     -3.2264702314010525539061647271983988409e-11_wp*scale(:)
  c(-11,:)=     -5.96264938781402337319841002642e-9_wp*scale(:)
  c(-10,:)=     -2.1656830629214041470164889350342e-7_wp*scale(:)
  c(-9 ,:)=      8.7969704055286288323596890609625e-7_wp*scale(:)
  c(-8 ,:)=     -0.00001133456724516819987751818232711775_wp*scale(:)
  c(-7 ,:)=      0.00021710795484646138591610188464622454_wp*scale(:)
  c(-6 ,:)=     -0.0021356291838797986414312219042358542_wp*scale(:)
  c(-5 ,:)=      0.00713761218453631422925717625758502986_wp*scale(:)
  c(-4 ,:)=     -0.0284696165863973422636410524436931061_wp*scale(:)
  c(-3 ,:)=      0.14327329352510759457155821037742893841_wp*scale(:)
  c(-2 ,:)=     -0.42498050943780130143385739554118569733_wp*scale(:)
  c(-1 ,:)=      0.65703074007121357894896358254040272157_wp*scale(:)
  c( 0 ,:)=     -0.42081655293724308770919536332797729898_wp*scale(:)
  c( 1 ,:)=     -0.21716117505137104371463587747283267899_wp*scale(:)
  c( 2 ,:)=      0.63457035267892488185929915286969303251_wp*scale(:)
  c( 3 ,:)=     -0.53298223962800395684936080758073568406_wp*scale(:)
  c( 4 ,:)=      0.23370490631751294307619384973520033236_wp*scale(:)
  c( 5 ,:)=     -0.05657736973328755112051544344507997075_wp*scale(:)
  c( 6 ,:)=      0.0080872029411844780634067667008050127_wp*scale(:)
  c( 7 ,:)=     -0.00093423623304808664741804536808932984_wp*scale(:)
  c( 8 ,:)=      0.00005075807947289728306309081261461095_wp*scale(:)
  c( 9 ,:)=     -4.62561497463184262755416490048242e-6_wp*scale(:)
  c( 10,:)=      6.3919128513793415587294752371778e-7_wp*scale(:)
  c( 11,:)=      1.87909235155149902916133888931e-8_wp*scale(:)
  c( 12,:)=      1.04757345962781829480207861447155543883e-10_wp*scale(:)
  c( 13,:)=     -4.84665690596158959648731537084025836e-13_wp*scale(:)
  c( 14,:)=      1.2392629629188986192855777620877e-17_wp*scale(:)
  !  <psi|D^2|phi_i>
  do i=-14,14
     b(i,:)=c(-i,:)
  enddo
  !<psi|D^2|psi_i>
  e(0,:)=   -24.875846029392331358907766562_wp*scale(:)
  e(1,:)=   -7.1440597663471719869313377994_wp*scale(:)
  e(2,:)=   -0.04251705323669172315864542163525830944_wp*scale(:)
  e(3,:)=   -0.26995931336279126953587091167128839196_wp*scale(:)
  e(4,:)=    0.08207454169225172612513390763444496516_wp*scale(:)
  e(5,:)=   -0.02207327034586634477996701627614752761_wp*scale(:)
  e(6,:)=    0.00409765642831595181639002667514310145_wp*scale(:)
  e(7,:)=   -0.00045167920287507774929432548999880117_wp*scale(:)
  e(8,:)=    0.00002398228524507599670405555359023135_wp*scale(:)
  e(9,:)=   -2.0904234952920365957922889447361e-6_wp*scale(:)
  e(10,:)=   3.7230763047369275848791496973044e-7_wp*scale(:)
  e(11,:)=   1.05857055496741470373494132287e-8_wp*scale(:)
  e(12,:)=   5.8138798302825405479592506674648873655e-11_wp*scale(:)
  e(13,:)=  -2.70800493626319438269856689037647576e-13_wp*scale(:)
  e(14,:)=   6.924474940639200152025730585882e-18_wp*scale(:)
  do i=1,14
     e(-i,:)=e(i,:)
  enddo
END SUBROUTINE prepare_sdc_wire
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! New conv routine taken out from convolut_kinetic_slab_T. These routines 
! significantly speedup the version included in convolut_kinetic_slab_T.
subroutine conv_kin_x_new_wire(n1,x,y,ndat,lowfil,lupfil,fil,ekin)!,mod_arr1)
  use module_defs, only: wp
  implicit none
  integer,intent(in)::n1,ndat,lowfil,lupfil
  real(wp), intent(in), dimension(lowfil:lupfil,3) :: fil
  real(wp),intent(in):: x(0:n1,ndat)
  real(wp),intent(inout)::y(0:n1,ndat),ekin
!  integer,intent(in) :: mod_arr1(lowfil:n1+lupfil)   
  real(wp) :: tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12,ekin_tmp
  integer :: i,i1,l,j

  ekin_tmp=0.0_wp
  !$omp parallel do default(private)&
  !$omp shared(ndat,n1,lowfil,lupfil,x,y,fil)&
  !$omp reduction(+:ekin_tmp)
  do i=0,ndat/12-1
    do i1=0,n1
      tt1=0.e0_wp
      tt2=0.e0_wp
      tt3=0.e0_wp
      tt4=0.e0_wp
      tt5=0.e0_wp
      tt6=0.e0_wp
      tt7=0.e0_wp
      tt8=0.e0_wp
      tt9 =0.e0_wp
      tt10=0.e0_wp
      tt11=0.e0_wp
      tt12=0.e0_wp
      do l=max(lowfil,-i1),min(lupfil,n1-i1)!lowfil,lupfil
         j=i1+l !j=mod_arr1(i1+l)
         tt1=tt1+x(j,i*12+1)*fil(l,1)
         tt2=tt2+x(j,i*12+2)*fil(l,1)
         tt3=tt3+x(j,i*12+3)*fil(l,1)
         tt4=tt4+x(j,i*12+4)*fil(l,1)
         tt5=tt5+x(j,i*12+5)*fil(l,1)
         tt6=tt6+x(j,i*12+6)*fil(l,1)
         tt7=tt7+x(j,i*12+7)*fil(l,1)
         tt8=tt8+x(j,i*12+8)*fil(l,1)
         tt9 =tt9 +x(j,i*12+9 )*fil(l,1)
         tt10=tt10+x(j,i*12+10)*fil(l,1)
         tt11=tt11+x(j,i*12+11)*fil(l,1)
         tt12=tt12+x(j,i*12+12)*fil(l,1)
      enddo
      y(i1,i*12+1)=y(i1,i*12+1)+tt1;    ekin_tmp=ekin_tmp+tt1*x(i1,i*12+1)
      y(i1,i*12+2)=y(i1,i*12+2)+tt2;    ekin_tmp=ekin_tmp+tt2*x(i1,i*12+2)
      y(i1,i*12+3)=y(i1,i*12+3)+tt3;    ekin_tmp=ekin_tmp+tt3*x(i1,i*12+3)
      y(i1,i*12+4)=y(i1,i*12+4)+tt4;    ekin_tmp=ekin_tmp+tt4*x(i1,i*12+4)
      y(i1,i*12+5)=y(i1,i*12+5)+tt5;    ekin_tmp=ekin_tmp+tt5*x(i1,i*12+5)
      y(i1,i*12+6)=y(i1,i*12+6)+tt6;    ekin_tmp=ekin_tmp+tt6*x(i1,i*12+6)
      y(i1,i*12+7)=y(i1,i*12+7)+tt7;    ekin_tmp=ekin_tmp+tt7*x(i1,i*12+7)
      y(i1,i*12+8)=y(i1,i*12+8)+tt8;    ekin_tmp=ekin_tmp+tt8*x(i1,i*12+8)
      y(i1,i*12+9 )=y(i1,i*12+9 )+tt9 ; ekin_tmp=ekin_tmp+tt9 *x(i1,i*12+9 )
      y(i1,i*12+10)=y(i1,i*12+10)+tt10; ekin_tmp=ekin_tmp+tt10*x(i1,i*12+10)
      y(i1,i*12+11)=y(i1,i*12+11)+tt11; ekin_tmp=ekin_tmp+tt11*x(i1,i*12+11)
      y(i1,i*12+12)=y(i1,i*12+12)+tt12; ekin_tmp=ekin_tmp+tt12*x(i1,i*12+12)
    enddo
  enddo
  !$omp end parallel do

  ekin=ekin+ekin_tmp

  ekin_tmp=0.0_wp
  !$omp parallel do default(private) &
  !$omp shared(ndat,n1,lowfil,lupfil,x,y,fil) &
  !$omp reduction(+:ekin_tmp)
  do i=(ndat/12)*12+1,ndat
    do i1=0,n1
      tt=0.e0_wp
      !do l=lowfil,lupfil
      !  j=mod_arr1(i1+l)
      do l=max(lowfil,-i1),min(lupfil,n1-i1)!lowfil,lupfil
         j=i1+l !j=mod_arr1(i1+l)
         tt=tt+x(j   ,i)*fil(l,1)
      enddo
      y(i1,i)=y(i1,i)+tt ; ekin_tmp=ekin_tmp+tt*x(i1,i)
    enddo
  enddo
  !$omp end parallel do

  ekin=ekin+ekin_tmp
END SUBROUTINE conv_kin_x_new_wire
