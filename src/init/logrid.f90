!cleaned version of the logrid_old.f90 in the unused directory (with newmethod=.true.)

subroutine make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     ibxy_c,ibzzx_c,ibyyzz_c,ibxy_f,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibyz_c,ibzxx_c,ibxxyy_c,ibyz_f,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
  !    creates complicated ib arrays    
  use module_base
  implicit none
  integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer i1,i2,i3,nt,m1,m2,m3,i_stat,i_all

  integer,intent(in):: ibyz_c(2,0:n2,0:n3),ibxy_c(2,0:n1,0:n2)
  integer,intent(in):: ibyz_f(2,0:n2,0:n3),ibxy_f(2,0:n1,0:n2)

  !    for shrink:    
  integer,intent(out):: ibzzx_c(2,-14:2*n3+16,0:n1) 
  integer,intent(out):: ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

  integer,intent(out):: ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
  integer,intent(out):: ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
  integer,intent(out):: ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

  !    for grow:    
  integer,intent(out):: ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
  integer,intent(out):: ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

  integer,intent(out):: ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer,intent(out):: ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer,intent(out):: ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

  logical,allocatable:: logrid_big(:)

  !    for real space:
  integer,intent(out):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)

  allocate(logrid_big((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
  call memocc(i_stat,product(shape(logrid_big))*kind(logrid_big),'logrid_big','make_all_ib')

  m1=nfu1-nfl1
  m2=nfu2-nfl2
  m3=nfu3-nfl3

  !   (0:n3,-14:2*n1+16,-14:2*n2+16) from grow
  !   (-14:2*n2+16,-14:2*n3+16,0:n1) from shrink
  !   (-14:2*n1+16,-14:2*n2+16,-14:2*n3+16) from real space

  !	for shrink:
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        ibxy_ff(:,i1,i2)=ibxy_f(:,i1,i2)
     enddo
  enddo

  call make_ib_inv(logrid_big, ibxy_c,ibzzx_c,ibyyzz_c,0,n1,0,n2,0,n3)
  call make_ib_inv(logrid_big,ibxy_ff,ibzzx_f,ibyyzz_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)

  !	for realspace:
  !-14:2*n2+16,-14:2*n3+16
  do i3=-14,2*n3+16
     do i2=-14,2*n2+16
        if (ibyyzz_c(1,i2,i3).ne.1000) then
           ibyyzz_r(1,i2,i3)=2*ibyyzz_c(1,i2,i3)
           ibyyzz_r(2,i2,i3)=2*ibyyzz_c(2,i2,i3)+30
        else
           ibyyzz_r(1,i2,i3)=1000
           ibyyzz_r(2,i2,i3)=-1000
        endif
     enddo
  enddo
  call squares(ibyyzz_r,2*n2+30,2*n3+30)

  !    for grow:

  do i2=nfl2,nfu2
     do i3=nfl3,nfu3
        ibyz_ff(:,i2,i3)=ibyz_f(:,i2,i3)
     enddo
  enddo

  call make_ib_c(logrid_big,ibyz_c,ibzxx_c,ibxxyy_c,n1,n2,n3)

  call make_ib(logrid_big,ibyz_ff,ibzxx_f,ibxxyy_f,&
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)

  call squares_1d(ibxxyy_f,2*nfl1-14,2*nfu1+16,2*nfl2-14,2*nfu2+16)

  i_all=-product(shape(logrid_big))*kind(logrid_big)
  deallocate(logrid_big,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_big','make_all_ib')

end subroutine make_all_ib


subroutine make_ib_inv(logrid_big,ibxy,ibzzx,ibyyzz,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
  !    This subroutine mimics the comb_grow_f one
  implicit none
  integer nt,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer,intent(out):: ibxy(2,nfl1:nfu1,nfl2:nfu2)
  integer,intent(out):: ibzzx(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
  integer,intent(out):: ibyyzz(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)
  logical logrid_big(nfl3:nfu3,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)! work array

  ! I3,i1,i2 -> i1,i2,i3 
  nt=(nfu1-nfl1+1)*(nfu2-nfl2+1)
  call ib_to_logrid_inv(ibxy,logrid_big,nfl3,nfu3,nt)

  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*(nfu3-nfl3)+31)*(nfu1-nfl1+1)
  call ib_from_logrid_inv(ibzzx,logrid_big,nfl2,nfu2,nt)
  call ib_to_logrid_inv(ibzzx,logrid_big,nfl2,nfu2,nt)

  ! I1,I2,I3  -> I2,I3,i1
  nt=(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31)
  call ib_from_logrid_inv( ibyyzz,logrid_big,nfl1,nfu1,nt)

end subroutine make_ib_inv


subroutine ib_to_logrid_inv(ib,logrid,nfl,nfu,ndat)
  ! This one mimics the comb_rot_grow_f_loc
  implicit none
  integer ndat,nfl,nfu,l,i
  integer ib(2,ndat)! input
  logical logrid(-14+2*nfl:2*nfu+16,ndat)! output

  logrid=.false.

  do l=1,ndat
     do i = 2*ib(1,l)-14 , 2*ib(2,l)+16
        logrid(i,l)=.true.
     enddo
  enddo

end subroutine ib_to_logrid_inv

subroutine ib_from_logrid_inv(ib,logrid,ml1,mu1,ndat)
  ! mimics the bounds subroutine    
  implicit none
  integer i,i1
  integer ml1,mu1,ndat
  integer ib(2,ndat)
  logical logrid(ndat,ml1:mu1)

  do i=1,ndat
     ib(1,i)= 1000
     ib(2,i)=-1000

     inner1:do i1=ml1,mu1
        if (logrid(i,i1)) then 
           ib(1,i)=i1
           exit inner1
        endif
     enddo inner1

     inner2:do i1=mu1,ml1,-1
        if (logrid(i,i1)) then 
           ib(2,i)=i1
           exit inner2
        endif
     enddo inner2
  enddo

end subroutine ib_from_logrid_inv


subroutine make_ib_c(logrid_big,ibyz,ibzxx,ibxxyy,n1,n2,n3)
  !    This subroutine mimics the comb_grow_f one
  implicit none
  integer nt,n1,n2,n3
  integer ibyz(2,0:n2,0:n3)! input
  integer ibzxx(2,0:n3,-14:2*n1+16)!output
  integer ibxxyy(2,-14:2*n1+16,-14:2*n2+16)!output
  logical logrid_big(0:n3,-14:2*n1+16,-14:2*n2+16)! work array

  call squares(ibyz,n2,n3)
  ! i1,i2,i3 -> i2,i3,I1
  nt=(n2+1)*(n3+1)
  call ib_to_logrid_rot(ibyz,logrid_big,0,n1,nt)

  ! i2,i3,I1 -> i3,I1,I2
  nt=(n3+1)*(2*n1+31)
  call ib_from_logrid(ibzxx,logrid_big,0,n2,nt)
  call squares(ibzxx,n3,2*n1+30)

  call ib_to_logrid_rot(ibzxx,logrid_big,0,n2,nt)

  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+31)*(2*n2+31)
  call ib_from_logrid(ibxxyy,logrid_big,0,n3,nt)
  call squares(ibxxyy,2*n1+30,2*n2+30)

end subroutine make_ib_c


subroutine make_ib(logrid_big,ibyz,ibzxx,ibxxyy,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
  !    This subroutine mimics the comb_grow_f one
  implicit none
  integer nt,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer ibyz(  2,nfl2:nfu2,nfl3:nfu3)! input
  integer ibzxx( 2,          nfl3:nfu3,2*nfl1-14:2*nfu1+16)!output
  integer ibxxyy(2,                    2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)!output
  logical logrid_big(           nfl3:nfu3,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)! work array

  ! i1,i2,i3 -> i2,i3,I1
  nt=(nfu2-nfl2+1)*(nfu3-nfl3+1)
  call ib_to_logrid_rot(  ibyz,logrid_big,nfl1,nfu1,nt)

  ! i2,i3,I1 -> i3,I1,I2
  nt=(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)
  call ib_from_logrid( ibzxx,logrid_big,nfl2,nfu2,nt)
  call ib_to_logrid_rot( ibzxx,logrid_big,nfl2,nfu2,nt)

  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31)
  call ib_from_logrid(ibxxyy,logrid_big,nfl3,nfu3,nt)

end subroutine make_ib

subroutine ib_to_logrid_rot(ib,logrid,nfl,nfu,ndat)
  ! This one mimics the comb_rot_grow_f_loc
  implicit none
  integer ndat,nfl,nfu,l,i
  integer ib(2,ndat)! input
  logical logrid(ndat,-14+2*nfl:2*nfu+16)! output

  logrid=.false.

  do l=1,ndat
     do i = 2*ib(1,l)-14 , 2*ib(2,l)+16
        logrid(l,i)=.true.
     enddo
  enddo

end subroutine ib_to_logrid_rot

subroutine ib_from_logrid(ib,logrid,ml1,mu1,ndat)
  ! mimics the bounds subroutine    
  implicit none
  integer i,i1
  integer ml1,mu1,ndat
  integer ib(2,ndat)
  logical logrid(ml1:mu1,ndat)

  do i=1,ndat
     ib(1,i)= 1000
     ib(2,i)=-1000

     inner1:do i1=ml1,mu1
        if (logrid(i1,i)) then 
           ib(1,i)=i1
           exit inner1
        endif
     enddo inner1

     inner2:do i1=mu1,ml1,-1
        if (logrid(i1,i)) then 
           ib(2,i)=i1
           exit inner2
        endif
     enddo inner2
  enddo

end subroutine ib_from_logrid

subroutine squares_1d(ib,nfl2,nfu2,nfl3,nfu3)
  ! modifies the ib array 
  ! so that it is made up of blocks of size 2
  ! the localization region is enlarged as a result

  ! works for even nfl2 only
  implicit none
  integer,intent(in)::nfl2,nfu2,nfl3,nfu3
  integer,intent(inout)::ib(2,nfl2:nfu2,nfl3:nfu3)

  integer i2,i3,ii2,ii3,ibmin,ibmax

  do i3=nfl3,nfu3
     do i2=nfl2/2,(nfu2-1)/2
        ii2=2*i2
        ibmin=min(ib(1,ii2,i3),ib(1,ii2+1,i3))

        ib(1,ii2,i3)=ibmin
        ib(1,ii2+1,i3)=ibmin

        ibmax=max(ib(2,ii2,i3),ib(2,ii2+1,i3))

        ib(2,ii2,i3)=ibmax
        ib(2,ii2+1,i3)=ibmax
     enddo
  enddo
end subroutine squares_1d


subroutine squares(ib,n2,n3)
  ! modifies the ib array 
  ! so that it is made up of squares 2x2
  ! the localization region is enlarged as a result
  implicit none
  integer,intent(in)::n2,n3
  integer,intent(inout)::ib(2,0:n2,0:n3)

  integer i2,i3,ii2,ii3,ibmin,ibmax

  do i3=0,(n3-1)/2
     do i2=0,(n2-1)/2
        ii3=2*i3
        ii2=2*i2
        ibmin=min(ib(1,ii2,ii3),ib(1,ii2+1,ii3),&
             ib(1,ii2,ii3+1),ib(1,ii2+1,ii3+1))

        ib(1,ii2,ii3)=ibmin
        ib(1,ii2+1,ii3)=ibmin
        ib(1,ii2,ii3+1)=ibmin
        ib(1,ii2+1,ii3+1)=ibmin

        ibmax=max(ib(2,ii2,ii3),ib(2,ii2+1,ii3),&
             ib(2,ii2,ii3+1),ib(2,ii2+1,ii3+1))

        ib(2,ii2,ii3)=ibmax
        ib(2,ii2+1,ii3)=ibmax
        ib(2,ii2,ii3+1)=ibmax
        ib(2,ii2+1,ii3+1)=ibmax
     enddo
  enddo
end subroutine squares
