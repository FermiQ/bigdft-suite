subroutine abi_md_isothermal(acell, acell_next, amass, bmass, dtion, etotal, etotal0, &
     & fcart, iatfix, itime, ktemp, mditemp, me, mttk_vars, natom, nnos, optcell, &
     & qmass, rprim, rprimd, rprim_next, rprimd_next, strten, strtarget, ucvol, &
     & ucvol_next, vel, vel_nexthalf, vmass, xcart, xcart_next, xred_next)

  use abi_defs_basis
  use abi_interfaces_lowlevel
  use abi_interfaces_geometry

  implicit none

  integer, intent(in) :: natom, itime, nnos, optcell, me
  real(dp), intent(in) :: mditemp, dtion, ktemp, bmass, vmass
  real(dp), intent(out) :: ucvol_next, etotal
  real(dp), intent(inout) :: etotal0, ucvol
  integer, intent(in) :: iatfix(3, natom)
  real(dp), intent(in) :: amass(natom), qmass(nnos)
  real(dp), intent(in) :: acell(3), rprim(3,3), rprimd(3,3)
  real(dp), intent(out) :: acell_next(3), rprim_next(3,3), rprimd_next(3,3)
  real(dp), intent(in) :: strtarget(6), strten(6)
  real(dp), intent(inout) :: vel(3,natom), fcart(3, natom)
  real(dp), intent(out) :: vel_nexthalf(3, natom)
  real(dp), intent(in) :: xcart(3, natom)
  real(dp), intent(out) :: xred_next(3,natom), xcart_next(3,natom)
  type(abi_mttk_type), intent(inout) :: mttk_vars

  integer,parameter :: lwork=8
  real(dp),parameter :: esh2=one/six,esh4=esh2/20._dp,esh6=esh4/42._dp
  real(dp),parameter :: esh8=esh6/72._dp,nosetol=tol10,v2tol=tol8
  integer :: iatom, idim, ierr
  real(dp) :: ekin, mttk_aloc, mttk_aloc2, mttk_bloc, polysh, vlogv
  real(dp) :: gmet(3,3), gprimd(3,3), rmet(3,3), mttk_alc(3)
  real(dp) :: mttk_alc2(3),mttk_blc(3),mttk_psh(3),mttk_tv(3,3),mttk_ubox(3,3)
  real(dp) :: mttk_uu(3),mttk_uv(3),mttk_veig(3),mttk_vt(3,3), work(lwork)
  real(dp), allocatable :: fcart_m(:,:)
  character(len=500) :: message

  !   write(6,*)'moldyn', nnos, qmass(:), bmass, vmass
  if(itime==0) then
     allocate(mttk_vars%glogs(nnos),mttk_vars%vlogs(nnos),&
          &    mttk_vars%xlogs(nnos))
     mttk_vars%glogs(:)=zero; mttk_vars%vlogs(:)=zero; mttk_vars%xlogs(:)=zero
     mttk_vars%vboxg(:,:)=zero
     vlogv=zero

     call abi_md_isokinetic_init(amass, mditemp, natom, vel)
  end if

  allocate(fcart_m(3,natom))
  !  XG070613 : Do not take away the following line , seems needed for the pathscale compiler
  if (me == 0) write(6,*)'moldyn',mttk_vars%vboxg(:,:)
  !  ##### sub  case optcell==0 Isothermal Ensemble ###########
  if(optcell==0) then
     !   There is no evolution of cell
     acell_next(:)=acell(:)
     ucvol_next=ucvol
     rprim_next(:,:)=rprim(:,:)
     rprimd_next(:,:)=rprimd(:,:)
     !   Update Thermostat variables and scale velocitie
     call abi_isotemp(amass,dtion,ekin,iatfix,ktemp,mttk_vars,natom,nnos,qmass,vel)
     !   Half velocity step
     do idim=1,3
        fcart_m(idim,:)=fcart(idim,:)/amass(:)
     end do
     vel_nexthalf(:,:)=vel(:,:)+dtion/two*fcart_m(:,:)
     !   New positions
     xcart_next(:,:)=xcart(:,:)+vel_nexthalf(:,:)*dtion
     !   Convert back to xred (reduced coordinates)
     call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)
     !   Computation of the forces for the new positions
     !   Compute LDA forces (big loop), fcart_m is used as dummy argument for fred
     call scfloop_main(acell, etotal, fcart, fcart_m, itime, me, natom, rprimd, xred_next)
     !   Next Half velocity step
     do idim=1,3
        fcart_m(idim,:)=fcart(idim,:)/amass(:)
     end do
     vel(:,:)=vel_nexthalf(:,:)+dtion/two*fcart_m(:,:)
     !   Update Thermostat variables and velocity
     call abi_isotemp(amass,dtion,ekin,iatfix,ktemp,mttk_vars,natom,nnos,qmass,vel)
     if(itime==0) etotal0=ekin+etotal
     write(message, '(a,es18.10,a,es18.10,a,es18.10)' )&
          & ' dtotal=',(ekin+etotal)-etotal0, &
          & ', ekin=',ekin,', epot=', etotal
     call abi_wrtout(std_out,message,'COLL')
     !   End of sub case optcell=0
     !   ##### sub  case optcell==1 Isothermal-Isenthalpic Ensemble (homogeneous cell deformation)##########
  else if (optcell==1) then
     !   Only homogeneous evolution of cell
     !   Evolution of cell we keep rprim constant
     rprim_next(:,:)=rprim(:,:)
     !   Update Thermostat variables and velocity
     call isopress(amass,dtion,ekin,iatfix,ktemp,natom,nnos,qmass,strten,strtarget,ucvol,mttk_vars,vel,vlogv,vmass)
     !   Half velocity step
     do idim=1,3
        fcart_m(idim,:)=fcart(idim,:)/amass(:)
     end do
     vel_nexthalf(:,:)=vel(:,:)+dtion/two*fcart_m(:,:)
     !   New positions
     mttk_aloc=exp(dtion/two*vlogv)
     mttk_aloc2=(vlogv*dtion/two)**2
     polysh=(((esh8*mttk_aloc2+esh6)*mttk_aloc2+esh4)*mttk_aloc2+esh2)*mttk_aloc2+one
     mttk_bloc=mttk_aloc*polysh*dtion
     xcart_next(:,:)=xcart(:,:)*mttk_aloc**2+vel_nexthalf(:,:)*mttk_bloc
     !   Update the volume and related quantities
     acell_next(:)=acell(:)*exp(dtion*vlogv)
     !   ucvol=ucvol*exp(dtion*vlogv)
     call abi_mkrdim(acell_next,rprim,rprimd_next)
     call abi_metric(gmet,gprimd,-1,rmet,rprimd_next,ucvol_next)
     !   Convert back to xred (reduced coordinates)
     call abi_xredxcart(natom,-1,rprimd_next,xcart_next,xred_next)
     !   Compute LDA forces (big loop), fcart_m is used as dummy argument for fred
     call scfloop_main(acell_next, etotal, fcart, fcart_m, itime, me, natom, &
          & rprimd_next, xred_next)
     !   Next Half velocity step
     do idim=1,3
        fcart_m(idim,:)=fcart(idim,:)/amass(:)
     end do
     vel(:,:)=vel_nexthalf(:,:)+dtion/two*fcart_m(:,:)
     !   Update Thermostat variables and velocity
     call isopress(amass,dtion,ekin,iatfix,ktemp,natom,nnos,qmass,strten,strtarget,ucvol,mttk_vars,vel,vlogv,vmass)
     if(itime==0) etotal0=ekin+etotal
     write(message, '(a,es18.10,a,es18.10,a,es18.10)' )&
          & ' dtotal=',(ekin+etotal)-etotal0, &
          & ', ekin=',ekin,', epot=', etotal
     call abi_wrtout(std_out,message,'COLL')
     !   End of sub case optcell = 1
     !   ##### sub  case optcell==2 Isothermal-Isenthalpic Ensemble (full cell deformation)##########
  else if (optcell==2) then
     !   Fisrt half step for extended variables
     call abi_isostress(amass,bmass,dtion,ekin,iatfix,ktemp, natom, nnos, qmass,strten,strtarget,ucvol,vel,mttk_vars)
     !   Half velocity step
     do idim=1,3
        fcart_m(idim,:)=fcart(idim,:)/amass(:)
     end do
     vel_nexthalf(:,:)=vel(:,:)+dtion/two*fcart_m(:,:)
     !   New positions
     mttk_vt(:,:)=mttk_vars%vboxg(:,:)
     call dsyev('V','U',3,mttk_vt,3,mttk_veig,work,lwork,ierr)
     mttk_tv(:,:)=transpose(mttk_vt)
     mttk_alc(:)=exp(dtion/two*mttk_veig(:))
     mttk_alc2(:)=(mttk_veig(:)*dtion/two)**2
     mttk_psh(:)=(((esh8*mttk_alc2(:)+esh6)*mttk_alc2(:)+esh4)*mttk_alc2(:)+esh2)*mttk_alc2(:)+one
     mttk_blc(:)=mttk_alc(:)*mttk_psh(:)*dtion
     !   Update the positions
     do iatom=1,natom
        mttk_uu(:)=matmul(mttk_tv,xcart(:,iatom))
        mttk_uv(:)=matmul(mttk_tv,vel_nexthalf(:,iatom))
        mttk_uu(:)=mttk_uu(:)*mttk_alc(:)**2+mttk_uv(:)*mttk_blc(:)
        xcart_next(:,iatom)=matmul(mttk_vt,mttk_uu)
     end do
     !   Update the box (rprimd and rprim)
     mttk_ubox(:,:)=matmul(mttk_tv,rprimd)
     do idim=1,3
        mttk_ubox(:,idim)=mttk_ubox(:,idim)*mttk_alc(:)**2
     end do
     rprimd_next(:,:)=matmul(mttk_vt,mttk_ubox)
     do idim=1,3
        rprim_next(idim,:)=rprimd_next(idim,:)/acell(:)
     end do
     !   Update the volume
     call abi_metric(gmet,gprimd,-1,rmet,rprimd_next,ucvol)
     !   Convert back to xred (reduced coordinates)
     call abi_xredxcart(natom,-1,rprimd_next,xcart_next,xred_next)

     !   Compute LDA forces (big loop), fcart_m is used as dummy argument for fred
     call scfloop_main(acell_next, etotal, fcart, fcart_m, itime, me, natom, &
          & rprimd_next, xred_next)
     !   Next Half velocity step
     do idim=1,3
        fcart_m(idim,:)=fcart(idim,:)/amass(:)
     end do
     vel(:,:)=vel_nexthalf(:,:)+dtion/two*fcart_m(:,:)
     !   Next half step for extended variables
     call abi_isostress(amass,bmass,dtion,ekin,iatfix,ktemp, natom, nnos, qmass,strten,strtarget,ucvol,vel,mttk_vars)
     if(itime==0) etotal0=ekin+etotal
     write(message, '(a,es18.10,a,es18.10,a,es18.10)' )&
          & ' dtotal=',(ekin+etotal)-etotal0, &
          & ', ekin=',ekin,', epot=', etotal
     call abi_wrtout(std_out,message,'COLL')
     !   Evolution of cell and volumr
     acell_next(:)=acell(:)
     ucvol_next=ucvol
     !   End of sub case optcell=2
  else
     write(message, '(a,a,a,a,i12,a,a)' ) ch10,&
          &    ' moldyn : BUG -',ch10,&
          &    '  Disallowed value for optcell=',optcell,ch10,&
          &    '  Allowed values with ionmov=13 : 0 to 2.'
     call abi_wrtout(std_out,message,'COLL')
     call abi_leave_new('COLL')
  end if
  deallocate(fcart_m)

end subroutine abi_md_isothermal
