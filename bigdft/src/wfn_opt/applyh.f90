!> @file
!!  Routine to calculate the action of the hamiltonian
!! @author
!!   Copyright (C) 2005-2013 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Calculate the action of the local hamiltonian on the orbitals
subroutine local_hamiltonian_old(iproc,nproc,npsidim_orbs,orbs,Lzd,hx,hy,hz,&
     ipotmethod,confdatarr,pot,psi,hpsi,pkernel,xc,alphaSIC,ekin_sum,epot_sum,eSIC_DC)!,&
!     dpbox,potential,comgp)
  use module_base
  use module_dpbox, only: denspot_distribution
  use module_types
  use module_xc
  use locreg_operations
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,npsidim_orbs
  integer, intent(in) :: ipotmethod !< Method which has to be chosen for applying the potential to the wavefunctions in the real space form:
                                    !! 0 is the traditional potential application
                                    !! 1 is the application of the exact exchange (which has to be precomputed and stored in the potential array)
                                    !! 2 is the application of the Perdew-Zunger SIC
                                    !! 3 is the application of the Non-Koopman's correction SIC
  real(gp), intent(in) :: hx,hy,hz,alphaSIC
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
  type(xc_info), intent(in) :: xc
  real(wp), dimension(npsidim_orbs), intent(in) :: psi              !< This dimension will be modified
  real(wp), dimension(*) :: pot                             !< the potential, with the dimension compatible with the ipotmethod flag
  real(gp), intent(out) :: ekin_sum,epot_sum,eSIC_DC
  real(wp), dimension(npsidim_orbs), intent(inout) :: hpsi
  type(coulomb_operator), intent(in) :: pkernel                     !< the PSolver kernel which should be associated for the SIC schemes
!!$  type(denspot_distribution),intent(in),optional :: dpbox
!!$  real(wp), dimension(*), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
!!$  type(p2pComms),intent(inout), optional:: comgp
  !!real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
  !!real(wp), dimension(max(dpbox%ndimrhopot,orbs%nspin)), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian'
  logical :: dosome
  integer :: iorb,npot,ispot,ispsi,ilr,ilr_orb!,jproc,ierr
  real(wp) :: exctXcoeff
  real(gp) :: ekin,epot,kx,ky,kz,eSICi,eSIC_DCi !n(c) etest
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable :: vsicpsir
  real(wp), dimension(:,:), allocatable :: psir
  !!write(*,*) 'condition',(present(dpbox) .and. present(potential) .and. present(comgp))

  call f_routine(id='local_hamiltonian')

  epot=0.d0
  ekin=0.d0

  !some checks
  exctXcoeff=xc_exctXfac(xc)

  if (exctXcoeff /= 0.0_gp .neqv. ipotmethod ==1) then
     if (iproc==0) write(*,*)&
          'ERROR (local_hamiltonian): potential method not compatible with exact exchange'
     stop
  end if

  if (.not.(associated(pkernel%kernel) .and. alphaSIC /=0.0_gp) .and. ipotmethod == 2) then
     if (iproc==0) write(*,*)&
          'ERROR (local_hamiltonian): potential method not compatible with SIC'
     stop
  end if

  ekin_sum=0.0_gp
  epot_sum=0.0_gp
  eSIC_DC=0.0_gp
!!$do jproc=0,nproc-1
!!$call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
!!$if (jproc==iproc) then
  !loop on the localisation regions (so to create one work array set per lr)
  loop_lr: do ilr=1,Lzd%nlr
    !check if this localisation region is used by one of the orbitals
    dosome=.false.
    do iorb=1,orbs%norbp
      dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
      if (dosome) exit
    end do
    if (.not. dosome) cycle loop_lr
      
    !components of the potential (four or one, depending on the spin)
    npot=orbs%nspinor
    if (orbs%nspinor == 2) npot=1
   
    ! Wavefunction in real space
    psir = f_malloc0((/ Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i, orbs%nspinor /),id='psir')

    call initialize_work_arrays_locham(Lzd%Llr(ilr),orbs%nspinor,.true.,wrk_lh)  
  
    ! wavefunction after application of the self-interaction potential
    if (ipotmethod == 2 .or. ipotmethod == 3) then
      vsicpsir = f_malloc((/ Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i, orbs%nspinor /),id='vsicpsir')
    end if

    ispsi=1
    loop_orbs: do iorb=1,orbs%norbp
      ilr_orb=orbs%inwhichlocreg(iorb+orbs%isorb)
      if (ilr_orb /= ilr) then
        ispsi=ispsi+&
             (Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor
        cycle loop_orbs
      end if
      
!!$      print *,'iorb+orbs%isorb,BEFORE',iorb+orbs%isorb,&
!!$                sum(psi(ispsi:&
!!$                ispsi+(Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor-1))

        
      call daub_to_isf_locham(orbs%nspinor,Lzd%Llr(ilr),wrk_lh,psi(ispsi),psir(1,1))

      !calculate the ODP, to be added to VPsi array
   
      !Perdew-Zunger SIC scheme
      eSIC_DCi=0.0_gp
      if (ipotmethod == 2) then
         !in this scheme the application of the potential is already done
         call PZ_SIC_potential(iorb,Lzd%Llr(ilr),orbs,xc,&
              0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psir,vsicpsir,eSICi,eSIC_DCi)
      !NonKoopmans' correction scheme
      else if (ipotmethod == 3) then 
         !in this scheme first we have calculated the potential then we apply it
         call vcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,&
              psir(1,1),1,vsicpsir(1,1),1)
         !for the moment the ODP is supposed to be valid only with one lr
         call psir_to_vpsi(npot,orbs%nspinor,Lzd%Llr(ilr),&
              pot(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspin+&
              (iorb-1)*Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor+1),&
              vsicpsir,eSICi)
      end if
   
      call psir_to_vpsi(npot,orbs%nspinor,Lzd%llr(ilr),&
           pot(orbs%ispot(iorb)),psir(1,1),epot,confdata=confdatarr(iorb))

!!$      !this ispot has to be better defined inside denspot structure
!!$      print *,'orbs, epot',orbs%isorb+iorb,epot,orbs%ispot(iorb),&
!!$           sum(pot(orbs%ispot(iorb):orbs%ispot(iorb):Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i-1))
   
      !ODP treatment (valid only for the nlr=1 case)
      if (ipotmethod==1) then !Exact Exchange
         ispot=1+Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*(orbs%nspin+iorb-1)
         !add to the psir function the part of the potential coming from the exact exchange
         call axpy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,exctXcoeff,pot(ispot),1,psir(1,1),1)
      else if (ipotmethod == 2) then !PZ scheme
         !subtract the sic potential from the vpsi function
         call axpy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,-alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
         !add the SIC correction to the potential energy
         epot=epot-alphaSIC*eSICi
         !accumulate the Double-Counted SIC energy
         eSIC_DC=eSIC_DC+alphaSIC*eSIC_DCi
      else if (ipotmethod == 3) then !NK scheme
         !add the sic potential from the vpsi function
         call axpy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
         epot=epot+alphaSIC*eSICi
         !accumulate the Double-Counted SIC energy
         eSIC_DC=eSIC_DC+alphaSIC*orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eSICi
      end if
   
      !apply the kinetic term, sum with the potential and transform back to Daubechies basis
      !k-point values, if present
      kx=orbs%kpts(1,orbs%iokpt(iorb))
      ky=orbs%kpts(2,orbs%iokpt(iorb))
      kz=orbs%kpts(3,orbs%iokpt(iorb))
 
      call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,Lzd%Llr(ilr),wrk_lh,&
           psir(1,1),hpsi(ispsi),ekin)

!!$      print *,'iorb+orbs%isorb,AFTER',iorb+orbs%isorb,&
!!$                sum(hpsi(ispsi:&
!!$                ispsi+(Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor-1))
    

      ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
      epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot
           !print *,'iorb+orbs%isorb',iorb+orbs%isorb,ekin,epot
      ispsi=ispsi+&
           (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
      !print *,'iorb,epot',orbs%isorb+iorb,epot

    enddo loop_orbs
   
    !deallocations of work arrays
    call f_free(psir)

    if (ipotmethod == 2 .or. ipotmethod ==3) then
       call f_free(vsicpsir)
    end if
    call deallocate_work_arrays_locham(wrk_lh)
   
  end do loop_lr
!!$end if
!!$end do

  call f_release_routine()

END SUBROUTINE local_hamiltonian_old

!> Calculate the action of the local potential on the orbitals
!! @param ipotmethod Indicates the method which has to be chosen for applying the potential to the wavefunctions in the 
!!                   real space form:
!!                   0 is the traditional potential application
!!                   1 is the application of the exact exchange (which has to be precomputed and stored in the potential array)
!!                   2 is the application of the Perdew-Zunger SIC
!!                   3 is the application of the Non-Koopman's correction SIC
subroutine psi_to_vlocpsi(iproc,npsidim_orbs,orbs,Lzd,&
     ipotmethod,confdatarr,pot,psi,vpsi,pkernel,xc,alphaSIC,epot_sum,evSIC,comgp,vpsi_noconf,econf_sum)
  use module_base
  use module_types
  use module_xc
  use locreg_operations
  use rhopotential, only: extract_potential
  implicit none
  integer, intent(in) :: iproc,ipotmethod,npsidim_orbs
  real(gp), intent(in) :: alphaSIC
  type(xc_info), intent(in) :: xc
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
  real(wp), dimension(npsidim_orbs), intent(in) :: psi !this dimension will be modified
  real(wp), dimension(*) :: pot !< the potential, with the dimension compatible with the ipotmethod flag
  real(gp), intent(out) :: epot_sum,evSIC
  real(wp), dimension(npsidim_orbs), intent(inout) :: vpsi
  type(coulomb_operator), intent(in) :: pkernel !< the PSolver kernel which should be associated for the SIC schemes
  type(p2pcomms),intent(in) :: comgp
  real(wp), dimension(npsidim_orbs), intent(inout),optional :: vpsi_noconf
  real(gp),intent(out),optional :: econf_sum
  !local variables
  character(len=*), parameter :: subname='psi_to_vlocpsi'
  logical :: dosome
  integer :: iorb,npot,ispot,ispsi,ilr,ilr_orb,nbox,nvctr,ispinor,ispin,size_lpot
  real(wp) :: exctXcoeff
  real(gp) :: epot,eSICi,eSIC_DCi,econf !n(c) etest
  type(workarr_sumrho) :: w
  real(wp), dimension(:,:), allocatable :: psir,vsicpsir,psir_noconf

  call f_routine(id='psi_to_vlocpsi')


  !some checks
  exctXcoeff=xc_exctXfac(xc)

  if (exctXcoeff /= 0.0_gp .neqv. ipotmethod ==1) then
     if (iproc==0) write(*,*)&
          'ERROR (local_hamiltonian): potential method not compatible with exact exchange'
     stop
  end if

  if (.not.(associated(pkernel%kernel) .and. alphaSIC /=0.0_gp) .and. ipotmethod == 2) then
     if (iproc==0) write(*,*)&
          'ERROR (local_hamiltonian): potential method not compatible with SIC'
     stop
  end if

  epot_sum=0.0_gp
  evSIC=0.0_gp
  if (present(econf_sum)) then
      econf_sum=0.0_gp
  end if

  call initialize_work_arrays_sumrho(lzd%nlr,lzd%llr,.true.,w)

  !loop on the localisation regions (so to create one work array set per lr)
  loop_lr: do ilr=1,Lzd%nlr
     !check if this localisation region is used by one of the orbitals
     dosome=.false.
     do iorb=1,orbs%norbp
        dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
        if (dosome) then
            exit
        end if
     end do
     if (.not. dosome) cycle loop_lr

     !initialise the work arrays
     call initialize_work_arrays_sumrho(lzd%llr(ilr),.false.,w)

     !box elements size
     nbox=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i

     !components of the potential (four or one, depending on the spin)
     npot=orbs%nspinor
     if (orbs%nspinor == 2) npot=1

     ! Wavefunction in real space
     psir = f_malloc0((/ nbox, orbs%nspinor /),id='psir')

     if (present(vpsi_noconf)) then
         psir_noconf = f_malloc((/ nbox, orbs%nspinor /),id='psir_noconf')
     end if

     !call to_zero(nbox*orbs%nspinor,psir(1,1))

     ! wavefunction after application of the self-interaction potential
     if (ipotmethod == 2 .or. ipotmethod == 3) then
        vsicpsir = f_malloc((/ nbox, orbs%nspinor /),id='vsicpsir')
     end if

  !n(c) etest=0.0_gp

  ispsi=1
  loop_orbs: do iorb=1,orbs%norbp
     ilr_orb=orbs%inwhichlocreg(iorb+orbs%isorb)
     nvctr=Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f
     if (ilr_orb /= ilr) then
        ispsi=ispsi+nvctr*orbs%nspinor
        cycle loop_orbs
     end if
     
     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form
     do ispinor=1,orbs%nspinor
        call daub_to_isf(Lzd%Llr(ilr),w,psi(ispsi+nvctr*(ispinor-1)),psir(1,ispinor))
     end do

     !ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspin+iorb-1)
     !etest=etest+dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,pot(ispot),1,psir(1,1),1)
     !print *,'epot, iorb,iproc,norbp',iproc,orbs%norbp,iorb

     !calculate the ODP, to be added to VPsi array

     !Perdew-Zunger SIC scheme
     eSIC_DCi=0.0_gp
     if (ipotmethod == 2) then
        !in this scheme the application of the potential is already done
        call PZ_SIC_potential(iorb,Lzd%Llr(ilr),orbs,xc,&
             0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
             pkernel,psir,vsicpsir,eSICi,eSIC_DCi)
     !NonKoopmans' correction scheme
     else if (ipotmethod == 3) then 
        !in this scheme first we have calculated the potential then we apply it
        call vcopy(nbox*orbs%nspinor,psir(1,1),1,vsicpsir(1,1),1)
        !for the moment the ODP is supposed to be valid only with one lr
        call psir_to_vpsi(npot,orbs%nspinor,Lzd%Llr(ilr),&
             pot(nbox*(orbs%nspin+(iorb-1)*orbs%nspinor)+1),&
             vsicpsir,eSICi)
     end if

     ! For the linear scaling case
     if (orbs%spinsgn(iorb+orbs%isorb)>0.d0) then
         ispin=1
     else
         ispin=2
     end if
     size_Lpot = Lzd%Llr(ilr_orb)%d%n1i*Lzd%Llr(ilr_orb)%d%n2i*Lzd%Llr(ilr_orb)%d%n3i

     call extract_potential(ispin, ilr_orb, size_lpot, lzd, pot, comgp)

     !apply the potential to the psir wavefunction and calculate potential energy
     if (present(vpsi_noconf)) then
         if (.not.present(econf_sum)) then
             call f_err_throw('econf must be present when psir_noconf is present')
         end if
         call vcopy(nbox*orbs%nspinor, psir(1,1), 1, psir_noconf(1,1), 1)
         call psir_to_vpsi(npot,orbs%nspinor,Lzd%Llr(ilr),&
              pot,psir,epot,confdata=confdatarr(iorb),vpsir_noconf=psir_noconf,econf=econf)
     else
         call psir_to_vpsi(npot,orbs%nspinor,Lzd%Llr(ilr),&
              pot,psir,epot,confdata=confdatarr(iorb))
     end if
     !!do i_stat=1,lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
     !!    write(1000+ilr_orb,*) orbs%ispot(iorb)+i_stat-1, pot(orbs%ispot(iorb)+i_stat-1)
     !!end do
     !this ispot has to be better defined inside denspot structure

     !ODP treatment (valid only for the nlr=1 case)
     if (ipotmethod==1) then !Exact Exchange
        ispot=1+nbox*(orbs%nspin+iorb-1)
        !add to the psir function the part of the potential coming from the exact exchange
        call axpy(nbox,exctXcoeff,pot(ispot),1,psir(1,1),1)
     else if (ipotmethod == 2) then !PZ scheme
        !subtract the sic potential from the vpsi function
        call axpy(nbox*orbs%nspinor,-alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
        !add the SIC correction to the potential energy
        epot=epot-alphaSIC*eSICi
        !accumulate the Double-Counted SIC energy
        evSIC=evSIC+alphaSIC*eSIC_DCi
     else if (ipotmethod == 3) then !NK scheme
        !add the sic potential from the vpsi function
        call axpy(nbox*orbs%nspinor,alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
        epot=epot+alphaSIC*eSICi
        !accumulate the Double-Counted SIC energy
        evSIC=evSIC+alphaSIC*orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eSICi
     end if

     do ispinor=1,orbs%nspinor
         call isf_to_daub(Lzd%Llr(ilr),w,psir(1,ispinor),vpsi(ispsi+nvctr*(ispinor-1)))
         if (present(vpsi_noconf)) then
             call isf_to_daub(Lzd%Llr(ilr),w,psir_noconf(1,ispinor),vpsi_noconf(ispsi+nvctr*(ispinor-1)))
         end if
     end do

     epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot
     if (present(econf_sum)) then
         econf_sum=econf_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*econf
     end if
     ispsi=ispsi+nvctr*orbs%nspinor
  enddo loop_orbs

  !deallocations of work arrays
  call f_free(psir)
  if (present(vpsi_noconf)) then
      call f_free(psir_noconf)
  end if
  if (ipotmethod == 2 .or. ipotmethod ==3) then
     call f_free(vsicpsir)
  end if

end do loop_lr

call deallocate_work_arrays_sumrho(w)

call f_release_routine()

END SUBROUTINE psi_to_vlocpsi


subroutine psi_to_kinpsi(iproc,npsidim_orbs,orbs,lzd,psi,hpsi,ekin_sum)
  use module_base
  use module_types
  use locreg_operations
  implicit none
  integer, intent(in) :: iproc,npsidim_orbs
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  real(wp), dimension(npsidim_orbs), intent(in) :: psi
  real(gp), intent(out) :: ekin_sum
  real(wp), dimension(npsidim_orbs), intent(inout) :: hpsi

  !local variables
  character(len=*), parameter :: subname='psi_to_kinpsi'
  logical :: dosome
  integer :: iorb,ispsi,ilr,ilr_orb
  real(gp) :: ekin
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable :: psir
  real(gp) :: kx,ky,kz


  ekin=0.d0
  ekin_sum=0.0_gp

  call initialize_work_arrays_locham(lzd%nlr,lzd%llr,orbs%nspinor,.true.,wrk_lh)  

  !loop on the localisation regions (so to create one work array set per lr)
  loop_lr: do ilr=1,Lzd%nlr
    !check if this localisation region is used by one of the orbitals
    dosome=.false.
    do iorb=1,orbs%norbp
      dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
      if (dosome) then
          exit
      end if
    end do
    if (.not. dosome) cycle loop_lr
   
    ! Wavefunction in real space
    psir = f_malloc0((/ Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i, orbs%nspinor /),id='psir')

    !initialise the work arrays
    call initialize_work_arrays_locham(Lzd%Llr(ilr),orbs%nspinor,.false.,wrk_lh)  

   
    ispsi=1
    loop_orbs: do iorb=1,orbs%norbp
      ilr_orb=orbs%inwhichlocreg(iorb+orbs%isorb)
      if (ilr_orb /= ilr) then
        ispsi=ispsi+&
             (Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor
        cycle loop_orbs
      end if
        
      !call daub_to_isf_locham(orbs%nspinor,Lzd%Llr(ilr),wrk_lh,psi(ispsi),psir(1,1))

      kx=orbs%kpts(1,orbs%iokpt(iorb))
      ky=orbs%kpts(2,orbs%iokpt(iorb))
      kz=orbs%kpts(3,orbs%iokpt(iorb))

      !call isf_to_daub_kinetic(lzd%hgrids(1),lzd%hgrids(2),lzd%hgrids(3),kx,ky,kz,orbs%nspinor,Lzd%Llr(ilr),wrk_lh,&
      !      psir(1,1),hpsi(ispsi),ekin)
      call psi_to_tpsi(lzd%hgrids,orbs%kpts(1,orbs%iokpt(iorb)),orbs%nspinor,&
           Lzd%Llr(ilr),psi(ispsi),wrk_lh,hpsi(ispsi),ekin)
   
      ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin

      ispsi=ispsi+&
           (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

    enddo loop_orbs

    call f_free(psir)

  end do loop_lr

  call deallocate_work_arrays_locham(wrk_lh)


end subroutine psi_to_kinpsi



!>   Transpose the wavefunction into a real and imaginary part to be treated with k-points
!!   to be used only when nspinor=2 or 4
!!   here the dimensions are n1->n1+1
subroutine transpose_for_kpoints(nspinor,n1,n2,n3,x,ww,direct)
  use module_base
  implicit none
  logical, intent(in) :: direct
  integer, intent(in) :: nspinor,n1,n2,n3
  real(wp), dimension(nspinor*n1*n2*n3), intent(inout) :: x,ww
  !local variables
  integer :: i1,i2,i3,idx,id,id2,id3,isd,ispinor,it

  !k-points also admitted in non-collinear case
  if (direct) then
     do ispinor=1,nspinor/2
        isd=(ispinor-1)*2*n1*n2*n3
        do idx=1,2
           do i3=1,n3
              id3=(i3-1)*n1*n2
              do i2=1,n2
                 id2=(i2-1)*n1
                 do i1=1,n1
                    id=i1+id2+id3+(idx-1)*n1*n2*n3+isd
                    it=idx+2*(i1-1)+2*id2+2*id3+isd
                    ww(it)=x(id)
                 end do
              end do
           end do
        end do
     end do
  else
     do ispinor=1,nspinor/2
        isd=(ispinor-1)*2*n1*n2*n3
        do idx=1,2
           do i3=1,n3
              id3=(i3-1)*n1*n2
              do i2=1,n2
                 id2=(i2-1)*n1
                 do i1=1,n1
                    id=i1+id2+id3+(idx-1)*n1*n2*n3+isd
                    it=idx+2*(i1-1)+2*id2+2*id3+isd
                    ww(id)=x(it)
                 end do
              end do
           end do
        end do
     end do
  end if
  
  !for mixed precision code it should be changed
  call vcopy(nspinor*n1*n2*n3,ww(1),1,x(1),1)
END SUBROUTINE transpose_for_kpoints


!>   routine for applying the local potentials
!!   supports the non-collinear case, the buffer for tails and different Boundary Conditions
!!   Optimal also for the complex wavefuntion case
!!   might generalize the buffers for two different localization regions, provided that the potential lies in a bigger region
subroutine apply_potential(n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot,psir,pot,epot,&
     ibyyzz_r) !optional
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot
  real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
  real(wp), dimension(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
       -14*nl3:2*n3+1+15*nl3-4*nbuf,npot), intent(in) :: pot
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
  real(gp), intent(out) :: epot
  !local variables
  integer :: i1,i2,i3,i1s,i1e,ispinor
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
  real(gp) :: epot_p
 
  !the Tail treatment is allowed only in the Free BC case
  if (nbuf /= 0 .and. nl1*nl2*nl3 == 0) stop 'NONSENSE: nbuf/=0 only for Free BC'

  epot=0.0_wp
!$omp parallel default(private)&
!$omp shared(pot,psir,n1,n2,n3,epot,ibyyzz_r,nl1,nl2,nl3,nbuf,nspinor)
  !case without bounds
  i1s=-14*nl1
  i1e=2*n1+1+15*nl1
  epot_p=0._gp
!$omp do
  do i3=-14*nl3,2*n3+1+15*nl3
     if (i3 >= -14+2*nbuf .and. i3 <= 2*n3+16-2*nbuf) then !check for the nbuf case
        do i2=-14*nl2,2*n2+1+15*nl2
           if (i2 >= -14+2*nbuf .and. i2 <= 2*n2+16-2*nbuf) then !check for the nbuf case
              !this if statement is inserted here for avoiding code duplication
              !it is to be seen whether the code results to be too much unoptimised
              if (present(ibyyzz_r)) then
                 !in this case we are surely in Free BC
                 !the min is to avoid to calculate for no bounds
                 do i1=-14+2*nbuf,min(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14-1
                    psir(i1,i2,i3,:)=0.0_wp
                 enddo
                 i1s=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf)
                 i1e=min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)
              end if
              !here we put the branchments wrt to the spin
              if (nspinor == 4) then
                 do i1=i1s,i1e
                    !wavefunctions
                    psir1=psir(i1,i2,i3,1)
                    psir2=psir(i1,i2,i3,2)
                    psir3=psir(i1,i2,i3,3)
                    psir4=psir(i1,i2,i3,4)
                    !potentials
                    pot1=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)
                    pot2=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,2)
                    pot3=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,3)
                    pot4=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,4)

                    !diagonal terms
                    tt11=pot1*psir1 !p1
                    tt22=pot1*psir2 !p2
                    tt33=pot4*psir3 !p3
                    tt44=pot4*psir4 !p4
                    !Rab*Rb
                    tt13=pot2*psir3 !p1
                    !Iab*Ib
                    tt14=pot3*psir4 !p1
                    !Rab*Ib
                    tt23=pot2*psir4 !p2
                    !Iab*Rb
                    tt24=pot3*psir3 !p2
                    !Rab*Ra
                    tt31=pot2*psir1 !p3
                    !Iab*Ia
                    tt32=pot3*psir2 !p3
                    !Rab*Ia
                    tt41=pot2*psir2 !p4
                    !Iab*Ra
                    tt42=pot3*psir1 !p4

                    ! Change epot later
                    epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
                         2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3

                    !wavefunction update
                    !p1=h1p1+h2p3-h3p4
                    !p2=h1p2+h2p4+h3p3
                    !p3=h2p1+h3p2+h4p3
                    !p4=h2p2-h3p1+h4p4
                    psir(i1,i2,i3,1)=tt11+tt13-tt14
                    psir(i1,i2,i3,2)=tt22+tt23+tt24
                    psir(i1,i2,i3,3)=tt33+tt31+tt32
                    psir(i1,i2,i3,4)=tt44+tt41-tt42
                 end do
              else
                 do ispinor=1,nspinor
                    do i1=i1s,i1e
                       !the local potential is always real
                       tt=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)*psir(i1,i2,i3,ispinor)
                       epot_p=epot_p+real(tt*psir(i1,i2,i3,ispinor),gp)
                       psir(i1,i2,i3,ispinor)=tt
                    end do
                 end do
              end if
              
              if (present(ibyyzz_r)) then
                 !the max is to avoid the calculation for no bounds
                 do i1=max(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14+1,2*n1+16-2*nbuf
                    psir(i1,i2,i3,:)=0.0_wp
                 enddo
              end if

           else
              do i1=-14,2*n1+16
                 psir(i1,i2,i3,:)=0.0_wp
              enddo
           endif
        enddo
     else
        do i2=-14,2*n2+16
           do i1=-14,2*n1+16
              psir(i1,i2,i3,:)=0.0_wp
           enddo
        enddo
     endif
  enddo
!$omp end do

!$omp critical
  epot=epot+epot_p
!$omp end critical

!$omp end parallel

END SUBROUTINE apply_potential

!> In this routine cprj will be communicated to all processors
subroutine gather_cprj(orbs, paw)
  use module_defs, only: dp
  use module_base, only: bigdft_mpi
  use module_types, only: orbitals_data, paw_objects
  use dynamic_memory
  use wrapper_MPI
  implicit none
  type(orbitals_data), intent(in) :: orbs
  type(paw_objects), intent(inout) :: paw
  ! Local variables:
  integer::iatom,ilmn,iorb,ikpts,jproc
  integer, parameter :: nspinor_ = 1
  ! Tabulated data to be send/received for mpi
  integer,allocatable,dimension(:):: ndsplt
  integer,allocatable,dimension(:):: ncntt 
  ! auxiliar arrays
  real(dp),allocatable,dimension(:,:,:,:)::raux  

  if (bigdft_mpi%nproc > 1) then

     !   Allocate temporary arrays
     ndsplt = f_malloc(0.to.bigdft_mpi%nproc-1,id='ndsplt')
     ncntt = f_malloc(0.to.bigdft_mpi%nproc-1,id='ncntt')
     raux = f_malloc0((/ 2, paw%lmnmax, paw%natom, orbs%norb * orbs%nkpts * nspinor_ /),id='raux')

     !   Set tables for mpi operations:
     !   receive buffer:
     do jproc=0,bigdft_mpi%nproc-1
        ncntt(jproc)=0
        do ikpts=1,orbs%nkpts
           ncntt(jproc)=ncntt(jproc)+&
                2*paw%lmnmax*paw%natom*orbs%norb_par(jproc,ikpts)*nspinor_
        end do
     end do
     !   Displacements table:
     !    ndspld(0)=0
     !    do jproc=1,bigdft_mpi%nproc-1
     !      ndspld(jproc)=ndspld(jproc-1)+ncntd(jproc-1)
     !    end do
     ndsplt(0)=0
     do jproc=1,bigdft_mpi%nproc-1
        ndsplt(jproc)=ndsplt(jproc-1)+ncntt(jproc-1)
     end do

     !   Transfer cprj to raux:
     do iorb=orbs%isorb * nspinor_+1,(orbs%isorb+orbs%norbp) * nspinor_
        do iatom=1,paw%natom
           do ilmn=1,paw%cprj(iatom,iorb)%nlmn
              raux(:,ilmn,iatom,iorb)=&
                   & paw%cprj(iatom,iorb)%cp(:,ilmn)
           end do
        end do
     end do
     !
     !    sendcnt=2*lmnmax*natom*orbs%norbp !N. of data to send
     !    recvcnt=2*lmnmax*natom*orbs%norb  !N. of data to receive
     !   
     !    call MPI_ALLGATHER(raux,sendcnt,MPI_DOUBLE_PRECISION,&
     !&     raux2,recvcnt,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call fmpi_allgather(sendbuf = raux(1,1,1,1), recvcounts = ncntt, &
          & displs = ndsplt, comm = bigdft_mpi%mpi_comm)
!!$     call MPI_ALLGATHERV(raux,ncntd,mpidtypw,&
!!$          &     raux2,ncntt,ndsplt,mpidtypw,MPI_COMM_WORLD,ierr)
     !
     !   Transfer back, raux to cprj:
     do iorb=1,orbs%norb * orbs%nkpts * nspinor_
        do iatom=1,paw%natom
           do ilmn=1,paw%cprj(iatom,iorb)%nlmn
              paw%cprj(iatom,iorb)%cp(:,ilmn)=raux(:,ilmn,iatom,iorb)
           end do
        end do
     end do
     !   Deallocate arrays:
     call f_free(ndsplt)
     call f_free(ncntt)
     call f_free(raux)
  end if


end subroutine gather_cprj

!> Find the starting and ending orbital for kpoint ikpt, and the corresponding nspinor
subroutine orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ikpt
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: isorb,ieorb,nspinor
  !local variables
  integer :: iorb

  !disable starting and ending points for the case no orbitals on a given processor
  if (orbs%norbp == 0) then
     isorb=1
     ieorb=0
  end if

  !find starting orbital
  isorb=1 !default if orbs%norbp==0
  do iorb=1,orbs%norbp
     if (orbs%iokpt(iorb)==ikpt) then
        isorb=iorb
        exit
     end if
  end do

  !find ending orbital
  ieorb=0 !default if orbs%norbp==0
  do iorb=orbs%norbp,1,-1
     if (orbs%iokpt(iorb)==ikpt) then
        ieorb=iorb
        exit
     end if
  end do

  !nspinor for this k-point
  nspinor=orbs%nspinor

END SUBROUTINE orbs_in_kpt


!>   Determine whether the k-point is complex of real
!!   Find the starting and ending orbital for kpoint ikpt, and the corresponding nspinor
subroutine ncplx_kpt(ikpt,orbs,ncplx)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ikpt
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: ncplx
  !local variables
  real(gp) :: kx,ky,kz

  !features of the k-point ikpt
  kx=orbs%kpts(1,ikpt)
  ky=orbs%kpts(2,ikpt)
  kz=orbs%kpts(3,ikpt)

  !evaluate the complexity of the k-point
  if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
     ncplx=1
  else
     ncplx=2
  end if

END SUBROUTINE ncplx_kpt


!!!> Calculate the action of the local hamiltonian on the orbitals
!!subroutine local_hamiltonianParabola(iproc,orbs,lr,hx,hy,hz,&
!!     nspin,pot,psi,hpsi,ekin_sum,epot_sum, nat, rxyz, onWhichAtom, at)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  use libxc_functionals
!!  implicit none
!!  integer, intent(in) :: iproc,nspin
!!  real(gp), intent(in) :: hx,hy,hz
!!  type(orbitals_data), intent(in) :: orbs
!!  type(locreg_descriptors), intent(in) :: lr
!!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi
!!  real(wp), dimension(*) :: pot
!!  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
!!  real(gp), intent(out) :: ekin_sum,epot_sum
!!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: hpsi
!!integer:: nat
!!real(8),dimension(3,nat):: rxyz 
!!integer,dimension(orbs%norbp):: onWhichAtom
!!type(atoms_data), intent(in) :: at
!!  !local variables
!!  character(len=*), parameter :: subname='local_hamiltonian'
!!  integer :: i_all,i_stat,iorb,npot,nsoffset,oidx,ispot
!!  real(wp) :: exctXcoeff
!!  real(gp) :: ekin,epot,kx,ky,kz,etest
!!  type(workarr_locham) :: wrk_lh
!!  real(wp), dimension(:,:), allocatable :: psir
!!real(8):: hxh, hyh, hzh
!!real(8),dimension(3):: rxyzShifted
!!
!!  exctXcoeff=libxc_functionals_exctXfac()
!!
!!  !initialise the work arrays
!!  call initialize_work_arrays_locham(lr,orbs%nspinor,wrk_lh)  
!!
!!  !components of the potential
!!  npot=orbs%nspinor
!!  if (orbs%nspinor == 2) npot=1
!!
!!  ! Wavefunction in real space
!!  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!  call memocc(i_stat,psir,'psir',subname)
!!
!!  call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,psir)
!!
!!  ekin_sum=0.0_gp
!!  epot_sum=0.0_gp
!!
!!  etest=0.0_gp
!!
!!  hxh=hx*.5d0
!!  hyh=hy*.5d0
!!  hzh=hz*.5d0
!!
!!  do iorb=1,orbs%norbp
!!
!!
!!     if(orbs%spinsgn(iorb+orbs%isorb)>0.0_gp .or. nspin == 1 .or. nspin == 4 ) then
!!        nsoffset=1
!!     else
!!        nsoffset=lr%d%n1i*lr%d%n2i*lr%d%n3i+1
!!     end if
!!
!!     oidx=(iorb-1)*orbs%nspinor+1
!!
!!     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
!!     !the psir wavefunction is given in the spinorial form
!!     call daub_to_isf_locham(orbs%nspinor,lr,wrk_lh,psi(1,oidx),psir)
!!
!!     !ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspin+iorb-1)
!!     !etest=etest+dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,pot(ispot),1,psir(1,1),1)
!!     !print *,'epot, iorb,iproc,norbp',iproc,orbs%norbp,iorb,etest
!!
!!     !apply the potential to the psir wavefunction and calculate potential energy
!!     select case(lr%geocode)
!!     case('F')
!!
!!   ! ATTENTION: WITH SHIFTED PARABOLA
!!    rxyzShifted(1)=rxyz(1,onWhichAtom(iorb))+orbs%parabolaShift(1,iorb)
!!    rxyzShifted(2)=rxyz(2,onWhichAtom(iorb))+orbs%parabolaShift(2,iorb)
!!    rxyzShifted(3)=rxyz(3,onWhichAtom(iorb))+orbs%parabolaShift(3,iorb)
!!        call apply_potentialParabola(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!             pot(nsoffset),epot, rxyzShifted, hxh, hyh, hzh, orbs%parabPrefacArr(at%astruct%iatype(onWhichAtom(iorb))), orbs%power, &
!!             lr%bounds%ibyyzz_r) !optional
!!   ! THIS WAS THE ORIGINAL
!!        !call apply_potentialParabola(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!        !     pot(nsoffset),epot, rxyz(1,onWhichAtom(iorb)), hxh, hyh, hzh, orbs%parabPrefacArr(at%astruct%iatype(onWhichAtom(iorb))),  &
!!        !     lr%bounds%ibyyzz_r) !optional
!!
!!        !call apply_potentialParabola(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!        !     pot(nsoffset),epot, rxyz(1,onWhichAtom(iorb)), hxh, hyh, hzh, orbs%parabPrefac,  &
!!        !     lr%bounds%ibyyzz_r) !optional
!!        !call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!        !     pot(nsoffset),epot,&
!!        !     lr%bounds%ibyyzz_r) !optional
!!          
!!     case('P') 
!!        !here the hybrid BC act the same way
!!        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,0,0,0,0,orbs%nspinor,npot,psir,&
!!             pot(nsoffset),epot)
!!
!!     case('S')
!!
!!        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,0,1,0,0,orbs%nspinor,npot,psir,&
!!             pot(nsoffset),epot)
!!     end select
!!
!!     !k-point values, if present
!!     kx=orbs%kpts(1,orbs%iokpt(iorb))
!!     ky=orbs%kpts(2,orbs%iokpt(iorb))
!!     kz=orbs%kpts(3,orbs%iokpt(iorb))
!!
!!     if (exctXcoeff /= 0.0_gp) then
!!        ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspin+iorb-1)
!!        !add to the psir function the part of the potential coming from the exact exchange
!!        call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i,exctXcoeff,pot(ispot),1,psir(1,1),1)
!!     end if
!!
!!     !apply the kinetic term, sum with the potential and transform back to Daubechies basis
!!     call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,lr,wrk_lh,&
!!          psir,hpsi(1,oidx),ekin)
!!
!!     ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
!!     epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot
!!
!!  enddo
!!
!!  !print *,'iproc,etest',etest
!!
!!  !deallocations of work arrays
!!  i_all=-product(shape(psir))*kind(psir)
!!  deallocate(psir,stat=i_stat)
!!  call memocc(i_stat,i_all,'psir',subname)
!!
!!  call deallocate_work_arrays_locham(lr,wrk_lh)
!!
!!END SUBROUTINE local_hamiltonianParabola
!!
!!
!!!> routine for applying the local potentials
!!!! supports the non-collinear case, the buffer for tails and different Boundary Conditions
!!!! Optimal also for the complex wavefuntion case
!!subroutine apply_potentialParabola(n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot,psir,pot,epot, rxyzParab, &
!!     hxh, hyh, hzh, parabPrefac, power, &
!!     ibyyzz_r) !optional
!!  use module_base
!!  implicit none
!!  integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot
!!  real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
!!  real(wp), dimension(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
!!       -14*nl3:2*n3+1+15*nl3-4*nbuf,npot), intent(in) :: pot
!!  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
!!  real(gp), intent(out) :: epot
!!real(8),dimension(3):: rxyzParab
!!real(8):: hxh, hyh, hzh, parabPrefac
!!integer:: power
!!  !local variables
!!  integer :: i1,i2,i3,i1s,i1e,ispinor
!!  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
!!  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
!!  real(gp) :: epot_p
!!real(8):: hxh2, hyh2, hzh2
!!integer:: ix0, iy0, iz0, istat, ipot
!!real(8),dimension(:,:,:,:),allocatable:: potCopy
!!
!!!write(*,'(a,2i7)') '-14*nl3, 2*n3+1+15*nl3', -14*nl3, 2*n3+1+15*nl3
!!!write(*,'(a,2i7)') '-14*nl2, 2*n2+1+15*nl2', -14*nl2, 2*n2+1+15*nl2
!!!write(*,'(a,2i7)') 'max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf), min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)', max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf), min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)
!!allocate(potCopy(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
!!       -14*nl3:2*n3+1+15*nl3-4*nbuf,npot), stat=istat)
!!potCopy=0.d0
!!  
!!  !the Tail treatment is allowed only in the Free BC case
!!  if (nbuf /= 0 .and. nl1*nl2*nl3 == 0) stop 'NONSENSE: nbuf/=0 only for Free BC'
!!
!!  epot=0.0_wp
!!  ! Copy the potential
!!!write(*,*) 'size(pot)', size(pot)
!!  !potCopy=pot
!!  do ipot=1,npot
!!      do i3=-14*nl3,2*n3+1+15*nl3-4*nbuf
!!          do i2=-14*nl2,2*n2+1+15*nl2-4*nbuf
!!              do i1=-14*nl1,2*n1+1+15*nl1-4*nbuf
!! !write(1,'(a,4i8)') 'i1, i2, i3, ipot', i1, i2, i3, ipot
!!                  potCopy(i1,i2,i3,ipot)=pot(i1,i2,i3,ipot)
!!              end do
!!          end do
!!      end do
!!  end do
!!  !potCopy(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
!!  !     -14*nl3:2*n3+1+15*nl3-4*nbuf,1:npot) &
!!  !  =pot(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
!!  !     -14*nl3:2*n3+1+15*nl3-4*nbuf,1:npot)
!!
!!   ix0=nint(rxyzParab(1)/hxh)
!!   iy0=nint(rxyzParab(2)/hyh)
!!   iz0=nint(rxyzParab(3)/hzh)
!!   hxh2=hxh**2
!!   hyh2=hyh**2
!!   hzh2=hzh**2
!!
!!!$omp parallel default(private)&
!!!$omp shared(pot,psir,n1,n2,n3,epot,ibyyzz_r,nl1,nl2,nl3,nbuf,nspinor)
!!  !case without bounds
!!  i1s=-14*nl1
!!  i1e=2*n1+1+15*nl1
!!  epot_p=0._gp
!!!$omp do
!!  do i3=-14*nl3,2*n3+1+15*nl3
!!     if (i3 >= -14+2*nbuf .and. i3 <= 2*n3+16-2*nbuf) then !check for the nbuf case
!!        do i2=-14*nl2,2*n2+1+15*nl2
!!           if (i2 >= -14+2*nbuf .and. i2 <= 2*n2+16-2*nbuf) then !check for the nbuf case
!!              !this if statement is inserted here for avoiding code duplication
!!              !it is to be seen whether the code results to be too much unoptimised
!!              if (present(ibyyzz_r)) then
!!                 !in this case we are surely in Free BC
!!                 !the min is to avoid to calculate for no bounds
!!                 do i1=-14+2*nbuf,min(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14-1
!!                    psir(i1,i2,i3,:)=0.0_wp
!!                 enddo
!!                 i1s=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf)
!!                 i1e=min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)
!!              end if
!!              
!!              !here we put the branchments wrt to the spin
!!              if (nspinor == 4) then
!!                 do i1=i1s,i1e
!!                    !wavefunctions
!!                    psir1=psir(i1,i2,i3,1)
!!                    psir2=psir(i1,i2,i3,2)
!!                    psir3=psir(i1,i2,i3,3)
!!                    psir4=psir(i1,i2,i3,4)
!!                    !potentials
!!                    pot1=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)
!!                    pot2=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,2)
!!                    pot3=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,3)
!!                    pot4=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,4)
!!
!!                    !diagonal terms
!!                    tt11=pot1*psir1 !p1
!!                    tt22=pot1*psir2 !p2
!!                    tt33=pot4*psir3 !p3
!!                    tt44=pot4*psir4 !p4
!!                    !Rab*Rb
!!                    tt13=pot2*psir3 !p1
!!                    !Iab*Ib
!!                    tt14=pot3*psir4 !p1
!!                    !Rab*Ib
!!                    tt23=pot2*psir4 !p2
!!                    !Iab*Rb
!!                    tt24=pot3*psir3 !p2
!!                    !Rab*Ra
!!                    tt31=pot2*psir1 !p3
!!                    !Iab*Ia
!!                    tt32=pot3*psir2 !p3
!!                    !Rab*Ia
!!                    tt41=pot2*psir2 !p4
!!                    !Iab*Ra
!!                    tt42=pot3*psir1 !p4
!!
!!                    ! Change epot later
!!                    epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
!!                         2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3
!!
!!                    !wavefunction update
!!                    !p1=h1p1+h2p3-h3p4
!!                    !p2=h1p2+h2p4+h3p3
!!                    !p3=h2p1+h3p2+h4p3
!!                    !p4=h2p2-h3p1+h4p4
!!                    psir(i1,i2,i3,1)=tt11+tt13-tt14
!!                    psir(i1,i2,i3,2)=tt22+tt23+tt24
!!                    psir(i1,i2,i3,3)=tt33+tt31+tt32
!!                    psir(i1,i2,i3,4)=tt44+tt41-tt42
!!                 end do
!!              else
!!                 do ispinor=1,nspinor
!!                    do i1=i1s,i1e
!!                       !the local potential is always real
!!                       ! Add the parabola to the potential
!!                       !tt=hxh**2*dble(i1-ix0)**2 + hyh**2*dble(i2-iy0)**2 + hzh**2*dble(i3-iz0)**2
!!                       !tt=hxh2*dble(i1-ix0)**2 + hyh2*dble(i2-iy0)**2 + hzh2*dble(i3-iz0)**2
!!                       tt=(hxh*dble(i1)-rxyzParab(1))**2 + (hyh*dble(i2)-rxyzParab(2))**2 + (hzh*dble(i3)-rxyzParab(3))**2
!!                       if(power==2) then
!!                           tt=parabPrefac*tt
!!                       else if(power==4) then
!!                           tt=parabPrefac*tt**2
!!                       else
!!                           write(*,'(a,i0)') "'power' must be 2 or 4, but we found ", power
!!                           stop
!!                       end if
!!                       potCopy(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)=potCopy(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)+tt
!!                       tt=potCopy(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)*psir(i1,i2,i3,ispinor)
!!                       epot_p=epot_p+real(tt*psir(i1,i2,i3,ispinor),gp)
!!                       psir(i1,i2,i3,ispinor)=tt
!!                    end do
!!                 end do
!!              end if
!!              
!!              if (present(ibyyzz_r)) then
!!                 !the max is to avoid the calculation for no bounds
!!                 do i1=max(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14+1,2*n1+16-2*nbuf
!!                    psir(i1,i2,i3,:)=0.0_wp
!!                 enddo
!!              end if
!!
!!           else
!!              do i1=-14,2*n1+16
!!                 psir(i1,i2,i3,:)=0.0_wp
!!              enddo
!!           endif
!!        enddo
!!     else
!!        do i2=-14,2*n2+16
!!           do i1=-14,2*n1+16
!!              psir(i1,i2,i3,:)=0.0_wp
!!           enddo
!!        enddo
!!     endif
!!  enddo
!!!$omp end do
!!
!!!$omp critical
!!  epot=epot+epot_p
!!!$omp end critical
!!
!!!$omp end parallel
!!
!!END SUBROUTINE apply_potentialParabola
