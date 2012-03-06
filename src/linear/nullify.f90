subroutine nullify_linearParameters(lin)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_linearParameters
  implicit none

  ! Calling arguments
  type(linearParameters),intent(out):: lin

  nullify(lin%potentialPrefac)
  nullify(lin%locrad)
  !nullify(lin%lphiRestart)
  !nullify(lin%lphiold)
  !nullify(lin%lhphiold)
  !nullify(lin%hamold)
  call nullify_orbitals_data(lin%orbs)
  call nullify_orbitals_data(lin%gorbs)
  call nullify_communications_arrays(lin%comms)
  call nullify_communications_arrays(lin%gcomms)
  nullify(lin%norbsPerType)
  !call nullify_p2pCommsSumrho(lin%comsr)
  call nullify_p2pComms(lin%comsr)
  !call nullify_p2pCommsGatherPot(lin%comgp)
  call nullify_p2pComms(lin%comgp)
  call nullify_largeBasis(lin%lb)
  call nullify_local_zone_descriptors(lin%lzd)
  !call nullify_p2pCommsOrthonormality(lin%comon)
  call nullify_p2pComms(lin%comon)
  call nullify_overlapParameters(lin%op)
  call nullify_matrixDescriptors(lin%mad)

end subroutine nullify_linearParameters

subroutine nullifyInputLinparameters(lin)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  type(linearInputParameters),intent(inout):: lin

  nullify(lin%locrad)
  nullify(lin%potentialPrefac)
  nullify(lin%potentialPrefac_lowaccuracy)
  nullify(lin%potentialPrefac_highaccuracy)
  nullify(lin%norbsPerType)

end subroutine nullifyInputLinparameters


subroutine nullify_p2pComms(p2pcomm)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_p2pComms
  implicit none

  ! Calling argument
  type(p2pComms),intent(out):: p2pcomm

  nullify(p2pcomm%noverlaps)
  nullify(p2pcomm%overlaps)
  nullify(p2pcomm%istarr)
  nullify(p2pcomm%istrarr)
  nullify(p2pcomm%sendBuf)
  nullify(p2pcomm%recvBuf)
  nullify(p2pcomm%auxarray)
  nullify(p2pcomm%comarr)
  nullify(p2pcomm%communComplete)
  nullify(p2pcomm%computComplete)
  nullify(p2pcomm%startingindex)
  nullify(p2pcomm%ise3)
  nullify(p2pcomm%requests)
end subroutine nullify_p2pComms


!!subroutine nullify_p2pCommsSumrho(comsr)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => nullify_p2pCommsSumrho
!!  implicit none
!!
!!  ! Calling argument
!!  type(p2pCommsSumrho),intent(out):: comsr
!!
!!  nullify(comsr%noverlaps)
!!  nullify(comsr%overlaps)
!!  nullify(comsr%istarr)
!!  nullify(comsr%istrarr)
!!  nullify(comsr%sendBuf)
!!  nullify(comsr%recvBuf)
!!  nullify(comsr%comarr)
!!  nullify(comsr%communComplete)
!!  nullify(comsr%computComplete)
!!  nullify(comsr%startingindex)
!!end subroutine nullify_p2pCommsSumrho
!!
!!
!!subroutine nullify_p2pCommsGatherPot(comgp)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => nullify_p2pCommsGatherPot
!!  implicit none
!!
!!  ! Calling argument
!!  type(p2pCommsGatherPot),intent(out):: comgp
!!
!!  nullify(comgp%noverlaps)
!!  nullify(comgp%overlaps)
!!  nullify(comgp%ise3)
!!  nullify(comgp%comarr)
!!  nullify(comgp%recvBuf)
!!  nullify(comgp%communComplete)
!!
!!end subroutine nullify_p2pCommsGatherPot


subroutine nullify_largeBasis(lb)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_largeBasis
  implicit none

  ! Calling argument
  type(largeBasis),intent(out):: lb
  call nullify_communications_arrays(lb%comms)
  call nullify_communications_arrays(lb%gcomms)
  call nullify_orbitals_data(lb%orbs)
  call nullify_orbitals_data(lb%gorbs)
  !call nullify_p2pCommsRepartition(lb%comrp)
  call nullify_p2pComms(lb%comrp)
  !call nullify_p2pCommsOrthonormality(lb%comon)
  call nullify_p2pComms(lb%comon)
  call nullify_overlapParameters(lb%op)
  !call nullify_p2pCommsGatherPot(lb%comgp)
  call nullify_p2pComms(lb%comgp)

end subroutine nullify_largeBasis


!!!subroutine nullify_p2pCommsRepartition(comrp)
!!!  use module_base
!!!  use module_types
!!!  use module_interfaces, exceptThisOne => nullify_p2pCommsRepartition
!!!  implicit none
!!!
!!!  ! Calling arguments
!!!  type(p2pCommsRepartition),intent(out):: comrp
!!!
!!!  nullify(comrp%comarr)
!!!  nullify(comrp%communComplete)
!!!  nullify(comrp%requests)
!!!
!!!end subroutine nullify_p2pCommsRepartition


!!subroutine nullify_p2pCommsOrthonormality(comon)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => nullify_p2pCommsOrthonormality
!!  implicit none
!!
!!  ! Calling argument
!!  type(p2pCommsOrthonormality),intent(out):: comon
!!
!!  nullify(comon%noverlaps)
!!  !!nullify(comon%overlaps)
!!  nullify(comon%comarr)
!!  nullify(comon%sendBuf)
!!  nullify(comon%recvBuf)
!!  nullify(comon%communComplete)
!!  nullify(comon%requests)
!!
!!end subroutine nullify_p2pCommsOrthonormality


subroutine nullify_overlapParameters(op)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_overlapParameters
  implicit none

  ! Calling argument
  type(overlapParameters),intent(out):: op

  nullify(op%noverlaps)
!  nullify(op%indexExpand)
!  nullify(op%indexExtract)
  nullify(op%overlaps)
  nullify(op%indexInRecvBuf)
  nullify(op%indexInSendBuf)
  !nullify(op%olr)
  nullify(op%wfd_overlap)
!  nullify(op%expseg)

end subroutine nullify_overlapParameters


!!subroutine nullify_expansionSegments(expseg)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => nullify_expansionSegments
!!  implicit none
!!
!!  ! Calling argument
!!  type(expansionSegments),intent(out):: expseg
!!
!!  nullify(expseg%segborders)
!!
!!end subroutine nullify_expansionSegments



subroutine nullify_linearInputGuess(lig)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_linearInputGuess
  implicit none

  ! Calling argument
  type(linearInputGuess),intent(out):: lig

  call nullify_local_zone_descriptors(lig%lzdig)
  call nullify_local_zone_descriptors(lig%lzdGauss)
  call nullify_orbitals_data(lig%orbsig)
  call nullify_orbitals_data(lig%orbsGauss)
  !call nullify_p2pCommsOrthonormality(lig%comon)
  call nullify_p2pComms(lig%comon)
  call nullify_overlapParameters(lig%op)
  !call nullify_p2pCommsGatherPot(lig%comgp)
  call nullify_p2pComms(lig%comgp)
  call nullify_matrixDescriptors(lig%mad)

end subroutine nullify_linearInputGuess


subroutine nullify_matrixDescriptors(mad)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_matrixDescriptors
  implicit none

  ! Calling argument
  type(matrixDescriptors),intent(out):: mad

  nullify(mad%keyv)
  nullify(mad%keyvmatmul)
  nullify(mad%nsegline)
  nullify(mad%keyg)
  nullify(mad%keygmatmul)
  nullify(mad%keygline)

end subroutine nullify_matrixDescriptors



subroutine nullify_local_zone_descriptors(lzd)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_local_zone_descriptors
  implicit none

  ! Calling arguments
  type(local_zone_descriptors),intent(out):: lzd
 
  nullify(lzd%rxyz)
  call nullify_locreg_descriptors(lzd%glr)
  nullify(lzd%llr)
  nullify(lzd%doHamAppl)
  !!nullify(lzd%cutoffweight) 
  !!nullify(lzd%cutoffweight)
 
end subroutine nullify_local_zone_descriptors



subroutine nullify_orbitals_data(orbs)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data),intent(out):: orbs
  
  nullify(orbs%norb_par)
  nullify(orbs%iokpt)
  nullify(orbs%ikptproc)
  nullify(orbs%inwhichlocreg)
  !nullify(orbs%inWhichLocregP)
  nullify(orbs%onWhichMPI)
  nullify(orbs%isorb_par)
  nullify(orbs%eval)
  nullify(orbs%occup)
  nullify(orbs%spinsgn)
  nullify(orbs%kwgts)
  nullify(orbs%kpts)
  nullify(orbs%ispot)
  orbs%npsidim_orbs=1
  orbs%npsidim_comp=1

end subroutine nullify_orbitals_data


subroutine nullify_communications_arrays(comms)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(communications_arrays),intent(out):: comms

  nullify(comms%ncntd)
  nullify(comms%ncntt)
  nullify(comms%ndspld)
  nullify(comms%ndsplt)
  nullify(comms%nvctr_par)
  
end subroutine nullify_communications_arrays


subroutine nullify_locreg_descriptors(lr)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_locreg_descriptors
  implicit none

  ! Calling arguments
  type(locreg_descriptors),intent(out):: lr


  !if(associated(lr%projflg)) then
     nullify(lr%projflg)
  !end if

  call nullify_wavefunctions_descriptors(lr%wfd)
  call nullify_convolutions_bounds(lr%bounds)

end subroutine nullify_locreg_descriptors


subroutine nullify_wavefunctions_descriptors(wfd)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(wavefunctions_descriptors),intent(out):: wfd

  !if(associated(wfd%keyg)) then
     nullify(wfd%keygloc)
     nullify(wfd%keyglob)
     nullify(wfd%keyvloc)
     nullify(wfd%keyvglob)
  !end if
end subroutine nullify_wavefunctions_descriptors


subroutine nullify_convolutions_bounds(bounds)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_convolutions_bounds
  implicit none

  ! Calling arguments
  type(convolutions_bounds),intent(out):: bounds

  call nullify_kinetic_bounds(bounds%kb)
  call nullify_shrink_bounds(bounds%sb)
  call nullify_grow_bounds(bounds%gb)
  !if(associated(bounds%ibyyzz_r)) then
     nullify(bounds%ibyyzz_r)
  !end if
end subroutine nullify_convolutions_bounds



subroutine nullify_kinetic_bounds(kb)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(kinetic_bounds),intent(out):: kb

  !if(associated(kb%ibyz_c))then
     nullify(kb%ibyz_c)
     nullify(kb%ibxz_c)
     nullify(kb%ibxy_c)
     nullify(kb%ibyz_f)
     nullify(kb%ibxz_f)
     nullify(kb%ibxy_f)
  !end if
end subroutine nullify_kinetic_bounds



subroutine nullify_shrink_bounds(sb)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(shrink_bounds),intent(out):: sb

  !if(associated(sb%ibzzx_c)) then
     nullify(sb%ibzzx_c)
     nullify(sb%ibyyzz_c)
     nullify(sb%ibxy_ff)
     nullify(sb%ibzzx_f)
     nullify(sb%ibyyzz_f)
  !end if

end subroutine nullify_shrink_bounds



subroutine nullify_grow_bounds(gb)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(grow_bounds),intent(out):: gb

  !if(associated(gb%ibzxx_c)) then
     nullify(gb%ibzxx_c)
     nullify(gb%ibxxyy_c)
     nullify(gb%ibyz_ff)
     nullify(gb%ibzxx_f)
     nullify(gb%ibxxyy_f)
  !end if
end subroutine nullify_grow_bounds



subroutine nullify_nonlocal_psp_descriptors(nlpspd)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(nonlocal_psp_descriptors),intent(out):: nlpspd

  nlpspd%natoms=0
  nullify(nlpspd%plr)

!!$  nullify(nlpspd%nvctr_p)
!!$  nullify(nlpspd%nseg_p)
!!$  nullify(nlpspd%keyv_p)
!!$  nullify(nlpspd%keyg_p)
!!$  nullify(nlpspd%nboxp_c)
!!$  nullify(nlpspd%nboxp_f)

end subroutine nullify_nonlocal_psp_descriptors



subroutine nullify_matrixMinimization(matmin)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(matrixMinimization),intent(out):: matmin

  nullify(matmin%mlr)
  nullify(matmin%inWhichLocregExtracted)
  nullify(matmin%inWhichLocregOnMPI)
  nullify(matmin%indexInLocreg)

end subroutine nullify_matrixMinimization


subroutine nullify_matrixLocalizationRegion(mlr)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(matrixLocalizationRegion),intent(out):: mlr

  nullify(mlr%indexInGlobal)

end subroutine nullify_matrixLocalizationRegion
