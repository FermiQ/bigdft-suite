!> @file
!!  Time-Dependent DFT ai la Casida
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Calculate the coupling matrix for the TD-DFT a la Casida
subroutine tddft_casida(iproc,nproc,dir_output,atoms,rxyz,n3p,n3parr,Glr,tddft_approach,orbs,&
     orbsv,fxc,pkernelseq,psi,psiv,exc_fac,bitp)
  use module_base
  use module_types
  use locregs
  use box
  implicit none
  integer, intent(in) :: iproc,nproc,n3p!,i3s
!  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: atoms
  type(orbitals_data), intent(in) :: orbs,orbsv
  type(locreg_descriptors), intent(in) :: Glr
  character(len=4), intent(in) :: tddft_approach
  character(len=*), intent(in) :: dir_output
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  real(wp), dimension(Glr%d%n1i,Glr%d%n2i,n3p), intent(in) :: fxc
  type(coulomb_operator), intent(inout) :: pkernelseq
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  real(wp), dimension(orbsv%npsidim_orbs), intent(in) :: psiv
  real(gp), intent(in) :: exc_fac
  type(box_iterator) :: bitp
  !local variables
  real(gp), dimension(3) :: chargec
  real(wp), dimension(:), allocatable :: psirocc,psirvirt

  !temporary call to the coupling matrix calculation
  psirocc = f_malloc(max(max(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*orbs%norbp, n3parr(0)*orbs%norb), 1),id='psirocc')
  psirvirt = f_malloc(max(max(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*orbsv%norbp, n3parr(0)*orbsv%norb), 1),id='psirvirt')

  call prepare_psirocc(iproc,nproc,Glr,orbs,n3p,n3parr,psi,psirocc)

  call prepare_psirocc(iproc,nproc,Glr,orbsv,n3p,n3parr,psiv,psirvirt)

  call center_of_charge(atoms,rxyz,chargec)

  call calculate_coupling_matrix(iproc,nproc,dir_output,bitp,tddft_approach,orbs%nspin,Glr%d%n1i*Glr%d%n2i*n3p,orbs,orbsv,&
       chargec,pkernelseq,fxc,psirocc,psirvirt)

!!$  call coupling_matrix_prelim(iproc,nproc,atoms%astruct%geocode,tddft_approach,orbs%nspin,Glr,orbs,orbsv,&
!!$       i3s,n3p,hxh,hyh,hzh,chargec,pkernelseq,fxc,psirocc,psirvirt,exc_fac)

  call f_free(psirocc)

  call f_free(psirvirt)

end subroutine tddft_casida



!!$
!!$!test the calculation of the coupling matrix with the OP2P embedded scheme
!!$if (orbs%nspin==2) then
!!$   sfac=1.0_gp
!!$   ngroup=2
!!$else
!!$   sfac=0.5_gp
!!$   ngroup=1
!!$end if
!!$!construct the OP2P scheme and test it
!!$!use temporaryly tyhe nvrct_par array
!!$nobj_occ = f_malloc((/ 0.to.nproc-1, 1.to.ngroup /),id='nobj_par')
!!$nobj_virt = f_malloc((/ 0.to.nproc-1, 1.to.ngroup /),id='nobj_par')
!!$
!!$call fill_nobj_par_for_OP2P(nproc,ngroup,orbs,nobj_occ)
!!$call fill_nobj_par_for_OP2P(nproc,ngroup,orbsv,nobj_virt)
!!$
!!$ndim=(wfd%nvctr_c+7*wfd_nvctr_f)*orbs%nspinor
!!$
!!$call initialize_OP2P_data(OP2P_occ,bigdft_mpi%mpi_comm,iproc,nproc,ngroup,ndim,nobj_occ,gpudirect,.false.)
!!$
!!$!start iterating on the occupied wavefunctions
!!$call set_OP2P_iterator(iproc,OP2P,iter_occ,orbs%norbp,psi_occ,res_fake) !res_fake should have the same size of psi_occ
!!$
!!$OP2P_occ: do
!!$   call OP2P_communication_step(iproc,OP2P_occ,iter_occ)
!!$   if(igpu==1) call synchronize() !to be removed in future version
!!$   if (iter_occ%event == OP2P_EXIT) exit OP2P_occ
!!$
!!$   !now we have available a set of the occupied orbitals
!!$   !let us now exchange all the virtual orbitals between the processors
!!$   call initialize_OP2P_data(OP2P_virt,bigdft_mpi%mpi_comm,iproc,nproc,ngroup,ndim,nobj_virt,gpudirect,.false.)
!!$   !iterating on the virtual wavefunctions
!!$   call set_OP2P_iterator(iproc,OP2P_virt,iter_virt,orbsv%norbp,psi_virt,res_virt) !res_virt should have the same size of psi_virt
!!$
!!$   OP2P_virt: do
!!$      call OP2P_communication_step(iproc,OP2P_virt,iter_virt)
!!$      if(igpu==1) call synchronize() !to be removed in future version
!!$      if (iter_occ%event == OP2P_EXIT) exit OP2P_virt
!!$      
!!$      !here the coupling matrix can be calculated
!!$      !we can calculate rho_{p,alpha}
!!$      !we can apply the poisson solver -> V_{p,alpha}
!!$
!!$      !calculate all the lines of the coupling matrix V_{palpha} rho_alphap in a distributed way
!!$      !by performing on-the-fly-construction of the density
!!$
!!$
!!$   end do OP2P_virt
!!$
!!$   call free_OP2P_data(OP2P_occ)
!!$
!!$
!!$   !inform
!!$   if (iproc == 0) then
!!$      call OP2P_info(iter_occ,OP2P_occ,prc,tel,trm)
!!$      call yaml_comment('Exact exchange calculation: '+prc**'(i3)'+&
!!$           '%; Time (s): Elapsed='+tel**'(1pg12.2)'&
!!$           +', Remaining='+trm**'(1pg12.2)')
!!$   end if
!!$end do OP2P_occ
!!$
!!$call free_OP2P_data(OP2P_occ)
!!$
!!$call f_free(nobj_occ)
!!$call f_free(nobj_virt)
