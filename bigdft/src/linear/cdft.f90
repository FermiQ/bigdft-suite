!> @file
!>  Contains files for constrained DFT
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> CDFT: calculates the weight matrix w_ab via the expression S^1/2PS^1/2, where S is the overlap of the whole system
!! CDFT: and P is a projector matrix onto the tmbs of the desired fragment
!! CDFT: for standalone CDFT calculations, assuming one charged fragment, for transfer integrals assuming two fragments
!! CDFT: where we constrain the difference - should later generalize this
subroutine calculate_weight_matrix_lowdin_wrapper(cdft,tmb,input,ref_frags,calculate_overlap_matrix,meth_overlap)
  use module_base
  use module_types
  use module_fragments
  use constrained_dft
  implicit none
  integer, intent(in) :: meth_overlap
  type(cdft_data), intent(inout) :: cdft
  type(input_variables),intent(in) :: input
  type(dft_wavefunction), intent(inout) :: tmb
  logical, intent(in) :: calculate_overlap_matrix
  type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags

  integer :: nfrag_charged
  !real(kind=gp), dimension(:,:), pointer :: ovrlp_half
  character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'

  call f_routine('calculate_weight_matrix_lowdin_wrapper')

  ! wrapper here so we can modify charge and avoid restructuring the code as much whilst still using the routine elsewhere
  if ((.not. input%lin%calc_transfer_integrals) .and. (input%lin%cdft_orb(1)==0)) then
     cdft%charge=ref_frags(input%frag%frag_index(cdft%ifrag_charged(1)))%nelec-input%frag%charge(cdft%ifrag_charged(1))
  end if

  if (input%lin%calc_transfer_integrals) then
     nfrag_charged=2
  else
     nfrag_charged=1
  end if

  !ovrlp_half=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/), id='ovrlp_half')
  call calculate_weight_matrix_lowdin(cdft%weight_matrix,cdft%weight_matrix_,nfrag_charged,cdft%ifrag_charged,tmb,input%frag,&
       ref_frags,calculate_overlap_matrix,.true.,meth_overlap,input%lin%cdft_orb)
  !call f_free_ptr(ovrlp_half)

  call f_release_routine()

end subroutine calculate_weight_matrix_lowdin_wrapper


subroutine calculate_weight_matrix_lowdin(weight_matrix,weight_matrix_,nfrag_charged,ifrag_charged,tmb,input_frag,&
     ref_frags,calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap,cdft_orb)
  use module_base
  use module_types
  use module_fragments
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized
  use sparsematrix_base, only: matrices, sparse_matrix, sparsematrix_malloc_ptr, &
                               DENSE_FULL, assignment(=), &
                               allocate_matrices, deallocate_matrices
  use sparsematrix, only: compress_matrix, uncompress_matrix, &
                          gather_matrix_from_taskgroups_inplace, uncompress_matrix2
  use transposed_operations, only: calculate_overlap_transposed
  use matrix_operations, only: overlapPowerGeneral
  implicit none
  type(sparse_matrix), intent(in) :: weight_matrix
  type(matrices), intent(inout) :: weight_matrix_
  type(fragmentInputParameters),intent(in) :: input_frag
  type(dft_wavefunction), intent(inout) :: tmb
  logical, intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  integer, intent(in) :: nfrag_charged, meth_overlap
  integer, dimension(2), intent(in) :: ifrag_charged
  integer, dimension(2), intent(in) :: cdft_orb
  !local variables
  integer :: ifrag,iorb,ifrag_ref,isforb,ierr,cdft_iorb1,cdft_iorb2,jorb
  integer,dimension(1) :: power
  real(kind=gp), allocatable, dimension(:,:) :: proj_mat, proj_ovrlp_half, weight_matrixp
  character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'
  real(kind=gp) :: max_error, mean_error
  type(matrices),dimension(1) :: inv_ovrlp

  call f_routine(id='calculate_weight_matrix_lowdin')

  if (cdft_orb(1) == 0) call allocate_matrices(tmb%linmat%smat(3), allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))

  if (calculate_overlap_matrix) then
     if(.not.tmb%can_use_transposed) then
         !if(.not.associated(tmb%psit_c)) then
         !    tmb%psit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='tmb%psit_c')
         !end if
         !if(.not.associated(tmb%psit_f)) then
         !    tmb%psit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='tmb%psit_f')
         !end if
         call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
              TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
         tmb%can_use_transposed=.true.
     end if

     call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
          tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%smat(1), tmb%linmat%auxs, tmb%linmat%ovrlp_)
     !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%smat(1), tmb%linmat%ovrlp_)
  end if   

  if (calculate_ovrlp_half) then
     tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(1), iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
     call uncompress_matrix2(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
          tmb%linmat%smat(1), tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
     ! Maybe not clean here to use twice tmb%linmat%smat(1), but it should not
     ! matter as dense is used
     if (cdft_orb(1) == 0) then
        power(1)=2
        call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc,bigdft_mpi%mpi_comm,&
             meth_overlap, 1, power, &
             tmb%orthpar%blocksize_pdsyev, &
             imode=2, ovrlp_smat=tmb%linmat%smat(1), inv_ovrlp_smat=tmb%linmat%smat(3), &
             ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, &
             check_accur=.false., max_error=max_error, mean_error=mean_error)
             !check_accur=.true., max_error=max_error, mean_error=mean_error)
        call f_free_ptr(tmb%linmat%ovrlp_%matrix)
     end if
  end if

  proj_mat=f_malloc0((/tmb%orbs%norb,tmb%orbs%norb/),id='proj_mat')
  !call f_zero(tmb%orbs%norb**2,proj_mat(1,1))

  ! spatial constraint - assume that if the first fragment is a spatial constraint, the second one must be also
  if (cdft_orb(1) == 0) then
     isforb=0
     do ifrag=1,input_frag%nfrag
        ifrag_ref=input_frag%frag_index(ifrag)
        if (ifrag==ifrag_charged(1)) then
           do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
              proj_mat(iorb+isforb,iorb+isforb)=1.0_gp
           end do
        end if
        if (nfrag_charged==2) then
           if (ifrag==ifrag_charged(2)) then
              do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
                 proj_mat(iorb+isforb,iorb+isforb)=-1.0_gp
              end do
           end if
        end if
        isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
     end do

  ! orbital constraint
  else
     ! find the orbital index for fragment 1
     call get_cdft_iorb(ifrag_charged(1),cdft_orb(1),cdft_iorb1,input_frag,ref_frags)
     if (nfrag_charged==2) call get_cdft_iorb(ifrag_charged(2),cdft_orb(2),cdft_iorb2,input_frag,ref_frags)

     isforb=0
     do ifrag=1,input_frag%nfrag
        ifrag_ref=input_frag%frag_index(ifrag)
        if (ifrag==ifrag_charged(1)) then
           do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
              do jorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
                 proj_mat(iorb+isforb,jorb+isforb) = ref_frags(ifrag_ref)%coeff(iorb,cdft_iorb1) &
                      * ref_frags(ifrag_ref)%coeff(jorb,cdft_iorb1)
              end do
           end do
        end if
        if (nfrag_charged==2) then
           if (ifrag==ifrag_charged(2)) then
              do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
                 do jorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
                    proj_mat(iorb+isforb,jorb+isforb) = ref_frags(ifrag_ref)%coeff(iorb,cdft_iorb2) &
                         * ref_frags(ifrag_ref)%coeff(jorb,cdft_iorb2)
                 end do
              end do
           end if
        end if
        isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
     end do
  end if

  proj_ovrlp_half=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/),id='proj_ovrlp_half')
  if (tmb%orbs%norbp>0) then
     if (cdft_orb(1) == 0) then
        call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, &
            tmb%orbs%norb, 1.d0, &
            proj_mat(1,1), tmb%orbs%norb, &
            inv_ovrlp(1)%matrix(1,tmb%orbs%isorb+1,1), tmb%orbs%norb, 0.d0, &
            proj_ovrlp_half(1,1), tmb%orbs%norb)
     else
        call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, &
            tmb%orbs%norb, 1.d0, &
            proj_mat(1,1), tmb%orbs%norb, &
            tmb%linmat%ovrlp_%matrix(1,tmb%orbs%isorb+1,1), tmb%orbs%norb, 0.d0, &
            proj_ovrlp_half(1,1), tmb%orbs%norb)

     end if
  end if
  call f_free(proj_mat)
  weight_matrixp=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='weight_matrixp')
  if (tmb%orbs%norbp>0) then
     if (cdft_orb(1) == 0) then
        call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, & 
          tmb%orbs%norb, 1.d0, &
          inv_ovrlp(1)%matrix(1,1,1), tmb%orbs%norb, &
          proj_ovrlp_half(1,1), tmb%orbs%norb, 0.d0, &
          weight_matrixp(1,1), tmb%orbs%norb)
     else
        call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, & 
          tmb%orbs%norb, 1.d0, &
          tmb%linmat%ovrlp_%matrix(1,1,1), tmb%orbs%norb, &
          proj_ovrlp_half(1,1), tmb%orbs%norb, 0.d0, &
          weight_matrixp(1,1), tmb%orbs%norb)
     end if
  end if
  call f_free(proj_ovrlp_half)
  weight_matrix_%matrix = sparsematrix_malloc_ptr(weight_matrix,iaction=DENSE_FULL,id='weight_matrix_%matrix')
  if (bigdft_mpi%nproc>1) then
     call mpi_allgatherv(weight_matrixp, tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, weight_matrix_%matrix, &
          tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, &
          mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
      if (weight_matrix%nspin/=1) then
          stop 'NEED TO FIX THE SPIN HERE: calculate_weight_matrix_lowdin'
      end if
     call vcopy(tmb%orbs%norb*tmb%orbs%norb*weight_matrix%nspin,weight_matrixp(1,1),1,weight_matrix_%matrix(1,1,1),1)
  end if
  call f_free(weight_matrixp)
  call compress_matrix(bigdft_mpi%iproc,bigdft_mpi%nproc,weight_matrix,weight_matrix_%matrix,weight_matrix_%matrix_compr)
  call f_free_ptr(weight_matrix_%matrix)
  if (cdft_orb(1) == 0) call deallocate_matrices(inv_ovrlp(1))
  if (calculate_ovrlp_half .and. cdft_orb(1) /= 0) call f_free_ptr(tmb%linmat%ovrlp_%matrix)
  call f_release_routine()

end subroutine calculate_weight_matrix_lowdin

     subroutine get_cdft_iorb(ifrag,cdft_orb,iorb,input_frag,ref_frags)
        use module_fragments
        implicit none
        integer, intent(in) :: ifrag, cdft_orb
        integer, intent(out) :: iorb
        type(fragmentInputParameters),intent(in) :: input_frag
        type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
 
        integer :: ifrag_ref, ihomo
        
        ifrag_ref = input_frag%frag_index(ifrag)
        ihomo = ceiling(ref_frags(ifrag_ref)%nelec/2.0d0)
        ! HOMO, HOMO-1...
        if (cdft_orb < 0) then
           iorb = ihomo + cdft_orb + 1
        ! LUMO, LUMO+1...
        else
           iorb = ihomo + cdft_orb
        end if

     end subroutine get_cdft_iorb



!> CDFT: calculates the weight matrix w_ab given w(r)
!! for the moment putting densities in global box and ignoring parallelization
subroutine calculate_weight_matrix_using_density(iproc,nproc,cdft,tmb,at,input,GPU,denspot)
  use module_base
  use module_types
  use constrained_dft, only: cdft_data
  use module_interfaces, only: LocalHamiltonianApplication
  use module_fragments
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized, start_onesided_communication
  use sparsematrix_base, only : matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix, only: gather_matrix_from_taskgroups_inplace
  use transposed_operations, only: calculate_overlap_transposed
  use rhopotential, only: full_local_potential
  use locregs_init, only: small_to_large_locreg
  use locreg_operations, only: confpot_data
  implicit none
  integer,intent(in) :: iproc, nproc
  type(cdft_data), intent(inout) :: cdft
  type(atoms_data), intent(in) :: at
  type(input_variables),intent(in) :: input
  type(dft_wavefunction), intent(inout) :: tmb
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers),intent(inout) :: GPU

  !integer :: iorb, jorb
  integer :: jproc
  real(kind=gp),dimension(:),allocatable :: hpsit_c, hpsit_f
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  type(energy_terms) :: energs
  character(len=*),parameter :: subname='calculate_weight_matrix_using_density'
  type(matrices) :: weight_
  integer,dimension(:),allocatable :: n3p

  energs=energy_terms_null()
  call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
  n3p = f_malloc(0.to.nproc-1,id='n3p')
  do jproc=0,nproc-1
      n3p(jproc) = max(denspot%dpbox%nscatterarr(jproc,2),1)
  end do
  call start_onesided_communication(bigdft_mpi%iproc,bigdft_mpi%nproc,&
       denspot%dpbox%mesh%ndims(1),denspot%dpbox%mesh%ndims(2),n3p,cdft%weight_function, &
       tmb%ham_descr%comgp%nrecvbuf,tmb%ham_descr%comgp%recvbuf,tmb%ham_descr%comgp,tmb%ham_descr%lzd)
  call f_free(n3p)

  allocate(confdatarrtmp(tmb%orbs%norbp))
  call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)

  call small_to_large_locreg(bigdft_mpi%iproc, tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%isorb, tmb%orbs%inwhichlocreg, &
       tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
       tmb%psi, tmb%ham_descr%psi)

  if (tmb%ham_descr%npsidim_orbs > 0) call f_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))

  call full_local_potential(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%orbs,tmb%ham_descr%lzd,2, &
       denspot%dpbox,denspot%xc,cdft%weight_function,denspot%pot_work,tmb%ham_descr%comgp)
  call LocalHamiltonianApplication(bigdft_mpi%iproc,bigdft_mpi%nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
       tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,&
       energs,input%SIC,GPU,2,denspot%xc,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
       potential=cdft%weight_function,comgp=tmb%ham_descr%comgp)

  deallocate(confdatarrtmp)

  print *,'CDFT Ekin,Epot,Eproj,Eh,Exc,Evxc',energs%ekin,energs%epot,energs%eproj,energs%eh,energs%exc,energs%evxc

  call f_free_ptr(denspot%pot_work)


  ! calculate w_ab
  ! Calculate the matrix elements <phi|H|phi>.
  if(.not.tmb%ham_descr%can_use_transposed) then
      !!if(associated(tmb%ham_descr%psit_c)) then
      !!    call f_free_ptr(tmb%ham_descr%psit_c)
      !!end if
      !!if(associated(tmb%ham_descr%psit_f)) then
      !!    call f_free_ptr(tmb%ham_descr%psit_f)
      !!end if

      !!tmb%ham_descr%psit_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='tmb%ham_descr%psit_c')
      !!tmb%ham_descr%psit_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='tmb%ham_descr%psit_f')
      call transpose_localized(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs, &
           tmb%ham_descr%collcom,TRANSPOSE_FULL,tmb%ham_descr%psi,tmb%ham_descr%psit_c,tmb%ham_descr%psit_f,tmb%ham_descr%lzd)
      tmb%ham_descr%can_use_transposed=.true.
  end if

  hpsit_c = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c')
  hpsit_f = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f')
  call transpose_localized(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,  &
       tmb%ham_descr%collcom,TRANSPOSE_FULL,tmb%hpsi,hpsit_c,hpsit_f,tmb%ham_descr%lzd)

  weight_ = matrices_null()
  call allocate_matrices(tmb%linmat%smat(2), allocate_full=.false., &
       matname='weight_', mat=weight_)

  call calculate_overlap_transposed(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%orbs,tmb%ham_descr%collcom, &
       tmb%ham_descr%psit_c,hpsit_c,tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%smat(2), tmb%linmat%auxm, weight_)
  call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc,bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
       tmb%linmat%smat(2), weight_)
  ! This can then be deleted if the transition to the new type has been completed.
  cdft%weight_matrix_%matrix_compr=weight_%matrix_compr
  call deallocate_matrices(weight_)

  call f_free(hpsit_c)
  call f_free(hpsit_f)

  ! debug
  !allocate(cdft%weight_matrix%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  !call memocc(istat, cdft%weight_matrix%matrix, 'cdft%weight_matrix%matrix', subname)
  !call uncompress_matrix(bigdft_mpi%iproc,cdft%weight_matrix)
  !do iorb=1,tmb%orbs%norb
  !   do jorb=1,tmb%orbs%norb
  !      write(87,*) iorb,jorb,cdft%weight_matrix%matrix(iorb,jorb)
  !   end do
  !end do
  !deallocate(cdft%weight_matrix%matrix, stat=istat)
  !call memocc(istat, iall, 'cdft%weight_matrix%matrix', subname)
  ! end debug

end subroutine calculate_weight_matrix_using_density


!> CDFT: calculates the weight function w(r)
!! for the moment putting densities in global box and ignoring parallelization
subroutine calculate_weight_function(in,ref_frags,cdft,ndimrho_all_fragments,rho_all_fragments,tmb,atoms,rxyz,denspot)
  use module_base
  use module_types
  use module_fragments
  use constrained_dft
  use io, only: plot_density
  implicit none
  type(input_variables), intent(in) :: in
  type(system_fragment), dimension(in%frag%nfrag_ref), intent(inout) :: ref_frags
  type(cdft_data), intent(inout) :: cdft
  integer, intent(in) :: ndimrho_all_fragments
  real(kind=wp),dimension(ndimrho_all_fragments),intent(in) :: rho_all_fragments
  type(dft_wavefunction), intent(in) :: tmb
  type(atoms_data), intent(in) :: atoms ! just for plotting
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz ! just for plotting
  type(DFT_local_fields), intent(in) :: denspot ! just for plotting

  integer :: ifrag, ifrag_charged, ref_frag_charged, iorb_start, ipt

  ! unecessary recalculation here, but needs generalizing anyway
  iorb_start=1
  do ifrag=1,in%frag%nfrag
     if (in%frag%charge(ifrag)/=0) then
        ifrag_charged=ifrag
        exit
     end if
     iorb_start=iorb_start+ref_frags(in%frag%frag_index(ifrag))%fbasis%forbs%norb
  end do
  ref_frag_charged=in%frag%frag_index(ifrag_charged)

  ! check both densities are the same size (for now calculating fragment density in global box)
  if (ndimrho_all_fragments /= cdft%ndim_dens) stop 'Fragment density dimension does not match global density dimension'

  ! only need to calculate the fragment density for the fragment with a charge (stored in frag%fbasis%density)
  ! ideally would use already calculated kernel here, but charges could vary so would need an actual fragment
  ! rather than a reference fragment - could point to original fragment but modify nks...
  ! also converting input charge to integer which might not be the case
  call calculate_fragment_density(ref_frags(ref_frag_charged),ndimrho_all_fragments,tmb,iorb_start,&
       int(in%frag%charge(ifrag_charged)),atoms,rxyz,denspot)

  ! the weight function becomes the fragment density of ifrag_charged divided by the sum of all fragment densities
  ! using a cutoff of 1.e-6 for fragment density and 1.e-12 for denominator
  do ipt=1,ndimrho_all_fragments
      if (denspot%rhov(ipt) >= 1.0e-12) then
         cdft%weight_function(ipt)=ref_frags(ref_frag_charged)%fbasis%density(ipt)/denspot%rhov(ipt)
      else
         cdft%weight_function(ipt)=0.0_gp
      end if
  end do

  ! COME BACK TO THIS
  cdft%charge=ref_frags(ref_frag_charged)%nelec-in%frag%charge(ref_frag_charged)

  ! plot the weight function, fragment density and initial total density (the denominator) to check

  call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'fragment_density.cube', &
       atoms,rxyz,denspot%pkernel,1,ref_frags(ref_frag_charged)%fbasis%density)

  call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'weight_function.cube', &
       atoms,rxyz,denspot%pkernel,1,cdft%weight_function)

  call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'initial_density.cube', &
       atoms,rxyz,denspot%pkernel,1,denspot%rhov)

  ! deallocate fragment density here as we don't need it any more
  if (associated(ref_frags(ref_frag_charged)%fbasis%density)) call f_free_ptr(ref_frags(ref_frag_charged)%fbasis%density)

end subroutine calculate_weight_function
