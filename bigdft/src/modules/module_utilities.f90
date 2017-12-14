module module_utilities

  private

  public :: calculate_fragment_multipoles_fn

  contains

    subroutine calculate_fragment_multipoles_fn(iproc, nproc, comm, &
               matrix_format, metadata_file, fragment_file, &
               overlap_file, kernel_file, kernel_matmul_file, lmax, multipoles_files)
           use futile
           use wrapper_mpi
           use wrapper_linalg
           use module_input_dicts, only: merge_input_file_to_dict
           use sparsematrix_base, only: sparse_matrix, matrices, matrices_null, assignment(=), &
                                        SPARSE_FULL, DENSE_FULL, DENSE_PARALLEL, SPARSE_TASKGROUP, &
                                        sparsematrix_malloc_ptr, sparsematrix_malloc0_ptr,&
                                        deallocate_sparse_matrix, deallocate_matrices, &
                                        sparse_matrix_metadata, deallocate_sparse_matrix_metadata
           use sparsematrix_init, only: bigdft_to_sparsebigdft, distribute_columns_on_processes_simple, &
                                        write_sparsematrix_info, init_matrix_taskgroups_wrapper
           use sparsematrix_io, only: read_sparse_matrix, write_sparse_matrix, read_linear_coefficients
           use sparsematrix, only: uncompress_matrix, uncompress_matrix_distributed2, resize_matrix_to_taskgroup, &
                                   transform_sparse_matrix, matrix_matrix_mult_wrapper
           use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_bigdft, &
                                             sparse_matrix_and_matrices_init_from_file_ccs, &
                                             sparse_matrix_metadata_init_from_file, &
                                             matrices_init_from_file_bigdft, &
                                             ccs_data_from_sparse_matrix, &
                                             ccs_matrix_write, &
                                             matrices_init, &
                                             get_selected_eigenvalues_from_FOE, &
                                             matrix_chebyshev_expansion
           use multipole, only: init_extract_matrix_lookup, extract_matrix, correct_multipole_origin
           implicit none
    
           ! Calling arguments
           integer,intent(in) :: iproc, nproc, comm, lmax
           character(len=*),intent(in) :: matrix_format, metadata_file, fragment_file
           character(len=*),intent(in) :: overlap_file, kernel_file, kernel_matmul_file
           character(len=*),dimension(-lmax:lmax,0:lmax),intent(in) :: multipoles_files
    
           ! Local variables
           external :: gather_timings
           type(sparse_matrix_metadata) :: smmd
           integer,dimension(:,:),allocatable :: fragment_atom_id, fragment_atom_type
           integer,dimension(:),allocatable :: fragment_nat
           logical,dimension(-lmax:lmax) :: file_present
           integer,dimension(:,:),allocatable :: lookup_kernel
           logical :: file_exists, perx, pery, perz
           integer :: ispin, iat, iiat, nat_frag, iat_frag, isf, mm, iitype, ifrag
           integer :: ist1, ist2, nmpmat, nsf_frag, ist, l, ll, m, nfrag, nat_frag_max
           type(matrices) :: ovrlp_mat, kernel_mat,  ovrlp_large, kernel_eff
           type(matrices),dimension(1) :: inv_ovrlp
           type(matrices),dimension(:,:),allocatable :: multipoles_matrices
           type(sparse_matrix),dimension(2) :: smat
           type(dictionary), pointer :: dict_timing_info
           character(len=32),dimension(:),allocatable :: fragment_atom_name
           real(kind=8) :: tr, spin_factor
           logical,dimension(:),allocatable :: supfun_in_fragment
           type(dictionary),pointer :: fragment_dict, fragment_item, fragment_atom
           type(mpi_environment) :: mpi_env
           real(kind=8),dimension(:,:,:),allocatable :: mpmat, kqmat
           real(kind=8),dimension(:,:),allocatable :: kmat, fragment_multipoles
           real(kind=8),dimension(:),allocatable :: fragment_charge
           real(kind=8),dimension(3) :: fragment_center

           integer :: i
    
           mpi_env = mpi_environment_null()
           call mpi_environment_set(mpi_env, iproc, nproc, comm, nproc)
           inquire(file=trim(fragment_file), exist=file_exists)
           if (file_exists) then
               fragment_dict => null()
               call merge_input_file_to_dict(fragment_dict, trim(fragment_file), mpi_env)
           else
               call f_err_throw('file'//trim(fragment_file)//'not present')
           end if
    
           call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
           if (iproc==0) then
               call yaml_mapping_open('Atomic System Properties')
               call yaml_map('Types of atoms',smmd%atomnames)
               call yaml_mapping_close()
           end if
    
           call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
                iproc, nproc, comm, smat(1), ovrlp_mat, &
                init_matmul=.false.)
    
           call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(kernel_file), &
                iproc, nproc, comm, smat(2), kernel_mat, &
                init_matmul=.true., filename_mult=trim(kernel_matmul_file))
    
           call init_matrix_taskgroups_wrapper(iproc, nproc, comm, .false., 2, smat)
    
           call resize_matrix_to_taskgroup(smat(1), ovrlp_mat)
           call resize_matrix_to_taskgroup(smat(2), kernel_mat)
    
    
           if (iproc==0) then
               call yaml_mapping_open('Matrix properties')
               call write_sparsematrix_info(smat(1), 'Overlap matrix')
               call write_sparsematrix_info(smat(2), 'Density kernel')
               call yaml_mapping_close()
           end if
    
           ! Check which multipole matrices are present
           ll = -1
           do l=0,lmax
               file_present(:) = .false.
               do m=-l,l
                   inquire(file=trim(multipoles_files(m,l)), exist=file_exists)
                   if (file_exists) then
                       file_present(m) = .true.
                   end if
               end do
               if (any(file_present(-l:l))) then
                   if (.not.all(file_present(-l:l))) then
                       call f_err_throw('for a given shell all matrices must be present', &
                            err_name='BIGDFT_RUNTIME_ERROR')
                   end if
                   ll = l
               end if
           end do
    
           allocate(multipoles_matrices(-ll:ll,0:ll))
           do l=0,ll
               do m=-l,l
                   call matrices_init_from_file_bigdft(matrix_format, trim(multipoles_files(m,l)), &
                        iproc, nproc, comm, smat(1), multipoles_matrices(m,l))
                   call resize_matrix_to_taskgroup(smat(1), multipoles_matrices(m,l))
               end do
           end do
    
           call timing(comm,'INIT','PR')
    
           mpi_env = mpi_environment_null()
           call mpi_environment_set(mpi_env, iproc, nproc, comm, nproc)
           inquire(file=trim(fragment_file), exist=file_exists)
           if (file_exists) then
               call merge_input_file_to_dict(fragment_dict, trim(fragment_file), mpi_env)
           else
               call f_err_throw('file '//trim(fragment_file)//' not present')
           end if
    
           supfun_in_fragment = f_malloc(smat(1)%nfvctr,id='supfun_in_fragment')
    
           ! Calculate the matrix KS
           kernel_eff = matrices_null()
           kernel_eff%matrix_compr = sparsematrix_malloc0_ptr(smat(2),iaction=SPARSE_TASKGROUP,id='kernel_eff%matrix_compr')
           ovrlp_large = matrices_null()
           ovrlp_large%matrix_compr = sparsematrix_malloc_ptr(smat(2), SPARSE_TASKGROUP, id='ovrlp_large%matrix_compr')
           call transform_sparse_matrix(iproc, smat(1), smat(2), SPARSE_TASKGROUP, 'small_to_large', &
                smat_in=ovrlp_mat%matrix_compr, lmat_out=ovrlp_large%matrix_compr)
           ! Should use the highlevel wrapper...
           do ispin=1,smat(2)%nspin
               ist=(ispin-1)*smat(2)%nvctrp_tg+1
               call matrix_matrix_mult_wrapper(iproc, nproc, smat(2), &
                    kernel_mat%matrix_compr(ist:), ovrlp_large%matrix_compr(ist:), kernel_eff%matrix_compr(ist:))
           end do
    
           ! Calculate the matrix S^-1
           inv_ovrlp(1) = matrices_null()
           inv_ovrlp(1)%matrix_compr = sparsematrix_malloc_ptr(smat(2), SPARSE_TASKGROUP, id='inv_ovrlp%matrix_compr')
           call matrix_chebyshev_expansion(iproc, nproc, comm, 1, [-1.d0], &
                smat(1), smat(2), ovrlp_mat, inv_ovrlp)
    
           ! Calculate the matrix S^-1*P, where P are the multipole matrices
           do l=0,ll
               do m=-l,l
                   call transform_sparse_matrix(iproc, smat(1), smat(2), SPARSE_TASKGROUP, 'small_to_large', &
                        smat_in=multipoles_matrices(m,l)%matrix_compr, lmat_out=ovrlp_large%matrix_compr)
                   call f_free_ptr(multipoles_matrices(m,l)%matrix_compr)
                   multipoles_matrices(m,l)%matrix_compr = sparsematrix_malloc_ptr(smat(2),iaction=SPARSE_TASKGROUP,&
                                                                   id='matrix_compr')
                   ! Should use the highlevel wrapper...
                   do ispin=1,smat(2)%nspin
                       ist=(ispin-1)*smat(2)%nvctrp_tg+1
                       call matrix_matrix_mult_wrapper(iproc, nproc, smat(2), &
                            inv_ovrlp(1)%matrix_compr(ist:), ovrlp_large%matrix_compr(ist:), &
                            multipoles_matrices(m,l)%matrix_compr(ist:))
                   end do
               end do
           end do
    
           call deallocate_matrices(ovrlp_large)
           call deallocate_matrices(inv_ovrlp(1))
    
    
    
           ! Periodicity in the three direction
           perx=(smmd%geocode /= 'F')
           pery=(smmd%geocode == 'P')
           perz=(smmd%geocode /= 'F')
    
    

           ! The number of fragments
           nfrag = dict_len(fragment_dict)

           ! The maximal number of atoms in a fragment
           nat_frag_max = 0
           fragment_item => null()
           do while (iterating(fragment_item,fragment_dict))
               nat_frag = dict_len(fragment_item)
               if (nat_frag>nat_frag_max) then
                   nat_frag_max = nat_frag
               end if
           end do

           ! Get the number of multipoles to be calculated
           nmpmat = 0
           do l=0,ll
               nmpmat = nmpmat + 2*l + 1
           end do

           fragment_atom_id = f_malloc0([nat_frag_max,nfrag],id='fragment_atom_id')
           fragment_atom_name = f_malloc_str(len(fragment_atom_name),nat_frag_max,id='fragment_atom_name')
           fragment_atom_type = f_malloc0([nat_frag_max,nfrag],id='fragment_atom_type')
           fragment_charge = f_malloc0(nfrag,id='fragment_charge')
           fragment_multipoles = f_malloc0([0.to.nmpmat,1.to.nfrag],id='fragment_multipoles')
           fragment_nat = f_malloc0(nfrag,id='fragment_nat')
    
           fragment_item => null()
           fragment_loop: do while (iterating(fragment_item,fragment_dict))

               ifrag = fragment_item%data%item + 1

               if (mod(ifrag-1,nproc)/=iproc) then
                   cycle fragment_loop
               end if

               nat_frag = dict_len(fragment_item)
               fragment_nat(ifrag) = nat_frag

               fragment_atom => dict_iter(fragment_item)
               iat_frag = 0
               fragment_charge(ifrag) = 0.d0
               fragment_center(1:3) = 0.d0
               do while (associated(fragment_atom))
                   iat_frag = iat_frag + 1
                   iiat = fragment_atom
                   iitype = smmd%iatype(iiat)
                   fragment_atom_id(iat_frag,ifrag) = iiat
                   !fragment_atom_name(iat_frag,ifrag) = trim(smmd%atomnames(iitype))
                   fragment_atom_type(iat_frag,ifrag) = iitype
                   fragment_charge(ifrag) = fragment_charge(ifrag) + smmd%nelpsp(iitype)
                   fragment_center(1:3) = fragment_center(1:3) + smmd%rxyz(1:3,iiat)*smmd%nelpsp(iitype)
                   fragment_atom => dict_next(fragment_atom)
               end do
               fragment_center(1:3) = fragment_center(1:3)/fragment_charge(ifrag)
    
               ! Count how many support functions belong to the fragment
               supfun_in_fragment(:) = .false.
               nsf_frag = 0
               outerloop: do isf=1,smat(1)%nfvctr
                   iiat = smmd%on_which_atom(isf)
                   do iat=1,nat_frag
                   if (iiat==fragment_atom_id(iat,ifrag)) then
                           supfun_in_fragment(isf) = .true.
                           nsf_frag = nsf_frag + 1
                           cycle outerloop
                       end if
                   end do
               end do outerloop
    
    
               lookup_kernel = f_malloc([nsf_frag,nsf_frag],id='lookup_kernel')
               call init_extract_matrix_lookup(smat(2), supfun_in_fragment, nsf_frag, lookup_kernel)
               kqmat = f_malloc([nsf_frag,nsf_frag,2],id='kqmat')
               kmat = f_malloc0([nsf_frag,nsf_frag],id='kmat')
               mpmat = f_malloc0([nsf_frag,nsf_frag,nmpmat],id='mpmat')
    
               do ispin=1,smat(1)%nspin
                   !!ist1 = (ispin-1)*smat(1)%nvctrp_tg+1
                   ist2 = (ispin-1)*smat(2)%nvctrp_tg+1
    
                   !call extract_matrix(smat(2), kernel_mat%matrix_compr(ist2:), nsf_frag, lookup_kernel, kmat(1,1))
                   call extract_matrix(smat(2), kernel_eff%matrix_compr(ist2:), nsf_frag, lookup_kernel, kmat(1,1))
    
                   mm = 0
                   do l=0,ll
                       do m=-l,l
                           mm = mm + 1
                           call extract_matrix(smat(2), multipoles_matrices(m,l)%matrix_compr(ist2:), &
                                nsf_frag, lookup_kernel, mpmat(1,1,mm))
                           call correct_multipole_origin(ispin, smmd%nat, l, m, nsf_frag, &
                                smmd, smat(2), fragment_center, smmd%rxyz, supfun_in_fragment, &
                                perx, pery, perz, smmd%cell_dim, &
                                multipoles_matrices(-1:1,0:1), lookup_kernel, mpmat(1,1,mm))
                       end do
                   end do
    
                   ! Now calculate KQ, where K is the kernel and Q the multipole matrix.
                   ! Take the trace of KQ, which is the multipole moment of the fragment.
                   ! Then calculate (KQ)^2 - KQ and take the trace, which gives the purity indicator.
                   mm = 0
                   do l=0,ll
                       do m=-l,l
                           mm = mm + 1
                           call gemm('n', 'n', nsf_frag, nsf_frag, nsf_frag, 1.d0, kmat(1,1), nsf_frag, &
                                mpmat(1,1,mm), nsf_frag, 0.d0, kqmat(1,1,1), nsf_frag)
                           tr = 0.d0
                           do isf=1,nsf_frag
                               tr = tr + kqmat(isf,isf,1)
                           end do
                           fragment_multipoles(mm,ifrag) = fragment_multipoles(mm,ifrag) + tr
                           if (l==0) then
                               if (smat(1)%nspin==1) then
                                   spin_factor = 2.0d0
                               else
                                   spin_factor = 1.0d0
                               end if
                               call dscal(nsf_frag**2, 1.d0/spin_factor, kqmat(1,1,1), 1)
                               call gemm('n', 'n', nsf_frag, nsf_frag, nsf_frag, 1.d0, kqmat(1,1,1), nsf_frag, &
                                    kqmat(1,1,1), nsf_frag, 0.d0, kqmat(1,1,2), nsf_frag)
                               call axpy(nsf_frag**2, -1.d0, kqmat(1,1,1), 1, kqmat(1,1,2), 1)
                               tr = 0.d0
                               do isf=1,nsf_frag
                                   tr = tr + kqmat(isf,isf,2)
                               end do
                               fragment_multipoles(0,ifrag) = fragment_multipoles(0,ifrag) + tr*spin_factor
                           end if
                       end do
                   end do
                
               end do
               fragment_multipoles(0,ifrag) = fragment_multipoles(0,ifrag)/fragment_charge(ifrag)
    
               !!if (iproc==0) then
               !!    call yaml_sequence(advance='no')
               !!    call yaml_map('Atom IDs',fragment_atom_id)
               !!    call yaml_map('Atom names',fragment_atom_name)
               !!    call yaml_map('Neutral fragment charge',fragment_charge(ifrag),fmt='(1es13.6)')
               !!    call yaml_map('Purity indicator',fragment_multipoles(0,ifrag),fmt='(1es13.6)')
               !!    mm = 0
               !!    do l=0,ll
               !!    call yaml_map('q'//adjustl(trim(yaml_toa(l))),fragment_multipoles(mm+1:mm+2*l+1,ifrag),fmt='(1es13.6)')
               !!        mm = mm + 2*l+1
               !!    end do
               !!end if
    
               call f_free(mpmat)
               call f_free(kmat)
               call f_free(kqmat)
               call f_free(lookup_kernel)
    
           end do fragment_loop


           ! Gather together the results
           call fmpi_allreduce(fragment_atom_type, FMPI_SUM, comm=comm)
           call fmpi_allreduce(fragment_atom_id, FMPI_SUM, comm=comm)
           call fmpi_allreduce(fragment_charge, FMPI_SUM, comm=comm)
           call fmpi_allreduce(fragment_multipoles, FMPI_SUM, comm=comm)
           call fmpi_allreduce(fragment_nat, FMPI_SUM, comm=comm)


           if (iproc==0) then
               call yaml_sequence_open('Fragment multipoles')
               fragment_item => null()
               do while (iterating(fragment_item,fragment_dict))
                   call yaml_sequence(advance='no')
                   ifrag = fragment_item%data%item + 1
                   nat_frag = fragment_nat(ifrag)
                   do iat=1,nat_frag
                       iitype = fragment_atom_type(iat,ifrag)
                       fragment_atom_name(iat) =  trim(smmd%atomnames(iitype))
                   end do
                   call yaml_map('Atom IDs',fragment_atom_id(1:nat_frag,ifrag))
                   call yaml_map('Atom names',fragment_atom_name(1:nat_frag))
                   call yaml_map('Neutral fragment charge',fragment_charge(ifrag),fmt='(1es13.6)')
                   call yaml_map('Purity indicator',fragment_multipoles(0,ifrag),fmt='(1es13.6)')
                   mm = 0
                   do l=0,ll
                   call yaml_map('q'//adjustl(trim(yaml_toa(l))),fragment_multipoles(mm+1:mm+2*l+1,ifrag),fmt='(1es13.6)')
                       mm = mm + 2*l+1
                   end do
               end do
               call yaml_sequence_close()
           end if
           call f_free(supfun_in_fragment)
    
    
           call f_free(fragment_charge)
           call f_free(fragment_nat)
           call f_free(fragment_multipoles)
           call f_free_str(len(fragment_atom_name),fragment_atom_name)
           call f_free(fragment_atom_type)
           call f_free(fragment_atom_id)
           call deallocate_sparse_matrix_metadata(smmd)
           call deallocate_sparse_matrix(smat(1))
           call deallocate_sparse_matrix(smat(2))
           call deallocate_matrices(ovrlp_mat)
           call deallocate_matrices(kernel_mat)
           call deallocate_matrices(kernel_eff)
           call dict_free(fragment_dict)
           do l=0,ll
               do m=-l,l
                   call deallocate_matrices(multipoles_matrices(m,l))
               end do
           end do
           deallocate(multipoles_matrices)
    
           call f_timing_checkpoint(ctr_name='LAST',mpi_comm=comm,nproc=mpisize(),&
                gather_routine=gather_timings)
    
    end subroutine calculate_fragment_multipoles_fn

end module module_utilities
