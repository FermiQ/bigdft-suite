!> @file
!!   Program to perform pre-/post-processing.
!! @author
!!   Copyright (C) 2015-2016 BigDFT group (SM)
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS


!> Test the input files and estimates the memory occupation versus the number
!! of processors
program utilities

   use module_base
   use yaml_output
   use module_types, only: local_zone_descriptors, local_zone_descriptors_null
   use module_types, only: bigdft_init_errors, bigdft_init_timing_categories, orbitals_data
   use public_enums, only: LINEAR_PARTITION_NONE
   use module_atoms, only: atoms_data, atoms_data_null, deallocate_atoms_data
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
   use postprocessing_linear, only: CHARGE_ANALYSIS_LOEWDIN, CHARGE_ANALYSIS_MULLIKEN, build_ks_orbitals_postprocessing
   use multipole, only: multipole_analysis_driver_new, init_extract_matrix_lookup, extract_matrix, correct_multipole_origin
   use multipole_base, only: lmax
   use io, only: io_read_descr_linear, read_psi_compress
   use bigdft_run, only: bigdft_init
   use locregs, only: nullify_locreg_descriptors
   ! It's very bad practice to oblige the use of module_interfaces.
   ! If it is not done, it crashes with segfault due to the optional arguments
   use module_interfaces, only: orbitals_descriptors, open_filename_of_iorb
   use module_input_dicts, only: merge_input_file_to_dict
   implicit none
   external :: gather_timings
   character(len=*), parameter :: subname='utilities'
   character(len=1) :: geocode
   character(len=3) :: do_ortho
   character(len=30) :: tatonam, radical, colorname, linestart, lineend, cname, methodc
   character(len=128) :: method_name, overlap_file, hamiltonian_file, kernel_file, coeff_file, pdos_file, metadata_file
   character(len=128) :: line, cc, output_pdos, conversion, infile, outfile, iev_min_, iev_max_, fscale_, matrix_basis
   character(len=128) :: kernel_matmul_file, istart_ks_, iend_ks_, matrix_format, fragment_file
   character(len=128),dimension(-lmax:lmax,0:lmax) :: multipoles_files
   logical :: multipole_analysis = .false.
   logical :: solve_eigensystem = .false.
   logical :: calculate_pdos = .false.
   logical :: convert_matrix_format = .false.
   logical :: calculate_selected_eigenvalues = .false.
   logical :: build_KS_orbitals = .false.
   logical :: calculate_fragment_multipoles = .false.
   type(atoms_data) :: at
   type(sparse_matrix_metadata) :: smmd
   integer :: istat, i_arg, ierr, nspin, icount, nthread, method, ntypes
   integer :: nfvctr_m, nseg_m, nvctr_m
   integer :: nfvctr_l, nseg_l, nvctr_l
   !integer :: nfvctrp_l, isfvctr_l, nfvctrp_m, isfvctr_m, nfvctrp_s, isfvctr_s
   integer :: iconv, iev, iev_min, iev_max, l, ll, m
   integer,dimension(:),pointer :: on_which_atom
   integer,dimension(:),pointer :: keyv_s, keyv_m, keyv_l, on_which_atom_s, on_which_atom_m, on_which_atom_l
   integer,dimension(:),pointer :: iatype, nzatom, nelpsp
   integer,dimension(:),pointer :: col_ptr, row_ind
   integer,dimension(:,:,:),pointer :: keyg_s, keyg_m, keyg_l
   integer,dimension(:),allocatable :: in_which_locreg, fragment_atom_id
   logical,dimension(-lmax:lmax) :: file_present
   real(kind=8),dimension(:),pointer :: matrix_compr, eval_ptr
   real(kind=8),dimension(:,:),pointer :: rxyz, coeff_ptr
   real(kind=8),dimension(:),allocatable :: eval, energy_arr, occups, phi, phi_tmp
   real(kind=8),dimension(:,:),allocatable :: denskernel, pdos, occup_arr
   logical,dimension(:,:),allocatable :: calc_array
   integer,dimension(:,:),allocatable :: lookup_kernel
   logical :: file_exists, found
   type(matrices) :: ovrlp_mat, hamiltonian_mat, kernel_mat, mat, ovrlp_large, kernel_eff
   type(matrices),dimension(1) :: ovrlp_minus_one_half, inv_ovrlp
   type(matrices),dimension(:,:),allocatable :: multipoles_matrices
   type(sparse_matrix) :: smat_s, smat_m, smat_l
   type(sparse_matrix),dimension(2) :: smat
   type(dictionary), pointer :: dict_timing_info
   integer :: iunit, nat, iat, iat_prev, ii, iitype, iorb, itmb, itype, ival, ios, ipdos, ispin
   integer :: jtmb, norbks, npdos, npt, ntmb, jjtmb, istart_ks, iend_ks, nat_frag, iat_frag
   integer :: isf, iiat, nsf_frag, mm, nmpmat, ist1, ist2
   character(len=20),dimension(:),pointer :: atomnames
   character(len=30),dimension(:),allocatable :: pdos_name
   character(len=32),dimension(:),allocatable :: fragment_atom_name
   real(kind=8),dimension(3) :: cell_dim
   character(len=2) :: backslash, num
   real(kind=8) :: energy, occup, occup_pdos, total_occup, fscale
   type(f_progress_bar) :: bar
   type(orbitals_data) :: orbs
   integer,parameter :: nkpt = 1
   real(kind=8),dimension(3,nkpt),parameter :: kpt = reshape((/0.d0,0.d0,0.d0/),(/3,nkpt/))
   real(kind=8),dimension(nkpt),parameter :: wkpt = (/1.d0/)
   type(local_zone_descriptors) :: lzd
   integer :: confpotorder, ilr, iiorb, iiorb_out, ispinor, nspinor, onwhichatom_tmp, npsidim_orbs, nsize, ist
   real(kind=8) :: confpotprefac, eval_tmp, tr, fragment_charge
   character(len=256) :: error, filename, KSgrid_file
   logical :: lstat
   logical,dimension(:),allocatable :: supfun_in_fragment
   integer,parameter :: ncolors = 12
   type(dictionary),pointer :: fragment_dict, fragment_item, fragment_atom
   type(mpi_environment) :: mpi_env
   real(kind=8),dimension(:,:,:),allocatable :: mpmat, kqmat
   real(kind=8),dimension(:,:),allocatable :: kmat
   real(kind=8),dimension(:),allocatable :: fragment_multipoles
   real(kind=8),dimension(3) :: fragment_center
   logical :: perx, pery, perz
   ! Presumably well suited colorschemes from colorbrewer2.org
   character(len=20),dimension(ncolors),parameter :: colors=(/'#a6cee3', &
                                                              '#1f78b4', &
                                                              '#b2df8a', &
                                                              '#33a02c', &
                                                              '#fb9a99', &
                                                              '#e31a1c', &
                                                              '#fdbf6f', &
                                                              '#ff7f00', &
                                                              '#cab2d6', &
                                                              '#6a3d9a', &
                                                              '#ffff99', &
                                                              '#b15928'/)
   !$ integer :: omp_get_max_threads

   call f_lib_initialize()

   ! Initialize MPI
   !call bigdft_mpi_init(ierr)
   call bigdft_init()

    if (bigdft_mpi%iproc==0) then
        call yaml_new_document()
        call print_logo()
    end if

   !Time initialization
   call f_timing_reset(filename='time.yaml',master=(bigdft_mpi%iproc==0),verbose_mode=.false.)

   if (bigdft_mpi%iproc==0) then
       call yaml_scalar('',hfill='~')
       call yaml_scalar('BIGDFT PRE-/POST-PROCESSING',hfill='~')
   end if

   if (bigdft_mpi%iproc==0) then
       call yaml_mapping_open('Parallel environment')
       call yaml_map('MPI tasks',bigdft_mpi%nproc)
       nthread = 1
       !$ nthread = omp_get_max_threads()
       call yaml_map('OpenMP threads',nthread)
       call yaml_mapping_close()
   end if

   ! Get arguments
   call get_command_argument(1, value = tatonam, status = istat)

   write(radical, "(A)") "input"
   if(trim(tatonam)=='' .or. istat>0) then
      write(*,'(1x,a)')&
         &   'Usage: ./utilities -a [option]'
      write(*,'(1x,a)')&
         &   '[option] can be the following: '
      write(*,'(1x,a)')&
           &   '"multipoles-analysis"" ' 
      write(*,'(1x,a)')&
           & 'perform a charge analysis (Loewdin or Mulliken)'

      stop
   else
      i_arg = 1
      loop_getargs: do
         call get_command_argument(i_arg, value = tatonam, status = istat)
         !call getarg(i_arg,tatonam)
         if(trim(tatonam)=='' .or. istat > 0) then
            exit loop_getargs
         else if (trim(tatonam)=='multipole-analysis') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = method_name)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_basis)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = hamiltonian_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_matmul_file)
            do l=0,lmax
                do m=-l,l
                    i_arg = i_arg + 1
                    call get_command_argument(i_arg, value = multipoles_files(m,l))
                end do
            end do
            !write(*,'(1x,2a)')&
            !   &   'perform a Loewdin charge analysis'
            multipole_analysis = .true.
            exit loop_getargs
         !!else if (trim(tatonam)=='solve-eigensystem') then
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = metadata_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = hamiltonian_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = overlap_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = coeff_file)
         !!   !write(*,'(1x,2a)')&
         !!   !   &   'perform a Loewdin charge analysis'
         !!   solve_eigensystem = .true.
         !!   exit loop_getargs
         !!else if (trim(tatonam)=='pdos') then
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = metadata_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = coeff_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = hamiltonian_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = overlap_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = pdos_file)
         !!   calculate_pdos = .true.
         !!else if (trim(tatonam)=='convert-matrix-format') then
         !!    i_arg = i_arg + 1
         !!    call get_command_argument(i_arg, value = conversion)
         !!    i_arg = i_arg + 1
         !!    call get_command_argument(i_arg, value = infile)
         !!    i_arg = i_arg + 1
         !!    call get_command_argument(i_arg, value = outfile)
         !!    convert_matrix_format = .true.
         !!else if (trim(tatonam)=='calculate-selected-eigenvalues') then
         !!    i_arg = i_arg + 1
         !!    call get_command_argument(i_arg, value = metadata_file)
         !!    i_arg = i_arg + 1
         !!    call get_command_argument(i_arg, value = overlap_file)
         !!    i_arg = i_arg + 1
         !!    call get_command_argument(i_arg, value = hamiltonian_file)
         !!    i_arg = i_arg + 1
         !!    call get_command_argument(i_arg, value = kernel_file)
         !!    i_arg = i_arg + 1
         !!    call get_command_argument(i_arg, value = iev_min_)
         !!    read(iev_min_,fmt=*,iostat=ierr) iev_min
         !!    i_arg = i_arg + 1
         !!    call get_command_argument(i_arg, value = iev_max_)
         !!    read(iev_max_,fmt=*,iostat=ierr) iev_max
         !!    i_arg = i_arg + 1
         !!    call get_command_argument(i_arg, value = fscale_)
         !!    read(fscale_,fmt=*,iostat=ierr) fscale
         !!    calculate_selected_eigenvalues = .true.
         else if (trim(tatonam)=='build-KS-orbitals') then
             i_arg = i_arg + 1
             call get_command_argument(i_arg, value = matrix_format)
             i_arg = i_arg + 1
             call get_command_argument(i_arg, value = metadata_file)
             i_arg = i_arg + 1
             call get_command_argument(i_arg, value = coeff_file)
             i_arg = i_arg + 1
             call get_command_argument(i_arg, value = istart_ks_)
             read(istart_ks_,fmt=*,iostat=ierr) istart_ks
             i_arg = i_arg + 1
             call get_command_argument(i_arg, value = iend_ks_)
             read(iend_ks_,fmt=*,iostat=ierr) iend_ks
             build_KS_orbitals = .true.
         else if (trim(tatonam)=='fragment-multipoles') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_format)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fragment_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_matmul_file)
            do l=0,lmax
                do m=-l,l
                    i_arg = i_arg + 1
                    call get_command_argument(i_arg, value = multipoles_files(m,l))
                end do
            end do
            !write(*,'(1x,2a)')&
            !   &   'perform a Loewdin charge analysis'
            calculate_fragment_multipoles = .true.
            exit loop_getargs
         end if
         i_arg = i_arg + 1
      end do loop_getargs
   end if


   if (multipole_analysis) then
       if (bigdft_mpi%iproc==0) then
           call yaml_comment('Multipole analysis',hfill='-')
       end if

       ! Determine the method
       select case(trim(method_name))
       case ('loewdin','LOEWDIN')
           method = CHARGE_ANALYSIS_LOEWDIN
       case ('mulliken','MULLIKEN')
           method = CHARGE_ANALYSIS_MULLIKEN
       !!case ('projector','PROJECTOR')
       !!    method = CHARGE_ANALYSIS_PROJECTOR
       case default
           call f_err_throw('Unknown Method for the multipole analysis',err_name='BIGDFT_INPUT_VARIABLES_ERROR')
       end select

       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       if (bigdft_mpi%iproc==0) then
           call yaml_mapping_open('Atomic System Properties')
           call yaml_map('Types of atoms',smmd%atomnames)
           call yaml_mapping_close()
       end if

       call sparse_matrix_and_matrices_init_from_file_bigdft('serial_text', trim(overlap_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_s, ovrlp_mat, &
            init_matmul=.false.)

       call sparse_matrix_and_matrices_init_from_file_bigdft('serial_text', trim(hamiltonian_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_m, hamiltonian_mat, &
            init_matmul=.false.)

       call sparse_matrix_and_matrices_init_from_file_bigdft('serial_text', trim(kernel_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_l, kernel_mat, &
            init_matmul=.true., filename_mult=trim(kernel_matmul_file))

       if (bigdft_mpi%iproc==0) then
           call yaml_mapping_open('Matrix properties')
           call write_sparsematrix_info(smat_s, 'Overlap matrix')
           call write_sparsematrix_info(smat_m, 'Hamiltonian matrix')
           call write_sparsematrix_info(smat_l, 'Density kernel')
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
               call matrices_init_from_file_bigdft('serial_text', trim(multipoles_files(m,l)), &
                    bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_s, multipoles_matrices(m,l))
           end do
       end do

       call timing(bigdft_mpi%mpi_comm,'INIT','PR')

       select case(method)
       case (CHARGE_ANALYSIS_MULLIKEN)
           methodc='loewdin'
           do_ortho='no'
       case (CHARGE_ANALYSIS_LOEWDIN)
           methodc='loewdin'
           do_ortho='yes'
       !!case (CHARGE_ANALYSIS_PROJECTOR)
       !!    methodc='projector'
       !!    do_ortho='no'
       case default
           call f_err_throw('wrong method',err_name='BIGDFT_RUNTIME_ERROR')
       end select

       ! Determine whether the multipoles matrices to be used are calculated analytically or on the grid.
       select case(trim(matrix_basis))
       case ('wavelet','WAVELET')
           call multipole_analysis_driver_new(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, 0, 11, &
                smmd, smat_s, smat_m, smat_l, &
                ovrlp_mat, hamiltonian_mat, kernel_mat, smmd%rxyz, &
                methodc, do_ortho=trim(do_ortho), projectormode='simple', &
                do_check=.false., calculate_multipole_matrices=.false., centers_auto=.true., &
                write_multipole_matrices_mode=0, &
                multipole_matrix_in=(/(/ovrlp_mat/)/))
       case ('realspace','REALSPACE')
           call multipole_analysis_driver_new(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, ll, 11, &
                smmd, smat_s, smat_m, smat_l, &
                ovrlp_mat, hamiltonian_mat, kernel_mat, smmd%rxyz, &
                methodc, do_ortho=trim(do_ortho), projectormode='simple', &
                do_check=.false., calculate_multipole_matrices=.false., centers_auto=.true., &
                write_multipole_matrices_mode=0, &
                multipole_matrix_in=multipoles_matrices)
       case default
           call f_err_throw('wrong value for matrix_basis',err_name='BIGDFT_RUNTIME_ERROR')
       end select

       call timing(bigdft_mpi%mpi_comm,'CALC','PR')

       call deallocate_sparse_matrix_metadata(smmd)
       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_l)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(kernel_mat)
       call deallocate_sparse_matrix(smat_m)
       call deallocate_matrices(hamiltonian_mat)
       do l=0,ll
           do m=-l,l
               call deallocate_matrices(multipoles_matrices(m,l))
           end do
       end do
       deallocate(multipoles_matrices)

       if (bigdft_mpi%iproc==0) then
           call yaml_comment('done',hfill='-')
       end if

   end if


   if (build_KS_orbitals) then


       if (bigdft_mpi%iproc==0) call yaml_comment('Reading from file '//trim(metadata_file),hfill='~')
       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)

       nspin = 1

       if (nspin/=1) then
           call f_err_throw('KS orbitals nor implemented for nspin/=1')
       end if

       iunit = 99
       if (bigdft_mpi%iproc==0) call yaml_comment('Reading from file '//trim(coeff_file),hfill='~')
       call f_open_file(iunit, file=trim(coeff_file), binary=.false.)
       call read_linear_coefficients(matrix_format, bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
            trim(coeff_file), nspin, ntmb, norbks, coeff_ptr)
       call f_close(iunit)

       in_which_locreg = f_malloc(ntmb,id='in_which_locreg')
       do iorb=1,ntmb
           in_which_locreg(iorb) = iorb
       end do
       !!write(*,*) 'HAC: MANUAL MODIFICATION OF IN_WHICH_LOCREG'
       !!in_which_locreg = (/5,1,2,3,4,6/)


       call orbitals_descriptors(bigdft_mpi%iproc, bigdft_mpi%nproc, ntmb, ntmb, 0, nspin, 1,&
            nkpt, kpt, wkpt, orbs, linear_partition=LINEAR_PARTITION_NONE)
       !!call init_linear_orbs(LINEAR_PARTITION_SIMPLE)


       nspinor = 1

       rxyz = f_malloc_ptr((/3,smmd%nat/),id='rxyz')

       lzd = local_zone_descriptors_null()
       allocate(lzd%llr(orbs%norb),stat=istat)
       do iorb=1,orbs%norb
           call nullify_locreg_descriptors(lzd%llr(iorb))
       end do

       ! Read the global grid
       KSgrid_file = 'KSgrid.dat'
       if (bigdft_mpi%iproc==0) call yaml_comment('Reading from file '//trim(KSgrid_file),hfill='~')
       call f_open_file(iunit, file=trim(KSgrid_file), binary=.false.)
       call io_read_descr_linear(iunit, .true., iiorb, eval_tmp, &
            lzd%glr%d%n1, lzd%glr%d%n2, lzd%glr%d%n3, &
            lzd%glr%ns1, lzd%glr%ns2, lzd%glr%ns3, lzd%hgrids, &
            lstat, error, onwhichatom_tmp, lzd%glr%locrad, lzd%glr%locregCenter, &
            confpotorder, confpotprefac, &
            nvctr_c=lzd%glr%wfd%nvctr_c, &
            nvctr_f=lzd%glr%wfd%nvctr_f, &
            nseg_c=lzd%glr%wfd%nseg_c, &
            nseg_f=lzd%glr%wfd%nseg_f, &
            keygloc=lzd%glr%wfd%keygloc, &
            keyglob=lzd%glr%wfd%keyglob, &
            keyvloc=lzd%glr%wfd%keyvloc, &
            keyvglob=lzd%glr%wfd%keyvglob)
       call f_close(iunit)

       ! THIS MUST BE MADE CLEANEE !!!!!!!!!
       !starting point of the region for interpolating functions grid
       lzd%glr%nsi1= 2 * lzd%glr%ns1 !- (Lnbl1 - Gnbl1)
       lzd%glr%nsi2= 2 * lzd%glr%ns2 !- (Lnbl2 - Gnbl2)
       lzd%glr%nsi3= 2 * lzd%glr%ns3 !- (Lnbl3 - Gnbl3)

       !dimensions of the fine grid inside the localisation region
       lzd%glr%d%nfl1=0
       lzd%glr%d%nfl2=0
       lzd%glr%d%nfl3=0

       !NOTE: This will not work with symmetries (must change it)
       lzd%glr%d%nfu1=lzd%glr%d%n1
       lzd%glr%d%nfu2=lzd%glr%d%n2
       lzd%glr%d%nfu3=lzd%glr%d%n3

       !dimensions of the interpolating scaling functions grid (reduce to +2 for periodic)
       lzd%glr%d%n1i=2*lzd%glr%d%n1+31
       lzd%glr%d%n2i=2*lzd%glr%d%n2+31
       lzd%glr%d%n3i=2*lzd%glr%d%n3+31



       filename='minBasis'
       if (bigdft_mpi%iproc==0) call yaml_comment('Reading from file '//trim(filename)//'*',hfill='~')

       npsidim_orbs = 0
       phi = f_malloc(npsidim_orbs,id='phi')
       !do iat = 1,smmd%nat
          do iorb=1,orbs%norbp
             !if(iat == smmd%on_which_atom(iorb+orbs%isorb)) then
                !shift = 1
                !do jorb = 1, iorb-1
                !   jlr = smmd%on_which_atom(jorb+isorb)
                !   shift = shift + Lzd%Llr(jlr)%wfd%nvctr_c+7*Lzd%Llr(jlr)%wfd%nvctr_f
                !end do
                ilr = in_which_locreg(iorb+orbs%isorb)
                do ispinor=1,nspinor
                   call open_filename_of_iorb(iunit, .false., trim(filename), &
                        orbs, iorb, ispinor, iiorb_out)

                   call io_read_descr_linear(iunit, .true., iiorb, eval_tmp, &
                        lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
                        lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, lzd%hgrids, &
                        lstat, error, onwhichatom_tmp, lzd%llr(ilr)%locrad, lzd%llr(ilr)%locregCenter, &
                        confpotorder, confpotprefac, &
                        nvctr_c=lzd%llr(ilr)%wfd%nvctr_c, &
                        nvctr_f=lzd%llr(ilr)%wfd%nvctr_f, &
                        nseg_c=lzd%llr(ilr)%wfd%nseg_c, &
                        nseg_f=lzd%llr(ilr)%wfd%nseg_f, &
                        keygloc=lzd%llr(ilr)%wfd%keygloc, &
                        keyglob=lzd%llr(ilr)%wfd%keyglob, &
                        keyvloc=lzd%llr(ilr)%wfd%keyvloc, &
                        keyvglob=lzd%llr(ilr)%wfd%keyvglob)

                    ! THIS MUST BE MADE CLEANER !!!!!!!!!
                    !starting point of the region for interpolating functions grid
                    lzd%llr(ilr)%nsi1= 2 * lzd%llr(ilr)%ns1 !- (Lnbl1 - Gnbl1)
                    lzd%llr(ilr)%nsi2= 2 * lzd%llr(ilr)%ns2 !- (Lnbl2 - Gnbl2)
                    lzd%llr(ilr)%nsi3= 2 * lzd%llr(ilr)%ns3 !- (Lnbl3 - Gnbl3)

                    !dimensions of the fine grid inside the localisation region
                    lzd%llr(ilr)%d%nfl1=max(lzd%llr(ilr)%ns1,lzd%glr%d%nfl1)-lzd%llr(ilr)%ns1 ! should we really substract isx (probably because the routines are coded with 0 as origin)?
                    lzd%llr(ilr)%d%nfl2=max(lzd%llr(ilr)%ns2,lzd%glr%d%nfl2)-lzd%llr(ilr)%ns2
                    lzd%llr(ilr)%d%nfl3=max(lzd%llr(ilr)%ns3,lzd%glr%d%nfl3)-lzd%llr(ilr)%ns3

                    !NOTE: This will not work with symmetries (must change it)
                    lzd%llr(ilr)%d%nfu1=min(lzd%llr(ilr)%ns1+lzd%llr(ilr)%d%n1,lzd%glr%d%nfu1)-lzd%llr(ilr)%ns1
                    lzd%llr(ilr)%d%nfu2=min(lzd%llr(ilr)%ns2+lzd%llr(ilr)%d%n2,lzd%glr%d%nfu2)-lzd%llr(ilr)%ns2
                    lzd%llr(ilr)%d%nfu3=min(lzd%llr(ilr)%ns3+lzd%llr(ilr)%d%n3,lzd%glr%d%nfu3)-lzd%llr(ilr)%ns3

                    !dimensions of the interpolating scaling functions grid (reduce to +2 for periodic)
                    lzd%llr(ilr)%d%n1i=2*lzd%llr(ilr)%d%n1+31
                    lzd%llr(ilr)%d%n2i=2*lzd%llr(ilr)%d%n2+31
                    lzd%llr(ilr)%d%n3i=2*lzd%llr(ilr)%d%n3+31


                   phi_tmp = f_malloc(npsidim_orbs,id='phi_tmp')
                   call f_memcpy(src=phi, dest=phi_tmp)
                   call f_free(phi)
                   nsize = lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
                   phi = f_malloc(npsidim_orbs+nsize,id='phi')
                   call f_memcpy(src=phi_tmp, dest=phi)
                   call f_free(phi_tmp)

                   call read_psi_compress(iunit, .true., lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%nvctr_f, &
                        phi(npsidim_orbs+1:npsidim_orbs+nsize), lstat, error)

                   npsidim_orbs = npsidim_orbs + nsize

                   call f_close(iunit)
                end do
             !end if
          enddo
       !end do



       ! have to copy the structures...
       at = atoms_data_null()
       at%astruct%iatype = f_malloc_ptr(smmd%nat,id='at%astruct%iatype')
       at%nzatom = f_malloc_ptr(smmd%ntypes,id='at%astruct%nzatom')
       at%nelpsp = f_malloc_ptr(smmd%ntypes,id='at%astruct%nelpsp')
       call f_memcpy(src=smmd%iatype, dest=at%astruct%iatype)
       call f_memcpy(src=smmd%nzatom, dest=at%nzatom)
       call f_memcpy(src=smmd%nelpsp, dest=at%nelpsp)

       !read(101,*) coeff_ptr
       
       !!rxyz(1:3,1) = (/7.576417, 5.980391, 5.963297/)
       !!rxyz(1:3,2) = (/5.771699, 5.702495, 5.771181/)
       !!rxyz(1:3,3) = (/7.728301, 7.297505, 7.228819/)
       !!write(*,*) 'rxyz',rxyz
       !!write(*,*) 'smmd%rxyz', smmd%rxyz
       do iat=1,smmd%nat
           rxyz(:3,iat) = smmd%rxyz(1:3,iat) + smmd%shift(1:3)
       end do
       if (istart_ks<1) then
           call f_err_throw('istart_ks<1')
       end if
       if (iend_ks>orbs%norb) then
           call f_err_throw('iend_ks>orbs%norb')
       end if
       call build_ks_orbitals_postprocessing(bigdft_mpi%iproc, bigdft_mpi%nproc, &
            orbs%norb, istart_ks, iend_ks, &
            nspin, nspinor, nkpt, kpt, wkpt, in_which_locreg, at, lzd, rxyz, &
            npsidim_orbs, phi, coeff_ptr(:,istart_ks:iend_ks))

       !call deallocate_atoms_data(at)
       call deallocate_local_zone_descriptors(lzd)
       call deallocate_orbitals_data(orbs)
       call f_free(in_which_locreg)
       call f_free_ptr(rxyz)
       call f_free(phi)
       call f_free_ptr(coeff_ptr)
       call deallocate_sparse_matrix_metadata(smmd)
       call f_free_ptr(at%astruct%iatype)
       call f_free_ptr(at%nzatom)
       call f_free_ptr(at%nelpsp)


   end if


   if (calculate_fragment_multipoles) then

       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       if (bigdft_mpi%iproc==0) then
           call yaml_mapping_open('Atomic System Properties')
           call yaml_map('Types of atoms',smmd%atomnames)
           call yaml_mapping_close()
       end if

       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, mpiworld(), smat(1), ovrlp_mat, &
            init_matmul=.false.)

       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(kernel_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, mpiworld(), smat(2), kernel_mat, &
            init_matmul=.true., filename_mult=trim(kernel_matmul_file))

       call init_matrix_taskgroups_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, mpiworld(), .true., 2, smat)

       call resize_matrix_to_taskgroup(smat(1), ovrlp_mat)
       call resize_matrix_to_taskgroup(smat(2), kernel_mat)


       if (bigdft_mpi%iproc==0) then
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
                    bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat(1), multipoles_matrices(m,l))
               call resize_matrix_to_taskgroup(smat(1), multipoles_matrices(m,l))
           end do
       end do

       call timing(bigdft_mpi%mpi_comm,'INIT','PR')

       mpi_env = mpi_environment_null()
       call mpi_environment_set(mpi_env, bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, bigdft_mpi%nproc)
       inquire(file=trim(fragment_file), exist=file_exists)
       if (file_exists) then
           call merge_input_file_to_dict(fragment_dict, fragment_file, mpi_env)
       else
           call f_err_throw('file'//trim(fragment_file)//'not present')
       end if

       supfun_in_fragment = f_malloc(smat(1)%nfvctr,id='supfun_in_fragment')

       ! Calculate the matrix KS
       kernel_eff = matrices_null()
       kernel_eff%matrix_compr = sparsematrix_malloc0_ptr(smat(2),iaction=SPARSE_TASKGROUP,id='kernel_eff%matrix_compr')
       ovrlp_large = matrices_null()
       ovrlp_large%matrix_compr = sparsematrix_malloc_ptr(smat(2), SPARSE_TASKGROUP, id='ovrlp_large%matrix_compr')
       call transform_sparse_matrix(bigdft_mpi%iproc, smat(1), smat(2), SPARSE_TASKGROUP, 'small_to_large', &
            smat_in=ovrlp_mat%matrix_compr, lmat_out=ovrlp_large%matrix_compr)
       ! Should use the highlevel wrapper...
       do ispin=1,smat(2)%nspin
           ist=(ispin-1)*smat(2)%nvctrp_tg+1
           call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smat(2), &
                kernel_mat%matrix_compr(ist:), ovrlp_large%matrix_compr(ist:), kernel_eff%matrix_compr(ist:))
       end do

       ! Calculate the matrix S^-1
       inv_ovrlp(1) = matrices_null()
       inv_ovrlp(1)%matrix_compr = sparsematrix_malloc_ptr(smat(2), SPARSE_TASKGROUP, id='inv_ovrlp%matrix_compr')
       call matrix_chebyshev_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, mpiworld(), 1, [-1.d0], &
            smat(1), smat(2), ovrlp_mat, inv_ovrlp)

       ! Calculate the matrix S^-1*P, where P are the multipole matrices
       do l=0,ll
           do m=-l,l
               call transform_sparse_matrix(bigdft_mpi%iproc, smat(1), smat(2), SPARSE_TASKGROUP, 'small_to_large', &
                    smat_in=multipoles_matrices(m,l)%matrix_compr, lmat_out=ovrlp_large%matrix_compr)
               call f_free_ptr(multipoles_matrices(m,l)%matrix_compr)
               multipoles_matrices(m,l)%matrix_compr = sparsematrix_malloc_ptr(smat(2),iaction=SPARSE_TASKGROUP,&
                                                               id='matrix_compr')
               ! Should use the highlevel wrapper...
               do ispin=1,smat(2)%nspin
                   ist=(ispin-1)*smat(2)%nvctrp_tg+1
                   call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smat(2), &
                        inv_ovrlp(1)%matrix_compr(ist:), ovrlp_large%matrix_compr(ist:), &
                        multipoles_matrices(m,l)%matrix_compr(ist:))
               end do
           end do
       end do

       call deallocate_matrices(ovrlp_large)
       call deallocate_matrices(inv_ovrlp(1))


       nmpmat = 0
       do l=0,ll
           nmpmat = nmpmat + 2*l + 1
       end do
       fragment_multipoles = f_malloc(0.to.nmpmat,id='fragment_multipoles')

       ! Periodicity in the three direction
       perx=(smmd%geocode /= 'F')
       pery=(smmd%geocode == 'P')
       perz=(smmd%geocode /= 'F')


       if (bigdft_mpi%iproc==0) then
           call yaml_sequence_open('Fragment multipoles')
       end if

       fragment_loop: do while (iterating(fragment_item,fragment_dict))

           nat_frag = dict_len(fragment_item)
           fragment_atom_id = f_malloc(nat_frag,id='fragment_atom_id')
           fragment_atom_name = f_malloc_str(len(fragment_atom_name),nat_frag,id='fragment_atom_name')
           fragment_atom => dict_iter(fragment_item)
           iat_frag = 0
           fragment_charge = 0.d0
           fragment_center(1:3) = 0.d0
           do while (associated(fragment_atom))
               iat_frag = iat_frag + 1
               iiat = fragment_atom
               iitype = smmd%iatype(iiat)
               fragment_atom_id(iat_frag) = iiat
               fragment_atom_name(iat_frag) = trim(smmd%atomnames(iitype))
               fragment_charge = fragment_charge + smmd%nelpsp(iitype)
               fragment_center(1:3) = fragment_center(1:3) + smmd%rxyz(1:3,iiat)*smmd%nelpsp(iitype)
               fragment_atom => dict_next(fragment_atom)
           end do
           fragment_center(1:3) = fragment_center(1:3)/fragment_charge

           ! Count how many support functions belong to the fragment
           supfun_in_fragment(:) = .false.
           nsf_frag = 0
           outerloop: do isf=1,smat(1)%nfvctr
               iiat = smmd%on_which_atom(isf)
               do iat=1,nat_frag
                   if (iiat==fragment_atom_id(iat)) then
                       supfun_in_fragment(isf) = .true.
                       nsf_frag = nsf_frag + 1
                       cycle outerloop
                   end if
               end do
           end do outerloop


           lookup_kernel = f_malloc([nsf_frag,nsf_frag],id='lookup_kernel')
           call init_extract_matrix_lookup(smat(2), supfun_in_fragment, nsf_frag, lookup_kernel)
           kqmat = f_malloc([nsf_frag,nsf_frag,2],id='kqmat')
           kmat = f_malloc([nsf_frag,nsf_frag],id='kmat')
           mpmat = f_malloc([nsf_frag,nsf_frag,nmpmat],id='mpmat')

           call f_zero(fragment_multipoles)
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
                       fragment_multipoles(mm) = fragment_multipoles(mm) + tr
                       if (l==0) then
                           if (smat(1)%nspin==1) then
                               call dscal(nsf_frag**2, 0.5d0, kqmat(1,1,1), 1)
                           end if
                           call gemm('n', 'n', nsf_frag, nsf_frag, nsf_frag, 1.d0, kqmat(1,1,1), nsf_frag, &
                                kqmat(1,1,1), nsf_frag, 0.d0, kqmat(1,1,2), nsf_frag)
                           call axpy(nsf_frag**2, -1.d0, kqmat(1,1,1), 1, kqmat(1,1,2), 1)
                           tr = 0.d0
                           do isf=1,nsf_frag
                               tr = tr + kqmat(isf,isf,2)
                           end do
                           fragment_multipoles(0) = fragment_multipoles(0) + tr/fragment_charge
                       end if
                   end do
               end do
            
           end do

           if (bigdft_mpi%iproc==0) then
               call yaml_sequence(advance='no')
               call yaml_map('Atom IDs',fragment_atom_id)
               call yaml_map('Atom names',fragment_atom_name)
               call yaml_map('Neutral fragment charge',fragment_charge,fmt='(1es13.6)')
               call yaml_map('Purity indicator',fragment_multipoles(0),fmt='(1es13.6)')
               mm = 0
               do l=0,ll
                   call yaml_map('q'//adjustl(trim(yaml_toa(l))),fragment_multipoles(mm+1:mm+2*l+1),fmt='(1es13.6)')
                   mm = mm + 2*l+1
               end do
           end if

           call f_free(mpmat)
           call f_free(kmat)
           call f_free(kqmat)
           call f_free(lookup_kernel)
           call f_free(fragment_atom_id)
           call f_free_str(len(fragment_atom_name),fragment_atom_name)

       end do fragment_loop

       if (bigdft_mpi%iproc==0) then
           call yaml_sequence_close()
       end if
       call f_free(supfun_in_fragment)


       call f_free(fragment_multipoles)
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

       call f_timing_checkpoint(ctr_name='LAST',mpi_comm=mpiworld(),nproc=mpisize(),&
            gather_routine=gather_timings)

   end if


   !!!if (solve_eigensystem) then

   !!!    call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
   !!!         bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_s, ovrlp_mat, &
   !!!         init_matmul=.false.)!, nat=nat, rxyz=rxyz, iatype=iatype, ntypes=ntypes, &
   !!!         !nzatom=nzatom, nelpsp=nelpsp, atomnames=atomnames)
   !!!    call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
   !!!    call sparse_matrix_and_matrices_init_from_file_bigdft(trim(hamiltonian_file), &
   !!!         bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_m, hamiltonian_mat, &
   !!!         init_matmul=.false.)

   !!!    ovrlp_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='ovrlp_mat%matrix')
   !!!    call uncompress_matrix(bigdft_mpi%iproc, bigdft_mpi%nproc, &
   !!!         smat_s, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_mat%matrix)
   !!!    hamiltonian_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='hamiltonian_mat%matrix')
   !!!    call uncompress_matrix(bigdft_mpi%iproc, bigdft_mpi%nproc, &
   !!!         smat_m, inmat=hamiltonian_mat%matrix_compr, outmat=hamiltonian_mat%matrix)
   !!!    eval = f_malloc(smat_s%nfvctr,id='eval')

   !!!    if (bigdft_mpi%iproc==0) then
   !!!        call yaml_comment('Diagonalizing the matrix',hfill='~')
   !!!    end if
   !!!    call diagonalizeHamiltonian2(bigdft_mpi%iproc, smat_s%nfvctr, hamiltonian_mat%matrix, ovrlp_mat%matrix, eval)
   !!!    if (bigdft_mpi%iproc==0) then
   !!!        call yaml_comment('Matrix succesfully diagonalized',hfill='~')
   !!!    end if
   !!!    iunit=99
   !!!    call f_open_file(iunit, file=trim(coeff_file), binary=.false.)
   !!!    !call writeLinearCoefficients(iunit, .true., nat, rxyz, smat_s%nfvctr, smat_s%nfvctr, &
   !!!    !     smat_s%nfvctr, hamiltonian_mat%matrix, eval)
   !!!    call write_linear_coefficients(0, trim(coeff_file), smmd%nat, smmd%rxyz, &
   !!!         smmd%iatype, smmd%ntypes, smmd%nzatom, &
   !!!         smmd%nelpsp, smmd%atomnames, smat_s%nfvctr, &
   !!!         smat_s%nfvctr, smat_s%nspin, hamiltonian_mat%matrix, eval)
   !!!    call f_close(iunit)

   !!!    call f_free(eval)
   !!!    call deallocate_matrices(ovrlp_mat)
   !!!    call deallocate_matrices(hamiltonian_mat)
   !!!    call deallocate_sparse_matrix(smat_s)
   !!!    call deallocate_sparse_matrix(smat_m)
   !!!    call deallocate_sparse_matrix_metadata(smmd)
   !!!    !call f_free_ptr(rxyz)
   !!!    !call f_free_ptr(iatype)
   !!!    !call f_free_ptr(nzatom)
   !!!    !call f_free_ptr(nelpsp)
   !!!    !call f_free_str_ptr(len(atomnames),atomnames)
   !!!end if


   !!!if (calculate_pdos) then
   !!!    iunit = 99
   !!!    if (bigdft_mpi%iproc==0) call yaml_comment('Reading from file '//trim(coeff_file),hfill='~')
   !!!    call f_open_file(iunit, file=trim(coeff_file), binary=.false.)
   !!!    call read_linear_coefficients(trim(coeff_file), nspin, ntmb, norbks, coeff_ptr, &
   !!!         eval=eval_ptr)
   !!!    call f_close(iunit)

   !!!    if (bigdft_mpi%iproc==0) call yaml_comment('Reading from file '//trim(overlap_file),hfill='~')
   !!!    call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
   !!!         bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_s, ovrlp_mat, &
   !!!         init_matmul=.false.)!, iatype=iatype, ntypes=ntypes, atomnames=atomnames, &
   !!!         !on_which_atom=on_which_atom)
   !!!    call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
   !!!    if (ntmb/=smat_s%nfvctr) call f_err_throw('ntmb/=smat_s%nfvctr')
   !!!    ovrlp_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_PARALLEL, id='ovrlp_mat%matrix')
   !!!    !call uncompress_matrix(bigdft_mpi%iproc, smat_s, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_mat%matrix)
   !!!    call uncompress_matrix_distributed2(bigdft_mpi%iproc, smat_s, DENSE_PARALLEL, &
   !!!         ovrlp_mat%matrix_compr, ovrlp_mat%matrix(1:,1:,1))

   !!!    if (bigdft_mpi%iproc==0) call yaml_comment('Reading from file '//trim(hamiltonian_file),hfill='~')
   !!!    call sparse_matrix_and_matrices_init_from_file_bigdft(trim(hamiltonian_file), &
   !!!         bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_m, hamiltonian_mat, &
   !!!         init_matmul=.false.)
   !!!    hamiltonian_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_PARALLEL, id='hamiltonian_mat%matrix')
   !!!    !call uncompress_matrix(bigdft_mpi%iproc, smat_m, inmat=hamiltonian_mat%matrix_compr, outmat=hamiltonian_mat%matrix)
   !!!    call uncompress_matrix_distributed2(bigdft_mpi%iproc, smat_m, DENSE_PARALLEL, &
   !!!         hamiltonian_mat%matrix_compr, hamiltonian_mat%matrix(1:,1:,1))

   !!!    iunit=99
   !!!    call f_open_file(iunit, file=pdos_file, binary=.false.)
   !!!    read(iunit,*) npdos

   !!!    calc_array = f_malloc((/ntmb,npdos/),id='calc_array')
   !!!    pdos_name = f_malloc_str(len(pdos_name),npdos,id='pdos_name')


   !!!    do ipdos=1,npdos
   !!!        do itmb=1,ntmb
   !!!            calc_array(itmb,ipdos) = .false.
   !!!        end do
   !!!    end do
   !!!    ipdos = 0
   !!!    !npdos_loop: do !ipdos=1,npdos
   !!!    if (bigdft_mpi%iproc==0) then
   !!!        call yaml_sequence_open('Atoms and support functions to be taken into account for each partial density of states')
   !!!    end if
   !!!        do 
   !!!            !read(iunit01,*,iostat=ios) cc, ival
   !!!            read(iunit,'(a128)',iostat=ios) line
   !!!            if (ios/=0) exit
   !!!            !write(*,*) 'line',line
   !!!            read(line,*,iostat=ios) cc, cname
   !!!            if (cc=='#') then
   !!!                ipdos = ipdos + 1
   !!!                pdos_name(ipdos) = trim(cname)
   !!!                if (bigdft_mpi%iproc==0) then
   !!!                    if (ipdos>1) then
   !!!                        call yaml_mapping_close()
   !!!                    end if
   !!!                    call yaml_sequence(advance='no')
   !!!                    call yaml_mapping_open(trim(pdos_name(ipdos)))
   !!!                end if
   !!!                cycle 
   !!!            end if
   !!!            read(line,*,iostat=ios) cc, ival
   !!!            if (bigdft_mpi%iproc==0) then
   !!!                call yaml_map(trim(cc),ival)
   !!!            end if
   !!!            found = .false.
   !!!            do itype=1,smmd%ntypes
   !!!                if (trim(smmd%atomnames(itype))==trim(cc)) then
   !!!                    iitype = itype
   !!!                    found = .true.
   !!!                    exit
   !!!                end if
   !!!            end do
   !!!            if (.not.found) then
   !!!                call f_err_throw('Atom type'//trim(cc)//' is unknown', &
   !!!                     err_name='BIGDFT_RUNTIME_ERROR')
   !!!            end if
   !!!            iat_prev = -1
   !!!            do itmb=1,ntmb
   !!!                iat = smmd%on_which_atom(itmb)
   !!!                if (iat/=iat_prev) then
   !!!                    ii = 0
   !!!                end if
   !!!                iat_prev = iat
   !!!                itype = smmd%iatype(iat)
   !!!                ii = ii + 1
   !!!                if (itype==iitype .and. ii==ival) then
   !!!                    if (calc_array(itmb,ipdos)) then
   !!!                        call f_err_throw('calc_array(itmb) must not be true here', &
   !!!                             err_name='BIGDFT_RUNTIME_ERROR')
   !!!                    end if
   !!!                    calc_array(itmb,ipdos) = .true.
   !!!                end if
   !!!            end do
   !!!        end do
   !!!    if (bigdft_mpi%iproc==0) then
   !!!        call yaml_mapping_close()
   !!!        call yaml_sequence_close()
   !!!    end if

   !!!    !end do npdos_loop
   !!!    call f_close(iunit)


   !!!    energy_arr = f_malloc0(norbks,id='energy_arr')
   !!!    occup_arr = f_malloc0((/npdos,norbks/),id='occup_arr')
   !!!    occups = f_malloc0(npdos,id='occups')

   !!!    denskernel = f_malloc((/ntmb,smat_s%nfvctrp/),id='denskernel')
   !!!    pdos = f_malloc0((/npt,npdos/),id='pdos')
   !!!    if (bigdft_mpi%iproc==0) then
   !!!        call yaml_comment('PDoS calculation',hfill='~')
   !!!        call yaml_mapping_open('Calculating PDoS')
   !!!    end if
   !!!    ! Calculate a partial kernel for each KS orbital
   !!!    !do ipdos=1,npdos
   !!!        !if (bigdft_mpi%iproc==0) call yaml_map('PDoS number',ipdos)
   !!!        !if (bigdft_mpi%iproc==0) call yaml_map('start, increment',(/ipdos,npdos/))
   !!!        !write(num,fmt='(i2.2)') ipdos
   !!!        if (bigdft_mpi%iproc==0) bar=f_progress_bar_new(nstep=norbks)
   !!!        do iorb=1,norbks
   !!!            if (bigdft_mpi%iproc==0) call yaml_map('orbital being processed',iorb)
   !!!            call gemm('n', 't', ntmb, smat_s%nfvctrp, 1, 1.d0, coeff_ptr(1,iorb), ntmb, &
   !!!                 coeff_ptr(smat_s%isfvctr+1,iorb), ntmb, 0.d0, denskernel(1,1), ntmb)
   !!!            energy = 0.d0
   !!!            call f_zero(occups)
   !!!            do ispin=1,nspin
   !!!                !$omp parallel default(none) &
   !!!                !$omp shared(ispin,smat_s,ntmb,denskernel,hamiltonian_mat,energy) &
   !!!                !$omp private(itmb,jtmb,jjtmb)
   !!!                !$omp do reduction(+:energy)
   !!!                do jtmb=1,smat_s%nfvctrp
   !!!                    jjtmb = smat_s%isfvctr + jtmb
   !!!                    do itmb=1,ntmb
   !!!                        energy = energy + denskernel(itmb,jtmb)*hamiltonian_mat%matrix(itmb,jtmb,ispin)
   !!!                    end do
   !!!                end do
   !!!                !$omp end do
   !!!                !$omp end parallel
   !!!                energy_arr(iorb) = energy
   !!!            end do
   !!!            !call mpiallred(energy, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
   !!!            do ispin=1,nspin
   !!!                !$omp parallel default(none) &
   !!!                !$omp shared(ispin,smat_s,ntmb,denskernel,ovrlp_mat,calc_array,npdos,occups) &
   !!!                !$omp private(itmb,jtmb,jjtmb,occup,ipdos)
   !!!                !$omp do reduction(+:occups)
   !!!                do jtmb=1,smat_s%nfvctrp
   !!!                    jjtmb = smat_s%isfvctr + jtmb
   !!!                    !if (calc_array(jjtmb,ipdos)) then
   !!!                        do itmb=1,ntmb!ipdos,ntmb,npdos
   !!!                            !if (calc_array(itmb,ipdos)) then
   !!!                                occup = denskernel(itmb,jtmb)*ovrlp_mat%matrix(itmb,jtmb,ispin)
   !!!                                do ipdos=1,npdos
   !!!                                    if (calc_array(itmb,ipdos) .and. calc_array(jjtmb,ipdos)) then
   !!!                                        occups(ipdos) = occups(ipdos) + occup
   !!!                                    end if
   !!!                                end do
   !!!                                !occup = occup + denskernel(itmb,jtmb)*ovrlp_mat%matrix(itmb,jtmb,ispin)
   !!!                            !end if
   !!!                        end do
   !!!                    !end if
   !!!                end do
   !!!                !$omp end do
   !!!                !$omp end parallel
   !!!                do ipdos=1,npdos
   !!!                    occup_arr(ipdos,iorb) = occups(ipdos)
   !!!                end do
   !!!            end do
   !!!            !occup_pdos = occup_pdos + occup
   !!!            !!if (iorb<norbks) then
   !!!            !!    write(iunit,'(2(a,es16.9),a)') '  ',occup,'*exp(-(x-',energy,')**2/(2*sigma**2)) + '//trim(backslash)
   !!!            !!else
   !!!            !!    write(iunit,'(2(a,es16.9),a)') '  ',occup,'*exp(-(x-',energy,')**2/(2*sigma**2))'
   !!!            !!end if
   !!!            !!!write(*,'(a,i6,3es16.8)')'iorb, eval(iorb), energy, occup', iorb, eval(iorb), energy, occup
   !!!            if (bigdft_mpi%iproc==0) call dump_progress_bar(bar,step=iorb)
   !!!        end do
   !!!        !total_occup = total_occup + occup_pdos
   !!!        !!if (ipdos==1) then
   !!!        !!    write(iunit,'(a,i0,a)') "plot f",ipdos,"(x) lc rgb 'color' lt 1 lw 2 w l title 'name'"
   !!!        !!else
   !!!        !!    write(iunit,'(a,i0,a)') "replot f",ipdos,"(x) lc rgb 'color' lt 1 lw 2 w l title 'name'"
   !!!        !!end if
   !!!        !!if (bigdft_mpi%iproc==0) call yaml_map('sum of PDoS',occup_pdos)
   !!!        !!output_pdos='PDoS_'//num//'.dat'
   !!!        !!call yaml_map('output file',trim(output_pdos))
   !!!        !!call f_open_file(iunit01, file=trim(output_pdos), binary=.false.)
   !!!        !!write(iunit01,'(a)') '#             energy                pdos'
   !!!        !!do ipt=1,npt
   !!!        !!    write(iunit01,'(2es20.12)') energy_bins(1,ipt), pdos(ipt,ipdos)
   !!!        !!end do
   !!!        !!call f_close(iunit01)
   !!!    !end do
   !!!    call mpiallred(occup_arr, mpi_sum, comm=bigdft_mpi%mpi_comm)
   !!!    call mpiallred(energy_arr, mpi_sum, comm=bigdft_mpi%mpi_comm)
   !!!    if (bigdft_mpi%iproc==0) call yaml_mapping_close()

   !!!    if (bigdft_mpi%iproc==0) then
   !!!        call yaml_comment('Calculation complete',hfill='=')
   !!!        output_pdos='PDoS.gp'
   !!!        call yaml_map('output file',trim(output_pdos))
   !!!        iunit = 99
   !!!        call f_open_file(iunit, file=trim(output_pdos), binary=.false.)
   !!!        write(iunit,'(a)') '# plot the DOS as a sum of Gaussians'
   !!!        write(iunit,'(a)') 'set samples 1000'
   !!!        write(iunit,'(a,2(es12.5,a))') 'set xrange[',eval_ptr(1),':',eval_ptr(ntmb),']'
   !!!        write(iunit,'(a)') 'sigma=0.01'
   !!!        write(backslash,'(a)') '\ '
   !!!        total_occup = 0.d0
   !!!        do ipdos=1,npdos
   !!!            occup_pdos = 0.d0
   !!!            call yaml_map('PDoS number',ipdos)
   !!!            call yaml_map('start, increment',(/ipdos,npdos/))
   !!!            write(num,fmt='(i2.2)') ipdos
   !!!            write(iunit,'(a,i0,a)') 'f',ipdos,'(x) = '//trim(backslash)
   !!!            do iorb=1,norbks
   !!!                if (iorb<norbks) then
   !!!                    write(iunit,'(2(a,es16.9),a)') '  ',occup_arr(ipdos,iorb),&
   !!!                        &'*exp(-(x-',energy_arr(iorb),')**2/(2*sigma**2)) + '//trim(backslash)
   !!!                else
   !!!                    write(iunit,'(2(a,es16.9),a)') '  ',occup_arr(ipdos,iorb),&
   !!!                        &'*exp(-(x-',energy_arr(iorb),')**2/(2*sigma**2))'
   !!!                end if
   !!!                occup_pdos = occup_pdos + occup_arr(ipdos,iorb)
   !!!            end do
   !!!            total_occup = total_occup + occup_pdos
   !!!            call yaml_map('sum of PDoS',occup_pdos)
   !!!        end do
   !!!        do ipdos=1,npdos
   !!!            if (ipdos<ncolors) then
   !!!                colorname = colors(ipdos)
   !!!            else
   !!!                colorname = 'color'
   !!!            end if
   !!!            if (ipdos<npdos) then
   !!!                lineend = ' ,'//trim(backslash)
   !!!            else
   !!!                lineend = ''
   !!!            end if
   !!!            if (ipdos==1) then
   !!!                linestart = 'plot'
   !!!            else
   !!!                linestart = '    '
   !!!            end if
   !!!            write(iunit,'(a,i0,a)') trim(linestart)//" f",ipdos,"(x) lc rgb '"//trim(colorname)//&
   !!!                &"' lt 1 lw 2 w l title '"//trim(pdos_name(ipdos))//"'"//trim(lineend)
   !!!        end do
   !!!        call f_close(iunit)
   !!!    end if
   !!!    if (bigdft_mpi%iproc==0) call yaml_map('sum of total DoS',total_occup)
   !!!    call mpibarrier(bigdft_mpi%mpi_comm)

   !!!    call f_free(pdos)
   !!!    call f_free(denskernel)
   !!!    call f_free_ptr(eval_ptr)
   !!!    call deallocate_matrices(ovrlp_mat)
   !!!    call deallocate_matrices(hamiltonian_mat)
   !!!    call deallocate_sparse_matrix(smat_s)
   !!!    call deallocate_sparse_matrix(smat_m)
   !!!    call deallocate_sparse_matrix_metadata(smmd)
   !!!    !call f_free_ptr(iatype)
   !!!    !call f_free_str_ptr(len(atomnames),atomnames)
   !!!    call f_free_str(len(pdos_name),pdos_name)
   !!!    call f_free(calc_array)
   !!!    !call f_free_ptr(on_which_atom)
   !!!    call f_free_ptr(coeff_ptr)
   !!!    call f_free(energy_arr)
   !!!    call f_free(occup_arr)
   !!!    call f_free(occups)
   !!!end if

   !!!if (convert_matrix_format) then
   !!!    select case (trim(conversion))
   !!!    case ('bigdft_to_ccs')
   !!!        iconv = 1
   !!!    case ('ccs_to_bigdft')
   !!!        iconv = 2
   !!!    case default
   !!!        call f_err_throw("wrong value for conversion; possible are 'bigdft_to_ccs' and 'ccs_to_bigdft'")
   !!!    end select

   !!!    if (iconv==1) then
   !!!        call sparse_matrix_and_matrices_init_from_file_bigdft(trim(infile), &
   !!!             bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
   !!!             smat, mat, init_matmul=.false.)
   !!!        row_ind = f_malloc_ptr(smat%nvctr,id='row_ind')
   !!!        col_ptr = f_malloc_ptr(smat%nfvctr,id='col_ptr')
   !!!        call ccs_data_from_sparse_matrix(smat, row_ind, col_ptr)
   !!!        if (bigdft_mpi%iproc==0) call ccs_matrix_write(trim(outfile), smat, row_ind, col_ptr, mat)
   !!!        call f_free_ptr(row_ind)
   !!!        call f_free_ptr(col_ptr)
   !!!    else if (iconv==2) then
   !!!        call sparse_matrix_and_matrices_init_from_file_ccs(trim(infile), bigdft_mpi%iproc, bigdft_mpi%nproc, &
   !!!             bigdft_mpi%mpi_comm, smat, mat, init_matmul=.false.)
   !!!        call write_sparse_matrix(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
   !!!             smat, mat, trim(outfile))
   !!!    end if

   !!!    call deallocate_sparse_matrix(smat)
   !!!    call deallocate_matrices(mat)
   !!!end if


   !!!if (calculate_selected_eigenvalues) then
   !!!    call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
   !!!    call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
   !!!         bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_s, ovrlp_mat, &
   !!!         init_matmul=.false.)
   !!!    call sparse_matrix_and_matrices_init_from_file_bigdft(trim(hamiltonian_file), &
   !!!         bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_m, hamiltonian_mat, &
   !!!         init_matmul=.false.)
   !!!    call sparse_matrix_and_matrices_init_from_file_bigdft(trim(kernel_file), &
   !!!         bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smat_l, kernel_mat, &
   !!!         init_matmul=.true.)
   !!!    call matrices_init(smat_l, ovrlp_minus_one_half(1))

   !!!    call timing(mpi_comm_world,'INIT','PR')
   !!!    if (iev_min<1 .or. iev_min>smat_s%nfvctr .or. iev_max>smat_s%nfvctr .or. iev_max<1) then
   !!!        if (bigdft_mpi%iproc==0) then
   !!!            call yaml_warning('The required eigenvalues are outside of the possible range, automatic ajustment')
   !!!        end if
   !!!    end if
   !!!    iev_min = max(iev_min,1)
   !!!    iev_min = min(iev_min,smat_s%nfvctr)
   !!!    iev_max = min(iev_max,smat_s%nfvctr)
   !!!    iev_max = max(iev_max,1)
   !!!    eval = f_malloc(iev_min.to.iev_max,id='eval')
   !!!    if (bigdft_mpi%iproc==0) then
   !!!        call yaml_mapping_open('Calculating eigenvalues using FOE')
   !!!    end if
   !!!    call get_selected_eigenvalues_from_FOE(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
   !!!         iev_min, iev_max, smat_s, smat_m, smat_l, ovrlp_mat, hamiltonian_mat, &
   !!!         ovrlp_minus_one_half, eval, fscale)

   !!!    call timing(mpi_comm_world,'CALC','PR')

   !!!    if (bigdft_mpi%iproc==0) then
   !!!        call yaml_sequence_open('values')
   !!!        do iev=iev_min,iev_max
   !!!            call yaml_sequence(advance='no')
   !!!            call yaml_mapping_open(flow=.true.)
   !!!            call yaml_map('ID',iev,fmt='(i6.6)')
   !!!            call yaml_map('eval',eval(iev),fmt='(es12.5)')
   !!!            call yaml_mapping_close()
   !!!        end do
   !!!        call yaml_sequence_close()
   !!!        call yaml_mapping_close()
   !!!    end if


   !!!    call deallocate_sparse_matrix(smat_s)
   !!!    call deallocate_sparse_matrix(smat_m)
   !!!    call deallocate_sparse_matrix(smat_l)
   !!!    call deallocate_matrices(ovrlp_mat)
   !!!    call deallocate_matrices(hamiltonian_mat)
   !!!    call deallocate_matrices(kernel_mat)
   !!!    call deallocate_matrices(ovrlp_minus_one_half(1))
   !!!    call deallocate_sparse_matrix_metadata(smmd)
   !!!    call f_free(eval)

   !!!    call timing(mpi_comm_world,'LAST','PR')

   !!!end if

   call build_dict_info(dict_timing_info)
   call f_timing_stop(mpi_comm=bigdft_mpi%mpi_comm,nproc=bigdft_mpi%nproc,&
        gather_routine=gather_timings,dict_info=dict_timing_info)
   call dict_free(dict_timing_info)

   if (bigdft_mpi%iproc==0) then
       call yaml_release_document()
   end if


   call bigdft_finalize(ierr)
   call f_lib_finalize()


  contains


    !> construct the dictionary needed for the timing information
    !! SM: This routine should go to a module
    subroutine build_dict_info(dict_info)
      use wrapper_MPI
      use dynamic_memory
      use dictionaries
      implicit none

      type(dictionary), pointer :: dict_info
      !local variables
      integer :: ierr,namelen,nthreads
      character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
      character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename
      type(dictionary), pointer :: dict_tmp
      !$ integer :: omp_get_max_threads

      call dict_init(dict_info)
!  bastian: comment out 4 followinf lines for debug purposes (7.12.2014)
      !if (DoLastRunThings) then
         call f_malloc_dump_status(dict_summary=dict_tmp)
         call set(dict_info//'Routines timing and number of calls',dict_tmp)
      !end if
      nthreads = 0
      !$  nthreads=omp_get_max_threads()
      call set(dict_info//'CPU parallelism'//'MPI tasks',bigdft_mpi%nproc)
      if (nthreads /= 0) call set(dict_info//'CPU parallelism'//'OMP threads',&
           nthreads)

      nodename=f_malloc0_str(MPI_MAX_PROCESSOR_NAME,0.to.bigdft_mpi%nproc-1,id='nodename')
      if (bigdft_mpi%nproc>1) then
         call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
         !gather the result between all the process
         call MPI_GATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
              nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,&
              bigdft_mpi%mpi_comm,ierr)
         if (bigdft_mpi%iproc==0) call set(dict_info//'Hostnames',&
                 list_new(.item. nodename))
      end if
      call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)

    end subroutine build_dict_info


end program utilities

!!$!> extract the different wavefunctions to verify if the completeness relation is satisfied
!!$subroutine completeness_relation
!!$
!!$  call wfn_filename(filename_out,radical,binary,ikpt,nspinor,nspin,ispinor,spin,iorb)
!!$
!!$  !loop that has to be done for each of the wavefunctions
!!$  unitwf=99
!!$  call f_open_file(unit=unitwf,file=filename_out)
!!$  call readonewave(unitwf,.not. binary,iorb,bigdft_mpi%iproc,&
!!$       it%lr%d%n1,it%lr%d%n2,it%lr%d%n3, &
!!$       Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
!!$       at,it%lr%wfd,rxyz_old,rxyz,&
!!$       psi_ptr,eval,psifscf)
!!$  call f_close(unitwf)
!!$
!!$end subroutine completeness_relation
