module fragment_functions
  use dynamic_memory
  use f_utils
  use dictionaries

  ! Data
  private

  integer :: N, nel, maxfrag, ifound
  real(kind=8),dimension(:,:),pointer :: K, xold

  public :: set_N
  public :: set_nel
  public :: set_maxfrag
  public :: allocate_K
  public :: allocate_xold
  public :: set_K
  public :: fcn_fragments
  public :: get_subkernel_purity
  public :: deallocate_K
  public :: deallocate_xold
  public :: orthogonalize_inputguess
  public :: set_ifound
  public :: set_xold


  contains

  subroutine set_N(N_in)
    implicit none
    integer,intent(in) :: N_in
    N = N_in
  end subroutine set_N

  subroutine set_nel(nel_in)
    implicit none
    integer,intent(in) :: nel_in
    nel = nel_in
  end subroutine set_nel

  subroutine set_maxfrag(maxfrag_in)
    implicit none
    integer,intent(in) :: maxfrag_in
    maxfrag = maxfrag_in
  end subroutine set_maxfrag

  subroutine set_ifound(ifound_in)
    implicit none
    integer,intent(in) :: ifound_in
    ifound = ifound_in
  end subroutine set_ifound

  subroutine allocate_K()
    implicit none
    K = f_malloc_ptr([N,N],id='K')
  end subroutine allocate_K

  subroutine allocate_xold()
    implicit none
    xold = f_malloc_ptr([N,maxfrag],id='xold')
  end subroutine allocate_xold

  subroutine deallocate_K()
    implicit none
    call f_free_ptr(K)
  end subroutine deallocate_K

  subroutine deallocate_xold()
    implicit none
    call f_free_ptr(xold)
  end subroutine deallocate_xold

  subroutine set_K(K_in)
    implicit none
    real(kind=8),dimension(:,:),intent(in) :: K_in
    if (size(K_in,1)/=N) call f_err_throw('wrong 1st dimension of K_in')
    if (size(K_in,2)/=N) call f_err_throw('wrong 2nd dimension of K_in')
    if (.not.associated(K)) then
        call f_err_throw('K is not associated')
    end if
    if (size(K,1)/=N) call f_err_throw('wrong 1st dimension of K')
    if (size(K,2)/=N) call f_err_throw('wrong 2nd dimension of K')
    call f_memcpy(src=K_in, dest=K)
  end subroutine set_K


  subroutine set_xold(it, weights)
    implicit none
    integer,intent(in) :: it
    real(kind=8),dimension(:),intent(in) :: weights
    if (size(weights)/=N) then
        call f_err_throw('wrong 1st dimension of weights')
    end if
    call f_memcpy(src=weights, dest=xold(1:N,it))
  end subroutine set_xold




subroutine fcn_fragments(m,n,x,fvec,iflag)
integer m,n,iflag
double precision x(n),fvec(m)
!real(kind=8),dimension(:,:),allocatable :: tmpmat
!real(kind=8),dimension(12) :: xold=[1,1,1,1,1,1,0,0,0,0,0,0]

  if (m/=n+2+ifound) then
      write(*,*) 'm, n, ifound', m, n, ifound
      call f_err_throw('m/=n+2+ifound')
  end if

  do i=1,n
      fvec(i) = 0.d0
      do j=1,n
          fvec(i) = fvec(i) + K(j,i)**2*x(j)
      end do
      fvec(i) = fvec(i) - x(i)*K(i,i)
  end do

  fvec(n+1) = 0.d0
  do j=1,n
      !fvec(n+1) = fvec(n+1) + x(j)*(x(j)-1.d0)
      fvec(n+1) = fvec(n+1) + (x(j)*(x(j)-1.d0))**2
      !fvec(n+1) = fvec(n+1) + x(j)**4 -2.d0*x(j)**3 + x(j)**2
      !write(*,*) 'new, total', x(j)**4 -2.d0*x(j)**3 + x(j)**2, fvec(n+1)
  end do

  fvec(n+2) = 0.d0
  do j=1,n
      !fvec(n+2) = fvec(n+2) + x(j)**2
      fvec(n+2) = fvec(n+2) + x(j)**1
  end do
  fvec(n+2) = fvec(n+2) - real(nel,kind=8)

  call f_memcpy(src=x(1:n), dest=xold(1:n,m-(n+1)))

  ! This will only be done if m>=n+3, meaning that we have previous solutions
  do i=n+3,m
      !write(*,*) 'additional constraint: i', i
      fvec(i) = 0.d0
      do j=1,n
          fvec(i) = fvec(i) + x(j)*xold(j,i-(n+2))
      end do
  end do
  
  !write(*,'(a,12(f8.3),es12.4)') 'fvec**2',fvec(1:12)**2, sum(fvec(1:12)**2)
  !write(*,'(a,es12.4)') 'sum(fvec**2)', sum(fvec**2)


!----------
!calculate the functions at x and
!return this vector in fvec.
!----------
return
end subroutine fcn_fragments


function get_subkernel_purity(n, K, in_frag, charge) result(trace)
implicit none

! Calling arguments
integer,intent(in) :: n
real(kind=8),dimension(n,n),intent(in) :: K
logical,dimension(n),intent(in) :: in_frag
real(kind=8),intent(in) :: charge
real(kind=8) :: trace

! Local variables
integer :: i, j
real(kind=8) :: tt

  trace = 0.d0
  do i=1,n
      if (in_frag(i)) then
          tt = 0.d0
          do j=1,n
              if (in_frag(j)) then
                  tt = tt + K(j,i)**2
              end if
          end do
          !write(*,*) 'tt, K(i,i)', tt, K(i,i)
          tt = tt - K(i,i)
          trace = trace + tt
      end if
  end do
  trace = trace/charge

end function get_subkernel_purity




subroutine orthogonalize_inputguess(m, n, x)
use wrapper_linalg
implicit none
integer,intent(in) :: n, m
real(kind=8),dimension(n),intent(inout) :: x
real(kind=8) :: tt, ttmax
integer :: i, iortho
  ! This will only be done if m>=n+3, meaning that we have previous solutions
  do iortho=1,10
      do i=n+3,m
          tt = dot(n, x(1), 1, xold(1,i-(n+2)), 1)
          !tt = tt/sqrt(dot(n, x(1), 1, x(1), 1))
          tt = tt/dot(n, xold(1,i-(n+2)), 1, xold(1,i-(n+2)), 1)
          call axpy(n, -tt, xold(1,i-(n+2)), 1, x(1), 1)
      end do
      ttmax = 0.d0
      do i=n+3,m
          tt = dot(n, x(1), 1, xold(1,i-(n+2)), 1)
          if (abs(tt)>ttmax) ttmax = abs(tt)
          !write(*,*) 'i, tt', i, tt
      end do
      if (ttmax<1.d-2) then
          !write(*,*) 'exit after iortho, ttmax',iortho, ttmax
          exit
      end if
  end do


end subroutine orthogonalize_inputguess

end module fragment_functions







!> @file
!!   Auxiliary program to perform pre-/post-processing
!! @author
!!   Copyright (C) 2016 CheSS developers
!!
!!   This file is part of CheSS.
!!   
!!   CheSS is free software: you can redistribute it and/or modify
!!   it under the terms of the GNU Lesser General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.
!!   
!!   CheSS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU Lesser General Public License for more details.
!!   
!!   You should have received a copy of the GNU Lesser General Public License
!!   along with CheSS.  If not, see <http://www.gnu.org/licenses/>.


program chess_toolbox

   !use module_base
  !use yaml_output
  use yaml_output
  use dictionaries
  use time_profiling
  use f_utils
  use wrapper_mpi
  use dynamic_memory
  use wrapper_linalg
   !!use module_types, only: bigdft_init_errors, bigdft_init_timing_categories
   !!use module_atoms, only: atoms_data, atoms_data_null, deallocate_atoms_data
   use sparsematrix_base
   use sparsematrix_init, only: bigdft_to_sparsebigdft, distribute_columns_on_processes_simple, &
                                write_sparsematrix_info, init_matrix_taskgroups_wrapper
   use sparsematrix_io, only: read_sparse_matrix, write_sparse_matrix, write_dense_matrix, read_dense_matrix
   use sparsematrix, only: uncompress_matrix, uncompress_matrix_distributed2, diagonalizeHamiltonian2, &
                           transform_sparse_matrix, compress_matrix, get_minmax_eigenvalues, &
                           resize_matrix_to_taskgroup, uncompress_matrix2
   use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_bigdft, &
                                     sparse_matrix_and_matrices_init_from_file_ccs, &
                                     sparse_matrix_metadata_init_from_file, &
                                     matrices_init_from_file_bigdft, &
                                     ccs_data_from_sparse_matrix, &
                                     ccs_matrix_write, &
                                     matrices_init, &
                                     matrix_matrix_multiplication, &
                                     matrix_chebyshev_expansion
   !!use postprocessing_linear, only: CHARGE_ANALYSIS_LOEWDIN, CHARGE_ANALYSIS_MULLIKEN
   !!use multipole, only: multipole_analysis_driver_new
   !!use multipole_base, only: lmax
   use sparsematrix_io, only: write_linear_coefficients, read_linear_coefficients
   !!use bigdft_run, only: bigdft_init
   use matrix_operations, only: matrix_for_orthonormal_basis
   use parallel_linalg, only: dgemm_parallel
   use f_random, only: f_random_number
   use highlevel_wrappers, only: calculate_eigenvalues, solve_eigensystem_lapack
   use foe, only: overlap_plusminus_onehalf, calculate_entropy_term
   use fragment_functions
   use minpack
   implicit none
   external :: gather_timings
   character(len=*), parameter :: subname='utilities'
   !character(len=1) :: geocode
   !character(len=3) :: do_ortho
   character(len=30) :: tatonam, radical, colorname, linestart, lineend, cname
   character(len=128) :: overlap_file, hamiltonian_file, hamiltonian_manipulated_file
   character(len=128) :: kernel_file, coeff_file, eval_file, pdos_file, metadata_file, output_file, output_bins_file
   character(len=128) :: line, cc, output_pdos, conversion, infile, outfile, iev_min_, iev_max_, fscale_
   character(len=128) :: ihomo_state_, homo_value_, lumo_value_, smallest_value_, largest_value_, scalapack_blocksize_, kT_
   character(len=128) :: accuracy_entropy_, nbin_, itype_, only_evals_, only_binned_values_
   character(len=128) :: frag_elements_start_, frag_elements_end_
   !!character(len=128),dimension(-lmax:lmax,0:lmax) :: multipoles_files
   character(len=128) :: kernel_matmul_file, fragment_file, manipulation_mode, diag_algorithm
   !logical :: multipole_analysis = .false.
   logical :: solve_eigensystem = .false.
   logical :: calculate_pdos = .false.
   logical :: convert_matrix_format = .false.
   logical :: calculate_selected_eigenvalues = .false.
   logical :: kernel_purity = .false.
   logical :: manipulate_eigenvalue_spectrum = .false.
   logical :: calculate_entropy = .false.
   logical :: analyze_density_matrix_dense = .false.
   logical :: find_fragments = .false.
   !!type(atoms_data) :: at
   type(sparse_matrix_metadata) :: smmd
   integer :: iproc, nproc
   integer :: istat, i_arg, ierr, nspin, nthread
   integer :: nfvctr_m, nseg_m, nvctr_m
   integer :: nfvctr_l, nseg_l, nvctr_l
   !integer :: nfvctrp_l, isfvctr_l, nfvctrp_m, isfvctr_m, nfvctrp_s, isfvctr_s
   integer :: iconv, iev_min, iev_max, jat, jat_start, jtype
   !integer,dimension(:),pointer :: on_which_atom
   !integer,dimension(:),pointer :: keyv_s, keyv_m, keyv_l
   !integer,dimension(:),pointer :: on_which_atom_s, on_which_atom_m, on_which_atom_l
   !integer,dimension(:),pointer :: nzatom
   integer,dimension(:),pointer :: col_ptr, row_ind
   !integer,dimension(:,:,:),pointer :: keyg_s, keyg_m, keyg_l
   integer,dimension(:),allocatable :: fragment_atom_id, fragment_supfun_id
   !!logical,dimension(-lmax:lmax) :: file_present
   real(kind=8),dimension(:),pointer :: eval_ptr
   real(kind=8),dimension(:,:),pointer :: coeff_ptr
   real(kind=8),dimension(:),allocatable :: eval, energy_arr, occups, darr, valarr
   real(kind=8),dimension(:,:),allocatable :: denskernel, pdos, occup_arr, hamiltonian_tmp, ovrlp_tmp, matrix_tmp
   real(kind=8),dimension(:,:),allocatable :: kernel_fragment, overlap_fragment, ksk_fragment, tmpmat
   logical,dimension(:,:),allocatable :: calc_array
   logical :: file_exists, found, found_a_fragment, found_icol, found_irow, only_evals, only_binned_values
   logical,dimension(3) :: periodic
   type(matrices) :: ovrlp_mat, hamiltonian_mat, kernel_mat, mat, ovrlp_large, KS_large, kernel_ortho
   type(matrices),dimension(1) :: ovrlp_minus_one_half
   !type(matrices),dimension(:,:),allocatable :: multipoles_matrices
   type(sparse_matrix) :: smat_s, smat_m, smat_l
   type(sparse_matrix),dimension(1) :: smat_arr1
   type(sparse_matrix),dimension(2) :: smat
   type(dictionary), pointer :: dict_timing_info
   integer :: iunit, iat, iat_prev, ii, iitype, iorb, itmb, itype, ival, ios, ipdos, ispin
   integer :: jtmb, norbks, npdos, npt, ntmb, jjtmb, nat_frag, nfvctr_frag, i, iiat
   integer :: icol, irow, icol_atom, irow_atom, iseg, iirow, iicol, j, ifrag, index_dot, ihomo_state, ieval
   integer :: scalapack_blocksize, nbin
   integer(f_long) :: iil, jjl, nnl
   character(len=128),dimension(:),allocatable :: pdos_name, fragment_atomnames
   real(kind=8),dimension(3) ::  ri, rj
   character(len=2) :: backslash, num
   real(kind=8) :: energy, occup, occup_pdos, total_occup, fscale, factor, scale_value, shift_value
   real(kind=8) :: maxdiff, meandiff, tt, tracediff, totdiff
   real(kind=8) :: homo_value, lumo_value, smallest_value, largest_value, gap, gap_target, actual_eval
   real(kind=8) :: mult_factor, add_shift, kT, eTS, accuracy_entropy, d, get_minimal_distance, val, val_mean, binwidth
   real(kind=8) :: dmin, dmax, dstart, dend
   real(mp),dimension(:),allocatable :: eval_min, eval_max
   type(f_progress_bar) :: bar
   integer,parameter :: ncolors = 12
   character(len=1024) :: outfile_base, outfile_extension, matrix_format
   real(kind=8),dimension(:),allocatable :: wa1, wa2, wa3, wa4, qtfa, diag, fvec, weights, qtf, purity_indicator_list
   real(kind=8),dimension(:,:),allocatable :: fjac
   integer,dimension(:),allocatable :: ipvt, natlist, natfrag
   integer,dimension(:,:),allocatable :: fraglist
   real(kind=8) :: epsfcn, gtol, xtol, ftol, charge, purity_indicator
   logical,dimension(:),allocatable :: in_frag, nat_in_frag, at_assigned
   integer :: maxfev, nfev, info, iflag, maxfrag, iel, ifound, nfrag, nfragmax, nfragfound
   integer :: frag_elements_start, frag_elements_end, l
   logical :: weights_integer, atom_complete
   !external :: fcn_fragments
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
!!real(kind=8),dimension(12,12) :: K = reshape([ &
!!   6.938149198647E-01*0.5d0, & !  6.5054268884712541E-01, &
!!   3.272295963759E-01*0.5d0, & !  3.2949680805585824E-01, &
!!  -8.232598245090E-01*0.5d0, & ! -7.9744818038650966E-01, &
!!  -3.374415636414E-01*0.5d0, & ! -3.3662077658064471E-01, &
!!  -5.623089185219E-02*0.5d0, & ! -7.5980028127545063E-02, &
!!  -6.528366322909E-02*0.5d0, & ! -7.3684217368715610E-02, &
!!   7.242523061010E-04*0.5d0, & ! -1.2512277871330202E-03, &
!!   2.872021888738E-03*0.5d0, & !  2.2569388373053366E-03, &
!!   1.778886557592E-03*0.5d0, & !  4.8311346900689316E-03, &
!!   3.269024530123E-03*0.5d0, & !  3.3868168991838957E-03, &
!!  -9.921706166357E-03*0.5d0, & ! -7.1606382026126979E-03, &
!!  -1.733865002582E-03*0.5d0, & ! -8.5501092318784650E-04, &
!!   3.272295963759E-01*0.5d0, & !  3.2949680805585824E-01, &
!!   1.828160565676E+00*0.5d0, & !  1.8142766035698217E+00, &
!!   1.531491157888E-01*0.5d0, & !  1.2573342964474479E-01, &
!!   9.590253229407E-02*0.5d0, & !  8.3026504183287608E-02, &
!!  -2.062748255364E-01*0.5d0, & ! -1.6221890701920783E-01, &
!!   3.618556568157E-01*0.5d0, & !  3.7054721394408252E-01, &
!!   1.201561640500E-04*0.5d0, & ! -1.2194754811245149E-04, &
!!  -8.656222443658E-03*0.5d0, & ! -4.9981248726740143E-03, &
!!  -5.384348451516E-03*0.5d0, & ! -4.7672118220079603E-03, &
!!  -8.418414296216E-03*0.5d0, & ! -7.0690308032024703E-03, &
!!   2.692708235747E-02*0.5d0, & !  2.1147900609865499E-02, &
!!   5.498002569627E-03*0.5d0, & !  3.7789743870511642E-03, &
!!  -8.232598245090E-01*0.5d0, & ! -7.9744818038650966E-01, &
!!   1.531491157888E-01*0.5d0, & !  1.2573342964474479E-01, &
!!   1.449739573264E+00*0.5d0, & !  1.4030428147089755E+00, &
!!  -2.059717562278E-01*0.5d0, & ! -1.3654876198174021E-01, &
!!  -1.656261339265E-01*0.5d0, & ! -1.0440414552220369E-01, &
!!   1.629979868882E-01*0.5d0, & !  1.8099981938571919E-01, &
!!  -1.660077632940E-03*0.5d0, & !  2.1053816511540126E-03, &
!!  -2.090230780194E-03*0.5d0, & ! -5.6860832890118511E-04, &
!!  -3.222215909655E-03*0.5d0, & ! -1.0051982793567579E-02, &
!!  -2.527136514363E-03*0.5d0, & !  1.0962495609968344E-03, &
!!   8.556425891157E-03*0.5d0, & !  2.0296568024845965E-03, &
!!   1.046493878568E-03*0.5d0, & ! -2.2742381910128176E-06, &
!!  -3.374415636414E-01*0.5d0, & ! -3.3662077658064471E-01, &
!!   9.590253229407E-02*0.5d0, & !  8.3026504183287608E-02, &
!!  -2.059717562278E-01*0.5d0, & ! -1.3654876198174021E-01, &
!!   1.911367448688E+00*0.5d0, & !  1.8964816824154958E+00, &
!!   1.340436977048E-02*0.5d0, & !  7.4380996259336262E-03, &
!!  -6.055808797439E-02*0.5d0, & ! -4.9894477672102464E-02, &
!!  -1.375815893670E-03*0.5d0, & !  1.4217061236529552E-04, &
!!   3.377611075509E-03*0.5d0, & !  6.0156170474603143E-03, &
!!   3.708146810445E-04*0.5d0, & !  2.6008616037529937E-03, &
!!   3.666160354128E-03*0.5d0, & ! -6.5022284714939069E-03, &
!!  -3.192439992136E-03*0.5d0, & ! -5.5945806147258441E-03, &
!!  -6.072425610887E-03*0.5d0, & ! -7.8622994567272425E-03, &
!!  -5.623089185219E-02*0.5d0, & ! -7.5980028127545063E-02, &
!!  -2.062748255364E-01*0.5d0, & ! -1.6221890701920783E-01, &
!!  -1.656261339265E-01*0.5d0, & ! -1.0440414552220369E-01, &
!!   1.340436977048E-02*0.5d0, & !  7.4380996259336262E-03, &
!!   1.457006798880E+00*0.5d0, & !  1.4183912699351939E+00, &
!!   8.444350854373E-01*0.5d0, & !  8.1878263380436023E-01, &
!!   8.661995700227E-03*0.5d0, & !  1.3683480430607777E-03, &
!!  -2.484291264832E-02*0.5d0, & ! -4.5540309103186849E-02, &
!!  -7.721767373712E-03*0.5d0, & ! -8.6840733871503956E-03, &
!!  -2.301031396063E-02*0.5d0, & ! -2.9176672522041612E-02, &
!!   5.265768990670E-02*0.5d0, & !  9.3406154473348435E-02, &
!!   2.640320090742E-02*0.5d0, & !  2.9912070463652879E-02, &
!!  -6.528366322909E-02*0.5d0, & ! -7.3684217368715610E-02, &
!!   3.618556568157E-01*0.5d0, & !  3.7054721394408252E-01, &
!!   1.629979868882E-01*0.5d0, & !  1.8099981938571919E-01, &
!!  -6.055808797439E-02*0.5d0, & ! -4.9894477672102464E-02, &
!!   8.444350854373E-01*0.5d0, & !  8.1878263380436023E-01, &
!!   6.683336108225E-01*0.5d0, & !  6.3247478790462075E-01, &
!!  -3.532808259313E-03*0.5d0, & ! -5.6303402248811217E-03, &
!!   3.067540132742E-02*0.5d0, & !  1.0711374355425686E-05, &
!!   1.740960169147E-02*0.5d0, & !  8.5993223276595553E-03, &
!!   2.786817284489E-02*0.5d0, & !  1.1108721075642548E-02, &
!!  -9.570937630381E-02*0.5d0, & ! -2.8349562302968374E-02, &
!!  -1.721101354245E-02*0.5d0, & ! -3.0935575444983382E-03, &
!!   7.242523061010E-04*0.5d0, & ! -1.2512277871330202E-03, &
!!   1.201561640500E-04*0.5d0, & ! -1.2194754811245149E-04, &
!!  -1.660077632940E-03*0.5d0, & !  2.1053816511540126E-03, &
!!  -1.375815893670E-03*0.5d0, & !  1.4217061236529552E-04, &
!!   8.661995700227E-03*0.5d0, & !  1.3683480430607777E-03, &
!!  -3.532808259313E-03*0.5d0, & ! -5.6303402248811217E-03, &
!!   6.775452402912E-01*0.5d0, & !  6.3891534399201377E-01, &
!!   3.477296941647E-01*0.5d0, & !  3.5455479907524345E-01, &
!!  -8.439439333169E-01*0.5d0, & ! -8.1310158325536031E-01, &
!!  -2.421804759560E-01*0.5d0, & ! -2.5307523261966697E-01, &
!!  -2.473013933785E-03*0.5d0, & ! -2.0451934147435610E-02, &
!!  -6.415160018399E-02*0.5d0, & ! -6.8451734329208674E-02, &
!!   2.872021888738E-03*0.5d0, & !  2.2569388373053366E-03, &
!!  -8.656222443658E-03*0.5d0, & ! -4.9981248726740143E-03, &
!!  -2.090230780194E-03*0.5d0, & ! -5.6860832890118511E-04, &
!!   3.377611075509E-03*0.5d0, & !  6.0156170474603143E-03, &
!!  -2.484291264832E-02*0.5d0, & ! -4.5540309103186849E-02, &
!!   3.067540132742E-02*0.5d0, & !  1.0711374355425686E-05, &
!!   3.477296941647E-01*0.5d0, & !  3.5455479907524345E-01, &
!!   1.824505308778E+00*0.5d0, & !  1.8119342627076809E+00, &
!!   2.038691374048E-01*0.5d0, & !  1.6381227320617225E-01, &
!!  -8.191619821176E-02*0.5d0, & ! -6.9703698558694530E-02, &
!!  -1.641304758174E-01*0.5d0, & ! -1.2993163626536569E-01, &
!!   3.498660418079E-01*0.5d0, & !  3.5846246321685293E-01, &
!!   1.778886557592E-03*0.5d0, & !  4.8311346900689316E-03, &
!!  -5.384348451516E-03*0.5d0, & ! -4.7672118220079603E-03, &
!!  -3.222215909655E-03*0.5d0, & ! -1.0051982793567579E-02, &
!!   3.708146810445E-04*0.5d0, & !  2.6008616037529937E-03, &
!!  -7.721767373712E-03*0.5d0, & ! -8.6840733871503956E-03, &
!!   1.740960169147E-02*0.5d0, & !  8.5993223276595553E-03, &
!!  -8.439439333169E-01*0.5d0, & ! -8.1310158325536031E-01, &
!!   2.038691374048E-01*0.5d0, & !  1.6381227320617225E-01, &
!!   1.457419959568E+00*0.5d0, & !  1.4132133557845523E+00, &
!!  -1.855646470249E-01*0.5d0, & ! -1.2198143352104893E-01, &
!!  -3.539189941822E-02*0.5d0, & ! -2.3593222798839068E-02, &
!!   2.950885123566E-02*0.5d0, & !  5.1252005652155611E-02, &
!!   3.269024530123E-03*0.5d0, & !  3.3868168991838957E-03, &
!!  -8.418414296216E-03*0.5d0, & ! -7.0690308032024703E-03, &
!!  -2.527136514363E-03*0.5d0, & !  1.0962495609968344E-03, &
!!   3.666160354128E-03*0.5d0, & ! -6.5022284714939069E-03, &
!!  -2.301031396063E-02*0.5d0, & ! -2.9176672522041612E-02, &
!!   2.786817284489E-02*0.5d0, & !  1.1108721075642548E-02, &
!!  -2.421804759560E-01*0.5d0, & ! -2.5307523261966697E-01, &
!!  -8.191619821176E-02*0.5d0, & ! -6.9703698558694530E-02, &
!!  -1.855646470249E-01*0.5d0, & ! -1.2198143352104893E-01, &
!!   1.703074846814E+00*0.5d0, & !  1.6677627145616152E+00, &
!!  -2.885435063655E-01*0.5d0, & ! -1.8701748632621470E-01, &
!!   5.667671718127E-01*0.5d0, & !  5.6902619024084800E-01, &
!!  -9.921706166357E-03*0.5d0, & ! -7.1606382026126979E-03, &
!!   2.692708235747E-02*0.5d0, & !  2.1147900609865499E-02, &
!!   8.556425891157E-03*0.5d0, & !  2.0296568024845965E-03, &
!!  -3.192439992136E-03*0.5d0, & ! -5.5945806147258441E-03, &
!!   5.265768990669E-02*0.5d0, & !  9.3406154473348435E-02, &
!!  -9.570937630381E-02*0.5d0, & ! -2.8349562302968374E-02, &
!!  -2.473013933785E-03*0.5d0, & ! -2.0451934147435610E-02, &
!!  -1.641304758174E-01*0.5d0, & ! -1.2993163626536569E-01, &
!!  -3.539189941822E-02*0.5d0, & ! -2.3593222798839068E-02, &
!!  -2.885435063655E-01*0.5d0, & ! -1.8701748632621470E-01, &
!!   1.658977490784E+00*0.5d0, & !  1.6195295556214508E+00, &
!!   6.644198725245E-01*0.5d0, & !  6.5469534647209005E-01, &
!!  -1.733865002582E-03*0.5d0, & ! -8.5501092318784650E-04, &
!!   5.498002569627E-03*0.5d0, & !  3.7789743870511642E-03, &
!!   1.046493878568E-03*0.5d0, & ! -2.2742381910128176E-06, &
!!  -6.072425610887E-03*0.5d0, & ! -7.8622994567272425E-03, &
!!   2.640320090742E-02*0.5d0, & !  2.9912070463652879E-02, &
!!  -1.721101354245E-02*0.5d0, & ! -3.0935575444983382E-03, &
!!  -6.415160018399E-02*0.5d0, & ! -6.8451734329208674E-02, &
!!   3.498660418079E-01*0.5d0, & !  3.5846246321685293E-01, &
!!   2.950885123566E-02*0.5d0, & !  5.1252005652155611E-02, &
!!   5.667671718127E-01*0.5d0, & !  5.6902619024084800E-01, &
!!   6.644198725245E-01*0.5d0, & !  6.5469534647209005E-01, &
!!   6.700542365811E-01*0.5d0],& !  6.3058157710140161E-01],&
!!  [12,12])

   call f_lib_initialize()

   ! Initialize MPI
   !call bigdft_mpi_init(ierr)
   !!call bigdft_init()
   ! MPI initialization; we have:
   ! iproc is the task ID
   ! nproc is the total number of tasks
   call mpiinit()
   iproc=mpirank()
   nproc=mpisize()

   call f_malloc_set_status(iproc=iproc)

   ! Initialize the sparse matrix errors and timings.
   call sparsematrix_init_errors
   call sparsematrix_initialize_timing_categories()


    if (iproc==0) then
        call yaml_new_document()
        !call print_logo()
    end if

   !Time initialization
   call f_timing_reset(filename='time.yaml',master=(iproc==0),verbose_mode=.false.)

   if (iproc==0) then
       call yaml_scalar('',hfill='~')
       call yaml_scalar('CHESS TOOLBOX',hfill='~')
   end if

   if (iproc==0) then
       call yaml_map('Timestamp of the run',yaml_date_and_time_toa())
       call yaml_mapping_open('Parallel environment')
       call yaml_map('MPI tasks',nproc)
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
         !!else if (trim(tatonam)=='multipole-analysis') then
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = method_name)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = matrix_basis)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = metadata_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = overlap_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = kernel_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = hamiltonian_file)
         !!   do l=0,lmax
         !!       do m=-l,l
         !!           i_arg = i_arg + 1
         !!           call get_command_argument(i_arg, value = multipoles_files(m,l))
         !!       end do
         !!   end do
         !!   !write(*,'(1x,2a)')&
         !!   !   &   'perform a Loewdin charge analysis'
         !!   multipole_analysis = .true.
         !!   exit loop_getargs
         else if (trim(tatonam)=='solve-eigensystem') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = itype_)
            read(itype_,fmt=*,iostat=ierr) itype
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_format)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = hamiltonian_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = coeff_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = eval_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = only_evals_)
            read(only_evals_,fmt=*,iostat=ierr) only_evals
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = scalapack_blocksize_)
            read(scalapack_blocksize_,fmt=*,iostat=ierr) scalapack_blocksize
            !write(*,'(1x,2a)')&
            !   &   'perform a Loewdin charge analysis'
            solve_eigensystem = .true.
            exit loop_getargs
         else if (trim(tatonam)=='pdos') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_format)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = coeff_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = hamiltonian_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = pdos_file)
            calculate_pdos = .true.
        else if (trim(tatonam)=='convert-matrix-format') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = conversion)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = infile)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = outfile)
            convert_matrix_format = .true.
        else if (trim(tatonam)=='calculate-selected-eigenvalues') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_format)
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
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = itype_)
            read(itype_,fmt=*,iostat=ierr) itype
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = iev_min_)
            read(iev_min_,fmt=*,iostat=ierr) iev_min
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = iev_max_)
            read(iev_max_,fmt=*,iostat=ierr) iev_max
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fscale_)
            read(fscale_,fmt=*,iostat=ierr) fscale
            calculate_selected_eigenvalues = .true.
        else if (trim(tatonam)=='kernel-purity') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_format)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_matmul_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fragment_file)
            kernel_purity = .true.
        else if (trim(tatonam)=='manipulate-eigenvalue-spectrum') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = manipulation_mode)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_format)
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
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = hamiltonian_manipulated_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ihomo_state_)
            read(ihomo_state_,fmt=*,iostat=ierr) ihomo_state
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = homo_value_)
            read(homo_value_,fmt=*,iostat=ierr) homo_value
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = lumo_value_)
            read(lumo_value_,fmt=*,iostat=ierr) lumo_value
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = smallest_value_)
            read(smallest_value_,fmt=*,iostat=ierr) smallest_value
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = largest_value_)
            read(largest_value_,fmt=*,iostat=ierr) largest_value
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = diag_algorithm)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = scalapack_blocksize_)
            read(scalapack_blocksize_,fmt=*,iostat=ierr) scalapack_blocksize
            manipulate_eigenvalue_spectrum = .true.
        else if (trim(tatonam)=='calculate-entropy') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_format)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_matmul_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kT_)
            read(kT_,fmt=*,iostat=ierr) kT
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = accuracy_entropy_)
            read(accuracy_entropy_,fmt=*,iostat=ierr) accuracy_entropy
            calculate_entropy = .true.
        else if (trim(tatonam)=='analyze-density-matrix-dense') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_format)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = nbin_)
            read(nbin_,fmt=*,iostat=ierr) nbin
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = output_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = output_bins_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = only_binned_values_)
            read(only_binned_values_,fmt=*,iostat=ierr) only_binned_values
            analyze_density_matrix_dense = .true.
        else if (trim(tatonam)=='find-fragments') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_format)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_matmul_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = frag_elements_start_)
            read(frag_elements_start_,fmt=*,iostat=ierr) frag_elements_start
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = frag_elements_end_)
            read(frag_elements_end_,fmt=*,iostat=ierr) frag_elements_end
            find_fragments = .true.
         end if
         i_arg = i_arg + 1
      end do loop_getargs
   end if


   !!if (multipole_analysis) then
   !!    if (iproc==0) then
   !!        call yaml_comment('Multipole analysis',hfill='-')
   !!    end if

   !!    ! Determine the method
   !!    select case(trim(method_name))
   !!    case ('loewdin','LOEWDIN')
   !!        method = CHARGE_ANALYSIS_LOEWDIN
   !!    case ('mulliken','MULLIKEN')
   !!        method = CHARGE_ANALYSIS_MULLIKEN
   !!    !!case ('projector','PROJECTOR')
   !!    !!    method = CHARGE_ANALYSIS_PROJECTOR
   !!    case default
   !!        call f_err_throw('Unknown Method for the multipole analysis',err_name='BIGDFT_INPUT_VARIABLES_ERROR')
   !!    end select

   !!    call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
   !!    if (iproc==0) then
   !!        call yaml_mapping_open('Atomic System Properties')
   !!        call yaml_map('Types of atoms',smmd%atomnames)
   !!        call yaml_mapping_close()
   !!    end if

   !!    call sparse_matrix_and_matrices_init_from_file_bigdft('serial_text', trim(overlap_file), &
   !!         iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
   !!         init_matmul=.true.)

   !!    call sparse_matrix_and_matrices_init_from_file_bigdft('serial_text', trim(kernel_file), &
   !!         iproc, nproc, mpiworld(), smat_l, kernel_mat, &
   !!         init_matmul=.true.)

   !!    call sparse_matrix_and_matrices_init_from_file_bigdft('serial_text', trim(hamiltonian_file), &
   !!         iproc, nproc, mpiworld(), smat_m, hamiltonian_mat, &
   !!         init_matmul=.true.)

   !!    ! Check which multipole matrices are present
   !!    do l=0,lmax
   !!        file_present(:) = .false.
   !!        do m=-l,l
   !!            inquire(file=trim(multipoles_files(m,l)), exist=file_exists)
   !!            if (file_exists) then
   !!                file_present(m) = .true.
   !!            end if
   !!        end do
   !!        if (any(file_present(-l:l))) then
   !!            if (.not.all(file_present(-l:l))) then
   !!                call f_err_throw('for a given shell all matrices must be present', &
   !!                     err_name='BIGDFT_RUNTIME_ERROR')
   !!            end if
   !!            ll = l
   !!        end if
   !!    end do

   !!    allocate(multipoles_matrices(-ll:ll,0:ll))
   !!    do l=0,ll
   !!        do m=-l,l
   !!            call matrices_init_from_file_bigdft(trim(multipoles_files(m,l)), &
   !!                 iproc, nproc, mpiworld(), smat_s, multipoles_matrices(m,l))
   !!        end do
   !!    end do

   !!    call timing(mpiworld(),'INIT','PR')

   !!    select case(method)
   !!    case (CHARGE_ANALYSIS_MULLIKEN)
   !!        methodc='loewdin'
   !!        do_ortho='no'
   !!    case (CHARGE_ANALYSIS_LOEWDIN)
   !!        methodc='loewdin'
   !!        do_ortho='yes'
   !!    !!case (CHARGE_ANALYSIS_PROJECTOR)
   !!    !!    methodc='projector'
   !!    !!    do_ortho='no'
   !!    case default
   !!        call f_err_throw('wrong method',err_name='BIGDFT_RUNTIME_ERROR')
   !!    end select

   !!    ! Determine whether the multipoles matrices to be used are calculated analytically or on the grid.
   !!    ! For the analytic case, only the overlap matrix is possible.
   !!    select case(trim(matrix_basis))
   !!    case ('wavelet','WAVELET')
   !!        call multipole_analysis_driver_new(iproc, nproc, 0, 11, &
   !!             smmd, smat_s, smat_m, smat_l, &
   !!             ovrlp_mat, hamiltonian_mat, kernel_mat, smmd%rxyz, &
   !!             methodc, do_ortho=trim(do_ortho), projectormode='simple', &
   !!             calculate_multipole_matrices=.false., do_check=.false., &
   !!             write_multipole_matrices=.false., &
   !!             multipole_matrix_in=(/(/ovrlp_mat/)/))
   !!    case ('realspace','REALSPACE')
   !!        call multipole_analysis_driver_new(iproc, nproc, ll, 11, &
   !!             smmd, smat_s, smat_m, smat_l, &
   !!             ovrlp_mat, hamiltonian_mat, kernel_mat, smmd%rxyz, &
   !!             methodc, do_ortho=trim(do_ortho), projectormode='simple', &
   !!             calculate_multipole_matrices=.false., do_check=.false., &
   !!             write_multipole_matrices=.false., &
   !!             multipole_matrix_in=multipoles_matrices)
   !!    case default
   !!        call f_err_throw('wrong value for matrix_basis',err_name='BIGDFT_RUNTIME_ERROR')
   !!    end select

   !!    call timing(mpiworld(),'CALC','PR')

   !!    call deallocate_sparse_matrix_metadata(smmd)
   !!    call deallocate_sparse_matrix(smat_s)
   !!    call deallocate_sparse_matrix(smat_l)
   !!    call deallocate_matrices(ovrlp_mat)
   !!    call deallocate_matrices(kernel_mat)
   !!    call deallocate_sparse_matrix(smat_m)
   !!    call deallocate_matrices(hamiltonian_mat)
   !!    do l=0,ll
   !!        do m=-l,l
   !!            call deallocate_matrices(multipoles_matrices(m,l))
   !!        end do
   !!    end do
   !!    deallocate(multipoles_matrices)

   !!    if (iproc==0) then
   !!        call yaml_comment('done',hfill='-')
   !!    end if

   !!end if


   if (solve_eigensystem) then
       if (iproc==0) then
           if (only_evals) then
               i = 1
           else
               i = 0
           end if
       end if
       call fmpi_bcast(i)!, root=0, comm=mpi_comm_world)
       if (i==1) then
           only_evals = .true.
       else
           only_evals = .false.
       end if
       call solve_eigensystem_lapack(iproc, nproc, mpi_comm_world, itype, matrix_format, metadata_file, &
            overlap_file, hamiltonian_file, scalapack_blocksize, write_coeff=.not.only_evals, write_eval=only_evals, &
            coeff_file=trim(coeff_file), eval_file=trim(eval_file))

       !!!if (iproc==0) call yaml_comment('Reading from file '//trim(overlap_file),hfill='~')
       !!call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
       !!     iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
       !!     init_matmul=.false.)!, nat=nat, rxyz=rxyz, iatype=iatype, ntypes=ntypes, &
       !!     !nzatom=nzatom, nelpsp=nelpsp, atomnames=atomnames)
       !!call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       !!!if (iproc==0) call yaml_comment('Reading from file '//trim(hamiltonian_file),hfill='~')
       !!call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(hamiltonian_file), &
       !!     iproc, nproc, mpiworld(), smat_m, hamiltonian_mat, &
       !!     init_matmul=.false.)

       !!ovrlp_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='ovrlp_mat%matrix')
       !!call uncompress_matrix(iproc, nproc, &
       !!     smat_s, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_mat%matrix)
       !!hamiltonian_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='hamiltonian_mat%matrix')
       !!call uncompress_matrix(iproc, nproc, &
       !!     smat_m, inmat=hamiltonian_mat%matrix_compr, outmat=hamiltonian_mat%matrix)
       !!eval = f_malloc(smat_s%nfvctr,id='eval')

       !!if (iproc==0) then
       !!    call yaml_comment('Diagonalizing the matrix',hfill='~')
       !!end if
       !!call diagonalizeHamiltonian2(iproc, nproc, mpiworld(), scalapack_blocksize, &
       !!     smat_s%nfvctr, hamiltonian_mat%matrix, ovrlp_mat%matrix, eval)
       !!if (iproc==0) then
       !!    call yaml_comment('Matrix successfully diagonalized',hfill='~')
       !!end if
       !!iunit=99
       !!call f_open_file(iunit, file=trim(coeff_file), binary=.false.)
       !!!call writeLinearCoefficients(iunit, .true., nat, rxyz, smat_s%nfvctr, smat_s%nfvctr, &
       !!!     smat_s%nfvctr, hamiltonian_mat%matrix, eval)
       !!call write_linear_coefficients(iproc, 0, trim(coeff_file), 2, smmd%nat, smmd%rxyz, &
       !!     smmd%iatype, smmd%ntypes, smmd%nzatom, &
       !!     smmd%nelpsp, smmd%atomnames, smat_s%nfvctr, &
       !!     smat_s%nfvctr, smat_s%nspin, hamiltonian_mat%matrix, eval)
       !!call f_close(iunit)

       !!call f_free(eval)
       !!call deallocate_matrices(ovrlp_mat)
       !!call deallocate_matrices(hamiltonian_mat)
       !!call deallocate_sparse_matrix(smat_s)
       !!call deallocate_sparse_matrix(smat_m)
       !!call deallocate_sparse_matrix_metadata(smmd)
       !call f_free_ptr(rxyz)
       !call f_free_ptr(iatype)
       !call f_free_ptr(nzatom)
       !call f_free_ptr(nelpsp)
       !call f_free_str_ptr(len(atomnames),atomnames)
   end if


   if (calculate_pdos) then
       iunit = 99
       !if (iproc==0) call yaml_comment('Reading from file '//trim(coeff_file),hfill='~')
       call f_open_file(iunit, file=trim(coeff_file), binary=.false.)
       call read_linear_coefficients('serial_text', iproc, nproc, mpiworld(), &
            trim(coeff_file), nspin, ntmb, norbks, coeff_ptr, &
            eval=eval_ptr)
       call f_close(iunit)

       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)

       !if (iproc==0) call yaml_comment('Reading from file '//trim(overlap_file),hfill='~')
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
            iproc, nproc, mpiworld(), smat(1), ovrlp_mat, &
            init_matmul=.false.)!, iatype=iatype, ntypes=ntypes, atomnames=atomnames, &
            !on_which_atom=on_which_atom)

       !if (iproc==0) call yaml_comment('Reading from file '//trim(hamiltonian_file),hfill='~')
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(hamiltonian_file), &
            iproc, nproc, mpiworld(), smat(2), hamiltonian_mat, &
            init_matmul=.false.)

       call init_matrix_taskgroups_wrapper(iproc, nproc, mpi_comm_world, .true., 2, smat)

       call resize_matrix_to_taskgroup(smat(1), ovrlp_mat)
       call resize_matrix_to_taskgroup(smat(2), hamiltonian_mat)

       if (iproc==0) then
           call yaml_mapping_open('Matrix properties')
           call write_sparsematrix_info(smat(1), 'Overlap matrix')
           call write_sparsematrix_info(smat(2), 'Hamiltonian matrix')
           call yaml_mapping_close()
       end if



       if (ntmb/=smat(1)%nfvctr) call f_err_throw('ntmb/=smat(1)%nfvctr')
       ovrlp_mat%matrix = sparsematrix_malloc_ptr(smat(1), iaction=DENSE_PARALLEL, id='ovrlp_mat%matrix')
       hamiltonian_mat%matrix = sparsematrix_malloc_ptr(smat(1), iaction=DENSE_PARALLEL, id='hamiltonian_mat%matrix')
       call uncompress_matrix_distributed2(iproc, smat(1), DENSE_PARALLEL, &
            ovrlp_mat%matrix_compr, ovrlp_mat%matrix(1:,1:,1))
       call uncompress_matrix_distributed2(iproc, smat(2), DENSE_PARALLEL, &
            hamiltonian_mat%matrix_compr, hamiltonian_mat%matrix(1:,1:,1))


       if (iproc==0) then
           call yaml_mapping_open('Matrix properties')
           call write_sparsematrix_info(smat(1), 'Overlap matrix')
           call write_sparsematrix_info(smat_l, 'Density kernel')
           call yaml_mapping_close()
       end if

       iunit=99
       call f_open_file(iunit, file=pdos_file, binary=.false.)
       read(iunit,*) npdos

       calc_array = f_malloc((/ntmb,npdos/),id='calc_array')
       pdos_name = f_malloc_str(len(pdos_name),npdos,id='pdos_name')


       do ipdos=1,npdos
           do itmb=1,ntmb
               calc_array(itmb,ipdos) = .false.
           end do
       end do
       ipdos = 0
       !npdos_loop: do !ipdos=1,npdos
       if (iproc==0) then
           call yaml_sequence_open('Atoms and support functions to be taken into account for each partial density of states')
       end if
           do 
               !read(iunit01,*,iostat=ios) cc, ival
               read(iunit,'(a128)',iostat=ios) line
               if (ios/=0) exit
               !write(*,*) 'line',line
               read(line,*,iostat=ios) cc, cname
               if (cc=='#') then
                   ipdos = ipdos + 1
                   pdos_name(ipdos) = trim(cname)
                   if (iproc==0) then
                       if (ipdos>1) then
                           call yaml_mapping_close()
                       end if
                       call yaml_sequence(advance='no')
                       call yaml_mapping_open(trim(pdos_name(ipdos)))
                   end if
                   cycle 
               end if
               read(line,*,iostat=ios) cc, ival
               if (iproc==0) then
                   call yaml_map(trim(cc),ival)
               end if
               found = .false.
               do itype=1,smmd%ntypes
                   if (trim(smmd%atomnames(itype))==trim(cc)) then
                       iitype = itype
                       found = .true.
                       exit
                   end if
               end do
               if (.not.found) then
                   call f_err_throw("Atom type '"//trim(cc)//"' is unknown", &
                        err_name='SPARSEMATRIX_RUNTIME_ERROR')
               end if
               iat_prev = -1
               do itmb=1,ntmb
                   iat = smmd%on_which_atom(itmb)
                   if (iat/=iat_prev) then
                       ii = 0
                   end if
                   iat_prev = iat
                   itype = smmd%iatype(iat)
                   ii = ii + 1
                   if (itype==iitype .and. ii==ival) then
                       if (calc_array(itmb,ipdos)) then
                           call f_err_throw('calc_array(itmb) must not be true here', &
                                err_name='SPARSEMATRIX_RUNTIME_ERROR')
                       end if
                       calc_array(itmb,ipdos) = .true.
                   end if
               end do
           end do
       if (iproc==0) then
           call yaml_mapping_close()
           call yaml_sequence_close()
       end if

       !end do npdos_loop
       call f_close(iunit)


       energy_arr = f_malloc0(norbks,id='energy_arr')
       occup_arr = f_malloc0((/npdos,norbks/),id='occup_arr')
       occups = f_malloc0(npdos,id='occups')

       denskernel = f_malloc((/ntmb,smat(1)%nfvctrp/),id='denskernel')
       pdos = f_malloc0((/npt,npdos/),id='pdos')
       if (iproc==0) then
           call yaml_comment('PDoS calculation',hfill='~')
           call yaml_mapping_open('Calculating PDoS')
       end if
       ! Calculate a partial kernel for each KS orbital
       !do ipdos=1,npdos
           !if (iproc==0) call yaml_map('PDoS number',ipdos)
           !if (iproc==0) call yaml_map('start, increment',(/ipdos,npdos/))
           !write(num,fmt='(i2.2)') ipdos
           if (iproc==0) bar=f_progress_bar_new(nstep=norbks)
           do iorb=1,norbks
               if (iproc==0) call yaml_map('orbital being processed',iorb)
               call gemm('n', 't', ntmb, smat(1)%nfvctrp, 1, 1.d0, coeff_ptr(1,iorb), ntmb, &
                    coeff_ptr(smat(1)%isfvctr+1,iorb), ntmb, 0.d0, denskernel(1,1), ntmb)
               energy = 0.d0
               call f_zero(occups)
               do ispin=1,nspin
                   !$omp parallel default(none) &
                   !$omp shared(ispin,smat,ntmb,denskernel,hamiltonian_mat,energy) &
                   !$omp private(itmb,jtmb,jjtmb)
                   !$omp do reduction( + : energy)
                   do jtmb=1,smat(1)%nfvctrp
                       jjtmb = smat(1)%isfvctr + jtmb
                       do itmb=1,ntmb
                           energy = energy + denskernel(itmb,jtmb)*hamiltonian_mat%matrix(itmb,jtmb,ispin)
                       end do
                   end do
                   !$omp end do
                   !$omp end parallel
                   energy_arr(iorb) = energy
               end do
               !call fmpi_allreduce(energy, 1, FMPI_SUM, comm=mpiworld())
               do ispin=1,nspin
!                   !$omp parallel default(none) &
!                   !$omp shared(ispin,smat,ntmb,denskernel,ovrlp_mat,calc_array,npdos,occups) &
!                   !$omp private(itmb,jtmb,jjtmb,occup,ipdos)
!                   !$omp do reduction(+:occups)
                   do jtmb=1,smat(1)%nfvctrp
                       jjtmb = smat(1)%isfvctr + jtmb
                       !if (calc_array(jjtmb,ipdos)) then
                           do itmb=1,ntmb!ipdos,ntmb,npdos
                               !if (calc_array(itmb,ipdos)) then
                                   occup = denskernel(itmb,jtmb)*ovrlp_mat%matrix(itmb,jtmb,ispin)
                                   do ipdos=1,npdos
                                       if (calc_array(itmb,ipdos) .and. calc_array(jjtmb,ipdos)) then
                                           occups(ipdos) = occups(ipdos) + occup
                                       end if
                                   end do
                                   !occup = occup + denskernel(itmb,jtmb)*ovrlp_mat%matrix(itmb,jtmb,ispin)
                               !end if
                           end do
                       !end if
                   end do
!                   !$omp end do
!                   !$omp end parallel
                   do ipdos=1,npdos
                       occup_arr(ipdos,iorb) = occups(ipdos)
                   end do
               end do
               !occup_pdos = occup_pdos + occup
               !!if (iorb<norbks) then
               !!    write(iunit,'(2(a,es16.9),a)') '  ',occup,'*exp(-(x-',energy,')**2/(2*sigma**2)) + '//trim(backslash)
               !!else
               !!    write(iunit,'(2(a,es16.9),a)') '  ',occup,'*exp(-(x-',energy,')**2/(2*sigma**2))'
               !!end if
               !!!write(*,'(a,i6,3es16.8)')'iorb, eval(iorb), energy, occup', iorb, eval(iorb), energy, occup
               if (iproc==0) call dump_progress_bar(bar,step=iorb)
           end do
           !total_occup = total_occup + occup_pdos
           !!if (ipdos==1) then
           !!    write(iunit,'(a,i0,a)') "plot f",ipdos,"(x) lc rgb 'color' lt 1 lw 2 w l title 'name'"
           !!else
           !!    write(iunit,'(a,i0,a)') "replot f",ipdos,"(x) lc rgb 'color' lt 1 lw 2 w l title 'name'"
           !!end if
           !!if (iproc==0) call yaml_map('sum of PDoS',occup_pdos)
           !!output_pdos='PDoS_'//num//'.dat'
           !!call yaml_map('output file',trim(output_pdos))
           !!call f_open_file(iunit01, file=trim(output_pdos), binary=.false.)
           !!write(iunit01,'(a)') '#             energy                pdos'
           !!do ipt=1,npt
           !!    write(iunit01,'(2es20.12)') energy_bins(1,ipt), pdos(ipt,ipdos)
           !!end do
           !!call f_close(iunit01)
       !end do
       call fmpi_allreduce(occup_arr, FMPI_SUM, comm=mpiworld())
       call fmpi_allreduce(energy_arr, FMPI_SUM, comm=mpiworld())
       if (iproc==0) call yaml_mapping_close()

       if (iproc==0) then
           call yaml_comment('Calculation complete',hfill='=')
           output_pdos='PDoS.gp'
           call yaml_map('output file',trim(output_pdos))
           iunit = 99
           call f_open_file(iunit, file=trim(output_pdos), binary=.false.)
           write(iunit,'(a)') '# plot the DOS as a sum of Gaussians'
           write(iunit,'(a)') 'set samples 1000'
           write(iunit,'(a,2(es12.5,a))') 'set xrange[',eval_ptr(1),':',eval_ptr(ntmb),']'
           write(iunit,'(a)') 'sigma=0.01'
           write(backslash,'(a)') '\ '
           total_occup = 0.d0
           do ipdos=1,npdos
               occup_pdos = 0.d0
               call yaml_map('PDoS number',ipdos)
               call yaml_map('start, increment',(/ipdos,npdos/))
               write(num,fmt='(i2.2)') ipdos
               write(iunit,'(a,i0,a)') 'f',ipdos,'(x) = '//trim(backslash)
               do iorb=1,norbks
                   if (iorb<norbks) then
                       write(iunit,'(2(a,es16.9),a)') '  ',occup_arr(ipdos,iorb),&
                           &'*exp(-(x-',energy_arr(iorb),')**2/(2*sigma**2)) + '//trim(backslash)
                   else
                       write(iunit,'(2(a,es16.9),a)') '  ',occup_arr(ipdos,iorb),&
                           &'*exp(-(x-',energy_arr(iorb),')**2/(2*sigma**2))'
                   end if
                   occup_pdos = occup_pdos + occup_arr(ipdos,iorb)
               end do
               total_occup = total_occup + occup_pdos
               call yaml_map('sum of PDoS',occup_pdos)
           end do
           do ipdos=1,npdos
               if (ipdos<ncolors) then
                   colorname = colors(ipdos)
               else
                   colorname = 'color'
               end if
               if (ipdos<npdos) then
                   lineend = ' ,'//trim(backslash)
               else
                   lineend = ''
               end if
               if (ipdos==1) then
                   linestart = 'plot'
               else
                   linestart = '    '
               end if
               write(iunit,'(a,i0,a)') trim(linestart)//" f",ipdos,"(x) lc rgb '"//trim(colorname)//&
                   &"' lt 1 lw 2 w l title '"//trim(pdos_name(ipdos))//"'"//trim(lineend)
           end do
           call f_close(iunit)
       end if
       if (iproc==0) call yaml_map('sum of total DoS',total_occup)
       call fmpi_barrier()

       call f_free(pdos)
       call f_free(denskernel)
       call f_free_ptr(eval_ptr)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(hamiltonian_mat)
       call deallocate_sparse_matrix(smat(1))
       call deallocate_sparse_matrix(smat(2))
       call deallocate_sparse_matrix_metadata(smmd)
       !call f_free_ptr(iatype)
       !call f_free_str_ptr(len(atomnames),atomnames)
       call f_free_str(len(pdos_name),pdos_name)
       call f_free(calc_array)
       !call f_free_ptr(on_which_atom)
       call f_free_ptr(coeff_ptr)
       call f_free(energy_arr)
       call f_free(occup_arr)
       call f_free(occups)
   end if

   if (convert_matrix_format) then
       select case (trim(conversion))
       case ('bigdft_to_ccs')
           iconv = 1
       case ('ccs_to_bigdft')
           iconv = 2
       case ('bigdft_to_dense')
           iconv = 3
       case ('binary_to_bigdft')
           iconv = 4
       case ('bigdft_to_binary')
           iconv = 5
       case default
           call f_err_throw("wrong value for conversion; possible are &
               &'bigdft_to_ccs',&
               &'ccs_to_bigdft',&
               &'bigdft_to_dense',&
               &'bigdft_to_binary',&
               &'binary_to_bigdft'")
       end select

       select case (iconv)
       case (1,3,5)
           !if (iproc==0) call yaml_comment('Reading from file '//trim(infile),hfill='~')
           call sparse_matrix_and_matrices_init_from_file_bigdft('serial_text', trim(infile), &
                iproc, nproc, mpiworld(), &
                smat_arr1(1), mat, init_matmul=.false.)
       case (2)
           !if (iproc==0) call yaml_comment('Reading from file '//trim(infile),hfill='~')
           call sparse_matrix_and_matrices_init_from_file_ccs(trim(infile), iproc, nproc, &
                mpiworld(), smat_arr1(1), mat, init_matmul=.false.)
       case(4)
           !if (iproc==0) call yaml_comment('Reading from file '//trim(infile),hfill='~')
           call sparse_matrix_and_matrices_init_from_file_bigdft('parallel_mpi-native', trim(infile), &
                iproc, nproc, mpiworld(), &
                smat_arr1(1), mat, init_matmul=.false.)
       end select
       call init_matrix_taskgroups_wrapper(iproc, nproc, mpiworld(), .false., 1, smat_arr1)
       call resize_matrix_to_taskgroup(smat_arr1(1), mat)
       if (iproc==0) then
           call yaml_mapping_open('Matrix properties')
           call write_sparsematrix_info(smat_arr1(1), 'Input matrix')
           call yaml_mapping_close()
       end if

       select case (iconv)
       case (1)
           row_ind = f_malloc_ptr(smat_arr1(1)%nvctr,id='row_ind')
           col_ptr = f_malloc_ptr(smat_arr1(1)%nfvctr,id='col_ptr')
           call ccs_data_from_sparse_matrix(smat_arr1(1), row_ind, col_ptr)
           if (iproc==0) call ccs_matrix_write(trim(outfile), smat_arr1(1), row_ind, col_ptr, mat)
           call f_free_ptr(row_ind)
           call f_free_ptr(col_ptr)
       case (2,4,5)
           !!call sparse_matrix_and_matrices_init_from_file_ccs(trim(infile), iproc, nproc, &
           !!     mpiworld(), smat_arr1, mat, init_matmul=.false.)
           if (len(outfile)>1024) then
               call f_err_throw('filename is too long')
           end if

           index_dot = index(outfile,'.',back=.true.)
           outfile_base = outfile(1:index_dot-1)
           outfile_extension = outfile(index_dot:)
           select case(iconv)
           case (2,4)
               if (trim(outfile_extension)/='.txt') then
                   call f_err_throw('Wrong file extension; must be .txt, but found '//trim(outfile_extension))
               end if
               call write_sparse_matrix('serial_text', iproc, nproc, mpiworld(), &
                    smat_arr1(1), mat, trim(outfile_base))
           case (5)
               if (trim(outfile_extension)/='.mpi') then
                   call f_err_throw('Wrong file extension; must be .mpi, but found '//trim(outfile_extension))
               end if
               call write_sparse_matrix('parallel_mpi-native', iproc, nproc, mpiworld(), &
                    smat_arr1(1), mat, trim(outfile_base))
           end select
       case (3)
           call write_dense_matrix(iproc, nproc, mpiworld(), smat_arr1(1), mat, &
                uncompress=.true., filename=trim(outfile), binary=.false.)

       end select

       call deallocate_sparse_matrix(smat_arr1(1))
       call deallocate_matrices(mat)
   end if


   if (calculate_selected_eigenvalues) then
       call calculate_eigenvalues(iproc, nproc, matrix_format, metadata_file, &
            overlap_file, hamiltonian_file, kernel_file, kernel_matmul_file, &
            itype, iev_min, iev_max, fscale)
       !!!call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       !!!call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
       !!!     iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
       !!!     init_matmul=.false.)
       !!!call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(hamiltonian_file), &
       !!!     iproc, nproc, mpiworld(), smat_m, hamiltonian_mat, &
       !!!     init_matmul=.false.)
       !!!call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(kernel_file), &
       !!!     iproc, nproc, mpiworld(), smat_l, kernel_mat, &
       !!!     init_matmul=.true., filename_mult=trim(kernel_matmul_file))
       !!!call matrices_init(smat_l, ovrlp_minus_one_half(1))

       !!!!call timing(mpiworld(),'INIT','PR')
       !!!call f_timing_checkpoint(ctr_name='INIT',mpi_comm=mpiworld(),nproc=mpisize(),&
       !!!             gather_routine=gather_timings)
       !!!if (iev_min<1 .or. iev_min>smat_s%nfvctr .or. iev_max>smat_s%nfvctr .or. iev_max<1) then
       !!!    if (iproc==0) then
       !!!        call yaml_warning('The required eigenvalues are outside of the possible range, automatic ajustment')
       !!!    end if
       !!!end if
       !!!iev_min = max(iev_min,1)
       !!!iev_min = min(iev_min,smat_s%nfvctr)
       !!!iev_max = min(iev_max,smat_s%nfvctr)
       !!!iev_max = max(iev_max,1)
       !!!eval = f_malloc(iev_min.to.iev_max,id='eval')
       !!!if (iproc==0) then
       !!!    call yaml_mapping_open('Calculating eigenvalues using FOE')
       !!!end if
       !!!call get_selected_eigenvalues_from_FOE(iproc, nproc, mpiworld(), &
       !!!     iev_min, iev_max, smat_s, smat_m, smat_l, ovrlp_mat, hamiltonian_mat, &
       !!!     ovrlp_minus_one_half, eval, fscale)

       !!!!!call timing(mpiworld(),'CALC','PR')
       !!!call f_timing_checkpoint(ctr_name='CALC',mpi_comm=mpiworld(),nproc=mpisize(),&
       !!!             gather_routine=gather_timings)

       !!!if (iproc==0) then
       !!!    call yaml_sequence_open('values')
       !!!    do iev=iev_min,iev_max
       !!!        call yaml_sequence(advance='no')
       !!!        call yaml_mapping_open(flow=.true.)
       !!!        call yaml_map('ID',iev,fmt='(i6.6)')
       !!!        call yaml_map('eval',eval(iev),fmt='(es12.5)')
       !!!        call yaml_mapping_close()
       !!!    end do
       !!!    call yaml_sequence_close()
       !!!    call yaml_mapping_close()
       !!!end if


       !!!call deallocate_sparse_matrix(smat_s)
       !!!call deallocate_sparse_matrix(smat_m)
       !!!call deallocate_sparse_matrix(smat_l)
       !!!call deallocate_matrices(ovrlp_mat)
       !!!call deallocate_matrices(hamiltonian_mat)
       !!!call deallocate_matrices(kernel_mat)
       !!!call deallocate_matrices(ovrlp_minus_one_half(1))
       !!!call deallocate_sparse_matrix_metadata(smmd)
       !!!call f_free(eval)

       !!call timing(mpiworld(),'LAST','PR')
       call f_timing_checkpoint(ctr_name='LAST',mpi_comm=mpiworld(),nproc=mpisize(),&
                    gather_routine=gather_timings)

   end if



   if (kernel_purity) then
       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
            iproc, nproc, mpiworld(), smat(1), ovrlp_mat, &
            init_matmul=.false.)
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(kernel_file), &
            iproc, nproc, mpiworld(), smat(2), kernel_mat, &
            init_matmul=.true., filename_mult=trim(kernel_matmul_file), &
            matmul_matrix=MATMUL_REPLICATE_MATRIX)
       !!if (iproc==0) then
       !!    write(*,*) 'FIRST: sum(abs(ovrlp_mat%matrix_compr))',sum(abs(ovrlp_mat%matrix_compr))
       !!    write(*,*) 'FIRST: sum(abs(ovrlp_large%matrix_compr))',sum(abs(ovrlp_large%matrix_compr))
       !!end if
       call init_matrix_taskgroups_wrapper(iproc, nproc, mpiworld(), .false., 2, smat)
       call resize_matrix_to_taskgroup(smat(1), ovrlp_mat)
       call resize_matrix_to_taskgroup(smat(2), kernel_mat)

       if (iproc==0) then
           call yaml_mapping_open('Matrix properties')
           call write_sparsematrix_info(smat(1), 'Overlap matrix')
           call write_sparsematrix_info(smat(2), 'Density kernel')
           call yaml_mapping_close()
       end if

       iunit=99
       call f_open_file(iunit, file=fragment_file, binary=.false.)
       read(iunit,*) nat_frag
       fragment_atomnames = f_malloc_str(len(fragment_atomnames),nat_frag,id='nat_frag')
       do iat=1,nat_frag
           read(iunit,*) fragment_atomnames(iat)
       end do
       call f_close(iunit)

       ! Calculate K'=S^1/2*K*S^1/2
       kernel_ortho = matrices_null()
       kernel_ortho%matrix_compr = sparsematrix_malloc_ptr(smat(2), iaction=SPARSE_TASKGROUP, id='kernel_ortho%matrix_compr')
       call matrix_for_orthonormal_basis(iproc, nproc, mpiworld(), &
            1020, smat(1), smat(2), &
            ovrlp_mat, kernel_mat, 'plus', kernel_ortho%matrix_compr)


       ovrlp_large = matrices_null()
       ovrlp_large%matrix_compr = sparsematrix_malloc_ptr(smat(2), iaction=SPARSE_FULL, id='ovrlp_large%matrix_compr')
       KS_large = matrices_null()
       KS_large%matrix_compr = sparsematrix_malloc_ptr(smat(2), iaction=SPARSE_FULL, id='KS_large%matrix_compr')
       fragment_atom_id = f_malloc(nat_frag,id='fragment_atom_id')
       fragment_supfun_id = f_malloc(smat(2)%nfvctr,id='fragment_supfun_id')
       jat_start = 1
       found_a_fragment = .false.
       ifrag = 0
       if (iproc==0) then
           call yaml_map('Fragment composition',fragment_atomnames)
           call yaml_comment('Starting fragment purity analysis',hfill='~')
           call yaml_sequence_open('Purity analysis of fragment')
       end if
       search_fragments: do
           ifrag = ifrag + 1
           fragment_loop: do iat=1,nat_frag
               found = .false.
               search_loop: do jat=jat_start,smmd%nat
                   jtype = smmd%iatype(jat)
                   if (trim(fragment_atomnames(iat))==trim(smmd%atomnames(jtype))) then
                       ! Found an atom beloning to the fragment
                       fragment_atom_id(iat) = jat
                       jat_start = jat + 1
                       found = .true.
                       exit search_loop
                   end if
               end do search_loop
               if (.not.found) then
                   jat_start = smmd%nat + 1
                   exit fragment_loop
               end if
           end do fragment_loop
           if (found) then
               found_a_fragment = .true.
           else
               exit search_fragments
           end if
           
           ! Count how many support functions belong to the fragment
           nfvctr_frag = 0
           do i=1,smat(2)%nfvctr
               iiat = smmd%on_which_atom(i)
               do iat=1,nat_frag
                   if (iiat==fragment_atom_id(iat)) then
                       nfvctr_frag = nfvctr_frag + 1
                       exit
                   end if
               end do
               fragment_supfun_id(i) = nfvctr_frag
           end do



           kernel_fragment = f_malloc0((/nfvctr_frag,nfvctr_frag/),id='kernel_fragment')
           overlap_fragment = f_malloc0((/nfvctr_frag,nfvctr_frag/),id='overlap_fragment')
           ksk_fragment = f_malloc((/nfvctr_frag,nfvctr_frag/),id='ksk_fragment')
           tmpmat = f_malloc((/nfvctr_frag,nfvctr_frag/),id='tmpmat')

           if (smat(2)%nspin==1) then
               factor = 0.5_mp
           else if (smat(2)%nspin==2) then
               factor = 1.0_mp
           end if

           call transform_sparse_matrix(iproc, smat(1), smat(2), SPARSE_FULL, 'small_to_large', &
                                        smat_in=ovrlp_mat%matrix_compr, lmat_out=ovrlp_large%matrix_compr)
           !!if (iproc==0) then
           !!    write(*,*) 'sum(abs(ovrlp_mat%matrix_compr))',sum(abs(ovrlp_mat%matrix_compr))
           !!    write(*,*) 'sum(abs(ovrlp_large%matrix_compr))',sum(abs(ovrlp_large%matrix_compr))
           !!end if

           if (iproc==0) then
               call yaml_comment('Fragment number'//trim(yaml_toa(ifrag)),hfill='-')
               call yaml_sequence(advance='no')
               call yaml_map('Atom IDs',fragment_atom_id)
               !call yaml_map('Sup Fun IDs',fragment_supfun_id)
               call yaml_map('Submatrix size',nfvctr_frag)
           end if

           ! Calculate KSK-K for the submatrices
           !!if (iproc==0) write(*,*) 'extracting ovrlp_large'
           call extract_fragment_submatrix(smmd, smat(2), ovrlp_large, nat_frag, nfvctr_frag, &
                fragment_atom_id, fragment_supfun_id, overlap_fragment)
           !!if (iproc==0) write(*,*) 'extracting kernel_mat'
           call extract_fragment_submatrix(smmd, smat(2), kernel_mat, nat_frag, nfvctr_frag, &
                fragment_atom_id, fragment_supfun_id, kernel_fragment)
           !!if (iproc==0) then
           !!    write(*,*) 'overlap_fragment',overlap_fragment
           !!    write(*,*) 'kernel_fragment',kernel_fragment
           !!end if
           call gemm('n', 'n', nfvctr_frag, nfvctr_frag, nfvctr_frag, factor, kernel_fragment(1,1), nfvctr_frag, &
                overlap_fragment(1,1), nfvctr_frag, 0.d0, tmpmat(1,1), nfvctr_frag)
           !!if (iproc==0) then
           !!    write(*,*) 'tmpmat',tmpmat
           !!end if
           call gemm('n', 'n', nfvctr_frag, nfvctr_frag, nfvctr_frag, 1.d0, tmpmat(1,1), nfvctr_frag, &
                kernel_fragment(1,1), nfvctr_frag, 0.d0, ksk_fragment(1,1), nfvctr_frag)
           !!if (iproc==0) then
           !!    write(*,*) 'ksk_fragment',ksk_fragment
           !!end if
           maxdiff = 0.0_mp
           meandiff = 0.0_mp
           tracediff = 0.0_mp
           do i=1,nfvctr_frag
               do j=1,nfvctr_frag
                   tt = kernel_fragment(j,i)-ksk_fragment(j,i)
                   maxdiff = max(maxdiff,abs(tt))
                   meandiff = meandiff + abs(tt)
                   if (i==j) then
                       tracediff = tracediff + tt
                   end if
               end do
           end do
           totdiff = meandiff
           meandiff = meandiff/real(nfvctr_frag,kind=mp)**2
           if (iproc==0) then
               call yaml_mapping_open('analyzing KSK-K')
               call yaml_map('maximal deviation',maxdiff,fmt='(es10.3)')
               call yaml_map('total deviation',totdiff,fmt='(es10.3)')
               call yaml_map('average deviation',meandiff,fmt='(es10.3)')
               call yaml_map('trace difference',tracediff,fmt='(es10.3)')
               call yaml_mapping_close()
           end if

           ! Calculate KS, extract the submatrices and calculate (KS)^2-(KS)
           call matrix_matrix_multiplication(iproc, nproc, smat(2), kernel_mat, ovrlp_large, KS_large)
           call extract_fragment_submatrix(smmd, smat(2), KS_large, nat_frag, nfvctr_frag, &
                fragment_atom_id, fragment_supfun_id, kernel_fragment)
           call gemm('n', 'n', nfvctr_frag, nfvctr_frag, nfvctr_frag, factor, kernel_fragment(1,1), nfvctr_frag, &
                kernel_fragment(1,1), nfvctr_frag, 0.d0, tmpmat(1,1), nfvctr_frag)

           maxdiff = 0.0_mp
           meandiff = 0.0_mp
           tracediff = 0.0_mp
           do i=1,nfvctr_frag
               do j=1,nfvctr_frag
                   tt = kernel_fragment(j,i)-tmpmat(j,i)
                   maxdiff = max(maxdiff,abs(tt))
                   meandiff = meandiff + abs(tt)
                   if (i==j) then
                       tracediff = tracediff + tt
                   end if
               end do
           end do
           totdiff = meandiff
           meandiff = meandiff/real(nfvctr_frag,kind=mp)**2

           if (iproc==0) then
               call yaml_mapping_open('analyzing (KS)^2-(KS)')
               call yaml_map('maximal deviation',maxdiff,fmt='(es10.3)')
               call yaml_map('total deviation',totdiff,fmt='(es10.3)')
               call yaml_map('average deviation',meandiff,fmt='(es10.3)')
               call yaml_map('trace difference',tracediff,fmt='(es10.3)')
               call yaml_mapping_close()
           end if

           call extract_fragment_submatrix(smmd, smat(2), kernel_ortho, nat_frag, nfvctr_frag, &
                fragment_atom_id, fragment_supfun_id, kernel_fragment)
           call gemm('n', 'n', nfvctr_frag, nfvctr_frag, nfvctr_frag, factor, kernel_fragment(1,1), nfvctr_frag, &
                kernel_fragment(1,1), nfvctr_frag, 0.d0, tmpmat(1,1), nfvctr_frag)
           maxdiff = 0.0_mp
           meandiff = 0.0_mp
           tracediff = 0.0_mp
           do i=1,nfvctr_frag
               do j=1,nfvctr_frag
                   tt = kernel_fragment(j,i)-tmpmat(j,i)
                   maxdiff = max(maxdiff,abs(tt))
                   meandiff = meandiff + abs(tt)
                   if (i==j) then
                       tracediff = tracediff + tt
                   end if
               end do
           end do
           totdiff = meandiff
           meandiff = meandiff/real(nfvctr_frag,kind=mp)**2
           if (iproc==0) then
               call yaml_mapping_open('analyzing (S^1/2KS^1/2)^2 - S^1/2KS^1/2')
               call yaml_map('maximal deviation',maxdiff,fmt='(es10.3)')
               call yaml_map('total deviation',totdiff,fmt='(es10.3)')
               call yaml_map('average deviation',meandiff,fmt='(es10.3)')
               call yaml_map('trace difference',tracediff,fmt='(es10.3)')
               call yaml_mapping_close()
           end if

           call f_free(kernel_fragment)
           call f_free(overlap_fragment)
           call f_free(ksk_fragment)
           call f_free(tmpmat)

       end do search_fragments
       if (iproc==0) then
           call yaml_sequence_close()
           call yaml_comment('Fragment purity analysis completed',hfill='~')
       end if
       if (.not.found_a_fragment) then
           call f_err_throw('The specified fragment is not present')
       end if

       call deallocate_sparse_matrix_metadata(smmd)
       call deallocate_sparse_matrix(smat(1))
       call deallocate_sparse_matrix(smat(2))
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(kernel_mat)
       call deallocate_matrices(ovrlp_large)
       call deallocate_matrices(KS_large)
       call deallocate_matrices(kernel_ortho)
       call f_free_str(len(fragment_atomnames),fragment_atomnames)
       call f_free(fragment_atom_id)
       call f_free(fragment_supfun_id)

       call f_timing_checkpoint(ctr_name='LAST',mpi_comm=mpiworld(),nproc=mpisize(),&
                    gather_routine=gather_timings)


   end if


   if (manipulate_eigenvalue_spectrum) then

       if (iproc==0) then
           call yaml_mapping_open('Input parameters')
           call yaml_map('Manipulation mode',trim(manipulation_mode))
           call yaml_map('Matrix format',trim(matrix_format))
           call yaml_map('Matrix metadata file',trim(metadata_file))
           call yaml_map('Input hamiltonian file',trim(hamiltonian_file))
           call yaml_map('Input overlap file',trim(overlap_file))
           call yaml_map('Output manipulated hamiltonian file',trim(hamiltonian_manipulated_file))
           call yaml_map('HOMO state',ihomo_state)
           call yaml_map('Target HOMO eigenvalue',homo_value)
           call yaml_map('Target LUMO eigenvalue',lumo_value)
           call yaml_map('Target lowest eigenvalue',smallest_value)
           call yaml_map('Target highest eigenvalue',largest_value)
           call yaml_mapping_close()
       end if

       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
            iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
            init_matmul=.false.)
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(hamiltonian_file), &
            iproc, nproc, mpiworld(), smat_m, hamiltonian_mat, &
            init_matmul=.false.)
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(kernel_file), &
            iproc, nproc, mpiworld(), smat_l, ovrlp_minus_one_half(1), &
            init_matmul=.true., filename_mult=trim(kernel_matmul_file), &
            matmul_matrix=MATMUL_REPLICATE_MATRIX)

       if (iproc==0) then
           call yaml_mapping_open('Matrix properties')
           call write_sparsematrix_info(smat_s, 'Overlap matrix')
           call write_sparsematrix_info(smat_m, 'Hamiltonian')
           call write_sparsematrix_info(smat_l, 'Kernel')
           call yaml_mapping_close()
       end if

       ovrlp_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='ovrlp_mat%matrix')
       call uncompress_matrix(iproc, nproc, &
            smat_s, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_mat%matrix)
       hamiltonian_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='hamiltonian_mat%matrix')
       call uncompress_matrix(iproc, nproc, &
            smat_m, inmat=hamiltonian_mat%matrix_compr, outmat=hamiltonian_mat%matrix)
       eval = f_malloc(smat_s%nfvctr,id='eval')

       !!! Diagonalize the original Hamiltonian matrix
       !!if (iproc==0) then
       !!    call yaml_comment('Diagonalizing the matrix',hfill='~')
       !!end if
       hamiltonian_tmp = f_malloc((/smat_s%nfvctr,smat_s%nfvctr/),id='hamiltonian_tmp')
       ovrlp_tmp = f_malloc((/smat_s%nfvctr,smat_s%nfvctr/),id='ovrlp_tmp')
       !!call f_memcpy(src=hamiltonian_mat%matrix,dest=hamiltonian_tmp)
       !!call f_memcpy(src=ovrlp_mat%matrix,dest=ovrlp_tmp)
       !!call diagonalizeHamiltonian2(iproc, nproc, mpiworld(), scalapack_blocksize, &
       !!     smat_s%nfvctr, hamiltonian_tmp, ovrlp_tmp, eval)
       !!if (iproc==0) then
       !!    call yaml_comment('Matrix succesfully diagonalized',hfill='~')
       !!end if

       !!if (iproc==0) then
       !!    call yaml_mapping_open('Eigenvalue spectrum')
       !!    call yaml_map('Smallest value',eval(1))
       !!    call yaml_map('Largest value',eval(smat_s%nfvctr))
       !!    call yaml_map('HOMO value',eval(ihomo_state))
       !!    call yaml_map('LUMO value',eval(ihomo_state+1))
       !!    call yaml_map('HOMO-LUMO gap',eval(ihomo_state+1)-eval(ihomo_state))
       !!    call yaml_mapping_close()
       !!end if

       ! Manipulate the Hamiltonian matrix such that it has the desired spectral properteis
       if (iproc==0) then
           call yaml_mapping_open('Manipulate the Hamiltonian spectrum')
       end if

       mode_if: if (manipulation_mode=='diagonal') then
           ! Set the matrix to zero
           !call f_zero(hamiltonian_mat%matrix)
           call f_zero(ovrlp_mat%matrix)
           call f_zero(hamiltonian_tmp)

           ! Set the lowest eigenvalue
           hamiltonian_tmp(1,1) = smallest_value

           ! Set the highest eigenvalue
           hamiltonian_tmp(smat_s%nfvctr,smat_s%nfvctr) = largest_value

           ! Set the HOMO eigenvalue
           hamiltonian_tmp(ihomo_state,ihomo_state) = homo_value

           ! Set the LUMO eigenvalue
           hamiltonian_tmp(ihomo_state+1,ihomo_state+1) = lumo_value

           ! Set the remaining values at random
           call f_random_number(tt, reset=.true.)
           mult_factor = homo_value-smallest_value
           add_shift = smallest_value
           do i=2,ihomo_state-1
               call f_random_number(tt)
               tt = tt*mult_factor+add_shift
               hamiltonian_tmp(i,i) = tt
           end do
           mult_factor = largest_value-lumo_value
           add_shift = lumo_value
           do i=ihomo_state+2,smat_s%nfvctr-1
               call f_random_number(tt)
               tt = tt*mult_factor+add_shift
               hamiltonian_tmp(i,i) = tt
           end do
           ! Calculate S^1/2
           !call matrices_init(smat_l, ovrlp_minus_one_half(1))
           call matrix_chebyshev_expansion(iproc, nproc, mpiworld(), 1, (/0.5_mp/), &
                smat_s, smat_l, ovrlp_mat, ovrlp_minus_one_half)
           call uncompress_matrix(iproc, nproc, &
                smat_l, inmat=ovrlp_minus_one_half(1)%matrix_compr, outmat=ovrlp_tmp)

           call dgemm_parallel(iproc, nproc, scalapack_blocksize, mpi_comm_world, &
                'n', 'n', smat_m%nfvctr, smat_m%nfvctr, smat_m%nfvctr, &
                1.0_mp, ovrlp_tmp(1:,1:), smat_m%nfvctr, &
                hamiltonian_tmp(1:,1:), smat_m%nfvctr, 0.0_mp, ovrlp_mat%matrix(1:,1:,1), smat_m%nfvctr)
           call f_free(hamiltonian_tmp)
           call dgemm_parallel(iproc, nproc, scalapack_blocksize, mpi_comm_world, &
                'n', 'n', smat_m%nfvctr, smat_m%nfvctr, smat_m%nfvctr, &
                1.0_mp, ovrlp_mat%matrix(1:,1:,1), smat_m%nfvctr, &
                ovrlp_tmp(1:,1:), smat_m%nfvctr, 0.0_mp, hamiltonian_mat%matrix(1:,1:,1), smat_m%nfvctr)
           do i=1,smat_s%nfvctr
               ovrlp_mat%matrix(i,i,1) = 1.0_mp
           end do
           !call yaml_map('H',ovrlp_tmp(:,:))
       else if (manipulation_mode=='full') then mode_if

           ! Diagonalize the original Hamiltonian matrix
           if (iproc==0) then
               call yaml_comment('Diagonalizing the matrix',hfill='~')
           end if
           !!hamiltonian_tmp = f_malloc((/smat_s%nfvctr,smat_s%nfvctr/),id='hamiltonian_tmp')
           !!ovrlp_tmp = f_malloc((/smat_s%nfvctr,smat_s%nfvctr/),id='ovrlp_tmp')
           call f_memcpy(src=hamiltonian_mat%matrix,dest=hamiltonian_tmp)
           call f_memcpy(src=ovrlp_mat%matrix,dest=ovrlp_tmp)
           call diagonalizeHamiltonian2(iproc, nproc, mpiworld(), 1, scalapack_blocksize, &
                smat_s%nfvctr, hamiltonian_tmp, ovrlp_tmp, eval)
           if (iproc==0) then
               call yaml_comment('Matrix succesfully diagonalized',hfill='~')
           end if

           if (iproc==0) then
               call yaml_mapping_open('Eigenvalue spectrum')
               call yaml_map('Smallest value',eval(1))
               call yaml_map('Largest value',eval(smat_s%nfvctr))
               call yaml_map('HOMO value',eval(ihomo_state))
               call yaml_map('LUMO value',eval(ihomo_state+1))
               call yaml_map('HOMO-LUMO gap',eval(ihomo_state+1)-eval(ihomo_state))
               call yaml_mapping_close()
           end if


           ! Scale the gap to the desired value
           gap = eval(ihomo_state+1)-eval(ihomo_state)
           gap_target = lumo_value-homo_value
           scale_value = gap_target/gap
           if (iproc==0) then
               call yaml_map('Scaling factor',scale_value)
           end if
           call vscal(smat_s%nfvctr**2, scale_value, hamiltonian_mat%matrix(1,1,1), 1)

           ! Move the HOMO level to the desired value
           call f_zero(ovrlp_tmp)
           shift_value=homo_value-scale_value*eval(ihomo_state)
           !if (iproc==0) then
           !    call yaml_map('shift_value',shift_value)
           !end if
           if (iproc==0) then
               call yaml_map('Shift value',shift_value)
           end if
           do i=1,smat_s%nfvctr
               do j=1,smat_s%nfvctr
                   !ovrlp_tmp(i,i)=shift_value*ovrlp_mat%matrix(i,i,1)
                   hamiltonian_mat%matrix(j,i,1) = hamiltonian_mat%matrix(j,i,1) + shift_value*ovrlp_mat%matrix(j,i,1)
               end do
           end do
           call axpy(smat_s%nfvctr**2, 1.0_mp, ovrlp_tmp(1,1), 1, hamiltonian_mat%matrix(1,1,1), 1)

           matrix_tmp = f_malloc((/smat_s%nfvctr,smat_s%nfvctr/),id='matrix_tmp')
           !call gemm('n', 't', smat_s%nfvctr, smat_s%nfvctr, 1, &
           !     1.0_mp, hamiltonian_tmp(1,1), smat_s%nfvctr, &
           !     hamiltonian_tmp(1,1), smat_s%nfvctr, 0.0_mp, ovrlp_tmp(1,1), smat_s%nfvctr)
           !call gemm('n', 'n', smat_s%nfvctr, smat_s%nfvctr, smat_s%nfvctr, &
           !     1.0_mp, ovrlp_mat%matrix(1,1,1), smat_s%nfvctr, &
           !     ovrlp_tmp(1,1), smat_s%nfvctr, 0.0_mp, matrix_tmp(1,1), smat_s%nfvctr)
           !call gemm('n', 'n', smat_s%nfvctr, smat_s%nfvctr, smat_s%nfvctr, &
           !     1.0_mp, matrix_tmp(1,1), smat_s%nfvctr, &
           !     ovrlp_mat%matrix(1,1,1), smat_s%nfvctr, 0.0_mp, ovrlp_tmp(1,1), smat_s%nfvctr)
           !call axpy(smat_s%nfvctr**2, smallest_value-(scale_value*eval(1)+shift_value), &
           !     ovrlp_tmp(1,1), 1, hamiltonian_mat%matrix(1,1,1), 1)
           ! Move the lowest eigenvalue to the desired value.
           ! This is also necessary for all eigevalues that are smaller than the new target value.
           if (iproc==0) then
               call yaml_sequence_open('Moving lower eigenvalues')
           end if
           do ieval=1,smat_s%nfvctr
               actual_eval = scale_value*eval(ieval)+shift_value
               if (ieval==1 .or. actual_eval<smallest_value) then
                   if (iproc==0) then
                       call yaml_sequence(advance='no')
                       call yaml_mapping_open(flow=.true.)
                       call yaml_map('eigenvalue number',ieval)
                       call yaml_map('Shift value',smallest_value-actual_eval)
                       call yaml_mapping_close()
                   end if
                   !call gemm('n', 't', smat_s%nfvctr, smat_s%nfvctr, 1, &
                   !     1.0_mp, hamiltonian_tmp(1,ieval), smat_s%nfvctr, &
                   !     hamiltonian_tmp(1,ieval), smat_s%nfvctr, 0.0_mp, ovrlp_tmp(1,1), smat_s%nfvctr)
                   !call gemm('n', 'n', smat_s%nfvctr, smat_s%nfvctr, smat_s%nfvctr, &
                   !     1.0_mp, ovrlp_mat%matrix(1,1,1), smat_s%nfvctr, &
                   !     ovrlp_tmp(1,1), smat_s%nfvctr, 0.0_mp, matrix_tmp(1,1), smat_s%nfvctr)
                   !call dgemm('n', 'n', smat_s%nfvctr, smat_s%nfvctr, smat_s%nfvctr, &
                   !     1.0_mp, matrix_tmp(1,1), smat_s%nfvctr, &
                   !     ovrlp_mat%matrix(1,1,1), smat_s%nfvctr, 0.0_mp, ovrlp_tmp(1,1), smat_s%nfvctr)
                   call dgemm_parallel(iproc, nproc, scalapack_blocksize, mpi_comm_world, &
                        'n', 't', smat_s%nfvctr, smat_s%nfvctr, 1, &
                        1.0_mp, hamiltonian_tmp(1:,ieval:ieval), smat_s%nfvctr, &
                        hamiltonian_tmp(1:,ieval:ieval), smat_s%nfvctr, 0.0_mp, ovrlp_tmp, smat_s%nfvctr)
                   call dgemm_parallel(iproc, nproc, scalapack_blocksize, mpi_comm_world, &
                        'n', 'n', smat_s%nfvctr, smat_s%nfvctr, smat_s%nfvctr, &
                        1.0_mp, ovrlp_mat%matrix(1:,1:,1), smat_s%nfvctr, &
                        ovrlp_tmp, smat_s%nfvctr, 0.0_mp, matrix_tmp, smat_s%nfvctr)
                   call dgemm_parallel(iproc, nproc, scalapack_blocksize, mpi_comm_world, &
                        'n', 'n', smat_s%nfvctr, smat_s%nfvctr, smat_s%nfvctr, &
                        1.0_mp, matrix_tmp, smat_s%nfvctr, &
                        ovrlp_mat%matrix(1:,1:,1), smat_s%nfvctr, 0.0_mp, ovrlp_tmp, smat_s%nfvctr)
                   call axpy(smat_s%nfvctr**2, smallest_value-actual_eval, &
                        ovrlp_tmp(1,1), 1, hamiltonian_mat%matrix(1,1,1), 1)
               end if
           end do
           if (iproc==0) then
               call yaml_sequence_close()
           end if

           ! Move the higest eigenvalue to the desired value.
           ! This is also necessary for all eigevalues that are bigger than the new target value.
           if (iproc==0) then
               call yaml_sequence_open('Moving upper eigenvalues')
           end if
           do ieval=1,smat_s%nfvctr
               actual_eval = scale_value*eval(ieval)+shift_value
               if (ieval==smat_s%nfvctr .or. actual_eval>largest_value) then
                   if (iproc==0) then
                       call yaml_sequence(advance='no')
                       call yaml_mapping_open(flow=.true.)
                       call yaml_map('eigenvalue number',ieval)
                       call yaml_map('Shift value',smallest_value-actual_eval)
                       call yaml_mapping_close()
                   end if
                   !call gemm('n', 't', smat_s%nfvctr, smat_s%nfvctr, 1, &
                   !     1.0_mp, hamiltonian_tmp(1,ieval), smat_s%nfvctr, &
                   !     hamiltonian_tmp(1,ieval), smat_s%nfvctr, 0.0_mp, ovrlp_tmp(1,1), smat_s%nfvctr)
                   !call gemm('n', 'n', smat_s%nfvctr, smat_s%nfvctr, smat_s%nfvctr, &
                   !     1.0_mp, ovrlp_mat%matrix(1,1,1), smat_s%nfvctr, &
                   !     ovrlp_tmp(1,1), smat_s%nfvctr, 0.0_mp, matrix_tmp(1,1), smat_s%nfvctr)
                   !call gemm('n', 'n', smat_s%nfvctr, smat_s%nfvctr, smat_s%nfvctr, &
                   !     1.0_mp, matrix_tmp(1,1), smat_s%nfvctr, &
                   !     ovrlp_mat%matrix(1,1,1), smat_s%nfvctr, 0.0_mp, ovrlp_tmp(1,1), smat_s%nfvctr)
                   call dgemm_parallel(iproc, nproc, scalapack_blocksize, mpi_comm_world, &
                        'n', 't', smat_s%nfvctr, smat_s%nfvctr, 1, &
                        1.0_mp, hamiltonian_tmp(1:,ieval:ieval), smat_s%nfvctr, &
                        hamiltonian_tmp(1:,ieval:ieval), smat_s%nfvctr, 0.0_mp, ovrlp_tmp, smat_s%nfvctr)
                   call f_free(hamiltonian_tmp)
                   call dgemm_parallel(iproc, nproc, scalapack_blocksize, mpi_comm_world, &
                        'n', 'n', smat_s%nfvctr, smat_s%nfvctr, smat_s%nfvctr, &
                        1.0_mp, ovrlp_mat%matrix(1:,1:,1), smat_s%nfvctr, &
                        ovrlp_tmp, smat_s%nfvctr, 0.0_mp, matrix_tmp, smat_s%nfvctr)
                   call dgemm_parallel(iproc, nproc, scalapack_blocksize, mpi_comm_world, &
                        'n', 'n', smat_s%nfvctr, smat_s%nfvctr, smat_s%nfvctr, &
                        1.0_mp, matrix_tmp, smat_s%nfvctr, &
                        ovrlp_mat%matrix(1:,1:,1), smat_s%nfvctr, 0.0_mp, ovrlp_tmp, smat_s%nfvctr)
                   call axpy(smat_s%nfvctr**2, largest_value-actual_eval, &
                        ovrlp_tmp(1,1), 1, hamiltonian_mat%matrix(1,1,1), 1)
               end if
           end do
           if (iproc==0) then
               call yaml_sequence_close()
               call yaml_mapping_close()
           end if
           !call axpy(smat_s%nfvctr**2, smallest_value-(scale_value*eval(smat_s%nfvctr)+shift_value), &
           !     ovrlp_tmp(1,1), 1, hamiltonian_mat%matrix(1,1,1), 1)

           call f_free(matrix_tmp)

       else mode_if
           call f_err_throw("wrong manipulation mode; must be 'diagonal' or 'full'")
       end if mode_if

       ! Compress the matrix
       call compress_matrix(iproc, nproc, smat_m, hamiltonian_mat%matrix, hamiltonian_mat%matrix_compr)

       ! Diagonalize the modified Hamiltonian matrix
       if (iproc==0) then
           call yaml_comment('Diagonalizing the matrix',hfill='~')
       end if
       !!hamiltonian_tmp = f_malloc((/smat_s%nfvctr,smat_s%nfvctr/),id='hamiltonian_tmp')
       !!ovrlp_tmp = f_malloc((/smat_s%nfvctr,smat_s%nfvctr/),id='ovrlp_tmp')
       !!call f_memcpy(src=hamiltonian_mat%matrix,dest=hamiltonian_tmp)
       !!call f_memcpy(src=ovrlp_mat%matrix,dest=ovrlp_tmp)
       !!call diagonalizeHamiltonian2(iproc, nproc, mpiworld(), scalapack_blocksize, &
       !!     smat_s%nfvctr, hamiltonian_tmp, ovrlp_tmp, eval)
       eval_min = f_malloc(smat_m%nspin,id='eval_min')
       eval_max = f_malloc(smat_m%nspin,id='eval_max')
       call get_minmax_eigenvalues(iproc, nproc, mpiworld(), 'generalized', scalapack_blocksize, &
            smat_m, hamiltonian_mat, eval_min, eval_max, &
            diag_algorithm, quiet=.true., smat2=smat_s, mat2=ovrlp_mat, evals=eval)
       !!do i=1,smat_m%nfvctr
       !!    write(*,*) 'i',eval(i)
       !!end do
       call f_free(eval_min)
       call f_free(eval_max)

       if (iproc==0) then
           call yaml_comment('Matrix succesfully diagonalized',hfill='~')
       end if

       if (iproc==0) then
           call yaml_mapping_open('Eigenvalue spectrum')
           !call yaml_map('EVALS',eval)
           call yaml_map('Smallest value',eval(1))
           call yaml_map('Largest value',eval(smat_s%nfvctr))
           call yaml_map('HOMO value',eval(ihomo_state))
           call yaml_map('LUMO value',eval(ihomo_state+1))
           call yaml_map('HOMO-LUMO gap',eval(ihomo_state+1)-eval(ihomo_state))
           call yaml_mapping_close()
       end if

       ! Write the manipulated Hamiltonian matrix
       index_dot = index(hamiltonian_manipulated_file,'.',back=.true.)
       outfile_base = hamiltonian_manipulated_file(1:index_dot-1)
       outfile_extension = hamiltonian_manipulated_file(index_dot:)
       if (trim(matrix_format)=='serial_text') then
           if (trim(outfile_extension)/='.txt') then
               call f_err_throw('Wrong file extension; must be .txt, but found '//trim(outfile_extension))
           end if
       else if (trim(matrix_format)=='parallel_mpi-native') then
           if (trim(outfile_extension)/='.mpi') then
               call f_err_throw('Wrong file extension; must be .mpi, but found '//trim(outfile_extension))
           end if
       else
           call f_err_throw('Wrong matrix format')
       end if
       !call compress_matrix(iproc, nproc, smat_m, hamiltonian_mat%matrix, hamiltonian_mat%matrix_compr)
       call write_sparse_matrix(matrix_format, iproc, nproc, mpiworld(), &
            smat_m, hamiltonian_mat, trim(outfile_base))

       !!!! Write the manipulated overlap matrix
       !!!index_dot = index(overlap_manipulated_file,'.',back=.true.)
       !!!outfile_base = overlap_manipulated_file(1:index_dot-1)
       !!!outfile_extension = overlap_manipulated_file(index_dot:)
       !!!if (trim(matrix_format)=='serial_text') then
       !!!    if (trim(outfile_extension)/='.txt') then
       !!!        call f_err_throw('Wrong file extension; must be .txt, but found '//trim(outfile_extension))
       !!!    end if
       !!!else if (trim(matrix_format)=='parallel_mpi-native') then
       !!!    if (trim(outfile_extension)/='.mpi') then
       !!!        call f_err_throw('Wrong file extension; must be .mpi, but found '//trim(outfile_extension))
       !!!    end if
       !!!else
       !!!    call f_err_throw('Wrong matrix format')
       !!!end if
       !!!call compress_matrix(iproc, nproc, smat_s, ovrlp_mat%matrix, ovrlp_mat%matrix_compr)
       !!!call write_sparse_matrix(matrix_format, iproc, nproc, mpiworld(), &
       !!!     smat_s, ovrlp_mat, trim(outfile_base))

       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_m)
       call deallocate_sparse_matrix(smat_l)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(hamiltonian_mat)
       call deallocate_matrices(ovrlp_minus_one_half(1))
       call deallocate_sparse_matrix_metadata(smmd)
       call f_free(eval)
       call f_free(ovrlp_tmp)

       !!call timing(mpiworld(),'LAST','PR')
       call f_timing_checkpoint(ctr_name='LAST',mpi_comm=mpiworld(),nproc=mpisize(),&
                    gather_routine=gather_timings)


   end if



   if (calculate_entropy) then
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
            iproc, nproc, mpiworld(), smat(1), ovrlp_mat, &
            init_matmul=.false.)
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(kernel_file), &
            iproc, nproc, mpiworld(), smat(2), kernel_mat, &
            init_matmul=.true., filename_mult=trim(kernel_matmul_file), &
            matmul_matrix=MATMUL_REPLICATE_MATRIX)
       call init_matrix_taskgroups_wrapper(iproc, nproc, mpiworld(), .false., 2, smat)
       call resize_matrix_to_taskgroup(smat(1), ovrlp_mat)
       call resize_matrix_to_taskgroup(smat(2), kernel_mat)

       if (iproc==0) then
           call yaml_mapping_open('Matrix properties')
           call write_sparsematrix_info(smat(1), 'Overlap matrix')
           call write_sparsematrix_info(smat(2), 'Density kernel')
           call yaml_mapping_close()
       end if

       ovrlp_minus_one_half(1) = matrices_null()
       ovrlp_minus_one_half(1)%matrix_compr = &
           sparsematrix_malloc_ptr(smat(2),iaction=SPARSE_TASKGROUP,id='ovrlp_minus_one_half(1)%matrix_compr')

       if (iproc==0) then
           call yaml_mapping_open('Calculate S^-1/2')
       end if
       call overlap_plusminus_onehalf('minus', 'ICE', iproc, nproc, mpiworld(), &
            smat(1), smat(2), ovrlp_mat, ovrlp_minus_one_half(1))
       if (iproc==0) then
           call yaml_mapping_close()
       end if

       if (iproc==0) then
           call yaml_mapping_open('Calculate entropy')
       end if
       call calculate_entropy_term(iproc, nproc, mpiworld(), &
            kT, smat(1), smat(2), &
            ovrlp_mat, kernel_mat, ovrlp_minus_one_half(1), &
            accuracy_entropy, eTS)
       if (iproc==0) then
           call yaml_map('eTS',eTS)
           call yaml_mapping_close()
       end if

       call f_timing_checkpoint(ctr_name='LAST',mpi_comm=mpiworld(),nproc=mpisize(),&
                    gather_routine=gather_timings)

       call deallocate_sparse_matrix(smat(1))
       call deallocate_sparse_matrix(smat(2))
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(kernel_mat)
       call deallocate_matrices(ovrlp_minus_one_half(1))

   end if




   if (analyze_density_matrix_dense) then

       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       kernel_mat = matrices_null()
       call read_dense_matrix(.not.matrix_format=='serial_text', trim(kernel_file), &
            iproc, nproc, mpiworld(), nfvctr_l, nspin, kernel_mat%matrix)

       if (nspin/=1) then
           call f_err_throw('dense density kernel analysis at the moment only possible for nspin=1')
       end if

       select case (smmd%geocode)
       case ('F')
           periodic = [.false.,.false.,.false.]
       case ('P')
           periodic = [.true.,.true.,.true.]
       case default
           call f_err_throw('Unknown boundary conditions')
       end select

       nnl = int(int(nfvctr_l,f_long)*int(nfvctr_l,f_long),f_long)
       darr = f_malloc(nnl,id='darr')
       valarr = f_malloc(nnl,id='valarr')

       ! Calculate the distances between the atoms to which the support functions belong
       dmin = huge(0._mp)
       dmax = 0._mp
       iil = int(0,f_long)
       if (iproc==0) then
           call yaml_comment('Calculating values',hfill='-')
           bar=f_progress_bar_new(nstep=nfvctr_l)
       end if
       do i=1,nfvctr_l
           iat = smmd%on_which_atom(i)
           ri = smmd%rxyz(1:3,iat)
           do j=1,nfvctr_l
               jat = smmd%on_which_atom(j)
               rj = smmd%rxyz(1:3,jat)
               iil = iil + int(1,f_long)
               d = get_minimal_distance(ri, rj, smmd%cell_dim, periodic)
               val = kernel_mat%matrix(j,i,1)
               if (d < dmin) then
                   dmin = d
               end if
               if (d > dmax) then
                   dmax = d
               end if
               darr(iil) = d
               valarr(iil) = val
           end do
           if ((mod(i,2)==0 .or. i==nfvctr_l) .and. iproc==0) call dump_progress_bar(bar,step=i)
       end do
       if (iproc==0) then
           call yaml_comment('Done',hfill='-')
       end if
       if (iil/=nnl) then
           call f_err_throw('iil/=nnl')
       end if

       if (.not.only_binned_values .and. iproc==0) then
           call yaml_comment('Writing values',hfill='-')
           call yaml_map('Output file for raw values',trim(output_file))
           iunit = 99
           call f_open_file(iunit, file=trim(output_file), binary=.false.)
           write(iunit,'(a)') '# distance      value'
           do iil=1,nnl
               write(iunit,'(2es20.12)') darr(iil), valarr(iil)
           end do
           call f_close(iunit)
       end if

       ! Calculate the mean value within bins.
       ! Make the maximal value slightly larger to take into account values which have the maximal values
       if (iproc==0) then
           call yaml_map('Output file for binned values',trim(output_bins_file))
           iunit = 99
           call f_open_file(iunit, file=trim(output_bins_file)//'', binary=.false.)
           write(iunit,'(a)') '# distance      value'

           dmax = dmax + 1.e-10
           nbin = 100
           binwidth = (dmax-dmin)/real(nbin,f_double)
           do i=0,nbin-1
               dstart = dmin + real(i,f_double)*binwidth
               dend = dmin + real(i+1,f_double)*binwidth
               val_mean = 0._mp
               jjl = int(0,f_long)
               !write(*,*) 'i, dstart, dend, dmin, dmax', i, dstart, dend, dmin, dmax
               do iil=1,nnl
                   if (darr(iil)>=dstart .and. darr(iil)<=dend) then
                       val_mean = val_mean + valarr(iil)
                       jjl = jjl + int(1,f_long)
                   end if 
               end do
               if (jjl>0) then
                   write(iunit,'(2es20.12)') 0.5_mp*(dstart+dend), val_mean/real(jjl,f_double)
               end if
           end do
           call yaml_comment('Done',hfill='-')
           call f_close(iunit)
       end if

       call f_free(darr)
       call f_free(valarr)

       call deallocate_matrices(kernel_mat)
       call deallocate_sparse_matrix_metadata(smmd)

   end if




   if (find_fragments) then

       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
            iproc, nproc, mpiworld(), smat(1), ovrlp_mat, &
            init_matmul=.false.)
       call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(kernel_file), &
            iproc, nproc, mpiworld(), smat(2), kernel_mat, &
            init_matmul=.true., filename_mult=trim(kernel_matmul_file), &
            matmul_matrix=MATMUL_REPLICATE_MATRIX)
       call init_matrix_taskgroups_wrapper(iproc, nproc, mpiworld(), .false., 2, smat)
       call resize_matrix_to_taskgroup(smat(1), ovrlp_mat)
       call resize_matrix_to_taskgroup(smat(2), kernel_mat)

       if (iproc==0) then
           call yaml_mapping_open('Matrix properties')
           call write_sparsematrix_info(smat(1), 'Overlap matrix')
           call write_sparsematrix_info(smat(2), 'Density kernel')
           call yaml_mapping_close()
       end if

       ! Calculate K'=S^1/2*K*S^1/2
       kernel_ortho = matrices_null()
       kernel_ortho%matrix_compr = sparsematrix_malloc_ptr(smat(2), iaction=SPARSE_TASKGROUP, id='kernel_ortho%matrix_compr')
       call matrix_for_orthonormal_basis(iproc, nproc, mpiworld(), &
            1020, smat(1), smat(2), &
            ovrlp_mat, kernel_mat, 'plus', kernel_ortho%matrix_compr)
       kernel_ortho%matrix = sparsematrix_malloc_ptr(smat(2), iaction=DENSE_FULL, id='kernel_ortho%matrix')
       call uncompress_matrix2(iproc, nproc, mpiworld(), smat(2), kernel_ortho%matrix_compr, kernel_ortho%matrix)

       call vscal(smat(1)%nfvctr**2, 0.5d0, kernel_ortho%matrix(1,1,1), 1)

       tmpmat = f_malloc([smat(1)%nfvctr,smat(1)%nfvctr],id='tmpmat')
       tmpmat = matmul(kernel_ortho%matrix(:,:,1),kernel_ortho%matrix(:,:,1))
       !write(*,*) 'kernel_ortho%matrix(1,1,1), tmpmat(1,1)', kernel_ortho%matrix(1,1,1), tmpmat(1,1)

       ! Maybe we should check here the consistency of nfvctr etc...

       ! Solve the overdetermined system in a least-square sense

       maxfrag = 4 !12

       call set_N(smat(1)%nfvctr)
       call set_maxfrag(maxfrag)
       call allocate_K()
       call set_K(kernel_ortho%matrix(:,:,1))
       call allocate_xold()

       in_frag = f_malloc(smat(1)%nfvctr,id='in_frag')
       nat_in_frag = f_malloc(smmd%nat,id='nat_in_frag')

       ftol = 1.e-5_mp
       xtol = 1.e-5_mp
       gtol = 1.e-5_mp
       maxfev = 10000
       epsfcn = 1.e-10_mp
       diag = f_malloc(smat(1)%nfvctr, id='diag')
       ipvt = f_malloc(smat(1)%nfvctr, id='ipvt')
       qtf = f_malloc(smat(1)%nfvctr, id='qtf')
       fvec = f_malloc(smat(1)%nfvctr+maxfrag+1, id='fvec')
       fjac = f_malloc([smat(1)%nfvctr+maxfrag+1,smat(1)%nfvctr], id='fjac')
       wa1 = f_malloc(smat(1)%nfvctr, id='wa1')
       wa2 = f_malloc(smat(1)%nfvctr, id='wa2')
       wa3 = f_malloc(smat(1)%nfvctr, id='wa3')
       wa4 = f_malloc(smat(1)%nfvctr+maxfrag+1, id='wa4')
       weights = f_malloc(smat(1)%nfvctr,id='weights')
       natlist = f_malloc(smmd%nat,id='natlist')
       at_assigned = f_malloc(smmd%nat,id='at_assigned')
       nfragmax = 100
       fraglist = f_malloc([smmd%nat,nfragmax],id='fraglist')
       natfrag = f_malloc(nfragmax,id='natfrag')
       purity_indicator_list = f_malloc(nfragmax, id='purity_indicator_list')

       if (iproc==0) then
           call yaml_sequence_open('Searching fragments')
       end if

       nfragfound = 0

       if (frag_elements_start<1) then
           frag_elements_start = 1
           if (iproc==0) then
               call yaml_warning('frag_elements_start is to small, adjusting to: '//trim(yaml_toa(frag_elements_start)))
           end if
       end if
       if (frag_elements_end>smat(1)%nfvctr) then
           frag_elements_end = smat(1)%nfvctr
           if (iproc==0) then
               call yaml_warning('frag_elements_end is to large, adjusting to: '//trim(yaml_toa(frag_elements_end)))
           end if
       end if

       do iel=frag_elements_start,frag_elements_end
           call set_nel(iel)
           ifound = 0
           call set_ifound(ifound)
           call yaml_comment('Searching fragment with '//trim(yaml_toa(iel))//' elements')
           do i=1,maxfrag
               !ii = i + 1
               ii = ifound + 2
               !weights(:) = 1.d0
               weights(:) = 0.d0
               !weights(13:18) = 0.7
               call random_number(weights)
               call orthogonalize_inputguess(smat(1)%nfvctr+ii, smat(1)%nfvctr, weights)
               !weights = 0.5d0
               !weights(1:6) = 0.8d0
               !weights(7:12) = 0.d0
               !write(*,*) 'call with smat(1)%nfvctr+ii, ii, ifound',smat(1)%nfvctr+ii, ii, ifound
               call fcn_fragments(smat(1)%nfvctr+ii, smat(1)%nfvctr, weights, fvec, iflag)
               !write(*,*) 'weights', weights
               !write(*,*) 'fvec', fvec
               call lmdif_wrapper(fcn_fragments, smat(1)%nfvctr+ii, smat(1)%nfvctr, weights, fvec, ftol, xtol, gtol, &
                    maxfev, epsfcn, diag, 1, 100._mp, 1, info, nfev, &
                    fjac, smat(1)%nfvctr+maxfrag+1, ipvt, qtf, &
                    wa1, wa2, wa3, wa4)
               !write(*,*) 'info', info

               ! Determine whether the found solutions corresponds to roughly integer weights.
               !write(*,'(a,2i5,12(f7.3))') 'iel, i, weights', iel, i, weights
               weights_integer = .true.
               do l=1,smat(1)%nfvctr
                   if (abs(weights(l))<0.1d0) then
                       in_frag(l) = .false.
                   else if (abs(weights(l)-1.d0)<0.1d0) then
                       in_frag(l) = .true.
                   else
                       !write(*,*) 'non-integer:',l, weights(l)
                       weights_integer = .false.
                   end if
               end do
               if (.not.weights_integer .or. all(.not.in_frag)) cycle

               ! Determine to which atoms the support functions belong.
               ! Also check that all support functions of a given atom belong to the selected fragment.
               at_assigned(:) = .false.
               atom_complete = .true.
               do l=1,smat(1)%nfvctr
                   iiat = smmd%on_which_atom(l)
                   if (in_frag(l)) then
                       if (at_assigned(iiat) .and. .not.nat_in_frag(iiat)) then
                           ! Not all support functions of this atom belong to the fragment.
                           atom_complete = .false.
                       end if
                       nat_in_frag(iiat) = .true.
                   else
                       if (at_assigned(iiat) .and. nat_in_frag(iiat)) then
                           ! Not all support functions of this atom belong to the fragment.
                           atom_complete = .false.
                       end if
                       nat_in_frag(iiat) = .false.
                   end if
                   at_assigned(iiat) = .true.
               end do

               if (.not.atom_complete) then
                   !write(*,*) 'cycle due no incomplete atom'
                   cycle
               end if

               ifound = ifound + 1
               call set_ifound(ifound)

               ! Set the xold array
               ! Convert the weights into proper 0 and 1
               weights(:) = 0.d0
               do l=1,smat(1)%nfvctr
                   if (in_frag(l)) then
                       weights(l) = 1.d0
                   else
                       weights(l) = 0.d0
                   end if
               end do
               call set_xold(ifound, weights)

               ! Check whether the fragment is indeed pure.
               nat_frag = 0
               do iat=1,smmd%nat
                   if (nat_in_frag(iat)) then
                       nat_frag = nat_frag + 1
                       natlist(nat_frag) = iat
                   end if
               end do

               charge = 0.d0
               do l=1,smmd%nat
                   if (nat_in_frag(l)) then
                       charge = charge + real(smmd%nelpsp(smmd%iatype(l)),kind=8)
                   end if
               end do
               !write(*,'(a,2i5,12(f7.3))') 'iel, i, weights', iel, i, weights
               !write(*,*) nat_in_frag, charge
               !write(*,*) 'purity indincator',  get_subkernel_purity(smat(1)%nfvctr, kernel_ortho%matrix, in_frag, charge)
               purity_indicator = get_subkernel_purity(smat(1)%nfvctr, kernel_ortho%matrix, in_frag, charge)
               if (iproc==0) then
                   call yaml_sequence(advance='no')
                   call yaml_mapping_open('frag')
                   call yaml_map('number of target elements',iel)
                   call yaml_map('fragment search',i)
                   call yaml_map('nat_in_frag',nat_frag)
                   call yaml_map('atom list',natlist(1:nat_frag))
                   call yaml_map('fragment charge',charge)
                   call yaml_map('purity indicator',purity_indicator)
                   call yaml_mapping_close()
               end if

               ! Check whether this fragment has already been found before.
               do ifrag=1,nfragfound
                   found = .true.
                   if (nat_frag /= natfrag(ifrag)) then
                       ! For sure a new fragment as the number of atoms is different
                       found = .false.
                       cycle
                   end if
                   do l=1,nat_frag
                       if (natlist(l) /= fraglist(l,ifrag)) then
                           found = .false.
                           cycle
                       end if
                   end do
                   if (found) then
                       ! This fragment has already been found
                       exit
                   end if
               end do
               if (.not.found) then
                   write(*,*) 'new fragment'
                   nfragfound = nfragfound + 1
                   if (nfragfound>nfragmax) then
                       call f_err_throw('maximal number of fragments found')
                   end if
                   natfrag(nfragfound) = nat_frag
                   call f_memcpy(src=natlist(1:nat_frag), dest=fraglist(1:nat_frag,ifrag))
                   purity_indicator_list(nfragfound) = purity_indicator
               else
                   write(*,*) 'FOUND fragment'
               end if
           end do
       end do

       if (iproc==0) then
           call yaml_sequence_close()
       end if

   call build_dict_info(dict_timing_info)
       if (iproc==0) then
           call yaml_sequence_open('List of calculated fragments')
           do ifrag=1,nfragfound
               call yaml_sequence(advance='no')
               call yaml_mapping_open('Fragment')
               call yaml_map('Number of atoms in fragment',natfrag(ifrag))
               call yaml_map('IDs of atom in fragment',fraglist(1:natfrag(ifrag),ifrag))
               call yaml_map('purity indicator',purity_indicator_list(ifrag))
               call yaml_mapping_close()
           end do
       end if

       call f_free(at_assigned)
       call f_free(in_frag)
       call f_free(nat_in_frag)
       call f_free(weights)
       call f_free(fvec)
       call f_free(diag)
       call f_free(fjac)
       call f_free(ipvt)
       call f_free(qtf)
       call f_free(wa1)
       call f_free(wa2)
       call f_free(wa3)
       call f_free(wa4)
       call deallocate_K()
       call deallocate_xold()
       call f_free(natlist)
       call f_free(natfrag)
       call f_free(fraglist)
       call f_free(purity_indicator_list)

       call deallocate_sparse_matrix(smat(1))
       call deallocate_sparse_matrix(smat(2))
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(kernel_mat)
       call deallocate_matrices(kernel_ortho)
       call deallocate_sparse_matrix_metadata(smmd)

   end if


   call f_timing_stop(mpi_comm=mpiworld(),nproc=nproc,&
        gather_routine=gather_timings,dict_info=dict_timing_info)
   call dict_free(dict_timing_info)

   if (iproc==0) then
       call yaml_release_document()
   end if


   !!call bigdft_finalize(ierr)
   call mpifinalize()

   call f_lib_finalize()




end program chess_toolbox

!!$!> extract the different wavefunctions to verify if the completeness relation is satisfied
!!$subroutine completeness_relation
!!$
!!$  call wfn_filename(filename_out,radical,binary,ikpt,nspinor,nspin,ispinor,spin,iorb)
!!$
!!$  !loop that has to be done for each of the wavefunctions
!!$  unitwf=99
!!$  call f_open_file(unit=unitwf,file=filename_out)
!!$  call readonewave(unitwf,.not. binary,iorb,iproc,&
!!$       it%lr%d%n1,it%lr%d%n2,it%lr%d%n3, &
!!$       Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
!!$       at,it%lr%wfd,rxyz_old,rxyz,&
!!$       psi_ptr,eval,psifscf)
!!$  call f_close(unitwf)
!!$
!!$end subroutine completeness_relation






subroutine extract_fragment_submatrix(smmd, smat, mat, nat_frag, nfvctr_frag, &
           fragment_atom_id, fragment_supfun_id, mat_frag)
  use sparsematrix_base
  !!use wrapper_mpi
  implicit none

  ! Calling arguments
  type(sparse_matrix_metadata),intent(in) :: smmd
  type(sparse_matrix),intent(in) :: smat
  type(matrices),intent(in) :: mat
  integer,intent(in) :: nat_frag, nfvctr_frag
  integer,dimension(nat_frag),intent(in) :: fragment_atom_id
  integer,dimension(smat%nfvctr),intent(in) :: fragment_supfun_id
  real(mp),dimension(nfvctr_frag,nfvctr_frag),intent(out) :: mat_frag

  ! Local variables
  integer :: iicol, icol, icol_atom, iat, iiat, ii, iirow, i, irow, irow_atom, iseg, icol_old
  logical :: found_icol, found_irow

  iicol = 0
  icol_old = -1
  seg_loop: do iseg=1,smat%nseg
      icol = smat%keyg(1,2,iseg)
      icol_atom = smmd%on_which_atom(icol)
      ! Search whether this column belongs to the fragment
      found_icol = .false.
      do iat=1,nat_frag
          iiat = fragment_atom_id(iat)
          if (icol_atom==iiat) then
              found_icol = .true.
          end if
      end do
      if (found_icol) then
          if (icol/=icol_old) then
              iicol = iicol + 1
              icol_old = icol
          end if
      else
          cycle seg_loop
      end if
      ii=smat%keyv(iseg)
      iirow = 0
      do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg) 
          irow = i
          irow_atom = smmd%on_which_atom(irow)
          ! Search whether this column belongs to the fragment
          found_irow = .false.
          do iat=1,nat_frag
              iiat = fragment_atom_id(iat)
              if (irow_atom==iiat) then
                  found_irow = .true.
              end if
          end do
          if (found_irow .and. found_icol) then
              iirow = iirow + 1
              !mat_frag(iirow,iicol) = mat%matrix_compr(ii)
              iirow = fragment_supfun_id(irow)
              iicol = fragment_supfun_id(icol)
              !!if (mpirank()==0) write(*,'(a,5i8,es16.6)') 'irow, icol, iirow, iicol, ii, val', irow, icol, iirow, iicol, ii, mat%matrix_compr(ii)
              mat_frag(iirow,iicol) = mat%matrix_compr(ii-smat%isvctrp_tg)
          end if
          !write(*,*) 'iirow, iicol, ii, mat_frag, mat_compr', iirow, iicol, ii, mat_frag(iirow,iicol), mat%matrix_compr(ii)
          ii = ii + 1
      end do
  end do seg_loop

end subroutine extract_fragment_submatrix




function get_minimal_distance(ra, rb, cell, periodic) result(d)
  use sparsematrix_base
  implicit none

  ! Calling arguments
  real(mp),dimension(3),intent(in) :: ra, rb, cell
  logical,dimension(3),intent(in) :: periodic
  real(mp) :: d

  ! Local variables
  integer :: i
  real(mp) :: di, shift

  d = 0._mp
  do i=1,3
      ! Distance in reduced coordinates
      di = abs(ra(i)-rb(i))/cell(i)
      !write(*,*) 'i, di', i, di
      if (periodic(i)) then
          !If the distance is larger than half the box size 
          !(which is 1 due to the reduced coordinates), 
          !we have to make a periodic wrap around
          shift=real(floor(di+0.5_mp),kind=mp)
      else
          shift = 0._mp
      end if
      ! Convert back to the corect cell size and add it up
      d = d + (cell(i)*(di-shift))**2
  end do
  d = sqrt(d)

end function get_minimal_distance



