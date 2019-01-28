module bigdft_matrices
  implicit none

  private

  public :: check_local_matrix_extents
  !!public :: get_modulo_array
  public :: init_matrixindex_in_compressed_fortransposed
  public :: init_bigdft_matrices


  contains

    subroutine check_local_matrix_extents(iproc, nproc, collcom, collcom_sr, orbs, &
               smmd, smat, aux, ind_min, ind_max)
          use module_base
          use sparsematrix_base, only: sparse_matrix, sparse_matrix_metadata
          use sparsematrix_init, only: get_sparsematrix_local_rows_columns, &
                                       check_projector_charge_analysis
          use communications_base, only: comms_linear
          use module_types, only: linmat_auxiliary, orbitals_data
          implicit none
    
          ! Caling arguments
          integer,intent(in) :: iproc, nproc
          type(comms_linear),intent(in) :: collcom, collcom_sr
          type(orbitals_data),intent(in) :: orbs
          type(sparse_matrix_metadata),intent(in) :: smmd
          type(sparse_matrix),intent(in) :: smat
          type(linmat_auxiliary),intent(in) :: aux
          integer,intent(out) :: ind_min, ind_max
    
          ! Local variables
          integer :: i, ii_ref, iorb, jorb, ii, iseg
          logical :: found

          real(kind=4) :: tr0, tr1, trt0, trt1
          real(kind=8) :: time0, time1, time2, time3, time4, time5, ttime
          logical, parameter :: extra_timing=.false.
          integer,dimension(:),pointer :: moduloarray
                        
    
          call timing(iproc,'matrix_extents','ON')
          call f_routine(id='check_local_matrix_extents')
          if (extra_timing) call cpu_time(trt0)  

          ind_min = smat%nvctr
          ind_max = 0

          !!call get_sparsematrix_local_extent(iproc, nproc, smmd, smat, ind_min, ind_max)

          if (extra_timing) call cpu_time(tr0)
          ! The operations done in the transposed wavefunction layout
          !call check_transposed_layout()
          !call get_modulo_array(smat%nfvctr, aux%mat_ind_compr2, moduloarray)
          !write(*,*) 'moduloarray',moduloarray
          !write(1000+iproc,*) 'calling find_minmax_transposed'
          call find_minmax_transposed(aux%mat_ind_compr2,collcom,smat%nfvctr,ind_min,ind_max)

          !write(*,'(a,2i8)') 'after check_transposed_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr1)
          if (extra_timing) time0=real(tr1-tr0,kind=8)    


          ! Now check the compress_distributed layout
          !call check_compress_distributed_layout()
          !call check_compress_distributed_layout(smat,ind_min,ind_max)

          !write(*,'(a,2i8)') 'after check_compress_distributed_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr0)
          if (extra_timing) time1=real(tr0-tr1,kind=8)        

          ! Now check the matrix matrix multiplications layout
          !!if (smat%smatmul_initialized) then
          !!    call check_matmul_layout(smat%smmm%nseq,smat%smmm%indices_extract_sequential,ind_min,ind_max)
          !!end if
          !write(*,'(a,2i8)') 'after check_matmul_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr1)
          if (extra_timing) time2=real(tr1-tr0,kind=8)    
    
          ! Now check the sumrho operations
          !call check_sumrho_layout()
          !write(1000+iproc,*) 'calling check_sumrho_layout'
          call check_sumrho_layout(collcom_sr,smat%nfvctr,aux%mat_ind_compr2,ind_min,ind_max)

          !call f_free_ptr(moduloarray)
          !write(*,'(a,2i8)') 'after check_sumrho_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr0)
          if (extra_timing) time3=real(tr0-tr1,kind=8)    
    
          !!!! Now check the pseudo-exact orthonormalization during the input guess
          !call check_ortho_inguess()
          !!call check_ortho_inguess(smat,ind_min,ind_max)
          !write(*,'(a,2i8)') 'after check_ortho_inguess: ind_min, ind_max', ind_min, ind_max
          !write(1000+iproc,*) 'calling check_projector_charge_analysis'
          call check_projector_charge_analysis(iproc, nproc, smmd, smat, ind_min, ind_max)
          if (extra_timing) call cpu_time(tr1)
          if (extra_timing) time4=real(tr1-tr0,kind=8)        

          ! Now check the submatrix extraction for the projector charge analysis
          !!call check_projector_charge_analysis(iproc, nproc, smmd, smat, ind_min, ind_max)

          !!write(*,'(a,3i8)') 'after check_local_matrix_extents: iproc, ind_min, ind_max', iproc, ind_min, ind_max

          !!call get_sparsematrix_local_rows_columns(smat, ind_min, ind_max, irow, icol)

          !!! Get the global indices of ind_min and ind_max
          !!do i=1,2
          !!    if (i==1) then
          !!        ii_ref = ind_min
          !!    else
          !!        ii_ref = ind_max
          !!    end if
          !!    ! Search the indices iorb,jorb corresponding to ii_ref
          !!    found=.false.

          !!    ! not sure if OpenMP is really worth it here
          !!    !$omp parallel default(none) &
          !!    !$omp private(iseg,ii,iorb,jorb) &
          !!    !$omp shared(smat,ii_ref,irow,icol,found,i)
          !!    !$omp do
          !!    outloop: do iseg=1,smat%nseg
          !!        if (.not. found) then
          !!           iorb = smat%keyg(1,2,iseg)
          !!           do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
          !!               ii = matrixindex_in_compressed(smat, jorb, iorb)
          !!               !if (iproc==0) write(*,'(a,5i9)') 'i, ii_ref, ii, iorb, jorb', i, ii_ref, ii, iorb, jorb
          !!               if (ii==ii_ref) then
          !!                   irow(i) = jorb
          !!                   icol(i) = iorb
          !!                   !exit outloop
          !!                   !SM: I think one should do this within a critical section since it is shared, just to be sure...
          !!                   !$omp critical
          !!                   found=.true.
          !!                   !$omp end critical
          !!               end if
          !!           end do
          !!        end if
          !!    end do outloop
          !!    !$omp end do
          !!    !$omp end parallel

          !!end do
          if (extra_timing) call cpu_time(tr0)
          if (extra_timing) time5=real(tr0-tr1,kind=8)    
          if (extra_timing) call cpu_time(trt1)  
          if (extra_timing) ttime=real(trt1-trt0,kind=8)  

          !write(1000+iproc,*) 'calling check_orbital_matrix_distribution'
          call check_orbital_matrix_distribution(orbs, smat, ind_min, ind_max)

          if (extra_timing.and.iproc==0) print*,'matextent',time0,time1,time2,time3,time4,time5,&
               time0+time1+time2+time3+time4+time5,ttime

          call f_release_routine()
          call timing(iproc,'matrix_extents','OF')    
    
    end subroutine check_local_matrix_extents


    subroutine find_minmax_transposed(mat_ind_compr2,collcom,nfvctr,ind_min,ind_max)
      use futile
      use communications_base, only: comms_linear
      use module_types, only: matrixindex_in_compressed_fortransposed2
      implicit none
      integer, intent(in) :: nfvctr
      type(comms_linear),intent(in) :: collcom
      !integer, dimension(:,:), intent(in) :: matrixindex_in_compressed_fortransposed
      type(matrixindex_in_compressed_fortransposed2),dimension(nfvctr),intent(in) :: mat_ind_compr2
      !integer, dimension(nfvctr), intent(in) :: moduloarray
      integer, intent(inout) :: ind_min,ind_max
      !local variables
      integer :: ipt, ii, i0, i, i0i, iiorb, j, i0j, jjorb, ind, iorb, jorb, ia, ib

      call f_routine(id='find_minmax_transposed')

      !$omp parallel default(none) &
      !$omp private(ipt,ii,i0,i,i0i,iiorb,iorb,j,i0j,jjorb,jorb,ind,ia,ib) &
      !$omp shared(collcom,ind_min,ind_max,mat_ind_compr2,nfvctr)
      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do ipt=1,collcom%nptsp_c
         ii=collcom%norb_per_gridpoint_c(ipt)
         i0 = collcom%isptsp_c(ipt)
         do i=1,ii
            i0i=i0+i
            iiorb=collcom%indexrecvorbital_c(i0i)
            !iorb=moduloarray(iiorb)
            iorb = modulo(iiorb-mat_ind_compr2(iiorb)%offset_compr,nfvctr)+1
            do j=1,ii
               i0j=i0+j
               jjorb=collcom%indexrecvorbital_c(i0j)
               !jorb=moduloarray(jjorb)
               !jorb = modulo(jjorb-mat_ind_compr(iiorb)%offset_compr,nfvctr)+1
               !ind = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
               !ind = matrixindex_in_compressed_fortransposed(jorb,iorb)
               !!if (jorb>ubound(mat_ind_compr(iiorb)%ind_compr,1)) then
               !!    write(*,*) 'ipt, ii, i, j, iiorb, iorb, jjorb, jorb, lbound, ubound', &
               !!                ipt, ii, i, j, iiorb, iorb, jjorb, jorb, &
               !!                lbound(mat_ind_compr(iiorb)%ind_compr,1), &
               !!                ubound(mat_ind_compr(iiorb)%ind_compr,1)
               !!end if
               !ind = mat_ind_compr(iiorb)%ind_compr(jorb)
               ia = jjorb-mat_ind_compr2(iiorb)%offset_compr
               ib = sign(1,ia)
               ind = mat_ind_compr2(iiorb)%section(ib)%ind_compr(jjorb)
               ind_min = min(ind_min,ind)
               ind_max = max(ind_max,ind)
            end do
         end do
      end do
      !$omp end do

      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do ipt=1,collcom%nptsp_f
         ii=collcom%norb_per_gridpoint_f(ipt)
         i0 = collcom%isptsp_f(ipt)
         do i=1,ii
            i0i=i0+i
            iiorb=collcom%indexrecvorbital_f(i0i)
            !iorb=moduloarray(iiorb)
            do j=1,ii
               i0j=i0+j
               jjorb=collcom%indexrecvorbital_f(i0j)
               !jorb=moduloarray(jjorb)
               !!jorb = modulo(jjorb-mat_ind_compr(iiorb)%offset_compr,nfvctr)+1
               !ind = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
               !ind = matrixindex_in_compressed_fortransposed(jorb,iorb)
               !!ind = mat_ind_compr(iiorb)%ind_compr(jorb)
               ia = jjorb-mat_ind_compr2(iiorb)%offset_compr
               ib = sign(1,ia)
               ind = mat_ind_compr2(iiorb)%section(ib)%ind_compr(jjorb)
               ind_min = min(ind_min,ind)
               ind_min = min(ind_min,ind)
               ind_max = max(ind_max,ind)
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel

      call f_release_routine()

    end subroutine find_minmax_transposed






    subroutine check_sumrho_layout(collcom_sr,nfvctr,mat_ind_compr2,ind_min,ind_max)
      use futile
      use communications_base, only: comms_linear
      use sparsematrix_init, only: matrixindex_in_compressed
      use module_types, only: matrixindex_in_compressed_fortransposed2
      implicit none
      integer, intent(in) :: nfvctr
      type(comms_linear),intent(in) :: collcom_sr
      !integer, dimension(:,:), intent(in) :: matrixindex_in_compressed_fortransposed
      type(matrixindex_in_compressed_fortransposed2),dimension(nfvctr),intent(in) :: mat_ind_compr2
      !integer, dimension(nfvctr), intent(in) :: moduloarray
      integer, intent(inout) :: ind_min,ind_max
      !local variables
      integer :: ipt, ii, i0, i, iiorb, ind, iorb, ia, ib

      call f_routine(id='check_sumrho_layout')

      !$omp parallel default(none) &
      !$omp private(ipt,ii,i0,iiorb,iorb,ind,i,ia,ib) &
      !$omp shared(collcom_sr,mat_ind_compr2,ind_min,ind_max)
      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do ipt=1,collcom_sr%nptsp_c
         ii=collcom_sr%norb_per_gridpoint_c(ipt)
         i0=collcom_sr%isptsp_c(ipt)
         do i=1,ii
            iiorb=collcom_sr%indexrecvorbital_c(i0+i)
            !iorb=moduloarray(iiorb)
            ia = iiorb-mat_ind_compr2(iiorb)%offset_compr
            ib = sign(1,ia)
            ind = mat_ind_compr2(iiorb)%section(ib)%ind_compr(iiorb)
            !ind=smat%matrixindex_in_compressed_fortransposed(iiorb,iiorb)
            !ind=matrixindex_in_compressed_fortransposed(iorb,iorb)
            !ind = mat_ind_compr(iiorb)%ind_compr(iorb)
            !ind=get_transposed_index(smat,iiorb,iiorb)
            ind_min = min(ind_min,ind)
            ind_max = max(ind_max,ind)
         end do
      end do
      !$omp end do
      !$omp end parallel

      call f_release_routine()

    end subroutine check_sumrho_layout



    !!subroutine get_modulo_array(nfvctr, mat_ind_compr, moduloarray)
    !!  use module_base
    !!  use module_types, only: matrixindex_in_compressed_fortransposed2
    !!  implicit none
    !!  ! Calling arguments
    !!  integer,intent(in) :: nfvctr
    !!  type(matrixindex_in_compressed_fortransposed),dimension(nfvctr),intent(in) :: mat_ind_compr
    !!  integer,dimension(:),pointer :: moduloarray
    !!  ! Local variables
    !!  integer :: i
    !!  moduloarray = f_malloc_ptr(nfvctr,id='moduloarray')
    !!  !$omp parallel default(none) &
    !!  !$omp shared(moduloarray,nfvctr, mat_ind_compr) &
    !!  !$omp private(i)
    !!  !$omp do
    !!  do i=1,nfvctr
    !!      moduloarray(i) = modulo(i-mat_ind_compr(i)%offset_compr,nfvctr)+1
    !!  end do
    !!  !$omp end do
    !!  !$omp end parallel
    !!end subroutine get_modulo_array


    subroutine init_matrixindex_in_compressed_fortransposed(iproc, nproc, collcom, collcom_shamop, &
               collcom_sr, sparsemat, aux)
      use module_base
      use sparsematrix_base, only: sparse_matrix
      use communications_base, only: comms_linear
      use sparsematrix_init, only: matrixindex_in_compressed
      use module_types, only: linmat_auxiliary, matrixindex_in_compressed_fortransposed2_null, linmat_auxiliary_null
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(comms_linear),intent(in) :: collcom, collcom_shamop, collcom_sr
      type(sparse_matrix),intent(inout) :: sparsemat
      type(linmat_auxiliary),intent(inout) :: aux
      
      ! Local variables
      integer :: iorb, jorb, istat, imin, imax, nmiddle, iiorb, jjorb
      integer :: ii, i, nlen, j, ifvctr
      integer,dimension(:),allocatable :: imin_old, imax_old, imin_new, imax_new
      !integer :: kproc,jproc,jjorbold,jjorb,isend,irecv,ilr,ijorb,iiorb,ind,ierr, irow, irowold, iseg
      !integer :: compressed_index
    !  integer,dimension(:,:),allocatable :: sendbuf, requests, iminmaxarr
    
      call f_routine(id='init_matrixindex_in_compressed_fortransposed')
      call timing(iproc,'init_matrCompr','ON')
    
      ! for the calculation of overlaps and the charge density
      !imin=minval(collcom%indexrecvorbital_c)
      !imin=min(imin,minval(collcom%indexrecvorbital_f))
      !imin=min(imin,minval(collcom_shamop%indexrecvorbital_c))
      !imin=min(imin,minval(collcom_shamop%indexrecvorbital_f))
      !imin=min(imin,minval(collcom_sr%indexrecvorbital_c))
      !imax=maxval(collcom%indexrecvorbital_c)
      !imax=max(imax,maxval(collcom%indexrecvorbital_f))
      !imax=max(imax,maxval(collcom_shamop%indexrecvorbital_c))
      !imax=max(imax,maxval(collcom_shamop%indexrecvorbital_f))
      !imax=max(imax,maxval(collcom_sr%indexrecvorbital_c))

      aux = linmat_auxiliary_null()

      imin_old = f_malloc(sparsemat%nfvctr,id='imin_old')
      imax_old = f_malloc(sparsemat%nfvctr,id='imax_old')
      imin_new = f_malloc(sparsemat%nfvctr,id='imin_new')
      imax_new = f_malloc(sparsemat%nfvctr,id='imax_new')
    
      nmiddle = sparsemat%nfvctr/2 + 1
    
      imin_old(:) = huge(1)
      imax_old(:) = 0
      imin_new(:) = huge(1)
      imax_new(:) = 0


      call get_minmax_lines(collcom, sparsemat, 'c', imin_old, imax_old, imin_new, imax_new)
      call get_minmax_lines(collcom, sparsemat, 'f', imin_old, imax_old, imin_new, imax_new)
      call get_minmax_lines(collcom_shamop, sparsemat, 'c', imin_old, imax_old, imin_new, imax_new)
      call get_minmax_lines(collcom_shamop, sparsemat, 'f', imin_old, imax_old, imin_new, imax_new)
      call get_minmax_lines(collcom_sr, sparsemat, 'c', imin_old, imax_old, imin_new, imax_new)

      !!allocate(aux%mat_ind_compr(sparsemat%nfvctr))
      allocate(aux%mat_ind_compr2(sparsemat%nfvctr))
      do ifvctr=1,sparsemat%nfvctr
          ! Determine with which size the array should be allocated
          if (imax_new(ifvctr)-imin_new(ifvctr)<0) then
              ! everything in either first or second half
              imin = imin_old(ifvctr)
              imax = imax_old(ifvctr)
          else
              ! in both half
              if (imax_old(ifvctr)-imin_old(ifvctr)>imax_new(ifvctr)-imin_new(ifvctr)) then
                  ! wrap around
                  imin = imin_new(ifvctr)
                  imax = imax_new(ifvctr)
              else
                  ! no wrap around
                  imin = imin_old(ifvctr)
                  imax = imax_old(ifvctr)
              end if
          end if
          !!aux%mat_ind_compr(ifvctr) = matrixindex_in_compressed_fortransposed_null()
          aux%mat_ind_compr2(ifvctr) = matrixindex_in_compressed_fortransposed2_null()
          !!aux%mat_ind_compr(ifvctr)%offset_compr = imin
          aux%mat_ind_compr2(ifvctr)%offset_compr = imin
          nlen = imax - imin + 1
          !!write(*,*) 'iproc, ifvctr, imin_old(ifvctr), imax_old(ifvctr), imin_new(ifvctr), imax_new(ifvctr), imin, imax, nlen',&
          !!            iproc, ifvctr, imin_old(ifvctr), imax_old(ifvctr), imin_new(ifvctr), imax_new(ifvctr), imin, imax, nlen
          !!aux%mat_ind_compr(ifvctr)%ind_compr = f_malloc_ptr(nlen,id='aux%linmat%mat_ind_compr%ind_compr')
          if (min(imax,sparsemat%nfvctr)>=imin) then
              aux%mat_ind_compr2(ifvctr)%section(1)%ind_compr = &
                  f_malloc_ptr(imin.to.min(imax,sparsemat%nfvctr),id='aux%linmat%mat_ind_compr%ind_compr')
          end if
          if (imax-sparsemat%nfvctr>=1) then
              aux%mat_ind_compr2(ifvctr)%section(-1)%ind_compr = &
                  f_malloc_ptr(1.to.imax-sparsemat%nfvctr,id='aux%linmat%mat_ind_compr%ind_compr')
          end if
          !$omp parallel do default(private) shared(sparsemat,aux,imin,imax,ifvctr)
          do jorb=imin,imax
              j = jorb - imin + 1
              jjorb = mod(jorb-1,sparsemat%nfvctr)+1
              !aux%mat_ind_compr(ifvctr)%ind_compr(j)=matrixindex_in_compressed(sparsemat, ifvctr, jjorb)
              !!aux%mat_ind_compr(ifvctr)%ind_compr(j)=matrixindex_in_compressed(sparsemat, jjorb, ifvctr)
              if (jorb<=sparsemat%nfvctr) then
                  aux%mat_ind_compr2(ifvctr)%section(1)%ind_compr(jorb) = &
                      matrixindex_in_compressed(sparsemat, jjorb, ifvctr)
              else
                  aux%mat_ind_compr2(ifvctr)%section(-1)%ind_compr(jorb-sparsemat%nfvctr) = &
                      matrixindex_in_compressed(sparsemat, jjorb, ifvctr)
              end if
          end do
          !$omp end parallel do
      end do

      call f_free(imin_old)
      call f_free(imax_old)
      call f_free(imin_new)
      call f_free(imax_new)

      !!!do i=1,size(collcom%indexrecvorbital_c)
      !!!    ii = mod(collcom%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
      !!!    imin_old = min(imin_old,ii)
      !!!    imax_old = max(imax_old,ii)
      !!!    if (ii>nmiddle) then
      !!!        imin_new = min(imin_new,ii)
      !!!    else
      !!!        imax_new = max(imax_new,ii+sparsemat%nfvctr)
      !!!    end if
      !!!end do
      !!!do i=1,size(collcom%indexrecvorbital_f)
      !!!    ii = mod(collcom%indexrecvorbital_f(i)-1,sparsemat%nfvctr)+1
      !!!    imin_old = min(imin_old,ii)
      !!!    imax_old = max(imax_old,ii)
      !!!    if (ii>nmiddle) then
      !!!        imin_new = min(imin_new,ii)
      !!!    else
      !!!        imax_new = max(imax_new,ii+sparsemat%nfvctr)
      !!!    end if
      !!!end do
      !!!do i=1,size(collcom_shamop%indexrecvorbital_c)
      !!!    ii = mod(collcom_shamop%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
      !!!    imin_old = min(imin_old,ii)
      !!!    imax_old = max(imax_old,ii)
      !!!    if (ii>nmiddle) then
      !!!        imin_new = min(imin_new,ii)
      !!!    else
      !!!        imax_new = max(imax_new,ii+sparsemat%nfvctr)
      !!!    end if
      !!!end do
      !!!do i=1,size(collcom_shamop%indexrecvorbital_f)
      !!!    ii = mod(collcom_shamop%indexrecvorbital_f(i)-1,sparsemat%nfvctr)+1
      !!!    imin_old = min(imin_old,ii)
      !!!    imax_old = max(imax_old,ii)
      !!!    if (ii>nmiddle) then
      !!!        imin_new = min(imin_new,ii)
      !!!    else
      !!!        imax_new = max(imax_new,ii+sparsemat%nfvctr)
      !!!    end if
      !!!end do
      !!!do i=1,size(collcom_sr%indexrecvorbital_c)
      !!!    ii = mod(collcom_sr%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
      !!!    imin_old = min(imin_old,ii)
      !!!    imax_old = max(imax_old,ii)
      !!!    if (ii>nmiddle) then
      !!!        imin_new = min(imin_new,ii)
      !!!    else
      !!!        imax_new = max(imax_new,ii+sparsemat%nfvctr)
      !!!    end if
      !!!end do
    
    
      !!write(*,*) 'iproc, imin_old, imax_old', iproc, imin_old, imax_old
      !!write(*,*) 'iproc, imin_new, imax_new', iproc, imin_new, imax_new
    
      !! values regardless of the spin
      !imin=mod(imin-1,sparsemat%nfvctr)+1
      !imax=mod(imax-1,sparsemat%nfvctr)+1
    
    
      !!! Determine with which size the array should be allocated
      !!if (imax_new-imin_new<0) then
      !!    ! everything in either first or second half
      !!    imin = imin_old
      !!    imax = imax_old
      !!    !sparsemat%offset_matrixindex_in_compressed_fortransposed = 1
      !!else
      !!    ! in both half
      !!    if (imax_old-imin_old>imax_new-imin_new) then
      !!        ! wrap around
      !!        imin = imin_new
      !!        imax = imax_new
      !!        !sparsemat%offset_matrixindex_in_compressed_fortransposed = imin_new
      !!    else
      !!        ! no wrap around
      !!        imin = imin_old
      !!        imax = imax_old
      !!        !sparsemat%offset_matrixindex_in_compressed_fortransposed = 1
      !!    end if
      !!end if
    
      !!! Check
      !!if (sparsemat%offset_matrixindex_in_compressed_fortransposed<sparsemat%nfvctr/2+1) then
      !!    stop 'sparsemat%offset_matrixindex_in_compressed_fortransposed<sparsemat%nfvctr/2+1'
      !!end if
    
      !!nlen = imax - imin + 1
      !!aux%offset_matrixindex_in_compressed_fortransposed = imin
    
      !!! This is a temporary solution for spin polarized systems
      !!imax=min(imax,orbs%norbu)
    
    
    
      !!allocate(sparsemat%matrixindex_in_compressed_fortransposed(imin:imax,imin:imax), stat=istat)
      !!call memocc(istat, sparsemat%matrixindex_in_compressed_fortransposed, &
      !sparsemat%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/imin.to.imax,imin.to.imax/),&
      !    id='sparsemat%matrixindex_in_compressed_fortransposed')
      !!sparsemat%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/nlen,nlen/),&
      !!    id='sparsemat%matrixindex_in_compressed_fortransposed')
      !!aux%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/nlen,nlen/),id='aux%matrixindex_in_compressed_fortransposed')
    
      !!!$omp parallel do default(private) shared(sparsemat,aux,imin,imax)
      !!do iorb=imin,imax
      !!    i = iorb - imin + 1
      !!    do jorb=imin,imax
      !!        j = jorb - imin + 1
      !!        !@ii=(jorb-1)*sparsemat%nfvctr+iorb
      !!        !@ispin=(ii-1)/sparsemat%nfvctr+1 !integer division to get the spin (1 for spin up (or non polarized), 2 for spin down)
      !!        !@iiorb=mod(iorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
      !!        !@jjorb=mod(jorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
      !!        !sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=compressed_index(iiorb,jjorb,orbs%norbu,sparsemat)
      !!        iiorb = mod(iorb-1,sparsemat%nfvctr)+1
      !!        jjorb = mod(jorb-1,sparsemat%nfvctr)+1
      !!        !sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=matrixindex_in_compressed(sparsemat, iiorb, jjorb)
      !!        !sparsemat%matrixindex_in_compressed_fortransposed(i,j)=matrixindex_in_compressed(sparsemat, iiorb, jjorb)
      !!        aux%matrixindex_in_compressed_fortransposed(i,j)=matrixindex_in_compressed(sparsemat, iiorb, jjorb)
      !!        !sendbuf(jorb,iorb)=compressed_index(jorb,iiorb,orbs%norb,sparsemat)
      !!        !sendbuf(iorb,jorb)=compressed_index(iiorb,jorb,orbs%norb,sparsemat)
      !!    end do
      !!end do
      !!!$omp end parallel do
    
      !@! Add the spin shift (i.e. the index is in the spin polarized matrix which is at the end)
      !@if (ispin==2) then
      !@    matrixindex_in_compressed = matrixindex_in_compressed + sparsemat%nvctr
      !@end if
    
      call timing(iproc,'init_matrCompr','OF')
      call f_release_routine()
    
    end subroutine init_matrixindex_in_compressed_fortransposed


    subroutine init_bigdft_matrices(iproc, nproc, atoms, orbs, &
               collcom_s, collcom_s_sr, collcom_m, lzd_s, lzd_m, in, linmat)
      use module_base
      use module_types
      use sparsematrix_wrappers, only: init_sparse_matrix_wrapper, check_kernel_cutoff
      use sparsematrix_init, only: sparse_matrix_metadata_init, init_matrix_taskgroups_wrapper
      use yaml_output
      use sparsematrix_memory, only: sparse_matrix_null
      use at_domain, only: domain_geocode
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(atoms_data),intent(in) :: atoms
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom_s, collcom_s_sr, collcom_m
      type(local_zone_descriptors),intent(inout) :: lzd_s
      type(local_zone_descriptors),intent(in) :: lzd_m
      type(input_variables), intent(inout) :: in
      type(linear_matrices),intent(inout),target :: linmat

      ! Local variables
      integer,dimension(2) :: irow, icol, iirow, iicol
      integer :: ind_min_s, ind_max_s
      integer :: ind_min_m, ind_max_m
      integer :: ind_min_l, ind_max_l
      type(sparse_matrix) :: smat_test
      integer,dimension(:,:),allocatable :: ind_minmax

      call f_routine(id='init_bigdft_matrices')

      call sparse_matrix_metadata_init(domain_geocode(atoms%astruct%dom), atoms%astruct%cell_dim, orbs%norb, &
           atoms%astruct%nat, atoms%astruct%ntypes, atoms%astruct%units, &           
           atoms%nzatom, atoms%nelpsp, atoms%astruct%atomnames, atoms%astruct%iatype, &
           atoms%astruct%rxyz, orbs%onwhichatom, linmat%smmd)

      ! Do not initialize the matrix multiplication to save memory. The multiplications
      ! are always done with the linmat%smat(3) type.
      call init_sparse_matrix_wrapper(iproc, nproc, &
           in%nspin, orbs, lzd_s, atoms%astruct, &
           in%store_index, init_matmul=.false., matmul_optimize_load_balancing=.false., &
           imode=1, smat=linmat%smat(1))
      call init_matrixindex_in_compressed_fortransposed(iproc, nproc, &
           collcom_s, collcom_m, collcom_s_sr, linmat%smat(1), &
           linmat%auxs)

      ! Do not initialize the matrix multiplication to save memory. The multiplications
      ! are always done with the linmat%smat(3) type.
      call init_sparse_matrix_wrapper(iproc, nproc, &
           in%nspin, orbs, lzd_m, atoms%astruct, &
           in%store_index, init_matmul=.false., matmul_optimize_load_balancing=.false., &
           imode=1, smat=linmat%smat(2))
      call init_matrixindex_in_compressed_fortransposed(iproc, nproc, &
           collcom_s, collcom_m, collcom_s_sr, linmat%smat(2), &
           linmat%auxm)

      ! check the extent of the kernel cutoff (must be at least shamop radius)
      if (iproc==0) then
          call yaml_comment('Sparse matrix initialization',hfill='-')
      end if
      call check_kernel_cutoff(iproc, orbs, atoms, in%hamapp_radius_incr, lzd_s, only_check=.true.)
      call init_sparse_matrix_wrapper(iproc, nproc, &
           in%nspin, orbs, lzd_s, atoms%astruct, &
           in%store_index, init_matmul=.true., matmul_optimize_load_balancing=in%cp%foe%matmul_optimize_load_balancing, &
           imode=2, smat=linmat%smat(3), smat_ref=linmat%smat(2))
      !write(1000+iproc,*) 'calling init_matrixindex_in_compressed_fortransposed'
      call init_matrixindex_in_compressed_fortransposed(iproc, nproc, &
           collcom_s, collcom_m, collcom_s_sr, linmat%smat(3), &
           linmat%auxl)


      !write(1000+iproc,*) 'calling check_local_matrix_extents s'
      call check_local_matrix_extents(iproc, nproc, collcom_s, &
           collcom_s_sr, orbs, linmat%smmd, linmat%smat(1), linmat%auxs, &
           ind_min_s, ind_max_s)
      !write(1000+iproc,*) 'calling check_local_matrix_extents m'
      call check_local_matrix_extents(iproc, nproc, collcom_m, &
           collcom_s_sr, orbs, linmat%smmd, linmat%smat(2), linmat%auxm, &
           ind_min_m, ind_max_m)
      !write(1000+iproc,*) 'calling check_local_matrix_extents l'
      call check_local_matrix_extents(iproc, nproc, collcom_m, &
           collcom_s_sr, orbs, linmat%smmd, linmat%smat(3), linmat%auxl, &
           ind_min_l, ind_max_l)

      ind_minmax = f_malloc([2,3],id='ind_minmax')
      ind_minmax(1:2,1) = [ind_min_s,ind_max_s]
      ind_minmax(1:2,2) = [ind_min_m,ind_max_m]
      ind_minmax(1:2,3) = [ind_min_l,ind_max_l]
      !write(1000+iproc,*) 'calling init_matrix_taskgroups_wrapper'
      call init_matrix_taskgroups_wrapper(iproc, nproc, bigdft_mpi%mpi_comm, in%enable_matrix_taskgroups, &
           3, linmat%smat, ind_minmax)
      !write(1000+iproc,*) 'after init_matrix_taskgroups_wrapper'
      call f_free(ind_minmax)

      call f_release_routine()

    end subroutine init_bigdft_matrices

    


    subroutine get_minmax_lines(collcom, sparsemat, coarse_fine, imin_old, imax_old, imin_new, imax_new)
      use futile
      use sparsematrix_base, only: sparse_matrix
      use module_types, only: comms_linear
      implicit none

      ! Calling arguments
      type(comms_linear),intent(in) :: collcom
      type(sparse_matrix),intent(in) :: sparsemat
      character(len=1),intent(in) :: coarse_fine
      integer,dimension(sparsemat%nfvctr),intent(inout) :: imin_old, imax_old, imin_new, imax_new

      ! Local variables
      integer :: ipt, ii, i0, iiorb, i, i0i, j, i0j, jj, nmiddle

      nmiddle = sparsemat%nfvctr/2 + 1

      select case (coarse_fine)
      case ('c')
          do ipt=1,collcom%nptsp_c
              ii = collcom%norb_per_gridpoint_c(ipt)
              i0 = collcom%isptsp_c(ipt)
              do i=1,ii
                  i0i = i0 + i
                  iiorb = collcom%indexrecvorbital_c(i0i)
                  do j=1,ii
                      i0j = i0 + j
                      !ii = mod(collcom%indexrecvorbital_c(i0i)-1,sparsemat%nfvctr)+1
                      !write(*,*) 'i, j, i0i, i0j', i, j, i0i, i0j
                      jj = collcom%indexrecvorbital_c(i0j)
                      imin_old(iiorb) = min(imin_old(iiorb),jj)
                      imax_old(iiorb) = max(imax_old(iiorb),jj)
                      if (jj>nmiddle) then
                          imin_new(iiorb) = min(imin_new(iiorb),jj)
                      else
                          imax_new(iiorb) = max(imax_new(iiorb),jj+sparsemat%nfvctr)
                      end if
                  end do
              end do
          end do
      case ('f')
          do ipt=1,collcom%nptsp_f
              ii = collcom%norb_per_gridpoint_f(ipt)
              i0 = collcom%isptsp_f(ipt)
              do i=1,ii
                  i0i = i0 + i
                  iiorb = collcom%indexrecvorbital_f(i0i)
                  do j=1,ii
                      i0j = i0 + j
                      !ii = mod(collcom%indexrecvorbital_f(i0i)-1,sparsemat%nfvctr)+1
                      jj = collcom%indexrecvorbital_f(i0j)
                      imin_old(iiorb) = min(imin_old(iiorb),jj)
                      imax_old(iiorb) = max(imax_old(iiorb),jj)
                      if (jj>nmiddle) then
                          imin_new(iiorb) = min(imin_new(iiorb),jj)
                      else
                          imax_new(iiorb) = max(imax_new(iiorb),jj+sparsemat%nfvctr)
                      end if
                  end do
              end do
          end do
      case default
          call f_err_throw("Wrong value for coarse_fine, must be 'c' or 'f'")
      end select

    end subroutine get_minmax_lines



    subroutine check_orbital_matrix_distribution(orbs, smat, ind_min, ind_max)
      use futile
      use module_types, only: orbitals_data
      use sparsematrix_base, only: sparse_matrix
      implicit none

      ! Calling arguments:
      type(orbitals_data),intent(in) :: orbs
      type(sparse_matrix),intent(in) :: smat
      integer,intent(inout) :: ind_min, ind_max

      ! Local variables
      integer :: iorb, ispin, iiorb, ishift, isegstart, isegend, iseg, i, ii

      call f_routine(id='check_orbital_matrix_distribution')

      do iorb=orbs%isorb+1,orbs%isorb+orbs%norbp
          if (orbs%spinsgn(iorb)>0.d0) then
              ispin=1
          else
              ispin=2
          end if
          iiorb = mod(iorb-1,smat%nfvctr)+1 ! spin-independent index
          !ishift=(ispin-1)*smat%nvctr!-smat%isvctrp_tg
          isegstart = smat%istsegline(iiorb)
          isegend = smat%istsegline(iiorb) + smat%nsegline(iiorb) - 1
          do iseg=isegstart,isegend
              ii=smat%keyv(iseg)
              do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                  !ind_min = min(ii+ishift,ind_min)
                  !ind_max = max(ii+ishift,ind_max)
                  ind_min = min(ii,ind_min)
                  ind_max = max(ii,ind_max)
                  ii=ii+1
              end do
          end do
      end do

      !write(*,*) 'end sub: ind_min, ind_max', ind_min, ind_max

      call f_release_routine()

    end subroutine check_orbital_matrix_distribution

end module bigdft_matrices
