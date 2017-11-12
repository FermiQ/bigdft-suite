!> @file
!!  Module to handle the rotation and translation of a scalar field
!!  according to a transformation
!!
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module handling the rotation and trasnlation of a scalar field as wavefunctions
module reformatting
  use module_defs, only: gp
  use f_enums
  implicit none
  
  private

  integer, parameter :: reformatting_reasons=7
  character(len=*), parameter :: HGRIDS=           'Grid spacing has changed'
  character(len=*), parameter :: UNIFORM=          'Box size has changed'
  character(len=*), parameter :: COMPRESSED_COARSE='Number of coarse grid points has changed'
  character(len=*), parameter :: COMPRESSED_FINE=  'Number of fine grid points has changed'
  character(len=*), parameter :: SHIFT=            'Molecule was shifted'
  character(len=*), parameter :: ROTATE=           'Molecule was rotated'
  character(len=*), parameter :: NO_REFORMAT=      'No reformatting required'
  character(len=*), parameter :: WRAP=             'Wrapping/unwrapping'

  type(f_enumerator), public :: REFORMAT_COPY=f_enumerator('REFORMAT_COPY',-123,null())
  type(f_enumerator), public :: REFORMAT_FULL=f_enumerator('REFORMAT_FULL',-122,null())
  type(f_enumerator), public :: REFORMAT_WRAP=f_enumerator('REFORMAT_WRAP',-121,null())

  
  public :: inspect_rototranslation,reformat_one_supportfunction,print_reformat_summary
  public :: reformatting_init_info

contains


  subroutine reformatting_init_info(info)
    use dictionaries, dict_set => set
    implicit none
    type(dictionary), pointer :: info
    call dict_init(info)

    call dict_set(info//NO_REFORMAT,0)
    call dict_set(info//HGRIDS,0)
    call dict_set(info//UNIFORM,0)
    call dict_set(info//COMPRESSED_COARSE,0)
    call dict_set(info//COMPRESSED_FINE,0)
    call dict_set(info//SHIFT,0)
    call dict_set(info//ROTATE,0)
    call dict_set(info//WRAP,0)
    
  end subroutine reformatting_init_info
  
  !> Print information about the reformatting due to restart
  subroutine print_reformat_summary(info)
    use module_base
    use yaml_output
    implicit none
    type(dictionary), pointer :: info
    !local variables
    integer :: i
    type(dictionary), pointer :: iter
    integer, dimension(:), allocatable :: reformat_reason ! array giving reasons for reformatting

    reformat_reason=f_malloc(0.to.dict_size(info)-1,id='reformat_reason')
    
    nullify(iter)
    i=0
    do while(iterating(iter,on=info))
       reformat_reason(i)=iter
       i=i+1
    end do
    
    if (mpisize(bigdft_mpi%mpi_comm) > 1) &
         call fmpi_allreduce(reformat_reason, FMPI_SUM, comm=bigdft_mpi%mpi_comm)

    if (mpirank(bigdft_mpi%mpi_comm)==0) then
       call yaml_mapping_open('Overview of the reformatting (several categories may apply)')
       nullify(iter)
       i=0
       do while(iterating(iter,on=info))
          call yaml_map(trim(dict_key(iter)),reformat_reason(i))
          i=i+1
       end do
!!$
!!$       call yaml_map('No reformatting required', reformat_reason(0))
!!$       call yaml_map('Grid spacing has changed', reformat_reason(1))
!!$       call yaml_map('Number of coarse grid points has changed', reformat_reason(2))
!!$       call yaml_map('Number of fine grid points has changed', reformat_reason(3))
!!$       call yaml_map('Box size has changed', reformat_reason(4))
!!$       call yaml_map('Molecule was shifted', reformat_reason(5))
!!$       call yaml_map('Molecule was rotated', reformat_reason(6))
!!$       call yaml_map('Wrapping/unwrapping', reformat_reason(7))
       call yaml_mapping_close()
    end if

    call f_free(reformat_reason)
    call dict_free(info)
  end subroutine print_reformat_summary

  subroutine rototranslations_shifts(rt,llr_src,llr_dest,&
       centre_src,centre_dest,da)
    use rototranslations
    use locregs
    implicit none
    type(rototranslation), intent(in) :: rt
    type(locreg_descriptors), intent(in) :: llr_src,llr_dest
    real(gp), dimension(3), intent(out) :: centre_src,centre_dest,da
    !local variables
    real(gp), dimension(3) :: oxyz_src,oxyz_dest

    call calculate_origins(llr_dest,llr_src,oxyz_src,oxyz_dest)
    call calculate_shifts(rt,llr_src%mesh_coarse,&
         llr_dest%mesh_coarse,oxyz_src,oxyz_dest,&
         centre_src,centre_dest,da)

  end subroutine rototranslations_shifts
  
  !> Checks whether reformatting is needed based on various criteria and returns final shift and centres needed for reformat
  function inspect_rototranslation(frag_trans,tol,dest_llr,src_llr,dest_glr,src_glr,info) result(en)
    use rototranslations
    use locregs
    use dictionaries, dict_set=>set
    use box
    implicit none
    real(gp), intent(in) :: tol ! tolerance for rotations and shifts
    type(locreg_descriptors), intent(in) :: dest_llr, src_llr,dest_glr,src_glr
    type(dictionary), pointer :: info
    type(rototranslation), intent(in) :: frag_trans ! includes centres of rotation in global coordinates, shift and angle
    type(f_enumerator) :: en
    ! local variables
    logical :: reformat_needed ! logical telling whether reformat is needed
    logical :: wrap_around ! periodic case where no reformatting is needed but tmb needs 'unwrapping'
    logical, dimension(3) :: per

    real(gp) :: displ !, mindist
    integer :: i
    real(gp), dimension(3) :: centre_src,centre_dest,da

    call rototranslations_shifts(frag_trans,src_llr,dest_llr,&
         centre_src,centre_dest,da)
!!$    call calculate_origins(dest_llr,src_llr,oxyz_src,oxyz_dest)
!!$    call calculate_shifts(frag_trans,src_llr%mesh_coarse,dest_llr%mesh_coarse,oxyz_src,oxyz_dest,&
!!$         centre_src,centre_dest,da)
    displ=square_gd(dest_llr%mesh_coarse,da)

    reformat_needed=.false.
    if (any(dest_llr%mesh_coarse%hgrids /= src_llr%mesh_coarse%hgrids)) then
       reformat_needed=.true.
       call increment(info//HGRIDS)
    end if
    if (any(dest_llr%mesh_coarse%ndims /= src_llr%mesh_coarse%ndims)) then
       reformat_needed=.true.
       call increment(info//UNIFORM)
    end if
    if (dest_llr%wfd%nvctr_c /= src_llr%wfd%nvctr_c) then
       reformat_needed=.true.
       call increment(info//COMPRESSED_COARSE)
    end if
    if (dest_llr%wfd%nvctr_f /= src_llr%wfd%nvctr_f) then
       reformat_needed=.true.
       call increment(info//COMPRESSED_FINE)
    end if
    if (abs(displ) > tol) then
       reformat_needed=.true.
       call increment(info//SHIFT)
    end if
    if (abs(frag_trans%theta) > tol) then
       reformat_needed=.true.
       call increment(info//ROTATE)
    end if
    if (.not. reformat_needed) then
       call increment(info//NO_REFORMAT)
    end if
    ! check to make sure we don't need to 'unwrap' (or wrap) tmb in periodic case
    wrap_around=.false.

    do i=1,3
       if (dimension_outside(i,dest_llr,dest_glr) .or. dimension_outside(i,src_llr,src_glr)) then
          if (different_starts(i,dest_llr,src_llr)) wrap_around = .true.
       end if
    end do
    if (wrap_around) then
       call increment(info//WRAP)
    end if

    !result
    if (reformat_needed) then
       en=REFORMAT_FULL
    else if (wrap_around) then
       en=REFORMAT_WRAP
    else
       en=REFORMAT_COPY
    end if

  contains

    subroutine increment(dict)
      implicit none
      type(dictionary), pointer :: dict
      !local variables
      integer :: idum
      idum=dict
      call dict_set(dict,idum+1)
    end subroutine increment

    pure function different_starts(idim,llr,llr_old) result(yes)
      implicit none
      integer, intent(in) :: idim
      type(locreg_descriptors), intent(in) :: llr,llr_old
      logical :: yes
      yes=.false.
      select case(idim)
      case(1)
         yes= llr%ns1/=llr_old%ns1
      case(2)
         yes= llr%ns2/=llr_old%ns2
      case(3)
         yes= llr%ns3/=llr_old%ns3
      end select

    end function different_starts

    
    pure function dimension_outside(idim,llr,glr) result(yes)
      implicit none
      integer, intent(in) :: idim
      type(locreg_descriptors), intent(in) :: llr,glr
      logical :: yes
      !local variables
      logical, dimension(3) :: per
      per=cell_periodic_dims(glr%mesh_coarse)
      yes=.false.
      select case(idim)
      case(1)
         yes=tmb_wrap(per(1),llr%ns1,llr%mesh_coarse%ndims(1),glr%mesh_coarse%ndims(1))
      case(2)
         yes=tmb_wrap(per(2),llr%ns2,llr%mesh_coarse%ndims(2),glr%mesh_coarse%ndims(2))
      case(3)
         yes=tmb_wrap(per(3),llr%ns3,llr%mesh_coarse%ndims(3),glr%mesh_coarse%ndims(3))
      end select
    end function dimension_outside
    
    !> Checks to see if a tmb wraps around the the unit cell
    !! knowing that there could have been a modulo operation
    pure function tmb_wrap(periodic,ns,n,nglr)
      use module_base
      implicit none
      logical, intent(in) :: periodic
      integer, intent(in) :: ns, n, nglr
      logical :: tmb_wrap

      tmb_wrap = .false.
      if (periodic) then
         !<=?
         if (ns<0 .or. ns+n>nglr) then
            tmb_wrap = .true.
         end if
      end if

    end function tmb_wrap
  end function inspect_rototranslation
  
  subroutine calculate_origins(llr,llr_old,oxyz_src,oxyz_dest)
    use locregs
    use bounds, only: locreg_mesh_coarse_origin,locreg_mesh_origin
    implicit none
    type(locreg_descriptors), intent(in) :: llr, llr_old
    real(gp), dimension(3), intent(out) :: oxyz_dest,oxyz_src
    !local variables
    integer, dimension(3) :: ns,ns_old
        
    ns_old(1)=llr_old%ns1
    ns_old(2)=llr_old%ns2
    ns_old(3)=llr_old%ns3
    ns(1)=llr%ns1
    ns(2)=llr%ns2
    ns(3)=llr%ns3

    oxyz_src=true_origin(llr_old)
    oxyz_dest=true_origin(llr)
!!$    oxyz_src=llr_old%mesh_coarse%hgrids*ns_old-0.5_gp*locreg_mesh_coarse_origin(llr_old%mesh_coarse)
!!$    oxyz_dest=*llr%mesh_coarse%hgrids*ns-0.5_gp*locreg_mesh_coarse_origin(llr%mesh_coarse)

  end subroutine calculate_origins

  !> this might be moved in locregs methods
  function true_origin(lr) result(oxyz)
    use locregs, only: locreg_descriptors
    use bounds, only: locreg_mesh_origin
    use box, only: cell_r
    implicit none
    type(locreg_descriptors), intent(in) :: lr
    real(gp), dimension(3) :: oxyz

    oxyz=-locreg_mesh_origin(lr%mesh)
    oxyz(1)=oxyz(1)+cell_r(lr%mesh,lr%nsi1+1,dim=1)
    oxyz(2)=oxyz(2)+cell_r(lr%mesh,lr%nsi2+1,dim=2)
    oxyz(3)=oxyz(3)+cell_r(lr%mesh,lr%nsi3+1,dim=3)

  end function true_origin

  !> Make frag_trans the argument so can eliminate need for interface
  subroutine reformat_one_supportfunction(llr,llr_old,&!geocode,hgrids_old,&
       n_old,psigold,&
       !hgrids,&
       n,&
       !centre_old,centre_new,da,&
       frag_trans,psi,psirold,tag)
    use module_base
    use locregs
    use rototranslations
    !use yaml_output
    use locreg_operations
    use bounds, only: ext_buffers_coarse,locreg_mesh_coarse_origin
    use box
    implicit none
    integer, dimension(3), intent(in) :: n,n_old
!    real(gp), dimension(3), intent(in) :: hgrids,hgrids_old
    !type(wavefunctions_descriptors), intent(in) :: wfd
    type(locreg_descriptors), intent(in) :: llr, llr_old
!    real(gp), dimension(3), intent(in) :: centre_old,centre_new,da
    type(rototranslation), intent(in) :: frag_trans
    real(wp), dimension(0:n_old(1),2,0:n_old(2),2,0:n_old(3),2), intent(in) :: psigold
    real(wp), dimension(llr%wfd%nvctr_c+7*llr%wfd%nvctr_f), intent(out) :: psi
    real(wp), dimension(llr_old%d%n1i,llr_old%d%n2i,llr_old%d%n3i), optional, intent(in) :: psirold
    integer, optional, intent(in) :: tag ! filename for printing functions, used for debugging only

    !local variables
    character(len=*), parameter :: subname='reformatonesupportfunction'
    character(len=1) :: geocode !< @copydoc poisson_solver::doc::geocode
    logical, dimension(3) :: per
    integer, dimension(3) :: nb
!!$  integer, dimension(3) :: ndims_tmp
    real(wp), external :: dnrm2
    real(gp), dimension(3) :: centre_src,centre_dest,da
    real(gp), dimension(3) :: oxyz_src,oxyz_dest
    real(wp), dimension(:), allocatable :: ww,wwold
    real(wp), dimension(:,:,:,:,:,:), allocatable :: psig
    real(wp), dimension(:,:,:), allocatable :: psir
    real(wp), dimension(:,:,:), pointer :: psifscfold, psifscf !no reason for pointers
    integer :: i,j,k
    ! isf version
    type(workarr_sumrho) :: w

    call f_routine(id=subname)

    geocode=cell_geocode(llr%mesh)
    
    ! old reformatting, otherwise in ISF (routines below should be incapsulated in a s1s0 method (libconv)
    if (.not. present(psirold)) then
       !conditions for periodicity in the three directions
       per(1)=(geocode /= 'F')
       per(2)=(geocode == 'P')
       per(3)=(geocode /= 'F')

       !buffers related to periodicity
       !WARNING: the boundary conditions are not assumed to change between new and old
       call ext_buffers_coarse(per(1),nb(1))
       call ext_buffers_coarse(per(2),nb(2))
       call ext_buffers_coarse(per(3),nb(3))

       psifscf = f_malloc_ptr(-nb.to.2*n+1+nb,id='psifscf')
       psifscfold = f_malloc_ptr(-nb.to.2*n_old+1+nb,id='psifscfold')

       wwold = f_malloc((2*n_old(1)+2+2*nb(1))*(2*n_old(2)+2+2*nb(2))*(2*n_old(3)+2+2*nb(3)),id='wwold')

       if (geocode=='F') then
          call synthese_grow(n_old(1),n_old(2),n_old(3),wwold,psigold,psifscfold)
       else if (geocode=='S') then
          call synthese_slab(n_old(1),n_old(2),n_old(3),wwold,psigold,psifscfold)
       else if (geocode=='P') then
          call synthese_per(n_old(1),n_old(2),n_old(3),wwold,psigold,psifscfold)
       end if

       call f_free(wwold)
    end if

    !if (present(tag)) then
    !   open(tag+10000)
    !   do i=-nb(1),2*n_old(1)+1+nb(1)
    !   do j=-nb(2),2*n_old(2)+1+nb(2)
    !   do k=-nb(3),2*n_old(3)+1+nb(3)
    !      write(tag+10000,'(3(I6,1x),1x,2(F12.6,1x))') i,j,k,psifscfold(i,j,k),&
    !           dnrm2((2*n_old(1)+2+2*nb(1))*(2*n_old(2)+2+2*nb(1))*(2*n_old(3)+2+2*nb(1)),psifscfold,1)
    !   end do
    !   end do
    !   end do
    !   close(tag+10000)
    !end if


!!$  !print the suggested order
!!$  call yaml_map('Suggested order for the transformation',irp)

!!$    if (.not. present(psirold)) then
!!$       ndims_old=(2*n_old+2+2*nb)
!!$       ndims_new=(2*n+2+2*nb)
!!$    else
!!$       ndims_old=[llr_old%d%n1i,llr_old%d%n2i,llr_old%d%n3i]
!!$       ndims_new=[llr%d%n1i,llr%d%n2i,llr%d%n3i]
!!$    end if

    !call calculate_origins(llr,llr_old,oxyz_src,oxyz_dest)
    call rototranslations_shifts(frag_trans,llr_old,llr,&
         centre_src,centre_dest,da)
    if (.not. present(psirold)) then
       call apply_rototranslation(frag_trans,.true.,&
            llr_old%mesh_fine,llr%mesh_fine,& 
            centre_src,centre_dest,da,&
            psifscfold,psifscf)
!!$       call field_rototranslation3D(nd+1,nrange,y_phi,frag_trans%Rmat,da,&
!!$            centre_old,centre_new,irp,&
!!$            hgridsh_old,ndims_old,psifscfold,&
!!$            hgridsh,ndims_new,psifscf)
    else
       psir=f_malloc((/llr%d%n1i,llr%d%n2i,llr%d%n3i/),id='psir')
       call apply_rototranslation(frag_trans,.false.,&
            llr_old%mesh,llr%mesh,&
            centre_src,centre_dest,da,&
            psirold,psir)
!!$       call field_rototranslation3D(nd+1,nrange,y_phi,frag_trans%Rmat,da,&
!!$            centre_old,centre_new,irp,&
!!$            hgridsh_old,ndims_old,psirold,&
!!$            hgridsh,ndims_new,psir)
    end if

    
    if (present(tag)) then
       open(tag+20000)
       do i=-nb(1),2*n(1)+1+nb(1)
          do j=-nb(2),2*n(2)+1+nb(2)
             do k=-nb(3),2*n(3)+1+nb(3)
                write(tag+20000,'(3(I6,1x),1x,2(F12.6,1x))') i,j,k,psifscf(i,j,k),&
                     dnrm2((2*n(1)+2+2*nb(1))*(2*n(2)+2+2*nb(1))*(2*n(3)+2+2*nb(1)),psifscf,1)
             end do
          end do
       end do
       close(tag+20000)
    end if


    !  print*, 'norm of psifscf ',dnrm2((2*n(1)+2+2*nb(1))*(2*n(2)+2+2*nb(1))*(2*n(3)+2+2*nb(1)),psifscf,1)
    if (.not. present(psirold)) then
       call f_free_ptr(psifscfold)
       psig = f_malloc((/ 0.to.n(1), 1.to.2, 0.to.n(2), 1.to.2, 0.to.n(3), 1.to.2 /),id='psig')
       ww = f_malloc((2*n(1)+2+2*nb(1))*(2*n(2)+2+2*nb(2))*(2*n(3)+2+2*nb(3)),id='ww')

       if (geocode=='F') then
          call analyse_shrink(n(1),n(2),n(3),ww,psifscf,psig)
       else if (geocode == 'S') then
          call analyse_slab(n(1),n(2),n(3),ww,psifscf,psig)
       else if (geocode == 'P') then
          call analyse_per(n(1),n(2),n(3),ww,psifscf,psig)
       end if

       call f_free_ptr(psifscf)


       if (present(tag)) then
          open(tag)
          do i=1,n(1)
             do j=1,n(2)
                do k=1,n(3)
                   write(tag,'(3(I6,1x),1x,2(F12.6,1x))') i,j,k,psig(i,1,j,1,k,1),&
                        dnrm2(8*(n(1)+1)*(n(2)+1)*(n(3)+1),psig,1)
                end do
             end do
          end do
          close(tag)
       end if


!!$    print*, 'norm new psig ',dnrm2(8*(n(1)+1)*(n(2)+1)*(n(3)+1),psig,1),n(1),n(2),n(3)
       call compress_plain(n(1),n(2),0,n(1),0,n(2),0,n(3),  &
            llr%wfd%nseg_c,llr%wfd%nvctr_c,llr%wfd%keygloc(1,1),llr%wfd%keyvloc(1),   &
            llr%wfd%nseg_f,llr%wfd%nvctr_f,&
            llr%wfd%keygloc(1,llr%wfd%nseg_c+min(1,llr%wfd%nseg_f)),&
            llr%wfd%keyvloc(llr%wfd%nseg_c+min(1,llr%wfd%nseg_f)),   &
            psig,psi(1),psi(llr%wfd%nvctr_c+min(1,llr%wfd%nvctr_f)))

       call f_free(psig)
       call f_free(ww)
    else
       call initialize_work_arrays_sumrho(llr,.true.,w)
       call f_zero(psi)
!!$     write(*,*) 'iproc,norm psirnew ',dnrm2(llr%d%n1i*llr%d%n2i*llr%d%n3i,psir,1),llr%d%n1i,llr%d%n2i,llr%d%n3i
       call isf_to_daub(llr,w,psir,psi)
       call deallocate_work_arrays_sumrho(w)
       call f_free(psir)
    end if

!!$  print*, 'norm of reformatted psi ',dnrm2(llr%wfd%nvctr_c+7*llr%wfd%nvctr_f,psi,1),llr%wfd%nvctr_c,llr%wfd%nvctr_f
!!$  print*, 'norm of reformatted psic ',dnrm2(llr%wfd%nvctr_c,psi,1)
!!$  print*, 'norm of reformatted psif ',dnrm2(llr%wfd%nvctr_f*7,psi(llr%wfd%nvctr_c+min(1,llr%wfd%nvctr_f)),1)

    call f_release_routine()

  contains

    !>determinant of a 3x3 matrix
    pure function det_33(a) result(det)
      implicit none
      real(gp), dimension(3,3), intent(in) :: a
      real(gp) :: det

      det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
           + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3))  &
           + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
    end function det_33

  END SUBROUTINE reformat_one_supportfunction

  subroutine calculate_shifts(rt,mesh_src,mesh_dest,oxyz_src,oxyz_dest,&
       centre_src,centre_dest,da)
    use box
    use rototranslations
    implicit none
    type(rototranslation), intent(in) :: rt
    type(cell), intent(in) :: mesh_src,mesh_dest
    real(gp), dimension(3), intent(in) :: oxyz_src,oxyz_dest
    real(gp), dimension(3), intent(out) :: centre_src,centre_dest,da

        !here for periodic BC we should check the correctness of the centre definition
    !in particular it appears weird that mesh_dest is always used for closest r
    !as well as the shift of hgrid/2 that have to be imposed for the detection of da array
    centre_src=closest_r(mesh_dest,rt%rot_center_src,oxyz_src)
    centre_dest=closest_r(mesh_dest,rt%rot_center_dest,oxyz_dest)

    da = centre_dest-centre_src-(mesh_src%hgrids-mesh_dest%hgrids)*0.5_gp
    
  end subroutine calculate_shifts
  
  subroutine apply_rototranslation(rt,hr,mesh_src,mesh_dest,&
       centre_src,centre_dest,da,psi_src,psi_dest)
    use module_defs
    use bounds, only: locreg_mesh_shape
    use box
    use rototranslations
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_strings
    implicit none
    logical, intent(in) :: hr !< boolean for high resolution approach (or the opposite?)
    type(rototranslation), intent(in) :: rt
    type(cell), intent(in) :: mesh_src,mesh_dest
    !> origin of the coordinate system in the two reference frames
    real(gp), dimension(3), intent(in) :: centre_src,centre_dest,da
    real(wp), dimension(mesh_src%ndim), intent(in) :: psi_src
    real(wp), dimension(mesh_dest%ndim), intent(out) :: psi_dest
    !local variables
!!$    integer, dimension(3) :: ndims_old,ndims_new
!!$    real(gp), dimension(3) :: hgridsh,hgridsh_old,centre_src,centre_dest,da
    real(wp), dimension(:), allocatable :: x_phi
    real(wp), dimension(:,:), allocatable :: y_phi
    integer :: itype, nd, nrange
!!$  real(gp), dimension(3) :: rrow
    !  real(gp), dimension(3,3) :: rmat !< rotation matrix
    !  real(gp) :: sint,cost,onemc,ux,uy,uz
    integer, dimension(3) :: irp

!!$    ndims_old=locreg_mesh_shape(mesh_src,hr)
!!$    ndims_new=locreg_mesh_shape(mesh_dest,hr)
!!$    call calculate_shifts(rt,mesh_src,mesh_dest,oxyz_src,oxyz_dest,&
!!$         centre_src,centre_dest,da)

!!$    ! transform to new structure
!!$    hgridsh=mesh_dest%hgrids
!!$    hgridsh_old=mesh_src%hgrids

    !create the scaling function array
    !use lots of points (to optimize one can determine how many points are needed at max)
    itype=16
    nd=2**20

    x_phi = f_malloc(0.to.nd,id='x_phi')
    y_phi = f_malloc((/0.to.nd,1.to.2/),id='y_phi')

    call my_scaling_function4b2B(itype,nd,nrange,x_phi,y_phi)
    !such check is rather a debug check, it might be removed
    if( abs(y_phi(nd/2,1)-1)>1.0e-10 ) then
       stop " wrong scaling function 4b2B: not a centered one "
    endif

    call f_free(x_phi)

    !call field_rototranslation(nd,nrange,y_phi,da,rt%rot_axis,centre_old,centre_new,rt%theta,&
    !     hgridsh_old,ndims_tmp,psifscf_tmp,hgridsh,(2*n+2+2*nb),psifscf)

!!$  sint=sin(rt%theta)
!!$  cost=cos(rt%theta)
!!$  onemc=1.0_gp-cost
!!$  ux=rt%rot_axis(1)
!!$  uy=rt%rot_axis(2)
!!$  uz=rt%rot_axis(3)

!!$  call yaml_sequence_open('Rotation matrix elements')
!!$  call yaml_sequence(trim(yaml_toa((/&
!!$       cost + onemc*ux**2   , ux*uy*onemc - uz*sint, ux*uz*onemc + uy*sint /),fmt='(1pg20.12)')))
!!$  call yaml_sequence(trim(yaml_toa((/&
!!$       ux*uy*onemc +uz*sint , cost + onemc*uy**2   , uy*uz*onemc - ux*sint /),fmt='(1pg20.12)')))
!!$  call yaml_sequence(trim(yaml_toa((/&
!!$       ux*uz*onemc -uy*sint , uy*uz*onemc + ux*sint, cost + onemc*uz**2    /),fmt='(1pg20.12)')))
!!$  call yaml_sequence_close()


!!$  !identify the rotation matrix elements
!!$  rmat=reshape([&
!!$       cost+onemc*ux**2    , ux*uy*onemc-uz*sint , ux*uz*onemc+uy*sint, &    !first row (xp)
!!$       ux*uy*onemc+uz*sint , cost+onemc*uy**2    , uy*uz*onemc-ux*sint, &   !second row (yp)
!!$       ux*uz*onemc-uy*sint , uy*uz*onemc+ux*sint , cost+onemc*uz**2], [3,3]) !third row (zp)

!!$  rmat(:,1)=(/cost+onemc*ux**2  ,ux*uy*onemc-uz*sint ,ux*uz*onemc+uy*sint/)
!!$  !second row (yp)
!!$  rmat(:,2)=(/ux*uy*onemc+uz*sint,cost+onemc*uy**2   ,uy*uz*onemc-ux*sint/)
!!$  !third row (zp)
!!$  rmat(:,3)=(/ux*uz*onemc-uy*sint,uy*uz*onemc+ux*sint,cost+onemc*uz**2   /)

!!$  !!write some output on the screen
!!$  !!print matrix elements, to be moved at the moment of identification of the transformation
!!$  call yaml_map('Rotation axis',rt%rot_axis,fmt='(1pg20.12)')
!!$  call yaml_map('Rotation angle (deg)',rt%theta*180.0_gp/pi_param,fmt='(1pg20.12)')
!!$  call yaml_map('Translation vector',da,fmt='(1pg20.12)')
!!$  call yaml_map('Rotation matrix ',rt%Rmat,fmt='(1pg20.12)')
!!$  call yaml_map('Rotation matrix',rmat,fmt='(1pg20.12)')
!!$  call yaml_map('Determinants',[det_33(rt%Rmat),det_33(rmat)])


!!$  !pay attention to what happens if two values are identical
!!$  !from where xp should be determined
!!$  rrow=abs(rmat(:,1))
!!$  irp(1)=maxloc(rrow,1)
!!$  !form where zp should be determined (note that the third line has been used)
!!$  rrow=abs(rmat(:,3))
!!$  !exclude of course the previously found direction
!!$  rrow(irp(1))=0.0_gp
!!$  irp(3)=maxloc(rrow,1)
!!$  !then the last dimension, which is determined by exclusion
!!$  rrow=1.0_gp
!!$  rrow(irp(1))=0.d0
!!$  rrow(irp(3))=0.d0
!!$  irp(2)=maxloc(rrow,1)

    !try different solutions, one of these should always work
    irp=selection(rt%Rmat)
    !otherwise we have a problem
    if (repeated(abs(irp))) then
       call f_err_throw('Determination of the best array failed, irp='//&
            trim(yaml_toa(irp,fmt='(i5)')),err_name='BIGDFT_RUNTIME_ERROR')
       return
    end if

    if (.not. hr) irp(:)=abs(irp)
    call field_rototranslation3D(nd+1,nrange,y_phi,rt%Rmat,da,&
         centre_src,centre_dest,irp,&
         mesh_src%hgrids,mesh_src%ndims,psi_src,&
         mesh_dest%hgrids,mesh_dest%ndims,psi_dest)

    call f_free(y_phi)

  end subroutine apply_rototranslation

  
  !> Select the best possible rotation sequence by
  !! considering the values of the coefficients of the
  !! rotation matrix
  pure function selection(rmat) result(irp)

    implicit none
    real(gp), dimension(3,3), intent(in) :: rmat !< rotation matrix
    integer, dimension(3) :: irp
    !local variables
    integer :: i
    integer, dimension(3) :: ib1,ib3
!!$      integer :: isgn
!!$      real(gp), dimension(3) :: rrow

    !determine ideal sequence for rotation, for important rows
    ib1=reorder(rmat(:,1),1)
    ib3=reorder(rmat(:,3),3)

    !verify if either one or three have multiple choices
    if (equabs(rmat(ib1(1),1),rmat(ib1(2),1)) .and. .not. equabs(rmat(ib3(1),3),rmat(ib3(2),3))) then
       !only ib1 has multiple choices, pick the one which is closest to cyclic permutation (if present)
       if (modulo(ib3(1),3) + 1 == ib1(2) .or. ib1(1)==ib3(1)) then
          !swap
          i=ib1(1)
          ib1(1)=ib1(2)
          ib1(2)=i
       end if
    else if (.not. equabs(rmat(ib1(1),1),rmat(ib1(2),1)) .and. equabs(rmat(ib3(1),3),rmat(ib3(2),3))) then
       !only ib3 has multiple choices
       if (modulo(ib3(2),3) + 1 == ib1(1) .or. ib3(1)==ib1(1)) then
          !swap
          i=ib3(1)
          ib3(1)=ib3(2)
          ib3(2)=i
       end if
    else if (equabs(rmat(ib1(1),1),rmat(ib1(2),1)) .and. equabs(rmat(ib3(1),3),rmat(ib3(2),3))) then
       !both of the row has multiple choices, therefore at least cyclic permutation must be present.
       !both of them are cyclic, choose the last one
       if (modulo(ib3(2),3) + 1 == ib1(1)) then
          !swap
          i=ib3(1)
          ib3(1)=ib3(2)
          ib3(2)=i
       else if (modulo(ib3(1),3) + 1 == ib1(2)) then
          !swap
          i=ib1(1)
          ib1(1)=ib1(2)
          ib1(2)=i
       else if (ib3(1) == ib1(1)) then
          !otherwise just ensure that the two are not equal
          !swap
          i=ib3(1)
          ib3(1)=ib3(2)
          ib3(2)=i
       end if
    else if (ib3(1) == ib1(1)) then
       !swap ib1,instead of ib3
       i=ib1(1)
       ib1(1)=ib1(2)
       ib1(2)=i
    end if
    !then assign the rotations
    irp(1)=ib1(1)
    irp(3)=ib3(1)

    !define the best for the second
    ib1=1
    ib1(irp(1))=0
    ib1(irp(3))=0
    irp(2)=maxloc(ib1,1)

!!$      irp(1)=ibest(1,1)
!!$      irp(3)=ibest(1,3)
!!$      if (ibest(1,3)==irp(1)) then
!!$         if (abs(abs(rmat(ibest(1,1),1))-abs(rmat(ibest(2,1),1)))< 1.d-12) then
!!$            irp(1)=ibest(2,1)
!!$         else if (abs(abs(rmat(ibest(1,3),3))-abs(rmat(ibest(2,3),3)))< 1.d-12) then
!!$            irp(3)=ibest(2,3)
!!$         else !better to preserve the first choice for the last
!!$            irp(1)=ibest(2,1)
!!$         end if
!!$      end if
!!$
!!$      !define the best for the second
!!$      ibest(:,2)=1
!!$      ibest(irp(1),2)=0
!!$      ibest(irp(3),2)=0
!!$      irp(2)=maxloc(ibest(:,2),1)

!!$      rrow=abs(rmat(:,1))
!!$      irp(1)=maxloc(rrow,1)
!!$      !form where zp should be determined (note that the third line has been used)
!!$      rrow=abs(rmat(:,3))
!!$      !exclude of course the previously found direction
!!$      rrow(irp(1))=0.0_gp
!!$      irp(3)=maxloc(rrow,1)
!!$      !then the last dimension, which is determined by exclusion
!!$      rrow=1.0_gp
!!$      rrow(irp(1))=0.d0
!!$      rrow(irp(3))=0.d0
!!$      irp(2)=maxloc(rrow,1)

!!$      !add to the transformations the sign of the axis of the chosen reference
!!$      !coordinate
!!$      !the second element has the sign which is the ratio of the previous two,
!!$      !plus a sign which is given by the fact that the order is a cyclic permutation
!!$      isgn=int(sign(1.0e0,real(rmat(irp(1),1)/rmat(irp(3),3))))
!!$      if (modulo(irp(1),3)+1 /= irp(2)) isgn=-isgn !cyclic permutation
!!$      irp(2)=isgn*irp(2)
!!$
!!$      !for the first and the third the sign is determined from the matrix element
!!$      isgn=int(sign(1.0e0,real(rmat(irp(1),1))))
!!$      irp(1)=isgn*irp(1)
!!$      isgn=int(sign(1.0e0,real(rmat(irp(3),3))))
!!$      irp(3)=isgn*irp(3)

  end function selection

  pure function repeated(ivec)
    implicit none
    integer, dimension(3), intent(in) :: ivec
    logical :: repeated
    repeated = ivec(1)==ivec(2) .or. ivec(2)==ivec(3) .or. ivec(1)==ivec(3)
  end function repeated

  !check if two objects are equal in absolute value modulo a given tolerance
  pure function equabs(a,b)
    implicit none
    real(gp), intent(in) :: a,b
    logical :: equabs
    real(gp), parameter :: tol=1.e-12_gp
    equabs=abs(abs(a)-abs(b)) < tol
  end function equabs

  !> defines the criterion for which one is better than two
  pure function better(idim,vec,one,two)
    implicit none
    integer, intent(in) :: idim,one,two
    real(gp), dimension(3), intent(in) :: vec
    logical :: better
    !local variables
    real(gp) :: vec1,vec2

    vec1=vec(one)
    vec2=vec(two)

    better=.false.

    !first criterion, most important: absolute value (clear separation)
    if (.not. equabs(vec1,vec2)) then
       better = abs(vec1)>abs(vec2)
    else
       !the two values are even. First choose the one which is positive
       if (sign(vec1,vec2) == vec1) then
          !the two objects have same sign and same absolute value
          if (one==idim .or. two==idim) then
             !better the one of the dimension
             better = one==idim
          else
             better = one<two .eqv. idim<=2
          end if
       else
          better = sign(1.0_gp,vec1)==1.0_gp
       end if
    end if

  end function better

  !> order the dimensions in terms of the maximum
  pure function reorder(vec,idim) result(imax)
    implicit none
    integer, intent(in) :: idim
    real(gp), dimension(3), intent(in) :: vec
    integer, dimension(3) :: imax
    !local variables
!!$      integer, dimension(3,3) :: ibest

    !initialization
    imax(1)=1
    imax(2)=2
    imax(3)=3
    if (better(idim,vec,2,1)) then
       if (better(idim,vec,3,1)) then
          if (better(idim,vec,2,3)) then
             !other worst case 2<3<1
             imax(1)=2
             imax(2)=3
             imax(3)=1
          else
             !  1>3<2, but 2<1 => 3<2<1
             imax(1)=3
             imax(3)=1
          end if
       else
          !2<1 and 3>1 => 2<1<3
          imax(1)=2
          imax(2)=1
       end if
    else
       if (better(idim,vec,3,2)) then
          if (better(idim,vec,3,1)) then
             !worst case, 3<1<2
             imax(1)=3
             imax(2)=1
             imax(3)=2
          else
             ! 1<3<2
             imax(2)=3
             imax(3)=2
          end if
       end if
    end if

    !once ordered preserve only the choices which are equal
!!$      if (abs(abs(vec(imax(2)))-abs(vec(imax(3)))) > 1.d-12 ) imax(3)=imax(2)
!!$      if (abs(abs(vec(imax(1)))-abs(vec(imax(2)))) > 1.d-12 ) imax(2:3)=imax(1)

  end function reorder

  
  !> Routine which directly applies the 3D transformation of the rototranslation
  !this routine has to be cleaned and the allocations of the work arrays have to be deplaced
  subroutine field_rototranslation3D(n_phi,nrange_phi,phi_ISF,Rmat,da,&
       centre_old,centre_new,&
       !newz,cost,sint,onemc,&
       iorder,hgrids_old,ndims_old,f_old,&
       hgrids_new,ndims_new,f_new)
    use module_base
!    use yaml_output
    implicit none
    integer, intent(in) :: n_phi,nrange_phi !< number of points of ISF array and real-space range
    integer, dimension(3), intent(in) :: iorder
    real(gp), dimension(3), intent(in) :: da !<coordinates of rigid shift vector
    real(gp), dimension(3), intent(in) :: centre_old,centre_new !<centre of rotation
    real(gp), dimension(3), intent(in) :: hgrids_old,hgrids_new !<dimension of old and new box
    integer, dimension(3), intent(in) :: ndims_old,ndims_new !<dimension of old and new box
    !real(gp), dimension(3) :: newz
    !real(gp), intent(in) ::  cost,sint,onemc
    real(gp), dimension(3,3), intent(in) :: Rmat !<rotation matrix
    real(gp), dimension(n_phi,2), intent(in) :: phi_ISF
    real(gp), dimension(ndims_old(1),ndims_old(2),ndims_old(3)), intent(in) :: f_old
    real(gp), dimension(ndims_new(1),ndims_new(2),ndims_new(3)), intent(out) :: f_new
    !local variables
    integer :: m_isf,k1,i,j,k,me,ms
    real(gp) :: dt,norm!,ux,uy,uz
    !$ real(gp) :: tt
    real(gp) :: scal !<scaling factor
    !real(gp), dimension(3,3) :: Rmat !<rotation matrix
    integer, dimension(3) :: isign,irp
    real(gp), dimension(:), allocatable :: shf
    real(gp), dimension(:), allocatable :: work,work2

    !print *,'3d'
    call f_routine(id='field_rototranslation3D')
    work =f_malloc(ndims_new(1)*(maxval(ndims_old))**2,id='work')
    work2=f_malloc(ndims_new(1)*ndims_new(2)*maxval(ndims_old),id='work2')

    m_isf=nrange_phi/2
    !shf=f_malloc(-m_isf .to. m_isf,id='shf')
    !for each of the dimensions build the interpolating vector which is needed

!!$    !identify the rotation matrix elements
!!$    ux=newz(1)
!!$    uy=newz(2)
!!$    uz=newz(3)
!!$    !first row (xp)
!!$    rmat(:,1)=(/cost+onemc*ux**2   ,ux*uy*onemc-uz*sint,ux*uz*onemc+uy*sint/)
!!$    !second row (yp)
!!$    rmat(:,2)=(/ux*uy*onemc+uz*sint,cost+onemc*uy**2   ,uy*uz*onemc-ux*sint/)
!!$    !third row (zp)
!!$    rmat(:,3)=(/ux*uz*onemc-uy*sint,uy*uz*onemc+ux*sint,cost+onemc*uz**2   /)


    !first step: determine xn from a coordinate n13o=xo or yo or zo
    !f_old (nxo,nyo,nzo) -> work(n11o,n12o,nxn) !n11o and n12o are the remaining dimensions
    !second step: determine yn from n22o=n11o or n12o
    !work(n11o,n12o,nxn) -> work2(n21o,nxn,nyn)
    !third step: determine zn from n21o
    !work2(n21o,nxn,nyn) -> f_new(xn,yn,zn)

    isign=1
    if (iorder(1)<0) isign(1)=2
    if (iorder(2)<0) isign(2)=2
    if (iorder(3)<0) isign(3)=2
    irp=abs(iorder)

    norm=0.0_gp

    !first step
    select case(irp(1))
    case(1) !xn is derived from xo 
       !$omp parallel default(shared) private(scal,k,j,i,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do k=1,ndims_old(3)
          do j=1,ndims_old(2)
             do i=1,ndims_new(1)
                call shift_and_start(irp(1),1,2,3,i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(1)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                !print *,'filter',sum(shf),dt,sum(shf**2),'vals',shf
                !work(j,k+(i-1)*ndims_old(3))
                work(j+ind(2,3,k,i))=convolve(irp(1),k1,j,k,ms,me,&
                     m_isf,shf,ndims_old(1),ndims_old(2),ndims_old(3),f_old,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(2) !xn is derived from yo
       !$omp parallel default(shared) private(scal,k,j,i,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do k=1,ndims_old(3)
          do j=1,ndims_old(1)
             do i=1,ndims_new(1)
                call shift_and_start(irp(1),1,1,3,i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(1)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                !work(j,k+(i-1)*ndims_old(3))
                work(j+ind(1,3,k,i))=convolve(irp(1),j,k1,k,ms,me,&
                     m_isf,shf,ndims_old(1),ndims_old(2),ndims_old(3),f_old,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(3) !xn is derived from zo
       !$omp parallel default(shared) private(scal,k,j,i,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do k=1,ndims_old(2)
          do j=1,ndims_old(1)
             do i=1,ndims_new(1)
                call shift_and_start(irp(1),1,1,2,i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(1)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                !work(k,j+(i-1)*ndims_old(2))
                work(j+ind(1,2,k,i))=convolve(irp(1),j,k,k1,ms,me,&
                     m_isf,shf,ndims_old(1),ndims_old(2),ndims_old(3),f_old,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    end select

    !second step
    select case(irp(1)*10+irp(2))
    case(21) !yp is derived from xo (and xp has been derived from y)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)
                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                work2(k+ind2(irp(3),i,j))=convolve(1,k1,k,i,ms,me,m_isf,shf,&
                     ndims_old(1),ndims_old(3),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(23) !yp is derived from zo (and xp has been derived from y)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                work2(k+ind2(irp(3),i,j))=convolve(2,k,k1,i,ms,me,m_isf,shf,&
                     ndims_old(1),ndims_old(3),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(12) !yp is derived from yo (and xp has been derived from x)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)
                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                !work2(k,i+(j-1)*ndims_new(1))
                work2(k+ind2(irp(3),i,j))=convolve(1,k1,k,i,ms,me,m_isf,shf,&
                     ndims_old(2),ndims_old(3),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(13) !yp is derived from zo (and xp has been derived from x)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)

                !              print *,'value fouund',dt,k1,j

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                !work2(k,i+(j-1)*ndims_new(1))
                work2(k+ind2(irp(3),i,j))=convolve(2,k,k1,i,ms,me,m_isf,shf,&
                     ndims_old(2),ndims_old(3),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(32) !yp is derived from yo (and xp has been derived from z)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                work2(k+ind2(irp(3),i,j))=convolve(2,k,k1,i,ms,me,m_isf,shf,&
                     ndims_old(1),ndims_old(2),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(31) !yp is derived from xo (and xp has been derived from z)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                work2(k+ind2(irp(3),i,j))=convolve(1,k1,k,i,ms,me,m_isf,shf,&
                     ndims_old(1),ndims_old(2),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    end select

    !third step
    !$omp parallel default(shared) private(scal,j,i,k,dt,k1,ms,me,shf,tt)
    allocate(shf(-m_isf:m_isf)) !the usage of f_malloc here can be considered
    !$omp do !!reduction(+:norm)
    do j=1,ndims_new(2)
       do i=1,ndims_new(1)
          do k=1,ndims_new(3)
             call shift_and_start(irp(3),3,2,3,i,j,k,&
                  dt,k1,ms,me,scal)

             call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(3)),shf)
             call redirect(m_isf,ms,me,shf,scal)
             f_new(i,j,k)=convolve(1,k1,i,j,ms,me,m_isf,shf,&
                  ndims_old(irp(3)),ndims_new(1),ndims_new(2),work2,scal)
             !tt=f_new(i,j,k)
             !tt=tt*tt
             !norm=norm+tt
          end do
       end do
    end do
    !$omp end do
    deallocate(shf)
    !$omp end parallel

    !print *,'norm of the function',sqrt(norm)

    call f_free(work)
    call f_free(work2)
    !call f_free(shf)
    call f_release_routine()

  contains

    !index of work array for step 1
    pure function ind(jc2,jc3,i2,i3)
      implicit none
      integer, intent(in) :: jc2,jc3,i2,i3
      integer :: ind

      ind=ndims_old(jc2)*(i2-1)+ndims_old(jc2)*ndims_old(jc3)*(i3-1)

    end function ind

    pure function ind2(jc3,i2,i3)
      implicit none
      integer, intent(in) :: jc3,i2,i3
      integer :: ind2

      ind2=ndims_old(jc3)*(i2-1)+ndims_old(jc3)*ndims_new(1)*(i3-1)

    end function ind2
    
    pure subroutine shift_and_start(ntr,istep,i2,i3,j1,j2,j3,&
         dt,istart,ms,me,alpha)
      use module_base
      implicit none
      integer, intent(in) :: ntr !< id of the dimension to be transformed
      integer, intent(in) :: istep,i2,i3
      integer, intent(in) :: j1,j2,j3
      integer, intent(out) :: istart,ms,me
      real(gp), intent(out) :: dt
      real(gp), intent(out) :: alpha !< scaling to preserve the normalization
      !local variables
      integer :: ivars
!!!! integer(kind=8) :: istart_shift
      real(gp), dimension(3) :: t
      real(gp) :: coord_old
!!!! real(gp) :: tt

      !define the coordinates in the reference frame, which depends on the transformed variables
      t(1)=-centre_new(1)+real(j1-1,gp)*hgrids_new(1) !the first step is always the same
      if (istep >=2) then
         t(2)=-centre_new(2)+real(j2-1,gp)*hgrids_new(2)
      else
         t(2)=-centre_old(i2)+real(j2-1,gp)*hgrids_old(i2)
      end if
      if (istep ==3) then
         t(3)=-centre_new(3)+real(j3-1,gp)*hgrids_new(3)
      else
         t(3)=-centre_old(i3)+real(j3-1,gp)*hgrids_old(i3)
      end if

      !code for the coords
      ivars=1000*istep+100+10*i2+i3

      !define the value of the shift of the variable we are going to transform
      !coordinate that has to be found in the old box, including the shift
      call old_coord(ntr,ivars,rmat,t(1),t(2),t(3),coord_old,alpha)
      coord_old=coord_old-da(ntr)

      !scaling factor to preserve the normalization of the convolution
      alpha=alpha*hgrids_new(ntr)/hgrids_old(ntr)

      !now finds ibar and delta
      call ibar_and_delta(coord_old,centre_old(ntr),hgrids_old(ntr),ndims_old(ntr),&
           centre_new(ntr),hgrids_new(ntr),m_isf,istart,dt,ms,me)      

      !if (sign(1.0_gp,alpha) < 0.0_gp) dt=-dt
!!!!seems that coord_old can be multiplied by  hold/(hnew*a), where a is the scaling between xnew and xold
!!!!central point of the convolution rounded to the grid points
!!!tt=(coord_old+centre_old(ntr)+hgrids_old(ntr))/hgrids_old(ntr)
!!!
!!!!after this points has been identified, it seems we might start finding istart
!!!!and the corresponding delta, followed by the step sizes
!!!if (tt > real(ndims_old(ntr),gp)) then
!!!   istart=ndims_old(ntr)
!!!else
!!!   istart=min(max(1,nint(tt)),ndims_old(ntr))
!!!end if
!!!
!!!!this shift brings the old point in the new reference frame and it defines the local translation.
!!!dt=real(istart,gp)-(coord_old+centre_new(ntr)+hgrids_new(ntr))/hgrids_old(ntr)
!!!
!!!!test: defines the local translation without new reference:
!!!!this shift brings the old point in the new reference frame and it defines the local translation.
!!!!dt=real(istart,gp)-tt
!!!
!!!
!!!!purify the shift to be lower than a multiple of the grid spacing
!!!istart_shift=nint(dt,kind=8)
!!!dt=dt-real(istart_shift,gp)
!!!   !!$istart=istart-istart_shift
!!!   !!$!identify extremes for the convolution
!!!   !!$ms=-min(m_isf,istart-1)
!!!   !!$me=min(m_isf,ndims_old(ntr)-istart)
!!!
!!!istart_shift=int(istart,kind=8)-istart_shift
!!!!identify extremes for the convolution, here a scaling factor should be added
!!!if (istart_shift - 1 > int(m_isf,kind=8)) then
!!!   ms=-m_isf
!!!else if (1 - istart_shift > int(m_isf,kind=8)) then !starting point is too far
!!!   ms=m_isf+1 !loop is disabled
!!!else
!!!   ms= 1 - int(istart_shift)
!!!end if
!!!! check the other extreme
!!!if (int(ndims_old(ntr),kind=8)-istart_shift > int(m_isf,kind=8)) then
!!!   me=m_isf
!!!else if (int(ndims_old(ntr),kind=8)-istart_shift < -int(m_isf,kind=8)) then
!!!   me=-m_isf-1 !disable loop
!!!else
!!!   me=ndims_old(ntr)-int(istart_shift)
!!!end if
!!!!check ms <= me ==> 1-istart <= ndims_old - istart ==> -2*istart <= ndims_old -1 ==> 2*istart >= 1 - ndims_old
!!!!reconvert to default integer value
!!!istart=int(istart_shift)
    end subroutine shift_and_start

    
  end subroutine field_rototranslation3D


  !> Correct the filter for the interpolation if the direction is inverted
  !! after this check put the absolute value in scal
  pure subroutine redirect(m_isf,ms,me,shf,scal)
    implicit none
    integer, intent(in) :: m_isf
    integer, intent(inout) :: ms,me
    real(gp), intent(inout) :: scal
    real(gp), dimension(-m_isf:m_isf), intent(inout) :: shf
    !local variables
!!!! integer :: ishf
!!!! real(gp) :: tt

    !do nothing otherwise
!!$      if (sign(1.0_gp,scal) < 0.0_gp) then
!!$         !swap array 
!!$         do ishf=1,m_isf
!!$            tt=shf(ishf)
!!$            shf(ishf)=shf(-ishf)
!!$            shf(-ishf)=tt
!!$         end do
!!$         !invert convolutions
!!$         ishf=me
!!$         me=-ms
!!$         ms=-ishf
!!$      end if

    scal=sqrt(abs(scal))
  end subroutine redirect



  !> identify the point from which the old data has to be interpolated
  !! and the shift with respect to this point
  !! identify also the start and end points of the convolution
  pure subroutine ibar_and_delta(xold,cold,hold,nold,cnew,hnew,m_isf,ibar,delta,ms,me)
    implicit none
    integer, intent(in) :: nold  !<number of points in of the old direction
    integer, intent(in) :: m_isf !<(half of the) size of the ISF support 
    real(gp), intent(in) :: xold !<previous coordinate in the old reference frame
    real(gp), intent(in) :: cold !<coordinate of the center of rotation in the old reference
    real(gp), intent(in) :: hold !<grid spacing of the old reference frame
    real(gp), intent(in) :: cnew !<coordinate of the center of rotation in the new reference
    real(gp), intent(in) :: hnew !<grid spacing of the new reference frame
    integer, intent(out) :: ibar !< point in the old frame center of the convolution 
    integer, intent(out) :: ms !< starting point of the convolution
    integer, intent(out) :: me !< ending point of the convolution
    real(gp), intent(out) :: delta !< shift of the convolution, defining the filter
    !local variables
    integer(kind=8) :: istart_shift
    real(gp) :: tt

    !seems that xold can be multiplied by  hold/(hnew*a), where a is the scaling between xnew and xold
    !central point of the convolution rounded to the grid points
    tt=(xold+cold+hold)/hold

    !after this point has been identified, we might start finding ibar
    !and the corresponding delta, followed by the step sizes
    if (tt > real(nold,gp)) then
       ibar=nold
    else
       ibar=min(max(1,nint(tt)),nold)
    end if

    !this shift brings the old point in the new reference frame and it defines the local translation.
    delta=real(ibar,gp)-(xold+cnew+hnew)/hold

    !purify the shift to be lower than a multiple of the grid spacing
    istart_shift=nint(delta,kind=8)
    delta=delta-real(istart_shift,gp)

    istart_shift=int(ibar,kind=8)-istart_shift
    !identify extremes for the convolution, here a scaling factor should be added
    if (istart_shift - 1 > int(m_isf,kind=8)) then
       ms=-m_isf
    else if (1 - istart_shift > int(m_isf,kind=8)) then !starting point is too far
       ms=m_isf+1 !loop is disabled
       delta=0.0_gp
    else
       ms= 1 - int(istart_shift)
    end if
    ! check the other extreme
    if (int(nold,kind=8)-istart_shift > int(m_isf,kind=8)) then
       me=m_isf
    else if (int(nold,kind=8)-istart_shift < -int(m_isf,kind=8)) then
       me=-m_isf-1 !disable loop
       delta=0.0_gp
    else
       me=nold-int(istart_shift)
    end if
    !check ms <= me ==> 1-ibar <= nold - ibar ==> -2*ibar <= nold -1 ==> 2*ibar >= 1 - nold
    !reconvert to default integer value. If the loop is disabled this conversion is useless
    ibar=int(istart_shift)

  end subroutine ibar_and_delta


  pure subroutine old_coord(icrd,ivars,rmat,x,y,z,coord,alpha)
    implicit none
    integer, intent(in) :: icrd !<id of the old coordinate to be retrieved
    integer, intent(in) :: ivars !< order of the variables in terms of 1000*istep+first*100+second*10+third
    real(gp), dimension(3,3), intent(in) :: rmat !< rotation matrix entries
    real(gp), intent(in) :: x,y,z !<coordinates to be used for the mapping
    real(gp), intent(out) :: coord
    real(gp), intent(out) :: alpha !<scaling between the new coordinate and the old one xold=alpha*xnew, xnew=x(istep),xold=x(icrd)
    !local variables


    !the possible cases are
    !123 (1,2,3)
    !122 (1
    !113 (2
    !121 (2
    !112 (3)

    coord=0.0_gp
    alpha=1.0_gp
    select case(icrd)
    case(1) !x coordinate
       select case(ivars)
       case(1123)!'xnyozo')
          coord=(x-rmat(2,1)*y-rmat(3,1)*z)/rmat(1,1)
          alpha=1.0_gp/rmat(1,1)
       case(2123)!'xnynzo')
          coord=(rmat(2,2)*x-rmat(2,1)*y+rmat(1,3)*z)/rmat(3,3)
          alpha=-rmat(2,1)/rmat(3,3)
       case(3123)!'xnynzn')
          coord=rmat(1,1)*x + rmat(1,2)*y  + rmat(1,3)*z
          alpha=rmat(1,3)
       case(2122)!'xnynyo')
          coord=(-rmat(3,2)*x + rmat(3,1)*y + rmat(1,3)*z)/rmat(2,3)
          alpha=rmat(3,1)/rmat(2,3)
       end select
    case(2) !y coordinate
       select case(ivars)
       case(1113)!'xnxozo')
          coord=(x-rmat(1,1)*y-rmat(3,1)*z)/rmat(2,1)
          alpha=1.0_gp/rmat(2,1)
       case(2121)!'xnynxo')
          coord=(rmat(3,2)*x-rmat(3,1)*y +rmat(2,3)*z)/rmat(1,3)
          alpha=-rmat(3,1)/rmat(1,3)
       case(2123)!'xnynzo')
          coord=(-rmat(1,2)*x +rmat(1,1)*y + rmat(2,3)*z )/rmat(3,3)
          alpha=rmat(1,1)/rmat(3,3)
       case(3123)!'xnynzn')
          coord=rmat(2,1)*x +rmat(2,2)*y +rmat(2,3)*z
          alpha=rmat(2,3)
       end select
    case(3) !z coordinate
       select case(ivars)
       case(1112)!'xnxoyo')
          coord=(x-rmat(1,1)*y - rmat(2,1)*z)/rmat(3,1)
          alpha=1.0_gp/rmat(3,1)
       case(2121)!'xnynxo')
          coord=(rmat(3,3)*z - rmat(2,2)*x + rmat(2,1)*y)/rmat(1,3)
          alpha=rmat(2,1)/rmat(1,3)
       case(2122)!'xnynyo')
          coord=(rmat(1,2)*x+rmat(3,3)*z-rmat(1,1)*y)/rmat(2,3)
          alpha=-rmat(1,1)/rmat(2,3)
       case(3123)!'xnynzn')
          coord=rmat(3,3)*z + rmat(3,2)*y + rmat(3,1)*x
          alpha=rmat(3,3)
       end select
    end select
  end subroutine old_coord

  pure function convolve(idim,i,j,k,ms,me,m_isf,shf,n1,n2,n3,f_in,alpha)
    implicit none
    integer, intent(in) :: idim !<dimension to be convolved
    integer, intent(in) :: n1,n2,n3,m_isf
    integer, intent(in) :: i,j,k !< starting point of the convolution
    integer, intent(in) :: ms,me !< extremes for the shift
    real(gp), dimension(-m_isf:m_isf), intent(in) :: shf
    real(gp), dimension(n1,n2,n3), intent(in) :: f_in
    real(gp), intent(in) :: alpha !<scaling factor to be applied after convolution
    real(gp) :: convolve
    !local variables
    integer :: l
    real(gp) :: tt

    tt=0.0_gp
    select case(idim)
    case(1)
       do l=ms,me
          tt=tt+shf(l)*f_in(i+l,j,k)
       end do
    case(2)
       do l=ms,me
          tt=tt+shf(l)*f_in(i,j+l,k)
       end do
    case(3)
       do l=ms,me
          tt=tt+shf(l)*f_in(i,j,k+l)
       end do
    end select

    !end of interpolate coefficient
    convolve=alpha*tt

  end function convolve

  pure subroutine define_filter(dt,nrange,nphi,phi,shf)
    use module_base
    implicit none
    integer, intent(in) :: nphi !< number of sampling points of the ISF function (multiple of nrange)
    integer, intent(in) :: nrange !< extension of the ISF domain in dimensionless units (even number)
    real(gp), intent(in) :: dt
    real(gp), dimension(nphi), intent(in) :: phi !< interpolating scaling function array
    real(gp), dimension(-nrange/2:nrange/2), intent(out) :: shf !< interpolating filter to be applied
    !local variables
    integer :: nunit,ish,ipos,m_isf,l,jisf

    m_isf=nrange/2
    !number of points for a unit displacement
    nunit=nphi/nrange 

    !evaluate the shift
    ish=nint(real(nunit,gp)*dt)

    !if (ish /= 0) print *,'dt',dt,ish

    !starting point in the filter definition
    ipos=ish+1
    if (ish<= 0) then
       jisf=-(abs(ish))/nunit-1
    else if (ish > 0) then
       jisf=ish/nunit+1
    else
       jisf=0
    end if
    jisf=jisf-m_isf

    !fill the filters in its nonzero coefficients
    do l=-m_isf,m_isf
       if (jisf >= -m_isf .and. jisf <= m_isf) then
          shf(l)=phi(ipos)
       else
          shf(l)=0.0_gp
       end if
       jisf=jisf+1
       ipos=ipos+nunit
    end do

    !if the sign is negative the filter can be inverted

  end subroutine define_filter



  !> Backward wavelet transform
  !! gives the anti-correlation
  subroutine back_trans_16_reversed(nd,nt,x,y)
    implicit none
    !Arguments
    integer, intent(in) :: nd !< Length of data set
    integer, intent(in) :: nt !< Length of data in data set to be transformed
    real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
    real(kind=8), intent(out) :: y(0:nd-1) !< Output data
    !Local variables
    integer :: i,j,ind

    integer, parameter :: m=18
    real(kind=8), dimension(-m:m) :: ch=(/0d0, 0d0,&
         3.571912260328699082d-6, -1.1450094552100700164d-6, &
         -0.00005642629040127758254d0, 0.00002345539585568117642d0, &
         0.0004069961892884996228d0, -0.0002465534369237166607d0, &
         -0.001634776719899382798d0, 0.00259729967896342247d0, &
         0.006477427625463336123d0, -0.01262044842878062896d0, &
         -0.02535252967734825372d0, 0.02966399618206407251d0, &
         0.06485097060728547963d0, -0.0289320622117497406d0, &
         0.0185085845718848147d0, 0.5048199552943667001d0, &
         0.970046711566057329d0, 0.7212353426722887695d0, &
         0.0294258861485558961d0, -0.2797722999367705543d0, &
         -0.0990303522418633099d0, 0.07410630821538452139d0,&
         0.04680637576666147908d0, -0.011843799423550127927d0, &
         -0.0122154536585793166d0, 0.0010521128108874154748d0, &
         0.00196569149666800115d0, -0.00008582923667387588177d0, &
         -0.0002141180336992365887d0, 3.667434093271785533d-6,&
         0.000011440737665613076119d0, 0d0, 0d0, 0d0, 0d0/)
    real(kind=8), dimension(-m:m) :: cg,cht,cgt

    !******** coefficients for wavelet transform *********************
    do i=-m,m
       cht(i)=0.d0
       cg(i)=0.d0
       cgt(i)=0.d0
    enddo

    ! the normalization is chosen such that a constant function remains the same constant 
    ! on each level of the transform

    cht( 0)=1.D0

    ! g coefficients from h coefficients
    do i=-m,m-1
       cg(i+1)=cht(-i)*(-1.d0)**(i+1)
       cgt(i+1)=ch(-i)*(-1.d0)**(i+1)
    enddo


    do i=0,nt/2-1
       y(2*i+0)=0.d0
       y(2*i+1)=0.d0

       do j=-m/2,m/2-1

          ! periodically wrap index if necessary
          ind=i-j
          loop99: do
             if (ind.lt.0) then 
                ind=ind+nt/2
                cycle loop99
             end if
             if (ind.ge.nt/2) then 
                ind=ind-nt/2
                cycle loop99
             end if
             exit loop99
          end do loop99

          y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
          y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
       end do
    end do

  END SUBROUTINE back_trans_16_reversed

  subroutine my_scaling_function4b2B(itype,nd,nrange,a,x)
    use module_base
    implicit none
    !Arguments
    !Type of interpolating functions
    integer, intent(in) :: itype
    !Number of points: must be 2**nex
    integer, intent(in) :: nd
    integer, intent(out) :: nrange
    real(kind=8), dimension(0:nd), intent(out) :: a
    real(kind=8), dimension(0:nd,2), intent(out) :: x
    !Local variables
    character(len=*), parameter :: subname='scaling_function4b2B'
    real(kind=8), dimension(:), allocatable :: y
    integer :: i,nt,ni

    call f_routine(id=subname)

    !Only itype=8,14,16,20,24,30,40,50,60,100
    select case(itype)
    case(16)
       !O.K.
    case default
       !print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
       !stop
       call f_err_throw('"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100, used:' &
          & // trim(yaml_toa(itype))//'"', err_name='BIGDFT_RUNTIME_ERROR')
    end select
!!$  write(unit=*,fmt="(1x,a,i0,a)") &
!!$       "Use interpolating scaling functions of ",itype," order"

    !Give the range of the scaling function
    !from -itype to itype
    ni=2*itype
    nrange = ni

    y = f_malloc(0.to.nd,id='y')

    ! plot scaling function
    call zero(nd+1,x(0,1))
    call zero(nd+1,y)
    nt=ni
    x(nt/2,1)=1.d0
    loop1: do
       nt=2*nt
       call back_trans_16(nd,nt,x(0,1),y)
       do i=0,nt-1
          x(i,1)=y(i)
       end do
       if (nt.eq.nd) then
          exit loop1
       end if
    end do loop1

    ! plot reversed scaling function
    call zero(nd+1,x(0,2))
    call zero(nd+1,y)
    nt=ni
    x(nt/2,2)=1.d0
    loop2: do
       nt=2*nt
       call back_trans_16_reversed(nd,nt,x(0,2),y)
       do i=0,nt-1
          x(i,2)=y(i)
       end do
       if (nt.eq.nd) then
          exit loop2
       end if
    end do loop2


    !open (unit=1,file='scfunction',status='unknown')
    do i=0,nd
       a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
       !write(1,*) a(i),x(i)
    end do
    !close(1)

    call f_free(y)
    call f_release_routine()
  END SUBROUTINE my_scaling_function4b2B

  
end module reformatting
