subroutine release_acceleration_OCL(GPU)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU
  call ocl_clean_command_queue(GPU%queue)
  call ocl_clean(GPU%context)
END SUBROUTINE release_acceleration_OCL

subroutine init_acceleration_OCL(GPU)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU

  call ocl_create_gpu_context(GPU%context)
  !call ocl_create_command_queue(GPU%queue,GPU%context)
  call ocl_build_programs(GPU%context)
  call ocl_create_command_queue_id(GPU%queue,GPU%context,GPU%id_proc)
  call init_event_list
END SUBROUTINE init_acceleration_OCL

subroutine allocate_data_OCL(n1,n2,n3,geocode,nspin,hx,hy,hz,wfd,orbs,GPU)
  use module_base
  use module_types
  implicit none
  character(len=1), intent (in) :: geocode
  integer, intent(in) :: n1,n2,n3,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  character(len=*), parameter :: subname='allocate_data_OCL'
  integer :: i_stat,iorb
  integer :: n1b, n2b, n3b
  integer, dimension(3) :: periodic

  if (geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif 
  if (geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif

  n1b = (n1+1) * 2
  n2b = (n2+1) * 2
  n3b = (n3+1) * 2

  if (periodic(1)==0) then
    n1b = n1b + 2*7 + 15
  endif
  if (periodic(2)==0) then
    n2b = n2b + 2*7 + 15
  endif
  if (periodic(3)==0) then
    n3b = n3b + 2*7 + 15
  endif


  !allocate the number of GPU pointers for the wavefunctions
  !allocate(GPU%psi(orbs%norbp+ndebug),stat=i_stat)
  !call memocc(i_stat,GPU%psi,'GPU%psi',subname)

  !allocate space on the card
  !allocate the compressed wavefunctions such as to be used as workspace
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c);
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f);
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work1)
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work2)
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work3)
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%d)

  if ( orbs%nspinor == 2) then
    call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_i);
    call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_i);
    call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work1_i)
    call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work2_i)
    call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work3_i)
    call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%d_i)
  end if
  !here spin value should be taken into account
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%rhopot_up)
  if( nspin == 2 ) then
    call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%rhopot_down)
  end if



  !allocate and copy the compression-decompression keys
  call ocl_create_read_buffer(GPU%context,wfd%nseg_c*4*2,GPU%keyg_c)
  call ocl_create_read_buffer(GPU%context,wfd%nseg_c*4,GPU%keyv_c)
  call ocl_create_read_buffer(GPU%context,wfd%nseg_f*4*2,GPU%keyg_f)
  call ocl_create_read_buffer(GPU%context,wfd%nseg_f*4,GPU%keyv_f)
  call ocl_enqueue_write_buffer(GPU%queue,GPU%keyg_c,wfd%nseg_c*2*4,wfd%keyg)
  call ocl_enqueue_write_buffer(GPU%queue,GPU%keyv_c,wfd%nseg_c*4,wfd%keyv)
  if (wfd%nseg_f > 0) then
     call ocl_enqueue_write_buffer(GPU%queue,GPU%keyg_f,wfd%nseg_f*2*4,wfd%keyg(1,wfd%nseg_c+1))
     call ocl_enqueue_write_buffer(GPU%queue,GPU%keyv_f,wfd%nseg_f*4,wfd%keyv(wfd%nseg_c+1))
  end if

  !for preconditioner
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_r);
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_r);
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_b);
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_b);
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_d);
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_d);
  if ( orbs%nspinor == 2) then
    call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_r_i);
    call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_r_i);
    call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_b_i);
    call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_b_i);
    call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_d_i);
    call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_d_i);
  end if
  !full_locham stategy (always true for the moment)
  GPU%full_locham=.true.

END SUBROUTINE allocate_data_OCL


subroutine free_gpu_OCL(GPU,orbs,nspin)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspin
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  character(len=*), parameter :: subname='free_gpu_OCL'
  integer :: i_stat,iorb,norbp,i_all
  

  call ocl_release_mem_object(GPU%d)
  call ocl_release_mem_object(GPU%work1)
  call ocl_release_mem_object(GPU%work2)
  call ocl_release_mem_object(GPU%work3)
  call ocl_release_mem_object(GPU%rhopot_up)
  if ( nspin == 2 ) then
    call ocl_release_mem_object(GPU%rhopot_down)
  endif
  call ocl_release_mem_object(GPU%keyg_c)
  call ocl_release_mem_object(GPU%keyv_c)
  call ocl_release_mem_object(GPU%keyg_f)
  call ocl_release_mem_object(GPU%keyv_f)
  call ocl_release_mem_object(GPU%psi_c)
  call ocl_release_mem_object(GPU%psi_f)
  if ( orbs%nspinor == 2) then
    call ocl_release_mem_object(GPU%psi_c_i)
    call ocl_release_mem_object(GPU%psi_f_i)
    call ocl_release_mem_object(GPU%work1_i)
    call ocl_release_mem_object(GPU%work2_i)
    call ocl_release_mem_object(GPU%work3_i)
    call ocl_release_mem_object(GPU%d_i)
  endif
  !for preconditioner
  call ocl_release_mem_object(GPU%psi_c_r)
  call ocl_release_mem_object(GPU%psi_f_r)
  call ocl_release_mem_object(GPU%psi_c_b)
  call ocl_release_mem_object(GPU%psi_f_b)
  call ocl_release_mem_object(GPU%psi_c_d)
  call ocl_release_mem_object(GPU%psi_f_d)
  if ( orbs%nspinor == 2) then
    call ocl_release_mem_object(GPU%psi_c_r_i)
    call ocl_release_mem_object(GPU%psi_f_r_i)
    call ocl_release_mem_object(GPU%psi_c_b_i)
    call ocl_release_mem_object(GPU%psi_f_b_i)
    call ocl_release_mem_object(GPU%psi_c_d_i)
    call ocl_release_mem_object(GPU%psi_f_d_i)
  endif

END SUBROUTINE free_gpu_OCL

subroutine daub_to_isf_OCL(orbs,lr,psi,psi_r,GPU)
  use module_base
  use module_types
  
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  type(GPU_pointers), intent(inout) :: GPU
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)), intent(in) :: psi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(out) :: psi_r
  
  integer, dimension(3) :: periodic
  integer :: isf

  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif
  if (lr%geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif 

  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if


  call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,&
          psi(1))
  call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,&
          psi(isf))

  call ocl_daub_to_isf(GPU%queue,(/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
          periodic,&
          lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
          GPU%psi_c,GPU%psi_f,&
          GPU%work1,GPU%work2,GPU%work3,GPU%d)

  call ocl_enqueue_read_buffer(GPU%queue,GPU%work2,lr%d%n1i*lr%d%n2i*lr%d%n3i*8,&
          psi_r(1))

END SUBROUTINE daub_to_isf_OCL

subroutine isf_to_daub_OCL(orbs,lr,psi_r,psi,GPU)
  use module_base
  use module_types
  
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  type(GPU_pointers), intent(inout) :: GPU
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)), intent(out) :: psi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(in) :: psi_r
  
  integer, dimension(3) :: periodic
  integer :: isf

  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif
  if (lr%geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif 

  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if


  call ocl_enqueue_write_buffer(GPU%queue,GPU%work1,lr%d%n1i*lr%d%n2i*lr%d%n3i*8,&
          psi_r(1))

  call ocl_isf_to_daub(GPU%queue,(/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
          periodic,&
          lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
          GPU%psi_c,GPU%psi_f,&
          GPU%work1,GPU%work2,GPU%work3,GPU%d)

  call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,&
          psi(1))
  call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,&
          psi(isf))

END SUBROUTINE isf_to_daub_OCL




subroutine local_hamiltonian_OCL(iproc,orbs,lr,hx,hy,hz,&
     nspin,pot,psi,hpsi,ekin_sum,epot_sum,GPU,ekin,epot)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor,orbs%norbp), intent(inout) :: psi
  real(wp), dimension(lr%d%n1i,lr%d%n2i,lr%d%n3i,nspin) :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor,orbs%norbp), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  real(gp), dimension(2,orbs%norbp), intent(out) :: ekin
  real(gp), dimension(2,orbs%norbp), intent(out) :: epot
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian_OCL'
  integer :: i_stat,iorb,isf,i
  real(gp), dimension(3) :: hgrids
  integer, dimension(3) :: periodic
  !stream ptr array
  real(kind=8), dimension(orbs%norbp) :: tab_stream_ptr
  real(kind=8) :: stream_ptr_first_trsf
  real(kind=8) :: rhopot
  integer :: n1, n2, n3

  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif 
  if (lr%geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif

  
  n1 = (lr%d%n1+1) * 2
  n2 = (lr%d%n2+1) * 2
  n3 = (lr%d%n3+1) * 2
  if (periodic(1)==0) then
    n1 = n1 + 2*7 + 15
  endif
  if (periodic(2)==0) then
    n2 = n2 + 2*7 + 15
  endif
  if (periodic(3)==0) then
    n3 = n3 + 2*7 + 15
  endif

  call ocl_enqueue_write_buffer(GPU%queue,GPU%rhopot_up,n1*n2*n3*8,pot) 
  if( nspin == 2 ) then
    call ocl_enqueue_write_buffer(GPU%queue,GPU%rhopot_down,n1*n2*n3*8,pot(1,1,1,2)) 
  end if
 
  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if

  hgrids(1)=0.5_gp*hx
  hgrids(2)=0.5_gp*hy
  hgrids(3)=0.5_gp*hz

  epot_sum=0.0_gp
  ekin_sum=0.0_gp

  do iorb=1,orbs%norbp

     if (orbs%spinsgn(orbs%isorb+iorb) > 0.0) then
       rhopot = GPU%rhopot_up
     else
       rhopot = GPU%rhopot_down
     endif
     !if orbs%nspinor /= 1 this implementation should be rediscussed
     if (.not. GPU%full_locham) then
        stop 'ONLY FULL LOCHAM IS IMPLEMENTED!'
     end if

     call ocl_enqueue_write_buffer_async(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,&
          psi(1,iorb)) 
     call ocl_enqueue_write_buffer_async(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,&
          psi(isf,iorb))
     if (orbs%nspinor == 2) then
       call ocl_enqueue_write_buffer_async(GPU%queue,GPU%psi_c_i,lr%wfd%nvctr_c*8,&
            psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+1,iorb)) 
       call ocl_enqueue_write_buffer_async(GPU%queue,GPU%psi_f_i,7*lr%wfd%nvctr_f*8,&
            psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+isf,iorb))
     end if
     !calculate the local hamiltonian
     !WARNING: the difference between full_locham and normal locham is inside
     call ocl_fulllocham_generic_k(GPU%queue,(/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
          (/periodic(1),periodic(2),periodic(3)/),&
          hgrids,&
          (/orbs%kpts(1,orbs%iokpt(iorb)),orbs%kpts(2,orbs%iokpt(iorb)),orbs%kpts(3,orbs%iokpt(iorb))/),&
          lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,& 
          lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,& 
          GPU%psi_c,GPU%psi_f,&
          GPU%psi_c_i,GPU%psi_f_i,&
          rhopot,&
          GPU%work1,GPU%work2,GPU%work3,&
          GPU%work1_i,GPU%work2_i,GPU%work3_i,&
          GPU%d,&
          GPU%d_i,&
          orbs%nspinor,&
          epot(1,iorb),ekin(1,iorb))
!,&
!          epot(iorb,2),ekin(iorb,2))

     call ocl_enqueue_read_buffer_async(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,hpsi(1,iorb))
     call ocl_enqueue_read_buffer_async(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,hpsi(isf,iorb))
     if (orbs%nspinor == 2) then
       call ocl_enqueue_read_buffer_async(GPU%queue,GPU%psi_c_i,lr%wfd%nvctr_c*8,hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+1,iorb))
       call ocl_enqueue_read_buffer_async(GPU%queue,GPU%psi_f_i,7*lr%wfd%nvctr_f*8,hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+isf,iorb))
     end if
     
  end do
  if (.not. ASYNCconv) then
     call ocl_finish(GPU%queue)
     do iorb=1,orbs%norbp
       ekin_sum = ekin_sum + orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*((ekin(1,iorb)+ekin(2,iorb))&
                  - (epot(1,iorb)+epot(2,iorb)))
       epot_sum = epot_sum + orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*(epot(1,iorb)+epot(2,iorb))
     end do
  endif
  
END SUBROUTINE local_hamiltonian_OCL


subroutine finish_hamiltonian_OCL(orbs,ekin_sum,epot_sum,GPU,ekin,epot)
  use module_base
  use module_types
  implicit none
  type(orbitals_data), intent(in) :: orbs
  real(gp), intent(out) :: ekin_sum,epot_sum
  type(GPU_pointers), intent(inout) :: GPU
  real(gp), dimension(2,orbs%norbp), intent(in) :: ekin
  real(gp), dimension(2,orbs%norbp), intent(in) :: epot

  integer :: iorb

  call ocl_finish(GPU%queue)
  ekin_sum=0
  epot_sum=0
  do iorb=1,orbs%norbp
    ekin_sum = ekin_sum + orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*((ekin(1,iorb)+ekin(2,iorb))&
                 - (epot(1,iorb)+epot(2,iorb)))
    epot_sum = epot_sum + orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*(epot(1,iorb)+epot(2,iorb))
  end do
END SUBROUTINE finish_hamiltonian_OCL

subroutine preconditionall_OCL(iproc,nproc,orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero,GPU)
  use module_base
  use module_types
  implicit none
  type(orbitals_data), intent(in) :: orbs
  integer, intent(in) :: iproc,nproc,ncong
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  real(dp), intent(out) :: gnrm,gnrm_zero
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(inout) :: hpsi
  !local variables
  character(len=*), parameter :: subname='preconditionall_OCL'
  integer ::  ierr,iorb,jorb,i_stat,ncplx,i_all,inds,isf,ikpt
  real(wp) :: scpr
  real(gp) :: cprecr,eval_zero,evalmax
  type(GPU_pointers), intent(inout) :: GPU
  type(workarr_precond) :: w
  integer, dimension(3) :: periodic
  real(wp), dimension(:,:), allocatable :: b
  real(gp), dimension(0:7) :: scal
  !stream ptr array
  real(kind=8), dimension(orbs%norbp) :: tab_stream_ptr

  !the eval array contains all the values
  !take the max for all k-points
  !one may think to take the max per k-point
  
  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif 
  if (lr%geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif


  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if
 
     !arrays for the CG procedure
     allocate(b(orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),orbs%norbp+ndebug),stat=i_stat)
     call memocc(i_stat,b,'b',subname)

     gnrm=0.0_dp
     gnrm_zero=0.0_dp
  call allocate_work_arrays(lr%geocode,lr%hybrid_on,orbs%nspinor,lr%d,w)
  if (orbs%norbp >0) ikpt=orbs%iokpt(1)
  do iorb=1,orbs%norbp
     !if it is the first orbital or the k-point has changed calculate the max
     if (orbs%iokpt(iorb) /= ikpt .or. iorb == 1) then
        !the eval array contains all the values
        !take the max for all k-points
        !one may think to take the max per k-point
        evalmax=orbs%eval((orbs%iokpt(iorb)-1)*orbs%norb+1)
        do jorb=1,orbs%norb
           evalmax=max(orbs%eval((orbs%iokpt(iorb)-1)*orbs%norb+jorb),evalmax)
        enddo
        eval_zero=evalmax
        ikpt=orbs%iokpt(iorb)
     end if


       if (orbs%kpts(1,orbs%iokpt(iorb))**2+orbs%kpts(2,orbs%iokpt(iorb))**2+&
           orbs%kpts(3,orbs%iokpt(iorb))**2 > 0.0_gp .or. orbs%nspinor==2 ) then
         ncplx=2
       else
         ncplx=1
       end if
       
        call cprecr_from_eval(lr%geocode,eval_zero,orbs%eval(orbs%isorb+iorb),cprecr)          

        do inds=1,orbs%nspinor,ncplx !the streams should be more if nspinor>1
           !the nrm2 function can be replaced here by ddot
           scpr=nrm2(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),hpsi(1,inds,iorb),1)
           if (orbs%occup(orbs%isorb+iorb) == 0.0_gp) then
              gnrm_zero=gnrm_zero+orbs%kwgts(orbs%iokpt(iorb))*scpr**2
           else
              !write(17,*)'iorb,gnrm',orbs%isorb+iorb,scpr**2
              gnrm=gnrm+orbs%kwgts(orbs%iokpt(iorb))*scpr**2
           end if
           call precondition_preconditioner(lr,ncplx,hx,hy,hz,scal,cprecr,w,&
                hpsi(1,inds,iorb),b(1,iorb))
           
           call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,&
                hpsi(1,inds,iorb))
           call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,&
                hpsi(isf,inds,iorb))
           
           call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c_b,lr%wfd%nvctr_c*8,&
                b(1,iorb))
           call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f_b,7*lr%wfd%nvctr_f*8,&
                b(isf,iorb))

           if(ncplx == 2) then
             call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c_i,lr%wfd%nvctr_c*8,&
                  hpsi(1,inds+1,iorb))
             call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f_i,7*lr%wfd%nvctr_f*8,&
                  hpsi(isf,inds+1,iorb))
           
             call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c_b_i,lr%wfd%nvctr_c*8,&
                  b(1+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),iorb))
             call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f_b_i,7*lr%wfd%nvctr_f*8,&
                  b(isf+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),iorb))
           endif 
           call ocl_preconditioner_generic_k(GPU%queue,&
                (/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
                (/periodic(1),periodic(2),periodic(3)/),&
                (/0.5_gp*hx,0.5_gp*hy,0.5_gp*hz/),&
                (/orbs%kpts(1,orbs%iokpt(iorb)),orbs%kpts(2,orbs%iokpt(iorb)),orbs%kpts(3,orbs%iokpt(iorb))/),&    
                cprecr,&
                ncong,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
                lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
                GPU%psi_c,GPU%psi_f,&
                GPU%psi_c_i,GPU%psi_f_i,&
                GPU%psi_c_r,GPU%psi_f_r,&
                GPU%psi_c_r_i,GPU%psi_f_r_i,&
                GPU%psi_c_b,GPU%psi_f_b,&
                GPU%psi_c_b_i,GPU%psi_f_b_i,&
                GPU%psi_c_d,GPU%psi_f_d,&
                GPU%psi_c_d_i,GPU%psi_f_d_i,&
                GPU%d,GPU%work1,GPU%work2,GPU%work3,&
                GPU%d_i,GPU%work1_i,GPU%work2_i,GPU%work3_i,&
                ncplx)
           

           call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,hpsi(1,inds,iorb))
           call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,hpsi(isf,inds,iorb))
           if ( ncplx == 2 ) then
             call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_c_i,lr%wfd%nvctr_c*8,hpsi(1,inds+1,iorb))
             call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_f_i,7*lr%wfd%nvctr_f*8,hpsi(isf,inds+1,iorb))
           endif

        end do
     end do

  call deallocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,w)

     !end of dynamic repartition
 

  i_all=-product(shape(b))*kind(b)
  deallocate(b,stat=i_stat)
  call memocc(i_stat,i_all,'b',subname)



END SUBROUTINE preconditionall_OCL

subroutine local_partial_density_OCL(iproc,nproc,orbs,&
     nrhotot,lr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc,nrhotot
  type(orbitals_data), intent(in) :: orbs
  integer, intent(in) :: nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
 
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
  real(dp), dimension(lr%d%n1i,lr%d%n2i,nrhotot,nspin), intent(inout) :: rho_p
  type(GPU_pointers), intent(inout) :: GPU
  
  integer:: iorb,iorb_r,i_stat,isf,iaddjmp
  real(kind=8) :: stream_ptr
  real(gp) :: hfac
  integer, dimension(3) :: periodic
  real(kind=8) :: rhopot

  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif 
  if (lr%geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif

  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if

  call set_d(GPU%queue, lr%d%n1i*lr%d%n2i*lr%d%n3i , 1.d-20,  GPU%rhopot_up)
  if ( nspin == 2 ) then
    call set_d(GPU%queue, lr%d%n1i*lr%d%n2i*lr%d%n3i , 1.d-20,  GPU%rhopot_down)
  end if
  !copy the wavefunctions on GPU
  do iorb=1,orbs%norbp*orbs%nspinor
     iorb_r = (iorb-1)/orbs%nspinor + 1
     
    call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,&
             psi(1,iorb)) 
    call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,&
             psi(isf,iorb))
 
    hfac=orbs%kwgts(orbs%iokpt(iorb_r))*orbs%occup(orbs%isorb+iorb_r)/(hxh*hyh*hzh);
    if (orbs%spinsgn(orbs%isorb+iorb_r) > 0.0) then
       rhopot = GPU%rhopot_up
    else
       rhopot = GPU%rhopot_down
    endif
  !calculate the density
   call ocl_locden_generic(GPU%queue, (/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
                          (/periodic(1),periodic(2),periodic(3)/),&
                          hfac,&
                          lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
                          lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
                          GPU%psi_c,GPU%psi_f,&
                          GPU%work1,GPU%work2,GPU%work3,&
                          rhopot)

  
  end do
  !copy back the results and leave the uncompressed wavefunctions on the card
  

  call ocl_enqueue_read_buffer(GPU%queue,GPU%rhopot_up,lr%d%n1i*lr%d%n2i*lr%d%n3i*8,rho_p)
  if( nspin == 2 ) then
    call ocl_enqueue_read_buffer(GPU%queue,GPU%rhopot_down,lr%d%n1i*lr%d%n2i*lr%d%n3i*8,rho_p(1,1,1,2))
  endif

END SUBROUTINE local_partial_density_OCL


