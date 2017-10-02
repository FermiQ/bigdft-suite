!> @file
!!  Routines to do frequencies calculation by finite difference
!! @author
!!    Copyright (C) 2010-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! @todo
!!  - Add higher order for finite difference
!!  - Maybe possibility to use Lanczos to determine lowest frequencies
!!  - Indicate correct formulae for entropy
!!  - Use the directory data (creat_dir_output) for the files hessian.dat and dynamical.dat


!> Calculate vibrational frequencies by frozen phonon approximation.
!! Use a file 'frequencies.res' to restart calculations.
program frequencies

   use module_base
!!$   use module_types
!!$   use module_interfaces
!!$   use m_ab6_symmetry
   use yaml_output
   use bigdft_run
   use dictionaries, only: f_err_throw
   use module_Atoms, only: move_this_coordinate
   implicit none

   !Parameters
   character(len=*), parameter :: subname='frequencies'
   !if |freq1-freq2|<tol_freq, freq1 and freq2 are identical
   real(gp), parameter :: tol_freq=1.d-11
   real(gp), parameter :: Temperature=300.0_gp !< Temperature (300K)
   character(len=*), dimension(3), parameter :: cc = (/ 'x', 'y', 'z' /)
   !File unit
   integer, parameter :: u_hessian=20, u_dynamical=21, u_freq=15, u_hess=35
   real(gp) :: alat,dd,rmass
   character(len=60) :: run_id
   !Input variables
   type(run_objects) :: runObj
   type(state_properties) :: outs
   !Atomic coordinates, forces
   real(gp), dimension(:,:), allocatable :: rxyz0      !< Atomic position of the reference configuration
   real(gp), dimension(:,:), allocatable :: fpos       !< Atomic forces used for the calculation of the Hessian
   real(gp), dimension(:,:), allocatable :: hessian    !< Hessian matrix
   real(gp), dimension(:,:), allocatable :: dynamical  !< Dynamical matrix
   real(gp), dimension(:,:), allocatable :: vectors, vectors0, vectors_old    !< Eigenvectors
   real(gp), dimension(:), allocatable :: eigens, eigens0, eigens_old     !< Real eigenvalues
   real(gp), dimension(:), allocatable :: sort_work    !< To sort the eigenvalues in ascending order
   integer, dimension(:), allocatable :: iperm         !< Array to sort eigenvalues
   integer, dimension(:), allocatable :: kmoves        !< Array which indicates moves to calculate for a given direction
   logical, dimension(:,:), allocatable :: moves       !< logical: .true. if already calculated
   real(gp), dimension(:,:), allocatable :: energies   !< Total energies for all moves
   real(gp), dimension(:,:,:), allocatable :: forces   !< Atomic forces for all moves
   real(gp), dimension(:), allocatable :: fpm, fmm    !< Hessian matrix

   !Function used to determine if the coordinate of the given atom is frozen

   character(len=max_field_length) :: prefix
   integer, dimension(:), allocatable :: ifrztyp0 !< To avoid to freeze the atoms for bigdft_state
   real(gp), dimension(3) :: freq_step
   real(gp) :: zpenergy,freq_exp,freq2_exp,vibrational_entropy,vibrational_energy,total_energy,tij,tji,dsym
   integer :: k,km,ii,jj,ik,imoves,order,n_order
   !integer :: iproc,nproc,igroup,ngroups
   integer :: iat,jat,i,j,ierr,infocode,ity,nconfig,nfree,istart,ios
!NN
   logical :: eigvbasis, eigvbasis_restart
   integer :: nvec, iter, maxiter, jstart, irestart
   real(gp):: fp, fm, alpha, ediff, conv_cut
   real(gp), external :: ddot
   logical :: converged
!
   logical :: exists
   integer :: FREQUENCIES_RUNTIME_ERROR
   !integer, dimension(4) :: mpi_info
   type(dictionary), pointer :: options


   call f_lib_initialize()
   !-finds the number of taskgroup size
   !-initializes the mpi_environment for each group
   !-decides the radical name for each run

   call bigdft_command_line_options(options)

   !-finds the number of taskgroup size
   !-initializes the mpi_environment for each group
   !-decides the radical name for each run
   call bigdft_init(options)

   if (bigdft_nruns(options) > 1) call f_err_throw('runs-file not supported for frequencies executable')

   call f_routine(id=subname)

   inquire(file='eigvbasis',exist=exists)
   eigvbasis=.false.
   if(exists)eigvbasis=.true.
   inquire(file='eigvbasis_restart',exist=exists)
   eigvbasis_restart=.false.
   if(exists)eigvbasis_restart=.true.
  

   !Define the errors for frequencies
   call f_err_define('FREQUENCIES_INPUT_ERROR',&
        'The input file for frequencies is missing.',&
        FREQUENCIES_RUNTIME_ERROR,&
        err_action='Please create it!')
   call f_err_define('FREQUENCIES_ORDER_ERROR',&
        'An invalid value for the order of the finite difference was given.',&
        FREQUENCIES_RUNTIME_ERROR,&
        err_action='Contact the developers')

   !print *,'iconfig,arr_radical(iconfig),arr_posinp(iconfig)',arr_radical(iconfig),arr_posinp(iconfig),iconfig,igroup
   ! Read all input files. This should be the sole routine which is called to initialize the run.
   call run_objects_init(runObj,options//'BigDFT'//0)! trim(run_id), 'posinp')

   ! Read all input files.
   call bigdft_get_run_properties(options//'BigDFT'//0, input_id = prefix)
   inquire(file=trim(prefix)//'.freq',exist=exists)
   if (.not. exists) call f_err_throw('(F) The input file "'//trim(prefix)//'.freq does not exist',&
                          err_name='FREQUENCIES_INPUT_ERROR')
   call frequencies_input_variables_new(bigdft_mpi%iproc,.true.,trim(prefix)//'.freq',runObj%inputs)

   !Order of the finite difference scheme
   order = runObj%inputs%freq_order
   if (order == -1) then
      n_order = 1
      kmoves = f_malloc(src=(/ -1 /),id='kmoves')
   else if (order == 1) then
      n_order = 1
      kmoves = f_malloc(src=(/ 1 /),id='kmoves')
   else if (order == 2) then
      n_order = 2
      kmoves = f_malloc(src=(/ -1, 1 /),id='kmoves')
   else if (order == 3) then
      n_order = 4
      kmoves = f_malloc(src=(/ -2, -1, 1, 2 /),id='kmoves')
   else
      call f_err_throw('(F) Frequencies: This order '//trim(yaml_toa(order))//' is not implemented!',&
           err_name='FREQUENCIES_ORDER_ERROR')
   end if

   ! Allocations
   call init_state_properties(outs, runObj%atoms%astruct%nat)
   rxyz0 = f_malloc((/ 3, runObj%atoms%astruct%nat /), id = 'rxyz0')
   ifrztyp0 = f_malloc(runObj%atoms%astruct%nat, id = 'ifrztyp0')
   moves = f_malloc((/ 1.to.n_order, 0.to.3*runObj%atoms%astruct%nat /),id='moves')
   energies = f_malloc((/ 1.to.n_order, 0.to.3*runObj%atoms%astruct%nat /),id='energies')
   forces = f_malloc((/ 1.to.3*runObj%atoms%astruct%nat, 1.to.n_order, 0.to.3*runObj%atoms%astruct%nat /),id='forces')
   fpos = f_malloc((/ 3*runObj%atoms%astruct%nat, n_order /),id='fpos')
   hessian = f_malloc((/ 3*runObj%atoms%astruct%nat, 3*runObj%atoms%astruct%nat /),id='hessian')
   dynamical = f_malloc((/ 3*runObj%atoms%astruct%nat, 3*runObj%atoms%astruct%nat /),id='dynamical')
!   print *, 'all allocations done me=', bigdft_mpi%iproc

   ! Initialize the Hessian and the dynamical matrix
   hessian = 0.d0
   dynamical = 0.d0
   ! Initialize freq_step (step to move atoms)
   freq_step(1) = runObj%inputs%freq_alpha*runObj%inputs%hx
   freq_step(2) = runObj%inputs%freq_alpha*runObj%inputs%hy
   freq_step(3) = runObj%inputs%freq_alpha*runObj%inputs%hz
!   freq_step(1)=0.001
!   freq_step(2)=0.001
!   freq_step(3)=0.001

   ! Reference positions.
   call vcopy(3*runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, rxyz0(1,1), 1)
   ! Remove frozen atoms in order to have the full atomic forces from bigdft_state
   ! If we want the Hessian and Dynamical matrices only for the freedom degrees not useful
   ! but permit to restart with more degrees of freedom (less frozen atoms)
   ifrztyp0 = runObj%atoms%astruct%ifrztyp
   runObj%atoms%astruct%ifrztyp = 0


!if(eigvbasis_restart)goto 999

   !Initialize the moves using a restart file if present
   !Regenerate it if trouble and indicate if all calculations are done
!   print *, 'to read restart =', bigdft_mpi%iproc
   call frequencies_check_restart(runObj%atoms%astruct%nat,n_order,imoves,moves,energies,forces,freq_step,runObj%atoms%amu,infocode)
   !We ignore at this stage the infocode
!   print *, 'reading restart done =', bigdft_mpi%iproc
!NN
call mpi_barrier(bigdft_mpi%mpi_comm,ierr)

   !Message
   if (bigdft_mpi%iproc == 0) then
      call yaml_map('(F) Frequency moves already calculated',imoves)
      call yaml_map('(F) Total Frequency moves',n_order*3*runObj%atoms%astruct%nat)
   end if

   !Reference state
   if (moves(1,0)) then
      call vcopy(3*outs%fdim, forces(1,1,0), 1, outs%fxyz(1,1), 1)
      outs%energy = energies(1,0)
      infocode=0
   else
      if (bigdft_mpi%iproc == 0) call yaml_comment('(F) Reference state calculation',hfill='=')
      call bigdft_state(runObj,outs,infocode)
      call frequencies_write_restart(0,0,0,runObj%atoms%astruct%rxyz,outs%energy,outs%fxyz, &
           & n_order=n_order,freq_step=freq_step,amu=runObj%atoms%amu)
      moves(:,0) = .true.
      call restart_inputs(runObj%inputs)
   end if
!   print *, 'setting reference state done =', bigdft_mpi%iproc
call mpi_barrier(bigdft_mpi%mpi_comm,ierr)

   if (bigdft_mpi%iproc == 0) then
      call yaml_map('(F) Exit signal for Wavefunction Optimization Finished',infocode)
      call yaml_comment('(F) Start Frequencies calculation',hfill='=')

      !This file contains the Hessian for post-processing: it is regenerated each time.
      call yaml_set_stream(unit=u_hessian,filename=trim(runObj%inputs%dir_output)//'hessian.yaml',&
             position='rewind',record_length=92,istat=ierr,setdefault=.false.,tabbing=0)

      call yaml_map('Step',freq_step,unit=u_hessian)
      call yaml_map('nat',runObj%atoms%astruct%nat,unit=u_hessian)
      call yaml_map('Energy',outs%energy,unit=u_hessian)
      call yaml_map('Forces',outs%fxyz,unit=u_hessian)

      !This file contains the dynamical matrix for post-processing: it is regenerated each time.
      call yaml_set_stream(unit=u_dynamical,filename=trim(runObj%inputs%dir_output)//'dynamical.yaml',&
             position='rewind',record_length=92,istat=ierr,setdefault=.false.,tabbing=0)
      call yaml_map('Step',freq_step,unit=u_dynamical)
      call yaml_map('nat',runObj%atoms%astruct%nat,unit=u_dynamical)
      call yaml_map('Energy',outs%energy,unit=u_dynamical)
      call yaml_map('Forces',outs%fxyz,unit=u_dynamical)
   end if

   !Number of considered degrees of freedom
   nfree = 0
   ! Loop over the atoms for the degrees of freedom
   
!   print *, 'start the basic loop 0 me=', bigdft_mpi%iproc
!call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
   do iat=1,runObj%atoms%astruct%nat

!      print *, "iat =", iat, " me=",bigdft_mpi%iproc, " part -1"
      ! Loop over x, y and z
      do i=1,3
!      print *, "i= ", i, "iat =", iat, " me=",bigdft_mpi%iproc, " part 0.1"
         if (.not.move_this_coordinate(ifrztyp0(iat),i)) then
            if (bigdft_mpi%iproc == 0) call yaml_comment( &
               '(F) The direction '// trim(yaml_toa(i)) // ' of the atom ' // trim(yaml_toa(iat)) // ' is frozen.')
            cycle
         end if
!      print *, "i= ", i, "iat =", iat, " me=",bigdft_mpi%iproc, " part 1"
!call mpi_barrier(bigdft_mpi%mpi_comm,ierr)

         ii = i+3*(iat-1)
         !One more degree of freedom
         nfree = nfree + 1
         if (i==1) then
            !Along x axis
            alat=runObj%atoms%astruct%cell_dim(1)
         else if (i==2) then
            !Along y axis
            alat=runObj%atoms%astruct%cell_dim(2)
         else
            !Along z axis
            alat=runObj%atoms%astruct%cell_dim(3)
         end if
!      print *, "i= ", i, "iat =", iat, " me=",bigdft_mpi%iproc, " part 2.0.0"
         km = 0
         do ik=1,n_order
            k = kmoves(ik)
            !-1-> 1, 1 -> 2, y = ( x + 3 ) / 2
            km = km + 1
            if (moves(km,ii)) then
               !This move is already done. We use the values from the restart file.
               fpos(:,km) = forces(:,km,ii)
!      print *, "km=", km, "ii= ", ii, "ik =", ik, "i= ", i, "iat =", iat, " me=",bigdft_mpi%iproc, " part 2.1.0"
               cycle
            end if
            !Displacement
            dd=real(k,gp)*freq_step(i)
            !We copy atomic positions
            call vcopy(3*runObj%atoms%astruct%nat, rxyz0(1,1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
            if (bigdft_mpi%iproc == 0) then
               call yaml_mapping_open('(F) Move',flow=.true.)
                  call yaml_map('atom',      iat)
                  call yaml_map('direction', k)
                  call yaml_map('axis',      cc(i))
                  call yaml_map('displacement (Bohr)', dd,fmt='(1pe20.10)')
               call yaml_mapping_close()
            end if
            if (runObj%atoms%astruct%geocode == 'P') then
               runObj%atoms%astruct%rxyz(i,iat)=modulo(rxyz0(i,iat)+dd,alat)
            else if (runObj%atoms%astruct%geocode == 'S') then
               runObj%atoms%astruct%rxyz(i,iat)=modulo(rxyz0(i,iat)+dd,alat)
            else
               runObj%atoms%astruct%rxyz(i,iat)=rxyz0(i,iat)+dd
            end if
            call bigdft_state(runObj,outs,infocode)
            call frequencies_write_restart(km,i,iat,runObj%atoms%astruct%rxyz,outs%energy,outs%fxyz)
            call vcopy(3*outs%fdim, outs%fxyz(1,1), 1, fpos(1,km), 1)
            moves(km,ii) = .true.
            call restart_inputs(runObj%inputs)
         end do
!call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!      print *, "i= ", i, "iat =", iat, " me=",bigdft_mpi%iproc, " part 5"
         ! Build the Hessian and the dynamical matrix
         do jat=1,runObj%atoms%astruct%nat
            rmass = amu_emass*sqrt(runObj%atoms%amu(runObj%atoms%astruct%iatype(iat))* &
                 & runObj%atoms%amu(runObj%atoms%astruct%iatype(jat)))
            do j=1,3
               jj = j+3*(jat-1)
               !Force is -dE/dR
               select case(order)
               case(-1)
                  dd = - (outs%fxyz(j,jat) - fpos(jj,1))/freq_step(i)
               case(1)
                  dd = - (fpos(jj,1) - outs%fxyz(j,jat))/freq_step(i)
               case(2)
                  dd = - (fpos(jj,2) - fpos(jj,1))/(2.d0*freq_step(i))
               case(3)
                  dd = - (fpos(jj,4) + fpos(jj,3) - fpos(jj,2) - fpos(jj,1))/(6.d0*freq_step(i))
               case default
                  call f_err_throw('(F) Frequencies: This order '//trim(yaml_toa(order))//' is not allowed!',&
                       err_name='FREQUENCIES_ORDER_ERROR')
               end select
               if (move_this_coordinate(ifrztyp0(jat),j)) hessian(jj,ii) = dd
               dynamical(jj,ii) = dd/rmass
            end do
         end do
!call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!      print *, "jat =", jat, "i= ", i, "iat =", iat, " me=",bigdft_mpi%iproc, " part 6"
         

         if (bigdft_mpi%iproc == 0) then
            call yaml_map('Atom'//trim(yaml_toa(iat))//' Coord.'//trim(yaml_toa(i)),hessian(:,ii),unit=u_hessian)
            call yaml_map('Atom'//trim(yaml_toa(iat))//' Coord.'//trim(yaml_toa(i)),dynamical(:,ii),unit=u_dynamical)
         end if

!call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!      print *,  "i= ", i, "iat =", iat, " me=",bigdft_mpi%iproc, " part 7"

      end do
   end do
!call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!   print *, 'end the basic loop 0 me=', bigdft_mpi%iproc
!   if (bigdft_mpi%iproc == 0) print *, 'end dynamical construction  0 '

   if (bigdft_mpi%iproc == 0) then
      ! Close the files
      call yaml_close_stream(unit=u_hessian)
      call yaml_close_stream(unit=u_dynamical)
   end if

   !Deallocations
   call f_free(fpos)
   call f_free(kmoves)

   !Symmetrization of the dynamical matrix
   !Even if we can calculate more second derivatives, we have only nfree diagonal terms
   dsym = 0.d0
   do i=1,3*runObj%atoms%astruct%nat
      do j=i+1,3*runObj%atoms%astruct%nat
         tij = dynamical(i,j)
         tji = dynamical(j,i)
         !We symmetrize
         dynamical(j,i) = 0.5d0 * (tij+tji)
         dynamical(i,j) = dynamical(j,i)
         dsym = dsym + (tij-tji)**2
      end do
   end do
   !Symmetrization of the hessian
   do i=1,3*runObj%atoms%astruct%nat
      do j=i+1,3*runObj%atoms%astruct%nat
         tij = hessian(i,j)
         tji = hessian(j,i)
         !We symmetrize
         hessian(j,i) = 0.5d0 * (tij+tji)
         hessian(i,j) = hessian(j,i)
      end do
   end do

!   if (bigdft_mpi%iproc == 0) print *, 'end hessian symm  0 '
    !write symmetrized hessian to file
!    open(unit=u_hess,file='hessian_symmetrized.dat')
!    do i=1,3*runObj%atoms%astruct%nat
!        write(u_hess,'(60(1x,es24.17))')(hessian(i,j),j=1,3*runObj%atoms%astruct%nat)
!    enddo
!    close(u_hess)
!   if (bigdft_mpi%iproc == 0) print *, 'end hessian construction  0 '

!999 continue

   !Allocations
!   print *, 'allocating eigens0 me=',bigdft_mpi%iproc
   eigens0    = f_malloc(3*runObj%atoms%astruct%nat,id='eigens0')
!   print *, 'allocating eigens me=',bigdft_mpi%iproc
   eigens    = f_malloc(3*runObj%atoms%astruct%nat,id='eigens')
!   print *, 'allocating vectors0 me=',bigdft_mpi%iproc
   vectors0   = f_malloc((/ 3*runObj%atoms%astruct%nat, 3*runObj%atoms%astruct%nat /),id='vectors0')
!   print *, 'allocating vectors_old me=',bigdft_mpi%iproc
   vectors_old   = f_malloc((/ 3*runObj%atoms%astruct%nat, 3*runObj%atoms%astruct%nat /),id='vectors_old')
!   print *, 'allocating vectors me=',bigdft_mpi%iproc
   vectors   = f_malloc((/ 3*runObj%atoms%astruct%nat, 3*runObj%atoms%astruct%nat /),id='vectors')
!   print *, 'allocating sort_work me=',bigdft_mpi%iproc
   sort_work = f_malloc(3*runObj%atoms%astruct%nat,id='sort_work')
!   print *, 'allocating iperm me=',bigdft_mpi%iproc
   iperm     = f_malloc(3*runObj%atoms%astruct%nat,id='iperm')
!   if (bigdft_mpi%iproc == 0) print *, 'allocated'

!if(eigvbasis_restart)goto 9999


   !print *, 'to solve hessian'
   !Diagonalise the hessian
   call solve(hessian,3*runObj%atoms%astruct%nat,eigens0,vectors0)

!   print *, 'solve hessian completed ', bigdft_mpi%iproc

!   if (bigdft_mpi%iproc == 0) print *, 'solve - hessian done'
!      call HDIAG(hessian,3*runObj%atoms%astruct%nat,3*runObj%atoms%astruct%nat,0,eigens0,vectors0,i)
   !print *, 'done solve hessian'
   !Sort eigenvalues in descending order (use abinit routine sort_dp)
   sort_work=eigens0
   do i=1,3*runObj%atoms%astruct%nat
      iperm(i)=i
   end do
!   print *, 'sort_dp  completed ', bigdft_mpi%iproc
   call abi_sort_dp(3*runObj%atoms%astruct%nat,sort_work,iperm,tol_freq)  !sagar
   if (bigdft_mpi%iproc == 0) then
      call yaml_comment('(F) Hessian results',hfill='=')
      call yaml_map('(F) Full Hessian Matrix Calculation',nfree == 3*runObj%atoms%astruct%nat)
      call yaml_map('(F) Number of calculated degrees of freedom',nfree)
      call yaml_map('(F) Hessian Eigenvalues',eigens0(iperm(3*runObj%atoms%astruct%nat:1:-1)),fmt='(1pe20.10)')
!NN-testing
!      if(eigvbasis)call write_restart_eigvbasis(-1,1,1,3*runObj%atoms%astruct%nat,hessian(1,1),eigens0(1),vectors0(1,1)) !note: hessian is not printed
   endif

   call f_free(hessian)
!   print *, 'f_free  completed ', bigdft_mpi%iproc


   !print *, 'to solve dynamical'
   !Diagonalise the dynamical matrix
   call solve(dynamical,3*runObj%atoms%astruct%nat,eigens,vectors)
!   print *, 'solve - dynamical done me=', bigdft_mpi%iproc
!     i=0
!      call HDIAG(dynamical,3*runObj%atoms%astruct%nat,3*runObj%atoms%astruct%nat,0,eigens,vectors,i)

   !print *, 'done solve dynamical'
   !Sort eigenvalues in descending order (use abinit routine sort_dp)
   sort_work=eigens
   do i=1,3*runObj%atoms%astruct%nat
      iperm(i)=i
   end do
   call abi_sort_dp(3*runObj%atoms%astruct%nat,sort_work,iperm,tol_freq)  !sagar
!   print *, 'sort_dp2 after dynamical done me=', bigdft_mpi%iproc

   if (bigdft_mpi%iproc == 0) then
      call yaml_comment('(F) Frequencies results',hfill='=')
      call yaml_map('(F) Full Dynamical Matrix Calculation',nfree == 3*runObj%atoms%astruct%nat)
      call yaml_map('(F) Number of calculated degrees of freedom',nfree)
      if (nfree == 3*runObj%atoms%astruct%nat) call yaml_map('(F) Dynamical matrix symmetrization',dsym)
      call yaml_map('(F) Eigenvalues',eigens(iperm(3*runObj%atoms%astruct%nat:1:-1)),fmt='(1pe20.10)')
      do i=1,3*runObj%atoms%astruct%nat
         if (eigens(i)<0.0_dp) then
            eigens(i)=-sqrt(-eigens(i))
         else
            eigens(i)= sqrt( eigens(i))
         end if
      end do
      call yaml_map('(F) Frequencies (Hartree)', eigens(iperm(3*runObj%atoms%astruct%nat:1:-1)),fmt='(1pe20.10)')
      call yaml_map('(F) Frequencies (cm-1)',    eigens(iperm(3*runObj%atoms%astruct%nat:1:-1))*Ha_cmm1,fmt='(f13.2)')
      call yaml_map('(F) Frequencies (THz)',     eigens(iperm(3*runObj%atoms%astruct%nat:1:-1))*Ha_THz,fmt='(f13.2)')
      ! Build frequencies.xyz in descending order. Use the v_sim format
      open(unit=u_freq,file='frequencies.xyz',status="unknown")
      do i=3*runObj%atoms%astruct%nat,1,-1
         write(u_freq,'(1x,i0,1x,1pe20.10,a)') runobj%atoms%astruct%nat
         write(u_freq,'(1x,a,i0,a,1pe20.10,a,0pf13.2,a,f13.2,a)') 'Mode ',i,': freq=', &
            & eigens(iperm(i)),' Ha,',eigens(iperm(i))*Ha_cmm1,' cm-1,',eigens(iperm(i))*Ha_THz,' Thz'
         ! Build the vector of the associated phonon
         do iat=1,runobj%atoms%astruct%nat
            ity=runobj%atoms%astruct%iatype(iat)
            write(u_freq,'(1x,a,1x,100(1pe20.10))') &
               &   trim(runobj%atoms%astruct%atomnames(ity)),rxyz0(:,iat),(vectors(3*(iat-1)+j,iperm(i)),j=1,3)
         end do
      end do
      close(unit=u_freq)

!      if(eigvbasis)call write_restart_eigvbasis(-1,1,1,3*runObj%atoms%astruct%nat,dynamical(1,1),eigens(1),vectors(1,1)) !note: hessian is not printed
!
!      !Vibrational entropy of the molecule
!      ! See : http://www.codessa-pro.com/descriptors/thermodynamic/entropy.htm)
!      !       http://www.ncsu.edu/chemistry/franzen/public_html/CH795N/lecture/XIV/XIV.html
!      !Zero-point energy
!      zpenergy = 0.0_gp
!      vibrational_energy=0.0_gp
!      vibrational_entropy=0.0_gp
!      !iperm: ascending order
!      !Remove zero frequencies:
!      if (nfree == 3*runobj%atoms%astruct%nat) then
!         if (runobj%atoms%astruct%nat == 2) then
!            istart = 6
!         else
!            istart = 7
!         end if
!      else
!         istart=3*runObj%atoms%astruct%nat-nfree+1
!      end if
!      do i=istart,3*runObj%atoms%astruct%nat
!         freq_exp=exp(eigens(iperm(i))*Ha_K/Temperature)
!         freq2_exp=exp(-eigens(iperm(i))*Ha_K/(2.0_gp*Temperature))
!         zpenergy=zpenergy+0.5_gp*eigens(iperm(i))
!         vibrational_energy=vibrational_entropy+eigens(iperm(i))*(0.5_gp+1.0_gp/(freq_exp-1.0_gp))
!         vibrational_entropy=vibrational_entropy + eigens(iperm(i))*freq2_exp/(1.0_gp-freq2_exp) - log(1.0_gp-freq2_exp)
!      end do
!      !Multiply by 1/kT
!      vibrational_entropy=vibrational_entropy*Ha_K/Temperature
!      total_energy=energies(1,0)+vibrational_energy
!      call yaml_map('(F) Zero-point energy (cm-1 and Hartree)', (/ zpenergy*Ha_cmm1, zpenergy /),fmt='(1pe20.10)')
!      call yaml_map('(F) Considered Temperature (Kelvin)',      Temperature,fmt='(f5.1)')
!      call yaml_map('(F) Vibrational entropy',                  vibrational_entropy,fmt='(1pe22.10)')
!      call yaml_map('(F) Vibrational Energy (cm-1 and Hartree)', &
!           &                     (/ vibrational_energy*Ha_cmm1, vibrational_energy/),fmt='(1pe20.10)')
!      call yaml_map('(F) Total Energy (Hartree)',               total_energy,fmt='(1pe22.10)')
   end if
!   print *, 'printing from master node done me=', bigdft_mpi%iproc

!9999 continue
   if(eigvbasis_restart)call f_free(hessian)
!   print *, 'freeing hessian done', bigdft_mpi%iproc
  

!NN
  if(eigvbasis)then
!  print *, 'starting eigvbasis me = ', bigdft_mpi%iproc 
  conv_cut=1.e-10
  maxiter=2
  ediff=1.0e-5
  nvec=3*runObj%atoms%astruct%nat
  call dcopy(nvec*nvec,vectors0,1,vectors_old,1)
!  print *, 'starting eigvbasis allocating eigens_old me = ', bigdft_mpi%iproc 
  eigens_old = f_malloc(nvec,id='eigens_old')
  call dcopy(nvec,eigens0,1,eigens_old,1)
!  print *, 'starting eigvbasis allocating hessian me= ', bigdft_mpi%iproc 
  hessian = f_malloc((/ nvec, nvec /),id='hessian') 
  fpm = f_malloc((/ nvec/),id='fpm') 
  fmm = f_malloc((/ nvec/),id='fmm') 
!  print *, 'starting eigvbasis allocated fpm/fmm me= ', bigdft_mpi%iproc 
  iter=1
  jstart=1
  hessian=0.d0
  dynamical=0.d0
  irestart=0
!  print *, 'starting reading eigvbasis new restart me= ', bigdft_mpi%iproc 
  if (bigdft_mpi%iproc == 0) &
     print *, 'EigVBasis | maxiter = ', maxiter
  inquire(file='eigvbasis_restart_new',exist=exists)
  if(exists)then
    open(109,file='eigvbasis_restart_new',status='old')
     irestart=1
!    read(109,*)irestart
     if (bigdft_mpi%iproc == 0) &
     print *, 'EigVBasis - will read restart eigvbasis_restart_new'
  end if
!  print *, 'done reading eigvbasis new restart me= ', bigdft_mpi%iproc 
!  call read_restart_eigvbasis(iter,jstart,nvec,hessian,eigens0,vectors0)
  if (bigdft_mpi%iproc == 0) &
     print *, 'after restart| iter =', iter, 'jstart =',jstart, 'nvec =', nvec
  optloop: do 
     if(iter.gt.maxiter)exit optloop
!calculate the derivative dE/dR
     do j=jstart,nvec
        if (bigdft_mpi%iproc == 0) & 
          print *, ' EigVBasis| iter =', iter, '/',maxiter,'  j=',j,'/',nvec
!        if(abs(eigens0(j))<5.e-3) then
!debug
   if(iter.gt.2)then
        alpha=dsqrt(ediff/abs(eigens0(j)))
        if (bigdft_mpi%iproc == 0) & 
          print *, ' EigVBasis| e= ',eigens0(j),' alpha =', alpha
          alpha=max(0.001d0,alpha)
          alpha=min(0.5d0,alpha)
!        else
!          alpha=0.001
!        end if
   else
        alpha=0.001d0
   end if
        if (bigdft_mpi%iproc == 0) & 
          call yaml_map('(F) Eigen Vector Basis: alpha',               alpha,fmt='(1pe22.10)')
        !displace along vector j
        if(irestart.eq.1)then
          do i=1,runObj%atoms%astruct%nat
            read(109,*,iostat=ios)outs%fxyz(1:3,i)
          end do
          if(ios.ne.0)then 
              irestart=0
              close(109)
          end if
          if(ios.eq.0)then 
              if(bigdft_mpi%iproc == 0) & 
                 print *, ' EigVBasis| j = ', j, ' Force (Delta) is RESTARTED' 
          end if
        end if
        if(irestart.eq.0)then
          if(bigdft_mpi%iproc == 0) & 
            print *, ' EigVBasis| j = ', j, ' Force (Delta) is  COMPUTED' 
          call displace(nvec,alpha,vectors0(1,j),rxyz0,runObj%atoms%astruct%rxyz)
          call bigdft_state(runObj,outs,infocode)
          call write_eigvbasis_restart_new(bigdft_mpi%iproc,109,runObj%atoms%astruct%nat,outs%fxyz)
        end if
        !force along vector i => (dE/dVector_i) for changes along Vector_j
        do i=1,nvec
          fpm(i)=ddot(nvec,outs%fxyz,1,vectors0(1,i),1)
        end do
        if(irestart.eq.1)then
          do i=1,runObj%atoms%astruct%nat
            read(109,*,iostat=ios)outs%fxyz(1:3,i)
          end do
          if(ios.ne.0)then 
              irestart=0
              close(109)
          end if
          if(ios.eq.0)then 
              if(bigdft_mpi%iproc == 0) & 
                 print *, ' EigVBasis| j = ', j, ' Force (-Delta) is RESTARTED' 
          end if
        end if
        if(irestart.eq.0)then
          if(bigdft_mpi%iproc == 0) & 
            print *, ' EigVBasis| j = ', j, ' Force (-Delta) is  COMPUTED' 
          call displace(nvec,-alpha,vectors0(1,j),rxyz0,runObj%atoms%astruct%rxyz)
          call bigdft_state(runObj,outs,infocode)
          call write_eigvbasis_restart_new(bigdft_mpi%iproc,109,runObj%atoms%astruct%nat,outs%fxyz)
        end if
        do i=1,nvec
          fmm(i)=ddot(nvec,outs%fxyz,1,vectors0(1,i),1)
        end do
        do i=1,nvec
          hessian(i,j)=-(fpm(i)-fmm(i))/(2.d0*alpha)
        end do
!        if (bigdft_mpi%iproc == 0) &
!         call write_restart_eigvbasis(1,iter,j,nvec,hessian(1,j),eigens0(1),vectors0(1,1))
     end do
!     irestart=0 ! this is done not to read the restart file again.
     iter=iter+1
     jstart=1

     fm=0.d0
     do i=1,nvec
       do j=1,nvec
          tij=hessian(i,j)
          tji=hessian(j,i)
          fm=max(fm,dabs(tij-tji))
          hessian(i,j)=(tij+tji)/2.d0
          hessian(j,i)=(tij+tji)/2.d0
          dynamical(i,j)=(tij+tji)/2.d0
          dynamical(j,i)=(tij+tji)/2.d0
       end do
     end do


     if (bigdft_mpi%iproc == 0) &
         call yaml_map('(F) Max Hessian Off Diag Variation ',fm,fmt='(1pe22.10)')
   
!call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!    print *, 'going to solve Hessian me=',bigdft_mpi%iproc

   !Diagonalise the hessian
     call solve(hessian,3*runObj%atoms%astruct%nat,eigens0,vectors0)

!call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!    print *, 'done solve Hessian me=',bigdft_mpi%iproc

     if (bigdft_mpi%iproc == 0)then 
        open(unit=u_freq,file='frequencies_Hessian_Eigvbasis.xyz',status="unknown")
        do i=3*runObj%atoms%astruct%nat,1,-1
           write(u_freq,'(1x,i0,1x,1pe20.10,a)') runobj%atoms%astruct%nat
           write(u_freq,'(1x,a,i0,a,1pe20.10)') 'Mode ',i,': freq=', eigens0(i)
           ! Build the vector of the associated phonon
           do iat=1,runobj%atoms%astruct%nat
              ity=runobj%atoms%astruct%iatype(iat)
              write(u_freq,'(1x,a,1x,100(1pe20.10))') &
                 trim(runobj%atoms%astruct%atomnames(ity)),rxyz0(:,iat),&
                   (vectors0(3*(iat-1)+j,i),j=1,3)
           end do
        end do
        close(unit=u_freq)
     end if



    !print *, 'Hessian in Q basis:'
    !do i=1,nvec
    !write(*,'(6f16.8)') hessian(i,1:6)
    !end do
    !print *, 'Eigen of Hessian in Q basis:'
    !write(*,'(6f16.8)') eigens0(1:6)


!NN-TODO
!      CALL HDIAG(DE2MAT,NQ,NQ,0,EIGWER,EIGVEC,NROT)
!      i=0
!      call HDIAG(hessian,3*runObj%atoms%astruct%nat,3*runObj%atoms%astruct%nat,0,eigens0,vectors0,i)

!     if (bigdft_mpi%iproc == 0) &
!        call write_restart_eigvbasis(-1,iter,nvec,nvec,hessian(1,1),eigens0(1),vectors0(1,1))
     !Sort eigenvalues in descending order (use abinit routine sort_dp)
!     sort_work=eigens0
!     do i=1,3*runObj%atoms%astruct%nat
!       iperm(i)=i
!     end do
!
  call reconstruct_dynamical(nvec,dynamical,vectors0,vectors_old,eigens0)
!call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!    print *, 'done reconstruct dynamical me=',bigdft_mpi%iproc
     do i=1,nvec
       do j=1,nvec
          tij=dynamical(i,j)
          tji=dynamical(j,i)
          dynamical(i,j)=(tij+tji)/2.d0
          dynamical(j,i)=(tij+tji)/2.d0
       end do
     end do
     if (bigdft_mpi%iproc == 0)then 
        open(unit=u_freq,file='frequencies_Hessian_Eigvbasis_Transform.xyz',status="unknown")
        do i=3*runObj%atoms%astruct%nat,1,-1
           write(u_freq,'(1x,i0,1x,1pe20.10,a)') runobj%atoms%astruct%nat
           write(u_freq,'(1x,a,i0,a,1pe20.10)') 'Mode ',i,': freq=', eigens0(i)
           ! Build the vector of the associated phonon
           do iat=1,runobj%atoms%astruct%nat
              ity=runobj%atoms%astruct%iatype(iat)
              write(u_freq,'(1x,a,1x,100(1pe20.10))') &
                 trim(runobj%atoms%astruct%atomnames(ity)),rxyz0(:,iat),&
                   (vectors0(3*(iat-1)+j,i),j=1,3)
           end do
        end do
        close(unit=u_freq)
     end if



!NN------------------------------------------------------------------------------------
     !constructing the mass weighted hessian matrix
     do iat=1,runObj%atoms%astruct%nat
       do i=1,3
          ii = i+3*(iat-1)
          do jat=1,runObj%atoms%astruct%nat
            rmass = amu_emass*sqrt(runObj%atoms%amu(runObj%atoms%astruct%iatype(iat))* &
                                   runObj%atoms%amu(runObj%atoms%astruct%iatype(jat)))
            do j=1,3
               jj = j+3*(jat-1)
               dynamical(jj,ii)=dynamical(jj,ii)/rmass
            end do
          end do
       end do
     end do 
    !print *, 'Dynamical in Q->R basis:'
    !do i=1,nvec
    !write(*,'(6f16.8)') dynamical(i,1:6)
    !end do

!    Diagonalise the dynamical
    eigens=0.d0
    vectors=0.d0
    call solve(dynamical,3*runObj%atoms%astruct%nat,eigens,vectors)
    if (bigdft_mpi%iproc == 0)then 
      print *, '*********Eigenvalues of Dynamical in Q->R basis:********'
      write(*,'(6f16.8)') eigens(1:3*runObj%atoms%astruct%nat)
    end if
    !print *, 'Eigenvectors of Dynamical in Q->R basis:'
    !do i=1,nvec
    !write(*,'(6f16.8)') vectors(i,1:6)
    !end do

!  square root of eigen values
   do i=1,3*runObj%atoms%astruct%nat
      if (eigens(i)<0.0_dp) then
         eigens(i)=-sqrt(-eigens(i))
      else
         eigens(i)= sqrt( eigens(i))
      end if
   end do
!NN------------------------------------------------------------------------------------
     !print *, 'finished eigns sqrt'
!
!     call sort_dp(3*runObj%atoms%astruct%nat,sort_work,iperm,tol_freq)
     if (bigdft_mpi%iproc == 0) then
        call yaml_comment('(F) EIGV: Hessian results',hfill='=')
        call yaml_map('(F) EIGV| Optimization Loop',iter)
!        call yaml_map('(F) EIGV| Hessian Eigenvalues',eigens0(iperm(3*runObj%atoms%astruct%nat:1:-1)),fmt='(1pe20.10)')

        call yaml_map('(F) EIGV| Dynamical Sqrt(Eigenvalues)',eigens(1:nvec),fmt='(1pe20.10)')
        call yaml_map('(F) EIGV| Frequencies (cm-1)', eigens(1:nvec)*Ha_cmm1,fmt='(f13.2)')
!
        open(unit=u_freq,file='frequencies_EVBasis.xyz',status="unknown")
        do i=3*runObj%atoms%astruct%nat,1,-1
           write(u_freq,'(1x,i0,1x,1pe20.10,a)') runobj%atoms%astruct%nat
           write(u_freq,'(1x,a,i0,a,1pe20.10)') 'Mode ',i,': freq=', eigens(i)
           ! Build the vector of the associated phonon
           do iat=1,runobj%atoms%astruct%nat
              ity=runobj%atoms%astruct%iatype(iat)
              write(u_freq,'(1x,a,1x,100(1pe20.10))') &
                 trim(runobj%atoms%astruct%atomnames(ity)),rxyz0(:,iat),&
                   (vectors(3*(iat-1)+j,i),j=1,3)
           end do
        end do
        close(unit=u_freq)
     endif
!     stop
     call check_conv_eigenvalues(nvec,eigens_old,eigens,conv_cut,converged,fm)
     if(bigdft_mpi%iproc == 0)call yaml_map('(F) Optimization: Change in eigenvalue',fm,fmt='(1pe20.10)')
     if((bigdft_mpi%iproc == 0).and.converged)&
         call yaml_map('(F) Optimization is converged at iteration',iter)
     if(converged)exit optloop
!
  end do optloop

!     do iat=1,runObj%atoms%astruct%nat
!        do i=1,3
!           ii= i+3*(iat-1)
!           do jat=1,runObj%atoms%astruct%nat
!              rmass = amu_emass*sqrt(runObj%atoms%amu(runObj%atoms%astruct%iatype(iat))* &
!                      runObj%atoms%amu(runObj%atoms%astruct%iatype(jat)))
!              do j=1,3
!                 jj = j+3*(jat-1)
!                 dynamical(jj,ii)=dynamical(jj,ii)/rmass  
!              end do
!           end do
!        end do
!     end do
!    !Diagonalise the dynamical matrix
!    call solve(dynamical,3*runObj%atoms%astruct%nat,eigens,vectors)
!   !Sort eigenvalues in descending order (use abinit routine sort_dp)
!   sort_work=eigens
!   do i=1,3*runObj%atoms%astruct%nat
!      iperm(i)=i
!   end do
!   call sort_dp(3*runObj%atoms%astruct%nat,sort_work,iperm,tol_freq)
!   call yaml_map('(F) EIGVBasis  Eigenvalues',eigens(iperm(3*runObj%atoms%astruct%nat:1:-1)),fmt='(1pe20.10)')
!   do i=1,3*runObj%atoms%astruct%nat
!      if (eigens(i)<0.0_dp) then
!         eigens(i)=-sqrt(-eigens(i))
!      else
!         eigens(i)= sqrt( eigens(i))
!      end if
!   end do
!   call yaml_map('(F) EIGVBasis Frequencies (Hartree)', eigens(iperm(3*runObj%atoms%astruct%nat:1:-1)),fmt='(1pe20.10)')
!   call yaml_map('(F) EIGVBasis Frequencies (cm-1)', eigens(iperm(3*runObj%atoms%astruct%nat:1:-1))*Ha_cmm1,fmt='(f13.2)')
!
     if (bigdft_mpi%iproc == 0) &
         call yaml_map('(F) EigVBasis Optimization Completed',.true.)
     call f_free(fpm)
     call f_free(fmm)
     call f_free(eigens_old)
     call f_free(hessian)
  end if
!NN deallocations

   ! De-allocations
   call f_free(vectors0) !NN
   call f_free(vectors_old) !NN
   call f_free(eigens0)  !NN

   call f_free(rxyz0)
   call f_free(ifrztyp0)

   call deallocate_state_properties(outs)

   call f_free(dynamical)
   call f_free(eigens)
   call f_free(vectors)

   call f_free(iperm)
   call f_free(sort_work)

   call f_free(moves)
   call f_free(energies)
   call f_free(forces)

   call f_release_routine()

   call free_run_objects(runObj)
   call dict_free(options)
   call bigdft_finalize(ierr)

   call f_lib_finalize()


contains
!  subroutine reconstruct_dynamical2(nvec,hessian,vectors,vectors0)
!  implicit none
!  integer :: nvec
!  real(gp), intent (inout) :: hessian(nvec,nvec)
!
!  real(gp) :: tmp
!  real(gp), allocatable :: hessian_new(:,:), vectors_new(:,:)
!  integer :: i, j, l, m
!
!  hessian_new = f_malloc((/ nvec, nvec /),id='vectors_new')
!  vectors_new = f_malloc((/ nvec, nvec /),id='u')
!
!  do i=1,nvec
!     tmp=0.d0
!     do l=1,nvec
!       tmp=tmp+vectors(i,l)*vectors0(i,l)
!     end do
!  end do 
!
!
!  do i=1,nvec
!    do j=1,nvec
!        hessian(i,j)=hessian_new(i,j)
!    end do
!  end do
!
!  call f_free(hessian_new) 
!  call f_free(u) 
!!
!  end subroutine reconstruct_dynamical2
  subroutine reconstruct_dynamical(nvec,hessian,vectors,vectors0,eigens)
!SImilarity transform to change H in Q basis to H' in Q' basis using unitary
!transform U=Q.Q'
! TODO: use dsptrf
  implicit none
  integer :: nvec
  real(gp), intent (inout) :: hessian(nvec,nvec)
  real(gp), intent (inout) ::  vectors(nvec,nvec),vectors0(nvec,nvec), eigens(nvec)

  real(gp) :: tmp
  real(gp), allocatable :: hessian_new(:,:), u(:,:)
  integer :: i, j, l, m

  hessian_new = f_malloc((/ nvec, nvec /),id='vectors_new')
  u = f_malloc((/ nvec, nvec /),id='u')
! print *, 'in reco =1'

!  do iat=1,natoms
!   do ii=1,3 
!     i=3*(iat-1)+ii
!     do jat=1,natoms
!       do jj=1,3
!         j=3*(iat-1)+ii
!         tmp=0.d0
!         do l=1,nvec
!            do m=1,nvec
!               tmp=tmp+hessian(l,m)*vectors(l,i)*vectors(m,j)
!            end do
!         end do
!         hessian_new(i,j)=tmp
!       end do
!     end do
!   end do
!  end do

!finding Unitary transform
   u = 0.d0
   do i=1,nvec
    do j=1,nvec
      do l=1,nvec
         u(l,i) = u(l,i)+ vectors(j,i) *vectors0(l,j)
      enddo
    enddo
   enddo
   do i=1,nvec
    do j=1,nvec
      vectors(j,i)=u(j,i)
   end do
  end do

!   do i=1,nvec
!     do j=1,nvec
!      u(i,j)=dot_product(vectors(i,1:nvec),vectors0(j,1:nvec))
!     end do
!   end do
! print *, 'in reco =2'

!Similarity Transform H' = U H U-1 = U H I U-1 = U H Q Q-1 U-1 = U E U-1 = U E U !(as U-1 = U* = U for real)
   do i=1,nvec
     do j=1,nvec
        hessian_new(i,j)=0.d0
        do l=1,nvec
        hessian_new(i,j)=hessian_new(i,j)+u(i,l)*eigens(l)*u(j,l)
     end do
     end do
   end do
! print *, 'in reco =3'

!   do i=1,nvec
!     do j=1,nvec
!        tmp=0.d0
!        do l=1,nvec
!           do m=1,nvec
!               !tmp=tmp+hessian(l,m)*vectors(l,i)*vectors(m,j)
!               tmp=tmp+hessian(l,m)*u(l,i)*u(m,j)
!           end do
!        end do
!        hessian_new(i,j)=tmp
!     end do
!   end do
      
  do i=1,nvec
    do j=1,nvec
        hessian(i,j)=hessian_new(i,j)
    end do
  end do
! print *, 'in reco =4'

  call f_free(hessian_new) 
  call f_free(u) 
!
  end subroutine reconstruct_dynamical

   subroutine check_conv_eigenvalues(nvec,eigens_old,eigens0,conv_cut,converged,mag)
   implicit none
   integer, intent(in) :: nvec
   real(gp), intent(in) :: eigens_old(nvec),eigens0(nvec), conv_cut 
   real(gp), intent(out) :: mag
   logical, intent(out) :: converged

   converged=.false.
   mag=0.d0
   do i=1,nvec
     mag=max(mag,dabs(eigens_old(nvec)-eigens0(nvec)))
   end do
   if(mag<conv_cut)converged=.true.
   end subroutine check_conv_eigenvalues

   !> Make displacements along the eigen vector
   subroutine displace(nvec,alpha,vector,rxyz0,rxyz)
   implicit none
   integer, intent(in) :: nvec
   real(gp), intent(in) :: alpha, vector(nvec), rxyz0(nvec)
   real(gp), intent(out) ::  rxyz(nvec)

   integer :: i
   do i=1,nvec
     rxyz(i) = rxyz0(i)+alpha*vector(i)
   end do
   end subroutine displace

   !> Solve the dynamical matrix
   subroutine solve(dynamical,n,eigens,vectors)
      implicit none
      integer, intent(in) :: n
      real(gp), intent(inout) :: dynamical(n,n)
      real(gp), intent(out) :: eigens(n),vectors(n,n)
      !Local variables
      character(len=*), parameter :: subname = "solve"
      integer :: info,lwork
      real(gp), dimension(:), allocatable :: work

      call f_routine(id=subname)
      lwork=100*n
      work=f_malloc(lwork,id='work')

!    print *, 'vectors inside solve:'
!    print *, 'n =', n
!
!    do i=1,n
!    write(*,'(6f16.8)') dynamical(i,1:6)
!    end do

      call dsyev('V','U',n,dynamical,n,eigens,work,lwork,info)
!    print *, 'eig vec after dsyev'
!    do i=1,n
!    write(*,'(6f16.8)') dynamical(i,1:6)
!    end do
!    print *, 'eig value after dsyev'
!    write(*,'(6f16.8)') eigens(1:6)

      vectors = dynamical

      if (info /= 0) then
         call yaml_warning('(F) Error from the routine dgeev: info=' // trim(yaml_toa(info)))
      end if

      !Put to zero if < 1.d-16
      do i=1,n
         do j=1,n
            if (abs(vectors(j,i)) < 1.d-16) vectors(j,i)=0.d0
         end do
      end do
      !de-allocation
      call f_free(work)

   END SUBROUTINE solve

   subroutine write_eigvbasis_restart_new(iproc,iunit,nat,fxyz)
   implicit none
   integer :: iproc, iunit, nat
   real(gp) :: fxyz(3,nat)
   integer :: i
   logical ::isopen

   inquire(unit=iunit,opened=isopen)
   if(isopen)close(iunit)

   if(iproc.ne.0)return

   open(unit=iunit,file='eigvbasis_restart_new',position='append')
   do i=1,nat
      write(iunit,'(3e20.10)')fxyz(1:3,i)
   end do
   close(iunit)
   print *, 'EigVBasis RESTART written'
   end subroutine write_eigvbasis_restart_new

   !> Check the restart file
   subroutine frequencies_check_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,ierror)
      implicit none
      !Arguments
      integer, intent(in) :: nat     !< Number of atoms
      integer, intent(in) :: n_order !< Order of the finite difference
      logical, dimension(n_order,0:3*nat), intent(out) :: moves            !< Contains moves already done
      real(gp), dimension(n_order,0:3*nat), intent(out) :: energies        !< Energies of the already moves
      real(gp), dimension(3*nat,n_order,0:3*nat), intent(out) :: forces    !< Forces of the already moves
      real(gp), dimension(3), intent(in) :: freq_step     !< Frequency step in each direction
      integer, intent(out) :: imoves                      !< Number of frequency already calculated   
      real(gp), dimension(:), intent(out) :: amu          !< Atomic masses
      integer, intent(out) :: ierror                      !< 0 means all calculations are done
      !Local variables
      character(len=*), parameter :: subname = "frequencies_check_restart"
      !We read the file
      call frequencies_read_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,ierror)
      !if (ierror /= 0) then
      !   !If error, we write a new file
      !   call frequencies_write_new_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,ierror)
      !end if
      !Check also if all calculations are done.
   end subroutine frequencies_check_restart


   !> Read the restart file associated to the calculation of the frequencies
   subroutine frequencies_read_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,ierror)
      implicit none
      !Arguments
      integer, intent(in) :: nat     !< Number of atoms
      integer, intent(in) :: n_order !< Order of the finite difference
      logical, dimension(n_order,0:3*nat), intent(out) :: moves            !< Contains moves already done
      real(gp), dimension(n_order,0:3*nat), intent(out) :: energies        !< Energies of the already moves
      real(gp), dimension(3*nat,n_order,0:3*nat), intent(out) :: forces    !< Forces of the already moves
      real(gp), dimension(3), intent(in) :: freq_step     !< Frequency step in each direction
      integer, intent(out) :: imoves                      !< Number of frequency already calculated
      real(gp), dimension(:), intent(out) :: amu          !< Atomic masses
      integer, intent(out) :: ierror                      !< Error when reading the file
      !Local variables
      character(len=*), parameter :: subname = "frequencies_read_restart"
      logical :: exists
      integer, parameter :: iunit = 15
      real(gp), dimension(3) :: steps
      real(gp), dimension(:), allocatable :: rxyz,fxyz
      real(gp) :: etot
      integer :: km,i,iat,ii,i_order

      call f_routine(id=subname)
      !Initialize by default to false
      imoves=0
      moves = .false.

      !Test if the file does exist.
      inquire(file='frequencies.res', exist=exists)
      if (.not.exists) then
         !There is no restart file.
         call f_zero(energies)
         if (bigdft_mpi%iproc == 0) call yaml_map('(F) File "frequencies.res" present',.false.)
         !Code error non zero
         ierror = -1
         return
      else
         if (bigdft_mpi%iproc == 0) call yaml_map('(F) File "frequencies.res" present',.true.)
      end if

      !We read the file
      open(unit=iunit,file='frequencies.res',status='old',form='unformatted',iostat=ierror)
      !First line is data for coherency of the calculation
      read(unit=iunit,iostat=ierror) i_order,steps,amu
      if (ierror /= 0) then
         !Read error, we exit
         if (bigdft_mpi%iproc == 0) then
            close(unit=iunit)
            call yaml_warning('(F) Error when reading the first line of "frequencies.res"')
         end if
         call f_release_routine()
         return
      else
         if (steps(1) /= freq_step(1) .or. steps(2) /= freq_step(2) .or. steps(3) /= freq_step(3)) then
            if (bigdft_mpi%iproc == 0) call yaml_warning('(F) The step to calculate frequencies is not the same: stop.')
            stop
         end if

         if (i_order > n_order) then
            if (bigdft_mpi%iproc == 0) then 
               call yaml_warning('(F) The number of points per direction is bigger in the "frequencies.res" file.')
               call yaml_warning('(F) Increase the order of the finite difference scheme')
            end if
            stop
         end if
      end if
      
      !Allocations
      rxyz=f_malloc(3*nat,id='rxyz')
      fxyz=f_malloc(3*nat,id='fxyz')

      !Read the reference state
      read(unit=iunit,iostat=ierror) iat,etot,rxyz,fxyz
      if (ierror /= 0 .or. iat /= 0) then
         !Read error, we assume that it is not calculated
         if (bigdft_mpi%iproc == 0) call yaml_map('(F) Reference state calculated in the "frequencies.res" file',.false.)
      else
         if (bigdft_mpi%iproc == 0) call yaml_map('(F) Reference state calculated in the "frequencies.res" file',.true.)
         energies(:,0) = etot
         forces(:,1,0) = fxyz
         moves(:,0) = .true.
      end if
      do
         read(unit=iunit,iostat=ierror) km,i,iat,rxyz,etot,fxyz
         if (ierror /= 0) then
            !Read error, we exit
            if (bigdft_mpi%iproc == 0) then
               close(unit=iunit)
               !Error if all moves are not read
               if (imoves < 3*nat+1) call yaml_warning('(F) The file "frequencies.res" is incomplete!')
            end if
            exit
         end if
         ii = i + 3*(iat-1)
         imoves = imoves + 1
         energies(km,ii) = etot
         forces(:,km,ii) = fxyz
         moves(km,ii) = .true.
      end do
      close(unit=iunit)

      !Deallocations
      call f_free(rxyz)
      call f_free(fxyz)

      call f_release_routine()

   END SUBROUTINE frequencies_read_restart



   !> write the full restart file
   !subroutine frequencies_write_new_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,ierror)
   !   implicit none
   !   !arguments
   !   integer, intent(in) :: nat     !< number of atoms
   !   integer, intent(in) :: n_order !< order of the finite difference
   !   logical, dimension(n_order,0:3*nat), intent(in) :: moves         !< contains moves already done
   !   real(gp), dimension(n_order,0:3*nat), intent(in) :: energies     !< energies of the already moves
   !   real(gp), dimension(3*nat,n_order,0:3*nat), intent(in) :: forces !< forces of the already moves
   !   real(gp), dimension(3), intent(in) :: freq_step    !< frequency step in each direction
   !   integer, intent(out) :: imoves                     !< number of frequency already calculated
   !   real(gp), dimension(:), intent(in) :: amu          !< atomic masses
   !   integer, intent(out) :: ierror                     !< error when reading the file
   !   !local variables
   !   integer, parameter :: iunit = 15

   !   if (bigdft_mpi%iproc ==0 ) then
   !      !this file is used as a restart
   !      open(unit=iunit,file='frequencies.res',status="unknown",form="unformatted")
   !      write(unit=iunit) n_order,freq_step,amu

   !      write(unit=iunit) 0,outs%energy,rxyz0,outs%fxyz
   !      do iat=1,runobj%atoms%astruct%nat
   !         if (ifrztyp0(iat) == 1) then
   !            if (bigdft_mpi%iproc == 0) call yaml_comment('(F) the atom ' // trim(yaml_toa(iat)) // ' is frozen.')
   !            cycle
   !         end if
   !         do i=1,3
   !            ii = i+3*(iat-1)
   !            km = 0
   !            do ik=1,n_order
   !               km = km + 1
   !               if (moves(km,ii)) then
   !                  write(unit=iunit) km,i,iat,outs%energy,rxyz0,outs%fxyz
   !               end if
   !            end do
   !         end do
   !      end do
   !      close(unit=iunit)
   !   end if
   !end subroutine frequencies_write_new_restart


   !> Write one move in the file restart (only moves==.true.)
   subroutine frequencies_write_restart(km,i,iat,rxyz,etot,fxyz,n_order,freq_step,amu)
      implicit none
      !Arguments
      integer, intent(in) :: km,i,iat
      real(gp), dimension(:,:), intent(in) :: rxyz
      real(gp), intent(in) :: etot
      real(gp), dimension(:,:), intent(in) :: fxyz
      integer, intent(in), optional :: n_order
      real(gp), intent(in), optional :: freq_step(3)
      real(gp), dimension(:), intent(in), optional :: amu
      !Local variables
      integer, parameter :: iunit = 15

      if (km == 0 .and. .not.(present(n_order).and.present(freq_step).and.present(amu))) then
         if (bigdft_mpi%iproc == 0) write(*,*) "Bug for use of frequencies_write_restart"
         if (bigdft_mpi%iproc == 0) call yaml_warning("(F) Bug for use of frequencies_write_restart")
         stop
      end if

      if (bigdft_mpi%iproc ==0 ) then
         !This file is used as a restart
         open(unit=iunit,file='frequencies.res',status="unknown",form="unformatted",position="append")
         if (km == 0) then
            write(unit=iunit) n_order,freq_step,amu
            write(unit=iunit) 0,etot,rxyz,fxyz
         else
            write(unit=iunit) km,i,iat,etot,rxyz,fxyz
         end if
         close(unit=iunit)
      end if
   END SUBROUTINE frequencies_write_restart


   subroutine restart_inputs(inputs)
     use module_input_keys
     use public_enums, only: ENUM_MEMORY
     implicit none
     !Argument
     type(input_variables), intent(inout) :: inputs
     !inputs%inputPsiId=1
     
     call inputpsiid_set_policy(ENUM_MEMORY,inputs%inputPsiId)

   END SUBROUTINE restart_inputs


   !> Integrate forces (not used)
   subroutine integrate_forces(iproc,n_moves) !n(c) energies,forces (arg:2,3)
   
      use module_base
   
      implicit none
      !Arguments
      integer, intent(in) :: iproc,n_moves
      !n(c) real(gp), intent(in) :: energies(n_moves)
      !n(c) real(gp), intent(in) :: forces(3*nat,n_moves)
      !Local variables
      character(len=*), parameter :: subname = "integrate_forces"
      real(gp), dimension(:), allocatable :: weight
      !n(c) real(gp) :: path
      integer :: i
   
      !Allocation
      weight = f_malloc(n_moves,id='weight')
   
      !Prepare the array of the correct weights of the iteration steps
      if (mod(n_moves,2).ne.1) then
         if (iproc == 0) write(*,*) 'the number of iteration steps has to be odd'
         stop
      end if
      weight(1)=1.d0/3.d0
      weight(2)=4.d0/3.d0
      do i=3,n_moves-2,2
         weight(i)=2.d0/3.d0
         weight(i+1)=4.d0/3.d0
      enddo
      weight(n_moves)=1.d0/3.d0
   
      !Start integration
      !n(c) path = 0_gp
      do i=1,n_moves
      end do
   
      !De-allocation
      call f_free(weight)
   
   END SUBROUTINE integrate_forces

!NN
  subroutine read_restart_eigvbasis(iter,jstart,nvec,hessian,eigens0,vectors0)
  implicit none
  integer  :: iter, jstart
  integer  ::  nvec
  real(gp) :: hessian(nvec,nvec),eigens0(nvec),vectors0(nvec*nvec)

  integer :: ios,iunit=15
  logical :: exists
  integer :: icall=0
  save    :: icall
   
  iter=1
  jstart=1

  icall = icall + 1
  if(icall.gt.1)return !only read the frequencies once


  open(iunit,file='frequencies_eigvbasis_1.res',status='unknown',form='unformatted')  
  inquire(file='frequencies_eigvbasis_1.res', exist=exists)
  if (.not.exists) then
    !There is no restart file.
    if (bigdft_mpi%iproc == 0) call yaml_map('(F) File "frequencies_eigvbasis_1.res" present',.false.)
    return
  else
    if (bigdft_mpi%iproc == 0) call yaml_map('(F) File "frequencies_eigvbasis_1.res" present',.true.)
  end if
  read(iunit,iostat=ios) (eigens0(i),i=1,nvec)
  read(iunit,iostat=ios) (vectors0(i),i=1,nvec*nvec)
  if((ios.ne.0).and.(bigdft_mpi%iproc == 0)) &
     call yaml_map('(F) File "frequencies_eigvbasis_1.res" present',.true.)
  close(iunit)
!
  open(iunit,file='frequencies_eigvbasis_2.res',status='unknown',form='unformatted')  
  inquire(file='frequencies_eigvbasis_2.res', exist=exists)
  if (.not.exists) then
    !There is no restart file.
    if (bigdft_mpi%iproc == 0) call yaml_map('(F) File "frequencies_eigvbasis_2.res" present',.false.)
    return
  else
    if (bigdft_mpi%iproc == 0) call yaml_map('(F) File "frequencies_eigvbasis_2.res" present',.true.)
  end if
  do i=1,nvec
    read(iunit,iostat=ierr) iter, jstart, hessian(i,1:nvec)
    if(ierr.ne.0)exit 
    if(jstart.eq.nvec)iter=iter+1 ! go to next opt. iteration, if j-loop is complete, i.e. jstart == nvec
  end do
  close(iunit)

  end subroutine read_restart_eigvbasis

  subroutine write_restart_eigvbasis(icontrol,iter,jstart,nvec,hessian,eigens0,vectors0)
  implicit none
  integer :: iter, jstart
  integer :: icontrol,nvec
  real(gp):: hessian(nvec),eigens0(nvec),vectors0(nvec*nvec)

  integer :: iunit=15
  logical :: exists
   
  if(icontrol.eq.-1)then
     open(iunit,file='frequencies_eigvbasis_1.res',status='unknown',form='unformatted')  
     write(iunit) eigens0(1:nvec)
     write(iunit) (vectors0(i),i=1,nvec*nvec)
     close(iunit)
     print *, 'wrote eigvbasis restart with : ONLY vectors and eigens'
  else if(icontrol.eq.0)then
     stop 'unknown icontrol'
  else
     if(jstart.eq.1)then
        open(iunit,file='frequencies_eigvbasis_2.res',status='unknown',form='unformatted') 
     else
        open(iunit,file='frequencies_eigvbasis_2.res',status='unknown',form='unformatted',position='append') !Append the "j-th hessian column  
     end if
     write(iunit) iter, jstart, hessian(1:nvec) !hessian is just some numbers; they will not be used as jstart=1
     close(iunit)
     print *, 'wrote eigvbasis restart with 1| jstart=',jstart, 'iter=',iter
  end if
  end subroutine write_restart_eigvbasis

END PROGRAM frequencies


!> Read the input variables needed for the frequencies calculation.
!! Every argument should be considered as mandatory.
subroutine frequencies_input_variables_new(iproc,dump,filename,in)
  use module_base
  use module_types
  use module_input
  implicit none
  !Arguments
  type(input_variables), intent(inout) :: in
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc
  logical, intent(in) :: dump
  !Local variables
  logical :: exists
  !n(c) integer, parameter :: iunit=111

  !Frequencies parameters
  call input_set_file(iproc,dump,trim(filename),exists,'Frequencies Parameters')  
  !if (exists) in%files = in%files + INPUTS_FREQ
  !call the variable, its default value, the line ends if there is a comment

  !Read in%freq_alpha (possible 1/64)
  call input_var(in%freq_alpha,'1/64',ranges=(/0.0_gp,1.0_gp/),&
       comment="Step size factor (alpha*hgrid)")
  !Read the order of finite difference scheme

  call input_var(in%freq_order,'2',exclusive=(/-1,1,2,3/),&
       comment="Order of the difference scheme")
  !Read the index of the method

  call input_var(in%freq_method,'1',exclusive=(/1/),&
       comment="Method used (only possible value=1)")
  call input_free((iproc == 0) .and. dump)

END SUBROUTINE frequencies_input_variables_new

