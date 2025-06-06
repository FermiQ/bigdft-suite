!> @file
!! BigDFT package performing ab initio calculation based on wavelets
!! @author
!!    Copyright (C) 2007-2017 BigDFT group<br/>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Main program to calculate electronic structure
program BigDFT
  use module_base
  use bigdft_run
  use public_keys, only: SKIP_RUN
  implicit none
  !input variables
  type(run_objects) :: runObj
  !output variables
  integer :: ierr
  character(len=60) :: posinp_id
  type(state_properties) :: outs
  type(dictionary), pointer :: run,options,runs

  call f_lib_initialize()

  call bigdft_command_line_options(options)

  !-finds the number of taskgroup size
  !-initializes the mpi_environment for each group
  !-decides the radical name for each run
  call bigdft_init(options)

  !case with parser information
  !this key will contain the runs which are associated to the current BigDFT instance
  !run => dict_iter(options .get. 'BigDFT')
  runs = options .get. 'BigDFT'
  nullify(run)
  do while(iterating(run,on=runs)) !associated(run))
     !here a loop on the documents of the input file starts
     !this loop is useful if we want to restart a run without saving the wavefunctions on the disk
     !number of atoms and number of orbitals in the run have to be the same
     !     do while(valid_dataset(runObj,on=run))

     call run_objects_init(runObj,run)
     !skip this run if after initialisation if told to be so
     !in the SKIP_RUN key
     if (dict_get(run,SKIP_RUN,.false.)) cycle

     call init_state_properties(outs,bigdft_nat(runObj))
     !central routine, in bigdft_run, needed also for QM/MM approaches
     call bigdft_get_run_properties(run, posinp_id = posinp_id)
     call process_run(posinp_id,runObj,outs)
     ! Deallocations.
     call deallocate_state_properties(outs)
     !     end do
     call free_run_objects(runObj)
     !run => dict_next(run)
  end do !loop over iconfig

  call dict_free(options)
  call bigdft_finalize(ierr)

  call f_lib_finalize()

END PROGRAM BigDFT
