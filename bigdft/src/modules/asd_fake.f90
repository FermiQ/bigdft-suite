!> @file
!! Fake routines to be replaced by the asd routine in the case of Spin-dynamics calculation
!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

subroutine asd_pred()
  use dictionaries, only: f_err_throw
  call f_err_throw('You should not enter in the routine asd_pred without'//&
       ' the correct ASD module in the distribution',err_name='BIGDFT_RUNTIME_ERROR')
end subroutine asd_pred

subroutine asd_corr()
  use dictionaries, only: f_err_throw
  call f_err_throw('You should not enter in the routine asd_corr without'//&
       ' the correct ASD module in the distribution',err_name='BIGDFT_RUNTIME_ERROR')
end subroutine asd_corr

subroutine allocate_asd()
  use dictionaries, only: f_err_throw
  call f_err_throw('You should not enter in the routine allocate_asd without'//&
       ' the correct ASD module in the distribution',err_name='BIGDFT_RUNTIME_ERROR')
end subroutine allocate_asd


