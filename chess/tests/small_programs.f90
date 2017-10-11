!> @file
!!   Flexible test of the matrix power expansion
!! @author
!!   Copyright (C) 2017 CheSS developers
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
!> program checking calls to small operations which are done within the package
program small_things
  use futile
  use fermi_level
  implicit none
  integer, parameter :: UP_=1,DOWN_=2
  integer :: occopt,nkpt
  real(f_double) :: kT,efermi,eTS
  integer, dimension(2) :: norbud
  real(f_double), dimension(:), allocatable :: wgts,eval,occup
  integer, dimension(:), allocatable :: ipiv
  real(f_double), dimension(:), pointer :: p_up,p_dw
  type(dictionary), pointer :: options

  call f_lib_initialize()
  !input variables
  call yaml_argparse(options,&
       '- {name: norb, shortname: n, default: 10, help_string: norb (tuple)}'//f_cr//&
       '- {name: nkpt, shortname: k, default: 1, help_string: k-points}'//f_cr//&
       '- {name: wgts, shortname: w, default: 1.0, help_string: weigths}'//f_cr//&
       '- {name: occopt, shortname: o, default: 1, help_string: smearing method}'//f_cr//&
       '- {name: tel, shortname: t, default: 1.e-3, help_string: temperature}')
  
  norbud=0
  norbud=options//'norb'
  nkpt=options//'nkpt'
  wgts=f_malloc(nkpt,id='wgts')
  wgts=options//'wgts'
  occopt=options//'occopt'
  kT=options//'tel'
  call dict_free(options)
  !determine fermi level starting from a random set of energies
  eval=f_malloc(sum(norbud),id='eval')
  occup=f_malloc0(sum(norbud),id='occup')
  ipiv=f_malloc(sum(norbud),id='ipiv')
  call f_random_number(eval)

  !manipulate separately spin up and spin down
  p_up=>f_subptr(occup,1.to.norbud(UP_))
  p_dw=>f_subptr(occup,from=norbud(UP_)+1,size=norbud(DOWN_))
  !determine the charge
  p_up(1:norbud(UP_)/2)=1.0_f_double
  p_dw(1:norbud(DOWN_)/2)=1.0_f_double

  !now sort the energies
  p_up=>f_subptr(eval,1.to.norbud(UP_))
  p_dw=>f_subptr(eval,from=norbud(UP_)+1,size=norbud(DOWN_))

  call sort_positions(norbud(UP_),p_up,ipiv)
  p_up=p_up(ipiv(norbud(UP_):1:-1))
  call sort_positions(norbud(DOWN_),p_dw,ipiv)
  p_dw=p_dw(ipiv(norbud(DOWN_):1:-1))

  call yaml_new_document()
  call yaml_map('Test energies',eval)
  call yaml_map('Electronic temperature',kT)

  call eval_to_occ(0,1,norbud(1),norbud(2),sum(norbud), nkpt, wgts, &
       eval, occup, .false., .true., kT, occopt, efermi, eTS, &
       norbud(1),norbud(2))

  call yaml_map('Fermi level',efermi)
  call yaml_map('Occupation numbers',occup)
  call yaml_map('TS',eTS)

  call yaml_release_document()
  call f_free(wgts,eval,occup)
  call f_free(ipiv)

  call f_lib_finalize()
end program small_things
