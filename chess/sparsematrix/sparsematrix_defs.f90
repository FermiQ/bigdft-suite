!> @file
!!   Basic file which allows to use the sparse matrices
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


!> Module defining the basic structures and routines related to the sparse matrix library
module sparsematrix_defs
  implicit none

  ! This module is public, such that all other modules using this one inherit all modules used in here

  ! Old and new version of the sparse matrix matrix multiplication
  integer,parameter :: MATMUL_NEW = 101
  integer,parameter :: MATMUL_OLD = 102
  integer,parameter :: matmul_version = MATMUL_NEW 

  ! Keywords for the onesided communications
  integer,parameter :: ONESIDED_POST   = 201
  integer,parameter :: ONESIDED_GATHER = 202
  integer,parameter :: ONESIDED_FULL   = 203

  ! Indicate whether the matrix to be applied during the matrix multiplications shall be replicated 
  ! to improve memory acces (at the cost of a larger memory footprint) or not.
  integer,parameter :: MATMUL_REPLICATE_MATRIX = 301
  integer,parameter :: MATMUL_ORIGINAL_MATRIX = 302
  integer,parameter :: matmul_matrix = MATMUL_REPLICATE_MATRIX ! MATMUL_ORIGINAL_MATRIX

end module sparsematrix_defs
