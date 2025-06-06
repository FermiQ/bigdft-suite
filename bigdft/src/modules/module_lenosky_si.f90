!> @file
!! Empirical potential for Si type Lenosky
!! @author
!!    Copyright (C) 2014-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!  Copyright (C) 2001-2002 Stefan Goedecker, CEA Grenoble
!!  This file is distributed under the terms of the
!!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .


!> Module defining quantities for Lensoky Si potential
module module_lenosky_si
use module_defs, only:gp
private

public lenosky_si_shift,lenosky_si
public init_lensic

logical, save :: initialized=.false.
real(gp), save :: alat_int(3)

contains

subroutine init_lensic(paramset,paramfile,geocode,units)
    use module_base
    use yaml_output
    implicit none
    !parameters
    character(len=*), intent(in) :: paramset
    character(len=*), intent(in) :: paramfile
    character(len=*), intent(in) :: geocode
    character(len=*), intent(in) :: units
    !local
    logical :: exists
    integer :: u
    call yaml_comment('Initializing lensic',hfill='-')
    call yaml_mapping_open('Reference Papers')
    call yaml_scalar('T. Lenosky, Modelling. Simul. Mater. Sci. Eng. 8, 825 (2000)')
    call yaml_scalar('S. Goedecker, Comp. Phys. Commun. 148, 124 (2002)')
    call yaml_mapping_close()


    initialized=.false.

    if(trim(geocode)/='F')then
        call f_err_throw('lensic only works with free '//&
             'boundary conditions. Specified boundary conditions are: '//&
             trim(adjustl(geocode)))
    endif

    if(trim(paramfile)/='none')then
        !lensic is only an ad-hoc wrapper for the preiodic lensoky
        !routine. We therefore need to specify cell legnths.
        !We don't use the cell values from the structure files
        !(xyz, ascii) because they are zeroed by the parsing routine
        !if bc=free. Setting the bc to periodic in the structure files
        !is also no solution, since sepcific parts of the code
        !are not getting executed in the periodic case (like the elimination
        !of roations).
        inquire(file=trim(adjustl(paramfile)),exist=exists)
        if(.not.exists)then
            call f_err_throw('Parameter file '//trim(adjustl(paramfile))//&
                 ' not found.')
        endif
        u=f_get_free_unit()
        open(unit=u,file=trim(adjustl(paramfile)))
            read(u,*) alat_int(1), alat_int(2), alat_int(3)
        close(u)
        if (units=='angstroem' .or. units=='angstroemd0') then
           ! if Angstroem convert to Bohr                                    
           alat_int(1)=alat_int(1)/Bohr_Ang                              
           alat_int(2)=alat_int(2)/Bohr_Ang                              
           alat_int(3)=alat_int(3)/Bohr_Ang                              
        else if  (units=='atomic' .or. units=='bohr'  .or.&  
             units== 'atomicd0' .or. units== 'bohrd0') then 
             !do nothing 
        else if (units == 'reduced') then                            
           !assume that for reduced coordinates cell size is in bohr         
        else                                                                 
           call f_err_throw('Length units in input file unrecognized.' // &  
                'Recognized units are angstroem or atomic = bohr',err_id=BIGDFT_INPUT_VARIABLES_ERROR)
           return                                                            
           !write(*,*) 'length units in input file unrecognized'             
           !write(*,*) 'recognized units are angstroem or atomic = bohr'     
           !stop                                                             
        endif  
        call yaml_map('Read cell parameters from file',trim(adjustl(paramfile)))
        initialized=.true.
    else
        select case(trim(paramset))
        case('default')
            call f_err_throw('No "default" parameter set for lensic defined.')
        case default
            call f_err_throw('Following parameter set for lensic force field '//&
                'is unknown: '//trim(adjustl(paramset)))
        end select
    endif
end subroutine init_lensic

subroutine lenosky_si(nat,alat,rat,fat,epot)
    implicit none
    integer::nat
    real(8) :: alat(3)
    real(8), intent(in)::rat(3,nat)
    real(8) :: fat(3,nat),epot
    real(8) :: coord=0.d0,ener_var=0.d0,coord_var=0.d0,count=0.d0

    call lenosky(nat,alat,rat,fat,epot,coord,ener_var,coord_var,count)
end subroutine lenosky_si

subroutine lenosky_si_shift(nat,rat,fat,epot)
    use module_base
    implicit none
    integer::nat,iat
    real(8), intent(in)::rat(3,nat)
    real(8) :: fat(3,nat),epot
    real(8) :: coord=0.d0,ener_var=0.d0,coord_var=0.d0,count=0.d0
    real(8) :: pos(3,nat),cmx,cmy,cmz

    if(.not.initialized)then
        call f_err_throw('Potential "lensic" not initialized',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif

    pos=rat
    !shift to center of cell:
        ! center of mass
    cmx=0.d0 ; cmy=0.d0 ; cmz=0.d0
    do iat=1,nat
        cmx=cmx+pos(1,iat)
        cmy=cmy+pos(2,iat)
        cmz=cmz+pos(3,iat)
    enddo
    cmx=cmx/nat ; cmy=cmy/nat ; cmz=cmz/nat
    pos(1,:)=pos(1,:)-cmx+0.5d0*alat_int(1)
    pos(2,:)=pos(2,:)-cmy+0.5d0*alat_int(2)
    pos(3,:)=pos(3,:)-cmz+0.5d0*alat_int(3)
    call lenosky(nat,alat_int,pos,fat,epot,coord,ener_var,coord_var,count)
end subroutine lenosky_si_shift


    subroutine lenosky(nat,alat,rxyz0,fxyz,ener,coord,ener_var,coord_var,count)
    use module_defs, only:gp
!     Evaluates the LENOSKY silicon potential with linear scaling
!     If publishable results are obtained with this program, citing the
!     following 2 references is very much appreciated:
!     1.  T. Lenosky, Modelling. Simul. Mater. Sci. Eng. 8, 825 (2000)
!     2.  S. Goedecker, Comp. Phys. Commun. 148, 124 (2002)
!
!     Parallelized using OpenMP
! Good Compiler options (last option only if paralellization with OpenMp desired)
! IBM Power3
!  xlf90_r -O2 -qarch=pwr3 -qtune=pwr3 -qmaxmem=-1 -qsmp=omp
! Dec Alpha
! f90 -O2 -arch ev67 -pipeline -omp
! Intel Pentium
!  ifc -w -xW -O2 -openmp


!  Copyright (C) 2001-2002 Stefan Goedecker, CEA Grenoble
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!
!     input: - "nat": number of atoms
!            - "alat": lattice constants of the orthorombic box
!               containing the particles
!            - "rxyz0": atomic positions in Angstroem.
!               If an atom is outside the box the program will bring it back
!               into the box by translations through alat
!     output:- "fxyz": forces in eV/A
!            - "ener": total energy in eV
!            - "coord": average coordination number
!            - "ener_var": variance of the energy/atom
!            - "coord_var": variance of the coordination number
!     I/Oput:- "count": is increased by one per call, has to be initialized
!               to 0.0_gp before first call of lenosky
        implicit real(gp) (a-h,o-z)
!$      interface
!$        integer ( kind=4 ) function omp_get_num_threads ( )
!$        end function omp_get_num_threads
!$      end interface
!$      interface
!$        integer ( kind=4 ) function omp_get_thread_num ( )
!$        end function omp_get_thread_num
!$      end interface

        dimension rxyz0(3,nat),fxyz(3,nat),alat(3)
        real(gp), ALLOCATABLE, DIMENSION(:,:) :: rxyz
        integer, ALLOCATABLE, DIMENSION(:,:) :: lsta
        integer, ALLOCATABLE, DIMENSION(:) :: lstb
        integer, ALLOCATABLE, DIMENSION(:) :: lay
        integer, ALLOCATABLE, DIMENSION(:,:,:,:) :: icell
        real(gp), ALLOCATABLE, DIMENSION(:,:) :: rel
        real(gp), ALLOCATABLE, DIMENSION(:,:) :: txyz
        real(gp), ALLOCATABLE, DIMENSION(:,:) :: f2ij,f3ij,f3ik

!        tmax_phi= 0.4500000d+01
!        cut=tmax_phi
        cut= 0.4500000d+01

        if (count.eq.0)  open(unit=10,file='lenosky.mon',status='unknown')
        count=count+1.0_gp

! linear scaling calculation of verlet list
        ll1=int(alat(1)/cut)
        if (ll1.lt.1) stop 'alat(1) too small'
        ll2=int(alat(2)/cut)
        if (ll2.lt.1) stop 'alat(2) too small'
        ll3=int(alat(3)/cut)
        if (ll3.lt.1) stop 'alat(3) too small'

! determine number of threads
        npr=1
!$omp parallel private(iam)  shared (npr)
!$       iam=omp_get_thread_num()
!$       if (iam.eq.0) npr=omp_get_num_threads()
!$omp end parallel

! linear scaling calculation of verlet list

     if (npr.le.1) then !serial if too few processors to gain by parallelizing

! set ncx for serial case, ncx for parallel case set below
        ncx=8
1234    ncx=ncx*2
        allocate(icell(0:ncx,-1:ll1,-1:ll2,-1:ll3))
        do 984,l3=-1,ll3
        do 984,l2=-1,ll2
        do 984,l1=-1,ll1
984     icell(0,l1,l2,l3)=0
        rlc1i=ll1/alat(1)
        rlc2i=ll2/alat(2)
        rlc3i=ll3/alat(3)

        do 983,iat=1,nat
        rxyz0(1,iat)=modulo(modulo(rxyz0(1,iat),alat(1)),alat(1))
        rxyz0(2,iat)=modulo(modulo(rxyz0(2,iat),alat(2)),alat(2))
        rxyz0(3,iat)=modulo(modulo(rxyz0(3,iat),alat(3)),alat(3))
        l1=int(rxyz0(1,iat)*rlc1i)
        l2=int(rxyz0(2,iat)*rlc2i)
        l3=int(rxyz0(3,iat)*rlc3i)

        ii=icell(0,l1,l2,l3)
        ii=ii+1
        icell(0,l1,l2,l3)=ii
        if (ii.gt.ncx) then
        write(10,*) count,'NCX too small',ncx
        deallocate(icell)
        goto 1234
        endif
        icell(ii,l1,l2,l3)=iat
983     continue

     else  ! parallel case

! periodization of particles can be done in parallel
!$omp parallel do shared (alat,nat,rxyz0) private(iat)
        do 5983,iat=1,nat
        rxyz0(1,iat)=modulo(modulo(rxyz0(1,iat),alat(1)),alat(1))
        rxyz0(2,iat)=modulo(modulo(rxyz0(2,iat),alat(2)),alat(2))
        rxyz0(3,iat)=modulo(modulo(rxyz0(3,iat),alat(3)),alat(3))
5983    continue
!$omp end parallel do

! assignment to cell is done serially
! set ncx for parallel case, ncx for serial case set above
        ncx=8
4321    ncx=ncx*2
        allocate(icell(0:ncx,-1:ll1,-1:ll2,-1:ll3))
        do 3984,l3=-1,ll3
        do 3984,l2=-1,ll2
        do 3984,l1=-1,ll1
3984    icell(0,l1,l2,l3)=0

        rlc1i=ll1/alat(1)
        rlc2i=ll2/alat(2)
        rlc3i=ll3/alat(3)

        do 6983,iat=1,nat
        l1=int(rxyz0(1,iat)*rlc1i)
        l2=int(rxyz0(2,iat)*rlc2i)
        l3=int(rxyz0(3,iat)*rlc3i)
        ii=icell(0,l1,l2,l3)
        ii=ii+1
        icell(0,l1,l2,l3)=ii
        if (ii.gt.ncx) then
        write(10,*) count,'NCX too small',ncx
        deallocate(icell)
        goto 4321
        endif
        icell(ii,l1,l2,l3)=iat
6983    continue

    endif


! duplicate all atoms within boundary layer
        laymx=ncx*(2*ll1*ll2+2*ll1*ll3+2*ll2*ll3+4*ll1+4*ll2+4*ll3+8)
        nn=nat+laymx
        allocate(rxyz(3,nn),lay(nn))
        do  iat=1,nat
        lay(iat)=iat
        rxyz(1,iat)=rxyz0(1,iat)
        rxyz(2,iat)=rxyz0(2,iat)
        rxyz(3,iat)=rxyz0(3,iat)
        enddo
        il=nat
! xy plane
        do l2=0,ll2-1
        do l1=0,ll1-1

        in=icell(0,l1,l2,0)
        icell(0,l1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,l1,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,l2,ll3-1)
        icell(0,l1,l2,-1)=in
        do ii=1,in
        i=icell(ii,l1,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo
        enddo


! yz plane
        do l3=0,ll3-1
        do l2=0,ll2-1

        in=icell(0,0,l2,l3)
        icell(0,ll1,l2,l3)=in
        do ii=1,in
        i=icell(ii,0,l2,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,l2,l3)
        icell(0,-1,l2,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo
        enddo


! xz plane
        do l3=0,ll3-1
        do l1=0,ll1-1

        in=icell(0,l1,0,l3)
        icell(0,l1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,l1,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,l1,ll2-1,l3)
        icell(0,l1,-1,l3)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo
        enddo


! x axis
        do l1=0,ll1-1

        in=icell(0,l1,0,0)
        icell(0,l1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,l1,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,0,ll3-1)
        icell(0,l1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,l1,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,l1,ll2-1,0)
        icell(0,l1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,ll2-1,ll3-1)
        icell(0,l1,-1,-1)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo


! y axis
        do l2=0,ll2-1

        in=icell(0,0,l2,0)
        icell(0,ll1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,0,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,l2,ll3-1)
        icell(0,ll1,l2,-1)=in
        do ii=1,in
        i=icell(ii,0,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,ll1-1,l2,0)
        icell(0,-1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,l2,ll3-1)
        icell(0,-1,l2,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo


! z axis
        do l3=0,ll3-1

        in=icell(0,0,0,l3)
        icell(0,ll1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,0,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,0,l3)
        icell(0,-1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,0,ll2-1,l3)
        icell(0,ll1,-1,l3)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,ll2-1,l3)
        icell(0,-1,-1,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo


! corners
        in=icell(0,0,0,0)
        icell(0,ll1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,0,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,0,0)
        icell(0,-1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,ll2-1,0)
        icell(0,ll1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,ll2-1,0)
        icell(0,-1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,0,ll3-1)
        icell(0,ll1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,0,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,ll1-1,0,ll3-1)
        icell(0,-1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,0,ll2-1,ll3-1)
        icell(0,ll1,-1,-1)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,ll1-1,ll2-1,ll3-1)
        icell(0,-1,-1,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        allocate(lsta(2,nat))
        nnbrx=24
2345    nnbrx=3*nnbrx/2
        allocate(lstb(nnbrx*nat),rel(5,nnbrx*nat))

        indlstx=0

!$omp parallel  &
!$omp private(iat,cut2,iam,ii,indlst,l1,l2,l3,myspace,npr) &
!$omp shared (indlstx,nat,nn,nnbrx,ncx,ll1,ll2,ll3,icell,lsta,lstb,lay, &
!$omp rel,rxyz,cut,myspaceout)


        npr=1
!$       npr=omp_get_num_threads()
        iam=0
!$       iam=omp_get_thread_num()

        cut2=cut**2
! assign contiguous portions of the arrays lstb and rel to the threads
        myspace=(nat*nnbrx)/npr
        if (iam.eq.0) myspaceout=myspace
! Verlet list, relative positions
        indlst=0
      do 6000,l3=0,ll3-1
      do 6000,l2=0,ll2-1
      do 6000,l1=0,ll1-1
      do 6600,ii=1,icell(0,l1,l2,l3)
        iat=icell(ii,l1,l2,l3)
        if ( ((iat-1)*npr)/nat .eq. iam) then
!       write(6,*) 'sublstiat:iam,iat',iam,iat
        lsta(1,iat)=iam*myspace+indlst+1
        call sublstiat_l(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, &
             rxyz,icell,lstb(iam*myspace+1),lay,rel(1,iam*myspace+1),cut2,indlst)
        lsta(2,iat)=iam*myspace+indlst
!        write(6,'(a,4(x,i3),100(x,i2))') &
!               'iam,iat,lsta',iam,iat,lsta(1,iat),lsta(2,iat), &
!                    (lstb(j),j=lsta(1,iat),lsta(2,iat))
        endif

6600    continue
6000    continue
!$omp critical
        indlstx=max(indlstx,indlst)
!$omp end critical
!$omp end parallel

           if (indlstx.ge.myspaceout) then
               write(10,*) count,'NNBRX too  small', nnbrx
               deallocate(lstb,rel)
               goto 2345
           endif

!$omp parallel  &
!$omp private(iam,npr,iat,iat1,iat2,lot,istop,tcoord,tcoord2, &
!$omp tener,tener2,txyz,f2ij,f3ij,f3ik,npjx,npjkx) &
!$omp shared (nat,nnbrx,lsta,lstb,rel,ener,ener2,fxyz,coord,coord2,istopg)

        npr=1
!$       npr=omp_get_num_threads()
        iam=0
!$       iam=omp_get_thread_num()

        npjx=300 ; npjkx=6000
        istopg=0

        if (npr.ne.1) then
! PARALLEL CASE
! create temporary private scalars for reduction sum on energies and
!        temporary private array for reduction sum on forces
!$omp critical
        allocate(txyz(3,nat),f2ij(3,npjx),f3ij(3,npjkx),f3ik(3,npjkx))
!$omp end critical
        if (iam.eq.0) then
        ener=0.0_gp
        ener2=0.0_gp
        coord=0.0_gp
        coord2=0.0_gp
        endif
!$omp do
        do 121,iat=1,nat
        fxyz(1,iat)=0.0_gp
        fxyz(2,iat)=0.0_gp
121     fxyz(3,iat)=0.0_gp
!$omp barrier

! Each thread treats at most lot atoms
        lot=int(float(nat)/float(npr)+.999999999999d0)
        iat1=iam*lot+1
        iat2=min((iam+1)*lot,nat)
!       write(6,*) 'subfeniat:iat1,iat2,iam',iat1,iat2,iam
        call subfeniat_l(iat1,iat2,nat,lsta,lstb,rel,tener,tener2, &
               tcoord,tcoord2,nnbrx,txyz,f2ij,npjx,f3ij,npjkx,f3ik,istop)
!$omp critical
        ener=ener+tener
        ener2=ener2+tener2
        coord=coord+tcoord
        coord2=coord2+tcoord2
        istopg=istopg+istop
        do 8093,iat=1,nat
        fxyz(1,iat)=fxyz(1,iat)+txyz(1,iat)
        fxyz(2,iat)=fxyz(2,iat)+txyz(2,iat)
        fxyz(3,iat)=fxyz(3,iat)+txyz(3,iat)
8093    continue
        deallocate(txyz,f2ij,f3ij,f3ik)
!$omp end critical

        else
! SERIAL CASE
        iat1=1
        iat2=nat
        allocate(f2ij(3,npjx),f3ij(3,npjkx),f3ik(3,npjkx))
        call subfeniat_l(iat1,iat2,nat,lsta,lstb,rel,ener,ener2, &
               coord,coord2,nnbrx,fxyz,f2ij,npjx,f3ij,npjkx,f3ik,istopg)
        deallocate(f2ij,f3ij,f3ik)

        endif
!$omp end parallel

!         write(6,*) 'ener,norm force', &
!                    ener,DNRM2(3*nat,fxyz,1)
        if (istopg.gt.0) stop 'DIMENSION ERROR (see WARNING above)'
        ener_var=ener2/nat-(ener/nat)**2
        coord=coord/nat
        coord_var=coord2/nat-coord**2

        deallocate(rxyz,icell,lay,lsta,lstb,rel)

        end subroutine


        subroutine subfeniat_l(iat1,iat2,nat,lsta,lstb,rel,tener,tener2, &
               tcoord,tcoord2,nnbrx,txyz,f2ij,npjx,f3ij,npjkx,f3ik,istop)
        use module_defs, only: gp
! for a subset of atoms iat1 to iat2 the routine calculates the (partial) forces
! txyz acting on these atoms as well as on the atoms (jat, kat) interacting
! with them and their contribution to the energy (tener).
! In addition the coordination number tcoord and the second moment of the
! local energy tener2 and coordination number tcoord2 are returned
        implicit real(gp) (a-h,o-z)
        dimension lsta(2,nat),lstb(nnbrx*nat),rel(5,nnbrx*nat),txyz(3,nat)
        dimension f2ij(3,npjx),f3ij(3,npjkx),f3ik(3,npjkx)
       real(gp) :: tmin_phi= 0.1500000d+01
       real(gp) :: tmax_phi= 0.4500000d+01
       real(gp) :: hi_phi= 3.00000000000d0
       real(gp) :: hsixth_phi=5.55555555555556D-002
       real(gp) :: h2sixth_phi=1.85185185185185D-002
       real(gp), dimension(0:9) :: cof_phi  = &
                    (/ 0.69299400000000d+01, -0.43995000000000d+00, &
                      -0.17012300000000d+01, -0.16247300000000d+01, &
                      -0.99696000000000d+00, -0.27391000000000d+00, &
                      -0.24990000000000d-01, -0.17840000000000d-01, &
                      -0.96100000000000d-02,  0.00000000000000d+00 /)
       real(gp), dimension(0:9) :: dof_phi  = &
                    (/ 0.16533229480429d+03,  0.39415410391417d+02, &
                       0.68710036300407d+01,  0.53406950884203d+01, &
                       0.15347960162782d+01, -0.63347591535331d+01, &
                      -0.17987794021458d+01,  0.47429676211617d+00, &
                      -0.40087646318907d-01, -0.23942617684055d+00 /)
       real(gp) :: tmin_rho= 0.1500000d+01
       real(gp) :: tmax_rho= 0.3500000d+01
       real(gp) :: hi_rho= 5.00000000000d0
       real(gp) :: hsixth_rho=3.33333333333333D-002
       real(gp) :: h2sixth_rho=6.66666666666667D-003
       real(gp), dimension(0:10) :: cof_rho =  &
                    (/ 0.13747000000000d+00, -0.14831000000000d+00, &
                      -0.55972000000000d+00, -0.73110000000000d+00, &
                      -0.76283000000000d+00, -0.72918000000000d+00, &
                      -0.66620000000000d+00, -0.57328000000000d+00, &
                      -0.40690000000000d+00, -0.16662000000000d+00, &
                       0.00000000000000d+00 /)
       real(gp), dimension(0:10) :: dof_rho =  &
                   (/ -0.32275496741918d+01, -0.64119006516165d+01, &
                       0.10030652280658d+02,  0.22937915289857d+01, &
                       0.17416816033995d+01,  0.54648205741626d+00, &
                       0.47189016693543d+00,  0.20569572748420d+01, &
                       0.23192807336964d+01, -0.24908020962757d+00, &
                      -0.12371959895186d+02 /)
       real(gp) :: tmin_fff= 0.1500000d+01
       real(gp) :: tmax_fff= 0.3500000d+01
       real(gp) :: hi_fff= 4.50000000000d0
       real(gp) :: hsixth_fff=3.70370370370370D-002
       real(gp) :: h2sixth_fff=8.23045267489712D-003
       real(gp), dimension(0:9) :: cof_fff  = &
                    (/ 0.12503100000000d+01,  0.86821000000000d+00, &
                       0.60846000000000d+00,  0.48756000000000d+00, &
                       0.44163000000000d+00,  0.37610000000000d+00, &
                       0.27145000000000d+00,  0.14814000000000d+00, &
                       0.48550000000000d-01,  0.00000000000000d+00 /)
       real(gp), dimension(0:9) :: dof_fff  = &
                    (/ 0.27904652711432d+02, -0.45230754228635d+01, &
                       0.50531739800222d+01,  0.11806545027747d+01, &
                      -0.66693699112098d+00, -0.89430653829079d+00, &
                      -0.50891685571587d+00,  0.66278396115427d+00, &
                       0.73976101109878d+00,  0.25795319944506d+01 /)
       real(gp) :: tmin_uuu= -0.1770930d+01
       real(gp) :: tmax_uuu= 0.7908520d+01
       real(gp) :: hi_uuu= 0.723181585730594d0
       real(gp) :: hsixth_uuu=0.230463095238095d0
       real(gp) :: h2sixth_uuu=0.318679429600340d0
       real(gp), dimension(0:7) :: cof_uuu  = &
                   (/ -0.10749300000000d+01, -0.20045000000000d+00, &
                       0.41422000000000d+00,  0.87939000000000d+00, &
                       0.12668900000000d+01,  0.16299800000000d+01, &
                       0.19773800000000d+01,  0.23961800000000d+01 /)
       real(gp), dimension(0:7) :: dof_uuu  = &
                   (/ -0.14827125747284d+00, -0.14922155328475d+00, &
                      -0.70113224223509d-01, -0.39449020349230d-01, &
                      -0.15815242579643d-01,  0.26112640061855d-01, &
                      -0.13786974745095d+00,  0.74941595372657d+00 /)
       real(gp) :: tmin_ggg= -0.1000000d+01
       real(gp) :: tmax_ggg= 0.8001400d+00
       real(gp) :: hi_ggg= 3.88858644327663d0
       real(gp) :: hsixth_ggg=4.28604761904762D-002
       real(gp) :: h2sixth_ggg=1.10221225156463D-002
       real(gp), dimension(0:7) :: cof_ggg  = &
                    (/ 0.52541600000000d+01,  0.23591500000000d+01, &
                       0.11959500000000d+01,  0.12299500000000d+01, &
                       0.20356500000000d+01,  0.34247400000000d+01, &
                       0.49485900000000d+01,  0.56179900000000d+01 /)
       real(gp), dimension(0:7) :: dof_ggg  = &
                    (/ 0.15826876132396d+02,  0.31176239377907d+02, &
                       0.16589446539683d+02,  0.11083892500520d+02, &
                       0.90887216383860d+01,  0.54902279653967d+01, &
                      -0.18823313223755d+02, -0.77183416481005d+01 /)

! initialize temporary private scalars for reduction sum on energies and
! private workarray txyz for forces forces
        tener=0.0_gp
        tener2=0.0_gp
        tcoord=0.0_gp
        tcoord2=0.0_gp
        istop=0
        do 121,iat=1,nat
        txyz(1,iat)=0.0_gp
        txyz(2,iat)=0.0_gp
121     txyz(3,iat)=0.0_gp


! calculation of forces, energy

        do 1000,iat=iat1,iat2

        dens2=0.0_gp
        dens3=0.0_gp
        jcnt=0
        jkcnt=0
        coord_iat=0.0_gp
        ener_iat=0.0_gp
        do 2000,jbr=lsta(1,iat),lsta(2,iat)
        jat=lstb(jbr)
        jcnt=jcnt+1
        if (jcnt.gt.npjx) then
            write(6,*) 'WARNING: enlarge npjx'
            istop=1
            return
        endif

        fxij=rel(1,jbr)
        fyij=rel(2,jbr)
        fzij=rel(3,jbr)
        rij=rel(4,jbr)
        sij=rel(5,jbr)

! coordination number calculated with soft cutoff between first
! nearest neighbor and midpoint of first and second nearest neighbor
        if (rij.le.2.36d0) then
        coord_iat=coord_iat+1.0_gp
        else if (rij.ge.3.12d0) then
        else
        xarg=(rij-2.36d0)*(1.0_gp/(3.12d0-2.36d0))
        coord_iat=coord_iat+(2*xarg+1.0_gp)*(xarg-1.0_gp)**2
        endif

! pairpotential term
        call splint(cof_phi,dof_phi,tmin_phi,tmax_phi, &
           hsixth_phi,h2sixth_phi,hi_phi,10,rij,e_phi,ep_phi)
        ener_iat=ener_iat+(e_phi*.5d0)
        txyz(1,iat)=txyz(1,iat)-fxij*(ep_phi*.5d0)
        txyz(2,iat)=txyz(2,iat)-fyij*(ep_phi*.5d0)
        txyz(3,iat)=txyz(3,iat)-fzij*(ep_phi*.5d0)
        txyz(1,jat)=txyz(1,jat)+fxij*(ep_phi*.5d0)
        txyz(2,jat)=txyz(2,jat)+fyij*(ep_phi*.5d0)
        txyz(3,jat)=txyz(3,jat)+fzij*(ep_phi*.5d0)

! 2 body embedding term
        call splint(cof_rho,dof_rho,tmin_rho,tmax_rho, &
             hsixth_rho,h2sixth_rho,hi_rho,11,rij,rho,rhop)
        dens2=dens2+rho
        f2ij(1,jcnt)=fxij*rhop
        f2ij(2,jcnt)=fyij*rhop
        f2ij(3,jcnt)=fzij*rhop

! 3 body embedding term
        call splint(cof_fff,dof_fff,tmin_fff,tmax_fff, &
             hsixth_fff,h2sixth_fff,hi_fff,10,rij,fij,fijp)

        do 3000,kbr=lsta(1,iat),lsta(2,iat)
        kat=lstb(kbr)
        if (kat.lt.jat) then
        jkcnt=jkcnt+1
        if (jkcnt.gt.npjkx) then
            write(6,*) 'WARNING: enlarge npjkx',npjkx
            istop=1
            return
        endif

! begin unoptimized original version:
!        fxik=rel(1,kbr)
!        fyik=rel(2,kbr)
!        fzik=rel(3,kbr)
!        rik=rel(4,kbr)
!        sik=rel(5,kbr)
!
!        call splint(cof_fff,dof_fff,tmin_fff,tmax_fff, &
!             hsixth_fff,h2sixth_fff,hi_fff,10,rik,fik,fikp)
!        costheta=fxij*fxik+fyij*fyik+fzij*fzik
!        call splint(cof_ggg,dof_ggg,tmin_ggg,tmax_ggg, &
!             hsixth_ggg,h2sixth_ggg,hi_ggg,8,costheta,gjik,gjikp)
! end unoptimized original version:

! begin optimized version
        rik=rel(4,kbr)
      if (rik.gt.tmax_fff) then
        fikp=0.0_gp ; fik=0.0_gp
        gjik=0.0_gp ;  gjikp=0.0_gp ; sik=0.0_gp
        costheta=0.0_gp ; fxik=0.0_gp ; fyik=0.0_gp ; fzik=0.0_gp
      else if (rik.lt.tmin_fff) then
        fxik=rel(1,kbr)
        fyik=rel(2,kbr)
        fzik=rel(3,kbr)
        costheta=fxij*fxik+fyij*fyik+fzij*fzik
        sik=rel(5,kbr)
        fikp=hi_fff*(cof_fff(1)-cof_fff(0)) -  &
             ( dof_fff(1)+2.0_gp*dof_fff(0) )*hsixth_fff
        fik=cof_fff(0) + (rik-tmin_fff)*fikp
        tt_ggg=(costheta-tmin_ggg)*hi_ggg
        if (costheta.gt.tmax_ggg) then
           gjikp=hi_ggg*(cof_ggg(8-1)-cof_ggg(8-2)) + &
                ( 2.0_gp*dof_ggg(8-1)+dof_ggg(8-2) )*hsixth_ggg
                 gjik=cof_ggg(8-1) + (costheta-tmax_ggg)*gjikp
        else
           klo_ggg=int(tt_ggg,kind=4)
           khi_ggg=klo_ggg+1
           cof_ggg_klo=cof_ggg(klo_ggg)
           dof_ggg_klo=dof_ggg(klo_ggg)
           b_ggg=tt_ggg-klo_ggg
           a_ggg=1.0_gp-b_ggg
           cof_ggg_khi=cof_ggg(khi_ggg)
           dof_ggg_khi=dof_ggg(khi_ggg)
           b2_ggg=b_ggg*b_ggg
           gjik=a_ggg*cof_ggg_klo
           gjikp=cof_ggg_khi-cof_ggg_klo
           a2_ggg=a_ggg*a_ggg
           cof1_ggg=a2_ggg-1.0_gp
           cof2_ggg=b2_ggg-1.0_gp
           gjik=gjik+b_ggg*cof_ggg_khi
           gjikp=hi_ggg*gjikp
           cof3_ggg=3.0_gp*b2_ggg
           cof4_ggg=3.0_gp*a2_ggg
           cof1_ggg=a_ggg*cof1_ggg
           cof2_ggg=b_ggg*cof2_ggg
           cof3_ggg=cof3_ggg-1.0_gp
           cof4_ggg=cof4_ggg-1.0_gp
           yt1_ggg=cof1_ggg*dof_ggg_klo
           yt2_ggg=cof2_ggg*dof_ggg_khi
           ypt1_ggg=cof3_ggg*dof_ggg_khi
           ypt2_ggg=cof4_ggg*dof_ggg_klo
           gjik=gjik + (yt1_ggg+yt2_ggg)*h2sixth_ggg
           gjikp=gjikp + ( ypt1_ggg - ypt2_ggg )*hsixth_ggg
        endif
      else
        fxik=rel(1,kbr)
        tt_fff=rik-tmin_fff
        costheta=fxij*fxik
        fyik=rel(2,kbr)
        tt_fff=tt_fff*hi_fff
        costheta=costheta+fyij*fyik
        fzik=rel(3,kbr)
        klo_fff=int(tt_fff,kind=4)
        costheta=costheta+fzij*fzik
        sik=rel(5,kbr)
        tt_ggg=(costheta-tmin_ggg)*hi_ggg
        if (costheta.gt.tmax_ggg) then
          gjikp=hi_ggg*(cof_ggg(8-1)-cof_ggg(8-2)) + &
                ( 2.0_gp*dof_ggg(8-1)+dof_ggg(8-2) )*hsixth_ggg
                gjik=cof_ggg(8-1) + (costheta-tmax_ggg)*gjikp
          khi_fff=klo_fff+1
          cof_fff_klo=cof_fff(klo_fff)
          dof_fff_klo=dof_fff(klo_fff)
          b_fff=tt_fff-klo_fff
          a_fff=1.0_gp-b_fff
          cof_fff_khi=cof_fff(khi_fff)
          dof_fff_khi=dof_fff(khi_fff)
          b2_fff=b_fff*b_fff
          fik=a_fff*cof_fff_klo
          fikp=cof_fff_khi-cof_fff_klo
          a2_fff=a_fff*a_fff
          cof1_fff=a2_fff-1.0_gp
          cof2_fff=b2_fff-1.0_gp
          fik=fik+b_fff*cof_fff_khi
          fikp=hi_fff*fikp
          cof3_fff=3.0_gp*b2_fff
          cof4_fff=3.0_gp*a2_fff
          cof1_fff=a_fff*cof1_fff
          cof2_fff=b_fff*cof2_fff
          cof3_fff=cof3_fff-1.0_gp
          cof4_fff=cof4_fff-1.0_gp
          yt1_fff=cof1_fff*dof_fff_klo
          yt2_fff=cof2_fff*dof_fff_khi
          ypt1_fff=cof3_fff*dof_fff_khi
          ypt2_fff=cof4_fff*dof_fff_klo
          fik=fik + (yt1_fff+yt2_fff)*h2sixth_fff
          fikp=fikp + ( ypt1_fff - ypt2_fff )*hsixth_fff
         else
              klo_ggg=int(tt_ggg,kind=4)
              khi_ggg=klo_ggg+1
           khi_fff=klo_fff+1
              cof_ggg_klo=cof_ggg(klo_ggg)
           cof_fff_klo=cof_fff(klo_fff)
              dof_ggg_klo=dof_ggg(klo_ggg)
           dof_fff_klo=dof_fff(klo_fff)
              b_ggg=tt_ggg-klo_ggg
           b_fff=tt_fff-klo_fff
              a_ggg=1.0_gp-b_ggg
           a_fff=1.0_gp-b_fff
              cof_ggg_khi=cof_ggg(khi_ggg)
           cof_fff_khi=cof_fff(khi_fff)
              dof_ggg_khi=dof_ggg(khi_ggg)
           dof_fff_khi=dof_fff(khi_fff)
              b2_ggg=b_ggg*b_ggg
           b2_fff=b_fff*b_fff
              gjik=a_ggg*cof_ggg_klo
           fik=a_fff*cof_fff_klo
              gjikp=cof_ggg_khi-cof_ggg_klo
           fikp=cof_fff_khi-cof_fff_klo
              a2_ggg=a_ggg*a_ggg
           a2_fff=a_fff*a_fff
              cof1_ggg=a2_ggg-1.0_gp
           cof1_fff=a2_fff-1.0_gp
              cof2_ggg=b2_ggg-1.0_gp
           cof2_fff=b2_fff-1.0_gp
              gjik=gjik+b_ggg*cof_ggg_khi
           fik=fik+b_fff*cof_fff_khi
              gjikp=hi_ggg*gjikp
           fikp=hi_fff*fikp
              cof3_ggg=3.0_gp*b2_ggg
           cof3_fff=3.0_gp*b2_fff
              cof4_ggg=3.0_gp*a2_ggg
           cof4_fff=3.0_gp*a2_fff
              cof1_ggg=a_ggg*cof1_ggg
           cof1_fff=a_fff*cof1_fff
              cof2_ggg=b_ggg*cof2_ggg
           cof2_fff=b_fff*cof2_fff
              cof3_ggg=cof3_ggg-1.0_gp
           cof3_fff=cof3_fff-1.0_gp
              cof4_ggg=cof4_ggg-1.0_gp
           cof4_fff=cof4_fff-1.0_gp
              yt1_ggg=cof1_ggg*dof_ggg_klo
           yt1_fff=cof1_fff*dof_fff_klo
              yt2_ggg=cof2_ggg*dof_ggg_khi
           yt2_fff=cof2_fff*dof_fff_khi
              ypt1_ggg=cof3_ggg*dof_ggg_khi
           ypt1_fff=cof3_fff*dof_fff_khi
              ypt2_ggg=cof4_ggg*dof_ggg_klo
           ypt2_fff=cof4_fff*dof_fff_klo
              gjik=gjik + (yt1_ggg+yt2_ggg)*h2sixth_ggg
           fik=fik + (yt1_fff+yt2_fff)*h2sixth_fff
              gjikp=gjikp + ( ypt1_ggg - ypt2_ggg )*hsixth_ggg
           fikp=fikp + ( ypt1_fff - ypt2_fff )*hsixth_fff
         endif
      endif
! end optimized version

        tt=fij*fik
        dens3=dens3+tt*gjik

        t1=fijp*fik*gjik
        t2=sij*(tt*gjikp)
        f3ij(1,jkcnt)=fxij*t1 + (fxik-fxij*costheta)*t2
        f3ij(2,jkcnt)=fyij*t1 + (fyik-fyij*costheta)*t2
        f3ij(3,jkcnt)=fzij*t1 + (fzik-fzij*costheta)*t2

        t3=fikp*fij*gjik
        t4=sik*(tt*gjikp)
        f3ik(1,jkcnt)=fxik*t3 + (fxij-fxik*costheta)*t4
        f3ik(2,jkcnt)=fyik*t3 + (fyij-fyik*costheta)*t4
        f3ik(3,jkcnt)=fzik*t3 + (fzij-fzik*costheta)*t4
        endif
3000        continue
2000        continue
        dens=dens2+dens3
        call splint(cof_uuu,dof_uuu,tmin_uuu,tmax_uuu, &
             hsixth_uuu,h2sixth_uuu,hi_uuu,8,dens,e_uuu,ep_uuu)
        ener_iat=ener_iat+e_uuu

! Only now ep_uu is known and the forces can be calculated, lets loop again
        jcnt=0
        jkcnt=0
        do 2200,jbr=lsta(1,iat),lsta(2,iat)
        jat=lstb(jbr)
        jcnt=jcnt+1
        txyz(1,iat)=txyz(1,iat)-ep_uuu*f2ij(1,jcnt)
        txyz(2,iat)=txyz(2,iat)-ep_uuu*f2ij(2,jcnt)
        txyz(3,iat)=txyz(3,iat)-ep_uuu*f2ij(3,jcnt)
        txyz(1,jat)=txyz(1,jat)+ep_uuu*f2ij(1,jcnt)
        txyz(2,jat)=txyz(2,jat)+ep_uuu*f2ij(2,jcnt)
        txyz(3,jat)=txyz(3,jat)+ep_uuu*f2ij(3,jcnt)

! 3 body embedding term
        do 3300,kbr=lsta(1,iat),lsta(2,iat)
        kat=lstb(kbr)
        if (kat.lt.jat) then
        jkcnt=jkcnt+1

        txyz(1,iat)=txyz(1,iat)-ep_uuu*(f3ij(1,jkcnt)+f3ik(1,jkcnt))
        txyz(2,iat)=txyz(2,iat)-ep_uuu*(f3ij(2,jkcnt)+f3ik(2,jkcnt))
        txyz(3,iat)=txyz(3,iat)-ep_uuu*(f3ij(3,jkcnt)+f3ik(3,jkcnt))
        txyz(1,jat)=txyz(1,jat)+ep_uuu*f3ij(1,jkcnt)
        txyz(2,jat)=txyz(2,jat)+ep_uuu*f3ij(2,jkcnt)
        txyz(3,jat)=txyz(3,jat)+ep_uuu*f3ij(3,jkcnt)
        txyz(1,kat)=txyz(1,kat)+ep_uuu*f3ik(1,jkcnt)
        txyz(2,kat)=txyz(2,kat)+ep_uuu*f3ik(2,jkcnt)
        txyz(3,kat)=txyz(3,kat)+ep_uuu*f3ik(3,jkcnt)


        endif
3300        continue
2200        continue

!        write(6,'(a,i4,x,e19.12,x,e10.3)') 'iat,ener_iat,coord_iat', &
!                                       iat,ener_iat,coord_iat
        tener=tener+ener_iat
        tener2=tener2+ener_iat**2
        tcoord=tcoord+coord_iat
        tcoord2=tcoord2+coord_iat**2

1000        continue

        return
        end subroutine


        subroutine sublstiat_l(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, &
                   rxyz,icell,lstb,lay,rel,cut2,indlst)
        use module_defs, only: gp
! finds the neighbours of atom iat (specified by lsta and lstb) and and
! the relative position rel of iat with respect to these neighbours
        implicit real(gp) (a-h,o-z)
        dimension rxyz(3,nn),lay(nn),icell(0:ncx,-1:ll1,-1:ll2,-1:ll3), &
                  lstb(0:myspace-1),rel(5,0:myspace-1)

        do 6363,k3=l3-1,l3+1
        do 6363,k2=l2-1,l2+1
        do 6363,k1=l1-1,l1+1
        do 6363,jj=1,icell(0,k1,k2,k3)
          jat=icell(jj,k1,k2,k3)
          if (jat.eq.iat) goto 6363
          xrel= rxyz(1,iat)-rxyz(1,jat)
          yrel= rxyz(2,iat)-rxyz(2,jat)
          zrel= rxyz(3,iat)-rxyz(3,jat)
          rr2=xrel**2 + yrel**2 + zrel**2
          if ( rr2 .le. cut2 ) then
           indlst=min(indlst,myspace-1)
           lstb(indlst)=lay(jat)
!        write(6,*) 'iat,indlst,lay(jat)',iat,indlst,lay(jat)
           tt=sqrt(rr2)
           tti=1.0_gp/tt
           rel(1,indlst)=xrel*tti
           rel(2,indlst)=yrel*tti
           rel(3,indlst)=zrel*tti
           rel(4,indlst)=tt
           rel(5,indlst)=tti
           indlst= indlst+1
          endif
6363        continue

        return
        end subroutine


        subroutine splint(ya,y2a,tmin,tmax,hsixth,h2sixth,hi,n,x,y,yp)
        use module_defs, only: gp
        implicit none
        !implicit real(gp) (a-h,o-z)
        !Arguments
        integer, intent(in) :: n
        real(gp), dimension(0:n-1), intent(in) :: y2a,ya
        real(gp), intent(in) :: tmin,tmax,hsixth,h2sixth,hi,x
        real(gp), intent(out) :: y,yp
        !Local variables
        integer :: khi,klo
        real(gp) :: tt
        real(gp) :: a,b,ya_khi,ya_klo,b2,a2,y2a_khi,y2a_klo
        real(gp) :: cof1,cof2,cof3,cof4,yt1,yt2,ypt1,ypt2
        ! interpolate if the argument is outside the cubic spline interval [tmin,tmax]
        tt=(x-tmin)*hi
        if (x.lt.tmin) then
          yp=hi*(ya(1)-ya(0)) -  &
          ( y2a(1)+2.0_gp*y2a(0) )*hsixth
          y=ya(0) + (x-tmin)*yp
        else if (x.gt.tmax) then
          yp=hi*(ya(n-1)-ya(n-2)) +  &
          ( 2.0_gp*y2a(n-1)+y2a(n-2) )*hsixth
          y=ya(n-1) + (x-tmax)*yp
        else
          ! otherwise evaluate cubic spline
          klo=int(tt,kind=4)
          khi=klo+1
          ya_klo=ya(klo)
          y2a_klo=y2a(klo)
          b=tt-klo
          a=1.0_gp-b
          ya_khi=ya(khi)
          y2a_khi=y2a(khi)
          b2=b*b
          y=a*ya_klo
          yp=ya_khi-ya_klo
          a2=a*a
          cof1=a2-1.0_gp
          cof2=b2-1.0_gp
          y=y+b*ya_khi
          yp=hi*yp
          cof3=3.0_gp*b2
          cof4=3.0_gp*a2
          cof1=a*cof1
          cof2=b*cof2
          cof3=cof3-1.0_gp
          cof4=cof4-1.0_gp
          yt1=cof1*y2a_klo
          yt2=cof2*y2a_khi
          ypt1=cof3*y2a_khi
          ypt2=cof4*y2a_klo
          y=y + (yt1+yt2)*h2sixth
          yp=yp + ( ypt1 - ypt2 )*hsixth
        endif
      end subroutine


end module
