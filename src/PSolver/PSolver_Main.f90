!!****f* H_potential/H_potential
!! FUNCTION
!!    Calculate the Hartree potential by solving Poisson equation 
!!       $\nabla^2 V(x,y,z)=-4 \pi \rho(x,y,z)$
!!    from a given $\rho$, 
!!    for different boundary conditions an for different data distributions.
!!    Following the boundary conditions, it applies the Poisson Kernel previously calculated.
!!
!! COPYRIGHT
!!    Copyright (C) 2002-2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SYNOPSIS
!!    geocode  Indicates the boundary conditions (BC) of the problem:
!!            'F' free BC, isolated systems.
!!                The program calculates the solution as if the given density is
!!                "alone" in R^3 space.
!!            'S' surface BC, isolated in y direction, periodic in xz plane                
!!                The given density is supposed to be periodic in the xz plane,
!!                so the dimensions in these direction mus be compatible with the FFT
!!                Beware of the fact that the isolated direction is y!
!!            'P' periodic BC.
!!                The density is supposed to be periodic in all the three directions,
!!                then all the dimensions must be compatible with the FFT.
!!                No need for setting up the kernel.
!!    datacode Indicates the distribution of the data of the input/output array:
!!            'G' global data. Each process has the whole array of the density 
!!                which will be overwritten with the whole array of the potential
!!            'D' distributed data. Each process has only the needed part of the density
!!                and of the potential. The data distribution is such that each processor
!!                has the xy planes needed for the calculation AND for the evaluation of the 
!!                gradient, needed for XC part, and for the White-Bird correction, which
!!                may lead up to 8 planes more on each side. Due to this fact, the information
!!                between the processors may overlap.
!!    nproc       number of processors
!!    iproc       label of the process,from 0 to nproc-1
!!    n01,n02,n03 global dimension in the three directions. They are the same no matter if the 
!!                datacode is in 'G' or in 'D' position.
!!    hx,hy,hz    grid spacings. For the isolated BC case for the moment they are supposed to 
!!                be equal in the three directions
!!    rhopot      main input/output array.
!!                On input, it represents the density values on the grid points
!!                On output, it is the Hartree potential
!!    karray      kernel of the poisson equation. It is provided in distributed case, with
!!                dimensions that are related to the output of the PS_dim4allocation routine
!!                it MUST be created by following the same geocode as the Poisson Solver.
!!    pot_ion     additional external potential that is added to the output, 
!!                when the XC parameter ixc/=0 and sumpion=.true.
!!                When sumpion=.true., it is always provided in the distributed form,
!!                clearly without the overlapping terms which are needed only for the XC part
!!    eh          Hartree energy
!!    offset      Total integral on the supercell of the final potential on output
!!                To be used only in the periodic case, ignored for other boundary conditions.
!!    sumpion     logical value which states whether to sum pot_ion to the final result or not
!!                if sumpion==.true. rhopot will be the Hartree potential + pot_ion+vxci
!!                                   pot_ion will be untouched
!!                if sumpion==.false. rhopot will be only the Hartree potential
!!                                    pot_ion will be the XC potential vxci
!! WARNING
!!    The dimensions of the arrays must be compatible with geocode, datacode, nproc, 
!!    ixc and iproc. Since the arguments of these routines are indicated with the *, it
!!    is IMPERATIVE to use the PS_dim4allocation routine for calculation arrays sizes.
!!    Moreover, for the cases with the exchange and correlation the density must be initialised
!!    to 10^-20 and not to zero.
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!! 
subroutine H_potential(geocode,datacode,iproc,nproc,n01,n02,n03,hx,hy,hz,&
     rhopot,karray,pot_ion,eh,offset,sumpion,&
     quiet) !optional argument
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=1), intent(in) :: datacode
  logical, intent(in) :: sumpion
  integer, intent(in) :: iproc,nproc,n01,n02,n03
  real(gp), intent(in) :: hx,hy,hz
  real(dp), intent(in) :: offset
  real(dp), dimension(*), intent(in) :: karray
  real(gp), intent(out) :: eh
  real(dp), dimension(*), intent(inout) :: rhopot
  real(wp), dimension(*), intent(inout) :: pot_ion
  character(len=3), intent(in), optional :: quiet
  !local variables
  character(len=*), parameter :: subname='H_potential'
  logical :: wrtmsg
  integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3,i3s_fake,i3xcsh_fake
  integer :: i_all,i_stat,ierr,ind,ind2,ind3,indp,ind2p,ind3p,i,j
  integer :: i1,i2,i3,j2,istart,iend,i3start,jend,jproc,i3xcsh,is_step,ind2nd
  integer :: nxc,nwbl,nwbr,nxt,nwb,nxcl,nxcr,nlim,ispin,istden,istglo
  real(dp) :: scal,ehartreeLOC,pot
  real(dp), dimension(:,:,:), allocatable :: zf
  integer, dimension(:,:), allocatable :: gather_arr

  call timing(iproc,'PSolv_comput  ','ON')

  !do not write anything on screen if quiet is set to yes
  if (present(quiet)) then
     if(quiet == 'yes' .or. quiet == 'YES') then
        wrtmsg=.false.
     else if(trim(quiet) == 'no' .or. trim(quiet) == 'NO') then
        wrtmsg=.true.
     else
        write(*,*)'ERROR: Unrecognised value for "quiet" option:',quiet
        stop
     end if
  else
     wrtmsg=.true.
  end if

  !calculate the dimensions wrt the geocode
  if (geocode == 'P') then
     if (iproc==0 .and. wrtmsg) &
          write(*,'(1x,a,3(i5),a,i5,a,i7,a)',advance='no')&
          'PSolver, periodic BC, dimensions: ',n01,n02,n03,'   proc',nproc,'  ixc:',0,' ... '
     call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'S') then
     if (iproc==0 .and. wrtmsg) &
          write(*,'(1x,a,3(i5),a,i5,a,i7,a)',advance='no')&
          'PSolver, surfaces BC, dimensions: ',n01,n02,n03,'   proc',nproc,'  ixc:',0,' ... '
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'F') then
     if (iproc==0 .and. wrtmsg) &
          write(*,'(1x,a,3(i5),a,i5,a,i7,a)',advance='no')&
          'PSolver, free  BC, dimensions: ',n01,n02,n03,'   proc',nproc,'  ixc:',0,' ... '
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else
     stop 'PSolver: geometry code not admitted'
  end if

  !array allocations
  allocate(zf(md1,md3,md2/nproc+ndebug),stat=i_stat)
  call memocc(i_stat,zf,'zf',subname)
  !initalise to zero the zf array
  zf=0.0_dp
  !call razero(md1*md3*md2/nproc,zf)

  !dimension for exchange-correlation (different in the global or distributed case)
  !let us calculate the dimension of the portion of the rhopot array to be passed 
  !to the xc routine
  !this portion will depend on the need of calculating the gradient or not, 
  !and whether the White-Bird correction must be inserted or not 
  !(absent only in the LB ixc=13 case)
  
  !nxc is the effective part of the third dimension that is being processed
  !nxt is the dimension of the part of rhopot that must be passed to the gradient routine
  !nwb is the dimension of the part of rhopot in the wb-postprocessing routine
  !note: nxc <= nwb <= nxt
  !the dimension are related by the values of nwbl and nwbr
  !      nxc+nxcl+nxcr-2 = nwb
  !      nwb+nwbl+nwbr = nxt
  istart=iproc*(md2/nproc)
  iend=min((iproc+1)*md2/nproc,m2)
  if (istart <= m2-1) then
     nxc=iend-istart
  else
     nxc=0
  end if

  if (datacode=='G') then
     !starting address of rhopot in the case of global i/o
     i3start=istart+1
  else if (datacode == 'D') then
     !distributed i/o
     i3start=1
  else
     stop 'PSolver: datacode not admitted'
  end if

  !this routine builds the values for each process of the potential (zf), multiplying by scal 

  !fill the array with the values of the charge density
  !no more overlap between planes
  !still the complex case should be defined
  do i3=1,nxc
     do i2=1,m3
        do i1=1,m1
           i=i1+(i2-1)*m1+(i3+i3start-2)*m1*m3
           zf(i1,i2,i3)=rhopot(i)
        end do
     end do
  end do

  if(geocode == 'P') then
     !no powers of hgrid because they are incorporated in the plane wave treatment
     scal=1.0_dp/real(n1*n2*n3,dp)
  else if (geocode == 'S') then
     !only one power of hgrid 
     !factor of -4*pi for the definition of the Poisson equation
     scal=-16.0_dp*atan(1.0_dp)*real(hy,dp)/real(n1*n2*n3,dp)
  else if (geocode == 'F') then
     !hgrid=max(hx,hy,hz)
     scal=hx*hy*hz/real(n1*n2*n3,dp)
  end if
  !here the case ncplx/= 1 should be added
  !eventually one may avoid zf array
  call timing(iproc,'PSolv_comput  ','OF')
  call G_PoissonSolver(geocode,iproc,nproc,1,n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,karray,&
       zf(1,1,1),&
       scal,hx,hy,hz,offset)
  call timing(iproc,'PSolv_comput  ','ON')

  
  !the value of the shift depends on the distributed i/o or not
  if (datacode=='G') then
     i3xcsh=istart !beware on the fact that this is not what represents its name!!!
     is_step=n01*n02*n03
  else if (datacode=='D') then
     i3xcsh=0 !shift not needed anymore
     is_step=m1*m3*nxt
  end if
 
  !if (iproc == 0) print *,'n03,nxt,nxc,geocode,datacode',n03,nxt,nxc,geocode,datacode

  ehartreeLOC=0.0_dp
  !recollect the final data
  !this part can be eventually removed once the zf disappears
  if (sumpion) then
     do j2=1,nxc
        i2=j2+i3xcsh 
        ind3=(i2-1)*n01*n02
        ind3p=(j2-1)*n01*n02
        do i3=1,m3
           ind2=(i3-1)*n01+ind3
           ind2p=(i3-1)*n01+ind3p
           do i1=1,m1
              ind=i1+ind2
              indp=i1+ind2p
              pot=zf(i1,i3,j2)
              ehartreeLOC=ehartreeLOC+rhopot(ind)*pot
              rhopot(ind)=real(pot,wp)+real(pot_ion(indp),wp)
           end do
        end do
     end do
  else
     do j2=1,nxc
        i2=j2+i3xcsh 
        ind3=(i2-1)*n01*n02
        do i3=1,m3
           ind2=(i3-1)*n01+ind3
           do i1=1,m1
              ind=i1+ind2
              pot=zf(i1,i3,j2)
              ehartreeLOC=ehartreeLOC+rhopot(ind)*pot
              rhopot(ind)=real(pot,wp)
           end do
        end do
     end do
  end if
  ehartreeLOC=ehartreeLOC*0.5_dp*hx*hy*hz
  
  i_all=-product(shape(zf))*kind(zf)
  deallocate(zf,stat=i_stat)
  call memocc(i_stat,i_all,'zf',subname)

  call timing(iproc,'PSolv_comput  ','OF')

  !gathering the data to obtain the distribution array
  !evaluating the total ehartree
  eh=real(ehartreeLOC,gp)
  if (nproc > 1) then
     call timing(iproc,'PSolv_commun  ','ON')

     eh=ehartreeLOC
     call mpiallred(eh,1,MPI_SUM,MPI_COMM_WORLD,ierr)

     call timing(iproc,'PSolv_commun  ','OF')

     if (datacode == 'G') then
        !building the array of the data to be sent from each process
        !and the array of the displacement

        call timing(iproc,'PSolv_comput  ','ON')
        allocate(gather_arr(0:nproc-1,2+ndebug),stat=i_stat)
        call memocc(i_stat,gather_arr,'gather_arr',subname)
        do jproc=0,nproc-1
           istart=min(jproc*(md2/nproc),m2-1)
           jend=max(min(md2/nproc,m2-md2/nproc*jproc),0)
           gather_arr(jproc,1)=m1*m3*jend
           gather_arr(jproc,2)=m1*m3*istart
        end do

        !gather all the results in the same rhopot array
        istart=min(iproc*(md2/nproc),m2-1)

        call timing(iproc,'PSolv_comput  ','OF')
        call timing(iproc,'PSolv_commun  ','ON')
        istden=1+n01*n02*istart
        istglo=1
        call MPI_ALLGATHERV(rhopot(istden),gather_arr(iproc,1),mpidtypw,&
             rhopot(istglo),gather_arr(0,1),gather_arr(0,2),mpidtypw,&
             MPI_COMM_WORLD,ierr)
        call timing(iproc,'PSolv_commun  ','OF')
        call timing(iproc,'PSolv_comput  ','ON')

        i_all=-product(shape(gather_arr))*kind(gather_arr)
        deallocate(gather_arr,stat=i_stat)
        call memocc(i_stat,i_all,'gather_arr',subname)

        call timing(iproc,'PSolv_comput  ','OF')

     end if
  end if

  !if(nspin==1 .and. ixc /= 0) eh=eh*2.0_gp
  if (iproc==0  .and. wrtmsg) write(*,'(a)')'done.'

end subroutine H_potential
!!***


!!****f* PSolver/PSolver
!! FUNCTION
!!    Calculate the Poisson equation $\nabla^2 V(x,y,z)=-4 \pi \rho(x,y,z)$
!!    from a given $\rho$, for different boundary conditions an for different data distributions.
!!    Following the boundary conditions, it applies the Poisson Kernel previously calculated.
!!
!! COPYRIGHT
!!    Copyright (C) 2002-2007 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SYNOPSIS
!!    geocode  Indicates the boundary conditions (BC) of the problem:
!!            'F' free BC, isolated systems.
!!                The program calculates the solution as if the given density is
!!                "alone" in R^3 space.
!!            'S' surface BC, isolated in y direction, periodic in xz plane                
!!                The given density is supposed to be periodic in the xz plane,
!!                so the dimensions in these direction mus be compatible with the FFT
!!                Beware of the fact that the isolated direction is y!
!!            'P' periodic BC.
!!                The density is supposed to be periodic in all the three directions,
!!                then all the dimensions must be compatible with the FFT.
!!                No need for setting up the kernel.
!!    datacode Indicates the distribution of the data of the input/output array:
!!            'G' global data. Each process has the whole array of the density 
!!                which will be overwritten with the whole array of the potential
!!            'D' distributed data. Each process has only the needed part of the density
!!                and of the potential. The data distribution is such that each processor
!!                has the xy planes needed for the calculation AND for the evaluation of the 
!!                gradient, needed for XC part, and for the White-Bird correction, which
!!                may lead up to 8 planes more on each side. Due to this fact, the information
!!                between the processors may overlap.
!!    nproc       number of processors
!!    iproc       label of the process,from 0 to nproc-1
!!    n01,n02,n03 global dimension in the three directions. They are the same no matter if the 
!!                datacode is in 'G' or in 'D' position.
!!    ixc         eXchange-Correlation code. Indicates the XC functional to be used 
!!                for calculating XC energies and potential. 
!!                ixc=0 indicates that no XC terms are computed. The XC functional codes follow
!!                the ABINIT convention.
!!    hx,hy,hz    grid spacings. For the isolated BC case for the moment they are supposed to 
!!                be equal in the three directions
!!    rhopot      main input/output array.
!!                On input, it represents the density values on the grid points
!!                On output, it is the Hartree potential, namely the solution of the Poisson 
!!                equation PLUS (when ixc/=0 sumpion=.true.) the XC potential 
!!                PLUS (again for ixc/=0 and sumpion=.true.) the pot_ion array. 
!!                The output is non overlapping, in the sense that it does not
!!                consider the points that are related to gradient and WB calculation
!!    karray      kernel of the poisson equation. It is provided in distributed case, with
!!                dimensions that are related to the output of the PS_dim4allocation routine
!!                it MUST be created by following the same geocode as the Poisson Solver.
!!    pot_ion     additional external potential that is added to the output, 
!!                when the XC parameter ixc/=0 and sumpion=.true., otherwise it corresponds 
!!                to the XC potential Vxc.
!!                When sumpion=.true., it is always provided in the distributed form,
!!                clearly without the overlapping terms which are needed only for the XC part
!!                When sumpion=.false. it is the XC potential and therefore it has 
!!                the same distribution of the data as the potential
!!                Ignored when ixc=0.
!!    eh,exc,vxc  Hartree energy, XC energy and integral of $\rho V_{xc}$ respectively
!!    offset      value of the potential at the point 1,1,1 of the grid.
!!                To be used only in the periodic case, ignored for other boundary conditions.
!!    sumpion     logical value which states whether to sum pot_ion to the final result or not
!!                if sumpion==.true. rhopot will be the Hartree potential + pot_ion+vxci
!!                                   pot_ion will be untouched
!!                if sumpion==.false. rhopot will be only the Hartree potential
!!                                    pot_ion will be the XC potential vxci
!!                this value is ignored when ixc=0. In that case pot_ion is untouched
!! WARNING
!!    The dimensions of the arrays must be compatible with geocode, datacode, nproc, 
!!    ixc and iproc. Since the arguments of these routines are indicated with the *, it
!!    is IMPERATIVE to use the PS_dim4allocation routine for calculation arrays sizes.
!!    Moreover, for the cases with the exchange and correlation the density must be initialised
!!    to 10^-20 and not to zero.
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!! 
subroutine PSolver(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
     rhopot,karray,pot_ion,eh,exc,vxc,offset,sumpion,nspin,&
     quiet) !optional argument
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=1), intent(in) :: datacode
  logical, intent(in) :: sumpion
  integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc,nspin
  real(gp), intent(in) :: hx,hy,hz
  real(dp), intent(in) :: offset
  real(dp), dimension(*), intent(in) :: karray
  real(gp), intent(out) :: eh,exc,vxc
  real(dp), dimension(*), intent(inout) :: rhopot
  real(wp), dimension(*), intent(inout) :: pot_ion
  character(len=3), intent(in), optional :: quiet
  !local variables
  character(len=*), parameter :: subname='PSolver'
  logical :: wrtmsg
  integer, parameter :: nordgr=4 !the order of the finite-difference gradient (fixed)
  integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3,i3s_fake,i3xcsh_fake
  integer :: i_all,i_stat,ierr,ind,ind2,ind3,ind4,ind4sh,i,j
  integer :: i1,i2,i3,j2,istart,iend,i3start,jend,jproc,i3xcsh,is_step,ind2nd
  integer :: nxc,nwbl,nwbr,nxt,nwb,nxcl,nxcr,nlim,ispin,istden,istglo
  real(dp) :: scal,ehartreeLOC,eexcuLOC,vexcuLOC,pot
  real(wp), dimension(:,:,:,:), allocatable :: zfionxc
  real(dp), dimension(:,:,:), allocatable :: zf
  integer, dimension(:,:), allocatable :: gather_arr
  real(dp), dimension(:), allocatable :: rhopot_G
  real(gp), dimension(:), allocatable :: energies_mpi

  call timing(iproc,'Exchangecorr  ','ON')


  !do not write anything on screen if quiet is set to yes
  if (present(quiet)) then
     if(quiet == 'yes' .or. quiet == 'YES') then
        wrtmsg=.false.
     else if(trim(quiet) == 'no' .or. trim(quiet) == 'NO') then
        wrtmsg=.true.
     else
        write(*,*)'ERROR: Unrecognised value for "quiet" option:',quiet
        stop
     end if
  else
     wrtmsg=.true.
  end if

  !calculate the dimensions wrt the geocode
  if (geocode == 'P') then
     if (iproc==0 .and. wrtmsg) &
          write(*,'(1x,a,3(i5),a,i5,a,i7,a)',advance='no')&
          'PSolver, periodic BC, dimensions: ',n01,n02,n03,'   proc',nproc,'  ixc:',ixc,' ... '
     call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'S') then
     if (iproc==0 .and. wrtmsg) &
          write(*,'(1x,a,3(i5),a,i5,a,i7,a)',advance='no')&
          'PSolver, surfaces BC, dimensions: ',n01,n02,n03,'   proc',nproc,'  ixc:',ixc,' ... '
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'F') then
     if (iproc==0 .and. wrtmsg) &
          write(*,'(1x,a,3(i5),a,i5,a,i7,a)',advance='no')&
          'PSolver, free  BC, dimensions: ',n01,n02,n03,'   proc',nproc,'  ixc:',ixc,' ... '
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else
     stop 'PSolver: geometry code not admitted'
  end if

  !array allocations
  allocate(zf(md1,md3,md2/nproc+ndebug),stat=i_stat)
  call memocc(i_stat,zf,'zf',subname)
  allocate(zfionxc(md1,md3,md2/nproc,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,zfionxc,'zfionxc',subname)

  !dimension for exchange-correlation (different in the global or distributed case)
  !let us calculate the dimension of the portion of the rhopot array to be passed 
  !to the xc routine
  !this portion will depend on the need of calculating the gradient or not, 
  !and whether the White-Bird correction must be inserted or not 
  !(absent only in the LB ixc=13 case)
  
  !nxc is the effective part of the third dimension that is being processed
  !nxt is the dimension of the part of rhopot that must be passed to the gradient routine
  !nwb is the dimension of the part of rhopot in the wb-postprocessing routine
  !note: nxc <= nwb <= nxt
  !the dimension are related by the values of nwbl and nwbr
  !      nxc+nxcl+nxcr-2 = nwb
  !      nwb+nwbl+nwbr = nxt
  istart=iproc*(md2/nproc)
  iend=min((iproc+1)*md2/nproc,m2)

  call xc_dimensions(geocode,ixc,istart,iend,m2,nxc,nxcl,nxcr,nwbl,nwbr,i3s_fake,i3xcsh_fake)
  nwb=nxcl+nxc+nxcr-2
  nxt=nwbr+nwb+nwbl

  if (datacode=='G') then
     !starting address of rhopot in the case of global i/o
     i3start=istart+2-nxcl-nwbl
     if((nspin==2 .and. nproc>1) .or. i3start <=0 .or. i3start+nxt-1 > n03 ) then
        !allocation of an auxiliary array for avoiding the shift of the density
        allocate(rhopot_g(m1*m3*nxt*2+ndebug),stat=i_stat)
        call memocc(i_stat,rhopot_g,'rhopot_g',subname)
        !here we should put the modulo of the results for the non-isolated GGA
        do ispin=1,nspin
           do i3=1,nxt
              do i2=1,m3
                 do i1=1,m1
                    i=i1+(i2-1)*m1+(i3-1)*m1*m3+(ispin-1)*m1*m3*nxt
                    j=i1+(i2-1)*n01+(modulo(i3start+i3-2,n03))*n01*n02+(ispin-1)*n01*n02*n03
                    rhopot_G(i)=rhopot(j)
                 end do
              end do
           end do
        end do
     end if
  else if (datacode == 'D') then
     !distributed i/o
     i3start=1
  else
     stop 'PSolver: datacode not admitted'
  end if

  !calculate the actual limit of the array for the zero padded FFT
  if (geocode == 'P') then
     nlim=n2
  else if (geocode == 'S') then
     nlim=n2
  else if (geocode == 'F') then
     nlim=n2/2
  end if

!!  print *,'density must go from',min(istart+1,m2),'to',iend,'with n2/2=',n2/2
!!  print *,'        it goes from',i3start+nwbl+nxcl-1,'to',i3start+nxc-1

  if (istart+1 <= m2) then 
     if(datacode=='G' .and. &
          ((nspin==2 .and. nproc > 1) .or. i3start <=0 .or. i3start+nxt-1 > n03 )) then
        !allocation of an auxiliary array for avoiding the shift 
        call xc_energy(geocode,m1,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,nxcl,nxcr,&
             ixc,hx,hy,hz,rhopot_G,pot_ion,sumpion,zf,zfionxc,&
             eexcuLOC,vexcuLOC,nproc,nspin)
        do ispin=1,nspin
           do i3=1,nxt
              do i2=1,m3
                 do i1=1,m1
                    i=i1+(i2-1)*m1+(i3-1)*m1*m3+(ispin-1)*m1*m3*nxt
                    j=i1+(i2-1)*n01+(modulo(i3start+i3-2,n03))*n01*n02+(ispin-1)*n01*n02*n03
                    rhopot(j)=rhopot_G(i)
                 end do
              end do
           end do
        end do
        !!          do i1=1,m1*m3*nxt
        !!             rhopot(n01*n02*(i3start-1)+i1)=rhopot_G(i1)
        !!          end do
        !!          do i1=1,m1*m3*nxt
        !!             rhopot(n01*n02*(i3start-1)+i1+n01*n02*n03)=rhopot_G(i1+m1*m3*nxt)
        !!          end do
        i_all=-product(shape(rhopot_G))*kind(rhopot_G)
        deallocate(rhopot_G,stat=i_stat)
        call memocc(i_stat,i_all,'rhopot_g',subname)
     else
        call xc_energy(geocode,m1,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,nxcl,nxcr,&
             ixc,hx,hy,hz,rhopot(1+n01*n02*(i3start-1)),pot_ion,sumpion,zf,zfionxc,&
             eexcuLOC,vexcuLOC,nproc,nspin)
     end if
  else if (istart+1 <= nlim) then !this condition ensures we have performed good zero padding
     do i2=istart+1,min(nlim,istart+md2/nproc)
        j2=i2-istart
        do i3=1,md3
           do i1=1,md1
              zf(i1,i3,j2)=0.0_dp
           end do
        end do
     end do
     eexcuLOC=0.0_dp
     vexcuLOC=0.0_dp
  else
     eexcuLOC=0.0_dp
     vexcuLOC=0.0_dp
  end if

  call timing(iproc,'Exchangecorr  ','OF')
  !this routine builds the values for each process of the potential (zf), multiplying by scal 
  if(geocode == 'P') then
     !no powers of hgrid because they are incorporated in the plane wave treatment
     scal=1.0_dp/real(n1*n2*n3,dp)
     !call P_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,karray,zf(1,1,1),&
     !     scal,hx,hy,hz,offset)
  else if (geocode == 'S') then
     !only one power of hgrid 
     !factor of -4*pi for the definition of the Poisson equation
     scal=-16.0_dp*atan(1.0_dp)*real(hy,dp)/real(n1*n2*n3,dp)
     !call S_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,karray,zf(1,1,1),&
     !     scal) !,hx,hy,hz,ehartreeLOC)
  else if (geocode == 'F') then
     !hgrid=max(hx,hy,hz)
     scal=hx*hy*hz/real(n1*n2*n3,dp)
     !call F_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,karray,zf(1,1,1),&
     !     scal)!,hgrid)!,ehartreeLOC)
  end if
  !here the case ncplx/= 1 should be added
  call G_PoissonSolver(geocode,iproc,nproc,1,n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,karray,zf(1,1,1),&
       scal,hx,hy,hz,offset)
  
  call timing(iproc,'PSolv_comput  ','ON')
  
  !the value of the shift depends on the distributed i/o or not
  if (datacode=='G') then
     i3xcsh=istart !beware on the fact that this is not what represents its name!!!
     is_step=n01*n02*n03
  else if (datacode=='D') then
     i3xcsh=nxcl+nwbl-1
     is_step=m1*m3*nxt
  end if
 
  !if (iproc == 0) print *,'n03,nxt,nxc,geocode,datacode',n03,nxt,nxc,geocode,datacode

  ehartreeLOC=0.0_dp
  !recollect the final data
  if (ixc==0) then !without XC the spin does not exist
     do j2=1,nxc
        i2=j2+i3xcsh !in this case the shift is always zero for a parallel run
        ind3=(i2-1)*n01*n02
        do i3=1,m3
           ind2=(i3-1)*n01+ind3
           do i1=1,m1
              ind=i1+ind2
              pot=zf(i1,i3,j2)
              ehartreeLOC=ehartreeLOC+rhopot(ind)*pot
              rhopot(ind)=real(pot,wp)
           end do
        end do
     end do
  else if (sumpion) then
     do j2=1,nxc
        i2=j2+i3xcsh
        ind3=(i2-1)*n01*n02
        do i3=1,m3
           ind2=(i3-1)*n01+ind3
           do i1=1,m1
              ind=i1+ind2
              pot=zf(i1,i3,j2)
              ehartreeLOC=ehartreeLOC+rhopot(ind)*pot
              rhopot(ind)=real(pot,wp)+real(zfionxc(i1,i3,j2,1),wp)
           end do
        end do
     end do
     !in the spin-polarised case the potential is given contiguously
     if (nspin==2) then
        !this start the count in the other component of the global array
        if (datacode=='G') ind=i3xcsh*n01*n02+n01*n02*n03
        do j2=1,nxc
           i2=j2+i3xcsh
           ind3=(i2-1)*n01*n02
           do i3=1,m3
              ind2=(i3-1)*n01+ind3
              do i1=1,m1
                 ind2nd=i1+ind2+is_step
                 ind=ind+1
                 pot=zf(i1,i3,j2)
                 ehartreeLOC=ehartreeLOC+rhopot(ind2nd)*pot
                 rhopot(ind)=real(pot,wp)+real(zfionxc(i1,i3,j2,2),wp)
              end do
           end do
        end do
     end if
  else
     if (datacode == 'G') then
        ind4sh=n01*n02*i3xcsh
     else
        ind4sh=0
     end if
     do j2=1,nxc
        i2=j2+i3xcsh
        ind3=(i2-1)*n01*n02
        do i3=1,m3
           ind2=(i3-1)*n01+ind3
           do i1=1,m1
              ind=i1+(i3-1)*n01+(i2-1)*n01*n02
              ind4=i1+(i3-1)*n01+(j2-1)*n01*n02+ind4sh
              pot=zf(i1,i3,j2)
              ehartreeLOC=ehartreeLOC+rhopot(ind)*pot
              rhopot(ind)=real(pot,wp)
              pot_ion(ind4)=zfionxc(i1,i3,j2,1)
           end do
        end do
     end do
     !in the spin-polarised (distributed) case the potential is given contiguously
     if (nspin==2) then
        if (datacode == 'D') ind4sh=-n01*n02*(nxt-nxc)
        do j2=1,nxc
           i2=j2+i3xcsh
           ind3=(i2-1)*n01*n02
           do i3=1,m3
              ind2=(i3-1)*n01+ind3
              do i1=1,m1
                 ind=i1+(i3-1)*n01+(i2-1)*n01*n02+is_step
                 ind4=i1+(i3-1)*n01+(j2-1)*n01*n02+is_step+ind4sh
                 pot=zf(i1,i3,j2)
                 ehartreeLOC=ehartreeLOC+rhopot(ind)*pot
                 pot_ion(ind4)=zfionxc(i1,i3,j2,2)
              end do
           end do
        end do
     end if
  end if
  ehartreeLOC=ehartreeLOC*0.5_dp*hx*hy*hz
  
  i_all=-product(shape(zf))*kind(zf)
  deallocate(zf,stat=i_stat)
  call memocc(i_stat,i_all,'zf',subname)
  i_all=-product(shape(zfionxc))*kind(zfionxc)
  deallocate(zfionxc,stat=i_stat)
  call memocc(i_stat,i_all,'zfionxc',subname)

  call timing(iproc,'PSolv_comput  ','OF')

  !gathering the data to obtain the distribution array
  !evaluating the total ehartree,eexcu,vexcu
  if (nproc > 1) then

     call timing(iproc,'PSolv_commun  ','ON')
     allocate(energies_mpi(6+ndebug),stat=i_stat)
     call memocc(i_stat,energies_mpi,'energies_mpi',subname)

     energies_mpi(1)=ehartreeLOC
     energies_mpi(2)=eexcuLOC
     energies_mpi(3)=vexcuLOC
     call MPI_ALLREDUCE(energies_mpi(1),energies_mpi(4),3,mpidtypd,MPI_SUM,MPI_COMM_WORLD,ierr)
     eh=energies_mpi(4)
     exc=energies_mpi(5)
     vxc=energies_mpi(6)

     i_all=-product(shape(energies_mpi))*kind(energies_mpi)
     deallocate(energies_mpi,stat=i_stat)
     call memocc(i_stat,i_all,'energies_mpi',subname)
     call timing(iproc,'PSolv_commun  ','OF')

     if (datacode == 'G') then
        !building the array of the data to be sent from each process
        !and the array of the displacement

        call timing(iproc,'PSolv_comput  ','ON')
        allocate(gather_arr(0:nproc-1,2+ndebug),stat=i_stat)
        call memocc(i_stat,gather_arr,'gather_arr',subname)
        do jproc=0,nproc-1
           istart=min(jproc*(md2/nproc),m2-1)
           jend=max(min(md2/nproc,m2-md2/nproc*jproc),0)
           gather_arr(jproc,1)=m1*m3*jend
           gather_arr(jproc,2)=m1*m3*istart
        end do

        !gather all the results in the same rhopot array
        istart=min(iproc*(md2/nproc),m2-1)

        call timing(iproc,'PSolv_comput  ','OF')
        call timing(iproc,'PSolv_commun  ','ON')
        istden=1+n01*n02*istart
        istglo=1
        do ispin=1,nspin
           if (ispin==2) then
              istden=istden+n01*n02*n03
              istglo=istglo+n01*n02*n03
           end if
           call MPI_ALLGATHERV(rhopot(istden),gather_arr(iproc,1),mpidtypw,&
                rhopot(istglo),gather_arr(0,1),gather_arr(0,2),mpidtypw,&
                MPI_COMM_WORLD,ierr)
           !if it is the case gather also the results of the XC potential
           if (ixc /=0 .and. .not. sumpion) then
              call MPI_ALLGATHERV(pot_ion(istden),gather_arr(iproc,1),&
                   mpidtypw,pot_ion(istglo),gather_arr(0,1),gather_arr(0,2),&
                   mpidtypw,MPI_COMM_WORLD,ierr)
           end if
        end do
        call timing(iproc,'PSolv_commun  ','OF')
        call timing(iproc,'PSolv_comput  ','ON')

        i_all=-product(shape(gather_arr))*kind(gather_arr)
        deallocate(gather_arr,stat=i_stat)
        call memocc(i_stat,i_all,'gather_arr',subname)

        call timing(iproc,'PSolv_comput  ','OF')

     end if

  else
     eh=real(ehartreeLOC,gp)
     exc=real(eexcuLOC,gp)
     vxc=real(vexcuLOC,gp)
  end if

  if(nspin==1 .and. ixc /= 0) eh=eh*2.0_gp
  if (iproc==0  .and. wrtmsg) write(*,'(a)')'done.'

end subroutine PSolver
!!***


!!****f* PSolver/PSolverNC
!! FUNCTION
!!    Transforms a generalized spin density into a pointwise collinear spin density which is
!!    then passed to the Poisson Solver (PSolver). 
!!
!! COPYRIGHT
!!    Copyright (C) 2002-2007 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SYNOPSIS
!!    geocode  Indicates the boundary conditions (BC) of the problem:
!!            'F' free BC, isolated systems.
!!                The program calculates the solution as if the given density is
!!                "alone" in R^3 space.
!!            'S' surface BC, isolated in y direction, periodic in xz plane                
!!                The given density is supposed to be periodic in the xz plane,
!!                so the dimensions in these direction mus be compatible with the FFT
!!                Beware of the fact that the isolated direction is y!
!!            'P' periodic BC.
!!                The density is supposed to be periodic in all the three directions,
!!                then all the dimensions must be compatible with the FFT.
!!                No need for setting up the kernel.
!!    datacode Indicates the distribution of the data of the input/output array:
!!            'G' global data. Each process has the whole array of the density 
!!                which will be overwritten with the whole array of the potential
!!            'D' distributed data. Each process has only the needed part of the density
!!                and of the potential. The data distribution is such that each processor
!!                has the xy planes needed for the calculation AND for the evaluation of the 
!!                gradient, needed for XC part, and for the White-Bird correction, which
!!                may lead up to 8 planes more on each side. Due to this fact, the information
!!                between the processors may overlap.
!!    nproc       number of processors
!!    iproc       label of the process,from 0 to nproc-1
!!    n01,n02,n03 global dimension in the three directions. They are the same no matter if the 
!!                datacode is in 'G' or in 'D' position.
!!    n3d         third dimension of the density. For distributed data, it takes into account 
!!                the enlarging needed for calculating the XC functionals.
!!                For global data it is simply equal to n03. 
!!                When there are too many processes and there is no room for the density n3d=0
!!    ixc         eXchange-Correlation code. Indicates the XC functional to be used 
!!                for calculating XC energies and potential. 
!!                ixc=0 indicates that no XC terms are computed. The XC functional codes follow
!!                the ABINIT convention.
!!    hx,hy,hz    grid spacings. For the isolated BC case for the moment they are supposed to 
!!                be equal in the three directions
!!    rhopot      main input/output array.
!!                On input, it represents the density values on the grid points
!!                On output, it is the Hartree potential, namely the solution of the Poisson 
!!                equation PLUS (when ixc/=0 sumpion=.true.) the XC potential 
!!                PLUS (again for ixc/=0 and sumpion=.true.) the pot_ion array. 
!!                The output is non overlapping, in the sense that it does not
!!                consider the points that are related to gradient and WB calculation
!!    karray      kernel of the poisson equation. It is provided in distributed case, with
!!                dimensions that are related to the output of the PS_dim4allocation routine
!!                it MUST be created by following the same geocode as the Poisson Solver.
!!    pot_ion     additional external potential that is added to the output, 
!!                when the XC parameter ixc/=0 and sumpion=.true., otherwise it corresponds 
!!                to the XC potential Vxc.
!!                When sumpion=.true., it is always provided in the distributed form,
!!                clearly without the overlapping terms which are needed only for the XC part
!!                When sumpion=.false. it is the XC potential and therefore it has 
!!                the same distribution of the data as the potential
!!                Ignored when ixc=0.
!!    eh,exc,vxc  Hartree energy, XC energy and integral of $\rho V_{xc}$ respectively
!!    offset      value of the potential at the point 1,1,1 of the grid.
!!                To be used only in the periodic case, ignored for other boundary conditions.
!!    sumpion     logical value which states whether to sum pot_ion to the final result or not
!!                if sumpion==.true. rhopot will be the Hartree potential + pot_ion+vxci
!!                                   pot_ion will be untouched
!!                if sumpion==.false. rhopot will be only the Hartree potential
!!                                    pot_ion will be the XC potential vxci
!!                this value is ignored when ixc=0. In that case pot_ion is untouched
!! WARNING
!!    The dimensions of the arrays must be compatible with geocode, datacode, nproc, 
!!    ixc and iproc. Since the arguments of these routines are indicated with the *, it
!!    is IMPERATIVE to use the PS_dim4allocation routine for calculation arrays sizes.
!!    Moreover, for the cases with the exchange and correlation the density must be initialised
!!    to 10^-20 and not to zero.
!!
!! AUTHOR
!!    Anders Bergman
!! CREATION DATE
!!    March 2008
!!
!! SOURCE
!! 
subroutine PSolverNC(geocode,datacode,iproc,nproc,n01,n02,n03,n3d,ixc,hx,hy,hz,&
     rhopot,karray,pot_ion,eh,exc,vxc,offset,sumpion,nspin)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=1), intent(in) :: datacode
  logical, intent(in) :: sumpion
  integer, intent(in) :: iproc,nproc,n01,n02,n03,n3d,ixc,nspin
  real(gp), intent(in) :: hx,hy,hz
  real(dp), intent(in) :: offset
  real(dp), dimension(*), intent(in) :: karray
  real(gp), intent(out) :: eh,exc,vxc
  real(wp), dimension(*), intent(inout) :: pot_ion
  real(dp), dimension(*), intent(inout) :: rhopot
  !local variables
  character(len=*), parameter :: subname='PSolverNC'
  real(dp) :: rhon,rhos,factor
  integer :: i_all,i_stat,i1,i2,i3,idx,offs
  real(dp), dimension(:,:,:), allocatable :: m_norm
  real(dp), dimension(:,:,:,:), allocatable :: rho_diag
  
  
  if(nspin==4) then
     !Allocate diagonal spin-density in real space
     if (n3d >0) then
        allocate(rho_diag(n01,n02,n3d,2+ndebug),stat=i_stat)
        call memocc(i_stat,rho_diag,'rho_diag',subname)
        allocate(m_norm(n01,n02,n3d+ndebug),stat=i_stat)
        call memocc(i_stat,m_norm,'m_norm',subname)
        !           print *,'Rho Dims',shape(rhopot),shape(rho_diag)
        idx=1
        offs=n01*n02*n3d 
        do i3=1,n3d
           do i2=1,n02
              do i1=1,n01
                 !rho_diag(i1,i2,i3,1)=rhopot(i1,i2,i3,1)
                 m_norm(i1,i2,i3)=&
                      sqrt(rhopot(idx+offs)**2+rhopot(idx+2*offs)**2+rhopot(idx+3*offs)**2)
                 rho_diag(i1,i2,i3,1)=(rhopot(idx)+m_norm(i1,i2,i3))*0.5_dp+1.00e-20
                 rho_diag(i1,i2,i3,2)=(rhopot(idx)-m_norm(i1,i2,i3))*0.5_dp+1.00e-20
                 idx=idx+1
                 !m_norm(i1,i2,i3)=sqrt(rhopot(i1,i2,i3,2)**2+rhopot(i1,i2,i3,3)**2+rhopot(i1,i2,i3,4)**2)
                 !rho_diag(i1,i2,i3,1)=(rhopot(i1,i2,i3,1)+m_norm(i1,i2,i3))*0.5d0!+1.00e-20
                 !rho_diag(i1,i2,i3,2)=(rhopot(i1,i2,i3,1)-m_norm(i1,i2,i3))*0.5d0!+1.00e-20
              end do
           end do
        end do
     else
        allocate(rho_diag(1,1,1,2+ndebug),stat=i_stat)
        call memocc(i_stat,rho_diag,'rho_diag',subname)
        allocate(m_norm(1,1,1+ndebug),stat=i_stat)
        call memocc(i_stat,m_norm,'m_norm',subname)
        rho_diag=0.0_dp
        m_norm=0.0_dp
     end if
     
     call PSolver(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
          rho_diag,karray,pot_ion,eh,exc,vxc,offset,sumpion,2)
     !print *,'Psolver R',eh,exc,vxc
     !open(17)
     idx=1
     do i3=1,n3d
        do i2=1,n02
           do i1=1,n01
              rhon=(rho_diag(i1,i2,i3,1)+rho_diag(i1,i2,i3,2))*0.5_dp
              rhos=(rho_diag(i1,i2,i3,1)-rho_diag(i1,i2,i3,2))*0.5_dp
              if(m_norm(i1,i2,i3)>rhopot(idx)*4.0e-20_dp)then
                 !                 if(m_norm(i1,i2,i3)>rhopot(i1,i2,i3,1)*4.0e-20)then
                 factor=rhos/m_norm(i1,i2,i3)
              else
                 factor=0.0_dp
              end if
              !write(17,'(3(i0,1x),5(1pe12.5))')i1,i2,i3,rhon,rhos,rho_diag(i1,i2,i3,1),rho_diag(i1,i2,i3,2),factor
              rhopot(idx)=rhon+rhopot(idx+3*offs)*factor
              rhopot(idx+offs)=rhopot(idx+offs)*factor
              rhopot(idx+2*offs)=-rhopot(idx+2*offs)*factor
              rhopot(idx+3*offs)=rhon-rhopot(idx+3*offs)*factor
              idx=idx+1
              !                 rhopot(i1,i2,i3,1)=rhon+rhopot(i1,i2,i3,4)*factor
              !                 rhopot(i1,i2,i3,2)=rhopot(i1,i2,i3,2)*factor
              !                 rhopot(i1,i2,i3,3)=-rhopot(i1,i2,i3,3)*factor
              !                 rhopot(i1,i2,i3,4)=rhon-rhopot(i1,i2,i3,4)*factor
           end do
        end do
     end do
     !close(17)
     i_all=-product(shape(rho_diag))*kind(rho_diag)
     deallocate(rho_diag,stat=i_stat)
     call memocc(i_stat,i_all,'rho_diag',subname)
     i_all=-product(shape(m_norm))*kind(m_norm)
     deallocate(m_norm,stat=i_stat)
     call memocc(i_stat,i_all,'m_norm',subname)
  else
     call PSolver(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
          rhopot,karray,pot_ion,eh,exc,vxc,offset,sumpion,nspin)
  end if

end subroutine PSolverNC
!!***


!!****f* PSolver/PS_dim4allocation
!! FUNCTION
!!    Calculate the dimensions needed for the allocation of the arrays 
!!    related to the Poisson Solver
!!
!! SYNOPSIS
!!    geocode  Indicates the boundary conditions (BC) of the problem:
!!            'F' free BC, isolated systems.
!!                The program calculates the solution as if the given density is
!!                "alone" in R^3 space.
!!            'S' surface BC, isolated in y direction, periodic in xz plane                
!!                The given density is supposed to be periodic in the xz plane,
!!                so the dimensions in these direction mus be compatible with the FFT
!!                Beware of the fact that the isolated direction is y!
!!            'P' periodic BC.
!!                The density is supposed to be periodic in all the three directions,
!!                then all the dimensions must be compatible with the FFT.
!!                No need for setting up the kernel.
!!    datacode Indicates the distribution of the data of the input/output array:
!!            'G' global data. Each process has the whole array of the density 
!!                which will be overwritten with the whole array of the potential
!!            'D' distributed data. Each process has only the needed part of the density
!!                and of the potential. The data distribution is such that each processor
!!                has the xy planes needed for the calculation AND for the evaluation of the 
!!                gradient, needed for XC part, and for the White-Bird correction, which
!!                may lead up to 8 planes more on each side. Due to this fact, the information
!!                between the processors may overlap.
!!    iproc,nproc number of process, number of processes
!!    n01,n02,n03 dimensions of the real space grid to be hit with the Poisson Solver
!!    ixc         eXchange-Correlation code. Indicates the XC functional to be used 
!!                for calculating XC energies and potential. 
!!                ixc=0 indicates that no XC terms are computed. The XC functional codes follow
!!                the ABINIT convention.
!!    n3d         third dimension of the density. For distributed data, it takes into account 
!!                the enlarging needed for calculating the XC functionals.
!!                For global data it is simply equal to n03. 
!!                When there are too many processes and there is no room for the density n3d=0
!!    n3p         third dimension for the potential. The same as n3d, but without 
!!                taking into account the enlargment for the XC part. For non-GGA XC, n3p=n3d.
!!    n3pi        Dimension of the pot_ion array, always with distributed data. 
!!                For distributed data n3pi=n3p
!!    i3xcsh      Shift of the density that must be performed to enter in the 
!!                non-overlapping region. Useful for recovering the values of the potential
!!                when using GGA XC functionals. If the density starts from rhopot(1,1,1),
!!                the potential starts from rhopot(1,1,i3xcsh+1). 
!!                For non-GGA XCs and for global distribution data i3xcsh=0
!!    i3s         Starting point of the density effectively treated by each processor 
!!                in the third direction.
!!                It takes into account also the XC enlarging. The array rhopot will correspond
!!                To the planes of third coordinate from i3s to i3s+n3d-1. 
!!                The potential to the planes from i3s+i3xcsh to i3s+i3xcsh+n3p-1
!!                The array pot_ion to the planes from i3s+i3xcsh to i3s+i3xcsh+n3pi-1
!!                For global disposition i3s is equal to distributed case with i3xcsh=0.
!!
!!
!! WARNING
!!    The XC enlarging due to GGA part is not present for surfaces and 
!!    periodic boundary condition. This is related to the fact that the calculation of the
!!    gradient and the White-Bird correction are not yet implemented for non-isolated systems
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!!
subroutine PS_dim4allocation(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,&
     n3d,n3p,n3pi,i3xcsh,i3s)
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=1), intent(in) :: datacode
  integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc
  integer, intent(out) :: n3d,n3p,n3pi,i3xcsh,i3s
  !local variables
  integer, parameter :: nordgr=4
  integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
  integer :: istart,iend,nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr


  !calculate the dimensions wrt the geocode
  if (geocode == 'P') then
     call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'S') then
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'F') then
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else
     stop 'PS_dim4allocation: geometry code not admitted'
  end if

  !formal start and end of the slice
  istart=iproc*(md2/nproc)
  iend=min((iproc+1)*md2/nproc,m2)

  if (datacode == 'D') then
     call xc_dimensions(geocode,ixc,istart,iend,m2,nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh)

     nwb=nxcl+nxc+nxcr-2
     nxt=nwbr+nwb+nwbl

     n3p=nxc
     n3d=nxt
     n3pi=n3p
  else if (datacode == 'G') then
     n3d=n03
     n3p=n03
     i3xcsh=0
     i3s=min(istart,m2-1)+1
     n3pi=max(iend-istart,0)
  else
     print *,datacode
     stop 'PS_dim4allocation: data code not admitted'
  end if

!!  print *,'P4,iproc',iproc,'nxc,ncxl,ncxr,nwbl,nwbr',nxc,nxcl,nxcr,nwbl,nwbr,&
!!       'ixc,n3d,n3p,i3xcsh,i3s',ixc,n3d,n3p,i3xcsh,i3s

end subroutine PS_dim4allocation
!!***


!!****f* PSolver/xc_dimensions
!! FUNCTION
!!   Calculate the dimensions to be used for the XC part, taking into account also
!!   the White-bird correction which should be made for some GGA functionals
!!
!! SYNOPSIS
!!    geocode   choice of the boundary conditions
!!
!!    ixc       XC id
!!
!!    m2        dimension to be parallelised
!!   
!!    nxc       size of the parallelised XC potential
!!
!!    ncxl,ncxr left and right buffers for calculating the WB correction after call drivexc
!!    nwbl,nwbr left and right buffers for calculating the gradient to pass to drivexc
!!    
!!    i3s       starting addres of the distributed dimension
!!
!!    i3xcsh    shift to be applied to i3s for having the striting address of the potential
!!
!! WARNING
!!    It is imperative that iend <=m2
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    May 2008
!!
!! SOURCE
!!
subroutine xc_dimensions(geocode,ixc,istart,iend,m2,nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh)
  use libxc_functionals

  implicit none

  character(len=1), intent(in) :: geocode
  integer, intent(in) :: ixc,istart,iend,m2
  integer, intent(out) :: nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh
  !local variables
  integer, parameter :: nordgr=4 !the order of the finite-difference gradient (fixed)
  logical :: use_gradient

  if (istart <= m2-1) then
     use_gradient = (ixc >= 11 .and. ixc <= 16) .or. &
          & (ixc < 0 .and. libxc_functionals_isgga())
     nxc=iend-istart
     if (use_gradient .and. geocode == 'F') then
        if (ixc==13) then
           !now the dimension of the part required for the gradient
           nwbl=min(istart,nordgr)
           nwbr=min(m2-iend,nordgr) !always m2 < iend
           nxcl=1
           nxcr=1
        else
           !now the dimension of the part required for the gradient
           if(istart<=nordgr) then
              nxcl=istart+1
              nwbl=0
           else
              nxcl=nordgr+1
              nwbl=min(nordgr,istart-nordgr)
           end if
           if(iend>=m2-nordgr+1) then
              nxcr=m2-iend+1
              nwbr=0
           else
              nxcr=nordgr+1
              nwbr=min(nordgr,m2-nordgr-iend)
           end if
        end if
     else if (geocode /= 'F' .and. use_gradient .and. nxc /= m2) then
        if (ixc==13) then
           !now the dimension of the part required for the gradient
           nwbl=nordgr
           nwbr=nordgr
           nxcl=1
           nxcr=1
        else
           nxcl=nordgr+1
           nwbl=nordgr
           nxcr=nordgr+1
           nwbr=nordgr
        end if
     !this case is also considered below
     !else if (geocode /= 'F' .and. use_gradient .and. nxc == m2) then
     else 
        nwbl=0
        nwbr=0
        nxcl=1
        nxcr=1
     end if
     i3xcsh=nxcl+nwbl-1
     i3s=istart+1-i3xcsh
  else
     nwbl=0
     nwbr=0
     nxcl=1
     nxcr=1
     nxc=0
     i3xcsh=0
     i3s=m2
  end if
end subroutine xc_dimensions
!!***


!!****f* PSolver/P_FFT_dimensions
!! FUNCTION
!!    Calculate four sets of dimension needed for the calculation of the
!!    convolution for the periodic system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with m2 and m3 exchanged
!!
!!    n1,n2,n3 the first FFT dimensions, for the moment supposed to be even
!!
!!    md1,md2,md3 the n1,n2,n3 dimensions. They contain the real unpadded space,
!!                properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! WARNING
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the other geometries of the Poisson Solver.
!!    The dimensions 2 and 3 are exchanged.
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    October 2006
!!
!! SOURCE
!!
subroutine P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
 implicit none
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 l3=m3 !beware of the half dimension
    call fourier_dim(l1,n1)
    if (n1 == m1) then
    else
       print *,'the FFT in the x direction is not allowed'
       print *,'n01 dimension',n01
       stop
    end if
    l1=l1+1
    call fourier_dim(l2,n2)
    if (n2 == m2) then
    else
       print *,'the FFT in the z direction is not allowed'
       print *,'n03 dimension',n03
       stop
    end if
    call fourier_dim(l3,n3)
    if (n3 == m3) then
    else
       print *,'the FFT in the y direction is not allowed'
       print *,'n02 dimension',n02
       stop
    end if

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3
151 if (nproc*(md2/nproc) < n2) then
    md2=md2+1
    goto 151
 endif


 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1 
 nd2=n2/2+1
 nd3=n3/2+1
250 if (modulo(nd3,nproc) /= 0) then
    nd3=nd3+1
    goto 250
 endif

end subroutine P_FFT_dimensions
!!***


!!****f* PSolver/S_FFT_dimensions
!! FUNCTION
!!    Calculate four sets of dimension needed for the calculation of the
!!    convolution for the surface system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with 2 and 3 exchanged
!!
!!    n1,n2 the first FFT dimensions, for the moment supposed to be even
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2     the n1,n2 dimensions. 
!!    md3         half of n3 dimension. They contain the real unpadded space,
!!                properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! WARNING
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the Poisson Solver with other geometries.
!!    Dimensions n02 and n03 were exchanged
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    October 2006
!!
!! SOURCE
!!
subroutine S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
 implicit none
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 l3=m3 !beware of the half dimension
    call fourier_dim(l1,n1)
    if (n1 == m1) then
    else
       print *,'the FFT in the x direction is not allowed'
       print *,'n01 dimension',n01
       stop
    end if
    l1=l1+1
    call fourier_dim(l2,n2)
    if (n2 == m2) then
    else
       print *,'the FFT in the z direction is not allowed'
       print *,'n03 dimension',n03
       stop
    end if
 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3/2
151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    goto 151
 endif


 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc

 !these two dimensions are like that since they are even
 nd1=n1/2+1
 nd2=n2/2+1

 nd3=n3/2+1
250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    goto 250
 endif

end subroutine S_FFT_dimensions
!!***


!!****f* PSolver/F_FFT_dimensions
!! FUNCTION
!!    Calculate four sets of dimension needed for the calculation of the
!!    zero-padded convolution
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension with the dimension 2 and 3 exchanged
!!
!!    n1,n2 the first FFT even dimensions greater that 2*m1, 2*m2
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2,md3 half of n1,n2,n3 dimension. They contain the real unpadded space,
!!                properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! WARNING
!!    The dimension m2 and m3 correspond to n03 and n02 respectively
!!    this is needed since the convolution routine manage arrays of dimension
!!    (md1,md3,md2/nproc)
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2006
!!
!! SOURCE
!!
subroutine F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
 implicit none
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space, inverted for convenience
 m1=n01
 m2=n03
 m3=n02
 ! real space grid dimension (suitable for number of processors)
 l1=2*m1
 l2=2*m2
 l3=m3 !beware of the half dimension
 do
    call fourier_dim(l1,n1)
    if (modulo(n1,2) == 0) then
       exit
    end if
    l1=l1+1
 end do
 do
    call fourier_dim(l2,n2)
    if (modulo(n2,2) == 0) then
       exit
    end if
    l2=l2+1
 end do
 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2/2
 md3=n3/2
151 if (nproc*(md2/nproc).lt.n2/2) then
    md2=md2+1
    goto 151
 endif

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1

250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    goto 250
 endif

end subroutine F_FFT_dimensions
!!***
