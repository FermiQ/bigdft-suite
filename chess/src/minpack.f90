module minpack

  implicit none

  private

  public :: lmdif_wrapper

  contains

    subroutine lmdif_wrapper(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn, &
                 diag,mode,factor,nprint,info,nfev,fjac,ldfjac, &
                 ipvt,qtf,wa1,wa2,wa3,wa4)
      integer,intent(in) :: m,n,maxfev,mode,nprint,ldfjac
      integer,intent(out) :: info,nfev
      integer,dimension(n),intent(out) :: ipvt
      real(kind=8),intent(in) :: ftol,xtol,gtol,epsfcn,factor
      real(kind=8),dimension(n),intent(inout)  :: x, diag, wa1, wa2, wa3
      real(kind=8),dimension(m),intent(out) :: fvec
      real(kind=8),dimension(ldfjac,n),intent(inout) :: fjac
      real(kind=8),dimension(n),intent(out) :: qtf
      real(kind=8),dimension(m),intent(inout) :: wa4
      external :: fcn

      call lmdif_wrapper(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn, &
           diag,mode,factor,nprint,info,nfev,fjac,ldfjac, &
           ipvt,qtf,wa1,wa2,wa3,wa4)

    end subroutine lmdif_wrapper

end module minpack
