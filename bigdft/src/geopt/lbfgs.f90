!> @file
!>   This file contains the LBFGS algorithm and supporting routines
!! @author
!!  Copyright (C) 2010-2011 BigDFT group 
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~/COPYING file
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the list of contributors, see ~/AUTHORS


!> Module which define the type parameterminimization
module minpar
   implicit none

   type parameterminimization
      !> general parameters for all methods
      character (len=10) :: approach
      integer :: iter=0
      integer :: iflag=0
      integer :: history
      logical ::converged
      !>parameters for print information
      integer :: verbosity
      integer :: MSAVE
      integer :: MP
      integer :: LP
      integer :: MAXFEV
      integer :: FINSTEP
      double precision :: ALPHA 
      double precision :: GTOL
      double precision :: XTOL
      double precision :: FTOL
      double precision :: STPMIN
      double precision :: STPMAX
      logical :: DIAGCO
      logical :: IWRITE
   end type parameterminimization

   !LG: this is a bad practice, the local variables should never go inside the modules
   type(parameterminimization), save :: parmin

end module minpar

!     ----------------------------------------------------------------------
!
!     ****************
!     LBFGS SUBROUTINE
!     ****************
!
!      SUBROUTINE LBFGS(IPROC,IN,N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
!      SUBROUTINE LBFGS(IPROC,IN,PAR,N,M,X,F,G,DIAG,IPRINT,EPS,W,IFLAG)
      SUBROUTINE LBFGS(IPROC,IN,PARMIN,N,M,X,F,G,DIAG,IPRINT,EPS,W,IFLAG)
      use module_types
!
!      use par_driver , only:driverparameters 
      use minpar, only: parameterminimization
      IMPLICIT NONE
!      type(driverparameters)::par
      type(parameterminimization)::parmin
      INTEGER N,M,IPRINT(2),IFLAG,IPROC
      DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
      DOUBLE PRECISION F,EPS
!      DOUBLE PRECISION F,EPS,XTOL
!      LOGICAL DIAGCO
!
!        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!                          JORGE NOCEDAL
!                        *** July 1990 ***
!
! 
!     This subroutine solves the unconstrained minimization problem
! 
!                      min F(x),    x= (x1,x2,...,xN),
!
!      using the limited memory BFGS method. The routine is especially
!      effective on problems involving a large number of variables. In
!      a typical iteration of this method an approximation Hk to the
!      inverse of the Hessian is obtained by applying M BFGS updates to
!      a diagonal matrix Hk0, using information from the previous M steps.
!      The user specifies the number M, which determines the amount of
!      storage required by the routine. The user may also provide the
!      diagonal matrices Hk0 if not satisfied with the default choice.
!      The algorithm is described in "On the limited memory BFGS method
!      for large scale optimization", by D. Liu and J. Nocedal,
!      Mathematical Programming B 45 (1989) 503-528.
! 
!      The user is required to calculate the function value F and its
!      gradient G. In order to allow the user complete control over
!      these computations, reverse  communication is used. The routine
!      must be called repeatedly under the control of the parameter
!      IFLAG. 
!
!      The steplength is determined at each iteration by means of the
!      line search routine MCVSRCH, which is a slight modification of
!      the routine CSRCH written by More' and Thuente.
! 
!      The calling statement is 
! 
!          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
! 
!      where
! 
!     N       is an INTEGER variable that must be set by the user to the
!             number of variables. It is not altered by the routine.
!             Restriction: N>0.
! 
!     M       is an INTEGER variable that must be set by the user to
!             the number of corrections used in the BFGS update. It
!             is not altered by the routine. Values of M less than 3 are
!             not recommended; large values of M will result in excessive
!             computing time. 3<= M <=7 is recommended. Restriction: M>0.
! 
!     X       is a DOUBLE PRECISION array of length N. On initial entry
!             it must be set by the user to the values of the initial
!             estimate of the solution vector. On exit with IFLAG=0, it
!             contains the values of the variables at the best point
!             found (usually a solution).
! 
!     F       is a DOUBLE PRECISION variable. Before initial entry and on
!             a re-entry with IFLAG=1, it must be set by the user to
!             contain the value of the function F at the point X.
! 
!     G       is a DOUBLE PRECISION array of length N. Before initial
!             entry and on a re-entry with IFLAG=1, it must be set by
!             the user to contain the components of the gradient G at
!             the point X.
! 
!     DIAGCO  is a LOGICAL variable that must be set to .TRUE. if the
!             user  wishes to provide the diagonal matrix Hk0 at each
!             iteration. Otherwise it should be set to .FALSE., in which
!             case  LBFGS will use a default value described below. If
!             DIAGCO is set to .TRUE. the routine will return at each
!             iteration of the algorithm with IFLAG=2, and the diagonal
!              matrix Hk0  must be provided in the array DIAG.
! 
! 
!     DIAG    is a DOUBLE PRECISION array of length N. If DIAGCO=.TRUE.,
!             then on initial entry or on re-entry with IFLAG=2, DIAG
!             it must be set by the user to contain the values of the 
!             diagonal matrix Hk0.  Restriction: all elements of DIAG
!             must be positive.
! 
!     IPRINT  is an INTEGER array of length two which must be set by the
!             user.
! 
!             IPRINT(1) specifies the frequency of the output:
!                IPRINT(1) < 0 : no output is generated,
!                IPRINT(1) = 0 : output only at first and last iteration,
!                IPRINT(1) > 0 : output every IPRINT(1) iterations.
! 
!             IPRINT(2) specifies the type of output generated:
!                IPRINT(2) = 0 : iteration count, number of function 
!                                evaluations, function value, norm of the
!                                gradient, and steplength,
!                IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
!                                variables and  gradient vector at the
!                                initial point,
!                IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
!                                variables,
!                IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.
! 
! 
!     EPS     is a positive DOUBLE PRECISION variable that must be set by
!             the user, and determines the accuracy with which the solution
!             is to be found. The subroutine terminates when
!
!                         ||G|| < EPS max(1,||X||),
!
!             where ||.|| denotes the Euclidean norm.
! 
!     XTOL    is a  positive DOUBLE PRECISION variable that must be set by
!             the user to an estimate of the machine precision (e.g.
!             10**(-16) on a SUN station 3/60). The line search routine will
!             terminate if the relative width of the interval of uncertainty
!             is less than XTOL.
! 
!     W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as
!             workspace for LBFGS. This array must not be altered by the
!             user.
! 
!     IFLAG   is an INTEGER variable that must be set to 0 on initial entry
!             to the subroutine. A return with IFLAG<0 indicates an error,
!             and IFLAG=0 indicates that the routine has terminated without
!             detecting errors. On a return with IFLAG=1, the user must
!             evaluate the function F and gradient G. On a return with
!             IFLAG=2, the user must provide the diagonal matrix Hk0.
! 
!             The following negative values of IFLAG, detecting an error,
!             are possible:
! 
!              IFLAG=-1  The line search routine MCSRCH failed. The
!                        parameter INFO provides more detailed information
!                        (see also the documentation of MCSRCH):
!
!                       INFO = 0  IMPROPER INPUT PARAMETERS.
!
!                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
!                                 UNCERTAINTY IS AT MOST XTOL.
!
!                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
!                                 REQUIRED AT THE PRESENT ITERATION.
!
!                       INFO = 4  THE STEP IS TOO SMALL.
!
!                       INFO = 5  THE STEP IS TOO LARGE.
!
!                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. 
!                                 THERE MAY NOT BE A STEP WHICH SATISFIES
!                                 THE SUFFICIENT DECREASE AND CURVATURE
!                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL.
!
! 
!              IFLAG=-2  The i-th diagonal element of the diagonal inverse
!                        Hessian approximation, given in DIAG, is not
!                        positive.
!           
!              IFLAG=-3  Improper input parameters for LBFGS (N or M are
!                        not positive).
! 
!
!
!    ON THE DRIVER:
!
!    The program that calls LBFGS must contain the declaration:
!
!                       EXTERNAL LB2
!
!    LB2 is a BLOCK DATA that defines the default values of several
!    parameters described in the COMMON section. 
!
! 
! 
!    COMMON:
! 
!     The subroutine contains one common area, which the user may wish to
!    reference:
! 
!         COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
! 
!    MP  is an INTEGER variable with default value 6. It is used as the
!        unit number for the printing of the monitoring information
!        controlled by IPRINT.
! 
!    LP  is an INTEGER variable with default value 6. It is used as the
!        unit number for the printing of error messages. This printing
!        may be suppressed by setting LP to a non-positive value.
! 
!    GTOL is a DOUBLE PRECISION variable with default value 0.9, which
!        controls the accuracy of the line search routine MCSRCH. If the
!        function and gradient evaluations are inexpensive with respect
!        to the cost of the iteration (which is sometimes the case when
!        solving very large problems) it may be advantageous to set GTOL
!        to a small value. A typical small value is 0.1.  Restriction:
!        GTOL should be greater than 1.D-04.
! 
!    STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which
!        specify lower and uper bounds for the step in the line search.
!        Their default values are 1.D-20 and 1.D+20, respectively. These
!        values need not be modified unless the exponents are too large
!        for the machine being used, or unless the problem is extremely
!        badly scaled (in which case the exponents should be increased).
! 
!
!  MACHINE DEPENDENCIES
!
!        The only variables that are machine-dependent are XTOL,
!        STPMIN and STPMAX.
! 
!
!  GENERAL INFORMATION
! 
!    Other routines called directly:  DAXPY, DDOT, LB1, MCSRCH
! 
!    Input/Output  :  No input; diagnostic messages on unit MP and
!                     error messages on unit LP.
! 
! 
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!      DOUBLE PRECISION GTOL,ONE,ZERO,GNORM,DDOT,STP1,FTOL,STPMIN,STPMAX,STP,YS,YY,SQ,YR,BETA,XNORM
      DOUBLE PRECISION ONE,ZERO,GNORM,DDOT,STP1,STP,YS,YY,SQ,YR,BETA,XNORM
!      INTEGER MP,LP,ITER,NFUN,POINT,ISPT,IYPT,MAXFEV,INFO,BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN
      INTEGER ITER,NFUN,POINT,ISPT,IYPT,INFO,BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN
      LOGICAL FINISH
      type(input_variables), intent(inout) :: IN 

!
      SAVE
      DATA ONE,ZERO/1.0D+0,0.0D+0/
!
!     INITIALIZE
!     ----------
!
      IF(IFLAG.EQ.0) GO TO 10
      GO TO (172,100) IFLAG
  10  ITER= 0
      IF(N.LE.0.OR.M.LE.0) GO TO 196
      IF(parmin%GTOL.LE.1.D-04) THEN
        IF(parmin%LP.GT.0.AND.IPROC==0) WRITE(parmin%LP,245)
        parmin%GTOL=9.D-01
      ENDIF
      NFUN= 1
      POINT= 0
      FINISH= .FALSE.
      IF(parmin%DIAGCO) THEN
         DO 30 I=1,N
 30      IF (DIAG(I).LE.ZERO) GO TO 195
      ELSE
         DO 40 I=1,N
 40      DIAG(I)= 1.0D0
      ENDIF
!
!     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
!     ---------------------------------------
!     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
!         OTHER TEMPORARY INFORMATION.
!     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
!     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
!         IN THE FORMULA THAT COMPUTES H*G.
!     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
!         STEPS.
!     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
!         GRADIENT DIFFERENCES.
!
!     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
!     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
!
      ISPT= N+2*M
      IYPT= ISPT+N*M     
      DO 50 I=1,N
 50   W(ISPT+I)= -G(I)*DIAG(I)
      GNORM= DSQRT(DDOT(N,G,1,G,1))
!      STP1= ONE/GNORM
!      STP11=2.d-2/GNORM  !original convention
      STP1=in%betax
!     PARAMETERS FOR LINE SEARCH ROUTINE
!     They are now set at initialization, see module par_driver
!      FTOL= 1.0D-4
!      FTOL= 1.0D-6  
!      MAXFEV= 20
!      MAXFEV= 10
!
!      IF(IPRINT(1).GE.0) CALL LB1(IPROC,PAR,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
      IF(IPRINT(1).GE.0) CALL LB1(IPROC,PARMIN,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
!
!    --------------------
!     MAIN ITERATION LOOP
!    --------------------
!
 80   ITER= ITER+1
      INFO=0
      BOUND=ITER-1
      IF(ITER.EQ.1) GO TO 165
      IF (ITER .GT. M)BOUND=M
!
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
      IF(.NOT.parmin%DIAGCO) THEN
         YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
         DO 90 I=1,N
   90    DIAG(I)= YS/YY
      ELSE
         IFLAG=2
         RETURN
      ENDIF
 100  CONTINUE
      IF(parmin%DIAGCO) THEN
        DO 110 I=1,N
 110    IF (DIAG(I).LE.ZERO) GO TO 195
      ENDIF
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!
      CP= POINT
      IF (POINT.EQ.0) CP=M
      W(N+CP)= ONE/YS
      DO 112 I=1,N
 112  W(I)= -G(I)
      CP= POINT
      DO 125 I= 1,BOUND
         CP=CP-1
         IF (CP.EQ. -1)CP=M-1
         SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
         INMC=N+M+CP+1
         IYCN=IYPT+CP*N
         W(INMC)= W(N+CP+1)*SQ
         CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
 125  CONTINUE
!
      DO 130 I=1,N
 130  W(I)=DIAG(I)*W(I)
!
      DO 145 I=1,BOUND
         YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
         BETA= W(N+CP+1)*YR
         INMC=N+M+CP+1
         BETA= W(INMC)-BETA
         ISCN=ISPT+CP*N
         CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
         CP=CP+1
         IF (CP.EQ.M)CP=0
 145  CONTINUE
!
!     STORE THE NEW SEARCH DIRECTION
!     ------------------------------
!
       DO 160 I=1,N
 160   W(ISPT+POINT*N+I)= W(I)
!
!     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION 
!     BY USING THE LINE SEARCH ROUTINE MCSRCH
!     ----------------------------------------------------
 165  NFEV=0
      STP=ONE
      IF (ITER.EQ.1) STP=STP1
      DO 170 I=1,N
 170  W(I)=G(I)
 172  CONTINUE
!      CALL MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,XTOL,MAXFEV,INFO,NFEV,DIAG)
!      CALL MCSRCH(PAR,N,X,F,G,W(ISPT+POINT*N+1),STP,INFO,NFEV,DIAG)
      CALL MCSRCH(PARMIN,N,X,F,G,W(ISPT+POINT*N+1),STP,INFO,NFEV,DIAG)
      IF (INFO .EQ. -1) THEN
        IFLAG=1
        RETURN
      ENDIF
      IF (INFO .NE. 1) GO TO 190
      NFUN= NFUN + NFEV
!
!     COMPUTE THE NEW STEP AND GRADIENT CHANGE 
!     -----------------------------------------
!
      NPT=POINT*N
      DO 175 I=1,N
      W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
 175  W(IYPT+NPT+I)= G(I)-W(I)
      POINT=POINT+1
      IF (POINT.EQ.M)POINT=0
!
!     TERMINATION TEST
!     ----------------
!
      GNORM= DSQRT(DDOT(N,G,1,G,1))
      XNORM= DSQRT(DDOT(N,X,1,X,1))
      XNORM= DMAX1(1.0D0,XNORM)
!      IF (GNORM/XNORM .LE. EPS) FINISH=.TRUE.
      IF (GNORM .LE. EPS) FINISH=.TRUE.
!
!      IF(IPRINT(1).GE.0) CALL LB1(IPROC,PAR,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
      IF(IPRINT(1).GE.0) CALL LB1(IPROC,PARMIN,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
      IF (FINISH) THEN
         IFLAG=0
         RETURN
      ENDIF
      GO TO 80
!
!     ------------------------------------------------------------
!     END OF MAIN ITERATION LOOP. ERROR EXITS.
!     ------------------------------------------------------------
!
 190  IFLAG=-1
      IF(parmin%LP.GT.0.AND.IPROC==0) WRITE(parmin%LP,200) INFO
      RETURN
 195  IFLAG=-2
      IF(parmin%LP.GT.0.AND.IPROC==0) WRITE(parmin%LP,235) I
      RETURN
 196  IFLAG= -3
      IF(parmin%LP.GT.0.AND.IPROC==0) WRITE(parmin%LP,240)
!
!     FORMATS
!     -------
!
 200  FORMAT(/' IFLAG= -1 ',/' LINE SEARCH FAILED. SEE',/&
               ' DOCUMENTATION OF ROUTINE MCSRCH',/' ERROR RETURN',/&
               ' OF LINE SEARCH: INFO= ',I2,/&
               ' POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT',/&
               ' OR INCORRECT TOLERANCES')
 235  FORMAT(/' IFLAG= -2',/' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,&
            ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
 240  FORMAT(/' IFLAG= -3',/' IMPROPER INPUT PARAMETERS (N OR M',&
            ' ARE NOT POSITIVE)')
 245  FORMAT(/'  GTOL IS LESS THAN OR EQUAL TO 1.D-04',&
            / ' IT HAS BEEN RESET TO 9.D-01')
      RETURN
      END
!
!     LAST LINE OF SUBROUTINE LBFGS
!
!
!      SUBROUTINE LB1(IPROC,PAR,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
      SUBROUTINE LB1(IPROC,PARMIN,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
!      use par_driver , only:driverparameters
      use minpar, only: parameterminimization
!
!     -------------------------------------------------------------
!     THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND
!     AMOUNT OF OUTPUT ARE CONTROLLED BY IPRINT.
!     -------------------------------------------------------------
      IMPLICIT NONE
!      type(driverparameters)::par
      type(parameterminimization)::parmin
!      INTEGER IPRINT(2),ITER,NFUN,LP,MP,N,M,IPROC
      INTEGER IPRINT(2),ITER,NFUN,N,M,IPROC,I
      DOUBLE PRECISION X(N),G(N),F,GNORM,STP!,GTOL,STPMIN,STPMAX
      LOGICAL FINISH
!      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
!
      IF (ITER.EQ.0)THEN
!      IF(IPROC==0)     WRITE(parmin%MP,10)
      IF(IPROC==0)     WRITE(parmin%MP,20) N,M
      IF(IPROC==0)     WRITE(parmin%MP,30)F,GNORM
                 IF (IPRINT(2).GE.1)THEN
                IF(IPROC==0)     WRITE(parmin%MP,40)
                IF(IPROC==0)     WRITE(parmin%MP,50) (X(I),I=1,N)
                IF(IPROC==0)     WRITE(parmin%MP,60)
                IF(IPROC==0)     WRITE(parmin%MP,50) (G(I),I=1,N)
                  ENDIF
!      IF(IPROC==0)     WRITE(parmin%MP,10)
!      IF(IPROC==0)     WRITE(parmin%MP,70)
      ELSE
          IF ((IPRINT(1).EQ.0).AND.(ITER.NE.1.AND..NOT.FINISH))RETURN
              IF (IPRINT(1).NE.0)THEN
                   IF(MOD(ITER-1,IPRINT(1)).EQ.0.OR.FINISH)THEN
!                         IF(IPRINT(2).GT.1.AND.ITER.GT.1.AND.IPROC==0) WRITE(parmin%MP,70)
!                         IF(IPROC==0)   WRITE(parmin%MP,80)ITER,NFUN,F,GNORM,STP
!                         IF(IPROC==0)   write(parmin%MP,'(i5,1x,e12.5,1x,e21.14,a,i5,1x,e12.5)') ITER,gnorm,f,' GEOPT BFGS ', nfun,STP
                         IF(IPROC==0)   write(parmin%MP,'(a,2(i5,1x),1x,1pe21.14,1x,1pe12.5,1x,1E12.5)')  & 
                                        ' MIN ',ITER,nfun,f,gnorm,STP
                         parmin%FINSTEP=ITER
                         parmin%ALPHA=STP
                   parmin%IWRITE=.true.
                   ELSE
                         RETURN
                   ENDIF
              ELSE
                   IF( IPRINT(2).GT.1.AND.FINISH.AND.IPROC==0) WRITE(parmin%MP,70)
                   IF(IPROC==0)   WRITE(parmin%MP,80)ITER,NFUN,F,GNORM,STP

              ENDIF
              IF (IPRINT(2).EQ.2.OR.IPRINT(2).EQ.3)THEN
                    IF (FINISH)THEN
                    IF(IPROC==0)    WRITE(parmin%MP,90)
                    ELSE
                    IF(IPROC==0)    WRITE(parmin%MP,40)
                    ENDIF
                    IF(IPROC==0)  WRITE(parmin%MP,50)(X(I),I=1,N)
                  IF (IPRINT(2).EQ.3)THEN
                    IF(IPROC==0)  WRITE(parmin%MP,60)
                    IF(IPROC==0)  WRITE(parmin%MP,50)(G(I),I=1,N)
                  ENDIF
              ENDIF
            IF (FINISH) WRITE(parmin%MP,100)
      ENDIF
!
 20   FORMAT('  N=',I5,'   NUMBER OF CORRECTIONS=',I2,/,  '       INITIAL VALUES')
 30   FORMAT('  F= ',1PD10.3,'   GNORM= ',1PD10.3)
 40   FORMAT(' VECTOR X= ')
 50   FORMAT(6(2X,1PD10.3))
 60   FORMAT(' GRADIENT VECTOR G= ')
 70   FORMAT(/'   I   NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)
 80   FORMAT(2(I4,1X),3X,3(1PD10.3,2X))
 90   FORMAT(' FINAL POINT X= ')
 100  FORMAT(/' THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.',/' IFLAG = 0')
!
      RETURN
      END
!     ******
!
!
!   ----------------------------------------------------------
!     DATA 
!   ----------------------------------------------------------
!
      BLOCK DATA LB2
      INTEGER LP,MP
      DOUBLE PRECISION GTOL,STPMIN,STPMAX
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
!      DATA MP,LP,GTOL,STPMIN,STPMAX/16,16,9.0D-01,1.0D-20,1.0D+20/
      DATA MP,LP,GTOL,STPMIN,STPMAX/16,16,9.0D-01,1.0D-20,20.d0/
      END
!
!
!!   ----------------------------------------------------------
!!
!      subroutine daxpy(n,da,dx,incx,dy,incy)
!!
!!     constant times a vector plus a vector.
!!     uses unrolled loops for increments equal to one.
!!     jack dongarra, linpack, 3/11/78.
!!
!      double precision dx(1),dy(1),da
!      integer i,incx,incy,ix,iy,m,mp1,n
!!
!      if(n.le.0)return
!      if (da .eq. 0.0d0) return
!      if(incx.eq.1.and.incy.eq.1)go to 20
!!
!!        code for unequal increments or equal increments
!!          not equal to 1
!!
!      ix = 1
!      iy = 1
!      if(incx.lt.0)ix = (-n+1)*incx + 1
!      if(incy.lt.0)iy = (-n+1)*incy + 1
!      do 10 i = 1,n
!        dy(iy) = dy(iy) + da*dx(ix)
!        ix = ix + incx
!        iy = iy + incy
!   10 continue
!      return
!!
!!        code for both increments equal to 1
!!        clean-up loop
!!
!   20 m = mod(n,4)
!      if( m .eq. 0 ) go to 40
!      do 30 i = 1,m
!        dy(i) = dy(i) + da*dx(i)
!   30 continue
!      if( n .lt. 4 ) return
!   40 mp1 = m + 1
!      do 50 i = mp1,n,4
!        dy(i) = dy(i) + da*dx(i)
!        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
!        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
!        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
!   50 continue
!      return
!      end
!
!
!!   ----------------------------------------------------------
!!
!      double precision function ddot(n,dx,incx,dy,incy)
!!
!!     forms the dot product of two vectors.
!!     uses unrolled loops for increments equal to one.
!!     jack dongarra, linpack, 3/11/78.
!!
!      double precision dx(1),dy(1),dtemp
!      integer i,incx,incy,ix,iy,m,mp1,n
!!
!      ddot = 0.0d0
!      dtemp = 0.0d0
!      if(n.le.0)return
!      if(incx.eq.1.and.incy.eq.1)go to 20
!!
!!        code for unequal increments or equal increments
!!          not equal to 1
!!
!      ix = 1
!      iy = 1
!      if(incx.lt.0)ix = (-n+1)*incx + 1
!      if(incy.lt.0)iy = (-n+1)*incy + 1
!      do 10 i = 1,n
!        dtemp = dtemp + dx(ix)*dy(iy)
!        ix = ix + incx
!        iy = iy + incy
!   10 continue
!      ddot = dtemp
!      return
!!
!!        code for both increments equal to 1
!!        clean-up loop
!!
!   20 m = mod(n,5)
!      if( m .eq. 0 ) go to 40
!      do 30 i = 1,m
!        dtemp = dtemp + dx(i)*dy(i)
!   30 continue
!      if( n .lt. 5 ) go to 60
!   40 mp1 = m + 1
!      do 50 i = mp1,n,5
!        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
!   50 continue
!   60 ddot = dtemp
!      return
!      end
!    ------------------------------------------------------------------
!
!     **************************
!     LINE SEARCH ROUTINE MCSRCH
!     **************************
!
!      SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL,MAXFEV,INFO,NFEV,WA)
!      SUBROUTINE MCSRCH(PAR,N,X,F,G,S,STP,INFO,NFEV,WA)
      SUBROUTINE MCSRCH(PARMIN,N,X,F,G,S,STP,INFO,NFEV,WA)
!      use par_driver , only:driverparameters
      use minpar, only: parameterminimization
      IMPLICIT NONE
!      type(driverparameters)::par
      type(parameterminimization)::parmin
      INTEGER N,INFO,NFEV!,MAXFEV
      DOUBLE PRECISION F,STP!,FTOL,GTOL,XTOL,STPMIN,STPMAX
      DOUBLE PRECISION X(N),G(N),S(N),WA(N)
!      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      SAVE
!
!                     SUBROUTINE MCSRCH
!                
!     A slight modification of the subroutine CSRCH of More' and Thuente.
!     The changes are to allow reverse communication, and do not affect
!     the performance of the routine. 
!
!     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
!     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
!
!     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
!     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
!     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
!     MINIMIZER OF THE MODIFIED FUNCTION
!
!          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
!
!     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
!     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
!     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
!     CONTAINS A MINIMIZER OF F(X+STP*S).
!
!     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
!     THE SUFFICIENT DECREASE CONDITION
!
!           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),
!
!     AND THE CURVATURE CONDITION
!
!           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).
!
!     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
!     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
!     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
!     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
!     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
!     SATISFIES THE SUFFICIENT DECREASE CONDITION.
!
!     THE SUBROUTINE STATEMENT IS
!
!        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF VARIABLES.
!
!       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
!         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
!         X + STP*S.
!
!       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
!         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
!
!       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
!         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
!         OF F AT X + STP*S.
!
!       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
!         SEARCH DIRECTION.
!
!       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
!         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
!         STP CONTAINS THE FINAL ESTIMATE.
!
!       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
!         communication implementation GTOL is defined in a COMMON
!         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
!         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
!         SATISFIED.
!
!       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
!         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
!         IS AT MOST XTOL.
!
!       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
!         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
!         communication implementatin they are defined in a COMMON
!         statement).
!
!       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
!         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
!         MAXFEV BY THE END OF AN ITERATION.
!
!       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
!
!         INFO = 0  IMPROPER INPUT PARAMETERS.
!
!         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
!
!         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
!                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
!
!         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
!                   IS AT MOST XTOL.
!
!         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
!
!         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
!
!         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
!
!         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
!                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
!                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
!                   TOLERANCES MAY BE TOO SMALL.
!
!       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
!         CALLS TO FCN.
!
!       WA IS A WORK ARRAY OF LENGTH N.
!
!     SUBPROGRAMS CALLED
!
!       MCSTEP
!
!       FORTRAN-SUPPLIED...ABS,MAX,MIN
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!     JORGE J. MORE', DAVID J. THUENTE
!
!     **********
      INTEGER INFOC,J
      LOGICAL BRACKT,STAGE1
      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,&
            FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,&
            STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO
      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/
      IF(INFO.EQ.-1) GO TO 45
      INFOC = 1
!
!     CHECK THE INPUT PARAMETERS FOR ERRORS.
!
      IF (N .LE. 0 .OR. STP .LE. ZERO .OR. parmin%FTOL .LT. ZERO .OR.&
         parmin%GTOL .LT. ZERO .OR. parmin%XTOL .LT. ZERO .OR. parmin%STPMIN .LT. ZERO &
         .OR. parmin%STPMAX .LT. parmin%STPMIN .OR. parmin%MAXFEV .LE. 0) RETURN
!
!     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
!     AND CHECK THAT S IS A DESCENT DIRECTION.
!
      DGINIT = ZERO
      DO 10 J = 1, N
         DGINIT = DGINIT + G(J)*S(J)
   10    CONTINUE
      IF (DGINIT .GE. ZERO) then
         write(parmin%LP,15)
   15    FORMAT(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
         RETURN
         ENDIF
!
!     INITIALIZE LOCAL VARIABLES.
!
      BRACKT = .FALSE.
      STAGE1 = .TRUE.
      NFEV = 0
      FINIT = F
      DGTEST = parmin%FTOL*DGINIT
      WIDTH = parmin%STPMAX - parmin%STPMIN
      WIDTH1 = WIDTH/P5
      DO 20 J = 1, N
         WA(J) = X(J)
   20    CONTINUE
!
!     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
!     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
!     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
!     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
!     THE INTERVAL OF UNCERTAINTY.
!     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
!     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
!
      STX = ZERO
      FX = FINIT
      DGX = DGINIT
      STY = ZERO
      FY = FINIT
      DGY = DGINIT
!
!     START OF ITERATION.
!
   30 CONTINUE
!
!        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
!        TO THE PRESENT INTERVAL OF UNCERTAINTY.
!
         IF (BRACKT) THEN
            STMIN = MIN(STX,STY)
            STMAX = MAX(STX,STY)
         ELSE
            STMIN = STX
            STMAX = STP + XTRAPF*(STP - STX)
            END IF
!
!        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
!
         STP = MAX(STP,parmin%STPMIN)
         STP = MIN(STP,parmin%STPMAX)
!
!        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
!        STP BE THE LOWEST POINT OBTAINED SO FAR.
!
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))&
           .OR. NFEV .GE. parmin%MAXFEV-1 .OR. INFOC .EQ. 0&
           .OR. (BRACKT .AND. STMAX-STMIN .LE. parmin%XTOL*STMAX)) STP = STX
!
!        EVALUATE THE FUNCTION AND GRADIENT AT STP
!        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
!        We return to main program to obtain F and G.
!
         DO 40 J = 1, N
            X(J) = WA(J) + STP*S(J)
   40       CONTINUE
         INFO=-1
         RETURN
!
   45    INFO=0
         NFEV = NFEV + 1
         DG = ZERO
         DO 50 J = 1, N
            DG = DG + G(J)*S(J)
   50       CONTINUE
         FTEST1 = FINIT + STP*DGTEST
!
!        TEST FOR CONVERGENCE.
!
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))&
           .OR. INFOC .EQ. 0) INFO = 6
         IF (STP .EQ. parmin%STPMAX .AND.&
            F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5
         IF (STP .EQ. parmin%STPMIN .AND.&
            (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4
         IF (NFEV .GE. parmin%MAXFEV) INFO = 3
         IF (BRACKT .AND. STMAX-STMIN .LE. parmin%XTOL*STMAX) INFO = 2
         IF (F .LE. FTEST1 .AND. ABS(DG) .LE. parmin%GTOL*(-DGINIT)) INFO = 1
!
!        CHECK FOR TERMINATION.
!
         IF (INFO .NE. 0) RETURN
!
!        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
!        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
!
         IF (STAGE1 .AND. F .LE. FTEST1 .AND.&
            DG .GE. MIN(parmin%FTOL,parmin%GTOL)*DGINIT) STAGE1 = .FALSE.
!
!        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
!        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
!        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
!        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
!        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
!
         IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN
!
!           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
!
            FM = F - STP*DGTEST
            FXM = FX - STX*DGTEST
            FYM = FY - STY*DGTEST
            DGM = DG - DGTEST
            DGXM = DGX - DGTEST
            DGYM = DGY - DGTEST
!
!           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!           AND TO COMPUTE THE NEW STEP.
!
            CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,&
                      BRACKT,STMIN,STMAX,INFOC)
!
!           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
!
            FX = FXM + STX*DGTEST
            FY = FYM + STY*DGTEST
            DGX = DGXM + DGTEST
            DGY = DGYM + DGTEST
         ELSE
!
!           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!           AND TO COMPUTE THE NEW STEP.
!
            CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,&
                      BRACKT,STMIN,STMAX,INFOC)
            END IF
!
!        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
!        INTERVAL OF UNCERTAINTY.
!
         IF (BRACKT) THEN
            IF (ABS(STY-STX) .GE. P66*WIDTH1)&
              STP = STX + P5*(STY - STX)
            WIDTH1 = WIDTH
            WIDTH = ABS(STY-STX)
            END IF
!
!        END OF ITERATION.
!
         GO TO 30
!
!     LAST LINE OF SUBROUTINE MCSRCH.
!
      END
      SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,STPMIN,STPMAX,INFO)
      IMPLICIT NONE
      INTEGER INFO
      DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
      LOGICAL BRACKT,BOUND
!
!     SUBROUTINE MCSTEP
!
!     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
!     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
!     A MINIMIZER OF THE FUNCTION.
!
!     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
!     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
!     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
!     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
!     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
!     WITH ENDPOINTS STX AND STY.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
!                        STPMIN,STPMAX,INFO)
!
!     WHERE
!
!       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
!         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
!         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
!         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
!         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
!
!       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
!         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
!         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
!         UPDATED APPROPRIATELY.
!
!       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
!         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
!         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
!         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
!
!       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
!         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
!         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
!         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
!
!       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
!         AND UPPER BOUNDS FOR THE STEP.
!
!       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
!         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
!         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
!         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!     JORGE J. MORE', DAVID J. THUENTE
!
      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
      INFO = 0
!
!     CHECK THE INPUT PARAMETERS FOR ERRORS.
!
      IF ((BRACKT .AND. (STP .LE. MIN(STX,STY) .OR.&
         STP .GE. MAX(STX,STY))) .OR.&
         DX*(STP-STX) .GE. 0.0 .OR. STPMAX .LT. STPMIN) RETURN
!
!     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
!
      SGND = DP*(DX/ABS(DX))
!
!     FIRST CASE. A HIGHER FUNCTION VALUE.
!     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
!     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
!     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
!
      IF (FP .GT. FX) THEN
         INFO = 1
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .LT. STX) GAMMA = -GAMMA
         P = (GAMMA - DX) + THETA
         Q = ((GAMMA - DX) + GAMMA) + DP
         R = P/Q
         STPC = STX + R*(STP - STX)
         STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
         IF (ABS(STPC-STX) .LT. ABS(STPQ-STX)) THEN
            STPF = STPC
         ELSE
           STPF = STPC + (STPQ - STPC)/2
           END IF
         BRACKT = .TRUE.
!
!     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
!     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
!     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
!     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
!
      ELSE IF (SGND .LT. 0.0) THEN
         INFO = 2
         BOUND = .FALSE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = ((GAMMA - DP) + GAMMA) + DX
         R = P/Q
         STPC = STP + R*(STX - STP)
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
            STPF = STPC
         ELSE
            STPF = STPQ
            END IF
         BRACKT = .TRUE.
!
!     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
!     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
!     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
!     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
!     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
!     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
!     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
!     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
!
      ELSE IF (ABS(DP) .LT. ABS(DX)) THEN
         INFO = 3
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
!
!        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
!        TO INFINITY IN THE DIRECTION OF THE STEP.
!
         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = (GAMMA + (DX - DP)) + GAMMA
         R = P/Q
         IF (R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN
            STPC = STP + R*(STX - STP)
         ELSE IF (STP .GT. STX) THEN
            STPC = STPMAX
         ELSE
            STPC = STPMIN
            END IF
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (BRACKT) THEN
            IF (ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
               END IF
         ELSE
            IF (ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
               END IF
            END IF
!
!     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
!     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
!     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
!     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
!
      ELSE
         INFO = 4
         BOUND = .FALSE.
         IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP .GT. STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
         ELSE IF (STP .GT. STX) THEN
            STPF = STPMAX
         ELSE
            STPF = STPMIN
            END IF
         END IF
!
!     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
!     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
!
      IF (FP .GT. FX) THEN
         STY = STP
         FY = FP
         DY = DP
      ELSE
         IF (SGND .LT. 0.0) THEN
            STY = STX
            FY = FX
            DY = DX
            END IF
         STX = STP
         FX = FP
         DX = DP
         END IF
!
!     COMPUTE THE NEW STEP AND SAFEGUARD IT.
!
      STPF = MIN(STPMAX,STPF)
      STPF = MAX(STPMIN,STPF)
      STP = STPF
      IF (BRACKT .AND. BOUND) THEN
         IF (STY .GT. STX) THEN
            STP = MIN(STX+0.66*(STY-STX),STP)
         ELSE
            STP = MAX(STX+0.66*(STY-STX),STP)
            END IF
         END IF
      RETURN
!
!     LAST LINE OF SUBROUTINE MCSTEP.
!
      END
