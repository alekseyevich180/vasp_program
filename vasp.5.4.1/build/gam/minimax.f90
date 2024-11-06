# 1 "minimax.F"
# 1 "./symbol.inc" 1 
!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define 1            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define 1             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define REAL_to_DBLE        convert REAL() to DBLE()
!#define 1                 compile for parallel machine with 1
!------------- end of user part --------------------------------
# 59


!
!   charge density: half grid mode Z direction
!










# 101


!
!   charge density real
!




# 118



!
!   wavefunctions: half grid mode for Z direction
!

# 134


!
!   wavefunctions real (gamma only)
!




























# 198

!
!   common definitions
!







# 211


!
!   mpi parallel macros
!














# 248

# 254

# 259

# 266

# 273



# 282







# 297







# 306


# 319










# 336

# 2 "minimax.F" 2 
!*****************************************************************
! This module contains routines to calculate the fourier transform
! of non-equally spaced grid points. The idea behind this is to
! apply Chebyshev's minimax theorem and applying the corresponding
! Remez algorithm.
!*****************************************************************

MODULE minimax
   USE prec

   INTEGER, PARAMETER :: qd=16
   REAL(qd),PARAMETER :: PIQ=3.1415926535897932384626433832795028841971693993751_qd
!parameter for Remez algorithm
   INTEGER, SAVE      :: NMINALT = 50
   INTEGER, SAVE      :: NMAXALT = 100
   INTEGER, SAVE      :: NALTCON = 200
!tweak minimization interval ?
   LOGICAL, SAVE      :: LTWEAK=.TRUE.

!some private parameter
!accuracy in minimization routines
   REAL(qd), PRIVATE, PARAMETER :: ACC=1.E-16_qd  
!truncation threshold for SVD and accuracy in FT routine
   REAL(qd),PARAMETER  :: ACC2=1.E-12_qd   
!tiny number for QDLUDCMP
   REAL (qd), PARAMETER :: VERYTINY=1.E-20_qd
!maximum iterations in Newton-Raphson algorithm
   INTEGER, PRIVATE, PARAMETER :: MAXNEWTON = 300
!damping factor for Newton-Raphson root finder
   REAL(qd), PRIVATE, PARAMETER :: DAMPNEWT=0.35_qd
!maximum iterations in Newton-Root finder
   INTEGER, PRIVATE, PARAMETER :: MAXNEWT = 100
!maximum iterations in Remez algorithm
   INTEGER, PRIVATE, SAVE :: MAXREMEZ = 50
!maximum bisections for FINDMAX and FINDMAX_LINEAR
   INTEGER, PRIVATE, PARAMETER :: MAXBISEC=5000
!parameter for LSQ search algorithm
   INTEGER, PARAMETER :: ITMAXLSQ = 1500
   INTEGER, PARAMETER :: NDATALSQ = 100
   REAL(qd), PARAMETER :: RANDLSQ =  50._qd
   INTEGER, PARAMETER :: IGUESS_MAX = 10

CONTAINS

!---------------------------------------------------------------
! function which should be interpolated
! with the minimax contition
!---------------------------------------------------------------

   FUNCTION INTERPOL(X)
      USE prec
      IMPLICIT NONE
      REAL(qd) X
      REAL(qd) INTERPOL
       
      IF (ABS(X)<1E-9_qd) THEN
         WRITE(*,*) ' INTERPOL: X too close to 0, returning 0',X
         INTERPOL=0._qd
      ELSE
         INTERPOL=1._qd/X
         ENDIF
   END FUNCTION INTERPOL

!----------------------LOGPARTITION-----------------------------
! logarithmic partition of inter a,b
!---------------------------------------------------------------

   FUNCTION LOGPARTITION(I,N,A,B)
      USE prec 
      IMPLICIT NONE
      REAL(qd) LOGPARTITION
      REAL(qd) A,B
      INTEGER I,N
      
      IF ( N<=1 ) THEN
         WRITE(*,*)'ERROR in LOGPARTITION, too less partitions chosen',N
         CALL M_exit(); stop
      ENDIF

      IF (I<1) THEN
         WRITE(*,*)'ERROR in LOGPARTITION, I range must be between 1 and N',I,N
         CALL M_exit(); stop
      ENDIF
      
      LOGPARTITION=(B*Log(2._qd)*Log(2._qd-I+N)+Log(1._qd+N)*((A-B)*Log(2._qd)-A*Log(2._qd-I+N)))/&
         ((Log(2._qd)-Log(1._qd+N))*Log(2._qd-I+N))
 
     RETURN
  ENDFUNCTION LOGPARTITION 

!****************************************************************************
! Main routine, caluclates imaginary time and frequency grid
! and corresponding forward Fouier transformation matrix,i.e.
! transformation from imaginary time to imanginary frequency
! ocptionally the sin transform is calculated
!****************************************************************************

   SUBROUTINE CALCULATE_MINIMAX_GRIDS(X1, X2, NTAU, TAU, TAUWEIGHT, NNU, &
      NU_COS, NU_COSWEIGHT, FTCOS, MAXITER_FT, IO, TYP, LFREQ, LTIME, ERRORS, &
      NU_SIN, NU_SINWEIGHT, FTSIN )
      USE prec 
      USE base
      IMPLICIT NONE
      REAL(q)             :: X1, X2            !energy interval
      INTEGER             :: NTAU              !number of time points
      REAL(q)             :: TAU(NTAU)         !imaginary times
      REAL(q)             :: TAUWEIGHT(NTAU)   !imaginary time weight
      INTEGER             :: NNU               !number of frequency points
      REAL(q)             :: NU_COS(NNU)       !imaginary frequencies (cos related)
      REAL(q)             :: NU_COSWEIGHT(NNU) !imaginary frequency weight (cos related)
      REAL(q)             :: FTCOS(NNU,NTAU)   !Fourier transformation matrix (cos related)
      INTEGER             :: MAXITER_FT        !iterations in sloppy remez algorithm
      TYPE(in_struct)     :: IO                !in_struct
      INTEGER             :: TYP               !type of quadrature (lsq or mm)
      REAL(q)             :: ERRORS(3+2*NNU)   !errors of approximations
      LOGICAL             :: LFREQ, LTIME      !flags
      REAL(q),OPTIONAL    :: NU_SIN(NNU)       !imaginary frequencies  (sin related)
      REAL(q),OPTIONAL    :: NU_SINWEIGHT(NNU) !imaginary frequency weight (sin related)
      REAL(q),OPTIONAL    :: FTSIN(NNU,NTAU)   !Fourier transformation matrix (sin related)
!local
      REAL(qd)            :: TIMES(2*NTAU)
      REAL(qd)            :: FREQ(2*NNU)
      REAL(qd),ALLOCATABLE:: FREQ_SIN(:) 
      REAL(qd)            :: TRANSFORMATION(NNU,NTAU)
      REAL(q)             :: R,A,B
      INTEGER             :: I,J,ITER
      REAL(q)             :: WORKFT(NNU)
      REAL(q)             :: ERR, ERRS(NNU)
      INTEGER             :: IPIV(NNU),INF
      LOGICAL             :: LSIN
      REAL(q)             :: RMAX
      INTEGER             :: INU
      INTEGER             :: ICOS, ISIN
   
      INU=IO%IU6
!check if Sin transform should be calculated
      IF (PRESENT(NU_SIN).AND.PRESENT(NU_SINWEIGHT).AND.PRESENT(FTSIN) )THEN
         LSIN=.TRUE.
         ALLOCATE(FREQ_SIN(2*NNU))
      ELSE
         LSIN=.FALSE.
      ENDIF
   
! first determine R
      R=X2/X1
!find maximum R, depending on NNU
      CALL FIND_RNMAX(R,NNU,IO%IU0)  
      A=MAX(0.01_q,X1)
      B=R*A

!set MAXREMEZ
      MAXREMEZ=MAXITER_FT
     
      IF ( NNU /= NTAU) THEN 
         IF (INU>=0 ) WRITE(*,*) 'Currently only 1-1 Fourier transform implemented, Sorry'
         NNU =-1
      ENDIF
     
!find out which grids to compute
      IF ( TYP == 141 ) THEN  
! ABS grid
         ICOS=4
         ISIN=4
      ELSE IF ( TYP == 142 ) THEN
! SIN grid
         ICOS=3
         ISIN=3
      ELSE IF (TYP == 143 ) THEN 
! COS + SIN grid
         ICOS=2
         ISIN=3
      ELSE
! COS grid
         ICOS=2
         ISIN=2
      ENDIF
      
      IF ( NNU > 32) THEN
         NNU=20
         NTAU=20          !default value for high accuracy
         IF ( INU >= 0 ) WRITE(*,'(" Number of grid points forced to ",I4)'),NNU
      ENDIF
   
!determine, which coefficients to use
      IF (INU>=0) WRITE(INU,'(" Interval for error minimization: [",E10.3,",",E10.3,"]")')A,B
      IF (INU>=0) WRITE(*,'(" Interval for error minimization: [",E10.3,",",E10.3,"]")')A,B
!decide if we start from a higher RMAX(N) and converge to actual R
      CALL RMAX2CURRENT_R(R,NNU,RMAX,IO%IU0)
      IF ( .NOT.LTWEAK) THEN
         R=RMAX
         B=RMAX*A
         IF (INU>=0) WRITE(INU,'(" Used interval: [",E10.3,",",E10.3,"]")')A,B
         IF (INU>=0) WRITE(*,'(" Used interval: [",E10.3,",",E10.3,"]")')A,B
      ENDIF

!frequency grid
      IF (LFREQ) THEN
!cos-nu points
!lsq coefficients may have to be found by converging from RMAX to R
         IF ( LTWEAK ) THEN
            CALL CONVERGE_FROM_RMAX2R_LSQ(ICOS,RMAX,R,NNU,FREQ,IO)
         ELSE 
            CALL NONLINEAR_LSQ_FIT(ICOS,REAL(R,KIND=qd),NNU,FREQ,IO)
         ENDIF 
!this calculates optimal nu (cos related) points
         CALL REMEZ_NL(ICOS,NNU,REAL(R,KIND=qd),ACC,FREQ,INU,ERR)  
!the scaled coefficients are obtained by multipliciation of A
         ERRORS(1) = ERR/A
         DO I=1,NNU
                NU_COS(I)=ABS(REAL(FREQ(I)*A,q))
          NU_COSWEIGHT(I)=REAL(FREQ(I+NNU)*A,q)
         ENDDO

!sin-nu points, if desired
         IF (LSIN) THEN
!sin only necessary if ISIN==3
            IF ( ISIN==3 ) THEN
!lsq coefficients may have to be found by converging from RMAX to R
               IF ( LTWEAK ) THEN
                  CALL CONVERGE_FROM_RMAX2R_LSQ(ISIN,RMAX,R,NNU,FREQ_SIN,IO)
               ELSE 
                  CALL NONLINEAR_LSQ_FIT(ISIN,REAL(R,KIND=qd),NNU,FREQ_SIN,IO)
               ENDIF 
!this calculates optimal nu (cos related) points
               CALL REMEZ_NL(ISIN,NNU,REAL(R,KIND=qd),ACC,FREQ_SIN,INU,ERR)  
            ELSE 
!otherwise use cos points
               FREQ_SIN(:)=FREQ(:)
            ENDIF
!the scaled coefficients are obtained by multipliciation of A
            ERRORS(2) = ERR/A
            DO I=1,NNU
               NU_SIN(I)=ABS(REAL(FREQ_SIN(I)*A,q))
               NU_SINWEIGHT(I)=REAL(FREQ_SIN(I+NNU)*A,q)
            ENDDO
         ENDIF
      ENDIF
  
!frequency grid
      IF (LTIME) THEN
         IF ( LTWEAK ) THEN
            CALL CONVERGE_FROM_RMAX2R_LSQ(1,RMAX,R,NTAU,TIMES,IO)
         ELSE
            CALL NONLINEAR_LSQ_FIT(1,REAL(R,KIND=qd),NTAU,TIMES,IO)
         ENDIF
         CALL REMEZ_NL(1,NTAU,REAL(R,KIND=qd),ACC,TIMES,INU,ERR)  
!scaled time points obtained by division by 2A
         TIMES=TIMES/2     
         ERRORS(3) = ERR/A 
         DO I=1,NTAU
            TAU(I)=ABS(REAL(TIMES(I)/A,q))
            TAUWEIGHT(I)=REAL(TIMES(I+NTAU)/A,q)
         ENDDO
      ENDIF

!----------------------------------------------------------
! cos fourier matrix necessary only for ACFDTRK, or ACFDTR
      IF (LTIME .AND. (.NOT.LSIN)) THEN
!finally the cos tranformation matrix
         CALL CALC_MINIMAX_FT(.TRUE., FREQ, NNU, TIMES, NTAU, R , &
                               0,TRANSFORMATION,INU,ERRS,A)
!scale onto proper interval
         TRANSFORMATION = TRANSFORMATION/A
         DO I=1,NNU
            DO J=1,NTAU
               FTCOS(I,J)=REAL(TRANSFORMATION(I,J)*COS(FREQ(I)*TIMES(J)),q)
            ENDDO
            ERRORS(3+I)=ERRS(I)
         ENDDO
      ENDIF
   ENDSUBROUTINE CALCULATE_MINIMAX_GRIDS

!****************************************************************************
! Remez algorithm for calculating the minimax approximation
! for a given function INTERPOL in the intervall A,B by a
! sum of weightend nonlinear rationals in NU
! This gives the optimal frequency or time grid
!
! COEFF(2*N) ... on entry : guess for minimax coefficients in [1,R]
! COEFF(2*N) ... on exit  : exact minimax coefficients in [1,R]
! if ITYPE=1 ... exp time grid
! if ITYPE=3 ... sin related coefficients are calculated
! if ITYPE=4 ... abs(=cos^2+sin^2) related coefficients are calculated
! otherwise  ... cos related coefficients are calculated
!
!****************************************************************************

   SUBROUTINE REMEZ_NL(ITYPE,N,R,ACC,COEFF,INU,ERROR,LCONV)
      USE prec
      IMPLICIT NONE
      INTEGER          :: ITYPE      !kind of error function?
      INTEGER          :: N          !order of interpolation
      REAL(qd)         :: R          !right interval boundary
      REAL(qd)         :: ACC        !desired accuracy
      REAL(qd)         :: COEFF(2*N) !coefficients
      INTEGER          :: INU        !writing unit
      REAL(q)          :: ERROR      !maximum error
      LOGICAL,OPTIONAL :: LCONV
!local
      INTEGER               :: I,INFO
      INTEGER               :: ITER
      REAL(qd)              :: Z,ZP(2*N),B,A
      REAL(qd)              :: LAMB(2*N+1),LAMBTEMP(2*N+1)
      REAL(qd)              :: MAXERROR2,MAXERROR1
      REAL(qd)              :: X0(2*N+1)
      REAL(qd), ALLOCATABLE :: RAN(:)
      LOGICAL               :: LCONV_
  
      LCONV_=.FALSE.
      IF ( PRESENT(LCONV))LCONV_=LCONV
   
      IF ( SIZE(COEFF) < 2*N ) THEN
         WRITE(*,*)'Error in REMEZ_NL, Sizes incompatible:',SIZE(COEFF),2*N
         CALL M_exit(); stop
      ENDIF
   
!this is the minimization interval
      A=REAL(1,KIND=qd)
      B=REAL(R,KIND=qd)
      
!write preamble
      IF ( INU >=0 ) THEN 
         WRITE(INU,196)
         WRITE(INU,197)
         IF (ITYPE==1) THEN
            WRITE(INU,201)
         ELSEIF (ITYPE==3) THEN
            WRITE(INU,203)
         ELSEIF(ITYPE==4) THEN
            WRITE(INU,204)
         ELSE
            WRITE(INU,202)
         ENDIF 
      ENDIF 
     
! this is the starting guess, obtained from varpro
      DO I=1, 2*N
         LAMB(I)=COEFF(I)
      ENDDO  
   
! The minimum of the L_\infty norm of the error function is:
      IF (ITYPE==1) THEN
         LAMB(2*N+1)=FINDMAX(N,LAMB,TAU_ERROR_FUNCTION,A,B)   
      ELSEIF (ITYPE==3) THEN
         LAMB(2*N+1)=FINDMAX(N,LAMB,NUSIN_ERROR_FUNCTION,A,B)   
      ELSEIF(ITYPE==4) THEN
         LAMB(2*N+1)=FINDMAX(N,LAMB,NUABS_ERROR_FUNCTION,A,B)   
      ELSE
         LAMB(2*N+1)=FINDMAX(N,LAMB,NUCOS_ERROR_FUNCTION,A,B)   
      ENDIF
   
!==========================================================================
      remez: DO ITER=1, MAXREMEZ  !start the remez algorithm
!==========================================================================
       
         IF(ITYPE==1) THEN 
            MAXERROR1=FINDMAX(N,LAMB,TAU_ERROR_FUNCTION,A,B)
         ELSEIF(ITYPE==3) THEN 
            MAXERROR1=FINDMAX(N,LAMB,NUSIN_ERROR_FUNCTION,A,B)
         ELSEIF(ITYPE==4) THEN
            MAXERROR1=FINDMAX(N,LAMB,NUABS_ERROR_FUNCTION,A,B)
         ELSE
            MAXERROR1=FINDMAX(N,LAMB,NUCOS_ERROR_FUNCTION,A,B)
         ENDIF 
!----------------------------------------------------------------------
! find an alternant of error function ( all local extrema)
         IF ( N < 27 ) THEN
            IF (ITYPE==1) THEN
               CALL FIND_ALTERNANT_EXTREMA(TAU_ERROR_FUNCTION,N,A,B,LAMB,X0,NMINALT,INFO)
               IF (INFO<=-1) CALL FIND_ALTERNANT_EXTREMA(TAU_ERROR_FUNCTION,N,A,B,LAMB,X0,NMAXALT,INFO)
            ELSEIF (ITYPE==3) THEN
               CALL FIND_ALTERNANT_EXTREMA(NUSIN_ERROR_FUNCTION,N,A,B,LAMB,X0,NMINALT,INFO)
               IF (INFO<=-1) CALL FIND_ALTERNANT_EXTREMA(NUSIN_ERROR_FUNCTION,N,A,B,LAMB,X0,NMAXALT,INFO)
            ELSEIF(ITYPE==4) THEN
               CALL FIND_ALTERNANT_EXTREMA(NUABS_ERROR_FUNCTION,N,A,B,LAMB,X0,NMINALT,INFO)
               IF (INFO<=-1) CALL FIND_ALTERNANT_EXTREMA(NUABS_ERROR_FUNCTION,N,A,B,LAMB,X0,NMAXALT,INFO)
            ELSE
               CALL FIND_ALTERNANT_EXTREMA(NUCOS_ERROR_FUNCTION,N,A,B,LAMB,X0,NMINALT,INFO)
               IF (INFO<=-1) CALL FIND_ALTERNANT_EXTREMA(NUCOS_ERROR_FUNCTION,N,A,B,LAMB,X0,NMAXALT,INFO)
            ENDIF
         ELSE
            IF (ITYPE==1) THEN
               CALL FIND_ALTERNANT_EXTREMA(TAU_ERROR_FUNCTION,N,A,B,LAMB,X0,NALTCON,INFO)
            ELSEIF (ITYPE==3) THEN
               CALL FIND_ALTERNANT_EXTREMA(NUSIN_ERROR_FUNCTION,N,A,B,LAMB,X0,NALTCON,INFO)
            ELSEIF(ITYPE==4) THEN
               CALL FIND_ALTERNANT_EXTREMA(NUABS_ERROR_FUNCTION,N,A,B,LAMB,X0,NALTCON,INFO)
            ELSE
               CALL FIND_ALTERNANT_EXTREMA(NUCOS_ERROR_FUNCTION,N,A,B,LAMB,X0,NALTCON,INFO)
            ENDIF
         ENDIF 
         IF (INFO <=-1 ) THEN
           IF (LCONV_)THEN
               IF (ITYPE==1) THEN
                  CALL FIND_ALTERNANT_EXTREMA(TAU_ERROR_FUNCTION,N,A,B,LAMB,X0,NALTCON,INFO)
               ELSEIF (ITYPE==3) THEN
                  CALL FIND_ALTERNANT_EXTREMA(NUSIN_ERROR_FUNCTION,N,A,B,LAMB,X0,NALTCON,INFO)
               ELSEIF(ITYPE==4) THEN
                  CALL FIND_ALTERNANT_EXTREMA(NUABS_ERROR_FUNCTION,N,A,B,LAMB,X0,NALTCON,INFO)
               ELSE
                  CALL FIND_ALTERNANT_EXTREMA(NUCOS_ERROR_FUNCTION,N,A,B,LAMB,X0,NALTCON,INFO)
               ENDIF
               IF (INU>=0 ) WRITE(*,100)NALTCON
               IF (INU>=0 ) WRITE(*,'(F45.32)')X0(:)
               CALL M_exit(); stop
            ELSE
               IF (INU>=0 ) WRITE(*,100)NMAXALT
               IF (INU>=0 ) WRITE(*,'(F45.32)')X0(:)
               CALL M_exit(); stop
            ENDIF
         ENDIF 
!----------------------------------------------------------------------
         
!----------------------------------------------------------------------
! solve non linear system for new coefficients
         LAMBTEMP(1:2*N+1)=LAMB(1:2*N+1)
         IF (ITYPE==1)THEN
            CALL SOLVE_NONLS(TAU_ERROR_FUNCTION_GRAD,N,LAMBTEMP,X0,ACC,INFO)
         ELSEIF (ITYPE==3)THEN
            CALL SOLVE_NONLS(NUSIN_ERROR_FUNCTION_GRAD,N,LAMBTEMP,X0,ACC,INFO)
         ELSEIF(ITYPE==4)THEN
            CALL SOLVE_NONLS(NUABS_ERROR_FUNCTION_GRAD,N,LAMBTEMP,X0,ACC,INFO)
         ELSE
            CALL SOLVE_NONLS(NUCOS_ERROR_FUNCTION_GRAD,N,LAMBTEMP,X0,ACC,INFO)
         ENDIF
!----------------------------------------------------------------------
      
!check sucess
         IF (INFO < 0 ) THEN 
            IF (INU>=0) WRITE(*,*)'ERROR in REMEZ_NL, no solution found for ITER=',ITER
            CALL M_exit(); stop
         ELSEIF ( INFO == 0 ) THEN 
            LAMB(1:2*N+1)=LAMBTEMP(1:2*N+1)
         ENDIF 
           
!now LAMB contains the new set of coefficients
!calculate error difference w.r.t. to previous step
         IF (ITYPE==1) THEN
            MAXERROR2=FINDMAX(N,LAMB,TAU_ERROR_FUNCTION,A,B)
         ELSEIF (ITYPE==3) THEN
            MAXERROR2=FINDMAX(N,LAMB,NUSIN_ERROR_FUNCTION,A,B)
         ELSEIF (ITYPE==4) THEN
            MAXERROR2=FINDMAX(N,LAMB,NUABS_ERROR_FUNCTION,A,B)
         ELSE
            MAXERROR2=FINDMAX(N,LAMB,NUCOS_ERROR_FUNCTION,A,B)
         ENDIF
   
         IF (INU>=0) WRITE(INU,198)ITER,MAXERROR1,ABS(MAXERROR2-MAXERROR1)
!----------------------------------------------------------------------
! if accuracy reached, get out!
!----------------------------------------------------------------------
         IF ( ABS(MAXERROR2-MAXERROR1) < ACC ) THEN
            EXIT remez
         ELSE
            CYCLE remez
         ENDIF
   
!=========================================================================
      ENDDO remez
!=========================================================================
   
!sort coefficients
      COEFF(1:2*N)=LAMB(1:2*N)
      CALL PIKSRT(N,COEFF)         
!determine the maximum error
      IF (ITYPE==1) THEN
         ERROR=FINDMAX(N,LAMB,TAU_ERROR_FUNCTION,A,B)
      ELSEIF (ITYPE==3) THEN
         ERROR=FINDMAX(N,LAMB,NUSIN_ERROR_FUNCTION,A,B)
      ELSEIF(ITYPE==4) THEN
         ERROR=FINDMAX(N,LAMB,NUABS_ERROR_FUNCTION,A,B)
      ELSE
         ERROR=FINDMAX(N,LAMB,NUCOS_ERROR_FUNCTION,A,B)
      ENDIF
       
      IF ( ITER > MAXREMEZ ) THEN
!at this stage the Remez algorithm did not converge after MAXREMEZ steps.
!write an appropriate Warning, but don't terminate the program
         IF( INU >=0 ) WRITE(INU,199)MAXREMEZ,ERROR
         IF( INU >=0 ) WRITE(*,199)MAXREMEZ,ERROR
      ELSE
! success!
         IF ( INU >=0 ) THEN
            WRITE(INU,200)ITER,ERROR
            IF (ITYPE==1) THEN     !exp grid
               WRITE(*,211)N,ERROR 
            ELSEIF (ITYPE==3) THEN !sin grid
               WRITE(*,213)N,ERROR
            ELSEIF (ITYPE==4) THEN !abs grid
               WRITE(*,214)N,ERROR
            ELSE                   !cos grid
               WRITE(*,212)N,ERROR
            ENDIF
         ENDIF
      ENDIF
      IF(INU>=0)WRITE(INU,196)
      RETURN
   
100   FORMAT('ERROR in REMEZ_NL: no alternant found!',/,&
             'Hint: try to increase sample points for alternant NMAXALT.',/,&
             '      Used parameter for this run:',I10)
196   FORMAT('  ',76('-'))
197   FORMAT(' |',15(' '),'Calculation minimax points for imaginary axes',16(' '),'|')
198   FORMAT(' |','  Remez',I4,':    maximum error= ',E15.7,',    dE= ',E15.7'      |')
199   FORMAT('  WARNING: Remez not finished after ',I2,' steps, maximum error:',F14.8,/,&
             '  Calculated quadrature is not the minimax approximation! ')
200   FORMAT(' |  finished after ',I2,' steps, maximum error:',E20.10,15(' '),'|')

201   FORMAT(' |',26(' '),'(for exp transformation)',26(' '),'|')
202   FORMAT(' |',26(' '),'(for cos transformation)',26(' '),'|')
203   FORMAT(' |',26(' '),'(for sin transformation)',26(' '),'|')
204   FORMAT(' |',26(' '),'(for abs transformation)',26(' '),'|')

211   FORMAT(' MM quadrature (exp) found (',I2,'th order), error=',E20.10)
212   FORMAT(' MM quadrature (cos) found (',I2,'th order), error=',E20.10)
213   FORMAT(' MM quadrature (sin) found (',I2,'th order), error=',E20.10)
214   FORMAT(' MM quadrature (abs) found (',I2,'th order), error=',E20.10)

   ENDSUBROUTINE REMEZ_NL

!****************************************************************************
! Performs a linear least square fit of the form
!
! \sum g(i,k) Cos(\nu_k \tau_i) Exp(-x \tau_i) = x/(x^2+\nu_k^2)
! if  LCOS=.TRUE. or
!
! \sum g(i,k) Sin(\nu_k \tau_i) Exp(-x \tau_i) = \nu_k/(x^2+\nu_k^2)
! if LCOS=.FALSE.
!
! The arrays NU and TAU containing the points \a_k and \tau_i
! must be given in addition to the interval (X1,X2) in which X may vary
! On output GAM(NNU,NTAU) is the optimal forward Fourier matrix for
! the RPA polarizability, i.e.
! transformation of imaginary time to imaginary frequency domain
!****************************************************************************

   SUBROUTINE CALC_MINIMAX_FT( LCOS, NU, NNU, TAU, NTAU, RD, MAXITER, GAM,INU,ERRORS,FACT)
      USE prec
      IMPLICIT NONE
      LOGICAL          :: LCOS               !cos or sin transform?
      REAL(qd)         :: NU(NNU),TAU(NTAU)  !frequency and time points
      INTEGER          :: NNU, NTAU          !number of points, respectively
      REAL(q)          :: RD                 !minimization interval
      INTEGER          :: MAXITER            !maximal number of Remez iterations
      REAL(qd)         :: GAM(NNU,NTAU)      !Transformation matrix
      INTEGER          :: INU                !for output
      REAL(q)          :: ERRORS(NNU)        !resulting errors for each frequency
      REAL(q),OPTIONAL :: FACT               !factor to multiply error
!local
      INTEGER,  PARAMETER   :: NDATA=100      !number of data points
      INTEGER,  PARAMETER   :: NPRINT=5000
      INTEGER               :: O,T,J, I       !loop variables
      REAL                  :: X1,X2              !minimization interval
      REAL(qd)              :: X(NDATA)       !data points
      REAL(qd)              :: E,ET,FP_(2)    !error function and its derivative values
      REAL(qd), ALLOCATABLE :: A(:,:)         !design matrix of fitting problem
      REAL(qd), ALLOCATABLE :: U(:,:)         !left singluar vectors
      REAL(qd), ALLOCATABLE :: S(:)           !singular values
      REAL(qd), ALLOCATABLE :: VT(:,:)        !right singular vectors
      REAL(qd), ALLOCATABLE :: VAR(:)         !variances
      REAL(qd)              :: B(NDATA)       !function values
      REAL(qd), ALLOCATABLE :: G(:)           !temporary fourier matrix
      REAL(qd)              :: W              !weight
      REAL(qd)              :: R              !interval length
      REAL(qd)              :: MAD            !mean absolute deviation
      INTEGER               :: SIGNS(NDATA)   !signs of error for each data point
      REAL(qd)              :: XI(NDATA)      !alternant
      INTEGER               :: INFO,ITER      !loop variables
      INTEGER               :: NEXTREMA       !number of consideres extrema
      REAL(qd)              :: F(NDATA,NTAU+1)!matrix of minimax coefficients
      REAL(qd), ALLOCATABLE :: GTEMP(:)       !temporary coefficients
      REAL(qd)              :: MAXERROR1      !some errors
      REAL(qd)              :: MAXERROR2      !some errors
      REAL(qd)              :: VARTOT1,VARTOT2!auxilary
      REAL(q)               :: XX(NPRINT),ERROR
   
!dump preamble
      IF (INU>=0) THEN              
         WRITE(INU,100)
         IF (LCOS) THEN
            WRITE(INU,101) 
         ELSE
            WRITE(INU,102) 
         ENDIF 
      ENDIF

!use conventional Chebyshev nodes as Data points
      X1=1.
      X2=REAL(RD)
!minimization interval length
      R=REAL(RD,KIND=qd)
!use exponentially distributed Chebyshev nodes
!decreases the maximum error for x->R by a factor of ~ 0.75
      IF ( ABS(X1-1._qd) > ACC2 ) THEN
         WRITE(*,*)'ERROR in CALC_MINIMAX_FT, X1/=1',X1
         CALL M_exit(); stop
      ENDIF 
      DO I=1, NDATA
         X(I)=LOG(R)+LOG(R)*COS(((2*(2*NDATA-I)+1)*PIQ)/(4*NDATA))
         X(I)=EXP(X(I))
      ENDDO  
   
!allocate memory for design matrix A and fourier vector for frequency
      ALLOCATE( A(NTAU, NDATA) ) 
      ALLOCATE( G(NTAU+1),GTEMP(NTAU+1) ) 
      ALLOCATE( VAR( NTAU+1)) 
      A=0._qd
      G=0._qd
      GTEMP=0._qd
      VAR = 0._qd
   
!---------------------------------------------------
      frequency: DO O=1,NNU                   !nu loop
!---------------------------------------------------
         time: DO T=1,NTAU                    !tau loop
            IF (LCOS) THEN
               data_cos: DO J=1,NDATA         !data points
!set up design matrix A and vector B
                  A(T,J)=COS( NU(O)*TAU(T) )*EXP(-X(J)*TAU(T))
   
!B contains frequency dependence of polarizability
                  B(J)=X(J)/(X(J)*X(J)+NU(O)*NU(O))
               ENDDO data_cos
            ELSE 
               data_sin: DO J=1,NDATA         !data points
!set up design matrix A and vector B
                  A(T,J)=SIN( NU(O)*TAU(T) )*EXP(-X(J)*TAU(T))
!B contains frequency dependence of polarizability
                  B(J)=NU(O)/(X(J)*X(J)+NU(O)*NU(O))
               ENDDO data_sin
            ENDIF
         ENDDO time

!     WRITE(*,'(12F18.10)')B(1:NDATA)
      
!the optimal fourier matrix transforming from
!imaginary time points tau to imaginary frequency
!nu is the solution of a linear fiting problem
!The least square solution is simply the solution
!vector a of A^T . A . a = A^T .B
!we solve this using SVD of A=U.S.V^T
         ALLOCATE( U(NTAU, NTAU) ) ; U=0._qd
         ALLOCATE( S(NDATA) )      ; S=0._qd
         ALLOCATE(VT(NDATA,NDATA)) ; VT=0._qd
       
!compute SVD of A
         CALL CALC_SVD( A, NTAU, NDATA, U, S, VT)
!caluclate least square coefficients
!initialize solution vector
         G=0
         VAR=0
         VARTOT1=0
         DO T=1,NTAU
            DO J=1,NTAU  
!kill large inverse singular values
               IF (ABS(S(J)) < ACC2 ) THEN
                  W=0._qd
               ELSE
                  W=1._qd/S(J)
               ENDIF
!(taken from NR chapter 14 for
!linear least square fiting problems)
               G(T)=G(T)+W*DOT_PRODUCT(VT(J,1:NDATA),B(1:NDATA))*U(T,J)
!variances
               VAR(T)= VAR(T) + W*W*U(T,J)*U(T,J)
            ENDDO
            VARTOT1 = VARTOT1 + VAR(T)
         ENDDO
   
!save least square coefficients
         GTEMP(1:NTAU)=G(1:NTAU)
!find the maximum error for either the cos or sin transform
         IF (LCOS) THEN
            MAXERROR1=FINDMAX_LINEAR(NTAU,GTEMP,TAU_NU_COS_ERROR,X1,X2,NU(O),TAU)
         ELSE
            MAXERROR1=FINDMAX_LINEAR(NTAU,GTEMP,TAU_NU_SIN_ERROR,X1,X2,NU(O),TAU)
         ENDIF
!clean memory for S, V, U
         DEALLOCATE( VT, S, U )
   
!now G is least square solution vector of fitting problem
!next we need to minimize the error function of fitting problem
!using the Remez algorithm.
!This is almost trivial, since the problem is linear!
!However in general the alternant of the error function has more than
!NTAU+1 points, this means we acutally have to solve an overdetermined
!system of equations. This is 1._q best using the SVD of the coupling
!matrix
      
!save X to XI, this is advantageous if we use the MAD and not
!the sloppy REMEZ algorithm described above
         NEXTREMA=NDATA
         XI(1:NDATA) = X(1:NDATA) 
!==========================================================================
         remez:DO ITER=1,MAXITER
!==========================================================================
            IF (LCOS) THEN
               MAXERROR1=FINDMAX_LINEAR(NTAU,GTEMP,TAU_NU_COS_ERROR,X1,X2,NU(O),TAU)
               CALL FIND_MAD(TAU_NU_COS_ERROR,NTAU,TAU,NU(O),GTEMP,X1,X2,X,NDATA,MAD,SIGNS, INFO)  
            ELSE
               MAXERROR1=FINDMAX_LINEAR(NTAU,GTEMP,TAU_NU_SIN_ERROR,X1,X2,NU(O),TAU)
               CALL FIND_MAD(TAU_NU_SIN_ERROR,NTAU,TAU,NU(O),GTEMP,X1,X2,X,NDATA,MAD,SIGNS, INFO)  
            ENDIF
       
!build Coefficient matrix and r.h.s
            DO J=1,NEXTREMA
               DO T=1,NTAU+1
                  IF (T<=NTAU ) THEN
                     IF (LCOS) THEN
                        F(J,T)=COS(NU(O)*TAU(T))*EXP(-XI(J)*TAU(T))
                     ELSE
                        F(J,T)=SIN(NU(O)*TAU(T))*EXP(-XI(J)*TAU(T))
                     ENDIF
                  ELSE
!F(J,T)=(-1._qd)**J
                     F(J,T)=SIGNS(J)*MAD
                  ENDIF
               ENDDO   
!b is now vector of r.h.s at extrema Xi
               IF (LCOS) THEN
                  B(J)=XI(J)/(XI(J)*XI(J)+NU(O)*NU(O))   
               ELSE
                  B(J)=NU(O)/(XI(J)*XI(J)+NU(O)*NU(O))   
               ENDIF
!WRITE(*,'(11F10.6," | ",F10.6)')F(J,1:NTAU+1),B(J)
            ENDDO
!allocate memory for U,S,VT
            ALLOCATE( U(NEXTREMA, NEXTREMA) ) 
            ALLOCATE( S(NTAU+1) )           
            ALLOCATE(VT(NTAU+1,NTAU+1))     
            U=0._qd
            S=0._qd
            VT=0._qd
   
!solve overdetermined system using SVD of coupling matrix F
            CALL CALC_SVD( F(1:NEXTREMA,1:NTAU+1), NEXTREMA, NTAU+1, U, S, VT)
        
!calculate new coefficients
            G=0
            VAR=0  
            VARTOT2=0
            DO T=1,NTAU+1
               DO J=1,NTAU+1  
!kill large inverse singular values
                  IF (ABS(S(J)) < ACC2) THEN
                     W=0._qd
                  ELSE
                     W=1._qd/S(J)
                  ENDIF
                  G(T)=G(T)+W*DOT_PRODUCT(U(1:NEXTREMA,J),B(1:NEXTREMA))*VT(J,T)
!variances
                  VAR(T)= VAR(T) + W*W*VT(J,T)*VT(J,T)
               ENDDO
               VARTOT2=VARTOT2 + VAR(T)
            ENDDO
   
!kill U,S,VT
            DEALLOCATE(U,S,VT)
       
!check convergence
            IF (LCOS) THEN
               MAXERROR2=FINDMAX_LINEAR(NTAU,G,TAU_NU_COS_ERROR,X1,X2,NU(O),TAU)
            ELSE
               MAXERROR2=FINDMAX_LINEAR(NTAU,G,TAU_NU_SIN_ERROR,X1,X2,NU(O),TAU)
            ENDIF
!breaking condition
            IF ( ABS(MAXERROR2-MAXERROR1)>ACC ) THEN 
               GTEMP(1:NTAU+1)=G(1:NTAU+1)
            ELSEIF ( ABS(MAXERROR2-MAXERROR1) < ACC .AND. ITER <  5 ) THEN
               GTEMP(1:NTAU+1)=G(1:NTAU+1)
            ELSEIF ( ABS(MAXERROR2-MAXERROR1) < ACC .AND. ITER >= 5 ) THEN
               EXIT remez
            ENDIF       
         
!==========================================================================
          ENDDO remez
!==========================================================================
   
         IF (INU>=0) THEN
            IF ( MAXITER >= 0 ) THEN 
               IF ( PRESENT(FACT) ) THEN
                  WRITE(INU,103)O,FACT*NU(O),MAXERROR1/FACT,ITER !,ABS(MAXERROR2-MAXERROR1)
               ELSE
                  WRITE(INU,103)O,NU(O),MAXERROR1,ITER !,ABS(MAXERROR2-MAXERROR1)
               ENDIF
            ENDIF
         ENDIF
!save error of current frequency point O
         IF ( PRESENT(FACT) ) THEN
            ERRORS(O)=MAXERROR2/FACT
         ELSE
            ERRORS(O)=MAXERROR2
         ENDIF
   
!save fourier matrix
         DO T=1,NTAU
            GAM(O,T)=G(T)
         ENDDO
   
!print error function
# 824


!---------------------------------------------------
      ENDDO frequency
!---------------------------------------------------
 
      IF (INU>=0) THEN
         WRITE(INU,110)
# 836

      ENDIF
      DEALLOCATE(A,G,GTEMP,VAR)
      RETURN
    
100   FORMAT('  ',76('-'))
101   FORMAT(' |',18(' '),'Calculation of cos transformation matrix',18(' '),'|')
102   FORMAT(' |',18(' '),'Calculation of sin transformation matrix',18(' '),'|')
103   FORMAT(' |    nu_',I2,'=',E15.7,' ERR= ',E15.7,' finished after ',I3,' steps     |')
110   FORMAT('  ',76('-'))
111   FORMAT(32E12.4)
   
   ENDSUBROUTINE CALC_MINIMAX_FT

!****************************************************************************
! Gauss-Chebyshev quadrature for the interval [X1,X2]
! quadruple version
!****************************************************************************

   SUBROUTINE GAUSS_CHEBYSHEV1Q(X1,X2,X,W,N)
      USE prec
      IMPLICIT NONE
      INTEGER N
      REAL(qd) X1,X2
      REAL(qd) X(N),W(N)
!local variables
      REAL(qd) XX(N),WW(N)
      INTEGER I,J,M
       
      M=N 
      DO I=1,M
!for [-1,1] interval
         XX(I)=COS(REAL( 2*I-1 ,KIND=qd)*PIQ/REAL(2*M,KIND=qd))
         WW(I)=PIQ/REAL(M,KIND=qd)
!scaled onto [X1,X2]
         X(I)=(X1+X2)/2._qd + (X2-X1)/2._qd*XX(I)
         W(I)=WW(I)*X1
      ENDDO
   
      RETURN
   ENDSUBROUTINE GAUSS_CHEBYSHEV1Q

!****************************RTSAFEL*****************************************
! simple root finder for linear error function (used in minimax_ft)
!****************************************************************************

   SUBROUTINE RTSAFEL(ERRORFUNC,N,TAU,NU,C,X1,X2,XACC,ROOT,INFO,GUESS)
      USE prec
      IMPLICIT NONE
      EXTERNAL ERRORFUNC
      INTEGER N,INFO
      REAL(qd)  X1,X2
      REAL(qd)  XACC
      REAL(qd)  C(N)
      REAL(qd)  NU, TAU(N)
      REAL(qd)  ROOT
      REAL(qd)  GUESS
      INTEGER  I
      REAL(qd)  F,FP(2),DX,Z
      
      Z=GUESS
      
      DO I=1,MAXNEWTON
         CALL ERRORFUNC(N,C,Z,F,FP,NU,TAU)
         DX=FP(1)/FP(2)
! WRITE(*,'(4F40.25)')GUESS,Z,DX,FP(1)
         Z=Z-DX
         IF (ABS(DX)<XACC) THEN
            INFO=0
            ROOT = Z
            RETURN
         ENDIF
      ENDDO
       
      INFO=1
      RETURN
   ENDSUBROUTINE RTSAFEL
   
!****************************************************************************
! sorts coefficients C small size N (<100) in ascending order (pick sort)
!****************************************************************************

   SUBROUTINE PIKSRT(N,C)
      USE prec
      REAL(qd) ARR(N) 
      REAL(qd) A,C(2*N)
    
      IF (N>100) THEN 
        WRITE(*,*)' Warning ing PIKSRT: This will take long since N=',N
      ENDIF
     
!abscissas
      ARR(1:N)=C(1:N)
      DO J=2, N
         A=ARR(J)
         DO I=J-1,1,-1
            IF (ARR(I)<=A) GOTO 10
            ARR(I+1)=ARR(I)
         ENDDO
         I=0
10       ARR(I+1)=A
      ENDDO
!store back to original array
      C(1:N)=ARR(1:N)
   
!weights
      ARR(1:N)=C(N+1:2*N)
      DO J=2, N
         A=ARR(J)
         DO I=J-1,1,-1
            IF (ARR(I)<=A) GOTO 11
            ARR(I+1)=ARR(I)
         ENDDO
         I=0
11       ARR(I+1)=A
      ENDDO
!store back to original array
      C(N+1:2*N)=ARR(1:N)
   
      RETURN
   ENDSUBROUTINE PIKSRT

!****************************************************************************
! This routine calculates the Singular value decomposition A=U.S.V^T
! of a given matrix  A(M,N)
!****************************************************************************

   SUBROUTINE CALC_SVD( AA, M, N, UU, SS, VVT )
      USE prec
      IMPLICIT NONE
      REAL(qd)  :: AA(M,N)  !matrix to be decomposed
      INTEGER   :: M,N      !number of rows and columns
      REAL(qd)  :: UU(M,M)  !left singular eigenvectors
      REAL(qd)  :: SS(N)    !signular values
      REAL(qd)  :: VVT(N,N) !right singular eigenvectors
!local
      INTEGER, PARAMETER   :: LWMAX=20000
      INTEGER              :: LDA, LDU, LDVT
      INTEGER              :: INFO, LWORK 
      INTEGER              :: I,J
      INTEGER              :: IWORK( 8*N )
      REAL(q), ALLOCATABLE ::  WORK(:)  
!work with double precision
      REAL(q) A(M,N), U(M,M),S(N),VT(N,N)  

!first define leading dimensions
      LDA = M
      LDU = M
      LDVT = N

      ALLOCATE( WORK( LWMAX ) ); WORK = 0._qd
 
!save quadruple precision into double precision variables
      DO I=1,M
         DO J=1,N
            A(I,J)=REAL(AA(I,J),q)
         ENDDO
      ENDDO
      
      DO J=1,M
         DO I=1,M
            U(J,I)=REAL(UU(J,I),q)
         ENDDO
      ENDDO
       
      DO J=1,N
         DO I=1,N
            VT(J,I)=REAL(VVT(J,I),q)
         ENDDO
         S(J)=REAL(SS(J),q)
      ENDDO

!Query the optimal workspace
      LWORK = -1
      CALL DGESDD( 'Singular vectors', M, N, A, LDA, S, U, LDU, VT,&
                  LDVT, WORK, LWORK, IWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      
!Compute SVD
      CALL DGESDD( 'Singular vectors', M, N, A, LDA, S, U, LDU, VT,&
                  LDVT, WORK, LWORK, IWORK, INFO )

!check for convergence
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'CAL_SVD Error: DGESDD failed to converge.',INFO
         CALL M_exit(); stop
      ENDIF

      DEALLOCATE(WORK ) 

!save quadruple precision into double precision variables
      DO J=1,N
         DO I=1,M
            AA(I,J)=REAL(A(I,J),KIND=qd)
         ENDDO
         SS(J)=REAL(S(J),KIND=qd)
      ENDDO
      
      DO J=1,M
         DO I=1,M
            UU(J,I)=REAL(U(J,I),KIND=qd)
         ENDDO
!IF(N==11)   WRITE(*,'(19F10.6)')U(J,1:M)
      ENDDO
       
      DO J=1,N
         DO I=1,N
            VVT(J,I)=REAL(VT(J,I),KIND=qd)
         ENDDO
      ENDDO

      RETURN
   ENDSUBROUTINE CALC_SVD

!***********************SOLVE_NONLS_FREQ*************************************
!    This routine solve the non linear system
!
!         y(x_i,L_i)+(-1)^i-1 L_{2N+1} = INTERPOL(x_i)
!
! on entry:  X(1:2*N+1) .... extrema of nonlinear error function
!            L(1:2*N)   .... coefficients of fit
!            L(2*N+1)   .... Maximum of errorfunction in [A,B]
!
! on exit:   L(1:2*N+1) .... new coefficients
!****************************************************************************

   SUBROUTINE SOLVE_NONLS(NLFITERGRAD,N,L,X,ACCURACY,INF)
      USE prec
      IMPLICIT NONE
      EXTERNAL               :: NLFITERGRAD !error function and its gradient
      INTEGER                :: N           !order of error function
      REAL(qd),INTENT(INOUT) :: L(2*N+1)    !coefficients of error function
      REAL(qd)               :: X(2*N+1)    !Alternant of error function
      REAL(qd)               :: ACCURACY    !accuracy
      INTEGER                :: INF         !information about convergence
!local
      INTEGER              :: NL            !number of variables
      INTEGER              :: ITER          !iteration loop variable
      INTEGER              :: I,K 
      INTEGER              :: INFO          !for QDLUDCMP
      INTEGER              :: LWORK         !for QDLUDCMP
      INTEGER,ALLOCATABLE  :: IPIV(:)       !for QDLUDCMP
      REAL(qd)             :: WORK(1000)    !for QDLUDCMP
      REAL(qd)             :: F             !error value at X
      REAL(qd)             :: DF(1:2*N)     !derivative of error at X
      REAL(qd)             :: VARIANCE      !variance of current solution to previous 1._q
      REAL(qd)             :: D             
      REAL(qd),ALLOCATABLE :: J(:,:)        !Jacobian
      REAL(qd),ALLOCATABLE :: LAMB(:)       !current solution
      REAL(qd),ALLOCATABLE :: Y(:)          !r.h.s. of system of equations
      REAL(qd),ALLOCATABLE :: DY(:)         !new direction
    
!allocate some stuff with size of total number of cofficients
      NL=2*N+1   
      ALLOCATE ( IPIV( NL ), J( NL, NL ), LAMB( NL ), Y( NL ), DY( NL ) )
    
      DO I = 1, NL
         LAMB(I)=L(I)
      ENDDO
   
!---------------------------------------------------------------------------
      conv: DO ITER=1,MAXNEWTON
!---------------------------------------------------------------------------
! we solve the nonlinear system iteratively using the Newton-Raphson
! method. For this we need the Jacobian w.r.t. the parameter L
         DO I=1,NL
!error function and its gradient
            CALL NLFITERGRAD(N,LAMB,X(I),F,DF)
!IF(ITER==1)WRITE(*,'(I4,F11.6,"|",F20.10,"|",10F20.10)')I,X(I),F,DF(1:MIN(NL-1,10))
!r.h.s.
            Y(I)=-(F-LAMB(NL)*(-1)**(I))   
!l.h.s.
            DO K=1,NL-1
            J(I,K)=DF(K)
            ENDDO
            J(I,NL)=-(-1)**(I)
         ENDDO
! LU decomposition
         CALL QDLUDCMP(J,NL,NL,IPIV,D,INFO)
         IF (INFO /= 0 ) THEN
            WRITE(*,*)'ERROR in SOLVE_NONLS_FREQ: QDLUDCMP failed with:',INFO
            CALL M_exit(); stop
         ENDIF 
      
!numerically more stable to use DGETRS than
         CALL QDLUBKSB( J, NL, NL, IPIV, Y)
      
         IF (ITER>0) THEN    !damp here
            DO I=1, NL
!weights must be positive!
               IF ( (LAMB(I)+DAMPNEWT*Y(I)<0._qd) .AND. (I>N .AND. I/=NL) ) THEN 
                  WRITE(*,*)ITER,'Error in SOLVE_NONLS, weight is negative:',I
                  WRITE(*,120)NMAXALT
                  CALL M_exit(); stop
               ELSE
                  LAMB(I)=LAMB(I)+DAMPNEWT*Y(I)
               ENDIF
            ENDDO  
         ELSE                
            DO I=1, NL
               IF ( (LAMB(I)+Y(I)<0._qd) .AND. (I>N .AND. I/=NL) ) THEN  !weights must be positive!
                  WRITE(*,*)ITER,'Error in SOLVE_NONLS_FREQ, weight is negative:',I
                  WRITE(*,120)NMAXALT
                  CALL M_exit(); stop
               ELSE
                  LAMB(I)=LAMB(I)+Y(I)
               ENDIF
            ENDDO  
         ENDIF
         
         VARIANCE=SQRT(DOT_PRODUCT(Y,Y))
          
         IF (VARIANCE<ACCURACY) THEN
            EXIT conv       !solution found
         ENDIF    
      
!---------------------------------------------------------------------------
      ENDDO  conv 
!---------------------------------------------------------------------------

      IF ( ITER > MAXNEWTON ) THEN
         WRITE(*,100)MAXNEWTON
         CALL M_exit(); stop
      ELSE 
!solution found
!store back to original array
         L(1:NL)=LAMB(1:NL)
         INF=0                   
      ENDIF 

      DEALLOCATE( IPIV, J, LAMB, Y, DY ) 
      
120   FORMAT('             Hint: try to increase sample points for alternant NMAXALT.',/,&
             '                   Used parameter for this run:',I10)
100   FORMAT('Error in SOLVE_NONLS: Newton-Raphson not converged after',I4,&
             ' iterations')   
   
   ENDSUBROUTINE SOLVE_NONLS

!****************************************************************************
! calculates error function:
!
!  \sum_i \gamma_i COS(\nu \tau_i)Exp(-x\tau_i) -x/(x^2-\nu^2)
!
! this is need for the optimal fourier transformation from \tau into
! \nu domain
!
!****************************************************************************

   SUBROUTINE TAU_NU_COS_ERROR( N, L, XX, Z, ZP,NU,TAU )
      USE prec
      IMPLICIT NONE
      INTEGER   :: N      ! order of error function
      REAL(qd)  :: L(N)   ! parameters
      REAL(qd)  :: XX     ! argument X of error
      REAL(qd)  :: Z      ! function value at X
      REAL(qd)  :: ZP(2)  ! derivatives of Y w.r.t. X
      REAL(qd)  :: NU     ! imaginary frequency \nu
      REAL(qd)  :: TAU(N) ! imaginary time points
!local
      INTEGER I,J,M
      REAL(qd) ZTMP,ZTMP_
      REAL(qd),ALLOCATABLE :: ZTMPP(:)
      REAL(qd)  X

      IF ( SIZE(L) < N) THEN
         WRITE(*,*)'Error in TAU_NU_COS_ERROR: Size of L incompatible:',SIZE(L),N
         CALL M_exit(); stop
      ENDIF
      
      X=REAL(XX,KIND=qd)
     
!initialize
      ZTMP=XX/(XX*XX+NU*NU)
    
      M=SIZE( ZP ) 
      ALLOCATE( ZTMPP( M ) ) 
! derivative of 1/X
      ZTMPP(1)=(-XX*XX+NU*NU)/((XX*XX+NU*NU)**2)
      DO J=2,M
      ZTMPP(J)=2._qd*(XX*XX*XX-3*XX*NU*NU)/((XX*XX+NU*NU)**3)
      ENDDO
     
!error function
      DO I=1,N 
        ZTMP = ZTMP- L(I)*COS(NU*TAU(I))*EXP(-X*TAU(I) )
      ENDDO
     
! Jth derivative
      DO J=1,M
         DO I=1,N
         ZTMPP(J) = ZTMPP(J) - L(I)*COS(NU*TAU(I))*EXP( -TAU(I)*X )*( -TAU(I) )**J 
         ENDDO
      ENDDO
     
      Z=-ZTMP
      ZP(1:M)=-ZTMPP(1:M)
      DEALLOCATE(ZTMPP)

      RETURN
   ENDSUBROUTINE TAU_NU_COS_ERROR

!****************************************************************************
! calculates error function:
!
!  \sum_i \gamma_i SIN(\nu \tau_i)Exp(-x\tau_i) -\nu/(x^2-\nu^2)
!
! this is need for the optimal fourier transformation from \tau into
! \nu domain
!
!****************************************************************************

   SUBROUTINE TAU_NU_SIN_ERROR( N, L, XX, Z, ZP,NU,TAU )
      USE prec
      IMPLICIT NONE
      INTEGER         :: N          ! order of error function
      REAL(qd)        :: L(N)       ! parameter of error function
      REAL(qd)        :: XX         ! argument X of error function
      REAL(qd)        :: Z          ! error at X
      REAL(qd)        :: ZP(2)      ! derivatives of error w.r.t. X
      REAL(qd)        :: NU         ! imaginary frequency \nu
      REAL(qd)        :: TAU(N)     ! imaginary time points
!local
      INTEGER I,J,M
      REAL(qd) ZTMP,ZTMP_
      REAL(qd),ALLOCATABLE :: ZTMPP(:)
      REAL(qd)  X

      IF ( SIZE(L) < N) THEN
         WRITE(*,*)'Error in TAU_NU_SIN_ERROR: Size of L incompatible:',SIZE(L),N
         CALL M_exit(); stop
      ENDIF
      
      X=REAL(XX,KIND=qd)
   
!initialize
      ZTMP=NU/(XX*XX+NU*NU)
      M=SIZE( ZP ) 
      ALLOCATE( ZTMPP( M ) ) 

! derivative of 1/X
      ZTMPP(1)=(-2._qd*XX*NU)/((XX*XX+NU*NU)**2)
      DO J=2,M
         ZTMPP(J)=-2._qd*NU*(-3._qd*XX*XX+NU*NU)/((XX*XX+NU*NU)**3)
      ENDDO
   
!error function
      DO I=1,N 
         ZTMP = ZTMP- L(I)*SIN(NU*TAU(I))*EXP(-X*TAU(I) )
      ENDDO
   
!Jth derivative
      DO J=1,M
         DO I=1,N
            ZTMPP(J) = ZTMPP(J) - L(I)*SIN(NU*TAU(I))*EXP( -TAU(I)*X )*( -TAU(I) )**J 
         ENDDO
      ENDDO
   
      Z=-ZTMP
      ZP(1:M)=-ZTMPP(1:M)
      DEALLOCATE(ZTMPP)
      RETURN
   ENDSUBROUTINE TAU_NU_SIN_ERROR

!****************************************************************************
!  Calculates the nonlinear error function at X
!
!   E(X,\lambda)=INTERPOL(X)-\sum_I=1^N \lambda(I+N)Exp(-\lambda(I)X)
!
! entry: N...Order
!        L...parameter array \lambda
!
! exit : Z...............Y(X,\lambda)
!        ZP(1:SIZE(ZP))..Derivatives up to SIZE(ZP) of Y(X,\lambda) w.r.t. X
!
! if INFO
!****************************************************************************

  SUBROUTINE TAU_ERROR_FUNCTION( N, L, X, Z, ZP )
     USE prec
     IMPLICIT NONE
     INTEGER            :: N        ! order of error
     REAL(qd)           :: L(2*N)   ! parameter
     REAL(qd)           :: X       ! argument X of error
     REAL(qd)           :: Z        ! eror at X
     REAL(qd)           :: ZP(2)    ! derivatives of Y w.r.t. X
!local
     INTEGER I,J,M
     REAL(qd) ZTMP,ZTMP_
     REAL(qd),ALLOCATABLE :: ZTMPP(:)
  
     IF ( SIZE(L) < 2*N) THEN
        WRITE(*,*)'Error in TAU_ERROR_FUNCTION: Size of L incompatible:',SIZE(L),2*N
        CALL M_exit(); stop
     ENDIF
     
!initialize
     ZTMP=INTERPOL(X)
     M=SIZE( ZP ) 
     ALLOCATE( ZTMPP( M ) ) 
! derivative of 1/X
     ZTMPP(1)=-INTERPOL(X)/X

     DO J=2,M
         ZTMPP(J)=(-J)*ZTMPP(J-1)/X
     ENDDO
  
     DO I=1,N 
        ZTMP = ZTMP - L(I+N)*EXP(-L(I)*X )
     ENDDO
  
! Jth derivative
     DO J=1,M
        DO I=1,N
           ZTMPP(J) = ZTMPP(J) - L(I+N)*EXP( -L(I)*X )*( -L(I) )**J
        ENDDO
     ENDDO
     
     Z=ZTMP
     ZP(1:M)=ZTMPP(1:M)
     DEALLOCATE(ZTMPP)
  
     RETURN
  ENDSUBROUTINE TAU_ERROR_FUNCTION  

!****************************************************************************
!  Calculates the nonlinear error function at X
!
!   E(X,\lambda)=INTERPOL(X) - 1/Pi \sum_I=1^N \lambda(I+N) 4 X^2/( X^2+\lambda(I)^2)^2
!
! entry: N...Order
!        L...parameter array \lambda
!
! exit : Z...............E(X,\lambda)
!        ZP(1:SIZE(ZP))..Derivatives up to SIZE(ZP) of E(X,\lambda)
!                        w.r.t. X
!****************************************************************************

   SUBROUTINE NUCOS_ERROR_FUNCTION(N,L,X,Z,ZP)
      USE PREC
      IMPLICIT NONE
      INTEGER        :: N       ! order of error
      REAL(qd)       :: L(2*N)  ! parameter array
      REAL(qd)       :: X       ! argument X of error
      REAL(qd)       :: Z       ! error at X
      REAL(qd)       :: ZP(2)   ! derivatives of error w.r.t. X
!local
      INTEGER I,J,M
      REAL(qd) ZTMP,ZTMP_
      REAL(qd),ALLOCATABLE :: ZTMPP(:)
      
      IF ( SIZE(L) < 2*N) THEN
         WRITE(*,*)'Error in NUCOS_ERROR_FUNCTION: Size of L incompatible:',SIZE(L),2*N
         CALL M_exit(); stop
      ENDIF

      ZTMP=INTERPOL(X)
      M=SIZE( ZP ) 
      ALLOCATE( ZTMPP( M ) ) 
! derivative of 1/X
      ZTMPP(1)=-INTERPOL(X)/X
      DO J=2,M
         ZTMPP(J)=(-J)*ZTMPP(J-1)/X
      ENDDO
     
      DO I=1,N
!error
         ZTMP = ZTMP - L(I+N)*( REAL(4,KIND=qd)*(X**2) )/( ( X*X+L(I)*L(I) )**2 )/PIQ
! 1st derivative
         ZTMPP(1)= ZTMPP(1)+((REAL(16,KIND=qd)*(X**3)*L(I+N))/((X*X+L(I)*L(I))**3) )/PIQ-&
                            (( REAL(8,KIND=qd)*X*L(I+N)/((X*X+L(I)*L(I))**2)) )/PIQ
! 2nd derivative
         ZTMPP(2)= ZTMPP(2)-(( REAL(96,KIND=qd)*(X**4)*L(I+N) )/( ( X*X + L(I)*L(I) )**4)/PIQ ) + &
                            (( REAL(80,KIND=qd)*(X**2)*L(I+N) )/( ( X*X + L(I)*L(I) )**3)/PIQ ) -&
                            (( REAL(8,KIND=qd)*L(I+N) )/( ( X*X + L(I)*L(I) )**2 )/PIQ )
      ENDDO
      
      Z=ZTMP
      ZP(1:M)=ZTMPP(1:M)
      DEALLOCATE(ZTMPP)
      RETURN
   ENDSUBROUTINE NUCOS_ERROR_FUNCTION

!****************************************************************************
!  Calculates the nonlinear error function at X
!
!   E(X,\lambda)=INTERPOL(X) - 1/Pi \sum_I=1^N \lambda(I+N) 4 X^2/( X^2+\lambda(I)^2)^2
!
! entry: N...Order
!        L...parameter array \lambda
!
! exit : Z...............E(X,\lambda)
!        ZP(1:SIZE(ZP))..Derivatives up to SIZE(ZP) of E(X,\lambda)
!                        w.r.t. LAMBDA
! (used for solving non-linear fiting problem)
!****************************************************************************

   SUBROUTINE NUCOS_ERROR_FUNCTION_GRAD(N,L,X,Z,ZP)
      USE prec
      IMPLICIT NONE
      INTEGER                :: N        ! order of fit
      REAL(qd)               :: L(2*N+1) ! parameter array
      REAL(qd)               :: X        ! argument X of Y
      REAL(qd)               :: Z        ! fit at X
      REAL(qd),INTENT(INOUT) :: ZP(2*N)  ! derivatives of Y w.r.t. L
!local
      INTEGER I,J,M
      REAL(qd) ZTMP
      REAL(qd),ALLOCATABLE :: ZTMPP(:)
      LOGICAL :: LDER
    
      IF ( SIZE(L) < 2*N) THEN
         WRITE(*,*)'Error in NUCOS_ERROR_FUNCTION_GRAD: Size of L incompatible:',SIZE(L),2*N
         CALL M_exit(); stop
      ENDIF
    
      IF (SIZE(ZP) < 2*N) THEN
         WRITE(*,*)'Error in NUCOS_ERROR_FUNCTION_GRAD, Size of ZP incompatible:',SIZE(ZP),2*N
         CALL M_exit(); stop
      ENDIF
     
!initialize
      ZTMP=INTERPOL(X)
      M=2*N 
      ALLOCATE( ZTMPP( M ) ) 
    
      DO I=1,N
         ZTMP      = ZTMP - L(I+N)*( 4._qd*(X**2) )/( ( X*X+L(I)*L(I) )**2 )/PIQ
         ZTMPP(I)  = +L(I)*L(I+N)*( 16._qd*(X**2) )/( ( X*X+L(I)*L(I) )**3 )/PIQ
         ZTMPP(I+N)= -( 4._qd*(X**2) )/( ( X*X+L(I)*L(I) )**2 )/PIQ
      ENDDO
    
      Z=ZTMP
      ZP(1:M)=ZTMPP(1:M)
     
      DEALLOCATE(ZTMPP)
      RETURN
   ENDSUBROUTINE NUCOS_ERROR_FUNCTION_GRAD
   
!****************************************************************************
!  Calculates the nonlinear error function at X
!
!   E(X,\lambda)=INTERPOL(X) - 1/Pi \sum_I=1^N \lambda(I+N) 4 \lambda(I)/( X^2+\lambda(I)^2)^2
!
! entry: N...Order
!        L...parameter array \lambda
!
! exit : Z...............E(X,\lambda)
!        ZP(1:SIZE(ZP))..Derivatives up to SIZE(ZP) of E(X,\lambda)
!                        w.r.t. X
!****************************************************************************

   SUBROUTINE NUSIN_ERROR_FUNCTION(N,L,X,Z,ZP)
      USE PREC
      IMPLICIT NONE
      INTEGER        :: N       ! order of error
      REAL(qd)       :: L(2*N)  ! parameter array
      REAL(qd)       :: X       ! argument X of error
      REAL(qd)       :: Z       ! error at X
      REAL(qd)       :: ZP(2)   ! derivatives of error w.r.t. X
!local
      INTEGER I,J,M
      REAL(qd) ZTMP,ZTMP_
      REAL(qd),ALLOCATABLE :: ZTMPP(:)
      LOGICAL :: LDER
     
      IF ( SIZE(L) < 2*N) THEN
         WRITE(*,*)'Error in NUSIN_ERROR_FUNCTION: Size of L incompatible:',SIZE(L),2*N
         CALL M_exit(); stop
      ENDIF
                                    
                 
      ZTMP=INTERPOL(X)
      M=SIZE( ZP ) 
      ALLOCATE( ZTMPP( M ) ) 
! derivative of 1/X
      ZTMPP(1)=-INTERPOL(X)/X
      DO J=2,M
         ZTMPP(J)=(-J)*ZTMPP(J-1)/X
      ENDDO
     
      DO I=1,N
         ZTMP = ZTMP - L(I+N)*( REAL(4,KIND=qd)*L(I)*L(I) )/( ( X*X+L(I)*L(I) )**2 )/PIQ
! 1st derivative
         ZTMPP(1) = ZTMPP(1)+((L(I)*L(I)*REAL(16,KIND=qd)*X*L(I+N))/((X*X+L(I)*L(I))**3))/PIQ
! 2nd derivative
         ZTMPP(2) = ZTMPP(2)-((L(I)*L(I)*REAL(96,KIND=qd)*X*X*L(I+N))/((X*X+L(I)*L(I))**4))/PIQ  &
                            +((L(I)*L(I)*REAL(16,KIND=qd)*L(I+N))/((X*X+L(I)*L(I))**3))/PIQ 
      ENDDO
     
      Z=ZTMP
      ZP(1:M)=ZTMPP(1:M)
      DEALLOCATE(ZTMPP)
      RETURN
   ENDSUBROUTINE NUSIN_ERROR_FUNCTION

!****************************************************************************
!  Calculates the nonlinear error function at X
!
!   E(X,\lambda)=INTERPOL(X) - 2/Pi \sum_I=1^N \lambda(I+N)  1/( X^2+\lambda(I)^2)^2
!
! entry: N...Order
!        L...parameter array \lambda
!
! exit : Z...............E(X,\lambda)
!        ZP(1:SIZE(ZP))..Derivatives up to SIZE(ZP) of E(X,\lambda)
!                        w.r.t. X
!****************************************************************************

   SUBROUTINE NUABS_ERROR_FUNCTION(N,L,X,Z,ZP)
      USE PREC
      IMPLICIT NONE
      INTEGER        :: N       ! order of error
      REAL(qd)       :: L(2*N)  ! parameter array
      REAL(qd)       :: X       ! argument X of error
      REAL(qd)       :: Z       ! error at X
      REAL(qd)       :: ZP(2)   ! derivatives of error w.r.t. X
!local
      INTEGER I,J,M
      REAL(qd) ZTMP,ZTMP_
      REAL(qd),ALLOCATABLE :: ZTMPP(:)
   
      IF ( SIZE(L) < 2*N) THEN
         WRITE(*,*)'Error in NUABS_ERROR_FUNCTION: Size of L incompatible:',SIZE(L),2*N
         CALL M_exit(); stop
      ENDIF
                 
      ZTMP=INTERPOL(X)
      M=SIZE( ZP ) 
      ALLOCATE( ZTMPP( M ) ) 
! derivative of 1/X
      ZTMPP(1)=-INTERPOL(X)/X
      DO J=2,M
        ZTMPP(J)=(-J)*ZTMPP(J-1)/X
      ENDDO
   
      DO I=1,N
         ZTMP = ZTMP - L(I+N)*(2._qd/PIQ)*1._qd/( X*X+L(I)*L(I))
! 1st derivative
         ZTMPP(1) = ZTMPP(1) + (2._qd/PIQ)*( 2._qd*X*L(I+N))/((X*X+L(I)*L(I))**2)
! 2nd derivative
         ZTMPP(2) = ZTMPP(2) - (2._qd/PIQ)*L(I+N)*(6._qd*X**2-2._qd*L(I)**2)/((X*X+L(I)*L(I))**3)
      ENDDO

      Z=ZTMP
      ZP(1:M)=ZTMPP(1:M)
      DEALLOCATE(ZTMPP)
      RETURN
   ENDSUBROUTINE NUABS_ERROR_FUNCTION

!****************************************************************************
!  Calculates the nonlinear error function at X
!
!   E(X,\lambda)=INTERPOL(X) - 2/Pi \sum_I=1^N \lambda(I+N) 1 \lambda(I)/( X^2+\lambda(I)^2)
!
! entry: N...Order
!        L...parameter array \lambda
!
! exit : Z...............E(X,\lambda)
!        ZP(1:SIZE(ZP))..Derivatives up to SIZE(ZP) of E(X,\lambda)
!                        w.r.t. LAMBDA
! (used for solving non-linear fiting problem)
!****************************************************************************

   SUBROUTINE NUABS_ERROR_FUNCTION_GRAD(N,L,X,Z,ZP)
      USE prec
      IMPLICIT NONE
      INTEGER                :: N       ! order of error function
      REAL(qd)               :: L(2*N+1)! parameter array
      REAL(qd)               :: X       ! argument X of Y
      REAL(qd)               :: Z       ! fit at X
      REAL(qd),INTENT(INOUT) :: ZP(2*N) ! derivatives of Y w.r.t. L
!local
      INTEGER                :: I,J,M
      REAL(qd)               :: ZTMP
      REAL(qd),ALLOCATABLE   :: ZTMPP(:)
   
      IF ( SIZE(L) < 2*N) THEN
         WRITE(*,*)'Error in NUABS_ERROR_FUNCTION_GRAD: Size of L incompatible:',SIZE(L),2*N
         CALL M_exit(); stop
      ENDIF
   
      IF (SIZE(ZP) < 2*N) THEN
         WRITE(*,*)'Error in NUABS_ERROR_FUNCTION_GRAD, Size of ZP incompatible:',SIZE(ZP),2*N
         CALL M_exit(); stop
      ENDIF
     
!initialize
      ZTMP=INTERPOL(X)
      M=2*N 
      ALLOCATE( ZTMPP( M ) ) 
   
      DO I=1,N
!function value
         ZTMP=ZTMP-(2._qd/PIQ)*L(I+N)/(X**2+L(I)**2)
!nonlinear coeficients derivative
         ZTMPP(I)=(2._qd/PIQ)*L(I+N)*(2._qd*L(I))/(( X*X+L(I)*L(I) )**2)
!linear coefficients derivative
         ZTMPP(I+N)=-(2._qd/PIQ)/(X**2+L(I)**2)
      ENDDO
   
      Z=ZTMP
      ZP(1:M)=ZTMPP(1:M)
     
      DEALLOCATE(ZTMPP)
      RETURN
   ENDSUBROUTINE NUABS_ERROR_FUNCTION_GRAD

!****************************************************************************
!  Calculates the nonlinear error function at X
!
!   E(X,\lambda)=INTERPOL(X) - 1/Pi \sum_I=1^N \lambda(I+N) 4 \lambda(I)/( X^2+\lambda(I)^2)^2
!
! entry: N...Order
!        L...parameter array \lambda
!
! exit : Z...............E(X,\lambda)
!        ZP(1:SIZE(ZP))..Derivatives up to SIZE(ZP) of E(X,\lambda)
!                        w.r.t. LAMBDA
! (used for solving non-linear fiting problem)
!****************************************************************************

   SUBROUTINE NUSIN_ERROR_FUNCTION_GRAD(N,L,X,Z,ZP)
      IMPLICIT NONE
      INTEGER                :: N       ! order of error function
      REAL(qd)               :: L(2*N+1)! parameter array
      REAL(qd)               :: X       ! argument X of Y
      REAL(qd)               :: Z       ! fit at X
      REAL(qd),INTENT(INOUT) :: ZP(2*N) ! derivatives of Y w.r.t. L
!local
      INTEGER                :: I,J,M
      REAL(qd)               :: ZTMP
      REAL(qd),ALLOCATABLE   :: ZTMPP(:)

      IF ( SIZE(L) < 2*N) THEN
         WRITE(*,*)'Error in NUSIN_ERROR_FUNCTION_GRAD: Size of L incompatible:',SIZE(L),2*N
         CALL M_exit(); stop
      ENDIF

      IF (SIZE(ZP) < 2*N) THEN
         WRITE(*,*)'Error in NUSIN_ERROR_FUNCTION_GRAD, Size of ZP incompatible:',SIZE(ZP),2*N
         CALL M_exit(); stop
      ENDIF
 
!initialize
      ZTMP=INTERPOL(X)
      M=2*N 
      ALLOCATE( ZTMPP( M ) ) 

      DO I=1,N
         ZTMP      = ZTMP - L(I+N)*( 4._qd*(L(I)**2) )/( ( X*X+L(I)*L(I) )**2 )/PIQ
!derivative w.r.t. abscissa
         ZTMPP(I)  =      + L(I+N)*(16._qd*(L(I)**3) )/( ( X*X+L(I)*L(I) )**3 )/PIQ &
                          - L(I+N)*( 8._qd*(L(I)   ) )/( ( X*X+L(I)*L(I) )**2 )/PIQ
!dervative w.r.t. weight
         ZTMPP(I+N)= -( 4._qd*(L(I)**2) )/( ( X*X+L(I)*L(I) )**2 )/PIQ
      ENDDO

      Z=ZTMP
      ZP(1:M)=ZTMPP(1:M)
      DEALLOCATE(ZTMPP)
      RETURN
   ENDSUBROUTINE NUSIN_ERROR_FUNCTION_GRAD

!****************************************************************************
!  Calculates the nonlinear error function at X
!
!   E(X,\lambda)=INTERPOL(X) - \sum_I=1^N \lambda(I+N)Exp(-\lambda(I)X)
!
! entry: N...Order
!        L...parameter array \lambda
!
! exit : Z...............E(X,\lambda)
!        ZP(1:SIZE(ZP))..Derivatives up to SIZE(ZP) of E(X,\lambda)
!                        w.r.t. LAMB
!****************************************************************************

   SUBROUTINE TAU_ERROR_FUNCTION_GRAD( N, L, X, Z, ZP )
      USE prec
      IMPLICIT NONE
      INTEGER                :: N       ! order of fit
      REAL(qd)               :: L(2*N)  ! parameter array
      REAL(qd)               :: X       ! argument X of Y
      REAL(qd)               :: Z       ! fit at X
      REAL(qd),INTENT(INOUT) :: ZP(2*N) ! derivatives of Y w.r.t. L
!local
      INTEGER J,M
      REAL(qd) ZTMP
      REAL(qd),ALLOCATABLE :: ZTMPP(:)
     
      IF ( SIZE(L) < 2*N) THEN
         WRITE(*,*)'Error in TAU_ERROR_FUNCTION_GRAD: Size of L incompatible:',SIZE(L),2*N
         CALL M_exit(); stop
      ENDIF
      
      IF (SIZE(ZP) < 2*N) THEN
         WRITE(*,*)'Error in TAU_ERROR_FUNCTION_GRAD, Size of ZP incompatible:',SIZE(ZP),2*N
         CALL M_exit(); stop
      ENDIF
   
!initialize
      ZTMP=INTERPOL(X)
      M=2*N
      ALLOCATE( ZTMPP( M ) )
      ZTMPP=0._qd
     
      DO J=1,N
!error function
         ZTMP = ZTMP - L(J+N)*EXP(-L(J)*X )
! derivative w.r.t. to exponent
         ZTMPP(J) = X*L(J+N)*EXP( -L(J)*X )
! derivative w.r.t. to weight
         ZTMPP(J+N) = -EXP( -L(J)*X )
      ENDDO
      
      Z=ZTMP
      ZP(1:M)=ZTMPP(1:M)
      DEALLOCATE(ZTMPP)
      
      RETURN
   ENDSUBROUTINE TAU_ERROR_FUNCTION_GRAD  
   
!****************************************************************************
! finds the mean absolute deviation
!****************************************************************************

   SUBROUTINE FIND_MAD(ERRORFUNC,N,TAU,NU,C,XX1,XX2,X,NDATA,MAD,SIGNS,INFO)
      USE prec
      USE m_unirnk
      IMPLICIT NONE
      EXTERNAL :: ERRORFUNC
      INTEGER  :: N                           ! fit order
      REAL(qd) :: TAU(N)                      ! imaginary time points
      REAL(qd) :: NU                          ! imaninary frequency point
      REAL(qd) :: C(N)
      REAL     :: XX1,XX2                     ! coefficients and interval
      REAL(qd) :: X(NDATA)                    ! data points
      INTEGER  :: SIGNS(NDATA)                ! tabulated signs of error for each data point
      INTEGER  :: NDATA                       ! number of data points
      REAL(qd) :: MAD                         ! mean absolute deviatation
      INTEGER  :: INFO
!local
      INTEGER :: I,J
      REAL(q) :: FI(NDATA)
      REAL    :: X00(NDATA)
      INTEGER :: LX(NDATA),K
      REAL(qd):: FA,DF1(2),X1,X2
      INTEGER :: INF
      
!first step is to tabulate signs of error function for
!each data point X(I) and its value
      DO I = 1 , NDATA
         CALL ERRORFUNC(N,C,X(I),FA,DF1,NU,TAU)
         IF ( FA > 0 ) SIGNS(I)=1
         IF ( FA < 0 ) SIGNS(I)=-1
         FI(I) = ABS( REAL( FA , q ) )   
      ENDDO 
       
!now sort the array FI
      CALL UNIRNK( FI, LX, K ) 
       
!caluclate Median of FI,
!for odd  K
      IF ( MOD ( K , 2 ) /= 0  ) THEN
!middle value is at K/2+1
         MAD = FI ( LX (  K/2 + 1  ) )  
!for even K
      ELSE
         MAD = ( FI ( LX ( K/2 ) ) +  FI ( LX (  K/2 + 1  ) ) )/2 
      ENDIF
   
!finally the MAD
      DO I =1 , NDATA
        FI ( I ) = ABS ( FI(I) - MAD )
      ENDDO

!now sort the array FI
      CALL UNIRNK( FI, LX, K ) 
    
!caluclate Median deviation of FI,
!for odd  K
      IF ( MOD ( K , 2 ) /= 0  ) THEN
!middle value is at K/2+1
         MAD = FI ( LX (  K/2 + 1  ) )  
!for even K
      ELSE
         MAD = ( FI ( LX ( K/2 ) ) +  FI ( LX (  K/2 + 1  ) ) )/2 
      ENDIF

      RETURN                                     !exit if all extrema are found
   ENDSUBROUTINE FIND_MAD

!****************************************************************************
! find alternant using extrema
!****************************************************************************

   SUBROUTINE FIND_ALTERNANT_EXTREMA(NLFIT,N,XX1,XX2,C,X,NDAT,INFO)
      USE prec
      USE m_unirnk
      IMPLICIT NONE
      EXTERNAL :: NLFIT
      INTEGER  :: N                      ! fit order
      REAL(qd) :: C(2*N)
      REAL(qd) :: XX1,XX2                ! coefficients and interval
      REAL(qd) :: X(2*N+1)               ! extrema
      REAL(qd) :: ACC                    ! accuracy
      INTEGER :: NDAT
      INTEGER :: INFO
!local
      INTEGER  :: I,J,K
      REAL(qd),ALLOCATABLE :: XI(:),WI(:),X0(:)
      REAL,ALLOCATABLE     :: X00(:)
      INTEGER,ALLOCATABLE  :: LX(:)
      REAL(qd) :: FA,DF2(2),DF1(2),FP_(2),X1,X2
      INTEGER  :: INF,ISTART,IFINAL,NDATA
!needed for rounding
      INTEGER(kind=8)    :: II
      REAL(qd),PARAMETER :: R=1.E+6_qd !8 significant digits
   
      INFO = 0
      
      X1=REAL(XX1,KIND=qd)
      X2=REAL(XX2,KIND=qd)  
     
!first of all find abscissas being in [x1,x2]
      DO I=1, N
        IF ( X1-C(I) < 0 ) EXIT
      ENDDO 
      ISTART=I
      DO I=1,N
        IF (X2-C(I) < 0 ) EXIT
      ENDDO
      
      IF ( I > N ) THEN
         IFINAL=N
      ELSE
         IFINAL=I
      ENDIF
   
      NDATA=(IFINAL-ISTART)*NDAT+2*NDAT
      ALLOCATE(XI(NDATA),WI(NDATA),X0(NDATA))
      ALLOCATE(X00(NDATA))
      ALLOCATE(LX(NDATA))
     
!setup sample nodes
      K=0
      DO I=1,NDAT
         K=K+1
         XI(K)=X1+(I-1)*(C(ISTART)-X1)/(NDAT-1._qd)
      ENDDO
      DO I=ISTART,IFINAL-1
         DO J=1,NDAT
            K=K+1
            XI(K)=C(I)+(J-1)*(C(I+1)-C(I))/(NDAT-1._qd)
         ENDDO
      ENDDO
      DO I=1,NDAT
         K=K+1
         XI(K)=C(IFINAL)+(I-1)*(X2-C(IFINAL))/(NDAT-1._qd)
      ENDDO
   
! first create a lot of chebyshev nodes
!      CALL GAUSS_CHEBYSHEV1Q(X1,X2,XI,WI,NDATA)
   
!----------------------------------------------
      node: DO I = 1, NDATA
!----------------------------------------------
         CALL RTSAFENL2(NLFIT,N,C,X1,X2,1.E-16_qd,X0(I),INF,XI(I))
         IF ( INF == -1 .OR. INF == 1) THEN   
            X0(I)=-100._qd
         ENDIF
   
!root must not be negative
         IF ( X0(I) < 0._qd ) THEN
            X0(I)=-100._qd
         ENDIF
         
!root must not be out of interval
         IF (X0(I)>0 .AND. (X0(I)<X1 .OR. X0(I)>X2) )THEN
            X0(I)=-100._qd
         ENDIF 
       
!backup for comparision
         X00(I)=REAL(X0(I))
!----------------------------------------------
      ENDDO node
!----------------------------------------------
 
!remove double entries
      CALL UNIRNK(X00,LX,K)
   
!write warning  if not enough points found
      IF ( K < 2*N-1  ) THEN 
         INFO = -1
         RETURN
      ENDIF
   
!most probably the first entry of X0 is -100
!however if this is not the case we need to use
!this point as first alternant point after X1
      IF ( .NOT. ( X00( LX(1) ) < 0._qd ) ) THEN
         X(2) = X0( LX(1) )
!save the remaining points
         DO I=2,2*N
            X( I + 1 ) = X0( LX(I) )
         ENDDO
      ELSE
         DO I=2, 2*N
            X( I ) = X0( LX(I) ) 
         ENDDO
      ENDIF
   
!CALL PRINT_ERROR_FUNCTION(NUCOS_ERROR_FUNCTION,N,C,X1,X2,1000,.TRUE.)
!add the end points
      K=K+2
      X(1)=X1  
      X(2*N+1)=X2
   
!A last check is necessary
!for high order N~20, still two or more points can coincide
      DO I=1,2*N+1
         DO J=1, 2*N+1
            IF ( I==J ) CYCLE
            IF ( ABS(X(I)-X(J)) < 1._qd/R ) THEN
               INFO = -I
               IF ( NDATA >1000 )WRITE(*,'(I5,2I3,3F45.32)')NDATA,I,J,X(I),X(J)
            ENDIF  
         ENDDO
      ENDDO

      IF ( INFO <= -1 ) RETURN
   
      INFO = 0
      RETURN                                     !exit if all extrema are found
   ENDSUBROUTINE FIND_ALTERNANT_EXTREMA

!****************************RTSAFENL2***************************************
! simple root finder for nonlinear error function
!****************************************************************************

   SUBROUTINE RTSAFENL2(NLFIT,N,C,X1,X2,XACC,ROOT,INFO,GUESS)
      USE prec
      EXTERNAL  NLFIT  !error function
      INTEGER   N,INFO
      REAL(qd)  X1,X2
      REAL(qd)  XACC   !accuracy
      REAL(qd)  C(2*N)
      REAL(qd)  ROOT
      REAL(qd)  GUESS
      INTEGER  J
      REAL(qd)  F,FP(2),DX,Z
   
      Z=GUESS
   
      DO I=1,MAXNEWT
         CALL NLFIT(N,C,Z,F,FP)
         DX=FP(1)/FP(2)
         Z=Z-DX
         IF (ABS(DX)<XACC) THEN
            INFO=0
            ROOT = Z
            RETURN
         ENDIF
      ENDDO
      
      INFO=1
      RETURN
   ENDSUBROUTINE RTSAFENL2

!********************************************************
! Calculated the maximum of L1 norm in the interval A,B
! for a linear error function
!********************************************************

   FUNCTION FINDMAX_LINEAR(N,L,FIT,X1,X2,NU,TAU)
      USE prec
      IMPLICIT NONE
      INTEGER N
      REAL(qd) L(1:N)
      REAL(qd) FINDMAX_LINEAR
      EXTERNAL FIT
      REAL X1,X2
      REAL(qd) NU
      REAL(qd) TAU(N)
!localz
      INTEGER I
      REAL(qd) E,ET,FP_(2)
      REAL(qd) X,XT,A,B
     
      IF (SIZE(L)<N) THEN 
         WRITE(*,*)' Error in FINDMAX_LINEAR: Sizes incompatible in L',SIZE(L),N
         CALL M_exit(); stop
      ENDIF
   
      A=REAL(X1,KIND=qd)
      B=REAL(X2,KIND=qd)
   
      CALL FIT(N,L,A,ET,FP_,NU,TAU)
      ET=ABS(ET)
      DO I=1,MAXBISEC
         X=(A+B)/2._qd+(B-A)/2._qd*COS((2*(2*MAXBISEC-I)+1)*PIQ/(4*MAXBISEC))
         CALL FIT(N,L,X,E,FP_,NU,TAU)
         ET=MAX(ABS(ET),ABS(E))
      ENDDO
   
      FINDMAX_LINEAR=ET
      RETURN
   ENDFUNCTION FINDMAX_LINEAR

!********************************************************
! prints the nonlinear error function and its first two
! derivatives w.r.t. argument X
!********************************************************

  SUBROUTINE PRINT_ERROR_FUNCTION(NLFIT,N,L,A,B,IMAX,LLOG,IT)
     USE prec
     IMPLICIT NONE
     EXTERNAL NLFIT
     INTEGER N,IMAX
     REAL(qd) L(2*N)
     REAL(qd) A,B
     LOGICAL,INTENT(IN),OPTIONAL :: LLOG
     INTEGER,INTENT(IN),OPTIONAL :: IT 
     INTEGER I,IFIN
     REAL(qd) X
     REAL(qd) F,DF(2)
     LOGICAL LLN
     CHARACTER(LEN=4) :: APP 
  
     IFIN=IMAX
    
     IF (PRESENT(LLOG)) THEN
        LLN=LLOG
     ELSE
        LLN=.FALSE.
     ENDIF
     
     IF (PRESENT(IT)) THEN
        WRITE(APP,'(I4)')IT
        APP=ADJUSTL(APP)
        OPEN(UNIT=90,FILE='erf'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
     ELSE
        OPEN(UNIT=90,FILE='erf.dat',ACTION='WRITE',STATUS='REPLACE')
     ENDIF
  
!choose abscissas logarithmic
     IF (LLN) THEN
        DO I=1,IFIN
!we use a logarithmic partition
           X=LOGPARTITION(I,IFIN,A,B)
           CALL NLFIT(N,L,X,F,DF)
           WRITE(90,'(4F40.30)')X,F,DF(1),DF(2)
        ENDDO
     ELSE
        DO I=1,IFIN
!equidistant abscissas
           X=A+(I-1)*(B-A)/REAL( IFIN,KIND=qd)
           CALL NLFIT(N,L,X,F,DF)
           WRITE(90,'(4F40.30)')X,F,DF(1),DF(2)
        ENDDO
     ENDIF
    
     CLOSE(90)
     RETURN
  ENDSUBROUTINE PRINT_ERROR_FUNCTION

!********************************************************
! Calculates the maximum of the L1 norm in the interval
! A,B of the exp error function
!********************************************************

   FUNCTION FINDMAX(N,L,NLFIT,AA,BB)
      USE prec
      INTEGER N
      REAL(qd) L(1:2*N)
      REAL(qd) FINDMAX
      EXTERNAL NLFIT
      REAL(qd) :: AA,BB
!localz
      INTEGER I
      REAL(qd) E,ET,FP_(2)
      REAL(qd) X,XT,A,B
      REAL(qd),ALLOCATABLE :: XI(:)
      INTEGER  :: ISTART
     
      IF (SIZE(L)<2*N) THEN 
         WRITE(*,*)' Error in FINDMAX: Sizes incompatible in L',SIZE(L),2*N
         CALL M_exit(); stop
      ENDIF

!smart sampling of minimization interval:
!first of all find first abscissa in [x1,x2]
      DO I=1, N
        IF ( AA-L(I) < 0 ) EXIT
      ENDDO 
      ISTART=I
      
      IF (ISTART+1 <= N)THEN       
         A=AA
         B=MIN(L(ISTART),BB)
      ELSE
         A=AA
         B=2*L(ISTART)
      ENDIF

!obtain the error at A
      CALL NLFIT(N,L,A,ET,FP_)
   
      ET=ABS(ET)
      DO I=1,MAXBISEC
         X=A+(I-1)*(B-A)/REAL(MAXBISEC-1,KIND=qd) 
         CALL NLFIT(N,L,X,E,FP_)
         E=ABS(E)
         IF ( ET < E ) THEN
            ET=E
            XT=X
         ENDIF
      ENDDO
   
      FINDMAX=ET
      RETURN
   ENDFUNCTION FINDMAX

!*********************************************************
!
! LU decomposition for quadruple precision  (NR)
!
! Given a matrix a(1:n,1:n) , with physical dimension
! np by np , this routine replaces it by the
! LU decomposition of a rowwise permutation of itself.
! a and n are input. a is output, arranged as in
! equation (2.3.14) of NR above; indx(1:n) is an output vector
! that records the row permutation effected by the partial
! pivoting; d is output as ±1 depending on whether the number
! of row interchanges was even or odd, respectively.
! This routine is used in combination with qdlubksb to solve
! linear equations
!
!*********************************************************

   SUBROUTINE QDLUDCMP(A,N,NP,INDX,D, INFO)
      IMPLICIT NONE
      REAL(qd) :: A(NP,NP)
      INTEGER  :: N
      INTEGER  :: NP
      INTEGER  :: INDX(N)
      REAL(qd) :: D
      INTEGER  :: INFO 
!local
      INTEGER, PARAMETER  :: NMAX=500
      INTEGER I,IMAX,J,K
      REAL(qd) AAMAX,DUM,SUM,VV(NMAX) 
 
!vv stores the implicit scaling of each row.
      D=1._qd
        
!no row interchanges yet.
!loop over rows to get the implicit scaling information.
      DO I=1,N
         AAMAX=0._qd
         DO J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
         ENDDO
         IF ( AAMAX .EQ. 0.) THEN
!singular matrix in qdludcmp no nonzero largest element.
            INFO=-I 
            RETURN
         ENDIF  
!save the scaling.
         VV(I)=1._qd/AAMAX
      ENDDO

!this is the loop over columns of crout’s method.
      DO J=1,N
!this is equation (2.3.12) except for i = j.
         DO I=1,J-1
            SUM=A(I,J)
            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            ENDDO
            A(I,J)=SUM
         ENDDO
!initialize for the search for largest pivot element.
         AAMAX=0._qd

!this is i = j of equation (2.3.12) and i = j + 1 . . . n
!of equation (2.3.13).
         DO I=J,N
            SUM=A(I,J)
            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            ENDDO
            A(I,J)=SUM
!figure of merit for the pivot.
            DUM=VV(I)*ABS(SUM)
!is it better than the best so far?
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            ENDIF
         ENDDO

!do we need to interchange rows?
         IF (J.NE.IMAX)THEN
!YES, DO SO...
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            ENDDO
!...and change the parity of d.
            D=-D
!also interchange the scale factor.
            VV(IMAX)=VV(J)
         ENDIF

         INDX(J)=IMAX
!if the pivot element is 0._q the matrix is singular (at least to the precision of the al-
!gorithm). for some applications on singular matrices, it is desirable to substitute tiny
!for 0._q.
         IF( A(J,J) .EQ. 0. )A(J,J)=VERYTINY

!now, finally, divide by the pivot element.
         IF(J.NE.N)THEN
            DUM=1._qd/A(J,J)
            DO I=J+1,N
                A(I,J)=A(I,J)*DUM
            ENDDO
         ENDIF
      ENDDO
 
      INFO=0
      RETURN
      ENDSUBROUTINE QDLUDCMP
 
!*********************************************************
! quadruple precision of DGETRS
!
! Solves the set of n linear equations A · X = B. Here a is
! input, not as the matrix A but rather as its LU
! decomposition, determined by the routine ludcmp. indx is
! input as the permutation vector returned by ludcmp.
! b(1:n) is input as the right-hand side vector B, and
! returns with the solution vector X. a , n , np, and
! indx are not modified by this routine and can be left
! in place for successive calls with different right-hand
! sides b. This routine takes into account the possibility
! that b will begin with many 0._q elements, so it is
! efficient for use in matrix inversion.
!*********************************************************

   SUBROUTINE QDLUBKSB(A,N,NP,INDX,B)
      INTEGER  N
      INTEGER  NP
      INTEGER  INDX(N)
      REAL(qd) A(NP,NP)
      REAL(qd) B(N)
!local
      INTEGER I,II,J,LL
      REAL(qd) SUM

      II=0
!when ii is set to a positive value, it will become the in-
!dex of the first nonvanishing element of b. we now do
!the forward substitution, equation (2.3.6). the only new
!wrinkle is to unscramble the permutation as we go.
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II.NE.0)THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            ENDDO
!a nonzero element was encountered, so from now on we will
!have to do the sums in the loop above.
         ELSEIF (SUM.NE.0._qd) THEN
            II=I
         ENDIF
         B(I)=SUM
      ENDDO

!now we do the backsubstitution, equation (2.3.7).
      DO I=N,1,-1
         SUM=B(I)
         DO J=I+1,N
            SUM=SUM-A(I,J)*B(J)
         ENDDO 
!store a component of the solution vector x.
         B(I)=SUM/A(I,I)
      ENDDO
!all 1._q!
      RETURN
   ENDSUBROUTINE QDLUBKSB
        
!============================================================================
! Least Square Routines

!****************************************************************************
! calling sequence for varpro.
! performs a nonlinear least square fit of order NORDER starting with the
! first order
!****************************************************************************

   SUBROUTINE NONLINEAR_LSQ_FIT(ITYPE,R,NORDER,C,IO)
      USE base
      USE separable_leastsq
      IMPLICIT NONE
      INTEGER         :: ITYPE      !base function selection
      REAL (qd)       :: R          !interval
      INTEGER         :: NORDER     !order of fit
      REAL(qd)        :: C(2*NORDER)!coefficients
      TYPE(in_struct) :: IO         !in_struct
!local
      INTEGER, PARAMETER :: NMAX=300
      REAL(qd),PARAMETER :: EPS=1.E-14_qd
      REAL(qd)  :: WEIGHTSCALE     !scaling factor for weights
      REAL (qd),ALLOCATABLE :: Y(:), T(:,:),A(:,:), W(:)
      REAL (qd),ALLOCATABLE :: ALF(:), BETA(:)
      REAL (qd),ALLOCATABLE :: RAN(:)
      REAL (qd),ALLOCATABLE :: ALFINIT(:), BETAINIT(:)
      INTEGER   :: I, IERR, IM1, IPRINT, IV, J, L, LPNL, N, NL, P
      INTEGER   :: IGUESS,ERROR
      LOGICAL   :: LRETRY, LEXIST
      CHARACTER(LEN=1)           :: CHK
      CHARACTER(LEN=100)         :: DUM
      INTEGER       :: NSTART, ITYP,NC
      INTEGER, SAVE :: IDEL=-1
      INTEGER       :: INU

!standard output is OUTCAR
      INU=IO%IU6 
    
      IF ( IO%NWRITE ==1 ) THEN
!surpress output of varpro, not recomended!!
         IPRINT = -1 
      ELSEIF (IO%NWRITE == 2 .OR. IO%NWRITE == 3  ) THEN
!normal output of varpro
         IPRINT = 0 
      ELSE
!print complete information
         IPRINT = 1
      ENDIF 
 
! determine weighting factor for different error functions
      IF ( ITYPE == 1 ) THEN
! time error function
         WEIGHTSCALE=1._qd
      ELSEIF( ITYPE==2 .OR. ITYPE==3 ) THEN   
!cos error function
         WEIGHTSCALE=PIQ/4._qd
      ELSEIF( ITYPE == 4 ) THEN
!abs error function
         WEIGHTSCALE=PIQ/2._qd
      ELSE
         IF ( INU>=0) WRITE(*,*)'ERROR in NONLINEAR_LSQ_FIT, type of error function is unkown',ITYPE
         CALL M_exit(); stop
      ENDIF 

      ALLOCATE(Y(NMAX),T(NMAX,1),A(NMAX,2*AMAX+2),W(NMAX))
      ALLOCATE(ALF(AMAX), BETA(BMAX))
      ALLOCATE(ALFINIT(AMAX), BETAINIT(BMAX))
 
!initialize saved coefficients
      ALFINIT=0._qd 
      BETAINIT=0._qd 
      ALF=0._qd
      BETA=0._qd

      IV=1  !# of linear independent variables
!starting order
      NSTART=1

!write to OUTCAR only if IO%NWRITE is higher as default
      IF (INU>=0 .AND. IO%NWRITE > 2 ) THEN
         WRITE (INU,110) NDATALSQ, L, NL, P, IV
         IF ( ITYPE == 2 ) THEN
            WRITE(INU,997)NORDER
         ELSEIF(ITYPE==3) THEN
            WRITE(INU,998)NORDER
         ELSEIF(ITYPE==4) THEN
            WRITE(INU,999)NORDER
         ELSE
            WRITE(INU,996)NORDER
         ENDIF
      ENDIF
!exponentially dirstributed Chebysev nodes yield smaller errors for x->1
      DO I =1, NDATALSQ
         T(I,1)=LOG(R)+LOG(R)*COS(((2*(2*NDATALSQ-I)+1)*PIQ)/(4._qd*NDATALSQ))
         T(I,1)=EXP(T(I,1))
         W(I) = 1.0_qd !WEIGHTS ARE SET TO 1
         IF ( T(I,1) > ACC2 ) THEN
            Y(I)=1._qd/(T(I,1))
         ELSE
            IF ( INU>=0) WRITE(INU,*)'ERROR IN NONLINEAR_LSQ_FIT: POINT',I,' TOO SMALL',T(I,1)
            CALL M_exit(); stop
         ENDIF 
!print sample points only for IO%NWRITE >2
         IF ( INU>=0 .AND.  IO%NWRITE > 2 ) WRITE (INU,210) I, T(I,1:IV), Y(I)
      ENDDO
 
      LRETRY=.FALSE.
      NC=0
      DO N =NSTART, NORDER
         NL=N  !# of current non-linear parameters
         L=NL  !# of linear parameters
         P=NL  !# of non-vanishing non-linear partial derivatives
!use as starting gues coefficients of previous order
         IF ( N>1)THEN
            ALF=0
            BETA=0
            ALF(1:N-1)=ALFINIT(1:N-1)
            BETA(1:N-1)=BETAINIT(1:N-1)
         ENDIF
! we need initial guesses for the parameters
         guess: DO IGUESS =1, IGUESS_MAX
            ALF(N)=0._qd+(IGUESS-1)*RANDLSQ/REAL(IGUESS_MAX-1,qd)
            BETA(N)=0._qd+(IGUESS-1)*RANDLSQ/REAL(IGUESS_MAX-1,qd)

!print initial guess of non-linear coefficients
            IF ( N<NORDER .AND. INU>=0 .AND. IO%NWRITE>2 ) WRITE (INU,130) ALF(1:N)
            IF ( INU>=0 .AND. IO%NWRITE>2) WRITE (INU,140)
            CALL VARPRO(N,N,NDATALSQ,NMAX,L+P+2,IV,T,Y,W,ADA,A,IPRINT,ALF,BETA,ITMAXLSQ,&
               ITYPE,IERR,INU,IO%NWRITE)
!varpro converged
            IF (IERR>0) THEN
               IF (LRETRY) THEN
                  IF ( IO%NWRITE>2 .AND. INU>=0 ) WRITE(INU,101)IGUESS
                  LRETRY=.FALSE.
               ENDIF 
               IF ( IO%NWRITE>2 .AND. INU>=0 ) WRITE(INU,102)NL,IERR
               IF ( IO%NWRITE>2 .AND. INU>=0 ) WRITE(*,102)NL,IERR
!converged -> save coefficients for next higher order
               ALFINIT=0
               BETAINIT=0
               ALFINIT(1:N)=ALF(1:N)
               BETAINIT(1:N)=BETA(1:N)
               EXIT guess
!varpro not converged
            ELSE
               IF (INU>=0 .AND. IO%NWRITE>2 ) WRITE (INU,'("                IERR =",2I7)') &
                  IERR,IGUESS
               IF ( .NOT. LRETRY ) THEN
                  IF (INU>=0 .AND. IO%NWRITE>2 ) WRITE(INU,103)N
               ENDIF 
               LRETRY=.TRUE.
            ENDIF  
         ENDDO guess

!truly converged, if less iterations than IGUESS_MAX needed
         IF ( IGUESS <= IGUESS_MAX) THEN
!optionally, dump coefficients to OUTCAR
            IF ( (N==NORDER .AND. IO%NWRITE>2) .AND. INU>=0) THEN
               IF (ITYPE==1)THEN
                  WRITE(INU,115)N
               ELSEIF(ITYPE==2)THEN
                  WRITE(INU,116)N
               ELSEIF(ITYPE==3)THEN
                  WRITE(INU,117)N
               ELSEIF(ITYPE==4)THEN
                  WRITE(INU,118)N
               ENDIF
               WRITE(INU,106)(BETA(I),ALF(I),I=1,N)
            ENDIF
!save coefficients
            DO I=1,N
               C(I)=ABS(ALF(I))
               C(I+N)=BETA(I)*WEIGHTSCALE
            ENDDO
         ELSE
            IF (INU>=0) WRITE(*,1001) 
            CALL M_exit(); stop
         ENDIF
      ENDDO

   IF ( INU>=0 .AND. IO%NWRITE>1 ) WRITE(*,1000)NORDER

110 FORMAT ('1   NON-LINEAR LEAST SQUARES PROBLEM'/  &
    ' NUMBER OF OBSERVATIONS =',i5, '   NUMBER OF LINEAR PARAMETERS =',i4/  &
    ' NUMBER OF NONLINEAR PARAMETERS =',i4,  &
    ' NUMBER OF NONVANISHING PARTIAL DERIVATIVES =',i4/  &
    ' NUMBER OF INDEPENDENT VARIABLES =',i4/ '    I    T(I)            Y(I)')
130 FORMAT ('0 INITIAL NONLINEAR PARAMETERS'/(4E20.7))
140 FORMAT ('0',50('*'))
210 FORMAT (i5,7E16.7)

101 FORMAT(/,'      success after ',I4,' attempts.')
102 FORMAT('      Order',I3' done after ',I6,' iterations')
103 FORMAT('          using random numbers for order',I3,' .')
115 FORMAT(' LSQ quadrature (exp) found (',I2,'th order):')
116 FORMAT(' LSQ quadrature (cos) found (',I2,'th order):')
117 FORMAT(' LSQ quadrature (sin) found (',I2,'th order):')
118 FORMAT(' LSQ quadrature (abs) found (',I2,'th order):')
106 FORMAT(2F45.32)

996  FORMAT(/,' ',68('*'),/,'  Looking for exp least square quadrature of order ',I3)
997  FORMAT(/,' ',68('*'),/,'  Looking for cos least square quadrature of order ',I3)
998  FORMAT(/,' ',68('*'),/,'  Looking for sin least square quadrature of order ',I3)
999  FORMAT(/,' ',68('*'),/,'  Looking for abs least square quadrature of order ',I3)
1000 FORMAT(' Found non-linear least square fit of order ',I2,&
            ', set NWRITE>2 for more details')
1001 FORMAT(/,'       maximum number of iterations reached, calculation terminated!',/,&
            '       Increasing the minimization interval might help!')
ENDSUBROUTINE NONLINEAR_LSQ_FIT

!****************************************************************************
!
! evaluate the functions phi and their partial derivatives
! d phi(j)/d alf(k), at the sample points t(i).
! isel = 1  means compute both functions and derivatives and
!           initialize inc and constant phi's
!      = 2  means compute only the nonconstant functions phi
!      = 3  means compute only the derivatives.
!
! this particular routine is for fitting two exponentials and
! 1._q constant term:
!
! eta(c, alf; t) = c  + c  * exp(-alf *t) + c * exp(-alf * t)
!                   1    2           1       3          2
!
! where c is the vector of linear parameters to be determined,
! and alf is the vector of nonlinear parameters to be determined.
! here n = 33, l = 3, nl = p = 2, ncfun = 1, phi(1) = 1, phi(2) =
! exp(-alf(1)*t), and phi(3) = exp(-alf(2)*t).  this example and
! test data are taken from reference 3 (see varpro).
!
!****************************************************************************

   SUBROUTINE ADA (LP1, NL, N, NMAX, LPP2, IV, A, INC, T, ALF, ISEL, ITYPE)
      USE SEPARABLE_LEASTSQ
      IMPLICIT NONE
      INTEGER, INTENT(IN)     :: LP1
      INTEGER, INTENT(IN)     :: NL
      INTEGER, INTENT(IN)     :: N
      INTEGER, INTENT(IN)     :: NMAX
      INTEGER, INTENT(IN)     :: LPP2
      INTEGER, INTENT(IN)     :: IV
      REAL (DP), INTENT(OUT)  :: A(:,:)        ! A(NMAX,LPP2)
      INTEGER, INTENT(OUT)    :: INC(:,:)      ! INC(12,LP1)
      REAL (DP), INTENT(IN)   :: T(:,:)        ! T(NMAX)
      REAL (DP), INTENT(IN)   :: ALF(:)        ! ALF(NL)
      INTEGER, INTENT(IN)     :: ISEL
      INTEGER, INTENT(IN)     :: ITYPE
!local
      INTEGER :: I,J,K
      
      SELECT CASE ( ISEL )
         CASE (1)
            GO TO 10
         CASE (2)
            GO TO 30
         CASE (3)
            GO TO 50
      END SELECT
      
!              ISEL = 1, in this case the incidence matrix inc is:
!                        diagonal !!
      
!                 1  0       <- ALF(1)
!                 0  1       <- ALF(2)
      
10    INC = 0
      DO I=1, NL
         INC(I,I) = 1
      ENDDO
!              ISEL = 1, SET CONSTANT FUNCTION PHI(1) (REMOVE IF NO
!              CONSTANT FUNCTIONS)
      
! A(1:N,1) = 1.0_DP
      
!              ISEL = 1 OR 2, COMPUTE NONCONSTANT FUNCTIONS PHI
30    CONTINUE
        
      IF ( ITYPE == 2 ) THEN
!cos related base function
         DO  I=1,N
            DO J=1,NL
                A(I,J) = (T(I,1)/(T(I,1)**2+ALF(J)**2))**2
            ENDDO
         END DO
      ELSEIF ( ITYPE == 3 ) THEN
!sin related base function
         DO  I=1,N
            DO J=1,NL
                A(I,J) = (ALF(J)/(T(I,1)**2+ALF(J)**2))**2
            ENDDO
         END DO
      ELSEIF ( ITYPE == 4 ) THEN
!cos**2 + sin**2 = abs related base function
         DO  I=1,N
            DO J=1,NL
                A(I,J) = 1._DP/(T(I,1)**2+ALF(J)**2)
            ENDDO
         END DO
      ELSE 
!exp base functin is default
         DO  I=1,N
            DO J=1,NL
                A(I,J) = EXP(-ALF(J)*T(I,1))
            ENDDO
         END DO
      ENDIF 
!             COLUMN L+1 = 3 IS LEFT FOR WORKSPACE.
      IF (ISEL == 2) GO TO 70
      
!              ISEL = 1 OR 3, COMPUTE DERIVATIVES DPHI(I) / D ALF(J)
50    CONTINUE
      IF ( ITYPE == 2 ) THEN   !cosine grid
!DERIVATIVE W.R.T. T
         DO  I=1,N
            DO J=NL+2,2*NL+1
               K=J-NL-1
               A(I,J) = -4._DP*T(I,1)**2*ALF(K)/((T(I,1)**2+ALF(K)**2)**3)
            ENDDO
         END DO
      ELSEIF ( ITYPE == 3) THEN  !sine grid
         DO  I=1,N
            DO J=NL+2,2*NL+1
               K=J-NL-1
               A(I,J) = 2._DP*ALF(K)*(T(I,1)**2-(ALF(K)**2))/((T(I,1)**2+ALF(K)**2)**3)
            ENDDO
         END DO
      ELSEIF ( ITYPE == 4 ) THEN !abs grid
         DO  I=1,N
            DO J=NL+2,2*NL+1
               K=J-NL-1
               A(I,J) = -2._DP*ALF(K)/((T(I,1)**2+ALF(K)**2)**2)
            ENDDO
         END DO
      ELSE                       !exp time grid
         DO  I=1,N
            DO J=NL+2,2*NL+1
               K=J-NL-1
               A(I,J) = -T(I,1)*EXP(-ALF(K)*T(I,1))
            ENDDO
         END DO
      ENDIF 
      
70    RETURN
   END SUBROUTINE ADA

!****************************************************************************
!
! Find LSQ grids by converging from RMAX to R
!
!****************************************************************************

   SUBROUTINE CONVERGE_FROM_RMAX2R_LSQ(ITYPE,RMAX,R,N,COF,IO)
      USE base
      USE separable_leastsq
      INTEGER   :: ITYPE           !base function selection
      REAL (q)  :: RMAX            !maximum interval
      REAL (q)  :: R               !desired interval
      INTEGER   :: N               !order of fit
      REAL(qd)  :: COF(:)          !coefficients
      TYPE(in_struct) :: IO         !in_struct
!local
      REAL(qd)  :: WEIGHTSCALE     !scaling factor for weights
      REAL(qd),ALLOCATABLE:: C(:)
      REAL (qd),ALLOCATABLE :: ALF(:), BETA(:)
      REAL (qd),ALLOCATABLE :: Y(:), T(:,:),A(:,:), W(:)
      INTEGER, PARAMETER :: NMAX=300
      REAL(q)   :: ERR
      INTEGER   :: ITER
      REAL(qd)  :: X1,X2
      REAL(qd)  :: X0(2*N+1)
      INTEGER   :: INFO
      INTEGER   :: I, L, P, NL, IERR,IPRINT, IV
      INTEGER   :: INU             !writing output

      INU=IO%IU6

      X1=1._qd
!initialize current minimization interval
      X2=REAL(RMAX,KIND=qd)

      IF ( SIZE(COF) /= 2*N ) THEN
         IF(INU>=0)WRITE(INU,*)'Error in CONVERGE_FROM_RMAX2R_LSQ, Size of COF inconsistent',&
            SIZE(COF),2*N
         CALL M_exit(); stop
      ENDIF
      ALLOCATE(C(2*N))
      C(1:2*N)=COF(1:2*N) 

!first find optimum coefficients
      CALL NONLINEAR_LSQ_FIT(ITYPE,X2,N,C,IO)

!now approach desired R from above
      IF ( INU>=0 ) WRITE(INU,'(A,F10.2,A,F10.2)')' Approaching R=',R,' from ',X2
      IF ( INU>=0 )WRITE(*,'(A,F10.2,A,F10.2)')' Approaching R=',R,' from ',X2

! determine weighting factor for different error functions
      IF ( ITYPE == 1 ) THEN
! time error function
         WEIGHTSCALE=1._qd
      ELSEIF( ITYPE==2 .OR. ITYPE==3 ) THEN   
!cos error function
         WEIGHTSCALE=PIQ/4._qd
      ELSEIF( ITYPE == 4 ) THEN
!abs error function
         WEIGHTSCALE=PIQ/2._qd
      ELSE
         IF ( INU>=0) WRITE(*,*)'ERROR in NONLINEAR_LSQ_FIT, type of error function is unkown',ITYPE
         CALL M_exit(); stop
      ENDIF 
   
!use these variables for VARPRO
      ALLOCATE(ALF(AMAX), BETA(BMAX))
      ALF=0._qd      
      BETA=0._qd  
      ALF(1:N)=C(1:N)    
      BETA(1:N)=C(N+1:2*N)/WEIGHTSCALE

!initialize sample points  and weights
      ALLOCATE(Y(NMAX),T(NMAX,1),A(NMAX,2*AMAX+2),W(NMAX))
!print info
      IF ( IO%NWRITE ==1 ) THEN
!surpress output of varpro, not recomended!!
         IPRINT = -1 
      ELSEIF (IO%NWRITE == 2 .OR. IO%NWRITE == 3  ) THEN
!normal output of varpro
         IPRINT = 0 
      ELSE
!print complete information
         IPRINT = 1
      ENDIF 

      NL=N  !# of current non-linear parameters
      L=NL  !# of linear parameters
      P=NL  !# of non-vanishing non-linear partial derivatives
      IV=1  !# of linear independent variables

      DO I =1, NDATALSQ
         T(I,1)=LOG(X2)+LOG(X2)*COS(((2*(2*NDATALSQ-I)+1)*PIQ)/(4._qd*NDATALSQ))
         T(I,1)=EXP(T(I,1))
         W(I) = 1.0_qd !WEIGHTS ARE SET TO 1
         IF ( T(I,1) > ACC ) THEN
            Y(I)=1._qd/(T(I,1))
         ELSE
            IF (INU>=0)WRITE(INU,*)'Error in CONVERGE_FROM_RMAX2R_LSQ: ',I,' too small',T(I,1)
            CALL M_exit(); stop
         ENDIF 
      ENDDO

      converge: DO ITER = 1, MAXNEWTON
! find an alternant of error function ( all local extrema)
         CALL VARPRO(N,N,NDATALSQ,NMAX,L+P+2,IV,T,Y,W,ADA,A,IPRINT,ALF,BETA,ITMAXLSQ,&
              ITYPE,IERR,INU,IO%NWRITE)
         IF (IERR<0) THEN
            WRITE(*,*)'ERROR in CONVERGE_FROM_RMAX2R_LSQ',ITER,IERR
            CALL M_exit(); stop
         ENDIF
         
!store new coefficients
         C(1:N)=ALF(1:N)
         C(N+1:2*N)=BETA(1:N)*WEIGHTSCALE         
         IF (ITYPE==3) THEN
            CALL FIND_ALTERNANT_EXTREMA(NUSIN_ERROR_FUNCTION,N,X1,X2,C,X0,NMINALT,INFO)
            IF (INFO<=-1)CALL FIND_ALTERNANT_EXTREMA(NUSIN_ERROR_FUNCTION,N,X1,X2,C,X0,NMAXALT,INFO)
            IF (INFO<=-1)CALL FIND_ALTERNANT_EXTREMA(NUSIN_ERROR_FUNCTION,N,X1,X2,C,X0,NALTCON,INFO)
            IF (INU>=0 .AND. IO%NWRITE>2)WRITE(INU,101)&
               X2,IERR,FINDMAX(N,C,NUSIN_ERROR_FUNCTION,X1,X2)
            IF (INU>=0 .AND. IO%NWRITE>2)WRITE(*,101)&
               X2,IERR,FINDMAX(N,C,NUSIN_ERROR_FUNCTION,X1,X2)
         ELSEIF(ITYPE==4) THEN
            CALL FIND_ALTERNANT_EXTREMA(NUABS_ERROR_FUNCTION,N,X1,X2,C,X0,NMINALT,INFO)
            IF (INFO<=-1)CALL FIND_ALTERNANT_EXTREMA(NUABS_ERROR_FUNCTION,N,X1,X2,C,X0,NMAXALT,INFO)
            IF (INFO<=-1)CALL FIND_ALTERNANT_EXTREMA(NUABS_ERROR_FUNCTION,N,X1,X2,C,X0,NALTCON,INFO)
            IF (INU>=0 .AND. IO%NWRITE>2)WRITE(INU,101)&
               X2,IERR,FINDMAX(N,C,NUABS_ERROR_FUNCTION,X1,X2)
            IF (INU>=0 .AND. IO%NWRITE>2)WRITE(*,101)&
               X2,IERR,FINDMAX(N,C,NUABS_ERROR_FUNCTION,X1,X2)
         ELSEIF(ITYPE==2) THEN
            CALL FIND_ALTERNANT_EXTREMA(NUCOS_ERROR_FUNCTION,N,X1,X2,C,X0,NMINALT,INFO)
            IF (INFO<=-1)CALL FIND_ALTERNANT_EXTREMA(NUCOS_ERROR_FUNCTION,N,X1,X2,C,X0,NMAXALT,INFO)
            IF (INFO<=-1)CALL FIND_ALTERNANT_EXTREMA(NUCOS_ERROR_FUNCTION,N,X1,X2,C,X0,NALTCON,INFO)
            IF (INU>=0 .AND. IO%NWRITE>2)WRITE(INU,101)&
               X2,IERR,FINDMAX(N,C,NUCOS_ERROR_FUNCTION,X1,X2)
            IF (INU>=0 .AND. IO%NWRITE>2)WRITE(*,101)&
               X2,IERR,FINDMAX(N,C,NUCOS_ERROR_FUNCTION,X1,X2)
         ELSE
            CALL FIND_ALTERNANT_EXTREMA(TAU_ERROR_FUNCTION,N,X1,X2,C,X0,NMINALT,INFO)
            IF (INFO<=-1)CALL FIND_ALTERNANT_EXTREMA(TAU_ERROR_FUNCTION,N,X1,X2,C,X0,NMAXALT,INFO)
            IF (INFO<=-1)CALL FIND_ALTERNANT_EXTREMA(TAU_ERROR_FUNCTION,N,X1,X2,C,X0,NALTCON,INFO)
            IF (INU>=0 .AND. IO%NWRITE>2)WRITE(INU,101)&
               X2,IERR,FINDMAX(N,C,TAU_ERROR_FUNCTION,X1,X2)
            IF (INU>=0 .AND. IO%NWRITE>2)WRITE(*,101)&
               X2,IERR,FINDMAX(N,C,TAU_ERROR_FUNCTION,X1,X2)
         ENDIF
         IF (INFO <=-1 ) THEN
            IF (INU>=0 ) WRITE(*,100)NALTCON
            IF (INU>=0 ) WRITE(*,'(F45.32)')X0(:)
            CALL M_exit(); stop
         ENDIF 
!if R is reached close enough, do a final minimization for optimum R and get out
         IF ( REAL(X0(2*N),q) < R ) THEN
            IF ( INU>=0 .AND. IO%NWRITE>2) THEN
               IF (ITYPE==1)THEN
                  WRITE(INU,115)N
               ELSEIF(ITYPE==2)THEN
                  WRITE(INU,116)N
               ELSEIF(ITYPE==3)THEN
                  WRITE(INU,117)N
               ELSEIF(ITYPE==4)THEN
                  WRITE(INU,118)N
               ENDIF
               WRITE(INU,106)(BETA(I),ALF(I),I=1,N)
            ENDIF
            EXIT converge
         ELSE
!set new minimization interval
            X2=X0(2*N)+0.01_qd
         ENDIF 

         DO I =1, NDATALSQ
            T(I,1)=LOG(X2)+LOG(X2)*COS(((2*(2*NDATALSQ-I)+1)*PIQ)/(4._qd*NDATALSQ))
            T(I,1)=EXP(T(I,1))
            W(I) = 1.0_qd !WEIGHTS ARE SET TO 1
            IF ( T(I,1) > EPS ) THEN
               Y(I)=1._qd/(T(I,1))
            ELSE
               IF (INU>=0)WRITE(INU,*)'Error in CONVERGE_FROM_RMAX2R_LSQ: ',I,' too small',T(I,1)
               CALL M_exit(); stop
            ENDIF 
         ENDDO
      ENDDO converge 

      IF ( ITER >  MAXNEWTON ) THEN
         IF ( INU>=0) WRITE(INU,*)'ERROR in CONVERGE_FROM_RMAX2R_LSQ, could not reach desired R'
         CALL M_exit(); stop
      ELSE
         COF(1:2*N)=C(1:2*N) 
      ENDIF 
      DEALLOCATE(C)

100 FORMAT('ERROR in CONVERGE_FROM_RMAX2R_LSQ: no alternant found!',/,&
           'Hint: try to increase sample points for alternant NMAXALT.',/,&
           '      Used parameter for this run:',I10)
101 FORMAT(' R=',F10.2,' done after',I5,' iterations, Error=',E14.6)
115 FORMAT(' LSQ quadrature (exp) found (',I2,'th order):')
116 FORMAT(' LSQ quadrature (cos) found (',I2,'th order):')
117 FORMAT(' LSQ quadrature (sin) found (',I2,'th order):')
118 FORMAT(' LSQ quadrature (abs) found (',I2,'th order):')
106 FORMAT(2F45.32)
   ENDSUBROUTINE CONVERGE_FROM_RMAX2R_LSQ

!****************************************************************************
!
! For each N there's a maximum interval [1,RNMAX].
! Coefficients of [1,R>RNMAX] coincide with the coefficients of [1,RNMAX].
! So we minimize on [1,RNMAX]
!
!****************************************************************************

   SUBROUTINE FIND_RNMAX(R,N,IU0)
      REAL(q) :: R            !current minimization interval
      INTEGER :: N            !order
      INTEGER :: IU0          
!local
      REAL(q) :: RNMAX(32)    !maximum of R for each N
      INTEGER :: I
!maximum of interval for each order N
      DATA RNMAX(1:32)/&      
      8.667,41.54, 146.8,436.0,1153,2807,6373,13749,&
      28387,56502,1.089E+5,2.042E+6,3.737E+5,6.691E+5,1.175E+6,2.027E+6,&
      3.440E+6,5.753E+6,9.491E+6,1.546E+7,2.491E+7,3.969E+7,6.258E+7,9.776E+7,&
      1.513E+8,2.325E+8,3.540E+8,5.353E+8,8.036E+8,1.198E+9,1.775E+9,2.614E+9/

!For a given N all R>R(N)_max yield identical MM and LSQ coefficients
!So we restrict the current interval, if necessary, to RNMAX
      R = MIN(R,RNMAX(N))

   ENDSUBROUTINE FIND_RNMAX

!****************************************************************************
!decide if we start from a higher RMAX(N) and converge to actual R
!****************************************************************************

   SUBROUTINE RMAX2CURRENT_R(R,N,RMAX,IU0)
      REAL(q) :: R            !current minimization interval
      INTEGER :: N            !order
      REAL(q) :: RMAX         
      INTEGER :: IU0          
!local
      REAL(q) :: CHECKPT(9)   !subdivisions
      INTEGER :: IVAL(9)      !supported orders of each subdivision
      INTEGER :: I, J         
!interval intersection (found empirically)
      DATA CHECKPT(1:9)/&
      25._q, 50._q, 106._q, 221._q, 1005._q ,&     
      10005._q, 20005._q, 50005._q, 100000._q/      
!maximum number of grid points supported for intersection
      DATA IVAL(1:9)/&
      12, 14, 16, 20, 24, &
      28, 32, 32, 32/

!adjust RMAX according to R, use the next larger interval
!find correct RMAX
      DO I = 1, SIZE(CHECKPT)
         IF (R-CHECKPT(I)<0) EXIT
      ENDDO
   
!if the interval is not supported, write a warning
      IF ( I>SIZE(CHECKPT) )THEN
         IF ( IU0 >=0 ) WRITE(IU0,'(A,2E14.6)')&
            ' Error in RMAX2CURRENT_R: interval too large',R,CHECKPT(SIZE(CHECKPT))
         CALL M_exit(); stop
      ELSE
!otherwise use the next larger interval as starting point for calculation
         RMAX=CHECKPT(I)
      ENDIF

! check if chosen RMAX supports given N,
      IF ( N > IVAL(I) .AND. I<SIZE(CHECKPT) ) THEN
! if not go to next larger interval, which supports this N
         rmax_n: DO J=I+1, SIZE(CHECKPT)
            IF ( N<= IVAL(J) ) THEN
               RMAX=CHECKPT(J)
               EXIT rmax_n
            ENDIF
         ENDDO rmax_n
      ENDIF 
    
   ENDSUBROUTINE RMAX2CURRENT_R

!============================================================================


ENDMODULE minimax
