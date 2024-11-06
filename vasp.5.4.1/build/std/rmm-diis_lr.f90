# 1 "rmm-diis_lr.F"
!#define debug
# 1 "./symbol.inc" 1 
!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

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


# 130

!
!   wavefunctions: full grid mode
!

# 167

!
!   wavefunctions complex
!




























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

# 3 "rmm-diis_lr.F" 2 
MODULE rmm_diis_lr
  USE prec
CONTAINS
!************************ SUBROUTINE LINEAR_RESPONSE_DIIS **************
!
! this subroutine solves the linear response equation
!    ( H(0) - e(0) S(0) ) |phi(1)> = - |xi>
! where xi is usually calculated to be
!    |xi> = ( H(1) - e(0) S(1) ) |phi(1)> - e(1) S(0) |phi(0)>
! i.e. the perturbation resulting from a change of the Hamiltonian
!
! in principle there is a related variational principle that reads
! < phi(1) | xi > + < xi | phi(1) > + <phi(1)| H(0) - e(0) S(0) |phi(1)>
! which could be optimised as well, but this requires to constrain
! the wavefunctions phi(1) to observe certain orthogonality constraints
!
! in the present implementation an inverse iteration like algorithm
! is the basic step in the linear response solver
! the routine is a variant of the rmm-diis.F routine
!
!  INFO%IALGO   determine type of preconditioning and the algorithm
!    8    TAP preconditioning
!    9    Jacobi like preconditioning
!    (TAP Teter Alan Payne is presently hardcoded)
!  WEIMIN  treshhold for total energy minimisation
!    is the fermiweight of a band < WEIMIN,
!    minimisation will break after a maximum of two iterations
!  EBREAK  absolut break condition
!    intra-band minimisation is stopped if DE is < EBREAK
!  DEPER   intra-band break condition (see below)
!  ICOUEV  number of intraband evalue minimisations
!  DESUM   total change of the variational quantity
!  RMS     norm of residual vector
!  LRESET  reset the wavefunction array entirely
!
!***********************************************************************

  SUBROUTINE LINEAR_RESPONSE_DIIS(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WXI,W0,WDES, &
       LMDIM,CDIJ,CQIJ, RMS,DESUM,ICOUEV, SV, CSHIFT, IU6, IU0, LRESET, IERROR)
    USE prec

    USE wave
    USE wave_high
    USE base
    USE lattice
    USE mpimy
    USE mgrid

    USE nonl_high
    USE hamil
    USE constant
    USE wave_mpi

    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    TYPE (info_struct) INFO
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W             ! LR of orbitals   ( H(0) - e(0) S(0) ) |phi(1)> = - |xi>
    TYPE (wavespin)    WXI           ! |xi>
    TYPE (wavespin)    W0            ! original, unpeturbed orbitals
    TYPE (wavedes)     WDES

    REAL(q)   SV(GRID%MPLWV*2,WDES%NCDIJ) ! local potential
    INTEGER LMDIM  
    REAL(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    REAL(q) CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    REAL(q) DESUM                    ! total change of e(1) related to phi(1)
    REAL(q) RMS                      ! magnitude of the residual vector
    INTEGER ICOUEV                   ! number of H | phi> evaluations
    REAL(q) CSHIFT                   ! complex shift
    INTEGER IU0, IU6                 ! units for output
    LOGICAL LRESET                   ! reset W0
    INTEGER, OPTIONAL :: IERROR      ! return error code
!----- local work arrays
    TYPE (wavedes1)    WDES1         ! descriptor for (1._q,0._q) k-point
    TYPE (wavefun1)    W1(WDES%NSIM) ! current wavefunction
    TYPE (wavefun1)    WTMP          ! temporary
! work arrays
    TYPE(wavefuna) :: W_INI, WOPT
    REAL(q),ALLOCATABLE:: PRECON(:,:)
    COMPLEX(q),ALLOCATABLE::    CWORK1(:)
    COMPLEX(q),ALLOCATABLE:: CHAM(:,:),B(:),CHAM_(:,:),B_(:)
    INTEGER,ALLOCATABLE :: IPIV(:)
    INTEGER :: NSIM,NRES             ! number of bands treated simultaneously
    INTEGER :: LD                    ! leading dimension of CF array
    INTEGER :: NITER                 ! maximum iteration count
    INTEGER :: NODE_ME, IONODE
    INTEGER :: NB(WDES%NSIM)         ! contains a list of bands currently optimized
    REAL(q) :: EVALUE0(WDES%NSIM)    ! eigenvalue e(0)
    COMPLEX(q) :: EVALUE0_C(WDES%NSIM)! version for complex shift
    REAL(q) :: FBREAK(WDES%NSIM)     ! relative break criterion for that band
    REAL(q) :: IT(WDES%NSIM)         ! current iteration for this band
    REAL(q) :: FNORM_                ! norm of residual vector
    REAL(q) :: FNORM(WDES%NSIM)      ! norm of residual vector for each band
    REAL(q) :: ORTH                  ! orthogonality condition
    REAL(q) :: EVAR_                 ! variational quantity
    REAL(q) :: EVAR(WDES%NSIM)       ! variational quantity for each band
    REAL(q) :: SLOCAL                ! average local potential
    REAL(q) :: DE_ATT                ! 1/4 of the total bandwidth
    REAL(q) :: EKIN                  ! kinetic energy
    REAL(q) :: TRIAL                 ! trial step
    REAL(q) :: OTRIAL                ! optimal trial step
    LOGICAL :: LSTOP
    INTEGER :: NP, ISP, NK, NB_DONE, N, IDUMP, ISPINOR, NPRO, M, MM, IND
    INTEGER :: I, ITER, IFAIL, N1, N2
    REAL(q) :: X, X2
    REAL(q) :: ESTART
    COMPLEX(q) :: C

    INFO%IALGO=8

    NODE_ME=0
    IONODE =0

    NODE_ME=WDES%COMM%NODE_ME
    IONODE =WDES%COMM%IONODE

!=======================================================================
!  INITIALISATION:
! maximum  number of iterations
!=======================================================================
    IF (PRESENT(IERROR)) THEN
       IERROR=0
    ENDIF
    NSIM=WDES%NSIM
! at least 6 iterations are required for save convergence
! since there is no other backup algorithm, safety first
    NITER=MAX(INFO%NDAV,6)
    NRES =NITER

    RMS   =0
    DESUM =0
    ESTART=0
    ICOUEV=0
    SLOCAL=MINVAL(REAL(W0%CELTOT(1,1:WDES%NKPTS,1:WDES%ISPIN),q))

    TRIAL = 0.3

    ALLOCATE(PRECON(WDES%NRPLWV,NSIM), &
         CWORK1(NRES), &
         CHAM(NRES,NRES),B(NRES),CHAM_(NRES,NRES),B_(NRES),IPIV(NRES))
    LD=WDES%NRPLWV*NRES*2

    CALL SETWDES(WDES,WDES1,0)
    CALL NEWWAVA(WOPT,  WDES1, NRES*2, NSIM)
    DO NP=1,NSIM
       CALL NEWWAV_R(W1(NP),WDES1)
    ENDDO
    CALL NEWWAV(WTMP,WDES1, .FALSE.)
!=======================================================================
    spin:    DO ISP=1,WDES%ISPIN
    kpoints: DO NK=1,WDES%NKPTS

    IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

    CALL SETWDES(WDES,WDES1,NK)
!=======================================================================
    DE_ATT=ABS(W0%CELTOT(WDES%NB_TOT,NK,ISP)-W0%CELTOT(1,NK,ISP))/2

    IF (INFO%LREAL) THEN
       CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
    ELSE
       CALL PHASE(WDES,NONL_S,NK)
    ENDIF

    NB=0          ! empty the list of bands, which are optimized currently
    NB_DONE=0     ! index the bands already optimised
!=======================================================================
    bands: DO
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
       newband: DO NP=1,NSIM
          IF (NB(NP)==0 .AND.  NB_DONE < WDES%NBANDS ) THEN
             NB_DONE=NB_DONE+1
             N     =NB_DONE
             NB(NP)=NB_DONE
             FBREAK(NP)=0
             IT(NP)  =0

             IDUMP=0
# 190


             IF (NODE_ME /= IONODE) IDUMP=0

             IF (IDUMP>=2) WRITE(*,*)
             IF (IDUMP>=2) WRITE(*,'(I3,1X)',ADVANCE='NO') N

             EVALUE0(NP) =W0%CELEN(N,NK,ISP)
             EVALUE0_C(NP)=EVALUE0(NP) +CMPLX(0.0_q,2.0_q*CSHIFT,q)

!   calculate the preconditioning matrix
! copy eigen energy from CELEN

             CALL SETUP_PRECOND( ELEMENT(W0, WDES1, N, ISP), 8,  IDUMP, PRECON(1,NP), & 
                  EVALUE0(NP)-SLOCAL, DE_ATT )

             IF (LRESET) THEN
                CALL APPLY_PRECOND( ELEMENT( WXI, WDES1, N, ISP), ELEMENT( W, WDES1, N, ISP), &
                     PRECON(1,NP), -1.0_q)
             ENDIF

             CALL SETWAV(W,W1(NP),WDES1,N,ISP)
!-----------------------------------------------------------------------
! FFT of the current trial wave function
!-----------------------------------------------------------------------
             CALL FFTWAV_W1(W1(NP))
             IF (LRESET) THEN
                IF ( INFO%LREAL ) THEN
                   CALL RPRO1(NONLR_S,WDES1,W1(NP))
                ELSE
                   CALL PROJ1(NONL_S,WDES1,W1(NP))
                ENDIF
             ENDIF
          ENDIF
       ENDDO newband
!=======================================================================
! if the NB list is now empty end the bands DO loop
!=======================================================================
       LSTOP=.TRUE.
       W1%LDO  =.FALSE.
       DO NP=1,NSIM
          IF ( NB(NP) /= 0 ) THEN
             LSTOP  =.FALSE.
             W1(NP)%LDO=.TRUE.     ! band not finished yet
             IT(NP) =IT(NP)+1   ! increase iteration count
          ENDIF
       ENDDO
       IF (LSTOP) THEN
          IF (IDUMP>=2) WRITE(*,*)
          EXIT bands
       ENDIF
!=======================================================================
! intra-band minimisation
!=======================================================================
!-----------------------------------------------------------------------
! calculate the vector (H(0)-e(0) S(0)) |phi(1)_opt >
!-----------------------------------------------------------------------
!  residual vector temporarily in upmost storage position (2*NRES)
!  to have uniform stride for result array
       IF (CSHIFT==0) THEN
          CALL HAMILTMU(WDES1, W1, NONLR_S, NONL_S, EVALUE0, &
               &     CDIJ, CQIJ, SV, ISP,  ELEMENTS_SECOND(WOPT, 2*NRES))
       ELSE
          CALL HAMILTMU_C(WDES1, W1, NONLR_S, NONL_S, EVALUE0_C, &
               &     CDIJ, CQIJ, SV, ISP, ELEMENTS_SECOND(WOPT, 2*NRES))

       ENDIF
       i2: DO NP=1,NSIM
          N=NB(NP); ITER=IT(NP); IF (.NOT. W1(NP)%LDO) CYCLE i2

          CALL TRUNCATE_HIGH_FREQUENCY_W1( ELEMENT( WOPT, 2*NRES, NP), .FALSE., INFO%ENINI)

          FNORM_ =0   ! norm of residual
          ORTH   =0   !
          EVAR_  =0
          IND=FLAT_INDEX(WOPT,2*NRES,NP)
          DO ISPINOR=0,WDES%NRSPINORS-1
             DO M=1,WDES1%NGVECTOR
                MM=M+ISPINOR*WDES1%NGVECTOR

!  |R> = H(0)-epsilon S(0) |phi(1)> + | xi >
                C=WOPT%CPTWFP(MM,IND)+WXI%CPTWFP(MM,N,NK,ISP)
!   <R|R>
                FNORM_ =FNORM_+C*CONJG(C)
!   <phi(0)| H(0)-e(0) S(0) |phi(1)> +  <phi(0)| xi >
!   since xi is orthogonal to phi(0), and <phi(0)| H(0)-e(0) S(0)
!   is (0._q,0._q) as well, orth should be (0._q,0._q)
                ORTH   =ORTH+C*CONJG(W0%CPTWFP(MM,N,NK,ISP))
!   variational quantity
!   <phi(1)|xi> + c.c + <phi(1)| H(0)-e(0) S(0)|phi(1)>
                EVAR_  =EVAR_+2*W%CPTWFP(MM,N,NK,ISP)*CONJG(WXI%CPTWFP(MM,N,NK,ISP)) & 
                     +W%CPTWFP(MM,N,NK,ISP)*CONJG(WOPT%CPTWFP(MM,IND))
             ENDDO
          ENDDO

          CALL M_sum_s(WDES%COMM_INB, 3, FNORM_, ORTH, EVAR_, 0._q)

          FNORM(NP)=FNORM_
          IF (IDUMP>=2) WRITE(*,'(E9.2,"R")',ADVANCE='NO') SQRT(ABS(FNORM_))
          IF (IDUMP>=2) WRITE(*,'(E9.2,"O")',ADVANCE='NO') ORTH
          IF (IDUMP>=2) WRITE(*,'(E9.2,"E")',ADVANCE='NO') EVAR_

          IF (ITER==1) THEN
! total norm of error vector at start
             RMS=RMS+WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)* &
                  &      SQRT(ABS(FNORM_))/WDES%NB_TOT
             ESTART=ESTART+WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)*EVAR_
          ELSE
             DESUM =DESUM +WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)*(EVAR_-EVAR(NP))
          ENDIF
          EVAR(NP)=EVAR_
! store variational quantity
          W%CELEN(N,NK,ISP)=EVAR_

! norm of total error vector before start
! norm smaller than EBREAK stop |e -e(app)| < | Residuum |

          IF (ABS(FNORM_)<INFO%EBREAK/10) THEN
             IF (IDUMP>=2) WRITE(*,'("X")',ADVANCE='NO')
             W1(NP)%LDO=.FALSE.
             CYCLE i2
          ENDIF

! break now before filling WOPT%CPTWFP if ITER > NITER
          IF (ITER>NITER) W1(NP)%LDO=.FALSE.
!-----------------------------------------------------------------------
! fill current wavefunctions into work array WOPT%CPTWFP at position ITER
!-----------------------------------------------------------------------
          IF (.NOT. W1(NP)%LDO) CYCLE i2

          CALL W1_COPY(W1(NP), ELEMENT(WOPT, ITER, NP))
          CALL W1_COPY(ELEMENT(WOPT, 2*NRES,NP), ELEMENT(WOPT, NRES+ITER,NP))

          IF (ITER > 1) THEN
! better conditioning for search
! w(iter-1)=w(iter)-w(iter-1)
             CALL W1_DSCAL( ELEMENT( WOPT, ITER-1, NP), -1.0_q)
             CALL W1_DAXPY( ELEMENT( WOPT, ITER, NP), 1.0_q, ELEMENT( WOPT, ITER-1, NP)) 

! gradient(iter-1)=gradient(iter)-gradient(iter-1)
             CALL W1_DSCAL( ELEMENT( WOPT, NRES+ITER-1, NP), -1.0_q)
             CALL W1_DAXPY( ELEMENT( WOPT, NRES+ITER, NP), 1.0_q, ELEMENT( WOPT, NRES+ITER-1, NP)) 
          ENDIF
!***********************************************************************
! inverse interation step
! minimize
!    | ( H - e S) | phi(1) > + | xi > |^ 2  -> min
! in the yet available subspace spanned by the wavefunction stored in CF
! if (1._q,0._q) denotes these wavefunctions as phi(1)_j, and R_j=  (H - e S) phi(1)_j
! the following equation is obtained:
!  sum_ij  b_i* < R_i | R_j > b_j + sum_i b_i* <R_i | xi > + c.c. -> min
! or equivalently
!  sum_j  < R_i | R_j > b_j  = - <R_i | xi >
! the new optimized wavefunction is given by solving this linear
! equation for b
! the scalar product < | > can be evaluated with any metric
!***********************************************************************
          CHAM=0
          B =0

!    A(n2,n1)=    < phi_n2 |  ( H - e S) ( H - e S)  | phi_n1 >
          builda: DO N1=1,ITER
             CALL W1_GEMV( (1._q,0._q), ELEMENTS( WOPT, NRES+N1, NRES+ITER, NP),  ELEMENT( WOPT, NRES+N1, NP), &
                  (0._q,0._q), CWORK1, 1)

             DO N2=N1,ITER
                CHAM(N2,N1)=       (CWORK1(N2-N1+1))
                CHAM(N1,N2)=(CONJG(CWORK1(N2-N1+1)))
             ENDDO
          ENDDO builda

!     B(n1) =   - <R_n1 | xi >= - < phi_n1 | ( H - e S) | xi >
          CALL W1_GEMV( (1._q,0._q), ELEMENTS( WOPT, NRES+1, NRES+ITER, NP),  &
                             ELEMENT( WXI, WDES1, N, ISP), (0._q,0._q), B(1), 1)

          DO N1=1,ITER
             B(N1)=     -(B(N1))
          ENDDO

          IF (ABS(CHAM(1,1))<1E-15) THEN
             IF (PRESENT(IERROR)) THEN
                IERROR=IERROR+1
                W1(NP)%LDO=.FALSE.
                CYCLE i2
             ELSE
                WRITE(0,*) 'internal ERROR: LINEAR_RESPONSE_DIIS matrix is zero, try to call with LRESET',N,NK,ITER,CHAM(1:ITER,1:ITER),B(1:ITER)
             CALL M_exit(); stop
             ENDIF
          ENDIF

          CHAM_=CHAM
          B_ =B
! calculate the solution of sum_j CHAM(i,j) * X(j) = B(i)
! overwrite B by X
          CALL ZGETRF( ITER, ITER, CHAM, NRES, IPIV, IFAIL )
          IF (IFAIL ==0) &
               CALL ZGETRS('N', ITER, 1, CHAM, NRES, IPIV, B, NRES, IFAIL)

          IF (.FALSE.) THEN
! dump the matrix and the solution vector
             IF (NODE_ME==IONODE) THEN
             N2=MIN(10,ITER)
             WRITE(6,*)
             DO N1=1,N2
                WRITE(*,'("m",I3,8E14.7)')N1, CHAM_(N1,1:N2)
             ENDDO
             WRITE(*,'(A4,8E14.7)') 'b', B_(1:N2)

             WRITE(*,*)
             WRITE(*,'(A4,8E14.7)') 'e', B (1:N2)
             ENDIF
          ENDIF

          IF (IFAIL/=0) THEN
             IF (IU6>=0) &
                  WRITE(IU6,219) IFAIL,ITER,N
             IF (IU0>=0) &
                  WRITE(IU0,219) IFAIL,ITER,N
!  try to save things somehow, goto next band
             W1(NP)%LDO=.FALSE.
             CYCLE i2
219          FORMAT('WARNING in EDDRMM_LR: call to GGETRF failed, returncode =',I4,I2,I2)
          ENDIF

          IF (ITER==2 .AND. IDUMP==2) THEN
! write out 'optimal trial step' i.e step which would have minimized
! the residuum
             IF (ITER==2) THEN
                OTRIAL= REAL( 1+B(1)/B(2) ,KIND=q)
                WRITE(*,'(1X,F7.4,"o")',ADVANCE='NO') OTRIAL
             ENDIF
          ENDIF

          IF (IDUMP >= 3) THEN
! set CWORK1(1) to < xi | xi >
             C=W1_DOT( ELEMENT(WXI, WDES1, N, ISP) , ELEMENT(WXI, WDES1, N, ISP))

             DO N1=1,ITER
                DO N2=1,ITER
                   C=C+CONJG(B(N2))*CHAM_(N2,N1)*B(N1)
                ENDDO
                C=C-B_(N1)*CONJG(B(N1))-CONJG(B_(N1))*B(N1)
             ENDDO
! residual after the step
             WRITE(*,'(1X,E9.2,"rs")',ADVANCE='NO') SQRT(ABS(C))
          ENDIF
!=======================================================================
! now performe the trial step (default TRIAL)
!=======================================================================
! W1=0
          CALL W1_DSCAL( W1(NP), 0.0_q)

! W1=W1 + B(I,1) *WOPT(I,NP)
          DO I=1,ITER
             CALL W1_GAXPY( ELEMENT(WOPT, I,NP), B(I), W1(NP))
          ENDDO

! trial step on wavefunction moving from the yet optimised wavefunction
! along the residual vector for that wavefunction
!      -  b_i { ( H(0) - e(0) S(0)) |phi(1)_i> + xi }
! this is somewhat dangerous in the very last step
          DO I=1,ITER
             CALL APPLY_PRECOND( ELEMENT(WOPT, NRES+I, NP), WTMP, PRECON(1,NP))
             CALL W1_GAXPY( WTMP, (-TRIAL*B(I)), W1(NP))
          ENDDO
          CALL ADD_PRECOND( ELEMENT(WXI, WDES1, N, ISP), W1(NP), PRECON(1,NP), -TRIAL)
! transform the wave-function to real space
          CALL FFTWAV_W1(W1(NP))
       ENDDO i2

! calculate results of projection operatores
       CALL W1_PROJALL(WDES1, W1, NONLR_S, NONL_S, NSIM)
!=======================================================================
! break of intra-band-minimisation
!=======================================================================
       i3: DO NP=1,NSIM
          N=NB(NP); ITER=IT(NP); IF (.NOT. W1(NP)%LDO) CYCLE i3

          ICOUEV=ICOUEV+1

          IF (ABS(FNORM(NP))<FBREAK(NP)) W1(NP)%LDO=.FALSE.
          IF (ITER==1) THEN
             FBREAK(NP)=ABS(FNORM(NP))*INFO%DEPER
          ENDIF
! evaluate residual vector in last step as well (just for testing)
          IF (ITER == NITER .AND. .NOT. IDUMP==2) W1(NP)%LDO=.FALSE.
       ENDDO i3

! (1._q,0._q) band just finished ?, set NB(NP) also to 0 and finish everything
       DO NP=1,NSIM
          N=NB(NP)
          IF (.NOT. W1(NP)%LDO .AND. N /=0 ) THEN
             NB(NP)=0
             IF (IDUMP==10) WRITE(*,*)
          ENDIF
       ENDDO
!=======================================================================
! move onto the next Band
!=======================================================================
    ENDDO bands
!=======================================================================
    ENDDO kpoints
    ENDDO spin
!=======================================================================
    IF (PRESENT(IERROR)) THEN
       CALL M_sum_i(WDES%COMM_INTER, IERROR ,1)
       CALL M_sum_i(WDES%COMM_KINTER, IERROR ,1)
    ENDIF
    CALL M_sum_d(WDES%COMM_INTER, RMS, 1)
    CALL M_sum_d(WDES%COMM_KINTER, RMS, 1)

    CALL M_sum_d(WDES%COMM_INTER, DESUM, 1)
    CALL M_sum_d(WDES%COMM_KINTER, DESUM, 1)

    CALL M_sum_i(WDES%COMM_INTER, ICOUEV ,1)
    CALL M_sum_i(WDES%COMM_KINTER, ICOUEV ,1)

    DO NP=1,NSIM
       CALL DELWAV_R(W1(NP))
    ENDDO
    CALL DELWAV(WTMP,.FALSE.)

    CALL DELWAVA(WOPT)
    DEALLOCATE(PRECON,CWORK1,CHAM,B,IPIV,CHAM_,B_)

! WRITE(*,*) 'start energy',ESTART

    RETURN
  END SUBROUTINE LINEAR_RESPONSE_DIIS
END MODULE RMM_DIIS_LR
