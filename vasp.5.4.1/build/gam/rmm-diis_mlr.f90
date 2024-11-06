# 1 "rmm-diis_mlr.F"
!#define debug


! #define dotiming

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

# 7 "rmm-diis_mlr.F" 2 
MODULE rmm_diis_mlr
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
       LMDIM,CDIJ,CQIJ,RMS,DESUM,ICOUEV,SV,CSHIFT,IU6,IU0,LRESET,IERROR,SW0STORE)
    USE prec

    USE wave
    USE wave_high
    USE base
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE constant
    USE wave_mpi

    USE dfast
    USE hamil, ONLY : APPLY_PRECOND,ADD_PRECOND,HAMILTMU,HAMILTMU_C
# 62

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
    COMPLEX(q), POINTER, OPTIONAL :: SW0STORE(:,:,:,:)
!----- local work arrays
    TYPE (wavedes1)    WDES1         ! descriptor for 1._q k-point
    TYPE (wavefun1)    W1(WDES%NSIM) ! current wavefunction
    TYPE (wavefun1)    WTMP          ! temporary
! work arrays
    TYPE(wavefuna) :: W_INI, WOPT
    REAL(q),ALLOCATABLE:: PRECON(:,:)
    REAL(q),ALLOCATABLE::    CWORK1(:)
    REAL(q),ALLOCATABLE:: CHAM(:,:),B(:),CHAM_(:,:),B_(:)
    INTEGER,ALLOCATABLE :: IPIV(:)
    INTEGER :: NSIM,NRES             ! number of bands treated simultaneously
    INTEGER :: NITER                 ! maximum iteration count
    INTEGER :: NODE_ME,IONODE,NCPU
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
    REAL(q) :: TRIAL                 ! trial step
    REAL(q) :: OTRIAL                ! optimal trial step
    INTEGER :: NP, ISP, NK, NB_DONE, N, IDUMP, ISPINOR, NPRO, M, MM, IND
    INTEGER :: I, ITER, IFAIL, N1, N2
    INTEGER :: NDO
    REAL(q) :: ESTART
    COMPLEX(q) :: C


    TYPE (wavefuna) WACC
    TYPE (wavefun1), ALLOCATABLE :: WSTRIP(:)

    COMPLEX(q), ALLOCATABLE, TARGET :: CW1(:,:)
    COMPLEX(q), POINTER :: CW1_RED(:,:)

    COMPLEX(q), POINTER :: SW0(:,:)
    COMPLEX(q), POINTER :: SW0_RED(:,:)
    LOGICAL LCALCSW0

    COMPLEX(q), ALLOCATABLE :: CWORKN(:,:)
    REAL(q), ALLOCATABLE :: DUMMY(:)
    LOGICAL, ALLOCATABLE :: LDO(:)

    REAL(q), ALLOCATABLE, TARGET :: COVL(:,:)

    REAL(q) EKIN,FAKT,X,X2

    INTEGER NRPLWV_RED,NB_START,NSTRIP,NSTRIP_ACT
    INTEGER NSTRPB,NPOSB1,NPOSB2,NSTRPL,NPOSPL
    LOGICAL DO_REDIS

    INTEGER ierr

# 150


    INFO%IALGO=8


    NODE_ME=WDES%COMM%NODE_ME
    IONODE =WDES%COMM%IONODE
    NCPU   =WDES%COMM_INTER%NCPU
# 162


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

    RMS   =0._q
    DESUM =0._q
    ESTART=0._q
    ICOUEV=0

    SLOCAL=MINVAL(REAL(W0%CELTOT(1,1:WDES%NKPTS,1:WDES%ISPIN)))

    TRIAL = 0.3

    ALLOCATE(PRECON(WDES%NRPLWV,NSIM),CWORK1(NRES),CHAM(NRES,NRES),B(NRES),CHAM_(NRES,NRES),B_(NRES),IPIV(NRES))

    CALL SETWDES(WDES,WDES1,0)

    CALL NEWWAVA(WOPT,WDES1,NRES*2,NSIM)
    DO NP=1,NSIM
       CALL NEWWAV_R(W1(NP),WDES1)
    ENDDO

    CALL NEWWAV(WTMP,WDES1, .FALSE.)


    IF (NCPU/=1) THEN
       DO_REDIS=.TRUE.
       NRPLWV_RED=WDES%NRPLWV/NCPU
    ELSE
       DO_REDIS=.FALSE.
       NRPLWV_RED=WDES%NRPLWV
    ENDIF

    LCALCSW0=.TRUE.
    IF (PRESENT(SW0STORE)) THEN
       IF (ASSOCIATED(SW0STORE)) THEN
! check dimensions
          IF (SIZE(SW0STORE,1)/=WDES%NRPLWV.OR.SIZE(SW0STORE,2)/=WDES%NBANDS &
         &   .OR.SIZE(SW0STORE,3)/=WDES%NKPTS.OR.SIZE(SW0STORE,4)/=WDES%ISPIN) THEN
             WRITE(*,*) 'ERROR: LINEAR_RESPONSE_DIIS: SW0STORE incorrectly dimensioned.'
             CALL M_exit(); stop
          ENDIF
          LCALCSW0=.FALSE.
       ELSE
          ALLOCATE(SW0STORE(WDES%NRPLWV,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN))
       ENDIF
    ELSE
       ALLOCATE(SW0(WDES%NRPLWV,WDES%NBANDS))
       IF (NCPU/=1) THEN
! get pointers for redistribution to "over-plane-wave-coefficients"
          CALL SET_WPOINTER( SW0_RED, NRPLWV_RED, WDES%NB_TOT, SW0(1,1))
       ELSE
          SW0_RED => SW0(:,:)
       ENDIF
    ENDIF

!   NSTRIP=NSTRIP_STANDARD_GLOBAL
    NSTRIP=MIN(4,WDES%NBANDS)

    IF (LCALCSW0) THEN
       ALLOCATE(WSTRIP(NSTRIP),DUMMY(NSTRIP),LDO(NSTRIP))
       DO NP=1,NSTRIP
          CALL NEWWAV(WSTRIP(NP),WDES1,INFO%LREAL)
       ENDDO
       IF (INFO%LREAL) THEN
          ALLOCATE(CWORKN(GRID%MPLWV*WDES%NRSPINORS,NSTRIP))
       ENDIF
    ENDIF

    ALLOCATE(CW1(WDES%NRPLWV,NSIM))
    IF (NCPU/=1) THEN
! get pointers for redistribution to "over-plane-wave-coefficients"
       CALL SET_WPOINTER( CW1_RED, NRPLWV_RED, NSIM*NCPU, CW1(1,1))
    ELSE
       CW1_RED => CW1(:,:)
    ENDIF   

    CALL NEWWAVA(WACC,WDES1,NSIM)

    ALLOCATE(COVL(WDES%NB_TOT,NSIM*NCPU))
# 258



    NONL_S%NK=-1
!=======================================================================
    spin:    DO ISP=1,WDES%ISPIN
    kpoints: DO NK=1,WDES%NKPTS


    IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

!=======================================================================
    CALL SETWDES(WDES,WDES1,NK)

    DE_ATT=ABS(W0%CELTOT(WDES%NB_TOT,NK,ISP)-W0%CELTOT(1,NK,ISP))/2

    IF (INFO%LREAL) THEN
       CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
    ELSE
       CALL PHASE(WDES,NONL_S,NK)
    ENDIF


    IF (PRESENT(SW0STORE)) THEN
       SW0 => SW0STORE(:,:,NK,ISP)
       IF (NCPU/=1) THEN
! get pointers for redistribution to "over-plane-wave-coefficients"
          CALL SET_WPOINTER( SW0_RED, NRPLWV_RED, WDES1%NB_TOT, SW0(1,1))
       ELSE
          SW0_RED => SW0(:,:)
       ENDIF
    ENDIF

    IF (LRESET) THEN
! if the bands distributed "over-plane-wave-coefficients", redistribute to "over-bands"
# 295

       IF (W0%OVER_BAND.AND.DO_REDIS) CALL REDIS_PW(WDES1,WDES1%NBANDS, W0%CPTWFP(1,1,NK,ISP))
# 301

! precompute a quantity that is needed to set up the preconditioner.
! doing it here will save us a lot of communication afterwards in
! steps with LRESET=.FALSE.
       DO NP=1,WDES1%NBANDS
          EKIN=0._q
          DO ISPINOR=0,WDES1%NRSPINORS-1
             DO M=1,WDES1%NGVECTOR
                MM=M+ISPINOR*WDES1%NGVECTOR
                C=W0%CPTWFP(MM,NP,NK,ISP)
                EKIN=EKIN+REAL( C*CONJG(C) ,KIND=q)*WDES1%DATAKE(M,ISPINOR+1)
             ENDDO
          ENDDO
          CALL M_sum_d(WDES1%COMM_INB,EKIN,1)
          IF (EKIN<2.0_q) EKIN=2.0_q
          EKIN=EKIN*1.5_q
! store in AUX
          W0%AUX(NP,NK,ISP)=EKIN
       ENDDO
# 322

! and back to "over-plane-wave-coefficients
       IF (W0%OVER_BAND.AND.DO_REDIS) CALL REDIS_PW(WDES1,WDES1%NBANDS, W0%CPTWFP(1,1,NK,ISP))
# 329

    ENDIF

    IF (LCALCSW0) THEN
! if the bands distributed "over-plane-wave-coefficients", redistribute to "over-bands"
# 336

       IF (W0%OVER_BAND.AND.DO_REDIS) CALL REDIS_PW(WDES1,WDES1%NBANDS, W0%CPTWFP(1,1,NK,ISP))
# 342

! compute S|phi(0)> and store it in SW0
       SW0=(0._q,0._q)
       DO NB_START=1,WDES1%NBANDS,NSTRIP
          NSTRIP_ACT=MIN(NSTRIP,WDES1%NBANDS-NB_START+1)
       
          DO NP=1,NSTRIP_ACT
             CALL W1_COPY(ELEMENT(W0,WDES1,NB_START+NP-1,ISP),WSTRIP(NP)); IF (INFO%LREAL) CALL FFTWAV_W1(WSTRIP(NP))
             SW0(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB_START+NP-1)=WSTRIP(NP)%CPTWFP(1:WDES1%NRSPINORS*WDES1%NGVECTOR)
          ENDDO
       
          IF (INFO%LREAL) THEN
             CWORKN=(0._q,0._q); LDO=.FALSE.; LDO(1:NSTRIP_ACT)=.TRUE.; DUMMY=0._q
             CALL RACCMU(NONLR_S,WDES1,WSTRIP,LMDIM,CQIJ,CQIJ,DUMMY,CWORKN,GRID%MPLWV*WDES1%NRSPINORS,NSTRIP,LDO)
             DO NP=1,NSTRIP_ACT
                DO ISPINOR=0,WDES1%NRSPINORS-1
                   CALL FFTEXT_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORKN(1+ISPINOR*WDES1%GRID%MPLWV,NP),SW0(1+ISPINOR*WDES1%NGVECTOR,NB_START+NP-1),GRID,.TRUE.)
                ENDDO
             ENDDO
          ELSE
             DO NP=1,NSTRIP_ACT 
                CALL VNLACC_ADD(NONL_S,WSTRIP(NP),CQIJ,CQIJ,1,0._q,SW0(:,NB_START+NP-1))
             ENDDO
          ENDIF 
       ENDDO

! redistribute from "over-bands" to "over-plane-wave-coefficients"
# 371

       IF (DO_REDIS) THEN
          CALL REDIS_PW(WDES1,WDES1%NBANDS, W0%CPTWFP(1,1,NK,ISP))
          CALL REDIS_PW(WDES1,WDES1%NBANDS,   SW0(1,1))
       ENDIF
# 380

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
          IF (NB(NP)==0 .AND. NB_DONE<WDES%NBANDS) THEN
             NB_DONE=NB_DONE+1
             N     =NB_DONE
             NB(NP)=NB_DONE
             FBREAK(NP)=0
             IT(NP)  =0

             IDUMP=0
# 403


             IF (NODE_ME /= IONODE) IDUMP=0

             IF (IDUMP>=2) WRITE(*,*)
             IF (IDUMP>=2) WRITE(*,'(I3,1X)',ADVANCE='NO') N

             EVALUE0(NP) =W0%CELEN(N,NK,ISP)
             EVALUE0_C(NP)=EVALUE0(NP) +CMPLX(0.0_q,2.0_q*CSHIFT,q)
# 417

             FAKT=2._q/W0%AUX(N,NK,ISP)
             DO ISPINOR=0,WDES1%NRSPINORS-1
                DO M=1,WDES1%NGVECTOR
                   MM=M+ISPINOR*WDES1%NGVECTOR
                   X=WDES1%DATAKE(M,ISPINOR+1)/W0%AUX(N,NK,ISP)
                   X2=27._q+X*(18._q+X*(12._q+8._q*X))
                   PRECON(MM,NP)=X2/(X2+16._q*X*X*X*X)*FAKT
                ENDDO
             ENDDO

             IF (LRESET) THEN
                CALL APPLY_PRECOND( ELEMENT( WXI, WDES1, N, ISP), ELEMENT( W, WDES1, N, ISP), &
                     PRECON(1,NP), -1.0_q)
             ENDIF

             CALL SETWAV(W,W1(NP),WDES1,N,ISP)

             CALL ZCOPY(WDES1%NRPLWV,W1(NP)%CPTWFP(1),1,CW1(1,NP),1)
             W1(NP)%CPTWFP=>CW1(:,NP)
# 449

          ENDIF
       ENDDO newband
!=======================================================================
! if the NB list is now empty end the bands DO loop
!=======================================================================
       NDO=0
       W1%LDO  =.FALSE.
       DO NP=1,NSIM
          IF ( NB(NP) /= 0 ) THEN
             NDO=NDO+1
             W1(NP)%LDO=.TRUE.     ! band not finished yet
             IT(NP) =IT(NP)+1   ! increase iteration count
          ENDIF
       ENDDO
       CALL M_sum_i(WDES%COMM_INTER,NDO,1)
       IF (NDO==0) THEN
          IF (IDUMP>=2) WRITE(*,*)
          EXIT bands
       ENDIF


       reset: IF (LRESET) THEN
       CALL MPI_barrier(WDES1%COMM_INTER%MPI_COMM,ierr) 
! redistribute from "over-bands" to "over-plane-wave-coefficients"
# 476

       IF (DO_REDIS) CALL REDIS_PW(WDES1,NSIM,CW1(1,1))
# 482

!COVL(i,j)=<phi(0)_i|S|w1_j>
       COVL=0._q
       IF (WDES1%NPL_RED/=0) THEN
          NSTRPB=WDES1%NB_TOT
          NDO=MIN(NSIM*NCPU,WDES1%NB_TOT)
          DO NPOSB1=1,NDO,NSTRPB; DO NPOSB2=1,WDES1%NB_TOT,NSTRPB
             CALL DGEMM('C','N',MIN(NSTRPB,WDES1%NB_TOT-NPOSB2+1),MIN(NSTRPB,NDO-NPOSB1+1),2* WDES1%NPL_RED, &
            &   1._q,SW0_RED(1,NPOSB2),2* NRPLWV_RED,CW1_RED(1,NPOSB1),2* NRPLWV_RED,0._q,COVL(NPOSB2,NPOSB1),WDES1%NB_TOT)
          ENDDO; ENDDO
# 494

       ENDIF
# 500

       CALL M_sum_d(WDES1%COMM_KIN,COVL(1,1),NSIM*NCPU*WDES1%NB_TOT)
# 506

! set COVL(i,j)=0 if |phi(0)_i> is an empty state.
       DO NP=1,WDES1%NB_TOT
          IF (W0%FERTOT(NP,NK,ISP)<1E-4_q) COVL(NP,:)=0._q
       ENDDO
# 513

! subtract |phi(0)_j><phi(0)_j|S|phi(1)_i> from |phi(1)_i>
       IF (WDES1%NPL_RED/=0) THEN
          NDO=MIN(NSIM*NCPU,WDES1%NB_TOT)
          NSTRPB=NDO; NSTRPL=2* WDES1%NPL_RED
          DO NPOSB1=1,NDO,NSTRPB; DO NPOSPL=1,2* WDES1%NPL_RED,NSTRPL
             CALL DGEMM('N','N',MIN(NSTRPL,2* WDES1%NPL_RED-NPOSPL+1),MIN(NSTRPB,NDO-NPOSB1+1),WDES1%NB_TOT, &
            &   -1._q,W0%CPTWFP(NPOSPL,1,NK,ISP),2* NRPLWV_RED,COVL(1,NPOSB1),WDES1%NB_TOT,1._q,CW1_RED(NPOSPL,NPOSB1),2* NRPLWV_RED)
          ENDDO; ENDDO
# 524

       ENDIF
# 530

! redistribute from "over-plane-wave-coefficients" to "over-bands"
       IF (DO_REDIS) CALL REDIS_PW(WDES1,NSIM,CW1(1,1))
# 537

       ENDIF reset
! test_
!   WRITE(*,*) 'pass lr.1.4',NK,ISP
! test_
!-----------------------------------------------------------------------
! FFT of the current trial wave functions
!-----------------------------------------------------------------------
       DO NP=1,NSIM
          IF (.NOT.W1(NP)%LDO) CYCLE
          CALL FFTWAV_W1(W1(NP))
          IF (LRESET) THEN
             IF ( INFO%LREAL ) THEN
                CALL RPRO1(NONLR_S,WDES1,W1(NP))
             ELSE
                CALL PROJ1(NONL_S,WDES1,W1(NP))
             ENDIF
          ENDIF
       ENDDO


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
# 570

         &   CDIJ, CQIJ, SV, ISP, WACC)

       ELSE
          CALL HAMILTMU_C(WDES1, W1, NONLR_S, NONL_S, EVALUE0_C, &
# 577

         &   CDIJ, CQIJ, SV, ISP, WACC)

       ENDIF

! test_
!   WRITE(*,*) 'pass lr.1.6',NK,ISP
! test_
       CALL MPI_barrier(WDES1%COMM_INTER%MPI_COMM,ierr)
# 588

! redistribute from "over-bands" to "over-plane-wave-coefficients"
       IF (DO_REDIS) CALL REDIS_PW(WDES1,NSIM,WACC%CPTWFP(1,1))
# 595

!COVL(i,j)=<phi(0)_i|wacc_j>
       COVL=0._q
       IF (WDES1%NPL_RED/=0) THEN
          NSTRPB=WDES1%NB_TOT
          NDO=MIN(NSIM*NCPU,WDES1%NB_TOT)
          DO NPOSB1=1,NDO,NSTRPB; DO NPOSB2=1,WDES1%NB_TOT,NSTRPB
             CALL DGEMM('C','N',MIN(NSTRPB,WDES1%NB_TOT-NPOSB2+1),MIN(NSTRPB,NDO-NPOSB1+1),2* WDES1%NPL_RED, &
            &   1._q,W0%CPTWFP(1,NPOSB2,NK,ISP),2* NRPLWV_RED,WACC%CW_RED(1,NPOSB1),2* NRPLWV_RED,0._q,COVL(NPOSB2,NPOSB1),WDES1%NB_TOT)
          ENDDO; ENDDO
# 607

       ENDIF
# 613

       CALL M_sum_d(WDES1%COMM_KIN,COVL(1,1),NSIM*NCPU*WDES1%NB_TOT)
# 619

! set COVL(i,j)=0 if |phi(0)_i> is an empty state.
       DO NP=1,WDES%NB_TOT
          IF (W0%FERTOT(NP,NK,ISP)<1E-4_q) COVL(NP,:)=0._q
       ENDDO
# 626

! subtract S|phi(0)_j><phi(0)_j|wacc_i> from |wacc_i>
       IF (WDES1%NPL_RED/=0) THEN
          NDO=MIN(NSIM*NCPU,WDES1%NB_TOT)
          NSTRPB=NDO; NSTRPL=2* WDES1%NPL_RED
          DO NPOSB1=1,NDO,NSTRPB; DO NPOSPL=1,2* WDES1%NPL_RED,NSTRPL
             CALL DGEMM('N','N',MIN(NSTRPL,2* WDES1%NPL_RED-NPOSPL+1),MIN(NSTRPB,NDO-NPOSB1+1),WDES1%NB_TOT, &
            &   -1._q,SW0_RED(NPOSPL,1),2* NRPLWV_RED,COVL(1,NPOSB1),WDES1%NB_TOT,1._q,WACC%CW_RED(NPOSPL,NPOSB1),2* NRPLWV_RED)
          ENDDO; ENDDO
# 637

       ENDIF
# 643

! redistribute from "over-plane-wave-coefficients" to "over-bands"
       IF (DO_REDIS) CALL REDIS_PW(WDES1,NSIM,WACC%CPTWFP(1,1))
# 650


       i2: DO NP=1,NSIM
          N=NB(NP); ITER=IT(NP); IF (.NOT. W1(NP)%LDO) CYCLE i2

          CALL ZCOPY(WDES%NRPLWV,WACC%CPTWFP(1,NP),1,WOPT%CPTWFP(1,2*NRES*NP),1)

!         CALL TRUNCATE_HIGH_FREQUENCY_W1( ELEMENT( WOPT, 2*NRES, NP), .FALSE., INFO%ENINI)

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
!   is 0._q as well, orth should be 0._q
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
! if 1._q denotes these wavefunctions as phi(1)_j, and R_j=  (H - e S) phi(1)_j
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
             CALL W1_GEMV( 1._q, ELEMENTS( WOPT, NRES+N1, NRES+ITER, NP),  ELEMENT( WOPT, NRES+N1, NP), &
                  0._q, CWORK1, 1)

             DO N2=N1,ITER
                CHAM(N2,N1)=       REAL(CWORK1(N2-N1+1),KIND=q)
                CHAM(N1,N2)=REAL((CWORK1(N2-N1+1)),KIND=q)
             ENDDO
          ENDDO builda

!     B(n1) =   - <R_n1 | xi >= - < phi_n1 | ( H - e S) | xi >
          CALL W1_GEMV( 1._q, ELEMENTS( WOPT, NRES+1, NRES+ITER, NP),  &
                             ELEMENT( WXI, WDES1, N, ISP), 0._q, B(1), 1)

          DO N1=1,ITER
             B(N1)=     -REAL(B(N1),KIND=q)
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
          CALL DGETRF( ITER, ITER, CHAM, NRES, IPIV, IFAIL )
          IF (IFAIL ==0) &
               CALL DGETRS('N', ITER, 1, CHAM, NRES, IPIV, B, NRES, IFAIL)

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
                   C=C+(B(N2))*CHAM_(N2,N1)*B(N1)
                ENDDO
                C=C-B_(N1)*(B(N1))-(B_(N1))*B(N1)
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


       CALL MPI_barrier(WDES1%COMM_INTER%MPI_COMM,ierr) 
# 862

! redistribute from "over-bands" to "over-plane-wave-coefficients"
       IF (DO_REDIS) CALL REDIS_PW(WDES1,NSIM,CW1(1,1))
# 869

!COVL(i,j)=<phi(0)_i|S|w1_j>
       COVL=0._q
       IF (WDES1%NPL_RED/=0) THEN
          NSTRPB=WDES1%NB_TOT
          NDO=MIN(NSIM*NCPU,WDES1%NB_TOT)
          DO NPOSB1=1,NDO,NSTRPB; DO NPOSB2=1,WDES1%NB_TOT,NSTRPB
             CALL DGEMM('C','N',MIN(NSTRPB,WDES1%NB_TOT-NPOSB2+1),MIN(NSTRPB,NDO-NPOSB1+1),2* WDES1%NPL_RED, &
            &   1._q,SW0_RED(1,NPOSB2),2* NRPLWV_RED,CW1_RED(1,NPOSB1),2* NRPLWV_RED,0._q,COVL(NPOSB2,NPOSB1),WDES1%NB_TOT)
          ENDDO; ENDDO
# 881

       ENDIF
# 887

       CALL M_sum_d(WDES1%COMM_KIN,COVL(1,1),NSIM*NCPU*WDES1%NB_TOT)
# 893

! set COVL(i,j)=0 if |phi(0)_i> is an empty state.
       DO NP=1,WDES1%NB_TOT
          IF (W0%FERTOT(NP,NK,ISP)<1E-4_q) COVL(NP,:)=0._q
       ENDDO
# 900

! subtract |phi(0)_j><phi(0)_j|S|phi(1)_i> from |phi(1)_i>
       IF (WDES1%NPL_RED/=0) THEN
          NDO=MIN(NSIM*NCPU,WDES1%NB_TOT)
          NSTRPB=NDO; NSTRPL=2* WDES1%NPL_RED
          DO NPOSB1=1,NDO,NSTRPB; DO NPOSPL=1,2* WDES1%NPL_RED,NSTRPL
             CALL DGEMM('N','N',MIN(NSTRPL,2* WDES1%NPL_RED-NPOSPL+1),MIN(NSTRPB,NDO-NPOSB1+1),WDES1%NB_TOT, &
            &   -1._q,W0%CPTWFP(NPOSPL,1,NK,ISP),2* NRPLWV_RED,COVL(1,NPOSB1),WDES1%NB_TOT,1._q,CW1_RED(NPOSPL,NPOSB1),2* NRPLWV_RED)
          ENDDO; ENDDO
# 911

       ENDIF
# 917

! redistribute from "over-plane-wave-coefficients" to "over-bands"
       IF (DO_REDIS) CALL REDIS_PW(WDES1,NSIM,CW1(1,1))
# 924

       DO NP=1,NSIM
          IF (W1(NP)%LDO) THEN
! transform the wave-functions to real space
             CALL FFTWAV_W1(W1(NP))
! store W1(NP)%CPTWFP back into W
             W%CPTWFP(1:WDES1%NRPLWV,NB(NP),NK,ISP)=CW1(1:WDES1%NRPLWV,NP)
          ENDIF
       ENDDO

! project onto projection operators
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

! 1._q band just finished ?, set NB(NP) also to 0 and finish everything
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


! after the above, all functions W0 are distributed "over-plane-wave-coefficients",
! paradoxically this means we have to set OVER_BAND=.TRUE. (this is really stupid)
! OVER_BAND=.TRUE. indicates we have to redistribute the wave functions in order to
! obtain the "over-bands" distribution (see REDIS_PW_OVER_BANDS)
    W0%OVER_BAND=.TRUE.
    IF (.NOT.PRESENT(SW0STORE)) THEN
       DEALLOCATE(SW0)
! redistribute from "over-plane-wave-coefficients" to "over-bands"
       CALL REDIS_PW_OVER_BANDS(WDES,W0)
! after this call OVER_BAND=.FALSE. to indicate that we do not have to redistribute the
! wave functions in order to reach "over-bands" distribution (super confusing)
    ENDIF
    NULLIFY(SW0,SW0_RED)
    IF (LCALCSW0) THEN
       DO NP=1,NSTRIP
          CALL DELWAV(WSTRIP(NP),INFO%LREAL)
       ENDDO
       DEALLOCATE(DUMMY,LDO)
       IF (INFO%LREAL) DEALLOCATE(CWORKN)
    ENDIF
    DEALLOCATE(CW1)
    NULLIFY(CW1_RED)
    CALL DELWAVA(WACC)
    DEALLOCATE(COVL)
# 1021



    RETURN
  END SUBROUTINE LINEAR_RESPONSE_DIIS

END MODULE rmm_diis_mlr
