# 1 "subrot.F"
!#define dotiming
# 1 "./symbol.inc" 1 
!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define NGZhalf             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define REAL_to_DBLE        convert REAL() to DBLE()
!#define 1                 compile for parallel machine with 1
!------------- end of user part --------------------------------
# 59

# 88

!
!   charge density: full grid mode
!










# 110

!
!   charge density complex
!






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











# 319










# 336

# 3 "subrot.F" 2 


MODULE subrot
  USE prec
  USE pead
  USE dfast
  USE hamil_high
CONTAINS
!************************ SUBROUTINE EDDIAG ****************************
! RCS:  $Id: subrot.F,v 1.10 2003/06/27 13:22:23 kresse Exp kresse $
!
! this subroutine calculates the electronic eigenvalues and
! optionally performes a  sub-space diagonalisation
! i.e. unitary transforms the wavefunctions so that the Hamiltonian
!  becomes diagonal in the subspace spanned by the wavefunctions
! IFLAG:
!  0 only eigenvalues (without  diagonalisation no sub-space matrix)
!  1 only eigenvalues and sub-space matrix (no diagonalisation)
!  2 eigenvalues using diagonalisation of sub-space matrix
!    do not rotate wavefunctions
! 12 eigenvalues using diagonalisation of sub-space matrix
!    no kinetic energy
!  3 eigenvalues and sub-space diagonalisation rotate wavefunctions
! 13 (for 13 no Jacoby algorithm is allowed)
!  4 eigenvalues and sub-space diagonalisation rotate wavefunctions
!    using Loewdin pertubation theory (conserves ordering)
!  5 eigenvalues and sub-space diagonalisation + orthogonalization
!    unfortunately this option turns out to slow down the
!    convergence of IALGO=48
!
! notes on the matrix CHAM used internally
! it is defined as
!  CHAM(n2,n1)   = < phi_n2 | H | phi_n1 >
! however the returned matrix CHAMHF is defined as the complex conjugated
!  CHAMHF(n2,n1) = < phi_n2 | H | phi_n1 >* = < phi_n1 | H | phi_n2 >
! since this is what rot.F uses
!
!  written by gK
!  last update Sep 24 2004 massive cleanup
!  optional arguments:
!  NBANDS_MAX   maximum number of bands for which eigenvalues are
!               calculated (only supported for IFLAG=0)
!  TODO: Martijn new arguments
!
!
!***********************************************************************

  SUBROUTINE EDDIAG(HAMILTONIAN, &
       GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
       LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IU0,EXHF, & 
       CHAMHF,LFIRST,LLAST,NKSTART,NKSTOP,NBANDS_MAX,EXHF_ACFDT)
    USE prec
    USE wave_high
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE hamil_high
    USE constant
    USE jacobi
    USE scala
    USE main_mpi
    USE fock
    USE pseudo
    USE ini
    USE sym_grad
    USE fileio

    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (symmetry) :: SYMM
    INTEGER LMDIM
    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    INTEGER            IFLAG            ! determines mode of diagonalisation
    COMPLEX(q)   SV(GRID%MPLWV,WDES%NCDIJ) ! local potential
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    INTEGER IU0, IU6
    REAL(q) EXHF
! if CHAMHF is present, CHAMHF is added to the subspace matrix at each iteration
! in addition calculation of Fock part is bypassed
    COMPLEX(q), OPTIONAL :: CHAMHF(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
! if LFIRST is supplied CHAMHF is set to CONJG(CHAMHF)-CHAM
! if LLAST  is supplied CHAMHF is set to CONJG(CHAM)
    LOGICAL, OPTIONAL :: LFIRST, LLAST
    INTEGER, OPTIONAL :: NKSTART, NKSTOP     ! start k-point
    INTEGER, OPTIONAL :: NBANDS_MAX  ! maximum band index
    LOGICAL :: LSCAAWARE_LOCAL
    REAL(q), OPTIONAL :: EXHF_ACFDT
    
! local
! work arrays for ZHEEV (blocksize times number of bands)
    INTEGER, PARAMETER :: LWORK=32
    COMPLEX(q)       CWRK(LWORK*WDES%NB_TOT)
    REAL(q)    R(WDES%NB_TOT)
# 106

    REAL(q)    RWORK(7*WDES%NB_TOT), ABSTOL, VL, VU
    INTEGER IWORK(5*WDES%NB_TOT), INFO(WDES%NB_TOT),IL, IU, NB_CALC

! work arrays (do max of 16 strips simultaneously)

    TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point
    TYPE (wavefun1)    W1             ! current wavefunction
    TYPE (wavefuna)    WA             ! array to store wavefunction
    TYPE (wavespin)    WFOCK          ! array to store the Fock (exchange contribution)
    TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character
    TYPE (wavefuna)    WOVL           ! array to hold non local part Q * wave function character
    TYPE (wavefuna)    WHAM           ! array to store accelerations for a selected block
    COMPLEX(q)         CDCHF          ! HF double counting energy
    INTEGER :: ICALL=0                ! number of calls
    COMPLEX(q),ALLOCATABLE,TARGET::  CHAM(:,:),COVL(:,:) ! Hamiltonian and overlap matrix
    TYPE (wavefun1)    WTMP(NSTRIP_STANDARD)

! array for asyncronous data redistribution
    TYPE (REDIS_PW_CTR),POINTER :: H_PW1, H_PW2

    LOGICAL, EXTERNAL :: USEFOCK_CONTRIBUTION

    INTEGER :: NB_TOT, NBANDS, NSTRIP, ISP, NK, N, NP, NPOS, NSTRIP_ACT, &
         NPOS_RED, NSTRIP_RED, IFAIL, MY_NKSTART, MY_NKSTOP, NSIM_LOCAL


    LscaAWARE_LOCAL=LSCAAWARE.AND.(.NOT.(IFLAG==4))
!#define timing
# 139

    IF (PRESENT(CHAMHF) .OR. IFLAG==1) LSCAAWARE_LOCAL=.FALSE.
# 143



    CDCHF=0         ! double counting HF
    IF (PRESENT(EXHF_ACFDT)) EXHF_ACFDT=0
    ICALL=ICALL+1

    NB_TOT=WDES%NB_TOT
    NBANDS=WDES%NBANDS
    IF (PRESENT(NBANDS_MAX) .AND. IFLAG == 0) NBANDS=NBANDS_MAX

    NSTRIP=NSTRIP_STANDARD
 
! allocate work space
    ALLOCATE(W1%CR(GRID%MPLWV*WDES%NRSPINORS))

    IF (.NOT. LscaAWARE_LOCAL) THEN
       ALLOCATE(CHAM(NB_TOT,NB_TOT))
    ELSE
       CALL INIT_scala(WDES%COMM_KIN, NB_TOT)
       ALLOCATE(CHAM(SCALA_NP(),SCALA_NQ()))
    ENDIF

    CALL SETWDES(WDES,WDES1,0)

    CALL NEWWAVA(WHAM, WDES1, NSTRIP)
    CALL NEWWAVA_PROJ(WNONL, WDES1)

    IF (IFLAG==5) THEN
       ALLOCATE(COVL(NB_TOT,NB_TOT))
       CALL NEWWAVA_PROJ(WOVL, WDES1)
    ENDIF

    IF (PRESENT(NKSTART)) THEN
       MY_NKSTART=NKSTART
    ELSE
       MY_NKSTART=1
    ENDIF
    IF (PRESENT(NKSTOP)) THEN
       MY_NKSTOP=NKSTOP
    ELSE
       MY_NKSTOP=WDES%NKPTS
    ENDIF

!=======================================================================
! start with HF part and store results WFOCK
!=======================================================================
    IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN

       CALL ALLOCW(WDES,WFOCK)
       NSIM_LOCAL=(W%WDES%NSIM*2+W%WDES%NB_PAR-1)/W%WDES%NB_PAR
       IF (LPEAD_NO_SCF()) THEN
          DO N=1,NSIM_LOCAL
             CALL NEWWAV(WTMP(N) , WDES1, .FALSE.)
          ENDDO
       ENDIF

       DO ISP=1,WDES%ISPIN
       DO NK=MY_NKSTART,MY_NKSTOP

          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

          CALL SETWDES(WDES,WDES1,NK)

          IF (NONLR_S%LREAL) THEN
             CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
          ELSE
             CALL PHASE(WDES,NONL_S,NK)
          ENDIF

          DO NPOS=1,NBANDS,NSIM_LOCAL
             NSTRIP_ACT=MIN(NBANDS+1-NPOS,NSIM_LOCAL)
             IF (LPEAD_NO_SCF()) THEN
! this is a little bit tricky/dirty
! if local field effects are not taken into account in the pead,
! (i.e. Hamiltonian constructed from initial wavefunctions)
! W_STORE (original wavefunctions) must be passed to FOCK_ACC and WTMP
! (wavefunction on which the action is calculated)
! must be supplied as an additional argument
                DO N=NPOS,NPOS+NSTRIP_ACT-1
                   NP=N-NPOS+1
                   CALL W1_COPY(ELEMENT(W,WDES1,N,ISP),WTMP(NP))
                ENDDO
                CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W_STORE,   &
                     NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                     WFOCK%CPTWFP(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, WTMP, LSYMGRAD=LSYMGRAD)
             ELSE
                IF (PRESENT(EXHF_ACFDT)) THEN
                   CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,  &
                        NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                        WFOCK%CPTWFP(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF,  LSYMGRAD=LSYMGRAD, EXHF_ACFDT=EXHF_ACFDT )
                ELSE
                   CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,  &
                        NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                        WFOCK%CPTWFP(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF,  LSYMGRAD=LSYMGRAD)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       ENDDO
       IF (LSYMGRAD) &
            CALL APPLY_SMALL_SPACE_GROUP_OP( W, WFOCK, NONLR_S, NONL_S,P, T_INFO%NIONS, LATT_CUR, SYMM, CQIJ, .FALSE. , -1, MY_NKSTART)
    ENDIF
!=======================================================================
    spin:  DO ISP=1,WDES%ISPIN
    kpoint: DO NK=MY_NKSTART,MY_NKSTOP

       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

!=======================================================================
       IF (LscaAWARE_LOCAL) CALL INIT_scala(WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))

       CALL SETWDES(WDES,WDES1,NK)
!=======================================================================
!  IFLAG=0 calculate eigenvalues  only
!=======================================================================
       IF (IFLAG==0) THEN
          W%CELEN(:,NK,ISP)=0
          DO N=1,NBANDS
! transform wavefunction to real space
! and calculate eigenvalues calling ECCP, no redistribution !
             CALL SETWAV(W, W1, WDES1, N, ISP) ! allocation for W1%CR 1._q above
             CALL FFTWAV_W1(W1)
!            IF (ASSOCIATED(HAMILTONIAN%AVEC)) THEN
!               CALL ECCP_VEC(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),HAMILTONIAN%AVEC, W1%CELEN)
             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL ECCP_TAU(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),LATT_CUR,HAMILTONIAN%MU(:,ISP),W1%CELEN)    
             ELSE
                CALL ECCP(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W1%CELEN)
             ENDIF
             IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN
                W1%CELEN=W1%CELEN+ &
                     W1_DOT( ELEMENT( W, WDES1, N, ISP), ELEMENT (WFOCK, WDES1, N, ISP))
             ENDIF

             W%CELEN(N,NK,ISP)=W1%CELEN
          ENDDO
! test
          CALL PEAD_EIGENVALUES(W,NK,ISP)
! test
          CYCLE kpoint
       ENDIF

       WA=ELEMENTS(W, WDES1, ISP)
!=======================================================================
!  IFLAG /= 0 calculate Hamiltonian CHAM
!=======================================================================
!  caclulate D |cfin_n> (D = non local strength of PP)
       IF (WDES%DO_REDIS .AND. LASYNC) THEN
          CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW1)
          CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW2)
          DO NPOS=1,NSTRIP
             CALL REDIS_PW_START(WDES, WA%CPTWFP(1,NPOS), NPOS, H_PW1)
          ENDDO
       ENDIF

       CALL OVERL(WDES1, .TRUE.,LMDIM,CDIJ(1,1,1,ISP), WA%CPROJ(1,1),WNONL%CPROJ(1,1))
       DO N=1,NBANDS
          CALL PEAD_ACC_ADD_CPROJ(WNONL%CPROJ(:,N),N,NK,ISP)
       ENDDO

       IF (IFLAG==5) THEN
          CALL OVERL(WDES1, WDES1%LOVERL,LMDIM,CQIJ(1,1,1,ISP), WA%CPROJ(1,1),WOVL%CPROJ(1,1))
       ENDIF

! redistribute the wavefunction characters
       CALL REDISTRIBUTE_PROJ(WA)
       CALL REDISTRIBUTE_PROJ(WNONL)
       IF (IFLAG==5) CALL REDISTRIBUTE_PROJ(WOVL)

       CHAM=0

       strip: DO NPOS=1,NBANDS,NSTRIP
          NSTRIP_ACT=MIN(NBANDS+1-NPOS,NSTRIP)

!  calculate V_{local} |phi> + T | phi >
!  for a block containing NSTRIP wavefunctions

! set Fock contribution
          IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN
             WHAM%CPTWFP(:,1:NSTRIP_ACT)=WFOCK%CPTWFP(:,NPOS:NPOS+NSTRIP_ACT-1,NK,ISP)
          ENDIF

          DO N=NPOS,NPOS+NSTRIP_ACT-1
             NP=N-NPOS+1
             CALL SETWAV(W, W1, WDES1, N, ISP)
             CALL FFTWAV_W1(W1)
             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL HAMILT_LOCAL_TAU(W1, SV, LATT_CUR, HAMILTONIAN%MU, ISP, WHAM%CPTWFP(:,NP), &
               &   USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF)), IFLAG/=12) 
             ELSE
                CALL HAMILT_LOCAL(W1, SV, ISP, WHAM%CPTWFP(:,NP), USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF)), IFLAG/=12)
             ENDIF
             CALL PEAD_ACC_ADD_PW(WHAM%CPTWFP(:,NP),N,NK,ISP)
             IF (WDES%DO_REDIS.AND. LASYNC) CALL REDIS_PW_START(WDES, WHAM%CPTWFP(1,NP), N, H_PW2)
          ENDDO
! redistribute wavefunctions
! after this redistributed up to and including 1...NPOS+NSTRIP_ACT
          IF (WDES%DO_REDIS) THEN
             IF (LASYNC) THEN
                DO N=NPOS,NPOS+NSTRIP_ACT-1
                   NP=N-NPOS+1
                   CALL REDIS_PW_STOP (WDES, WA%CPTWFP(1,N), N, H_PW1)
                   IF (N+NSTRIP<=NBANDS) &
                        CALL REDIS_PW_START(WDES, WA%CPTWFP(1,N+NSTRIP), N+NSTRIP, H_PW1)
                   CALL REDIS_PW_STOP (WDES, WHAM%CPTWFP(1,NP), N, H_PW2)
                ENDDO
             ELSE
                CALL REDISTRIBUTE_PW( ELEMENTS( WA, NPOS, NPOS-1+NSTRIP_ACT))
                CALL REDISTRIBUTE_PW( ELEMENTS( WHAM, 1, NSTRIP_ACT))
             ENDIF
          ENDIF

          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR

          IF (.NOT. LscaAWARE_LOCAL) THEN
             CALL ORTH1('U', &
               WA%CW_RED(1,1),WHAM%CPTWFP(1,1),WA%CPROJ_RED(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1))
          ELSE
             CALL ORTH1_DISTRI('U', &
               WA%CW_RED(1,1),WHAM%CPTWFP(1,1),WA%CPROJ_RED(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1), & 
               WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
          ENDIF
       ENDDO strip
# 374


       IF (WDES%DO_REDIS .AND. LASYNC) THEN
          CALL REDIS_PW_DEALLOC(H_PW1)
          CALL REDIS_PW_DEALLOC(H_PW2)
       ENDIF

       
       IF (.NOT. LscaAWARE_LOCAL) THEN
          CALL M_sum_z(WDES%COMM_KIN,CHAM(1,1),NB_TOT*NB_TOT)
! add lower triangle
          DO N=1,NB_TOT
             DO NP=N+1,NB_TOT
                CHAM(NP,N)=CONJG(CHAM(N,NP))
             ENDDO
          ENDDO
       ENDIF

       IF (PRESENT(CHAMHF).AND.PRESENT(LFIRST)) THEN
          IF (LFIRST) THEN
             CHAMHF(:,:,NK,ISP)=CONJG(CHAMHF(:,:,NK,ISP))-CHAM(:,:)
          ENDIF
          CHAM(:,:)=CHAM(:,:)+CHAMHF(:,:,NK,ISP)
       ENDIF
# 400

!-----------------------------------------------------------------------
! calculate the overlap matrix
!-----------------------------------------------------------------------
       IF (IFLAG==5) THEN
          COVL=(0._q,0._q)

          DO NPOS=1,NB_TOT-NSTRIP_STANDARD_GLOBAL,NSTRIP_STANDARD_GLOBAL
             CALL ORTH1('U',WA%CW_RED(1,1),WA%CW_RED(1,NPOS),WA%CPROJ_RED(1,1), &
                  WOVL%CPROJ_RED(1,NPOS),NB_TOT, &
                  NPOS,NSTRIP_STANDARD_GLOBAL,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
          ENDDO

          CALL ORTH1('U',WA%CW_RED(1,1),WA%CW_RED(1,NPOS),WA%CPROJ_RED(1,1), &
               WOVL%CPROJ_RED(1,NPOS),NB_TOT, &
               NPOS,NB_TOT-NPOS+1,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
          CALL M_sum_z(WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT)
# 419

       ENDIF
# 423

!=======================================================================
! IFLAG =2
!  simply copy eigenvalues
!=======================================================================
       IF (IFLAG==12) THEN
          IF (IU0>=0) CALL DUMP_HAM( "Hamilton matrix",WDES, CHAM)
          CALL M_exit(); stop
       ENDIF

       IF (.NOT. LscaAWARE_LOCAL) THEN
          DO N=1,WDES%NB_TOTK(NK,ISP)
             W%CELTOT(N,NK,ISP)=CHAM(N,N)
          ENDDO
       ENDIF
!=======================================================================
! IFLAG =4 use Loewdin perturbation to get rotation matrix
! this preserves the ordering of the eigenvalues
! MIND: does not work for real matrices
!=======================================================================
       IF (IFLAG==4) THEN
          CALL LOEWDIN_DIAG(WDES%NB_TOTK(NK,ISP), NB_TOT, CHAM)

          CALL ORSP(WDES%NB_TOTK(NK,ISP), NB_TOT, NB_TOT, CHAM)
!test writegamma
!          IF (IU0>=0) CALL WRITEGAMMA(NK, WDES%NB_TOTK(NK,ISP), NB_TOT, CHAM, .TRUE.)
# 451

       ELSE
!=======================================================================
! IFLAG > 1 and IFLAG <4
! diagonalization of CHAM
! we have lots of choices for the parallel version
! this  makes things rather complicated
! to allow for reasonable simple programming, once the diagonalisation
! has been 1._q I jump to line 100
!=======================================================================
          IF (IFLAG==1) GOTO 1000


          DO N=1,NB_TOT
             IF (.NOT. LscaAWARE_LOCAL) THEN
                IF (ABS(AIMAG(CHAM(N,N)))>1E-2_q .AND. IU0>=0) THEN
                   WRITE(IU0,'(A,I5,E14.3)')'WARNING in EDDIAG: sub space matrix is not hermitian',N,AIMAG(CHAM(N,N))
                ENDIF
                CHAM(N,N)= REAL( CHAM(N,N) ,KIND=q)
             ELSE
                CALL BG_CHANGE_DIAGONALE(WDES%NB_TOTK(NK,ISP),CHAM(1,1),IU0)
             ENDIF
          ENDDO


!
! parallel versions
! if fast Jacobi method exists use it (T3D, T3E only)
! use the first line in that case
          

          IFAIL=0



          IF ( LJACOBI .AND. IFLAG /=13 ) THEN
             IF (IU0>=0) WRITE(IU0,*)'jacoby called'
             CALL jacDSSYEV(WDES%COMM_KIN, CHAM(1,1), R, NB_TOT)
             CALL M_sum_z(WDES%COMM_KIN, CHAM(1,1),NB_TOT*NB_TOT)
             CALL M_sum_z(WDES%COMM_KIN, R , NB_TOT)

             GOTO 100
          ENDIF

! use 1 if available in parallel version
          IF ( LscaLAPACK .AND. IFLAG /=5 ) THEN
             IF (.NOT. LscaAWARE_LOCAL) THEN
                CALL pDSSYEX_ZHEEVX(WDES%COMM_KIN, CHAM(1,1), R,  NB_TOT, WDES%NB_TOTK(NK,ISP))
                CALL M_sum_z(WDES%COMM_KIN, CHAM(1,1),NB_TOT*NB_TOT)
             ELSE
                CALL BG_pDSSYEX_ZHEEVX(WDES%COMM_KIN, CHAM(1,1), R,  WDES%NB_TOTK(NK,ISP))
             ENDIF
             
             GOTO 100
          ENDIF


!
!  seriell codes
!
# 534

          IF (IFLAG == 5) THEN
             CALL ZHEGV &
                  (1,'V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT,COVL(1,1),NB_TOT, &
                  R,CWRK,LWORK*NB_TOT, RWORK,  IFAIL)
          ELSE
# 544

             ABSTOL=1E-10_q
             VL=0 ; VU=0 ; IL=0 ; IU=0
             ALLOCATE(COVL(NB_TOT,NB_TOT))
             CALL ZHEEVX( 'V', 'A', 'U', WDES%NB_TOTK(NK,ISP), CHAM(1,1) , NB_TOT, VL, VU, IL, IU, &
                  ABSTOL , NB_CALC , R, COVL(1,1), NB_TOT, CWRK, &
                  LWORK*NB_TOT, RWORK, IWORK, INFO, IFAIL )         
             CHAM=COVL
             DEALLOCATE(COVL)

          ENDIF

! T3D uses a global sum which does not guarantee to give the same results on all nodes
! the following line is required to make the code waterproof (we had problems)
! since we now use a propritary sum (see mpi.F) we should not require
! this broadcast anymore
! 
! CALL M_bcast_z(WDES%COMM, CHAM(1,1), NB_TOT*NB_TOT)

100       CONTINUE
          

          IF (IFAIL/=0) THEN
             WRITE(0,*) 'ERROR in EDDIAG: call to ZHEEV/ZHEEVX/DSYEV/DSYEVX failed! '// &
                  &              'error code was ',IFAIL
             CALL M_exit(); stop
          ENDIF

          DO N=1,WDES%NB_TOTK(NK,ISP)
             W%CELTOT(N,NK,ISP)=R(N)
          ENDDO
       ENDIF
# 578

!=======================================================================
! IFLAG > 2
! rotate wavefunctions
!=======================================================================
       IF (IFLAG==2) GOTO 1000

       IF (.NOT. LscaAWARE_LOCAL) THEN
          CALL LINCOM('F',WA%CW_RED(:,:),WA%CPROJ_RED(:,:),CHAM(1,1), &
            WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), & 
            WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
            WA%CW_RED(:,:),WA%CPROJ_RED(:,:))
       ELSE
          CALL LINCOM_DISTRI('F',WA%CW_RED(1,1),WA%CPROJ_RED(1,1),CHAM(1,1), &
            WDES%NB_TOTK(NK,ISP), & 
            WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
            WDES%COMM_KIN, NBLK )
       ENDIF

# 599


1000   CONTINUE
       
!  back redistribution over bands
       IF (WDES%DO_REDIS) THEN
          CALL REDISTRIBUTE_PROJ( ELEMENTS( W, WDES1, ISP))
          IF (LASYNC) THEN
             W%OVER_BAND=.TRUE.
          ELSE
             CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))
          ENDIF
          ! "redis ok"
       ENDIF

! return Hamiltonian (IFLAG==1) or rotation matrix (IFLAG=3,4)
       IF (PRESENT(LLAST)) THEN       
          IF (LLAST) THEN
             IF (IFLAG==1) CHAM=CONJG(CHAM)
             CHAMHF(:,:,NK,ISP)=CHAM(:,:)
          ENDIF
       ENDIF
# 623


    ENDDO kpoint
    ENDDO spin


    CALL M_sum_z(WDES%COMM_KIN,CDCHF,1)
    CALL M_sum_z(WDES%COMM_KINTER,CDCHF,1)
    EXHF=CDCHF

    IF (PRESENT(EXHF_ACFDT)) THEN
       CALL M_sum_d(WDES%COMM_KIN,EXHF_ACFDT,1)
       CALL M_sum_d(WDES%COMM_KINTER,EXHF_ACFDT,1)
    ENDIF

! need to correct CELTOT at this point so that it is correct on all nodes
    IF (IFLAG==0) CALL MRG_CEL(WDES,W)

    IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
       CALL KPAR_SYNC_CELTOT(WDES,W)
    ENDIF

    IF (IFLAG>=3.AND.(LHFCALC.OR.LUSEPEAD()).AND.WDES%COMM_KINTER%NCPU.GT.1) THEN
       CALL KPAR_SYNC_ALL(WDES,W)
    ENDIF


! deallocation ...
    DEALLOCATE(CHAM,W1%CR)
    CALL DELWAVA(WHAM)
    CALL DELWAVA_PROJ(WNONL)
    IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN 
       CALL DEALLOCW(WFOCK)
       IF (LPEAD_NO_SCF()) THEN
          DO N=1,NSIM_LOCAL
             CALL DELWAV(WTMP(N), .FALSE.)
          ENDDO
       ENDIF
    ENDIF

    IF (IFLAG==5) THEN
       DEALLOCATE(COVL)
       CALL DELWAVA_PROJ(WOVL)
    ENDIF

    IF (PRESENT(CHAMHF).AND.PRESENT(LFIRST)) THEN
       LFIRST=.FALSE.
    ENDIF
    RETURN

  END SUBROUTINE EDDIAG

!***********************************************************************
!
!  for a selected k-points and spin component
! ) calculate onsite ((1._q,0._q) centre) contribution to the Hamilton matrix
!  if LOVERL is .TRUE.
!
!    CCORR(m,n) += <psi_k,n| beta_i> D_ij <beta_j | psi_k,m>
!
! ) if the local potential SV is passed  down the Hamilton matrix
!
!    CCORR(m,n) += <psi_k,n|  SV | psi_k,m>
!
!  is added as well
!  the calling routine must initialize CCORR to (0._q,0._q)
!
!***********************************************************************

  SUBROUTINE ONE_CENTER_BETWEEN_STATES(HAMILTONIAN, LATT_CUR, LOVERL, WDES, W, NK, ISP, LMDIM, &
       CDIJC, CCORR, SV)
    USE prec
    USE wave_high
    USE dfast
    USE lattice
    USE hamil_high

    TYPE (ham_handle)  HAMILTONIAN
    TYPE (latt)        LATT_CUR
    LOGICAL LOVERL
    INTEGER NK, ISP, LMDIM
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    COMPLEX(q) CDIJC(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)  ! (1._q,0._q) centre correction
    COMPLEX(q) ::  CCORR(WDES%NB_TOT,WDES%NB_TOT)
    COMPLEX(q), OPTIONAL  :: SV(WDES%GRID%MPLWV,WDES%NCDIJ)   ! local potential
! local
    COMPLEX(q)   , POINTER :: CPROJ_RED(:,:)
    INTEGER NPOS
    INTEGER NSTRIP, NSTRIP_ACT, NPOS_RED, NSTRIP_RED
    TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character
    TYPE (wavefuna)    WOVL           ! array to hold non local part
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1             ! current wavefunction
    TYPE (wavefuna)    WHAM           ! store Hamiltonian times wavefunction
    TYPE (wavefuna)    WA             ! array to store wavefunction

    NSTRIP=NSTRIP_STANDARD
!-----------------------------------------------------------------------
! non local part
!-----------------------------------------------------------------------
    IF (LOVERL) THEN
       CALL SETWDES(WDES,WDES1,NK)

       CALL NEWWAVA_PROJ(WNONL, WDES1)
       IF (WDES%DO_REDIS) THEN
          CALL NEWWAVA_PROJ(WOVL, WDES1)
          CALL WA_COPY_CPROJ(ELEMENTS(W, WDES1, ISP), WOVL)
       ELSE
          WOVL=ELEMENTS(W, WDES1, ISP)
       ENDIF
       
       CALL OVERL(WDES1, .TRUE., LMDIM, CDIJC(1,1,1,ISP), WOVL%CPROJ(1,1), WNONL%CPROJ(1,1))
       
       IF (WDES%DO_REDIS) THEN
          CALL REDIS_PROJ(WDES1, WDES%NBANDS, WNONL%CPROJ(1,1))
          CALL REDIS_PROJ(WDES1, WDES%NBANDS, WOVL%CPROJ (1,1))
       ENDIF
       
       DO NPOS=1,WDES%NBANDS,NSTRIP
          NSTRIP_ACT=MIN(WDES%NBANDS+1-NPOS,NSTRIP)
          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR
          
          CALL ORTH1('U', &
               W%CPTWFP(1,1,NK,ISP),W%CPTWFP(1,1,NK,ISP),WOVL%CPROJ(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),WDES%NB_TOT, &
               NPOS_RED, NSTRIP_RED, 0,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CCORR(1,1))
!             attention              -

          
       ENDDO
       CALL DELWAVA_PROJ(WNONL)
       IF (WDES%DO_REDIS) CALL DELWAVA_PROJ(WOVL)
    ENDIF
!-----------------------------------------------------------------------
! local part
!-----------------------------------------------------------------------
    IF (PRESENT(SV)) THEN
! allocate work space
       ALLOCATE(W1%CR(WDES%GRID%MPLWV*WDES%NRSPINORS))

       CALL SETWDES(WDES,WDES1,NK)
       CALL NEWWAVA(WHAM, WDES1, NSTRIP)

       WA=ELEMENTS(W, WDES1, ISP)

       strip: DO NPOS=1,WDES%NBANDS,NSTRIP
          NSTRIP_ACT=MIN(WDES%NBANDS+1-NPOS,NSTRIP)

!  calculate V_{local} |phi> + T | phi >
!  for a block containing NSTRIP wavefunctions
          DO N=NPOS,NPOS+NSTRIP_ACT-1
             NP=N-NPOS+1

             CALL SETWAV(W, W1, WDES1, N, ISP)
             CALL FFTWAV_W1(W1)
             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL HAMILT_LOCAL_TAU(W1, SV, LATT_CUR, HAMILTONIAN%MU, ISP,  WHAM%CPTWFP(:,NP), .FALSE., .FALSE.)
             ELSE
                CALL HAMILT_LOCAL(W1, SV, ISP,  WHAM%CPTWFP(:,NP), .FALSE., .FALSE.)
             ENDIF
          ENDDO
! redistribute wavefunctions
! after this redistributed up to and including 1...NPOS+NSTRIP_ACT
          IF (WDES%DO_REDIS) THEN
             CALL REDISTRIBUTE_PW( ELEMENTS( WA, NPOS, NPOS-1+NSTRIP_ACT))
             CALL REDISTRIBUTE_PW( ELEMENTS( WHAM, 1, NSTRIP_ACT))
          ENDIF

          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR

          CALL ORTH1('U', &
               WA%CW_RED(1,1),WHAM%CPTWFP(1,1),WA%CPROJ_RED(1,1), &
               WA%CPROJ_RED(1,NPOS_RED),WDES%NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CCORR(1,1))
!           attention                   ---

       ENDDO strip

       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))

! deallocation ...
       DEALLOCATE(W1%CR)
       CALL DELWAVA(WHAM)

    ENDIF

    CALL M_sum_z(WDES%COMM_KIN,CCORR(1,1),WDES%NB_TOT*WDES%NB_TOT)

  END SUBROUTINE ONE_CENTER_BETWEEN_STATES


!***********************************************************************
!
! read in the file GAMMA and add the density matrix to the
! diagonal density matrix F ((1._q,0._q)-electron occupancies)
! then diagonalize the resulting matrix and rotate the
! (1._q,0._q) electron orbitals
!
!***********************************************************************

  SUBROUTINE ADD_GAMMA_FROM_FILE( WDES, W, IO )
    USE prec
    USE wave_high
    USE dfast
    USE base
    USE fileio
    IMPLICIT NONE
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (in_struct)   IO
! local
    INTEGER NK, ISP, NB, LMDIM
    TYPE (wavedes1)    WDES1
    TYPE (wavefuna)    WA             ! array to store wavefunction
    COMPLEX(q), ALLOCATABLE :: CHAM(:,:)
    REAL(q)    R(WDES%NB_TOT)
! LAPACK
    INTEGER :: IFAIL
    INTEGER, PARAMETER :: LWORK=32
    COMPLEX(q)       CWRK(LWORK*WDES%NB_TOT)
    REAL(q)    RWORK(3*WDES%NB_TOT)

    CALL OPENGAMMA
    CALL READGAMMA_HEAD( WDES%NKPTS , WDES%NB_TOT, IO)

    ALLOCATE(CHAM(WDES%NB_TOT, WDES%NB_TOT))

    DO ISP=1,WDES%ISPIN
    DO NK=1,WDES%NKPTS

       CHAM=0
       CALL READGAMMA(NK, WDES%NB_TOTK(NK,ISP), SIZE(CHAM,1), CHAM, IO)

       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

! add (1._q,0._q)-particle occupancies
       DO NB=1, WDES%NB_TOTK(NK,ISP)
          CHAM(NB,NB)=CHAM(NB,NB)+W%FERTOT(NB,NK,ISP)
       ENDDO
# 867

! diagonalize the (1._q,0._q)-particle matrix
! do not use the X routines, iterative diagonalization is not save
! for density matrix
       IFAIL=0
       
! change sign of density matrix to sort occupied states as lowest states
       CHAM=-CHAM
# 879

       CALL ZHEEV &
            ('V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),WDES%NB_TOT, &
            R,CWRK,LWORK*WDES%NB_TOT, RWORK,  IFAIL)

# 886

! change sign of eigenvalues
       R=-R

       IF (IFAIL/=0) THEN
          WRITE(0,*) 'ERROR in ADD_GAMMA_FROM_FILE: call to ZHEEV/ DSYEV failed! '// &
               &              'error code was ',IFAIL
          CALL M_exit(); stop
       ENDIF

       CALL SETWDES(WDES,WDES1,NK)

       WA=ELEMENTS(W, WDES1, ISP)

!  distribution over plane wave coefficients
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PROJ( ELEMENTS( W, WDES1, ISP))
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))

       CALL LINCOM('F',WA%CW_RED(:,:),WA%CPROJ_RED(:,:),CHAM(1,1), &
            WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), & 
            WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,WDES%NB_TOT, &
            WA%CW_RED(:,:),WA%CPROJ_RED(:,:))

!  back redistribution over bands
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PROJ( ELEMENTS( W, WDES1, ISP))
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))

! updated (1._q,0._q)-electron occupancies
       W%FERTOT(1:WDES%NB_TOTK(NK,ISP),NK,ISP)=R(1:WDES%NB_TOTK(NK,ISP))
    ENDDO
    ENDDO

! just in case sync everything back to all nodes if KPAR is used
    CALL KPAR_SYNC_ALL(WDES,W)

    DEALLOCATE(CHAM)

  END SUBROUTINE ADD_GAMMA_FROM_FILE


!************************ SUBROUTINE EDDIAG_EXACT **********************
!
! this subroutine performs a full diagonalization of the Hamiltonian
! right now it is stupidly implemented since
! the number of bands is increased for all k-points
! than EDDIAG is called
! and finally the number of bands is set back to the original
! value
! this requires a lot of storage but is still convenient for
! GW and RPA calculations
!
!***********************************************************************

  SUBROUTINE EDDIAG_EXACT(HAMILTONIAN, &
       GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
       LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IU0,IU6,EXHF,EXHF_ACFDT)
    USE prec
    USE wave_high
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE hamil_high
    USE main_mpi
    USE pseudo
    USE poscar
    USE ini
    USE choleski
    USE fock
    USE scala
    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (symmetry) ::   SYMM      
    INTEGER LMDIM
    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    INTEGER            IFLAG            ! determines mode of diagonalisation
    COMPLEX(q)   SV(GRID%MPLWV,WDES%NCDIJ) ! local potential
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    INTEGER IU0, IU6
    REAL(q) EXHF
    REAL(q) EXHF_ACFDT
! local
    INTEGER NB_TOT    ! maximum number of plane wave coefficients = number of bands
    TYPE (wavedes)     WDES_TMP
    TYPE (wavespin)    W_TMP
    INTEGER NK, DEGREES_OF_FREEDOM

! just make sure that data distribution is over bands
    CALL REDIS_PW_OVER_BANDS(WDES, W)
! are all bands calculated anyway
    DEGREES_OF_FREEDOM=MAXVAL(WDES%NPLWKP_TOT)
    IF (WDES%LGAMMA) THEN
       DEGREES_OF_FREEDOM=DEGREES_OF_FREEDOM*2-1
    ENDIF

    IF (DEGREES_OF_FREEDOM<=WDES%NB_TOT) THEN
       IFLAG=3
       CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
            LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IU0,EXHF,EXHF_ACFDT=EXHF_ACFDT)
    ELSE

       NB_TOT=((DEGREES_OF_FREEDOM+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
    
       WDES_TMP=WDES
       WDES_TMP%NB_TOT=NB_TOT
       WDES_TMP%NBANDS=NB_TOT/WDES%NB_PAR
       CALL INIT_SCALAAWARE( WDES_TMP%NB_TOT, WDES_TMP%NRPLWV, WDES_TMP%COMM_KIN )
       
       NULLIFY(WDES_TMP%NB_TOTK)
       ALLOCATE(WDES_TMP%NB_TOTK(WDES%NKDIM,2))
! set the maximum number of bands k-point dependent
       DO NK=1,WDES_TMP%NKPTS
          IF (WDES_TMP%LGAMMA) THEN
             WDES_TMP%NB_TOTK(NK,:)=MIN(WDES_TMP%NB_TOT,WDES_TMP%NPLWKP_TOT(NK)*2-1)
          ELSE
             WDES_TMP%NB_TOTK(NK,:)=MIN(WDES_TMP%NB_TOT,WDES_TMP%NPLWKP_TOT(NK))
          ENDIF
       ENDDO
       CALL RESETUP_FOCK_WDES(WDES_TMP, LATT_CUR, LATT_CUR, -1)
       
       CALL ALLOCW(WDES_TMP,W_TMP)

       CALL DUMP_ALLOCATE(IU6)

       W_TMP%FERTOT(:,:,:)=0
       W_TMP%CELTOT(:,:,:)=0
! copy data back to work array
       W_TMP%CPTWFP(:,1:WDES%NBANDS,:,:)    =W%CPTWFP(:,1:WDES%NBANDS,:,:)
       W_TMP%CPROJ(:,1:WDES%NBANDS,:,:) =W%CPROJ(:,1:WDES%NBANDS,:,:)
       W_TMP%CELTOT(1:WDES%NB_TOT,:,:)=W%CELTOT(1:WDES%NB_TOT,:,:)
       W_TMP%FERTOT(1:WDES%NB_TOT,:,:)=W%FERTOT(1:WDES%NB_TOT,:,:)

! random initialization beyond WDES%NBANDS
       CALL WFINIT(WDES_TMP, W_TMP, 1E10_q, WDES%NB_TOT+1) ! ENINI=1E10 not cutoff restriction
! get characters
       CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W_TMP)
! orthogonalization
       CALL ORTHCH(WDES_TMP,W_TMP, WDES%LOVERL, LMDIM,CQIJ)
! and diagonalization
       IFLAG=3
       CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W_TMP,WDES_TMP,SYMM, &
            LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IU0,EXHF,EXHF_ACFDT=EXHF_ACFDT)
! copy data back to original array
       W%CPTWFP(:,1:WDES%NBANDS,:,:)    =W_TMP%CPTWFP(:,1:WDES%NBANDS,:,:)
       W%CPROJ(:,1:WDES%NBANDS,:,:) =W_TMP%CPROJ(:,1:WDES%NBANDS,:,:)
       W%CELTOT(1:WDES%NB_TOT,:,:)=W_TMP%CELTOT(1:WDES%NB_TOT,:,:)
       W%FERTOT(1:WDES%NB_TOT,:,:)=W_TMP%FERTOT(1:WDES%NB_TOT,:,:)
!
       CALL DEALLOCW(W_TMP)
       DEALLOCATE(WDES_TMP%NB_TOTK)

       CALL RESETUP_FOCK_WDES(WDES, LATT_CUR, LATT_CUR, -1)

    ENDIF

  END SUBROUTINE EDDIAG_EXACT

END MODULE subrot


!***********************************************************************
!
! dump a "Hamilton matrix" between the calculated states
!
!***********************************************************************

  
  SUBROUTINE DUMP_HAM( STRING, WDES, CHAM)
    USE wave
    CHARACTER (LEN=*) :: STRING
    TYPE (wavedes)     WDES
    COMPLEX(q) ::  CHAM(WDES%NB_TOT,WDES%NB_TOT)
    INTEGER N1, N2, NPL2
    INTEGER NB_TOT

    NB_TOT=WDES%NB_TOT

    WRITE(*,*) STRING
    NPL2=MIN(10,NB_TOT)
    DO N1=1,NPL2
       WRITE(*,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
    ENDDO
    WRITE(*,*)

    DO N1=1,NPL2
       WRITE(6,2)N1,(AIMAG( CHAM(N1,N2)),N2=1,NPL2)
    ENDDO
    WRITE(*,*)

1   FORMAT(1I2,3X,40F9.5)
!1   FORMAT(1I2,3X,40F14.9)
2   FORMAT(1I2,3X,40F9.5)

  END SUBROUTINE DUMP_HAM


!***********************************************************************
!
! dump a "Hamilton matrix" between the calculated states
! single precision version
!
!***********************************************************************

  
  SUBROUTINE DUMP_HAM_SINGLE( STRING, WDES, CHAM)
    USE wave
    CHARACTER (LEN=*) :: STRING
    TYPE (wavedes)     WDES
    COMPLEX(qs) ::  CHAM(WDES%NB_TOT,WDES%NB_TOT)
    INTEGER N1, N2, NPL2
    INTEGER NB_TOT

    NB_TOT=WDES%NB_TOT

    WRITE(*,*) STRING
    NPL2=MIN(10,NB_TOT)
    DO N1=1,NPL2
       WRITE(*,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
    ENDDO
    WRITE(*,*)

    DO N1=1,NPL2
       WRITE(6,2)N1,(AIMAG( CHAM(N1,N2)),N2=1,NPL2)
    ENDDO
    WRITE(*,*)

1   FORMAT(1I2,3X,40F9.5)
2   FORMAT(1I2,3X,40F9.5)
!2   FORMAT(1I2,3X,40E9.1)

  END SUBROUTINE DUMP_HAM_SINGLE


!***********************************************************************
!
! dump a "Hamilton matrix" between the calculated states
!
!***********************************************************************

  
  SUBROUTINE DUMP_HAM_SELECTED( STRING, WDES, CHAM, NDIM, NBANDS)
    USE wave
    CHARACTER (LEN=*) :: STRING
    TYPE (wavedes)     WDES
    INTEGER NDIM
    INTEGER NBANDS
    COMPLEX(q) ::  CHAM(NDIM,NDIM)
! local
    INTEGER N1, N2, NPL2
    INTEGER NB_TOT

    NB_TOT=WDES%NB_TOT

    WRITE(*,*) STRING
    NPL2=MIN(12,NBANDS)
    DO N1=1,NPL2
       WRITE(*,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
    ENDDO
    WRITE(*,*)

    DO N1=1,NPL2
       WRITE(6,2)N1,(AIMAG( CHAM(N1,N2)),N2=1,NPL2)
    ENDDO
    WRITE(*,*)

1   FORMAT(1I2,3X,24F9.5)
2   FORMAT(1I2,3X,24F9.5)

  END SUBROUTINE DUMP_HAM_SELECTED


!=======================================================================
!
! small routine to dump a distributed matrix
! the descriptor in DESCA must properly describe the matrix
! since RECON_SLICE calls check
!
!=======================================================================


  SUBROUTINE DUMP_HAM_DISTRI( STRING, WDES, CHAM_DISTRI, NB_TOT, DESCA, IU)
    USE wave
    USE scala
    IMPLICIT NONE
    CHARACTER (LEN=*) :: STRING            ! string to dump
    TYPE (wavedes)    :: WDES              ! wave function descriptor
    COMPLEX(q)              :: CHAM_DISTRI(*)    ! distributed matrix
    INTEGER           :: NB_TOT
    INTEGER           :: DESCA(*)          ! distributed matrix descriptor array
    INTEGER           :: IU                ! unit to write to (not dump for IU<0)
! local
    INTEGER, PARAMETER :: NDUMP=16
    COMPLEX(q),ALLOCATABLE  :: CHAM(:,:)
    INTEGER N1, N2
    INTEGER COLUMN_HIGH, COLUMN_LOW

    COLUMN_LOW=1

    COLUMN_HIGH=MIN(COLUMN_LOW+NDUMP-1,NB_TOT)

    ALLOCATE(CHAM(NB_TOT, COLUMN_HIGH-COLUMN_LOW+1))

    CALL RECON_SLICE(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  DESCA, COLUMN_LOW, COLUMN_HIGH)

    CALL M_sum_z(WDES%COMM_KIN, CHAM(1,1), SIZE(CHAM))

    IF (IU>=0) THEN
    WRITE(IU,*) STRING
    DO N1=1,NDUMP
       WRITE(IU,1)N1+COLUMN_LOW-1,(REAL( CHAM(N1,N2+COLUMN_LOW-1) ,KIND=q) ,N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)

    DO N1=1,NDUMP
       WRITE(IU,2)N1+COLUMN_LOW-1,(AIMAG( CHAM(N1,N2+COLUMN_LOW-1)),N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)

    ENDIF

    DEALLOCATE(CHAM)
!1   FORMAT(1I2,3X,40F9.5)
!2   FORMAT(1I2,3X,40F9.5)
1   FORMAT(1I2,3X,40F7.4)
2   FORMAT(1I2,3X,40F7.4)
!1   FORMAT(1I2,3X,40F14.9)

  END SUBROUTINE DUMP_HAM_DISTRI


!************************ SUBROUTINE ORSP   ****************************
!
! this subroutine perfomes a gram-schmidt orthogonalistion of a set
! of vectors (all elements on local node)
! the subroutine uses BLAS 3 calls
!
!***********************************************************************

  SUBROUTINE ORSP(NBANDS, NPL, NRPLWV, CPTWFP)
    USE prec
    IMPLICIT NONE

    INTEGER NBANDS
    INTEGER NPL
    INTEGER NRPLWV
    COMPLEX(q) CPTWFP(NRPLWV,NBANDS)
! local
    COMPLEX(q) CPRO(NBANDS)
    COMPLEX(q), EXTERNAL :: ZDOTC
    REAL(q), EXTERNAL ::  DDOT
    INTEGER I, N
    REAL(q) WFMAG

    IF (NBANDS> NRPLWV) THEN
       WRITE(*,*) 'internal error in ORSP: leading dimension of matrix too small'
       CALL M_exit(); stop
    ENDIF

    CPRO=0
    DO N=1,NBANDS

! normalise the vector

       WFMAG=ZDOTC(NPL,CPTWFP(1,N),1,CPTWFP(1,N),1)
       CALL ZDSCAL(NPL,1/SQRT(WFMAG),CPTWFP(1,N),1)

! now orthogonalise all higher vectors to the
! present vector

       IF (NBANDS/=N ) THEN
          CALL ZGEMV( 'C', NPL , NBANDS-N ,(1._q,0._q) , CPTWFP(1,N+1), &
               &             NRPLWV, CPTWFP(1,N), 1 , (0._q,0._q) ,  CPRO, 1)

          DO I=1,NBANDS
             CPRO(I)=CONJG(CPRO(I))
          ENDDO

          CALL ZGEMM( 'N', 'T' , NPL , NBANDS-N , 1 , -(1._q,0._q) , &
               &             CPTWFP(1,N), NRPLWV , CPRO , NBANDS , &
               &             (1._q,0._q) , CPTWFP(1,N+1) , NRPLWV )
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE ORSP


!***********************************************************************
!
! use Loewdin perturbation to determine a rotation matrix
! this preserves the ordering of the eigenvalues
! MIND: does not work for real matrices
!
!***********************************************************************

  SUBROUTINE LOEWDIN_DIAG(NB_TOT, NBDIM, CHAM)
    USE prec
    IMPLICIT NONE
    INTEGER NB_TOT, NBDIM
    COMPLEX(q) :: CHAM(NBDIM, NB_TOT)
! local
    REAL(q), PARAMETER  :: DIFMAX=0.001_q
    REAL(q) DIFCEL
    INTEGER N1, N2
    COMPLEX(q) :: CROT
    REAL(q) :: FAKT

    IF (NB_TOT>NBDIM) THEN
       WRITE(*,*) 'internal error in LOEWDIN_DIAG: leading dimension of matrix too small'
       CALL M_exit(); stop
    ENDIF


    DO N2=1,NB_TOT
       DO N1=1,N2-1
          DIFCEL= REAL( CHAM(N2,N2)-CHAM(N1,N1) ,KIND=q)
          IF (ABS(DIFCEL)<DIFMAX) THEN
             CROT  =0
          ELSE
             CROT  =CONJG(CHAM(N1,N2))/DIFCEL
             IF (ABS(CROT)>0.1_q) THEN
                FAKT= 0.1_q/ABS(CROT)
                CROT  = CROT*FAKT
             ENDIF
          ENDIF
          CHAM(N2,N1) =-CROT
          CHAM(N1,N2) =-CONJG(CROT)
       ENDDO
    ENDDO
    DO N1=1,NB_TOT
       CHAM(N1,N1)=1
    ENDDO
  END SUBROUTINE LOEWDIN_DIAG


!*******************************************************************
!  calculate the matrix elements of a local potential
!
!  CHAM(i,j)= <psi_i,k| V |psi_j,k> = int psi_i,k*(r) V(r) psi_j,k(r)
!
! between states
! the argument CVTOT must be in real space
! the result is retured in CHAM
!
!*******************************************************************

  SUBROUTINE LOCAL_BETWEEN_STATES( HAMILTONIAN, W, LATT_CUR, P, T_INFO, IRDMAX, LMDIM, &
       GRID_SOFT, GRIDC, GRIDUS, SOFT_TO_C, C_TO_US, CVTOT, CHAM)
    USE prec
    USE wave_high
    USE lattice
    USE poscar
    USE pseudo
    USE pot
    USE pawm
    USE subrot
    USE hamil_high
    USE us
    
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (wavespin)    W
    
    INTEGER  IRDMAX      ! allocation required for augmentation
    TYPE (latt)        LATT_CUR
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (grid_3d)     GRIDC                  ! grid for potentials/charge
    TYPE (grid_3d)     GRID_SOFT              ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDUS                 ! grid for augmentation
    TYPE (transit)     SOFT_TO_C              ! index table between GRID_SOFT and GRIDC
    TYPE (transit)     C_TO_US                ! index table between GRID_SOFT and GRIDC
    COMPLEX(q)  CVTOT(GRIDC%MPLWV,W%WDES%NCDIJ) ! local potential
    INTEGER LMDIM
    COMPLEX(q)       CHAM(W%WDES%NB_TOT,W%WDES%NB_TOT,W%WDES%NKPTS,W%WDES%ISPIN)
! local
    INTEGER ISP, NK
    COMPLEX(q) ::   SV(W%WDES%GRID%MPLWV,W%WDES%NCDIJ)   ! local potential
    COMPLEX(q) :: CDIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    COMPLEX(q) :: CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    INTEGER IRDMAA
    REAL(q)  DISPL(3,T_INFO%NIONS)


    IF (W%WDES%COMM_KINTER%NCPU.NE.1) THEN
       CALL M_stop('LOCAL_BETWEEN_STATES: KPAR>1 not tested (but seems ok), sorry.')
!PK Trivial but callers must be adapted
       CALL M_exit(); stop
    END IF


    DISPL=0

! get  the non local strenght parameters
    CALL SETDIJ_(W%WDES, GRIDC, GRIDUS, C_TO_US, LATT_CUR, P, T_INFO, W%WDES%LOVERL, &
         LMDIM, CDIJ, CQIJ, CVTOT, .FALSE., IRDMAA, IRDMAX, DISPL)

! transform CVTOT to reciprocal space (required by SET_SV)
    DO ISP=1,W%WDES%NCDIJ
       CALL FFT_RC_SCALE(CVTOT(1,ISP),CVTOT(1,ISP),GRIDC)
    ENDDO

! now set SV from CVTOT
    CALL SET_SV( W%WDES%GRID, GRIDC, GRID_SOFT, W%WDES%COMM_INTER, SOFT_TO_C, W%WDES%NCDIJ, SV, CVTOT)
    
    CHAM=0

    DO ISP=1,W%WDES%NCDIJ
       DO NK=1,W%WDES%NKPTS
          CALL ONE_CENTER_BETWEEN_STATES( HAMILTONIAN, LATT_CUR, W%WDES%LOVERL, W%WDES, W, NK, ISP, LMDIM, &
               CDIJ, CHAM(1,1,NK,ISP), SV)
!         CALL DUMP_HAM( "Hamiltonian", W%WDES, CHAM(1,1,NK,ISP))

       ENDDO
    ENDDO
       

  END SUBROUTINE LOCAL_BETWEEN_STATES
