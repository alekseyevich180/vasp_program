# 1 "scala.F"
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

# 2 "scala.F" 2 
!=======================================================================
! RCS:  $Id: scala.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!
! Module containing wrapper for 1
! written by Gilles de Wijs (gD) and Georg Kresse (gK)
! modified to run on any number of nodes by Dario Alfe
! BG_ routines were added by IBM fellow (bG)
!
! 10.09.2014   gK:
! - removed the T3D support (broken anyway most likely)
! - BG_ routines are now always compiled in
! - added some comments
! - added double check that BG_INIT_scala was called before calling
!    compute routines (simply not sure this is 1._q properly in VASP)
! - removed a lot of unused variables
! - some clean up, but a lot of off and creapy stuff is still around
!   specifically the redundance between BG_INIT_SCALA and INIT_SCALA
!   is odd
!
!=======================================================================

# 129


 MODULE scala
    USE prec
    USE mpimy
# 139

    LOGICAL, PUBLIC :: LELPA      = .FALSE.   ! use ELPA


    LOGICAL, PUBLIC :: LscaLAPACK = .TRUE.
    LOGICAL, PUBLIC :: LscaLU     = .FALSE.   ! use parallel LU  decomposition
    LOGICAL, PUBLIC :: LscaAWARE  = .FALSE.   ! use parallel LU  decomposition
    LOGICAL, PUBLIC :: LscaAWARE_read  = .FALSE. ! LscaAWARE was read from file

!
! customize if required
!
      INTEGER,PRIVATE,PARAMETER :: NB=16               ! blocking factor for distribution of matrices
! P4 optimal, larger values were even slightly better (160) but still
! slower on a Gigabit cluster than ZHEEVX

! customization of matrix diagonalization

      INTEGER,PRIVATE :: NCLUST=48            ! maximum cluster of eigenvector
      REAL(q),PRIVATE :: ABSTOL_DEF=1e-10     ! specifies eigenvector orthogonality tolerance
      REAL(q),PRIVATE :: ORFAC=-1.e0          ! controls reorthogonalisation of eigenvectors
      REAL(qs),PRIVATE:: ORFAC_single=-1.e0   ! controls reorthogonalisation of eigenvectors

!  see manpage for pSSYEVX pZHEEVX:
!    NCLUST is the number of eigenvalues in the largest cluster,
!    where a cluster is defined as a set of close eigenvalues:
!    {W(K),...,W(K+NCLUST-1)|W(J+1)<= W(J)+orfac*norm(A)}
!
! end customization
!
      INTEGER,PARAMETER,PRIVATE :: DTYPE_=1, &
                           CTXT_=2, M_=3, N_=4, MB_=5, NB_=6,  &
                           RSRC_=7, CSRC_ =8, LLD_=9, DLEN_=9

! DESCA( DTYPE_ )The descriptor type.  In most cases,
!                DTYPE_ = BLOCK_CYCLIC_2D = 1
! DESCA( CTXT_ ) The BLACS context handle, indicating
!                the BLACS process grid A is distribu-
!                ted over. The context itself is glo-
!                bal, but the handle (the integer
!                value) may vary.
! DESCA( M_ )    The number of rows in the global
!                array A.
! DESCA( N_ )    The number of columns in the global
!                array A.
! DESCA( MB_ )   The blocking factor used to distribute
!                the rows of the array.
! DESCA( NB_ )   The blocking factor used to distribute
!                the columns of the array.
! DESCA( RSRC_ ) The process row over which the first
!                row of the array A is distributed.
! DESCA( CSRC_ ) The process column over which the
!                first column of the array A is
!                distributed.
! DESCA( LLD_ )  The leading dimension of the local
!                array.  LLD_A >= MAX(1,LOCr(M_A)).
! DLEN_          dimension of DESCA

!
! INIT_scala initializes the following variables
!
! usually VASP uses a "light-weighted" implementation of scalapack, in the
! sense that the matrices are stored in the conventional seriel manner
! and only when the 1 routines are called, the matrices are
! stored distributed and then treated by 1
! this has some disatvantages, for instance the storage requirements
! are huge
! if the variable  LscaAWARE is set in the INCAR, however
! VASP uses distributed matrices from the outside
! however, even in this case, 1 is essentially completely
! "unvisible" to VASP, and all 1 related routines
! should go in here
!
! in several place (bse.F, chi_GG.F) that strategy did not work out,
! and the implementations is directly in the main VASP code
!
! most routines in this module use a global array descriptor DESCSTD
! which is set up by INIT_SCALA. INIT_SCALA can be called
! repeatedly if the matrix dimensions or the processor grid change
! however, all subroutines that have an argument DESCA
! are save and they can be used in different contexts
!
      TYPE scalapack_des
         INTEGER :: MAGIC          ! tells whether MODULE is initialized
         INTEGER :: N              ! dimension of matrix currently treated
! in all routines checked against the actual dimension
         INTEGER :: MPI_COMM       ! mpi communicator used for 1
         INTEGER :: NPROCS         ! number of PEs
         INTEGER :: ICTXT          ! context handle of grid (usually 1 context)
         INTEGER :: NP             ! number of rows on the processor
         INTEGER :: NQ             ! number of cols on the processor
         INTEGER :: NPROW,NPCOL    ! processor grid dimensions
         INTEGER :: MYROW,MYCOL    ! processor coordinates in ps. grid
         INTEGER :: LWWORK,LIWORK,LRWORK,MALLOCPQ ! sizes of scalapack workarrays
      END TYPE scalapack_des
      
! global sca-lapack descriptor
! used as "default" in VASP

      TYPE (scalapack_des), SAVE, PRIVATE :: GSD 
      INTEGER,SAVE   :: DESCSTD( DLEN_ ) ! distributed matrix descriptor array

! few notes: many of the variables are redundant with DESCSTS
!
! GSD%NP= NUMROC(GSD%N,DESCSTD(MB_),GSD%MYROW,0,GSD%NPROW)
! GSD%NQ= NUMROC(GSD%N,DESCSTD(NB_),GSD%MYCOL,0,GSD%NPCOL)
! CALL BLACS_GRIDINFO(DESCSTD(CTXT_),GSD%NPROW,GSD%NPCOL,GSD%MYROW,GSD%MYCOL) to determine GSD%NPROW,...GSD%MYCOL
!
! not sure how lightweighted the evaluation is
!

! size of arrays defined with COMPLEX(q)
# 253

      INTEGER, PARAMETER,PRIVATE :: MCOMP=2


 CONTAINS
!=======================================================================
!
! this subroutine sets LscaAWARE to a default value if not read from
! INCAR
!
!=======================================================================
    SUBROUTINE  INIT_SCALAAWARE( NB_TOT, NRPLWV, COMM )
      IMPLICIT NONE

      INTEGER :: NB_TOT, NRPLWV
      TYPE (communic) COMM
! local
      REAL(q) :: STORAGE_HAM, STORAGE_WAVE
      
# 275

      STORAGE_HAM =16._q*NB_TOT*NB_TOT
      STORAGE_WAVE=16._q*NB_TOT*NRPLWV/ COMM%NCPU

      IF (.NOT. LscaAWARE_read ) THEN
         IF (STORAGE_HAM>STORAGE_WAVE) THEN
            LscaAWARE=LscaLAPACK
         ELSE
            LscaAWARE=.FALSE.
         ENDIF
      ENDIF

    END SUBROUTINE

!=======================================================================
!
! function to determine the blocking factor used in 1
!
!=======================================================================

    SUBROUTINE  QUERRY_SCALA_NB( NB_QUERRY)
      NB_QUERRY=NB
    END SUBROUTINE QUERRY_SCALA_NB

!=======================================================================
!
! messy scalapack initialization, calculate all required workspace
!
! the first routine
!    SUBROUTINE INIT_scala(myCOMM, N )
! sets up a global structure for 1 communication
! this is sufficiently flexible for most parts of VASP
!
! However in the more complicated cases, such as GW, it might
! be necessary to use different contextes for instance
! subdividing the nodes into different communication universes
! in that case
!   SUBROUTINE INIT_scala_DESC(myCOMM, N , DESCA, GS)
! should be used with local DESCA and GS (scalapack_des) descriptors
!
!=======================================================================

    SUBROUTINE INIT_scala(myCOMM, N )
      TYPE(communic) :: myCOMM
      INTEGER N

      CALL INIT_scala_DESC(myCOMM, N , DESCSTD, GSD)

    END SUBROUTINE INIT_scala

!
! here is the context safe version
!

    SUBROUTINE INIT_scala_DESC(myCOMM, N , DESCA, GS)

      USE main_mpi            ! to get world communicator
      IMPLICIT NONE

      INTEGER              :: DESCA( DLEN_ ) ! distributed matrix descriptor array
      TYPE (scalapack_des) :: GS             ! descriptor that should be initialized by the routine
      TYPE(communic) :: myCOMM               ! VASP 1 communicator
      INTEGER N                              ! dimension of matrices to be handled
! local
      INTEGER INFO,NEIG
      INTEGER NP0,NQ0,NN
      INTEGER ierror

      INTEGER,INTRINSIC :: MAX,MIN
      INTEGER,EXTERNAL :: NUMROC,ICEIL
      INTEGER :: FACTORS(4), PWR(4)        
      DATA FACTORS /2, 3, 5, 7/   ! a few prime factors likely to be be used
      INTEGER :: NNPCOL, NNPROW, MR, NDUM, NMIN, I

      EXTERNAL BLACS_PINFO
      EXTERNAL DESCINIT
# 353


! if gs%magic is set and the communicator changes or
! the size N of the matrix changes, then reset the GS structure
      IF (GS%MAGIC == 1789231) THEN
         IF (myCOMM%MPI_COMM /= GS%MPI_COMM .OR. GS%N /=N) THEN
! exit the BLACS context first
            CALL BLACS_GRIDEXIT(GS%ICTXT)
! require re-initialization
            GS%MAGIC=0
         ENDIF
      ENDIF

! if GS%MAGIC is not set to the right value reinitialize the GS structure
      IF (GS%MAGIC /= 1789231) THEN

      GS%MAGIC=1789231
      GS%N=N
      GS%MPI_COMM=myCOMM%MPI_COMM

! number of nodes (for (1._q,0._q) image)
      GS%NPROCS = myCOMM%NCPU

! determine processor grid, according to naive guess
      PWR=0
      MR = GS%NPROCS
      GS%NPROW=1; GS%NPCOL=1
      DO I=1,4
         DO
           IF( MOD(MR,FACTORS(I)) == 0 ) THEN
              PWR(I) = PWR(I) + 1
              MR = MR/FACTORS(I)
           ELSE
              EXIT
           ENDIF
         ENDDO
         NNPROW = PWR(I)/2
         NNPCOL = PWR(I) - NNPROW
         IF(MOD(I,2)==0) THEN         ! to make the grid more uniformly distributed (hopefully)
           NDUM = NNPCOL
           NNPCOL = NNPROW
           NNPROW = NDUM
         ENDIF
         GS%NPROW = GS%NPROW*FACTORS(I)**NNPROW
         GS%NPCOL = GS%NPCOL*FACTORS(I)**NNPCOL
      ENDDO
      IF( GS%NPROW*GS%NPCOL < GS%NPROCS )THEN  ! final adjustment to the total number of procs
         NMIN = MIN(GS%NPROW,GS%NPCOL)
         NDUM = GS%NPROCS/GS%NPROW/GS%NPCOL
         IF(GS%NPROW==NMIN)THEN          ! again, hoping that in this way the grid is more uniform
           GS%NPROW = GS%NPROW * NDUM   ! you might think to a better algorithm ( I dont think it
         ELSE                             ! really matters )
           GS%NPCOL = GS%NPCOL * NDUM
         ENDIF
      ENDIF


! make processor grids
      CALL MPI_barrier( myCOMM%MPI_COMM, ierror )

      CALL PROCMAP( myCOMM, GS%ICTXT, 2, GS%NPROW, GS%NPCOL)

! calculate local size of matrices

      CALL BLACS_GRIDINFO( GS%ICTXT, GS%NPROW, GS%NPCOL, GS%MYROW, GS%MYCOL )

      GS%NP = NUMROC(N,NB,GS%MYROW,0,GS%NPROW)   ! get number of rows on proc
      GS%NQ = NUMROC(N,NB,GS%MYCOL,0,GS%NPCOL)   ! get number of cols on proc

      CALL DESCINIT(DESCA,N,N,NB,NB,0,0,GS%ICTXT,MAX( 1, GS%NP ),INFO ) ! setup descriptor
      IF (INFO.NE.0) THEN
        WRITE(*,*) 'internal error in INIT_SCALA: DESCA, DESCINIT, INFO: ', INFO
        CALL M_exit(); stop
      ENDIF

# 431


! calculate scalapack workspace and allocate, pDSSYEX_ZHEEVX only

! iwork
      GS%LIWORK=6*MAX(N,GS%NPROW*GS%NPCOL+1,4)
! wwork
      NN = MAX( N, NB, 2 )
      NEIG = N                ! number of eigenvectors requested
      IF ((DESCA(MB_).NE.NB).OR.(DESCA(NB_).NE.NB) .OR. &
          (DESCA(RSRC_).NE.0) .OR.(DESCA( CSRC_).NE.0)) THEN
        WRITE(*,*) 'internal error in INIT_SCALA: pSSYEZS_ETC, ERROR'
        CALL M_exit(); stop
      ENDIF
      NP0 = NUMROC( NN, NB, 0, 0, GS%NPROW )
      NQ0 = MAX( NUMROC( NEIG, NB, 0, 0, GS%NPCOL ), NB )
# 451

      GS%LRWORK=4*N+MAX( 5*NN, NP0 * NQ0 )+ICEIL( NEIG, GS%NPROW*GS%NPCOL)*NN+ &
         NCLUST*N
      GS%LWWORK=N + MAX((NP0 +NQ0 +NB ) * NB, 3)


      GS%MALLOCPQ=MAX(GS%NP*GS%NQ,1)

! gK: dump global variables to see what has been set
!      WRITE(*,*) GS%NPROCS         ,' number of PEs'
!      WRITE(*,*) GS%ICTXT          ,' context handle of grid'
!      WRITE(*,*) GS%NPROW,GS%NPCOL    ,' processor grid dimensions'
!      WRITE(*,*) GS%MYROW,GS%MYCOL    ,' processor coordinates in ps. grid'
!      WRITE(*,*) GS%NP             ,' number of rows on the processor'
!      WRITE(*,*) GS%NQ             ,' number of cols on the processor'
!      WRITE(*,*) GS%LWWORK,GS%LIWORK,GS%LRWORK,GS%MALLOCPQ ,' sizes of scalapack workarrays'

      ENDIF

    END SUBROUTINE INIT_scala_DESC


!=======================================================================
!
! any routine should either call INIT_scala or CHECK_scala
! CHECK_scala, gracefully aborts if the GSD%N is incompatible
! to the supplied value N
!
!=======================================================================

    SUBROUTINE CHECK_scala (N, STRING)
      INTEGER N
      CHARACTER (LEN=*) :: STRING
      IF (N /= GSD%N) THEN
         WRITE(0,*) 'internal error in VASP: BG_INIT_SCALA_BG was not called beforehand for this k-point'
         WRITE(0,*) '  CHECK_scala was called by ',STRING
         CALL M_exit(); stop
      ENDIF

    END SUBROUTINE CHECK_scala

!=======================================================================
!
! routine to querry the GSD%NP and GSD%NQ variables
!
!=======================================================================

    FUNCTION SCALA_NP()
      INTEGER SCALA_NP
      SCALA_NP=GSD%NP
    END FUNCTION SCALA_NP


    FUNCTION SCALA_NQ()
      INTEGER SCALA_NQ
      SCALA_NQ=GSD%NQ
    END FUNCTION SCALA_NQ

!=======================================================================
!
!     this subroutine maps processors onto a grid
!     and allows to performe 1 operations on a subset of nodes
!     taken from BLACS example code (written by Clint Whaley 7/26/94)
!     modified by gD and gK
!
!=======================================================================

    SUBROUTINE PROCMAP(VCOMM,CONTEXT, MAPPING, NPROW, NPCOL)
      USE main_mpi     ! to get 1 communication handle
      IMPLICIT NONE
      TYPE(communic) :: VCOMM

!     .. Scalar Arguments ..
      INTEGER :: CONTEXT, MAPPING, NPROW, NPCOL
      INTEGER :: IMAP(NPROW, NPCOL)
      INTEGER ierror
!     ..
!
!  Purpose
!  =======
!  PROCMAP maps NPROW*NPCOL onto an VASP-Communicator
!  It is a nice system independent implementation
!
!  Arguments
!  =========
!
!  CONTEXT      (output) INTEGER
!               This integer is used by the BLACS to indicate a context.
!               A context is a universe where messages exist and do not
!               interact with other contexts messages.  The context includes
!               the definition of a grid, and each processs coordinates in it.
!
!  MAPPING      (input) INTEGER
!               Way to map processes to grid.  Choices are:
!               1 : row-major natural ordering
!               2 : column-major natural ordering
!
!  NPROW        (input) INTEGER
!               The number of process rows the created grid
!               should have.
!
!  NPCOL        (input) INTEGER
!               The number of process columns the created grid
!               should have.
!
!  =====================================================================
!
      INTEGER,EXTERNAL :: BLACS_PNUM
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDMAP
      INTEGER TMPCONTXT, NPROCS, I, J, K
      INTEGER BLACS_IDS(VCOMM%NCPU), MY_BLACS_ID
      INTEGER NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP
!
!     See how many processes there are in the system
!     NPROCS here refers to all processors being used

      NPROCS=VCOMM%NCPU
!
!     Temporarily map all processes into 1 x NPROCS grid
!
      TMPCONTXT=VCOMM%MPI_COMM
      CALL MPI_barrier( VCOMM%MPI_COMM, ierror )
      CALL BLACS_GRIDINIT( TMPCONTXT, 'Row', 1, NPROCS )
      CALL MPI_barrier( VCOMM%MPI_COMM, ierror )
      
!
!     If we want a row-major natural ordering
!
      CALL BLACS_GRIDINFO( TMPCONTXT, NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP)

      MY_BLACS_ID= BLACS_PNUM(TMPCONTXT, MYROW_TMP, MYCOL_TMP)
      CALL MPI_ALLGATHER(MY_BLACS_ID,1,MPI_INTEGER, &
     &     BLACS_IDS(1),1,MPI_INTEGER,VCOMM%MPI_COMM,IERROR)
      IF (IERROR.NE.0) CALL M_stop_ierr('ERROR allgather in procmap',ierror)

      K=1
      IF (MAPPING .EQ. 1) THEN
         DO I = 1, NPROW
            DO J = 1, NPCOL
               IMAP(I, J) = BLACS_IDS( K )
               K = K + 1
            END DO
         END DO
!
!     If we want a column-major natural ordering
!
      ELSE IF (MAPPING .EQ. 2) THEN
         DO J = 1, NPCOL
            DO I = 1, NPROW
               IMAP(I, J) = BLACS_IDS( K )
               K = K + 1
            END DO
         END DO
      ELSE
         WRITE(*,*) 'INIT_SCALA, PROCMAP: unknown mapping, STOP'
         CALL M_exit(); stop
      END IF
!
!     Free temporary context
!
      CALL MPI_barrier( VCOMM%MPI_COMM, ierror )
      CALL BLACS_GRIDEXIT(TMPCONTXT)
!      WRITE (*,*) 'MAPPING',NPCOL,NPROW
!
!     Apply the new mapping to form desired context
!
      CONTEXT=VCOMM%MPI_COMM
      CALL BLACS_GRIDMAP( CONTEXT, IMAP, NPROW, NPROW, NPCOL )

!      WRITE(*,'("local map",16I3)') ((IMAP(I,J),I=1,NPROW),J=1,NPCOL)

      RETURN
    END SUBROUTINE PROCMAP

!=======================================================================
!
! the following routines rely on the global GSD communicator
! they however re-initialize GSD as well, so should be safe
! against changes of the matrix dimension
!
!=======================================================================
!
! this subroutine determines U^-1: A= U+ U and invert U
!   1. calls INIT_scala
!   2. calls the setup of distributed the distributed matrix
!   3. calls PDPOTRF and PDTRTRI (gamma) or PZPOTRF and PZTRTRI
!   4. calls reconstruction of distributed data into patched matrix
! note: the sum over processors of the patched matrix is not realised
!       in this module
!
!=======================================================================


    SUBROUTINE pPOTRF_TRTRI (COMM, AMATIN, NDIM, N)
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER NDIM           ! leading dimension of AMATIN
      INTEGER N              ! NxN matrix to be distributed
      COMPLEX(q)    AMATIN(NDIM,N) ! input/output matrix
! local variables
      INTEGER INFO

      COMPLEX(q), ALLOCATABLE ::  A(:)

      INTEGER,EXTERNAL :: NUMROC
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER,EXTERNAL :: BLACS_PNUM

      CALL INIT_scala(COMM, N)

!-----------------------------------------------------------------------
! allocate workarray (on T3D allocated in init_T3D)
!-----------------------------------------------------------------------
      ALLOCATE(A(GSD%MALLOCPQ))
!-----------------------------------------------------------------------
! do actual calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN
!       get parts of the matrix into the local arrays
        CALL DISTRI(AMATIN,NDIM,N,A,DESCSTD)

        INFO=0
# 677

        CALL PZPOTRF( 'U', N, A, 1, 1, DESCSTD, INFO )


        IF (INFO.NE.0) THEN
          WRITE(*,*) 'pPOTRF_TRTRI, POTRF, INFO:',INFO
          WRITE(*,*) 'STOP'
          CALL M_exit(); stop
        ENDIF
# 688

         CALL PZTRTRI &

     &    ('U', 'N', N, A, 1, 1, DESCSTD, INFO)
        IF (INFO.NE.0) THEN
          WRITE(*,*) 'pPOTRF_TRTRI, TRTRI, INFO:',INFO
          CALL M_exit(); stop
        ENDIF

      CALL RECON(AMATIN,NDIM,N,A,DESCSTD)
   ENDIF calc

!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(A)


      RETURN
    END SUBROUTINE pPOTRF_TRTRI


!=======================================================================
!
! this subroutine determines A^-1: A= U+ U and invert (U+ U)
!   1. calls INIT_scala
!   2. calls the setup of distributed the distributed matrix
!   3. calls PDPOTRF and PDTRTRI (gamma) or PZPOTRF and PZTRTRI
!   4. calls reconstruction of distributed data into patched matrix
! note: the sum over processors of the patched matrix is not realised
!       in this module
!
!=======================================================================


    SUBROUTINE pPOTRF_POTRI (COMM, AMATIN, NDIM, N)
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER NDIM           ! leading dimension of AMATIN
      INTEGER N              ! NxN matrix to be distributed
      COMPLEX(q)    AMATIN(NDIM,N) ! input/output matrix
! local variables
      INTEGER INFO

      COMPLEX(q), ALLOCATABLE ::  A(:)

      INTEGER,EXTERNAL :: NUMROC
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER,EXTERNAL :: BLACS_PNUM

      CALL INIT_scala(COMM, N)
!-----------------------------------------------------------------------
! allocate workarray
!-----------------------------------------------------------------------
      ALLOCATE(A(GSD%MALLOCPQ))
!-----------------------------------------------------------------------
! do actual calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN
!       get parts of the matrix into the local arrays
        CALL DISTRI(AMATIN,NDIM,N,A,DESCSTD)

        INFO=0
# 755

        CALL PZPOTRF( 'U', N, A, 1, 1, DESCSTD, INFO )


        IF (INFO.NE.0) THEN
          WRITE(*,*) 'pPOTRF_TRTRI, POTRF, INFO:',INFO
          WRITE(*,*) 'STOP'
          CALL M_exit(); stop
        ENDIF
# 766

         CALL PZPOTRI &

     &    ('U', N, A, 1, 1, DESCSTD, INFO)
        IF (INFO.NE.0) THEN
          WRITE(*,*) 'pPOTRF_POTRI, POTRI, INFO:',INFO
          CALL M_exit(); stop
        ENDIF

      CALL RECON(AMATIN,NDIM,N,A,DESCSTD)
   ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(A)

      RETURN
    END SUBROUTINE pPOTRF_POTRI


!=======================================================================
!
! call to pZHEEVX respectively pDSYEVX
! i.e. diagonalization of an symmetric or hermitian matrix
!
!=======================================================================

    SUBROUTINE pDSSYEX_ZHEEVX(COMM, AMATIN, W, NDIM, N, ABSTOL_, AMATIN_SCALA)
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N                          ! NxN matrix to be distributed
      INTEGER NDIM                       ! leading dimension of matrix
      COMPLEX(q)    AMATIN(NDIM,N)             ! input/output matrix
      REAL(q) W(N)                       ! eigenvalues
      REAL(q), OPTIONAL      ::  ABSTOL_ ! tolerance
      COMPLEX(q), OPTIONAL, TARGET ::  AMATIN_SCALA(GSD%MALLOCPQ) ! alternative input/output matrix (1 format)
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(q) :: GAP(COMM%NCPU)
      REAL(q) :: ZERO=0.0_q
      INTEGER INFO
      INTEGER M,NZ
      INTEGER I2
      REAL(q) :: ABSTOL

      COMPLEX(q), ALLOCATABLE,TARGET ::  A(:) ! A matrix to be diagonalized (allocatable)
      COMPLEX(q), POINTER  ::  AP(:)          ! A matrix to be diagonalized (pointer) either to A or AMATIN_SCALA
      COMPLEX(q), ALLOCATABLE    ::  Z(:)     ! Z eigenvector matrix
      REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
      COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array

      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE
# 826


      CALL INIT_scala(COMM, N)
      
      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE
      
      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
         WRITE(*,*) "pDSSYEX_ZHEEVX: too many processors ", &
              GSD%NPROW*GSD%NPCOL,COMM%NCPU
         CALL M_exit(); stop
      ENDIF
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      IF (PRESENT(AMATIN_SCALA)) THEN
         AP=>AMATIN_SCALA
      ELSE
         ALLOCATE(A(GSD%MALLOCPQ))
         AP=>A
      ENDIF
      ALLOCATE(Z(GSD%MALLOCPQ))
      ALLOCATE(IWORK(GSD%LIWORK))
      ALLOCATE(WWORK(GSD%LWWORK))
      ALLOCATE(RWORK(GSD%LRWORK))

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

!     get parts of the matrix into the local arrays
         IF (.NOT. PRESENT(AMATIN_SCALA)) THEN
            CALL DISTRI(AMATIN, NDIM, N,A(1),DESCSTD)
         ENDIF
         INFO=0

!     call scalapack routine

# 878


        IF (.NOT. LELPA) THEN
# 885

         CALL PZHEEVX( 'V', 'A', 'U', N, AP(1), 1, 1, DESCSTD, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, &
              WWORK, GSD%LWWORK,  RWORK, GSD%LRWORK, IWORK, GSD%LIWORK, &
              IFAIL, ICLUSTR, GAP, INFO )

         IF (M.NE.N) THEN
            WRITE (*,*) 'GSD%LWWORK',GSD%LWWORK,GSD%LRWORK,GSD%LIWORK,DESCSTD(M_)
            WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvalues found",M,N
            CALL M_exit(); stop
         ENDIF
         IF (NZ.NE.N) THEN
            WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvectors computed",M,N
            CALL M_exit(); stop
         ENDIF
         DO I2=1,N
            IF (IFAIL(I2).NE.0) THEN
               WRITE(*,*) "ERROR in subspace rotation PSSYEVX: I2,IFAIL= ",I2,IFAIL(I2)
               CALL M_exit(); stop
            ENDIF
         ENDDO
         IF (INFO<0) THEN
            WRITE(*,*) "ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= ",INFO
            CALL M_exit(); stop
         ELSE IF (INFO>0) THEN
            IF (MOD(INFO,2).NE.0) THEN
               WRITE(*,*) "ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= ",INFO
               CALL M_exit(); stop
            ELSEIF (MOD(INFO/2,2).NE.0) THEN
! this error condition happens very often, but does not imply that
! the result is not usefull
! WRITE(*,*) "WARNING eigenvector not reorthogonalized PDSYEVX/ PZHEEVX: IFAIL= ",INFO
            ELSE
               WRITE(*,*) "ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= ",INFO
               CALL M_exit(); stop
            ENDIF      
         ENDIF
        ENDIF

        IF (PRESENT(AMATIN_SCALA)) THEN
           AMATIN_SCALA=Z
        ELSE
           CALL RECON(AMATIN, NDIM, N,Z, DESCSTD)
        ENDIF

     ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      IF (.NOT.PRESENT(AMATIN_SCALA)) THEN
         DEALLOCATE(A)
      ENDIF
      DEALLOCATE(Z)
      DEALLOCATE(IWORK)
      DEALLOCATE(WWORK)
      DEALLOCATE(RWORK)

      RETURN
    END SUBROUTINE pDSSYEX_ZHEEVX


!=======================================================================
!
! call to pZHEEVX respectively pDSYEVX
! i.e. diagonalization of an symmetric or hermitian matrix
!
! a single precision matrix is passed to the routine
!
!=======================================================================

    SUBROUTINE pSSYEX_CHEEVX(COMM, AMATIN, W, NDIM, N, ABSTOL_, AMATIN_SCALA)
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               ! NxN matrix to be distributed
      INTEGER NDIM            ! leading dimension of matrix
      COMPLEX(qs)   AMATIN(NDIM,N)  ! input/output matrix
      REAL(q) W(N)            ! eigenvalues
      REAL(q), OPTIONAL ::  ABSTOL_ ! tolerance
      COMPLEX(q), OPTIONAL, TARGET ::  AMATIN_SCALA(GSD%MALLOCPQ) ! alternative input/output matrix (1 format)
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(q) :: GAP(COMM%NCPU)
      REAL(q) :: ZERO=0.0_q
      INTEGER INFO
      INTEGER M,NZ
      INTEGER I2
      REAL(q) :: ABSTOL
      COMPLEX(q), ALLOCATABLE,TARGET ::  A(:)     ! A matrix to be diagonalized (allocatable)
      COMPLEX(q), POINTER            ::  AP(:)    ! A matrix to be diagonalized (pointer) either A or AMATIN_SCALA
      COMPLEX(q), ALLOCATABLE    ::  Z(:)     ! Z eigenvector matrix
      REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
      COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array

      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE
# 987

      CALL INIT_scala(COMM, N)

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
         WRITE(*,*) "pDSSYEX_ZHEEVX: too many processors ", &
              GSD%NPROW*GSD%NPCOL,COMM%NCPU
         CALL M_exit(); stop
      ENDIF
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      IF (PRESENT(AMATIN_SCALA)) THEN
         AP=>AMATIN_SCALA
      ELSE
         ALLOCATE(A(GSD%MALLOCPQ))
         AP=>A
      ENDIF
      ALLOCATE(Z(GSD%MALLOCPQ))
      ALLOCATE(IWORK(GSD%LIWORK))
      ALLOCATE(WWORK(GSD%LWWORK))
      ALLOCATE(RWORK(GSD%LRWORK))

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

!     get parts of the matrix into the local arrays
         IF (.NOT. PRESENT(AMATIN_SCALA)) THEN
            CALL DISTRI_SINGLE(AMATIN, NDIM, N,A(1),DESCSTD)
         ENDIF
         INFO=0
         
!     call scalapack routine
# 1037


        IF (.NOT. LELPA) THEN
# 1044

         CALL PZHEEVX( 'V', 'A', 'U', N, AP(1), 1, 1, DESCSTD, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, &
                    WWORK, GSD%LWWORK,  RWORK, GSD%LRWORK, IWORK, GSD%LIWORK, &
                    IFAIL, ICLUSTR, GAP, INFO )

         IF (M.NE.N) THEN
            WRITE (*,*) 'GSD%LWWORK',GSD%LWWORK,GSD%LRWORK,GSD%LIWORK,DESCSTD(M_)
            WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvalues found",M,N
            CALL M_exit(); stop
         ENDIF
         IF (NZ.NE.N) THEN
            WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvectors computed",M,N
            CALL M_exit(); stop
         ENDIF
         DO I2=1,N
            IF (IFAIL(I2).NE.0) THEN
               WRITE(*,*) "ERROR in subspace rotation PSSYEVX: I2,IFAIL= ",I2,IFAIL(I2)
               CALL M_exit(); stop
            ENDIF
         ENDDO
         IF (INFO<0) THEN
            WRITE(*,*) "ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= ",INFO
            CALL M_exit(); stop
         ELSE IF (INFO>0) THEN

            IF (MOD(INFO,2).NE.0) THEN
               WRITE(*,*) "ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= ",INFO
               CALL M_exit(); stop
            ELSEIF (MOD(INFO/2,2).NE.0) THEN
! this error condition happens very often, but does not imply that
! the result is not usefull
! WRITE(*,*) "WARNING eigenvector not reorthogonalized PDSYEVX/ PZHEEVX: IFAIL= ",INFO
            ELSE
               WRITE(*,*) "ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= ",INFO
               CALL M_exit(); stop
            ENDIF
         ENDIF
        ENDIF

        IF (PRESENT(AMATIN_SCALA)) THEN
           AMATIN_SCALA=Z
        ELSE
           CALL RECON_SINGLE(AMATIN, NDIM, N,Z,DESCSTD)
        ENDIF
      ENDIF calc

!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      IF (.NOT.PRESENT(AMATIN_SCALA)) THEN
         DEALLOCATE(A)
      ENDIF
      DEALLOCATE(Z)
      DEALLOCATE(IWORK)
      DEALLOCATE(WWORK)
      DEALLOCATE(RWORK)
      RETURN
    END SUBROUTINE pSSYEX_CHEEVX

!=======================================================================
!
! clean the matrix elements along the diagonal (force them to be real)
! as it is 1._q in the seriel version
!
!=======================================================================


    SUBROUTINE BG_CHANGE_DIAGONALE(N, A, IU0)
      IMPLICIT NONE
      COMPLEX(q)    A(*) ! distributed matrix
      INTEGER N    ! dimension of the matrix
      INTEGER IU0  ! IO unit for error report
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'BG_CHANGE_DIAGONALE' )

      CALL BLACS_GRIDINFO(DESCSTD(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCSTD(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCSTD(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCSTD(NB_)
        J1RES=MIN(DESCSTD(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCSTD(MB_)
            I1RES=MIN(DESCSTD(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCSTD(MB_)*MYROW+NPROW*(I1-1)+I2) == (DESCSTD(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                 IF (ABS(AIMAG(A(IROW_JCOL)))>1E-2_q .AND. IU0>=0) THEN
                    WRITE(IU0,*)'WARNING: Sub-Space-Matrix is not hermitian subr', &
                    &                 AIMAG(A(IROW_JCOL)),DESCSTD(NB_)*MYCOL+NPCOL*(J1-1)+J2
                 ENDIF
                A(IROW_JCOL) = REAL( A(IROW_JCOL) ,KIND=q)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCSTD(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE BG_CHANGE_DIAGONALE


!=======================================================================
!
! call to pCHEEVX respectively PSSYEVX
! i.e. diagonalization of an symmetric or hermitian matrix
!
!=======================================================================

    SUBROUTINE pSSYEX_CHEEVX_single(COMM, AMATIN, W, NDIM, N, ABSTOL_, AMATIN_SCALA)
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               ! NxN matrix to be distributed
      INTEGER NDIM            ! leading dimension of matrix
      COMPLEX(qs)   AMATIN(NDIM,N)  ! input/output matrix
      REAL(q) W(N)            ! eigenvalues
      REAL(q), OPTIONAL ::  ABSTOL_ ! tolerance
      COMPLEX(qs), OPTIONAL, TARGET ::  AMATIN_SCALA(GSD%MALLOCPQ) ! alternative input/output matrix (1 format)
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(qs) :: GAP(COMM%NCPU)
      REAL(qs) :: ZERO=0.0_q
      INTEGER INFO
      INTEGER M,NZ
      INTEGER I2
      REAL(qs) :: ABSTOL
      REAL(qs) :: W_single(N)            ! eigenvalues
      COMPLEX(qs), ALLOCATABLE,TARGET ::  A(:)     ! A matrix to be diagonalized (allocatable)
      COMPLEX(qs), POINTER            ::  AP(:)    ! A matrix to be diagonalized (pointer) either A or AMATIN_SCALA
      COMPLEX(qs), ALLOCATABLE    ::  Z(:)     ! Z eigenvector matrix
      REAL(qs), ALLOCATABLE ::  RWORK(:) ! work array
      COMPLEX(qs), ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE

      CALL INIT_scala(COMM, N)

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
         WRITE(*,*) "pDSSYEX_ZHEEVX: too many processors ", &
              GSD%NPROW*GSD%NPCOL,COMM%NCPU
         CALL M_exit(); stop
      ENDIF
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      IF (PRESENT(AMATIN_SCALA)) THEN
         AP=>AMATIN_SCALA
      ELSE
         ALLOCATE(A(GSD%MALLOCPQ))
         AP=>A
      ENDIF
      ALLOCATE(Z(GSD%MALLOCPQ))
      ALLOCATE(IWORK(GSD%LIWORK))
      ALLOCATE(WWORK(GSD%LWWORK))
      ALLOCATE(RWORK(GSD%LRWORK))
      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

!     get parts of the matrix into the local arrays
         IF (.NOT. PRESENT(AMATIN_SCALA)) THEN
            CALL DISTRI_SINGLE_SINGLE(AMATIN, NDIM, N,A(1),DESCSTD)
         ENDIF
         INFO=0
         
!     call scalapack routine
# 1245

         CALL PCHEEVX( 'V', 'A', 'U', N, AP(1), 1, 1, DESCSTD, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W_single, ORFAC_single, Z, 1, 1, DESCSTD, &
                    WWORK, GSD%LWWORK,  RWORK, GSD%LRWORK, IWORK, GSD%LIWORK, &
                    IFAIL, ICLUSTR, GAP, INFO )

         IF (M.NE.N) THEN
            WRITE (*,*) 'GSD%LWWORK',GSD%LWWORK,GSD%LRWORK,GSD%LIWORK,DESCSTD(M_)
            WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvalues found",M,N
            CALL M_exit(); stop
         ENDIF
         IF (NZ.NE.N) THEN
            WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvectors computed",M,N
            CALL M_exit(); stop
         ENDIF
         DO I2=1,N
            IF (IFAIL(I2).NE.0) THEN
               WRITE(*,*) "ERROR in subspace rotation PSSYEVX: I2,IFAIL= ",I2,IFAIL(I2)
               CALL M_exit(); stop
            ENDIF
         ENDDO
! cluster check, gap check ?
!        IF (INFO.NE.0) THEN
!          IF (NODE_ME==IONODE) WRITE(0,*) "WARNING in subspace rotation, PSSYEVX: INFO=",INFO
!        ENDIF
!     reconstruct from the Z matrix

         IF (PRESENT(AMATIN_SCALA)) THEN
            AMATIN_SCALA=Z
         ELSE
            CALL RECON_SINGLE_SINGLE(AMATIN, NDIM, N,Z,DESCSTD)
         ENDIF
         W=W_single
      ENDIF calc

!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      IF (.NOT.PRESENT(AMATIN_SCALA)) THEN
         DEALLOCATE(A)
      ENDIF
      DEALLOCATE(Z)
      DEALLOCATE(IWORK)
      DEALLOCATE(WWORK)
      DEALLOCATE(RWORK)

      RETURN
    END SUBROUTINE pSSYEX_CHEEVX_single

!=======================================================================
!
! matrix diagonalization 1 aware version
! this version is virtually identical to pDSSYEX_ZHEEVX
! when called with the additional AMATIN_SCALA argument
! i.e.
!  BG_pDSSYEX_ZHEEVX(COMM, A, W, N) could be replaced by
!  pDSSYEX_ZHEEVX(COMM, ADUMMY, W, NDUMMY, N, AMATIN_SCALA=A)
! but lets leave it here for the time being
!
!=======================================================================

    SUBROUTINE BG_pDSSYEX_ZHEEVX(COMM, A, W, N)
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               ! NxN matrix to be distributed
      COMPLEX(q)    A(GSD%NP,GSD%NQ)  ! input/output matrix
      REAL(q) W(N)            ! eigenvalues
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(q) :: GAP(COMM%NCPU)
      REAL(q) :: ZERO=0.0_q
      INTEGER INFO
      INTEGER M,NZ
      INTEGER I2

      COMPLEX(q)    Z(GSD%NP,GSD%NQ)
      REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
      COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE
# 1333

      CALL INIT_scala(COMM, N )

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
        WRITE(*,*) "pDSSYEX_ZHEEVX: too many processors ", &
                     GSD%NPROW*GSD%NPCOL,COMM%NCPU
        CALL M_exit(); stop
      ENDIF
!-----------------------------------------------------------------------
! allocation (on T3D allocated in init_T3D)
!-----------------------------------------------------------------------
      ALLOCATE(IWORK(GSD%LIWORK))
      ALLOCATE(WWORK(GSD%LWWORK))
      ALLOCATE(RWORK(GSD%LRWORK))
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

      INFO=0

# 1365


!     call scalapack routine
      IF (.NOT. LELPA) THEN
# 1373

        CALL PZHEEVX( 'V', 'A', 'U', N, A, 1, 1, DESCSTD, ZERO, ZERO, 13,  &
                    -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, &
                    WWORK, GSD%LWWORK,  RWORK, GSD%LRWORK, IWORK, GSD%LIWORK, &
                    IFAIL, ICLUSTR, GAP, INFO )

        IF (M.NE.N) THEN
          WRITE (*,*) 'GSD%LWWORK',GSD%LWWORK,GSD%LRWORK,GSD%LIWORK,DESCSTD(M_)
          WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvalues found",M,N
          CALL M_exit(); stop
        ENDIF
        IF (NZ.NE.N) THEN
          WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvectors computed",M,N
          CALL M_exit(); stop
        ENDIF
        DO I2=1,N
          IF (IFAIL(I2).NE.0) THEN
            WRITE(*,*) "ERROR in subspace rotation PSSYEVX: I2,IFAIL= ",I2,IFAIL(I2)
            CALL M_exit(); stop
          ENDIF
        ENDDO
        IF (INFO<0) THEN
           WRITE(*,*) "ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= ",INFO
           CALL M_exit(); stop
        ELSE IF (INFO>0) THEN
           IF (MOD(INFO,2).NE.0) THEN
              WRITE(*,*) "ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= ",INFO
              CALL M_exit(); stop
           ELSEIF (MOD(INFO/2,2).NE.0) THEN
! this error condition happens very often, but does not imply that
! the result is not usefull
! WRITE(*,*) "WARNING eigenvector not reorthogonalized PDSYEVX/ PZHEEVX: IFAIL= ",INFO
           ELSE
              WRITE(*,*) "ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= ",INFO
              CALL M_exit(); stop
           ENDIF
        ENDIF
      ENDIF

      A=Z

   ENDIF calc

!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(IWORK)
      DEALLOCATE(WWORK)
      DEALLOCATE(RWORK)

      RETURN
    END SUBROUTINE BG_pDSSYEX_ZHEEVX


!=======================================================================
!
! matrix LU decompositions 1 aware version
!  Determine A = U+ U (POTRF) and invert U (TRTRI)
!
!=======================================================================

    SUBROUTINE BG_pPOTRF_TRTRI(A, N, INFO)
      IMPLICIT NONE

      INTEGER N                  ! NxN matrix to be distributed
      COMPLEX(q)    A(GSD%NP,GSD%NQ)     ! input/output matrix
      INTEGER INFO

      CALL CHECK_scala (N, 'BG_pPOTRF_TRTRI' )

# 1445

         CALL PZPOTRF &

          & ('U',N, A(1,1),1,1,DESCSTD,INFO)
# 1451

         CALL PZTRTRI &

          & ('U','N',N,A(1,1),1,1,DESCSTD,INFO)

      END SUBROUTINE BG_pPOTRF_TRTRI 

!=======================================================================
!
! matrix LU decompositions 1 aware version
!  Determine A = U+ U (POTRF) and invert A = U+U (POTRI)
!
!=======================================================================

    SUBROUTINE BG_pPOTRF_POTRI(A, N, INFO)
      IMPLICIT NONE

      INTEGER N               ! NxN matrix to be distributed
      COMPLEX(q)    A(GSD%NP,GSD%NQ)     ! input/output matrix
      INTEGER INFO

      CALL CHECK_scala (N, 'BG_pPOTRF_POTRI' )
# 1475

         CALL PZPOTRF &

          & ('U',N, A(1,1),1,1,DESCSTD,INFO)
# 1481

         CALL PZPOTRI &

          & ('U',N, A(1,1), 1, 1, DESCSTD, INFO)

     END SUBROUTINE BG_pPOTRF_POTRI



!=======================================================================
!
! the following routines are called with descriptors DESCA
! the routines can be called from any context or subgroup
! of 1 processed as long as DESCA is properly set up
! they also to not refer to GSD
! in that sense they are much safer to use, but matter
! of fact the caller must "know" scalapack
!
! for the time being the check however that the passed matrix
! dimensions are consistent with those set up in GSD
! by calling CHECK_scala
! this should be removed after more extensive testing
!
!=======================================================================

    SUBROUTINE DISTRI(AMATIN, NDIM, N, A,DESCA)
      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      COMPLEX(q)    AMATIN(NDIM,N),A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'DISTRI' )

      IF (N>NDIM) THEN
         WRITE(*,*) 'internal error in scala.F: leading dimension of matrix too small'
         CALL M_exit(); stop
      END IF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      ! "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI

! same as above using a single precision input array

    SUBROUTINE DISTRI_SINGLE(AMATIN, NDIM, N, A,DESCA)
      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      COMPLEX(qs)   AMATIN(NDIM,N)
      COMPLEX(q)    A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'DISTRI_SINGLE' )

      IF (N>NDIM) THEN
         WRITE(*,*) 'internal error in scala.F: leading dimension of matrix too small'
         CALL M_exit(); stop
      ENDIF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      ! "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI_SINGLE

! same as above using a single precision input array

    SUBROUTINE DISTRI_SINGLE_SINGLE(AMATIN, NDIM, N, A,DESCA)
      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      COMPLEX(qs)   AMATIN(NDIM,N)
      COMPLEX(qs)   A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'DISTRI_SINGLE_SINGLE' )


      IF (N>NDIM) THEN
         WRITE(*,*) 'internal error in scala.F: leading dimension of matrix too small'
         CALL M_exit(); stop
      ENDIF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      ! "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI_SINGLE_SINGLE


!=======================================================================
!
! merge distributed matrix into (1._q,0._q) large matrix
!
!=======================================================================

    SUBROUTINE RECON(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      COMPLEX(q)    AMATIN(NDIM, N),A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'RECON' )

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON

! same as above using a single precision output array

    SUBROUTINE RECON_SINGLE(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      COMPLEX(qs)   AMATIN(NDIM, N)
      COMPLEX(q)    A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'RECON_SINGLE' )

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SINGLE


! same as above using a single precision output array

    SUBROUTINE RECON_SINGLE_SINGLE(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      COMPLEX(qs)   AMATIN(NDIM, N)
      COMPLEX(qs)   A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'RECON_SINGLE_SINGLE' )

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SINGLE_SINGLE

!=======================================================================
!
! matrix diagonalization 1 aware version
! this version expects a proper 1 descriptor to
! note that this routine does stop when eigenvectors
!
!=======================================================================

     SUBROUTINE PDSSYEX_ZHEEVX_DESC(A, W, N, DESCA_, COMM, JOBIN)
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       COMPLEX(q)    A(:)            ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ ) ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER ICLUSTR(2*COMM%NCPU)
       INTEGER IFAIL(N)
       REAL(q) :: GAP(COMM%NCPU)
       REAL(q) :: ZERO=0.0_q
       INTEGER INFO
       INTEGER M,NZ
       INTEGER I2
       
       COMPLEX(q), ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A)))

# 1866

       CALL PZHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
            IFAIL, ICLUSTR, GAP, INFO )
       RWORK(1)=RWORK(1)+NCLUST*N

       LWWORK= MAX(1,INT(WWORK(1)))
       LIWORK= MAX(1,INT(IWORK(1)))
       LRWORK= MAX(1,INT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

# 1884

       CALL PZHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
            IFAIL, ICLUSTR, GAP, INFO )

       IF (M.NE.N) THEN
          WRITE (*,*) 'LWWORK',LWWORK,LRWORK,LIWORK,DESCA(M_)
          WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvalues found",M,N
          CALL M_exit(); stop
       ENDIF
       IF (NZ.NE.N) THEN
          WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvectors computed",M,N
          CALL M_exit(); stop
       ENDIF
       DO I2=1,N
          IF (IFAIL(I2).NE.0) THEN
             WRITE(*,*) "ERROR in subspace rotation PSSYEVX: I2,IFAIL= ",I2,IFAIL(I2)
             CALL M_exit(); stop
          ENDIF
       ENDDO
       IF (INFO<0) THEN
          WRITE(*,*) "ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= ",INFO
          CALL M_exit(); stop
       ELSE IF (INFO>0) THEN
          IF (MOD(INFO,2).NE.0) THEN
             WRITE(*,*) "ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= ",INFO
             CALL M_exit(); stop
          ELSEIF (MOD(INFO/2,2).NE.0) THEN
! this error condition happens very often, but does not imply that
! the result is not usefull
             WRITE(*,*) "WARNING eigenvector not reorthogonalized PDSYEVX/ PZHEEVX: IFAIL= ",INFO
          ELSE
             WRITE(*,*) "ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= ",INFO
             CALL M_exit(); stop
          ENDIF
       ENDIF
    
       A=Z
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

    END SUBROUTINE PDSSYEX_ZHEEVX_DESC



!=======================================================================
!
! matrix diagonalization using divide and conquere 1
! this version expects a proper 1 descriptor to
!
!=======================================================================

     SUBROUTINE PDSSYEV_ZHEEVD_DESC(A, W, N, DESCA_, COMM, JOBIN)
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       COMPLEX(q)    A(:)            ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ )  ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       
       COMPLEX(q), ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A)))

# 1980

       CALL PZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )

       LWWORK= MAX(1,INT(WWORK(1)))
       LIWORK= MAX(1,INT(IWORK(1)))
       LRWORK= MAX(1,INT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

# 1996

       CALL PZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )


       A=Z
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          WRITE(*,*) "ERROR in diagonalization PDSYEVD/ PZHEEVD: IFAIL= ",INFO
          CALL M_exit(); stop
       ELSE IF (INFO>0) THEN
          WRITE(*,*) "WARNING: PDSYEVD/ PZHEEVD failed using fallback PDSYEV/ PZHEEV"
          CALL PDSSYEV_ZHEEV_DESC(A, W, N, DESCA, COMM)
       ENDIF

    END SUBROUTINE PDSSYEV_ZHEEVD_DESC

!
! same version of above with two dimensional array A
!

     SUBROUTINE PDSSYEV_ZHEEVD_DESC2(A, W, N, DESCA_, COMM, JOBIN)
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       COMPLEX(q)    A(:,:)          ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ ) ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       
       COMPLEX(q), ALLOCATABLE ::     Z(:,:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,1), SIZE(A,2)))

# 2063

       CALL PZHEEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )

       LWWORK= MAX(1,INT(WWORK(1)))
       LIWORK= MAX(1,INT(IWORK(1)))
       LRWORK= MAX(1,INT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

# 2079

       CALL PZHEEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )


       A=Z
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          WRITE(*,*) "ERROR in diagonalization PDSYEVD/ PZHEEVD: IFAIL= ",INFO
          CALL M_exit(); stop
       ELSE IF (INFO>0) THEN
          WRITE(*,*) "WARNING: PDSYEVD/ PZHEEVD failed using fallback PDSYEV/ PZHEEV"
          CALL PDSSYEV_ZHEEV_DESC2(A, W, N, DESCA, COMM)
       ENDIF

    END SUBROUTINE PDSSYEV_ZHEEVD_DESC2



!=======================================================================
!
! matrix diagonalization using standard 1  diagonalization
! usually not very fast but robust
!
!=======================================================================

     SUBROUTINE PDSSYEV_ZHEEV_DESC(A, W, N, DESCA_, COMM, JOBIN)
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       COMPLEX(q)    A(:)            ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ ) ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       COMPLEX(q), ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A)))

# 2147

       CALL PZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )

       LWWORK= MAX(1,INT(WWORK(1)))
       LRWORK= MAX(1,INT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

# 2162

       CALL PZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )

       IF (INFO.NE.0) THEN
          WRITE(*,*) "ERROR in diagonalization PDSYEV/ PZHEEV: IFAIL= ",INFO
          CALL M_exit(); stop
       ENDIF
    
       A=Z
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PDSSYEV_ZHEEV_DESC


     SUBROUTINE PDSSYEV_ZHEEV_DESC2(A, W, N, DESCA_, COMM, JOBIN)
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       COMPLEX(q)    A(:,:)          ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ ) ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       COMPLEX(q), ALLOCATABLE ::     Z(:,:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,1), SIZE(A,2)))

# 2218

       CALL PZHEEV( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )

       LWWORK= MAX(1,INT(WWORK(1)))
       LRWORK= MAX(1,INT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

# 2233

       CALL PZHEEV( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )

       IF (INFO.NE.0) THEN
          WRITE(*,*) "ERROR in diagonalization PDSYEV/ PZHEEV: IFAIL= ",INFO
          CALL M_exit(); stop
       ENDIF
    
       A=Z
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PDSSYEV_ZHEEV_DESC2


!=======================================================================
!
! the routine makes an (approximately unitary) transformation matrix
! exactly unitary (or orthogonal)
!
! if the orbitals (observing <n|m>= delta_mn )
! are transformed by a matrix A, they will no longer be
! exactly orthonormal
!    |n'> = |n> A_nn'  (see LINCOM)
! but rather observe
!    <n'| m'> =  A^+_n'm A_mm' = (A+ A)_n'm'
!
! The routine calculate S=  A+ A then determines
!    S= U+ U  (U upper triangular matrix)
! inverts U and updates the transformation matrix to
!    A <-  A U^-1
! (1._q,0._q) can easily show that
!   (A U^-1)^+ (A U^-1) = 1
!
!=======================================================================

     SUBROUTINE MAKE_UNITARY_DESC(A, N, DESCA )
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       COMPLEX(q)    A(:)            ! input/output matrix
       INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
       INTEGER :: INFO 
! local variables
       COMPLEX(q), ALLOCATABLE ::     Z(:)

       ALLOCATE(Z(SIZE(A)))

! calculate S= A+ A
       CALL PZGEMM( 'C', 'N', N, N, N, (1._q,0._q),  &
          A(1), 1, 1, DESCA, &
          A(1), 1, 1, DESCA, (0._q,0._q),  &
          Z(1), 1, 1, DESCA )

! determine S= U+ U  (U upper triangular matrix)
# 2295

       CALL PZPOTRF &

          & ('U', N, Z(1), 1, 1, DESCA, INFO)

       IF (INFO/=0) THEN
          WRITE(*,*) "internal error in VASP: MAKE_UNITARY_DESC PZPOTRF failed",INFO
          CALL M_exit(); stop
       ENDIF
! invert U <- U^-1
# 2307

       CALL PZTRTRI &

          & ('U','N', N, Z(1), 1, 1, DESCA, INFO)
       IF (INFO/=0) THEN
          WRITE(*,*) "internal error in VASP: MAKE_UNITARY_DESC PZTRTRI failed",INFO
          CALL M_exit(); stop
       ENDIF
! triangular update of matrix A <- A U^-1
! note that the lower part of Z is set so DGEMM does not work here
# 2319

       CALL PZTRMM &

     &                ('R', 'U', 'N', 'N' , N, N, (1._q,0._q), &
     &                 Z, 1, 1, DESCA, A, 1, 1, DESCA)

       DEALLOCATE(Z)
       
     END SUBROUTINE MAKE_UNITARY_DESC

!=======================================================================
!
! add diagonal elements AD (usually eigenvalues) to the supplied
! distributed matrix
!
!=======================================================================


    SUBROUTINE ADD_TO_DIAGONALE(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q) A(:)        ! distributed matrix
      COMPLEX(q) AD(:)       ! diagonal elements to be added to A
      INTEGER N              ! dimension of the matrix
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == & 
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                A(IROW_JCOL) = A(IROW_JCOL)+AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE ADD_TO_DIAGONALE

!
! same as above for a real valued function AD
!

    SUBROUTINE ADD_TO_DIAGONALE_REAL(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q)    A(:)       ! distributed matrix
      REAL(q) AD(:)      ! diagonal elements to be subtracted from A
      INTEGER N          ! dimension of the matrix
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == & 
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                A(IROW_JCOL) = A(IROW_JCOL)+AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE ADD_TO_DIAGONALE_REAL


!=======================================================================
!
! determine diagonal elements AD of the supplied
! distributed matrix  A
! important note: global sum is required after the call !!
!
!=======================================================================


    SUBROUTINE DETERMINE_DIAGONALE(N, A, AD, DESCA)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q) A(:)    ! distributed matrix
      COMPLEX(q) AD(:)   ! on exit: diagonal elements
      INTEGER N          ! dimension of the matrix
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AD=0

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == & 
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=A(IROW_JCOL)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DETERMINE_DIAGONALE

!
! this version is for a COMPLEX(q) input array (REAL(q) for gammaonly)
!

    SUBROUTINE DETERMINE_DIAGONALE_GDEF(N, A, AD, DESCA)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q)    A(:)       ! distributed matrix
      COMPLEX(q) AD(:)   ! on exit: diagonal elements
      INTEGER N          ! dimension of the matrix
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AD=0

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == & 
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=A(IROW_JCOL)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DETERMINE_DIAGONALE_GDEF

!=======================================================================
!
! divide distribute matrix by eigenvalue difference
! this is usually required in standard perturbation theory
! the operation is only 1._q in the occupied-unoccupied part
! of the matrix
! the first version is used to construct the density matrix
! from a given perturbation
! the second version is used to construct an approximately
! unitary rotation matrix
!
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF(N, A, CELTOT, FERTOT, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q)    A(:)          ! distributed matrix
      COMPLEX(q) CELTOT(:)  ! eigenvalues
      REAL(q)    FERTOT(:)  ! occupancies
      INTEGER N             ! dimension of the matrix
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      REAL(q), PARAMETER  :: DIFMAX=0.001_q
      REAL(q) :: DIFCELL
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW_JCOL = IROW_JCOL+1
! identical occupancies, clear rotation matrix
               DIFCELL=(FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)- & 
                        FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
               IF (ABS(DIFCELL)<DIFMAX) THEN
                  A(IROW_JCOL) = 0
               ELSE
! determine inverse rotation matrix
                  DIFCELL= (REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)-   & 
                       CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)))/ DIFCELL
                  
! huge rotation element, should never occur since the
! occupancies should cancel those
                  IF (ABS(DIFCELL)<DIFMAX) THEN
                     A(IROW_JCOL) = 0
                  ELSE
                     A(IROW_JCOL) = A(IROW_JCOL)/DIFCELL
                     IF (ABS(A(IROW_JCOL))>0.4_q) THEN
! damp rotation matrix as in Loewdin case
                        A(IROW_JCOL)=A(IROW_JCOL)/ABS(ABS(A(IROW_JCOL)))*0.4_q
                     ENDIF
                  ENDIF
               END IF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF

!=======================================================================
!
! divide distribute matrix by eigenvalue difference
! this is usually required in standard perturbation theory
! the operation is only 1._q in the occupied-unoccupied part
! of the matrix
!
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF_SIGNED(N, A, CELTOT, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q)    A(:)          ! distributed matrix
      COMPLEX(q) CELTOT(:)  ! eigenvalues
      INTEGER N             ! dimension of the matrix
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      REAL(q), PARAMETER  :: DIFMAX=0.001_q
      REAL(q) :: DIFCELL
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              DIFCELL=REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)- & 
                           CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              IF (ABS(DIFCELL)<DIFMAX) THEN
                 A(IROW_JCOL) = 0
              ELSE
! the equations above also involve
! the difference between the occupancies, for the density matrix they would drop out
! in other words: rotations in the (un)occupied-(un)occupied  block are irrelevant
! for energies
! tests indicate that inclusion of DIFCELL makes indeed little difference for energies
                 A(IROW_JCOL) = A(IROW_JCOL)/DIFCELL
                 IF (ABS(A(IROW_JCOL))>0.4_q) THEN
! damp rotation matrix as in Loewdin case
                    A(IROW_JCOL)=A(IROW_JCOL)/ABS(ABS(A(IROW_JCOL)))*0.4_q
                 ENDIF
              ENDIF
! set diagonal to 1
              IF (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2==DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) THEN
                 A(IROW_JCOL)=1
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE
    
!=======================================================================
!
! calculate |H_ai|^2 / (epsilon_a-epsilon_i) *(f_a-f_i)
! both the matrix H_ai and the eigenvalues, and the occupancies are passed
! down
! the calculation is only performed for |epsilon_a - epsilon_i| >DIFFMAX
!
!=======================================================================

    SUBROUTINE SUM_SECOND_ORDER(N, A, CELTOT, FERTOT, CSUM, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q)    A(:)          ! distributed matrix
      REAL(q) CELTOT(:)     ! eigenvalues epsilon
      REAL(q) FERTOT(:)     ! occupancies
      COMPLEX(q) CSUM       ! final result |H_ai|^2 / (epsilon_a-epsilon_i)
      INTEGER N             ! dimension of the matrix
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      REAL(q), PARAMETER  :: DIFMAX=0.001_q
      REAL(q) :: DIFCELL
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CSUM=0
      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              DIFCELL=REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)- & 
                          CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              IF (ABS(DIFCELL)>=DIFMAX) THEN
                 CSUM=CSUM+A(IROW_JCOL)*CONJG(A(IROW_JCOL))/DIFCELL* & 
                 (FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)-FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SUM_SECOND_ORDER


!=======================================================================
!
! clear occupied-occupied and unoccupied-unoccupied block of
! a square matrix
!
!=======================================================================

    SUBROUTINE CLEAR_OCCOCC_UNOCCUNOCC(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q)    A(:)          ! distributed matrix
      INTEGER N             ! dimension of the matrix
      INTEGER LAST_OCCUPIED ! last band treated as occupied
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) >LAST_OCCUPIED .AND.  & 
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) .OR. &
                  ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<=LAST_OCCUPIED .AND. & 
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)>  LAST_OCCUPIED)) THEN
! do nothing
               ELSE
! remove the element from the perturbation
                  A(IROW_JCOL)=0
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCOCC_UNOCCUNOCC

!=======================================================================
!
! clear occupied-unoccupied and unoccupied-occupied blocks of
! a square matrix
!
!=======================================================================

    SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q)    A(:)          ! distributed matrix
      INTEGER N             ! dimension of the matrix
      INTEGER LAST_OCCUPIED ! last band treated as occupied
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) >LAST_OCCUPIED .AND.  & 
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) .OR. &
                  ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<=LAST_OCCUPIED .AND. & 
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)>  LAST_OCCUPIED)) THEN
! remove the element from the perturbation
                  A(IROW_JCOL)=0
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC
 
!=======================================================================
!
! set occupied-occupied block of a square matrix to identity matrix
!
!=======================================================================

    SUBROUTINE SET_OCCOCC_IDENTITY(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q)    A(:)          ! distributed matrix
      INTEGER N             ! dimension of the matrix
      INTEGER LAST_OCCUPIED ! last band treated as occupied
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<= LAST_OCCUPIED .AND.  & 
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) THEN
! remove the element from the perturbation
                  A(IROW_JCOL)=0._q
                  IF ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) .EQ. (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)) THEN
                     A(IROW_JCOL)=1._q
                  ENDIF
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SET_OCCOCC_IDENTITY

!=======================================================================
!
! conjugate the upper triangle of a matrix
! this routine constructs from the change of the density matrix
! the unitary transformation for the orbitals
! (1._q,0._q) can also simply take the lower triangle and change the sign
! in the upper
!
!=======================================================================

    SUBROUTINE MINUS_UPPER(N, A, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q)    A(:)          ! distributed matrix
      INTEGER N             ! dimension of the matrix
      INTEGER LAST_OCCUPIED ! last band treated as occupied
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) < &
                  (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                   A(IROW_JCOL)=-(A(IROW_JCOL))
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MINUS_UPPER


!=======================================================================
!
! multiply distribute matrix by eigenvalue difference
!   1/(iw - epsilon_i) and  1/(iw - epsilon_a)
! from the left and right
!
!=======================================================================


    SUBROUTINE MULTIPLY_BY_IVEIGENVALUE(N, A, CELTOTINV, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
      COMPLEX(q) A(:)         ! distributed matrix
      COMPLEX(q) CELTOTINV(:) ! multiplication factors e.g. 1/(iw-eigenvalue)
      INTEGER N               ! dimension of the matrix
      INTEGER LAST_OCCUPIED   ! last band treated as occupied
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*CELTOTINV(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) & 
                                           *CELTOTINV(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MULTIPLY_BY_IVEIGENVALUE

!
! use only real part of transformation
!
    SUBROUTINE MULTIPLY_BY_IVEIGENVALUE_REAL(N, A, CELTOTINV, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
      COMPLEX(q) A(:)         ! distributed matrix
      COMPLEX(q) CELTOTINV(:) ! multiplication factors e.g. 1/(iw-eigenvalue)
      INTEGER N               ! dimension of the matrix
      INTEGER LAST_OCCUPIED   ! last band treated as occupied
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*REAL(CELTOTINV(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) & 
                                           *CELTOTINV(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),q)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MULTIPLY_BY_IVEIGENVALUE_REAL


!=======================================================================
!
! right multiply matrix by a vector (most likely a BLACS routines
! exists to do this
!  A <- A R    or  A_ij = A_ij R_j
!
!=======================================================================
!
    SUBROUTINE RIGHT_MULTIPLY_BY_R(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
      COMPLEX(q) A(:)               ! distributed matrix
      REAL(q) R(:)            ! multiplication factors R
      INTEGER N               ! dimension of the matrix
      INTEGER LAST_OCCUPIED   ! last band treated as occupied
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RIGHT_MULTIPLY_BY_R


    SUBROUTINE RIGHT_MULTIPLY_BY_C(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
      COMPLEX(q) A(:)         ! distributed matrix
      COMPLEX(q) R(:)         ! multiplication factors R
      INTEGER N               ! dimension of the matrix
      INTEGER LAST_OCCUPIED   ! last band treated as occupied
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RIGHT_MULTIPLY_BY_C

!=======================================================================
!
! distribute slice of the global matrix AMATIN to a distributed matrix A
! as input the routine expects several columns of the matrix
! to be distributed
! the columns of the input matrix are COLUMN_LOW up to COLUMN_HIGH
! all rows are passed to the routine
!
!=======================================================================

    SUBROUTINE DISTRI_SLICE(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   ! first dimension of matrix
      INTEGER N                                      ! rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                ! starting and end row supplied in AMATIN
      COMPLEX(q)    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  ! input matrix (globally known)
      COMPLEX(q)    A(*)                                   ! on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         ! distributed matrix descriptor array
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      ! "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     WRITE(*,*) 'distri',NDIM, N

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN

             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   A(IROW+(JCOL-1)*DESCA(LLD_))= &
                   AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO
      ! "DISTRI",ITEST,MYCOL,MYROW
      RETURN
    END SUBROUTINE DISTRI_SLICE

!
! distribute matrix and add hermitian conjugated elements
!
!
    SUBROUTINE DISTRI_SLICE_HERM(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   ! first dimension of matrix
      INTEGER N                                      ! rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                ! starting and end row supplied in AMATIN
      COMPLEX(q)    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  ! input matrix (globally known)
      COMPLEX(q)    A(*)                                   ! on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         ! distributed matrix descriptor array
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index
      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      ! "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     set conjugated part first
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                       A(IROW+(JCOL-1)*DESCA(LLD_))=CONJG(AMATIN(COL, ROW-COLUMN_LOW+1))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     non conjugated part
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                     A(IROW+(JCOL-1)*DESCA(LLD_))=AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DISTRI_SLICE_HERM

!
! distribute matrix and add hermitian conjugated elements
! in this version it is assumed that the lower triangle of the input matrix is set
! this requires some extra precaution that previously properly set elements
! are not overwritten
! only the lower part of the input matrix is used as source
!
    SUBROUTINE DISTRI_SLICE_HERM_L(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   ! first dimension of matrix
      INTEGER N                                      ! rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                ! starting and end row supplied in AMATIN
      COMPLEX(q)    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  ! input matrix (globally known)
      COMPLEX(q)    A(*)                                   ! on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         ! distributed matrix descriptor array
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index
      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      ! "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     set conjugated part
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                     IF (COL>=ROW) A(IROW+(JCOL-1)*DESCA(LLD_))=CONJG(AMATIN(COL, ROW-COLUMN_LOW+1))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     non conjugated part
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   IF (ROW>=COL) A(IROW+(JCOL-1)*DESCA(LLD_))=AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DISTRI_SLICE_HERM_L

!=======================================================================
!
! reconstruct a slice of the originally distributed matrix AMATIN
!
!=======================================================================

    SUBROUTINE RECON_SLICE(AMATIN, NDIM, N, A, DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_ )             ! distributed matrix descriptor
      INTEGER COLUMN_LOW, COLUMN_HIGH    ! first and last column of the matrix to be merged
      INTEGER NDIM                       ! first dimension of matrix AMATIN
      INTEGER N                          ! rank of matrix (not used except by CHECK_scala)
      COMPLEX(q)    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
! on return: matrix stored in serial layout between column
!   COLUMN_LOW and COLUMN_HIGH
      COMPLEX(q)    A(*)                       ! distributed matrix
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   AMATIN(ROW,COL-COLUMN_LOW+1)= &
                     A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SLICE

!
! complex version of the previous routine
!

    SUBROUTINE RECON_SLICE_C(AMATIN, NDIM, N, A, DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_ )             ! distributed matrix descriptor
      INTEGER COLUMN_LOW, COLUMN_HIGH    ! first and last column of the matrix to be merged
      INTEGER NDIM                       ! first dimension of matrix AMATIN
      INTEGER N                          ! rank of matrix (not used except by CHECK_scala)
      COMPLEX(q) AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
! on return: matrix stored in serial layout between column
!   COLUMN_LOW and COLUMN_HIGH
      COMPLEX(q) A(*)                    ! distributed matrix
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   AMATIN(ROW,COL-COLUMN_LOW+1)= &
                     A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SLICE_C

!
! similar version as above but add complex conjugated elements
! in the lower triangle  in AMATIN
!
    SUBROUTINE RECON_SLICE_HERM(AMATIN,NDIM,N,A,DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_)              ! distributed matrix descriptor array
      INTEGER COLUMN_LOW, COLUMN_HIGH    ! first and last column of the matrix to be merged
      INTEGER NDIM                       ! first dimension of matrix AMATIN
      INTEGER N                          ! rank of matrix (not used except by CHECK_scala)
      COMPLEX(q)    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
! on return: matrix stored in serial layout between column
      COMPLEX(q)    A(*)                       ! distributed matrix
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'RECON_SLICE_HERM' )

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     set conjugated part first restrict to lower triangle
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                     IF (ROW<COL) & 
                     AMATIN(COL, ROW-COLUMN_LOW+1)= & 
                        CONJG(A(IROW+(JCOL-1)*DESCA(LLD_)))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     now the non conjugated part but restrict to upper triangle
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   IF (ROW<=COL) & 
                      AMATIN(ROW,COL-COLUMN_LOW+1)= &
                            A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SLICE_HERM
END MODULE




   
!========================== ORTH1_DISTRI ===============================
!
! ORTH1_DISTRI is essentially identical to ORTH1
! it calculates a slice or stripe of the matrix
!   O(I,J) =  <CPTWFP(I) | 1+ Q | CFW(J) >
!
! and THEN DISTRIBUTES THE SLICE in a 1 compatible manner
! using DISTRI_SLICE
! ORTH1 requires the storage of a global NB_TOT time NB_TOT matrix
! whereas ORTH1_DISTRI only the storage of the distributed matrix
!
! arguments are identical to ORTH1, except for additional communicator
!
!=======================================================================

     SUBROUTINE ORTH1_DISTRI(CSEL,CPTWFP,CFW,CPROJ,CPROW, & 
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK)
        USE prec
        USE scala
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        COMPLEX(q)      CPROJ(NPROD,NBANDS)
        COMPLEX(q)      CPROW(NPROD,NSTRIP)
        COMPLEX(q)      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI' )

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        CALL M_sum_z(MY_COMM,COVL(1,1),NBANDS*NSTRIP)

    
# 3660

        IF (.NOT. LELPA) THEN
           CALL DISTRI_SLICE(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF
      END SUBROUTINE ORTH1_DISTRI

!
! same as above, this version always adds the Hermitian conjugated elements
! to the matrix
! costs little extra time, and is more secure, in the sense that
! we can use the constructed matrix in places where the lower triangle is required as well
!
      
      SUBROUTINE ORTH1_DISTRI_HERM(CSEL,CPTWFP,CFW,CPROJ,CPROW, & 
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK)
        USE prec
        USE scala
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        COMPLEX(q)      CPROJ(NPROD,NBANDS)
        COMPLEX(q)      CPROW(NPROD,NSTRIP)
        COMPLEX(q)      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI_HERM' )

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        CALL M_sum_z(MY_COMM,COVL(1,1),NBANDS*NSTRIP)

        IF (CSEL=='U') THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ELSE
           CALL DISTRI_SLICE_HERM_L(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF
      END SUBROUTINE ORTH1_DISTRI_HERM


!========================== ORTH1_DISTRI_HERM_DESC =====================
!
! same as above, this version requires an additional argument, the descriptor DESCA
! of the matrix
!
!=======================================================================
      
      SUBROUTINE ORTH1_DISTRI_HERM_DESC(CSEL,CPTWFP,CFW,CPROJ,CPROW, & 
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK, DESCA)
        USE prec
        USE scala
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        COMPLEX(q)      CPROJ(NPROD,NBANDS)
        COMPLEX(q)      CPROW(NPROD,NSTRIP)
        COMPLEX(q)      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL
        INTEGER DESCA(*)

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        CALL M_sum_z(MY_COMM,COVL(1,1),NBANDS*NSTRIP)

        IF (CSEL=='U') THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCA,NPOS,NPOS+NSTRIP-1)
        ELSE
           CALL DISTRI_SLICE_HERM_L(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCA,NPOS,NPOS+NSTRIP-1)
        ENDIF
      END SUBROUTINE ORTH1_DISTRI_HERM_DESC


!========================== ORTH1_DISTRI_DAVIDSON ======================
!
! ORTH1_DISTRI_DAVIDSON is identical to ORTH1_DISTRI
! but also changes the intermediate COVL matrix as required in
! the davidson.F routine (requires vector R)
!
!=======================================================================

      SUBROUTINE ORTH1_DISTRI_DAVIDSON(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,R,MY_COMM, NBANDSK )
        USE prec
        USE scala

        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK, I
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        COMPLEX(q)      CPROJ(NPROD,NBANDS)
        COMPLEX(q)      CPROW(NPROD,NSTRIP)
        COMPLEX(q)      COVL(NBANDS,NSTRIP), COVL2(*)
        REAL(q) :: R(NSTRIP)
        CHARACTER*(*) CSEL

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI_DAVIDSON' )

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        COVL(NPOS:NPOS+NSTRIP-1,1:NSTRIP)=0
        IF (MY_COMM%NODE_ME==MY_COMM%IONODE) THEN
           DO I=1,NSTRIP
              COVL(NPOS-1+I,I)=R(I)
           ENDDO
        ENDIF

        CALL M_sum_z(MY_COMM,COVL(1,1),NBANDS*NSTRIP)

# 3785

        IF (.NOT. LELPA) THEN
           CALL DISTRI_SLICE(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF

      END SUBROUTINE ORTH1_DISTRI_DAVIDSON

!========================== LINCOM_DISTRI ==============================
!
! LINCOM_DISTRI does essentially the same as LINCOM but
! accepts a distributed matrix COVL2 instead of a global
! NB_TOT times NB_TOT matrix
! the required matrix COVL is reconstructed on the fly in memory
! slice by slice
! the disadvantage is that the transformed wavefunctions need
! to be stored in a temporary array
!
!=======================================================================

      SUBROUTINE LINCOM_DISTRI(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        COMPLEX(q)        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        COMPLEX(q)        COVL(*)
! local
        COMPLEX(q)        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

        CALL CHECK_scala (NIN, 'LINCOM_DISTRI' )

! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK
          CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)

          CALL M_sum_z(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT)
          
          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)
        
        CALL M_sum_z(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT)

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

        CF   =CF_RESULT
        CPROF=CPROF_RESULT

      END SUBROUTINE LINCOM_DISTRI

!
! similar version as above but add the Hermitian conjugated elements
! in the lower triangle
! can be used for "full" mode (MODE="F") if the transformation
! is Hermitian and only upper triangle is set (ZPOTRI does that)
!
      SUBROUTINE LINCOM_DISTRI_HERM(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        COMPLEX(q)        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        COMPLEX(q)        COVL(*)
! local
        COMPLEX(q)        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

        CALL CHECK_scala (NIN, 'LINCOM_DISTRI_HERM' )

! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK
          CALL RECON_SLICE_HERM(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)

          CALL M_sum_z(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT)

          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE_HERM(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)
      
        CALL M_sum_z(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT)

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

        CF   =CF_RESULT
        CPROF=CPROF_RESULT

      END SUBROUTINE LINCOM_DISTRI_HERM


!========================== LINCOM_DISTRI_DESC =========================
!
! LINCOM_DISTRI_DESC essentially the same as LINCOM_DISTRI
! but additionally requires the distributed matrix descriptor
!
!=======================================================================

      SUBROUTINE LINCOM_DISTRI_DESC(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK, DESCA)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        COMPLEX(q)        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        COMPLEX(q)        COVL(*)
        INTEGER DESCA(*)                               ! distributed matrix descriptor array
! local
        COMPLEX(q)        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK
          CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCA  ,NPOS,NPOS+NOUT-1)

          CALL M_sum_z(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT)
          
          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCA  ,NPOS,NPOS+NOUT-1)
        
        CALL M_sum_z(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT)

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

        CF   =CF_RESULT
        CPROF=CPROF_RESULT

      END SUBROUTINE LINCOM_DISTRI_DESC

!=======================================================================
!
! routine to call DISTRI_SLICE without explicitly checking
! the interface
!
!=======================================================================
      SUBROUTINE DISTRI_SLICE_NOINT(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
        USE prec
        USE scala
        IMPLICIT NONE
      
        INTEGER NDIM                                   ! first dimension of matrix
        INTEGER N                                      ! rank of matrix
        INTEGER COLUMN_LOW, COLUMN_HIGH                ! starting and end row supplied in AMATIN
        COMPLEX(q)    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  ! input matrix (globally known)
        COMPLEX(q)    A(*)                                   ! output distributed matrix
        INTEGER DESCA(*)                               ! distributed matrix descriptor array

        CALL DISTRI_SLICE(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)

      END SUBROUTINE DISTRI_SLICE_NOINT


      SUBROUTINE BG_pDSSYEX_ZHEEVX_NOINT(COMM,A,W,N)
      USE prec
      USE scala
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               ! NxN matrix to be distributed
      COMPLEX(q)    A(*)            ! input/output matrix
      REAL(q) W(N)            ! eigenvalues

      CALL BG_pDSSYEX_ZHEEVX(COMM,A(1),W,N)

      END SUBROUTINE BG_pDSSYEX_ZHEEVX_NOINT

!=======================================================================
!
! dummy routines if 1 is not compiled in
!
!=======================================================================
# 4084





