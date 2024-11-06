# 1 "dfast.F"
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

# 2 "dfast.F" 2 
MODULE dfast
  USE prec
  USE mpimy
!***********************************************************************
!
! this module contains a few high level BLAS 3 routines
! to calculate inproducts between sets of wave functions
! and to perform transformation (sub space rotations)
!
!***********************************************************************

!
! the overlap or Hamiltonmatrix between states is calculated blockwise
! if the upper of lower triangle is required
! NSTRIPD is the maximum value ever used, whereas NSTRIP_STANDARD
! is the typical value used
  INTEGER, PARAMETER :: NSTRIPD=16
  INTEGER, SAVE :: NSTRIP_STANDARD
  INTEGER, SAVE :: NSTRIP_STANDARD_GLOBAL
!
! when the wavefunctions are rotated (unitary transformations)
! this is 1._q in blocks to save storage for the transformed wavefunctions
  INTEGER :: NBLK=256

  INTERFACE
     SUBROUTINE ORTH1(CSEL,CPTWFP,CFW,GPROJ,CPROW,NBANDS, &
          &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
       USE prec
       IMPLICIT COMPLEX(q) (C)
       IMPLICIT REAL(q) (A-B,D-H,O-Z)
       
       REAL(q)      GPROJ
       REAL(q)      CPROW
       REAL(q)      COVL
       CHARACTER (LEN=*) CSEL
     END SUBROUTINE ORTH1
  END INTERFACE
  
  INTERFACE
     SUBROUTINE ORTH2(CPTWFP,CFW,GPROJ,CPROW,NBANDS, &
          &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
       USE prec
       IMPLICIT COMPLEX(q) (C)
       IMPLICIT REAL(q) (A-B,D-H,O-Z)

       REAL(q)      GPROJ
       REAL(q)      CPROW
       REAL(q)      COVL
     END SUBROUTINE ORTH2
  END INTERFACE

  INTERFACE
     SUBROUTINE ORTH3(CPTWFP,CFW,GPROJ,CPROW, &
          & NBANDS,NSTRIP,NDIM,NPOS, &
          & NPL,NPRO,NPLDIM,NPROD,COVL)
       USE prec
       IMPLICIT COMPLEX(q) (C)
       IMPLICIT REAL(q) (A-B,D-H,O-Z)

       REAL(q)      GPROJ
       REAL(q)      CPROW
       REAL(q)      COVL
     END SUBROUTINE ORTH3
  END INTERFACE


  TYPE parallel_gemm
     INTEGER :: N          ! dimension of matrix
     TYPE (communic) :: COMM
     INTEGER :: NPROC      ! number of cores
     INTEGER, ALLOCATABLE :: NCOL(:), OFFSET(:), NCTOT(:), OFFDATA(:)
  END TYPE

  CONTAINS
!***********************************************************************
!
! set default for the blocking
!
!***********************************************************************

    SUBROUTINE SET_NBLK_NSTRIP(WDES)
      USE wave
      TYPE (wavedes)  WDES
      
      INTEGER NCPU


      NCPU   =WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
# 92

! use typically 1._q quarter of plane wave coefficients but never more than 256
# 96

      IF (NBLK==-1) NBLK=MIN(256,MAX(32,(WDES%NRPLWV/256)*64))

      CALL M_max_i(WDES%COMM, NBLK, 1)
! set NSTRIP between [1 and NSTRIPD]
      NSTRIP_STANDARD=MAX(MIN(NSTRIPD,32/NCPU,WDES%NBANDS),1)
      NSTRIP_STANDARD_GLOBAL=NSTRIP_STANDARD*NCPU

    END SUBROUTINE SET_NBLK_NSTRIP


!************************ SUBROUTINE LINCOM ****************************
!
! build linear combinations of wavefunctions according to matrix CTRANS
! this subroutine performes implicitly a MATRIX x MATRIX multiplication,
! it is needed for the unitary transformation of the wavefunctions or
! for orthogonalisation routines and uses a blocked algorithm
! to save storage
!
!  COUT_n,k =   sum_kp CIN_n,kp CH_ kp,k
! where n  =  1 ... NPL    (leading dimension of array = NPLDIM)
!       k  =  1 ... NOUT
!       kp =  1 ... NIN    (NIN must be greater or equal NOUT)
! Important: on exit the input array will be overwritten by output data!
!
!  MODE determines the mode for the transformation
!  "U"   CTRANS   upper triangle set
!  "L"   CTRANS   lower triangle set
!  "F"   CTRANS   all components set
!  "C"   CTRANS   contains the conjugated transformation matrix
!  "T"   CTRANS   contains the transposed transformation matrix
!
!  "A"   used only for Davidson
!  "B"   used only for Davidson
!
! LINCOM_BLAS has a BLAS like calling sequence and does not
! allow for the Davidson like functionality
!
!***********************************************************************


    SUBROUTINE LINCOM(MODE,CF,CPROF,CTRANS,NIN,NOUT,NPL, &
       &           NPRO,NPLDIM,NPROD,LDTRAN,CFA,CPROFA)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      CHARACTER (1) MODE
      DIMENSION   CF(NPLDIM,NIN),CFA(NPLDIM,NIN)
      REAL(q)        CPROF(NPROD,NIN),CPROFA(NPROD,NIN)
      REAL(q)        CTRANS(LDTRAN,NIN)

! work array
      COMPLEX(q),ALLOCATABLE ::   CBLOCK(:,:)
      ALLOCATE(CBLOCK(NBLK,LDTRAN))

      CALL LINBAS(MODE,CF,CBLOCK,CTRANS,NIN,NOUT,2* NPL, &
     &            2* NPLDIM,LDTRAN,2* NBLK,CFA)
      IF (NPRO/=0) THEN
      CALL LINBAS(MODE,CPROF,CBLOCK,CTRANS,NIN,NOUT, NPRO, &
     &             NPROD,LDTRAN,2* NBLK,CPROFA)
      ENDIF
      DEALLOCATE(CBLOCK)

      RETURN
    END SUBROUTINE LINCOM


    SUBROUTINE LINCOM_BLAS(MODE, NPL, NPRO, NIN, NOUT, CF, NPLDIM, CPROF, NPROD, CTRANS, LDTRAN)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      CHARACTER (1) MODE
      DIMENSION   CF(NPLDIM,NIN)
      REAL(q)        CPROF(NPROD,NIN)
      REAL(q)        CTRANS(LDTRAN,NIN)

! work array
      COMPLEX(q),ALLOCATABLE ::   CBLOCK(:,:)
      ALLOCATE(CBLOCK(NBLK,LDTRAN))

      CALL LINBAS_BLAS(MODE, 2* NPL, NIN, NOUT, CF, 2* NPLDIM, CTRANS, LDTRAN, CBLOCK, 2* NBLK)
      IF (NPRO/=0) THEN
      CALL LINBAS_BLAS(MODE, NPRO,   NIN, NOUT, CPROF, NPROD,  CTRANS, LDTRAN, CBLOCK, 2* NBLK)
      ENDIF
      DEALLOCATE(CBLOCK)

      RETURN
    END SUBROUTINE LINCOM_BLAS

!************************ SUBROUTINE SETUP_PARALLEL_GEMM ***************
!
! this subroutine sets up a simple communicator structure
! to perform parallel DGEMM calls
! the source matrices need to be known on all nodes
! the final "destination" matrix is essembled using a
! allgatherv
!
!***********************************************************************

    SUBROUTINE SETUP_PARALLEL_GEMM(COMM, NB_TOT, PGEMM_HANDLE)
      USE mpimy
      IMPLICIT NONE
      TYPE(communic) :: COMM          ! used communicator
      INTEGER        :: NB_TOT        ! matrix dimension
      TYPE (parallel_gemm), POINTER ::  PGEMM_HANDLE
! local
      INTEGER        :: I


      ALLOCATE(PGEMM_HANDLE)

      PGEMM_HANDLE%NPROC =COMM%NCPU
      PGEMM_HANDLE%N     =NB_TOT
      PGEMM_HANDLE%COMM  =COMM

      ALLOCATE( PGEMM_HANDLE%NCOL(PGEMM_HANDLE%NPROC),  PGEMM_HANDLE%OFFSET(PGEMM_HANDLE%NPROC), & 
                PGEMM_HANDLE%NCTOT(PGEMM_HANDLE%NPROC), PGEMM_HANDLE%OFFDATA(PGEMM_HANDLE%NPROC) )

      CALL DISTRIBUTE_COLUMNS(NB_TOT, PGEMM_HANDLE%NPROC, PGEMM_HANDLE%NCOL)
      PGEMM_HANDLE%OFFSET(1) = 0
      PGEMM_HANDLE%OFFDATA(1) = 0
      PGEMM_HANDLE%NCTOT(1) = PGEMM_HANDLE%NCOL(1) * NB_TOT
      DO I=2,PGEMM_HANDLE%NPROC
         PGEMM_HANDLE%OFFSET(I)  =PGEMM_HANDLE%OFFSET(I-1) + PGEMM_HANDLE%NCOL(I-1)
         PGEMM_HANDLE%NCTOT(I)   = NB_TOT * PGEMM_HANDLE%NCOL(I)
         PGEMM_HANDLE%OFFDATA(I) = NB_TOT * PGEMM_HANDLE%OFFSET(I)
      ENDDO


      CONTAINS

      SUBROUTINE  DISTRIBUTE_COLUMNS( N, NPROC, NCOL )
!
! calculate number of columns per core
!

        INTEGER :: N_PER_PROC, REMAINDER, ME1, NPROC, N, NCOL(NPROC)
        INTEGER :: NCOL_MIN=16
        INTEGER :: NPROC_MAX

        NCOL = 0

! only NPROC_MAX cores take part in the communication
! otherwise communication gets too expensive
        NPROC_MAX=MIN(MAX(1,N/NCOL_MIN),NPROC)
        
        N_PER_PROC = N/NPROC_MAX
        REMAINDER   = N - NPROC_MAX * N_PER_PROC
        DO ME1 = 1,NPROC_MAX
           IF( ME1 <= REMAINDER )THEN
              NCOL(ME1) = N_PER_PROC + 1
           ELSE
              NCOL(ME1) = N_PER_PROC
           ENDIF
        ENDDO
        
        RETURN
      END SUBROUTINE DISTRIBUTE_COLUMNS
      
     
    END SUBROUTINE SETUP_PARALLEL_GEMM


    SUBROUTINE RELEASE_PARALLEL_GEMM(PGEMM_HANDLE)
      IMPLICIT NONE

      TYPE (parallel_gemm), POINTER ::  PGEMM_HANDLE

      IF (ASSOCIATED(PGEMM_HANDLE)) THEN
         DEALLOCATE( PGEMM_HANDLE%NCOL,  PGEMM_HANDLE%OFFSET, & 
                     PGEMM_HANDLE%NCTOT, PGEMM_HANDLE%OFFDATA )
         DEALLOCATE(PGEMM_HANDLE)
      ENDIF

    END SUBROUTINE RELEASE_PARALLEL_GEMM

!************************ SUBROUTINE LINBAS ****************************
!
! the call has the same interface as DGEMM except
! for the handle which needs to be passed along as well
! note that this version works in principle only for square matrices
! and is not intended for general purpose use
! I think also the leading dimensions must be correct and identical
!
!***********************************************************************

    SUBROUTINE PARALLEL_GGEMM( PGEMM_HANDLE, TRANSA , TRANSB, N1, N2, N3, ALPHA, A, LDA, &
         B, LDB, BETA, C, LDC)
      IMPLICIT NONE

      CHARACTER(LEN=1), INTENT(IN) :: TRANSA
      CHARACTER(LEN=1), INTENT(IN) :: TRANSB
      INTEGER N1, N2, N3
      INTEGER LDA, LDB, LDC
      REAL(q)             :: ALPHA, BETA
      REAL(q)             :: A(LDA,*),B(LDB,*),C(LDC,*)
      TYPE (parallel_gemm), POINTER ::  PGEMM_HANDLE
! local
      INTEGER MPIDATA
      INTEGER INFO

      REAL(q) :: CTMP(LDC,LDC)


      MPIDATA = MPI_double_precision
# 307


      IF (.NOT. ASSOCIATED(PGEMM_HANDLE))  THEN
         WRITE(*,*) 'internal error in PARALLEL_GGEMM: PGEMM_HANDLE not set up'
         CALL M_exit(); stop
      ENDIF

      IF (N1 /= N2 .OR. N2 /= N3 .OR. LDC /= N1) THEN
         WRITE(*,*) 'internal error in PARALLEL_GGEMM: not a square matrix',N1,N2,N3,LDC
         CALL M_exit(); stop
      ENDIF

      IF (N1 /= PGEMM_HANDLE%N) THEN
         WRITE(*,*) 'internal error in PARALLEL_GGEMM: not correctly set up',N1 /= PGEMM_HANDLE%N
         CALL M_exit(); stop
      ENDIF

      IF (TRANSB=='N' .OR. TRANSB=='n') THEN
         CALL DGEMM( TRANSA, TRANSB, N1,  PGEMM_HANDLE%NCOL( PGEMM_HANDLE%COMM%NODE_ME), N3, ALPHA, A, LDA, &
              B(1,1+PGEMM_HANDLE%OFFSET( PGEMM_HANDLE%COMM%NODE_ME)), LDB, BETA, C(1,1+PGEMM_HANDLE%OFFSET( PGEMM_HANDLE%COMM%NODE_ME)), LDC )

         CTMP(1:LDC,1:LDC)=C(1:LDC,1:LDC)
         CALL MPI_allgatherv (CTMP(1,1+PGEMM_HANDLE%OFFSET(PGEMM_HANDLE%COMM%NODE_ME)), PGEMM_HANDLE%NCTOT(PGEMM_HANDLE%COMM%NODE_ME), MPIDATA, &
              C, PGEMM_HANDLE%NCTOT, PGEMM_HANDLE%OFFDATA, MPIDATA, PGEMM_HANDLE%COMM%MPI_COMM, INFO)
      ELSE
         WRITE(*,*) 'internal error in PARALLEL_GGEMM: the second matrix needs to be stored "N"'
         CALL M_exit(); stop
      ENDIF
# 338


    END SUBROUTINE PARALLEL_GGEMM
END MODULE dfast


!************************ SUBROUTINE LINBAS ****************************
!
! build linear combinations of set of vectors according to matrix CTRANS
! this subroutine performes implicitly a MATRIX x MATRIX multiplication,
! it is needed for the unitary transformation of the wavefunctions or
! for orthogonalisation routines and uses a blocked algorithm
! to save storage
! LINBAS is only called from LINCOM
!
! LINBAS_BLAS  allows only in place transformation of a set of
! of orbitals without adding another set
! furthermore the calling sequence is more BLAS like
!
!***********************************************************************

    SUBROUTINE LINBAS(MODE,CF,CBLOCK,CTRANS,NIN,NOUT,NPL, &
     &           NPLDIM,LDTRAN,NBLK,CFA)
      USE prec

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

! Might be that ZGEMM performs much much faster than ZTRMM ... ??
# 369

      PARAMETER(IUSETR=1)


      CHARACTER (1) MODE
      LOGICAL     LTRI,LADD,LBOTH,LTRANS
      REAL(q)     CF(NPLDIM,NIN),CFA(NPLDIM,NIN)
      REAL(q)     CBLOCK(NBLK,LDTRAN)
      REAL(q)     CTRANS(LDTRAN,NIN)

      IF (NOUT>NIN) THEN
         WRITE(*,1)
 1       FORMAT('internal error in routine LINBAS: wrong arguments, NOUT>NIN')
         CALL M_exit(); stop
      ENDIF

      LTRI =(MODE=='U').OR.(MODE=='u').OR. &
     &      (MODE=='L').OR.(MODE=='l')
      LADD =(MODE=='A').OR.(MODE=='a')
      LBOTH=(MODE=='B').OR.(MODE=='b')
      LTRANS=(MODE=='T').OR.(MODE=='t').OR.(MODE=='C').OR.(MODE=='c')
      IF (LTRI.AND.(IUSETR==0)) THEN
         DO 4 N2=1,NIN
            IF ((MODE=='L').OR.(MODE=='l')) THEN
!DIR$ IVDEP
!OCL NOVREC
               DO N1=1,N2-1
               CTRANS(N1,N2)= 0._q
               ENDDO
            ELSE
!DIR$ IVDEP
!OCL NOVREC
               DO N1=N2+1,NIN
               CTRANS(N1,N2)= 0._q
               ENDDO
            ENDIF
    4    ENDDO
      ENDIF

! Try to get best load balance, maximum block size < NBLK ...
      NBLOCK=NBLK

      DO 70 IBLOCK=0,NPL-1,NBLOCK
         ILENPL=MIN(NBLOCK,NPL-IBLOCK)
         IADDPL=MIN(IBLOCK,NPL-1)
         ILENPL=MAX(ILENPL,0)

         IF (LTRI.AND.(IUSETR/=0)) THEN
! 'Triangular update':

            CALL DTRMM &
# 422

     &                ('R',MODE,'N','N', ILENPL,NOUT,1._q, &
     &                 CTRANS,LDTRAN,CF(1+IADDPL,1), NPLDIM)
         ELSE
! 'Full update':
            IF (LBOTH.OR.(.NOT.LADD)) THEN
               DO 30 N1=1,NIN
                  DO 10 M=1,ILENPL
                     CBLOCK(M,N1)=CF(M+IADDPL,N1)
   10             CONTINUE
   30          CONTINUE
               IF (LTRANS) THEN
               CALL DGEMM('N',MODE, ILENPL, NOUT, NIN, 1._q, &
     &               CBLOCK(1,1), NBLK, CTRANS(1,1), &
     &               LDTRAN, 0._q, CF(IADDPL+1,1), NPLDIM)
               ELSE
               CALL DGEMM('N', 'N', ILENPL, NOUT, NIN, 1._q, &
     &               CBLOCK(1,1), NBLK, CTRANS(1,1), &
     &               LDTRAN, 0._q, CF(IADDPL+1,1), NPLDIM)
               ENDIF
            ENDIF
            IF (LBOTH.OR.LADD) THEN
               IADDT=0
               IF (LBOTH) IADDT=NIN
               DO 60 N1=1,NIN
                  DO 40 M=1,ILENPL
                     CBLOCK(M,N1)=CFA(M+IADDPL,N1)
   40             CONTINUE
   60          CONTINUE
               CALL DGEMM('N', 'N', ILENPL, NOUT, NIN, 1._q, &
     &             CBLOCK(1,1), NBLK, CTRANS(1+IADDT,1), &
     &                   LDTRAN, 1._q, CF(IADDPL+1,1), NPLDIM)
            ENDIF
         ENDIF

   70 CONTINUE

      RETURN
    END SUBROUTINE


    SUBROUTINE LINBAS_BLAS(MODE, NPL, NIN, NOUT, CF,NPLDIM, CTRANS, LDTRAN, CBLOCK, NBLK)
      USE prec

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

! Might be that ZGEMM performs much much faster than ZTRMM ... ??
# 472

      PARAMETER(IUSETR=1)


      CHARACTER (1) MODE
      LOGICAL     LTRI,LADD,LBOTH,LTRANS
      REAL(q)     CF(NPLDIM,NIN)
      REAL(q)     CBLOCK(NBLK,LDTRAN)
      REAL(q)     CTRANS(LDTRAN,NIN)

      IF (NOUT>NIN) THEN
         WRITE(*,1)
 1       FORMAT('internal error in routine LINBAS: wrong arguments, NOUT>NIN')
         CALL M_exit(); stop
      ENDIF

      LTRI =(MODE=='U').OR.(MODE=='u').OR. &
     &      (MODE=='L').OR.(MODE=='l')
      LADD =(MODE=='A').OR.(MODE=='a')
      LTRANS=(MODE=='T').OR.(MODE=='t').OR.(MODE=='C').OR.(MODE=='c')
      IF (LTRI.AND.(IUSETR==0)) THEN
         DO 4 N2=1,NIN
            IF ((MODE=='L').OR.(MODE=='l')) THEN
!DIR$ IVDEP
!OCL NOVREC
               DO N1=1,N2-1
               CTRANS(N1,N2)= 0._q
               ENDDO
            ELSE
!DIR$ IVDEP
!OCL NOVREC
               DO N1=N2+1,NIN
               CTRANS(N1,N2)= 0._q
               ENDDO
            ENDIF
    4    ENDDO
      ENDIF

! Try to get best load balance, maximum block size < NBLK ...
      NBLOCK=NBLK

      DO 70 IBLOCK=0,NPL-1,NBLOCK
         ILENPL=MIN(NBLOCK,NPL-IBLOCK)
         IADDPL=MIN(IBLOCK,NPL-1)
         ILENPL=MAX(ILENPL,0)

         IF (LTRI.AND.(IUSETR/=0)) THEN
! 'Triangular update':

            CALL DTRMM &
# 524

     &                ('R',MODE,'N','N', ILENPL,NOUT,1._q, &
     &                 CTRANS,LDTRAN,CF(1+IADDPL,1), NPLDIM)
         ELSE
! 'Full update':
            IF (LBOTH.OR.(.NOT.LADD)) THEN
               DO 30 N1=1,NIN
                  DO 10 M=1,ILENPL
                     CBLOCK(M,N1)=CF(M+IADDPL,N1)
   10             CONTINUE
   30          CONTINUE
               IF (LTRANS) THEN
               CALL DGEMM('N',MODE, ILENPL, NOUT, NIN, 1._q, &
     &               CBLOCK(1,1), NBLK, CTRANS(1,1), &
     &               LDTRAN, 0._q, CF(IADDPL+1,1), NPLDIM)
               ELSE
               CALL DGEMM('N', 'N', ILENPL, NOUT, NIN, 1._q, &
     &               CBLOCK(1,1), NBLK, CTRANS(1,1), &
     &               LDTRAN, 0._q, CF(IADDPL+1,1), NPLDIM)
               ENDIF
            ENDIF
         ENDIF

   70 CONTINUE

      RETURN
    END SUBROUTINE

!************************ SUBROUTINE ORTH1 *****************************
! calculates a stripe of columns of the overlap between two set of vectors
!
!   O(I,J) =  <CPTWFP(I) |  CFW(J) > +  <GPROJ(I) | CPROW(J) >
!       J=NPOS ,..., NPOS+NSTRIP-1
!       I=1 ,..., NBANDS
!
! usually CWF might be either equal CPTWFP, and CPROW = Q | GPROJ>
! or it might hold
! H_kinetic + H_local CPTWFP, and CPROW = D | GPROJ>
!
! ORTH1 determines only the lower part of the resultant matrix
!  and is suitable if the resulting matrix is Hermitian
!
! ORTH2 calculates the entire matrix
!
!
!                ....    CFW(NPOS) CFW(NPOS+1) ... CFW(NPOS+NSTRIP-1) ....
!
! C(1)*             X    O2         O2       O2       O2              X
! ...               X    O2         O2       O2       O2              X
! C(NPOS)*          X    C          C        C        C               X
! ...               X    C          C        C        C               X
! C(NPOS+NSTRIP-1)* X    C          C        C        C               X
! ...               X    C          C        C        C               X
!
! in a single call the matrix elements marked with C are calculated
! by ORTH1, ORTH2 calculates additionally calculates those marked with O2
! those marked with X are not updated in a single call
!
!***********************************************************************

    SUBROUTINE ORTH1(CSEL,CPTWFP,CFW,GPROJ,CPROW,NBANDS, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      REAL(q)      GPROJ(NPROD,NBANDS)
      REAL(q)      CPROW(NPROD,NSTRIP)
      REAL(q)      COVL(NBANDS,NBANDS)
      CHARACTER (LEN=*) CSEL

      IF (NSTRIP+NPOS-1 > NBANDS) THEN
        WRITE(*,*)'internal error in ORTH1: dim=',NSTRIP+NPOS,NBANDS
        CALL M_exit(); stop
      ENDIF
!
! update of lower triangular part
!
    IF (CSEL(1:1) == 'L' .OR. CSEL(1:1) == 'l') THEN
      IF (NPL/=0) THEN
! for historic reasons 1._q can "unroll" the loop over-plane
! wave coefficients (not used at present)
      NBLOCK=2* NPL

      DO NPOSPL=1,2* NPL-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,1._q, &
              CPTWFP(NPOSPL,NPOS),2* NPLDIM,CFW(NPOSPL,1), &
              2* NPLDIM,1._q,COVL(NPOS,NPOS),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NBANDS-NPOS+1,NSTRIP,2* NPL-NPOSPL+1,1._q, &
              CPTWFP(NPOSPL,NPOS),2* NPLDIM,CFW(NPOSPL,1), &
              2* NPLDIM,1._q,COVL(NPOS,NPOS),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN
! for historic reasons 1._q can "unroll" the loop non-local
! wave coefficients (not used at present)
      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,1._q, &
              GPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
              NPROD,1._q,COVL(NPOS,NPOS),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NBANDS-NPOS+1,NSTRIP,NPRO-NPOSPR+1,1._q, &
              GPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
              NPROD,1._q,COVL(NPOS,NPOS),NBANDS)
      ENDIF
!
! update of upper triangular part
!
    ELSE IF (CSEL(1:1) == 'U' .OR. CSEL(1:1) == 'u') THEN
      IF (NPL/=0) THEN
      NBLOCK=2* NPL

      DO NPOSPL=1,2* NPL-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NSTRIP+NPOS-1,NSTRIP,NBLOCK,1._q, &
              CPTWFP(NPOSPL,1),2* NPLDIM,CFW(NPOSPL,1), &
              2* NPLDIM,1._q,COVL(1,NPOS),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NSTRIP+NPOS-1,NSTRIP,2* NPL-NPOSPL+1,1._q, &
              CPTWFP(NPOSPL,1),2* NPLDIM,CFW(NPOSPL,1), &
              2* NPLDIM,1._q,COVL(1,NPOS),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NSTRIP+NPOS-1,NSTRIP,NBLOCK,1._q, &
              GPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,1._q,COVL(1,NPOS),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NSTRIP+NPOS-1,NSTRIP,NPRO-NPOSPR+1,1._q, &
              GPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,1._q,COVL(1,NPOS),NBANDS)
      ENDIF

    ELSE
      WRITE(*,*)'internal error in ORTH1: CSEL=',CSEL
    ENDIF

      RETURN
    END SUBROUTINE ORTH1



!************************ SUBROUTINE ORTH2 *****************************
! calculates a stripe of columns of the overlap between two set of vectors
!
!   O(I,J) =  <CPTWFP(I) |  CFa(J) > +  <GPROJ(I) | CPROW(J) >
!       J=NPOS ,..., NPOS+NSTRIP-1
!       I=1 ,..., NBANDS
!
!
! for general case (i.e. resulting matrix must not be Hermitian)
! see comments in ORTH1
!***********************************************************************

    SUBROUTINE ORTH2(CPTWFP,CFW,GPROJ,CPROW,NBANDS, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      REAL(q)      GPROJ(NPROD,NBANDS)
      REAL(q)      CPROW(NPROD,NSTRIP)
      REAL(q)      COVL(NBANDS,NBANDS)

      IF (NPL/=0) THEN
! here external blocking can be 1._q, but pretty useless
      NBLOCK=2* NPL
      DO NPOSPL=1,2* NPL-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NBANDS,NSTRIP,NBLOCK,1._q, &
     &         CPTWFP(NPOSPL,1),2* NPLDIM,CFW(NPOSPL,1), &
     &         2* NPLDIM,1._q,COVL(1,NPOS),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NBANDS,NSTRIP,2* NPL-NPOSPL+1,1._q, &
     &         CPTWFP(NPOSPL,1),2* NPLDIM,CFW(NPOSPL,1), &
     &         2* NPLDIM,1._q,COVL(1,NPOS),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN
! here external blocking can be 1._q, but pretty useless
      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NBANDS,NSTRIP,NBLOCK,1._q, &
     &         GPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
     &          NPROD,1._q,COVL(1,NPOS),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NBANDS,NSTRIP,NPRO-NPOSPR+1,1._q, &
     &         GPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
     &          NPROD,1._q,COVL(1,NPOS),NBANDS)
      ENDIF

      RETURN
    END SUBROUTINE ORTH2

!
! similar to the previous subroutine
! but this 1._q allows for a different leading dimension NDIM
! for the matrix to be set up (COVL)
!
    SUBROUTINE ORTH3(CPTWFP,CFW,GPROJ,CPROW, &
         & NBANDS,NSTRIP,NDIM,NPOS, &
         & NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      REAL(q)      GPROJ(NPROD,NBANDS)
      REAL(q)      CPROW(NPROD,NSTRIP)
      REAL(q)      COVL(NDIM,*)

      IF (NPL/=0) THEN
! here external blocking can be 1._q, but pretty useless
         NBLOCK=2* NPL
         DO NPOSPL=1,2* NPL-NBLOCK,NBLOCK
            CALL DGEMM('T','N',NBANDS,NSTRIP,NBLOCK,1._q, &
                 &         CPTWFP(NPOSPL,1),2* NPLDIM,CFW(NPOSPL,1), &
                 &         2* NPLDIM,1._q,COVL(1,NPOS),NDIM)
         ENDDO
         CALL DGEMM('T','N',NBANDS,NSTRIP,2* NPL-NPOSPL+1,1._q, &
              &         CPTWFP(NPOSPL,1),2* NPLDIM,CFW(NPOSPL,1), &
              &         2* NPLDIM,1._q,COVL(1,NPOS),NDIM)
      ENDIF
      
      IF (NPRO/=0) THEN
! here external blocking can be 1._q, but pretty useless
         NBLOCK=NPRO
         DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
            CALL DGEMM('T','N',NBANDS,NSTRIP,NBLOCK,1._q, &
                 &         GPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
                 &          NPROD,1._q,COVL(1,NPOS),NDIM)
         ENDDO
         CALL DGEMM('T','N',NBANDS,NSTRIP,NPRO-NPOSPR+1,1._q, &
              &         GPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              &          NPROD,1._q,COVL(1,NPOS),NDIM)
      ENDIF
      
      RETURN
    END SUBROUTINE ORTH3

!************************ SUBROUTINE ORTH1_NOSUBINDEX ******************
!
!   O(I,J) =  <CPTWFP(I) |  CFa(J) > +  <GPROJ(I) | CPROW(J) >
!       J=NPOS ,..., NPOS+NSTRIP-1
!       I=1 ,..., NBANDS
!
! this routine is for use in 1 aware versions of subrot.F and
! choleski2.F
! the ORTH1 version accesses COVL at the storage position NPOS.
! the ORTH1_NOSUBINDEX version accesses COVL at storage position 1
! regardless of NPOS
!
!***********************************************************************


    SUBROUTINE ORTH1_NOSUBINDEX(CSEL,CPTWFP,CFW,GPROJ,CPROW,NBANDS, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      REAL(q)      GPROJ(NPROD,NBANDS)
      REAL(q)      CPROW(NPROD,NSTRIP)
      REAL(q)      COVL(NBANDS,NBANDS)
      CHARACTER*(*) CSEL

      IF (NSTRIP+NPOS-1 > NBANDS) THEN
        WRITE(*,*)'internal error in ORTH1_NOSUBINDEX: dim=',NSTRIP+NPOS,NBANDS
        CALL M_exit(); stop
      ENDIF
!
! update of lower triangular part
!
    IF (CSEL(1:1) == 'L' .OR. CSEL(1:1) == 'l') THEN
      IF (NPL/=0) THEN
      NBLOCK=2* NPL

      DO NPOSPL=1,2* NPL-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,1._q, &
              CPTWFP(NPOSPL,NPOS),2* NPLDIM,CFW(NPOSPL,1), &
              2* NPLDIM,1._q,COVL(NPOS,1),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NBANDS-NPOS+1,NSTRIP,2* NPL-NPOSPL+1,1._q, &
              CPTWFP(NPOSPL,NPOS),2* NPLDIM,CFW(NPOSPL,1), &
              2* NPLDIM,1._q,COVL(NPOS,1),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,1._q, &
              GPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
              NPROD,1._q,COVL(NPOS,1),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NBANDS-NPOS+1,NSTRIP,NPRO-NPOSPR+1,1._q, &
              GPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
              NPROD,1._q,COVL(NPOS,1),NBANDS)
      ENDIF
!
! update of upper triangular part
!
    ELSE IF (CSEL(1:1) == 'U' .OR. CSEL(1:1) == 'u') THEN
      IF (NPL/=0) THEN
      NBLOCK=2* NPL

      DO NPOSPL=1,2* NPL-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NSTRIP+NPOS-1,NSTRIP,NBLOCK,1._q, &
              CPTWFP(NPOSPL,1),2* NPLDIM,CFW(NPOSPL,1), &
              2* NPLDIM,1._q,COVL(1,1),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NSTRIP+NPOS-1,NSTRIP,2* NPL-NPOSPL+1,1._q, &
              CPTWFP(NPOSPL,1),2* NPLDIM,CFW(NPOSPL,1), &
              2* NPLDIM,1._q,COVL(1,1),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NSTRIP+NPOS-1,NSTRIP,NBLOCK,1._q, &
              GPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,1._q,COVL(1,1),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NSTRIP+NPOS-1,NSTRIP,NPRO-NPOSPR+1,1._q, &
              GPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,1._q,COVL(1,1),NBANDS)
      ENDIF

    ELSE
      WRITE(*,*)'internal error in ORTH1_NOSUBINDEX: CSEL=',CSEL
    ENDIF

    RETURN
    END SUBROUTINE ORTH1_NOSUBINDEX

!************************ SUBROUTINE ORTH2_NOSUBINDEX ******************
!   O(I,J) =  <CPTWFP(I) |  CFa(J) > +  <GPROJ(I) | CPROW(J) >
!       J=NPOS ,..., NPOS+NSTRIP-1
!       I=1 ,..., NBANDS
!
! this routine is for use in 1 aware versions of subrot.F and
! choleski2.F
! the ORTH2 version accesses COVL at the storage position NPOS.
! the ORTH2_NOSUBINDEX version accesses COVL at storage position 1
! regardless of NPOS
!
!***********************************************************************


    SUBROUTINE ORTH2_NOSUBINDEX(CPTWFP,CFW,GPROJ,CPROW,NBANDS, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      REAL(q)      GPROJ(NPROD,NBANDS)
      REAL(q)      CPROW(NPROD,NSTRIP)
      REAL(q)      COVL(NBANDS,NBANDS)
!
! update of upper triangular part
!
      IF (NPL/=0) THEN
      NBLOCK=2* NPL

      DO NPOSPL=1,2* NPL-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NBANDS,NSTRIP,NBLOCK,1._q, &
              CPTWFP(NPOSPL,1),2* NPLDIM,CFW(NPOSPL,1), &
              2* NPLDIM,1._q,COVL(1,1),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NBANDS,NSTRIP,2* NPL-NPOSPL+1,1._q, &
              CPTWFP(NPOSPL,1),2* NPLDIM,CFW(NPOSPL,1), &
              2* NPLDIM,1._q,COVL(1,1),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL DGEMM('T','N',NBANDS,NSTRIP,NBLOCK,1._q, &
              GPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,1._q,COVL(1,1),NBANDS)
      ENDDO
      CALL DGEMM('T','N',NBANDS,NSTRIP,NPRO-NPOSPR+1,1._q, &
              GPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,1._q,COVL(1,1),NBANDS)
      ENDIF

    RETURN
    END SUBROUTINE ORTH2_NOSUBINDEX

!************************ SUBROUTINE LINCOM_SLICE **********************
!
! operates on a block of the matrix CF calculates a linear combination
!  COUT_n,k =   sum_kp CIN_n,kp CH_ kp,k
! this is the version required for 1 aware operation
!
!***********************************************************************


    SUBROUTINE LINCOM_SLICE(MODE,CF,CPROF,CTRANS,NIN,NPOS,NOUT,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN,CF_RESULT,CPROF_RESULT)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      CHARACTER*1 MODE
      DIMENSION   CF(NPLDIM,NIN),CFA(NPLDIM,NIN)
      REAL(q)        CPROF(NPROD,NIN),CPROFA(NPROD,NIN)
      REAL(q)        CTRANS(LDTRAN,NOUT)
      INTEGER     NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN

      CALL LINBAS_SLICE(MODE,CF,CF_RESULT,CTRANS,NIN,NPOS,NOUT,2* NPL, &
     &            2* NPLDIM,LDTRAN)
      IF (NPRO/=0) THEN
      CALL LINBAS_SLICE(MODE,CPROF,CPROF_RESULT,CTRANS,NIN,NPOS,NOUT, NPRO, &
     &             NPROD,LDTRAN)
      ENDIF

      RETURN
    END SUBROUTINE LINCOM_SLICE
      


    SUBROUTINE LINBAS_SLICE(MODE,CF,CF_RESULT,CTRANS,NIN,NPOS,NOUT,NPL, &
     &           NPLDIM,LDTRAN)
      USE prec

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      INTEGER NIN,NPOS,NOUT,NPL,NPLDIM,LDTRAN

      CHARACTER*1 MODE
      LOGICAL     LTRI
      REAL(q)     CF(NPLDIM,NIN)
      REAL(q)     CF_RESULT(NPLDIM,NIN)
      REAL(q)     CTRANS(LDTRAN,NOUT)

      IF (NOUT>NIN) THEN
         WRITE(*,1)
 1       FORMAT('internal error in routine LINBAS: wrong arguments, NOUT>NIN')
         CALL M_exit(); stop
      ENDIF

      LTRI =(MODE=='U').OR.(MODE=='u').OR. &
     &      (MODE=='L').OR.(MODE=='l')
! this part makes sure that the upper or lower triangle is cleared
! if we use DGEMM instead of ZTRMM
      IF (LTRI) THEN
         DO 4 N2=1,NOUT
            IF ((MODE=='L').OR.(MODE=='l')) THEN
               WRITE(*,*) ' internal error: lower triangle is not yet supported'
               CALL M_exit(); stop
!DIR$ IVDEP
!OCL NOVREC
               DO N1=1,N2+NPOS-2
                  CTRANS(N1,N2)= 0._q
               ENDDO
            ELSE
!DIR$ IVDEP
!OCL NOVREC
               DO N1=N2+NPOS,NIN
                  CTRANS(N1,N2)= 0._q
               ENDDO
            ENDIF
    4    ENDDO
      ENDIF

      ILENPL=NPL

      IF (LTRI) THEN
! 'Triangular update':
         CALL DGEMM('N', 'N', ILENPL, NOUT, NOUT+NPOS-1, 1._q, &
     &        CF(1,1), NPLDIM, CTRANS(1,1), &
     &        LDTRAN, 0._q, CF_RESULT(1,1), NPLDIM)
      ELSE
! 'Full update':
         CALL DGEMM('N', 'N', ILENPL, NOUT, NIN, 1._q, &
     &        CF(1,1), NPLDIM, CTRANS(1,1), &
     &        LDTRAN, 0._q, CF_RESULT(1,1), NPLDIM)
      ENDIF

      RETURN
    END SUBROUTINE LINBAS_SLICE
