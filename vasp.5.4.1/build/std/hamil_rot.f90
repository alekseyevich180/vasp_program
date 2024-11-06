# 1 "hamil_rot.F"
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

# 2 "hamil_rot.F" 2 
!======================================================================
! RCS:  $Id: hamil_rot.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!  MODUL contains some routines only used for all bands simultaneous
!  update scheme implemented in rot.f
!  dont worry about this routines if you do not use rot.f
!
!======================================================================

!************************* SUBROUTINE HAMNNP ***************************
!
! this subroutine calculates the Hamiltonian enclosed between
! (1._q,0._q) wavefunction and a set of wavefunctions
!   H_n' n = <C_n | H | C_n'>
!***********************************************************************

      SUBROUTINE HAMNNP(GRID,NONLR_S,NONL_S,W,WDES,CR, LREAL,SV,N,NK,ISP, &
         LMDIM,CDIJ,CHAM, CACC)
      USE prec

      
      USE nonl_high
      USE mpimy
      USE mgrid
      USE wave

      IMPLICIT NONE

      TYPE (grid_3d)     GRID
      TYPE (nonlr_struct)NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES

      LOGICAL LREAL
      INTEGER N     !  band index
      INTEGER NK    !  k-point inded
      INTEGER ISP   ! spin component
      INTEGER LMDIM ! leading dimension of CDIJ and CQIJ arrays

      REAL(q)    CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q) CR(GRID%RL%NP*WDES%NRSPINORS)
      REAL(q)      SV(GRID%MPLWV*2,WDES%NRSPINORS*WDES%NRSPINORS) ! local potential
      COMPLEX(q)       CHAM(WDES%NBANDS,WDES%NBANDS)
      COMPLEX(q) CACC(WDES%NRPLWV)
! local work arrays
      TYPE (wavedes1)    WDES1
      COMPLEX(q)       DWORK(WDES%NBANDS)
      COMPLEX(q)       CPROF(WDES%NPROD)
      COMPLEX(q) CWORK1(GRID%MPLWV*WDES%NRSPINORS)
      INTEGER NPL      ! number of plane wave coefficient
      INTEGER NGVECTOR ! number of plane wave coefficients per spinor

      INTEGER ISPINOR,ISPINOR_,NPRO,NPRO_,NIS,NT,LMMAXC,NI,L,LP
      INTEGER M,MM,N1,N2

      NPL=WDES%NPLWKP(NK)
      NGVECTOR=WDES%NGVECTOR(NK)

      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,GRID)

!-----------------------------------------------------------------------
! calculate the local contribution
!-----------------------------------------------------------------------
      CWORK1=0
      CALL VHAMIL(WDES1,GRID,SV(1,1),CR(1),CWORK1(1)) 
!-----------------------------------------------------------------------
! transform back to reciprocal space and add kinetic energy
!-----------------------------------------------------------------------
      DO ISPINOR=0,WDES1%NRSPINORS-1
         CALL FFTEXT_MPI(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV),CACC(1+ISPINOR*NGVECTOR),GRID,.FALSE.)
         DO M=1,NGVECTOR
            MM=M+ISPINOR*NGVECTOR
            CACC(MM)=CONJG(CACC(MM)+W%CPTWFP(MM,N,NK,ISP)* WDES1%DATAKE(M,ISPINOR+1))
         ENDDO
      ENDDO
!=======================================================================
! calculate the non local part
!=======================================================================
      CPROF=0

      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
      DO ISPINOR_=0,WDES1%NRSPINORS-1

      NPRO =ISPINOR *(WDES1%NPRO/2)
      NPRO_=ISPINOR_*(WDES1%NPRO/2)

      NIS =1
      DO NT=1,WDES1%NTYP
        LMMAXC=WDES1%LMMAX(NT)
        IF (LMMAXC==0) GOTO 100
        DO NI=NIS,WDES1%NITYP(NT)+NIS-1

          DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
          DO LP=1,LMMAXC
            CPROF(L+NPRO)=CPROF(L+NPRO)+ CONJG(CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+ISP)*W%CPROJ(LP+NPRO_,N,NK,ISP))
          ENDDO; ENDDO
          NPRO = LMMAXC+NPRO
          NPRO_= LMMAXC+NPRO_
        ENDDO
 100    NIS = NIS+WDES1%NITYP(NT)
      ENDDO
      ENDDO
      ENDDO spinor
!=======================================================================
! calculate the Hamilton-matrix defined as
! CHAM(n2,n1) = < C(n1) | H | C(n2) >
!=======================================================================
      CALL ZGEMV( 'T', NPL , WDES%NBANDS-N+1 ,(1._q,0._q) ,W%CPTWFP(1,N,NK,ISP), &
     &              WDES%NRPLWV, CACC, 1 , (0._q,0._q) ,  DWORK, 1)

      CALL ZGEMV( 'T',  WDES%NPRO, WDES%NBANDS-N+1 ,(1._q,0._q) ,W%CPROJ(1,N,NK,ISP), &
     &             WDES%NPROD, CPROF, 1 , (1._q,0._q) ,  DWORK, 1)

      N1=N
      DO N2=N1,WDES%NBANDS
        CHAM(N2,N1)= (DWORK(N2-N1+1))
        CHAM(N1,N2)= (CONJG(DWORK(N2-N1+1)))
      ENDDO

      DO M=1,NPL
         CACC(M)=CONJG(CACC(M))
      ENDDO

      RETURN 
      END SUBROUTINE



!************************* SUBROUTINE HAMIL0 ***************************
!
! this subroutine calculates the gradient vector of the band N, if the
! eigenmatrix is not diagonal   i.e.
!   H | C_n >- H_n n' S | C_n' >
! especially  due to the overlap-operator this is a somewhat complicated
! operation
!
! CHAM contains the (possibly non diagonal) sub-space Hamilton matrix
! the  wavefunction must be given in reciprocal space C(..,N) and real
! space CR and the projection operators must be stored in CPROJ(..,N)
! the routine assumes that the local part and the kinetic energy part of
! the Hamiltonian are already stored in CACC
! the remaining contributions are added to CACC
!
!***********************************************************************

      SUBROUTINE HAMIL0(GRID,NONLR_S,NONL_S,W,WDES,CR, LREAL,SV,N,NK,ISP, &
         LMDIM,CQIJ,CDIJ,CHAM,  CACC)
      USE prec

      
      USE nonl_high
      USE mpimy
      USE mgrid
      USE wave
      IMPLICIT NONE

      TYPE (grid_3d)     GRID
      TYPE (nonlr_struct)NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavedes)     WDES
      TYPE (wavedes1)    WDES1
      TYPE (wavespin)    W

      LOGICAL LREAL

      INTEGER N     !  band index
      INTEGER NK    !  k-point inded
      INTEGER ISP   ! spin component
      INTEGER LMDIM ! leading dimension of CDIJ and CQIJ arrays

      REAL(q)    CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q) CR(GRID%RL%NP*WDES%NRSPINORS)
      REAL(q)      SV(GRID%MPLWV*2,WDES%NRSPINORS*WDES%NRSPINORS)
      COMPLEX(q)       CHAM(WDES%NBANDS,WDES%NBANDS)
      COMPLEX(q) CACC(WDES%NRPLWV)
! work arrays
      COMPLEX(q)       CWORK(WDES%NBANDS)

      COMPLEX(q)       CWORK2(WDES%NPROD),CWORK3(WDES%NPROD)
      COMPLEX(q) CPTACC(GRID%MPLWV*WDES%NRSPINORS)
      INTEGER NPL      ! number of plane wave coefficient
      INTEGER NGVECTOR ! number of plane wave coefficients per spinor
      INTEGER ISPINOR,ISPINOR_,NPRO,NPRO_,NIS,NT,LMMAXC,NI,L,LP
      INTEGER M,MM,NP


      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,GRID)
!=======================================================================
! calculate sum_i n' Q_ij j H_nn' C_in
!=======================================================================

      NPL=WDES%NPLWKP(NK)
      NGVECTOR=WDES%NGVECTOR(NK)

!---First store the elements of the Hamiltonmatrix in the correct order
      DO NP=1,WDES%NBANDS
        CWORK(NP)=CHAM(N,NP)
      ENDDO
!-----------------------------------------------------------------------
!  now sum_n  -H_nn'  C_i n'     and   sum_n -H_nn' C_Gn'
!-----------------------------------------------------------------------
      CALL ZGEMM( 'N', 'N' , NPL , 1 , WDES%NBANDS , -(1._q,0._q), &
     &             W%CPTWFP(1,1,NK,ISP), WDES%NRPLWV , CWORK , WDES%NBANDS , &
     &             (1._q,0._q) , CACC , WDES%NRPLWV )

      CALL ZGEMM( 'N', 'N' , WDES%NPRO , 1 , WDES%NBANDS  , -(1._q,0._q), &
     &             W%CPROJ(1,1,NK,ISP),  WDES%NPROD , CWORK , WDES%NBANDS , &
     &             (0._q,0._q) , CWORK3 , WDES%NPROD  )

!-----------------------------------------------------------------------
!  now sum_ij Q_ij  (sum_np -H_n np  C_i np)   + D_ij C_in
!-----------------------------------------------------------------------
      CWORK2=0
      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
      DO ISPINOR_=0,WDES1%NRSPINORS-1

      NPRO =ISPINOR *(WDES1%NPRO/2)
      NPRO_=ISPINOR_*(WDES1%NPRO/2)

      NIS =1
      DO NT=1,WDES1%NTYP
        LMMAXC=WDES1%LMMAX(NT)
        IF (LMMAXC==0) GOTO 100
        DO NI=NIS,WDES1%NITYP(NT)+NIS-1

          DO L =1,LMMAXC
          DO LP=1,LMMAXC
           CWORK2(L+NPRO)=CWORK2(L+NPRO)+CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+ISP)*W%CPROJ (LP+NPRO_,N,NK,ISP) &
                                        +CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+ISP)*CWORK3(LP+NPRO_)
          ENDDO; ENDDO
          NPRO = LMMAXC+NPRO
          NPRO_= LMMAXC+NPRO_
        ENDDO
 100    NIS = NIS+WDES1%NITYP(NT)
      ENDDO
      ENDDO
      ENDDO spinor

      CPTACC=0
!=======================================================================
! non-local contribution using real-space projection add to CPTACC
! and FFT, finally add to CACC
!=======================================================================
      IF (LREAL) THEN
         CALL RACC0(NONLR_S,WDES1,CWORK2(1),CPTACC(1))
         DO ISPINOR=0,WDES1%NRSPINORS-1
            CALL FFTEXT_MPI(NGVECTOR,WDES%NINDPW(1,NK),CPTACC(1+ISPINOR*WDES1%GRID%MPLWV),CACC(1+ISPINOR*NGVECTOR),GRID,.TRUE.)
         ENDDO
!=======================================================================
! calculate the non local contribution using reciprocal space add to CACC
!=======================================================================
      ELSE
        CALL VNLAC0(NONL_S,WDES1,CWORK2(1),CACC(1))
      ENDIF
      RETURN
      END SUBROUTINE
