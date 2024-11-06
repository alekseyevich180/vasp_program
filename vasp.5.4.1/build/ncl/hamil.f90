# 1 "hamil.F"
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

# 2 "hamil.F" 2 

! work around PGI compiler bug

MODULE hamil_helper
  INTERFACE
     SUBROUTINE PW_NORM_WITH_METRIC(WDES1, C, FNORM, FMETRIC, METRIC)
       USE prec
       USE wave
       REAL(q) FNORM
       TYPE (wavedes1)    WDES1
       COMPLEX(q) :: C
       REAL(q), OPTIONAL :: FMETRIC
       REAL(q), OPTIONAL :: METRIC(WDES1%NPL)
     END SUBROUTINE PW_NORM_WITH_METRIC
  END INTERFACE
END MODULE hamil_helper


MODULE hamil
  USE prec
  USE wave_high
!***********************************************************************
!
! this module implements most of the low and high level
! routines to calculated the action of a local Hamiltonian
! onto the wavefunctions
!
!***********************************************************************

  USE hamil_helper
# 44


  INTERFACE
     SUBROUTINE KINHAMIL( WDES1, GRID, CVR,  LADD, DATAKE, EVALUE, CPTWFP, CH)
       USE mgrid
       USE wave

       TYPE (wavedes1)    WDES1
       TYPE (grid_3d)     GRID
       COMPLEX(q) :: CVR
       LOGICAL    :: LADD
       REAL(q)    :: DATAKE
       REAL(q)    :: EVALUE
       COMPLEX(q) :: CPTWFP
       COMPLEX(q) :: CH
     END SUBROUTINE KINHAMIL
  END INTERFACE

  INTERFACE
     SUBROUTINE VHAMIL(WDES1,GRID,SV,CR,CVR)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       COMPLEX(q)   SV
       COMPLEX(q) :: CR
       COMPLEX(q) :: CVR
     END SUBROUTINE VHAMIL
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE(WDES1,  CHARGE, NDIM, CR1, CR2, WEIGHT)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       INTEGER NDIM
       COMPLEX(q)   CHARGE
       COMPLEX(q) :: CR1,CR2
       REAL(q) :: WEIGHT
     END SUBROUTINE PW_CHARGE
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_CMPLX(WDES1,  CHARGE, NDIM, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       INTEGER NDIM
       COMPLEX(q)   CHARGE
       COMPLEX(q) :: CR1,CR2
       REAL(q) :: WEIGHT
     END SUBROUTINE PW_CHARGE_CMPLX
  END INTERFACE

  INTERFACE
     SUBROUTINE VHAMIL_TRACE(WDES1, GRID, SV, CR, CVR, WEIGHT)
       USE prec
       USE mgrid
       USE wave

       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1

       COMPLEX(q)   SV
       COMPLEX(q) :: CR,CVR
       REAL(q)    :: WEIGHT
     END SUBROUTINE VHAMIL_TRACE
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_TRACE(WDES1,  CHARGE, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       COMPLEX(q)   CHARGE
       COMPLEX(q) :: CR1,CR2
     END SUBROUTINE PW_CHARGE_TRACE
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_TRACE_GDEF(WDES1,  CHARGE, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       COMPLEX(q)   CHARGE
       COMPLEX(q) :: CR1,CR2
     END SUBROUTINE PW_CHARGE_TRACE_GDEF
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_TRACE_NO_CONJG(WDES1,  CHARGE, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       COMPLEX(q)   CHARGE
       COMPLEX(q) :: CR1,CR2
     END SUBROUTINE PW_CHARGE_TRACE_NO_CONJG
  END INTERFACE

  INTERFACE
     SUBROUTINE ECCP_NL_(LMDIM,LMMAXC,CDIJ,CQIJ,EVALUE,CPROJ1,CPROJ2,CNL)
       USE prec
       IMPLICIT NONE
       COMPLEX(q)  CNL
       INTEGER LMMAXC, LMDIM
       COMPLEX(q) CDIJ,CQIJ
       REAL(q) EVALUE
       COMPLEX(q) CPROJ1,CPROJ2
     END SUBROUTINE ECCP_NL_
  END INTERFACE
CONTAINS

!************************* SUBROUTINE ECCP   ***************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
! this subroutine calculates the expectation value of <c|H|cp>
! where c and cp are two wavefunctions
!  FFT and non-local projections of wavefunctions st be supplied
!***********************************************************************

  SUBROUTINE ECCP(WDES1,W1,W2,LMDIM,CDIJ,GRID,SV, CE)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)


    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    INTEGER NGVECTOR, ISPINOR
    COMPLEX(q)      CNL
    COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
!=======================================================================
! calculate the local contribution
!=======================================================================
    CLOCAL=0
    NGVECTOR=WDES1%NGVECTOR

    DO ISPINOR =0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
          DO M=1,GRID%RL%NP
             MM =M+ISPINOR *GRID%MPLWV
             MM_=M+ISPINOR_*GRID%MPLWV
             CLOCAL=CLOCAL+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
          ENDDO
       ENDDO
    ENDDO

    CLOCAL=CLOCAL/GRID%NPLWV
!=======================================================================
! kinetic energy contribution
!=======================================================================
    CKIN=0
    DO ISPINOR=0,WDES1%NRSPINORS-1
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CKIN=CKIN+W1%CPTWFP(MM)*CONJG(W2%CPTWFP(MM))*WDES1%DATAKE(M,ISPINOR+1)
       ENDDO
    ENDDO
!=======================================================================
! non local contribution
!=======================================================================
    CNL =0
    NPRO=0
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1

          NPRO =ISPINOR *(WDES1%NPRO/2)
          NPRO_=ISPINOR_*(WDES1%NPRO/2)

          NIS =1
          DO NT=1,WDES1%NTYP
             LMMAXC=WDES1%LMMAX(NT)
             IF (LMMAXC==0) GOTO 310
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                CALL ECCP_NL(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CNL)
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
310          NIS = NIS+WDES1%NITYP(NT)
          ENDDO
       ENDDO
    ENDDO spinor

    CE=(CLOCAL+CKIN+CNL)
    CALL M_sum_z(WDES1%COMM_INB, CE, 1)

  END SUBROUTINE ECCP

!************************* SUBROUTINE ECCP   ***************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
! this subroutine calculates the expectation value of <c|H|cp>
! where c and cp are two wavefunctions
!  FFT and non-local projections of wavefunctions must be supplied
!***********************************************************************

  SUBROUTINE ECCP_VEC(WDES1,W1,W2,LMDIM,CDIJ,GRID,SV,AVEC, CE)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)


    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    INTEGER NGVECTOR, ISPINOR
    COMPLEX(q)      CNL
    COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q)   AVEC(:,:)

    COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK1(WDES1%NRPLWV)
    COMPLEX(q) :: CWORK2(WDES1%GRID%MPLWV)

    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
!=======================================================================
! calculate the local contribution
!=======================================================================
    CLOCAL=0
    NGVECTOR=WDES1%NGVECTOR

    DO ISPINOR =0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
          DO M=1,GRID%RL%NP
             MM =M+ISPINOR *GRID%MPLWV
             MM_=M+ISPINOR_*GRID%MPLWV
             CLOCAL=CLOCAL+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
          ENDDO
       ENDDO
    ENDDO

    CLOCAL=CLOCAL/GRID%NPLWV
!=======================================================================
! kinetic energy contribution
!=======================================================================
    CVR=0
    CALL KINHAMIL_VEC( WDES1, GRID, CVR,  .FALSE., & 
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
                  AVEC(1,1), CWORK1, CWORK2, 0.0_q, W1%CPTWFP(1), CWORK1)

    CKIN=0
    DO ISPINOR=0,WDES1%NRSPINORS-1
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CKIN=CKIN+CWORK1(MM)*CONJG(W2%CPTWFP(MM))
       ENDDO
    ENDDO
!=======================================================================
! non local contribution
!=======================================================================
    CNL =0
    NPRO=0

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1

          NPRO =ISPINOR *(WDES1%NPRO/2)
          NPRO_=ISPINOR_*(WDES1%NPRO/2)

          NIS =1
          DO NT=1,WDES1%NTYP
             LMMAXC=WDES1%LMMAX(NT)
             IF (LMMAXC==0) GOTO 310
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                CALL ECCP_NL(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CNL)
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
310          NIS = NIS+WDES1%NITYP(NT)
          ENDDO
       ENDDO
    ENDDO spinor

    CE=(CLOCAL+CKIN+CNL)
    CALL M_sum_z(WDES1%COMM_INB, CE, 1)

  END SUBROUTINE ECCP_VEC


!************************* SUBROUTINE ECCP   ***************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
! this subroutine calculates the expectation value of <c|H|cp>
! where c and cp are two wavefunctions
!  FFT and non-local projections of wavefunctions must be supplied
!***********************************************************************

  SUBROUTINE ECCP_TAU(WDES1,W1,W2,LMDIM,CDIJ,GRID,SV,LATT_CUR,MU,CE)
    USE prec
    USE mpimy
    USE mgrid
    USE lattice
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID
    TYPE (latt)     :: LATT_CUR

    INTEGER NGVECTOR, ISPINOR
    COMPLEX(q)      CNL
    COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q)   MU(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS)

    COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CH(WDES1%NRPLWV),CWORK1(WDES1%NRPLWV),CWORK2(WDES1%GRID%MPLWV)

    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
!=======================================================================
! calculate the local contribution
!=======================================================================
    CLOCAL=0
    NGVECTOR=WDES1%NGVECTOR

    DO ISPINOR =0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
          DO M=1,GRID%RL%NP
             MM =M+ISPINOR *GRID%MPLWV
             MM_=M+ISPINOR_*GRID%MPLWV
             CLOCAL=CLOCAL+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
          ENDDO
       ENDDO
    ENDDO

    CLOCAL=CLOCAL/GRID%NPLWV
!=======================================================================
! kinetic energy contribution
!=======================================================================
    CVR=0; CH=0
    CALL KINHAMIL_TAU( WDES1, GRID, CVR,  .FALSE., .TRUE., & 
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
                  LATT_CUR, MU(1,1), CWORK1, CWORK2, 0.0_q, W1%CPTWFP(1), CH(1))

    CKIN=0
    DO ISPINOR=0,WDES1%NRSPINORS-1
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CKIN=CKIN+CH(MM)*CONJG(W2%CPTWFP(MM))
       ENDDO
    ENDDO
!=======================================================================
! non local contribution
!=======================================================================
    CNL =0
    NPRO=0

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1

          NPRO =ISPINOR *(WDES1%NPRO/2)
          NPRO_=ISPINOR_*(WDES1%NPRO/2)

          NIS =1
          DO NT=1,WDES1%NTYP
             LMMAXC=WDES1%LMMAX(NT)
             IF (LMMAXC==0) GOTO 310
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                CALL ECCP_NL(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CNL)
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
310          NIS = NIS+WDES1%NITYP(NT)
          ENDDO
       ENDDO
    ENDDO spinor

    CE=(CLOCAL+CKIN+CNL)
    CALL M_sum_z(WDES1%COMM_INB, CE, 1)

  END SUBROUTINE ECCP_TAU


!*********************************************************************
!
! calculate the band structure energy
!
!*********************************************************************

  FUNCTION BANDSTRUCTURE_ENERGY(WDES, W) RESULT(E)
    IMPLICIT NONE
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W      ! wavefunction
    REAL(q)            E
    INTEGER            ISP, NK, NB

    E=0
    DO ISP=1,WDES%ISPIN
       DO NK=1,WDES%NKPTS

          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

          DO NB=1,WDES%NB_TOT
             E=E+WDES%RSPIN* REAL( W%CELTOT(NB,NK,ISP) ,KIND=q) *WDES%WTKPT(NK)*W%FERTOT(NB,NK,ISP)
          ENDDO
       ENDDO
    ENDDO
    CALL M_sum_d(WDES%COMM_KINTER, E, 1)
  END FUNCTION BANDSTRUCTURE_ENERGY


!**********************************************************************
!
! calculate the kinetic energy of each wavefunction
!
!**********************************************************************

  SUBROUTINE KINETIC_ENERGY(W )
    USE wave_high
    IMPLICIT NONE
    TYPE (wavespin), TARGET :: W
!    local
    TYPE (wavedes1)    WDES1
    TYPE (wavedes), POINTER :: WDES
    
    INTEGER ISP, NK, NB, ISPINOR, MM, M
    REAL(q) :: CKIN
    
    
    WDES=>W%WDES
    
    spin:  DO ISP=1,WDES%ISPIN
       kpoint: DO NK=1,WDES%NKPTS


          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
!PK MRG_AUX already zeros contributions from non-kpoint nodes, but be cautious and (0._q,0._q) workspace anyway
             DO NB=1,WDES%NBANDS
                W%AUX(NB,NK,ISP)=0.
             END DO
          ELSE

          CALL SETWDES(WDES,WDES1,NK)
          
          DO NB=1,WDES%NBANDS
             CKIN=0
             DO ISPINOR=0,WDES1%NRSPINORS-1
                DO M=1,WDES1%NGVECTOR
                   MM=M+ISPINOR*WDES1%NGVECTOR
                   CKIN=CKIN+W%CPTWFP(MM,NB,NK,ISP)*CONJG(W%CPTWFP(MM,NB,NK,ISP))*WDES1%DATAKE(M,ISPINOR+1)
                ENDDO
             ENDDO
             W%AUX(NB,NK,ISP)=CKIN
          ENDDO

          ENDIF

       END DO kpoint
    END DO spin
    
    CALL MRG_AUX(WDES,W)
  END SUBROUTINE KINETIC_ENERGY
  
!************************* SUBROUTINE SETUP_PRECOND ********************
!
! subroutine to set up a preconditioning matrix
! IALGO determines the type of preconditioning
!
!   0,6,8 TAP preconditioning
!    9    Jacobi like preconditioning
!   else  no preconditioning
!
! the Jacobi like preconditer required the eigenvalue and
! a complex shift
!
!***********************************************************************


  SUBROUTINE SETUP_PRECOND( W1, IALGO, IDUMP, PRECON, EVALUE, DE_ATT )
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    INTEGER IALGO         ! chosen algorithm
    INTEGER IDUMP         ! dump flag
    REAL(q) :: EVALUE     ! eigenvalue minus average local potential
    REAL(q) :: DE_ATT     ! complex shift
    REAL(q) :: PRECON(*)  ! return: preconditioning matrix
! local
    INTEGER M, MM, NK, ISPINOR
    REAL(q) :: EKIN, FAKT, X, X2
    COMPLEX(q) :: CPT

    IF (IALGO==0 .OR. IALGO==8 .OR. IALGO==6) THEN
       EKIN=0

       DO ISPINOR=0,W1%WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
          DO M=1,W1%WDES1%NGVECTOR
             MM=M+ISPINOR*W1%WDES1%NGVECTOR
             CPT=W1%CPTWFP(MM)
             EKIN =EKIN+ REAL( CPT*CONJG(CPT) ,KIND=q) * W1%WDES1%DATAKE(M,ISPINOR+1)
          ENDDO
       ENDDO
       CALL M_sum_d(W1%WDES1%COMM_INB, EKIN, 1)

       IF (EKIN<2.0_q) EKIN=2.0_q
       EKIN=EKIN*1.5_q
       IF (IDUMP==2)  WRITE(*,'(E9.2,"E")',ADVANCE='NO') EKIN

       FAKT=2._q/EKIN

       DO ISPINOR=0,W1%WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
          DO M=1,W1%WDES1%NGVECTOR
             MM=M+ISPINOR*W1%WDES1%NGVECTOR
             X=W1%WDES1%DATAKE(M,ISPINOR+1)/EKIN
             X2= 27+X*(18+X*(12+8*X))
             PRECON(MM)=X2/(X2+16*X*X*X*X)*FAKT
          ENDDO
       ENDDO
    ELSE IF (IALGO==9) THEN
       DO ISPINOR=0,W1%WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
          DO M=1,W1%WDES1%NGVECTOR
             MM=M+ISPINOR*W1%WDES1%NGVECTOR
             X=MAX(W1%WDES1%DATAKE(M,ISPINOR+1)-EVALUE,0._q)
             PRECON(MM)= REAL( 1._q/(X+ CMPLX( 0 , DE_ATT ,KIND=q) ) ,KIND=q) !new
          ENDDO
       ENDDO
    ELSE
       DO ISPINOR=0,W1%WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
          DO M=1,W1%WDES1%NGVECTOR
             MM=M+ISPINOR*W1%WDES1%NGVECTOR
             PRECON(MM)=1.0_q
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE SETUP_PRECOND


  SUBROUTINE APPLY_PRECOND( W1, W2, PRECON, MUL)
    USE gridq
    IMPLICIT NONE
    TYPE (wavefun1)    W1, W2
    REAL(q) :: PRECON(*)
    REAL(q), OPTIONAL :: MUL
    REAL(q) :: MUL_

    IF (PRESENT(MUL)) THEN
       MUL_=MUL
    ELSE
       MUL_=1.0
    ENDIF

    CALL REAL_CMPLX_CMPLX_MUL( W1%WDES1%NPL,  PRECON(1), W1%CPTWFP(1) , MUL_, W2%CPTWFP(1), 0.0_q, W2%CPTWFP(1))
    
  END SUBROUTINE APPLY_PRECOND

  SUBROUTINE ADD_PRECOND( W1, W2, PRECON, MUL)
    USE gridq
    IMPLICIT NONE
    TYPE (wavefun1)    W1, W2
    REAL(q) :: PRECON(*)
    REAL(q), OPTIONAL :: MUL
    REAL(q) :: MUL_

    IF (PRESENT(MUL)) THEN
       MUL_=MUL
    ELSE
       MUL_=1.0
    ENDIF

    CALL REAL_CMPLX_CMPLX_MUL( W1%WDES1%NPL,  PRECON(1), W1%CPTWFP(1) , MUL_, W2%CPTWFP(1), 1.0_q, W2%CPTWFP(1))
    
  END SUBROUTINE ADD_PRECOND


  SUBROUTINE SIMPLE_PRECOND( W1 )
    USE gridq
    IMPLICIT NONE
    TYPE (wavefun1)    W1

    REAL(q) :: EKIN, FAKT, X, X2
    INTEGER :: ISPINOR, M, MM

! old version
    EKIN=10
    FAKT=2._q/EKIN
    DO ISPINOR=0,W1%WDES1%NRSPINORS-1
       DO M=1,W1%WDES1%NGVECTOR
          MM=M+ISPINOR*W1%WDES1%NGVECTOR
          X=W1%WDES1%DATAKE(M,ISPINOR+1)/EKIN
          X2=27+X*(18+X*(12+8*X))
          W1%CPTWFP(MM)= W1%CPTWFP(MM)*X2/(X2+16*X*X*X*X)*FAKT
       ENDDO
    ENDDO

  END SUBROUTINE SIMPLE_PRECOND


!***********************************************************************
!
! truncate high frequency components if LDELAY is set
! or if use_eini_without_ldelay is defined
!
!***********************************************************************

  SUBROUTINE TRUNCATE_HIGH_FREQUENCY_W1( W1, LDELAY, ENINI)
    USE prec
    IMPLICIT NONE
    LOGICAL LDELAY     ! usually truncation only during delay phase (LDELAY set)
    REAL(q) ENINI      ! cutoff at which the wavefunction is truncated
    TYPE (wavefun1)    W1
! local
    INTEGER ISPINOR

    IF (LDELAY.OR.W1%WDES1%LSPIRAL) THEN
       DO ISPINOR=1,W1%WDES1%NRSPINORS
          CALL TRUNCATE_HIGH_FREQUENCY_ONE(W1%WDES1%NGVECTOR, &
               W1%CPTWFP(1+(ISPINOR-1)*W1%WDES1%NGVECTOR), W1%WDES1%DATAKE(1,ISPINOR), ENINI)
       ENDDO
    END IF
  END SUBROUTINE TRUNCATE_HIGH_FREQUENCY_W1


!***********************************************************************
!
! determine norm of a wavefunction and norm with a metric
!
!***********************************************************************

  SUBROUTINE PW_NORM_WITH_METRIC_W1(W1, FNORM, FMETRIC, METRIC)
    IMPLICIT NONE
    REAL(q) FNORM
    TYPE (wavefun1)    W1
    REAL(q), OPTIONAL :: FMETRIC
    REAL(q), OPTIONAL :: METRIC(W1%WDES1%NPL)

    IF (PRESENT(FMETRIC) .AND. PRESENT(METRIC)) THEN
       CALL PW_NORM_WITH_METRIC( W1%WDES1, W1%CPTWFP(1), FNORM, FMETRIC, METRIC)
    ELSE
       CALL PW_NORM_WITH_METRIC( W1%WDES1, W1%CPTWFP(1), FNORM)
    ENDIF
  END SUBROUTINE PW_NORM_WITH_METRIC_W1



!************************* SUBROUTINE HAMILT ***************************
!
! this subroutine calculates the H acting onto a wavefuntion
! the  wavefunction must be given in reciprocal space C and real
! space CR
! CH contains the result
!   H- EVALUE Q |phi>
! where Q is the nondiagonal part of the overlap matrix S
!***********************************************************************

  SUBROUTINE HAMILT( &
       &    W1, NONLR_S, NONL_S, EVALUE, CDIJ, CQIJ, SV, ISP, CH)
    USE mpimy
    USE mgrid
    
    USE nonl_high
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S

    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)      SV(:,:)
    INTEGER    ISP
    COMPLEX(q) CH(:)
    REAL(q)    EVALUE
! local variables
    COMPLEX(q) :: CWORK1(W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS)


    IF (NONLR_S%LREAL) THEN
! calculate the local contribution (result in CWORK1)
       CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CWORK1(1)) 
! non-local contribution in real-space

       CALL RACC(NONLR_S, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)

       CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), .FALSE., &
            W1%WDES1%DATAKE(1,1), 0.0_q, W1%CPTWFP(1), CH(1))
    ELSE
! calculate the local contribution (result in CWORK1)
       CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CWORK1(1)) 

! calculate the non local contribution in reciprocal space
       CALL VNLACC(NONL_S, W1, CDIJ, CQIJ, ISP, EVALUE, CH)
       CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), .TRUE., &
            W1%WDES1%DATAKE(1,1), 0.0_q, W1%CPTWFP(1), CH(1))

    ENDIF
    RETURN
  END SUBROUTINE HAMILT

!************************* SUBROUTINE HAMILTMU *************************
!
! this subroutine calculates the H acting onto a set of wavefuntions
! the  wavefunction must be given in reciprocal space C and real
! space CR
! CH contains the result
!   H- EVALUE S |phi>
! where S is the overlap matrix S
!
!***********************************************************************

  SUBROUTINE HAMILTMU( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, ISP, WRESULT)
    USE mpimy
    USE mgrid
    
    USE nonl_high
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S    
    REAL(q)    EVALUE(:)               ! eigenvalues
    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)      SV(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT
    
! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
    INTEGER NP
    
! calculate the local contribution (result in CWORK1)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN
          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO

! non-local contribution in real-space

    IF (NONLR_S%LREAL) THEN
       CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN
             CALL KINHAMIL( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO


    ELSE

! calculate the non local contribution in reciprocal space

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN
             CALL VNLACC(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CPTWFP(:,NP))
             CALL KINHAMIL( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO
    ENDIF

    RETURN
  END SUBROUTINE HAMILTMU


!************************* SUBROUTINE HAMILTMU *************************
!
! this subroutine calculates the H acting onto a set of wavefuntions
! the  wavefunction must be given in reciprocal space C and real
! space CR
! CH contains the result
!   H- EVALUE S |phi>
! where S is the overlap matrix S
! this version includes a vector potential A
!
!***********************************************************************

  SUBROUTINE HAMILTMU_VEC( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, AVEC, ISP, WRESULT)
    USE mpimy
    USE mgrid
    
    USE nonl_high
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S    
    REAL(q)    EVALUE(:)               ! eigenvalues
    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)      SV(:,:)
    COMPLEX(q)      AVEC(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT
    
! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
    COMPLEX(q) :: CWORK2(WDES1%NRPLWV)
    COMPLEX(q) :: CWORK3(WDES1%GRID%MPLWV)


    
    INTEGER NP
    
! calculate the local contribution (result in CWORK1)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN
          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO

! non-local contribution in real-space

    IF (NONLR_S%LREAL) THEN
       CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE,CWORK1)

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN
             CALL KINHAMIL_VEC( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), AVEC(1,1), &
                  CWORK2, CWORK3, EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO


    ELSE

! calculate the non local contribution in reciprocal space

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN
             CALL VNLACC(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CPTWFP(:,NP))
             CALL KINHAMIL_VEC( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), AVEC(1,1), &
                  CWORK2, CWORK3, EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO


    ENDIF

    RETURN
  END SUBROUTINE HAMILTMU_VEC


!************************* SUBROUTINE HAMILTMU_TAU *********************
!
! this subroutine calculates the H acting onto a set of wavefuntions
! for a kinetic energy functional
!
!***********************************************************************

  SUBROUTINE HAMILTMU_TAU( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, LATT_CUR, MU, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE lattice
    USE nonl_high
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (latt)        LATT_CUR    
    REAL(q)    EVALUE(:)               ! eigenvalues
    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)      SV(:,:)
    COMPLEX(q)      MU(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT
    
! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
    COMPLEX(q) :: CWORK2(WDES1%NRPLWV)
    COMPLEX(q) :: CWORK3(WDES1%GRID%MPLWV)
    
    INTEGER NP

! calculate the local contribution (result in CWORK1)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN
          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO

! non-local contribution in real-space

    IF (NONLR_S%LREAL) THEN
       CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE,CWORK1)

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN
             CALL KINHAMIL_TAU( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
                  CWORK2, CWORK3, EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO

    ELSE

! calculate the non local contribution in reciprocal space

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN
             CALL VNLACC(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CPTWFP(:,NP))
             CALL KINHAMIL_TAU( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
                  CWORK2, CWORK3, EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO

    ENDIF

    RETURN
  END SUBROUTINE HAMILTMU_TAU


!************************* SUBROUTINE HAMILTMU_C ***********************
!
! identical to previous subroutine but with a complex eigenvalue EVALUE
!
!***********************************************************************

  SUBROUTINE HAMILTMU_C( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, ISP, WRESULT)
    USE mpimy
    USE mgrid
    
    USE nonl_high
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S    
    COMPLEX(q)    EVALUE(:)               ! eigenvalues
    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)      SV(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT
    
! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
    INTEGER NP
    
! calculate the local contribution (result in CWORK1)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN
          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO

! non-local contribution in real-space

    IF (NONLR_S%LREAL) THEN
       CALL RACCMU_C_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE,CWORK1)

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN
             CALL KINHAMIL_C( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO

    ELSE

! calculate the non local contribution in reciprocal space

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN
             CALL VNLACC_C(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CPTWFP(:,NP))
             CALL KINHAMIL_C( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO
    ENDIF

    RETURN
  END SUBROUTINE HAMILTMU_C

!************************* SUBROUTINE HAMILT_LOCAL *********************
!
! this subroutine calculates the local and kinetic energy part of H
! acting onto a wavefuntion
! the  wavefunction must be given in reciprocal space C and real
! space CR
! LADD=.TRUE.
!  CH = CH + SV* W1 + T * W1
! LADD=.FALSE.
!  CH = SV* W1 + T * W1
! LKIN determines whether the kinetic energy contribution is added or not
! the default is add the kinetic energy contribution
!
!***********************************************************************

  SUBROUTINE HAMILT_LOCAL(W1, SV, ISP, CH, LADD, LKIN)
    USE mgrid
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    COMPLEX(q)   SV(:,:)
    INTEGER :: ISP
    COMPLEX(q) CH(:)
    LOGICAL LADD
    LOGICAL, OPTIONAL :: LKIN
! local variables
    COMPLEX(q) :: CWORK1(W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS)

    CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CWORK1(1)) 
    IF (.NOT. PRESENT(LKIN)) THEN
       CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
            W1%WDES1%DATAKE(1,1), 0.0_q, W1%CPTWFP(1), CH(1))
    ELSE
       IF (LKIN) THEN
          CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
               W1%WDES1%DATAKE(1,1), 0.0_q, W1%CPTWFP(1), CH(1))
       ELSE
          CALL FFTHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
               W1%WDES1%DATAKE(1,1), 0.0_q, W1%CPTWFP(1), CH(1))
       ENDIF
    ENDIF

  END SUBROUTINE HAMILT_LOCAL


  SUBROUTINE HAMILT_LOCAL_TAU(W1, SV, LATT_CUR, MU, ISP, CH, LADD, LKIN)
    USE mgrid
    USE lattice
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    TYPE (latt)        LATT_CUR
    COMPLEX(q)   SV(:,:)
    COMPLEX(q)   MU(:,:)
    INTEGER ISP
    COMPLEX(q) CH(:)
    LOGICAL LADD
    LOGICAL, OPTIONAL :: LKIN
! local variables
    COMPLEX(q) :: CVR(W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK1(W1%WDES1%NRPLWV)
    COMPLEX(q) :: CWORK2(W1%WDES1%GRID%MPLWV)
    LOGICAL :: LKINETIC=.TRUE.

    IF (PRESENT(LKIN)) LKINETIC=LKIN

    CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CVR(1))

    CALL KINHAMIL_TAU( W1%WDES1, W1%WDES1%GRID, CVR(1), LADD, LKINETIC, W1%WDES1%DATAKE(1,1), &
   &   W1%WDES1%IGX(1), W1%WDES1%IGY(1), W1%WDES1%IGZ(1), W1%WDES1%VKPT(1), LATT_CUR, MU(1,ISP), CWORK1, CWORK2, &
   &   0.0_q, W1%CPTWFP(1), CH(1))

  END SUBROUTINE HAMILT_LOCAL_TAU


!************************* SUBROUTINE ECCP_NL_ALL   ********************
!
! this subroutine calculates the non local contribution to the
! expectation value of sum_ij <c2|p_j> D_ij - e Q_ij <p_i|c1>
! where c and cp are two wavefunctions
!
!***********************************************************************

  SUBROUTINE ECCP_NL_ALL(WDES1,W1,W2,CDIJ,CQIJ,EVALUE,CE)
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    COMPLEX(q) CDIJ(:,:,:,:)
    COMPLEX(q) CQIJ(:,:,:,:)
    REAL(q) EVALUE
! return
    COMPLEX(q) CE
! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NI, NIS, NT, LMMAXC

    CE=0

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1

          NPRO =ISPINOR *(WDES1%NPRO/2)
          NPRO_=ISPINOR_*(WDES1%NPRO/2)

          NIS =1
          DO NT=1,WDES1%NTYP
             LMMAXC=WDES1%LMMAX(NT)
             IF (LMMAXC==0) GOTO 310
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                CALL ECCP_NL_(SIZE(CDIJ,1),LMMAXC, &
                     CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),EVALUE, &
                     W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CE)
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
310          NIS = NIS+WDES1%NITYP(NT)
          ENDDO
       ENDDO
    ENDDO spinor

    CALL M_sum_z(WDES1%COMM_INB, CE, 1)

    RETURN
  END SUBROUTINE ECCP_NL_ALL


!************************* SUBROUTINE ECCP_NL_ALL_ION ******************
!
! this subroutine calculates the non local contribution to the
! expectation value of sum_ij <c2|p_j> D_ij - e Q_ij <p_i|c1>
! where c and cp are two wavefunctions
!
!***********************************************************************

  SUBROUTINE ECCP_NL_ALL_ION(WDES1,W1,W2,CDIJ,CQIJ,EVALUE,CE)
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    COMPLEX(q) CDIJ(:,:,:,:)
    COMPLEX(q) CQIJ(:,:,:,:)
    REAL(q) EVALUE
! return
    COMPLEX(q) CE(*)
! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NI, NIS, NT, LMMAXC

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1

          NPRO =ISPINOR *(WDES1%NPRO/2)
          NPRO_=ISPINOR_*(WDES1%NPRO/2)

          NIS =1
          DO NT=1,WDES1%NTYP
             LMMAXC=WDES1%LMMAX(NT)
             IF (LMMAXC==0) GOTO 310
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                CALL ECCP_NL_(SIZE(CDIJ,1),LMMAXC, &
                     CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),EVALUE, &
                     W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CE(NI))
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
310          NIS = NIS+WDES1%NITYP(NT)
          ENDDO
       ENDDO
    ENDDO spinor

    RETURN
  END SUBROUTINE ECCP_NL_ALL_ION

END MODULE hamil

!************************* SUBROUTINE VHAMIL  **************************
!
! low level F77 to calculate the product of the local potential
! times a wavefunction stored in CR and returns the result in CVR
!   CVR= SV * CR
! this low level routine should be called only from within
! hamil.F
!
!***********************************************************************

   SUBROUTINE VHAMIL(WDES1,GRID,SV,CR,CVR)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1

     COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS),CVR(GRID%MPLWV*WDES1%NRSPINORS)
! local variables
     REAL(q) RINPLW
     INTEGER ISPINOR, ISPINOR_, M, MM, MM_

     RINPLW=1._q/GRID%NPLWV

     IF (WDES1%NRSPINORS==1) THEN
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           CVR(M)= SV(M,1) *CR(M)*RINPLW
        ENDDO
     ELSE
        CVR(1:GRID%MPLWV*2)=0
        DO ISPINOR =0,1
           DO ISPINOR_=0,1
!DIR$ IVDEP
!OCL NOVREL
              DO M=1,GRID%RL%NP
                 MM =M+ISPINOR *GRID%MPLWV
                 MM_=M+ISPINOR_*GRID%MPLWV
                 CVR(MM)= CVR(MM)+ SV(M,1+ISPINOR_+2*ISPINOR) *CR(MM_)*RINPLW
              ENDDO
           ENDDO
        ENDDO
     ENDIF
   END SUBROUTINE VHAMIL

   SUBROUTINE VHAMIL_TRACE(WDES1, GRID, SV, CR, CVR, WEIGHT)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1

     COMPLEX(q)   SV(GRID%MPLWV) ! local potential
     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS),CVR(GRID%MPLWV*WDES1%NRSPINORS)
     REAL(q)    :: WEIGHT
! local variables
     INTEGER ISPINOR, M, MM

!DIR$ IVDEP
!OCL NOVREL
     DO M=1,GRID%RL%NP
        CVR(M)= CVR(M)+SV(M) *CR(M)*WEIGHT
     ENDDO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           MM =M+GRID%MPLWV
           CVR(MM)= CVR(MM)+ SV(M) *CR(MM)*WEIGHT
        ENDDO
     ENDIF
   END SUBROUTINE VHAMIL_TRACE

!************************* SUBROUTINE PW_CHARGE ************************
!
! low level F77 to calculate charge density from two plane wave states
! and add the weighted spinor like density to an array
! this low level routine should be called only from within
! hamil.F
!  CHARGE = CR1*CONJG(CR2) * WEIGHT + CHARGE
!
!***********************************************************************

!
! this version yields an COMPLEX(q) array spinor density
!

   SUBROUTINE PW_CHARGE(WDES1,  CHARGE, NDIM, CR1, CR2, WEIGHT)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     INTEGER NDIM
     COMPLEX(q)   CHARGE(NDIM,WDES1%NRSPINORS*WDES1%NRSPINORS)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     REAL(q) :: WEIGHT
! local variables
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_

     IF (WDES1%NRSPINORS==1) THEN
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%GRID%RL%NP
           CHARGE(M,1)=CHARGE(M,1)+ & 
                CONJG(CR2(M))*CR1(M)*WEIGHT
        ENDDO
     ELSE
        spinor: DO ISPINOR=0,1
           DO ISPINOR_=0,1
!DIR$ IVDEP
!OCL NOVREL
              DO M=1,WDES1%GRID%RL%NP
                 MM =M+ISPINOR *WDES1%GRID%MPLWV
                 MM_=M+ISPINOR_*WDES1%GRID%MPLWV
                 CHARGE(M,1+ISPINOR_+2*ISPINOR)=CHARGE(M,1+ISPINOR_+2*ISPINOR)+ & 
                      CONJG(CR2(MM_))*CR1(MM)*WEIGHT
              ENDDO
           ENDDO
        ENDDO spinor
     ENDIF

   END SUBROUTINE PW_CHARGE

!
! this version yields a COMPLEX(q) charge array, spinor density
! (used by exchange routines)

   SUBROUTINE PW_CHARGE_CMPLX(WDES1,  CHARGE, NDIM, CR1, CR2 )
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     INTEGER NDIM
     COMPLEX(q)   CHARGE(NDIM,WDES1%NRSPINORS*WDES1%NRSPINORS)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     REAL(q) :: WEIGHT
! local variables
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_

     IF (WDES1%NRSPINORS==1) THEN
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%GRID%RL%NP
           CHARGE(M,1)=CONJG(CR2(M))*CR1(M)
        ENDDO
     ELSE
        CHARGE=0
        spinor: DO ISPINOR=0,1
           DO ISPINOR_=0,1
!DIR$ IVDEP
!OCL NOVREL
              DO M=1,WDES1%GRID%RL%NP
                 MM =M+ISPINOR *WDES1%GRID%MPLWV
                 MM_=M+ISPINOR_*WDES1%GRID%MPLWV
                 CHARGE(M,1+ISPINOR_+2*ISPINOR)=CHARGE(M,1+ISPINOR_+2*ISPINOR)+ & 
                      CONJG(CR2(MM_))*CR1(MM)
              ENDDO
           ENDDO
        ENDDO spinor
     ENDIF

   END SUBROUTINE PW_CHARGE_CMPLX

!
! this version traces of spinors to yield only the density
! results array is COMPLEX(q) (used by exchange routines)
   SUBROUTINE PW_CHARGE_TRACE(WDES1,  CHARGE, CR1, CR2)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     COMPLEX(q)   CHARGE(WDES1%GRID%MPLWV)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
! local variables
     REAL(q) RINPLW
     INTEGER M,MM

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%GRID%RL%NP
         CHARGE(M)=CONJG(CR2(M))*CR1(M)
      ENDDO
      IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%GRID%RL%NP
            MM =M+WDES1%GRID%MPLWV
            CHARGE(M)=CHARGE(M)+CONJG(CR2(MM))*CR1(MM)
         ENDDO
      ENDIF

   END SUBROUTINE PW_CHARGE_TRACE

!
! this version traces of spinors to yield only the density
! result array as well as orbital input is COMPLEX(q)

   SUBROUTINE PW_CHARGE_TRACE_GDEF(WDES1,  CHARGE, CR1, CR2)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     COMPLEX(q) :: CHARGE(WDES1%GRID%MPLWV)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
! local variables
     REAL(q) RINPLW
     INTEGER M,MM

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%GRID%RL%NP
         CHARGE(M)=CONJG(CR2(M))*CR1(M)
      ENDDO
      IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%GRID%RL%NP
            MM =M+WDES1%GRID%MPLWV
            CHARGE(M)=CHARGE(M)+CONJG(CR2(MM))*CR1(MM)
         ENDDO
      ENDIF

   END SUBROUTINE PW_CHARGE_TRACE_GDEF

   SUBROUTINE PW_CHARGE_TRACE_NO_CONJG(WDES1,  CHARGE, CR1, CR2)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     COMPLEX(q)   CHARGE(WDES1%GRID%MPLWV)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
! local variables
     REAL(q) RINPLW
     INTEGER M,MM

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%GRID%RL%NP
         CHARGE(M)=CR2(M)*CR1(M)
      ENDDO
      IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%GRID%RL%NP
            MM =M+WDES1%GRID%MPLWV
            CHARGE(M)=CHARGE(M)+CR2(MM)*CR1(MM)
         ENDDO
      ENDIF

   END SUBROUTINE PW_CHARGE_TRACE_NO_CONJG


!************************* SUBROUTINE KINHAMIL  ************************
!
! low level F77 to calculate
! the kinetic energy part of a Hamiltonian
! possibly substract the wavefunction times an eigenvalue
! and adds a component after FFT to real space
! this low level routine should be called only from within
! hamil.F
! LADD=.TRUE.
!  CH= CH + FFT(CVR) + DATAKE*CPTWFP - EVALUE*CPTWFP
! LADD=.FALSE.
!  CH= FFT(CVR) + DATAKE*CPTWFP - EVALUE*CPTWFP
!
! FFTHAMIL exclude the kinetic energy contribution
!
!***********************************************************************

   SUBROUTINE KINHAMIL( WDES1, GRID, CVR,  LADD, DATAKE, EVALUE, CPTWFP, CH)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) ! usually potential times wavefunction
     LOGICAL    :: LADD                            ! if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: EVALUE                          ! subtract EVALUE*wavefunction
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                ! wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                ! result
! local
     INTEGER ISPINOR, M, MM

     DO ISPINOR=0,WDES1%NRSPINORS-1
        CALL FFTEXT_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE+CPTWFP(MM)*DATAKE(M,ISPINOR+1)
        ENDDO
     ENDDO
   END SUBROUTINE KINHAMIL


   SUBROUTINE FFTHAMIL( WDES1, GRID, CVR,  LADD, DATAKE, EVALUE, CPTWFP, CH)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) ! usually potential times wavefunction
     LOGICAL    :: LADD                            ! if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: EVALUE                          ! subtract EVALUE*wavefunction
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                ! wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                ! result
! local
     INTEGER ISPINOR, M, MM

     DO ISPINOR=0,WDES1%NRSPINORS-1
        CALL FFTEXT_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE
        ENDDO
     ENDDO
   END SUBROUTINE FFTHAMIL

!************************* SUBROUTINE KINHAMIL_VEC *********************
!
! low level F77 to calculate
! similar to KINHAMIL_VEC
! but adds the action of the vector potential on the wavefunction
!
! -(e A p) e(iGr)
!  |e|A p  e(iGr) = (|e| A hbar) -i nabla e(iGr)
! (|e| A hbar) G e(iGr)
!
! the coefficients |e| A hbar are included in the stored AVEC
!
!
!***********************************************************************

   SUBROUTINE KINHAMIL_VEC( WDES1, GRID, CVR,  LADD, & 
                  DATAKE, IGX, IGY, IGZ, VKPT, AVEC, CWORK1, CWORK2, EVALUE, CPTWFP, CH)
     USE prec
     USE mgrid
     USE wave
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) ! usually potential times wavefunction
     LOGICAL    :: LADD                            ! if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: VKPT(3)                         ! k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                ! component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                ! component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                ! component of each G vector in direction of third rec lattice vector
     COMPLEX(q)      :: AVEC(GRID%MPLWV,3)     ! magnetic vector potential
     COMPLEX(q) :: CWORK1(WDES1%NRPLWV)
     COMPLEX(q) :: CWORK2(GRID%MPLWV)
     REAL(q)    :: EVALUE                          ! subtract EVALUE*wavefunction
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                ! wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                ! result
     COMPLEX(q), PARAMETER :: CI=(0._q,1._q)
! local
     INTEGER ISPINOR, M, MM
     REAL(q)    :: RINPLW

     RINPLW=2*PI/GRID%NPLWV
      
     DO ISPINOR=0,WDES1%NRSPINORS-1
! component along first axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGX(M)+VKPT(1))
        ENDDO
! FFT to real space
        CALL FFTWAV_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,1) *CWORK2(M) *RINPLW
        ENDDO

! component along second axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGY(M)+VKPT(2))
        ENDDO
! FFT to real space
        CALL FFTWAV_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,2) *CWORK2(M)*RINPLW
        ENDDO

! component along third axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGZ(M)+VKPT(3))
        ENDDO
! FFT to real space
        CALL FFTWAV_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,3) *CWORK2(M)*RINPLW
        ENDDO

        
        CALL FFTEXT_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE+CPTWFP(MM)*DATAKE(M,ISPINOR+1)
        ENDDO
     ENDDO
   END SUBROUTINE KINHAMIL_VEC


!************************* SUBROUTINE KINHAMIL_TAU *********************
!
! low level F77 to calculate
! similar to KINHAMIL_TAU
! but adds the action of   - nabla mu  nabla
!
! - nabla mu nable e(iGr)
!
!***********************************************************************

   SUBROUTINE KINHAMIL_TAU( WDES1, GRID, CVR, LADD, LKINETIC, DATAKE, & 
                  IGX, IGY, IGZ, VKPT, LATT_CUR, MU, CWORK1, CWORK2, EVALUE, CPTWFP, CH)
     USE prec
     USE mgrid
     USE wave
     USE lattice
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     TYPE (latt)        LATT_CUR
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) ! usually potential times wavefunction
     LOGICAL    :: LADD                            ! if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: VKPT(3)                         ! k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                ! component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                ! component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                ! component of each G vector in direction of third rec lattice vector
     COMPLEX(q)      :: MU(GRID%MPLWV*WDES1%NRSPINORS*WDES1%NRSPINORS) ! \mu = d f_xc / d \tau
     REAL(q)    :: EVALUE                          ! subtract EVALUE*wavefunction
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                ! wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                ! result
     LOGICAL    :: LKINETIC
! local
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_
     REAL(q)    :: RINPLW
     REAL(q)    :: GX,GY,GZ
     COMPLEX(q) :: CWORK1(WDES1%NRPLWV)
     COMPLEX(q) :: CWORK2(GRID%MPLWV)
     COMPLEX(q), PARAMETER :: CI=(0._q,1._q)
     REAL(q)    :: DCMU                            ! for diagnostics: <\psi| \nabla | \mu \nabla \psi >

     DCMU=0

     RINPLW=1.0_q/GRID%NPLWV
     DO ISPINOR=0,WDES1%NRSPINORS-1
        CALL FFTEXT_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
        IF (LKINETIC) THEN
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              MM=M+ISPINOR*WDES1%NGVECTOR
              CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE+CPTWFP(MM)*DATAKE(M,ISPINOR+1)
           ENDDO
        ELSE
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              MM=M+ISPINOR*WDES1%NGVECTOR
              CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE
           ENDDO
        ENDIF

! component along x-axis
        CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M)=CPTWFP(MM)*GX*CITPI
           ENDDO
! FFT to real space
           CALL FFTWAV_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*GRID%MPLWV
              CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M) *RINPLW
           ENDDO
        ENDDO

        CALL FFTEXT_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK1(M)*GX*CITPI
           DCMU=DCMU-CWORK1(M)*GX*CITPI*CONJG(CPTWFP(MM))
        ENDDO

! component along y-axis
        CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M)=CPTWFP(MM)*GY*CITPI
           ENDDO
! FFT to real space
           CALL FFTWAV_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*GRID%MPLWV
              CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M)*RINPLW
           ENDDO
        ENDDO

        CALL FFTEXT_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK1(M)*GY*CITPI
           DCMU=DCMU-CWORK1(M)*GY*CITPI*CONJG(CPTWFP(MM))
        ENDDO

! component along z-axis
        CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M)=CPTWFP(MM)*GZ*CITPI
           ENDDO
! FFT to real space
           CALL FFTWAV_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*GRID%MPLWV
              CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M)*RINPLW
           ENDDO
        ENDDO

        CALL FFTEXT_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK1(M)*GZ*CITPI
           DCMU=DCMU-CWORK1(M)*GZ*CITPI*CONJG(CPTWFP(MM))
        ENDDO
     ENDDO
   END SUBROUTINE KINHAMIL_TAU


!************************* SUBROUTINE PSCURRENT ************************
!
! gradient psi contribution to current: Im[psi^*(r) nabla psi(r)]
! density contribution to current: psi^*(r) * psi(r)
! for a supplied weighted k-point and psi
! with psi(r) = e(iGr)
! we get contributions proportional to (iG) e(iGr)
!
!***********************************************************************

   SUBROUTINE PSCURRENT( WDES1, GRID, CR, & 
                  IGX, IGY, IGZ, VKPT, WEIGHT, AVEC, CWORK1, CWORK2, CPTWFP, J, C)
     USE prec
     USE mgrid
     USE wave
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS)  ! wavefunction in real space
     REAL(q)    :: VKPT(3)                         ! k-point in reciprocal lattice
     REAL(q)    :: WEIGHT
     INTEGER    :: IGX(WDES1%NGDIM)                ! component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                ! component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                ! component of each G vector in direction of third rec lattice vector
     COMPLEX(q)      :: AVEC(GRID%MPLWV,3)     ! magnetic vector potential
     COMPLEX(q) :: CWORK1(WDES1%NRPLWV)
     COMPLEX(q) :: CWORK2(GRID%MPLWV)
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                ! wavefunction
     COMPLEX(q)      :: J(GRID%RL%NP,3)                 ! three component vector storing current density in real space
     COMPLEX(q)      :: C(GRID%RL%NP)                   ! charge density
! local
     INTEGER ISPINOR, M, MM

     DO ISPINOR=0,WDES1%NRSPINORS-1
! component along first axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGX(M)+VKPT(1))*CITPI
        ENDDO
! FFT to real space
        CALL FFTWAV_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           J(M,1)=J(M,1)+AIMAG(CWORK2(M)*CONJG(CR(MM)))*WEIGHT
        ENDDO

! component along second axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGY(M)+VKPT(2))*CITPI
        ENDDO
! FFT to real space
        CALL FFTWAV_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           J(M,2)=J(M,2)+AIMAG(CWORK2(M)*CONJG(CR(MM)))*WEIGHT
        ENDDO

! component along third axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGZ(M)+VKPT(3))*CITPI
        ENDDO
! FFT to real space
        CALL FFTWAV_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           J(M,3)=J(M,3)+AIMAG(CWORK2(M)*CONJG(CR(MM)))*WEIGHT
        ENDDO

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           C(M)=C(M)+CR(MM)*CONJG(CR(MM))*WEIGHT
        ENDDO
     ENDDO
   END SUBROUTINE PSCURRENT



!************************* SUBROUTINE KINHAMIL_C ***********************
!
! low level routine for complex eigenvalues
!
!***********************************************************************

   SUBROUTINE KINHAMIL_C( WDES1, GRID, CVR,  LADD, DATAKE, EVALUE, CPTWFP, CH)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) ! usually potential times wavefunction
     LOGICAL    :: LADD                            ! if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     COMPLEX(q) :: EVALUE                          ! subtract EVALUE*wavefunction
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                ! wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                ! result
! local
     INTEGER ISPINOR, M, MM

     DO ISPINOR=0,WDES1%NRSPINORS-1
        CALL FFTEXT_MPI(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE+CPTWFP(MM)*DATAKE(M,ISPINOR+1)
        ENDDO
     ENDDO
   END SUBROUTINE KINHAMIL_C

!***********************************************************************
!
! low level F77 version
! determine norm of a wavefunction and norm with a metric
! (if supplied)
!
!***********************************************************************

  SUBROUTINE PW_NORM_WITH_METRIC(WDES1, C, FNORM, FMETRIC, METRIC)
    USE wave
    IMPLICIT NONE
    REAL(q) FNORM
    TYPE (wavedes1)    WDES1
    COMPLEX(q) :: C(WDES1%NPL)
    REAL(q), OPTIONAL :: FMETRIC
    REAL(q), OPTIONAL :: METRIC(WDES1%NPL)

    INTEGER ISPINOR, M, MM

    FNORM=0
    IF (PRESENT (METRIC).AND. PRESENT(FMETRIC)) THEN
       FMETRIC=0
       DO ISPINOR=0,WDES1%NRSPINORS-1
          DO M=1,WDES1%NGVECTOR
             MM=M+ISPINOR*WDES1%NGVECTOR
             FNORM =FNORM+C(MM)*CONJG(C(MM))
             FMETRIC =FMETRIC+C(MM)*CONJG(C(MM))*METRIC(MM)
          ENDDO
       ENDDO
       CALL M_sum_s(WDES1%COMM_INB, 2, FNORM, FMETRIC, 0._q, 0._q)
    ELSE
       DO ISPINOR=0,WDES1%NRSPINORS-1
          DO M=1,WDES1%NGVECTOR
             MM=M+ISPINOR*WDES1%NGVECTOR
             FNORM =FNORM+C(MM)*CONJG(C(MM))
          ENDDO
       ENDDO
       CALL M_sum_s(WDES1%COMM_INB, 1, FNORM, 0._q, 0._q, 0._q)
    ENDIF

  END SUBROUTINE PW_NORM_WITH_METRIC

!***********************************************************************
!
! low level F77 routine
! truncate high frequency components if LDELAY is set
! or if use_eini_without_ldelay is defined
! first version is scheduled for removal
!
!***********************************************************************


  SUBROUTINE TRUNCATE_HIGH_FREQUENCY( WDES1, C, LDELAY, ENINI)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavedes1)    WDES1
    COMPLEX(q) C(WDES1%NPL)
    LOGICAL LDELAY     ! usually truncation only during delay phase (LDELAY set)
    REAL(q) ENINI      ! cutoff at which the wavefunction is truncated
! local
    INTEGER ISPINOR

    IF (LDELAY.OR.WDES1%LSPIRAL) THEN
       DO ISPINOR=1,WDES1%NRSPINORS
          CALL TRUNCATE_HIGH_FREQUENCY_ONE(WDES1%NGVECTOR, &
               C(1+(ISPINOR-1)*WDES1%NGVECTOR), WDES1%DATAKE(1,ISPINOR), ENINI)
       ENDDO
    END IF
  END SUBROUTINE TRUNCATE_HIGH_FREQUENCY


  SUBROUTINE TRUNCATE_HIGH_FREQUENCY_ONE(NP, CPTWFP, DATAKE, ENINI)
    USE prec
    IMPLICIT NONE
    INTEGER NP
    COMPLEX(q) :: CPTWFP(NP)
    REAL(q)    :: DATAKE(NP)
    REAL(q)    :: ENINI
! local
    INTEGER N
    
    DO N=1,NP
       IF (DATAKE(N)>ENINI) CPTWFP(N)=0
    ENDDO
    
  END SUBROUTINE TRUNCATE_HIGH_FREQUENCY_ONE


!************************* SUBROUTINE ECCP_NL_   ***********************
!
! this subroutine calculates the expectation value of
!    <c| D- evalue Q |cp>
! where c and cp are two wavefunctions; non local part only
! for (1._q,0._q) ion only
!
!***********************************************************************

  SUBROUTINE ECCP_NL_(LMDIM,LMMAXC,CDIJ,CQIJ,EVALUE,CPROJ1,CPROJ2,CNL)
    USE prec
    IMPLICIT NONE
    COMPLEX(q)  CNL
    INTEGER LMMAXC, LMDIM
    COMPLEX(q) CDIJ(LMDIM,LMDIM),CQIJ(LMDIM,LMDIM)
    REAL(q) EVALUE
    COMPLEX(q) CPROJ1(LMMAXC),CPROJ2(LMMAXC)

! local
    INTEGER L, LP

    DO L=1,LMMAXC
       DO LP=1,LMMAXC
          CNL=CNL+(CDIJ(LP,L)-EVALUE*CQIJ(LP,L))*CPROJ1(LP)*CONJG(CPROJ2(L))
       ENDDO
    ENDDO
  END SUBROUTINE ECCP_NL_
