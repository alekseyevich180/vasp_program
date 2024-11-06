# 1 "stufak.F"
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

# 2 "stufak.F" 2 
!************************ SUBROUTINE STUFAK ****************************
! RCS:  $Id: stufak.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
! this subroutine calculates the structure factor on the grid of
! reciprocal lattice vectors
! cstrf(g) = sum over ions (-exp(ig.r)) where r is the position of the
! ion
!***********************************************************************

      SUBROUTINE STUFAK(GRIDC,T_INFO,CSTRF)
      USE prec

      USE mpimy
      USE mgrid
      USE poscar
      USE constant
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO

      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

! loop over all types of atoms
      NIS=1
      typ: DO NT=1,T_INFO%NTYP
      CALL STUFAK_ONE(GRIDC,T_INFO%NITYP(NT),T_INFO%POSION(1,NIS),T_INFO%VCA(NT),CSTRF(1,NT))
      NIS=NIS+T_INFO%NITYP(NT)

      ENDDO typ

      RETURN
      END


!************************ SUBROUTINE STUFAK_ONE ************************
!
! this subroutine calculates the structure factor on the grid of
! for 1._q species (i.e. partial structure factor)
! cstrf(g) = sum over ions (-exp(ig.r)) where r is the position of the
! ion
!
!***********************************************************************

      SUBROUTINE STUFAK_ONE(GRIDC,NIONS,POSION,VCA,CSTRF)
      USE prec

      USE mpimy
      USE mgrid
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      REAL(q) POSION(3,NIONS)
      REAL(q) VCA

      TYPE (grid_3d)     GRIDC

      COMPLEX(q) CSTRF(GRIDC%RC%NP)
      CSTRF=0
!=======================================================================
! calculate partial structur-factor
!=======================================================================
      ion: DO NI=1,NIONS
!=======================================================================
! loop over all grid points
!=======================================================================
# 89

!-----------------------------------------------------------------------
! more envolved version which is faster on most (scalar) machines
! and includes support for parallel machines
!-----------------------------------------------------------------------
      CX =EXP(-CITPI*POSION(1,NI))
      G1 =POSION(1,NI)*(-(GRIDC%NGX/2-1))

      col: DO NC=1,GRIDC%RC%NCOL
        N=(NC-1)*GRIDC%RC%NROW+1

        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)
        G2=POSION(2,NI)*GRIDC%LPCTY(N2)
        G3=POSION(3,NI)*GRIDC%LPCTZ(N3)
        CE=EXP(-CITPI*(G3+G2+G1))*VCA
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
        DO N1P=0,GRIDC%RC%NROW-1

          N1=MOD(N1P+(-(GRIDC%NGX/2-1))+GRIDC%NGX,GRIDC%NGX)
          CSTRF(N+N1)=CSTRF(N+N1)+CE
          CE=CE*CX
        ENDDO
      ENDDO col

!-----------------------------------------------------------------------
!  next ion + next type
!-----------------------------------------------------------------------
      ENDDO ion

      RETURN
      END


