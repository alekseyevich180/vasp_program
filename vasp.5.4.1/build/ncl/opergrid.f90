# 1 "opergrid.F"
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

# 2 "opergrid.F" 2 

!***********************************************************************
! RCS:  $Id: opergrid.F,v 1.3 2003/06/27 13:22:21 kresse Exp kresse $
!
!  sum all points on the grid and dump onto screen
!
!***********************************************************************

      SUBROUTINE SUMGRD(CSRC,GRID)
      USE prec
      USE mgrid

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      DIMENSION CSRC(GRID%NGX_rd,GRID%NGY_rd,GRID%NGZ_rd)

      CSUM=0

      DO NZ=1,GRID%NGZ_rd
      DO NY=1,GRID%NGY_rd
      DO NX=1,GRID%NGX_rd
        CSUM=CSUM+CSRC(NX,NY,NZ)
      ENDDO
      ENDDO
      ENDDO
      WRITE(*,*) CSUM
      RETURN
      END

!***********************************************************************
!
!  sum all points on the real grid and dump onto screen
!
!***********************************************************************

      SUBROUTINE SUMRGR(CSRC, GRID)
      USE prec
      USE mgrid

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CSRC(GRID%NGX,GRID%NGY,GRID%NGZ)
      COMPLEX(q) CSUM
      CSUM=0

      DO NZ=1,GRID%NGZ
      DO NY=1,GRID%NGY
      DO NX=1,GRID%NGX
        CSUM=CSUM+CSRC(NX,NY,NZ)
      ENDDO
      ENDDO
      ENDDO
      WRITE(*,*) CSUM
      RETURN
      END

!***********************SUBROUTINE DUMP1******************************
!
!  This debug-routine dumps a compressed array for (1._q,0._q) K-point
!
!**********************************************************************

      SUBROUTINE DUMP1(NK,NRPLWV,IGX,IGY,IGZ,NPLWKP,C1)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION IGX(NRPLWV,1),IGY(NRPLWV,1),IGZ(NRPLWV,1)
      DIMENSION NPLWKP(1)
      REAL(q)      C1(1)

      DO 50 K=1,NPLWKP(NK)
        WRITE(*,10)K,IGX(K,NK),IGY(K,NK),IGZ(K,NK), &
     & REAL( C1(K) ,KIND=q)
   50 CONTINUE
   10 FORMAT(I3,2X,3I2,6E14.7)
      WRITE(*,*)
      RETURN
      END SUBROUTINE

!***********************SUBROUTINE DUMPA ******************************
!
!  This debug-routine dumps a uncompressed array for (1._q,0._q) K-Point
!
!**********************************************************************

      SUBROUTINE DUMPA(NK,NGX,NGY,NGZ,CA)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q)   CA(NGX,NGY,NGZ)

      NDUMP=MIN(NGX,NGY,NGZ,10)

      DO NZ=1,NDUMP
      WRITE(*,*)'NZ=',NZ
      DO NY=1,NDUMP
        WRITE(*,'(I2,24E10.3)') NY,(ABS(CA(NX,NY,NZ)),NX=1,NDUMP)
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

!*********************************************************************
!
! routine to dump a (NxN) Matrix
!
!*********************************************************************

      SUBROUTINE DDUMP(C,NBANDS)
      USE prec
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      COMPLEX(q) C(NBANDS,NBANDS)

      NBA=MIN(10,NBANDS)

      DO 830 N1=NBANDS-NBA+1,NBANDS
      WRITE(*,1)N1, &
     &         (REAL( C(N1,N2) ,KIND=q) ,N2=NBANDS-NBA+1,NBANDS)
 830  CONTINUE

      WRITE(*,*)
      DO 840 N1=NBANDS-NBA+1,NBANDS
      WRITE(*,2)N1, &
     &         (AIMAG(C(N1,N2)),N2=NBANDS-NBA+1,NBANDS)
 840  CONTINUE

      WRITE(*,*)

    1 FORMAT(1I2,3X,20F9.5)
    2 FORMAT(1I2,3X,20F9.5)
      RETURN
      END SUBROUTINE

