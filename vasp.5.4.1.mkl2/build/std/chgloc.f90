# 1 "chgloc.F"
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

# 2 "chgloc.F" 2 

!************************ SUBROUTINE CHGLOC ****************************
! RCS:  $Id: chgloc.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!  calculate local charge density n(r)
!
!***********************************************************************

      SUBROUTINE CHGLOC(NBANDS,NKDIM,LDIMP,NIONS,ISPIN, &
     &     PAR,FERWE)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION PAR(NBANDS,NKDIM,LDIMP,NIONS,ISPIN)
      DIMENSION FERWE(NBANDS,NKDIM,ISPIN)

      RETURN
      END
