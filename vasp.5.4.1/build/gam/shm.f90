# 1 "shm.F"
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

# 2 "shm.F" 2 
!************************************************************************
! RCS:  $Id: shm.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!  this module contains the routines required to handle
!  shmem on the T3D
!  it allocates the required amount and returns a position to the
!  workspace as required
!
!***********************************************************************

!=======================================================================
!
! the first subroutine calculates the maximum amount of workspace
! required
!
!=======================================================================

    SUBROUTINE SHM_MAX(WDES, MPLMAX, MALLOC)
      USE prec
      USE wave
      USE main_mpi
      USE scala
      IMPLICIT NONE
      INTEGER MPLMAX            ! maxmimum number of fft coefficients
      INTEGER MALLOC            ! amount of workspace need
      INTEGER MSCALA
      TYPE (wavedes)    WDES

      INTEGER, PARAMETER :: MCOMP_=1
# 33

      MALLOC=0

# 68


    END SUBROUTINE SHM_MAX

!=======================================================================
!
! the second subroutine allocates the required amount of workspace
!
!=======================================================================

    SUBROUTINE SHM_ALLOC(MALLOC)
      USE prec
      USE scala
      IMPLICIT NONE
      INTEGER MALLOC            ! amount of workspace need
      INTEGER MALLOC_DONE       ! amount of workspace allocated
      INTEGER INFO
# 92


      RETURN
    END SUBROUTINE SHM_ALLOC


!=======================================================================
!
! check whether sufficient workspace was allocated
!
!=======================================================================

    FUNCTION ISHM_CHECK(MALLOC)
      USE prec
      IMPLICIT NONE
      INTEGER ISHM_CHECK        ! returns 1 if workspace is sufficient
      INTEGER MALLOC            ! amount of workspace need
      INTEGER MALLOC_DONE       ! amount of workspace allocated

# 121

      ISHM_CHECK=0

      RETURN
    END FUNCTION ISHM_CHECK

