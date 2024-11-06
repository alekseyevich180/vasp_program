# 1 "smart_allocate.F"
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

# 2 "smart_allocate.F" 2 
!**********************************************************************
!
! this module can be used if time for allocate and
! deallocation turns out to be a problem
!
!**********************************************************************


MODULE smart_allocate
  USE prec
  IMPLICIT NONE
  CONTAINS

  SUBROUTINE SMART_ALLOCATE_REAL(A,N)
    REAL(q),POINTER :: A(:)
    INTEGER N

    IF (ASSOCIATED(A)) THEN
       IF (SIZE(A) < N) THEN
          DEALLOCATE(A)
       ENDIF
    ENDIF
    IF (.NOT. ASSOCIATED(A)) THEN
       ALLOCATE(A(N))
    ENDIF

  END SUBROUTINE SMART_ALLOCATE_REAL

  SUBROUTINE SMART_ALLOCATE_WAVE(A,N)
    COMPLEX(q),POINTER :: A(:)
    INTEGER N

    IF (ASSOCIATED(A)) THEN
       IF (SIZE(A) < N) THEN
          DEALLOCATE(A)
       ENDIF
    ENDIF
    IF (.NOT. ASSOCIATED(A)) THEN
       ALLOCATE(A(N))
    ENDIF

  END SUBROUTINE SMART_ALLOCATE_WAVE

  SUBROUTINE SMART_ALLOCATE_COMPLEX(A,N)
    COMPLEX(q),POINTER :: A(:)
    INTEGER N

    IF (ASSOCIATED(A)) THEN
       IF (SIZE(A) < N) THEN
          DEALLOCATE(A)
       ENDIF
    ENDIF
    IF (.NOT. ASSOCIATED(A)) THEN
       ALLOCATE(A(N))
    ENDIF

  END SUBROUTINE SMART_ALLOCATE_COMPLEX

END MODULE smart_allocate
      
