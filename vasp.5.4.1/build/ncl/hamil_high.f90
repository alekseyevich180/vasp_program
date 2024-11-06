# 1 "hamil_high.F"
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

# 2 "hamil_high.F" 2 
MODULE hamil_high
  USE prec
  USE hamil
! handle for orbital magnetization and metaGGAs
  TYPE ham_handle
     COMPLEX(q),POINTER :: AVTOT(:,:)    ! local vector magnetization potential
     COMPLEX(q)  ,POINTER    :: AVEC(:,:)     ! soft part of vector magnetization potential
     COMPLEX(q),POINTER :: MUTOT(:,:)    ! derivative of energy density with respect to kinetic energy density
     COMPLEX(q)  ,POINTER    :: MU(:,:)       ! same as MUTOT, but on GRID instead of GRIDC
  END TYPE ham_handle

!***********************************************************************
!
! this module implements high level routines to calculate the action
! of the Hamiltonian onto a wavefunction
! high level routines for manipulations wavefunctions are implemented
! as well
!
!***********************************************************************
END MODULE hamil_high

