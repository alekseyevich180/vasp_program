# 1 "dmft.F"
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

# 2 "dmft.F" 2 
MODULE dmft 
  USE prec
  USE fock
  USE chi_base
  USE wpot
  USE lattice
  USE xi
  USE mlwf
  USE constant
  IMPLICIT NONE

!**********************************************************************
!
! The routines in this module compute
!
! \int_BZ dk1 dk2 dk3 dk4
!  < w_n1 w'_n2 | W | w_n4 w'_n3 > =
!
! \int_BZ dk1 dk2 dk3 dk4
!  \int_d3r d3r' w*_k1,n1(r) w_k2+q,n4(r) w*_k2,n2(r') w_k1-q,n3(r') W(r',r)
!
! where W is the bare or a screened Coulomb interaction
!
! Caveat: the difference vectors between any two k-points must be
! included in the k-point set. This requires to use Gamma centered
! meshes.
!
! These routines are similar to the stuff in the local_field.F and
! ump2.F routines.
!
!**********************************************************************

!
  REAL(q),ALLOCATABLE :: WPOS(:,:)
! two electron four orbital integrals in a Wannier basis
  COMPLEX(q), ALLOCATABLE :: TWOE_4WANNIER(:,:,:,:,:,:,:)
 
! some parameter from ADD_XI routine, needed for PIJKL
  REAL(q), PARAMETER :: PIJKL_EMPTY_THRESHHOLD=0.00001
  REAL(q), PARAMETER :: DEGENERACY_THRESHOLD=1E-2_q
  LOGICAL,SAVE :: LWPOT
!
! wannier states for DMFT matrix elements
!
  INTEGER,SAVE :: NWLOW              
  INTEGER,SAVE :: NWHIGH              
  INTEGER,SAVE :: NWTOT
!print all matrix elements calculated
  LOGICAL,SAVE :: LPRINTALL

# 3570


END MODULE dmft
