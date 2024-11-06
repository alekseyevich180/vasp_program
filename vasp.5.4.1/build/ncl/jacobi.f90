# 1 "jacobi.F"
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

# 2 "jacobi.F" 2 
!=======================================================================
! RCS:  $Id: jacobi.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
! Module containing fast T3D/T3E specific implementation of
! Jacobis matrix diagonalization
!  written by I.J.Bush at Daresbury Laboratory in Feburary 1995.
! based on an algorithm  described by Littlefield (see below)
! put in MODULE and wrapper for VASP by gK
!=======================================================================

! should be compiled only if shmem-put is allowed (T3D_SMA)
! but the routine seems to work even on T3E with data streaming enabled
! thus (1._q,0._q) can always use it if F90_T3D is specified in the
! makefile

!#if defined(F90_T3D) && defined(gammareal)

# 1305

 MODULE jacobi
   LOGICAL :: LJACOBI=.FALSE.
   CONTAINS

      SUBROUTINE jacDSSYEV(COMM,AMATIN,W,N)
      USE prec
      USE mpimy
      TYPE (communic) COMM
      INTEGER N             ! order of the matrix
      COMPLEX(q)    AMATIN(N,N)   ! input/output matrix
      REAL(q) W(N)          ! eigenvalues
      WRITE(*,*) 'internal ERROR:  jacDSSYEV is not supported'
      CALL M_exit(); stop
      END SUBROUTINE

 END MODULE

