# 1 "cg.F"
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

# 2 "cg.F" 2 
!**********************************************************************
!
! Module which implementes a force-based Conjugate Gradient optimizer
!
!**********************************************************************

  MODULE cg
    USE prec
    USE lattice

    IMPLICIT NONE
    SAVE
    PRIVATE 
    PUBLIC :: cg_step,cg_init  ! call cg_init from opt_init

    INTEGER :: nions,iu0,iu6
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: Fold,R,g,g_unit

    REAL(q) :: finite_step,maxmove
    LOGICAL :: fdstep

!**********************************************************************
!
! Conjugate gradient method
!
!**********************************************************************
  CONTAINS
    SUBROUTINE cg_step(optflag,posion,toten,force,latt_a,latt_b)

      LOGICAL optflag
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)

      REAL(q) :: a1,a2,gam,step_size
      REAL(q) :: fp1,fp2,Favg,curvature

      R = posion

! Convert the position into Cartesian coordinates
      CALL dirkar(nions,R,latt_a)

      IF (fdstep) THEN
        IF (iu6>=0) WRITE(iu6,'(A5,A)') 'OPT:',' CG fdstep'
        fdstep = .false.
        optflag = .true.
        a1 = ABS(SUM(force*Fold))
        a2 = SUM(Fold*Fold)
        IF ((a1 .LE. 0.5_q*a2) .AND. (a2 .NE. 0.0_q)) THEN
          gam = SUM(force*(force-Fold))/a2
        ELSE
          gam = 0.0_q
        END IF
        IF (iu6>=0) WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' CG gam',gam
        g = force+g*gam
        g_unit = g/SQRT(SUM(g*g))
        Fold = force
! Move from the original configuration
        R = R+g_unit*finite_step
      ELSE
        IF (iu6>=0) WRITE(iu6,'(A5,A)') 'OPT:',' CG step'
        fdstep = .true.
        optflag = .false.
        fp1 = SUM(Fold*g_unit)
        fp2 = SUM(force*g_unit)
        curvature = (fp1-fp2)/finite_step
        IF (iu6>=0) WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' CG curvature',curvature
        IF (curvature .LT.  0.0_q) THEN
         step_size = maxmove
        ELSE
          Favg = 0.5_q*(fp1+fp2)
          step_size = Favg/curvature
          IF (ABS(step_size) .GT. maxmove) THEN
            step_size = SIGN(maxmove,step_size)-SIGN(finite_step,step_size)
          ELSE
            step_size = step_size-0.5_q*finite_step   ! (*)
          END IF
        END IF
        IF (iu6>=0) WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' CG step_size',step_size
! Move from the configuration after the fd_step, so (*) has a "-" sign
        R = R+g_unit*step_size
      END IF

! convert position back to direct coordinates
      CALL kardir(nions,R,latt_b)
! update posion
      posion = R

    END SUBROUTINE cg_step

!**********************************************************************
! Conjugate gradient initilizer
!**********************************************************************

    SUBROUTINE cg_init(T_INFO,IO)
      USE base
      USE poscar
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO

      INTEGER IDUM,IERR,N
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      nions=T_INFO%nions
      iu0=IO%IU0
      iu6=IO%IU6

! read in variables used for conjugate gradients

      finite_step=0.005_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FDSTEP','=','#',';','F', &
     &            IDUM,finite_step,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''FDSTEP'' from file INCAR.'
         CALL M_exit(); stop
      ENDIF

      maxmove=0.2_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'MAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''MAXMOVE'' from file INCAR.'
         CALL M_exit(); stop
      ENDIF

! initialize allocatable variables
 
      ALLOCATE(Fold(3,nions),R(3,nions))
      ALLOCATE(g(3,nions),g_unit(3,nions))

      Fold=0._q
      g=0._q

      IF (iu6>=0) THEN
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' CG, FDSTEP',finite_step
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' CG, MAXMOVE',maxmove
      END IF

! start optimization with a finite difference step
      fdstep=.true.

    END SUBROUTINE cg_init

  END MODULE cg
