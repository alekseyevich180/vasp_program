# 1 "chain.F"
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

# 2 "chain.F" 2 
!**********************************************************************
!
! Module which controls the running of four methods, the nudged
! elastic band (neb.F), the dimer (dimer.F), lanczos (lanczos.F), and
! the dynamical matrix method (dynmat.F).  The purpose of this module
! is to determine which of these methods should be run, and to do so at
! each ionic step.
!
! A set of force-based optimers are also included for use with normal
! optimizations, or transition state calculations.  The optimizers are
! steepest-descent, quick-min, conjugate-gradients, and LBFGS.
!
! NOTE: the vasp folks have now added their own dynamical matrix
! code.  The only advantage of this dynamical matrix implementation
! is that you can combine forces from multiple vasp runs.  This
! makes it easier to separate a large calculation into several jobs,
! and systematically check for convergence of the normal modes, or
! prefactors of reactions.
!
! For more information see: http://theory.cm.utexas.edu/vtsttools/
!
! Contributers:
!   Andri Arnaldsson
!   Sam Chill
!   Graeme Henkelman
!   Hannes Jonsson
!   Daniel Sheppard
!   Blas Uberuaga
!   Lijun Xu
!   Liang Zhang
!
! Email: henkelman@utexas.edu
!
!**********************************************************************

  MODULE chain
    USE prec
    USE main_mpi
    USE poscar
    USE lattice
    USE neb
    USE dynmat
    USE dimer
    USE lanczos
    USE bbm
    USE instanton
    USE opt

    IMPLICIT NONE
    SAVE
    PRIVATE
    PUBLIC :: chain_force, chain_init, parallel_tempering
    PUBLIC :: Sum_Chain, And_Chain, LHYPER_NUDGE
    INTEGER mpmd_client_rank
# 58

    INTEGER :: ICHAIN, IOPT, ISIF_local
    LOGICAL :: optflag, fconverge, ftot_flag, LINTERACT, LMPMD
    LOGICAL :: sconverge, cell_flag, twodim_flag
    REAL(q) :: EDIFFG_local, ftot_val, jacobian, PSTRESS_local 
    REAL(q),ALLOCATABLE :: Free(:,:)

# 72


!!! Tempering variables

    INTEGER :: nions_max

! tag VCAIMAGES allows to perform two MD's with e.g. different POTCAR files
! allowing to do force averaging
    LOGICAL :: LVCAIMAGES
    REAL(q) :: VCAIMAGES   ! weight of first images, weight of second image is 1-VCAIMAGES
! average forces over two images
! it is recommended to combine this with the VCA type calculations
! where the POSCAR/POTCAR/INCAR etc. are strictly identical
! and the VCA tag is set for (1._q,0._q) type to 0 and 1 respectively
! in the two POTCAR files
! forces and energies are averaged over the two nodes
! IMAGES in mpi_main.F is forced to two in this case


! the tag LTEMPER allows to perform parallel tempering runs
! IMAGES images of VASP are kicked up and the
! temperatures between the images (TEBEG) are swapped using
! a Monte-Carlo algorithm (with corresponding temperature scaling)
    LOGICAL :: LTEMPER
    INTEGER :: NTEMPER     ! replica exchange every NTEMPER steps

    INTEGER, ALLOCATABLE, SAVE :: ATTEMPTS(:)
    INTEGER, ALLOCATABLE, SAVE :: SUCCESS(:)
!!!


!**********************************************************************
!  General force routine for any method using the repeated image mode.
!  The variable ICHAIN determines which method to use (see chain_init)
!**********************************************************************

  CONTAINS
    SUBROUTINE chain_force(nions, posion, toten, force, stress, a, b, IU6)
!SUBROUTINE chain_force(nions,posion,toten,force,a,b,IU6)
      INTEGER :: nions,ni,nj,IU6,I
      REAL(q) :: ftot, frms, fmaxatom, ftemp, fmaxdim, toten
      REAL(q),DIMENSION(3,nions) :: posion,force
      REAL(q),DIMENSION(3,nions) :: posion_vasp, force_vasp, force_dimlan
      REAL(q),DIMENSION(3,3) :: a, b
      REAL(q),DIMENSION(3,3) :: stress, hstress, stress_vasp
      REAL(q),DIMENSION(3,3) :: sdA, sdB, sdA2
      REAL(q) :: stot, smaxdim
      REAL(q) :: omega
      LOGICAL :: stopcar_exists, newcar_exists
      INTEGER NIOND, NIONPD, NTYPPD, NTYPD, ierr
      TYPE (latt):: LATT_CUR
      TYPE (TYPE_info) :: T_I
      TYPE (dynamics) :: DYN
      TYPE (in_struct) :: IO

!!! try parallel temperating and VCA code here

      INTEGER node
      REAL(q),ALLOCATABLE,SAVE :: force_all(:,:,:)

!      IF (images==0) RETURN


      node = comm_chain%node_me

!======================================================================
! Parallel tempering return
!======================================================================
      IF (LTEMPER) THEN
         RETURN
!======================================================================
! VCA average forces over two images
!======================================================================
      ELSE IF (LVCAIMAGES) THEN
         ALLOCATE(force_all(3,nions,2))
         force_all(:,:,1:2) = 0
         force_all(:,:,node) = force
         CALL M_sum_d( comm_chain, force_all(1,1,1), nions*3*2)
         force = 0
         force(:,:) = force_all(:,:,1)*VCAIMAGES + force_all(:,:,2)*(1-VCAIMAGES)
         force = force
         DEALLOCATE(force_all)
         RETURN
      ENDIF

!!!

! Interactive Mode only have rank 0 process write data
      IF (LINTERACT .AND. IU6>0) THEN
! Write the force-energy file.
        WRITE(*,*) 'LINTERACT: Writing FU file.'
        OPEN(UNIT = 1, FILE = "FU")
        WRITE(1, *) toten
        DO ni=1,nions
          WRITE(1, *) force(1, ni), force(2, ni), force(3, ni)
        ENDDO
        CLOSE(1)
! Wait for the NEWCAR or STOPCAR.
        WRITE(*,*) 'LINTERACT: Waiting for NEWCAR.'
        DO
          INQUIRE(FILE = "STOPCAR", EXIST = stopcar_exists)
          INQUIRE(FILE = "NEWCAR", EXIST = newcar_exists)
          IF (stopcar_exists .OR. newcar_exists) exit
          CALL Sleep(1)
        ENDDO
! Read the NEWCAR in.
        IF (newcar_exists) THEN
          CALL RD_POSCAR_HEAD(LATT_CUR, T_I, NIOND, NIONPD, NTYPD, NTYPPD, IO%IU0, IO%IU6)
          CALL RD_POSCAR(LATT_CUR, T_I, DYN, NIOND, NIONPD, NTYPD, NTYPPD, IO%IU0, IO%IU6)
          posion = DYN%POSION
          CALL unlink("NEWCAR")
        ENDIF
      ENDIF


          CALL MPI_Barrier(comm_chain%mpi_comm, ierr)


! /Interactive Mode

! MPMD mode
# 198

      CALL MPI_Barrier(comm_chain%mpi_comm, ierr)
! End MPMD mode

! save the original stress tensor
      stress_vasp = stress
! if we are not using vasp optimizers we need to add
! the pressure to the stress tensor
      IF (IOPT .NE. 0) THEN
! find volume of cell
        omega =  a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) 
        omega = omega - a(1,2)*(a(2,1)*a(3,3) - a(3,1)*a(2,3)) 
        omega = omega + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2)) 
        DO I=1,3
          stress(I,I) = stress(I,I) - omega * PSTRESS_local 
        ENDDO
      ENDIF

! for dimer/lanczos, save the force and posion for the vasp stop criteria
      IF (ICHAIN==2 .OR. ICHAIN==3) THEN
        force_vasp = force
        posion_vasp = posion
      ENDIF

# 240


! optflag indicates who has control
!  true: optimizer is active
!  false: chain method is active (for dimer/lanczos)
      IF (IMAGES==0 .AND. ICHAIN==0) THEN
        optflag = .TRUE.
      ELSE
        IF (ICHAIN==0) THEN
          CALL neb_step(optflag,posion,toten,force,stress,a,b)
!CALL neb_step(optflag,posion,toten,force,a,b)
        ELSEIF (ICHAIN==1) THEN
          CALL dynmat_step(optflag,posion,toten,force,a,b)
        ELSEIF (ICHAIN==2) THEN
          CALL dimer_step(optflag,posion,toten,force,a,b)
          force_dimlan=force
        ELSEIF (ICHAIN==3) THEN
          CALL lanczos_step(optflag,posion,toten,force,a,b)
          force_dimlan=force         
        ELSEIF (ICHAIN==6) THEN
          CALL bbm_step(optflag,posion,toten,force,a,b) 
!INS_BEGIN
        ELSEIF (ICHAIN==4) THEN
          CALL instanton_step(optflag,posion,toten,force,a,b)
!INS_END
        ENDIF
      ENDIF

      IF (IU6>0) THEN
        WRITE(IU6,*) "stress matrix after NEB project (eV)"
        WRITE(IU6,'(3F13.5)') stress
      ENDIF

! 2D calculations (only in x,y)
      IF (twodim_flag) THEN
        stress(1,3) = 0._q
        stress(2,3) = 0._q
        stress(3,3) = 0._q
      ENDIF

! (0._q,0._q) out any added forces on frozen atoms
      force = force*Free

! for dimer/lanczos, check force criteria using the vasp force
      IF (ICHAIN==2 .OR. ICHAIN==3) THEN
        force_dimlan = force
        force = force_vasp
      ENDIF

      ftot = 0._q
      fmaxatom = 0._q
      fmaxdim = 0._q
      DO ni = 1,nions
        ftemp = 0._q
        DO nj = 1,3
          ftot = ftot + force(nj,ni)**2
          ftemp = ftemp + force(nj,ni)**2
          IF(fmaxdim.LT.ABS(force(nj,ni))) THEN
            fmaxdim = ABS(force(nj,ni))
          ENDIF
        ENDDO
        IF (ftemp .GT. fmaxatom) fmaxatom = ftemp
      ENDDO
      frms = SQRT(ftot/REAL(nions))
      fmaxatom = SQRT(fmaxatom)
      IF (IU6>=0) WRITE(IU6,4693) fmaxatom,frms
 4693 FORMAT(1x,' FORCES: max atom, RMS ',2f12.6)
      IF (IU6>=0) WRITE(IU6,4694) SQRT(ftot),fmaxdim
 4694 FORMAT(1x,' FORCE total and by dimension',2f12.6)

! added to monitor the stress vector
      stot = 0._q
      smaxdim = 0._q
      DO ni=1,3
        DO nj=1,3
          stot = stot + stress(nj,ni)**2
          IF(smaxdim .LT. ABS(stress(nj,ni))) THEN
            smaxdim = ABS(stress(nj,ni))
          ENDIF
        ENDDO
      ENDDO
 4695 FORMAT(1x,' Stress total and by dimension',2f12.6)
      IF (IU6>=0) WRITE(IU6,4695) SQRT(stot),smaxdim

! for dimer/lanczos, use projected force in optimizer
      IF (ICHAIN==2 .OR. ICHAIN==3) THEN
        force = force_dimlan
      ENDIF

! stops based on the Magnitude of the Force
      IF(ftot_flag) THEN
        fconverge = (SQRT(ftot) .LT. ftot_val)
        CALL and_chain(fconverge)
        IF (fconverge) THEN
          IF (IU6>=0) WRITE(IU6,*) 'CONVERGED based on Magnitude of Force'
          CALL M_exit(); stop
        ENDIF
      ENDIF

      IF (IOPT .NE. 0) THEN
        hstress = stress
        CALL sdotA(hstress,a)
! freeze out rotation
        hstress(2,1) = 0._q
        hstress(3,1) = 0._q
        hstress(3,2) = 0._q
        hstress = hstress/jacobian
        fconverge = (fmaxatom .LT. ABS(EDIFFG_local))
        CALL and_chain(fconverge)
        IF (.NOT. fconverge) THEN
! our own optimizers (optflag: do or do not optimize)
!CALL opt_step(optflag,posion,toten,force,a,b)
          CALL opt_step(optflag,posion,toten,force,hstress,a,b)
        ELSE
          IF (cell_flag) THEN
            sconverge = (smaxdim .LT. ABS(EDIFFG_local)*dble(nions))
            CALL and_chain(sconverge)
          ELSE 
            sconverge = .TRUE.
          ENDIF
          IF (.NOT. sconverge) THEN
            CALL opt_step(optflag,posion,toten,force,hstress,a,b)
          ELSE
            IF (IU6>=0) WRITE(IU6,*) 'OPT: skip step - force has converged'
! GH: remove any finite difference steps taken by dimer or lan
            IF (ICHAIN==2 .OR. ICHAIN==3)  posion=posion_vasp
            IF (ICHAIN==2) THEN
              CALL dimer_fin()
            ENDIF
          ENDIF
        ENDIF
      ENDIF

# 379


! for dimer/lanczos, return the true force to vasp
    IF (ICHAIN==2 .OR. ICHAIN==3) THEN
      force = force_vasp
    ENDIF
! return unprojected stress vasp
    stress = stress_vasp

    END SUBROUTINE chain_force


!**********************************************************************
! Initialize the chain (repeated image mode) and determine which of the
! three possible methods to use based on the ICHAIN variable:
!   ICHAIN==0: nudged elastic band (default)
!   ICHAIN==1: dynamical matrix
!   ICHAIN==2: dimer method
!   ICHAIN==3: lanczos method
!   ICHAIN==6: bond-boost method
!**********************************************************************

    SUBROUTINE chain_init(T_INFO, IO)
      USE base
      TYPE (in_struct) :: IO
      TYPE (type_info) :: T_INFO
      INTEGER :: NI,NJ,IU0,IU6

      INTEGER :: IERR,N,IDUM
      CHARACTER*1 :: CHARAC
      COMPLEX(q) :: CDUM 
      LOGICAL :: LDUM
      REAL(q) :: RDUM
      REAL(q) :: omega
      REAL, PARAMETER :: EVTOJ=1.60217733E-19_q
      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt):: LATT_CUR
      TYPE (TYPE_info) :: T_I
      TYPE (dynamics) :: DYN

# 421

      IU0 = IO%IU0
      IU6 = IO%IU6

! write the version number
      IF(IU6>=0) WRITE(IU6,'(/,A,/)') ' VTST: version 3.1, (03/28/14)'

      IF(IU6>=0) WRITE(IU6,*) 'CHAIN: initializing optimizer'

! initialize optimizer
      CALL opt_init(T_INFO, IO)
      optflag = .FALSE.

! initialize chain based method
      ICHAIN = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ICHAIN','=','#',';','I', &
     &            ICHAIN,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''ICHAIN'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF

      LINTERACT = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LINTERACT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LINTERACT,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LINTERACT'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF

      LMPMD = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LMPMD','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LMPMD,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LMPMD'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF

      IOPT = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'IOPT','=','#',';','I', &
     &            IOPT,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''IOPT'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF

      EDIFFG_local = 0.1_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'EDIFFG','=','#',';','F', &
     &            IDUM,EDIFFG_local,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''EDIFFG'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF

!     determines if cell should change in NEB
      cell_flag = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LNEBCELL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,cell_flag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LNEBCELL'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF

!     sets stress tensor to two dimensions
      twodim_flag = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LTWODIM','=','#',';','L', &
     &            IDUM,RDUM,CDUM,twodim_flag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LTWODIM'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF

!     USED to converge based on Magnitude of the force
      ftot_flag = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FMAGFLAG','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ftot_flag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''FMAGFLAG'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF

      ftot_val = 0.01_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FMAGVAL','=','#',';','F', &
     &            IDUM,ftot_val,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''FMAGVAL'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF

! make sure that convergence is force based when using IOPT
      IF((IOPT .NE. 0) .AND. (EDIFFG_local .GT. 0.0_q)) THEN
        IF(IU6>=0) WRITE(IU6,*) 'Must set  EDIFFG < 0 when using IOPT > 0'
        CALL M_exit(); stop
      ENDIF

! check that solid state NEB uses our Optimizers
      IF((cell_flag .EQV. .TRUE.) .AND. ((IOPT .NE. 3) .AND. (IOPT .NE. 7))) THEN
        IF(IU6>=0) WRITE(IU6,*) 'Must set  IOPT = 3 or 7 when using LNEBCELL=.TRUE.'
        CALL M_exit(); stop
      ENDIF

      PSTRESS_local = 0.0_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PSTRESS','=','#',';','F', &
     &            IDUM,PSTRESS_local,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''PSTRESS'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF
      PSTRESS_local = PSTRESS_local/(EVTOJ*1E22_q)

      ISIF_local = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ISIF','=','#',';','I', &
     &            ISIF_local,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''ISIF'' from file INCAR.'
        CALL M_exit(); stop
      ENDIF

! check that if ISIF is set, either vasp optimizers are used or IOPT=3,7
      IF((IOPT>0) .AND. (ISIF_local>2)) THEN
        IF (ISIF_local/=3) THEN
          IF(IU6>=0) WRITE(IU6,*) 'Only ISIF = 3 is supported with IOPT = 3 or 7'
          CALL M_exit(); stop
        ELSE IF((IOPT.NE.3) .AND. (IOPT.NE.7)) THEN
          IF(IU6>=0) WRITE(IU6,*) 'Must set  IOPT = 3 or 7 when using ISIF = 3'
          IF(IU6>=0) WRITE(IU6,*) 'Alternatively, use the vasp optimizers'
          CALL M_exit(); stop
        ENDIF
      ENDIF

      IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Read ICHAIN ',ICHAIN

      IF(IMAGES==0) THEN
          CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
               NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
          CALL RD_POSCAR(LATT_CUR, T_I, DYN, &
               NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)

          IF (T_I%NIONS /= T_INFO%NIONS) THEN
             IF (IU0>=0) WRITE(IU0,*)'ERROR: image mode number of ions wrong'
          CALL M_exit(); stop
          ENDIF
          omega = LATT_CUR%OMEGA
          jacobian = (omega/DBLE(T_INFO%nions))**(1._q/3._q)*sqrt(DBLE(T_INFO%nions))
      ENDIF

      IF(ICHAIN==0) THEN
        IF(IMAGES>0) THEN
          jacobian = 0.0_q
          IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Running the NEB'
        ENDIF
        CALL neb_init(T_INFO,IO,jacobian)
      ELSEIF(ICHAIN==1) THEN
        IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Dynamical Matrix'
        CALL dynmat_init(T_INFO,IO)
      ELSEIF(ICHAIN==2) THEN
        IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Dimer method'
        CALL dimer_init(T_INFO,IO)
        IF(iopt==0) THEN
          IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Must set IOPT>0 to use the Dimer method'
          CALL M_exit(); stop
        ENDIF
      ELSEIF(ICHAIN==3) THEN
        IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Lanczos method'
        CALL lanczos_init(T_INFO,IO)
        IF(iopt==0) THEN
          IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Must set IOPT>0 to use the Lanczos method'
          CALL M_exit(); stop
        ENDIF
      ELSEIF (ICHAIN==4) THEN
        IF (IO%IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Instanton method'
        CALL instanton_init(T_INFO,IO)
      ELSEIF(ICHAIN==6) THEN
      IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Running the bond-boost method'
      CALL bbm_init(T_INFO,IO)
      ENDIF

! Make vector to (0._q,0._q) out force on frozen atoms
      ALLOCATE(Free(3,T_INFO%nions))
      Free(:,:)=1._q
      IF(T_INFO%LSDYN) THEN
        DO NI=1,T_INFO%nions
          DO NJ=1,3
            IF (.NOT.T_INFO%LSFOR(NJ,NI)) Free(NJ,NI)=0._q
          ENDDO
        ENDDO
      ENDIF

# 614


# 631



!!! try putting the replica exchange code back here

! quick return, if we are not running in image mode
      IF (images==0) RETURN


! read the VCAIMAGES
      VCAIMAGES=-1

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'VCAIMAGES','=','#',';','F', &
           &        IDUM,VCAIMAGES,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) WRITE(IO%IU0,*)'Error reading item ''VCAIMAGES'' from file INCAR.'
         CALL M_exit(); stop
      ENDIF
      CALL XML_INCAR('VCAIMAGES','F',IDUM,VCAIMAGES,CDUM,LDUM,CHARAC,N)
      IF (VCAIMAGES==-1) THEN
         LVCAIMAGES=.FALSE.
      ELSE
         LVCAIMAGES=.TRUE.
      ENDIF


! LTEMPER -- use subspace diagonalization or not (default is TRUE):
      LTEMPER=.FALSE.
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LTEMPER','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LTEMPER,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
         WRITE(IO%IU0,*)'Error reading item ''LTEMPER'' from file INCAR.'
         CALL M_exit(); stop
      ENDIF
      CALL XML_INCAR('LTEMPER','L',IDUM,RDUM,CDUM,LTEMPER,CHARAC,N)

      IF (LTEMPER) THEN
! read NTEMPER
         NTEMPER=200
         CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'NTEMPER','=','#',';','I', &
              &            NTEMPER,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
         IF (((IERR/=0).AND.(IERR/=3)).OR. &
              &                    ((IERR==0).AND.(N<1))) THEN
            IF (IO%IU0>=0) &
                 WRITE(IO%IU0,*)'Error reading item ''NTEMPER'' from file INCAR.'
            CALL M_exit(); stop
         ENDIF

         CALL XML_INCAR('NTEMPER','I',NTEMPER,RDUM,CDUM,LDUM,CHARAC,N)

         ALLOCATE(ATTEMPTS(images-1))
         ALLOCATE(SUCCESS(images-1))
         ATTEMPTS=0
         SUCCESS=0
      ENDIF

      IF (LVCAIMAGES .OR. LTEMPER) THEN
         nions_max=T_INFO%NIONS
         CALL M_max_i(comm_chain, nions_max, 1 )
         IF (T_INFO%NIONS /= nions_max) THEN
            IF (IO%IU0>=0) WRITE(IO%IU0,*)'ERROR: image mode number of ions wrong'
            CALL M_exit(); stop
         ENDIF
      ENDIF

!!!

    END SUBROUTINE chain_init
        


# 765


!**********************************************************************
! Returns true if hyper nudged elastic band method is used
!**********************************************************************

      FUNCTION LHYPER_NUDGE()
      LOGICAL LHYPER_NUDGE
!      IF (images==0 .OR. spring /= 0 ) THEN
      IF (images==0 ) THEN
        LHYPER_NUDGE=.FALSE.
      ELSE
        LHYPER_NUDGE=.TRUE.
      ENDIF
      END FUNCTION LHYPER_NUDGE

!**********************************************************************
! 1 routines
!**********************************************************************

! Sum over elements

    SUBROUTINE sum_chain( value )
      REAL(q) :: value
      REAL(q) :: value_all(images)
      INTEGER node

      IF (images==0    ) RETURN
      IF (LTEMPER) RETURN

! VCAIMAGES returns average value for energies etc.
      IF (LVCAIMAGES) THEN
         node=comm_chain%node_me
         value_all=0
         value_all(node)=value
         CALL M_sum_d( comm_chain, value_all, 2)
         value=value_all(1)*VCAIMAGES+value_all(2)*(1-VCAIMAGES)
      ELSE
         CALL M_sum_d( comm_chain, value, 1 )
      ENDIF


    END SUBROUTINE sum_chain


!**********************************************************************
!
! also the logical break conditionions must be
!
!**********************************************************************

    SUBROUTINE and_chain( value )
      LOGICAL :: value
      REAL(q) :: sum

      IF (images==0) RETURN

! if (1._q,0._q) node is .FALSE., .FALSE. is returned on all nodes
      IF (value) THEN
         sum=0
      ELSE
         sum=1
      ENDIF
      CALL M_sum_d( comm_chain, sum, 1 )

      IF (sum>=1) THEN
         value=.FALSE.
      ELSE
         value=.TRUE.
      ENDIF
    END SUBROUTINE and_chain


!**************** SUBROUTINE s2ts ************************************
! transform stress tensor to force on vectors
! true_stress = (stress dot B)^T
! ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!***********************************************************************

      SUBROUTINE s2ts(V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER N
      DIMENSION V(3,3),BASIS(3,3)

      DO N=1,3
        V1 = V(N,1)*BASIS(1,1) + V(N,2)*BASIS(2,1) + V(N,3)*BASIS(3,1)
        V2 = V(N,1)*BASIS(1,2) + V(N,2)*BASIS(2,2) + V(N,3)*BASIS(3,2)
        V3 = V(N,1)*BASIS(1,3) + V(N,2)*BASIS(2,3) + V(N,3)*BASIS(3,3)
        V(N,1) = V1
        V(N,2) = V2
        V(N,3) = V3
      ENDDO

      RETURN
      END SUBROUTINE s2ts

!**************** SUBROUTINE ts2s ************************************
! transforms force on vectors to stress tensor
!  stress=true_stress^T dot A
!  direct lattice      (BASIS must be equal to A direct lattice)
! to cartesian coordinates
!***********************************************************************

      SUBROUTINE ts2s(V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER N
      DIMENSION V(3,3),BASIS(3,3)

      DO N=1,3
        V1 = V(N,1)*BASIS(1,1) + V(N,2)*BASIS(1,2) + V(N,3)*BASIS(1,3)
        V2 = V(N,1)*BASIS(2,1) + V(N,2)*BASIS(2,2) + V(N,3)*BASIS(2,3)
        V3 = V(N,1)*BASIS(3,1) + V(N,2)*BASIS(3,2) + V(N,3)*BASIS(3,3)
        V(N,1) = V1
        V(N,2) = V2
        V(N,3) = V3
      ENDDO

      RETURN
      END SUBROUTINE ts2s

      SUBROUTINE sdotA(V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER I,J,K
      DIMENSION V(3,3),BASIS(3,3),AC(3,3)
      AC=0
      DO J=1,3
         DO I=1,3
!A(I,J)=AC(I,J)
            DO K=1,3
               AC(I,J) = AC(I,J) + V(I,K)*BASIS(K,J)
            ENDDO
         ENDDO
      ENDDO
      V = AC

      RETURN
      END SUBROUTINE sdotA

      SUBROUTINE sdotB(V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER I,J,K
      DIMENSION V(3,3),BASIS(3,3),AC(3,3)
      AC=0
      DO J=1,3
         DO I=1,3
!A(I,J)=AC(I,J)
            DO K=1,3
               AC(I,J) = AC(I,J) + V(I,K)*BASIS(J,K)
            ENDDO
         ENDDO
      ENDDO
      V = AC

      RETURN
      END SUBROUTINE sdotB

!======================================================================
! Returns a unit vector along v1
!======================================================================
      FUNCTION return_unit(v1,nmax)
      INTEGER :: nmax
      REAL(q) :: v1(3,nmax)
      REAL(q),DIMENSION(3,nmax) :: return_unit
      IF (SUM(v1*v1) .NE. 0._q) THEN
        return_unit = v1*(1._q/SQRT(SUM(v1*v1)))
      ELSE
        return_unit = v1*0._q
      ENDIF
      END FUNCTION return_unit

# 985


      SUBROUTINE PARALLEL_TEMPERING(NSTEP, NIONS, POSION, VEL, TOTEN, FORCE, TEBEG, TEEND, &
        A, B, IU6)
      USE constant
      IMPLICIT NONE
      INTEGER :: NSTEP           ! step
      INTEGER :: NIONS           ! number of ions
      INTEGER :: IU6             ! std-out  (unit for OUTCAR file)
      REAL(q) :: POSION(3,nions) ! position of ions
      REAL(q) :: VEL(3,nions)    ! velocity of ions
      REAL(q) :: FORCE(3,nions)
      REAL(q) :: TOTEN           ! total energy
      REAL(q) :: TEBEG           ! temperature of ensemble
      REAL(q) :: TEEND           ! temperature of ensemble (must equal TEBEG)
      REAL(q) :: A(3,3),B(3,3)   ! lattice and reciprocal lattice
! local
      REAL(q),ALLOCATABLE :: TEBEG_old(:), TEBEG_new(:), TOTEN_all(:)
      REAL(q) :: E1, E2
      INTEGER,ALLOCATABLE :: ID(:)
      INTEGER :: SWAP
      INTEGER :: NODE, I, J
      REAL(q) :: TEBEG_new_local, E
      REAL(q) :: value


      IF (.NOT. LTEMPER) RETURN
      IF (IU6>=0) WRITE(IU6,*) 'parallel tempering routine entered NSTEP=, NTEMPER=', NSTEP,NTEMPER
      CALL RANDOM_NUMBER(value)

      IF (value <1.0_q/NTEMPER) THEN

         ALLOCATE(TEBEG_old(images), TEBEG_new(images), TOTEN_all(images), ID(images))

! my id, i.e. which subdir I am running in
         node = comm_chain%node_me

         TEBEG_old = 0
         TEBEG_old(node) = TEBEG
         TOTEN_all = 0
         TOTEN_all(node) = TOTEN
         ID = 0
         ID(node) = node

         CALL M_sum_d( comm_chain, TEBEG_old,  images)
         CALL M_sum_d( comm_chain, TOTEN_all,  images)
         CALL M_sum_i( comm_chain, ID,  images)

! sort the temperatures ascendingly and store the corresponding node id
         CALL SORT_ASC_REAL(images, TEBEG_old, ID)
!======================================================================
! now select swaps randomly
! images/2 swaps
! the array SWAP stores the list of images to be swapped
! 0 no swap
! 1 two lowest temperatures are swapped
! 2 next two are swapped and so on
!======================================================================
         CALL RANDOM_NUMBER(value)  ! values between [0,1[
         SWAP=(images-1)*value+1    ! create a random number [1,images-1]

         CALL M_bcast_i( comm_chain, SWAP,  1)
         IF (IU6>=0) WRITE(IU6,'(A,16I10)')   ' attempting swapping', SWAP

         TEBEG_new=TEBEG_old

         IF (SWAP/=0) THEN
            ATTEMPTS(SWAP)=ATTEMPTS(SWAP)+1

            E=(1.0_q/TEBEG_old(SWAP)-1.0_q/TEBEG_old(SWAP+1))/ BOLKEV * (TOTEN_all(ID(SWAP))-TOTEN_all(ID(SWAP+1)))
            E1=(1.0_q/TEBEG_old(SWAP)-1.0_q/TEBEG_old(SWAP+1))/ BOLKEV
            E2=(TOTEN_all(ID(SWAP))-TOTEN_all(ID(SWAP+1)))
            E=EXP(E)
            CALL RANDOM_NUMBER(value) ! values between [0,1[
            CALL M_bcast_d( comm_chain, value,  1)
            IF (value>E) THEN
! Metropolis forbids swap
               SWAP=0
            ELSE
               SUCCESS(SWAP)=SUCCESS(SWAP)+1
            ENDIF
         ENDIF
         IF (SWAP/=0) THEN
            TEBEG_new(SWAP)  =TEBEG_old(SWAP+1)
            TEBEG_new(SWAP+1)=TEBEG_old(SWAP)
         ENDIF
         IF (IU6>=0) WRITE(IU6,'(A,16F10.4)') '  1/T1-1/T2         ', E1
         IF (IU6>=0) WRITE(IU6,'(A,16F10.4)') '  E1  -E2           ', E2
         IF (IU6>=0) WRITE(IU6,'(A,16F10.7)')   '            random  ', value
         IF (IU6>=0) WRITE(IU6,'(A,16I10)')   '            swapping', SWAP

         IF (IU6>=0) WRITE(IU6,'(A,16F14.7)') ' parallel tempering old TOTEN ', (TOTEN_all(ID(I)),I=1,images)
         IF (IU6>=0) WRITE(IU6,'(A,16F14.7)') ' parallel tempering old TEBEG ', TEBEG_old
         IF (IU6>=0) WRITE(IU6,'(A,16F14.7)') ' parallel tempering new TEBEG ', TEBEG_new
         IF (IU6>=0) WRITE(IU6,'(A,16F14.7)') ' Acceptance ratio for swaps          ', REAL(SUCCESS,q)/MAX(1,ATTEMPTS)

         DO I=1,images
            IF ( ID(I)==node) THEN
               IF (TEBEG /= TEBEG_old(I)) THEN
                  WRITE(*,*) 'internal error in PARALELL_TEMPERING:', I,ID(I), TEBEG, TEBEG_old(I)
               ENDIF
               TEBEG_new_local=TEBEG_new(I)
            ENDIF
         END DO

         VEL(:,:)=VEL(:,:)*SQRT(TEBEG_new_local/TEBEG)

         TEBEG=TEBEG_new_local
         TEEND=TEBEG_new_local

         IF (IU6>=0) WRITE(IU6,*)
         IF (IU6>=0) WRITE(IU6,"('   TEBEG  = ',F6.1,';   TEEND  =',F6.1)") TEBEG,TEEND
         IF (IU6>=0) WRITE(IU6,*)

         DEALLOCATE(TEBEG_old, TEBEG_new, TOTEN_all, ID)
      ENDIF


    END SUBROUTINE PARALLEL_TEMPERING

!**********************************************************************
! sorts RA in descending order, and rearanges an index array RB
! seems to be a quicksort, but I am not sure
! subroutine writen by Florian Kirchhof
!**********************************************************************

  SUBROUTINE SORT_ASC_REAL(N, RA, RB)
    REAL(q) :: RA(N)
    INTEGER :: RB(N)
    REAL(q) :: RRA
    INTEGER :: RRB
    INTEGER :: N, L, IR, J, I

    IF (N<=1) RETURN

    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       RRA=RA(L)
       RRB=RB(L)
    ELSE
       RRA=RA(IR)
       RRB=RB(IR)
       RA(IR)=RA(1)
       RB(IR)=RB(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          RA(1)=RRA
          RB(1)=RRB
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(RA(J).GT.RA(J+1))J=J+1
       ENDIF
       IF(RRA.GT.RA(J))THEN
          RA(I)=RA(J)
          RB(I)=RB(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    RA(I)=RRA
    RB(I)=RRB
    GO TO 10
  END SUBROUTINE SORT_ASC_REAL

!!! end parallel tempering code

END MODULE chain
