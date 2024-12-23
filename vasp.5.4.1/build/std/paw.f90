# 1 "paw.F"
!#define debug
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

# 3 "paw.F" 2 
!*******************************************************************
! RCS:  $Id: paw.F,v 1.14 2003/06/27 13:22:22 kresse Exp kresse $
!
!  PAW-module
!  implements the top level PAW functions
!  most of the functionality is implemented in the MODULES radial
!
!  all routine written by Georg Kresse
!  (even the symmetrization routines :)
!*******************************************************************
  MODULE pawm
    USE prec
    USE pseudo
    USE paw
    CONTAINS

!*******************************************************************
!
!  start up procedure for PAW
!  checks internal consistency of the QPAW
!  and sets up the compensation charges
!  on the radial grid, and spline coefficients which are used
!  to interpolate compensation charges in us.F
!
!*******************************************************************

      SUBROUTINE SET_PAW_AUG(NTYPD, P, IU6, LMAX_CALC, LCOMPAT)
        USE radial
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        INTEGER NTYPD
        INTEGER IU6
        TYPE (potcar),TARGET :: P(NTYPD)
        TYPE (rgrid),POINTER :: R
        INTEGER LMAX_CALC  ! if -1 mimic US PP (see below)
! local variables
        INTEGER NTYP,CHANNELS,LMAX,L,LP,N,I
        INTEGER, PARAMETER :: NQ=2
        REAL(q)  QQ(NQ)
        REAL(q)  A(NQ),B,ALPHA
        REAL(q)  SUM,QR,BJ,STEP,X
        REAL(q),PARAMETER ::  TH=1E-6_q
        LOGICAL LCOMPAT


! if LMAX_CALC == -1 the PAW will run in a special mode
! resulting in essentially US-PP like behaviour
! i.e. all terms are linearized around the atomic reference configuration
! and (1._q,0._q) center terms are not evaluated

        IF (LMAX_CALC ==-1) THEN
           MIMIC_US=.TRUE.
        ELSE
           MIMIC_US=.FALSE.
        ENDIF

        LMAX_MIX=0

        typ: DO NTYP=1,NTYPD
           IF (.NOT. ASSOCIATED(P(NTYP)%QPAW)) CYCLE

! maximal L mixed in Broyden mixer
           LMAX_MIX=MAX(LMAX_MIX,P(NTYP)%LMAX_CALC)

           R => P(NTYP)%R
           CALL  RAD_ALIGN(R,R%RMAX) ! reallign RMAX with grid

           CALL RAD_CHECK_QPAW( R, P(NTYP)%LMAX, &
                P(NTYP)%WAE, P(NTYP)%WPS, P(NTYP)%QPAW, P(NTYP)%QTOT , P(NTYP)%LPS, P(NTYP)%RDEP )
           IF ((USELDApU().OR.LCALC_ORBITAL_MOMENT()).AND.INTEGRALS_LDApU()) THEN
              CALL OVERLAP_AE(R,P(NTYP)%RDEP,NTYP,NTYPD,P(NTYP)%LMAX,P(NTYP)%WAE,P(NTYP)%LPS)
           ENDIF

           CHANNELS=P(NTYP)%LMAX
           DO L=1, CHANNELS
              DO LP=1, P(NTYP)%LMAX
                 IF (P(NTYP)%LPS(L)==P(NTYP)%LPS(LP)) THEN
                    P(NTYP)%QION(L,LP)=P(NTYP)%QPAW(L,LP,0)
                 ENDIF
              ENDDO
           ENDDO

           LMAX=0
           DO I=1,CHANNELS
              LMAX=MAX( P(NTYP)%LPS(I),LMAX )
           ENDDO

           LMAX=LMAX*2+1              ! maximum l in augmentation charges
! to allow use of (1._q,0._q)-center dipole operators increase by (1._q,0._q)
           ALLOCATE(P(NTYP)%QDEP(NPSRNL,5,0:LMAX), &
                    P(NTYP)%AUG (R%NMAX,0:LMAX) )

!           IF (IU6>=0) WRITE(IU6,"(' L augmenation charges for type=',I4)") NTYP

           ll: DO L=0,LMAX

! find q values
              CALL AUG_SETQ(L,R,R%RMAX,QQ,A,LCOMPAT)
!              IF (IU6>=0) WRITE(IU6,2) L,QQ,A

! setup augmentation charge on radial grid rho(r) r^2

              DO N=1,R%NMAX
                 SUM=0
                 IF (R%R(N) <= R%RMAX) THEN
                    DO I=1,NQ
                       QR=QQ(I)*R%R(N)
                       CALL SBESSEL( QR, BJ, L)
                       SUM=SUM+BJ*A(I)*R%R(N)*R%R(N)
                    ENDDO
                 ENDIF
                 P(NTYP)%AUG(N,L)=SUM
              ENDDO

! setup spline for augmentation charge
! the spline ends at PSDMAX*(NPSRNL-1)/NPSRNL see SETDEP
              STEP= R%RMAX/(NPSRNL-1)
! this resets PSDMAX which is from now on no longer identical to R%RMAX
              P(NTYP)%PSDMAX=NPSRNL*STEP

              DO N=1,NPSRNL
                 X=STEP*(N-1)
                 SUM=0
                 DO I=1,NQ
                    QR=QQ(I)*X
                    CALL SBESSEL( QR, BJ, L)
                    SUM=SUM+BJ*A(I)
                 ENDDO
                 P(NTYP)%QDEP(N,1,L) = X
                 P(NTYP)%QDEP(N,2,L) = SUM
              ENDDO
! derivative at startpoint
              X=STEP/1000
              SUM=0
              DO I=1,NQ
                 QR=QQ(I)*X
                 CALL SBESSEL( QR, BJ, L)
                 SUM=SUM+BJ*A(I)
              ENDDO
              SUM=(SUM-P(NTYP)%QDEP(1,2,L))/X
! derivative is (0._q,0._q) for all but L=1
              IF (L/=1) THEN
                 SUM=0
              ENDIF
              CALL SPLCOF(P(NTYP)%QDEP(1,1,L),NPSRNL,NPSRNL,SUM)
           ENDDO ll
           P(NTYP)%AUG_SOFT =>P(NTYP)%AUG

        ENDDO typ

! mix only L=SET_LMAX_MIX_TO components in Broyden mixer
        LMAX_MIX=MIN(LMAX_MIX,SET_LMAX_MIX_TO)

! if US-PP are mimicked no need to mix any onsite components
        IF (MIMIC_US) LMAX_MIX=-1

      END SUBROUTINE SET_PAW_AUG


!*******************************************************************
!
! WRT_RHO_PAW write the PAW occupancies to a file specified by
! IU
!
!*******************************************************************

    SUBROUTINE WRT_RHO_PAW(P, T_INFO, LOVERL, RHOLM_STORE, COMM, IU )
      USE pseudo
      USE poscar
      USE wave
      USE constant
      IMPLICIT NONE

      INTEGER IU               ! io unit
      TYPE (type_info)   T_INFO
      TYPE (potcar),TARGET::  P(T_INFO%NTYP)
      LOGICAL  LOVERL          ! overlap matrix used ?
      REAL(q) RHOLM_STORE(:)   ! storage for the channel occupancies
      TYPE(communic) :: COMM

! local variables
      TYPE (potcar),POINTER :: PP
      INTEGER NT, NI, I
      INTEGER IBASE, IADD, NELEMENTS, LYMAX, LMMAX
      INTEGER, EXTERNAL :: MAXL_AUG
      INTEGER NODE_ME, IONODE

      REAL(q), ALLOCATABLE::  BUFFER(:)


      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE


!=======================================================================
! quick return if possible
!=======================================================================
      IF (.NOT.LOVERL .OR. MIMIC_US ) RETURN

      LYMAX =MAXL_AUG(T_INFO%NTYP,P)
      LMMAX=(LYMAX+1)**2

      ALLOCATE( BUFFER(LMMAX*LMMAX))

!=======================================================================
! cycle all ions and write the required elements
!=======================================================================
      IBASE=1

      ion: DO NI=1,T_INFO%NIONS
         BUFFER=0
         IADD  =0

         NT=T_INFO%ITYP(NI)
         PP=> P(NT)

         NELEMENTS =0 
         IF (DO_LOCAL(NI)) THEN
            CALL RETRIEVE_RHOLM( BUFFER, RHOLM_STORE(IBASE:), &
                 METRIC(IBASE:), IADD, PP, .TRUE., NELEMENTS)

            IF (NELEMENTS > LMMAX*LMMAX) THEN
               WRITE(*,*)'internal ERROR: WRT_RHO_PAW running out of buffer'
               CALL M_exit(); stop
            ENDIF
            
            IBASE=IBASE+IADD
         ENDIF

         CALL M_sum_i(COMM, NELEMENTS, 1)
         CALL M_sum_d(COMM, BUFFER, NELEMENTS)
         
         IF (NODE_ME==IONODE) THEN
         WRITE(IU,'("augmentation occupancies",2I4)') NI, NELEMENTS
         WRITE(IU,'(5E15.7)') (BUFFER(I),I=1,NELEMENTS)
         ENDIF

      ENDDO ion

# 274

      DEALLOCATE(BUFFER)
    END SUBROUTINE WRT_RHO_PAW

!*******************************************************************
!
! RD_RHO_PAW read the PAW occupancies from a file specified by
! IU
! only components up to LMAX_MIX are read in
! the remaining components are required by the routine
! SET_DD_PAW
! the current program logic is that these are determined from CRHODE
! which is initialised according to atomic occupancies
!
!*******************************************************************

    SUBROUTINE RD_RHO_PAW(P, T_INFO, LOVERL, RHOLM_STORE, COMM, IU, IERR )
      USE pseudo
      USE poscar
      USE wave
      USE constant
      IMPLICIT NONE

      INTEGER IU               ! io unit
      INTEGER IERR             ! error status on return
! 0 ok, 1 error occured
      TYPE (type_info)   T_INFO
      TYPE (potcar),TARGET::  P(T_INFO%NTYP)
      LOGICAL  LOVERL          ! overlap matrix used ?
      REAL(q) RHOLM_STORE(:)   ! storage for the channel occupancies
      TYPE(communic) :: COMM

! local variables
      TYPE (potcar),POINTER :: PP
      INTEGER NT, NI, I
      INTEGER IBASE, IADD, NELEMENTS, LYMAX, LMMAX
      INTEGER, EXTERNAL :: MAXL_AUG
      INTEGER NI_READ, NELEMENTS_READ
      INTEGER NODE_ME, IONODE
      INTEGER IFLAG

      REAL(q), ALLOCATABLE::  BUFFER(:)
      CHARACTER (1) CH


      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE


!=======================================================================
! quick return if possible
!=======================================================================
      IF (.NOT.LOVERL .OR. MIMIC_US ) RETURN

      LYMAX =MAXL_AUG(T_INFO%NTYP,P)
      LMMAX=(LYMAX+1)**2

      ALLOCATE( BUFFER(LMMAX*LMMAX))

!=======================================================================
! cycle all ions and write the required elements
!=======================================================================
      IBASE=1

      ion: DO NI=1,T_INFO%NIONS
         BUFFER=0
         IADD  =0

         NT=T_INFO%ITYP(NI)
         PP=> P(NT)

         IERR=0


         IF (NODE_ME==IONODE) THEN
         NI_READ=0
         NELEMENTS_READ=0

         READ(IU,'(24X,2I4)',IOSTAT=IERR) &
              NI_READ, NELEMENTS_READ

         IF (NELEMENTS_READ > LMMAX*LMMAX) THEN
            WRITE(*,*)'internal ERROR: RD_RHO_PAW running out of buffer'
            CALL M_exit(); stop
         ENDIF

         IF (IERR == 0 ) READ(IU,*,IOSTAT=IERR) (BUFFER(I),I=1,NELEMENTS_READ)
         IF (NI_READ /= NI .OR. IERR/=0 ) THEN
            IERR=1
         ENDIF
         ENDIF

         CALL M_bcast_i( COMM, NELEMENTS_READ , 1)
         CALL M_sum_d(COMM, BUFFER, NELEMENTS_READ)

         NELEMENTS=0
         IF (DO_LOCAL(NI)) THEN
            CALL STORE_RHOLM( BUFFER, RHOLM_STORE(IBASE:), &
                 METRIC(IBASE:), IADD, PP, .TRUE. , NELEMENTS)
            IBASE=IBASE+IADD
         ENDIF

         CALL MPI_barrier( COMM%MPI_COMM, I )
         CALL M_sum_i(COMM, NELEMENTS, 1)

         IF (NODE_ME==IONODE) THEN
         IF (NELEMENTS /= NELEMENTS_READ)  THEN
           IERR=1
         ENDIF
         ENDIF

         CALL M_bcast_i( COMM, IERR, 1)
         
         IF (IERR/=0) THEN
            IF (NODE_ME==IONODE) WRITE(*,*) 'RD_RHO_PAW: ion', NI,'data corrupt'
            DEALLOCATE (BUFFER)
            RETURN
         ENDIF
         
      ENDDO ion

! Now try to find out whether or not the imaginary part
! of the occupancies were written to this file as well
      IFLAG=0; CH=' '
      IF (NODE_ME==IONODE) THEN
      READ(IU,'(A)',ADVANCE='No',END=100) CH
      IF (CH=='a') IFLAG=1
      ENDIF
 100  CONTINUE 
      CALL M_sum_i(COMM, IFLAG, 1)

! if so then read the imaginary parts from file
      IF (IFLAG>0) THEN
      IF (NODE_ME==IONODE) WRITE(*,*) 'reading imaginary part of occupancies ...'

      IBASE=1

      again: DO NI=1,T_INFO%NIONS
         BUFFER=0
         IADD  =0

         NT=T_INFO%ITYP(NI)
         PP=> P(NT)

         IERR=0


         IF (NODE_ME==IONODE) THEN
         NI_READ=0
         NELEMENTS_READ=0

         IF (NI==1) THEN
            READ(IU,'(40X,2I4)',IOSTAT=IERR) &
                 NI_READ, NELEMENTS_READ
         ELSE
            READ(IU,'(41X,2I4)',IOSTAT=IERR) &
                 NI_READ, NELEMENTS_READ
         ENDIF

         IF (NELEMENTS_READ > LMMAX*LMMAX) THEN
            WRITE(*,*)'internal ERROR: RD_RHO_PAW running out of buffer'
            CALL M_exit(); stop
         ENDIF

         IF (IERR == 0 ) READ(IU,*,IOSTAT=IERR) (BUFFER(I),I=1,NELEMENTS_READ)
         IF (NI_READ /= NI .OR. IERR/=0 ) THEN
            IERR=1
         ENDIF
         ENDIF
# 470

      ENDDO again
      ENDIF

      DEALLOCATE(BUFFER)

    END SUBROUTINE RD_RHO_PAW

!*******************************************************************
!
! this subroutine calculates the average local potential
! inside the PAW sphere it is required for relaxed core calculations
!
!*******************************************************************

    SUBROUTINE SET_AVERAGEPOT(PP)
      USE prec
      USE constant
      USE radial
      USE pseudo
      
      IMPLICIT NONE
      
      TYPE (potcar) PP

      REAL(q) :: POT(PP%R%NMAX)

! local variables
      INTEGER I
      REAL(q) SUM
      REAL(q) SCALE
      REAL(q) DHARTREE

! scale is simply Y_00
      SCALE=1/(2*SQRT(PI))
! calculate the Hatree potential of the core electrons
      CALL RAD_POT_HAR(0,PP%R,POT, PP%RHOAE, DHARTREE)
!     WRITE(0,'(8F14.7)') APOT(R%NMAX-7:R%NMAX)* R%R(R%NMAX-7:R%NMAX)/RYTOEV/AUTOA*SCALE
! add the potential of the nucleus (delta like charge at origin)
      POT=POT*SCALE-FELECT/PP%R%R*(PP%ZCORE+PP%ZVALF_ORIG)

      DO I=1,PP%R%NMAX
         POT(I)=POT(I)+PP%POTPSC(PP%R%NMAX)-POT(PP%R%NMAX)
      ENDDO
! test
!      WRITE(*,*) 'SET_AVERAGEPOT: V_h[Nc](Rmax)+V_h[Z](Rmax):',POT(PP%R%NMAX)
! test
! add the valence contributions
      POT=POT-PP%POTAE
! test
!      WRITE(*,*) 'SET_AVERAGEPOT: VAERMAX:',POT(PP%R%NMAX)
! test

!-----------------------------------------------------------------------
! calculate \int V(r) aug(r)
!-----------------------------------------------------------------------
      SUM=0
      DO I=1,PP%R%NMAX
         SUM=SUM+POT(I)*PP%AUG(I,0)*PP%R%SI(I)
      ENDDO

      PP%AVERAGEPOT(1)=SUM
! test
!      WRITE(*,*) 'SET_AVERAGEPOT: PP%AVERAGEPOT(1)=',PP%AVERAGEPOT(1)
! test

    END SUBROUTINE SET_AVERAGEPOT

!*******************************************************************
!
! this subroutine calculates the (1._q,0._q)-centre (on site) potentials
! for atomic occupancies within the augmentation spheres;
! these are stored as the "reference" potentials in POTAE and POTPS
! the actual potentials are referenced to these potentials
!
! it is called to reduce numerical errors introduced by
! different implementations of the exchange correlation functional
! used in the pseudopotential generation code and VASP
!
! additionally it calculates an energy correction to the
! atomic reference energy if the type of the exchange correlation
! functional has changed
! TODO: Cu atom: energy using new updated potential is
!                very different for LMAXPAW = -1
!                PAW terms are identical to old version
!                probably energy has to be reconsidered
!
!*******************************************************************

    SUBROUTINE SET_PAW_ATOM_POT( P , T_INFO, LOVERL, LMDIM, EALLAT, LMETAGGA, IU6 )

      USE pseudo
      USE asa
      USE poscar
      USE wave
      USE constant
      USE radial
      USE setexm
      USE pawfock_inter
      USE pawkli
      USE meta
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar),TARGET::      P(T_INFO%NTYP)
      INTEGER LMDIM
      LOGICAL LMETAGGA
      INTEGER  IU6              ! io unit for OUTCAR
      REAL(q) EALLAT            ! atomic reference energy
      LOGICAL  LOVERL

! local variables
      TYPE (potcar),POINTER:: PP
      INTEGER NT,LYMAX,NI,NDIM,LMMAX,M,K
      INTEGER,PARAMETER :: ISPIN=1
      INTEGER, EXTERNAL :: MAXL_AUG,MAXL1,MAXL
      REAL(q)  RHOLM(LMDIM*LMDIM), DLM(LMDIM*LMDIM)
      REAL(q)  CRHODE(LMDIM,LMDIM,1),CHF(LMDIM,LMDIM,1)
      REAL(q),ALLOCATABLE :: POT(:,:,:), RHO(:,:,:),POTOLD(:),POTDIFF(:)
      REAL(q) :: DOUBLEC_AE,DOUBLEC_PS,EXCG_OLD,EXCG,VX
      REAL(q) :: DOUBLEPS,DOUBLEAE,DOUBLEC_HF,EX_HF,SCALE
      REAL(q) :: DOUBLEPS_OLD,DOUBLEAE_OLD
      REAL(q) :: Z
      REAL(q) :: DEXC,DEXC_OLD,DEXCM
      REAL(q) :: DEATOM             ! change of valence-valence and core-valence energy
! related to change of the exchange correlation functional
! variables required to store core wavefunctions
      INTEGER MAXNL
      REAL(q), ALLOCATABLE :: W(:,:), EIG(:)
      INTEGER, ALLOCATABLE :: N(:), LC(:)
      LOGICAL, EXTERNAL :: USEFOCK_ONECENTER, USEFOCK_AE_ONECENTER
      INTEGER, EXTERNAL :: ONE_CENTER_NMAX_FOCKAE
      INTEGER :: NAE

      INTEGER :: NTP, NTPP, NTMAX
      REAL(q) :: WTOTAL

      LOGICAL :: VERSION46=.FALSE.
!=======================================================================
! quick return and allocation of work space
!=======================================================================
      DOUBLEC_AE=0
      DOUBLEC_PS=0
      IF (.NOT.LOVERL) RETURN

      LYMAX =MAXL_AUG(T_INFO%NTYP,P)
      LMMAX=(LYMAX+1)**2

      CALL RAD_AUXILIARY_FUNCTIONS_METAGGA(MAXL(T_INFO%NTYP,P))

      NDIM=0
      DO NT=1,T_INFO%NTYP
         IF (ASSOCIATED(P(NT)%QPAW)) THEN
            NDIM=MAX(NDIM, P(NT)%R%NMAX)
         END IF
      ENDDO

      IF (NDIM == 0) RETURN
      ALLOCATE ( POT( NDIM, LMMAX, ISPIN ), RHO( NDIM, LMMAX, ISPIN ),POTOLD(NDIM), POTDIFF(NDIM))
!=======================================================================
! cycle all ions and recalculate the local potentials
! in VASP the local reference potentials are those for the atom
! if the XC potential has been changed on the fly, the calculation
! of the reference energy (atomic energy) is somewhat involved
!=======================================================================
      SCALE=2*SQRT(PI)

      DO NT=1,T_INFO%NTYP
         DEATOM=0
         PP=> P(NT)
         IF ( .NOT. ASSOCIATED(P(NT)%QPAW )) CYCLE

! calculate  total electronic core charge
! round up in all cases, since some electronic charge might be outside the
! augmentation sphere
         CALL SIMPI(PP%R, PP%RHOAE, Z)
         PP%ZCORE= AINT(Z*sqrt(4*PI)+0.9)

! calculate the LDA/GGA exchange correlation energy of the core
! the contribution from outside the core radius (R%NMAX)
! was read from the POTCAR file (entry DEXC) and stored in P%DEXCCORE
! (vasp.5.1: the latest potentials are constructed such that DEXC is essentially (0._q,0._q))
         IF (LUSE_THOMAS_FERMI) &
              CALL PUSH_XC_TYPE(LEXCH, LDAX, ALDAC, AGGAX, AGGAC, 0.0_q)

         CALL APPLY_ONE_CENTER_AEXX()
         CALL RAD_CORE_XC( PP%R, PP%RHOAE, DEXC)
         CALL RESTORE_ONE_CENTER_AEXX

         IF (LUSE_THOMAS_FERMI) CALL POP_XC_TYPE

! the same for metaGGA
         IF (LMETAGGA) THEN
            CALL RAD_CORE_META_XC( PP%R, PP%RHOAE, PP%TAUAE, DEXCM)
         ELSE
            DEXCM=0
         ENDIF
! DEXCCORE is later subtracted from the energy
         PP%DEXCCOREM=PP%DEXCCORE+DEXCM
         PP%DEXCCORE =PP%DEXCCORE+DEXC

! switch to original potential type
         IF(.NOT. VERSION46) CALL PUSH_XC_TYPE(PP%LEXCH,1.0_q,1.0_q,1.0_q,1.0_q, 0.0_q)
         CALL RAD_CORE_XC( PP%R, PP%RHOAE, DEXC_OLD)
         IF(.NOT. VERSION46) CALL POP_XC_TYPE

! negative change of core energy due to change of exchange correlation functional
         DEATOM=DEXC_OLD-DEXC

! calculate the energies and potential of the spherical atom

! set RHOLM and CRHODE to atomic occupancies
         CALL SET_CRHODE_ATOM(CRHODE(:,:,1), PP)
! transform CRHODE to llp,LM
         RHOLM=0
         CALL TRANS_RHOLM( CRHODE(:,:,1), RHOLM(:), PP )
! get the spin up density
         CRHODE=CRHODE/2

         LYMAX =MAXL1(PP)*2

!-----------------------------------------------------------------------
! update reference (1._q,0._q) centre AE potential
! PP%POTAE stores the "valence only contribution"
! (except for the DFT exchange correlation part)
!  V_H(rho_val) + V_xc(rho_val+rho_core)
!-----------------------------------------------------------------------
         CALL SET_AVERAGEPOT(PP)

         POTOLD=0
         POTOLD(1:SIZE(PP%POTAE))=PP%POTAE
         PP%POTAE=0
         RHO=0
! get atomic charge density
         CALL RAD_CHARGE( RHO(:,:,1), PP%R, RHOLM, PP%LMAX, PP%LPS, PP%WAE )

         IF (LUSE_THOMAS_FERMI) &
              CALL PUSH_XC_TYPE(LEXCH, LDAX, ALDAC, AGGAX, AGGAC, 0.0_q)
         CALL APPLY_ONE_CENTER_AEXX()
! recalculate  V_H(rho_val) + V_xc(rho_val+rho_core) using the present
! (1._q,0._q)-centre xc potential
         CALL RAD_POT( PP%R, ISPIN, 1, 1, .FALSE.,   &
                 RHO, PP%RHOAE, PP%POTAE,  POT, DOUBLEAE, EXCG)
! test
!        CALL RAD_POT_METAGGA( PP%R, ISPIN, 1,   &
!                RHO, PP%RHOAE, PP%POTAE,  POT, DOUBLEAE, EXCG)
! test
         CALL RESTORE_ONE_CENTER_AEXX
         IF (LUSE_THOMAS_FERMI) CALL POP_XC_TYPE

! WRITE(*,*) 'all electron core+ valence ',EXCG

! recalculate  V_H(rho_val) + V_xc(rho_val+rho_core) using the original
! (1._q,0._q)-centre xc potential
         IF(.NOT. VERSION46) CALL PUSH_XC_TYPE(PP%LEXCH,1.0_q,1.0_q,1.0_q,1.0_q, 0.0_q)
         CALL RAD_POT( PP%R, ISPIN, 1, 1, .FALSE.,  &
                 RHO, PP%RHOAE, PP%POTAE,  POT, DOUBLEAE_OLD, EXCG_OLD)
! test
!        CALL RAD_POT_METAGGA( PP%R, ISPIN, 1,  &
!                RHO, PP%RHOAE, PP%POTAE,  POT, DOUBLEAE_OLD, EXCG_OLD)
! test
! WRITE(0,*) 'all electron core+ valence old',EXCG_OLD
         IF(.NOT. VERSION46) CALL POP_XC_TYPE

! change of (core+val) energy due to change of exchange correlation functional
! minus change of core energy
         DEATOM=DEATOM-EXCG_OLD+EXCG

! "double counting" correction for PAW
! if the PAW (1._q,0._q) center terms are not updated:
! it is given by  \int dr (e_xc_new(r) - v_xc_old(r)) rho(r)
         VX=EXCG_OLD-DOUBLEAE_OLD
         DOUBLEAE=EXCG-VX

         DO K=1,PP%R%NMAX
           PP%POTAE(K)=-POT(K,1,1)/SCALE
           PP%POTAE_XCUPDATED(K)=-POT(K,1,1)/SCALE
         ENDDO


         PP%POTAE_XCUPDATED=0
! at this point POT contains the old "valence" only (1._q,0._q)-electron potential
! calculate (1._q,0._q) center contributions and store them in CHF
! CHF = -\int Q_ij(r) (V^H(rho_v) + V^xc(rho_v+rho_c))
! this term cancels the contribution stored in DION by the
! pseudo-potential generation code
         DLM=0
         CALL RAD_POT_WEIGHT( PP%R, 1, 0, POT)
         CALL RAD_PROJ(  POT(:,:,1)  , PP%R,1._q, DLM, PP%LMAX, PP%LPS, PP%WAE )
! transform them using Clebsch Gordan coefficients and store in CHF
         CHF=0
         CALL TRANS_DLM( CHF(:,:,1), DLM , PP )
         CHF=-CHF
!#define DFT_core_valence
# 780

! get new (1._q,0._q) center AE potential for atom ("valence only")
! this might be different from above if the exchange correlation functional
! has changed
!  V_H(rho_val) + V_xc_updated(rho_val+rho_core)
         CALL RAD_POT( PP%R, ISPIN, 1, 1, .FALSE.,  &
                 RHO, PP%RHOAE, PP%POTAE_XCUPDATED,  POT, DOUBLEAE_OLD, EXCG_OLD)

         DLM=0
         CALL RAD_POT_WEIGHT( PP%R, 1, 0, POT)
         CALL RAD_PROJ(  POT(:,:,1)  , PP%R,1._q, DLM, PP%LMAX, PP%LPS, PP%WAE )
! transform them using Clebsch Gordan coefficients and add to CHF
         CALL TRANS_DLM( CHF(:,:,1), DLM , PP )

! add those contributions to PP%DION (PP strength parameters)
         CALL ADD_CDIJ(CHF(:,:,1), PP)

! recalculate potential again (RAD_POT_WEIGHT has multiplied it by a weight)
         CALL RAD_POT( PP%R, ISPIN, 1, 1, .FALSE.,  &
                 RHO, PP%RHOAE, PP%POTAE_XCUPDATED,  POT, DOUBLEAE_OLD, EXCG_OLD)
         DO K=1,PP%R%NMAX
            PP%POTAE_XCUPDATED(K)=-POT(K,1,1)/SCALE
         ENDDO


         IF (USEFOCK_ONECENTER()) THEN
! exact exchange energy for (1._q,0._q)-center AE wave functions
            CALL SETUP_PAWFOCK_AE(NT, PP)
            CALL CALC_PAWFOCK(NT, PP, CRHODE, CHF, EX_HF)

! WRITE(*,*) 'all electron contribution HF',EX_HF
! valence-valence non-local exchange for atom
            DEATOM=DEATOM-EX_HF

! free the (1._q,0._q)-center terms so that AE-PS is calculated
! upon calling CALC_PAWFOCK
            IF (.NOT. USEFOCK_AE_ONECENTER()) CALL RELEASE_PAWFOCK

! difference between (1._q,0._q) center AE and PS Fock energies
!  E_x(AE)-E_x(PS) = - DOUBLEC_HF
            
            CALL CALC_PAWFOCK(NT, PP, CRHODE, CHF, DOUBLEC_HF)

            DOUBLEAE = DOUBLEAE- DOUBLEC_HF
! WRITE(*,*) 'doublec HF',DOUBLEC_HF
            CALL RELEASE_PAWFOCK


# 830

! for convenience the core-valence exchange is added directly to PP%QATOM
! this however is only allowed if the (1._q,0._q) center terms
! are added to CDIJ (or if the  switch is used)
! the main reason for this is that the DFT core-valence interaction
! needs to be subtracted which is 1._q by setting the reference AE
! potential to V_xc^DFT(rho_val + rho_core) and going through SET_DD_PAW
               EX_HF=0

!core-valence HF
! AEXX_OLD=AEXX
! HFSCREEN_OLD=HFSCREEN
! AEXX=1.0
! HFSCREEN=0.0
!end core-valence HF

               CALL SET_PAW_CORE_FOCK(PP, CRHODE(:,:,1), EX_HF)

!core-valence HF
! AEXX=AEXX_OLD
! HFSCREEN=HFSCREEN_OLD
! PP%RHOAE=0
!end core-valence HF

! valence-core non-local exchange for atom
                IF (.NOT. MIMIC_US) DEATOM=DEATOM+EX_HF
# 858


         ENDIF
!-----------------------------------------------------------------------
! update the local (1._q,0._q) centre pseudo potential
!-----------------------------------------------------------------------
! store the PS potential value at PP%R%NMAX
         PP%VPSRMAX=-PP%POTPS(PP%R%NMAX)+PP%POTPSC(PP%R%NMAX)

         POTOLD=0
         POTOLD(1:SIZE(PP%POTPS))=PP%POTPS
         PP%POTPS=0
         RHO=0
         CALL RAD_CHARGE( RHO(:,:,1), PP%R, RHOLM, PP%LMAX, PP%LPS, PP%WPS )
         CALL RAD_AUG_CHARGE(  RHO(:,:,1), PP%R, RHOLM, PP%LMAX, PP%LPS,  &
              LYMAX, PP%AUG, PP%QPAW )

! get double counting corrections
         IF (LUSE_THOMAS_FERMI) CALL PUSH_XC_TYPE(LEXCH, LDAX, ALDAC, AGGAX, AGGAC, 0.0_q)
         CALL RAD_POT( PP%R, ISPIN, 1, 1, .FALSE.,  &
                 RHO, PP%RHOPS, PP%POTPS, POT, DOUBLEPS, EXCG)
! test
!        CALL RAD_POT_METAGGA( PP%R, ISPIN, 1,  &
!                RHO, PP%RHOPS, PP%POTPS, POT, DOUBLEPS, EXCG)
! test
         IF (LUSE_THOMAS_FERMI) CALL POP_XC_TYPE

         IF(.NOT. VERSION46) CALL PUSH_XC_TYPE(PP%LEXCH,1.0_q,1.0_q,1.0_q,1.0_q, 0.0_q)
         CALL RAD_POT( PP%R, ISPIN, 1, 1, .FALSE.,  &
                 RHO, PP%RHOPS, PP%POTPS, POT, DOUBLEPS_OLD, EXCG_OLD)
! test
!        CALL RAD_POT_METAGGA( PP%R, ISPIN, 1,  &
!                RHO, PP%RHOPS, PP%POTPS, POT, DOUBLEPS_OLD, EXCG_OLD)
! test
         IF(.NOT. VERSION46) CALL POP_XC_TYPE

! double counting correction for PAW (1._q,0._q) centre pseudo terms
! as above \int dr (e_xc_new(r) - v_xc_old(r)) rho_PS(r)
         VX=EXCG_OLD-DOUBLEPS_OLD
         DOUBLEPS=EXCG-VX


! at this point POT contains the old valence only (1._q,0._q)-electron potential
! add the core potential (PP%POTPSC-PP%POTPS, PP%POTPS is (0._q,0._q) here)
         POT(1:SIZE(PP%POTPS),1,1)=POT(1:SIZE(PP%POTPS),1,1)+(PP%POTPSC-PP%POTPS)*SCALE
         POTOLD(1:SIZE(PP%POTPS))=POT(1:SIZE(PP%POTPS),1,1)

! calculate (1._q,0._q) center contributions and store them in CHF
         DLM=0
         CALL RAD_POT_WEIGHT( PP%R, 1, 0, POT)
         CALL RAD_PROJ(  POT(:,:,1)  , PP%R,-1._q, DLM, PP%LMAX, PP%LPS, PP%WPS )
         CALL RAD_AUG_PROJ( POT(:,:,1), PP%R, DLM, PP%LMAX, PP%LPS, &
                  0, PP%AUG, PP%QPAW )

! transform them using Clebsch Gordan coefficients and store in CHF
         CHF=0
         CALL TRANS_DLM( CHF(:,:,1), DLM , PP )
! add those contributions to PP%DION (PP strength parameters)
         CHF=-CHF

! get new (1._q,0._q) center potential for atom (valence only)
! this might be different from above if the exchange correlation potential
! or the way the augmentation is 1._q has changed
    
! if required do more accurate augmentation of charge density using
! accurate charge restoration
         IF (ONE_CENTER_NMAX_FOCKAE()>0) THEN
            CALL RAD_AUG_CHARGE_FOCK(  RHO(:,:,1), PP%R, RHOLM, PP%LMAX, PP%LPS,  &
                    LYMAX, PP%AUG_FOCK, PP%QPAW_FOCK)
         ENDIF


! redefine local ionic pseudopotential
! ) by adding the valence contribution for the conventional standard augmentation scheme
!   and the corresponding xc functional
! ) and subtracting the current valence contribution for the PP xc functional
         IF (ONE_CENTER_NMAX_FOCKAE()>0) THEN
! actually it is an open matter what to chose as local reference potential
! e.g.  the new functional or the old
! for OEP methods I lean towards PP generation potential
! for other cases towards new xc type
            CALL PUSH_XC_TYPE(PP%LEXCH,1.0_q,1.0_q,1.0_q,1.0_q, 0.0_q)
            CALL RAD_POT( PP%R, ISPIN, 1, 1, .FALSE.,  &
                 RHO, PP%RHOPS, PP%POTPS, POT, DOUBLEPS_OLD, EXCG_OLD)
            CALL POP_XC_TYPE
            POT(1:SIZE(PP%POTPS),1,1)=POT(1:SIZE(PP%POTPS),1,1)+(PP%POTPSC-PP%POTPS)*SCALE

            POTDIFF(1:SIZE(PP%POTPS))=(POT(1:SIZE(PP%POTPS),1,1)-POTOLD(1:SIZE(PP%POTPS)))/SCALE

! redefine core pseudo potential on plane wave grid
            CALL FOURPOT_TO_Q( PP%RDEP, POTDIFF, PP%PSP(:,2), SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)
            PP%PSCORE=PP%PSP(1,2)
            CALL SPLCOF(PP%PSP(1,1), SIZE(PP%PSP,1), SIZE(PP%PSP,1), 0._q)

! redefine core pseudo potential on radial grid
            PP%POTPSC=PP%POTPSC-POTDIFF

! sanity test compare fourier transfrom of new PP%PSP with new POTPSC
! WRITE(78,'(2E16.7)') (PP%R%R(K),PP%POTPSC(K),K=1,PP%R%NMAX)
! WRITE(78,*)
! CALL POTTORHO(  PP%ZVALF, SIZE(PP%PSP,1), PP%PSP(:,2), PP%PSGMAX/ SIZE(PP%PSP,1), &
!      .TRUE. , PP%R%NMAX, PP%R%R ,  PP%POTPSC )
! WRITE(78,'(2E16.7)') (PP%R%R(K),PP%POTPSC(K),K=1,PP%R%NMAX)

! recalculate the total local potential for the atom on the grid (core+ valence)
! should be identical to original POTLOK
         ENDIF

         CALL RAD_POT( PP%R, ISPIN, 1, 1, .FALSE.,  &
                 RHO, PP%RHOPS, PP%POTPS, POT, DOUBLEPS_OLD, EXCG_OLD)
         POT(1:SIZE(PP%POTPS),1,1)=POT(1:SIZE(PP%POTPS),1,1)+(PP%POTPSC-PP%POTPS)*SCALE

         DLM=0
         CALL RAD_POT_WEIGHT( PP%R, 1, 0, POT)
         CALL RAD_PROJ(  POT(:,:,1)  , PP%R,-1._q, DLM, PP%LMAX, PP%LPS, PP%WPS )
         CALL RAD_AUG_PROJ( POT(:,:,1), PP%R, DLM, PP%LMAX, PP%LPS, &
                  0, PP%AUG, PP%QPAW )

         IF (ONE_CENTER_NMAX_FOCKAE()>0) THEN
            CALL RAD_AUG_PROJ_FOCK( POT(:,:,1), PP%R, DLM, PP%LMAX, PP%LPS, &
                 0, PP%AUG_FOCK, PP%QPAW_FOCK )
         ENDIF
! transform them using Clebsch Gordan coefficients and add to CDIJ
         CALL TRANS_DLM( CHF(:,:,1), DLM , PP )

! add those contributions to PP%DION (PP strength parameters)
         CALL ADD_CDIJ(CHF(:,:,1), PP)

! recalculate reference potential again (to get rid of weights)
         CALL RAD_POT( PP%R, ISPIN, 1, 1, .FALSE.,  &
                 RHO, PP%RHOPS, PP%POTPS, POT, DOUBLEPS_OLD, EXCG_OLD)

! now update the local potential to the new potential
         DO K=1,PP%R%NMAX
           PP%POTPS(K)=-POT(K,1,1)/SCALE
! WRITE(77,'(I4,5F14.7)') K,PP%R%R(K)/AUTOA,RHO(K,1,1)*AUTOA*SCALE,PP%POTPS(K),PP%POTPS(K)-POTOLD(K)
         ENDDO

         DOUBLEC_PS= DOUBLEC_PS- DOUBLEPS*T_INFO%NITYP(NT)*T_INFO%VCA(NT)
         DOUBLEC_AE= DOUBLEC_AE+(DOUBLEAE-PP%DEXCCORE)*T_INFO%NITYP(NT)*T_INFO%VCA(NT)

! update energy of atom
         PP%EATOM=PP%EATOM-DEATOM                      ! atomic energy without convergence correction
         PP%EATOM_CORRECTED=PP%EATOM_CORRECTED-DEATOM  ! convergence corrected energy
! update total energy of all atoms
         EALLAT=EALLAT-DEATOM*T_INFO%NITYP(NT)*T_INFO%VCA(NT)

!core-valence HF
! IF (.NOT. USEFOCK_ONECENTER()) THEN
!   AEXX=1.0
!   HFSCREEN = 0.0
!   CALL SET_PAW_CORE_FOCK(PP, CRHODE(:,:,1), EX_HF)
!   AEXX=AEXX_OLD
!   HFSCREEN=HFSCREEN_OLD
!   PP%RHOAE=0      ! switch off DFT core/valence exchange and correlation
! ENDIF
!end core-valence HF

         IF (T_INFO%DARWIN_R(NT) /=0) THEN
            WRITE(*,*) 'pseudopotential strength before Darwin'
            DO M=1,PP%LMAX
               WRITE(*,'(10F14.7)') PP%DION(:,M)
            ENDDO
            CALL RAD_POT_DARWIN( PP%R, 1, POT , T_INFO%DARWIN_V(NT), T_INFO%DARWIN_R(NT))
! weight properly for simpson integration
            CALL RAD_POT_WEIGHT( PP%R, 1, 1, POT)
! extract DARWIN potential contribution into array RHOLM
            RHOLM=0
            CALL RAD_PROJ(  POT(:,:,1)  , PP%R,1._q, RHOLM, PP%LMAX, PP%LPS, PP%WAE )
! transform them using Clebsch Gordan coefficients and add to CDIJ
            CHF=0
            CALL TRANS_DLM( CHF(:,:,1), RHOLM , PP )
! add those contributions to PP%DION (PP strength parameters)
            CALL ADD_CDIJ(CHF(:,:,1), PP)
            WRITE(*,*) 'pseudopotential strength after Darwin'
            DO M=1,PP%LMAX
               WRITE(*,'(10F14.7)') PP%DION(:,M)
            ENDDO
         ENDIF

!        CALL SET_PAW_METAGGA(PP)

      ENDDO
      
!  store double counting corrections for the atomic case
!  applied only for MIMIC_US
      DOUBLEC_PS_ATOM=DOUBLEC_PS 
      DOUBLEC_AE_ATOM=DOUBLEC_AE

      DEALLOCATE( POT, RHO , POTOLD, POTDIFF )

      RETURN
    END SUBROUTINE SET_PAW_ATOM_POT

    SUBROUTINE SET_PAW_ATOM_POT_RHF(P,T_INFO,LOVERL,EALLAT,IO)

      USE prec
      USE base
      USE poscar
      USE rhfatm
      IMPLICIT NONE
      TYPE(type_info) T_INFO
      TYPE(potcar), TARGET :: P(T_INFO%NTYP)
      TYPE(in_struct) IO
      REAL(q) EALLAT
      LOGICAL LOVERL
! local variables
      TYPE (potcar), POINTER:: PP
      INTEGER NT
      REAL(q) EATOM_OLD,DEATOM

      IF (.NOT.LOVERL) RETURN

      CALL RHFATM_WRITER(IO)

      DO NT=1,T_INFO%NTYP
         PP=> P(NT)

         EATOM_OLD=PP%EATOM

         IF (.NOT.ASSOCIATED(P(NT)%QPAW )) CYCLE

         CALL RHFATM_SET_PP(PP,IO)

         DEATOM=PP%EATOM-EATOM_OLD

         EALLAT=EALLAT+T_INFO%NITYP(NT)*DEATOM
      ENDDO

      RETURN
    END SUBROUTINE SET_PAW_ATOM_POT_RHF

!*******************************************************************
!
!  SET_PAW_CORE_FOCK
! calculates the core-valence exchange contribution
! for Hartree-Fock
!
! sum_c int dr' phi_c(r') phi_v(r') int dr phi_c(r) phi_v(r)/(r-r')
!
! this is boundary condition independent since the
! valence and core are orthogonal int dr phi_c(r) phi_v(r)=0.
! The Coloumb potential
!  V(r') = int dr phi_c(r) phi_v(r) / (r-r')
! is therefore boundary condition independent!
!
!*******************************************************************

    SUBROUTINE SET_PAW_CORE_FOCK(PP, CRHODE, HF_ENERGY)
! variables required to store core wavefunctions
      USE pseudo
      USE pawfock
      USE prec
      USE cl
      IMPLICIT NONE

      TYPE (potcar),POINTER:: PP
      REAL(q)              :: CRHODE(:,:)   ! valence occupancy
      REAL(q)              :: HF_ENERGY     ! valence-core interaction energy
!local
      INTEGER CHANC,CHANNELS
      REAL(q), ALLOCATABLE :: W(:,:), A(:,:),B(:,:), EIG(:)
      INTEGER, ALLOCATABLE :: N(:), LC(:)
      REAL(q), ALLOCATABLE :: S(:,:,:,:,:)
      REAL(q) :: NORM
      REAL(q), ALLOCATABLE :: COCC(:,:),DFOCK(:,:),DHARTREE(:,:)
      INTEGER :: CH1,CH2, LYMAX, LMMAX, RNMAX, M, LM1, LM2
      INTEGER, EXTERNAL :: MAXL1

! number of grid points
      RNMAX =PP%R%NMAX

! determine the number of core channels
      CALL CL_INIT_CORE_CONF(PP,CHANC)
      CHANNELS=CHANC+PP%LMAX
      ALLOCATE( W(RNMAX,CHANNELS),A(RNMAX,CHANNELS), B(RNMAX,CHANNELS), N(CHANNELS), LC(CHANNELS), EIG(CHANNELS))

      CALL SET_CORE_WF( PP%RHOAE, PP%POTAE , PP%R, PP%ZCORE, PP%ZVALF_ORIG , &
           W, N, LC, EIG, A_=A, B_=B)

      DO CH1=1,PP%LMAX
! copy valence partial waves to A
         A(:,CHANC+CH1)=PP%WAE(:,CH1)
         LC(CHANC+CH1) =PP%LPS(CH1)
      END DO
# 1155

      LYMAX =3    ! no g orbitals at this point

      ALLOCATE(S(CHANNELS, CHANNELS, CHANNELS, CHANNELS, 0:LYMAX))

      S=0
      CALL COLOUMB_4TERM( A, LC, PP%R, CHANNELS, S, 8)
      DEALLOCATE(A, B, W)

      LMMAX=0
      DO CH1=1,CHANNELS
         LMMAX=LMMAX+(LC(CH1)*2+1)
      ENDDO
      ALLOCATE(COCC(LMMAX,LMMAX), DHARTREE(LMMAX,LMMAX), DFOCK(LMMAX,LMMAX))
!
! set diagonal components of the occupancy matrix (core only)
!
      COCC=0
      LMMAX=0
      DO CH1=1,CHANC
         DO M=1,LC(CH1)*2+1
            COCC(LMMAX+M,LMMAX+M)=1
         ENDDO
         LMMAX=LMMAX+(LC(CH1)*2+1)
      ENDDO

! partially include valence states in the setup
! this does not work properly but might be required if
! NBANDSGW_LOW is supported
!
! set atomic occupancies
! CALL SET_CRHODE_ATOM(COCC(LMMAX+1:,LMMAX+1:), PP)
! spin
! COCC(LMMAX+1:,LMMAX+1:)=COCC(LMMAX+1:,LMMAX+1:)/2
! selectively clear states
! COCC(8:10,8:10)=0
! COCC(12:,12:)=0


      CALL CALC_DHARTREE(S,COCC,CHANNELS,LC,DHARTREE,DFOCK)

!DO CH1=1,SIZE(COCC,1)
!   WRITE(*,'(40F7.2)') COCC(:,CH1)
!ENDDO
!DO CH1=1,SIZE(COCC,1)
!   WRITE(*,'(40F7.2)') DFOCK(:,CH1)
!ENDDO

      CALL ADD_CDIJ(DFOCK(LMMAX+1:,LMMAX+1:),PP)

      HF_ENERGY=0
      DO LM1=1,PP%LMMAX
         DO LM2=1,PP%LMMAX

            HF_ENERGY=HF_ENERGY+ &
                 DFOCK(LMMAX+LM1,LMMAX+LM2)*CRHODE(LM1,LM2)
# 1214

         ENDDO
      ENDDO
! double the energy since only up component is stored in CRHODE
      HF_ENERGY=HF_ENERGY*2

      DEALLOCATE(COCC,DHARTREE,DFOCK)
      DEALLOCATE(S)
      DEALLOCATE( N, LC, EIG)
      CALL CL_CLEAR_CORE_CONF()
    END SUBROUTINE SET_PAW_CORE_FOCK




!*******************************************************************
!
! this is the central routine of the PAW method
! it calculates the on site terms and adds them to D(I,J)
! at the same time double counting corrections are calculated
!
!*******************************************************************

    SUBROUTINE SET_DD_PAW(WDES, P , T_INFO, LOVERL, &
         ISPIN, LMDIM, CDIJ, RHOLM_STORE, CRHODE, &
         E, LMETA, LASPH, LCOREL )
      USE pseudo
      USE asa
      USE poscar
      USE wave
      USE constant
      USE radial
      USE base
      USE relativistic
      USE LDAPLUSU_MODULE
      USE cl
      USE pawfock_inter
      USE setexm
      USE egrad
      USE meta
      USE setxcmeta
      USE hyperfine
      IMPLICIT NONE

      TYPE (type_info) T_INFO
      TYPE (potcar),TARGET::  P(T_INFO%NTYP)
      TYPE (wavedes)  WDES
      TYPE (energy)   E
      INTEGER LMDIM, ISPIN
      REAL(q)  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      REAL(q)  CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      REAL(q)  RHOLM_STORE(:,:)
      LOGICAL  LOVERL
      LOGICAL  LMETA      ! calculate meta GGA contribution
      LOGICAL  LASPH      ! calculate aspherical corrections to potential
      LOGICAL  LCOREL     ! calculate accurate core level shifts
! local variables
      TYPE (potcar),POINTER:: PP
      INTEGER NT,LYMAX,NI,NDIM,LMMAX,NIP,ISP,IBASE,IADD,ISIZE,K,ITMP,NCDIJ,LMAX
      INTEGER, EXTERNAL :: MAXL_AUG,MAXL1
      LOGICAL, EXTERNAL :: USEFOCK_CONTRIBUTION, USEFOCK_AE_ONECENTER
! automatic arrays crash on IBM and SGI (thank's to Lucian Anton NAG, Zhengji Zhao SGI)
!      REAL(q) DDLM(LMDIM*LMDIM),RHOLM(LMDIM*LMDIM),RHOLM_(LMDIM*LMDIM,WDES%NCDIJ)
!      REAL(q) CTMP(LMDIM,LMDIM,MAX(2,WDES%NCDIJ)),CSO(LMDIM,LMDIM,WDES%NCDIJ), &
!              CHF(LMDIM,LMDIM,WDES%NCDIJ)
!      REAL(q) COCC(LMDIM,LMDIM,MAX(2,WDES%NCDIJ)),COCC_IM(LMDIM,LMDIM)
      REAL(q),ALLOCATABLE :: DDLM(:),RHOLM(:),RHOLM_(:,:)
      REAL(q),ALLOCATABLE :: CTMP(:,:,:),CSO(:,:,:),CHF(:,:,:)
      REAL(q),ALLOCATABLE :: COCC(:,:,:),COCC_IM(:,:)
      REAL(q),ALLOCATABLE :: CMETAGGA(:,:,:)
      REAL(q),ALLOCATABLE :: POT(:,:,:), RHO(:,:,:), POTAE(:,:,:), RHOAE(:,:,:)
# 1287




      REAL(q),ALLOCATABLE :: RHOCOL(:,:,:), KINDENSCOL(:,:,:),KINDENSWCOL(:,:,:)
! core level shifts
      REAL(q),ALLOCATABLE :: DRHOCORE(:)
      REAL(q) :: DOUBLEC_AE,DOUBLEC_PS
      REAL(q) :: DOUBLEPS,DOUBLEAE
      REAL(q) :: EXCM,EXCG,EXCGA
      REAL(q) :: DOUBLEC_LDAU, DOUBLEC_HF
! euler angles of the global spin quantisation axis
      REAL(q) :: ALPHA,BETA
! #define robdeb
      INTEGER :: LMAX_TAU, LMMAX_TAU
      REAL(q),ALLOCATABLE :: TAUAE(:,:,:),TAUPS(:,:,:),MUAE(:,:,:),MUPS(:,:,:)
! kinetic energy density (true and Weizsaecker)
      REAL(q),ALLOCATABLE :: KINDENSAE(:,:),KINDENSPS(:,:)
      REAL(q),ALLOCATABLE :: WKDAE(:,:),WKDPS(:,:)
      REAL(q),ALLOCATABLE :: RHOUPD(:,:,:),RHOLMUPD(:,:)
      REAL(q),POINTER :: NULPOINTER(:)
      REAL(q) :: SPI2,TMP
      INTEGER II,RNMAX, RNMAX_CL
      INTEGER, EXTERNAL :: ONE_CENTER_NMAX_FOCKAE
# 1313

! variables required to store core wavefunctions
      INTEGER MAXNL
      REAL(q), ALLOCATABLE :: W(:,:), EIG(:)
      INTEGER, ALLOCATABLE :: N(:), LC(:)
! needed to distribute over COMM_KINTER
      INTEGER IDONE
      LOGICAL LSKIP

!=======================================================================
! quick return and allocation of work space
!=======================================================================

      DOUBLEC_AE=0
      DOUBLEC_PS=0
      E%PAWPSM=0; E%PAWAEM=0
      E%PAWPSG=0; E%PAWAEG=0
      E%PAWCORE=0
      E%PAWCOREM=0

      CL_SHIFT= 0

      NULLIFY(NULPOINTER)

      IF (.NOT.LOVERL) RETURN

! mimic US-PP just set the double counting corrections correctly
      IF (MIMIC_US) THEN
         DOUBLEC_AE=DOUBLEC_AE_ATOM
         DOUBLEC_PS=DOUBLEC_PS_ATOM

         E%PAWAE=DOUBLEC_AE
         E%PAWPS=DOUBLEC_PS
         RETURN
      ENDIF

      SPI2= 2*SQRT(PI)

      LYMAX =MAXL_AUG(T_INFO%NTYP,P)

      NDIM=0
      DO NT=1,T_INFO%NTYP
         IF (ASSOCIATED(P(NT)%QPAW)) THEN
            NDIM=MAX(NDIM, P(NT)%R%NMAX)
         END IF
      ENDDO

      IF (NDIM == 0) RETURN

      IF (LUSE_THOMAS_FERMI) CALL PUSH_XC_TYPE(LEXCH, LDAX, ALDAC, AGGAX, AGGAC, 0.0_q)

      ALLOCATE(DDLM(LMDIM*LMDIM),RHOLM(LMDIM*LMDIM))
      ALLOCATE(RHOLM_(LMDIM*LMDIM,WDES%NCDIJ))
      ALLOCATE(CTMP(LMDIM,LMDIM,MAX(2,WDES%NCDIJ)),CSO(LMDIM,LMDIM,WDES%NCDIJ),CHF(LMDIM,LMDIM,WDES%NCDIJ))
      ALLOCATE(COCC(LMDIM,LMDIM,MAX(2,WDES%NCDIJ)),COCC_IM(LMDIM,LMDIM))

      LMMAX=(LYMAX+1)**2
      NCDIJ = WDES%NCDIJ

      ALLOCATE ( POT( NDIM, LMMAX, NCDIJ ), RHO( NDIM, LMMAX, NCDIJ ), &
     &   POTAE( NDIM, LMMAX, NCDIJ ), RHOAE( NDIM, LMMAX, NCDIJ), DRHOCORE(NDIM))
      ALLOCATE (RHOCOL( NDIM, LMMAX, NCDIJ ))
# 1377

! allocate kinetic energy density arrays for metagga
      ALLOCATE (KINDENSAE(NDIM,NCDIJ),KINDENSPS(NDIM,NCDIJ)) 
      ALLOCATE (WKDAE(NDIM,NCDIJ),WKDPS(NDIM,NCDIJ))
      ALLOCATE (RHOUPD(NDIM,1,NCDIJ),RHOLMUPD(LMDIM*LMDIM,NCDIJ))
      ALLOCATE (CMETAGGA(LMDIM,LMDIM,NCDIJ))

! for spin orbit coupling set the euler angles
      IF ( WDES%LSORBIT ) &
         CALL EULER(WDES%SAXIS, ALPHA, BETA)
!=======================================================================
! cycle all ions and add corrections to pseudopotential strength CDIJ
!=======================================================================
      IBASE=1; IDONE=0

      ion: DO NI=1,T_INFO%NIONS
         NIP=NI_LOCAL(NI, WDES%COMM_INB) ! local storage index
         NT=T_INFO%ITYP(NI)

         LSKIP=.FALSE.

! DO_LOCAL represents a distribution of the work on the
! (1._q,0._q)-center terms over the procs in COMM_INB and COMM_INTER (=COMM_KIN).
! The following allows an additional round-robin distribution over COMM_KINTER as well.
         IF (DO_LOCAL(NI)) THEN
            IDONE=IDONE+1; LSKIP=(MOD(IDONE,WDES%COMM_KINTER%NCPU)+1/=WDES%COMM_KINTER%NODE_ME)
         ENDIF

! if this element is not treated locally CYCLE
         IF (.NOT. DO_LOCAL(NI).OR.LSKIP) THEN
! for PAW, set CDIJ to (0._q,0._q) if it resides on local node
! and if the element is not treated locally
            IF (ASSOCIATED(P(NT)%QPAW)) THEN
               IF (NIP /= 0) THEN
                  DO ISP=1,NCDIJ
                     CDIJ(:,:,NIP,ISP)=0
                  ENDDO
               ENDIF
            ELSE
! US PP: initialize to (0._q,0._q) if we are not on first node in COMM_INTER
! (at the end, we use a global sum over COMM_INTER)

               IF (WDES%COMM_INTER%NODE_ME*WDES%COMM_KINTER%NODE_ME /=1 .AND. NIP /=0 ) THEN
                  DO ISP=1,NCDIJ
                     CDIJ(:,:,NIP,ISP)=0
                  ENDDO
               ENDIF

            ENDIF
! in case we skip work on this ion, we still have to advance IBASE
            IF (LSKIP) THEN
               PP=> PP_POINTER(P, NI, NT)
               IBASE=IBASE+COUNT_RHO_PAW_ELEMENTS(PP)
            ENDIF
            CYCLE ion
         ENDIF

! The CMBJ parameter may depend on the atomic type
! if so, it is set here to the right entry in CMBJ_TYP
         CALL SET_CMBJ_RAD(NT)

         PP=> PP_POINTER(P, NI, NT)
         CALL SET_RSGF_TYPE(NT)
!        CALL SET_RSGF_SIMPLE(PP)
         LYMAX =MAXL1(PP)*2
         RNMAX =PP%R%NMAX

         LMAX_TAU=LYMAX+2; LMMAX_TAU=(LMAX_TAU+1)**2
         ALLOCATE (TAUAE(NDIM,LMMAX_TAU,NCDIJ),TAUPS(NDIM,LMMAX_TAU,NCDIJ), &
        &   MUAE(NDIM,LMMAX_TAU,NCDIJ),MUPS(NDIM,LMMAX_TAU,NCDIJ),KINDENSCOL(NDIM,LMMAX_TAU,NCDIJ))
!-----------------------------------------------------------------------
! first set RHOLM (i.e. the on site occupancy matrix)
! and then the lm dependent charge densities RHO and RHOAE
! (excluding augmentation charges yet)
!-----------------------------------------------------------------------
         RHOLM_=0
!      WRITE(*,'("RHOMIX",6F10.6)') RHOLM_STORE
         COCC=0
         ISIZE=UBOUND(RHOLM_STORE,1)

         DO ISP=1,NCDIJ
! retrieve the (1._q,0._q) center on site charge densities to RHOLM_
            IF ( LMAX_MIX < PP%LMAX_CALC) &
                 CALL TRANS_RHOLM( CRHODE(:,:,NIP,ISP), RHOLM_(:,ISP), PP )
! retrieve mixed elements from RHOLM_STORE and overwrite them in RHOLM_
            CALL RETRIEVE_RHOLM( RHOLM_(:,ISP), RHOLM_STORE(IBASE:,ISP), &
                           METRIC(IBASE:), IADD, PP, .FALSE.,  ITMP)

! calculate the total radial angular decomposed charge distributions
            LMMAX=(LYMAX+1)**2
            RHOAE(:,:,ISP)=0; RHO(:,:,ISP)=0
            CALL RAD_CHARGE( RHOAE(:,:,ISP), PP%R,RHOLM_(:,ISP), PP%LMAX, PP%LPS, PP%WAE )
            CALL RAD_CHARGE( RHO(:,:,ISP), PP%R, RHOLM_(:,ISP), PP%LMAX, PP%LPS, PP%WPS )
# 1473

! for LDA+U or Hartree fock we need the mixed occupancy matrix COCC
            IF (USELDApU() .OR. LDO_METAGGA() .OR. & 
           &   USEFOCK_CONTRIBUTION() .OR. USEFOCK_AE_ONECENTER() .OR. &
           &   LHYPERFINE()) THEN
! calculate the occupancy matrix COCC from RHOLM_(:,ISP)
               CALL TRANS_RHOLMI( COCC(:,:,ISP), RHOLM_(:,ISP), PP )
# 1495

            ENDIF
!              WRITE(*,*) 'spin',ISP
!              CALL DUMP_DLLMM("cocc",COCC(:,:,ISP),PP)
!              CALL DUMP_DLLMM("crhode",CRHODE(:,:,NIP,ISP),PP)
         ENDDO

! bring COCC to (up,down)
         IF (WDES%ISPIN==2) CALL OCC_FLIP2(COCC,LMDIM)

         CALL RAD_KINEDEN(PP,PP%WAE,LMAX_TAU,NCDIJ,COCC,TAUAE)
         CALL RAD_KINEDEN(PP,PP%WPS,LMAX_TAU,NCDIJ,COCC,TAUPS)

! bring COCC to spinor representation
         IF (WDES%LNONCOLLINEAR) CALL OCC_FLIP4(COCC,LMDIM)

!-----------------------------------------------------------------------
! calculation and output of radial kinetic energy density
! (remember augmentation charges are still excluded)
!-----------------------------------------------------------------------
         IF ( LMETA ) THEN
! RHOLM_ is the (total, magnetization) representation
            DO ISP=1,NCDIJ
               RHOLMUPD(:,ISP)=RHOLM_(:,ISP)
            ENDDO

            KINDENSAE=0; WKDAE=0
            KINDENSPS=0; WKDPS=0
            IF (ISPIN==2) THEN
! go to (spin up, down) representation
               CALL FLIP_RAD(RHOLMUPD,RHOLMUPD,LMDIM*LMDIM)
               CALL FLIP_RAD(RHOAE(:,1,1:2),RHOUPD(:,1,1:2),RNMAX)
            ELSE
! for ISPIN =1 and ISPIN = 4 stay in total, magnetization representation
               RHOUPD(:,1,1)=RHOAE(:,1,1)
            ENDIF

            DO ISP=1,ISPIN
               IF (ISPIN ==4 .AND. ISP==1) THEN 
               CALL RAD_KINETIC_EDENS(KINDENSAE(:,ISP),WKDAE(:,ISP),PP%R,RHOLMUPD(:,ISP), &
                    PP%LMAX,PP%LPS, PP%WAE, PP%RHOAE, PP%TAUAE, RHOUPD(:,:,ISP), 1)
               ELSE
               CALL RAD_KINETIC_EDENS(KINDENSAE(:,ISP),WKDAE(:,ISP),PP%R,RHOLMUPD(:,ISP), &
                    PP%LMAX,PP%LPS, PP%WAE, PP%RHOAE, PP%TAUAE, RHOUPD(:,:,ISP), ISPIN)
               ENDIF
            ENDDO

            IF (ISPIN==2) THEN
               CALL FLIP_RAD(RHO(:,1,1:2),RHOUPD(:,1,1:2),RNMAX)
            ELSE
               RHOUPD(:,1,1)=RHO(:,1,1)
            ENDIF
            DO ISP=1,ISPIN
               CALL RAD_KINETIC_EDENS(KINDENSPS(:,ISP),WKDPS(:,ISP),PP%R,RHOLMUPD(:,ISP), &
                    PP%LMAX,PP%LPS, PP%WPS, NULPOINTER, NULPOINTER, RHOUPD(:,:,ISP),ISPIN)
            ENDDO
         ENDIF
!-----------------------------------------------------------------------
! add augmentation charges now
!-----------------------------------------------------------------------
         DO ISP=1,NCDIJ
            CALL RAD_AUG_CHARGE(  RHO(:,:,ISP), PP%R, RHOLM_(:,ISP), PP%LMAX, PP%LPS,  &
                  LYMAX, PP%AUG, PP%QPAW )
            IF (ONE_CENTER_NMAX_FOCKAE()>0) THEN
               CALL RAD_AUG_CHARGE_FOCK(  RHO(:,:,ISP), PP%R, RHOLM_(:,ISP), PP%LMAX, PP%LPS, &
                    LYMAX, PP%AUG_FOCK, PP%QPAW_FOCK)
            ENDIF
            CALL RAD_INT( PP%R,  LYMAX, RHO(:,:,ISP), RHOAE(:,:,ISP) )
         ENDDO
         IBASE=IBASE+IADD
!-----------------------------------------------------------------------
! now finish the meta GGA stuff
!-----------------------------------------------------------------------
            IF ( LMETA ) THEN
! output of kinetic energy density
# 1599


! calculate E(xc) for meta-GGA
!               CALL RAD_META_GGA(PP%R, ISPIN, LYMAX, PP%LMAX_CALC, &
!                    RHO, PP%RHOPS, POT,KINDENSPS,WKDPS,EXCM)
               CALL RAD_META_GGA_ASPH(PP%R, ISPIN, LYMAX, PP%LMAX_CALC, &
                    RHO, PP%RHOPS, POT,KINDENSPS,WKDPS,EXCM)
               E%PAWPSM=E%PAWPSM-EXCM
!               CALL RAD_META_GGA(PP%R, ISPIN, LYMAX, PP%LMAX_CALC, &
!                    RHOAE, PP%RHOAE, POTAE,KINDENSAE,WKDAE,EXCM)
               CALL RAD_META_GGA_ASPH(PP%R, ISPIN, LYMAX, PP%LMAX_CALC, &
                    RHOAE, PP%RHOAE, POTAE,KINDENSAE,WKDAE,EXCM)
               E%PAWAEM=E%PAWAEM+EXCM-PP%DEXCCOREM
               E%PAWCOREM=E%PAWCOREM+PP%DEXCCOREM
# 1618

            ENDIF


         CALL EGRAD_EFG_RAD_HAR_ONLY(T_INFO,NI,PP,RHO,RHOAE)
!-----------------------------------------------------------------------
! calculate the local radial potential
! mind in the non-collinear case the potential V(r) = d E(r) / d rho (r)
! and the potential vec mu(r) = d E(r) / d vec m (r) are stored in
! POT and POTAE (potentials need to be real), whereas
! in the collinear case the spin up and down potentials
! are stored in POT and POTAE
! probably (1._q,0._q) should rewrite this in the collinear case
!-----------------------------------------------------------------------
! initialise the spin orbit contributions to D_ij to 0
         CSO=0
! Hartree Fock contribution set to (0._q,0._q)
         CHF=0

         IF ( WDES%LNONCOLLINEAR ) THEN
! bring KINDENSPS from spinor to 2 component spin up and spin down presentation

            CALL RAD_MAG_DENSITY_KINDENS(RHO, RHOCOL, TAUPS, KINDENSCOL, LYMAX, LMAX_TAU, PP%R)
# 1643

! do LDA+U instead of LSDA+U (set magnetisation density to 0)
            IF (L_NO_LSDA()) RHOCOL(:,:,2:WDES%NCDIJ)=0
# 1648

            CALL RAD_POT( PP%R, 2, LYMAX, PP%LMAX_CALC, LASPH,   &
               RHOCOL, PP%RHOPS, PP%POTPS, POT, DOUBLEPS, EXCG)

            CALL SET_CMBJ_ONE_CENTER_FACT(-1._q)
            CALL RAD_POT_METAGGA( PP%R, 2, PP%LMAX_CALC, LMAX_TAU, LASPH,  &
               RHOCOL, RHOCOL, PP%RHOPS, PP%POTPS, KINDENSCOL, POT, POT, DOUBLEPS, EXCG, MUPS, PP%TAUPS)

            E%PAWPSG=E%PAWPSG-EXCG
            CALL RAD_MAG_DIRECTION( RHO, RHOCOL, POT, LYMAX, PP%R)
# 1660

            CALL RAD_MAG_DIRECTION_KINDENS( RHO, TAUPS, LMAX_TAU, POT, MUPS, PP%R)

! bring KINDENSAE from spinor to 2 component spin up and spin down presentation
            CALL RAD_MAG_DENSITY_KINDENS(RHOAE, RHOCOL, TAUAE, KINDENSCOL, LYMAX, LMAX_TAU, PP%R)

! do LDA+U instead of LSDA+U (set magnetisation density to 0)
            IF (L_NO_LSDA()) RHOCOL(:,:,2:WDES%NCDIJ)=0

            CALL APPLY_ONE_CENTER_AEXX()
            CALL RAD_POT( PP%R, 2, LYMAX, PP%LMAX_CALC, LASPH,  &
                 RHOCOL, PP%RHOAE, PP%POTAE_XCUPDATED,  POTAE, DOUBLEAE,EXCG)

            CALL SET_CMBJ_ONE_CENTER_FACT( 1._q)
            CALL RAD_POT_METAGGA( PP%R, 2, PP%LMAX_CALC, LMAX_TAU, LASPH, &
                 RHOCOL, RHOCOL, PP%RHOAE, PP%POTAE_XCUPDATED, KINDENSCOL, POTAE, POTAE, DOUBLEAE, EXCG, MUAE, PP%TAUAE)

            CALL RESTORE_ONE_CENTER_AEXX

            E%PAWAEG=E%PAWAEG+EXCG-PP%DEXCCORE
            E%PAWCORE=E%PAWCORE+PP%DEXCCORE
            CALL RAD_MAG_DIRECTION( RHOAE, RHOCOL, POTAE, LYMAX, PP%R)

            CALL RAD_MAG_DIRECTION_KINDENS( RHOAE, TAUAE, LMAX_TAU, POTAE, MUAE, PP%R)

            IF (WDES%LSORBIT) &
              CALL SPINORB_STRENGTH(POTAE(:,1,1), PP%RHOAE, PP%POTAE_XCUPDATED, PP%R, CSO, &
                PP%LMAX, PP%LPS ,PP%WAE, PP%ZCORE+PP%ZVALF_ORIG, THETA=BETA, PHI=ALPHA)

         ELSE
! collinear case
            DRHOCORE=0

! do LDA+U instead of LSDA+U (set magnetisation density to 0)
            RHOCOL=RHO
            IF (L_NO_LSDA()) RHOCOL(:,:,2)=0
# 1698

! cl shifts DRHOCORE is (0._q,0._q) here, since for the pseudo terms
! we do not include the core electron in the exchange correlation term
! but only in the Hartree term
            CALL RAD_POT( PP%R, ISPIN, LYMAX, PP%LMAX_CALC, LASPH,   &
               RHOCOL, PP%RHOPS-DRHOCORE(1:RNMAX), PP%POTPS, POT, DOUBLEPS, EXCG)

            CALL SET_CMBJ_ONE_CENTER_FACT(-1._q)
            CALL RAD_POT_METAGGA( PP%R, ISPIN, PP%LMAX_CALC, LMAX_TAU, LASPH, &
               RHOCOL, RHOCOL, PP%RHOPS-DRHOCORE(1:RNMAX), PP%POTPS, TAUPS, POT, POT, DOUBLEPS, EXCG, MUPS, PP%TAUPS)

            CALL SET_CL_DRHOCORE_PS(DRHOCORE, NT, PP%R, PP%AUG)
            CALL ADD_CL_HARTREE_POT(DRHOCORE, NT, ISPIN, POT, PP%R)
!CALL RAD_POT_HAR_ONLY( PP%R, ISPIN, LYMAX, PP%LMAX_CALC,RHOCOL,  POT, DOUBLEPS)

            E%PAWPSG=E%PAWPSG-EXCG

            DRHOCORE=0
            
! do LDA+U instead of LSDA+U (set magnetisation density to 0)
            RHOCOL=RHOAE
            IF (L_NO_LSDA()) RHOCOL(:,:,2)=0
 
            CALL SET_CL_DRHOCORE_AE(DRHOCORE, NT, PP%RHOAE, PP%POTAE , PP%R, PP%ZCORE, PP%ZVALF_ORIG )

            CALL APPLY_ONE_CENTER_AEXX()

            CALL RAD_POT( PP%R, ISPIN, LYMAX, PP%LMAX_CALC, LASPH,   &
               RHOCOL, PP%RHOAE-DRHOCORE(1:RNMAX), PP%POTAE_XCUPDATED,  POTAE, DOUBLEAE,EXCG)

            CALL SET_CMBJ_ONE_CENTER_FACT( 1._q)
            CALL RAD_POT_METAGGA( PP%R, ISPIN, PP%LMAX_CALC, LMAX_TAU, LASPH, &
               RHOCOL, RHOCOL, PP%RHOAE-DRHOCORE(1:RNMAX), PP%POTAE_XCUPDATED, TAUAE, POTAE, POTAE, DOUBLEAE, EXCG, MUAE, PP%TAUAE)

            CALL RESTORE_ONE_CENTER_AEXX

            CALL ADD_CL_HARTREE_POT(DRHOCORE, NT, ISPIN, POTAE, PP%R)

!CALL RAD_POT_HAR_ONLY( PP%R, ISPIN, LYMAX, PP%LMAX_CALC,RHOCOL,  POTAE, DOUBLEAE)


            E%PAWAEG=E%PAWAEG+EXCG-PP%DEXCCORE
            E%PAWCORE=E%PAWCORE+PP%DEXCCORE
         ENDIF

         DOUBLEC_PS= DOUBLEC_PS-DOUBLEPS*T_INFO%VCA(NT)
! the core-core exchange correlation energy is included
! up to this point in \int dr (e_xc(rho_c+rho_v) - v_xc(rho_c+rho_v) rho_v(r)
! subtract it now
         DOUBLEC_AE= DOUBLEC_AE+DOUBLEAE*T_INFO%VCA(NT)-PP%DEXCCORE*T_INFO%VCA(NT)

         CALL HYPERFINE_RAD(T_INFO,NI,PP,RHO,RHOAE,POTAE,COCC)
!-----------------------------------------------------------------------
! calculate core level shift for averaged up and down potential
! or total potential (in non collinear case stored in POTAE(..,..,1))
!-----------------------------------------------------------------------
         DO RNMAX_CL=1,RNMAX
           IF (PP%RDEP>0 .AND.  PP%R%R(RNMAX_CL)-PP%RDEP > -5E-3) EXIT
         ENDDO

         IF (.NOT. LCOREL) THEN
           CL_SHIFT(1,NI)=0
           IF (NCDIJ==2) THEN
              CALL RAD_CL_SHIFT( (POT(1:RNMAX,1,1)+POT(1:RNMAX,1,2))/SPI2/2, &
                      (POTAE(1:RNMAX,1,1)+POTAE(1:RNMAX,1,2))/SPI2/2, PP%R, CL_SHIFT(1,NI), PP%AUG(:,0))
           ELSE
              CALL RAD_CL_SHIFT( POT(1:RNMAX,1,1)/SPI2, POTAE(1:RNMAX,1,1)/SPI2, PP%R, CL_SHIFT(1,NI), PP%AUG(:,0))
           ENDIF
         ELSE

           CALL CL_INIT_CORE_CONF(PP,MAXNL)
           ALLOCATE( W(RNMAX,MAXNL), N(MAXNL), LC(MAXNL), EIG(MAXNL))

! first version for cl-shifts
! ---------------------------
! calculate core wavefunctions in PAW sphere for atomic reference potential

           CL_SHIFT(1:MAXNL,NI) =0
           CALL SET_CORE_WF( PP%RHOAE, PP%POTAE , PP%R, PP%ZCORE, PP%ZVALF_ORIG , &
             W, N, LC, EIG)
           IF (MAXNL > SIZE(CL_SHIFT,DIM=1)) THEN
              WRITE(0,*) 'internal error: increase CL_MAXNL in main.F',MAXNL,SIZE(CL_SHIFT,DIM=1)
              CALL M_exit(); stop
           ENDIF

! now calculate the first order change caused by the current potential
! and subtract the pseudo contribution
           IF (NCDIJ==2) THEN
             CALL RAD_CL_SHIFT_AE( (PP%POTPSC-PP%POTPS)+ (POT(1:RNMAX,1,1)+POT(1:RNMAX,1,2))/2/SPI2 , &
                   PP%R, CL_SHIFT(:,NI), &
                   PP%AUG(:,0),  MAXNL, W, EIG, (POTAE(:,1,1)+POTAE(:,1,2))/SPI2/2)
           ELSE
             CALL RAD_CL_SHIFT_AE( (PP%POTPSC-PP%POTPS)+ POT(1:RNMAX,1,1)/SPI2 , &
                   PP%R, CL_SHIFT(:,NI), &
                   PP%AUG(:,0),  MAXNL, W, EIG, POTAE(:,1,1)/SPI2)
           ENDIF
!  WRITE(0,*) CL_SHIFT(1:MAXNL,NI)

! CL_SHIFT(1:MAXNL,NI) =0
! version for cl-shifts that uses exact potential
! ------------------------------------------------
! solve radial Schroedinger equation for *current* potential
! IF (NCDIJ==2) THEN
!    CALL SET_CORE_WF( PP%RHOAE, PP%POTAE , PP%R, PP%ZCORE, PP%ZVALF_ORIG , &
!    W, N, LC, EIG, (POTAE(1:RNMAX,1,1)+POTAE(1:RNMAX,1,2))/2 , NMAX= RNMAX_CL)

! subtract only the pseudo contribution
!    CALL RAD_CL_SHIFT_AE( (PP%POTPSC-PP%POTPS)+ (POT(1:RNMAX,1,1)+POT(1:RNMAX,1,2))/2/SPI2 , &
!         PP%R, CL_SHIFT(:,NI), &
!         PP%AUG(:,0),  MAXNL, W, EIG)
! ELSE
!    CALL SET_CORE_WF( PP%RHOAE, PP%POTAE , PP%R, PP%ZCORE, PP%ZVALF_ORIG , &
!    W, N, LC, EIG, POTAE(1:RNMAX,1,1),  NMAX= RNMAX_CL)

! subtract only the pseudo contribution
!    CALL RAD_CL_SHIFT_AE( (PP%POTPSC-PP%POTPS)+ POT(1:RNMAX,1,1)/SPI2 , &
!         PP%R, CL_SHIFT(:,NI), &
!         PP%AUG(:,0),  MAXNL, W, EIG)
! ENDIF
! WRITE(0,*)  CL_SHIFT(1:MAXNL,NI)
           DEALLOCATE(W, N, LC, EIG)
           CALL CL_CLEAR_CORE_CONF
         ENDIF

!-----------------------------------------------------------------------
! calculate the PAW correction terms to the pseudopotential strength D
! I have defined the PAW contribution in a way that in the limit of
! atomic occupancies no contributions are added
!-----------------------------------------------------------------------

! multiply potentials by simpson weights
         CALL RAD_POT_WEIGHT( PP%R, NCDIJ, LYMAX, POTAE)
         CALL RAD_POT_WEIGHT( PP%R, NCDIJ, LYMAX, POT)
# 1833

         CTMP=0
         DO ISP=1,NCDIJ
            DDLM=0
            CALL RAD_PROJ(  POTAE(:,:,ISP), PP%R, 1._q, DDLM, PP%LMAX, PP%LPS, PP%WAE )
            CALL RAD_PROJ(  POT(:,:,ISP)  , PP%R,-1._q, DDLM, PP%LMAX, PP%LPS, PP%WPS )

            CALL RAD_AUG_PROJ( POT(:,:,ISP), PP%R, DDLM, PP%LMAX, PP%LPS, &
                  LYMAX, PP%AUG, PP%QPAW )
# 1850

            IF (ONE_CENTER_NMAX_FOCKAE()>0) THEN
               CALL RAD_AUG_PROJ_FOCK( POT(:,:,ISP), PP%R, DDLM, PP%LMAX, PP%LPS, &
                    LYMAX, PP%AUG_FOCK, PP%QPAW_FOCK )
            ENDIF
! transform them using Clebsch Gordan coefficients and add to CDIJ
            CALL TRANS_DLM( CTMP(:,:,ISP), DDLM , PP )
         ENDDO

         
         CMETAGGA=0
         CALL RAD_PROJ_METAGGA(PP,PP%WAE,LMAX_TAU,MUAE, 1._q,CMETAGGA)
         CALL RAD_PROJ_METAGGA(PP,PP%WPS,LMAX_TAU,MUPS,-1._q,CMETAGGA)

! non-collinear case: strength parameters need to go to the spinor presentation now
         IF (WDES%LNONCOLLINEAR) THEN
            CALL DIJ_FLIP(CTMP,LMDIM)
            CALL DIJ_FLIP(CMETAGGA,LMDIM)
         ENDIF

         IF (USELDApU() .OR. USEFOCK_CONTRIBUTION() .OR. USEFOCK_AE_ONECENTER()) THEN
            IF (WDES%ISPIN==1.AND.(.NOT.WDES%LNONCOLLINEAR)) THEN
               COCC(:,:,1)=COCC(:,:,1)*0.5_q
! LDA+U requires up and down density
               COCC(:,:,2)=COCC(:,:,1)
            ENDIF

!---------------------------------------------------------------
            IF (USEFOCK_AE_ONECENTER()) THEN
               CALL SETUP_PAWFOCK_AE(NT, PP)  ! this initializes AE part only
               CALL CALC_PAWFOCK(NT, PP, COCC, CHF, DOUBLEC_HF)
               DOUBLEC_AE = DOUBLEC_AE + DOUBLEC_HF*T_INFO%VCA(NT)
            ELSE IF (USEFOCK_CONTRIBUTION()) THEN
               CALL CALC_PAWFOCK(NT, PP, COCC, CHF, DOUBLEC_HF)
               DOUBLEC_AE = DOUBLEC_AE + DOUBLEC_HF*T_INFO%VCA(NT)
            ENDIF
         
! correction terms from LDA+U
            IF (USELDApU()) THEN
!---------------------------------------------------------------
               CALL LDAPLUSU(LMDIM,NI,NT,COCC, CTMP, PP, DOUBLEC_LDAU)
               DOUBLEC_AE = DOUBLEC_AE + DOUBLEC_LDAU*T_INFO%VCA(NT)
            ENDIF
         ENDIF

         IF (LCALC_ORBITAL_MOMENT().AND.WDES%LNONCOLLINEAR) THEN
            COCC=CRHODE(:,:,NIP,:)
            CALL OCC_FLIP4(COCC,LMDIM) ! go to spinor representation
            CALL CALC_ORBITAL_MOMENT(LMDIM, NI, NT, COCC, PP, 0._q,0._q)
         ENDIF

         IF (WDES%LSORBIT) THEN
            COCC=CRHODE(:,:,NIP,:)
            CALL OCC_FLIP4(COCC,LMDIM) ! go to spinor representation
            CALL CALC_SPINORB_MATRIX_ELEMENTS(WDES,PP,T_INFO,NI,CSO,COCC)
         ENDIF

         CDIJ(:,:,NIP,:)=CDIJ(:,:,NIP,:)+(CTMP+CSO+CHF+CMETAGGA)*T_INFO%VCA(NT)

         DEALLOCATE(TAUAE,TAUPS,MUAE,MUPS,KINDENSCOL)
         CALL UNSET_RSGF_TYPE
! we can deallocate the calculated range-separated Greens functions here,
! which means that these kernels are recalculated for each ion
         CALL DEALLOCATE_RSGF
      ENDDO ion
! or we could deallocate here, which means that the Greens functions are
! only thrown away after SET_DD_PAW has finished completely
!     CALL DEALLOCATE_RSGF
!=======================================================================
! now distribute the DIJ to all nodes which hold DIJ (using global sum)
!=======================================================================

      CALL M_sum_d(WDES%COMM_INTER, CDIJ, LMDIM*LMDIM*WDES%NIONS*NCDIJ)
      CALL M_sum_d(WDES%COMM_KINTER,CDIJ, LMDIM*LMDIM*WDES%NIONS*NCDIJ)
# 1927

      CALL M_sum_d(WDES%COMM, CL_SHIFT, SIZE(CL_SHIFT))

      CALL M_sum_d(WDES%COMM, DOUBLEC_AE, 1)
      CALL M_sum_d(WDES%COMM, DOUBLEC_PS, 1)
      CALL M_sum_d(WDES%COMM, E%PAWPSG, 1)
      CALL M_sum_d(WDES%COMM, E%PAWAEG, 1)
      CALL M_sum_d(WDES%COMM, E%PAWCORE, 1)

      IF ( LMETA ) THEN
         CALL M_sum_d(WDES%COMM, E%PAWPSM, 1)
         CALL M_sum_d(WDES%COMM, E%PAWAEM, 1)
         CALL M_sum_d(WDES%COMM, E%PAWCOREM, 1)
      ENDIF
# 1952

      DEALLOCATE(KINDENSAE,WKDAE,WKDPS,KINDENSPS,RHOUPD,RHOLMUPD)
      DEALLOCATE(RHOCOL)
# 1957

      DEALLOCATE(POTAE, RHOAE, POT, RHO, DRHOCORE )
      DEALLOCATE(DDLM, RHOLM, RHOLM_)
      DEALLOCATE(CTMP,CSO,CHF,CMETAGGA)
      DEALLOCATE(COCC,COCC_IM)

      E%PAWAE=DOUBLEC_AE
      E%PAWPS=DOUBLEC_PS

      IF (LUSE_THOMAS_FERMI) CALL POP_XC_TYPE

      CALL RELEASE_PAWFOCK
    END SUBROUTINE SET_DD_PAW

  END MODULE pawm
