#include "symbol.inc"

MODULE mlr_main
  USE prec
  USE vaspxml
  USE lr_helper
  IMPLICIT NONE
CONTAINS

!*********************************************************************
!
! calculate second derivatives using linear response theory
!
! implemented by gK
! the main scheduler is based on the finite difference kernel of 
! Orest Dubay (at least the skeleton is essentially identical)
! it steps through all ions (at least those that are allowed to move)
! and calls the main linear response kerner to calculate the
! linear response of the wavefunction and the second derivatives
! with respect to the ionic positions
!
!*********************************************************************

  SUBROUTINE LR_SKELETON( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHAM, &
          IBRION,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,FORCE)

    USE base
    USE lattice
    USE finite_differences
    USE charge
    USE pseudo
    USE lattice
    
    USE nonl_high
    USE msymmetry
    USE mpimy
    USE mgrid
    USE mkpoints
    USE constant
    USE poscar
    USE wave
    USE pot
    USE subrot
    USE pawm
    USE rmm_diis
    USE choleski
    USE david
    USE wave_high
    USE mlrf_main
    USE lri_main
    USE subrot_cluster
    USE kpoints_change
    USE full_kpoints
    USE hamil_high
    USE pead
    USE fock
    USE morbitalmag
    USE meta
    USE mlwf, ONLY : WANNIER90
    USE wannier_interpolation
! solvation__
      USE solvation
! solvation__
    IMPLICIT NONE
!=======================================================================
!  structures
!=======================================================================
    TYPE (tau_handle)  KINEDEN
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavedes)     WDES
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W          ! wavefunction
    TYPE (latt)        LATT_CUR
    TYPE (dynamics)    DYN
    TYPE (info_struct) INFO
    TYPE (in_struct)   IO
    TYPE (mixing)      MIX
    TYPE (kpoints_struct) KPOINTS
    TYPE (symmetry)    SYMM
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
    TYPE (grid_3d)     GRIDB      ! Broyden grid
    TYPE (transit)     B_TO_C     ! index table between GRIDB and GRIDC
    TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
    TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
    TYPE (energy)      E
    TYPE (latt)        LATT_INI

    INTEGER IBRION
    INTEGER LMDIM,IRDMAX,IRDMAA,NEDOS
    REAL(q) TOTEN,EFERMI

    COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
    COMPLEX(q)  CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
    RGRID       DENCOR(GRIDC%RL%NP)           ! partial core
    COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
    COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

!  augmentation related quantities
    OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
             CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
             CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!  paw sphere charge density
    INTEGER N_MIX_PAW
    REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ),RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ)
!  charge-density and potential on soft grid
    COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
    RGRID       SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)
!  density of states
    REAL(q)    DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
!  Hamiltonian
    GDEF       CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
    REAL(q)   XCSIF(3,3)
    REAL(q) :: FORCE(3,T_INFO%NIONS)   ! forces in cartesian coordinates

! local variables related to finite difference code

    TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
    INTEGER :: NIONS
    INTEGER :: IU6                ! OUTCAR file
    INTEGER :: IU6K               ! OUTCAR k-points
    INTEGER :: IU0                ! stdout

    TYPE (type_info)   T_INFO_0
    INTEGER :: IUDYNMAT           ! DYNMAT file
    REAL(q) :: X
    REAL(q),ALLOCATABLE      :: INITIAL_FORCE(:,:)
    REAL(q),ALLOCATABLE      :: DISPL_FORCES(:,:,:)
    REAL(q),ALLOCATABLE      :: INT_STRAIN(:,:,:), PIEZO(:,:,:)
    REAL(q),ALLOCATABLE      :: BORN_CHARGES(:,:,:), BORN_CHARGES2(:,:,:)
    REAL(q),ALLOCATABLE      :: SECOND_DERIV(:,:)
    INTEGER                  :: DOF
    INTEGER                  :: PROCESSED_DISPL
    INTEGER                  :: I,J,K,M,N, IDIR
    REAL(q),ALLOCATABLE      :: WORK(:,:)
    REAL(q),ALLOCATABLE      :: EIGENVECTORS(:,:)
    REAL(q),ALLOCATABLE      :: EIGENVALUES(:)
    INTEGER                  :: IERROR
    REAL(q) :: EPSILON(3,3),  ELASTIC(3,3,3,3),  ELASTICP(6,3,3)
! variables for the reevaluation of energy
    REAL(q) DESUM, RMS       ! change
    REAL(q) TOTEN_           ! energy
    REAL(q) EDIFF            ! break condition
    INTEGER ICOUEV
    INTEGER NSIM
! linear response using symmetry
    INTEGER                  :: NKORIG
    REAL(q),ALLOCATABLE      :: D(:,:,:)
    INTEGER,ALLOCATABLE      :: ND(:), IDIRD(:,:)
    REAL(q)                  :: WORKD(3,3,T_INFO%NIONS)
    INTEGER                  :: IWORK(T_INFO%NIONS)
    CHARACTER(3)             :: IDIR_TEXT(3)=(/"x","y","z"/)
    LOGICAL                  :: LDO(T_INFO%NIONS)
    REAL(q)                  :: WDMAT(3,T_INFO%NIONS,3,T_INFO%NIONS),DMAT(3,3,T_INFO%NIONS,T_INFO%NIONS)
    REAL(q)                  :: ST(3,3,3,T_INFO%NIONS), AST(3,3,3,T_INFO%NIONS)
    REAL(q)                  :: FACT
! variables to store G [H,r] phi
    COMPLEX(qs), ALLOCATABLE :: RPHI(:,:,:,:,:)
    GDEFS, ALLOCATABLE ::  RPHI_CPROJ(:,:,:,:,:)
    TYPE (skpoints_trans)   :: KPOINTS_TRANS

! symmetry related quantities (common block)
    INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
    REAL(q)  GTRANS,AP
    COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
         GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

    ! NULLIFY to make sure a later "IF(ASSOCIATED())" works correctly.
    NULLIFY(KPOINTS_TRANS%CPHASE)

    IF (KINTER<0) THEN
       IF (LHFCALC) THEN
          CALL INTERPOLATE_BAND_STR(HAMILTONIAN, KPOINTS, GRID, LATT_CUR, LATT_INI, & 
            T_INFO, NONLR_S, NONL_S, W, NEDOS, DOS, DOSI, INFO%NELECT, INFO%NUP_DOWN,  &
            LMDIM, P, SV, CQIJ, CDIJ, SYMM, IO%IU0, IO%IU6)
       ELSE
          CALL INTERPOLATE_BAND_STR_DFT(HAMILTONIAN, E, KPOINTS, GRID, LATT_CUR, LATT_INI, & 
            T_INFO, NONLR_S, NONL_S, W, NEDOS, DOS, DOSI,  &
            LMDIM, P, SV, CQIJ, CDIJ, SYMM, INFO, IO, GRIDC, GRIDUS, C_TO_US, IRDMAX )
       ENDIF
       RETURN
    ENDIF

    IF (KINTER>0) THEN
#if defined(WANNIERINTERPOL) && ! defined(gammareal)
#ifdef VASP2WANNIER90
       IF (WANNIER90()) THEN
          CALL WANNIER_INTERPOLATE(HAMILTONIAN,KPOINTS,GRID,LATT_CUR,LATT_INI,&
               T_INFO,NONLR_S,NONL_S,W,LMDIM,P,SV,CQIJ,CDIJ,SYMM,INFO,IO)
       ENDIF
#endif
#endif
    ENDIF


    IF (LHFCALC) THEN
       CALL VTUTOR('E','LHFLINEARRESPONSE',0.0_q,1,0,1,(0.0_q,0.0_q),1,.TRUE.,1, &
            IO%IU6,3)
       CALL VTUTOR('S','LHFLINEARRESPONSE',0.0_q,1,0,1,(0.0_q,0.0_q),1,.TRUE.,1, &
            IO%IU0,3)
    ENDIF

    IF (LDO_METAGGA()) THEN
       CALL VTUTOR('E','METAGGARESPONSE',0.0_q,1,0,1,(0.0_q,0.0_q),1,.TRUE.,1, &
            IO%IU6,3)
       CALL VTUTOR('S','METAGGARESPONSE',0.0_q,1,0,1,(0.0_q,0.0_q),1,.TRUE.,1, &
            IO%IU0,3)
    ENDIF

    IF (INFO%LREAL) THEN
       CALL VTUTOR('W','LINEARRESPONSE LREAL',0.0_q,1,0,1,(0.0_q,0.0_q),1,.TRUE.,1, &
            IO%IU6,3)
       CALL VTUTOR('W','LINEARRESPONSE LREAL',0.0_q,1,0,1,(0.0_q,0.0_q),1,.TRUE.,1, &
            IO%IU0,3)
    ENDIF

    NIONS = T_INFO%NIONS
    IU6   = IO%IU6
    IU0   = IO%IU0

    IF (IO%NWRITE>=3) THEN
       IU6K  = IO%IU6
    ELSE
       IU6K  = -1
    ENDIF

    IF (IBRION>0) THEN
       IF (IBRION==8) THEN
!new
          ALLOCATE(D(3,3,T_INFO%NIONS), ND(T_INFO%NIONS),IDIRD(3,T_INFO%NIONS))

          CALL FREDOM(SYMM%ROTMAP,ISYMOP,INVMAP,NROTK,NPCELL,D,ND,1,T_INFO%NTYP,T_INFO%NIONS, &
             T_INFO%NITYP,WORKD,IWORK, & 
             LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3), &
             LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3),IDIRD)
          DOF=SUM(ND)

          IF (IU6>=0) THEN
             WRITE(IU6,*)
             WRITE(IU6,'(A,I5,A)') ' Found ',DOF,' degrees of freedom:'
             WRITE(IU6,*)' ----------------------------------------------'
             WRITE(IU6,*)
             DO J=1,NIONS
                IF (ND(J).GT.0) THEN
                   WRITE(IU6,'(A,I5,A,3A2)') '     directions for atom ',J,':  ',IDIR_TEXT(IDIRD(ND(J),J))
                ENDIF
             ENDDO
          ENDIF
!end new
       ELSE
          CALL COUNT_DOF(NIONS, T_INFO%LSFOR, T_INFO%LSDYN, DOF)
       ENDIF
    ELSE
       DOF=0
    ENDIF

    T_INFO_0=T_INFO
    NULLIFY(T_INFO_0%POSION)
    ALLOCATE(T_INFO_0%POSION(3,NIONS))
    ALLOCATE(INITIAL_FORCE(3,NIONS))
    ALLOCATE(DISPL_FORCES(DOF,3,NIONS),BORN_CHARGES(3,3,NIONS),BORN_CHARGES2(3,3,NIONS), & 
      PIEZO(3,3,3),INT_STRAIN(DOF,3,3))

    DISPL_FORCES=0
    BORN_CHARGES=0 ; BORN_CHARGES2=0
    PIEZO=0
    INT_STRAIN=0

    T_INFO_0%POSION             = T_INFO%POSION
    INITIAL_FORCE               = FORCE
    NKORIG                      = WDES%NKPTS



!=======================================================================
! reset the potential and recalculate the ground state wavefunctions
! with very high precision for this setup
!=======================================================================
    IF (IU0>=0) THEN
       WRITE (IU0,*) 'Linear response reoptimize wavefunctions to high precision'
       WRITE ( 17,*) 'Linear response reoptimize wavefunctions'
    ENDIF

    IF (IU6>=0) THEN
       WRITE (IU6,*) 'Linear response reoptimize wavefunctions'
    ENDIF

    IF (INFO%LREAL) THEN
       CALL RSPHER(GRID,NONLR_S,LATT_CUR)
    ENDIF

    CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES,  &
         INFO,P,T_INFO,E,LATT_CUR, &
         CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

!   CALL VECTORPOT(GRID, GRIDC, GRID_SOFT, SOFT_TO_C,  WDES%COMM_INTER, & 
!        LATT_CUR, T_INFO%POSION, HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT)
                  
    CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
         LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

    CALL SETDIJ_AVEC(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
         LMDIM,CDIJ,HAMILTONIAN%AVTOT, NONLR_S, NONL_S, IRDMAX)
    
    CALL SET_DD_MAGATOM(WDES, T_INFO, P, LMDIM, CDIJ) 
    CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
         E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )
    !
    ! only if the tolerance is very tight for the occupied states
    ! it is possible to attain a tight tolerance for the perturbed states
    ! with the present DIIS algorithm

    EDIFF = 1E-10
    NSIM=WDES%NSIM*2
#ifdef MPI
    NSIM=((WDES%NSIM*2+WDES%COMM_INTER%NCPU-1)/WDES%COMM_INTER%NCPU)*WDES%COMM_INTER%NCPU
#endif

    X=INFO%EBREAK

    ! for some non obvious reasons 
    IF (.NOT. LMAGBLOCH) THEN 
       INFO%EBREAK=0.25*EDIFF
    ELSE
       IF (IU0>=0) WRITE(IU0,*) 'WARNING: EBREAK remains at default, EDDAV has a problem for too tight convergence criteria'
       INFO%EBREAK=EDIFF
    ENDIF

    DO I=1,3
        CALL EDDAV(HAMILTONIAN,P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
             LMDIM,CDIJ,CQIJ, RMS,DESUM,ICOUEV, SV, E%EXHF, IO%IU6,IO%IU0, .FALSE., .TRUE., .FALSE.)

       E%EBANDSTR=BANDSTRUCTURE_ENERGY(WDES, W)
       TOTEN_=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF+Ediel_sol

       IF (IO%IU0>=0) THEN
          WRITE(IO%IU0,1000) I, TOTEN_, TOTEN_-TOTEN, DESUM, ICOUEV, RMS
          WRITE(17,1000)     I, TOTEN_, TOTEN_-TOTEN, DESUM, ICOUEV, RMS
       ENDIF
1000   FORMAT('DAV: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5,I6,'  ',E10.3)
       IF(ABS(DESUM)<EDIFF.AND.ABS(TOTEN_-TOTEN)<EDIFF) EXIT

       TOTEN=TOTEN_
    ENDDO
    CALL KPAR_SYNC_ALL(WDES,W)

    INFO%EBREAK=X    ! restore the break condition

    IF (IO%IU0>=0 .AND. IO%LOPEN) CALL WFORCE(17)

!=======================================================================
! determine G [H, r] |phi> = d/ dk | phi(k)>
!=======================================================================
    NULLIFY(DEG_CLUSTER)
    CALL FIND_DEG_CLUSTERS(WDES, W, DEG_CLUSTER)

    IF (LEPSILON .OR. KINTER/=0 .OR. LMAGBLOCH ) THEN

    CALL SET_NABIJ_AUG(P,T_INFO%NTYP)
    ALLOCATE(RPHI(WDES%NRPLWV,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN,3), &
             RPHI_CPROJ(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN,3))

    IF (LUSEPEAD()) THEN
       IF (LMAGBLOCH) THEN
          LPEAD_RETURN_Q_CPROJ=.FALSE.
       ENDIF
       RPHI=0
       RPHI_CPROJ=0
       CALL PEAD_DPSI_DK_ALL(W,KPOINTS,P,CQIJ,LATT_CUR,T_INFO,RPHI,RPHI_CPROJ)
       LPEAD_RETURN_Q_CPROJ=.TRUE.
    ELSE 
       RPHI=0
       RPHI_CPROJ=0

       DO IDIR=1,3
          IF (IU0>=0) THEN
             WRITE (IU0,*)'Linear response G [H, r] |phi>, progress :'
             WRITE (IU0,'(A,I3)') &
                  '  Direction: ',IDIR
             WRITE (17,*)'Linear response G [H, r] |phi>, progress :'
             WRITE (17,'(A,I3)') &
                  '  Direction: ',IDIR
          END IF

          IF (IU6>=0) THEN
             WRITE (IU6,*)'Linear response G [H, r] |phi>, progress :'
             WRITE (IU6,'(A,I3)') &
                  '  Direction: ',IDIR
          ENDIF

          CALL LRF_RPHI( &
             P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
             T_INFO,INFO,IO,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX, &
             CDIJ,CQIJ,SV,LMDIM,DEG_CLUSTER, IDIR, RPHI(:,:,:,:,IDIR), RPHI_CPROJ(:,:,:,:,IDIR), & 
             .TRUE.) !KINTER==0)

          IF (IO%LOPEN) CALL WFORCE(IO%IU6)
       ENDDO
    ENDIF

    IF (KINTER>0) THEN
!       CALL INTERPOLATE_BAND_STR(HAMILTONIAN, KPOINTS, GRID, LATT_CUR, LATT_INI, & 
!            T_INFO, NONLR_S, NONL_S, W, NEDOS, DOS, DOSI, INFO%NELECT, INFO%NUP_DOWN, &
!            LMDIM, P, SV, CQIJ, CDIJ, SYMM, IO%IU0, IO%IU6, RPHI)
       CALL INTERPOLATE_BAND_STR_DFT( HAMILTONIAN, E, KPOINTS, GRID, LATT_CUR, LATT_INI, & 
            T_INFO, NONLR_S, NONL_S, W, NEDOS, DOS, DOSI,  &
            LMDIM, P, SV, CQIJ, CDIJ, SYMM, INFO, IO, GRIDC, GRIDUS, C_TO_US, IRDMAX , RPHI)

       ! well this deallocates most of the quantities we have used, but not all
       DEALLOCATE(RPHI, RPHI_CPROJ)
       RETURN
    ENDIF

    IF (LMAGBLOCH) THEN

       CALL BLOCH_CURRENT( W, GRID_SOFT, GRIDC, GRIDUS, C_TO_US, SOFT_TO_C, P, LATT_CUR, & 
          HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT, CHTOT, NONLR_S, NONL_S, & 
          RPHI, RPHI_CPROJ, CDIJ, CQIJ, SV, EFERMI, &
          T_INFO, LMDIM, CRHODE, IRDMAX, IO%IU6, IO%IU0)

       ! well this deallocates most of the quantities we have used, but not all
       DEALLOCATE(RPHI, RPHI_CPROJ)
       RETURN
    ENDIF
!=======================================================================
! response to external excluding local field effects
! i.e. in the independent particle approximation
!=======================================================================
    DO IDIR=1,3

       IF (IU0>=0) THEN
          WRITE (IU0,*)'Linear response to external field (no local field effect), progress :'
          WRITE (IU0,'(A,I3)') &
               '  Direction: ',IDIR
          WRITE (17,*)'Linear response to external field (no local field effect), progress :'
          WRITE (17,'(A,I3)') &
               '  Direction: ',IDIR
       END IF
     
       IF (IU6>=0) THEN
          WRITE (IU6,*)'Linear response to external field (no local field effect), progress :'
          WRITE (IU6,'(A,I3)') &
               '  Direction: ',IDIR
       ENDIF

       IF (SYMM%ISYM>0) THEN
          DYN%VEL=0
          DYN%VEL(IDIR,:)=1.0
          CALL KARDIR(T_INFO%NIONS,DYN%VEL,LATT_CUR%B)
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
       ENDIF

       CALL LRF_MAIN( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,RPHI,RPHI_CPROJ, &
          LATT_CUR, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI, DEG_CLUSTER, KPOINTS_TRANS, EPSILON(IDIR,:), BORN_CHARGES(IDIR,:,:), PIEZO(IDIR,:,:),  & 
          IDIR, .FALSE.)

          IF (IO%LOPEN) CALL WFORCE(IO%IU6)
    ENDDO
    ! reinitialise symmetry
    IF (SYMM%ISYM>0) THEN
       DYN%VEL=0
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
    ENDIF

    IF (SYMM%ISYM>0) CALL TSYM(EPSILON,ISYMOP,NROTK,LATT_CUR%A)
    IF (IO%IU6>=0) THEN
       WRITE(IO%IU6,1110) '(INDEPENDENT PARTICLE, excluding Hartree and local field effects)',EPSILON
       WRITE(IO%IU6,130)
       CALL XML_TENSOR("epsilon_rpa",EPSILON)
    ENDIF

!=======================================================================
! response to external field including local field effects
! if do_kpoints_chage is not set the k-point set remains unchanged
! and symmetrization of the final tensor is performed
! results are not exact but errors are sometimes small (not recommended)
!=======================================================================
#define do_kpoints_change
#ifndef do_kpoints_change
    DO IDIR=1,3

       IF (IU0>=0) THEN
          WRITE (IU0,*)'Linear response to external field, progress :'
          WRITE (IU0,'(A,I3)') &
               '  Direction: ',IDIR
          WRITE (17,*)'Linear response to external field, progress :'
          WRITE (17,'(A,I3)') &
               '  Direction: ',IDIR
       END IF
     
       IF (IU6>=0) THEN
          WRITE (IU6,*)'Linear response to external field, progress :'
          WRITE (IU6,'(A,I3)') &
               '  Direction: ',IDIR
       ENDIF

       IF (SYMM%ISYM>0) THEN
          DYN%VEL=0
          DYN%VEL(IDIR,:)=1.0
          CALL  KARDIR(T_INFO%NIONS,DYN%VEL,LATT_CUR%B)
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
       ENDIF

       CALL LRF_MAIN( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,RPHI,RPHI_CPROJ, &
          LATT_CUR, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI, DEG_CLUSTER, KPOINTS_TRANS, EPSILON(IDIR,:), BORN_CHARGES(IDIR,:,:), PIEZO(IDIR,:,:),  & 
          IDIR, .TRUE.)

          IF (IO%LOPEN) CALL WFORCE(IO%IU6)
    ENDDO
    ! reinitialise symmetry
    IF (SYMM%ISYM>0) THEN
       DYN%VEL=0
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
    ENDIF

    IF (SYMM%ISYM>0) CALL TSYM(EPSILON,ISYMOP,NROTK,LATT_CUR%A)
    IF (IO%IU6>=0) THEN
       IF (LRPA) THEN
          WRITE(IO%IU6,1100) '(including local field effects in RPA (Hartree))',EPSILON
       ELSE
          WRITE(IO%IU6,1100) '(including local field effects in DFT)',EPSILON
       ENDIF
       WRITE(IO%IU6,130)
       CALL XML_TENSOR("epsilon",EPSILON)
    ENDIF

!=======================================================================
! response to external field including local field effects
! this version changes the k-point set on the fly
! and symmetrizes the polarisation vector
!=======================================================================
#else
    DO IDIR=1,3
       IF (IU0>=0) THEN
          WRITE (IU0,*)'Linear response to external field, progress :'
          WRITE (IU0,'(A,I3)') &
               '  Direction: ',IDIR
          WRITE (17,*)'Linear response to external field, progress :'
          WRITE (17,'(A,I3)') &
               '  Direction: ',IDIR
       END IF
     
       IF (IU6>=0) THEN
          WRITE (IU6,*)'Linear response to external field, progress :'
          WRITE (IU6,'(A,I3)') &
               '  Direction: ',IDIR
       ENDIF
       !
       ! reinitialise symmetry part for field in direction IDIR
       ! presently this is done by supplying a velocity field to the ions
       !
       IF (SYMM%ISYM>0) THEN
          DYN%VEL=0
          DYN%VEL(IDIR,:)=1.0
          CALL  KARDIR(T_INFO%NIONS,DYN%VEL,LATT_CUR%B)
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
#ifdef oldsym
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#else
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#endif
          CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR, KPOINTS_TRANS)

          ! Loewdin perturbation theory to improve states at added k-points
          CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
          &    LMDIM,CDIJ,CQIJ,4,SV,T_INFO,P,IO%IU0,DESUM, NKSTART=NKORIG+1)
          CALL KPAR_SYNC_ALL(WDES,W)

          CALL FIND_DEG_CLUSTERS(WDES, W, DEG_CLUSTER)
       ENDIF

       CALL LRF_MAIN( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,RPHI,RPHI_CPROJ, &
          LATT_CUR, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI, DEG_CLUSTER, KPOINTS_TRANS, EPSILON(IDIR,:), BORN_CHARGES(IDIR,:,:), PIEZO(IDIR,:,:), & 
          IDIR, .TRUE.)
       IF (IO%LOPEN) CALL WFORCE(IO%IU6)

       IF (SYMM%ISYM>0) CALL DEALLOCATE_KPOINTS_TRANS(KPOINTS_TRANS)
    ENDDO

    IF (SYMM%ISYM>0) THEN
       DYN%VEL=0
       CALL  KARDIR(T_INFO%NIONS,DYN%VEL,LATT_CUR%B)
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
       CALL TSYM(EPSILON,ISYMOP,NROTK,LATT_CUR%A)
    ENDIF

    IF (IO%IU6>=0) THEN
       IF (LRPA) THEN
          WRITE(IO%IU6,1100) '(including local field effects in RPA (Hartree))',EPSILON
       ELSE
          WRITE(IO%IU6,1100) '(including local field effects in DFT)',EPSILON
       ENDIF
       WRITE(IO%IU6,130)
       CALL XML_TENSOR("epsilon",EPSILON)
    ENDIF
#endif
1100   FORMAT(// &
            " MACROSCOPIC STATIC DIELECTRIC TENSOR ",A/, &
            " ------------------------------------------------------"/, &
            &       3(6X,3F13.6/), &
            " ------------------------------------------------------"/)

1110   FORMAT(// &
            " HEAD OF MICROSCOPIC STATIC DIELECTRIC TENSOR ",A/, &
            " ------------------------------------------------------"/, &
            &       3(6X,3F13.6/), &
            " ------------------------------------------------------"/)


    ENDIF
130 FORMAT (5X, //, &
     &'----------------------------------------------------', &
     &'----------------------------------------------------'//)

!=======================================================================
! ionic displacements
!=======================================================================

    IF (IU0>=0 .AND. DOF>0) THEN
       WRITE (IU0,*) 'Linear response DOF=',DOF
       WRITE (17,*) 'Linear response DOF=',DOF
    ENDIF

    IF (IU6>=0 .AND. DOF>0) THEN
       WRITE (IU6,*) 'Linear response:'
       WRITE (IU6,*) '  Degrees of freedom DOF   = ',DOF
    END IF

    DO PROCESSED_DISPL=1,DOF

       IF (IU0>=0) THEN
          WRITE (IU0,*)'Linear response progress:'
          WRITE (IU0,'(A,I3,A,I3)') &
               '  Degree of freedom: ',PROCESSED_DISPL,'/',DOF
          WRITE (17,*)'Linear response progress:'
          WRITE (17,'(A,I3,A,I3)') &
               '  Degree of freedom: ',PROCESSED_DISPL,'/',DOF
       END IF
     
       IF (IU6>=0) THEN
          WRITE (IU6,*)'Linear response progress:'
          WRITE (IU6,'(A,I3,A,I3)') &
               '  Degree of freedom: ',PROCESSED_DISPL,'/',DOF
       ENDIF
       T_INFO%POSION=T_INFO_0%POSION
       IF (DYN%IBRION==8) THEN
          CALL FIND_IJ_ID(.FALSE.,NIONS,PROCESSED_DISPL,ND,J,IDIR)
          IDIR=IDIRD(IDIR,J)
       ELSE
          CALL FIND_IJ(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,PROCESSED_DISPL,IDIR,J)
       ENDIF

       IF (SYMM%ISYM>0) THEN
          DYN%VEL=0
          DYN%VEL(IDIR,J)=1.0
          CALL  KARDIR(1,DYN%VEL(:,J),LATT_CUR%B)
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
#ifdef oldsym
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#else
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#endif
          CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR, KPOINTS_TRANS)

          ! Loewdin perturbation theory to improve states at added k-points
          CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
          &    LMDIM,CDIJ,CQIJ,4,SV,T_INFO,P,IO%IU0,DESUM, NKSTART=NKORIG+1)
          CALL KPAR_SYNC_ALL(WDES,W)

          ! or alternatively, but not more accurate 3 steps Davidson
          ! problem is that this might yield wrong piezoelectric tensors
!          X=INFO%EBREAK
!          INFO%EBREAK=0.25*EDIFF
!          DO I=1,3
!             CALL EDDAV(HAMILTONIAN,P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
!             LMDIM,CDIJ,CQIJ, RMS,DESUM,ICOUEV, SV, E%EXHF, IO%IU6,IO%IU0, .FALSE., .TRUE., .FALSE.)
!             
!             E%EBANDSTR=BANDSTRUCTURE_ENERGY(WDES, W)
!             TOTEN_=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF
!
!             IF (IO%IU0>=0) THEN
!                WRITE(IO%IU0,1000) I, TOTEN_, TOTEN_-TOTEN, DESUM, ICOUEV, RMS
!                WRITE(17,1000)     I, TOTEN_, TOTEN_-TOTEN, DESUM, ICOUEV, RMS
!             ENDIF
!             IF(ABS(DESUM)<EDIFF.AND.ABS(TOTEN_-TOTEN)<EDIFF) EXIT
!             TOTEN=TOTEN_
!          ENDDO
!          CALL KPAR_SYNC_ALL(WDES,W)
!          INFO%EBREAK=X         ! restore the break condition

          CALL FIND_DEG_CLUSTERS(WDES, W, DEG_CLUSTER)
       ENDIF

       CALL LR_MAIN( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
          T_INFO,T_INFO_0,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,DISPL_FORCES(PROCESSED_DISPL,:,:),INT_STRAIN(PROCESSED_DISPL,:,:), & 
          IDIR, J, BORN_CHARGES2(:,IDIR,J), &
          DEG_CLUSTER, KPOINTS_TRANS, RPHI, RPHI_CPROJ)
       IF (IO%LOPEN) CALL WFORCE(IO%IU6)

       IF (SYMM%ISYM>0) CALL DEALLOCATE_KPOINTS_TRANS(KPOINTS_TRANS)
    ENDDO

    CALL FREE_DEG_CLUSTERS(WDES,DEG_CLUSTER)

    IF (LEPSILON) THEN
       DEALLOCATE(RPHI, RPHI_CPROJ)
    ENDIF

    ! restore original symmetry
    IF (SYMM%ISYM>0) THEN
       DYN%VEL=0
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
#ifdef oldsym
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
            T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#else
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
            T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#endif
       CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
       CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
    ENDIF
!=======================================================================
!
! final processing and output 
! Born effective charges, piezoelectric tensors, and internal strain
!
!=======================================================================
    IF (DYN%IBRION==8) THEN
       !
       ! version for case symmetry was used to reduce the number of displacments
       !
       CALL PRINT_DYNMAT_ID(.FALSE.,NIONS,DOF,1.0_q,T_INFO%NTYP,T_INFO%NITYP,T_INFO%POMASS,DISPL_FORCES,D,ND,IU6)

       DO PROCESSED_DISPL=1,DOF
          CALL FIND_IJ_ID(.FALSE.,NIONS,PROCESSED_DISPL,ND,J,IDIR)
          DO N=1,T_INFO%NIONS
             DMAT(1:3,IDIR,N,J)=DISPL_FORCES(PROCESSED_DISPL,1:3,N)
          END DO
          ST(1:3,1:3,IDIR,J)=INT_STRAIN  (PROCESSED_DISPL,1:3,1:3)
       END DO

       CALL MKDMAT(SYMM%ROTMAP,ISYMOP,INVMAP,NROTK,NPCELL,D,DMAT,ST,ND,      &
            1,T_INFO%NTYP,T_INFO%NIONS,T_INFO%NITYP,WORKD,WDMAT,AST,IWORK,SYMM%TAU,SYMM%TAUROT,   &
            SYMM%WRKROT, &
            LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3), &
            LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3))

       LDO=.TRUE.
       CALL SYDMAT(DMAT,SYMM%ROTMAP,ISYMOP,NROTK,NPCELL,1,T_INFO%NTYP, &
            T_INFO%NIONS,T_INFO%NITYP,WDMAT,LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3),LDO)
       CALL STMAT(ST  ,SYMM%ROTMAP,ISYMOP,NROTK,NPCELL,1,T_INFO%NTYP, &
            T_INFO%NIONS,T_INFO%NITYP,AST  ,LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3),LDO)

       ! now reset DOF to NIONS*3
       ! and copy results from DMAT back to DISPL_FORCES

       CALL COUNT_DOF(NIONS, T_INFO%LSFOR, T_INFO%LSDYN, DOF)
       DEALLOCATE(DISPL_FORCES)
       ALLOCATE(DISPL_FORCES(DOF,3,NIONS))

       DO PROCESSED_DISPL=1,DOF
          CALL FIND_IJ(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,PROCESSED_DISPL,IDIR,J)
          DO N=1,NIONS
             DISPL_FORCES(PROCESSED_DISPL,1:3,N)=DMAT(1:3,IDIR,N,J)
          ENDDO
       ENDDO
    ELSE
       !
       ! version for case no symmetry was used
       !
       CALL PRINT_DYNMAT(NIONS,DOF,1.0_q,T_INFO%NTYP,T_INFO%NITYP,T_INFO%POMASS,DISPL_FORCES,T_INFO%LSDYN,T_INFO%LSFOR,IU6)

       DO PROCESSED_DISPL=1,DOF
          CALL FIND_IJ(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,PROCESSED_DISPL,IDIR,J)
          ST(1:3,1:3,IDIR,J)=INT_STRAIN  (PROCESSED_DISPL,1:3,1:3)
       END DO
    ENDIF

    IF (IU6>=0) THEN
       WRITE (IU6,130) 
       IF (LEPSILON) THEN
          ! for the sake of having everything at the very end print EPSILON again
          IF (LRPA) THEN
             WRITE(IO%IU6,1100) '(including local field effects in RPA (Hartree))',EPSILON
          ELSE
             WRITE(IO%IU6,1100) '(including local field effects in DFT)',EPSILON
          ENDIF

          FACT=EVTOJ*1E20_q/LATT_CUR%OMEGA
          
170       FORMAT(/ ' PIEZOELECTRIC TENSOR  for field in x, y, z        (e  Angst)',/ &
               10X,'XX', 10X,'YY', 10X,'ZZ',10X,'XY', 10X,'YZ', 10X,'ZX'/ &
               '  ----------------------------------------------------', &
               '----------------------------')

180       FORMAT(/ ' ',A,'  for field in x, y, z        (C/m^2)',/ &
               &        10X,'XX', 10X,'YY', 10X,'ZZ', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
               &        '  ----------------------------------------------------', &
               &        '----------------------------')
          IF (.NOT.LRPA) THEN
             WRITE (IU6,170)
             DO I =1,3
                WRITE (IU6,140) IDIR_TEXT(I),(PIEZO(I,J,J),J=1,3), & 
                     PIEZO(I,1,2),PIEZO(I,2,3),PIEZO(I,3,1)
             ENDDO
             
             CALL TSYM3(PIEZO,ISYMOP,NROTK,LATT_CUR%A)
          
             WRITE (IU6,180) 'PIEZOELECTRIC TENSOR'
             DO I =1,3
                WRITE (IU6,140) IDIR_TEXT(I),(PIEZO(I,J,J)*FACT,J=1,3), & 
                     PIEZO(I,1,2)*FACT,PIEZO(I,2,3)*FACT,PIEZO(I,3,1)*FACT
             ENDDO
140          FORMAT(2X,A1,6F12.5)
             
          
             WRITE (IU6,*)
             WRITE (IU6,*) 'BORN EFFECTIVE CHARGES (in e, cummulative output)'
             WRITE (IU6,*) '-------------------------------------------------'
             
             DO N=1,T_INFO%NIONS
                WRITE (IU6,'(" ion ",I4)') N
                DO IDIR =1,3
                   WRITE (IU6,'(I5,3F12.5)') IDIR,BORN_CHARGES(IDIR,:,N)
                ENDDO
             ENDDO

             CALL XML_BORN_CHARGES(BORN_CHARGES,T_INFO%NIONS)
             
          ENDIF
       ENDIF
    ENDIF
    IF (DOF>0 .AND. IU6>=0) THEN
       IF (DYN%ISIF>0) THEN
160       FORMAT(/ ' INTERNAL STRAIN TENSOR FOR ION ',I4,' for displacements in x,y,z  (eV/Angst):',/ &
               10X,'X', 11X,'Y', 11X,'Z', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
               '  ----------------------------------------------------', &
               '----------------------------')
          DO N=1,T_INFO%NIONS
             WRITE(IU6,160) N
             DO I =1,3
                WRITE (IU6,140) IDIR_TEXT(I),(ST(J,J,I,N),J=1,3),ST(1,2,I,N),ST(2,3,I,N),ST(3,1,I,N)
             ENDDO
          ENDDO
       END IF
       WRITE (IU6,130)
    ENDIF
!=======================================================================
!
! final processing and output
! vibrational frequencies
!
!=======================================================================
    IF (DOF>0 .AND. IU6>=0 ) THEN

       ALLOCATE(SECOND_DERIV(DOF,DOF))
       
       DO N=1,DOF
          CALL FIND_IJ(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,N,I,J)       
          DO M=1,DOF
             SECOND_DERIV(M,N)=DISPL_FORCES(M,I,J)
          END DO
       END DO
       
       IF (IU6>=0) THEN
          WRITE (IU6,*) 
          WRITE (IU6,*) 'SECOND DERIVATIVES (NOT SYMMETRIZED)'
          WRITE (IU6,*) '------------------------------------'
          CALL PRINT_SECOND_DERIV(NIONS,DOF,SECOND_DERIV,T_INFO%LSFOR,T_INFO%LSDYN,IU6)
       END IF
       
       DO N=1,DOF
          DO M=N+1,DOF
             X=0.5_q*(SECOND_DERIV(N,M)+SECOND_DERIV(M,N))
             SECOND_DERIV(N,M)=X
             SECOND_DERIV(M,N)=X
          END DO
       END DO

       ALLOCATE(WORK(DOF,32),EIGENVECTORS(DOF,DOF),EIGENVALUES(DOF))

       EIGENVECTORS=SECOND_DERIV
       N=1
       DO I=1,T_INFO%NTYP
          DO J=1,T_INFO%NITYP(I)
             DO K=1,3
                CALL FIND_DOF_INDEX(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,K,N,M)
                IF (M>0) EIGENVECTORS(:,M)=EIGENVECTORS(:,M)/SQRT(T_INFO%POMASS(I))
                IF (M>0) EIGENVECTORS(M,:)=EIGENVECTORS(M,:)/SQRT(T_INFO%POMASS(I))
             END DO
             N=N+1
          END DO
       END DO
       
       CALL XML_TAG("dynmat")
       CALL XML_VECARRAY("hessian")
       CALL XML_ARRAY_REAL(EIGENVECTORS)
       CALL XML_CLOSE_TAG

       CALL DSYEV &
            ('V','U',DOF,EIGENVECTORS,DOF, &
            EIGENVALUES,WORK,32*DOF, IERROR)
       IF (IERROR/=0) THEN
          IF (IU6>=0) THEN
             WRITE(IU6,*) "Error while diagonalisation DSYEV INFO=",IERROR
             WRITE(IU6,*) "Some of (or all) eigenvectors and eigenvalues are not correct !"
          END IF
       END IF

       CALL XML_VEC_REAL(EIGENVALUES,"eigenvalues",'(ES16.8)')
       CALL XML_VECARRAY("eigenvectors")
       CALL XML_ARRAY_REAL(EIGENVECTORS,'(ES16.8)')
       CALL XML_CLOSE_TAG
       CALL XML_CLOSE_TAG

       CALL PRINT_EIGENVECTORS(NIONS,DOF,T_INFO_0%POSION,LATT_CUR%A, &
            EIGENVECTORS,EIGENVALUES,       &
            T_INFO%LSFOR,T_INFO%LSDYN,IU6)

       N=1
       DO I=1,T_INFO%NTYP
          DO J=1,T_INFO%NITYP(I)
             DO K=1,3
                CALL FIND_DOF_INDEX(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,K,N,M)
                IF (M>0) EIGENVECTORS(M,:)=EIGENVECTORS(M,:)/SQRT(T_INFO%POMASS(I))
             END DO
             N=N+1
          END DO
       END DO

       IF (IU6>=0 .AND. IO%NWRITE>=3) THEN
          WRITE(IU6,*) "Eigenvectors after division by SQRT(mass)"
          CALL PRINT_EIGENVECTORS(NIONS,DOF,T_INFO_0%POSION,LATT_CUR%A, &
               EIGENVECTORS,EIGENVALUES,       &
               T_INFO%LSFOR,T_INFO%LSDYN,IU6)
       ENDIF
       DEALLOCATE(WORK, EIGENVECTORS, EIGENVALUES)
       
       IF (DYN%ISIF>0 .OR. LEPSILON) THEN
          WRITE(IO%IU6,130)
          ! invert the matrix of the second derivatives
          SECOND_DERIV=-SECOND_DERIV
          CALL INV_SECOND_DERIV(SECOND_DERIV, DOF, IU6 )

          ! ionic contribution to macroscopic dielectric tensor
          IF (LEPSILON .AND. IO%IU6>=0 .AND. .NOT. LRPA) THEN
             CALL EPSILON_ION( T_INFO, DOF, SECOND_DERIV, BORN_CHARGES, EPSILON )
             ! induced polariation -> field
             EPSILON=EPSILON*2*TPI/(LATT_CUR%OMEGA)*FELECT
             WRITE(IO%IU6,1100) 'IONIC CONTRIBUTION',EPSILON

             CALL XML_TENSOR("epsilon_ion",EPSILON)
          ENDIF

          IF (DYN%ISIF>0.AND.IO%IU6>=0) THEN
100       FORMAT(/ &
            A / &
            ' Direction', &
            4X,'XX', 10X,'YY', 10X,'ZZ', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
            ' --------------------------------------------------------------------------------'/ &
            ' XX     ',6F12.4/ &
            ' YY     ',6F12.4/ &
            ' ZZ     ',6F12.4/ &
            ' XY     ',6F12.4/ &
            ' YZ     ',6F12.4/ &
            ' ZX     ',6F12.4/ &
            ' --------------------------------------------------------------------------------'/)
             CALL ELASTIC_ION( T_INFO, DOF, SECOND_DERIV, ST, ELASTIC )

             ELASTICP(1,:,:)=ELASTIC(1,1,:,:)
             ELASTICP(2,:,:)=ELASTIC(2,2,:,:)
             ELASTICP(3,:,:)=ELASTIC(3,3,:,:)
             ELASTICP(4,:,:)=ELASTIC(1,2,:,:)
             ELASTICP(5,:,:)=ELASTIC(2,3,:,:)
             ELASTICP(6,:,:)=ELASTIC(3,1,:,:)

             FACT=EVTOJ*1E22_q/LATT_CUR%OMEGA

             WRITE(IU6,100) ' ELASTIC MODULI IONIC CONTR (kBar)', ( &
                  (ELASTICP(J,I,I)*FACT,I=1,3), &
                   ELASTICP(J,1,2)*FACT,ELASTICP(J,2,3)*FACT,ELASTICP(J,3,1)*FACT,J=1,6)
          END IF

          ! ionic contribution to piezoelectric tensor
          IF (LEPSILON .AND. DYN%ISIF>0.AND.IO%IU6>=0 .AND. .NOT. LRPA) THEN

             CALL PIEZO_ION( T_INFO, DOF, SECOND_DERIV, ST, BORN_CHARGES, PIEZO )
             
             WRITE (IU6,180) 'PIEZOELECTRIC TENSOR IONIC CONTR'
             FACT=EVTOJ*1E20_q/LATT_CUR%OMEGA
             DO I =1,3
                WRITE (IU6,140) IDIR_TEXT(I),(PIEZO(I,J,J)*FACT,J=1,3), & 
                     PIEZO(I,1,2)*FACT,PIEZO(I,2,3)*FACT,PIEZO(I,3,1)*FACT
             ENDDO
          ENDIF
       END IF
       DEALLOCATE(SECOND_DERIV)
    END IF
    IF (IO%IU6>=0) WRITE(IO%IU6,130)

    DEALLOCATE(T_INFO_0%POSION)
    DEALLOCATE(INITIAL_FORCE)
    DEALLOCATE(DISPL_FORCES,BORN_CHARGES,BORN_CHARGES2, &
      PIEZO,INT_STRAIN)

    IF (IBRION==8) THEN
      DEALLOCATE(D, ND, IDIRD)
    ENDIF


    IF (IU0>=0 .AND. DOF>0) THEN
       WRITE (IU0,*) 'Linear response finished'
       WRITE (17,*) 'Linear response finished'
    ENDIF

    IF (IO%LOPEN) CALL WFORCE(IO%IU6)

  END SUBROUTINE LR_SKELETON



!************************ SUBROUTINE INTERPOLATE_BAND_STR **************
!
! this routine interpolates the band structure to a dense k-point
! grid 
! the strategy is fairly simply
! first the first derivative with respect to k is calculated
!    (presently only the subspace of considered orbitals)
! next the k-point grid is partly shifted and a diagonalization
! in the subspace is performed
! very little data is stored, to keep the routine simple
! and concise
! 
!***********************************************************************

    SUBROUTINE INTERPOLATE_BAND_STR(HAMILTONIAN, KPOINTS, GRID, LATT_CUR, LATT_INI, & 
       &    T_INFO, NONLR_S, NONL_S, W, NEDOS, DOS, DOSI, NELECT, NUP_DOWN, &
       &    LMDIM, P, SV, CQIJ, CDIJ, SYMM, IU0, IU6, RPHI)
      USE prec
      USE wave_high
      USE lattice
      USE poscar
      USE mpimy
      USE mgrid
      USE nonl_high
      USE base
      USE pseudo
      USE kpoints_change
      USE constant
      USE choleski
      USE subrot
      USE mlrf_main
      USE hamil_high
      IMPLICIT NONE

      TYPE (ham_handle)  HAMILTONIAN
      TYPE (kpoints_struct) KPOINTS
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      TYPE (latt)        LATT_INI
      TYPE (type_info)   T_INFO
      TYPE (nonlr_struct)NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (potcar)      P(:)
      TYPE (wavespin), TARGET :: W
      INTEGER LMDIM
      OVERLAP         CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      OVERLAP         CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      TYPE (symmetry) SYMM
      RGRID           SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)
      INTEGER IU0, IU6
      COMPLEX(qs), OPTIONAL :: RPHI(:,:,:,:,:)
      INTEGER    NEDOS
      REAL(q)    DOS(NEDOS,W%WDES%ISPIN),DOSI(NEDOS,W%WDES%ISPIN), NELECT, NUP_DOWN
!    local
      INTEGER NKPTS_ORIG, I, ISP, NK, NPOS, NB
      INTEGER IKX, IKY, IKZ, IKXP, IKYP, IKZP, IWZ
      REAL(q) DISPL(3), DISPL_CART(3), DISPL_FOUND(3), DIST, DIST_FOUND, KSTART
      REAL(q) EFERMI, ENTROPY, PAR(1,1,1,1,W%WDES%NCDIJ),DOSPAR(1,1,1,W%WDES%NCDIJ)
      REAL(q)    DOS_TMP(NEDOS,W%WDES%ISPIN),DOSI_TMP(NEDOS,W%WDES%ISPIN), EFERMI_TMP
      INTEGER :: N
      REAL(q) :: EXHF
!-----------------------------------------------------------------------
! double k-points in new structure
! in the upper part of the wavefunction arrays the wavefunctions 
! corresponding to the new shifted k-points are stored
! in the lower part the original wavefunctions are stored
!-----------------------------------------------------------------------

      IF (KPOINTS%NKPX<1 .OR. KPOINTS%NKPY<1 .OR. KPOINTS%NKPZ<1 )  THEN
         ! here we need some warning
         ! the routine works only using regular gamma centered grids
         WRITE(*,*) ' INTERPOLATE_BAND_STR: requires a regular k-point grid '
         RETURN
      ENDIF

      CALL CHECK_FULL_KPOINTS

      ! original number of k-points
      NKPTS_ORIG=W%WDES%NKPTS
#ifdef oldsym
      CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LSORBIT.AND..NOT.W%WDES%LSPIRAL, &
           T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM, IU6, IU0, W%WDES%VKPT(:,1:NKPTS_ORIG))
#else
      CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LNONCOLLINEAR, &
           T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM, IU6, IU0, W%WDES%VKPT(:,1:NKPTS_ORIG))
#endif
      CALL KPAR_SYNC_ALL(W%WDES,W)
      CALL RE_GEN_LAYOUT( GRID, W%WDES, KPOINTS, LATT_CUR, LATT_INI, IU6, IU0)
      CALL REALLOCATE_WAVE( W, GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR)

      CALL XML_TAG("eigenvalues", comment="interpolated")
!-----------------------------------------------------------------------
! loop over new super k-point grid
! a  kind of Wigner Seitz procedure is used to generate the new supergrid
! this can be bypassed by setting IWZ to 0
! a few notes are in place here
! ) only gamma centered grids do not violate the symmetry
!   and are therefore recommended
! ) such super grids are created by setting KINTER to an odd value
! ) otherwise the super grids might spoil the symmetry
!-----------------------------------------------------------------------
      IWZ=2
      KSTART=-REAL(ABS(KINTER),q)/2+0.5_q
      DOS =0
      DOSI=0
      EFERMI=0

      DO IKX=0,ABS(KINTER)-1
      DO IKY=0,ABS(KINTER)-1
      DO IKZ=0,ABS(KINTER)-1

         ! search the equivalent k-point with the shortest length
         DIST_FOUND=1E6
         DO IKXP=-IWZ,IWZ
         DO IKYP=-IWZ,IWZ
         DO IKZP=-IWZ,IWZ
            DISPL(1)=((IKX+KSTART)/ABS(KINTER)+IKXP)/KPOINTS%NKPX
            DISPL(2)=((IKY+KSTART)/ABS(KINTER)+IKYP)/KPOINTS%NKPY
            DISPL(3)=((IKZ+KSTART)/ABS(KINTER)+IKZP)/KPOINTS%NKPZ
            ! k-point displacement in cartesian coordinates
            DISPL_CART=DISPL
            CALL DIRKAR(1, DISPL_CART(1), LATT_CUR%B)
            DIST=DISPL_CART(1)**2+DISPL_CART(2)**2+DISPL_CART(3)**2
            IF (DIST<DIST_FOUND) THEN
               DISPL_FOUND=DISPL
               DIST_FOUND =DIST
            ENDIF
         ENDDO
         ENDDO
         ENDDO
         
         DISPL=DISPL_FOUND
         DISPL_CART=DISPL
         CALL DIRKAR(1, DISPL_CART(1), LATT_CUR%B)

         DO NK=1,NKPTS_ORIG
            W%WDES%VKPT(:,NKPTS_ORIG+NK)=W%WDES%VKPT(:,NK)+DISPL
         ENDDO

         IF (PRESENT(RPHI)) THEN
            DO ISP=1,W%WDES%ISPIN
               DO NK=1,NKPTS_ORIG
                  DO N=1,W%WDES%NBANDS
                     W%CW(:,N,NKPTS_ORIG+NK,ISP)=W%CW(:,N,NK,ISP)+TPI*( &
                          RPHI(:,N,NK,ISP,1)*(0.0_q,1.0_q)*DISPL_CART(1)+ & 
                          RPHI(:,N,NK,ISP,2)*(0.0_q,1.0_q)*DISPL_CART(2)+ &
                          RPHI(:,N,NK,ISP,3)*(0.0_q,1.0_q)*DISPL_CART(3))
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO ISP=1,W%WDES%ISPIN
               DO NK=1,NKPTS_ORIG
                  DO N=1,W%WDES%NBANDS
                     W%CW(:,N,NKPTS_ORIG+NK,ISP)=W%CW(:,N,NK,ISP)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         CALL SET_DATAKE(W%WDES, LATT_CUR%B)
         IF (NONLR_S%LREAL) THEN
            CALL RSPHER(GRID,NONLR_S,LATT_CUR)
         ELSE
            CALL SPHER(GRID,NONL_S,P,W%WDES,LATT_CUR, 1)
         ENDIF
         CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)

         CALL ORTHCH(W%WDES,W, W%WDES%LOVERL, LMDIM, CQIJ)
         
         CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,W%WDES,SYMM, &
              LMDIM,CDIJ,CQIJ,2,SV,T_INFO,P,IU0,EXHF, & 
              NKSTART=NKPTS_ORIG+1)

         IF (KPOINTS%ISMEAR>=0) THEN
         ! well this is not elegant, we now simply set the weights for
         ! lower k-point set to 0, and that for the upper one to the
         ! proper weight
         KPOINTS%WTKPT(NKPTS_ORIG+1:NKPTS_ORIG*2)=KPOINTS%WTKPT(1:NKPTS_ORIG)
         KPOINTS%WTKPT(1:NKPTS_ORIG)=0

         ! sum dos
         CALL DENSTA( IU0, IU6, W%WDES, W, KPOINTS, NELECT, &
              NUP_DOWN, ENTROPY, EFERMI_TMP, KPOINTS%SIGMA, .TRUE., &
              NEDOS, 0, 0, DOS_TMP, DOSI_TMP, PAR, DOSPAR)

         ! restore weights
         KPOINTS%WTKPT(1:NKPTS_ORIG)=KPOINTS%WTKPT(NKPTS_ORIG+1:NKPTS_ORIG*2)
         KPOINTS%WTKPT(NKPTS_ORIG+1:NKPTS_ORIG*2)=0

         DOS =DOS +DOS_TMP *(1_q/REAL(KINTER,q)/REAL(KINTER,q)/REAL(KINTER,q))
         DOSI=DOSI+DOSI_TMP*(1_q/REAL(KINTER,q)/REAL(KINTER,q)/REAL(KINTER,q))
         EFERMI=EFERMI+EFERMI_TMP*(1_q/REAL(KINTER,q)/REAL(KINTER,q)/REAL(KINTER,q))
         ENDIF

         IF (IU6>=0) THEN
            
            DO ISP=1,W%WDES%ISPIN
               WRITE(IU6,'(/" k-point displacement            :",3F10.4 )') DISPL
               WRITE(IU6,'( " k-point displacement (cartesian):",3F10.4 )') DISPL_CART
               IF (W%WDES%ISPIN==2) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
               DO N=NKPTS_ORIG+1,W%WDES%NKPTS
                  WRITE(IU6,2201) N-NKPTS_ORIG, W%WDES%VKPT(1,N),W%WDES%VKPT(2,N),W%WDES%VKPT(3,N), &
                       &             (I,REAL( W%CELTOT(I,N,ISP) ,KIND=q) ,W%FERTOT(I,N,ISP)*W%WDES%RSPIN,I=1,W%WDES%NB_TOT)
               ENDDO
            ENDDO
2201        FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
                 &         '  band No.  band energies     occupation '/ &
                 &           (3X,I4,3X,F10.4,3X,F10.5))

            
         ENDIF

         CALL XML_VEC_REAL( DISPL, "displacement" )
         CALL XML_VEC_REAL( DISPL_CART, "displacement_cart" )

         CALL XML_EIGENVAL_NOHEAD( W%CELTOT(:,NKPTS_ORIG+1:W%WDES%NKPTS,:), W%FERTOT(:,NKPTS_ORIG+1:W%WDES%NKPTS,:), &
              W%WDES%NB_TOT, W%WDES%NKPTS-NKPTS_ORIG, W%WDES%ISPIN)
      
      ENDDO
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
! restore old data layout
!-----------------------------------------------------------------------
      CALL XML_CLOSE_TAG
#ifdef oldsym
      CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LSORBIT.AND..NOT.W%WDES%LSPIRAL, &
           T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,-1)
#else
      CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LNONCOLLINEAR, &
           T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,-1)
#endif
      CALL RE_GEN_LAYOUT( GRID, W%WDES, KPOINTS, LATT_CUR, LATT_INI, -1, -1)
      CALL REALLOCATE_WAVE( W, GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR)

      IF (KPOINTS%ISMEAR>=0) THEN
      CALL XML_DOS(EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, .FALSE., &
           DOS, DOSI, DOSPAR, NEDOS, 1, 1, W%WDES%NCDIJ, comment='interpolated')
      ENDIF


    END SUBROUTINE INTERPOLATE_BAND_STR




!************************ SUBROUTINE INTERPOLATE_BAND_STR_DFT **********
!
! This routine interpolates the band structure on a dense k-point
! grid. The strategy is fairly simple. First the first derivative with 
! respect to k is calculated (presently only the subspace of considered 
! orbitals).
! The routine is much more sophisticated than the routine above
! it first calculates the full k-point set and then performs
! the interpolation 
!
! This method currently only works with regular Gamma centered grids
! and the internal density of states calculation is only done for 
! ISMEAR/=-5
!
! TODO:
!  - collect all routines in a separate file (if possible)
!  - fix the density of state ISMEAR=-5 problem (if possible)
!  - consider if it is possible to map the shifted grids
!    into a new total CW array such that this could be used
!    for the optical routines (maybe all routines, by doing a
!    preconverger, then interpolation to obtain tighter wavefunctions,
!    then calling whatever routine based on this). Should give a general
!    speedup compared to the alternative (with loss of accuracy))
!  - add option to return arrays of velocities, eigenvalues and kpoints 
!  - add test to wether printout is needed (i.e. if values are returned
!    instead) 
! 
!
!***********************************************************************

    SUBROUTINE INTERPOLATE_BAND_STR_DFT(HAMILTONIAN, E, KPOINTS, GRID, LATT_CUR, LATT_INI, & 
       &    T_INFO, NONLR_S, NONL_S, W, NEDOS, DOS, DOSI, &
       &    LMDIM, P, SV, CQIJ, CDIJ, SYMM, INFO, IO, GRIDC, GRIDUS, C_TO_US, IRDMAX, RPHI)
      USE prec
      USE wave_high
      USE hamil_high
      USE lattice
      USE poscar
      USE mpimy
      USE mgrid
      USE nonl_high
      USE base
      USE pseudo
      USE kpoints_change
      USE constant
      USE choleski
      USE subrot
      USE subrot_cluster
      USE mlr_optic
      USE mlrf_main
      USE ini
      USE david
      IMPLICIT NONE

      TYPE (ham_handle)  HAMILTONIAN
      TYPE (energy)      E
      TYPE (kpoints_struct) KPOINTS
      TYPE (kpoints_struct) KPOINTS_INTER ! kpoint array to obtain IBZ of interpolated grid
      TYPE (skpoints_full),SAVE,POINTER ::  KPOINTS_INTER_FULL ! full interpolated grid
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      TYPE (latt)        LATT_INI
      TYPE (type_info)   T_INFO
      TYPE (nonlr_struct)NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (potcar)      P(:)
      TYPE (wavespin), TARGET :: W
      !TYPE (wavespin), TARGET :: W_OLD
      TYPE (wavespin)     :: W2     ! stores  the derivative of the wavefunction with respect to k_i
      TYPE (wavefun)      :: WTMP   ! temporary
      LOGICAL :: LDONE, LVELOCITY
      TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
      TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
      TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (info_struct) INFO
      TYPE (in_struct)   IO
      INTEGER LMDIM
      OVERLAP         CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      OVERLAP         CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      TYPE (symmetry) SYMM
      RGRID           SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)
      COMPLEX(qs), OPTIONAL :: RPHI(:,:,:,:,:)
      INTEGER    NEDOS
      REAL(q)    DOS(NEDOS,W%WDES%ISPIN),DOSI(NEDOS,W%WDES%ISPIN)
    ! local
      TYPE (wavedes)  :: WDES_INTER
      TYPE (wavespin) :: W_INTER

      INTEGER NKPTS_ORIG, I, L, K, J, ISP, NK, NPOS, NB, KPOINTS_INSIDE, DISPINDEX, IBZINDEX, IDIR, IRDMAX
      INTEGER IKX, IKY, IKZ, IKXP, IKYP, IKZP, IWZ, ISP_IRZ, IK
      REAL(q) DISPL(3), DISPL_CART(3), DIST, DIST_FOUND, KSTART, DISPL_FOUND(3,ABS(KINTER*KINTER*KINTER))
      REAL(q) EFERMI, ENTROPY, PAR(1,1,1,1,W%WDES%NCDIJ),DOSPAR(1,1,1,W%WDES%NCDIJ)
      INTEGER :: N, NSIM, NELM, ICOUEV
      REAL(q) :: EXHF, RMS, DESUM1, TOTEN
      LOGICAL, ALLOCATABLE :: K_FOUND_IN_GRID(:)
      REAL(q), ALLOCATABLE :: VELTEMP(:,:,:,:), VKPT_OLD(:,:), EINTERPOL_FULL(:,:,:,:),EINTERPOL(:,:,:,:), KINTERPOL(:,:) ! vectors to store eigenvalues, velocities and kpoints
      COMPLEX(q), ALLOCATABLE :: CW_OLD(:,:,:,:) ! temporary array for original W%CW
      REAL(q) :: S(3,3) ! mapping aray used to generate correct signs on the velocities
      REAL(q) DER(3) ! velocity temp array
!-----------------------------------------------------------------------
! Warnings and start up checks. This routine only works for regular
! Gamma centered grids. Implement more bulletproof tests ASAP!
!-----------------------------------------------------------------------
      IF (KPOINTS%NKPX<1 .OR. KPOINTS%NKPY<1 .OR. KPOINTS%NKPZ<1 )  THEN
         ! here we need some warning
         ! the routine works only using regular gamma centered grids
         WRITE(*,*) ' INTERPOLATE_BAND_STR: requires a regular k-point grid '
         RETURN
      ENDIF
      CALL CHECK_FULL_KPOINTS
      ! original number of k-points
      NKPTS_ORIG=W%WDES%NKPTS
      ! determine interpolated IBZ list for later exclusion of
      ! points outside IBZ, LINVERSION_IN is copied from the method call above
      ! as well as LNOSYM
      CALL START_TIMING("velo")
      CALL START_TIMING("timo")

      ! initialize things such as smearing etc.
      KPOINTS_INTER=KPOINTS
#ifdef oldsym
      CALL RD_KPOINTS_KINTER(KPOINTS_INTER, LATT_CUR, KINTER, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LSORBIT.AND..NOT.W%WDES%LSPIRAL, SYMM%ISYM<0, IO%IU6, IO%IU0)
#else
      CALL RD_KPOINTS_KINTER(KPOINTS_INTER, LATT_CUR, KINTER, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LNONCOLLINEAR, SYMM%ISYM<0, IO%IU6, IO%IU0)
#endif
      ! also, generate the full grid to get ISYMOP for later rotation
      ! of eigenvalues onto the full grid. This is not done unless 
      ! velocities are needed (saves significant computational time if
      ! only DOS is wanted)
      CALL STOP_TIMING("timo",IO%IU6,"POINTGEN_IBZ")

      CALL START_TIMING("timo")
      ALLOCATE(KPOINTS_INTER_FULL)
      CALL IBZKPT_HF(LATT_CUR, KPOINTS_INTER, KPOINTS_INTER_FULL, T_INFO%NIONS, SYMM%ROTMAP, SYMM%MAGROT,SYMM%ISYM, IO%IU6, IO%IU0)
      CALL STOP_TIMING("timo",IO%IU6,"POINTGEN_FULL")
!-----------------------------------------------------------------------
! precalculate all displacements
!-----------------------------------------------------------------------
      IWZ=2
      KSTART=-REAL(ABS(KINTER),q)/2+0.5_q
      ! for even meshes allow k-points at the boundary of the BZ
      IF ( MOD(KINTER,2)==0) THEN
         ! for even meshes there seems to be a problem
         IWZ=2
         KSTART=(-REAL(ABS(KINTER),q))/2
      ENDIF

      IK=0
      DO IKX=0,ABS(KINTER)-1
      DO IKY=0,ABS(KINTER)-1
      DO IKZ=0,ABS(KINTER)-1
         IK=IK+1
         ! search the equivalent k-point with the shortest distance
         DIST_FOUND=1E6
         DO IKXP=-IWZ,IWZ
         DO IKYP=-IWZ,IWZ
         DO IKZP=-IWZ,IWZ
            DISPL(1)=((IKX+KSTART)/ABS(KINTER)+IKXP)/KPOINTS%NKPX
            DISPL(2)=((IKY+KSTART)/ABS(KINTER)+IKYP)/KPOINTS%NKPY
            DISPL(3)=((IKZ+KSTART)/ABS(KINTER)+IKZP)/KPOINTS%NKPZ
            ! k-point displacement in cartesian coordinates
            DISPL_CART=DISPL
            CALL DIRKAR(1, DISPL_CART(1), LATT_CUR%B)
            DIST=DISPL_CART(1)**2+DISPL_CART(2)**2+DISPL_CART(3)**2
!test
            IF (DIST<=DIST_FOUND) THEN
!            IF (DIST<DIST_FOUND) THEN
               DISPL_FOUND(:,IK)=DISPL
               DIST_FOUND =DIST
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Create temporary vectors to store new kpoints and energies 
! symmetry layout is needed
!-----------------------------------------------------------------------
      ALLOCATE(EINTERPOL(4,W%WDES%NB_TOT,KPOINTS_INTER%NKPTS,W%WDES%ISPIN), & 
               KINTERPOL(3,KPOINTS_INTER%NKPTS), & 
               K_FOUND_IN_GRID(KPOINTS_INTER%NKPTS), &
               VKPT_OLD(3,NKPTS_ORIG), & 
               CW_OLD(SIZE(W%CW,1),W%WDES%NBANDS,NKPTS_ORIG, W%WDES%ISPIN), &
               VELTEMP(3,W%WDES%NB_TOT,NKPTS_ORIG,W%WDES%ISPIN))

      K_FOUND_IN_GRID=.FALSE.

      IF (LVEL) CALL ALLOCW(W%WDES,W2,WTMP,WTMP)

      ! Store old wave coefficients and kpoints
      DO ISP=1,W%WDES%ISPIN
         DO NK=1,NKPTS_ORIG
            IF (ISP==1) THEN
               VKPT_OLD(:,NK)=W%WDES%VKPT(:,NK)
            ENDIF
            DO N=1,W%WDES%NBANDS   
            CW_OLD(:,N,NK,ISP)=W%CW(:,N,NK,ISP)
            ENDDO
         ENDDO
      ENDDO
      ! initialize and set a couple of indices
      EINTERPOL=0.0
      KINTERPOL=0.0
      IBZINDEX=1
      DISPINDEX=1

!-----------------------------------------------------------------------
! now replace those k-points in the new dense grid (KPOINTS_INTER)
! that are not in the grid generated by shifting the original grid
! by symmetry equivalent k-points
!-----------------------------------------------------------------------
      ! first mark those k-points that are ok in KPOINTS_INTER
      DO IK=1,ABS(KINTER*KINTER*KINTER)
      DO N=1,NKPTS_ORIG
         NK=KPOINT_IN_FULL_GRID_KINTER(W%WDES%VKPT(:,N)+DISPL_FOUND(:,IK),KPOINTS_INTER)
         IF (NK/=-1) K_FOUND_IN_GRID(NK)=.TRUE.
      ENDDO
      ENDDO

      LDONE=.FALSE.

      ! now search for each k-point in the full BZ instead of IRZ
      ! and check whether this one maps onto a k-point in KPOINTS_INTER 
      ! which we have not yet found
      ! if so, replace that one in KPOINTS_INTER
      DO IK=1,ABS(KINTER*KINTER*KINTER)
      DO N=1,NKPTS_ORIG
         NK=KPOINT_IN_FULL_GRID_KINTER(W%WDES%VKPT(:,N)+DISPL_FOUND(:,IK),KPOINTS_INTER)
         IF (NK==-1) THEN
            ! must match to some k-point in the full BZ
            I=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,N)+DISPL_FOUND(:,IK), KPOINTS_INTER_FULL)
            ! equivalent one in IRZ
            NK=KPOINTS_INTER_FULL%NEQUIV(I)
            IF (.NOT. K_FOUND_IN_GRID(NK)) THEN
               ! replace that one in KPOINTS_INTER
               KPOINTS_INTER%VKPT(:,NK)=KPOINTS_INTER_FULL%VKPT(:,I)
               K_FOUND_IN_GRID(NK)=.TRUE.
               LDONE=.TRUE.
            ENDIF
         ENDIF
      ENDDO
      ENDDO

      ! final check: the shifted grid now hopefully generates all points in KPOINTS_INTER ...
      ! if not
      DO NK=1, KPOINTS_INTER%NKPTS
         IF (.NOT. K_FOUND_IN_GRID(NK) ) THEN
            IF (IO%IU0>=0) WRITE(IO%IU0,'(I6,3F14.7)') NK, KPOINTS_INTER%VKPT(:,NK)
            CALL VTUTOR('E','INTERPOLATE_K',0.0_q,1,0,1,(0.0_q,0.0_q),1,.TRUE.,1, &
               IO%IU0,3)
            CALL VTUTOR('S','INTERPOLATE_K',0.0_q,1,0,1,(0.0_q,0.0_q),1,.TRUE.,1, &
               IO%IU6,3)

            STOP
         ENDIF
      ENDDO

      ! unfortunately need to regenerate the IBZKPT_HF to get the proper equivalance list
      IF (LVEL .AND. LDONE) &
           CALL IBZKPT_HF(LATT_CUR, KPOINTS_INTER, KPOINTS_INTER_FULL, T_INFO%NIONS, SYMM%ROTMAP, SYMM%MAGROT,SYMM%ISYM, IO%IU6, IO%IU0)

      ! initialize WDES_INTER to fit KPOINTS_INTER
      WDES_INTER=W%WDES
      WDES_INTER%NKPTS= KPOINTS_INTER%NKPTS
      WDES_INTER%VKPT =>KPOINTS_INTER%VKPT
      WDES_INTER%WTKPT=>KPOINTS_INTER%WTKPT

      CALL ALLOCW_NOPLANEWAVE(WDES_INTER, W_INTER)
!-----------------------------------------------------------------------
! loop over new super k-point grid
! a  kind of Wigner Seitz procedure is used to generate the new supergrid
! this can be bypassed by setting IWZ to 0
! a few notes are in place here
! ) only gamma centered grids do not violate the symmetry
!   and are therefore recommended
! ) such super grids are created by setting KINTER to an odd value
! ) otherwise the super grids might spoil the symmetry
!-----------------------------------------------------------------------
      DOS =0
      DOSI=0
      EFERMI=0

      CALL START_TIMING("shifts")
      IF (IO%IU0>=0) WRITE(IO%IU0,*) ' '
      IF (IO%IU0>=0) WRITE(IO%IU0,*) '------------------------------------------------------'
      IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Starting shift of original grid'
      IF (IO%IU0>=0) WRITE(IO%IU0,*) '  original grid:', KPOINTS%NKPX, KPOINTS%NKPY, KPOINTS%NKPZ
      IF (IO%IU0>=0) WRITE(IO%IU0,*) '  new grid:     ', KPOINTS%NKPX*ABS(KINTER), KPOINTS%NKPY*ABS(KINTER), KPOINTS%NKPZ*ABS(KINTER)
      IF (IO%IU0>=0) WRITE(IO%IU0,*) '------------------------------------------------------'

      DO IK=1,ABS(KINTER*KINTER*KINTER)
         IF (IO%IU0>=0) THEN
            WRITE(IO%IU0,*) 'displacement', IK
            WRITE(17,*)     'displacement', IK
         ENDIF
         DISPL=DISPL_FOUND(:,IK)

         DISPL_CART=DISPL
         CALL DIRKAR(1, DISPL_CART(1), LATT_CUR%B)

         DO NK=1,NKPTS_ORIG
            W%WDES%VKPT(:,NK)=VKPT_OLD(:,NK)+DISPL
         ENDDO

         IF (KINTER>=0) THEN
            DO ISP=1,W%WDES%ISPIN
               DO NK=1,NKPTS_ORIG
                  DO N=1,W%WDES%NBANDS
                     W%CW(:,N,NK,ISP)=CW_OLD(:,N,NK,ISP)+TPI*( &
                          RPHI(:,N,NK,ISP,1)*(0.0_q,1.0_q)*DISPL_CART(1)+ & 
                          RPHI(:,N,NK,ISP,2)*(0.0_q,1.0_q)*DISPL_CART(2)+ &
                          RPHI(:,N,NK,ISP,3)*(0.0_q,1.0_q)*DISPL_CART(3))
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO ISP=1,W%WDES%ISPIN
               DO NK=1,NKPTS_ORIG
                  DO N=1,W%WDES%NBANDS
                     W%CW(:,N,NK,ISP)=CW_OLD(:,N,NK,ISP)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         CALL START_TIMING("timo")
         CALL SET_DATAKE(W%WDES, LATT_CUR%B)
         IF (NONLR_S%LREAL) THEN
            CALL RSPHER(GRID,NONLR_S,LATT_CUR)
         ELSE
            CALL SPHER(GRID,NONL_S,P,W%WDES,LATT_CUR, 1)
         ENDIF
         CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)

         CALL ORTHCH(W%WDES,W, W%WDES%LOVERL, LMDIM, CQIJ)        
!-----------------------------------------------------------------------
! gK: get accurate eigenvectors at least for LVEL
! but maybe always so maybe get rid of the (LVEL) ?
!-----------------------------------------------------------------------
         IF (.NOT. LINTERFAST) THEN
            ! Davidson 
            ! well I do not know whether this is required, but just in case
            NSIM=W%WDES%NSIM*2
#ifdef MPI
            NSIM=((W%WDES%NSIM*2+W%WDES%COMM_INTER%NCPU-1)/W%WDES%COMM_INTER%NCPU)*W%WDES%COMM_INTER%NCPU
#endif
            DO NELM=1,10
               CALL EDDAV(HAMILTONIAN,P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,W%WDES, NSIM, &
                    LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV,E%EXHF, IO%IU6,IO%IU0, .FALSE., .TRUE., .FALSE.)
               E%EBANDSTR=BANDSTRUCTURE_ENERGY(W%WDES, W)

               IF (IO%IU0>=0) WRITE(17, 200)      NELM,E%EBANDSTR,E%EBANDSTR-TOTEN,DESUM1,ICOUEV,RMS
               IF (IO%IU0>=0) WRITE(IO%IU0, 200)  NELM,E%EBANDSTR,E%EBANDSTR-TOTEN,DESUM1,ICOUEV,RMS

               TOTEN=E%EBANDSTR
           
200            FORMAT('DAV: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
                    &       I6,'  ',E10.3)
               
               IF (ABS(DESUM1) < ABS(INFO%EDIFF) ) EXIT
            ENDDO
         ELSE IF (LVEL) THEN
            CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,W%WDES,SYMM, &
                 LMDIM,CDIJ,CQIJ,3,SV,T_INFO,P,IO%IU0,EXHF)
         ELSE
            CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,W%WDES,SYMM, &
                 LMDIM,CDIJ,CQIJ,2,SV,T_INFO,P,IO%IU0,EXHF)
         ENDIF

         E%EBANDSTR=BANDSTRUCTURE_ENERGY(W%WDES, W)

        CALL STOP_TIMING("timo",IO%IU6,"DIAG")

!-----------------------------------------------------------------------       
! call routines to calculate the velocity
! ATTENTION: Not quite sure how ISPIN=2 works for these routines.
!            please use with caution
!-----------------------------------------------------------------------
         IF (LVEL) THEN
           IF (DISPINDEX==1) THEN
              IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Calculating velocities on the displaced fine grid:'
           ENDIF

           W2%CW   =0
           W2%CPROJ=0
           NULLIFY(DEG_CLUSTER)
           CALL FIND_DEG_CLUSTERS(W%WDES, W, DEG_CLUSTER)
           ! well is this really important ?
           ! IF (IO%IU0>=0) WRITE(IO%IU0,'(A,I4)')     'displacement ',DISPINDEX
           DO IDIR=1,3
              ! IF (IO%IU0>=0) WRITE(IO%IU0,'(A,I4)') '   direction ', IDIR
               ! these routines follow routines in linear_optics.F and 
               ! elinear_response.F and therein
               IF (LUSEPEAD()) THEN
                  W2%CW=0
                  W2%CPROJ=0
                  W2%CELTOT=0
                  CALL PEAD_DPSI_DK_IDIR(W,KPOINTS,P,CQIJ,LATT_CUR,T_INFO,IDIR,W2)
               ELSE
                  CALL FOCK_K_DER_ANALYT(KPOINTS, GRID, LATT_CUR, LATT_INI, T_INFO,  NONLR_S, NONL_S, W, W2, &
                       &   LMDIM, P, CQIJ, SYMM, IDIR, LDONE, IO%IU0, IO%IU6)
          
                  CALL LRF_RPHI0( &
                       P,NONLR_S,NONL_S,W,LATT_CUR, &
                       T_INFO,INFO,IO,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX, &
                       CDIJ,CQIJ,SV,LMDIM,DEG_CLUSTER, IDIR, W2, LDONE)
               ENDIF
               ! copy velocities into array
               DO ISP=1,W%WDES%ISPIN
                  DO NK=1,W%WDES%NKPTS
                     DO N=1,W%WDES%NB_TOT
                       VELTEMP(IDIR,N,NK,ISP)=REAL(W2%CELTOT(N,NK,ISP),KIND=q)
                    ENDDO
                 ENDDO
               ENDDO
           ENDDO
         ELSE
           IF (DISPINDEX==0) THEN
              IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Skipping velocity calculations'
           ENDIF
         ENDIF
!-----------------------------------------------------------------------
! Extract kpoints in the IBZ, all others are skipped. We have no use
! for points outside IBZ. Consider inlining this where the disp kpoints
! are calculated. This could be difficult in general (due to EDDIAG), 
! so this brute force method is probably the best 
!-----------------------------------------------------------------------
         DO N=1,NKPTS_ORIG
            NK=KPOINT_IN_FULL_GRID_KINTER(W%WDES%VKPT(:,N),KPOINTS_INTER)

            ! check if kpoints are inside IBZ, if not do not print out values
            IF (NK/=-1) THEN
               ! copy points to array so that a later printout in vasprun.xml is possible
               DO IDIR=1,3
                  KINTERPOL(IDIR,NK)=W%WDES%VKPT(IDIR,N)
               ENDDO
               DO ISP=1,W%WDES%ISPIN
                  DO NB=1,W%WDES%NB_TOT
                     EINTERPOL(1,NB,NK,ISP)=REAL(W%CELTOT(NB,N,ISP),KIND=q)
                     W_INTER%CELTOT(NB,NK,ISP)=W%CELTOT(NB,N,ISP)
                     DO IDIR=1,3
                        EINTERPOL(IDIR+1,NB,NK,ISP)=VELTEMP(IDIR,NB,N,ISP)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         DISPINDEX=DISPINDEX+1
      ENDDO
      CALL STOP_TIMING("shifts",IO%IU6,"SHIFTS")
!-----------------------------------------------------------------------
! Loop IBZ grid and print eigenvalues and velocities to the xml file
! consider to remove this in the future. Currently this is only for
! debuging
!-----------------------------------------------------------------------
      CALL XML_TAG("eigenvalues", comment="interpolated_ibz")
      IF (LVEL) CALL XML_TAG("velocities")
      CALL XML_TAG("kpoints_ibz")
      CALL XML_KPOINTS_LIST(KINTERPOL, KPOINTS_INTER%WTKPT)
      CALL XML_CLOSE_TAG("kpoints_ibz")
      IF (LVEL) THEN
         CALL XML_EIGENVALUES_EXT(EINTERPOL,4, W%WDES%NB_TOT, KPOINTS_INTER%NKPTS, W%WDES%ISPIN)
      ELSE
         CALL XML_EIGENVALUES_EXT(EINTERPOL,1, W%WDES%NB_TOT, KPOINTS_INTER%NKPTS, W%WDES%ISPIN)
      ENDIF
      IF (LVEL) CALL XML_CLOSE_TAG("velocities")
!-----------------------------------------------------------------------
! loop full grid and lay all eigenvalues and velocities on the 
! full grid and print to xml file. For the velocities the correct sign
! is needed. This is obtained by generating the transformation matrix S.
! In the future it should be sufficient to only write this sign array to 
! the xml file (this is integer and should occupy less memory, than the 
! full velocity grid) and then regenerate the sign during transport
! calculations.
!-----------------------------------------------------------------------
      IF (LVEL) THEN
        ALLOCATE(EINTERPOL_FULL(4,W%WDES%NB_TOT,KPOINTS_INTER_FULL%NKPTS,W%WDES%ISPIN))
        DO N=1,KPOINTS_INTER_FULL%NKPTS
           NK=KPOINTS_INTER_FULL%NEQUIV(N)
           ! generate transformation matrix to get proper sign on 
           ! the velocities
           S=0
           DO L=1,3
              DO K=1,3
                 DO J=1,3
                    DO I=1,3
                       S(L,I)=S(L,I)+LATT_CUR%A(L,K)*KPOINTS_INTER_FULL%ISYMOP(J,K,N)*LATT_CUR%B(I,J)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           DO ISP=1,W%WDES%ISPIN
              DO NB=1,W%WDES%NB_TOT
                 EINTERPOL_FULL(1,NB,N,ISP)=EINTERPOL(1,NB,NK,ISP)
                 DER=0
                 DO IDIR=1,3
                    DER(:)=DER(:)+S(:,IDIR)*EINTERPOL(IDIR+1,NB,NK,ISP)
                 ENDDO
                 DO IDIR=1,3
                    EINTERPOL_FULL(IDIR+1,NB,N,ISP)=DER(IDIR)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO 
        ! finally call xml routines to print values similar to old LVEL.
        ! should have a test wether these are wanted. Could be that 
        ! we would like to call this routine without printing velocities.
        ! Though, if this is the case, we should certainly return the
        ! velocities, eigenvalues and kpoints.
        CALL XML_TAG("electronvelocities")
        CALL XML_TAG("kpoints")
        CALL XML_KPOINTS_LIST(KPOINTS_INTER_FULL%VKPT, KPOINTS_INTER_FULL%WTKPT)
        CALL XML_CLOSE_TAG("kpoints")
        CALL XML_EIGENVALUES_EXT(EINTERPOL_FULL, 4, W%WDES%NB_TOT, KPOINTS_INTER_FULL%NKPTS, W%WDES%ISPIN)
        CALL XML_CLOSE_TAG("electronvelocities")
     ENDIF
!-----------------------------------------------------------------------
! cleanup, also cleanup KPOINTS_INTER, but it refuses as this
! does not have an ALLOCATE pointer
!----------------------------------------------------------------------

      CALL DENSTA( IO%IU0, IO%IU6, WDES_INTER, W_INTER, KPOINTS_INTER, INFO%NELECT, &
           INFO%NUP_DOWN, ENTROPY, EFERMI, KPOINTS_INTER%SIGMA, .TRUE., &
           NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
 
      E%EBANDSTR=BANDSTRUCTURE_ENERGY(WDES_INTER, W_INTER)

      TOTEN=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF

      IF (IO%IU0>=0) THEN
         WRITE(IO%IU0,1000) TOTEN
         WRITE(17,1000)     TOTEN
      ENDIF
1000  FORMAT('Energy from interpolated bandstructure ',E20.12)

      IF (LVEL) THEN
         DEALLOCATE(EINTERPOL,K_FOUND_IN_GRID, EINTERPOL_FULL, KINTERPOL, KPOINTS_INTER_FULL,VELTEMP)
      ELSE
         DEALLOCATE(EINTERPOL,K_FOUND_IN_GRID, KINTERPOL,VELTEMP)
      ENDIF

      CALL DEALLOCW(W_INTER)
!-----------------------------------------------------------------------
! restore old data layout
!-----------------------------------------------------------------------
      CALL XML_CLOSE_TAG
      ! restore old wave coefficients and kpoints
      DO ISP=1,W%WDES%ISPIN
         DO NK=1,NKPTS_ORIG
            IF (ISP==1) THEN
               W%WDES%VKPT(:,NK)=VKPT_OLD(:,NK)
            ENDIF
            DO N=1,W%WDES%NBANDS
               W%CW(:,N,NK,ISP)=CW_OLD(:,N,NK,ISP)
            ENDDO
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
! dump dos
!-----------------------------------------------------------------------
      CALL XML_DOS(EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, .FALSE., &
           DOS, DOSI, DOSPAR, NEDOS, 1, 1, W%WDES%NCDIJ, comment='interpolated')
      IF (IO%IU0>=0) WRITE(IO%IU0,*) '------------------------------------------------------'
      CALL STOP_TIMING("velo",IO%IU6,"KINTER")
    END SUBROUTINE INTERPOLATE_BAND_STR_DFT


END MODULE mlr_main
