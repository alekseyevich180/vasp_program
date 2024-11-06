# 1 "ilinear_response.F"
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

# 2 "ilinear_response.F" 2 

MODULE lri_main

CONTAINS
!*********************************************************************
!
! calculate second derivatives using linear response theory
! implemented by gK
! this is the main subroutine the direction and the index of the ion
! to be displaced must be supplied
! the second derivatives are returned in the array F
!
!*********************************************************************

  SUBROUTINE LR_MAIN( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W0,LATT_CUR, &
          T_INFO,T_INFO_0,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,F,SIF,IDIR, ION, BORN_CHARGE, &
          DEG_CLUSTER, KPOINTS_TRANS, RPHI, RPHI_CPROJ)

    USE base
    USE lattice
    USE finite_differences
    USE charge
    USE pseudo
    USE lattice
    USE force    
    USE nonl_high
    USE msymmetry
    USE mpimy
    USE mgrid
    USE mkpoints
    USE constant
    USE poscar
    USE wave
    USE pot
    USE ini
    USE pawm
    USE hamil_lr
    USE rmm_diis
    USE rmm_diis_lr
    USE subrot_lr
    USE charge
    USE us
    USE wave_high
    USE choleski
    USE broyden
    USE subrot
    USE ebs
    USE mlrf_main
    USE subrot_cluster
    USE kpoints_change
    USE hamil_high
    USE meta
! solvation__
      USE solvation
! solvation__
    IMPLICIT NONE
!=======================================================================
!  structures
!=======================================================================
    TYPE (tau_handle)  KINEDEN
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (type_info)   T_INFO, T_INFO_0
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavedes)     WDES
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W0         ! unperturbed wavefunctions
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
    
    INTEGER LMDIM,IRDMAX,NEDOS
    REAL(q) TOTEN,EFERMI

    COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
    COMPLEX(q)  CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
    REAL(q)       DENCOR(GRIDC%RL%NP)           ! partial core
    COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
    COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

!  augmentation related quantities
    REAL(q)  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
             CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
             CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!  paw sphere charge density
    INTEGER N_MIX_PAW
    REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ),RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ)
!  charge-density and potential on soft grid
    COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
    REAL(q)       SV(GRID%MPLWV*2,WDES%NCDIJ)
!  density of states
    REAL(q)    DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
! local l-projected wavefunction characters (not really used here)
    REAL(q)    PAR(1,1,1,1,WDES%NCDIJ),DOSPAR(1,1,1,WDES%NCDIJ)
!  Hamiltonian
    REAL(q) :: F(:,:)                    ! force constants
    REAL(q) :: SIF(:,:)                  ! internal strain tensor
    INTEGER :: IDIR, ION                 ! direction and ion index
    REAL(q) :: BORN_CHARGE(3)            ! Born effective charges
! G [H,r] phi
    COMPLEX(qs) :: RPHI(:,:,:,:,:)
    REAL(qs)       :: RPHI_CPROJ(:,:,:,:,:)
    TYPE (eigenf_cluster_pointer) :: DEG_CLUSTER(WDES%NKPTS,WDES%ISPIN)
    TYPE (skpoints_trans)   :: KPOINTS_TRANS
!   finite difference stencils
    REAL(q) :: STENCIL2_X(2)=(/0.5_q,-0.5_q/)
    REAL(q) :: STENCIL2_W(2)=(/1.0_q,-1.0_q/)
    REAL(q) :: STENCIL4_X(4)=(/1.0_q,0.5_q,-0.5_q,-1.0_q/)
    REAL(q) :: STENCIL4_W(4)=(/-2._q/12._q,16._q/12._q,-16._q/12._q,2._q/12._q/)
! local variables
    INTEGER :: IONODE, NODE_ME, NN
    INTEGER :: IRDMAA
    REAL(q) :: DESUM,DE                  ! change of energy
    REAL(q) :: CSHIFT=0.002              ! complex shift
! CSHIFT shifts the eigenvalue e(0) by a small complex shift c when calculating
! |xi> = H(1) - e(1) S(0) - (e(0) +c) S(1) |phi(0)>
! this implies that the change in the overlap operator is put into the
! "complex part" of |xi>
! it is the most elegant method to retrieve <phi(0)_k| S(1) |phi(0)_j>,
! the change of overlap between two states, in the
! limit of two degenerate eigenvalue pairs (see subrot_lr.F)
! alternatively 1._q must guarantee that degenerated initial  states |phi(0)>
! belong to different symmetry class
! this can be achived by calling EDDIAG_LR with RESOLVE_DEG=.TRUE. in the
! first few iterations
! presently both methods are applied simultaneously
! therefore setting CSHIFT to 0, should not change results
! another method would be to calculate  <phi(0)_k| S(1) |phi(0)_j> explicitly
! and store it seperately
! Note that for metals, degenerate states at the Fermi-level must be
! correctly resolved since their occupancies change
! gK: decreased CSHIFT to 0.002, oddly there seems to be issues
! with larger shifts. Most likely some odd problem in LINEAR_RESPONSE_DIIS
! when shifts are used
    REAL(q) :: RMS(INFO%NELM)            ! magnitude of residual vector
    INTEGER :: ICOUEV                    ! number of H | phi> evaluations
    INTEGER :: N                         ! main loop counter (electronic step)
    INTEGER :: NORDER                    ! order in the MP smearing
    REAL(q) :: TOTENL=0                  ! total energy in previous step
    COMPLEX(q),ALLOCATABLE :: CHTOT1(:,:)! derivative change of total charge
    COMPLEX(q),ALLOCATABLE :: CHTOT2(:,:)! derivative change of total charge
    COMPLEX(q),ALLOCATABLE :: CVTOT1(:,:)! derivative of local potential
    COMPLEX(q),ALLOCATABLE :: CWORK (:,:)! derivative of local potential
    COMPLEX(q),ALLOCATABLE:: CSTRF1(:,:) ! structure-factors for displaced ions
    COMPLEX(q),ALLOCATABLE:: CSTRF2(:,:) ! structure-factor
    COMPLEX(q),ALLOCATABLE:: CSTRF3(:,:) ! structure-factor
    COMPLEX(q),ALLOCATABLE:: CSTRF4(:,:) ! structure-factor
    REAL(q),ALLOCATABLE :: CDIJ1(:,:,:,:)! derivative of non local strength
    REAL(q),ALLOCATABLE::CRHODE1(:,:,:,:)! derivative of onsite occupancies
    REAL(q),ALLOCATABLE::CTMP(:,:,:,:)   ! temporary
    REAL(q),ALLOCATABLE::CRHODE2(:,:,:)  ! derivative for the ion moved
    REAL(q),ALLOCATABLE   :: SV1(:,:)      ! derivative of soft potential
    REAL(q),ALLOCATABLE   :: SV2(:,:)      ! work array
    COMPLEX(q),ALLOCATABLE:: CHDEN1(:,:) ! 1st order change of pseudo charge
    REAL(q)  RHOLM1(N_MIX_PAW,WDES%NCDIJ)! 1._q center occupancy matrix (mixed)

    TYPE (wavespin)     :: WXI           ! stores H(1) - epsilon S(1) | phi_0>
    TYPE (wavespin)     :: W1            ! first order change of wavefunction
    TYPE (wavefun)      :: WTMP          ! temporary
    TYPE (energy)          E1,E2         ! energy arrays
    INTEGER I                            ! loop counter for various cases
    INTEGER NIP
    INTEGER ISP
    REAL(q)             :: POTIM
    REAL(q)             :: RMST,RMSC,RMSP,WEIGHT,STEP ! mixer
    REAL(q)             :: Z(3)
    INTEGER             :: IERRBR
    REAL(q)             :: SCALE=0.02_q
    REAL(q), EXTERNAL :: TR_POT_CHARGE
    INTEGER :: ierror
! forces and stress tensors
    REAL(q)   XCSIF(3,3),EWSIF(3,3),TSIF(3,3),PRESS,VTMP(3)
    REAL(q)   EWIFOR(3,T_INFO%NIOND),TIFOR(3,T_INFO%NIOND), DISPL(3,T_INFO%NIONS)

    IONODE=0
    NODE_ME=0

    IONODE  = WDES%COMM%IONODE
    NODE_ME = WDES%COMM%NODE_ME

! finite differences for second derivative of potential
! too small values cause convergence problems (noise upon
! calculation of POT_DER)
! larger ones (1E-2) smooth the energy surface but increase frequencies
!   POTIM=5E-3_q
! gK 02.05.2012: SiO2 ELASTIC MODULI show convergence only around POTIM=2E-3
! added fourth order stencil in POT_DER4, which converges already
! at POTIM=1E-2_q


    POTIM=1E-2_q ! pretty anything between 1E-2 and 2E-3 will do
# 215

!=======================================================================
!  allocate all required work arrays
!=======================================================================
    CALL REINIT_DEG_CLUSTERS(WDES,DEG_CLUSTER)

    ALLOCATE(CHTOT1(GRIDC%MPLWV,WDES%NCDIJ), &
             CHTOT2(GRIDC%MPLWV,WDES%NCDIJ), &
             CVTOT1(GRIDC%MPLWV,WDES%NCDIJ), &
             CSTRF1(GRIDC%MPLWV,T_INFO%NTYP), &
             CSTRF2(GRIDC%MPLWV,T_INFO%NTYP), &
             CWORK(GRIDC%MPLWV,WDES%NCDIJ), &
             CDIJ1(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
             CRHODE1(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
             CTMP(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
             CRHODE2(LMDIM,LMDIM,WDES%NCDIJ), &
             SV1(GRID%MPLWV*2,WDES%NCDIJ), &
             SV2(GRID%MPLWV*2,WDES%NCDIJ), &
             CHDEN1(GRID_SOFT%MPLWV,WDES%NCDIJ))

    ALLOCATE(CSTRF3(GRIDC%MPLWV,T_INFO%NTYP), &
             CSTRF4(GRIDC%MPLWV,T_INFO%NTYP))


    CALL ALLOCW(WDES,WXI,WTMP,WTMP)
    CALL ALLOCW(WDES,W1,WTMP,WTMP)

    W1%CPTWFP   =0
    W1%GPROJ=0
!=======================================================================
! setup structure factors
!=======================================================================
    CALL DISPLACE_ATOM(T_INFO_0%POSION, T_INFO%POSION, ION, IDIR, LATT_CUR%A, LATT_CUR%B, (POTIM/2))
    CALL STUFAK(GRIDC,T_INFO,CSTRF1)

    CALL DISPLACE_ATOM(T_INFO_0%POSION, T_INFO%POSION, ION, IDIR, LATT_CUR%A, LATT_CUR%B, -(POTIM/2))
    CALL STUFAK(GRIDC,T_INFO,CSTRF2)


    CALL DISPLACE_ATOM(T_INFO_0%POSION, T_INFO%POSION, ION, IDIR, LATT_CUR%A, LATT_CUR%B, (POTIM))
    CALL STUFAK(GRIDC,T_INFO,CSTRF3)

    CALL DISPLACE_ATOM(T_INFO_0%POSION, T_INFO%POSION, ION, IDIR, LATT_CUR%A, LATT_CUR%B, -(POTIM))
    CALL STUFAK(GRIDC,T_INFO,CSTRF4)


    CALL DISPLACE_ATOM(T_INFO_0%POSION, T_INFO%POSION, ION, IDIR, LATT_CUR%A, LATT_CUR%B, 0._q)
!=======================================================================
!
!  calculate
!  ) the first order change of the 1._q-centre occupancies
!   (CRHODE2) and the charge density for fixed wavefunction coefficients
!   (CHTOT2) i.e. the change due to the update of the wavefunction
!   coefficients is not included
!  ) and the derivative of the wave function coefficients for a fixed
!    potential
!
!=======================================================================
    SV1=0
    CDIJ1=0

    CALL LR_HAMIL(GRID,INFO,LATT_CUR,NONLR_S,NONL_S, W0, WXI, WDES, &
         LMDIM,CDIJ1, CDIJ, CQIJ, CRHODE2, ION, IDIR,  SV1, CSHIFT, RMS(1), ICOUEV, .FALSE.)

    CALL MRG_FERWE(WDES,WXI)

    IF (INFO%LOVERL) THEN
! derivative of augmentation charge with respect to the moved ion
! rho(0)_ij  d Q_ij(r)/dR for the single displaced ion
       CALL AUG_DER( &
            WDES, GRIDC, GRIDUS, C_TO_US, &
            LATT_CUR, P, T_INFO, SYMM, INFO%LOVERL, &
            LMDIM, CRHODE, CHTOT2, IRDMAX, ION, IDIR)

       CHDEN1=0
       CRHODE1=0
! contribution
! < phi(0) | d/dR p_i> D(0)_ij - e(0) Q(0)_ij  <p_j | phi(0) >
       
       NIP=NI_LOCAL(ION, WDES%COMM_INB)  ! local storage index in CRHODE1
       IF (NIP/=0) CRHODE1(:,:,NIP,:)=CRHODE1(:,:,NIP,:)+CRHODE2
! change storage convention to (total, magnetization)
       
       CALL US_FLIP(WDES, LMDIM, CRHODE1, INFO%LOVERL, .FALSE.)
! add first order change of 1._q center occupancy matrix times
! augmentation charge distribution rho(1)_ij Q_ij(r)
       CALL DEPLE_ADD(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
            LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
            LMDIM, CRHODE1, CHTOT2, CHDEN1, IRDMAX)

       CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
            CRHODE1, RHOLM1)

! symmetrise total pseudo charge density CHTOT2 if required
       IF (SYMM%ISYM ==1) THEN
          IF (WDES%LNONCOLLINEAR) THEN
             CALL RHOSYM(CHTOT2(1,1),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,1)
             IF (.NOT.WDES%LSPIRAL) &
                  &   CALL SYMFIELD(CHTOT2(1,2),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,WDES%SAXIS,LATT_CUR)
          ELSE
             DO ISP=1,WDES%ISPIN
                CALL RHOSYM(CHTOT2(1,ISP),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,ISP)
             ENDDO
          ENDIF
       ENDIF
    ELSE
       CHTOT2=0
    ENDIF

!=======================================================================
! 1st derivative of charge in the first iteration is given by
! the derivative of the atomic charge density stored in CHTOT1
!=======================================================================
    CRHODE1=0
    RHOLM1 =0
    CHTOT1 =0
    CALL CHARGEDER(GRIDC,P,T_INFO,LATT_CUR,CHTOT1(1,1),.FALSE.,ION,IDIR)

# 348


! to make timing more sensefull syncronize now
    CALL MPI_barrier( WDES%COMM%MPI_COMM, ierror )
    CALL START_TIMING("LOOP")

    IF (IO%IU0>=0) THEN
       WRITE(IO%IU0,142)
       WRITE(17,142)
    ENDIF
142 FORMAT('       N       E                     dE             ' &
         ,'d eps       ncg     rms          rms(c)')

    INFO%LMIX=.FALSE.

130 FORMAT (5X, //, &
     &'----------------------------------------------------', &
     &'----------------------------------------------------'//)

140 FORMAT (5X, //, &
     &'----------------------------------------- Iteration ', &
     &I4,'(',I4,')  ---------------------------------------'//)
    ! 'electron entered'

!     CALL DIPOL_RESET()
!=======================================================================
!  main selfconsistent loop for linear response
!=======================================================================
    MIX%LRESET=.TRUE.
    TOTEN=0


    INFO%LCHCON=.FALSE.
electron: DO N=1,INFO%NELM

       IF(IO%IU6>=0) WRITE(IO%IU6,140) ION*3+IDIR,N
       CALL XML_TAG("scstep")

!=======================================================================
!
!  calculate first derivative of potential (CVTOT1, SV1)
!  and derivatives of non local strenght (CDIJ1)
!
!=======================================================================
       CALL START_TIMING("G")

! positive displacement
       CALL DISPLACE_ATOM(T_INFO_0%POSION, T_INFO%POSION, ION, IDIR, LATT_CUR%A, LATT_CUR%B, 0.0_q)
       DISPL=0 ; DISPL(IDIR,ION)=1


       CALL POT_DER4( WDES, GRID, GRIDC, GRIDUS, GRID_SOFT, C_TO_US, SOFT_TO_C, &
            INFO, T_INFO, P,  LATT_CUR, E1, E2, CSTRF, CSTRF1, CSTRF2, CSTRF3, CSTRF4, &
            CHTOT, CHTOT1, CHTOTL, DENCOR, CVTOT1, CWORK, SV1, SV2, & 
            LMDIM, CDIJ1, CRHODE, CRHODE1, CQIJ, CTMP, &
            N_MIX_PAW, RHOLM, RHOLM1, RHOLM_LAST,  &
            IRDMAX, DISPL, POTIM)
# 412


! second order change of double counting corrections
! E1%DENC = - Tr [ rho(1) V_H rho(1) ]/2
! see below for the definition of rho_c(1)
       E1%DENC =(E1%DENC + E2%DENC -2*E%DENC )*(1/POTIM)**2*2
       E1%XCENC=(E1%XCENC+ E2%XCENC-2*E%XCENC)*(1/POTIM)**2*2
       E1%PAWPS=(E1%PAWPS+ E2%PAWPS-2*E%PAWPS)*(1/POTIM)**2*2
       E1%PAWAE=(E1%PAWAE+ E2%PAWAE-2*E%PAWAE)*(1/POTIM)**2*2

       CALL STOP_TIMING("G",IO%IU6,"POT+DIJ")
!=======================================================================
!
!  calculate |xi> = H(1) - e(0) S(1) - e(1) S(0) | phi(0)>
!  where H(1) and S(1) are the first order change of the
!  Hamiltonian and the overlap operator respectively and
!  e(1) is the first order change of the eigenenergy
!  e(1) is evaluated as well and stored in W1%CELTOT
!
!=======================================================================

       CALL LR_HAMIL(GRID,INFO,LATT_CUR,NONLR_S,NONL_S, W0, WXI, WDES, &
             LMDIM,CDIJ1, CDIJ, CQIJ, CRHODE2, ION, IDIR,  SV1, CSHIFT, RMS(N), ICOUEV, .TRUE.)
    
       CALL MRG_FERWE(WDES,WXI)
       CALL MRG_CEL(WDES,WXI)

# 464


       CALL STOP_TIMING("G",IO%IU6,"HAMIL1")

!
! the double counting corrections are
!  - Tr [ rho(1)-rho_c(1) V_Hxc  rho(1)-rho_c(1) ]
! where rho(1) is the total 1st order change of the charge
! and rho_c(1) is the change of the charge density for fixed
! wavefunction coefficients
! the d.c calculated above are however  Tr [ rho(1) V_Hxc  rho(1) ]
! adding Tr [  (rho(1) V_Hxc rho_c(1)) ] should work, but seems to be
! too inaccurate
!       E1%DENC = E1%DENC+E1%XCENC+TR_POT_CHARGE( GRIDC, WDES%NCDIJ, CVTOT1, CHTOT2 )
! test
!       CALL POTHAR(GRIDC,LATT_CUR, CHTOT1,CWORK,E2%DENC)
!       WRITE(*,*) E2%DENC ,E1%DENC,E2%DENC-E1%DENC
!       CALL POTXC2(GRIDC,LATT_CUR, CHTOT, DENCOR, CHTOT1,CWORK,E2%XCENC)
!       WRITE(*,*) E2%XCENC,E1%XCENC,E2%XCENC-E1%XCENC
! Hartree and LDA contribution
       CALL POTHAR(GRIDC,LATT_CUR, CHTOT1-CHTOT2,CWORK,E1%DENC)
       CALL POTXC2(GRIDC,LATT_CUR, CHTOT, DENCOR, CHTOT1-CHTOT2,CWORK,E1%XCENC)
!=======================================================================
!
!  calculate the first order change of the wavefunction
!    H(0) - epsilon S(0) | phi_1> = - |xi>
!  the first order change of the energy is returned in WXI%CELEN
!  the first order change of the norm in WXI%FERWE
!
!=======================================================================
       CALL EDDIAG_LR(W1,W0,WXI,LMDIM,CQIJ,INFO%LOVERL,.TRUE.,CSHIFT,IO%IU0,DEG_CLUSTER)
       CALL STOP_TIMING("G",IO%IU6,"LRDIAG")

       CALL LINEAR_RESPONSE_DIIS(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W1,WXI,W0,WDES, &
            LMDIM,CDIJ,CQIJ, RMS(N),DESUM,ICOUEV, SV, CSHIFT, IO%IU6,IO%IU0, LRESET=(N==1))

       CALL MRG_CEL(WDES,W1)
       CALL STOP_TIMING("G",IO%IU6,"LRDIIS")

       IF (N<=RESOLVE_DEG_NTIMES) THEN
          CALL EDDIAG_LR(W1,W0,WXI,LMDIM,CQIJ,INFO%LOVERL,.TRUE.,CSHIFT,IO%IU0,DEG_CLUSTER,.TRUE.)
       ELSE
          CALL EDDIAG_LR(W1,W0,WXI,LMDIM,CQIJ,INFO%LOVERL,.TRUE.,CSHIFT,IO%IU0,DEG_CLUSTER)
       ENDIF

       CALL STOP_TIMING("G",IO%IU6,"LRDIAG")

! sum f(0) (<phi(1)|xi> + c.c + <phi(1)| H(0)|phi(1)>)
       W1%FERTOT=W0%FERTOT  ! set the occupancy matrix W1 to W0
       E1%EBANDSTR= BANDSTRUCTURE_ENERGY(WDES, W1)

! normalize phi(1) by adding - <phi(0)| S(1) |phi(0)> phi(0) /2
! not required if complex shift method is used
       IF (CSHIFT==0) THEN
          CALL  NORM_LR( W1, W0, WXI, WDES)
! subtract e(1) <phi(0)| S(1) |phi(0)> f(0)
          WXI%FERTOT= WXI%FERTOT*W0%FERTOT ! set WXI%FERWE to <phi(0)| S(1) |phi(0)> f(0)
          E1%EBANDSTR= E1%EBANDSTR - BANDSTRUCTURE_ENERGY(WDES, WXI)
       ENDIF
       
# 546

# 603


!=======================================================================
! calculate first derivative of 1._q electron occupancies
!=======================================================================

! occupancies (FERWE) use finite differences
! this could be 1._q analytically for some methods
! but not for all, so what the hell, precision is good enough

       W1%CELTOT=W0%CELTOT+WXI%CELTOT*(POTIM/2)
       CALL DENSTA( IO%IU0, IO%IU6, WDES, W1, KPOINTS, INFO%NELECT, &
            INFO%NUP_DOWN, E1%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
            NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
       WXI%FERTOT=W1%FERTOT

       W1%CELTOT=W0%CELTOT-WXI%CELTOT*(POTIM/2)
       CALL DENSTA( IO%IU0, IO%IU6, WDES, W1, KPOINTS, INFO%NELECT, &
           INFO%NUP_DOWN, E2%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
           NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
! finished, store the derivative in W1
! and set the eigenvalue array (until this point stored in WXI)
       W1%CELTOT=   WXI%CELTOT
       W1%FERTOT  =(WXI%FERTOT-W1%FERTOT)   *(1/POTIM)
       E1%EENTROPY=(E1%EENTROPY-E2%EENTROPY)*(1/POTIM)

       TOTENL=TOTEN
       TOTEN =E1%EBANDSTR+E1%DENC+E1%XCENC+E1%PAWPS+E1%PAWAE+E1%EENTROPY+Ediel_sol
       DE= (TOTEN-TOTENL)
       
       IF (IO%IU0>=0) THEN
          WRITE(IO%IU0,200,ADVANCE='NO') N, TOTEN*SCALE**2, DE*SCALE**2, DESUM*SCALE**2, ICOUEV, RMS(N)*SCALE
          WRITE(17,200,ADVANCE='NO')     N, TOTEN*SCALE**2, DE*SCALE**2, DESUM*SCALE**2, ICOUEV, RMS(N)*SCALE
       ENDIF
 200   FORMAT('RMM: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5,I6,'  ',E10.3)
!=======================================================================
! calculate first derivative of charge density
!=======================================================================

    chg:IF (.NOT. INFO%LCHCON .AND. N>=RESOLVE_DEG_NTIMES) THEN

! store CHTOTL for later use in the mixer
       CHTOTL=CHTOT1
       RHOLM_LAST=RHOLM1
! determine first order change of 1._q center occupancy matrix
       CALL DEPSUM1(W0, W1, WDES, LMDIM, CRHODE1, INFO%LOVERL)

! add the derivative of occupancy matrix related to the change of the projectors
! < phi(0) | d/dR p_i> D(0)_ij - e(0) Q(0)_ij  <p_j | phi(0) >
       NIP=NI_LOCAL(ION, WDES%COMM_INB)  ! local storage index in CRHODE1
       IF (NIP/=0) CRHODE1(:,:,NIP,:)=CRHODE1(:,:,NIP,:)+CRHODE2

! change storage convention to (total, magnetization)
       CALL US_FLIP(WDES, LMDIM, CRHODE1, INFO%LOVERL, .FALSE.)

! soft pseudo charge
       CALL SOFT_CHARGE1(GRID,GRID_SOFT,W0,W1,WDES, CHDEN1)
! change storage convention to (total, magnetization)
       CALL RC_FLIP(CHDEN1, GRID_SOFT, WDES%NCDIJ, .FALSE.)

! symmetrisation of soft pseudo charge
       IF (SYMM%ISYM ==2) THEN
          IF (WDES%LNONCOLLINEAR) THEN
             CALL RHOSYM(CHDEN1(1,1),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,1)
             IF (.NOT.WDES%LSPIRAL) &
                  CALL SYMFIELD(CHDEN1(1,2),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,WDES%SAXIS,LATT_CUR)
          ELSE
             DO ISP=1,WDES%ISPIN
                CALL RHOSYM(CHDEN1(1,ISP),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,ISP)
             ENDDO
          ENDIF
       ENDIF

! derivative of augmentation charge with respect to the moved ion
! rho(0)_ij  d Q_ij(r)/dR for the single displaced ion
       CALL AUG_DER( &
            WDES, GRIDC, GRIDUS, C_TO_US, &
            LATT_CUR, P, T_INFO, SYMM, INFO%LOVERL, &
            LMDIM, CRHODE, CHTOT1, IRDMAX, ION, IDIR) 

! add first order change of 1._q center occupancy matrix times
! augmentation charge distribution rho(1)_ij Q_ij(r)
       CALL DEPLE_ADD(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
                 LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
                 LMDIM, CRHODE1, CHTOT1, CHDEN1, IRDMAX)

! set the 1._q center occupancy matrix passed through mixer (RHOLM1)
       CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
            CRHODE1, RHOLM1)

! symmetrise total pseudo charge density CHTOT if required
       IF (SYMM%ISYM ==1) THEN
          IF (WDES%LNONCOLLINEAR) THEN
             CALL RHOSYM(CHTOT1(1,1),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,1)
             IF (.NOT.WDES%LSPIRAL) &
                  &   CALL SYMFIELD(CHTOT1(1,2),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,WDES%SAXIS,LATT_CUR)
          ELSE
             DO ISP=1,WDES%ISPIN
                CALL RHOSYM(CHTOT1(1,ISP),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,ISP)
             ENDDO
          ENDIF
       ENDIF

# 747

!=======================================================================
! mix charge density
!=======================================================================
       IF (MIX%IMIX/=0) THEN
          CALL START_TIMING("G")
          INFO%LMIX=.TRUE.

          CALL SET_RHO0(GRIDC, CHTOT1, 0.0_q)

          IF (MIX%IMIX==4) THEN
!  broyden mixing ... :
             CALL BRMIX(KINEDEN,GRIDB,GRIDC,IO,MIX,B_TO_C, &
                  (2*GRIDC%MPLWV),CHTOT1,CHTOTL,WDES%NCDIJ,LATT_CUR%B, &
                  LATT_CUR%OMEGA, N_MIX_PAW, RHOLM1, RHOLM_LAST, &
                  RMST,RMSC,RMSP,WEIGHT,.TRUE.,IERRBR)
             MIX%LRESET=.FALSE.
          ELSE
!  simple mixing ... :
             RMST=0
             CALL MIX_SIMPLE(GRIDC,MIX,WDES%NCDIJ, CHTOT1,CHTOTL, &
                  N_MIX_PAW, RHOLM1, RHOLM_LAST, LATT_CUR%B, LATT_CUR%OMEGA, RMST)
          ENDIF
          CALL STOP_TIMING("G",IO%IU6,"MIXING")
          ! "mixing is ok"

          IF (IO%IU0>=0) THEN
             WRITE(IO%IU0,300) RMST*SCALE
             WRITE(17,300) RMST*SCALE
          ENDIF
300       FORMAT('   ',E10.3)
       ELSE
          IF (IO%IU0>=0) THEN
             WRITE(IO%IU0,*)
             WRITE(17,*)
          ENDIF
       ENDIF

    ELSE  chg
       IF (IO%IU0>=0) THEN
          WRITE(IO%IU0,*)
          WRITE(17,*)
       ENDIF       
    ENDIF chg
!=======================================================================
! total time used for this step
!=======================================================================
       CALL SEPERATOR_TIMING(IO%IU6)
       CALL STOP_TIMING("LOOP",IO%IU6,XMLTAG='total')

       NORDER=0 ; IF (KPOINTS%ISMEAR>=0) NORDER=KPOINTS%ISMEAR

!======================== end of loop ENDLSC ===========================
! end of the selfconsistent calculation loop
! and write commands
!=======================================================================
       INFO%LABORT=.FALSE.

       IF(ABS(DESUM)*SCALE**2<INFO%EDIFF.AND.ABS(DE)*SCALE**2<INFO%EDIFF) INFO%LABORT=.TRUE.
! rms stable for three steps, ok stop, no way to improve
       IF(N>RESOLVE_DEG_NTIMES+3) THEN
          IF (ABS((RMS(N)-RMS(N-1))/RMS(N))<1E-1 .AND. ABS((RMS(N)-RMS(N-2))/RMS(N))<1E-1) INFO%LABORT=.TRUE.
       ENDIF
!-----do not stop before minimum number of iterations is reached
       IF (N < ABS(INFO%NELMIN)) INFO%LABORT=.FALSE.

       IF (NODE_ME==IONODE) THEN
! mixing related output
      IF ( INFO%LMIX .AND. MIX%IMIX==4 ) THEN
        IF (IERRBR/=0) THEN
          IF (IO%IU0>=0) &
          WRITE(IO%IU0,*) 'ERROR: Broyden mixing failed, tried ''simple '// &
                      'mixing'' now and reset mixing at next step!'
          WRITE(IO%IU6,*) 'ERROR: Broyden mixing failed, tried ''simple '// &
                       'mixing'' now and reset mixing at next step!'
        ENDIF

        WRITE(IO%IU6,2441) RMST,RMSC,RMSP,WEIGHT
        IF (ABS(RMST-RMSC)/RMST> 0.1_q) THEN
           WRITE(IO%IU6,*) ' WARNING: grid for Broyden might be to small'
        ENDIF
        IF (MIX%NEIG > 0) THEN
           WRITE(IO%IU6,2442) MIX%AMEAN,MIX%EIGENVAL(1:MIX%NEIG)
        ENDIF
      ENDIF

 2441 FORMAT(/ &
     &       ' Broyden mixing:'/ &
     &       '  rms(total) =',E12.5,'    rms(broyden)=',E12.5,/ &
     &       '  rms(prec ) =',E12.5/ &
     &       '  weight for this iteration ',F10.2)

 2442 FORMAT(/' eigenvalues of (default mixing * dielectric matrix)' / &
             '  average eigenvalue GAMMA= ',F8.4,/ (10F8.4))

      WRITE(IO%IU6,210) E1%DENC*SCALE**2,E1%XCENC*SCALE**2,E1%PAWPS*SCALE**2,E1%PAWAE*SCALE**2, &
           E1%EENTROPY*SCALE**2,E1%EBANDSTR*SCALE**2,TOTEN*SCALE**2,(TOTEN-E1%EENTROPY)*SCALE**2

210   FORMAT(/ &
              ' Free energy of the ion-electron system (eV)'/ &
     &        '  ---------------------------------------------------'/ &
     &        '  -Hartree energ DENC   = ',F18.8/ &
     &        '  -V(xc)+E(xc)   XCENC  = ',F18.8/ &
     &        '  PAW double counting   = ',2F18.8/ &
     &        '  entropy T*S    EENTRO = ',F18.8/ &
     &        '  eigenvalues    EBANDS = ',F18.8/ &
     &        '  ---------------------------------------------------'/ &
     &        '  free energy    TOTEN  = ',F18.8,' eV'// &
     &        '  energy without entropy =',F18.8)

      IF (IO%LOPEN) CALL WFORCE(IO%IU6)
      IF (IO%LOPEN) CALL WFORCE(17)

      ENDIF
!=======================================================================
!  xml related output
!=======================================================================
       CALL XML_TAG("energy")
       IF (INFO%LABORT .OR. N==1) THEN
          CALL XML_ENERGY(TOTEN, TOTEN-E1%EENTROPY, TOTEN-E1%EENTROPY/(2+NORDER))
       ELSE
          CALL XML_ENERGY(TOTEN, TOTEN-E1%EENTROPY, TOTEN-E1%EENTROPY/(2+NORDER))
       ENDIF
       CALL XML_CLOSE_TAG
       
       CALL XML_CLOSE_TAG("scstep")


       IF (INFO%LABORT) THEN
          EXIT electron
       ENDIF
       INFO%LSOFT=.FALSE.

    ENDDO electron

    IF (IO%IU6>=0) THEN
       
       DO ISP=1,WDES%ISPIN
          IF (WDES%ISPIN==2) WRITE(IO%IU6,'(/A,I1)') ' spin component ',ISP
          DO NN=1,WDES%NKPTS
             WRITE(IO%IU6,2201)NN,WDES%VKPT(1,NN),WDES%VKPT(2,NN),WDES%VKPT(3,NN), &
                  &             (I,REAL( W1%CELTOT(I,NN,ISP) ,KIND=q) ,W1%FERTOT(I,NN,ISP)*WDES%RSPIN,I=1,WDES%NB_TOT)
          ENDDO
       ENDDO
       
2201   FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
            &         '  band No.  band energies     occupation '/ &
            &           (3X,I4,3X,F10.4,3X,F10.5))
       IF (IO%IU6>=0) WRITE(IO%IU6,131)
131    FORMAT (5X, //, &
            &       '------------------------ aborting loop because EDIFF', &
            &       ' is reached ----------------------------------------'//)
       
    ENDIF
!======================== end of loop ENDLSC ===========================
!
! force and stress
!
!=======================================================================
! add the derivate of the wavefunction character w.r.t the
! projection operators
    W1%GPROJ=W1%GPROJ+WXI%GPROJ

! no correction to force
    INFO%LCORR =.FALSE.
    F  =0
    SIF=0

!-----------------------------------------------------------------------

    POTIM=1E-3_q
    DO I=1,4
    STEP  =STENCIL4_X(I)*POTIM
    WEIGHT=STENCIL4_W(I)/POTIM
# 929


    WXI%CPTWFP    =W0%CPTWFP    +STEP* W1%CPTWFP
    WXI%GPROJ =W0%GPROJ +STEP* W1%GPROJ
    WXI%CELTOT=W0%CELTOT+STEP* W1%CELTOT
    WXI%FERTOT=W0%FERTOT+STEP* W1%FERTOT

    CHTOT2 =CHTOT +CHTOT1*STEP
    CTMP   =CRHODE+CRHODE1*STEP
    RHOLM_LAST=RHOLM  +RHOLM1*STEP

    CALL DISPLACE_ATOM(T_INFO_0%POSION, T_INFO%POSION, ION, IDIR, LATT_CUR%A, LATT_CUR%B, STEP)
    DISPL=0 ;  DISPL(IDIR,ION)=STEP

    CALL STUFAK(GRIDC,T_INFO,CSTRF1)
    IF (.NOT. INFO%LREAL) THEN
       CALL PHASE(WDES,NONL_S,0)
    ENDIF


    IF (INFO%LCORE) CALL RHOPAR(GRIDC,T_INFO,INFO,LATT_CUR,P,CSTRF1,DENCOR,IO%IU6)
    CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
         INFO,P,T_INFO,E1,LATT_CUR, &
         CHTOT2,CSTRF1,CVTOT1,DENCOR,SV1, SOFT_TO_C,XCSIF)

    CALL SETDIJ_(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO_0,INFO%LOVERL, &
         LMDIM,CDIJ1,CQIJ,CVTOT1,.TRUE.,IRDMAA,IRDMAX, DISPL)

    CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ1(1,1,1,1), RHOLM_LAST , CTMP(1,1,1,1), &
         E1,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )

    CALL FEWALD(T_INFO%POSION,EWIFOR,LATT_CUR%A,LATT_CUR%B,LATT_CUR%ANORM,LATT_CUR%BNORM, &
         LATT_CUR%OMEGA,EWSIF,E%TEWEN,T_INFO%NTYP,P%ZVALF,T_INFO%VCA, & 
         T_INFO%NIONS,T_INFO%NIOND,T_INFO%ITYP,T_INFO%NITYP,IO%IU6)
    CALL FORCE_AND_STRESS( &
         KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,WXI,LATT_CUR, &
         T_INFO,T_INFO_0,DYN,INFO,IO,MIX,SYMM,GRID,GRID_SOFT, &
         GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C, &
         CHTOT2,CHTOTL,DENCOR,CVTOT1,CSTRF1, &
         CDIJ1,CQIJ,CTMP,N_MIX_PAW,RHOLM_LAST,RHOLM_LAST, &
         CHDEN,SV1, &
         LMDIM, IRDMAX, .FALSE.,DYN%ISIF/=0,.FALSE.,.TRUE.,  &
         XCSIF, EWSIF, TSIF, EWIFOR, TIFOR, PRESS, TOTEN )
    F  =F  +TIFOR*WEIGHT
    SIF=SIF+TSIF *WEIGHT
    ENDDO
! restore all quantities to their original value

    CALL DISPLACE_ATOM(T_INFO_0%POSION, T_INFO%POSION, ION, IDIR, LATT_CUR%A, LATT_CUR%B,0._q)

! real space projector need to recalculated since force routine destroys them
    IF (INFO%LREAL) THEN
       CALL RSPHER(GRID,NONLR_S,LATT_CUR)
    ELSE
       CALL PHASE(WDES,NONL_S,0)
    ENDIF
    IF (INFO%LCORE) CALL RHOPAR(GRIDC,T_INFO,INFO,LATT_CUR,P,CSTRF,DENCOR,IO%IU6)

! write internal strain tensor and final force constants and to OUTCAR
    IF (NODE_ME==IONODE) THEN

    IF (DYN%ISIF/=0) THEN
       WRITE(IO%IU6,160) ION,IDIR,SIF
    ENDIF

160 FORMAT(/ ' INTERNAL STRAIN TENSOR  FOR ION ',I4,'  DIRECTION ',I1,3X,'(eV/Angst):',/ &
         &   ' -----------------------------------------------------------------------------',/ &
         (3X,3F9.5))

    WRITE(IO%IU6,100) ION,IDIR

    DO I=1,T_INFO%NIONS
       VTMP=T_INFO%POSION(1:3,I)
       CALL  DIRKAR(1,VTMP,LATT_CUR%A)
       WRITE(IO%IU6,110) VTMP,F(:,I)
    ENDDO
    VTMP=SUM(F,2)

    WRITE(IO%IU6,120) VTMP

100 FORMAT(// ' POSITION    ',15X,'FORCE-CONSTANT FOR ION ',I4,' DIRECTION ',I1,' (eV/Angst/Angst)'/ &
         ' ----------------------------------------------', &
         '-------------------------------------')
110 FORMAT((3F13.5,3X,3F14.6))
120 FORMAT( ' ----------------------------------------------', &
         '-------------------------------------',/ &
         '    total drift:      ',20X,3F14.6//)

    WRITE(IO%IU0,"(' force on displaced ion    ',I4,' direction',I2,' :',3F8.3)") ION,IDIR,F(:,ION)

    ENDIF
!======================== end of loop ENDLSC ===========================
!
! change of polarisation
! at this point W1%CPTWFP contains the first order change of the wavefunc.
! and W1%GPROJ the corresponding first order change of the character
!  plus the change of the character w.r.t the projection operators
! the total derivative of the polarisation is given by the total
! derivative of Equ. (24) in
! M. Gaidos, K. Hummer, G. Kresse,
! Phys. Rev. B {\bf 73}, 045112 (2006).
! the first term yields <phi_nk| d/dR_i Q(r-R_i)  |phi_nk>
!        (rigid augmentation)
! otherwise contributions from the first term drop out since the
! wavefunction remain orthogonal
! the last 1._q drops out (a dipole moving in a field gives no contribution)
! then there are terms involving the change on the wavefunction on the
! plane wave  grid   < -i d/dk phi_nk |  d/dR_i phi_nk>
! and terms involving the total change of the wavefunction character
!     \sum_j beta_j* (<p_j | d/dR_i phi_nk> + <d/dR_i p_j | phi_nk >)
! beta_j is defined at the top of LRF_RPHI
!
!=======================================================================
    IF (LEPSILON) THEN
       DO I=1,3
          CALL KDER_WAVE (WXI, GRID, WDES, I, RPHI, RPHI_CPROJ,  &
               T_INFO, P, LATT_CUR, KPOINTS_TRANS)

! rotate the polarization vector in the same manner as W1
          CALL SUBROT_DEG_ALL(WDES, WXI, .FALSE., .TRUE., DEG_CLUSTER )
          WXI%FERTOT=W0%FERTOT
! calculate
! < d phi_nk/d_k |  d phi_nk/d R_i>   and
! \sum_j beta_j* (<p_j | d phi_nk/d R_i> + <d p_j/d R_i | phi_nk >)
          CALL SUM_INPROD(WDES, WXI, W1, Z(I), .TRUE., LMDIM)
       ENDDO
       VTMP=0
       VTMP(IDIR)=VTMP(IDIR)+P(T_INFO%ITYP(ION))%ZVALF*T_INFO%VCA(T_INFO%ITYP(ION))+ONE_CENTRE_CHARGE(WDES, T_INFO, P, ION, LMDIM, CQIJ, CRHODE)

       IF (SYMM%ISYM>0) CALL POLSYM(Z,LATT_CUR%A)
       IF (SYMM%ISYM>0) CALL POLSYM(VTMP,LATT_CUR%A)

       BORN_CHARGE=Z+VTMP
       IF (IO%IU0>=0) THEN
          WRITE(IO%IU0,"(' Born effective charge ion ',I4,' direction',I2,' :',3F8.3)") ION,IDIR,Z+VTMP
          WRITE(IO%IU6,"(' BORN EFFECTIVE CHARGE FOR ION ',I4,' DIRECTION',I2,' :',3F8.3)") ION,IDIR,Z+VTMP

          VTMP(IDIR)=VTMP(IDIR)-P(T_INFO%ITYP(ION))%ZVALF*T_INFO%VCA(T_INFO%ITYP(ION))
          WRITE(IO%IU6,"('                   rigid augmentation           :',3F8.3)") VTMP

          VTMP(IDIR)=P(T_INFO%ITYP(ION))%ZVALF*T_INFO%VCA(T_INFO%ITYP(ION))
          WRITE(IO%IU6,"('                   ionic contribution           :',3F8.3)") VTMP
          WRITE(IO%IU6,"('                   Berry contribution           :',3F8.3)") Z
       ENDIF
    END IF

    IF (IO%IU6>=0) WRITE(IO%IU6,130)

    CALL SUBROT_DEG_ALL(WDES, W0, .TRUE., .TRUE., DEG_CLUSTER )

    CALL DEALLOCW(WXI)
    CALL DEALLOCW(W1)
    DEALLOCATE(CHTOT1,CHTOT2,CVTOT1,CSTRF1,CSTRF2,CWORK,CDIJ1,CRHODE1,CRHODE2,CTMP,SV1,SV2,CHDEN1)

    DEALLOCATE(CSTRF3,CSTRF4)

  END SUBROUTINE LR_MAIN

END MODULE lri_main
