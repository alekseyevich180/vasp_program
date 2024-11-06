# 1 "fast_aug.F"
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

# 2 "fast_aug.F" 2 
!
! the following define statements are used to determine
! which routines are used do transform
!     lmn,l'm'n' to LMN and vice versa
! (Clebsch Gordon like transformations)
!
! if kernel_dgemm is define DGEMM is used to tranform
! if kernel_dgemv is define DGEMV is used to tranform
! if none is define DO loops are used
!
! kernel_dgemv should be always set, since this improves
! performance by a factor 4
! kernel_dgemm is only implemented in few versions
! and can be additionally set
! using mkl this is however often slower then the DGEMV kernels
!

!#define kernel_dgemm


!***********************************************************************
!
! this module implements all support routines required for a fast
! augmentation of pseudo chargedensities on a plane wave grid
!
! it is based massively on nonr.F
! the augmentation is 1._q in real space
!
! written by gK
!
!***********************************************************************

!***********************************************************************
!
! the first set of routines allows for the fast evaluation of
! augmentation charges on the course FFT grid (GRID)
! the routines rely on nonlr.F and mimics the functionality of us.F
!
!***********************************************************************


MODULE augfast
  USE prec
  USE poscar
  USE nonlr
  IMPLICIT none

  INTEGER, PARAMETER, PRIVATE :: LMAX_AUG_ACCURATE=4

! this handle can be used to quickly extract the 1._q center
! occupancies up to a certain main quantum number LMAX
! these 1._q center terms are stored in consecutive order
! and can be used to calculate e.g. 1._q center corrections
! in MP2
  TYPE one_center_handle
     INTEGER :: LMAX_ONE_CENTER             ! maximum L for 1._q center terms
     INTEGER :: TOTAL_ENTRIES               ! total number of entries in density related arrays
     INTEGER :: TOTAL_ENTRIES_CHI           ! total number of entries in response function arrays
! note that this is SUM (H%ENTRIES_FOR_TYPE(:)*WDES%NITYP(:))
     INTEGER, POINTER :: ENTRIES_FOR_TYPE(:)! number of Clebsch Gordan transformed entries for each type
     INTEGER, POINTER :: CHANNEL1(:,:)
     INTEGER, POINTER :: CHANNEL2(:,:)
     INTEGER, POINTER :: LM(:,:)
     INTEGER, POINTER :: L(:,:)
! this array stores for each combined index consisting of i=(channel1,channel2,L,M), type
! the corresponding CHANNEL number and the L quantum number
     INTEGER, POINTER :: NR_OF_CG(:,:)      !
! this array stores for each combined index i=(channel1,channel2,L,M) and type
! how many elements contribute after transformation using CG coefficients
! the typical loop constructs is
! one_center_charge = sum_(i=1, NR_OF_CG(:,:)-1)
!                          CG(i)*CRHODE(INDEX(i))
     REAL(q), POINTER :: CG(:)              ! array for Clebsch Gordan coefficients
     INTEGER, POINTER :: INDEX(:)           ! position at which contribution is stored
     REAL(q), POINTER :: POSION(:,:)        ! position of ions (link to T_INFO%POSION)
! this index is 1._q dimensional and given by first_index+ WDES(NT)%LMMAX*second_index
! where first_index and second_index are the usual indices into CRHODE or CDIJ
! like arrays
     REAL(q), POINTER :: POT(:)             ! kernel for four orbital integrals
! for each type a matrix with ENTRIES_FOR_TYPE x ENTRIES_FOR_TYPE entries are stored
  END TYPE one_center_handle
CONTAINS

!************************ SUBROUTINE SETUP_FASTAUG ********************
!
! the following subroutine sets up the arrays required to perform a fast
! augmentation using the non-local pseudopotential routines (in real
!  space)
! the  final FAST_AUG structure is a non local pseudopotential
! descriptor it can be handled by all routines in nonlr.F
! the corresponding deallocation routine is also available
!
!***********************************************************************


  SUBROUTINE SETUP_FASTAUG( T_INFO, P, WDES, GRID, LATT_CUR,  LMDIM, SZPRECFOCK, &
      FAST_AUG, TRANS_MATRIX, LMAX_AUG, LMAX_FOCKAE, NMAX_FOCKAE, QMAX_FOCKAE, EXXOEP, IU6, IU0 )
    USE poscar
    USE nonlr
    USE pseudo
    USE constant
    USE asa
    USE wave
    USE mgrid
    USE lattice
    USE ini
    IMPLICIT NONE

    TYPE (type_info) T_INFO
    TYPE (potcar),TARGET::  P(T_INFO%NTYP)
    TYPE (wavedes) :: WDES       ! wave function descriptor
    TYPE (grid_3d) :: GRID       ! 3d- grid descriptor
    TYPE (latt)    :: LATT_CUR   ! lattice parameters
    CHARACTER (6)  :: SZPRECFOCK

    INTEGER        :: LMAX_AUG   ! maximum l augmentation
    INTEGER        :: LMAX_FOCKAE! maximum l for precise restoration of AE charge density
    INTEGER        :: NMAX_FOCKAE! maximum number of functions
    REAL(q)        :: QMAX_FOCKAE(:)

    TYPE (nonlr_struct) :: FAST_AUG
    REAL(q), POINTER :: TRANS_MATRIX(:,:,:,:) 

    INTEGER        :: LMDIM      ! leading dimension of arrays like CDIJ
    INTEGER :: EXXOEP
    INTEGER :: IU6, IU0          ! io units
! local
    INTEGER :: NIS, NIP, N, NI, NT, NTP, LOCAL_IONS, NPRO, LMMAXC, L, CH1, CH2, MAXL_TMP, LMAX_FOCKAE_TYPE
    INTEGER :: LASTTYP, NALLOC, CHANNELS_WITHOUT_ACCURATE
    INTEGER :: LMMAX_AUG, LMAX_FOCKAE_TYPE_STORE(T_INFO%NTYP)
    LOGICAL :: LREALLOCATE
    REAL(q) :: RMAXL(0:10)       ! outermost point for augmentation charge for each L quantum number
    REAL(q) :: RMAX              ! outermost point for augmentation charges (max (RMAXL))
    REAL(q) :: QMAX              ! maximum allowed q values in augmentation charge
    LOGICAL, EXTERNAL :: USE_OEP_IN_GW
    TYPE (communic) COMM

    COMM = WDES%COMM_INB

! set number of ions and types to 0._q (indicates FAST_AUG is not set up)
    FAST_AUG%NTYP =0
    FAST_AUG%NIONS=0

    NULLIFY(TRANS_MATRIX)
    IF (LMAX_AUG < 0) RETURN
!
! set up the arrays required to perform the fast
! augmentation on the grid GRID
! the applied routines are equivalent to the nonl pseudopotentials
! in real space
!
    IF ( WDES%LSPIRAL) THEN
       IF (IU0>0 ) WRITE(IU0, *) 'SETUP_FASTAUG: spin spirals are presently not supported'
       CALL M_exit(); stop
    ENDIF
    CALL  NONLR_SETUP(FAST_AUG , T_INFO, P, .TRUE., WDES%LSPIRAL)

    IF (IU6>=0) THEN
       WRITE(IU6,*)
       WRITE(IU6,'(A)') ' Radii for the augmentation spheres in the non-local exchange'
    ENDIF


! initialize the various modified augmentation charges
! such as the soft augmentation charge used on the HF grid
! the augmentation charges that allow to restore the all-electron density
! on the plane wave grid etc.
    
    DO NT=1,T_INFO%NTYP
! maximum L quantum number for this species
       MAXL_TMP=0
       DO L=1,P(NT)%LMAX 
          MAXL_TMP=MAX(MAXL_TMP,P(NT)%LPS(L))
       ENDDO

       IF (.NOT.ASSOCIATED(P(NT)%QPAW)) THEN
          FAST_AUG%LMAX(NT)    = -1
          FAST_AUG%LMMAX(NT)   =  0
          FAST_AUG%CHANNELS(NT)=  0
          FAST_AUG%PSRMAX(NT)  = P(NT)%PSDMAX
          NULLIFY(P(NT)%QPAW_FOCK)
          LMAX_FOCKAE_TYPE_STORE(NT)=-1
       ELSE

          FAST_AUG%LMAX(NT)    = MIN(LMAX_AUG,2*MAXL_TMP)
          FAST_AUG%LMMAX(NT)   =(FAST_AUG%LMAX(NT)+1)*(FAST_AUG%LMAX(NT)+1)
          FAST_AUG%CHANNELS(NT)= FAST_AUG%LMAX(NT)+1

          IF (LMAX_FOCKAE>=0 .OR. EXXOEP/=0 .OR. USE_OEP_IN_GW()) THEN

! if LMAX_FOCKAE is set, do exactly the same as in the past
! in particular perform no optimization of augmentation charges
! (causing incompatibilities between density (DFT) part and Fock part
! this is required for all OEP treatments

             RMAX=P(NT)%R%RMAX
             RMAXL=RMAX
             QMAX=0

          ELSE IF (SZPRECFOCK(1:1)=='a' .OR. SZPRECFOCK(1:1)=='n') THEN
!
! soften the augmentation charges for exact exchange part
! by increasing the radius for the augmentation
! small softening
!
             RMAX=MIN(P(NT)%R%RMAX*1.25,P(NT)%R%REND)
             RMAXL=RMAX
             QMAX=2.0_q*SQRT(WDES%ENMAX/RYTOEV)/AUTOA

          ELSE IF (SZPRECFOCK(1:1)=='f' ) THEN

! soften p and d somewhat more
! this makes the energies more accurate,
! however restarting from a "normal" WAVECAR requires more iterations
!
             RMAX=MIN(P(NT)%R%RMAX*1.35,P(NT)%R%REND)
             RMAXL=RMAX
             RMAXL(0)=MIN(P(NT)%R%RMAX*1.25,P(NT)%R%REND)
             QMAX=2.0_q*SQRT(WDES%ENMAX/RYTOEV)/AUTOA

          ELSE IF (SZPRECFOCK(1:1)=='s' ) THEN

! soften the channels even more, this flag is not really recommended
! (testing only)
             RMAX=MIN(P(NT)%R%RMAX*1.45,P(NT)%R%REND)
             RMAXL=RMAX
             RMAXL(0)=MIN(P(NT)%R%RMAX*1.35,P(NT)%R%REND)
             QMAX=2.0_q*SQRT(WDES%ENMAX/RYTOEV)/AUTOA

          ELSE
! fallback: PRECFOCK = Medium
! use the same augmentation as for density functional part

             RMAX=P(NT)%R%RMAX
             RMAXL=RMAX
             QMAX=0

          ENDIF


          IF (IU6>=0) THEN
             WRITE(IU6,'(A,I3,A,F8.3,A,F8.3,A)') ' for species ',NT,' augmentation radius',RMAX,' (default was',P(NT)%R%RMAX,')'
             IF (QMAX/=0) &
             WRITE(IU6,'(A,2F8.1)') '       energy cutoff for augmentation ',(QMAX*AUTOA)**2*RYTOEV
          ENDIF

! PSRMAX determines the step size
! the outermost point is FAST_AUG%PSRMAX(NT)/NPSRNL*(NPSRNL-1) in nonlr.F
! this is brain damaging I know (gK)
          FAST_AUG%PSRMAX(NT)  =RMAX*NPSRNL/(NPSRNL-1)

          CHANNELS_WITHOUT_ACCURATE=FAST_AUG%CHANNELS(NT)

          IF (LMAX_FOCKAE>=0 ) THEN
! increase LMMAX and CHANNELS as required
             LMAX_FOCKAE_TYPE=MIN(FAST_AUG%LMAX(NT), LMAX_FOCKAE)
             
             FAST_AUG%LMMAX(NT)   =FAST_AUG%LMMAX(NT) + &
                  (LMAX_FOCKAE_TYPE+1)*(LMAX_FOCKAE_TYPE+1)*NMAX_FOCKAE
             FAST_AUG%CHANNELS(NT)=FAST_AUG%CHANNELS(NT)+ (LMAX_FOCKAE_TYPE+1)*NMAX_FOCKAE
             
! setup QDEP_FOCK and AUG_FOCK (i.e. the augmentation charges on
! a regular grid and the spherical grids)
! QDEP_FOCK contains the augmentation charges on a equally spaced grid
! this array is used by nonlr.F to quickly add the augmentation charges
! the array first contains the regular augmentations charges 0:LMAX_FOCKAE_TYPE
! and then 0-normalized augmentation charges to more accurately restore the
! shape of the charge density
! it is certain that QDEP and QDEP_FOCK are identical from 0_LMAX_FOCKAE_TYPE
             CALL SET_BESSEL_AUG(P(NT), RMAX, RMAXL, QMAX, CHANNELS_WITHOUT_ACCURATE-1, LMAX_FOCKAE_TYPE, NMAX_FOCKAE )
             
             ALLOCATE(P(NT)%QPAW_FOCK( SIZE(P(NT)%QPAW,1) , SIZE(P(NT)%QPAW,2), & 
                  0:LMAX_FOCKAE_TYPE, NMAX_FOCKAE))

! determine AUG_FOCK
             CALL RAD_QPAW_FOCK( P(NT)%R, P(NT)%LMAX, RMAXL, QMAX_FOCKAE(NT), &
                  P(NT)%WAE, P(NT)%WPS, P(NT)%AUG, P(NT)%AUG_FOCK, P(NT)%LPS, &
                  P(NT)%QPAW,P(NT)%QPAW_FOCK, LMAX_FOCKAE_TYPE, NMAX_FOCKAE , IU6)
             
             FAST_AUG%BETA(NT)%PSPRNL=>P(NT)%QDEP_FOCK(:,:,:)
          ELSE
             LMAX_FOCKAE_TYPE=LMAX_FOCKAE

! setup QDEP_FOCK and AUG_SOFT
             CALL SET_BESSEL_AUG(P(NT), RMAX, RMAXL, QMAX, FAST_AUG%LMAX(NT), -1, 0 )
             FAST_AUG%BETA(NT)%PSPRNL=>P(NT)%QDEP_FOCK(:,:,0:FAST_AUG%LMAX(NT))
          ENDIF
          
          LMAX_FOCKAE_TYPE_STORE(NT)=LMAX_FOCKAE_TYPE

          NULLIFY(FAST_AUG%BETA(NT)%LPS)
          ALLOCATE(FAST_AUG%BETA(NT)%LPS(FAST_AUG%CHANNELS(NT)))

          DO L=1,CHANNELS_WITHOUT_ACCURATE
             FAST_AUG%BETA(NT)%LPS(L)=L-1
          ENDDO
          
          DO N=1,NMAX_FOCKAE
             DO L=0,LMAX_FOCKAE_TYPE
                FAST_AUG%BETA(NT)%LPS(CHANNELS_WITHOUT_ACCURATE+L+1+(N-1)*(LMAX_FOCKAE_TYPE+1))=L
             ENDDO
          ENDDO
       END IF
    ENDDO

! determine the maximum possible number of augmentation channels
    LMMAX_AUG=0
    NIS   =1
    DO NT=1,WDES%NTYP
! determine global storage index for this ion
       NI=NI_GLOBAL(NIS, WDES%COMM_INB)
       NTP=T_INFO%ITYP(NI)
       IF (LMAX_FOCKAE>=0 ) THEN
          LMAX_FOCKAE_TYPE=MIN(FAST_AUG%LMAX(NTP), LMAX_FOCKAE)
       ELSE
          LMAX_FOCKAE_TYPE=-1
       ENDIF
! TODO: LMAX_FOCKAE_TYPE_STORE(NTP)   can be removed
!  after this has been carefully evaluated
       IF (LMAX_FOCKAE_TYPE_STORE(NTP) /= LMAX_FOCKAE_TYPE) THEN
          WRITE(*,*) 'internal error in SETUP_FASTAUG: LMAX_FOCKAE_TYPE incorrect',NT,NTP,LMAX_FOCKAE_TYPE_STORE(NTP),LMAX_FOCKAE_TYPE
          CALL M_exit(); stop
       ENDIF

       LMMAX_AUG=MAX(LMMAX_AUG, (FAST_AUG%LMAX(NTP)+1)**2 + & 
                  MIN(FAST_AUG%LMAX(NTP)+1,LMAX_FOCKAE+1)**2*NMAX_FOCKAE)
       NIS = NIS+WDES%NITYP(NT)
    ENDDO
! old statement and less tight bounds on allocation
!LMMAX_AUG=(LMAX_AUG+1)*(LMAX_AUG+1) + (LMAX_FOCKAE+1)*(LMAX_FOCKAE+1)*NMAX_FOCKAE

!
! loop over all types that are local on this node
! set up the transformation matrix
!
    ALLOCATE(TRANS_MATRIX( LMDIM, LMDIM, LMMAX_AUG, WDES%NTYP))
    CALL REGISTER_ALLOCATE(8._q*SIZE(TRANS_MATRIX), "HF")

    NIS   =1
    typ:  DO NT=1,WDES%NTYP
! determine global storage index for this ion
       NI=NI_GLOBAL(NIS, WDES%COMM_INB)
       NTP=T_INFO%ITYP(NI)
       IF (NTP /= WDES%NT_GLOBAL(NT)) THEN
          WRITE(*,*) 'internal error in SETUP_FASTAUG: WDES%NT_GLOBAL inconsistent with NTP',WDES%NT_GLOBAL(NT), NTP
       ENDIF
       IF (LMAX_FOCKAE>=0 ) THEN
          LMAX_FOCKAE_TYPE=MIN(FAST_AUG%LMAX(NTP), LMAX_FOCKAE)
       ELSE
          LMAX_FOCKAE_TYPE=-1
       ENDIF

       CALL SETUP_TRANS_MATRIX( P(NTP), TRANS_MATRIX(:, :, :, NT), FAST_AUG%LMAX(NTP), LMAX_FOCKAE_TYPE , NMAX_FOCKAE)
       IF (T_INFO%VCA(NTP)/=1.0) THEN
          TRANS_MATRIX(:, :, :, NT)=TRANS_MATRIX(:, :, :, NT)*T_INFO%VCA(NTP)
       ENDIF
       NIS = NIS+WDES%NITYP(NT)
    ENDDO typ

! smooth augmentation charges
!   FAST_AUG%RSMOOTH=LATT_CUR%ANORM(1)/GRID%NGX*(1.25-1E-4)

    CALL REAL_OPTLAY(GRID, LATT_CUR, FAST_AUG, .TRUE., LREALLOCATE, IU6, IU0)
    CALL NONLR_ALLOC(FAST_AUG)
!
! probably something most be 1._q here for spin spirals ...
!
    IF (ASSOCIATED(FAST_AUG%CRREXP)) THEN
       FAST_AUG%CRREXP=1
    ELSE
!       WRITE(0,*)'report SETUP_FASTAUG: FAST_AUG%CRREXP is not associated, bypassing set to 1'
    ENDIF

  END SUBROUTINE SETUP_FASTAUG


!************************ SUBROUTINE FASTAUG_SUPER      ****************
!
! FASTAUG_SUPER copies the FASTAUG structure and reinitialize
! those elements that might have changed when a super cell is used
!
!***********************************************************************


  SUBROUTINE FASTAUG_SUPER( T_INFO, P, WDES, GRID, LATT_CUR,  LMDIM, & 
       FAST_AUG, FAST_AUG_NEW, TRANS_MATRIX, LMAX_FOCKAE, NMAX_FOCKAE )
    USE poscar
    USE nonlr
    USE pseudo
    USE constant
    USE asa
    USE wave
    USE mgrid
    USE lattice
    USE ini
    IMPLICIT NONE

    TYPE (type_info) T_INFO
    TYPE (potcar),TARGET::  P(T_INFO%NTYP)
    TYPE (wavedes) :: WDES       ! wave function descriptor
    TYPE (grid_3d) :: GRID       ! 3d- grid descriptor
    TYPE (latt)    :: LATT_CUR   ! lattice parameters
    INTEGER        :: LMDIM      ! leading dimension of arrays like CDIJ
    TYPE (nonlr_struct) :: FAST_AUG
    TYPE (nonlr_struct) :: FAST_AUG_NEW
    REAL(q), POINTER :: TRANS_MATRIX(:,:,:,:) 
    INTEGER        :: LMAX_FOCKAE, NMAX_FOCKAE

! local
    TYPE (communic) COMM
    INTEGER :: LMMAX_AUG, NI, NIS, NT, NTP, LMAX_FOCKAE_TYPE
    LOGICAL :: LREALLOCATE

    COMM = WDES%COMM_INB

    NULLIFY(TRANS_MATRIX)

! start with a simple copy of FAST_AUG
    FAST_AUG_NEW=FAST_AUG

! now relink to proper ion number
    IF (T_INFO%NTYP /=FAST_AUG_NEW%NTYP ) THEN
       WRITE(0,*) 'internal error in FASTAUG_SUPER: number of types must not change'
    ENDIF

    FAST_AUG_NEW%NTYP =T_INFO%NTYP
    FAST_AUG_NEW%NIONS=T_INFO%NIONS
    FAST_AUG_NEW%NITYP  =>T_INFO%NITYP
    FAST_AUG_NEW%ITYP   =>T_INFO%ITYP
    FAST_AUG_NEW%POSION =>T_INFO%POSION

! loop over all types and determine maximum L quantum number
    NIS   =1
    LMMAX_AUG=0
    DO NT=1,WDES%NTYP
! determine global storage index for this ion
       NI=NI_GLOBAL(NIS, WDES%COMM_INB)
       NTP=T_INFO%ITYP(NI)
       IF (LMAX_FOCKAE>=0 ) THEN
          LMAX_FOCKAE_TYPE=MIN(FAST_AUG_NEW%LMAX(NTP), LMAX_FOCKAE)
       ELSE
          LMAX_FOCKAE_TYPE=-1
       ENDIF

       LMMAX_AUG=MAX(LMMAX_AUG,(FAST_AUG_NEW%LMAX(NTP)+1)**2 + & 
                  MIN(FAST_AUG_NEW%LMAX(NTP)+1,LMAX_FOCKAE+1)**2*NMAX_FOCKAE)
       NIS = NIS+WDES%NITYP(NT)
    ENDDO
!
! loop over all types that are local on this node
! setup the transformation matrix
!
    ALLOCATE(TRANS_MATRIX( LMDIM, LMDIM, LMMAX_AUG, WDES%NTYP))
    CALL REGISTER_ALLOCATE(8._q*SIZE(TRANS_MATRIX), "HF")

    NIS   =1

    typ:  DO NT=1,WDES%NTYP
! determine global storage index for this ion
       NI=NI_GLOBAL(NIS, WDES%COMM_INB)
       NTP=T_INFO%ITYP(NI)
       IF (NTP /= WDES%NT_GLOBAL(NT)) THEN
          WRITE(*,*) 'internal error in FASTAUG_SUPER: WDES%NT_GLOBAL inconsistent with NTP',WDES%NT_GLOBAL(NT), NIS, NI, NTP
          CALL M_exit(); stop
       ENDIF
       IF (LMAX_FOCKAE>=0 ) THEN
          LMAX_FOCKAE_TYPE=MIN(FAST_AUG_NEW%LMAX(NTP), LMAX_FOCKAE)
       ELSE
          LMAX_FOCKAE_TYPE=LMAX_FOCKAE
       ENDIF

       CALL SETUP_TRANS_MATRIX( P(NTP), TRANS_MATRIX(:, :, :, NT), FAST_AUG_NEW%LMAX(NTP), LMAX_FOCKAE_TYPE , NMAX_FOCKAE)
       IF (T_INFO%VCA(NTP)/=1.0) THEN
          TRANS_MATRIX(:, :, :, NT)=TRANS_MATRIX(:, :, :, NT)*T_INFO%VCA(NTP)
       ENDIF
       NIS = NIS+WDES%NITYP(NT)
    ENDDO typ

!  not sure whether reallocation in real space can be 1._q better not
!    CALL REAL_OPTLAY(GRID, LATT_CUR, FAST_AUG_NEW, .TRUE., LREALLOCATE, -1, -1)
    CALL REAL_OPTLAY(GRID, LATT_CUR, FAST_AUG_NEW, .FALSE., LREALLOCATE, -1, -1)
    CALL NONLR_ALLOC(FAST_AUG_NEW)
!
! probably something must be 1._q here for spin spirals ...
!
    IF (ASSOCIATED(FAST_AUG_NEW%CRREXP)) THEN
       FAST_AUG_NEW%CRREXP=1
    ELSE
!       WRITE(0,*)'report FASTAUG_SUPER: FAST_AUG_NEW%CRREXP is not associated, bypassing set to 1'
    ENDIF

  END SUBROUTINE FASTAUG_SUPER


  SUBROUTINE DEALLOC_FASTAUG( FAST_AUG, TRANS_MATRIX )
    USE nonlr
    IMPLICIT NONE

    TYPE (nonlr_struct) :: FAST_AUG
    REAL(q), POINTER :: TRANS_MATRIX(:,:,:,:) 
! local
    INTEGER NT

    IF (.NOT. ASSOCIATED(TRANS_MATRIX)) RETURN

    DEALLOCATE(TRANS_MATRIX)
    DO NT=1,FAST_AUG%NTYP
       DEALLOCATE(FAST_AUG%BETA(NT)%LPS)
    ENDDO
    CALL  NONLR_DEALLOC(FAST_AUG)

  END SUBROUTINE DEALLOC_FASTAUG


  SUBROUTINE COPY_FASTAUG( FAST_AUG1, FAST_AUG2)
    USE nonlr
    IMPLICIT NONE
    TYPE (nonlr_struct) :: FAST_AUG1, FAST_AUG2

    FAST_AUG2=FAST_AUG1
    NULLIFY(FAST_AUG2%NLIMAX,FAST_AUG2%NLI,FAST_AUG2%RPROJ,FAST_AUG2%CRREXP)
    CALL NONLR_ALLOC(FAST_AUG2)
    
  END SUBROUTINE COPY_FASTAUG

!************************ SUBROUTINE SETUP_AUG_DES  ********************
!
! in adddition to a non local pseudopotential descriptor the non local
! pseudopotential routines require also a descriptor for the layout
! of the site occupancy matrix (usually the wavefunction character
!  (wavedes1)
! this descriptor AUG_DES describes how the moments are stored and
! arranged
! this descriptor must be compatible to FAST_AUG and is for most
! of the entries only a copy of FAST_AUG
!
!***********************************************************************


  SUBROUTINE SETUP_AUG_DES(GRID, WDES, AUG_DES, FAST_AUG )
    USE prec
    USE wave
    USE mgrid
    USE nonlr
    IMPLICIT NONE

    TYPE (wavedes)      WDES
    TYPE (grid_3d)      GRID
    TYPE (wavedes1)     AUG_DES
    TYPE (nonlr_struct) FAST_AUG
! local
    TYPE (communic) COMM
    INTEGER :: LMMAX_AUG
    INTEGER NALLOC, LASTTYP, NIS, NPRO, LMMAXC, IONS_LOCAL, NT, NT_LOCAL, NI

    COMM = WDES%COMM_INB

! first copy all constituents from WDES
! includes NTYP, NIONS, NITYP, NT_GLOBAL (all correct!)
    CALL SETWDES(WDES,AUG_DES,1)
! now copy the grid information from GRID
    CALL SETWGRID_OLD(AUG_DES,GRID)
! finally correct those entries that are different
! essentially only the L quantum number related parts must be
! updated
    AUG_DES%LNONCOLLINEAR=.FALSE.   ! loop over spinors is 1._q externally
    AUG_DES%RINPL        =1         ! multiplicator used by RACC0 and RACCMU nonlr.F
    AUG_DES%NRSPINORS    =1         ! no spinors, since we have only charge

! the maximum pair of NLM quantum numbers is different and
! needs to be reset
    ALLOCATE(AUG_DES%LMMAX(SIZE(WDES%LMMAX)))
    IF (SIZE(AUG_DES%LMMAX) < AUG_DES%NTYP) THEN
       WRITE(*,*)'internal ERROR 1 in SETUP_AUG_DES: size insufficient',SIZE(AUG_DES%LMMAX), AUG_DES%NTYP
       CALL M_exit(); stop
    ENDIF
! seriel version simply copy to FAST_AUG%LMMAX
    AUG_DES%LMMAX(1:AUG_DES%NTYP)=FAST_AUG%LMMAX(1:AUG_DES%NTYP)
!
! parallel version needs some more work
! determine how occupancies are distributed (this is very similar to
! after the occupancies have been merged from all node (MRG_PROJ)
! (only required if NCORE > 1)

    ALLOCATE(AUG_DES%NPRO_POS(SIZE(WDES%NPRO_POS)))
    AUG_DES%NPRO_POS=0

    LASTTYP=0
    NIS    =1
    NPRO   =0

    IONS_LOCAL=0
    NT_LOCAL=0
    LASTTYP=0
    DO NT=1,FAST_AUG%NTYP
       LMMAXC=FAST_AUG%LMMAX(NT)
       DO NI=NIS,FAST_AUG%NITYP(NT)+NIS-1
! does this element reside on local node
          IF (NI_LOCAL(NI,COMM) /=0 ) THEN
             IF (NT /= LASTTYP) THEN
                NT_LOCAL=NT_LOCAL+1
                LASTTYP=NT
             ENDIF
             AUG_DES%LMMAX(NT_LOCAL)=LMMAXC
             IONS_LOCAL=IONS_LOCAL+1
             AUG_DES%NPRO_POS(IONS_LOCAL)=NPRO
          ENDIF
          NPRO= LMMAXC+NPRO
       ENDDO
100    NIS = NIS+FAST_AUG%NITYP(NT)
    ENDDO

    IF (IONS_LOCAL /= AUG_DES%NIONS) THEN
       WRITE(*,*)'internal ERROR 2 in SETUP_AUG_DES: local ions do not match ',IONS_LOCAL,AUG_DES%NIONS
       CALL M_exit(); stop
    ENDIF
    IF (NT_LOCAL /= AUG_DES%NTYP) THEN
       WRITE(*,*)'internal ERROR 3 in SETUP_AUG_DES: local types do not match ',NT_LOCAL,AUG_DES%NTYP
       CALL M_exit(); stop
    ENDIF


    AUG_DES%NPRO         =SUM(AUG_DES%LMMAX *AUG_DES%NITYP)
    AUG_DES%NPROD        =((AUG_DES%NPRO+AUG_DES%NB_PAR-1)/AUG_DES%NB_PAR)*AUG_DES%NB_PAR 
! total number of entries equals FAST_AUG
    AUG_DES%NPRO_TOT     =SUM(FAST_AUG%NITYP*FAST_AUG%LMMAX)

    NPRO=AUG_DES%NPRO
    CALL M_sum_i(COMM,NPRO ,1)
    IF (NPRO/= AUG_DES%NPRO_TOT) THEN
       WRITE(*,*)'internal ERROR 1 in SETUP_AUG_DES: NPRO does not match ',NPRO,AUG_DES%NPRO_TOT

       WRITE(*,*) 'NTYP',COMM%NODE_ME,AUG_DES%NTYP,AUG_DES%NPRO
       WRITE(*,*) 'NITYP',COMM%NODE_ME,AUG_DES%NITYP
       WRITE(*,*) 'LMMAX',COMM%NODE_ME,AUG_DES%LMMAX

       CALL M_exit(); stop
    ENDIF

  END SUBROUTINE SETUP_AUG_DES

  SUBROUTINE DEALLOC_AUG_DES(AUG_DES)
    USE wave
    TYPE (wavedes1)     AUG_DES

    DEALLOCATE(AUG_DES%LMMAX)

    DEALLOCATE(AUG_DES%NPRO_POS)

  END SUBROUTINE DEALLOC_AUG_DES


!*******************************************************************
!
!  this subroutine sets up an array of orthogonal functions that
!  can be used to restore the AE charge density on the plane wave
!  grid accurately
!  VASP usually restores all moments (monopole, dipole, quadrupole etc.
!  on the plane wave grid).
!  Sometimes this is not sufficient (GW, RPA correlation)
!  and 1._q would like to restore the AE charge density accurately
!  on the plane wave grid.
!  This is 1._q by developing the difference between the PS
!  and AE charge density into a set of mutual orthogonal functions
!  with 0._q moment l, where l is the angular quantum number
!
!*******************************************************************
  
  SUBROUTINE SET_BESSEL_AUG( P, RMAX, RMAXL, QMAX, LMAX, LMAX_BESSEL, NMAX_BESSEL )
    USE radial
    USE pseudo
    IMPLICIT NONE
    TYPE (potcar),TARGET :: P
    INTEGER :: IU6
    INTEGER :: LMAX_BESSEL, NMAX_BESSEL
! local variables
    TYPE (rgrid),POINTER :: R
    INTEGER  CHANNELS,LMAX,L,LP,N,I,NBESSEL,INDEX_BESSEL
    INTEGER  NQ
    REAL(q)  RMAX
    REAL(q)  RMAXL(0:)
    REAL(q)  QMAX
    INTEGER, PARAMETER :: NDIM=20
    REAL(q)  QQ(NDIM), A(NDIM,NMAX_BESSEL), AA(NDIM)
    REAL(q)  SUM,QR,BJ,STEP,X
    INTEGER :: LMAX_SOFT

    IF (NDIM<2+NMAX_BESSEL) THEN
       WRITE(0,*) 'internal error in SET_BESSEL_AUG: increase NDIM'
    ENDIF

! allocate QDEP_FOCK
    R=>P%R

    ALLOCATE(P%QDEP_FOCK(NPSRNL,5, 0:LMAX + (LMAX_BESSEL+1)*NMAX_BESSEL))

! AUG_SOFT must contain the same number of data as PP%AUG
! otherwise pawfock.F is in deep trouble
    LMAX_SOFT=0
    DO I=1,P%LMAX
       LMAX_SOFT=MAX( P%LPS(I),LMAX_SOFT )
    ENDDO

    LMAX_SOFT=MAX(LMAX_SOFT*2+1,LMAX) ! maximum l in augmentation charges as used in paw.F
    ALLOCATE(P%AUG_SOFT(R%NMAX,0:LMAX_SOFT) )

    ALLOCATE(P%AUG_FOCK (R%NMAX, 0:LMAX_BESSEL, NMAX_BESSEL))

    STEP= RMAX/(NPSRNL-1)

! now set up new entries

!
! lower panel 0:LMAX (these are the standard 1._q normalized augmentation charges)
! usually the P%AUG is equal to P%AUG_SOFT and P%QDEP_FOCK to P%QDEP
!
    DO L=0,LMAX_SOFT

       NQ=2
! find q values for simple scheme
       CALL AUG_SETQ(L,R,RMAXL(L),QQ(1:2),AA(1:2),LCOMPAT)

       IF (QMAX/=0) THEN
          CALL OPTREAL_AUG(L,NPSRNL,QMAX/NPSRNL/2, QMAX, RMAXL(L)*NPSNL/(NPSNL-1), NQ, QQ, AA(1))
!          WRITE(*,'(8F10.3)') QQ(1:NQ)
!          WRITE(*,'(8F10.3)') AA(1:NQ)
       ENDIF

! setup augmentation charge on radial grid rho(r) r^2
       DO N=1,R%NMAX
          SUM=0
          IF (R%R(N) <= RMAXL(L)) THEN
             DO I=1,NQ
                QR=QQ(I)*R%R(N)
                CALL SBESSEL( QR, BJ, L)
                SUM=SUM+BJ*AA(I)*R%R(N)*R%R(N)
             ENDDO
          ENDIF
          P%AUG_SOFT(N,L)=SUM
       ENDDO
! set up augmentation charge for spline interpolation on linear grid
       IF (L<=LMAX) THEN
       DO N=1,NPSRNL
          X=STEP*(N-1)
          SUM=0
          IF (X <= RMAXL(L)) THEN
             DO I=1,NQ
                QR=QQ(I)*X
                CALL SBESSEL( QR, BJ, L)
                SUM=SUM+BJ*AA(I)
             ENDDO
          ENDIF
          P%QDEP_FOCK(N,1,L) = X
          P%QDEP_FOCK(N,2,L) = SUM
       ENDDO
! derivative at startpoint
       X=STEP/1000
       SUM=0
       DO I=1,NQ
          QR=QQ(I)*X
          CALL SBESSEL( QR, BJ, L)
          SUM=SUM+BJ*AA(I)
       ENDDO
       SUM=(SUM-P%QDEP_FOCK(1,2,L))/X
! derivative is 0._q for all but L=1
       IF (L/=1) THEN
          SUM=0
       ENDIF
       CALL SPLCOF(P%QDEP_FOCK(1,1,L),NPSRNL,NPSRNL,SUM)
       ENDIF

    ENDDO

    NQ=2+NMAX_BESSEL

!
! upper panel: these are the 0._q normalized augmentation charges
!

    DO L=0,LMAX_BESSEL
       IF (QMAX>0) THEN
          WRITE(0,*) 'internal error SET_BESSEL_AUG: combination of LMAX_BESSEL and QMAX is not allowed', LMAX_BESSEL, QMAX
          CALL M_exit(); stop
       ENDIF
! find q values (either 1._q per function or two
! depending on whether the function is continous at outermost node or not
       CALL AUG_SETQ_NTH(L,R,RMAXL(L),QQ,A,NMAX_BESSEL,NQ)

       DO NBESSEL=1,NMAX_BESSEL
! setup augmentation charge on radial grid  rho(r) r^2
          DO N=1,R%NMAX
             SUM=0
             IF (R%R(N) <= RMAXL(L)) THEN
                DO I=1,NQ
                   QR=QQ(I)*R%R(N)
                   CALL SBESSEL( QR, BJ, L)
                   SUM=SUM+BJ*A(I,NBESSEL)*R%R(N)*R%R(N)
                ENDDO
             ENDIF
             P%AUG_FOCK(N, L, NBESSEL)=SUM
          ENDDO

! setup spline for augmentation charge
! this is stored consecutively in QDEP_FOCK
! after the original entries
          INDEX_BESSEL=LMAX+L+1+(NBESSEL-1)*(LMAX_BESSEL+1)

          DO N=1,NPSRNL
             X=STEP*(N-1)
             SUM=0
             IF (X <= RMAXL(L)) THEN
                DO I=1,NQ
                   QR=QQ(I)*X
                   CALL SBESSEL( QR, BJ, L)
                   SUM=SUM+BJ*A(I,NBESSEL)
                ENDDO
             ENDIF
             P%QDEP_FOCK(N,1,INDEX_BESSEL) = X
             P%QDEP_FOCK(N,2,INDEX_BESSEL) = SUM
          ENDDO

! derivative at startpoint
          X=STEP/1000
          SUM=0
          DO I=1,NQ
             QR=QQ(I)*X
             CALL SBESSEL( QR, BJ, L)
             SUM=SUM+BJ*A(I,NBESSEL)
          ENDDO
          SUM=(SUM-P%QDEP_FOCK(1,2,INDEX_BESSEL))/X
! derivative is 0._q for all but L=1
          IF (L/=1) THEN
             SUM=0
          ENDIF
          CALL SPLCOF(P%QDEP_FOCK(1,1,INDEX_BESSEL),NPSRNL,NPSRNL,SUM)
       ENDDO
    ENDDO

  END SUBROUTINE SET_BESSEL_AUG

!*******************************************************************
!
!  AUG_SETQ_NTQ
!  find a set of 1._q or two {q_i} and coefficients {A_i} such that
!     sum_i  A_i j_l(q_i,Rc) =0
!  and condition
!     sum_i  A_i d/d r j_l(q_i,Rc) =0
!  is fulfilled.
!  Furthermore the l moment of the functions is 1._q
!     sum_i int_0^Rc A_i j_l(q_i r) r^(l+2) dr = 1
!  and the functions are mutually orthogonal to each other
!  currently for up to two functions (NBESSEL) the coefficients
!  are returned in the array A_i
!  implementing more coefficients would be straight forward
!  in principle
!
!*******************************************************************

  SUBROUTINE AUG_SETQ_NTH(L,R,RMAX,QQ,A,NBESSEL,NQ)
    USE radial
    IMPLICIT NONE
    
    INTEGER  :: L
    TYPE (rgrid) R
    REAL(q)  :: RMAX     ! outermost radius for augmentation charges
    REAL(q)  :: QQ(:)
    REAL(q)  :: A (:,:)
    INTEGER  :: NBESSEL  ! "order" of required bessel functions
    INTEGER  :: NQ       ! number of required bessel functions
! local
    REAL(q)  AMAT(NQ,NQ), BMAT(NQ,NQ), NORM(NQ)
    INTEGER  IPIV(NQ)
! gaussian integration
    INTEGER, PARAMETER :: M=32
    INTEGER :: ITYPE, I, N, NB
    REAL(q) :: SUM,  INT, QR, BJ, BJP, BJPP
    REAL(q)  WR(M),RA(M)
    INTEGER IFAIL
    EXTERNAL GAUSSI2

    IF (NBESSEL>2) THEN
       WRITE(0,*) 'internal error in AUG_SETQ_NTH: presently only two functions are supported'
       CALL M_exit(); stop
    ENDIF
!-----------------------------------------------------------------------
! search for q so that j(q Rc)=0
!-----------------------------------------------------------------------
    CALL AUG_BEZERO(QQ,L,NQ)
    QQ=QQ/RMAX
!-----------------------------------------------------------------------
! set the matrix
!-----------------------------------------------------------------------
    ITYPE=0
    CALL GAUSSI(GAUSSI2,0.0_q,RMAX, ITYPE,M,WR,RA,IFAIL)
    AMAT=0
    DO I=1,NQ
! int j_l(qr) r^(2+l) dr
       SUM=0
       INT=0
       DO N=1,M
          CALL SBESSEL(QQ(I)*RA(N),BJ,L)
          SUM=SUM+BJ*RA(N)**(2+L)*WR(N)
          INT=INT+BJ*BJ*RA(N)**2*WR(N)
       ENDDO
! first entry moment
       AMAT(1,I) = SUM
       NORM(I)   = INT
       
! second entry derivative at boundary
       IF (NQ>1) THEN
          CALL SBESSE3( QQ(I)*RMAX, BJ, BJP, BJPP, L)
          BJP =BJP *QQ(I)
          BJPP=BJPP*QQ(I)*QQ(I)
          AMAT(2,I) = BJP
       ENDIF
    ENDDO
    BMAT=AMAT
!-----------------------------------------------------------------------
!  solve the linear equations  B_i = AMAT_ii' A_i'
!-----------------------------------------------------------------------
! derivative and r^(2+l) norm 0._q
    A(1,1)=0
    A(2,1)=0
! coefficient 3 equal 1
    AMAT(3,3)=1
    A(3,1)=1
    
    IFAIL=0
    CALL DGETRF( 3, 3, AMAT, NQ, IPIV, IFAIL )
    CALL DGETRS('N', 3, 1, AMAT, NQ, IPIV, A(1,1), NQ, IFAIL)
!-----------------------------------------------------------------------
!  second function
!-----------------------------------------------------------------------
    IF (NBESSEL>=2) THEN
       A(4,1)=0
       AMAT=BMAT
! derivative and r^(2+l) norm 0._q
       A(1,2)=0
       A(2,2)=0
! function orthogonal to previous function
       A(3,2)=0
       AMAT(3,:)=A(:,1)*NORM(:)
! coefficient 4 equal 1
       A(4,2)=1
       AMAT(4,4)=1
       
       IFAIL=0
       CALL DGETRF( 4, 4, AMAT, NQ, IPIV, IFAIL )
       CALL DGETRS('N', 4, 1, AMAT, NQ, IPIV, A(1,2), NQ, IFAIL)
    ENDIF
!-----------------------------------------------------------------------
!  normalize the functions properly
!-----------------------------------------------------------------------
! norm of squared function should be 1._q
! so that it can be used like a projector function
    DO NB=1,NBESSEL
       INT=0
       DO N=1,M
          SUM=0
          DO I=1,NQ
             CALL SBESSEL(QQ(I)*RA(N),BJ,L)
             SUM=SUM+BJ*A(I,NB)
          ENDDO
          INT=INT+SUM*SUM*RA(N)**2*WR(N)
       ENDDO
       A(:,NB)=A(:,NB)/SQRT(INT)
    ENDDO
# 1001

  END SUBROUTINE AUG_SETQ_NTH

!*******************************************************************
!
!  RAD_QPAW_FOCK
!  subroutine calculates the difference between the AE and PS
!  charge and develops this difference into a set of Besselfunctions
!  the result is stored in QPAW_FOCK
!  it allows to accurately restore the AE charge density on
!  the plane wave grid up to a certain cutoff
!
!*******************************************************************

  SUBROUTINE RAD_QPAW_FOCK( R, CHANNELS, RMAXL, QMAX_FOCKAE, WAE, WPS , AUG, AUG_FOCK, L, QPAW, & 
                                QPAW_FOCK, LMAX_BESSEL, NMAX_BESSEL, IU6 )
    USE radial
    IMPLICIT NONE

    TYPE (rgrid),INTENT(IN) :: R
    INTEGER,INTENT(IN) :: CHANNELS
    REAL(q),INTENT(IN) :: RMAXL(0:), QMAX_FOCKAE
    REAL(q),INTENT(IN) :: WAE(:,:),WPS(:,:)   ! AE and soft wavefunctions
    REAL(q),INTENT(IN) :: AUG(:,0:)           ! 1-normalized L-dep compensation charge
    REAL(q),INTENT(IN) :: AUG_FOCK(:,0:,:)    ! additional accurate shape restoring compensation charges
    REAL(q),INTENT(IN) :: QPAW(:,:,0:)        ! moments of compensation charge
    REAL(q),INTENT(OUT) :: QPAW_FOCK(:,:,0:,:) ! moments of compensation charge
    INTEGER,INTENT(IN) :: L(:)
    INTEGER,INTENT(IN) :: LMAX_BESSEL         ! maximum L quantum number of Bessel projection
    INTEGER,INTENT(IN) :: NMAX_BESSEL         ! maximum N quantum number of Bessel projection
! local variables
    REAL(q) :: RHOT(R%NMAX),RHOJ(R%NMAX),POT(R%NMAX)
    INTEGER CH1, CH2, I, NBESSEL
    REAL(q) :: RES, INT, INTPS, BJ
    INTEGER :: LL,LLP,LMIN,LMAX,LMAIN,N
! decides which version is used
    INTEGER :: IVERSION=1, IFIT, NFIT
    INTEGER :: IU6
    REAL(q) :: QFIT(NMAX_BESSEL) 
    REAL(q) :: AMAT(NMAX_BESSEL,NMAX_BESSEL),A(NMAX_BESSEL),ABEST(NMAX_BESSEL)
    REAL(q) :: ERROR,SUM_ERROR,TOTAL,MAX_ERROR, ERRORQ, MAX_ERRORQ
    INTEGER :: IPIV(NMAX_BESSEL)


    QPAW_FOCK=0
    SUM_ERROR=0
    MAX_ERROR=0
    MAX_ERRORQ=0
    TOTAL=0
    DO CH1=1,CHANNELS
       DO CH2=CH1,CHANNELS
! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         DO LMAIN=LMIN,LMAX,2
            IF (LMAIN > LMAX_BESSEL) CYCLE
            DO I=1,R%NMAX
               RHOT(I)=(WAE(I,CH1)*WAE(I,CH2)-WPS(I,CH1)*WPS(I,CH2)-QPAW(CH1,CH2,LMAIN)*AUG(I,LMAIN))*R%R(I)**LMAIN
            ENDDO
            CALL SIMPI(R, RHOT , RES)
! screwed if we do not have 0._q charge
            IF ( ABS(RES) > 1E-4 ) THEN
               WRITE(0,1) CH1,CH2,RES,QPAW(CH1,CH2,LMAIN)
 1             FORMAT('internal error RAD_QPAW_FOCK: QPAW is incorrect',/ &
               '      channels',2I3,' QPAW=',E20.10,' int=',10E20.10)
               CALL M_exit(); stop
            ENDIF
! calculate now without weight
            DO I=1,R%NMAX
               RHOT(I)=(WAE(I,CH1)*WAE(I,CH2)-WPS(I,CH1)*WPS(I,CH2)-QPAW(CH1,CH2,LMAIN)*AUG(I,LMAIN))
            ENDDO
!-----------------------------------------------------------------------
! develop RHOT into functions stored in AUG_FOCK
! VERSION I
! fit the Besseltransform of the AE-PS charge density at 1._q or
! two selected q values
! this might be changed to a least square fit, but 1._q would need
! to decide for a target weighting function. This could be based
! on a simple model for the dielectric function, but I hardly believe
! this is required
!-----------------------------------------------------------------------
            IF (IVERSION==1) THEN
            ERROR=100
            NFIT=-1
!            DO IFIT=0,(NMAX_BESSEL-1)*100
            DO IFIT=0,0
               IF (NMAX_BESSEL==1) THEN
! at 6 A-1 typically the 1-epsilon^-1 is smaller than 0.05
                  QFIT=6.0       !  6 A-1 =  140 eV
                  IF (QMAX_FOCKAE>0) THEN
                     QFIT=QMAX_FOCKAE
                  ENDIF
               ELSE
! at 10 A-1 there is usually hardly any RPA correlation left
                  QFIT(1)=5.0       !  5 A-1 = 95 eV
                  IF (QMAX_FOCKAE>0) THEN
                     QFIT(1)=QMAX_FOCKAE
                  ENDIF
                  QFIT(2)=QFIT(1)*2 ! 10 A-1 = 380 eV
               ENDIF

! function RHOT contains already r^2
               DO N=1,NMAX_BESSEL
                  DO I=1,R%NMAX
                     CALL SBESSEL(  R%R(I)*QFIT(N), BJ, LMAIN)
                     RHOJ(I)=RHOT(I)*BJ
                  ENDDO
                  CALL SIMPI(R, RHOJ , A(N))
                  DO NBESSEL=1,NMAX_BESSEL
                     DO I=1,R%NMAX
                        CALL SBESSEL(  R%R(I)*QFIT(N), BJ, LMAIN)
                        RHOJ(I)=AUG_FOCK(I, LMAIN, NBESSEL)*BJ
                     ENDDO
                     CALL SIMPI(R, RHOJ , AMAT(N, NBESSEL))
                  ENDDO
               ENDDO
               I=0
! for a single element that is the result
!              A(1)=A(1)/AMAT(1,1)
               CALL DGETRF( NMAX_BESSEL, NMAX_BESSEL, AMAT, NMAX_BESSEL, IPIV, I )
               CALL DGETRS('N', NMAX_BESSEL, 1, AMAT, NMAX_BESSEL, IPIV, A, NMAX_BESSEL, I)
               QPAW_FOCK(CH1, CH2, LMAIN, :)=A
               QPAW_FOCK(CH2, CH1, LMAIN, :)=A

               DO I=1,R%NMAX
                  RHOJ(I)=WAE(I,CH1)*WAE(I,CH2)
               ENDDO
               CALL RAD_POT_HAR(LMAIN,R,POT(:),RHOJ(:),INT)
               
               DO I=1,R%NMAX
                  RHOJ(I)=(WPS(I,CH1)*WPS(I,CH2)+QPAW(CH1,CH2,LMAIN)*AUG(I,LMAIN))+ & 
                  SUM(QPAW_FOCK(CH1, CH2, LMAIN, :)* &
                  AUG_FOCK(I, LMAIN, :))
               ENDDO
               CALL RAD_POT_HAR(LMAIN,R,POT(:),RHOJ(:),INTPS)
               IF (ABS(INT-INTPS)<ERROR) THEN
                  ABEST=A
                  NFIT=IFIT
                  ERROR=ABS(INT-INTPS)
               ENDIF
            ENDDO
            QPAW_FOCK(CH1, CH2, LMAIN, :)=ABEST
            QPAW_FOCK(CH2, CH1, LMAIN, :)=ABEST
# 1148

            SUM_ERROR=SUM_ERROR+ERROR
            MAX_ERROR=MAX(MAX_ERROR,ABS(ERROR))
            TOTAL=TOTAL+1
!-----------------------------------------------------------------------
! VERSION II
! the functions stored in AUG_FOCK have the property to be
! orthogonal w.r.t the weight  r^2
! int aug_n(r) aug_n'(r) r^2 dr = delta n n'
! use them as projector functions
! this might be more suitable for Hartree-Fock, but does not work so
! well at low energies (and that is what we need for GW
! and to a lesser extend RPA correlation energies)
!-----------------------------------------------------------------------
            ELSE
            DO NBESSEL=1,NMAX_BESSEL
! AUG_FOCK stores aug_n r^2
! the factor r^2 is already included in RHOT hence divide by r^2
               RHOJ(:)=RHOT(:)*AUG_FOCK(:, LMAIN, NBESSEL)/R%R(:)**2
               CALL SIMPI(R, RHOJ , RES)

               QPAW_FOCK(CH1, CH2, LMAIN, NBESSEL)=RES
               QPAW_FOCK(CH2, CH1, LMAIN, NBESSEL)=RES
! remove projected part (debugging below requires this)
               RHOT(:)=RHOT(:)-RES*AUG_FOCK(:, LMAIN, NBESSEL)
            ENDDO
            ENDIF

            DO I=1,R%NMAX
               RHOJ(I)=WAE(I,CH1)*WAE(I,CH2)- & 
                  ((WPS(I,CH1)*WPS(I,CH2)+QPAW(CH1,CH2,LMAIN)*AUG(I,LMAIN))+ & 
                    SUM(QPAW_FOCK(CH1, CH2, LMAIN, :)* &
                        AUG_FOCK(I, LMAIN, :)))
            ENDDO
            CALL RAD_QPAW_FOCK_ERROR( R, LMAIN, RHOJ, QFIT(NMAX_BESSEL), ERRORQ)


            MAX_ERRORQ=MAX(MAX_ERRORQ,ABS(ERRORQ))
# 1221


         ENDDO
       ENDDO
    ENDDO
    IF (IU6>=0) THEN
       WRITE(IU6,'( " error exchange integrals (PW-AE in eV) (mean, max)",2F6.3, " density error[0,qmax] ",F6.4)') SUM_ERROR/TOTAL,MAX_ERROR,MAX_ERRORQ
    ENDIF
  END SUBROUTINE RAD_QPAW_FOCK
      

! this small function develops a function into a set of spherical
! besselfunction
  SUBROUTINE RAD_QPAW_FOCK_HELPER( R, L, RHOT)
    USE radial
    IMPLICIT NONE

    TYPE (rgrid) R
    INTEGER  :: L
    REAL(q)  :: RHOT(R%NMAX)
! loc
    REAL(q) :: RHOJ(R%NMAX)
    INTEGER, PARAMETER :: NBESSEL=100
    REAL(q), PARAMETER :: QMAX=25.0
    REAL(q)  :: QQ(NBESSEL), INT, BJ
    INTEGER  :: N, I

    CALL AUG_BEZERO(QQ,L,NBESSEL)

    QQ=QQ/QQ(NBESSEL)*QMAX
    DO I=1,R%NMAX
       WRITE(78,*) R%R(I), RHOT(I)
    ENDDO
    WRITE(78,*)

    DO N=1, NBESSEL
! function RHOT contains already r^2
       DO I=1,R%NMAX
         CALL SBESSEL(  R%R(I)*QQ(N), BJ, L)
         RHOJ(I)=RHOT(I)*BJ
      ENDDO

       CALL SIMPI(R, RHOJ , INT)
       WRITE(77,*) QQ(N), INT
    END DO
    WRITE(77,*)
    
  END SUBROUTINE RAD_QPAW_FOCK_HELPER


  SUBROUTINE RAD_QPAW_FOCK_ERROR( R, L, RHOT, QMAX, ERROR)
    USE radial
    IMPLICIT NONE

    TYPE (rgrid) R
    INTEGER  :: L
    REAL(q)  :: RHOT(R%NMAX)
! loc
    REAL(q) :: RHOJ(R%NMAX)
    INTEGER, PARAMETER :: NBESSEL=100
    REAL(q)  :: QQ(NBESSEL), INT, BJ
    REAL(q)  :: QMAX, ERROR
    INTEGER  :: N, I

    CALL AUG_BEZERO(QQ,L,NBESSEL)

    QQ=QQ/QQ(NBESSEL)*QMAX
    
    ERROR=0

    DO N=1, NBESSEL
! function RHOT contains already r^2
       DO I=1,R%NMAX
         CALL SBESSEL(  R%R(I)*QQ(N), BJ, L)
         RHOJ(I)=RHOT(I)*BJ
      ENDDO

       CALL SIMPI(R, RHOJ , INT)
       ERROR=MAX(ERROR,ABS(INT))
    END DO
    
  END SUBROUTINE RAD_QPAW_FOCK_ERROR
      
!*******************************************************************
!
!  calculate the transformation matrix that
!  is required to go from the occupancies RHO(lm, l'm') to
!  the onsite L,M dependent occupancies RHO(LM)
!  where RHO(lm,l'm')  are the occupancies of each channel
!  and RHO(LM) is given by
!
!  RHO(LM) =  sum C(LM,ll',mm')  RHO(lm,l'm') QPAW(llp)
!
!  LM is a shorthand for the combined index (L+1)*(L+1)+M+1
!  C(LM,ll',mm') are the Clebsch Gordan coefficients
!  to avoid the expensive complicated loop constructs this
!  subroutines stores the explicit transformation matrix
!
!    TRANS_MATRIX(:,:,LM,ITYP)
!  so that
!  RHO(LM)= sum TRANS_MATRIX(lm,l'm',LM,type) RHO(lm,l'm')
!
!*******************************************************************

  SUBROUTINE SETUP_TRANS_MATRIX( P, TRANS_MATRIX, LMAX, LMAX_FOCKAE, NMAX_FOCKAE)
    USE pseudo
    USE asa
    IMPLICIT NONE
    TYPE (potcar) P        ! PP descriptor
    REAL(q) :: TRANS_MATRIX(:, :, :)
    INTEGER :: LMAX, LMAX_FOCKAE, NMAX_FOCKAE
! local varible
    INTEGER CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP, LMMAX
    INTEGER NBESSEL, DIM3

! initialize everything to 0
    LMMAX=(LMAX+1)*(LMAX+1)
    TRANS_MATRIX(:,:,:)=0

! short cut LMAX=-1
    IF (LMMAX==0) RETURN

    DIM3=SIZE( TRANS_MATRIX, 3)

! loop over all channels (l,epsilon)
    LM=1
    DO CH1=1,P%LMAX
    LMP=1
    DO CH2=1,P%LMAX

! quantum numbers l and lp of these two channels
      LL =P%LPS(CH1)
      LLP=P%LPS(CH2)

      CALL YLM3LOOKUP(LL,LLP,LMINDX)

      DO M =1,2*LL+1
      DO MP=1,2*LLP+1
         LMINDX=LMINDX+1

         ISTART=INDCG(LMINDX)
         IEND  =INDCG(LMINDX+1)

         DO  IC=ISTART,IEND-1
            IF (JS(IC) <= LMMAX ) THEN
             TRANS_MATRIX(LM+M-1,LMP+MP-1,JS(IC))= &
                   TRANS_MATRIX(LM+M-1,LMP+MP-1,JS(IC))+ &
                   YLM3(IC)*P%QPAW(CH1,CH2,JL(IC))
            ENDIF
! base index of JC is 1._q
            IF (JS(IC) <= (LMAX_FOCKAE+1)*(LMAX_FOCKAE+1)) THEN
               IF (JS(IC)+LMMAX+(LMAX_FOCKAE+1)*(LMAX_FOCKAE+1)*(NMAX_FOCKAE-1)>DIM3) THEN
                  WRITE(0,*) 'internal error in  SETUP_TRANS_MATRIX: size insufficient', & 
                       JS(IC)+LMMAX+(LMAX_FOCKAE+1)*(LMAX_FOCKAE+1)*(NMAX_FOCKAE-1),DIM3
               ENDIF
               DO NBESSEL=1,NMAX_FOCKAE
                 TRANS_MATRIX(LM+M-1,LMP+MP-1,JS(IC)+LMMAX+(LMAX_FOCKAE+1)*(LMAX_FOCKAE+1)*(NBESSEL-1))= &
                 TRANS_MATRIX(LM+M-1,LMP+MP-1,JS(IC)+LMMAX+(LMAX_FOCKAE+1)*(LMAX_FOCKAE+1)*(NBESSEL-1))+ &
                   YLM3(IC)*P%QPAW_FOCK(CH1,CH2, JL(IC), NBESSEL )
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      ENDDO

      LMP=LMP+2*LLP+1
    ENDDO
    LM =LM +2*LL +1
    ENDDO
!    WRITE(*,*) SIZE(TRANS_MATRIX,1),LMAX_FOCKAE,NMAX_FOCKAE
!    WRITE(*,*) 'TRANS_MATRIX dump'
!    WRITE(*,'(8(8F8.4/))')  TRANS_MATRIX(:,:, 1:LMMAX+(LMAX_FOCKAE+1)*(LMAX_FOCKAE+1)*NMAX_FOCKAE)

  END SUBROUTINE SETUP_TRANS_MATRIX

!*******************************************************************
!
! the 1._q center handle allows to retrieve selected
! 1._q center occupancies in order to treat the 1._q-center terms
! exactly
!
!*******************************************************************

  SUBROUTINE SET_UP_ONE_CENTER_H( WDES, P, T_INFO, LMAX_CALC, H)
    USE pseudo
    USE wave
    USE asa
    USE pawfock_inter
    IMPLICIT NONE
    TYPE (wavedes) :: WDES ! wave function descriptor
    TYPE (type_info) T_INFO
    TYPE (potcar), TARGET :: P(:)     ! P descriptors
    INTEGER :: LMAX_CALC   ! maximum L
    TYPE (one_center_handle), POINTER :: H
! local varible
    TYPE (potcar), POINTER :: PP      ! P descriptors
    INTEGER :: CH1,CH2,LL,LLP,LM,LMP,LMIN,LMAX,LMINDX,ISTART,IEND,IC,M,MP
    INTEGER :: LMAIN, LMMAIN, MMAIN
    INTEGER :: NT, NTP, NIS, NI, NALLOC, BASE_INDEX
    REAL(q), ALLOCATABLE :: CG(:)
    INTEGER, ALLOCATABLE :: INDEX(:)
    INTEGER :: INDEX_INTO_CG, ENTRIES_FOR_TYPE
! maximum possible index for this array (2L+1)*(2L+1)
! L is the maximum 1._q center angular momentum number
    IF (LMAX_CALC<0) RETURN

    ALLOCATE(H)
!-----------------------------------------------------------------------
! first part set up the required handle
!-----------------------------------------------------------------------
    NALLOC=MAXVAL(WDES%LMMAX(:))**2
! this is a little bit generous and certainly suffices
    ALLOCATE(H%ENTRIES_FOR_TYPE(WDES%NTYP), & 
             H%NR_OF_CG(NALLOC, WDES%NTYP), &
             H%CHANNEL1(NALLOC, WDES%NTYP), H%CHANNEL2(NALLOC, WDES%NTYP), &  
             H%L(NALLOC, WDES%NTYP), H%LM(NALLOC, WDES%NTYP))

    H%POSION=>T_INFO%POSION
    H%ENTRIES_FOR_TYPE=0
    H%NR_OF_CG        =0
    
! I think this should suffice once and for ever to temporary store CG and INDEX
    NALLOC=NALLOC*NALLOC*WDES%NTYP*(LMAX_CALC+1)*(LMAX_CALC+1)
    ALLOCATE(CG(NALLOC), INDEX(NALLOC))
   
    H%LMAX_ONE_CENTER=LMAX_CALC
    H%TOTAL_ENTRIES=0
    H%TOTAL_ENTRIES_CHI=0

    NIS=1
    INDEX_INTO_CG=0

    DO NT=1,WDES%NTYP
! determine global type index for this ion
      NTP=T_INFO%ITYP(NI_GLOBAL(NIS, WDES%COMM_INB))

      ENTRIES_FOR_TYPE=0

! loop over all channels CH1,CH2= (l,epsilon),(l,epsilon)'
      LM=1
      DO CH1=1,P(NTP)%LMAX
      LMP=LM
      DO CH2=CH1,P(NTP)%LMAX

! quantum numbers l and lp of these two channels
         LL =P(NTP)%LPS(CH1)
         LLP=P(NTP)%LPS(CH2)

! now collect contributions for linearly stored
! Clebsch Gordan transformed charge densities
         LMIN=ABS(LL-LLP) ; LMAX=MIN(LMAX_CALC,P(NTP)%LMAX_CALC,ABS(LL+LLP))
! due to sum rules only LMIN to LMAX in steps of 2 are allowed
         DO LMAIN=LMIN,LMAX,2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            ENTRIES_FOR_TYPE=ENTRIES_FOR_TYPE+1
            IF (ENTRIES_FOR_TYPE > SIZE(H%NR_OF_CG,1)) THEN
               WRITE(0,*) 'internal error in SET_UP_ONE_CENTER_H: SIZE(H%NR_OF_CG) not sufficient'
               CALL M_exit(); stop
            ENDIF
! store the channel numbers
            H%CHANNEL1(ENTRIES_FOR_TYPE,NT)=CH1
            H%CHANNEL2(ENTRIES_FOR_TYPE,NT)=CH2
            H%LM      (ENTRIES_FOR_TYPE,NT)=LMMAIN
            H%L       (ENTRIES_FOR_TYPE,NT)=LMAIN

! now loop over all M and MP indices
! and search for those contributions that contribute to
! this transformed entry
            CALL YLM3LOOKUP(LL,LLP,LMINDX)

            DO M =1,2*LL+1
            DO MP=1,2*LLP+1
               LMINDX=LMINDX+1

               ISTART=INDCG(LMINDX)
               IEND  =INDCG(LMINDX+1)

               DO  IC=ISTART,IEND-1
                  IF (JS(IC) == LMMAIN) THEN
                     H%NR_OF_CG(ENTRIES_FOR_TYPE,NT)=H%NR_OF_CG(ENTRIES_FOR_TYPE,NT)+1
                     INDEX_INTO_CG=INDEX_INTO_CG+1
                     IF (INDEX_INTO_CG > NALLOC) THEN
                        WRITE(0,*) 'internal error in SET_UP_ONE_CENTER_H: dimension of CG and INDEX not sufficient'
                        CALL M_exit(); stop
                     ENDIF
! store CG coefficient
                     CG(INDEX_INTO_CG)=YLM3(IC)
! convert index LM+M-1,LMP+MP-1 into 1._q dim index
                     INDEX(INDEX_INTO_CG)=LM+M-1+(LMP+MP-1-1)*WDES%LMMAX(NT)

! the interchanged channels contribute to same element
! after CG transformation
                     IF (CH1/=CH2) THEN
                        H%NR_OF_CG(ENTRIES_FOR_TYPE,NT)=H%NR_OF_CG(ENTRIES_FOR_TYPE,NT)+1
                        INDEX_INTO_CG=INDEX_INTO_CG+1
                        IF (INDEX_INTO_CG > NALLOC) THEN
                           WRITE(0,*) 'internal error in SET_UP_ONE_CENTER_H: dimension of CG and INDEX not sufficient'
                           CALL M_exit(); stop
                        ENDIF
! store CG coefficient
                        CG(INDEX_INTO_CG)=YLM3(IC)
! convert index LM+M-1,LMP+MP-1 into 1._q dim index
                        INDEX(INDEX_INTO_CG)=(LM+M-1-1)*WDES%LMMAX(NT)+LMP+MP-1
                     ENDIF

                  ENDIF
               ENDDO
            ENDDO
            ENDDO

         ENDDO
         ENDDO
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO

      H%ENTRIES_FOR_TYPE(NT)=ENTRIES_FOR_TYPE
      H%TOTAL_ENTRIES=H%TOTAL_ENTRIES+H%ENTRIES_FOR_TYPE(NT)*WDES%NITYP(NT)
      H%TOTAL_ENTRIES_CHI=H%TOTAL_ENTRIES_CHI+H%ENTRIES_FOR_TYPE(NT)*H%ENTRIES_FOR_TYPE(NT)*WDES%NITYP(NT)
! increase ion index
      NIS = NIS+WDES%NITYP(NT)
    ENDDO

! make TOTAL_ENTRIES_CHI dividable by 2
    H%TOTAL_ENTRIES_CHI=MOD(H%TOTAL_ENTRIES_CHI+1,2)
!-----------------------------------------------------------------------
! here a little comment applies
! the CG coefficients are all  1/sqrt(3.1415*4) and every element
! in the CRHODE matrix contributes only to a single element in
! the transformed density matrix
! 1._q could use this to speed up the code significantly
!-----------------------------------------------------------------------

    ALLOCATE(H%CG(INDEX_INTO_CG),H%INDEX(INDEX_INTO_CG))
    H%CG(1:INDEX_INTO_CG)   =CG   (1:INDEX_INTO_CG)
    H%INDEX(1:INDEX_INTO_CG)=INDEX(1:INDEX_INTO_CG)
! deallocate temporary array
    DEALLOCATE(CG, INDEX)
!-----------------------------------------------------------------------
! second part set up the matrix elements in the potential arrays
!-----------------------------------------------------------------------
! just in case
    CALL RELEASE_PAWFOCK

    NALLOC=SUM(H%ENTRIES_FOR_TYPE(:)*H%ENTRIES_FOR_TYPE(:))
    ALLOCATE(H%POT(NALLOC))
    H%POT=0
    
    NIS=1
    BASE_INDEX    =1
    typ:  DO NT=1,WDES%NTYP
! determine global type index for this ion type
       NTP=T_INFO%ITYP(NI_GLOBAL(NIS, WDES%COMM_INB))
       PP=>P(NTP)
       CALL SETUP_PAWFOCK(NTP, PP)

       CALL SETUP_PAWFOCK_MATRIX( H%ENTRIES_FOR_TYPE(NT), H%CHANNEL1(:,NT), H%CHANNEL2(:,NT), H%L(:,NT), H%LM(:,NT), & 
            PP,H%POT(BASE_INDEX:BASE_INDEX+H%ENTRIES_FOR_TYPE(NT)*H%ENTRIES_FOR_TYPE(NT)-1))
            
       BASE_INDEX=BASE_INDEX+H%ENTRIES_FOR_TYPE(NT)*H%ENTRIES_FOR_TYPE(NT)
! increase ion index
       NIS = NIS+WDES%NITYP(NT)
    ENDDO typ
    CALL RELEASE_PAWFOCK
  END SUBROUTINE SET_UP_ONE_CENTER_H

!*******************************************************************
!
! deallocate 1._q center handle(WDES%NTYP)
!
!*******************************************************************

  SUBROUTINE DEALLOCATE_ONE_CENTER_H( H)
    TYPE (one_center_handle), POINTER :: H
    
    IF (ASSOCIATED(H)) THEN
       DEALLOCATE(H%ENTRIES_FOR_TYPE, & 
             H%NR_OF_CG, &
             H%CHANNEL1, H%CHANNEL2, H%L, H%LM)

       DEALLOCATE(H%CG,H%INDEX)

       DEALLOCATE(H)
       
    ENDIF
  END SUBROUTINE DEALLOCATE_ONE_CENTER_H


!***********************************************************************
!
!  apply the phase factor e^(i q R_i) to a linear CG transformed
!  occupancy matrix.
!  This is required since the four orbital integrals
!  and the occupancies matrices might be set off by a reciprocal
!  lattice vector
!   e^(i G r)
!  hence a shift of e^(i G R_ion) is missing in CRHO
!  The algebra is fairly easy to work out. One simply needs
!  to work out what happens if 1._q of the wavefunctions shifts by
!  e^(i G r), and how this changes the corresponding GPROJ
!  To test this, shift all atoms by a constant vector and check
!  whether results are invariant under the shift
!
!***********************************************************************

  SUBROUTINE APPLY_PHASE_ONE_CENTER(WDES, H, CRHO, VKPT )
    USE prec
    USE constant
    IMPLICIT NONE

    TYPE (wavedes) WDES
    TYPE (one_center_handle) :: H
    REAL(q) :: CRHO(:,:)       ! 1._q center charge density for a number of vectors
    REAL(q) :: VKPT(3)
! local
    REAL(q),PARAMETER :: TINY=1E-4_q
    COMPLEX (q) CSHIFT
    INTEGER NIS, NRHO, NT, NI, NIP

    IF (ABS(VKPT(1))>TINY .OR. ABS(VKPT(2))>TINY .OR. ABS(VKPT(3))>TINY) THEN
       NRHO=1
       NIS=1
       DO NT=1,WDES%NTYP
       
          DO NI=NIS,WDES%NITYP(NT)+NIS-1
             NIP=NI_GLOBAL(NI,WDES%COMM_INB)
             
             CSHIFT=EXP(CITPI*(VKPT(1)*H%POSION(1,NIP)+ & 
                               VKPT(2)*H%POSION(2,NIP)+ &
                               VKPT(3)*H%POSION(3,NIP)))
             
             CRHO(NRHO:NRHO+H%ENTRIES_FOR_TYPE(NT)-1,:)=CRHO(NRHO:NRHO+H%ENTRIES_FOR_TYPE(NT)-1,:)*CSHIFT
             NRHO = NRHO+H%ENTRIES_FOR_TYPE(NT)
          ENDDO
          
          NIS = NIS+WDES%NITYP(NT)
       ENDDO
    ENDIF
  END SUBROUTINE APPLY_PHASE_ONE_CENTER


!*******************************************************************
!
! apply 1._q-center potential to a set of 1._q center occupancies
!
!*******************************************************************

  SUBROUTINE APPLY_ONE_CENTER_H( WDES, H, CRHO, CRHO_TMP, N)
    USE wave
    IMPLICIT NONE

    TYPE (wavedes) WDES
    TYPE (one_center_handle) :: H
    REAL(q) :: CRHO(:,:)       ! 1._q center charge density for a number of vectors
    REAL(q) :: CRHO_TMP(:,:)   ! temporary work array
    INTEGER :: N            ! number of vectors in CRHO and CRHO_TMP
! local
    INTEGER :: NIS, NT, NI
    INTEGER :: NPOT, NRHO

    NIS =1
    NRHO=1
    NPOT=1
    DO NT=1,WDES%NTYP
       
       DO NI=NIS,WDES%NITYP(NT)+NIS-1
          
          CALL DGEMM('N','N', H%ENTRIES_FOR_TYPE(NT), N, H%ENTRIES_FOR_TYPE(NT), 1._q, &
               H%POT(NPOT), H%ENTRIES_FOR_TYPE(NT), CRHO(NRHO,1), SIZE(CRHO,1), &
               0._q, CRHO_TMP(NRHO,1), SIZE(CRHO,1) )

          NRHO = NRHO+H%ENTRIES_FOR_TYPE(NT)
       ENDDO

    NPOT= NPOT+H%ENTRIES_FOR_TYPE(NT)*H%ENTRIES_FOR_TYPE(NT)
    NIS = NIS+WDES%NITYP(NT)
    ENDDO

    IF (WDES%LGAMMA) THEN
       CALL DCOPY( SIZE(CRHO,1)*N,  CRHO_TMP(1,1), 1,  CRHO(1,1), 1)
    ELSE
       CALL ZCOPY( SIZE(CRHO,1)*N,  CRHO_TMP(1,1), 1,  CRHO(1,1), 1)
    ENDIF

  END SUBROUTINE APPLY_ONE_CENTER_H


  SUBROUTINE APPLY_ONE_CENTER_H_NOOVERWRT( WDES, H, CRHO, CRHO_TMP, N)
    USE wave
    IMPLICIT NONE

    TYPE (wavedes) WDES
    TYPE (one_center_handle) :: H
    REAL(q) :: CRHO(:,:)       ! 1._q center charge density for a number of vectors
    REAL(q) :: CRHO_TMP(:,:)   ! temporary work array
    INTEGER :: N            ! number of vectors in CRHO and CRHO_TMP
! local
    INTEGER :: NIS, NT, NI
    INTEGER :: NPOT, NRHO

    NIS =1
    NRHO=1
    NPOT=1
    DO NT=1,WDES%NTYP
       
       DO NI=NIS,WDES%NITYP(NT)+NIS-1
          
          CALL DGEMM('N','N', H%ENTRIES_FOR_TYPE(NT), N, H%ENTRIES_FOR_TYPE(NT), 1._q, &
               H%POT(NPOT), H%ENTRIES_FOR_TYPE(NT), CRHO(NRHO,1), SIZE(CRHO,1), &
               0._q, CRHO_TMP(NRHO,1), SIZE(CRHO,1) )

          NRHO = NRHO+H%ENTRIES_FOR_TYPE(NT)
       ENDDO

    NPOT= NPOT+H%ENTRIES_FOR_TYPE(NT)*H%ENTRIES_FOR_TYPE(NT)
    NIS = NIS+WDES%NITYP(NT)
    ENDDO


  END SUBROUTINE APPLY_ONE_CENTER_H_NOOVERWRT
      

!************************* SUBROUTINE CALC_RHOLM_TRANS ****************
!
!  calculate net LM moment of the soft augmentation charges
!  (compensation charges in PAW language) RHO(LM) at all sites
!  this routine is similar to the routine CALC_RHOLM in the paw.F routine
!  but tuned for efficiency
!  by using the stored transformation matrix
!
!  RHO(LM)= sum TRANS_MATRIX(lm,l'm',LM,type) RHO(lm,l'm')
!
!**********************************************************************

    SUBROUTINE CALC_RHOLM_TRANS(WDES, AUG_DES, TRANS_MATRIX, CRHODE, CRHOLM,  FACT)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)  WDES
      TYPE (wavedes1) AUG_DES
      REAL(q) :: TRANS_MATRIX(:,:,:,:) 
      REAL(q) :: CRHODE(:,:,:,:)
      REAL(q) :: CRHOLM(:)
      REAL(q) :: FACT          ! arbitrary multiplication factor
! local
      INTEGER ISPINOR, ISPINOR_, NIS, NT, LMMAXC, NI, LM, L, LP, NPRO, NPRO_, NAUG

      CRHOLM=0 

      spinor: DO ISPINOR=0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      NAUG= (ISPINOR_+2*ISPINOR)*AUG_DES%NPRO+1

      NIS =1
      DO NT=1,WDES%NTYP
        LMMAXC=WDES%LMMAX(NT)
        IF (LMMAXC==0) GOTO 100
        DO NI=NIS,WDES%NITYP(NT)+NIS-1

          DO LM=1,AUG_DES%LMMAX(NT)
          DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
          DO LP=1,LMMAXC
             CRHOLM(NAUG+LM-1)=CRHOLM(NAUG+LM-1)+ &
                 CRHODE(LP,L,NI,ISPINOR_+2*ISPINOR+1)*TRANS_MATRIX(LP,L,LM,NT)*FACT
          ENDDO
          ENDDO
          ENDDO

          NAUG = NAUG+ AUG_DES%LMMAX(NT)
        ENDDO
 100    NIS = NIS+WDES%NITYP(NT)
      ENDDO

      ENDDO
      ENDDO spinor

      RETURN
    END SUBROUTINE CALC_RHOLM_TRANS

!************************ SUBROUTINE DEPSUM_TWO_BANDS_RHOLM ************
!
! this subroutine is a combination of DEPSUM_TWO_BANDS
! and CALC_RHOLM_TRANS
! for efficiency reason 1._q should use this subroutine
! this version returns the charge from the up and down channel
! seperately
! (not used presently in VASP)
!
!***********************************************************************

    SUBROUTINE DEPSUM_TWO_BANDS_RHOLM( CPROJ1, CPROJ2, WDES, AUG_DES, & 
         TRANS_MATRIX, CRHOLM, FACT, LOVERL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1) WDES
      TYPE (wavedes1) AUG_DES
      LOGICAL LOVERL
      REAL(q) :: TRANS_MATRIX(:, :, :, :)
      REAL(q)  :: CPROJ1(:)
      REAL(q)  :: CPROJ2(:)
      REAL(q)  :: CTMP
      REAL(q)  :: CRHOLM(:)
      REAL(q) :: FACT              ! arbitrary multiplication factor
! local
      INTEGER :: ISPINOR, LMBASE, NT, LMMAXC, NIS, NI, L, LP, LM, NAUG
      
      IF (.NOT.LOVERL) RETURN

      CRHOLM=0

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)+1
      NAUG=   ISPINOR *AUG_DES%NPRO+1

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 100

         ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1
            CALL RHOLM_KERNEL(CPROJ1(LMBASE), CPROJ2(LMBASE), &
                 TRANS_MATRIX(1,1,1,NT), SIZE(TRANS_MATRIX,1), &
                 CRHOLM(NAUG), FACT, LMMAXC, AUG_DES%LMMAX(NT))

            LMBASE = LMMAXC+LMBASE
            NAUG   = NAUG+AUG_DES%LMMAX(NT)

         ENDDO ion
# 1865

100      NIS = NIS+WDES%NITYP(NT)
      ENDDO typ

      ENDDO spinor

      RETURN
    END SUBROUTINE DEPSUM_TWO_BANDS_RHOLM

!****************** SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_TRACE ************
!
! this subroutine calculates the total charge from two bands
! it differs from the previous routine in the way spinors are handled
! this version returns only the spin-trace (total charge)
!
!***********************************************************************

    SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_TRACE( CPROJ1, CPROJ2, WDES, AUG_DES, & 
         TRANS_MATRIX, CRHOLM, FACT, LOVERL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1) WDES
      TYPE (wavedes1) AUG_DES
      LOGICAL LOVERL
      REAL(q) :: TRANS_MATRIX(:, :, :, :)
      REAL(q)  :: CPROJ1(:)
      REAL(q)  :: CPROJ2(:)
      REAL(q)  :: CTMP
      REAL(q)  :: CRHOLM(:)
      REAL(q) :: FACT              ! arbitrary multiplication factor
! local
      INTEGER :: ISPINOR, LMBASE, NT, LMMAXC, NIS, NI, L, LP, LM, NAUG
 
      IF (.NOT.LOVERL) RETURN

      CRHOLM=0

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)+1
      NAUG=   1

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 100

         ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

            CALL RHOLM_KERNEL(CPROJ1(LMBASE), CPROJ2(LMBASE), &
                 TRANS_MATRIX(1,1,1,NT), SIZE(TRANS_MATRIX,1), &
                 CRHOLM(NAUG), FACT, LMMAXC, AUG_DES%LMMAX(NT))

            LMBASE = LMMAXC+LMBASE
            NAUG   = NAUG+AUG_DES%LMMAX(NT)
            
         ENDDO ion
# 1931

100      NIS = NIS+WDES%NITYP(NT)
      ENDDO typ

      ENDDO spinor

      RETURN
    END SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_TRACE

    SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_TRACE_NOAE( CPROJ1, CPROJ2, WDES, AUG_DES, & 
         TRANS_MATRIX, CRHOLM, FACT, LOVERL, LMAX)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1) WDES
      TYPE (wavedes1) AUG_DES
      LOGICAL LOVERL
      REAL(q) :: TRANS_MATRIX(:, :, :, :)
      REAL(q)  :: CPROJ1(:)
      REAL(q)  :: CPROJ2(:)
      REAL(q)  :: CTMP
      REAL(q)  :: CRHOLM(:)
      REAL(q) :: FACT              ! arbitrary multiplication factor
      INTEGER :: LMAX(:)
! local
      INTEGER :: ISPINOR, LMBASE, NT, LMMAXC, NIS, NI, L, LP, LM, NAUG
      INTEGER LMMAX
      
      IF (.NOT.LOVERL) RETURN

      CRHOLM=0

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)+1
      NAUG=   1

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 100
         LMMAX=(LMAX(NT)+1)*(LMAX(NT)+1)
         IF (LMMAX > AUG_DES%LMMAX(NT)) THEN
            WRITE(*,*) 'internal error in DEPSUM_TWO_BANDS_RHOLM_TRACE_NOAE: LMMAX too large', LMMAX, AUG_DES%LMMAX(NT)
            CALL M_exit(); stop
         ENDIF

         ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

            CALL RHOLM_KERNEL(CPROJ1(LMBASE), CPROJ2(LMBASE), &
                 TRANS_MATRIX(1,1,1,NT), SIZE(TRANS_MATRIX,1), &
                 CRHOLM(NAUG), FACT, LMMAXC, LMMAX)

            LMBASE = LMMAXC+LMBASE
            NAUG   = NAUG+AUG_DES%LMMAX(NT)
            
         ENDDO ion

100      NIS = NIS+WDES%NITYP(NT)
      ENDDO typ

      ENDDO spinor

      RETURN
    END SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_TRACE_NOAE

!************************ SUBROUTINE DEPSUM_VECTOR  ******************
!
! this subroutine calculates  the occupancy matrix for 1._q ion
! as opposed to the previous routine the input array is a vector
!
!   G(r', N alpha) G*+(r', N alpha') T_alpha alpha'^NLM
!
! performance is much improved by blocking over the first index
!
!***********************************************************************

    SUBROUTINE DEPSUM_VECTOR( GO_R_PROJ, GU_R_PROJ, GWORK, NP, TRANS_MATRIX, LMMAXC, LMBASE, NLM_LMMAX )
      USE prec
      USE dfast

      REAL(q) :: GO_R_PROJ(:,:)
      REAL(q) :: GU_R_PROJ(:,:)
      REAL(q) :: GWORK(:,:)
      INTEGER :: NP        ! number of grid points
      REAL(q) :: TRANS_MATRIX(:, :, :)
      INTEGER :: LMMAXC    ! maximum number of orbital channels
      INTEGER :: LMBASE    ! base index
      INTEGER :: NLM_LMMAX ! maximum number of charge channels
! local
      INTEGER :: NR, NRP, L, LP, LM
      REAL(q) :: GO_GU_PROJ(NBLK)

! block over first index
      DO NR=1, NP, NBLK
         NRP=MIN(NR+NBLK-1, NP)
! loop over PAW alpha and alpha'
         DO L=1,LMMAXC
            DO LP=1,LMMAXC
               GO_GU_PROJ(1:NRP-NR+1)=GO_R_PROJ(NR:NRP, L+LMBASE)*(GU_R_PROJ(NR:NRP, LP+LMBASE))
! loop over charge index NLM
               DO LM=1,NLM_LMMAX
! index build from NI and LM
                  GWORK(NR:NRP,LM)=GWORK(NR:NRP,LM)+ &
                       GO_GU_PROJ(1:NRP-NR+1)*TRANS_MATRIX(LP,L,LM)

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
    END SUBROUTINE DEPSUM_VECTOR

!****************** DEPSUM_TWO_BANDS_RHOLM_NO_CONJG ********************
!
! this subroutine calculates the total charge from two bands
! it differs from the previous routine by the fact that the
! second wavefunction is not conjugated
!
!***********************************************************************

    SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_NO_CONJG( CPROJ1, CPROJ2, WDES, AUG_DES, & 
         TRANS_MATRIX, CRHOLM, FACT, LOVERL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1) WDES
      TYPE (wavedes1) AUG_DES
      LOGICAL LOVERL
      REAL(q) :: TRANS_MATRIX(:, :, :, :)
      REAL(q)  :: CPROJ1(:)
      REAL(q)  :: CPROJ2(:)
      REAL(q)  :: CTMP
      REAL(q)  :: CRHOLM(:)
      REAL(q) :: FACT              ! arbitrary multiplication factor
! local
      INTEGER :: ISPINOR, LMBASE, NT, LMMAXC, NIS, NI, L, LP, LM, NAUG
      
      IF (.NOT.LOVERL) RETURN

      CRHOLM=0

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)+1
      NAUG=   1

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 100

         ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

            CALL RHOLM_KERNEL_NO_CONJG(CPROJ1(LMBASE), CPROJ2(LMBASE), &
                 TRANS_MATRIX(1,1,1,NT), SIZE(TRANS_MATRIX,1), &
                 CRHOLM(NAUG), FACT, LMMAXC, AUG_DES%LMMAX(NT))

            LMBASE = LMMAXC+LMBASE
            NAUG   = NAUG+AUG_DES%LMMAX(NT)
            
         ENDDO ion

100      NIS = NIS+WDES%NITYP(NT)
      ENDDO typ

      ENDDO spinor

      RETURN
    END SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_NO_CONJG


!****************** SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_FULL *************
!
! this subroutine calculates the total charge from two bands
! it differs from the previous routine in the way spinors are handled
! this version is a full spinor version i.e. four components
! are returned (up,up), (up, down), (down, up) and (down,down)
!
!***********************************************************************

    SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_FULL(CPROJ1, CPROJ2, WDES, & 
         AUG_DES, TRANS_MATRIX, CRHOLM, LOVERL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1) WDES
      TYPE (wavedes1) AUG_DES
      LOGICAL LOVERL
      REAL(q) :: TRANS_MATRIX(:, :, :, :)
      REAL(q)  :: CPROJ1(:)
      REAL(q)  :: CPROJ2(:)
      REAL(q)  :: CTMP
      REAL(q)  :: CRHOLM(:)
! local
      INTEGER :: ISPINOR, ISPINOR_, LMBASE, LMBASE_, NT, LMMAXC, NIS, NI, L, LP, LM, NAUG
      INTEGER LMMAX
      
      IF (.NOT.LOVERL) RETURN

      CRHOLM=0

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)+1
      LMBASE_=ISPINOR_*(WDES%NPRO/2)+1

      NAUG= (ISPINOR_+2*ISPINOR)*AUG_DES%NPRO+1

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 100

         ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1
            CALL RHOLM_KERNEL(CPROJ1(LMBASE), CPROJ2(LMBASE), &
                 TRANS_MATRIX(1,1,1,NT), SIZE(TRANS_MATRIX,1), &
                 CRHOLM( NAUG ), 1.0_q, LMMAXC, AUG_DES%LMMAX(NT))
            LMBASE = LMMAXC+LMBASE
            LMBASE_= LMMAXC+LMBASE_
            NAUG = NAUG+AUG_DES%LMMAX(NT)

         ENDDO ion
# 2166

100      NIS = NIS+WDES%NITYP(NT)
      ENDDO typ

      ENDDO
      ENDDO spinor

      RETURN
    END SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_FULL

!
! this version is very similar but does not include the additional 1._q-center
! augmentation charges
! it is required in the the local Hartree Fock routines (pwlhf)
!

    SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_FULL_NOAE(CPROJ1, CPROJ2, WDES, & 
         AUG_DES, TRANS_MATRIX, CRHOLM, LOVERL, LMAX)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1) WDES
      TYPE (wavedes1) AUG_DES
      LOGICAL LOVERL
      REAL(q) :: TRANS_MATRIX(:, :, :, :)
      REAL(q)  :: CPROJ1(:)
      REAL(q)  :: CPROJ2(:)
      REAL(q)  :: CTMP
      REAL(q)  :: CRHOLM(:)
      INTEGER :: LMAX(:)
! local
      INTEGER :: ISPINOR, ISPINOR_, LMBASE, LMBASE_, NT, LMMAXC, NIS, NI, L, LP, LM, NAUG
      INTEGER LMMAX
      
      IF (.NOT.LOVERL) RETURN

      CRHOLM=0

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)+1
      LMBASE_=ISPINOR_*(WDES%NPRO/2)+1

      NAUG= (ISPINOR_+2*ISPINOR)*AUG_DES%NPRO+1

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 100
         LMMAX=(LMAX(NT)+1)*(LMAX(NT)+1)
         IF (LMMAX > AUG_DES%LMMAX(NT)) THEN
            WRITE(*,*) 'internal error in DEPSUM_TWO_BANDS_RHOLM_FULL_NOAE: LMMAX too large', LMMAX, AUG_DES%LMMAX(NT)
            CALL M_exit(); stop
         ENDIF

         ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1
            CALL RHOLM_KERNEL(CPROJ1(LMBASE), CPROJ2(LMBASE), &
                 TRANS_MATRIX(1,1,1,NT), SIZE(TRANS_MATRIX,1), &
                 CRHOLM( NAUG ), 1.0_q, LMMAXC, LMMAX)  ! AUG_DES%LMMAX(NT)
            LMBASE = LMMAXC+LMBASE
            LMBASE_= LMMAXC+LMBASE_
            NAUG = NAUG+AUG_DES%LMMAX(NT)

         ENDDO ion

100      NIS = NIS+WDES%NITYP(NT)
      ENDDO typ

      ENDDO
      ENDDO spinor

      RETURN
    END SUBROUTINE DEPSUM_TWO_BANDS_RHOLM_FULL_NOAE


!************************ SUBROUTINE DEPSUM_TWO_BANDS ******************
!
! this subroutine calculates  the occupancy matrix
! for two wavefunctions, for which the projected character is supplied
!
!  rho_\alpha,\beta :=   <phi_1 | p_beta> <p_alpha | phi_2>
!  CRHODE(LP,L,NI,ISPINOR_+2*ISPINOR)=
!             CPROJ1(L+LMBASE,N,NK)*(CPROJ2(LP+LMBASE_,N,NK))
!
! mind that alpha and beta are interchanged as everywhere in VASP
!
!***********************************************************************

    SUBROUTINE DEPSUM_TWO_BANDS( CPROJ1, CPROJ2, WDES, ISP, CRHODE, WEIGHT, LOVERL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)     WDES
      INTEGER ISP
      LOGICAL LOVERL
      REAL(q)  :: CPROJ1(WDES%NPRO)
      REAL(q)  :: CPROJ2(WDES%NPRO)
      REAL(q)  :: CRHODE(:,:,:,:)
      REAL(q)  :: WEIGHT
! local
      INTEGER :: ISPINOR, ISPINOR_, LMBASE, LMBASE_, NT, LMMAXC, NIS, NI, L, LP
      
      IF (.NOT.LOVERL) RETURN

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)
      LMBASE_=ISPINOR_*(WDES%NPRO/2)

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 100

         ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

!DIR$ IVDEP
!OCL NOVREC
            DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
               DO LP=1,LMMAXC
                  CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)= CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)+  &
                       CPROJ1(L+LMBASE)*(CPROJ2(LP+LMBASE_))*WEIGHT
               ENDDO
            ENDDO
      
            LMBASE = LMMAXC+LMBASE
            LMBASE_= LMMAXC+LMBASE_

         ENDDO ion
100      NIS = NIS+WDES%NITYP(NT)
      ENDDO typ

      ENDDO
      ENDDO spinor

      RETURN
    END SUBROUTINE DEPSUM_TWO_BANDS


!****************** DEPSUM_TWO_BANDS_ONE_CTR_TR ************************
!
! this subroutine calculates the total charge from two bands
! it also accumulates the 1._q center occupancies up to a maximum
! L quantum number
!
!***********************************************************************

    SUBROUTINE DEPSUM_TWO_BANDS_ONE_CTR_TR( CPROJ1, CPROJ2, WDES, AUG_DES, & 
         TRANS_MATRIX, CRHOLM, H, CRHO_ONE_CENTER, FACT, LOVERL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1) WDES
      TYPE (wavedes1) AUG_DES
      LOGICAL LOVERL
      REAL(q) :: TRANS_MATRIX(:, :, :, :)
      REAL(q)  :: CPROJ1(:)
      REAL(q)  :: CPROJ2(:)
      REAL(q)  :: CTMP
      REAL(q)  :: CRHOLM(:)
      TYPE (one_center_handle) :: H
      REAL(q)  :: CRHO_ONE_CENTER(H%TOTAL_ENTRIES)
      REAL(q) :: FACT              ! arbitrary multiplication factor
! local
      INTEGER :: ISPINOR, LMBASE, NT, LMMAXC, NIS, NI, L, LP, LM, NAUG
      INTEGER :: NRHO, INDEX_INTO_CG, INDEX_INTO_CG_TYPE
      
      IF (.NOT.LOVERL) RETURN

      CRHOLM=0
      CRHO_ONE_CENTER=0

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)+1
      NAUG=   1
      NRHO=   1
      INDEX_INTO_CG_TYPE=1

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 100
         
         ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1
            INDEX_INTO_CG=INDEX_INTO_CG_TYPE
            CALL RHOLM_ONE_CENTER_KERNEL( CPROJ1(LMBASE), CPROJ2(LMBASE), &
                 TRANS_MATRIX(1,1,1,NT), SIZE(TRANS_MATRIX,1), &
                 CRHOLM(NAUG), H%ENTRIES_FOR_TYPE(NT), CRHO_ONE_CENTER(NRHO),  &
                 INDEX_INTO_CG, H%NR_OF_CG(1,NT), H%CG(1), H%INDEX(1), FACT, LMMAXC, AUG_DES%LMMAX(NT))

            LMBASE = LMMAXC+LMBASE
            NAUG   = NAUG+AUG_DES%LMMAX(NT)
            NRHO   = NRHO+H%ENTRIES_FOR_TYPE(NT)
            
         ENDDO ion
         INDEX_INTO_CG_TYPE=INDEX_INTO_CG

100      NIS = NIS+WDES%NITYP(NT)
      ENDDO typ

      ENDDO spinor

      RETURN
    END SUBROUTINE DEPSUM_TWO_BANDS_ONE_CTR_TR


!************************* SUBROUTINE CALC_DLLMM_TRANS ****************
!
!  transform D(LM) to D(lm,l'm'), where D(LM) is defined as
!
!  D(LM) = \int V Y(L,M) Q(L)(r) d^3 r
!
!  and \int Q(L)(r) r^2 dr is integrating to 1. Result is added to CDIJ
!
!  this routine is similar to the routine CALC_DLLMM in the paw.F routine
!  but tuned for efficiency
!  and uses the stored transformation matrix
!  D(lm,l'm') = sum_LM TRANS_MATRIX(lm,l'm',LM,type) D(LM)
!
!**********************************************************************

    SUBROUTINE CALC_DLLMM_TRANS(WDES, AUG_DES, TRANS_MATRIX, CDIJ, CDLM)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes) WDES
      TYPE (wavedes1) AUG_DES
      REAL(q) :: TRANS_MATRIX(:,:,:,:) 
      REAL(q) CDIJ(:,:,:,:)
      REAL(q) CDLM(:)
      
! local
      INTEGER ISPINOR, NIS, NT, LMMAXC, NI, LM, L, LP, NPRO, NAUG

      CDIJ=0

      spinor: DO ISPINOR=0,WDES%NRSPINORS-1
         NIS =1
         NAUG=   ISPINOR *AUG_DES%NPRO+1

         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC/=0) THEN

               DO NI=NIS,WDES%NITYP(NT)+NIS-1
                  CALL DLLMM_KERNEL(CDIJ(1,1,NI,ISPINOR+1),SIZE(CDIJ,1), &
                       TRANS_MATRIX(1,1,1,NT), SIZE(TRANS_MATRIX,1), &
                       CDLM( NAUG), LMMAXC, AUG_DES%LMMAX(NT))
                  NAUG   = NAUG+AUG_DES%LMMAX(NT)
               ENDDO
# 2430

            ENDIF
            NIS = NIS+WDES%NITYP(NT)
         ENDDO
      ENDDO spinor

      RETURN
    END SUBROUTINE CALC_DLLMM_TRANS
!
! this version is very similar but does not include the additional 1._q-center
! augmentation charges
! it is required in the the local Hartree Fock routines (pwlhf)
!

    SUBROUTINE CALC_DLLMM_TRANS_NOAE(WDES, AUG_DES, TRANS_MATRIX, CDIJ, CDLM, LMAX)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes) WDES
      TYPE (wavedes1) AUG_DES
      REAL(q) :: TRANS_MATRIX(:,:,:,:) 
      REAL(q) CDIJ(:,:,:,:)
      REAL(q) CDLM(:)
      INTEGER :: LMAX(:)
      
! local
      INTEGER ISPINOR, NIS, NT, LMMAXC, NI, LM, L, LP, NPRO, NAUG
      INTEGER LMMAX

      CDIJ=0

      spinor: DO ISPINOR=0,WDES%NRSPINORS-1
         NIS =1
         NAUG=   ISPINOR *AUG_DES%NPRO+1

         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC/=0) THEN
               LMMAX=(LMAX(NT)+1)*(LMAX(NT)+1)
               DO NI=NIS,WDES%NITYP(NT)+NIS-1
                  CALL DLLMM_KERNEL(CDIJ(1,1,NI,ISPINOR+1),SIZE(CDIJ,1), &
                       TRANS_MATRIX(1,1,1,NT), SIZE(TRANS_MATRIX,1), &
                       CDLM( NAUG), LMMAXC, LMMAX)
                  NAUG   = NAUG+AUG_DES%LMMAX(NT)
               ENDDO
            ENDIF
            NIS = NIS+WDES%NITYP(NT)
         ENDDO
      ENDDO spinor

      RETURN
    END SUBROUTINE CALC_DLLMM_TRANS_NOAE

!************************ SUBROUTINE CALC_DLLMM_VECTOR  ****************
!
! this subroutine calculates
!
! G'(r', N alpha)=sum_beta,NLM W(r', NLM ) G(r', N beta) T_alpha beta^NLM
! G'(r', N alpha)=sum_beta G(r', N beta) sum_NLM T_alpha beta^NLM W(r', NLM)
!
! for a set of real space grid point
! performance is much improved by blocking over the first index
!
!***********************************************************************

    SUBROUTINE CALC_DLLMM_VECTOR( GO_R_PROJ, GI_R_PROJ, GWORK, NP, TRANS_MATRIX, LMMAXC, LMBASE, NLM_LMMAX )
      USE prec
      USE dfast

      REAL(q) :: GO_R_PROJ(:,:) ! out:  G'(r', N alpha)
      REAL(q) :: GI_R_PROJ(:,:) ! in:   G(r', N beta)
      REAL(q) :: GWORK(:,:)     ! in:   W(r',NLM)
      INTEGER :: NP        ! number of grid points
      REAL(q) :: TRANS_MATRIX(:, :, :)
      INTEGER :: LMMAXC    ! maximum number of orbital channels
      INTEGER :: LMBASE    ! base index
      INTEGER :: NLM_LMMAX ! maximum number of charge channels
! local
      INTEGER :: NR, NRP, L, LP, LM
      REAL(q) :: PROJ(NBLK)

! block over first real grid point
      DO NR=1, NP, NBLK
         NRP=MIN(NR+NBLK-1, NP)
! loop over PAW alpha and alpha'
         DO L=1,LMMAXC
            DO LP=1,LMMAXC
! p(r') = sum_NLM T_alpha beta^NLM W(r', NLM)
               PROJ=0
               DO LM=1,NLM_LMMAX
                  PROJ(1:NRP-NR+1)=PROJ(1:NRP-NR+1)+GWORK(NR:NRP,LM)*TRANS_MATRIX(LP,L,LM)
               ENDDO
! sum_beta G(r', N beta) p(r')
               GO_R_PROJ(NR:NRP, LMBASE+LP)=GO_R_PROJ(NR:NRP, LMBASE+LP)+ &
                       GI_R_PROJ(NR:NRP,LMBASE+L)*PROJ(1:NRP-NR+1)
            ENDDO
         ENDDO
      ENDDO
      
    END SUBROUTINE CALC_DLLMM_VECTOR

END MODULE augfast

!***********************************************************************
!
!  apply the phase factor e^(i q R_i) to a linear CG transformed
!  occupancy matrix.
!  this is the same subroutine as above but without an
!  explicit calling interfact
!
!***********************************************************************


  SUBROUTINE APPLY_PHASE_ONE_CENTER_NOINT(WDES, H, CRHO, N, VKPT )
    USE prec
    USE constant
    USE augfast
    IMPLICIT NONE

    TYPE (wavedes) WDES
    TYPE (one_center_handle) :: H
    INTEGER :: N
    REAL(q) :: CRHO(H%TOTAL_ENTRIES,N)    ! 1._q center charge density for a number of vectors
    REAL(q) :: VKPT(3)
! local
    REAL(q),PARAMETER :: TINY=1E-4_q
    COMPLEX (q) CSHIFT
    INTEGER NIS, NRHO, NT, NI, NIP

    IF (ABS(VKPT(1))>TINY .OR. ABS(VKPT(2))>TINY .OR. ABS(VKPT(3))>TINY) THEN
       NRHO=1
       NIS=1
       DO NT=1,WDES%NTYP
       
          DO NI=NIS,WDES%NITYP(NT)+NIS-1
             NIP=NI_GLOBAL(NI,WDES%COMM_INB)
             
             CSHIFT=EXP(CITPI*(VKPT(1)*H%POSION(1,NIP)+ & 
                               VKPT(2)*H%POSION(2,NIP)+ &
                               VKPT(3)*H%POSION(3,NIP)))
             
             CRHO(NRHO:NRHO+H%ENTRIES_FOR_TYPE(NT)-1,:)=CRHO(NRHO:NRHO+H%ENTRIES_FOR_TYPE(NT)-1,:)*CSHIFT
             NRHO = NRHO+H%ENTRIES_FOR_TYPE(NT)
          ENDDO
          
          NIS = NIS+WDES%NITYP(NT)
       ENDDO
    ENDIF
  END SUBROUTINE APPLY_PHASE_ONE_CENTER_NOINT


!*******************************************************************
!
! calculate the squared 1._q center occupancy matrix
!
!*******************************************************************

  SUBROUTINE GEMM_ONE_CENTER( LREALSTORE, LREAL, WDES, H, CRHO, CRHOW, N, & 
       RESPONSEONE, RESPONSEONER)
    USE wave
    USE augfast
    IMPLICIT NONE

    LOGICAL :: LREALSTORE            ! response function is real values
    LOGICAL :: LREAL                 ! CRHO and CRHOW are real vectors
    TYPE (wavedes) WDES
    TYPE (one_center_handle) :: H
    REAL(q) :: CRHO(H%TOTAL_ENTRIES,N)  ! 1._q center charge density for a number of vectors
    REAL(q) :: CRHOW(H%TOTAL_ENTRIES,N) ! 1._q center charge density for a number of vectors times weight
    INTEGER :: N                     ! number of vectors in CRHO and CRHOW
    COMPLEX(q) :: RESPONSEONE(*)     ! complex valued response function
    REAL(q)    :: RESPONSEONER(*)    ! real valued response function
! local
    INTEGER :: NIS, NT, NI
    INTEGER :: NPOT, NRHO, NRHOW

    NIS =1
    NRHO=1
    NRHOW=1
    NPOT=1
    DO NT=1,WDES%NTYP
       
       DO NI=NIS,WDES%NITYP(NT)+NIS-1
          IF (LREALSTORE) THEN
! CRHO and CRHOW hold H%ENTRIES_FOR_TYPE(NT) real entries
! resulting matrix is REAL as well
             CALL DGEMM('N','T', H%ENTRIES_FOR_TYPE(NT), H%ENTRIES_FOR_TYPE(NT), N, 1.0_q, &
                  CRHOW(1,1),  H%TOTAL_ENTRIES, CRHO(1,1),  H%TOTAL_ENTRIES, &
                  1.0_q, RESPONSEONER(NPOT), H%ENTRIES_FOR_TYPE(NT))
! eat up entries from NRHOW
             NRHOW= NRHOW+H%ENTRIES_FOR_TYPE(NT)
          ELSE IF (LREAL) THEN
! CRHOW holds H%ENTRIES_FOR_TYPE(NT) complex entries, CRHO real entries
! resulting matrix is has H%ENTRIES_FOR_TYPE(NT) x H%ENTRIES_FOR_TYPE(NT)
! complex entries
             CALL DGEMM('N','T', 2*H%ENTRIES_FOR_TYPE(NT), H%ENTRIES_FOR_TYPE(NT), N, 1.0_q, &
                  CRHOW(1,1),  2*H%TOTAL_ENTRIES, CRHO(1,1),  H%TOTAL_ENTRIES, &
                  1.0_q, RESPONSEONE(NPOT), 2*H%ENTRIES_FOR_TYPE(NT))
! eat up 2 real entries from CRHOW (1._q complex)
             NRHOW= NRHOW+2*H%ENTRIES_FOR_TYPE(NT)
          ELSE
! CRHOW and CRHO hold H%ENTRIES_FOR_TYPE(NT) complex entries
             CALL ZGEMM('N','C', H%ENTRIES_FOR_TYPE(NT), H%ENTRIES_FOR_TYPE(NT), N, (1.0_q,0.0_q), &
                  CRHOW(1,1),  H%TOTAL_ENTRIES, CRHO(1,1),  H%TOTAL_ENTRIES, &
                  (1.0_q,0.0_q), RESPONSEONE(NPOT), H%ENTRIES_FOR_TYPE(NT))
! eat up NRHOW complex entries from NRHOW
             NRHOW= NRHOW+H%ENTRIES_FOR_TYPE(NT)
          ENDIF

          NRHO= NRHO+H%ENTRIES_FOR_TYPE(NT)
          NPOT= NPOT+H%ENTRIES_FOR_TYPE(NT)*H%ENTRIES_FOR_TYPE(NT)
       ENDDO

    NIS = NIS+WDES%NITYP(NT)
    ENDDO

  END SUBROUTINE GEMM_ONE_CENTER


!**********************************************************************
!
! low level F77 subroutine to calculate RHOLM
!
!**********************************************************************

!
! version that maps onto DGEMV
!
  SUBROUTINE RHOLM_KERNEL( CPROJ1, CPROJ2, TRANS_MATRIX, LD, CRHOLM, FACT, LMMAXC, LMMAX)
    USE prec
    USE wave
    IMPLICIT NONE

    INTEGER :: LD
    INTEGER LMMAX
    REAL(q) :: TRANS_MATRIX(LD* LD, LMMAX)
    INTEGER LMMAXC
    REAL(q)  :: CPROJ1(LMMAXC)
    REAL(q)  :: CPROJ2(LMMAXC)
    REAL(q) :: CRHOLM(LMMAX)
    REAL(q) :: FACT              ! arbitrary multiplication factor
! local

    REAL(q)  :: CTMP(LD, LD)
# 2677


    CALL RHOLM_KERNEL_AUX( CPROJ1, CPROJ2, LD, CTMP, FACT, LMMAXC)

! CRHOLM(LM)=CRHOLM(LM)+TRANS_MATRIX(LP,L,LM)* CTMP(LP,L)

    CALL DGEMV( 'T' , LD*LD, LMMAX, 1._q , TRANS_MATRIX(1,1) , &
                  LD*LD, CTMP(1,1), 1 , 1._q ,  CRHOLM(1), 1)
# 2692

!    WRITE(*,'(8F14.7)') CRHOLM(1:8)

  END SUBROUTINE RHOLM_KERNEL

  SUBROUTINE RHOLM_KERNEL_AUX( CPROJ1, CPROJ2, LD, CTMP, FACT, LMMAXC)
    USE prec
    USE wave
    IMPLICIT NONE

    INTEGER LMMAXC
    INTEGER LD
    REAL(q)  :: CPROJ1(LMMAXC)
    REAL(q)  :: CPROJ2(LMMAXC)
    REAL(q) :: FACT              ! arbitrary multiplication factor
! local
    REAL(q)  :: CTMP(LD, LD)
    INTEGER L, LP, LM

    CTMP=0
!DIR$ IVDEP
!OCL NOVREC
    DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
       DO LP=1,LMMAXC
          CTMP(LP,L)=CPROJ1(L)*(CPROJ2(LP))*FACT
       ENDDO
    ENDDO
  END SUBROUTINE
# 2757

!
! version that maps onto DGEMM
!
  SUBROUTINE RHOLM_KERNEL_DGEMM( CPROJ1, CPROJ2, TRANS_MATRIX, LD, CRHOLM, FACT, LMMAXC, LMMAX, NIONS)
    USE prec
    USE wave
    IMPLICIT NONE

    INTEGER :: LD
    INTEGER :: LMMAX, LMMAXC
    INTEGER :: NIONS
    REAL(q) :: TRANS_MATRIX(LD* LD, LMMAX)
    REAL(q)    :: CPROJ1(LMMAXC,NIONS)
    REAL(q)    :: CPROJ2(LMMAXC,NIONS)

!   REAL(q)    :: CRHOLM(LMMAX,NIONS)
    REAL(q) :: CRHOLM( 2* LMMAX, NIONS)

    REAL(q) :: FACT              ! arbitrary multiplication factor
! local
    INTEGER :: L, LP, NI

    REAL(q)  :: CTMP(LD, LD, NIONS)
# 2784

    REAL(q)  :: CM

    CTMP=0
    DO NI=1,NIONS
!DIR$ IVDEP
!OCL NOVREC
       DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
          DO LP=1,LMMAXC
             CM=CPROJ1(L,NI)*(CPROJ2(LP,NI))*FACT
             CTMP(LP, L, NI)=REAL(CM, KIND=q)
# 2799

          ENDDO
       ENDDO
    ENDDO

! CRHOLM(LM)=CRHOLM(LM)+TRANS_MATRIX(LP,L,LM)* CTMP(LP,L)

    CALL DGEMM( 'T' , 'N', LMMAX, NIONS, LD*LD, 1._q , TRANS_MATRIX(1,1) , &
                  LD*LD, CTMP(1,1,1), LD*LD , 1._q ,  CRHOLM(1,1), LMMAX)
# 2816

!    WRITE(*,'(8F14.7)') CRHOLM(1:4)
  END SUBROUTINE RHOLM_KERNEL_DGEMM


  SUBROUTINE RHOLM_KERNEL_NO_CONJG( CPROJ1, CPROJ2, TRANS_MATRIX, LD, CRHOLM, FACT, LMMAXC, LMMAX)
    USE prec
    USE wave
    IMPLICIT NONE

    INTEGER :: LD
    INTEGER LMMAX
    REAL(q) :: TRANS_MATRIX(LD, LD, LMMAX)
    INTEGER LMMAXC
    REAL(q)  :: CPROJ1(LMMAXC)
    REAL(q)  :: CPROJ2(LMMAXC)
    REAL(q)  :: CRHOLM(LMMAX)
    REAL(q) :: FACT              ! arbitrary multiplication factor
! local
    REAL(q)  :: CTMP
    INTEGER L, LP, LM

    DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
       DO LP=1,LMMAXC
          CTMP=CPROJ1(L)*CPROJ2(LP)*FACT
!DIR$ IVDEP
!OCL NOVREC
          DO LM=1,LMMAX
             CRHOLM(LM)=CRHOLM(LM)+CTMP*TRANS_MATRIX(LP,L,LM)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE RHOLM_KERNEL_NO_CONJG


!**********************************************************************
!
! low level F77 subroutine to calculate RHOLM
! and contracted 1._q center occupancy matrix
!
!**********************************************************************


  SUBROUTINE RHOLM_ONE_CENTER_KERNEL( CPROJ1, CPROJ2, TRANS_MATRIX, LD, & 
       CRHOLM, ENTRIES_FOR_TYPE, CRHO_ONE_CENTER,  & 
       INDEX_INTO_CG, NR_OF_CG, CG, INDEX, FACT, LMMAXC, LMMAX)
    USE prec
    USE wave
    IMPLICIT NONE

    INTEGER :: LD
    INTEGER LMMAX
    REAL(q) :: TRANS_MATRIX(LD, LD, LMMAX)
    INTEGER LMMAXC
    REAL(q)  :: CPROJ1(LMMAXC)
    REAL(q)  :: CPROJ2(LMMAXC)
    REAL(q)  :: CRHOLM(LMMAX)
    INTEGER :: ENTRIES_FOR_TYPE
    REAL(q)    :: CRHO_ONE_CENTER(ENTRIES_FOR_TYPE)
    INTEGER :: INDEX_INTO_CG
    INTEGER :: NR_OF_CG(ENTRIES_FOR_TYPE)      ! number of Clebsch Gordan coeff for each entry
    REAL(q) :: CG(*)                           ! equivalent to CG array in H
    INTEGER :: INDEX(*)                        ! equivalent to INDEX array in H
    REAL(q) :: FACT                            ! arbitrary multiplication factor
! local
    REAL(q)  :: CTMP(LMMAXC* LMMAXC), CTMP_
    INTEGER L, LP, LM, I, J

    DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
       DO LP=1,LMMAXC
          CTMP_=CPROJ1(L)*(CPROJ2(LP))*FACT
          CTMP(LP+(L-1)*LMMAXC)=CTMP_
!DIR$ IVDEP
!OCL NOVREC
          DO LM=1,LMMAX
             CRHOLM(LM)=CRHOLM(LM)+CTMP_*TRANS_MATRIX(LP,L,LM)
          ENDDO
       ENDDO
    ENDDO

    DO I=1,ENTRIES_FOR_TYPE
       DO J=1,NR_OF_CG(I)
          CRHO_ONE_CENTER(I)=CRHO_ONE_CENTER(I)+CG(INDEX_INTO_CG)*CTMP(INDEX(INDEX_INTO_CG))
          INDEX_INTO_CG=INDEX_INTO_CG+1
       ENDDO
    ENDDO

  END SUBROUTINE RHOLM_ONE_CENTER_KERNEL



!**********************************************************************
!
! low level F77 subroutine to calculate CDLM
!
!**********************************************************************

!
! version that maps onto DGEMV
!
  SUBROUTINE DLLMM_KERNEL( CDIJ, LDCDIJ, TRANS_MATRIX, LD, CDLM, LMMAXC, LMMAX)
    USE prec
    USE wave
    IMPLICIT NONE

    INTEGER :: LD
    INTEGER LMMAX
    REAL(q) :: TRANS_MATRIX(LD*LD, LMMAX)
    INTEGER :: LDCDIJ

    REAL(q) CDLM( LMMAX)
    REAL(q) CDIJ( LDCDIJ*LDCDIJ)
# 2936

    INTEGER LMMAXC
! local
    INTEGER L, LP, LM

! CDIJ(L)= CDIJ(L)+TRANS_MATRIX(L,LM)*CDLM(LM),  LM=1,LMMAX, L=1,LD*LD

    CALL DGEMV( 'N' , LD*LD, LMMAX, 1._q , TRANS_MATRIX(1,1) , &
                  LD*LD, CDLM(1), 1 , 1._q ,  CDIJ(1), 1)
# 2952


  END SUBROUTINE DLLMM_KERNEL
# 2987

!
! version that maps onto DGEMM
!
  SUBROUTINE DLLMM_KERNEL_DGEMM( CDIJ, LDCDIJ, NIONS, TRANS_MATRIX, LD, CDLM, LMMAX)
    USE prec
    USE wave
    IMPLICIT NONE

    INTEGER :: LD
    INTEGER :: LMMAX
    REAL(q) :: TRANS_MATRIX(LD*LD, LMMAX)
    INTEGER :: LDCDIJ
    INTEGER :: NIONS

    REAL(q) CDLM( LMMAX,NIONS)
    REAL(q) CDIJ( LDCDIJ*LDCDIJ,NIONS)
# 3012

! local
    INTEGER L, LP, LM


    CALL DGEMM( 'N', 'N', LD*LD, NIONS, LMMAX, 1._q , TRANS_MATRIX(1,1), LD*LD, &
   &           CDLM(1,1), LMMAX, 1._q , CDIJ(1,1), LDCDIJ*LDCDIJ)
# 3035

  END SUBROUTINE DLLMM_KERNEL_DGEMM
