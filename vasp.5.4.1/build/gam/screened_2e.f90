# 1 "screened_2e.F"
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

# 2 "screened_2e.F" 2 

!*********************************************************************
!
! This module implements a set a modules to handle
! the screened two electron integrals
! the central quantity is the structure screened_2e_handle
! which allows to map a set of (k1,q) points onto (k1',k2)
! where k1' lies in the IRZ
! this is in fact much more complicated than 1._q would expect
! the routine also properly symmetrizes the two electron integrals
! handling in particular degeneracies between eigenvalue pairs
!
!*********************************************************************

MODULE screened_2e
  USE prec
  IMPLICIT NONE
  LOGICAL :: LSPECTRALGW=.FALSE.

  TYPE screened_2e_handle
     INTEGER :: NUMBER_OF_NQ
     INTEGER :: NUMBER_OF_NQ_FULL
     INTEGER, POINTER :: NQ(:)
     INTEGER, POINTER :: K1_IN_IRZ(:,:)
     INTEGER, POINTER :: K2_IN_IRZ(:,:)
     LOGICAL, POINTER :: REQUIRED(:,:)
     REAL(q), POINTER :: WTKPT(:,:)
     INTEGER, POINTER :: NUMBER_OF_REDUNDANT(:,:)
     INTEGER, POINTER :: K2_STORE_INDEX(:,:)
  END TYPE screened_2e_handle


  TYPE screened_2e_spline
     INTEGER :: NOMEGA
     REAL(q), POINTER :: REAL_PART(:,:)
     REAL(q), POINTER :: IMAG_PART(:,:)
  END TYPE screened_2e_spline

! the following variables determine the spacing and the
! maximum and minimum energies for the frequency dependent selfenergy
! see routines CALC_SELFENERGY
  INTEGER :: NOMEGARES=201
  REAL(q) :: OMEGAMINR=-100,OMEGAMAXR=100

CONTAINS
!*********************************************************************
!
! Bloch integrals of the type
!
! \sum_k1 \sum_k2 u_*k1(r) u*_k2(r) F(r,r') u_k1(r') u*_k2(r')
!
! can be performed
! ) either by restricting k1 to the IRZ and looping over all k2
! ) or by restricting q= k1-k2 to the IRZ and looping over all k1
!   (with k2=k1-q).
! If eigenvalues are required the first version is easier to use,
! but the SCREENED_TWO_ELECTRON_INTEGRAL and ADD_XI use the second
! variant
! the screened_2e_handle stores all required quantities to map loops
! of the first kind to the second kind
!
!   NUMBER_OF_NQ       number of q-points in IRZ
!   NUMBER_OF_NQ_FULL  total number of q-points in full BZ
!   K1_IN_IRZ(NQ,NK1)  stores the index NK1 mapped into the IRZ
!   K2_IN_IRZ(NQ,NK1)  stores the corresponding index NK2,
!                      when NK1 is mapped into the IRZ
!   REQUIRED(NQ,NK1)   determines whether this integral needs
!                      to be calculated, due to symmetry some
!                      integrals are redundant
!                      .T. integral needs to be calculated
!                      .F. integral does not need to be calculated
!   WTKPT(K1_IN_IRZ,K2_IN_IRZ)
!                      stores the weight of the two electron integral
!                      this takes care of the fact that some
!                      integrals are symmetry equivalent
!                      and properly weights the contributions
!   NUMBER_OF_REDUNDANT(K1_IN_IRZ,K2_IN_IRZ)
!                      counter for the number of in principle
!                      symmetry inequivalent two electron integrals
!                      usually 1 or 0
!   K2_STORE_INDEX(K1_IN_IRZ,K2_IN_IRZ)
!                      since not all K2_IN_IRZ indices might be required
!                      (only those with WTKPT(K1_IN_IRZ,K2_IN_IRZ)/=0
!                       are required)
!                      the compact index allows to use a more compact
!                      storage mode for screened to electron integrals
!                      it stores the actual storage position NK2_STORE
!                      for each K1, K2 pair
!
!*********************************************************************

  SUBROUTINE SETUP_IRZ_MAP(S2E, WDES, IU0, IU6)
    USE wave
    USE kpoints_change
    USE fock
    TYPE (wavedes)     WDES
    TYPE (screened_2e_handle) S2E
    INTEGER IU0, IU6
! local
! if LALLINT is set,
! all two electron integrals even if related by symmetry
! are calculated, however the way the loop is 1._q, it is yet impossible
! to calculate all symmetry equivalent (consider K1=0, K2 is limited to K1+Q
! where Q is in the IRZ
    LOGICAL,PARAMETER :: LALLINT=.FALSE.

    INTEGER NQ_INDEX_IN_FULL_GRID, NQ, NK1, NK2, K1_IN_IRZ, K2_IN_IRZ
    LOGICAL :: KPOINT_VISITED(KPOINTS_ORIG%NKPTS)
    REAL(q), PARAMETER :: TINY=1.E-6_q
    REAL(q) WEIGHT_SUM
    LOGICAL :: NQ_REQUIRED(KPOINTS_FULL_ORIG%NKPTS)
    INTEGER :: IERR
    KPOINT_VISITED=.FALSE.
    IF (WDES%NKPTS /= KPOINTS_FULL_ORIG%NKPTS) THEN
       WRITE(0,*)'internal error in SETUP_IRZ_MAP: WDES%NKPTS needs to be identical to KPOINTS_FULL_ORIG%NKPTS'
       WRITE(0,*)'  usually this requires to regenerate the k-point grid without symmetry'
       WRITE(0,*)'  WDES%NKPTS, KPOINTS_FULL_ORIG%NKPTS=', WDES%NKPTS, KPOINTS_FULL_ORIG%NKPTS
       CALL M_exit(); stop
    ENDIF

    CALL SET_KPOINT_MAP_INTO_BZ(KPOINTS_FULL)
    CALL SET_KPOINT_MAP_INTO_BZ(KPOINTS_FULL_ORIG)

    ALLOCATE(S2E%NQ(KPOINTS_ORIG%NKPTS), &
             S2E%K1_IN_IRZ(KPOINTS_ORIG%NKPTS, WDES%NKPTS), &
             S2E%K2_IN_IRZ(KPOINTS_ORIG%NKPTS, WDES%NKPTS), &
             S2E%REQUIRED(KPOINTS_ORIG%NKPTS, WDES%NKPTS), &
             S2E%WTKPT(KPOINTS_ORIG%NKPTS, WDES%NKPTS), &
             S2E%NUMBER_OF_REDUNDANT(KPOINTS_ORIG%NKPTS, WDES%NKPTS), &
             S2E%K2_STORE_INDEX(KPOINTS_ORIG%NKPTS, WDES%NKPTS))

    S2E%REQUIRED=.TRUE.
    S2E%WTKPT=0
    S2E%NUMBER_OF_REDUNDANT=0
    S2E%NUMBER_OF_NQ=0
    S2E%K2_STORE_INDEX=0

! first loop over all k-point pairs and determined all possible relative
! vectors, map them into the BZ and tag them in the array NQ_REQUIRED
    NQ_REQUIRED=.FALSE.
    DO NK1=1,1
       IF (KPOINTS_FULL_ORIG%WTKPT(NK1)==0) CYCLE
       DO NK2=1,KPOINTS_FULL_ORIG%NKPTS
          IF (KPOINTS_FULL_ORIG%WTKPT(NK2)==0) CYCLE

          NQ_INDEX_IN_FULL_GRID= & 
               KPOINT_IN_FULL_GRID(KPOINTS_FULL_ORIG%VKPT(:,NK1)-KPOINTS_FULL_ORIG%VKPT(:,NK2),KPOINTS_FULL_ORIG)
          NQ_REQUIRED(NQ_INDEX_IN_FULL_GRID)=.TRUE.
       ENDDO
    ENDDO

! loop over all k-points
    S2E%NUMBER_OF_NQ_FULL=0
    DO NQ_INDEX_IN_FULL_GRID=1,KPOINTS_FULL_ORIG%NKPTS
!
! test whether this k-point is required or has already been visited
!
       IF (.NOT. NQ_REQUIRED(NQ_INDEX_IN_FULL_GRID)) CYCLE
       IF (SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL_ORIG%VKPT(:,NQ_INDEX_IN_FULL_GRID))) CYCLE
       S2E%NUMBER_OF_NQ_FULL=S2E%NUMBER_OF_NQ_FULL+1
    ENDDO

! loop over all difference q-points
    qpoints: DO NQ_INDEX_IN_FULL_GRID=1,KPOINTS_FULL_ORIG%NKPTS
!
! test whether this k-point is required or has already been visited
!
       IF (.NOT. NQ_REQUIRED(NQ_INDEX_IN_FULL_GRID)) CYCLE
       IF (SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL_ORIG%VKPT(:,NQ_INDEX_IN_FULL_GRID))) CYCLE
       NQ=KPOINTS_FULL_ORIG%NEQUIV(NQ_INDEX_IN_FULL_GRID)  ! index in the IRZ
       IF (NQ > SIZE(KPOINT_VISITED)) THEN
          WRITE(0,*)'internal error in SETUP_IRZ_MAP: KPOINT_VISITED array too small'
          WRITE(0,*)' this usually means a point does not map into the IRZ',NQ
          CALL M_exit(); stop
       ENDIF
       IF (KPOINT_VISITED(NQ)) CYCLE
       KPOINT_VISITED(NQ)=.TRUE.

       S2E%NUMBER_OF_NQ=S2E%NUMBER_OF_NQ+1
       IF (S2E%NUMBER_OF_NQ>SIZE(S2E%NQ)) THEN
          WRITE(0,*)'internal error in SETUP_IRZ_MAP: NQ array too small'
          CALL M_exit(); stop
       ENDIF
       S2E%NQ(S2E%NUMBER_OF_NQ)=NQ
! loop over all k-points NK1
       DO NK1=1, WDES%NKPTS
! if the k-point has 0._q weight we do not need to include it
          IF (WDES%WTKPT(NK1)==0) THEN
             S2E%REQUIRED(NQ, NK1)=.FALSE.
             CYCLE
          ENDIF
! determine NK2 = NK1-NQ
          NK2=KPOINT_IN_FULL_GRID(WDES%VKPT(:,NK1)-KPOINTS_ORIG%VKPT(:,NQ),KPOINTS_FULL)
! now bring NK1 to IRZ, and rotate NK2 by the same operation
          K1_IN_IRZ=KPOINT_IN_FULL_GRID( &
               BRING_TO_IRZ( WDES%VKPT(:,NK1), WDES%VKPT(:,NK1), .TRUE.), &
               KPOINTS_FULL)
          K2_IN_IRZ=KPOINT_IN_FULL_GRID( &
               BRING_TO_IRZ( WDES%VKPT(:,NK1), KPOINTS_FULL%VKPT(:,NK2)), &
               KPOINTS_FULL)
          S2E%K1_IN_IRZ(NQ, NK1)=K1_IN_IRZ
          S2E%K2_IN_IRZ(NQ, NK1)=K2_IN_IRZ

! this particular two electron integral has already been calculated
! skip it
          IF (S2E%WTKPT(K1_IN_IRZ, K2_IN_IRZ)/=0) THEN
             S2E%REQUIRED(NQ, NK1)=.FALSE.
          ENDIF
! if LALLINT is set,
! all two electron integrals even if related by symmetry
! are calculated
          IF (LALLINT) THEN
             S2E%REQUIRED(NQ, NK1)=.TRUE.
          ENDIF

          IF ( KPOINTS_ORIG%WTKPT(NQ) /=0) THEN
! add to weight of q-point (difference vecor) of this k-point pair
             S2E%WTKPT(K1_IN_IRZ, K2_IN_IRZ)=S2E%WTKPT(K1_IN_IRZ, K2_IN_IRZ)+ & 
                  KPOINTS_ORIG%WTKPT(NQ)
          ELSE
! strange the weight was 0._q
! this can happen if auxiliary k-points with 0._q weights are used
             S2E%WTKPT(K1_IN_IRZ, K2_IN_IRZ)=S2E%WTKPT(K1_IN_IRZ, K2_IN_IRZ)+1._q/KPOINTS_FULL_ORIG%NKPTS_NON_ZERO
          ENDIF

          IF (S2E%REQUIRED(NQ, NK1)) THEN
             S2E%NUMBER_OF_REDUNDANT(K1_IN_IRZ, K2_IN_IRZ)= &
                  S2E%NUMBER_OF_REDUNDANT(K1_IN_IRZ, K2_IN_IRZ)+1
          ENDIF
       ENDDO

    ENDDO qpoints

    S2E%WTKPT=S2E%WTKPT/KPOINTS_FULL_ORIG%NKPTS_NON_ZERO
!
! the summed weights now correspond exactly to the original weight
! KPOINTS_ORIG%WTKPT(K1_IN_IRZ)==SUM(S2E%WTKPT(K1_IN_IRZ, :))
! this is the proper weight for total energy calculations
! we however want to calculate shifts of the eigenvalues
! divide by KPOINTS_ORIG%WTKPT(K1_IN_IRZ) such that SUM(S2E%WTKPT(K1_IN_IRZ, :))==1
!
    IERR=0
    DO K1_IN_IRZ=1, KPOINTS_ORIG%NKPTS
       IF (KPOINTS_ORIG%WTKPT(K1_IN_IRZ)/=0) THEN
          S2E%WTKPT(K1_IN_IRZ, :)=S2E%WTKPT(K1_IN_IRZ, :)/KPOINTS_ORIG%WTKPT(K1_IN_IRZ)*NKREDX*NKREDY*NKREDZ
          IF (ODDONLY)  S2E%WTKPT(K1_IN_IRZ, :)=S2E%WTKPT(K1_IN_IRZ, :)*2
          IF (EVENONLY) S2E%WTKPT(K1_IN_IRZ, :)=S2E%WTKPT(K1_IN_IRZ, :)*2

          WEIGHT_SUM=SUM(S2E%WTKPT(K1_IN_IRZ, :))

          IF (ABS(WEIGHT_SUM-1)>=TINY) THEN
             IERR=IERR+1
             S2E%WTKPT(K1_IN_IRZ, :)=S2E%WTKPT(K1_IN_IRZ, :)/WEIGHT_SUM
          ENDIF
       ELSE
          WEIGHT_SUM=SUM(S2E%WTKPT(K1_IN_IRZ, :))
          IF (WEIGHT_SUM/=0) THEN
             WRITE(0,*) 'internal error in VASP: SETUP_IRZ_MAP weight at k-point with zero WTKPT'
             CALL M_exit(); stop
          ENDIF
       ENDIF

! divide final weights by the number of integral evaluations
       DO NK2=1,WDES%NKPTS
          IF (S2E%NUMBER_OF_REDUNDANT(K1_IN_IRZ, NK2)/=0) THEN
             S2E%WTKPT(K1_IN_IRZ, NK2)=S2E%WTKPT(K1_IN_IRZ, NK2)/S2E%NUMBER_OF_REDUNDANT(K1_IN_IRZ, NK2)
          ENDIF
       ENDDO
    ENDDO

    IF (IERR>0) THEN
         CALL VTUTOR('E','GW kweight', &
     &               0.0_q,1,1,1,(1.0_q,0.0_q),1,.TRUE.,1,IU0,3)
         CALL VTUTOR('E','GW kweight', &
     &               0.0_q,1,1,1,(1.0_q,0.0_q),1,.TRUE.,1,IU6,3)
    ENDIF
!
!  now set up the compact index
!
    DO K1_IN_IRZ=1 , KPOINTS_ORIG%NKPTS
       NK2=0
       DO K2_IN_IRZ=1, WDES%NKPTS
          IF (S2E%NUMBER_OF_REDUNDANT(K1_IN_IRZ, K2_IN_IRZ)/=0) THEN
             NK2=NK2+1
             S2E%K2_STORE_INDEX(K1_IN_IRZ,K2_IN_IRZ)=NK2
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE SETUP_IRZ_MAP


  SUBROUTINE DEALLOCATE_IRZ_MAP(S2E)
   TYPE (screened_2e_handle) S2E

   DEALLOCATE(S2E%NQ, &
             S2E%K1_IN_IRZ, &
             S2E%K2_IN_IRZ, &
             S2E%REQUIRED, &
             S2E%WTKPT, &
             S2E%NUMBER_OF_REDUNDANT, &
             S2E%K2_STORE_INDEX)

  END SUBROUTINE
    
!*********************************************************************
!
! BRING_TO_IRZ determines to which k-point VKPT1 corresponds in
! the IRZ and applies the symmetry operation that brings VKPT1 into the IRZ
! to VKPT2
!
!*********************************************************************

  FUNCTION BRING_TO_IRZ( VKPT1, VKPT2, LTEST)
    USE wave
    USE kpoints_change
    REAL(q) :: VKPT1(3)
    REAL(q) :: VKPT2(3)
    REAL(q) :: BRING_TO_IRZ(3)
    LOGICAL, OPTIONAL :: LTEST
! local
    INTEGER NK1_FULL_ORIG   ! k-points index in the original k-point structure
    INTEGER NK1_IN_IRZ      ! k-points index in the irz
    INTEGER NK, NK_FOUND
    REAL(q) :: VT(3)
    REAL(q), PARAMETER :: TINY=1.E-6_q
    REAL(q) :: AMAT(3,3)
    INTEGER :: IPIV(3), INFO

! map k1 to IRZ -> k1_IRZ
    NK1_FULL_ORIG=KPOINT_IN_FULL_GRID(VKPT1,KPOINTS_FULL_ORIG)
    NK1_IN_IRZ=KPOINTS_FULL_ORIG%NEQUIV(NK1_FULL_ORIG)

! integer rotation matrix, needss to be transposed
    AMAT=TRANSPOSE(KPOINTS_FULL_ORIG%IROTOP(:,:,NK1_FULL_ORIG))
! LU decomposition
    CALL DGETRF( 3, 3, AMAT, 3, IPIV, INFO )
! solve IROTOP VT= VKPT2
    VT(:)=VKPT2(:)
    CALL DGETRS('N', 3, 1, AMAT, 3, IPIV, VT, 3, INFO)

! search for this k-point in KPOINTS_FULL_ORIG
    NK_FOUND=KPOINT_IN_FULL_GRID(VT,KPOINTS_FULL_ORIG)

! if LTEST is true compare found k-point with NK1_IN_IRZ
    IF (PRESENT(LTEST)) THEN
       IF (NK_FOUND/=NK1_IN_IRZ) THEN
          WRITE(0,*) 'internal error in BRING_TO_IRZ: NK_FOUND is not idential to NK1_FULL_ORIG ', NK_FOUND, NK1_IN_IRZ
          CALL M_exit(); stop
       ENDIF
    END IF
! return k-point vector
    BRING_TO_IRZ=KPOINTS_FULL_ORIG%VKPT(:,NK_FOUND)
    RETURN

! old version

! now loop over all k-points in KPOINTS_FULL_ORIG and determine the
! 1._q that maps onto the supplied k-point VKPT2 when the symmetry operation
! is applied that takes k1_IRZ to k1
! maybe it would be easier to invert IROTOP and apply that to k2
    NK_FOUND=-1
    DO NK=1,KPOINTS_FULL_ORIG%NKPTS
       VT(1)=KPOINTS_FULL_ORIG%VKPT(1,NK)*KPOINTS_FULL_ORIG%IROTOP(1,1,NK1_FULL_ORIG)+ &
             KPOINTS_FULL_ORIG%VKPT(2,NK)*KPOINTS_FULL_ORIG%IROTOP(2,1,NK1_FULL_ORIG)+ &
             KPOINTS_FULL_ORIG%VKPT(3,NK)*KPOINTS_FULL_ORIG%IROTOP(3,1,NK1_FULL_ORIG)
       VT(2)=KPOINTS_FULL_ORIG%VKPT(1,NK)*KPOINTS_FULL_ORIG%IROTOP(1,2,NK1_FULL_ORIG)+ &
             KPOINTS_FULL_ORIG%VKPT(2,NK)*KPOINTS_FULL_ORIG%IROTOP(2,2,NK1_FULL_ORIG)+ &
             KPOINTS_FULL_ORIG%VKPT(3,NK)*KPOINTS_FULL_ORIG%IROTOP(3,2,NK1_FULL_ORIG)
       VT(3)=KPOINTS_FULL_ORIG%VKPT(1,NK)*KPOINTS_FULL_ORIG%IROTOP(1,3,NK1_FULL_ORIG)+ &
             KPOINTS_FULL_ORIG%VKPT(2,NK)*KPOINTS_FULL_ORIG%IROTOP(2,3,NK1_FULL_ORIG)+ &
             KPOINTS_FULL_ORIG%VKPT(3,NK)*KPOINTS_FULL_ORIG%IROTOP(3,3,NK1_FULL_ORIG)
       IF(  ABS(MOD(VT(1)-VKPT2(1)+6.5_q,1._q)-0.5_q) <= TINY .AND. & 
            ABS(MOD(VT(2)-VKPT2(2)+6.5_q,1._q)-0.5_q) <= TINY .AND. & 
            ABS(MOD(VT(3)-VKPT2(3)+6.5_q,1._q)-0.5_q) <= TINY) THEN
          NK_FOUND=NK
          EXIT
       ENDIF
    ENDDO
    IF (NK_FOUND==-1) THEN
       WRITE(0,*) 'internal error in BRING_TO_IRZ: could not find ',VKPT2
       CALL M_exit(); stop
    ENDIF

  END FUNCTION BRING_TO_IRZ



!*********************************************************************
!
! ENTER_IN_SCREENED_2E enters an entry in the screened two electron
! integral table
! the routine is quite complicated since the storage arrangement
! is not that simple we have calculated
! int dr dr'
!   phi_k1,nb1(r) phi_k1-q,nb2(r)  W(r,r')  phi_k1-q,nb2(r') phi_k1,nb1(r')
! determine k2=k1-q and store entry
!
!*********************************************************************

  SUBROUTINE ENTER_IN_SCREENED_2E(S2E, NQ, NK1, NB1, NB2, SCREENED_TWO_ELECTRON_INTEGRALS, ZENTRY )
    USE kpoints_change

    TYPE (screened_2e_handle) S2E
    INTEGER :: NQ   ! index of q-points
    INTEGER :: NK1  ! index of k
    INTEGER :: NB1  ! band index corresponding to NK1
    INTEGER :: NB2  ! band index corresponding to NK2= NK1-NQ
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:)
    COMPLEX(q) ZENTRY
! local
    INTEGER K1_IN_IRZ
    INTEGER K2_IN_IRZ,  K2_IN_FULL_ORIG_IRZ,  K2_IN_FULL_ORIG, K2_STORE_INDEX
    REAL(q), PARAMETER :: TINY=1.E-3_q

    IF (NQ>SIZE(S2E%K1_IN_IRZ,1)) THEN
       WRITE(0,*) 'internal error in ENTER_IN_SCREENED_2E: NQ too large',NQ
       CALL M_exit(); stop
    ENDIF
    IF (NK1>SIZE(S2E%K1_IN_IRZ,2)) THEN
       WRITE(0,*) 'internal error in ENTER_IN_SCREENED_2E: NK1 too large',NK1
       CALL M_exit(); stop
    ENDIF

    K1_IN_IRZ=S2E%K1_IN_IRZ(NQ, NK1)
    K2_IN_IRZ=S2E%K2_IN_IRZ(NQ, NK1)

    IF (K1_IN_IRZ>SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,2)) THEN
       WRITE(0,*) 'internal error in ENTER_IN_SCREENED_2E: K1_IN_IRZ too large',K1_IN_IRZ
       CALL M_exit(); stop
    ENDIF


    IF (S2E%REQUIRED(NQ, NK1)) THEN
       K2_STORE_INDEX=S2E%K2_STORE_INDEX(K1_IN_IRZ, K2_IN_IRZ)
       IF (K2_STORE_INDEX>SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,4)) THEN
          WRITE(0,*) 'internal error in ENTER_IN_SCREENED_2E: K2_IN_IRZ too large',K2_IN_IRZ
          CALL M_exit(); stop
       ENDIF

       SCREENED_TWO_ELECTRON_INTEGRALS(NB1, K1_IN_IRZ, NB2, K2_STORE_INDEX)=ZENTRY+ &
       SCREENED_TWO_ELECTRON_INTEGRALS(NB1, K1_IN_IRZ, NB2, K2_STORE_INDEX)
!       WRITE(78,'(4I5,4F10.4)') K1_IN_IRZ,NB1,K2_IN_IRZ,NB2,ZENTRY
!test for comparison with HF routine
!      K2_IN_FULL_ORIG=KPOINT_IN_FULL_GRID( &
!           KPOINTS_FULL%VKPT(:,K2_IN_IRZ), &
!           KPOINTS_FULL_ORIG)
!      K2_IN_FULL_ORIG_IRZ=KPOINTS_FULL_ORIG%NEQUIV(K2_IN_FULL_ORIG)
!      IF (NK1==K1_IN_IRZ .AND. NB1<=8 .AND. NB2<=4) &
!      WRITE(*,'(4I5,4F10.4)') K1_IN_IRZ,NB1,K2_IN_FULL_ORIG,NB2, &
!           SCREENED_TWO_ELECTRON_INTEGRALS(NB1, K1_IN_IRZ, NB2, K2_STORE_INDEX)
!end test
    ENDIF
  END SUBROUTINE ENTER_IN_SCREENED_2E

!*********************************************************************
!
! REQUIRED_IN_SCREENED_2E returns true if a two electron integral
! needs to be calculated
! some integrals are redundant due to symmetry
!
!*********************************************************************


  FUNCTION REQUIRED_SCREENED_2E(S2E, NQ, NK1)
    TYPE (screened_2e_handle) S2E
    INTEGER :: NQ, NK1
    COMPLEX(q) ZENTRY
    LOGICAL REQUIRED_SCREENED_2E
! local
    INTEGER K1_IN_IRZ
    REAL(q), PARAMETER :: TINY=1.E-6_q

    IF (NQ>SIZE(S2E%K1_IN_IRZ,1)) THEN
       WRITE(0,*) 'internal error in REQUIRED_IN_SCREENED_2E: NQ too large',NQ
       CALL M_exit(); stop
    ENDIF
    IF (NK1>SIZE(S2E%K1_IN_IRZ,2)) THEN
       WRITE(0,*) 'internal error in REQUIRED_IN_SCREENED_2E: NK1 too large',NK1
       CALL M_exit(); stop
    ENDIF

    REQUIRED_SCREENED_2E=S2E%REQUIRED(NQ, NK1)
  END FUNCTION REQUIRED_SCREENED_2E


!*********************************************************************
!
! CLEANUP_SCREENED_2E divides the stored two electron integrals
! by the number of redundant integrals
! furthermore the integrals are averaged over degenerated eigenvalue
! pairs
! this is not particularly elegant but simple to implement
! alternative solutions:
! ) break the symmetry of the wavefunctions using a small perturbation
!  'induced' by NQ
! ) cleanest: calculate
!  V_(1,1'),2(omega)=
!      sum_G G' (u_1(G') u*_2(G'))* RESPONSEFUN(G',G) (u_1(G) u*_2(G))
!  where 1' is a state with the same eigenenergy as 1
!  this is obviously not straightforward ...
!
!*********************************************************************

  SUBROUTINE CLEANUP_SCREENED_2E(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, ISYM, LKSUM)
    USE wave
    USE kpoints_change
    USE subrot_cluster

    TYPE (wavespin)     W
    TYPE (screened_2e_handle) S2E
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:)
    INTEGER :: ISYM
    LOGICAL :: LKSUM
! local
    COMPLEX(q) :: SUM
    TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
    TYPE (eigenf_cluster),POINTER :: PDEG_CLUSTER

    INTEGER :: NK1, NK2, NK2P, NB1, NB2, ISP
    INTEGER :: NBANDSGW

    NBANDSGW=SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,1)

    CALL M_sum_single(W%WDES%COMM_KINTER, SCREENED_TWO_ELECTRON_INTEGRALS, SIZE(SCREENED_TWO_ELECTRON_INTEGRALS)*2)

!
! divide by symmetry related multiplicity
!
!    DO ISP=1, W%WDES%ISPIN
!       DO NK1=1,KPOINTS_ORIG%NKPTS
!          DO NK2=1, W%WDES%NKPTS
!             IF (S2E%NUMBER_OF_REDUNDANT(NK1, NK2)/=0) THEN
!                NK2P=S2E%K2_STORE_INDEX(NK1, NK2)
!                SCREENED_TWO_ELECTRON_INTEGRALS(:, NK1, :, NK2P, ISP)=SCREENED_TWO_ELECTRON_INTEGRALS(:, NK1, :, NK2P, ISP)/S2E%NUMBER_OF_REDUNDANT(NK1, NK2)
!             ENDIF
!          ENDDO
!       ENDDO
!    ENDDO

    IF (ISYM>0) THEN
    NULLIFY(DEG_CLUSTER)
    CALL FIND_DEG_CLUSTERS(W%WDES, W, DEG_CLUSTER, LDISTRIBUTED_IN = .FALSE.)
!
! average along first index
!
    DO ISP=1, W%WDES%ISPIN
       DO NK1=1,KPOINTS_ORIG%NKPTS

          PDEG_CLUSTER=>DEG_CLUSTER(NK1,ISP)%DEG_CLUSTER
          DO
             IF (.NOT.ASSOCIATED(PDEG_CLUSTER)) EXIT
             IF (PDEG_CLUSTER%BAND_START<NBANDSGW) THEN

             DO NK2=1, W%WDES%NKPTS
                IF (S2E%NUMBER_OF_REDUNDANT(NK1, NK2)/=0) THEN
                NK2P=S2E%K2_STORE_INDEX(NK1, NK2)
                DO NB2=1, W%WDES%NBANDS
                   SUM=0
                   DO NB1=PDEG_CLUSTER%BAND_START,MIN(PDEG_CLUSTER%BAND_END,NBANDSGW)
                      SUM=SUM+SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, NK2P, ISP)
                   ENDDO
                   SUM=SUM/(MIN(PDEG_CLUSTER%BAND_END,NBANDSGW)-PDEG_CLUSTER%BAND_START+1)

                   DO NB1=PDEG_CLUSTER%BAND_START,MIN(PDEG_CLUSTER%BAND_END,NBANDSGW)
                      SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, NK2P, ISP)=SUM
                   ENDDO
                ENDDO
                ENDIF
             ENDDO
             ENDIF

             PDEG_CLUSTER=>PDEG_CLUSTER%NEXT
          ENDDO
       ENDDO
    ENDDO

!
! average along second index
!
    DO ISP=1, W%WDES%ISPIN
       DO NK2=1, W%WDES%NKPTS

          PDEG_CLUSTER=>DEG_CLUSTER(NK2,ISP)%DEG_CLUSTER
          DO
             IF (.NOT.ASSOCIATED(PDEG_CLUSTER)) EXIT
             DO NK1=1,KPOINTS_ORIG%NKPTS
                IF (S2E%NUMBER_OF_REDUNDANT(NK1, NK2)/=0) THEN
                NK2P=S2E%K2_STORE_INDEX(NK1, NK2)
                DO NB1=1, NBANDSGW
                   SUM=0
                   DO NB2=PDEG_CLUSTER%BAND_START,PDEG_CLUSTER%BAND_END
! determine local storage index
                      IF ( MOD(NB2-1,W%WDES%NB_PAR)+1 == W%WDES%NB_LOW) THEN
                         SUM=SUM+SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, 1+(NB2-1)/W%WDES%NB_PAR, NK2P, ISP)
                      ENDIF
                   ENDDO
                   IF (LKSUM) THEN
                      CALL M_sum_z(W%WDES%COMM_INTER, SUM, 1)
                   ELSE
                      CALL M_sum_z(W%WDES%COMM_INTER, SUM, 1)
                   ENDIF

                   SUM=SUM/(PDEG_CLUSTER%BAND_END-PDEG_CLUSTER%BAND_START+1)
                   
                   DO NB2=PDEG_CLUSTER%BAND_START,PDEG_CLUSTER%BAND_END
! determine local storage index
                      IF ( MOD(NB2-1,W%WDES%NB_PAR)+1 == W%WDES%NB_LOW) THEN
                         SCREENED_TWO_ELECTRON_INTEGRALS( NB1, NK1, 1+(NB2-1)/W%WDES%NB_PAR, NK2P, ISP)=SUM
                      ENDIF
                   ENDDO
                ENDDO
                ENDIF
             ENDDO
             PDEG_CLUSTER=>PDEG_CLUSTER%NEXT
          ENDDO
       ENDDO
    ENDDO

    CALL FREE_DEG_CLUSTERS(W%WDES,DEG_CLUSTER)
# 639

    ENDIF
  END SUBROUTINE CLEANUP_SCREENED_2E

!*********************************************************************
!
! This subroutine sets degenerated eigenvalues
! to the average of the eigenvalues in the cluster
!
!*********************************************************************

  SUBROUTINE CLEANUP_CELEN(W)
    USE wave
    USE subrot_cluster

    TYPE (wavespin)     W
    COMPLEX(q) :: SUM
    TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
    TYPE (eigenf_cluster),POINTER :: PDEG_CLUSTER

    INTEGER :: NK, NB, ISP

    NULLIFY(DEG_CLUSTER)
    CALL FIND_DEG_CLUSTERS(W%WDES, W, DEG_CLUSTER, LDISTRIBUTED_IN = .FALSE.)

    DO ISP=1, W%WDES%ISPIN
       DO NK=1, W%WDES%NKPTS

          PDEG_CLUSTER=>DEG_CLUSTER(NK,ISP)%DEG_CLUSTER
          DO
             IF (.NOT.ASSOCIATED(PDEG_CLUSTER)) EXIT
             SUM=0
             DO NB=PDEG_CLUSTER%BAND_START,PDEG_CLUSTER%BAND_END
                SUM=SUM+W%CELTOT(NB, NK, ISP)
             ENDDO
             SUM=SUM/(PDEG_CLUSTER%BAND_END-PDEG_CLUSTER%BAND_START+1)
                   
             DO NB=PDEG_CLUSTER%BAND_START,PDEG_CLUSTER%BAND_END
                W%CELTOT(NB, NK, ISP)=SUM
             ENDDO

             SUM=0
             DO NB=PDEG_CLUSTER%BAND_START,PDEG_CLUSTER%BAND_END
                SUM=SUM+W%FERTOT(NB, NK, ISP)
             ENDDO
             SUM=SUM/(PDEG_CLUSTER%BAND_END-PDEG_CLUSTER%BAND_START+1)
                   
             DO NB=PDEG_CLUSTER%BAND_START,PDEG_CLUSTER%BAND_END
                W%FERTOT(NB, NK, ISP)=SUM
             ENDDO

             PDEG_CLUSTER=>PDEG_CLUSTER%NEXT
          ENDDO
       ENDDO
    ENDDO
    
    CALL FREE_DEG_CLUSTERS(W%WDES,DEG_CLUSTER)
    
  END SUBROUTINE CLEANUP_CELEN


!*********************************************************************
!
! This subroutine averages a complex helper array
! over degenerated  eigenvalue pairs
!
!*********************************************************************

  SUBROUTINE CLEANUP_CELEN_HELPER(W, NBANDSGW, CELTOT, FERTOT)
    USE wave
    USE subrot_cluster

    TYPE (wavespin)     W
    INTEGER    :: NBANDSGW
    COMPLEX(q), OPTIONAL :: CELTOT(:,:,:)
    REAL(q), OPTIONAL :: FERTOT(:,:,:)

    TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
    TYPE (eigenf_cluster),POINTER :: PDEG_CLUSTER
    COMPLEX(q) :: SUM

    INTEGER :: NK, NB, ISP

    NULLIFY(DEG_CLUSTER)
    CALL FIND_DEG_CLUSTERS(W%WDES, W, DEG_CLUSTER, LDISTRIBUTED_IN = .FALSE.)

    DO ISP=1, W%WDES%ISPIN
       DO NK=1, W%WDES%NKPTS

          PDEG_CLUSTER=>DEG_CLUSTER(NK,ISP)%DEG_CLUSTER
          DO
             IF (.NOT.ASSOCIATED(PDEG_CLUSTER)) EXIT
             IF (PDEG_CLUSTER%BAND_START<=NBANDSGW) THEN
             IF (PRESENT (CELTOT)) THEN
                SUM=0
                DO NB=PDEG_CLUSTER%BAND_START,MIN(PDEG_CLUSTER%BAND_END,NBANDSGW)
                   SUM=SUM+CELTOT(NB, NK, ISP)
                ENDDO
                SUM=SUM/(MIN(PDEG_CLUSTER%BAND_END,NBANDSGW)-PDEG_CLUSTER%BAND_START+1)
                
                DO NB=PDEG_CLUSTER%BAND_START,PDEG_CLUSTER%BAND_END
                   CELTOT(NB, NK, ISP)=SUM
                ENDDO
             ENDIF

             IF (PRESENT(FERTOT)) THEN
                SUM=0
                DO NB=PDEG_CLUSTER%BAND_START,MIN(PDEG_CLUSTER%BAND_END,NBANDSGW)
                   SUM=SUM+FERTOT(NB, NK, ISP)
                ENDDO
                SUM=SUM/(MIN(PDEG_CLUSTER%BAND_END,NBANDSGW)-PDEG_CLUSTER%BAND_START+1)
                
                DO NB=PDEG_CLUSTER%BAND_START,PDEG_CLUSTER%BAND_END
                   FERTOT(NB, NK, ISP)=SUM
                ENDDO
                
             ENDIF

             ENDIF

             PDEG_CLUSTER=>PDEG_CLUSTER%NEXT
          ENDDO
       ENDDO
    ENDDO
    
    CALL FREE_DEG_CLUSTERS(W%WDES,DEG_CLUSTER)
    
  END SUBROUTINE CLEANUP_CELEN_HELPER

!*********************************************************************
!
! SPLINE_2E determines a spline fit for the frequency dependent
! screened two electron integral along the freqency axis
! the value at the outermost grid point is subtracted and
! returned
!
!*********************************************************************


  SUBROUTINE SPLINE_2E(SCREENED_TWO_ELECTRON_INTEGRALS, OMEGA, S2E_SPLINE )
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:)
    REAL(q) :: OMEGA(:)
! local
    INTEGER :: NOMEGA
    TYPE (screened_2e_spline) S2E_SPLINE

    IF (SIZE(OMEGA) /= SIZE(SCREENED_TWO_ELECTRON_INTEGRALS)) THEN
       WRITE(0,*) 'internal error in SPLINE_2E: different boundaries',SIZE(OMEGA),SIZE(SCREENED_TWO_ELECTRON_INTEGRALS)
       CALL M_exit(); stop
    ENDIF

    NOMEGA=SIZE(OMEGA)

    IF ( ASSOCIATED(S2E_SPLINE%REAL_PART) .AND. S2E_SPLINE%NOMEGA/= NOMEGA) THEN
       DEALLOCATE(S2E_SPLINE%REAL_PART, S2E_SPLINE%IMAG_PART)
    ENDIF
    IF ( .NOT. ASSOCIATED(S2E_SPLINE%REAL_PART)) THEN
       ALLOCATE(S2E_SPLINE%REAL_PART(NOMEGA,5))
       ALLOCATE(S2E_SPLINE%IMAG_PART(NOMEGA,5))
    ENDIF
    
    S2E_SPLINE%NOMEGA=NOMEGA

! set the frequencies
    S2E_SPLINE%REAL_PART(:,1)= OMEGA
    S2E_SPLINE%IMAG_PART(:,1)= OMEGA

! set the energies
    S2E_SPLINE%REAL_PART(:,2)= REAL(SCREENED_TWO_ELECTRON_INTEGRALS,q)
    S2E_SPLINE%IMAG_PART(:,2)= AIMAG(SCREENED_TWO_ELECTRON_INTEGRALS)
    
    
    CALL SPLCOF(S2E_SPLINE%REAL_PART(1,1),NOMEGA,NOMEGA,0.0_q)
! for the imaginary part I am not certain
! even or odd
! use the natural boundary conditions for the time being
    CALL SPLCOF(S2E_SPLINE%IMAG_PART(1,1),NOMEGA,NOMEGA,1E30_q)

  END SUBROUTINE SPLINE_2E

!*********************************************************************
!
! SET_OMEGA_GRID sets the grid for the frequency integration
! if OMEGAWEIGHT is supplied Gauss weights are returned
! this is only sensible for integration along the complex frequency
! axis since Gauss Legendre integration fails for non analytical
! functions
!
!*********************************************************************

  SUBROUTINE SET_OMEGA_GRID(OMEGAMAX, BANDGAPE, OMEGATL, SIGMA, OMEGA, IU6,  OMEGAWEIGHT, VERSION)
    USE gauss_quad
    REAL(q) :: OMEGAMAX                  ! transition energy up to which dense grid is used
    REAL(q) :: OMEGAMAX_LOCAL
    REAL(q) :: BANDGAPE                  ! band gap
    REAL(q) :: OMEGATL                   ! trandition energy up to which coarse grid is used
    REAL(q) :: SIGMA                     ! broadening parameter
    REAL(q) :: OMEGA(:)                  ! frequencies
    REAL(q), OPTIONAL :: OMEGAWEIGHT(:)  ! Gauss weights
    INTEGER, OPTIONAL :: VERSION         ! type of grid
    INTEGER :: IU6
! local
    LOGICAL LGAUSS                       ! set gauss weights
    INTEGER NOMEGA, I, IFAIL
    REAL(q) :: O, A, B
    EXTERNAL GAUSSI2
    INTEGER :: GRIDVERSION               ! controls what kind of grid is used
! 1 linear grid x
! 2 linear x + x^2
! 3 linear x +switch/(1-x^2)
! 4 exponential  - ln (1-x) / alpha
! switch is set such that 1/4 of the grid points are used for the tail
    REAL(q) :: OMIN                      
    REAL(q) :: SWITCH

    IF (PRESENT(VERSION)) THEN
       GRIDVERSION=VERSION
    ELSE
       GRIDVERSION=5
    ENDIF

    NOMEGA=SIZE(OMEGA)
! quick return if no grid points
    IF (NOMEGA==0) RETURN

    IF (PRESENT(OMEGAWEIGHT)) THEN
       LGAUSS=.TRUE.
       IF (NOMEGA/=SIZE(OMEGAWEIGHT)) THEN
          WRITE(0,*) 'internal error in SET_OMEGA_GRID: OMEGAWEIGHT and OMEGA have a different size', & 
               SIZE(OMEGA),SIZE(OMEGAWEIGHT)
          CALL M_exit(); stop
       ENDIF
    ELSE
       LGAUSS=.FALSE.
    ENDIF

    IF (NOMEGA==1) THEN
       OMEGA(1)=0
    ELSE IF (LGAUSS) THEN
! linear to cubic
! GAUSSI supplies the lower bound as last value in OMEGA, which is
! not desired, invert boundaries
# 886

       CALL GAUSS_LEGENDRE(0.0_q, 1.0_q, OMEGA, OMEGAWEIGHT, NOMEGA)
! DO I = 1,NOMEGA
!    WRITE(*,'(4F14.6)')OMEGA(I),OMEGAWEIGHT(I)
! ENDDO

       
       IF (IFAIL/=0) THEN
          WRITE(0,*) 'internal error in SET_OMEGA_GRID: NOMEGA is not compatible with Gauss integration'
          CALL M_exit(); stop
       ENDIF

! linear version
       IF (GRIDVERSION==1) THEN
          DO I=1,NOMEGA
             OMEGA(I)=  OMEGAMAX*OMEGA(I)
             OMEGAWEIGHT(I)=OMEGAWEIGHT(I)*OMEGAMAX
          ENDDO
       ELSE IF (GRIDVERSION==2) THEN
! linear to cubic
          DO I=1,NOMEGA
             O=OMEGA(I)
             OMEGA(I)=  OMEGAMAX*(O+O**3)/2.0
! derivative of function
             OMEGAWEIGHT(I)=OMEGAWEIGHT(I)*OMEGAMAX*(1+3*O**2)/2.0
          ENDDO
       ELSE IF (GRIDVERSION==3) THEN
! linear to cubic
          DO I=1,NOMEGA
             O=OMEGA(I)
             OMEGA(I)=  OMEGAMAX*(O+O**2)/2.0
! derivative of function
             OMEGAWEIGHT(I)=OMEGAWEIGHT(I)*OMEGAMAX*(1+2*O)/2.0
          ENDDO
       ELSE IF (GRIDVERSION==41) THEN
! simple logarithmic version
! assuming the function to be integrated behaves like g(x)=e^(-alpha x)
! function is given by the inversion of G(x)=-1/alpha (e^(-alpha x)-1)
! x= - ln (1-y) / alpha with y=[0,1[
          A=-LOG(1-OMEGA(NOMEGA))/OMEGATL
          DO I=1,NOMEGA
             O=OMEGA(I)
             OMEGA(I)=  -LOG(1-O)/A
! derivative of function
             OMEGAWEIGHT(I)=OMEGAWEIGHT(I)*1.0_q/(1-O)/A
          ENDDO
       ELSE IF (GRIDVERSION==4) THEN
! default for AC-FDT
! logarithmic with modified behaviour for small y
! this maintains the original logarithmic form for y approaching 1 (x->infty)
! but adds or subtracts a quadratic function at small y
! for b=0 original version is obtained. Now x is given by
! x= - ln ( (1-b) (1-y) + b (1-y)**2 ) / alpha with y=[0,1[
! for negative b the number of grid points close to 0._q is increased
! this usually yields much improved convergence compared to the previous version
          B= -0.75
          O= OMEGA(NOMEGA)
          A=-LOG((1-B)*(1-O)+B*(1-O)**2)/OMEGATL
          DO I=1,NOMEGA
             O=OMEGA(I)
             OMEGA(I)=  -LOG((1-B)*(1-O)+B*(1-O)**2)/A
! derivative of function
             OMEGAWEIGHT(I)=OMEGAWEIGHT(I)*1.0_q/((1-B)*(1-O)+B*(1-O)**2)/A*((1-B)+2*B*(1-O))
          ENDDO
       ELSE IF (GRIDVERSION==40) THEN
! logarithmic with modified behaviour for small y
! this 1._q choses B automatically such that the smallest value on
! the grid matches SIGMA (the broadening width)
! as a rule of thumb:
!  1._q gets exact integrals if SIGMA equals the band gap / 5
          OMIN=SIGMA
! better to supply some default (you never know what the users do)
          IF (OMIN==0) OMIN=0.1
          
          B= 0
          DO
             O= OMEGA(NOMEGA)
             A=-LOG((1-B)*(1-O)+B*(1-O)**2)/OMEGATL
             O= OMEGA(1)
             O= -LOG((1-B)*(1-O)+B*(1-O)**2)/A
             IF (O<OMIN) EXIT
             B=B-.001
          ENDDO
          DO I=1,NOMEGA
             O=OMEGA(I)
             OMEGA(I)=  -LOG((1-B)*(1-O)+B*(1-O)**2)/A
! derivative of function
             OMEGAWEIGHT(I)=OMEGAWEIGHT(I)*1.0_q/((1-B)*(1-O)+B*(1-O)**2)/A*((1-B)+2*B*(1-O))
          ENDDO

       ELSEIF (GRIDVERSION==45) THEN
!integrate the complete interval [0,inf)
          CALL GAUSS_LEGENDRE(0.0_q, 1.0_q, OMEGA, OMEGAWEIGHT, NOMEGA)
!use nonlinear variable transformation of Gauss grid
          DO I = 1, NOMEGA 
!weights
             OMEGAWEIGHT(I)=OMEGAWEIGHT(I)/((1._q-OMEGA(I))**2)
!abscissas
             OMEGA(I)=OMEGA(I)/(1._q-OMEGA(I))
          ENDDO 
       ELSEIF (GRIDVERSION==46) THEN
!integrate the interval [0,OMEGATL)
          CALL GAUSS_LEGENDRE(0.0_q, 1.0_q, OMEGA, OMEGAWEIGHT, NOMEGA)
!use linear variable transformation of Gauss grid
          DO I = 1, NOMEGA 
!weights
             OMEGAWEIGHT(I)=OMEGAWEIGHT(I)*OMEGATL
!abscissas
             OMEGA(I)=OMEGA(I)*OMEGATL
          ENDDO 
       ELSEIF (GRIDVERSION==47) THEN
!integrate the interval [0,OMEGATL)
          CALL GAUSS_LEGENDRE(0.0_q, 1.0_q, OMEGA, OMEGAWEIGHT, NOMEGA)
!use linear variable transformation of Gauss grid
          DO I = 1, NOMEGA 
!weights
             OMEGAWEIGHT(I)=OMEGAWEIGHT(I)*OMEGATL
!abscissas
             OMEGA(I)=OMEGA(I)*OMEGATL
          ENDDO 
       ELSE IF (GRIDVERSION==44) THEN
! logarithmic
! assuming the function to be integrated behaves like g(x)=x^((1-B)/B) e^(-alpha x^(1/B))
! function is given by the inversion of G(x)=-(e^(-(alpha x)^(1/B))-1)
! x= (-ln (1-y))^B / alpha with y=[0,1[
! for B = 1.0 GRIDVERSION=4 is recovered
          B = 1.8_q
! search for optimal B, such that lowest frequency matches SIGMA
!DO
          A=(-LOG(1-OMEGA(NOMEGA)))**B/OMEGATL
!   O= (-LOG(1-OMEGA(1)))**B/A
!   IF (O>SIGMA) EXIT
!   B=B*0.999
!ENDDO
          
          DO I=1,NOMEGA
             O=OMEGA(I)
             OMEGA(I)= (-LOG(1-O))**B/A
! derivative of function
             OMEGAWEIGHT(I)=OMEGAWEIGHT(I)*B*(-LOG(1-O))**(B-1.0_q)*(1.0_q/(1-O))/A
          ENDDO
       ELSE IF (GRIDVERSION==50) THEN
! nonlinear
! w' = w/(1+w)  (0<= w <inf) -> \int_0^inf dw f(w) = \int_0^1 dw'/(1-w')^2 f((w'/(1-w'))
          DO I = 1, NOMEGA
             OMEGAWEIGHT(I) = OMEGAWEIGHT(I)/((1-OMEGA(I))**2)
             OMEGA(I) = OMEGA(I)/(1._q-OMEGA(I))
          ENDDO
       ENDIF
    ELSE
! simple case generate simple grids
! linear version
       IF (GRIDVERSION==1) THEN
          DO I=1,NOMEGA
             OMEGA(I)=(OMEGAMAX*(I-1))/MAX(NOMEGA-1,1)
          ENDDO
       ELSE IF (GRIDVERSION==2) THEN
! linear to cubic
          DO I=1,NOMEGA
             A=I
             OMEGA(I)=(OMEGAMAX)*((A-1)/(NOMEGA-1)+(A-1)**3/(NOMEGA-1)**3)/2
          ENDDO
       ELSE IF (GRIDVERSION==3) THEN
! linear to 1/quadratic
! 1/quadratic is the analytical behaviour at large w
          SWITCH=5
          DO
             O=1
             B=(O+SWITCH/(1-O**2/1.1)-SWITCH)
             O=REAL(NOMEGA*3._q/4,q)/(NOMEGA-1)
             O=(OMEGATL)*(O+SWITCH/(1-O**2/1.1)-SWITCH)/B
             IF (O>OMEGAMAX) EXIT
             SWITCH=SWITCH-0.01
          ENDDO
          
          DO I=1,NOMEGA
             O=REAL(I-1,q)/(NOMEGA-1)
             OMEGA(I)=(OMEGATL)*(O+SWITCH/(1-O**2/1.1)-SWITCH)/B
          ENDDO
       ELSE  IF (GRIDVERSION==4) THEN
! 1 /sqrt(w^2+1/(1-y))
          DO
             O=1
             B=(O+SWITCH/(1-O**2/1.1)-SWITCH)
             O=REAL(NOMEGA*3._q/4,q)/(NOMEGA-1)
             O=(OMEGATL)*(O+SWITCH/(1-O**2/1.1)-SWITCH)/B
             IF (O>OMEGAMAX) EXIT
             SWITCH=SWITCH-0.01
          ENDDO
          
          DO I=1,NOMEGA
             O=REAL(I-1,q)/(NOMEGA-1)
             OMEGA(I)=(OMEGATL)*(O+SWITCH/(1-O**2/1.1)-SWITCH)/B
          ENDDO
       ELSE  IF (GRIDVERSION==5) THEN
          OMEGAMAX_LOCAL=OMEGAMAX
          DO I=1,10
! assume g(x)=(1+x^2)/(1+x^4)
! function is given by the inversion of G(x)=\int=0,x g(x') dx'
             CALL SET_OMEGA_GRID_NUM(OMEGAMAX_LOCAL, OMEGATL, OMEGA, GRIDVERSION)
             IF (OMEGA(2) < MAX(BANDGAPE/1.1,0.05_q) ) EXIT
! can not achive target
             OMEGAMAX_LOCAL=OMEGAMAX_LOCAL/1.1
          ENDDO
       ENDIF

    ENDIF
  END SUBROUTINE SET_OMEGA_GRID

!*********************************************************************
!
! set the grid numerically
! assuming that the function to be integrated behaves like
!  g(x)
! the optimal grid is given by the inverse of the integral
!  G(x)=\int=0,x g(x') dx'
! we assume that
!  g(x) = (1+x^2)/ (1+x^4)
! The function has a maximum slightly larger than 1 at around
! 0.7 and decays like 1/x^2, which is the analytical behaviour
! at larger x
!
!*********************************************************************
  SUBROUTINE SET_OMEGA_GRID_NUM(OMEGAMAX, OMEGATL, OMEGA, VERSION)
    REAL(q) :: OMEGAMAX
    REAL(q) :: OMEGATL                   ! tail modelling
    REAL(q) :: OMEGA(:)                  ! frequencies
    INTEGER :: VERSION                   ! type of grid
    INTEGER :: IU6
! local
    INTEGER, PARAMETER ::  N=1000
    REAL(q) :: RHO(N), RHO_INT(10000), RHO_INT_MAX, RHO_INT_TARGET
    REAL(q) :: O, D, SUM, X
    INTEGER :: I, NO, NOMEGA
    
    NOMEGA=SIZE(OMEGA)
! set up g(x) = (1+x^2)/(1+x^4)
    D=OMEGATL/(N-1)
    DO I=1,N
       O=(I-1)*D
       X=O/OMEGAMAX
       RHO(I)= (1+X**2)/(1+X**4)
    ENDDO
! calculate the integral G(x)=\int=0,x g(x') dx' using trapezoidal rule
    SUM=0
    RHO_INT(1)=0
    DO I=2,N
       SUM=SUM+D*(RHO(I-1)+RHO(I))/2
       RHO_INT(I)=SUM
    ENDDO
! assume linear spacing in RHO_INT array
! map that onto an integerindex
! which than corresponds to the frequency
    RHO_INT_MAX=RHO_INT(N)

    OMEGA(1)=0
    DO NO=2,NOMEGA-1
       RHO_INT_TARGET=RHO_INT_MAX*(NO-1)/(NOMEGA-1)
       DO I=1,N
          IF (RHO_INT_TARGET<=RHO_INT(I)) THEN
! linear interpolation
             OMEGA(NO)=((RHO_INT_TARGET-RHO_INT(I-1))/(RHO_INT(I)-RHO_INT(I-1))+I-2)/(N-1)*OMEGATL
             EXIT
          ENDIF
       ENDDO
    ENDDO
    OMEGA(NOMEGA)=OMEGATL
  END SUBROUTINE SET_OMEGA_GRID_NUM



!*********************************************************************
! This routine sets the imaginary grids for the green function
! GRIDVERSIONS
! least square and minimax grids are available
!*********************************************************************
  SUBROUTINE SET_IMAG_GRIDS(BANDGAPE, MAXTRANSE,&
                            NTAU, TAU , TAUWEIGHT, &
                            NNU,  NU_COS, NU_COSWEIGHT, FTCOS, MAXITER_FT, &
                            IO, GRIDVERSION , ERRORS, LSPACETIME, &
                            NU_SIN, NU_SINWEIGHT, FTSIN )
    USE base
    USE gauss_quad
    USE minimax
    USE constant
    IMPLICIT NONE
    REAL(q) :: BANDGAPE                  !band gap
    REAL(q) :: MAXTRANSE                 !maximum transition considered in chi
    INTEGER :: NNU, NTAU                 ! number of tau and nu points
    REAL(q) :: TAU(:), TAUWEIGHT(:)      ! time grid and weighitng points
    REAL(q) :: NU_COS(:),NU_COSWEIGHT(:) ! frequency grid and weighting points (cos related)
    REAL(q) :: FTCOS(:,:)                ! Cos transformation matrix tau -> NU
    INTEGER :: MAXITER_FT 
    TYPE(in_struct) IO
    INTEGER :: GRIDVERSION               ! type of grid
    REAL(q) :: ERRORS(3+2*NNU)
    LOGICAL :: LSPACETIME                ! routine executed form spacetime GW or reciprocal code?
    REAL(q),OPTIONAL :: NU_SIN(:),NU_SINWEIGHT(:) ! frequency grid and weighting points (sin related)
    REAL(q),OPTIONAL :: FTSIN(:,:)       ! Sin transformation matrix tau -> NU
   
    IF (NNU /= NTAU ) THEN 
     WRITE(*,*)'Error in SET_IMAG_GRIDS, NNU and NTAU differ!',NNU,NTAU
     CALL M_exit(); stop
    ENDIF
  
!---------------------------------------------------------------------------
! minimax grids:
! for RPA and MP2 the cos and time grid with fourier transform is sufficient
    IF ( GRIDVERSION >= 140 .AND. GRIDVERSION < 150) THEN
!---------------------------------------------------------------------------
!the spacetime RPA code needs the time points
       IF ( LSPACETIME ) THEN
           IF ( PRESENT( NU_SIN) .AND. PRESENT(NU_SINWEIGHT) .AND. PRESENT(FTSIN) ) THEN
              CALL CALCULATE_MINIMAX_GRIDS (BANDGAPE, MAXTRANSE, NTAU, TAU , TAUWEIGHT, &
                 NNU, NU_COS , NU_COSWEIGHT, FTCOS,MAXITER_FT, IO, GRIDVERSION,.TRUE.,.TRUE.,ERRORS,&
                 NU_SIN, NU_SINWEIGHT, FTSIN)
           ELSE
              CALL CALCULATE_MINIMAX_GRIDS (BANDGAPE, MAXTRANSE, NTAU, TAU , TAUWEIGHT, &
                 NNU, NU_COS , NU_COSWEIGHT, FTCOS,MAXITER_FT, IO, GRIDVERSION,.TRUE.,.TRUE., ERRORS)
           ENDIF 
!the reciprocal RPA code does not need the time points
       ELSE
           CALL CALCULATE_MINIMAX_GRIDS (BANDGAPE, MAXTRANSE, NTAU, TAU , TAUWEIGHT, &
              NNU, NU_COS , NU_COSWEIGHT, FTCOS,MAXITER_FT, IO, GRIDVERSION,.TRUE.,.FALSE., ERRORS)
       ENDIF 
!---------------------------------------------------------------------------
! for testing use only frequency grid for R = 100
    ELSE IF ( GRIDVERSION == 199 ) THEN
!---------------------------------------------------------------------------
       IF ( PRESENT( NU_SIN) ) THEN
!          CALL GENERATE_MATSUBARA_GRIDS(NTAU, TAU, NU_COS, NU_SIN)
          CALL GENERATE_UNIFORM_GRIDS(MAXTRANSE,NTAU, TAU,TAUWEIGHT,&
            NU_COS, NU_COSWEIGHT, NU_SIN, NU_SINWEIGHT) 
       ELSE
          CALL CALCULATE_MINIMAX_GRIDS (1._q, 100._q , NTAU, TAU , TAUWEIGHT, &
            NNU, NU_COS , NU_COSWEIGHT, FTCOS, MAXITER_FT, IO, 1,.FALSE.,.TRUE.,ERRORS)
       ENDIF 
    ENDIF

  RETURN
ENDSUBROUTINE SET_IMAG_GRIDS


!*********************************************************************
!
!  Calculates uniform grids
!
!*********************************************************************
   
  SUBROUTINE GENERATE_UNIFORM_GRIDS(BANDMAX,NTAU,TAU,TAUW,&
     NU_COS,NU_COSW, NU_SIN, NU_SINW)
     USE prec    
     USE gauss_quad
     INTEGER :: NTAU
     REAL(q) :: BANDMAX 
     REAL(q) :: TAU(:), NU_COS(:), NU_SIN(:)
     REAL(q) :: TAUW(:), NU_COSW(:), NU_SINW(:)
!local
     INTEGER :: I, J 
     REAL(q),PARAMETER :: SIGMA=0.1_q
     REAL(q),PARAMETER :: BETA=1._q  !inverse temperature
   
     CALL GAUSS_LEGENDRE(0._q, BANDMAX, TAU,TAUW,NTAU)
     NU_COS=TAU
     NU_SIN=TAU
     NU_SINW=TAUW
     NU_COSW=TAUW
  ENDSUBROUTINE GENERATE_UNIFORM_GRIDS

!*********************************************************************
!
! INTEGRATE_W_2E_SPLINE
! integrate the dynamically screened two electron integrals
! over frequency to determine
!
!  sigma(w)=
!  i/(2 pi) int dw' W(w')/(w'+w-e_2+ i delta sign(e_2-mu)) [e (id w')]
!
! this version uses a spline fit, but since the spectral
! function is rather spiky, this it hardly an improvement
!
!*********************************************************************

  SUBROUTINE INTEGRATE_W_2E_SPLINE(SCREENED_TWO_ELECTRON_INTEGRALS, OMEGA, & 
       OMEGARES, SIGMA, EPSILON, ISIGN, LSHIFT_BOUNDARY_TO_ZERO, SHIFT)
    USE constant
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:)  ! screened two electron integrals
    REAL(q) :: OMEGA(:)                          ! frequencies for the screened two electron integrals
    REAL(q) :: OMEGARES(:)                       ! energies at which the selfenergy is required
    COMPLEX(q) :: SIGMA(:)                       ! self energies
    REAL(q) :: EPSILON                           ! eigenvalue
    REAL(q) :: ISIGN                             ! sign(epsilon-mu)
    LOGICAL :: LSHIFT_BOUNDARY_TO_ZERO           ! set boundary values to 0._q
    REAL(q) :: SHIFT                             ! complex shift
! local
    INTEGER :: NOMEGA_INTERPOLATE
    TYPE (screened_2e_spline) S2E_SPLINE
    REAL(q) :: OMEGAMAX, OMEGAP
    REAL(q) :: DELTA                             ! step width for integration
    REAL(q) :: FDER, FREAL, FIMAG
    COMPLEX(q) :: S2EMAX
    INTEGER :: I, I_R

    DELTA=SHIFT/4

    NULLIFY(S2E_SPLINE%REAL_PART)
    CALL SPLINE_2E(SCREENED_TWO_ELECTRON_INTEGRALS, OMEGA, S2E_SPLINE )

    IF (LSHIFT_BOUNDARY_TO_ZERO) THEN
       S2EMAX=SCREENED_TWO_ELECTRON_INTEGRALS( S2E_SPLINE%NOMEGA )
    ELSE
       S2EMAX=0
    ENDIF

    IF (SIZE(OMEGARES) /= SIZE(SIGMA)) THEN
       WRITE(0,*) 'internal error in SPLINE_2E: different boundaries', SIZE(OMEGARES), SIZE(SIGMA)
       CALL M_exit(); stop
    ENDIF
    
    OMEGAMAX=OMEGA(SIZE(OMEGA))

    NOMEGA_INTERPOLATE=OMEGAMAX/DELTA
    
    SIGMA=0
    DO I=0,NOMEGA_INTERPOLATE
       OMEGAP=I*DELTA
       CALL SPLVAL(OMEGAP, FREAL, FDER, S2E_SPLINE%REAL_PART(1,1), S2E_SPLINE%NOMEGA, S2E_SPLINE%NOMEGA)
       CALL SPLVAL(OMEGAP, FIMAG, FDER, S2E_SPLINE%IMAG_PART(1,1), S2E_SPLINE%NOMEGA, S2E_SPLINE%NOMEGA)

       IF (I/=0) THEN
       DO I_R=1,SIZE(OMEGARES)
          SIGMA(I_R)=SIGMA(I_R)+ &
               (CMPLX(FREAL,FIMAG,q)-S2EMAX)*( &
! contribution from negative frequencies
          1.0_q/(-OMEGAP+OMEGARES(I_R)-EPSILON+ISIGN*CMPLX(0._q,SHIFT,q))+ &
! contribution from positive frequencies
          1.0_q/(+OMEGAP+OMEGARES(I_R)-EPSILON+ISIGN*CMPLX(0._q,SHIFT,q)))
       ENDDO
       ELSE
       DO I_R=1,SIZE(OMEGARES)
          SIGMA(I_R)=SIGMA(I_R)+ &
               (CMPLX(FREAL,FIMAG,q)-S2EMAX)*( &
          1.0_q/(-OMEGAP+OMEGARES(I_R)-EPSILON+ISIGN*CMPLX(0._q,SHIFT,q)))
       ENDDO
       ENDIF
    ENDDO
! multiply by i/(2pi)
    SIGMA=SIGMA*CMPLX(0._q,0.5_q/PI*DELTA,q)

! add integral from a constant function
!  i/(2 pi) int dw' S2EMAX / (w'+w-e_2 + i delta sign(e_2-mu)) e (id w')
!  = - S2EMAX Theta(-e_2-mu) =  - S2EMAX Theta (-ISIGN)
    IF (ISIGN==-1.0_q) THEN
       SIGMA=SIGMA-S2EMAX
    ENDIF

    DEALLOCATE(S2E_SPLINE%REAL_PART, S2E_SPLINE%IMAG_PART)
    
  END SUBROUTINE INTEGRATE_W_2E_SPLINE


!*********************************************************************
!
! INTEGRATE_W_2E_SIMPLE
! integrate the dynamically screened two electron integrals
! over frequency to determine
!
!  sigma(w)=
!  i/(2 pi) int dw' W(w')/(w-w'-e_2+ i delta sign(e_2-mu))
!
! This is Equ. (10) Shishkin, Kresse, PRB 74, 035101
!
! this is the crudest version using a linear interpolation of W(w')
! the frequencies w are supplied in the array OMEGARES
! the result is returned in SIGMA
! the screened two electron integrals and the eigenvalue epsilon
! are supplied by the calling routine
! this routine is rather crude, but it seems to work reasonably;
! for a frequency independent W it should essentially yield -1 for
! occupied states and 0 for unoccupied states
! this is achieved by analytically performing the
! integration for the constant contribution
!
! INTEGRATE_W_2E_SPECTRAL is more "accurate" for spectral properties,
! although, QP shifts converge fast with NOMEGA using this routine
!
!*********************************************************************

  SUBROUTINE INTEGRATE_W_2E_SIMPLE(SCREENED_TWO_ELECTRON_INTEGRALS, OMEGA, & 
       OMEGARES, SIGMA, EPSILON, ISIGN, LSHIFT_BOUNDARY_TO_ZERO, SHIFT)
    USE constant
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:)  ! screened two electron integrals
    REAL(q) :: OMEGA(:)                          ! frequencies for the screened two electron integrals
    REAL(q) :: OMEGARES(:)                       ! energies at which the selfenergy is required
    COMPLEX(q) :: SIGMA(:)                       ! self energies
    REAL(q) :: EPSILON                           ! eigenvalue
    REAL(q) :: ISIGN                             ! sign(epsilon-mu)
    LOGICAL :: LSHIFT_BOUNDARY_TO_ZERO           ! set boundary values to 0._q, and integrate
! constant shift using an analytical formula
    REAL(q) :: SHIFT                             ! complex shift
! local
    INTEGER :: NOMEGA_INTERPOLATE
    REAL(q) :: OMEGAMAX, OMEGAP
    REAL(q) :: DELTA                             ! step width for integration
    COMPLEX(q) :: F
    COMPLEX(q) :: S2EMAX
    INTEGER :: I, I_R

    DELTA=SHIFT/2

    IF (LSHIFT_BOUNDARY_TO_ZERO) THEN
       S2EMAX=SCREENED_TWO_ELECTRON_INTEGRALS( SIZE(OMEGA) )
    ELSE
       S2EMAX=0
    ENDIF

    IF (SIZE(OMEGARES) /= SIZE(SIGMA)) THEN
       WRITE(0,*) 'internal error in INTEGRATE_W_2E_SIMPLE: different boundaries', SIZE(OMEGARES), SIZE(SIGMA)
       CALL M_exit(); stop
    ENDIF

    OMEGAMAX=OMEGA(SIZE(OMEGA))
    
    NOMEGA_INTERPOLATE=OMEGAMAX/DELTA

    SIGMA=0
    DO I=0,NOMEGA_INTERPOLATE
       OMEGAP=I*DELTA

       CALL LINEAR_INTERPOLATE(OMEGAP, F, OMEGA(:), SCREENED_TWO_ELECTRON_INTEGRALS(:) , SIZE(OMEGA))
       F=F-S2EMAX
       IF (I/=0) THEN
          DO I_R=1,SIZE(OMEGARES)
          
             SIGMA(I_R)=SIGMA(I_R)+ &
! contribution from negative frequencies
                  F/(+OMEGARES(I_R)-OMEGAP-EPSILON+ISIGN*CMPLX(0._q,SHIFT,q))+ &
! contribution from positive frequencies
                  F/(+OMEGARES(I_R)+OMEGAP-EPSILON+ISIGN*CMPLX(0._q,SHIFT,q))
          ENDDO
       ELSE
          DO I_R=1,SIZE(OMEGARES)
             SIGMA(I_R)=SIGMA(I_R)+ &
                  F/(OMEGARES(I_R)-EPSILON+ISIGN*CMPLX(0._q,SHIFT,q))
          ENDDO
       ENDIF
    ENDDO
! multiply by i/(2pi)
    SIGMA=SIGMA*CMPLX(0._q,0.5_q/PI*DELTA,q)

! add integral from a constant function
!  i/(2 pi) int dw' S2EMAX / (w'+w-e_2 + i delta sign(e_2-mu)) e (id w')
!  = - S2EMAX Theta(-e_2-mu) =  - S2EMAX Theta (-ISIGN)
    IF (ISIGN==-1.0_q) THEN
       SIGMA=SIGMA-S2EMAX
    ENDIF

  END SUBROUTINE INTEGRATE_W_2E_SIMPLE


!*********************************************************************
!
! INTEGRATE_W_2E_SPECTRAL
! integrate the dynamically screened two electron integrals
! over frequency to determine
!
!  sigma(w)=
!  i/(2 pi) int dw' W(w')/(w-w'-e+ i delta sign(e-mu)) [e (id w')]
!
! this version calculates the spectral function of sigma(w) as:
! ISIGN= -1/ energy e negative:  Imag sigma( w) =  Imag Wres(-w+e)
!                                Imag sigma(-w) =  Imag Wres( w+e)
! ISIGN= +1/ energy e positive:  Imag sigma( w) =  Imag Wres( w-e)
! and determines the real part by a Hilbert transformation
! in the equations above Wres yields only a contribution
! for positive arguments (resonant part of W)
!  sigma(w)=
!  -1/(pi) int dw' Imag W(w')/(w-w'+ i delta sign(e-mu)) [e (id w')]
!
! the relation can be shown by e.g. using a single
! plasmon pole at w_p  W(w')= 1/(w'-w_p+i eta) + 1/(-w'-w_p + i eta)
! and performing the contour integrations over w'
! (at this point in the code poles are below the real axis in the
!  matrix elements stored in W, see note-conjugation in chi_base.F)
!
! for spectral functions this version is an "improvement" of the SIMPLE
! version, however, convergence of QP shifts is worse with this version
!
!*********************************************************************

  SUBROUTINE INTEGRATE_W_2E_SPECTRAL(SCREENED_TWO_ELECTRON_INTEGRALS, OMEGA, & 
       OMEGARES, SIGMA, EPSILON, ISIGN, LSHIFT_BOUNDARY_TO_ZERO, SHIFT)
    USE constant
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:)  ! screened two electron integrals
    REAL(q) :: OMEGA(:)                          ! frequencies for the screened two electron integrals
    REAL(q) :: OMEGARES(:)                       ! energies at which the selfenergy is required
    COMPLEX(q) :: SIGMA(:)                       ! self energies
    REAL(q) :: EPSILON                           ! eigenvalue
    REAL(q) :: ISIGN                             ! sign(epsilon-mu)
    LOGICAL :: LSHIFT_BOUNDARY_TO_ZERO           ! set boundary values to 0._q, and integrate
! constant shift using an analytical formula
    REAL(q) :: SHIFT                             ! complex shift
! local
    INTEGER :: NOMEGA_INTERPOLATE
    REAL(q) :: OMEGAMAX, OMEGAP
    REAL(q) :: DELTA                             ! step width for integration
    COMPLEX(q) :: F
    COMPLEX(q) :: S2EMAX
    INTEGER :: I, I_R

    DELTA=SHIFT/2


    IF (LSHIFT_BOUNDARY_TO_ZERO) THEN
       S2EMAX=SCREENED_TWO_ELECTRON_INTEGRALS( SIZE(OMEGA) )
    ELSE
       S2EMAX=0
    ENDIF

    IF (SIZE(OMEGARES) /= SIZE(SIGMA)) THEN
       WRITE(0,*) 'internal error in INTEGRATE_W_2E_SPECTRAL: different boundaries', SIZE(OMEGARES), SIZE(SIGMA)
       CALL M_exit(); stop
    ENDIF
    
    OMEGAMAX=OMEGA(SIZE(OMEGA))
    
    NOMEGA_INTERPOLATE=OMEGAMAX/DELTA

    SIGMA=0
    DO I=0,NOMEGA_INTERPOLATE
       OMEGAP=I*DELTA

       CALL LINEAR_INTERPOLATE(OMEGAP, F, OMEGA(:), SCREENED_TWO_ELECTRON_INTEGRALS(:) , SIZE(OMEGA))
       F=AIMAG(F-S2EMAX) ! take imaginary part (spectral function)
       IF (ISIGN<0) THEN
          DO I_R=1,SIZE(OMEGARES)
! Hilbert transformation of spectral function for
! sigma(w') =  W(-w'+epsilon)
             SIGMA(I_R)=SIGMA(I_R)+ &
               F/(OMEGARES(I_R)+OMEGAP-EPSILON+ISIGN*CMPLX(0._q,SHIFT,q))
          ENDDO
       ELSE
          DO I_R=1,SIZE(OMEGARES)
! Hilbert transformation of spectral function for
! sigma(w') =  W(w'-epsilon)
             SIGMA(I_R)=SIGMA(I_R)+ &
! contribution from positive frequencies
               F/(OMEGARES(I_R)-OMEGAP-EPSILON+ISIGN*CMPLX(0._q,SHIFT,q))
          ENDDO
       ENDIF
    ENDDO
! multiply by 1/(pi)
    SIGMA=-SIGMA/PI*DELTA

! add integral from a constant function
!  i/(2 pi) int dw' S2EMAX / (w'+w-e_2 + i delta sign(e_2-mu)) e (id w')
!  = - S2EMAX Theta(-e_2-mu) =  - S2EMAX Theta (-ISIGN)
    IF (ISIGN==-1.0_q) THEN
       SIGMA=SIGMA-S2EMAX
    ENDIF

  END SUBROUTINE INTEGRATE_W_2E_SPECTRAL

!*********************************************************************
!
! INTEGRATE_W_2E_SPECTRAL_IMAG
! integrate the dynamically screened two electron integrals
! over real frequency to determine
!
!  sigma( i w)=
!  i/(2 pi) int dw' W(w')/(i w-w'-e_2 )
!
! along the imaginary axis from the poles of the selfenery at the real axis
! see above
!
!*********************************************************************

  SUBROUTINE INTEGRATE_W_2E_SPECTRAL_IMAG(SCREENED_TWO_ELECTRON_INTEGRALS, OMEGA, & 
       OMEGARES, SIGMA, EPSILON, ISIGN, LSHIFT_BOUNDARY_TO_ZERO, SHIFT)
    USE constant
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:)  ! screened two electron integrals
    REAL(q) :: OMEGA(:)                          ! frequencies for the screened two electron integrals
    REAL(q) :: OMEGARES(:)                       ! energies at which the selfenergy is required
    COMPLEX(q) :: SIGMA(:)                       ! self energies
    REAL(q) :: EPSILON                           ! eigenvalue
    REAL(q) :: ISIGN                             ! sign(epsilon-mu)
    LOGICAL :: LSHIFT_BOUNDARY_TO_ZERO           ! set boundary values to 0._q, and integrate
! constant shift using an analytical formula
    REAL(q) :: SHIFT                             ! complex shift
! local
    INTEGER :: NOMEGA_INTERPOLATE
    REAL(q) :: OMEGAMAX, OMEGAP
    REAL(q) :: DELTA                             ! step width for integration
    COMPLEX(q) :: F
    COMPLEX(q) :: S2EMAX
    INTEGER :: I, I_R

    DELTA=SHIFT/2
  
    IF (LSHIFT_BOUNDARY_TO_ZERO) THEN
       S2EMAX=SCREENED_TWO_ELECTRON_INTEGRALS( SIZE(OMEGA) )
    ELSE
       S2EMAX=0
    ENDIF

    IF (SIZE(OMEGARES) /= SIZE(SIGMA)) THEN
       WRITE(0,*) 'internal error in SPLINE_2E: different boundaries', SIZE(OMEGARES), SIZE(SIGMA)
       CALL M_exit(); stop
    ENDIF
    
    OMEGAMAX=OMEGA(SIZE(OMEGA))
    
    NOMEGA_INTERPOLATE=OMEGAMAX/DELTA

    SIGMA=0
    DO I=0,NOMEGA_INTERPOLATE
       OMEGAP=I*DELTA

       CALL LINEAR_INTERPOLATE(OMEGAP, F, OMEGA(:), SCREENED_TWO_ELECTRON_INTEGRALS(:) , SIZE(OMEGA))
       F=AIMAG(F-S2EMAX)
       IF (ISIGN<0) THEN
          DO I_R=1,SIZE(OMEGARES)
! Hilbert transformation of spectral function for
! sigma(w') =  W(-w'+epsilon)
             SIGMA(I_R)=SIGMA(I_R)+ &
               F/(CMPLX(0._q,OMEGARES(I_R),q)+OMEGAP-EPSILON)
          ENDDO
       ELSE
          DO I_R=1,SIZE(OMEGARES)
! Hilbert transformation of spectral function for
! sigma(w') =  W(w'-epsilon)
             SIGMA(I_R)=SIGMA(I_R)+ &
! contribution from positive frequencies
               F/(CMPLX(0._q,OMEGARES(I_R),q)-OMEGAP-EPSILON)
          ENDDO
       ENDIF
    ENDDO
! multiply by 1/(pi)
    SIGMA=-SIGMA*1.0/PI*DELTA

! add integral from a constant function
!  i/(2 pi) int dw' S2EMAX / (w'+w-e_2 + i delta sign(e_2-mu)) e (id w')
!  = - S2EMAX Theta(-e_2-mu) =  - S2EMAX Theta (-ISIGN)
    IF (ISIGN==-1.0_q) THEN
       SIGMA=SIGMA-S2EMAX
    ENDIF

  END SUBROUTINE INTEGRATE_W_2E_SPECTRAL_IMAG


!*********************************************************************
!
! small helper routine that interpolates values FP defined on
! a grid XP to a specified point X
! the result is returned in F
!
!*********************************************************************

  SUBROUTINE LINEAR_INTERPOLATE(X, F, XP, FP, NAC)
    USE prec
    USE main_mpi
    REAL(q), INTENT(IN)   :: X       ! point on which F should be evaluated
    COMPLEX(q),INTENT(OUT)::F        ! return: function value
    REAL(q), INTENT(IN)  :: XP(:)    ! grid
    COMPLEX(qs), INTENT(IN):: FP(:)  ! function values
    INTEGER, INTENT(IN) :: NAC       ! number of grid points
! local
    INTEGER I, J, K
    REAL(q) DX, DX0
!  interval bisectioning
    I=1
    J=NAC
    IF (X  <XP(I)) GO TO 60
    IF (X  <XP(J)) GO TO 70
    K=J-1
    GOTO 90
60  K=1
    GOTO 90
70  K=(I+J)/2
    IF(I==K) GOTO 90
    IF (X   <XP(K)) GO TO 80
    I=K
    GOTO 70
80  J=K
    GOTO 70
!
90  CONTINUE
    IF (K+1 > NAC) THEN
       F=0
    ELSE
       DX =X      -XP(K)
       DX0=XP(K+1)-XP(K)
       F   =FP(K+1)*(DX/DX0)+FP(K)*(1.0_q-DX/DX0)
    ENDIF
  END SUBROUTINE LINEAR_INTERPOLATE

!*********************************************************************
!
! INTEGRATE_W_2E_SIMPLE_C
! integrate the dynamically screened two electron integrals
! using complex contour integration
!
!*********************************************************************

  SUBROUTINE INTEGRATE_W_2E_SIMPLE_C(SCREENED_TWO_ELECTRON_INTEGRALS, OMEGA, OMEGAWEIGHT, &
       NOMEGA_REAL, OMEGARES, SIGMA, EPSILON, LSHIFT_BOUNDARY_TO_ZERO, ISIGN)
    USE constant
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:)  ! screened two electron integrals
    REAL(q) :: OMEGA(:)                          ! frequencies for the screened two electron integrals
    REAL(q) :: OMEGAWEIGHT(:)                    ! weights for Gauss integration along complex axis
    INTEGER :: NOMEGA_REAL                       ! number of real frequencies
    REAL(q) :: OMEGARES(:)                       ! energies at which the selfenergy is required
    COMPLEX(q) :: SIGMA(:)                       ! self energies
    REAL(q) :: EPSILON                           ! eigenvalue
    LOGICAL :: LSHIFT_BOUNDARY_TO_ZERO           ! set boundary values to 0._q, and integrate
! constant shift using an analytical formula
    REAL(q) :: ISIGN                             ! sign(epsilon-mu)
! local
    INTEGER :: NOMEGA_INTERPOLATE
    REAL(q) :: OMEGAMAX, DELTA, OMEGAP
    COMPLEX(q) :: F, S2EMAX
    INTEGER :: I, I_R
    REAL(q), PARAMETER :: DEGENERACY_THRESHOLD=1E-2_q

    IF (LSHIFT_BOUNDARY_TO_ZERO) THEN
       S2EMAX=SCREENED_TWO_ELECTRON_INTEGRALS( SIZE(OMEGA) )
    ELSE
       S2EMAX=0
    ENDIF

    IF (SIZE(OMEGARES) /= SIZE(SIGMA)) THEN
       WRITE(0,*) 'internal error in INTEGRATE_W_2E_SIMPLE_C: different boundaries', SIZE(OMEGARES), SIZE(SIGMA)
       CALL M_exit(); stop
    ENDIF
    SIGMA=0

! Gauss integration along imaginary axis
    IF (NOMEGA_REAL+1<SIZE(OMEGA)) THEN
       DO I=NOMEGA_REAL+1,SIZE(OMEGA)-1
          OMEGAP=OMEGA(I)
          F=(SCREENED_TWO_ELECTRON_INTEGRALS(I)-S2EMAX)
          DO I_R=1,SIZE(OMEGARES)
             SIGMA(I_R)=SIGMA(I_R)+OMEGAWEIGHT(I)*F* &
                 (OMEGARES(I_R)-EPSILON)/((OMEGARES(I_R)-EPSILON)**2+OMEGAP**2)
          ENDDO
       ENDDO
       SIGMA=-SIGMA/PI
    ENDIF

    DO I_R=1,SIZE(OMEGARES)
       CALL LINEAR_INTERPOLATE(ABS(OMEGARES(I_R)-EPSILON), F, & 
            OMEGA(1:NOMEGA_REAL), SCREENED_TWO_ELECTRON_INTEGRALS(1:NOMEGA_REAL) , NOMEGA_REAL)
       IF (ISIGN>=0) THEN
! empty state
          IF (ABS(OMEGARES(I_R)-EPSILON)<DEGENERACY_THRESHOLD) THEN
             SIGMA(I_R)=SIGMA(I_R)+(F-S2EMAX)/2
          ELSE IF (OMEGARES(I_R)>EPSILON) THEN
             SIGMA(I_R)=SIGMA(I_R)+(F-S2EMAX)
          ENDIF
       ELSE
! filled state
          IF (ABS(OMEGARES(I_R)-EPSILON)<DEGENERACY_THRESHOLD) THEN
             SIGMA(I_R)=SIGMA(I_R)-(F-S2EMAX)/2
          ELSE IF (OMEGARES(I_R)<EPSILON) THEN
             SIGMA(I_R)=SIGMA(I_R)-(F-S2EMAX)
          ENDIF
       ENDIF
    ENDDO

! add integral from a constant function
!  i/(2 pi) int dw' S2EMAX / (w'+w-e_2 + i delta sign(e_2-mu)) e (id w')
!  = - S2EMAX Theta(-e_2-mu) =  - S2EMAX Theta (-ISIGN)
    IF (ISIGN==-1.0_q) THEN
       SIGMA=SIGMA-S2EMAX
    ENDIF
  END SUBROUTINE INTEGRATE_W_2E_SIMPLE_C


!*********************************************************************
!
! INTEGRATE_W_2E_IMAG
! integrate the dynamically screened two electron integrals
! over the imaginary frequency to determine
!
!  sigma(i w)=
!  i/(2 pi) int diw' W(i w')/(i w'+i w-e_2+ i delta sign(e_2-mu)) [e (-d w')]
! =1/(2 pi) int d w' W(i w')/(w'+ w+i e_2+ delta sign(e_2-mu)) [e (-d w')]
!
! this version uses a spline fit through the data points
! and integrates using simple discretisation
!
!*********************************************************************

  SUBROUTINE INTEGRATE_W_2E_IMAG(SCREENED_TWO_ELECTRON_INTEGRALS, OMEGA, & 
       OMEGARES, SIGMA, EPSILON, ISIGN, LSHIFT_BOUNDARY_TO_ZERO)
    USE main_mpi
    USE constant
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:)  ! screened two electron integrals
    REAL(q) :: OMEGA(:)                          ! frequencies for the screened two electron integrals
    REAL(q) :: OMEGARES(:)                       ! energies at which the selfenergy is required
    COMPLEX(q) :: SIGMA(:)                       ! self energies
    REAL(q) :: EPSILON                           ! eigenvalue
    REAL(q) :: ISIGN                             ! sign(epsilon-mu)
    LOGICAL :: LSHIFT_BOUNDARY_TO_ZERO           ! set boundary values to 0._q, and integrate
! constant shift using an analytical formula
! local
    INTEGER :: NOMEGA_INTERPOLATE
    REAL(q) :: OMEGAMAX, OMEGAP, FREAL, FIMAG, FDER
    TYPE (screened_2e_spline) S2E_SPLINE
    REAL(q) :: DELTA                             ! step width for integration
    COMPLEX(q) :: F
    COMPLEX(q) :: S2EMAX
    INTEGER :: I, I_R

    NULLIFY(S2E_SPLINE%REAL_PART)
    CALL SPLINE_2E(SCREENED_TWO_ELECTRON_INTEGRALS, OMEGA, S2E_SPLINE )

    IF (LSHIFT_BOUNDARY_TO_ZERO) THEN
       S2EMAX=SCREENED_TWO_ELECTRON_INTEGRALS( SIZE(OMEGA) )
    ELSE
       S2EMAX=0
    ENDIF

    IF (SIZE(OMEGARES) /= SIZE(SIGMA)) THEN
       WRITE(0,*) 'internal error in INTEGRATE_W_2E_IMAG: different boundaries', SIZE(OMEGARES), SIZE(SIGMA)
       CALL M_exit(); stop
    ENDIF

    OMEGAMAX=OMEGA(SIZE(OMEGA))
! I am not sure what a sensible value is, this might require further investigation (gK)
    DELTA=(OMEGARES(SIZE(OMEGARES))-OMEGARES(1))/SIZE(OMEGARES)/10
        
    NOMEGA_INTERPOLATE=OMEGAMAX/DELTA

    SIGMA=0
    DO I=0,NOMEGA_INTERPOLATE
       OMEGAP=I*DELTA

       CALL SPLVAL(OMEGAP, FREAL, FDER, S2E_SPLINE%REAL_PART(1,1), S2E_SPLINE%NOMEGA, S2E_SPLINE%NOMEGA)
       CALL SPLVAL(OMEGAP, FIMAG, FDER, S2E_SPLINE%IMAG_PART(1,1), S2E_SPLINE%NOMEGA, S2E_SPLINE%NOMEGA)
       F=CMPLX(FREAL,FIMAG,q)-S2EMAX

       IF (I/=0) THEN
          DO I_R=1,SIZE(OMEGARES)
             SIGMA(I_R)=SIGMA(I_R)+ &
! contribution from negative frequencies
                  F/(-OMEGAP+OMEGARES(I_R)+CMPLX(0.0_q,EPSILON))+ &
! contribution from positive frequencies
                  F/(+OMEGAP+OMEGARES(I_R)+CMPLX(0.0_q,EPSILON))
          ENDDO
       ELSE
! central point only counted once
          DO I_R=1,SIZE(OMEGARES)
             SIGMA(I_R)=SIGMA(I_R)+ &
                  F/(OMEGARES(I_R)+CMPLX(0.0_q,EPSILON))
          ENDDO
       ENDIF
    ENDDO
    SIGMA=SIGMA*CMPLX(0._q,0.5_q/PI*DELTA,q)

! add integral from a constant function (occupied states only)
!  i/(2 pi) int dw' S2EMAX / (w'+w+i e_2 + delta sign(e_2-mu)) e (-d w')
!  = - S2EMAX Theta(-e_2-mu) =  - S2EMAX Theta (-ISIGN)
    IF (ISIGN==-1.0_q) THEN
       SIGMA=SIGMA-S2EMAX
    ENDIF

    DEALLOCATE(S2E_SPLINE%REAL_PART, S2E_SPLINE%IMAG_PART)

  END SUBROUTINE INTEGRATE_W_2E_IMAG


!*********************************************************************
!
! INTEGRATE_W_2E_IMAG_GAUSS
! integrate the dynamically screened two electron integrals
! along the imaginary axis using a Gauss integration
!
!*********************************************************************

  SUBROUTINE INTEGRATE_W_2E_IMAG_GAUSS(SCREENED_TWO_ELECTRON_INTEGRALS, OMEGA, & 
       OMEGAWEIGHT, OMEGARES, SIGMA, EPSILON, ISIGN, LSHIFT_BOUNDARY_TO_ZERO)
    USE main_mpi
    USE constant
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:)  ! screened two electron integrals
    REAL(q) :: OMEGA(:)                          ! frequencies for the screened two electron integrals
    REAL(q) :: OMEGAWEIGHT(:)                    ! omega weight
    REAL(q) :: OMEGARES(:)                       ! energies at which the selfenergy is required
    COMPLEX(q) :: SIGMA(:)                       ! self energies
    REAL(q) :: EPSILON                           ! eigenvalue
    REAL(q) :: ISIGN                             ! sign(epsilon-mu)
    LOGICAL :: LSHIFT_BOUNDARY_TO_ZERO           ! set boundary values to 0._q, and integrate
! constant shift using an analytical formula
! local
    REAL(q) :: OMEGAMAX, OMEGAP, FREAL, FIMAG, FDER
    COMPLEX(q) :: F
    INTEGER :: I, I_R, NOMEGAGAUSS

    NOMEGAGAUSS=SIZE(OMEGA)

    IF (SIZE(OMEGARES) /= SIZE(SIGMA)) THEN
       WRITE(0,*) 'internal error in INTEGRATE_W_2E_IMAG: different boundaries', SIZE(OMEGARES), SIZE(SIGMA)
       CALL M_exit(); stop
    ENDIF
    
    OMEGAMAX=OMEGA(SIZE(OMEGA))
    
    SIGMA=0
!frequency point of interest
    DO I_R=1,SIZE(OMEGARES)
!use Gaussian quadrature
       DO I=1,NOMEGAGAUSS
!scaled abscissa
          OMEGAP=OMEGA(I)
!kernel at OMEGAP
          F=OMEGAWEIGHT(I)*SCREENED_TWO_ELECTRON_INTEGRALS(I)
       
          SIGMA(I_R)=SIGMA(I_R)+ &
! contribution from negative frequencies
               F/(-OMEGAP+OMEGARES(I_R)+EPSILON*(0.0_q,1.0_q))+ &
! contribution from positive frequencies
               F/(+OMEGAP+OMEGARES(I_R)+EPSILON*(0.0_q,1.0_q))
       ENDDO
    ENDDO
    SIGMA=SIGMA*CMPLX(0._q,0.5_q/PI,q)

  END SUBROUTINE INTEGRATE_W_2E_IMAG_GAUSS

!***********************************************************************
!
! calculate selfenergy enclosed between all states
! this is 1._q by integrating the dynamically screened two electron
! integrals over frequency using a Gauss integration or preferably
! (because of the non analytic behaviour) a linear interpolation
!
!***********************************************************************

  SUBROUTINE CALC_SELFENERGY(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, &
       SHIFT, OMEGA, CELTOT_X, LFOCK_ADD, IU6, IU0 , EFERMI)
    USE wave
    USE vaspxml
    TYPE (wavespin)     W
    TYPE (screened_2e_handle) S2E
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,:)
! screened two electron integrals
    REAL(q) :: SHIFT              ! complex shift
    REAL(q) :: OMEGA(:)           ! energies for screened two electron inegrals
    COMPLEX(q) :: CELTOT_X(:,:,:) ! exact exchange contribution
    LOGICAL :: LFOCK_ADD          ! add exact exchange contribution to sigma
    INTEGER :: IU6                ! stdout
    INTEGER :: IU0                ! stdin
    REAL(q) :: EFERMI             ! Fermi level
! local
    INTEGER NOMEGA, NBANDSGW, I, ISPIN, NK1, NB1, NK2, NB2, ISP
    REAL(q) :: ISIGN, EPSILON, ESHIFT
    REAL(q),ALLOCATABLE :: OMEGARES(:)
    REAL(q), ALLOCATABLE :: A(:,:)
    COMPLEX(q),ALLOCATABLE :: SIGMA(:), SIGMA_IMAG(:), SIGMAFROMK(:)

    ESHIFT=0
!    CALL SHIFT_CELEN( W, ESHIFT)

    NOMEGA=SIZE(OMEGA)
    NBANDSGW=SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,1)

    ALLOCATE(OMEGARES(NOMEGARES), SIGMA(NOMEGARES), SIGMA_IMAG(NOMEGARES), SIGMAFROMK(NOMEGARES), A(3,NOMEGARES))

    DO I=1, NOMEGARES
       OMEGARES(I)=OMEGAMINR+(OMEGAMAXR-OMEGAMINR)/(NOMEGARES-1)*(I-1)
    ENDDO

    IF (IU0>=0) WRITE(IU0,'(A,2F7.2)') ' calculating selfenergy CALC_SELFENERGY between w=',OMEGARES(1),OMEGARES(NOMEGARES)
    IF (IU6>=0) WRITE(IU6,'(A,2F7.2)') ' calculating selfenergy CALC_SELFENERGY between w=',OMEGARES(1),OMEGARES(NOMEGARES)

    DO ISP=1, W%WDES%ISPIN
       IF (W%WDES%ISPIN==2 .AND. IU6>=0) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
! select gamma point only for the test calculations
       DO NK1=1,SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,2)
       IF (W%WDES%WTKPT(NK1)==0) CYCLE ! 0._q k-point weighted points are useless right now in the GW code

       IF (IU6>=0) WRITE(IU6,1) NK1,W%WDES%VKPT(1,NK1),W%WDES%VKPT(2,NK1),W%WDES%VKPT(3,NK1)
       DO NB1=1,MIN(W%WDES%NB_TOTK(NK1,ISP),NBANDSGW)
!          IF (NB1>1) THEN
!             IF (ABS(W%CELTOT(NB1, NK1, ISP)-W%CELTOT(NB1-1, NK1, ISP))<1E-2) CYCLE
!          ENDIF

          SIGMA=0
          SIGMA_IMAG=0
          DO NK2=1,W%WDES%NKPTS
             IF (S2E%NUMBER_OF_REDUNDANT(NK1, NK2)/=0) THEN
                DO NB2=1,W%WDES%NBANDS
                   IF (W%AUX(NB2, NK2, ISP)==0 ) THEN
                      CYCLE     ! bands that should not be included are bypassed
                   ENDIF

                   IF ( W%FERWE(NB2, NK2, ISP) > 0.5)  THEN
                      ISIGN=-1.0
                   ELSE
                      ISIGN= 1.0
                   ENDIF
                   EPSILON=W%CELEN(NB2, NK2, ISP)-EFERMI

                   SIGMAFROMK=0
                   IF (LSPECTRALGW) THEN
                      CALL INTEGRATE_W_2E_SPECTRAL( & 
                           SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,:), &
                           OMEGA, OMEGARES, SIGMAFROMK, EPSILON, ISIGN, .NOT. LFOCK_ADD, SHIFT)
                   ELSE
                      CALL INTEGRATE_W_2E_SIMPLE( &
                           SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,:), &
                           OMEGA, OMEGARES, SIGMAFROMK, EPSILON, ISIGN, .NOT. LFOCK_ADD, SHIFT)
                   ENDIF

! half occupied states should give no contribution use smooth interpolation
                   IF (ISIGN==-1) THEN
                      SIGMAFROMK=SIGMAFROMK*(2*W%FERWE(NB2, NK2, ISP)-1)
                   ELSE
                      SIGMAFROMK=SIGMAFROMK*(1-2*W%FERWE(NB2, NK2, ISP))
                   ENDIF
                   SIGMA=SIGMA+SIGMAFROMK*S2E%WTKPT(NK1, NK2)

                   CALL INTEGRATE_W_2E_SPECTRAL_IMAG( & 
                        SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,:), &
                        OMEGA, OMEGARES, SIGMAFROMK, EPSILON, ISIGN, .NOT. LFOCK_ADD, SHIFT)

! half occupied states should give no contribution use smooth interpolation
                   IF (ISIGN==-1) THEN
                      SIGMAFROMK=SIGMAFROMK*(2*W%FERWE(NB2, NK2, ISP)-1)
                   ELSE
                      SIGMAFROMK=SIGMAFROMK*(1-2*W%FERWE(NB2, NK2, ISP))
                   ENDIF
                   SIGMA_IMAG=SIGMA_IMAG+SIGMAFROMK*S2E%WTKPT(NK1, NK2)

                ENDDO
             ENDIF
          ENDDO

          CALL M_sum_z(W%WDES%COMM_INTER, SIGMA, SIZE(SIGMA))
          CALL M_sum_z(W%WDES%COMM_INTER, SIGMA_IMAG, SIZE(SIGMA_IMAG))

!          WRITE(200+NB1,'(3F14.7)') (OMEGARES(I),SIGMA(I),I=1,NOMEGARES)
!          WRITE(300+NB1,'(3F14.7)') (OMEGARES(I),SIGMA_IMAG(I),I=1,NOMEGARES)

! add HF contribution
          IF (LFOCK_ADD) THEN
             SIGMA=SIGMA+CELTOT_X(NB1, NK1, ISP)
             SIGMA_IMAG=SIGMA_IMAG+CELTOT_X(NB1, NK1, ISP)
          ENDIF

          SIGMA=SIGMA-ESHIFT
          SIGMA_IMAG=SIGMA_IMAG-ESHIFT

          IF (IU6>=0) THEN
               WRITE(IU6,'(1X,I6,3X,F10.4,3X,F10.5,A)') NB1,REAL( W%CELTOT(NB1,NK1,ISP) ,KIND=q) ,W%FERTOT(NB1,NK1,ISP)*W%WDES%RSPIN, & 
                    '  selfenergy along real axis'
               WRITE(IU6,'(3F14.7)') (OMEGARES(I),SIGMA(I),I=1,NOMEGARES,10)
               WRITE(IU6,'(1X,I6,3X,F10.4,3X,F10.5,A)') NB1,REAL( W%CELTOT(NB1,NK1,ISP) ,KIND=q) ,W%FERTOT(NB1,NK1,ISP)*W%WDES%RSPIN, & 
                    '  selfenergy along imaginary axis'
               WRITE(IU6,'(3F14.7)') (OMEGARES(I),SIGMA_IMAG(I),I=1,NOMEGARES,10)
               WRITE(IU6,*)
               A(1,:)=OMEGARES
               A(2,:)=REAL(SIGMA,q)
               A(3,:)=AIMAG(SIGMA)
               CALL XML_VECARRAY("selfenergy")
               CALL XML_ARRAY_REAL(A)
               CALL XML_CLOSE_TAG

               A(2,:)=REAL(SIGMA_IMAG,q)
               A(3,:)=AIMAG(SIGMA_IMAG)
               CALL XML_VECARRAY("selfenergy along imaginary axis")
               CALL XML_ARRAY_REAL(A)
               CALL XML_CLOSE_TAG
          ENDIF
       ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(OMEGARES, SIGMA, SIGMA_IMAG, SIGMAFROMK, A)

 1  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
          &         '  band No.  band energies     occupation '/)

  END SUBROUTINE CALC_SELFENERGY

!***********************************************************************
!
! calculate correlation contribution to the quasiparticle energy
! using a contour deformation
! from my point of view pretty useless
! information along real axis is anyway necessary, and introducing
! a complex contour does just add extra cost without improving
! accuracy (gK)
!
!***********************************************************************

  SUBROUTINE CALC_SELFENERGY_C(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, &
       SHIFT, OMEGA, OMEGAWEIGHT, NOMEGA_REAL, CELTOT_X, LFOCK_ADD, IU6, IU0, EFERMI)
    USE wave
    USE vaspxml
    TYPE (wavespin)     W
    TYPE (screened_2e_handle) S2E
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,:)
! screened two electron integrals
    REAL(q) :: SHIFT              ! complex shift
    REAL(q) :: OMEGA(:)           ! energies for screened two electron inegrals
    REAL(q) :: OMEGAWEIGHT(:)     ! weights for Gauss integration
    INTEGER :: NOMEGA_REAL        ! real frequencies between 1:NOMEGA_REAL
    COMPLEX(q) :: CELTOT_X(:,:,:) ! exact exchange contribution
    LOGICAL :: LFOCK_ADD          ! add exact exchange contribution to sigma
    INTEGER :: IU6                ! stdout
    INTEGER :: IU0                ! stdin
    REAL(q) :: EFERMI
! local
    INTEGER NOMEGA, NBANDSGW, I, ISPIN, NK1, NB1, NK2, NB2, ISP
    REAL(q) :: ISIGN, EPSILON
    REAL(q),ALLOCATABLE :: OMEGARES(:)
    REAL(q), ALLOCATABLE :: A(:,:)
    COMPLEX(q),ALLOCATABLE :: SIGMA(:), SIGMAFROMK(:)

    NOMEGA=SIZE(OMEGA)
    NBANDSGW=SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,1)

    ALLOCATE(OMEGARES(NOMEGARES), SIGMA(NOMEGARES), SIGMAFROMK(NOMEGARES), A(3,NOMEGARES))

    DO I=1, NOMEGARES
       OMEGARES(I)=OMEGAMINR+(OMEGAMAXR-OMEGAMINR)/(NOMEGARES-1)*(I-1)
    ENDDO

    IF (IU0>=0) WRITE(IU0,'(A,2F7.2)') ' calculating selfenergy between w=',OMEGARES(1),OMEGARES(NOMEGARES)
    IF (IU0>=0) WRITE(IU0,'(A,2F7.2)') ' complex contour integration'
    IF (IU6>=0) WRITE(IU6,'(A,2F7.2)') ' calculating selfenergy between w=',OMEGARES(1),OMEGARES(NOMEGARES)

    DO ISP=1, W%WDES%ISPIN
       IF (W%WDES%ISPIN==2 .AND. IU6>=0) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
! select gamma point only for the test calculations
       DO NK1=1,SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,2)
       IF (W%WDES%WTKPT(NK1)==0) CYCLE ! 0._q k-point weighted points are useless right now in the GW code

       IF (IU6>=0) WRITE(IU6,1) NK1,W%WDES%VKPT(1,NK1),W%WDES%VKPT(2,NK1),W%WDES%VKPT(3,NK1)
       DO NB1=1,MIN(W%WDES%NB_TOTK(NK1,ISP),NBANDSGW)

          SIGMA=0
          DO NK2=1,W%WDES%NKPTS
             IF (S2E%NUMBER_OF_REDUNDANT(NK1, NK2)/=0) THEN
                DO NB2=1,W%WDES%NBANDS
                   IF (W%AUX(NB2, NK2, ISP)==0 ) CYCLE ! bands that should not be included are bypassed

                   IF ( W%FERWE(NB2, NK2, ISP) > 0.5)  THEN
                      ISIGN=-1.0
                   ELSE
                      ISIGN= 1.0
                   ENDIF
                   EPSILON=W%CELEN(NB2, NK2, ISP)-EFERMI
                   CALL INTEGRATE_W_2E_SIMPLE_C( & 
                        SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,:), &
                        OMEGA, OMEGAWEIGHT, NOMEGA_REAL, OMEGARES, SIGMAFROMK, EPSILON,.NOT. LFOCK_ADD, ISIGN)
                   IF (ISIGN==-1) THEN
                      SIGMAFROMK=SIGMAFROMK*(2*W%FERWE(NB2, NK2, ISP)-1)
                   ELSE
                      SIGMAFROMK=SIGMAFROMK*(1-2*W%FERWE(NB2, NK2, ISP))
                   ENDIF
                   SIGMA=SIGMA+SIGMAFROMK*S2E%WTKPT(NK1, NK2)
                ENDDO
             ENDIF
          ENDDO
          CALL M_sum_z(W%WDES%COMM_INTER, SIGMA, SIZE(SIGMA))

          IF (LFOCK_ADD) THEN
             SIGMA=SIGMA+CELTOT_X(NB1, NK1, ISP)
          ENDIF

          IF (IU6>=0) THEN
               WRITE(IU6,'(1X,I6,3X,F10.4,3X,F10.5,A)') NB1,REAL( W%CELTOT(NB1,NK1,ISP) ,KIND=q) ,W%FERTOT(NB1,NK1,ISP)*W%WDES%RSPIN, & 
                    '  selfenergy along real axis'
             WRITE(IU6,'(3F14.7)') (OMEGARES(I),SIGMA(I),I=1,NOMEGARES,10)
             WRITE(IU6,*)
             A(1,:)=OMEGARES
             A(2,:)=REAL(SIGMA,q)
             A(3,:)=AIMAG(SIGMA)
             CALL XML_VECARRAY("selfenergy")
             CALL XML_ARRAY_REAL(A)
             CALL XML_CLOSE_TAG
          ENDIF
          
       ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(OMEGARES, SIGMA, SIGMAFROMK, A)

 1  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
          &         '  band No.  band energies     occupation '/)

  END SUBROUTINE CALC_SELFENERGY_C

!***********************************************************************
!
! calculate selfenergy enclosed between all states along the imaginary
! frequeny axis
!
!***********************************************************************

  SUBROUTINE CALC_SELFENERGY_IMAG(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, &
       SHIFT, OMEGA, OMEGAWEIGHT, NOMEGA_REAL, CELTOT_X, LFOCK_ADD, IU6, IU0, EFERMI)
    USE wave
    USE ratpolfit
    USE vaspxml
    TYPE (wavespin)     W
    TYPE (screened_2e_handle) S2E
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,:)
! screened two electron integrals
    REAL(q) :: SHIFT              ! complex shift
    REAL(q) :: EFERMI             ! fermi level
    REAL(q) :: OMEGA(:)           ! energies for screened two electron inegrals
    REAL(q) :: OMEGAWEIGHT(:)     ! weights for Gauss integration
    INTEGER :: NOMEGA_REAL        ! real frequencies between 1:NOMEGA_REAL
    COMPLEX(q) :: CELTOT_X(:,:,:) ! exact exchange contribution
    LOGICAL :: LFOCK_ADD          ! add exact exchange contribution to sigma
    INTEGER :: IU6                ! stdout
    INTEGER :: IU0                ! stdin
! local
    INTEGER NOMEGA, NBANDSGW, I, ISPIN, NK1, NB1, NK2, NB2, ISP
    REAL(q) :: ISIGN, EPSILON
    REAL(q),ALLOCATABLE :: OMEGARES(:)
    REAL(q), ALLOCATABLE :: A(:,:)
    COMPLEX(q),ALLOCATABLE :: SIGMA_IMAG(:), SIGMAFROMK_IMAG(:)

    NOMEGA=SIZE(OMEGA)
    NBANDSGW=SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,1)

    ALLOCATE(OMEGARES(NOMEGARES), SIGMA_IMAG(NOMEGARES), & 
      SIGMAFROMK_IMAG(NOMEGARES), A(3,NOMEGARES) )

    DO I=1, NOMEGARES
       OMEGARES(I)=OMEGAMINR+(OMEGAMAXR-OMEGAMINR)/(NOMEGARES-1)*(I-1)
    ENDDO

    IF (IU0>=0) WRITE(IU0,'(A,2F7.2)') ' calculating selfenergy CALC_SELFENERGY_IMAG between i w=',OMEGARES(1),OMEGARES(NOMEGARES)
    IF (IU6>=0) WRITE(IU6,'(A,2F7.2)') ' calculating selfenergy CALC_SELFENERGY_IMAG between i w=',OMEGARES(1),OMEGARES(NOMEGARES)

    DO ISP=1, W%WDES%ISPIN
       IF (W%WDES%ISPIN==2 .AND. IU6>=0) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
! select gamma point only for the test calculations
       DO NK1=1,SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,2)
       IF (W%WDES%WTKPT(NK1)==0) CYCLE ! 0._q k-point weighted points are useless right now in the GW code

       IF (IU6>=0) WRITE(IU6,1) NK1,W%WDES%VKPT(1,NK1),W%WDES%VKPT(2,NK1),W%WDES%VKPT(3,NK1)
       DO NB1=1,MIN(W%WDES%NB_TOTK(NK1,ISP),NBANDSGW)
!          IF (NB1>1) THEN
!             IF (ABS(W%CELTOT(NB1, NK1, ISP)-W%CELTOT(NB1-1, NK1, ISP))<1E-6) CYCLE
!          ENDIF

          SIGMA_IMAG=0
          DO NK2=1,W%WDES%NKPTS
             IF (S2E%NUMBER_OF_REDUNDANT(NK1, NK2)/=0) THEN
                DO NB2=1,W%WDES%NBANDS
                   IF (W%AUX(NB2, NK2, ISP)==0 ) CYCLE ! bands that should not be included are bypassed

                   IF ( W%FERWE(NB2, NK2, ISP) > 0.5)  THEN
                      ISIGN=-1.0
                   ELSE
                      ISIGN= 1.0
                   ENDIF
                   EPSILON=W%CELEN(NB2, NK2, ISP)-EFERMI

                   CALL INTEGRATE_W_2E_IMAG( & 
                        SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,NOMEGA_REAL+1:), &
                        OMEGA(NOMEGA_REAL+1:), OMEGARES, SIGMAFROMK_IMAG, & 
                        EPSILON, ISIGN, .NOT. LFOCK_ADD)
!                   CALL INTEGRATE_W_2E_IMAG_GAUSS( &
!                        SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,NOMEGA_REAL+1:), &
!                        OMEGAWEIGHT(NOMEGA_REAL+1:),OMEGA(NOMEGA_REAL+1:), OMEGARES, SIGMAFROMK_IMAG, &
!                        EPSILON, ISIGN, .NOT. LFOCK_ADD)

                   IF (ISIGN==-1) THEN
                      SIGMAFROMK_IMAG=SIGMAFROMK_IMAG*(2*W%FERWE(NB2, NK2, ISP)-1)
                   ELSE
                      SIGMAFROMK_IMAG=SIGMAFROMK_IMAG*(1-2*W%FERWE(NB2, NK2, ISP))
                   ENDIF

                   SIGMA_IMAG=SIGMA_IMAG+SIGMAFROMK_IMAG*S2E%WTKPT(NK1, NK2)
                ENDDO
             ENDIF
          ENDDO
          CALL M_sum_z(W%WDES%COMM_INTER, SIGMA_IMAG, SIZE(SIGMA_IMAG))

!          WRITE(300+NB1,'(3F14.7)') (OMEGARES(I),SIGMA_IMAG(I),I=1,NOMEGARES)

!          CALL FIT_SCREENEDTWOE_HALF( REAL(W%CELTOT(NB1, NK1, ISP),q), OMEGA(NOMEGA_REAL+1:), &
!               SIGMA_IMAG, OMEGAWEIGHT(NOMEGA_REAL+1:), OMEGARES, SIGMA)

! add HF contribution
          IF (LFOCK_ADD) THEN
             SIGMA_IMAG=SIGMA_IMAG+CELTOT_X(NB1, NK1, ISP)
          ENDIF

          IF (IU6>=0) THEN
               WRITE(IU6,'(1X,I6,3X,F10.4,3X,F10.5,A)') NB1,REAL( W%CELTOT(NB1,NK1,ISP) ,KIND=q) ,W%FERTOT(NB1,NK1,ISP)*W%WDES%RSPIN, & 
                    '  selfenergy along imaginary axis'
             WRITE(IU6,'(3F14.7)') (OMEGARES(I),SIGMA_IMAG(I),I=1,NOMEGARES,10)
             WRITE(IU6,*)
             A(1,:)=OMEGARES
             A(2,:)=REAL(SIGMA_IMAG,q)
             A(3,:)=AIMAG(SIGMA_IMAG)
             CALL XML_VECARRAY("selfenergy along imaginary axis")
             CALL XML_ARRAY_REAL(A)
             CALL XML_CLOSE_TAG
          ENDIF
       ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(OMEGARES, SIGMA_IMAG, SIGMAFROMK_IMAG, A)
    CALL M_exit(); stop

 1  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
          &         '  band No.  band energies     occupation '/)
  END SUBROUTINE CALC_SELFENERGY_IMAG

!*********************************************************************
! use a Gauss integration on the imaginary frequency axis
! to get the self-energy on the imaginary axis
!*********************************************************************

  SUBROUTINE CALC_SELFENERGY_IMAG_GAUSS(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, &
       SHIFT, OMEGA, OMEGAWEIGHT, NOMEGA_REAL, CELTOT_X, LFOCK_ADD, IU6, IU0,EFERMI)
    USE wave
    USE ratpolfit
    USE vaspxml
    TYPE (wavespin)     W
    TYPE (screened_2e_handle) S2E
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,:)
! screened two electron integrals
    REAL(q) :: SHIFT              ! complex shift
    REAL(q) :: EFERMI             ! fermi level
    REAL(q) :: OMEGA(:)           ! energies for screened two electron inegrals
    REAL(q) :: OMEGAWEIGHT(:)     ! weights for Gauss integration
    INTEGER :: NOMEGA_REAL        ! real frequencies between 1:NOMEGA_REAL
    COMPLEX(q) :: CELTOT_X(:,:,:) ! exact exchange contribution
    LOGICAL :: LFOCK_ADD          ! add exact exchange contribution to sigma
    INTEGER :: IU6                ! stdout
    INTEGER :: IU0                ! stdin
! local
    INTEGER NOMEGA, NBANDSGW, I, ISPIN, NK1, NB1, NK2, NB2, ISP
    REAL(q) :: ISIGN, EPSILON
    REAL(q),ALLOCATABLE :: OMEGARES(:)
    REAL(q), ALLOCATABLE :: A(:,:)
    COMPLEX(q),ALLOCATABLE :: SIGMA_IMAG(:), SIGMAFROMK_IMAG(:), SIGMA(:)


    IF (W%WDES%COMM_KINTER%NCPU.NE.1) THEN
       CALL M_stop('CALC_SELFENERGY_IMAG: KPAR>1 not implemented, sorry.')
!PK Look over all of screened_2e carefully. Likely just k-point skips and COMM_KINTER sums required.
       CALL M_exit(); stop
    END IF

    NOMEGA=SIZE(OMEGA)
    NBANDSGW=SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,1)

    ALLOCATE(OMEGARES(NOMEGARES), SIGMA_IMAG(NOMEGARES), & 
      SIGMAFROMK_IMAG(NOMEGARES), SIGMA(NOMEGARES), A(3,NOMEGARES) )

    WRITE(*,*)' Using Gauss-Legendre integration '


! don't shift here
!    CALL SHIFT_CELEN( W)

    DO I=1, NOMEGARES
       OMEGARES(I)=OMEGAMINR+(OMEGAMAXR-OMEGAMINR)/(NOMEGARES-1)*(I-1)
    ENDDO

    IF (IU0>=0) WRITE(IU0,'(A,2F7.2)') 'calculating selfenergy CALC_SELFENERGY_IMAG_GAUSS between w=',OMEGARES(1),OMEGARES(NOMEGARES)
    IF (IU0>=0) WRITE(IU0,'(A,2F7.2)') 'analytical continuation from complex frequencies'
    IF (IU6>=0) WRITE(IU6,'(A,2F7.2)') 'calculating selfenergy CALC_SELFENERGY_IMAG_GAUSS between w=',OMEGARES(1),OMEGARES(NOMEGARES)

    DO ISP=1, W%WDES%ISPIN
       IF (W%WDES%ISPIN==2 .AND. IU6>=0) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
! select gamma point only for the test calculations
       DO NK1=1,1  ! SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,2)
       IF (IU6>=0) WRITE(IU6,1) NK1,W%WDES%VKPT(1,NK1),W%WDES%VKPT(2,NK1),W%WDES%VKPT(3,NK1)
       DO NB1=1,NBANDSGW
          IF (NB1>1) THEN
             IF (ABS(W%CELTOT(NB1, NK1, ISP)-W%CELTOT(NB1-1, NK1, ISP))<1E-6) CYCLE
          ENDIF

          SIGMA_IMAG=0
          DO NK2=1,W%WDES%NKPTS
             IF (S2E%NUMBER_OF_REDUNDANT(NK1, NK2)/=0) THEN
                DO NB2=1,W%WDES%NBANDS
                   IF (W%AUX(NB2, NK2, ISP)==0 ) CYCLE ! bands that should not be included are bypassed

                   IF ( W%FERWE(NB2, NK2, ISP) > 0.5)  THEN
                      ISIGN=-1.0
                   ELSE
                      ISIGN= 1.0
                   ENDIF
                   EPSILON=W%CELEN(NB2, NK2, ISP) - EFERMI
                   CALL INTEGRATE_W_2E_IMAG_GAUSS( & 
                        SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,NOMEGA_REAL+1:), &
                        OMEGA(NOMEGA_REAL+1:),OMEGAWEIGHT, OMEGARES, SIGMAFROMK_IMAG, & 
                        EPSILON, ISIGN, .NOT. LFOCK_ADD)
                   IF (ISIGN==-1) THEN
                      SIGMAFROMK_IMAG=SIGMAFROMK_IMAG*(2*W%FERWE(NB2, NK2, ISP)-1)
                   ELSE
                      SIGMAFROMK_IMAG=SIGMAFROMK_IMAG*(1-2*W%FERWE(NB2, NK2, ISP))
                   ENDIF

                   SIGMA_IMAG=SIGMA_IMAG+SIGMAFROMK_IMAG*S2E%WTKPT(NK1, NK2)
                ENDDO
             ENDIF
          ENDDO
          CALL M_sum_z(W%WDES%COMM_INTER, SIGMA_IMAG, SIZE(SIGMA_IMAG))

          SIGMA=0
          WRITE(200+NB1,'(3F14.7)') (OMEGARES(NOMEGA_REAL+I),SIGMA_IMAG(I),I=1,NOMEGARES)

!          CALL FIT_SCREENEDTWOE_HALF( REAL(W%CELTOT(NB1, NK1, ISP),q), OMEGA(NOMEGA_REAL+1:), &
!               SIGMA_IMAG, OMEGAWEIGHT(NOMEGA_REAL+1:), OMEGARES, SIGMA)

          IF (LFOCK_ADD) THEN
             SIGMA=SIGMA+CELTOT_X(NB1, NK1, ISP)
          ENDIF
          IF (IU6>=0) THEN
             WRITE(IU6,'(3F14.7)') (OMEGARES(I),SIGMA(I),I=1,NOMEGARES)
             WRITE(IU6,*)
             A(1,:)=OMEGARES
             A(2,:)=REAL(SIGMA,q)
             A(3,:)=AIMAG(SIGMA)
             CALL XML_VECARRAY("selfenergy")
             CALL XML_ARRAY_REAL(A)
             CALL XML_CLOSE_TAG
          ENDIF
WRITE(*,*)'band ' ,NB1,'done',W%CELTOT(NB1,NK1,ISP)-EFERMI
       ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(OMEGARES, SIGMA_IMAG, SIGMAFROMK_IMAG, SIGMA, A)

 1  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
          &         '  band No.  band energies     occupation '/)
  END SUBROUTINE CALC_SELFENERGY_IMAG_GAUSS

!***********************************************************************
!
! calculate the selfenergy enclosed between all states
! in this version it is assumed that the SCREENED_TWO_ELECTRON_INTEGRALS
! contains already the correct Hilbert transformed W
!
!***********************************************************************

  SUBROUTINE CALC_SELFENERGY_LINEAR(W, S2E, &
       SCREENED_TWO_ELECTRON_INTEGRALS, &
       SHIFT, OMEGA, CELTOT_X, LFOCK_ADD, IU6, IU0, EFERMI)
    USE wave
    USE vaspxml
    TYPE (wavespin)     W
    TYPE (screened_2e_handle) S2E
    COMPLEX(qs):: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,:)
    REAL(q) :: SHIFT              ! complex shift
    REAL(q) :: OMEGA(:)           ! energies for screened two electron inegrals
    COMPLEX(q) :: CELTOT_X(:,:,:) ! exact exchange contribution
    LOGICAL :: LFOCK_ADD          ! add exact exchange contribution to sigma
    INTEGER :: IU6                ! stdout
    INTEGER :: IU0                ! stdin
    REAL(q) :: EFERMI             ! Fermi-level
! local
    COMPLEX(qs) :: SE_MINUS(SIZE(OMEGA))
    COMPLEX(qs) :: SE_PLUS(SIZE(OMEGA))
    INTEGER NOMEGA, NBANDSGW, I, ISPIN, NK1, NB1, NK2, NB2, ISP
    REAL(q) :: ISIGN, EPSILON
    REAL(q),ALLOCATABLE :: OMEGARES(:)
    REAL(q), ALLOCATABLE :: A(:,:)
    COMPLEX(q),ALLOCATABLE :: SIGMA(:), SIGMAFROMK(:)

    IF (.NOT.LFOCK_ADD) THEN
       WRITE(0,*) 'internal error in CALC_SELFENERGY_LINEAR: LFOCK_ADD should be .TRUE., assuming it'
    ENDIF

    NOMEGA=SIZE(OMEGA)
    NBANDSGW=SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,1)

    ALLOCATE(OMEGARES(NOMEGARES), SIGMA(NOMEGARES), SIGMAFROMK(NOMEGARES), A(3,NOMEGARES) )

    IF (NOMEGARES/=1) THEN 
       DO I=1, NOMEGARES
          OMEGARES(I)=OMEGAMINR+(OMEGAMAXR-OMEGAMINR)/(NOMEGARES-1)*(I-1)
       ENDDO
    ENDIF
    IF (IU0>=0) WRITE(IU0,'(A,2F7.2)') ' calculating selfenergy between w=',OMEGARES(1),OMEGARES(NOMEGARES)
    IF (IU6>=0) WRITE(IU6,'(A,2F7.2)') ' calculating selfenergy between w=',OMEGARES(1),OMEGARES(NOMEGARES)

    DO ISP=1, W%WDES%ISPIN
       IF (W%WDES%ISPIN==2 .AND. IU6>=0) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
! select gamma point only for the test calculations
       DO NK1=1,SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,2)
          IF (W%WDES%WTKPT(NK1)==0) CYCLE ! 0._q k-point weighted points are useless right now in the GW code

          IF (IU6>=0) WRITE(IU6,1) NK1,W%WDES%VKPT(1,NK1),W%WDES%VKPT(2,NK1),W%WDES%VKPT(3,NK1)
          DO NB1=1,MIN(W%WDES%NB_TOTK(NK1,ISP),NBANDSGW)
!             IF (NB1>1) THEN
!                IF (ABS(W%CELTOT(NB1, NK1, ISP)-W%CELTOT(NB1-1, NK1, ISP))<1E-2) CYCLE
!             ENDIF

             SIGMA=0
             DO NK2=1,W%WDES%NKPTS
                IF (S2E%NUMBER_OF_REDUNDANT(NK1, NK2)/=0) THEN
                   DO NB2=1,W%WDES%NBANDS
                      IF (W%AUX(NB2, NK2, ISP)==0 ) CYCLE ! bands that should not be included are bypassed

                      IF ( W%FERWE(NB2, NK2, ISP) > 0.5)  THEN
                         ISIGN=-1.0  ! occupied, below Fermi-level
                         SE_MINUS(:)= SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,1+NOMEGA:2*NOMEGA)   
                         SE_PLUS(:) =-SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,1:NOMEGA)
                      ELSE
                         ISIGN= 1.0  ! unoccupied, above Fermi-level
                         SE_MINUS(:)= SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,1:NOMEGA)

                         SE_PLUS(:) =-SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,1+NOMEGA:2*NOMEGA)
                      ENDIF
                      EPSILON=W%CELEN(NB2, NK2, ISP)-EFERMI

                      CALL INTERPOLATE_W_2E_SIMPLE( & 
                           SE_MINUS(:), SE_PLUS(:), &
                           OMEGA, OMEGARES, SIGMAFROMK, EPSILON)
                      SIGMA=SIGMA+SIGMAFROMK*S2E%WTKPT(NK1, NK2)
                   ENDDO
                ENDIF
             ENDDO
             CALL M_sum_z(W%WDES%COMM_INTER, SIGMA, SIZE(SIGMA))

!             WRITE(200+NB1,'(3F14.7)') (OMEGARES(I),SIGMA(I),I=1,NOMEGARES)

! add HF contribution
             SIGMA=SIGMA+CELTOT_X(NB1, NK1, ISP)

             IF (IU6>=0) THEN
                WRITE(IU6,'(1X,I6,3X,F10.4,3X,F10.5,A)') NB1,REAL( W%CELTOT(NB1,NK1,ISP) ,KIND=q) ,W%FERTOT(NB1,NK1,ISP)*W%WDES%RSPIN, & 
                    '  selfenergy along real axis'
                WRITE(IU6,'(3F14.7)') (OMEGARES(I),SIGMA(I),I=1,NOMEGARES,10)
                WRITE(IU6,*)
                A(1,:)=OMEGARES
                A(2,:)=REAL(SIGMA,q)
                A(3,:)=AIMAG(SIGMA)
                CALL XML_VECARRAY("selfenergy")
                CALL XML_ARRAY_REAL(A)
                CALL XML_CLOSE_TAG
             ENDIF

          ENDDO
       ENDDO
       CALL M_exit(); stop
    ENDDO
    DEALLOCATE(OMEGARES, SIGMA, SIGMAFROMK, A)

1   FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
         &         '  band No.  band energies     occupation '/)
2   FORMAT((3X,I4,3X,F10.4,3X,F10.5))

  END SUBROUTINE CALC_SELFENERGY_LINEAR

!*********************************************************************
!
! INTERPOLATE_W_2E_SIMPLE
! interpolate the dynamically screened two electron integrals
!
!*********************************************************************

 SUBROUTINE INTERPOLATE_W_2E_SIMPLE(  &
                           SE_MINUS, SE_PLUS, &
                           OMEGA, OMEGARES, SIGMA, EPSILON )
    USE constant
    COMPLEX(qs) :: SE_MINUS(:)             ! screened two electron integrals for negative frequencies
    COMPLEX(qs) :: SE_PLUS(:)              ! screened two electron integrals for positive frequencies
    REAL(q) :: OMEGA(:)                    ! frequencies for the screened two electron integrals
    REAL(q) :: OMEGARES(:)                 ! energies at which the selfenergy is required
    COMPLEX(q) :: SIGMA(:)                 ! self energies
    REAL(q) :: EPSILON                     ! eigenvalue
! local
    REAL(q) :: OMEGAP
    COMPLEX(q) :: F
    INTEGER :: I

    IF (SIZE(OMEGARES) /= SIZE(SIGMA)) THEN
       WRITE(0,*) 'internal error in SPLINE_2E: different boundaries', SIZE(OMEGARES), SIZE(SIGMA)
       CALL M_exit(); stop
    ENDIF
 
    DO I=1,SIZE(OMEGARES)
       OMEGAP = OMEGARES(I)-EPSILON
       IF (OMEGAP>0.0_q) THEN
          CALL LINEAR_INTERPOLATE( OMEGAP, F, OMEGA(:), SE_PLUS(:) , SIZE(OMEGA))
       ELSE 
          CALL LINEAR_INTERPOLATE(-OMEGAP, F, OMEGA(:), SE_MINUS(:) , SIZE(OMEGA))
       ENDIF
       
       SIGMA(I)=F
    ENDDO
 
  END SUBROUTINE INTERPOLATE_W_2E_SIMPLE


!***********************************************************************
!
! calculate correlation contribution to the quasiparticle energy
! this is 1._q by integrating the dynamically screened two electron
! integrals over frequency at three frequencies centered
! around the eigenvalues
!  sigma_1,2(e_1+-delta)=
!  i/(2 pi) int dw' V_1,2(w') / (w'+w-e_2 + i delta sign(e_2-mu)) e (id w')
! where
!  V_1,2(w')=
!      sum_G G' (u_1(G') u*_2(G'))* RESPONSEFUN(G',G,w') (u_1(G) u*_2(G))
!
!***********************************************************************


  SUBROUTINE QP_SHIFT(W, LSPECTRAL, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, &
      SHIFT, OMEGA, OMEGAWEIGHT, NOMEGA_REAL, CELTOT_HARTREE_KINETIC, CELTOT_X, &
      LFOCK_ADD, NELM, LRESTORE, LG0W0, IU6, IU0)
    USE wave
    USE ratpolfit
    USE kpoints_change

    TYPE (wavespin)     W
    LOGICAL LSPECTRAL
    TYPE (screened_2e_handle) S2E
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,:)  
! screened two electron integrals
    REAL(q) :: SHIFT              ! complex shift
    REAL(q) :: OMEGA(:)           ! energies for screened two electron inegrals
    REAL(q) :: OMEGAWEIGHT(:)     ! weights for Gauss integration
    INTEGER :: NOMEGA_REAL        ! number of points along real axis
    COMPLEX(q) :: CELTOT_HARTREE_KINETIC(:,:,:)
    COMPLEX(q) :: CELTOT_X(:,:,:) ! exact exchange contribution
    LOGICAL :: LFOCK_ADD          ! add exact exchange contribution to sigma
    INTEGER :: NELM               ! iteration counter (supplied by external routine)
    LOGICAL :: LRESTORE           ! restore old eigenvalues after calculation of QP shifts
    LOGICAL :: LG0W0              ! do not update the position of the poles of the 1._q electron Green function
    INTEGER :: IU6                ! stdout
    INTEGER :: IU0                ! stdin
! local
    INTEGER NOMEGA, NBANDSGW, I, ISPIN, NK1, NB1, NK2, NB2, ISP, ITER, NITER
    REAL(q) :: ISIGN, EPSILON, EDIFF, QPEPSILON, Z, ESHIFT
! to determine the Z factor central differences are used
    REAL(q) :: DELTAOMEGA=0.2
    REAL(q) :: OMEGARES(3), SIGMA(3)
    COMPLEX(q) :: SIGMAFROMK(3), SIGMA_CMPLX(3)
    COMPLEX(q), ALLOCATABLE :: SIGMAFROMK_IMAG(:), SIGMA_IMAG(:)
    REAL(q)    :: CELNEW(W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN) 
    REAL(q)    :: CELOLD(W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN)

    NOMEGA=SIZE(OMEGA)
    NBANDSGW=SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,1)

    ALLOCATE(SIGMA_IMAG(NOMEGA-NOMEGA_REAL), SIGMAFROMK_IMAG(NOMEGA-NOMEGA_REAL))

    ESHIFT=0
!    CALL SHIFT_CELEN( W, ESHIFT)
    CELNEW=W%CELTOT
    CELOLD=W%CELTOT

    IF (NOMEGA>2) THEN
       NITER=4
    ELSE
       NITER=1
    ENDIF

    IF (LSPECTRAL .OR. (NOMEGA>2 .AND. NOMEGA_REAL<=1)) THEN
       NITER=1
    ENDIF
! these iterations are only 1._q if LSPECTRAL=.FALSE.
  iteration: DO ITER=1,NITER
    IF (IU0>=0) WRITE(IU0,'(//A,I2)') ' calculate QP shifts <psi_nk| G(iteration)W_0 |psi_nk>: iteration',ITER+NELM-1
    IF (IU6>=0) THEN
       WRITE(IU6,'(//A,I2)') ' QP shifts <psi_nk| G(iteration)W_0 |psi_nk>: iteration',ITER+NELM-1
       IF (ITER+NELM-1==1) THEN
          WRITE(IU6,'(A)') ' for sc-GW calculations column KS-energies equals QP-energies in previous step'
          WRITE(IU6,'(A)') ' and V_xc(KS)=  KS-energies - (<T + V_ion + V_H > + <T+V_H+V_ion>^1  + <V_x>^1)'
       ENDIF
    ENDIF
       
    DO ISP=1, W%WDES%ISPIN
       IF (W%WDES%ISPIN==2 .AND. IU6>=0) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
       DO NK1=1,SIZE(SCREENED_TWO_ELECTRON_INTEGRALS,2)
       IF (W%WDES%WTKPT(NK1)==0) CYCLE ! 0._q k-point weighted points are useless right now in the GW code

       IF ((NOMEGA==2 .OR. NOMEGA==1) .AND.IU6>=0) THEN
          WRITE(IU6,3) & 
            NK1,W%WDES%VKPT(1,NK1),W%WDES%VKPT(2,NK1),W%WDES%VKPT(3,NK1)
       ELSE IF (ITER+NELM-1==1.AND.IU6>=0) THEN
          WRITE(IU6,1) & 
            NK1,W%WDES%VKPT(1,NK1),W%WDES%VKPT(2,NK1),W%WDES%VKPT(3,NK1)
       ELSE IF (IU6>=0) THEN
          WRITE(IU6,2) & 
            NK1,W%WDES%VKPT(1,NK1),W%WDES%VKPT(2,NK1),W%WDES%VKPT(3,NK1)
       ENDIF

       DO NB1=1,NBANDSGW
! OMEGARES are those frequencies at which we want to determine the
! self energy
          OMEGARES(2)=CELNEW(NB1, NK1, ISP)
          OMEGARES(1)=OMEGARES(2)-DELTAOMEGA
          OMEGARES(3)=OMEGARES(2)+DELTAOMEGA

          SIGMA=0
          SIGMA_CMPLX=0
          SIGMA_IMAG=0
          DO NK2=1,W%WDES%NKPTS
             IF (S2E%NUMBER_OF_REDUNDANT(NK1, NK2)/=0) THEN
                DO NB2=1,W%WDES%NBANDS
                   IF (W%AUX(NB2, NK2, ISP)==0 ) CYCLE ! bands that should not be included are bypassed

                   IF ( W%FERWE(NB2, NK2, ISP) > 0.5)  THEN
                      ISIGN=-1.0
                   ELSE
                      ISIGN= 1.0
                   ENDIF
! epsilon (E_n') is the pole of the 1._q-electron Green's function
! for state n'
! the structure of the the selfenergy is
! int w' W(w') G(w+w') = int w' W(w') |psi_n'><psi_n'| 1/(w+w'-E_n')
                   EPSILON=W%CELTOT((NB2-1)*W%WDES%NB_PAR+W%WDES%NB_LOW, NK2, ISP)

! eigenvalue difference (this is only correct if ITER=1)
                   EDIFF = W%CELTOT(NB1,NK1,ISP)-EPSILON

                   SIGMAFROMK=0
                   SIGMAFROMK_IMAG=0
                   IF (LSPECTRAL) THEN
! in the spectral method the first entry in SCREENED_TWO_ELECTRON_INTEGRALS
! stores the selfenergy contributions from the pair NK1 NK2 and the second
! 1._q the derivative
                      SIGMAFROMK(2)=SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 1)
                      SIGMAFROMK(1)=SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 1) & 
                        +DELTAOMEGA*SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 2)
                      SIGMAFROMK(3)=SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 1) & 
                        -DELTAOMEGA*SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 2)
                   ELSE IF (NOMEGA>2 .AND. NOMEGA_REAL<=1) THEN
! analytical continuation
! select if 1._q or less frequencies along real axis are known
                      WRITE(*,*) 'internal error: analytical continuation is not implemented, sorry stopping'
                      CALL M_exit(); stop
                      CALL INTEGRATE_W_2E_IMAG( & 
                           SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,NOMEGA_REAL+1:), &
                           OMEGA(NOMEGA_REAL+1:), OMEGA(NOMEGA_REAL+1:), SIGMAFROMK_IMAG, & 
                           EPSILON, ISIGN, .NOT. LFOCK_ADD)
! half occupied states should give no contribution use smooth interpolation
                      IF (ISIGN==-1) THEN
                         SIGMAFROMK     =SIGMAFROMK*(2*W%FERWE(NB2, NK2, ISP)-1)
                         SIGMAFROMK_IMAG=SIGMAFROMK_IMAG*(2*W%FERWE(NB2, NK2, ISP)-1)
                      ELSE
                         SIGMAFROMK     =SIGMAFROMK*(1-2*W%FERWE(NB2, NK2, ISP))
                         SIGMAFROMK_IMAG=SIGMAFROMK_IMAG*(1-2*W%FERWE(NB2, NK2, ISP))
                      ENDIF
                   ELSE IF (NOMEGA>2 .AND. NOMEGA_REAL/=NOMEGA) THEN
! dynamical selfenergy using complex contour
! select of real and imaginary frequencies are used
                      CALL INTEGRATE_W_2E_SIMPLE_C( & 
                           SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,:), &
                           OMEGA, OMEGAWEIGHT, NOMEGA_REAL, OMEGARES, SIGMAFROMK, EPSILON, .NOT. LFOCK_ADD, ISIGN)
! half occupied states should give no contribution use smooth interpolation
                      IF (ISIGN==-1) THEN
                         SIGMAFROMK=SIGMAFROMK*(2*W%FERWE(NB2, NK2, ISP)-1)
                      ELSE
                         SIGMAFROMK=SIGMAFROMK*(1-2*W%FERWE(NB2, NK2, ISP))
                      ENDIF
                   ELSE IF (NOMEGA>2) THEN
! dynamical selfenergy integral along real axis
! default
                      IF (LSPECTRALGW) THEN
                         CALL INTEGRATE_W_2E_SPECTRAL( & 
                           SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,1:NOMEGA_REAL), &
                           OMEGA(1:NOMEGA_REAL), OMEGARES, SIGMAFROMK, EPSILON, ISIGN, .NOT. LFOCK_ADD, SHIFT)
                      ELSE
                         CALL INTEGRATE_W_2E_SIMPLE( & 
                           SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP,1:NOMEGA_REAL), &
                           OMEGA(1:NOMEGA_REAL), OMEGARES, SIGMAFROMK, EPSILON, ISIGN, .NOT. LFOCK_ADD, SHIFT)
                      ENDIF
! half occupied states should give no contribution use smooth interpolation
                      IF (ISIGN==-1) THEN
                         SIGMAFROMK=SIGMAFROMK*(2*W%FERWE(NB2, NK2, ISP)-1)
                      ELSE
                         SIGMAFROMK=SIGMAFROMK*(1-2*W%FERWE(NB2, NK2, ISP))
                      ENDIF
                   ELSE IF (NOMEGA==2) THEN
! COHSEX
! two frequencies, higher 1._q is assumed to be exchange
! this converges slower then NOMEGA=1 and is mainly for testing
! SEX first
                      SIGMAFROMK(2)=-W%FERWE(NB2, NK2, ISP)* &
                         SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 1)
                      IF (LFOCK_ADD) THEN
! COH term (only correlation W-v)
                         SIGMAFROMK(1)=0.5_q* &
                              (SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 1))
                      ELSE
! COH term (first frequency point stores correlation, second point exchange)
                         SIGMAFROMK(1)=0.5_q* &
                              (SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 1)- &
                               SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 2))
                      ENDIF
                   ELSE
! static case, simple SEX  (W-v) *(-f)
                      SIGMAFROMK(2)=-W%FERWE(NB2, NK2, ISP)* &
                         SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 1)
! and COH (W-v) 0.5
                      SIGMAFROMK(1)=0.5* &
                         SCREENED_TWO_ELECTRON_INTEGRALS(NB1, NK1, NB2, S2E%K2_STORE_INDEX(NK1, NK2), ISP, 1)
                   ENDIF
                   SIGMA      =SIGMA      +SIGMAFROMK*S2E%WTKPT(NK1, NK2)
                   SIGMA_CMPLX=SIGMA_CMPLX+SIGMAFROMK*S2E%WTKPT(NK1, NK2)
                   SIGMA_IMAG =SIGMA_IMAG+SIGMAFROMK_IMAG*S2E%WTKPT(NK1, NK2)
                ENDDO
             ENDIF
          ENDDO
          CALL M_sum_d(W%WDES%COMM_INTER, SIGMA, SIZE(SIGMA))
          CALL M_sum_z(W%WDES%COMM_INTER, SIGMA_CMPLX, SIZE(SIGMA_CMPLX))
          CALL M_sum_z(W%WDES%COMM_INTER, SIGMA_IMAG, SIZE(SIGMA_IMAG))

          IF (NOMEGA>2 .AND. NOMEGA_REAL<=1) THEN
             SIGMAFROMK=0
             CALL FIT_SCREENEDTWOE_HALF( REAL(W%CELTOT(NB1, NK1, ISP),q), OMEGA(NOMEGA_REAL+1:), &
                  SIGMA_IMAG, OMEGAWEIGHT(NOMEGA_REAL+1:), OMEGARES, SIGMAFROMK)
             SIGMA=SIGMAFROMK
          ENDIF

! shift final result by the same value as eigenvalues have been shifted
          SIGMA=SIGMA-ESHIFT
          SIGMA_CMPLX=SIGMA_CMPLX-ESHIFT

! add Fock exchange term to QP energies
          IF (LFOCK_ADD) THEN
             IF (NOMEGA<=2) THEN
! COHSEX only add exchange to SEX term
                SIGMA(2)=SIGMA(2)+CELTOT_X(NB1, NK1, ISP)
                SIGMA_CMPLX(2)=SIGMA_CMPLX(2)+CELTOT_X(NB1, NK1, ISP)
             ELSE
                SIGMA=SIGMA+CELTOT_X(NB1, NK1, ISP)
                SIGMA_CMPLX=SIGMA_CMPLX+CELTOT_X(NB1, NK1, ISP)
             ENDIF
          ENDIF
          QPEPSILON=SIGMA(2)

          IF (NOMEGA>2) THEN
! Z-factor if full self energy has been evaluated
             Z=1/(1-(SIGMA(3)-SIGMA(1))/(2*DELTAOMEGA))
          ELSE
! Z-factor=1 in static case
             Z=1
          ENDIF
! eps_iter+1 = eps_iter  + Z ( sigma_c + sigma_x + <T+V_H+V_ion> +
!                             <T+V_H+V_ion>^1  + <V_x>^1 - eps_iter)
          QPEPSILON=REAL(CELNEW(NB1, NK1, ISP),q)+Z*( QPEPSILON+ &
               REAL(-CELNEW(NB1, NK1, ISP)+CELTOT_HARTREE_KINETIC(NB1, NK1, ISP),q))
          IF (NOMEGA<=2) THEN
             QPEPSILON=QPEPSILON+SIGMA(1)
          ENDIF

          IF ((NOMEGA<=2).AND. IU6>=0) THEN
! two frequencies: static COHSEX
             WRITE(IU6,10)& 
               NB1,REAL(CELNEW(NB1, NK1, ISP),q), QPEPSILON, SIGMA(2), & 
               REAL(CELNEW(NB1, NK1, ISP)-CELTOT_HARTREE_KINETIC(NB1, NK1, ISP), q), & 
               REAL(CELTOT_X(NB1, NK1, ISP), q), SIGMA(1), &
               W%FERTOT(NB1,NK1,ISP)*W%WDES%RSPIN
          ELSE IF (ITER+NELM-1==1 .AND. IU6>=0) THEN
             WRITE(IU6,10) & 
               NB1,REAL(CELNEW(NB1, NK1, ISP),q), QPEPSILON, SIGMA(2), & 
               REAL(CELNEW(NB1, NK1, ISP)-CELTOT_HARTREE_KINETIC(NB1, NK1, ISP), q), & 
               REAL(CELTOT_X(NB1, NK1, ISP), q), Z, &
               W%FERTOT(NB1,NK1,ISP)*W%WDES%RSPIN, IMAG(SIGMA_CMPLX(2))
          ELSE IF (IU6>=0)  THEN
             WRITE(IU6,10) & 
               NB1,REAL(CELNEW(NB1, NK1, ISP),q), QPEPSILON, SIGMA(2), & 
               REAL(CELTOT_HARTREE_KINETIC(NB1, NK1, ISP), q), & 
               REAL(CELTOT_X(NB1, NK1, ISP), q), Z, &
               W%FERTOT(NB1,NK1,ISP)*W%WDES%RSPIN, IMAG(SIGMA_CMPLX(2))
          ENDIF
          CELNEW(NB1, NK1, ISP)=QPEPSILON
       ENDDO
       ENDDO
    ENDDO
! update CELTOT
! CELTOT is used to determine the poles in the 1._q electron Green function
! if this update is *not 1._q*, the poles of the 1._q electron Green function
! are fixed, and only the energy at which the selfenergy is evaluated is updated
    IF (.NOT.LG0W0) THEN
       W%CELTOT=CELNEW
    ENDIF
! update eigenvalues for those k-points that are not in the IRZ
    DO ISP=1, W%WDES%ISPIN
       DO NK1=1, W%WDES%NKPTS
! determine index in original full grid and map into IRZ
          NK2=KPOINTS_FULL_ORIG%NEQUIV(KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,NK1),KPOINTS_FULL_ORIG))
          W%CELTOT(:,NK1,ISP)=W%CELTOT(:,NK2,ISP)
       ENDDO
    ENDDO

    ENDDO iteration

    W%CELTOT=CELNEW
! update eigenvalues for those k-points that are not in the IRZ
    DO ISP=1, W%WDES%ISPIN
       DO NK1=1, W%WDES%NKPTS
          IF (W%WDES%WTKPT(NK1)==0) CYCLE ! 0._q k-point weighted points are useless right now in the GW code
! determine index in original full grid and map into IRZ
          NK2=KPOINTS_FULL_ORIG%NEQUIV(KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,NK1),KPOINTS_FULL_ORIG))
          W%CELTOT(:,NK1,ISP)=W%CELTOT(:,NK2,ISP)
       ENDDO
    ENDDO

    IF (LRESTORE) THEN
       W%CELTOT=CELOLD
    ENDIF
    DEALLOCATE(SIGMA_IMAG, SIGMAFROMK_IMAG)
    

 1  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
          &         "  band No.  KS-energies  QP-energies   sigma(KS)   V_xc(KS)     V^pw_x(r,r')   Z            occupation"/)
 2  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
          &         "  band No. old QP-enery  QP-energies   sigma(KS)   T+V_ion+V_H  V^pw_x(r,r')   Z            occupation"/)
 3  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
          &         "  band No.  KS-energies  QP-energies   SEX         V_xc(KS)     V_x(r,r')      COH          occupation"/)
10  FORMAT((3X,I4,3X,8(F10.4,3X)))

  END SUBROUTINE QP_SHIFT


!************************ SUBROUTINE DETERMINE_SLOT ********************
!
! MAXIM before programming commenting and not afterwards
!
!***********************************************************************

  SUBROUTINE DETERMINE_SLOT( EDIFF, OMEGAA, MOMEGA1, MOMEGA2, WEIGHT1, WEIGHT2)

    USE constant
    IMPLICIT NONE
! MAXIM seperate variables that are passed to subroutine from
! those that are local
    REAL(q)  :: EDIFF                 ! difference between eigenvalues
    REAL(q)  :: OMEGAA(:)             ! frequencies
    INTEGER MOMEGA1, MOMEGA2          ! on return: MAXIM
    REAL(q), OPTIONAL :: WEIGHT1, WEIGHT2 ! on return: MAXIM
! local
! MAXIM: never use REAL always use REAL(q) instead of REAL
    INTEGER :: N, M, MOMEGA
    REAL(q) :: DIFF, OMEGA1, OMEGA2

    IF (EDIFF<OMEGAA(1).OR. EDIFF>OMEGAA(size(OMEGAA))) THEN
       WRITE(*,*) 'ERROR in DETERMINE_SLOT: the frequency range is not sufficient',EDIFF,REAL(OMEGAA(size(OMEGAA)))
       WRITE(*,*) ' please set OMEGATL manually in the INCAR file.'
       WRITE(*,*) ' Usually increasing it by 20 % with respect to the default solves this issue.' 
       CALL M_exit(); stop
    ENDIF

    DO MOMEGA=2,SIZE(OMEGAA)
       IF (OMEGAA(MOMEGA)>EDIFF) EXIT
    ENDDO

    OMEGA1 = OMEGAA(MOMEGA)
    OMEGA2 = OMEGAA(MOMEGA-1)
    MOMEGA1 = MOMEGA
    MOMEGA2 = MOMEGA - 1

! calculate the weights of the density (sum up to 1._q)
    IF (PRESENT(WEIGHT2) .AND. PRESENT(WEIGHT1)) THEN
! weights such that the sum to 1
       WEIGHT1 = (EDIFF-OMEGA2)/(OMEGA1-OMEGA2)
       WEIGHT2 = 1-WEIGHT1
! scale such that the integral is 1._q
       WEIGHT1 = WEIGHT1/(OMEGA1-OMEGA2)
       WEIGHT2 = WEIGHT2/(OMEGA1-OMEGA2)
    ENDIF
  END SUBROUTINE DETERMINE_SLOT


!************************ SUBROUTINE DETERMINE_SLOT_INTER ************
!
! similar to previous version, but returns two weight set of weights
! first  1._q allows to interpolate the value at EDIF [W_INTER(:,1)]
! second 1._q allows to determine the derivative at EDIF [W_INTER(:,2)]
!
!**********************************************************************

  SUBROUTINE DETERMINE_SLOT_INTER( EDIFF, OMEGAA, MOMEGA1, MOMEGA2, W_INTER)

    USE constant
    IMPLICIT NONE
! MAXIM seperate variables that are passed to subroutine from
! those that are local
    REAL(q)  :: EDIFF                 ! MAXIM: documentation
    REAL(q)  :: OMEGAA(:)             ! frequencies
    INTEGER MOMEGA1, MOMEGA2          ! MAXIM: documentation
    REAL(q)  :: W_INTER(2,2)          ! weights for interpolation
! local
! MAXIM: never use REAL always use REAL(q) instead of REAL
    REAL(q) :: ABSEDIFF
    INTEGER :: N, MOMEGA
    REAL(q) :: OMEGA1, OMEGA2

    ABSEDIFF=ABS(EDIFF)
    IF (ABSEDIFF<OMEGAA(1) .OR. ABSEDIFF>OMEGAA(size(OMEGAA))) THEN
       WRITE(*,*) 'ERROR in DETERMINE_SLOT_INTER: the frequency range is not sufficient',ABSEDIFF,REAL(OMEGAA(size(OMEGAA)))
       CALL M_exit(); stop
    ENDIF

    DO MOMEGA=2,SIZE(OMEGAA)
       IF (OMEGAA(MOMEGA)>ABSEDIFF) EXIT
    ENDDO

    OMEGA1 = OMEGAA(MOMEGA)
    OMEGA2 = OMEGAA(MOMEGA-1)
    MOMEGA1 = MOMEGA
    MOMEGA2 = MOMEGA - 1

! weights for interpolation
    W_INTER(1,1) = (ABSEDIFF-OMEGA2)/(OMEGA1-OMEGA2)
    W_INTER(2,1) = 1-W_INTER(1,1)
! weights for derivative
    W_INTER(1,2) =   1/(OMEGA1-OMEGA2)
    W_INTER(2,2) =  -1/(OMEGA1-OMEGA2)

    IF (EDIFF>=0) THEN
       W_INTER(:,1)=-W_INTER(:,1)
    ENDIF

  END SUBROUTINE DETERMINE_SLOT_INTER


  SUBROUTINE DETERMINE_SLOT_INTER_WEIGHT( EDIFF, OMEGAA, MOMEGA1, MOMEGA2, W_INTER, WEIGHT, IERROR)

    USE constant
    IMPLICIT NONE
! MAXIM seperate variables that are passed to subroutine from
! those that are local
    REAL(q)  :: EDIFF                 ! MAXIM: documentation
    REAL(q)  :: OMEGAA(:)             ! frequencies
    INTEGER MOMEGA1, MOMEGA2          ! MAXIM: documentation
    REAL(q)  :: W_INTER(2,2)          ! weights for interpolation
    REAL(q)  :: WEIGHT
    INTEGER  :: IERROR
! local
! MAXIM: never use REAL always use REAL(q) instead of REAL
    REAL(q) :: ABSEDIFF
    INTEGER :: N, MOMEGA
    REAL(q) :: OMEGA1, OMEGA2

    ABSEDIFF=ABS(EDIFF)
    IERROR=0
    IF (ABSEDIFF<OMEGAA(1) .OR. ABSEDIFF>OMEGAA(size(OMEGAA))) THEN
       WRITE(*,*) 'ERROR in DETERMINE_SLOT_INTER_WEIGHT: the frequency range is not sufficient',ABSEDIFF,REAL(OMEGAA(size(OMEGAA)))
       IERROR=1
       RETURN
    ENDIF

    DO MOMEGA=2,SIZE(OMEGAA)
       IF (OMEGAA(MOMEGA)>ABSEDIFF) EXIT
    ENDDO

    OMEGA1 = OMEGAA(MOMEGA)
    OMEGA2 = OMEGAA(MOMEGA-1)
    MOMEGA1 = MOMEGA
    MOMEGA2 = MOMEGA - 1

! weights for interpolation
    W_INTER(1,1) = (ABSEDIFF-OMEGA2)/(OMEGA1-OMEGA2)
    W_INTER(2,1) = 1-W_INTER(1,1)
! weights for derivative
    W_INTER(1,2) =   1/(OMEGA1-OMEGA2)
    W_INTER(2,2) =  -1/(OMEGA1-OMEGA2)

    W_INTER=W_INTER*WEIGHT

    IF (EDIFF>=0) THEN
       W_INTER(:,1)=-W_INTER(:,1)
    ENDIF

  END SUBROUTINE DETERMINE_SLOT_INTER_WEIGHT

!*******************  READ_SCREENED_2E_FILE  **************************
!
! rudimentary routine to read and write TWOEINT file
!
! READ_SCREENED_2E_FILE
! read the screened two electron integral if possible from the
! file
! IERR is 0._q if the operation was sucessfull and otherwise
! IERR is set to a value larger than 0
!
! WRITE_SCREENED_2E_FILE
! write the file
!
! particularly in the parallel mode, this routine is a real hack
! each node reads and writes to a different file
! if 1._q fails it returns an error
!
!**********************************************************************

  SUBROUTINE READ_SCREENED_2E_FILE(SCREENED_TWO_ELECTRON_INTEGRALS, IERR, IU0, & 
      NOMEGA, NOMEGAR, OMEGAMAX, OMEGATL, NODE)
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,:)
    INTEGER :: IERR, IU0        
    INTEGER, OPTIONAL :: NODE   ! append this integer to TWOEINT file name
! local
    INTEGER :: N(6), N_FILE(6)
    INTEGER :: IU=72
    CHARACTER (2) :: APP
    INTEGER :: NOMEGA, NOMEGAR, NOMEGA_FILE, NOMEGAR_FILE
    REAL(q) :: OMEGAMAX, OMEGATL, OMEGAMAX_FILE, OMEGATL_FILE

    N=SHAPE(SCREENED_TWO_ELECTRON_INTEGRALS)

    APP="  "
    IF (PRESENT(NODE)) THEN
       WRITE (APP  , "(I1,I1)") MOD(node/10,10),MOD(node,10)
       OPEN( UNIT=IU, FILE="TWOEINT"//APP, STATUS="OLD", IOSTAT=IERR)
    ELSE
       OPEN( UNIT=IU, FILE="TWOEINT", STATUS="OLD", IOSTAT=IERR)
    ENDIF

    IF (IERR==0) THEN
       IF (IU0>=0) WRITE(IU0,*) 'reading TWOEINT'//APP
          
       READ(IU,'(2I10,2E14.7)', IOSTAT=IERR) NOMEGA_FILE, NOMEGAR_FILE, OMEGAMAX_FILE, OMEGATL_FILE
       IF (NOMEGA_FILE /= NOMEGA .OR. NOMEGAR_FILE /= NOMEGAR .OR. ABS(OMEGAMAX_FILE-OMEGAMAX)>1E-3_q .OR. ABS(OMEGATL_FILE-OMEGATL)>1E-3_q) THEN
          IERR=1
       ENDIF
       IF (IERR==0) THEN
          READ(IU,*) N_FILE
          IF (SUM(ABS(N_FILE-N)) /= 0) THEN
             WRITE(*,*) 'error in READ_SCREENED_2E_FILE:  TWOEINT'//APP//' not compatible'
             CALL M_exit(); stop
          ENDIF
       ENDIF

       IF (IERR==0) THEN
          READ(IU,'(10E14.7)')  SCREENED_TWO_ELECTRON_INTEGRALS 
       ENDIF
       IF (IERR==0 .AND. IU0>=0) THEN
          WRITE(IU0,*) 'the TWOEINT file was read successfully'
       ENDIF
       CLOSE(IU)
    ELSE
       CLOSE(IU)
    ENDIF

    IF (IERR/=0 .AND. IU0>=0) THEN
       WRITE(IU0,*) 'the TWOEINT file was not read'
    ENDIF
  END SUBROUTINE READ_SCREENED_2E_FILE


  SUBROUTINE READ_SCREENED_2E_FILE_HEAD( IERR, IU0, NOMEGA, NOMEGAR, OMEGAMAX, OMEGATL, NODE)
    INTEGER :: IERR, IU0        
    INTEGER, OPTIONAL :: NODE   ! append this integer to TWOEINT file name
! local
    INTEGER :: IU=72
    CHARACTER (2) :: APP
    INTEGER :: NOMEGA, NOMEGAR, NOMEGA_FILE, NOMEGAR_FILE
    REAL(q) :: OMEGAMAX, OMEGATL, OMEGAMAX_FILE, OMEGATL_FILE

    APP="  "
    IF (PRESENT(NODE)) THEN
       WRITE (APP  , "(I1,I1)") MOD(node/10,10),MOD(node,10)
       OPEN( UNIT=IU, FILE="TWOEINT"//APP, STATUS="OLD", IOSTAT=IERR)
    ELSE
       OPEN( UNIT=IU, FILE="TWOEINT", STATUS="OLD", IOSTAT=IERR)
    ENDIF

    IF (IERR==0) THEN
       READ(IU,'(2I10,2E14.7)', IOSTAT=IERR) NOMEGA_FILE, NOMEGAR_FILE, OMEGAMAX_FILE, OMEGATL_FILE
       IF (IERR==0) THEN
          NOMEGA=NOMEGA_FILE
          NOMEGAR=NOMEGAR_FILE
          OMEGAMAX=OMEGAMAX_FILE
          OMEGATL =OMEGATL_FILE
       ENDIF
       CLOSE(IU)
    ELSE
       CLOSE(IU)
    ENDIF

    IF (IERR/=0 .AND. IU0>=0) THEN
       WRITE(IU0,*) 'the head of TWOEINT'//APP//' was not read sucessfully'
    ELSE IF (IU0>=0) THEN
       WRITE(IU0,*) 'read OMEGAMAX, OMEGATL, NOMEGA and NOMEGAR from TWOEINT'//APP
    ENDIF
  END SUBROUTINE READ_SCREENED_2E_FILE_HEAD


  SUBROUTINE WRITE_SCREENED_2E_FILE(SCREENED_TWO_ELECTRON_INTEGRALS, IU6, & 
      NOMEGA, NOMEGAR, OMEGAMAX, OMEGATL, NODE)
    COMPLEX(qs) :: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,:)
    INTEGER :: IU6
    INTEGER, OPTIONAL :: NODE   ! append this integer to TWOEINT file name
! local
    INTEGER :: N(6)
    INTEGER :: IU=72, IERR
    CHARACTER (2) :: APP
    INTEGER :: NOMEGA, NOMEGAR
    REAL(q) :: OMEGAMAX, OMEGATL

    N=SHAPE(SCREENED_TWO_ELECTRON_INTEGRALS)

    APP="  "
    IF (PRESENT(NODE)) THEN
       WRITE (APP  , "(I1,I1)") MOD(NODE/10,10),MOD(NODE,10)
       OPEN( UNIT=IU, FILE="TWOEINT"//APP, IOSTAT=IERR)
    ELSE
       OPEN( UNIT=IU, FILE="TWOEINT", IOSTAT=IERR)
    ENDIF
    IF (IU6>=0) WRITE(IU6,*) 'writing TWOEINT'//APP//' file'

    IF (IERR==0) THEN
       WRITE(IU,'(2I10,2E14.7)') NOMEGA, NOMEGAR, OMEGAMAX, OMEGATL
       WRITE(IU,'(7I10)') N
       WRITE(IU,'(10E14.7)')  SCREENED_TWO_ELECTRON_INTEGRALS 
       IF (IU6>=0) WRITE(IU6,*) 'TWOEINT file was written successfully'
       CLOSE( UNIT=IU)
    ELSE
       IF (IU6>=0) WRITE(IU6,*) 'TWOEINT file was not written'
    ENDIF

  END SUBROUTINE WRITE_SCREENED_2E_FILE

!*******************  SHIFT_CELEN *************************************
!
! determine the valence band maximum and the conduction band
! minimum and shift the eigenenergies such that the Fermi-level
! lies exactly at 0._q
!
!**********************************************************************

  SUBROUTINE SHIFT_CELEN( W, ESHIFT)
    USE wave
    TYPE (wavespin)     W
    REAL(q), OPTIONAL :: ESHIFT

    INTEGER :: ISP, NK, N
    REAL(q) :: VBM, CBM, EFERMI
   
    VBM=-1000
    CBM= 1000
    DO ISP=1, W%WDES%ISPIN
       DO NK=1, W%WDES%NKPTS
          DO N=1,W%WDES%NB_TOT

! filled band
             IF ( W%FERTOT(N, NK, ISP)>0.5) THEN
                VBM=MAX(VBM, REAL( W%CELTOT(N, NK, ISP), q))
             ELSE
                CBM=MIN(CBM, REAL( W%CELTOT(N, NK, ISP), q))
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    EFERMI= (VBM+CBM)/2

    W%CELTOT= W%CELTOT-EFERMI
    IF (PRESENT(ESHIFT)) THEN
       ESHIFT=EFERMI
    ENDIF
    
  END SUBROUTINE SHIFT_CELEN

!*******************  MEAN_CBM_VBM  ************************************
!
! determine the mean eigenvalue of the last state with an occupancy
! larger 0.5 (valence band or occupied)
! and the first state with an occupancy smaller 0.5 (conduction band
! or virtual state)
! the routine also introduces a weight measuring the deviation from
! the integer occupancy
! if there is only 1._q orbital with fraction occupancy, the
! Fermi-level will be placed slightly above or below that orbital
! depending on the occupancy
!
! the second routine checks the supplied Fermi-level and
! ajusts the Fermi-level if it lies below the CBM or above
! the VBM
!
!**********************************************************************

  SUBROUTINE MEAN_CBM_VBM( W, EFERMI)
    USE wave
    TYPE (wavespin)     W
    REAL(q) :: EFERMI
! local
    INTEGER :: ISP, NK, N
    REAL(q) :: VBM, CBM, VBM_WEIGHT, CBM_WEIGHT
   
    VBM=-1000
    CBM= 1000
    VBM_WEIGHT=0
    CBM_WEIGHT=0

    DO ISP=1, W%WDES%ISPIN
       DO NK=1, W%WDES%NKPTS
          DO N=1,W%WDES%NB_TOT

! filled band
             IF ( W%FERTOT(N, NK, ISP)>0.5) THEN
! new VBM
                IF (REAL( W%CELTOT(N, NK, ISP), q)> VBM) THEN
                   VBM=REAL( W%CELTOT(N, NK, ISP), q)
! weight for this state
! weight is larger for fractionally filled state
! if the occupancy is 1, the state is weighted by 0.01
                   VBM_WEIGHT=MAX((1-W%FERTOT(N, NK, ISP)),0.00001_q)
                ENDIF
             ELSE
! new CBM
                IF (REAL( W%CELTOT(N, NK, ISP), q)< CBM) THEN
                   CBM=REAL( W%CELTOT(N, NK, ISP), q)
! weight for this state
! weight is larger for fractionally filled state
! if the occupancy is 0, the state is weighted by 0.01
                   CBM_WEIGHT=MAX((W%FERTOT(N, NK, ISP)),0.00001_q)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    EFERMI= (VBM*VBM_WEIGHT+CBM* CBM_WEIGHT)/(VBM_WEIGHT+CBM_WEIGHT)


  END SUBROUTINE MEAN_CBM_VBM

!
! if NUP_DOWN is set in the INCAR, the Fermi-level
! is determined independently for up and down electrons
!

  SUBROUTINE MEAN_CBM_VBM_SPIN( W, EFERMI, NUP_DOWN)
    USE wave
    TYPE (wavespin)     W
    REAL(q) :: EFERMI(W%WDES%ISPIN)
    REAL(q) :: NUP_DOWN 
! local
    INTEGER :: ISP, NK, N
    REAL(q) :: VBM, CBM, VBM_WEIGHT, CBM_WEIGHT

    IF (NUP_DOWN >= 0 .AND. W%WDES%ISPIN>1 ) THEN
       DO ISP=1, W%WDES%ISPIN
          VBM=-1000
          CBM= 1000
          VBM_WEIGHT=0
          CBM_WEIGHT=0
          

          DO NK=1, W%WDES%NKPTS
             DO N=1,W%WDES%NB_TOT
! filled band
                IF ( W%FERTOT(N, NK, ISP)>0.5) THEN
! new VBM
                   IF (REAL( W%CELTOT(N, NK, ISP), q)> VBM) THEN
                      VBM=REAL( W%CELTOT(N, NK, ISP), q)
! weight for this state
! weight is larger for fractionally filled state
! if the occupancy is 1, the state is weighted by 0.01
                      VBM_WEIGHT=MAX((1-W%FERTOT(N, NK, ISP)),0.00001_q)
                   ENDIF
                ELSE
! new CBM
                   IF (REAL( W%CELTOT(N, NK, ISP), q)< CBM) THEN
                      CBM=REAL( W%CELTOT(N, NK, ISP), q)
! weight for this state
! weight is larger for fractionally filled state
! if the occupancy is 0, the state is weighted by 0.01
                      CBM_WEIGHT=MAX((W%FERTOT(N, NK, ISP)),0.00001_q)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          EFERMI(ISP)= (VBM*VBM_WEIGHT+CBM* CBM_WEIGHT)/(VBM_WEIGHT+CBM_WEIGHT) 
      ENDDO
    ELSE
! call standard routine
       CALL MEAN_CBM_VBM( W, EFERMI(1))
       IF (W%WDES%ISPIN>1) EFERMI(2)=EFERMI(1)
    ENDIF

  END SUBROUTINE MEAN_CBM_VBM_SPIN


END MODULE screened_2e
