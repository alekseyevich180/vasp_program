# 1 "pot.F"
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

# 2 "pot.F" 2 
      MODULE pot
      USE prec
      USE charge
      CONTAINS
!************************ SUBROUTINE POTLOK ****************************
! RCS:  $Id: pot.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!
! this subroutine calculates  the total local potential CVTOT
! which is the sum of the hartree potential, the exchange-correlation
! potential and the ionic local potential
! the routine also calculates the total local potential SV on the small
! grid
! on entry:
!  CHTOT(:,1)    density
!  CHTOT(:,2)    respectively CHTOT(:,2:4) contain the magnetization
! on return (LNONCOLLINEAR=.FALSE.):
!  CVTOT(:,1)    potential for up
!  CVTOT(:,2)    potential for down
! on return (LNONCOLLINEAR=.TRUE.):
!  CVTOT(:,1:4)  spinor representation of potential
!
!***********************************************************************

    SUBROUTINE POTLOK(GRID,GRIDC,GRID_SOFT, COMM_INTER, WDES,  &
                  INFO,P,T_INFO,E,LATT_CUR,  &
                  CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF )
      USE prec
      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE setexm
      USE base
      USE xcgrad
      USE wave
      USE mdipol
      USE meta
      USE Constrained_M_modular
      USE main_mpi, ONLY: COMM
! solvation__
      USE solvation
! solvation__

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID,GRIDC,GRID_SOFT
      TYPE (wavedes)     WDES
      TYPE (transit)     SOFT_TO_C
      TYPE (info_struct) INFO
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (energy)      E
      TYPE (latt)        LATT_CUR
      TYPE (communic)    COMM_INTER

      REAL(q)   SV(GRID%MPLWV*2, WDES%NCDIJ)
      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP), &
                 CHTOT(GRIDC%MPLWV, WDES%NCDIJ), CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      REAL(q)      DENCOR(GRIDC%RL%NP)
      REAL(q)    XCSIF(3,3),TMPSIF(3,3)
! work arrays (allocated after call to FEXCG)
      COMPLEX(q), ALLOCATABLE::  CWORK1(:),CWORK(:,:)
      REAL(q) ELECTROSTATIC
      LOGICAL, EXTERNAL :: L_NO_LSDA_GLOBAL

# 76

      
      MWORK1=MAX(GRIDC%MPLWV,GRID_SOFT%MPLWV)
      ALLOCATE(CWORK1(MWORK1),CWORK(GRIDC%MPLWV,WDES%NCDIJ))

!-----------------------------------------------------------------------
!
!  calculate the exchange correlation potential and the dc. correction
!
!-----------------------------------------------------------------------
      EXC     =0
      E%XCENC =0
      E%EXCG  =0
      E%CVZERO=0
      XCSIF   =0

      CVTOT=0
      
# 104


! test
! xc: IF (ISLDAXC()) THEN
  xc: IF (ISLDAXC().AND.(.NOT.LDO_METAGGA())) THEN
! test
! transform the charge density to real space
        EXCG  =0
        XCENCG=0
        CVZERG=0
        TMPSIF=0

        DO ISP=1,WDES%NCDIJ
           CALL FFT3D_MPI(CHTOT(1,ISP),GRIDC,1)
        ENDDO

# 130


        IF (WDES%ISPIN==2) THEN

! get the charge and the total magnetization
          CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)
! do LDA+U instead of LSDA+U
          IF (L_NO_LSDA_GLOBAL()) CWORK(:,2)=0
!
          IF (ISGGA()) THEN
! gradient corrections to LDA
! unfortunately FEXCGS requires (up,down) density
! instead of (rho,mag)
             CALL RL_FLIP(CWORK, GRIDC, 2, .TRUE.)
! GGA potential
             CALL FEXCGS(2, GRIDC, LATT_CUR, XCENCG, EXCG, CVZERG, TMPSIF, &
                  CWORK, CVTOT, DENCOR)
             CALL RL_FLIP(CWORK, GRIDC, 2, .FALSE.)
          ENDIF

! add LDA part of potential
          CALL FEXCF(GRIDC,LATT_CUR%OMEGA, &
             CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
             E%CVZERO,EXC,E%XCENC,XCSIF, .TRUE.)
!gk COH
! add Coulomb hole
          CALL COHSM1_RGRID(2, CWORK(1,1), CVTOT(1,1), DENCOR, GRIDC, LATT_CUR%OMEGA, .TRUE.)
!gK COHend
! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)

! get the proper direction vx = v0 + hat m delta v
          CALL MAG_DIRECTION(CHTOT(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
        ELSEIF (WDES%LNONCOLLINEAR) THEN
          IF (ISGGA()) THEN
! GGA potential
             CALL FEXCGS(4, GRIDC, LATT_CUR, XCENCG, EXCG, CVZERG, TMPSIF, &
                  CHTOT, CVTOT, DENCOR)
          ENDIF

! FEXCF requires (up,down) density instead of (rho,mag)
          CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)
! quick hack to do LDA+U instead of LSDA+U
          IF (L_NO_LSDA_GLOBAL()) CWORK(:,2)=0
! end of hack
! add LDA part of potential
          CALL FEXCF(GRIDC,LATT_CUR%OMEGA, &
             CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
             E%CVZERO,EXC,E%XCENC,XCSIF, .TRUE.)
!gk COH
! add Coulomb hole
          CALL COHSM1_RGRID(2, CWORK(1,1), CVTOT(1,1), DENCOR, GRIDC, LATT_CUR%OMEGA, .TRUE.)
!gK COHend
! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)
! get the proper direction vx = v0 + hat m delta v
                    
          CALL MAG_DIRECTION(CHTOT(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
       ELSE
          IF (ISGGA()) THEN
! gradient corrections to LDA
             CALL FEXCG(GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,TMPSIF, &
                  CHTOT,CVTOT,DENCOR)
          ENDIF
                
! LDA part of potential
          CALL FEXCP(GRIDC,LATT_CUR%OMEGA, &
               CHTOT,DENCOR,CVTOT,CWORK,E%CVZERO,EXC,E%XCENC,XCSIF,.TRUE.)
!gk COH
! add Coulomb hole
          CALL COHSM1_RGRID(1, CHTOT(1,1), CVTOT(1,1), DENCOR, GRIDC, LATT_CUR%OMEGA, .TRUE.)
!gK COHend
       ENDIF

# 248


       XCSIF=XCSIF+TMPSIF
       E%EXCG=EXC+EXCG
       E%XCENC=E%XCENC+XCENCG
       E%CVZERO=E%CVZERO+CVZERG

      ELSE xc
         DO ISP=1,WDES%NCDIJ
            CALL FFT3D_MPI(CHTOT(1,ISP),GRIDC,1)
         ENDDO
      ENDIF xc
!-MM- changes to accomodate constrained moments
!-----------------------------------------------------------------------
! add constraining potential
!-----------------------------------------------------------------------

# 272


!-MM- end of addition

!-----------------------------------------------------------------------
! calculate the total potential
!-----------------------------------------------------------------------
! add external electrostatic potential
      DIP%ECORR=0
      DIP%E_ION_EXTERN=0

      IF (DIP%LCOR_DIP) THEN
! get the total charge and store it in CWORK
          IF  ( WDES%NCDIJ > 1) THEN
             CALL MAG_DENSITY(CHTOT,CWORK, GRIDC, WDES%NCDIJ)
          ELSE
             CALL RL_ADD(CHTOT,1.0_q,CHTOT,0.0_q,CWORK,GRIDC)
          ENDIF

           CALL CDIPOL(GRIDC, LATT_CUR,P,T_INFO, &
             CWORK,CSTRF,CVTOT(1,1), WDES%NCDIJ, INFO%NELECT )

         CALL EXTERNAL_POT(GRIDC, LATT_CUR, CVTOT(1,1))
      ENDIF

      DO ISP=1,WDES%NCDIJ
         CALL FFT_RC_SCALE(CHTOT(1,ISP),CHTOT(1,ISP),GRIDC)
         CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
      ENDDO
!-----------------------------------------------------------------------
! FFT of the exchange-correlation potential to reciprocal space
!-----------------------------------------------------------------------
      RINPL=1._q/GRIDC%NPLWV
      DO  ISP=1,WDES%NCDIJ 
         CALL RL_ADD(CVTOT(1,ISP),RINPL,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
         CALL FFT3D_MPI(CVTOT(1,ISP),GRIDC,-1)
      ENDDO
!-----------------------------------------------------------------------
! add the hartree potential and the double counting corrections
!-----------------------------------------------------------------------
      CALL POTHAR(GRIDC, LATT_CUR, CHTOT, CWORK,E%DENC)
      DO I=1,GRIDC%RC%NP
         CVTOT(I,1)=CVTOT(I,1)+CWORK(I,1)
      ENDDO
! solvation__
!-----------------------------------------------------------------------
! add the dielectric corrections to CVTOT and the energy
!-----------------------------------------------------------------------
      CALL SOL_Vcorrection(INFO,T_INFO,LATT_CUR,P,WDES,GRIDC,CHTOT,CVTOT)
! solvation__
!-----------------------------------------------------------------------
!  add local pseudopotential potential
!-----------------------------------------------------------------------
      IF(INFO%TURBO==0)THEN
         CALL POTION(GRIDC,P,LATT_CUR,T_INFO,CWORK,CWORK1,CSTRF,E%PSCENC)
      ELSE
         CALL POTION_PARTICLE_MESH(GRIDC,P,LATT_CUR,T_INFO,CWORK,E%PSCENC,E%TEWEN)
      ENDIF

      ELECTROSTATIC=0
      NG=1
      col: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW
        FACTM=1
        IF (N3 /= 1) FACTM=2

        ELECTROSTATIC=ELECTROSTATIC+ FACTM* CWORK(NG,1)*CONJG(CHTOT(NG,1))
        NG=NG+1
      ENDDO row
      ENDDO col
      ELECTROSTATIC=ELECTROSTATIC+E%PSCENC-E%DENC+E%TEWEN

      E%PSCENC=E%PSCENC + DIP%ECORR + DIP%E_ION_EXTERN

      DO I=1,GRIDC%RC%NP
         CVTOT(I,1)=CVTOT(I,1)+CWORK(I,1)
      ENDDO
      CALL POT_FLIP(CVTOT, GRIDC,WDES%NCDIJ )
!=======================================================================
! if overlap is used :
! copy CVTOT to SV and set contribution of unbalanced lattice-vectors
! to (0._q,0._q),  then  FFT of SV and CVTOT to real space
!=======================================================================

      DO ISP=1,WDES%NCDIJ
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
         CALL CP_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CVTOT(1,ISP),CWORK1)
         CALL SETUNB(CWORK1,GRID_SOFT)
         CALL FFT3D_MPI(CWORK1,GRID_SOFT, 1)
         CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,SV(1,ISP),GRID_SOFT)

!  final result is only correct for first in-band-group
! (i.e. proc with nodeid 1 in COMM_INTER)
!  copy to other in-band-groups using COMM_INTER
! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)

         CALL M_bcast_d(COMM_INTER, SV(1,ISP), GRID%RL%NP)
# 373

         CALL FFT3D_MPI(CVTOT(1,ISP),GRIDC,1)
      ENDDO

      DEALLOCATE(CWORK1,CWORK)
      RETURN
    END SUBROUTINE POTLOK

!***********************************************************************
!
! small helper routine to set the local potential on the coarse
! plane wave grid from the full dense (augmentation) grid
! CVTOT must be supplied in reciprocal space (not usually the case)
! and is returned in real space
!
!***********************************************************************


    SUBROUTINE SET_SV( GRID, GRIDC, GRID_SOFT, COMM_INTER, SOFT_TO_C, NCDIJ, SV, CVTOT)

      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT NONE

      INTEGER NCDIJ
      TYPE (grid_3d)     GRID,GRIDC,GRID_SOFT
      TYPE (transit)     SOFT_TO_C

      REAL(q)   SV(GRID%MPLWV*2, NCDIJ)
      COMPLEX(q) CVTOT(GRIDC%MPLWV, NCDIJ)
      TYPE (communic)    COMM_INTER
! work arrays
      COMPLEX(q) ::  CWORK1(GRID_SOFT%MPLWV)
      INTEGER ISP


      DO ISP=1,NCDIJ
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
         CALL CP_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CVTOT(1,ISP),CWORK1)
         CALL SETUNB(CWORK1,GRID_SOFT)
         CALL FFT3D_MPI(CWORK1,GRID_SOFT, 1)
         CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,SV(1,ISP),GRID_SOFT)

!  final result is only correct for first in-band-group
! (i.e. proc with nodeid 1 in COMM_INTER)
!  copy to other in-band-groups using COMM_INTER
! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)

         CALL M_bcast_d(COMM_INTER, SV(1,ISP), GRID%RL%NP)
# 425

         CALL FFT3D_MPI(CVTOT(1,ISP),GRIDC,1)
      ENDDO
    END SUBROUTINE SET_SV

!************************ SUBROUTINE POTXC  ****************************
!
! this subroutine to calculate the XC-potential including gradient
! corrections
! this routine is required to calculate the partial core
! corrections to the forces
! on entry:
!  CHTOT(:,1)    density
!  CHTOT(:,2)    respectively CHTOT(:,2:4) contain the magnetization
! on return (LNONCOLLINEAR=.FALSE.):
!  CVTOT(:,1)    average potential
!  CVTOT(:,2:4)  magnetic field
!
!***********************************************************************

      SUBROUTINE POTXC(GRIDC, INFO, WDES, LATT_CUR, CVTOT,CHTOT,DENCOR)
      USE prec

      USE xcgrad
      USE setexm
      USE mpimy
      USE mgrid
      USE lattice
      USE base
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (wavedes)     WDES
      TYPE (info_struct) INFO
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CHTOT(GRIDC%MPLWV,WDES%NCDIJ),CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      REAL(q)      DENCOR(GRIDC%RL%NP)
! work arrays
      REAL(q)    XCSIF(3,3)
      COMPLEX(q), ALLOCATABLE:: CWORK(:,:)

      CVTOT = 0 
      ALLOCATE(CWORK(GRIDC%MPLWV,WDES%NCDIJ))

      DO ISP=1,WDES%NCDIJ
         CALL FFT3D_MPI(CHTOT(1,ISP),GRIDC,1)
      ENDDO
        IF (WDES%ISPIN==2) THEN

! get the charge and the total magnetization
          CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)

          IF (ISGGA()) THEN
! gradient corrections to LDA
! unfortunately FEXCGS requires (up,down) density
! instead of (rho,mag)
             CALL RL_FLIP(CWORK, GRIDC, 2, .TRUE.)
! GGA potential
             CALL FEXCGS(2, GRIDC, LATT_CUR, XCENCG, EXCG, CVZERG, XCSIF, &
                  CWORK, CVTOT, DENCOR)
             CALL RL_FLIP(CWORK, GRIDC, 2, .FALSE.)
          ENDIF

! add LDA part of potential
          CALL FEXCF(GRIDC,LATT_CUR%OMEGA, &
             CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
             CVZERO,EXC,XCENC,XCSIF, .TRUE.)
! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)

! get the proper direction vx = v0 + hat m delta v
          CALL MAG_DIRECTION(CHTOT(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
        ELSEIF (WDES%LNONCOLLINEAR) THEN
!-MM- gradient corrections in the noncollinear case are calculated
!     a bit differently than in the collinear case
          IF (ISGGA()) THEN
! GGA potential
             CALL FEXCGS(4, GRIDC, LATT_CUR, XCENCG, EXCG, CVZERG, XCSIF, &
                  CHTOT, CVTOT, DENCOR)
          ENDIF

! FEXCF requires (up,down) density instead of (rho,mag)
          CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)
! add LDA part of potential
          CALL FEXCF(GRIDC,LATT_CUR%OMEGA, &
             CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
             CVZERO,EXC,XCENC,XCSIF, .TRUE.)
! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)
! get the proper direction vx = v0 + hat m delta v
          CALL MAG_DIRECTION(CHTOT(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
!-MM- end of changes to calculation of gga in noncollinear case

       ELSE
          IF (ISGGA()) THEN
! gradient corrections to LDA
             CALL FEXCG(GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,XCSIF, &
                  CHTOT,CVTOT,DENCOR)
          ENDIF

! LDA part of potential
          CALL FEXCP(GRIDC,LATT_CUR%OMEGA, &
               CHTOT,DENCOR,CVTOT,CWORK,CVZERO,EXC,XCENC,XCSIF,.TRUE.)
       ENDIF

       DO ISP=1,WDES%NCDIJ
          CALL FFT_RC_SCALE(CHTOT(1,ISP),CHTOT(1,ISP),GRIDC)
          CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
       ENDDO
!-----------------------------------------------------------------------
! FFT of the exchange-correlation potential to reciprocal space
!-----------------------------------------------------------------------
      RINPL=1._q/GRIDC%NPLWV
      DO  ISP=1,WDES%NCDIJ 
         CALL RL_ADD(CVTOT(1,ISP),RINPL,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
         CALL FFT3D_MPI(CVTOT(1,ISP),GRIDC,-1)
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
      ENDDO
      DEALLOCATE(CWORK)

      END SUBROUTINE POTXC

!tb beg
!c this is basically the subroutine POTXC modified to work with the dDsC correction
!c scheme
      SUBROUTINE POTXC_ddsc(GRIDC, INFO, WDES, LATT_CUR, CVTOT,CHTOT,DENCOR,IVDW,XDM)
      USE prec

      USE xcgrad
      USE setexm
      USE mpimy
      USE mgrid
      USE lattice
      USE base
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (wavedes)     WDES
      TYPE (info_struct) INFO
      TYPE (latt)        LATT_CUR
      INTEGER IVDW
      COMPLEX(q) CHTOT(GRIDC%MPLWV,WDES%NCDIJ),CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      REAL(q)      DENCOR(GRIDC%RL%NP)
! work arrays
      REAL(q)    XCSIF(3,3),XDM(GRIDC%RL%NP)
      COMPLEX(q), ALLOCATABLE:: CWORK(:,:)

      CVTOT = 0
      ALLOCATE(CWORK(GRIDC%MPLWV,WDES%NCDIJ))

      DO ISP=1,WDES%NCDIJ
         CALL FFT3D_MPI(CHTOT(1,ISP),GRIDC,1)
      ENDDO
        IF (WDES%ISPIN==2) THEN

! get the charge and the total magnetization
          CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)
          IF (ISGGA()) THEN
! gradient corrections to LDA
! unfortunately FEXCGS requires (up,down) density
! instead of (rho,mag)
             CALL RL_FLIP(CWORK, GRIDC, 2, .TRUE.)
! GGA potential
             CALL FEXCGS_ddsc(2, GRIDC, LATT_CUR, XCENCG, EXCG, CVZERG, XCSIF, &
                  CWORK, CVTOT, DENCOR,XDM,IVDW)
             CALL RL_FLIP(CWORK, GRIDC, 2, .FALSE.)
          ENDIF

! add LDA part of potential
          CALL FEXCF(GRIDC,LATT_CUR%OMEGA, &
             CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
             CVZERO,EXC,XCENC,XCSIF, .TRUE.)
! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)

! get the proper direction vx = v0 + hat m delta v
          CALL MAG_DIRECTION(CHTOT(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
        ELSEIF (WDES%LNONCOLLINEAR) THEN
!-MM- gradient corrections in the noncollinear case are calculated
!     a bit differently than in the collinear case
          IF (ISGGA()) THEN
! GGA potential
             CALL FEXCGS_ddsc(4, GRIDC, LATT_CUR, XCENCG, EXCG, CVZERG, XCSIF, &
                  CHTOT, CVTOT, DENCOR,XDM,IVDW)
          ENDIF

! FEXCF requires (up,down) density instead of (rho,mag)
          CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)
! add LDA part of potential
          CALL FEXCF(GRIDC,LATT_CUR%OMEGA, &
             CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
             CVZERO,EXC,XCENC,XCSIF, .TRUE.)
! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)
! get the proper direction vx = v0 + hat m delta v
          CALL MAG_DIRECTION(CHTOT(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
!-MM- end of changes to calculation of gga in noncollinear case

       ELSE
          IF (ISGGA()) THEN
! gradient corrections to LDA
             CALL FEXCG_ddsc(GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,XCSIF, &
                  CHTOT,CVTOT,DENCOR,XDM,IVDW)
          ENDIF

! LDA part of potential
          CALL FEXCP(GRIDC,LATT_CUR%OMEGA, &
               CHTOT,DENCOR,CVTOT,CWORK,CVZERO,EXC,XCENC,XCSIF,.TRUE.)
       ENDIF

       DO ISP=1,WDES%NCDIJ
          CALL FFT_RC_SCALE(CHTOT(1,ISP),CHTOT(1,ISP),GRIDC)
          CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
       ENDDO
!-----------------------------------------------------------------------
! FFT of the exchange-correlation potential to reciprocal space
!-----------------------------------------------------------------------
      RINPL=1._q/GRIDC%NPLWV
      DO  ISP=1,WDES%NCDIJ
         CALL RL_ADD(CVTOT(1,ISP),RINPL,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
         CALL FFT3D_MPI(CVTOT(1,ISP),GRIDC,-1)
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
      ENDDO
      DEALLOCATE(CWORK)

      END SUBROUTINE POTXC_ddsc
!tb end

!************************ SUBROUTINE POTION_PARTICLE_MESH **************
!
!***********************************************************************

      SUBROUTINE POTION_PARTICLE_MESH(GRIDC,P, LATT_CUR, T_INFO, CVPS, PSCENC, TEWEN, FOR)
      USE ini
      USE prec
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE constant
      USE charge
      IMPLICIT NONE
      TYPE (grid_3d) GRIDC
      TYPE (type_info) T_INFO
      TYPE (potcar) P(T_INFO%NTYP)
      TYPE (latt) LATT_CUR
      COMPLEX(q) CVPS(GRIDC%RC%NP)
      REAL(q) PSCENC,TEWEN
      REAL(q), OPTIONAL :: FOR(:,:)
! local variables
      REAL(q) ZVSUM,G,GX,GY,GZ,TEWEN0
      INTEGER NT,NIS,NI,N,N1,N2,N3,NC,FACTM
      REAL(q), ALLOCATABLE :: QTOT(:)
      COMPLEX(q),ALLOCATABLE :: CHG(:)

      ZVSUM=0
      DO NT=1,T_INFO%NTYP
         ZVSUM=ZVSUM+P(NT)%ZVALF*T_INFO%NITYP(NT)*T_INFO%VCA(NT)
      ENDDO

      PSCENC=0
      DO NT=1,T_INFO%NTYP
         PSCENC=PSCENC+P(NT)%PSCORE*T_INFO%VCA(NT)*T_INFO%NITYP(NT)*(ZVSUM/LATT_CUR%OMEGA)
      ENDDO

      TEWEN=-PSCENC
      DO NT=1,T_INFO%NTYP
         TEWEN=TEWEN-P(NT)%ESELF*T_INFO%VCA(NT)*T_INFO%NITYP(NT)
      ENDDO

      CVPS =0._q
      ALLOCATE(CHG(GRIDC%MPLWV))
      CHG=0

! add pseudo chg.dens. and FFT to recipr.space
      DO NT=1,T_INFO%NTYP
         P(NT)%USESPL=>P(NT)%RHOSPL
         P(NT)%USEZ=-P(NT)%ZVALF
         P(NT)%USECUT=P(NT)%RCUTRHO
      END DO

      ALLOCATE(QTOT(T_INFO%NIONS))

      CALL RHOADD(T_INFO,LATT_CUR,P,GRIDC,CHG,QTOT)

! multiply by 4pi/G^2
! loop over recipr.lattice points
      N=0
      TEWEN0=0._q
      col: DO NC=1,GRIDC%RC%NCOL
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         row: DO N1=1,GRIDC%RC%NROW
            N=N+1
            FACTM=1
            IF (N3 /= 1) FACTM=2

!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
            GX= GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
            GY= GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
            GZ= GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)

            G=SQRT(GX**2+GY**2+GZ**2)*2*PI
            IF ( (GRIDC%LPCTX(N1)/=0) .OR. (GRIDC%LPCTY(N2)/=0) .OR. (GRIDC%LPCTZ(N3)/=0) ) THEN                
!              CVPS(N) = CVPS(N) + ((EDEPS/LATT_CUR%OMEGA)*CHG(N)/G**2)
               CVPS(N)=((EDEPS/LATT_CUR%OMEGA)*CHG(N)/G**2)
            ELSE
               CVPS(N)=0._q
            ENDIF
            TEWEN0=TEWEN0+FACTM* CVPS(N)*CONJG(CHG(N))*0.5_q

         ENDDO row
      ENDDO col
      CALL M_sum_d(GRIDC%COMM,TEWEN0,1)
      TEWEN=TEWEN+TEWEN0

      CALL SETUNB(CVPS,GRIDC)
      DEALLOCATE(CHG)

      IF (PRESENT(FOR)) THEN
! bring the potential to real space
         CALL FFT3D_MPI(CVPS(1),GRIDC,1)
! calculate force on each ion
!
! (dE/dR_i,a) = e^2 \int d^3 r (drho_i(|r-R_i|)/dR_i,a) V(r)
! i=1..NION
! a=1..3 (x,y,z)
!
! drho_i(|r-R_i|)/dR_i,a = -rho'(|r-R_i|) (r_a-R_i,a)/|r-R_i|
!
! Evaluate f_i(r)=d rho/dr*(r-R_i)/|r-R_i| on the grid for each ion
! and integrate F_i = \int dr f_i(r) V(r)
         NIS=1
         type: DO NT=1,T_INFO%NTYP
            IF (.NOT.ASSOCIATED(P(NT)%USESPL)) THEN
               NIS = NIS+T_INFO%NITYP(NT); CYCLE type
            ENDIF
            ions: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1
               CALL RHODER(T_INFO,LATT_CUR,P(NT),NI,GRIDC,QTOT,CVPS,FOR(1,NI))
            ENDDO ions
            NIS = NIS+T_INFO%NITYP(NT)
         ENDDO type
         FOR=FOR*LATT_CUR%OMEGA
! take the potential back to reciprocal space
         CALL RL_ADD(CVPS(1),1._q/GRIDC%NPLWV,CVPS(1),0.0_q,CVPS(1),GRIDC)
         CALL FFT3D_MPI(CVPS(1),GRIDC,-1)
         CALL SETUNB_COMPAT(CVPS(1),GRIDC)
      ENDIF

      DEALLOCATE(QTOT)

      RETURN
      END SUBROUTINE POTION_PARTICLE_MESH

      END MODULE

!************************ SUBROUTINE POTHAR ****************************
!
! this subroutine calculates the hartree potential from the electronic
! charge density. The correction to the
! total energy due to overcounting the hartree energy on summing the
! electronic eigenvalues is also computed (hartree contribution to the
! total energy = 0.5*sum (vh(g)*rho(-g)). the sum of eigenvalues gives
! sum (vh(g)*rho(-g)) where vh(g) is the hartree potential at wavevector
! g and rho(g) is the charge density at wavevector g)
!
!***********************************************************************

      SUBROUTINE POTHAR(GRIDC,LATT_CUR, CHTOT,CVD,DENC)
      USE prec
      USE mpimy
      USE mgrid
      USE lattice
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      COMPLEX(q) CVD(GRIDC%RC%NP),CHTOT(GRIDC%RC%NP)

      DENC=0._q
!=======================================================================
! scale the hartree potential by edeps divided by the volume of the unit
! cell
!=======================================================================
      SCALE=EDEPS/LATT_CUR%OMEGA/TPI**2
!=======================================================================
! calculate the hartree potential on the grid of reciprocal lattice
! vectors and the correction to the total energy
!=======================================================================
      NI=0
      col: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW

        NI=NI+1
        FACTM=1
        IF (N3 /= 1) FACTM=2

        GX= (GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
        GY= (GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
        GZ= (GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))

        GSQU=GX**2+GY**2+GZ**2
!=======================================================================
! since the G=0 coulomb contributions to the hartree, ewald and
! electron-ion energies are individually divergent but together sum to
! (0._q,0._q), set the hartree potential at G=0 to (0._q,0._q).
!=======================================================================
        IF ((GRIDC%LPCTX(N1)==0).AND.(GRIDC%LPCTY(N2)==0).AND.(GRIDC%LPCTZ(N3)==0)) &
     & THEN
          CVD(NI)=(0.0_q,0.0_q)
        ELSE
          CVD(NI)=CHTOT(NI)/GSQU*SCALE
        ENDIF
      ENDDO row
      ENDDO col
      CALL SETUNB(CVD,GRIDC)
!=======================================================================
! calculate the correction to the total energy
!=======================================================================
      NI=0
      col2: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row2: DO N1=1,GRIDC%RC%NROW

        NI=NI+1
        FACTM=1
        IF (N3 /= 1) FACTM=2

        DUM=FACTM* CVD(NI)*CONJG(CHTOT(NI))
        DENC=DENC+DUM
      ENDDO row2
      ENDDO col2
      DENC=-DENC/2
      CALL M_sum_d(GRIDC%COMM,DENC,1)

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE POTION ****************************
!
! this subroutine calculates the pseudopotential and its derivatives
! multiplied by the partial structur-factors
! on the grid of  reciprocal lattice vectors
!
!***********************************************************************

      SUBROUTINE POTION(GRIDC,P,LATT_CUR,T_INFO,CVPS,CDVPS,CSTRF,PSCENC)
      USE prec

      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      COMPLEX(q) CVPS(GRIDC%RC%NP),CDVPS(GRIDC%RC%NP)

!=======================================================================
! calculate the contribution to the total energy from the non-coulomb
! part of the g=0 component of the pseudopotential and the force on the
! unit cell due to the change in this energy as the size of the cell
! changes
!=======================================================================
      ZVSUM=0
      DO NT=1,T_INFO%NTYP
         ZVSUM=ZVSUM+P(NT)%ZVALF*T_INFO%NITYP(NT)*T_INFO%VCA(NT)
      ENDDO

      PSCENC=0
      DO NT=1,T_INFO%NTYP
         PSCENC=PSCENC+P(NT)%PSCORE*T_INFO%VCA(NT)*T_INFO%NITYP(NT)*(ZVSUM/LATT_CUR%OMEGA)
      ENDDO

      CVPS =0
      CDVPS=0
!=======================================================================
! loop over all types of atoms
! multiply structur factor by local pseudopotential
!=======================================================================
      typ: DO NT=1,T_INFO%NTYP

      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS
      ZZ=  -4*PI*P(NT)%ZVALF*FELECT

      N=0
      col: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW
        N=N+1
!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
        GX= GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
        GY= GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
        GZ= GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*2*PI
        IF ( ( (GRIDC%LPCTX(N1)/=0) .OR. (GRIDC%LPCTY(N2)/=0) .OR. &
     &         (GRIDC%LPCTZ(N3)/=0) ) .AND. (G<PSGMA2) ) THEN
!=======================================================================
! convert the magnitude of the reciprocal lattice vector to a position
! in the pseudopotential arrays and interpolate the pseudopotential and
! its derivative
!=======================================================================
        I  =INT(G*ARGSC)+1
        REM=G-P(NT)%PSP(I,1)
        VPST =(P(NT)%PSP(I,2)+REM*(P(NT)%PSP(I,3)+ &
     &                     REM*(P(NT)%PSP(I,4)  +REM*P(NT)%PSP(I,5))))
        DVPST= P(NT)%PSP(I,3)+REM*(P(NT)%PSP(I,4)*2+REM*P(NT)%PSP(I,5)*3)

        CVPS (N)=CVPS (N)+( VPST+ ZZ / G**2)  /LATT_CUR%OMEGA*CSTRF(N,NT)
        CDVPS(N)=CDVPS(N)+( DVPST- 2*ZZ /G**3)/LATT_CUR%OMEGA*CSTRF(N,NT)
        ELSE
        CVPS (N)=0._q
        CDVPS(N)=0._q
        ENDIF

      ENDDO row
      ENDDO col
      ENDDO typ
! set local potential to (0._q,0._q)
!      CVPS=0
!      CDVPS=0
      CALL SETUNB(CVPS,GRIDC)
      CALL SETUNB(CDVPS,GRIDC)

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE  MAG_DIRECTION  *******************
!
! on entry CVTOT must contain the v_xc(up) and v_xc(down)
! on return CVTOT contains
! collinear case:
!  CVTOT(:,1) =  (v_xc(up) + v_xc(down))/2
!  CVTOT(:,2) =  (v_xc(up) - v_xc(down))/2
! non collinear case:
!  CVTOT(:,1) =  (v_xc(up) + v_xc(down))/2
!  CVTOT(:,2) =  hat m_x (v_xc(up) - v_xc(down))/2
!  CVTOT(:,3) =  hat m_y (v_xc(up) - v_xc(down))/2
!  CVTOT(:,4) =  hat m_z (v_xc(up) - v_xc(down))/2
! where hat m is the unit vector of the local magnetization density
!
!***********************************************************************



      SUBROUTINE MAG_DIRECTION(CHTOT, CVTOT, GRID, NCDIJ)

      USE prec
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)     GRID
      INTEGER NCDIJ
      
      REAL(q) CHTOT(GRID%MPLWV*2, NCDIJ), &
            CVTOT(GRID%MPLWV*2, NCDIJ)
! local variables
      INTEGER K
      REAL(q) :: NORM2,DELTAV,V0
      
      IF (NCDIJ==2) THEN
         DO K=1,GRID%RL%NP
            V0    =(CVTOT(K,1)+CVTOT(K,2))/2
            DELTAV=(CVTOT(K,1)-CVTOT(K,2))/2
            CVTOT(K,1) = V0

            CVTOT(K,2) = DELTAV
# 1019

         ENDDO
      ELSE IF (NCDIJ==4) THEN
         DO K=1,GRID%RL%NP
            V0    =(CVTOT(K,1)+CVTOT(K,2))/2
            DELTAV=(CVTOT(K,1)-CVTOT(K,2))/2

            NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CHTOT(K,2)+ CHTOT(K,3)*CHTOT(K,3)+ CHTOT(K,4)*CHTOT(K,4))),1.E-20_q)
# 1030


            CVTOT(K,1) = V0
            CVTOT(K,2) = DELTAV * REAL(CHTOT(K,2),KIND=q) / NORM2
            CVTOT(K,3) = DELTAV * REAL(CHTOT(K,3),KIND=q) / NORM2
            CVTOT(K,4) = DELTAV * REAL(CHTOT(K,4),KIND=q) / NORM2
         ENDDO
         
      ELSE
         WRITE(*,*) 'internal error: MAG_DIRECTION called with NCDIJ=',NCDIJ
         CALL M_exit(); stop
      ENDIF
         
      END SUBROUTINE MAG_DIRECTION


      SUBROUTINE MAG_DIRECTION_KINDENS(CHTOT, TAU, MUTOT, CVTOT, GRID, NCDIJ, LATT_CUR)

      USE prec
      USE constant
      USE lattice
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      INTEGER NCDIJ
      
      REAL(q) CHTOT(GRID%MPLWV*2, NCDIJ), TAU(GRID%MPLWV*2, NCDIJ), &
            MUTOT(GRID%MPLWV*2, NCDIJ), CVTOT(GRID%MPLWV*2, NCDIJ)
! local variables
      INTEGER K
      REAL(q) :: NORM2,DELTAV,V0,TAUPROJ,FACT
!#define debug3
# 1067

      
      IF (NCDIJ==2) THEN
! local potential
         DO K=1,GRID%RL%NP
            V0    =(CVTOT(K,1)+CVTOT(K,2))/2
            DELTAV=(CVTOT(K,1)-CVTOT(K,2))/2
            CVTOT(K,1) = V0

            CVTOT(K,2) = DELTAV
# 1079

         ENDDO
! \mu = dE_xc / d\tau
         DO K=1,GRID%RL%NP
            V0    =(MUTOT(K,1)+MUTOT(K,2))/2
            DELTAV=(MUTOT(K,1)-MUTOT(K,2))/2
            MUTOT(K,1) = V0

            MUTOT(K,2) = DELTAV
# 1090

         ENDDO
      ELSE IF (NCDIJ==4) THEN
! local potential and \mu = dE_xc / d\tau
         DO K=1,GRID%RL%NP
            V0    =REAL((CVTOT(K,1)+CVTOT(K,2))/2,KIND=q)
            DELTAV=REAL((CVTOT(K,1)-CVTOT(K,2))/2,KIND=q)

            NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CHTOT(K,2)+ CHTOT(K,3)*CHTOT(K,3)+ CHTOT(K,4)*CHTOT(K,4))),1.E-20_q)
# 1104


            CVTOT(K,1) = V0
            CVTOT(K,2) = DELTAV * REAL(CHTOT(K,2),KIND=q) / NORM2
            CVTOT(K,3) = DELTAV * REAL(CHTOT(K,3),KIND=q) / NORM2
            CVTOT(K,4) = DELTAV * REAL(CHTOT(K,4),KIND=q) / NORM2
         ENDDO
! \mu = dE_xc / d\tau
         DO K=1,GRID%RL%NP
            V0    =REAL((MUTOT(K,1)+MUTOT(K,2))/2,KIND=q)
            DELTAV=REAL((MUTOT(K,1)-MUTOT(K,2))/2,KIND=q)

            NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CHTOT(K,2)+ CHTOT(K,3)*CHTOT(K,3)+ CHTOT(K,4)*CHTOT(K,4))),1.E-20_q)
# 1122


            MUTOT(K,1) = V0
            MUTOT(K,2) = DELTAV * REAL(CHTOT(K,2),KIND=q) / NORM2
            MUTOT(K,3) = DELTAV * REAL(CHTOT(K,3),KIND=q) / NORM2
            MUTOT(K,4) = DELTAV * REAL(CHTOT(K,4),KIND=q) / NORM2


            TAUPROJ = (TAU(K,2)*CHTOT(K,2)+TAU(K,3)*CHTOT(K,3)+TAU(K,4)*CHTOT(K,4)) /NORM2/NORM2
# 1135


            FACT=1._q/HSQDTM!*2._q
!           FACT=0._q
            IF (NORM2>1E-3_q) THEN
               CVTOT(K,2)=CVTOT(K,2)+DELTAV*REAL(TAU(K,2)-TAUPROJ*CHTOT(K,2),KIND=q)/NORM2*FACT
               CVTOT(K,3)=CVTOT(K,3)+DELTAV*REAL(TAU(K,3)-TAUPROJ*CHTOT(K,3),KIND=q)/NORM2*FACT
               CVTOT(K,4)=CVTOT(K,4)+DELTAV*REAL(TAU(K,4)-TAUPROJ*CHTOT(K,4),KIND=q)/NORM2*FACT
            ENDIF
# 1170

         ENDDO
# 1178

      ELSE
         WRITE(*,*) 'internal error: MAG_DIRECTION called with NCDIJ=',NCDIJ
         CALL M_exit(); stop
      ENDIF
 
      END SUBROUTINE MAG_DIRECTION_KINDENS


!************************ SUBROUTINE MAG_DENSITY ***********************
!
! this subroutine calculates the total charge density and the
! absolute magnitude of the magnetization density
! on entry:
!  CHTOT  rho, m_x, m_y, m_z
! on exit:
!  CWORK  rho, sqrt(m_x^2 + m_y^2 + m_z^2)
! in the collinear case, it this means a simple copy CHTOT to CWORK
!
!***********************************************************************

      SUBROUTINE MAG_DENSITY(CHTOT, CWORK, GRID, NCDIJ)

      USE prec
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)     GRID
      INTEGER NCDIJ
      
      REAL(q) CHTOT(GRID%MPLWV*2, NCDIJ), &
            CWORK(GRID%MPLWV*2, NCDIJ)
! local
      INTEGER K

      
      IF (NCDIJ==2) THEN
         DO K=1,GRID%RL%NP
            CWORK(K,1)=CHTOT(K,1)

            CWORK(K,2)=CHTOT(K,2)
# 1221

         ENDDO
      ELSE IF (NCDIJ==4) THEN
         DO K=1,GRID%RL%NP
            CWORK(K,1)=CHTOT(K,1)
            CWORK(K,2)=SQRT(ABS(CHTOT(K,2)*CHTOT(K,2)+ CHTOT(K,3)*CHTOT(K,3) + CHTOT(K,4)*CHTOT(K,4)))
         ENDDO
      ELSE
         WRITE(*,*) 'internal error: MAG_DENSITY called with NCDIJ=',NCDIJ
         CALL M_exit(); stop
      ENDIF
      END SUBROUTINE MAG_DENSITY


!************************ SUBROUTINE POT_FLIP **************************
!
!
! rearranges the storage mode for spin components of potentials:
! for the collinear case calculate:
!  v0 1 + v_z
!  v0 1 - v_z
! for the non collinear case calculate
!  v  = v0 1 + sigma_x v_x + simga_y v_y + sigma_z v_z
!
!***********************************************************************

      SUBROUTINE POT_FLIP(CVTOT, GRID, NCDIJ)
      USE prec
      USE mgrid
      IMPLICIT NONE
      INTEGER NCDIJ
      TYPE (grid_3d)     GRID
      COMPLEX(q) :: CVTOT(GRID%MPLWV, NCDIJ)

! local
      COMPLEX(q) :: C00,CX,CY,CZ
      REAL(q) :: FAC
      INTEGER K
      
      IF (NCDIJ==2) THEN
         FAC=1.0_q
         DO K=1,GRID%RC%NP
            C00=CVTOT(K,1)
            CZ =CVTOT(K,2)

            CVTOT(K,1)= (C00+CZ)*FAC           
            CVTOT(K,2)= (C00-CZ)*FAC           
         ENDDO
      ELSE IF (NCDIJ==4) THEN
         FAC=1.0_q
         DO K=1,GRID%RC%NP
            C00=CVTOT(K,1)
            CX =CVTOT(K,2)
            CY =CVTOT(K,3)
            CZ =CVTOT(K,4)

            CVTOT(K,1)= (C00+CZ)*FAC           
            CVTOT(K,2)= (CX-CY*(0._q,1._q))*FAC
            CVTOT(K,3)= (CX+CY*(0._q,1._q))*FAC
            CVTOT(K,4)= (C00-CZ)*FAC           
         ENDDO
      ELSE IF (NCDIJ==1) THEN
      ENDIF

    END SUBROUTINE POT_FLIP

    SUBROUTINE POT_FLIP_RL(CVTOT, GRID, NCDIJ)
      USE prec
      USE mgrid
      IMPLICIT NONE
      INTEGER NCDIJ
      TYPE (grid_3d)     GRID
      
      REAL(q) CVTOT(GRID%MPLWV*2, NCDIJ)
! local
      COMPLEX(q) :: C00,CX,CY,CZ
      REAL(q) :: FAC
      INTEGER K

      IF (NCDIJ==2) THEN
         FAC=1.0_q
         DO K=1,GRID%RL%NP
            C00=CVTOT(K,1)
            CZ =CVTOT(K,2)

            CVTOT(K,1)= (C00+CZ)*FAC
            CVTOT(K,2)= (C00-CZ)*FAC
         ENDDO
      ELSE IF (NCDIJ==4) THEN
         FAC=1.0_q
         DO K=1,GRID%RL%NP
            C00=CVTOT(K,1)
            CX =CVTOT(K,2)
            CY =CVTOT(K,3)
            CZ =CVTOT(K,4)

            CVTOT(K,1)= (C00+CZ)*FAC
            CVTOT(K,2)= (CX-CY*(0._q,1._q))*FAC
            CVTOT(K,3)= (CX+CY*(0._q,1._q))*FAC
            CVTOT(K,4)= (C00-CZ)*FAC
         ENDDO
      ELSE IF (NCDIJ==1) THEN
      ENDIF

    END SUBROUTINE POT_FLIP_RL


!************************ SUBROUTINE EXTERNAL_POT **********************
!
! this subroutine can be used to add an external potential
! the units of the potential are eV
!
!***********************************************************************

      SUBROUTINE EXTERNAL_POT(GRIDC, LATT_CUR, CVTOT)

      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      REAL(q)      CVTOT(GRIDC%MPLWV*2)

      RETURN

      NG=0

      IF (GRIDC%RL%NFAST==3) THEN
! mpi version: x-> N2, y-> N3, z-> N1
         N2MAX=GRIDC%NGX
         N3MAX=GRIDC%NGY
         N1MAX=GRIDC%NGZ

         DO NC=1,GRIDC%RL%NCOL
            N2= GRIDC%RL%I2(NC)
            N3= GRIDC%RL%I3(NC)
            DO N1=1,GRIDC%RL%NROW
               NG=NG+1
               CVTOT(NG)=CVTOT(NG)
            ENDDO
         ENDDO
      ELSE
! conventional version: x-> N1, y-> N2, z-> N3
         N1MAX=GRIDC%NGX
         N2MAX=GRIDC%NGY
         N3MAX=GRIDC%NGZ

         DO NC=1,GRIDC%RL%NCOL
            N2= GRIDC%RL%I2(NC)
            N3= GRIDC%RL%I3(NC)
            DO N1=1,GRIDC%RL%NROW
               NG=NG+1
               CVTOT(NG)=CVTOT(NG)
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE
