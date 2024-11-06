# 1 "xml_writer.F"
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

# 2 "xml_writer.F" 2 
!=======================================================================
!
! this routine writes out all the tags that are going to be used
! by vasp including defaults automactically set
! it should allow to redo calculations and might be important
! for some postprocessing programs
!
! after the tag SZNAM1 the calling interface is exactly equivalent to
! reader.F, even though some of the variables are not used in the
! xml_writer subroutine
!
!=======================================================================

      SUBROUTINE XML_WRITER( &
     &        NPAR, &
     &        SZNAM1,ISTART,IALGO,IMIX,MAXMIX,MREMOVE, &
     &        AMIX,BMIX,AMIX_MAG,BMIX_MAG,AMIN, &
     &        WC,INIMIX,MIXPRE,MIXFIRST,LFOUND,LDIAG,LSUBROT,LREAL,LREALD, &
     &        LPDENS,IBRION,ICHARG,INIWAV,NELM,NELMIN,NELMDL,EDIFF, &
     &        EDIFFG,NSW,ISIF,IWAVPR,ISYM,NBLOCK,KBLOCK,ENMAX,POTIM, &
     &        TEBEG,TEEND,NFREE, &
     &        NPACO,APACO,NTYPIN,NTYPD,SMASS,SCALEE,POMASS,DARWIN_V,DARWIN_R, &
     &        RWIGS,NELECT,NUP_DOWN,TIME,EMIN,EMAX,EFERMI,ISMEAR,SPACING,LGAMMA, & 
     &        PSTRESS,NDAV, &
     &        SIGMA,LTET,WEIMIN,EBREAK,DEPER,NWRITE,LCORR, &
     &        IDIOT,NIONS,NTYPP,lmusic,LOPTICS,STM, &
     &        ISPIN,ATOMOM,NIOND,LWAVE,LCHARG,LVTOT,LVHAR,SZPREC, &
     &        ENAUG,LORBIT,LELF,ROPT,ENINI, &
     &        NGX,NGY,NGZ,NGXF,NGYF,NGZF,NBANDS,NEDOS,NBLK,LATT_CUR, &
     &        LPLANE_WISE,LCOMPAT,LMAX_CALC,LMAX_MIX,NSIM,LPARD,LPAW,LADDGRID, &
     &        LNONCOLLINEAR,LSORBIT,SAXIS,LMETAGGA, &
     &        LSPIRAL,LZEROZ,QSPIRAL, &
     &        LASPH,LORBITALREAL, &
     &        TURBO,IRESTART,NREBOOT,NMIN,EREF, &
     &        NLSPLINE)

      USE prec
      USE sym_prec
      USE ini
      USE lattice
      USE scala
      USE wave_mpi
      USE constant
      USE pseudo   ! for subroutine EXTYP
      USE vaspxml
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (latt)        LATT_CUR

      CHARACTER (1)    CHARAC
      CHARACTER (40)   SZNAM1
      CHARACTER (80)   SZNAM
      CHARACTER (6)    SZPREC


      LOGICAL   LDUM,MIXFIRST,LFOUND,LDIAG,LSUBROT,LREAL,LREALD,LPDENS,LTET,LOPTICS, &
     &          LCORR,LOPEN,lmusic,LWAVE,LCHARG,LVTOT,LVHAR, &
     &          LORBIT_,LELF,LCOMPAT,LPARD,LPAW,LADDGRID, &
     &          LNONCOLLINEAR,LSORBIT,LMETAGGA,LPLANE_WISE,LASPH,LORBITALREAL
      INTEGER   IGGA2
      DIMENSION POMASS(NTYPD),RWIGS(NTYPP), &
     &          ROPT(NTYPD),DARWIN_V(NTYPD),DARWIN_R(NTYPD)
      DIMENSION ATOMOM(*)
      REAL(q)   SAXIS(3)
      REAL(q)   NELECT,NUP_DOWN
      REAL(q)   STM(7)
      INTEGER   TURBO,IRESTART,NREBOOT,NMIN
      REAL(q)   EREF
!-MM- Spin spiral stuff
      LOGICAL   LSPIRAL,LZEROZ
      REAL(q)   QSPIRAL(3)
      REAL(q)   SPACING
      LOGICAL   LGAMMA
      LOGICAL   NLSPLINE
!-MM- end of addition
      CALL XML_TAG("separator","general")
      CALL XML_INCAR('SYSTEM','S',IDUM,RDUM,CDUM,LDUM,SZNAM1,LEN(SZNAM1))
      CALL XML_INCAR('LCOMPAT','L',IDUM,RDUM,CDUM,LCOMPAT,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","electronic")
      CALL XML_INCAR('PREC','S',IDUM,RDUM,CDUM,LDUM,SZPREC,1)
      CALL XML_INCAR('ENMAX','F',IDUM,ENMAX,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('ENAUG','F',IDUM,ENAUG,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('EDIFF','F',IDUM,EDIFF,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('IALGO','I',IALGO,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('IWAVPR','I',IWAVPR,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NBANDS','I',NBANDS,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NELECT','F',IDUM,NELECT,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('TURBO','I',TURBO,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('IRESTART','I',IRESTART,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NREBOOT','I',NREBOOT,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NMIN','I',NMIN,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('EREF','F',IDUM,EREF,CDUM,LDUM,CHARAC,1)

      CALL XML_TAG("separator","electronic smearing")
      CALL XML_INCAR('ISMEAR','I',ISMEAR,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('SIGMA','F',IDUM,SIGMA,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('KSPACING','F',IDUM,SPACING,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('KGAMMA','L',IDUM,RDUM,CDUM,LGAMMA,CHARAC,1)
      CALL XML_CLOSE_TAG
      
      CALL XML_TAG("separator","electronic projectors")
      CALL XML_INCAR('LREAL','L',IDUM,RDUM,CDUM,LREAL,CHARAC,1)
      CALL XML_INCAR_V('ROPT','F',IDUM,ROPT,CDUM,LDUM,CHARAC,NTYPIN)
      CALL XML_INCAR('LMAXPAW','I',LMAX_CALC,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('LMAXMIX','I',LMAX_MIX,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NLSPLINE','L',IDUM,RDUM,CDUM,NLSPLINE,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","electronic startup")
      CALL XML_INCAR('ISTART','I',ISTART,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('ICHARG','I',ICHARG,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('INIWAV','I',INIWAV,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","electronic spin")
      CALL XML_INCAR('ISPIN','I',ISPIN,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('LNONCOLLINEAR','L',IDUM,RDUM,CDUM,LNONCOLLINEAR,CHARAC,1)
      NMAGMOM=NIONS
      IF (LNONCOLLINEAR) NMAGMOM=3*NIONS
      CALL XML_INCAR_V('MAGMOM','F',IDUM,ATOMOM,CDUM,LDUM,CHARAC,NMAGMOM)
      CALL XML_INCAR('NUPDOWN','F',IDUM,NUP_DOWN,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('LSORBIT','L',IDUM,RDUM,CDUM,LSORBIT,CHARAC,1)
      CALL XML_INCAR_V('SAXIS','F',IDUM,SAXIS,CDUM,LDUM,CHARAC,3)
      CALL XML_INCAR('LSPIRAL','L',IDUM,RDUM,CDUM,LSPIRAL,CHARAC,1)
      CALL XML_INCAR_V('QSPIRAL','F',IDUM,QSPIRAL,CDUM,LDUM,CHARAC,3)
      CALL XML_INCAR('LZEROZ','L',IDUM,RDUM,CDUM,LZEROZ,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","electronic exchange-correlation")
      CALL XML_INCAR('LASPH','L',IDUM,RDUM,CDUM,LASPH,CHARAC,1)
      CALL XML_INCAR('LMETAGGA','L',IDUM,RDUM,CDUM,LMETAGGA,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","electronic convergence")
      CALL XML_INCAR('NELM','I',NELM,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NELMDL','I',NELMDL,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NELMIN','I',NELMIN,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('ENINI','F',IDUM,ENINI,CDUM,LDUM,CHARAC,1)

      CALL XML_TAG("separator","electronic convergence detail")
      CALL XML_INCAR('LDIAG','L',IDUM,RDUM,CDUM,LDIAG,CHARAC,1)
      CALL XML_INCAR('LSUBROT','L',IDUM,RDUM,CDUM,LSUBROT,CHARAC,1)
      CALL XML_INCAR('WEIMIN','F',IDUM,WEIMIN,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('EBREAK','F',IDUM,EBREAK,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('DEPER','F',IDUM,DEPER,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NRMM','I',NDAV,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('TIME','F',IDUM,TIME,CDUM,LDUM,CHARAC,1)
      CALL XML_CLOSE_TAG
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","electronic mixer")
      CALL XML_INCAR('AMIX','F',IDUM,AMIX,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('BMIX','F',IDUM,BMIX,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('AMIN','F',IDUM,AMIN,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('AMIX_MAG','F',IDUM,AMIX_MAG,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('BMIX_MAG','F',IDUM,BMIX_MAG,CDUM,LDUM,CHARAC,1)

      CALL XML_TAG("separator","electronic mixer details")
      CALL XML_INCAR('IMIX','I',IMIX,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('MIXFIRST','L',IDUM,RDUM,CDUM,MIXFIRST,CHARAC,1)
      CALL XML_INCAR('MAXMIX','I',MAXMIX,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('WC','F',IDUM,WC,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('INIMIX','I',INIMIX,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('MIXPRE','I',MIXPRE,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('MREMOVE','I',MREMOVE,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_CLOSE_TAG
      CALL XML_CLOSE_TAG

      CALL XML_WRITE_EFIELD_
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","grids")
      CALL XML_INCAR('NGX','I',NGX,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NGY','I',NGY,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NGZ','I',NGZ,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NGXF','I',NGXF,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NGYF','I',NGYF,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NGZF','I',NGZF,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('ADDGRID','L',IDUM,RDUM,CDUM,LADDGRID,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","ionic")
      CALL XML_INCAR('NSW','I',NSW,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('IBRION','I',IBRION,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('ISIF','I',ISIF,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('PSTRESS','F',IDUM,PSTRESS,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('EDIFFG','F',IDUM,EDIFFG,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NFREE','I',NFREE,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('POTIM','F',IDUM,POTIM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('SMASS','F',IDUM,SMASS,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('SCALEE','F',IDUM,SCALEE,CDUM,LDUM,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","ionic md")
      CALL XML_INCAR('TEBEG','F',IDUM,TEBEG,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('TEEND','F',IDUM,TEEND,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NBLOCK','I',NBLOCK,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('KBLOCK','I',KBLOCK,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NPACO','I',NPACO,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('APACO','F',IDUM,APACO,CDUM,LDUM,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","symmetry")
      CALL XML_INCAR('ISYM','I',ISYM,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('SYMPREC','F',IDUM,TINY,CDUM,LDUM,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","dos")
      CALL XML_INCAR('LORBIT','L',IDUM,RDUM,CDUM,LORBIT_,CHARAC,1)
      CALL XML_INCAR_V('RWIGS','F',IDUM,RWIGS,CDUM,LDUM,CHARAC,NTYPP)
      CALL XML_INCAR('NEDOS','I',NEDOS,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('EMIN','F',IDUM,EMIN,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('EMAX','F',IDUM,EMAX,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('EFERMI','F',IDUM,EFERMI,CDUM,LDUM,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","writing")
      CALL XML_INCAR('NWRITE','I',NWRITE,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('LWAVE','L',IDUM,RDUM,CDUM,LWAVE,CHARAC,1)
      CALL XML_INCAR('LCHARG','L',IDUM,RDUM,CDUM,LCHARG,CHARAC,1)
      CALL XML_INCAR('LPARD','L',IDUM,RDUM,CDUM,LPARD,CHARAC,1)
      CALL XML_INCAR('LVTOT','L',IDUM,RDUM,CDUM,LVTOT,CHARAC,1)
      CALL XML_INCAR('LVHAR','L',IDUM,RDUM,CDUM,LVHAR,CHARAC,1)
      CALL XML_INCAR('LELF','L',IDUM,RDUM,CDUM,LELF,CHARAC,1)
      CALL XML_INCAR('LOPTICS','L',IDUM,RDUM,CDUM,LOPTICS,CHARAC,1)
      CALL XML_INCAR_V('STM','F',IDUM,STM,CDUM,LDUM,CHARAC,7)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","performance")
      CALL XML_INCAR('NPAR','I',NPAR,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NSIM','I',NSIM,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('NBLK','I',NBLK,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('LPLANE','L',IDUM,RDUM,CDUM,LPLANE_WISE,CHARAC,1)
      CALL XML_INCAR('LSCALAPACK','L',IDUM,RDUM,CDUM,LscaLAPACK,CHARAC,1)
      CALL XML_INCAR('LSCAAWARE','L',IDUM,RDUM,CDUM,LSCAAWARE,CHARAC,1)
      CALL XML_INCAR('LSCALU','L',IDUM,RDUM,CDUM,LSCALU,CHARAC,1)
      CALL XML_INCAR('LASYNC','L',IDUM,RDUM,CDUM,LASYNC,CHARAC,1)
      CALL XML_INCAR('LORBITALREAL','L',IDUM,RDUM,CDUM,LORBITALREAL,CHARAC,1)
      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","miscellaneous")
      CALL XML_INCAR('IDIOT','I',IDIOT,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('LMUSIC','L',IDUM,RDUM,CDUM,LMUSIC,CHARAC,1)
      CALL XML_INCAR_V('POMASS','F',IDUM,POMASS,CDUM,LDUM,CHARAC,NTYPIN)
      CALL XML_INCAR_V('DARWINR','F',IDUM,DARWIN_R,CDUM,LDUM,CHARAC,NTYPIN)
      CALL XML_INCAR_V('DARWINV','F',IDUM,DARWIN_V,CDUM,LDUM,CHARAC,NTYPIN)
      CALL XML_INCAR('LCORR','L',IDUM,RDUM,CDUM,LCORR,CHARAC,1) 
      CALL XML_CLOSE_TAG
      

    END SUBROUTINE XML_WRITER
