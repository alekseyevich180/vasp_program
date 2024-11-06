#include "symbol.inc"

module gpu
  implicit none
  logical, save :: use_gpu = .true.
end module gpu

!Hacene Mohamed
!
!The GPU data type

#ifndef GPU_DATA_INC
#define GPU_DATA_INC

module GPU_DATA

  USE prec
  
  implicit none
  type, public :: GPU_TYPE
      INTEGER(8) CPTWFP,CFW,CVR,CPTDUM,CH
      INTEGER(8) CW1,CW2,CR1,CR2,CR1_AL,CR,CWW
      INTEGER(8) CR1_ALL

      INTEGER(8)  CPROJ,CPROW,COVL,CF,CFA,CBLOCK,CTRANS,CPROJ1,CPROJ2,CHAM,CAP,CA
      !GDEF CPROJ,CPROW,COVL,CF,CFA,CBLOCK,CTRANS,CPROJ1,CPROJ2,CHAM,CAP,CA
      REAL(q)     DATAKE
      RGRID       SV,CDWORK
      INTEGER(8)  NINDPW
  end type GPU_TYPE

  type, public :: GPU_W
      COMPLEX(q)  CPTWFP,CFW,CR,CVR,CPTDUM
      GDEF        CPROJ,CPROW,COVL,CF,CFA,CBLOCK,CTRANS
      RGRID       SV,CDWORK
  end type GPU_W

  
end module GPU_DATA

#endif
