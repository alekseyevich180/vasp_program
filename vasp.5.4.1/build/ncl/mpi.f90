# 1 "mpi.F"
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

# 2 "mpi.F" 2 




!=======================================================================
!
! in most 1 implementations the collective communcation
! is slow, hence by default we avoid them
! if you want to use them define 1 (here or in the makefile)
!
! 1 use MPI_alltoall
! avoid_async    tries to post syncronised send and read operations
!                such that collisions are avoided, this is usually slower
! PROC_GROUP     does communication is block wise in a group of
!                roughly PROC_GROUP processors
!                this reduces collisions
!
! in addition our own implementation of the collective communication
! routines allow
!=======================================================================

!#define 1
!#define avoid_async





# 40

!
! OPENMPI seems to have a bcast bug
! whenever a node initiates a bcast the node master node seems
! to return immediately, whereas all the other nodes seem
! to wait until the master node initiates ANOTHER 1 call
! a simple work around is to initiate a barrier command after
! every single bcast
! to do this set MPI_bcast_with_barrier the makefile
! or undocument this line

!
! alternatively the main VASP code can initiate a barrier after
! a block of bcast calls
! this might be more efficient since less barrier are initiated
! this requires to define MPI_barrier_after_bcast
! in the makefile.
! The corresponding barriers have been inserted in
! wave_high.F and wave_mpi.F

!=======================================================================
!
! 1 communication routines for VASP
! all communication should be 1._q using this interface to allow
! adaption of other communication routines
! routines were entirely rewritten by Kresse Georg,
! but functionallity is similar to a module written by Peter
! Lockey at Daresbury
!
!======================================================================
      MODULE mpimy
      USE prec

!
! communication description include file
! this structure allows the use a several communicators
! it currently supports 1 and CRAY shmem calls
!
! Nbranch is the number of branches at each node in the gsum routines
! two should be always fine
!
      INTEGER Nbranch
      PARAMETER( Nbranch=2 )

      TYPE communic
!only COMM
        INTEGER MPI_COMM          ! MPI_Communicator
        INTEGER NODE_ME           ! node id starting from 1 ... NCPU
        INTEGER IONODE            ! node which has to do IO (set to 0 for no IO)
        INTEGER NCPU              ! total number of proc in this communicator
# 92

      END TYPE

! Standard 1 include file.
! I would like to have everything in the header but "freaking" SGI
! compiler can not handle this, thus I have to use an include file
      INCLUDE "mpif.h"

# 106

! There are no global local sum routines in 1, thus some workspace
! is required to store the results of the global sum
      INTEGER,PARAMETER ::  NZTMP=8000/2, NDTMP=8000, NITMP=8000
! workspace for integer, complex, and real
      COMPLEX(q),SAVE :: ZTMP_m(NZTMP)
      REAL(q),SAVE    :: DTMP_m(NDTMP)
      INTEGER,SAVE    :: ITMP_m(NITMP)

# 118



      CONTAINS

# 266

!----------------------------------------------------------------------
!
! M_init: initialise the basic communications
! (number of nodes, determine ionode)
!
!----------------------------------------------------------------------
!
# 276

      SUBROUTINE M_init(COMM)

      IMPLICIT NONE
      INCLUDE "pm.inc"

     TYPE(communic) COMM
      INTEGER i, ierror
# 286


      call MPI_init( ierror )
      IF ( ierror /= MPI_success ) THEN
         WRITE(*,*) 'Initpm: Error in MPI_init'
        STOP
      ENDIF
!
! initial communicator is world wide
! set only NCPU, NODE_ME and IONODE
! no internal setup 1._q at this point
!
# 306

      COMM%MPI_COMM= MPI_comm_world


      call MPI_comm_rank( COMM%MPI_COMM, COMM%NODE_ME, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Initpm: Error in MPI_comm_rank',ierror)
      COMM%NODE_ME= COMM%NODE_ME+1

      call MPI_comm_size( COMM%MPI_COMM, COMM%NCPU , ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Initpm: Error in MPI_comm_size',ierror)

      COMM%IONODE = 1

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_divide: creates a 2 dimensional cartesian topology
!  and a communicator along rows and columns of the process matrix
!
!----------------------------------------------------------------------

      SUBROUTINE M_divide( COMM, NPAR, COMM_INTER, COMM_INB, reorder)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM, COMM_INTER, COMM_INB, COMM_CART
      INTEGER NPAR,NPAR_2
      INTEGER, PARAMETER :: ndims=2
      INTEGER :: dims(ndims)
      LOGICAL :: periods(ndims), reorder, remain_dims(ndims)
      INTEGER :: ierror

      IF (NPAR >= COMM%NCPU) NPAR=COMM%NCPU
      dims(1)       = NPAR
      dims(2)       = COMM%NCPU/ NPAR
      IF (dims(1)*dims(2) /= COMM%NCPU ) THEN
         WRITE(0,*) 'M_divide: can not subdivide ',COMM%NCPU,'nodes by',NPAR
      ENDIF

      periods(ndims)=.FALSE.

      CALL MPI_Cart_create( COMM%MPI_COMM , ndims, dims, periods, reorder, &
                COMM_CART%MPI_COMM , ierror)
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Dividepm: Error in MPI_Cart_create', ierror)
! create the in-band communicator
      remain_dims(1)= .FALSE.
      remain_dims(2)= .TRUE.

      CALL MPI_Cart_sub( COMM_CART%MPI_COMM, remain_dims, COMM_INB%MPI_COMM, ierror )
      IF ( ierror /= MPI_success )&
         CALL M_stop_ierr('Dividepm: Error in MPI_Cart_sub (1) ', ierror)

! create the inter-band communicator
      remain_dims(1)= .TRUE.
      remain_dims(2)= .FALSE.

      CALL MPI_Cart_sub( COMM_CART%MPI_COMM, remain_dims, COMM_INTER%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Dividepm: Error in MPI_Cart_sub (2) ', ierror)
! overwrite initial communicator by new (1._q,0._q)
      COMM=COMM_CART


      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_divide_shmem:
!
!----------------------------------------------------------------------

      SUBROUTINE M_divide_shmem(COMM_WORLD,COMM,N,COMM_SHMEM)
      USE iso_c_binding
      IMPLICIT NONE
      INCLUDE "pm.inc"
      TYPE (communic) COMM_WORLD,COMM,COMM_SHMEM
      INTEGER N
! local variables
      CHARACTER*(MPI_MAX_PROCESSOR_NAME) myname,pname
      INTEGER :: COMM_group,COMM_SHMEM_group
      INTEGER :: group(0:COMM%NCPU-1)
      INTEGER :: resultlen,id_in_group
      INTEGER :: ierror
      INTEGER :: I,IGRP
! shmem variables for sanity check
      INTEGER :: shmid
      TYPE(c_ptr) :: address
      INTEGER*8 k

      INTEGER*4, POINTER :: IDSHMEM(:)
      INTEGER :: ID

      IF (N>0) THEN
! explicitly specified division of COMM
         I=COMM%NCPU/N
         IF (I*N/=COMM%NCPU) THEN
            WRITE(*,*) 'M_divide_shmem: ERROR: can not subdivide ',COMM%NCPU,' nodes by',N
            STOP
         ENDIF
         IGRP=0
         DO I=1,COMM%NCPU
            IF ((I-1)/N==(COMM%NODE_ME-1)/N) THEN
               group(IGRP)=I-1
               IGRP=IGRP+1
            ENDIF
         ENDDO
      ENDIF
      IF (N<0) THEN
! attempt automatic division of COMM
         IGRP=0
         CALL MPI_Get_processor_name(myname,resultlen,ierror)
         DO I=1,COMM%NCPU
            IF (I==COMM%NODE_ME) pname=myname
            CALL MPI_bcast(pname,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,I-1,COMM%MPI_COMM,ierror)
            IF (ierror/=MPI_success) &
               CALL M_stop_ierr('ERROR: MPI_bcast (M_divide_shmem:1) returns',ierror)

            CALL MPI_barrier(COMM%MPI_COMM,ierror)

            IF (pname(1:resultlen)==myname(1:resultlen)) THEN
               group(IGRP)=I-1
               IGRP=IGRP+1
            ENDIF
         ENDDO
      ENDIF

      CALL MPI_Comm_group(COMM%MPI_COMM,COMM_group,ierror)
      CALL MPI_Group_incl(COMM_group,IGRP,group,COMM_SHMEM_group,ierror)
      CALL MPI_Comm_create(COMM%MPI_COMM,COMM_SHMEM_group,COMM_SHMEM%MPI_COMM,ierror)

      CALL M_initc(COMM_SHMEM)

! sanity check
      k=1
      IF (COMM_SHMEM%NODE_ME==1) CALL getshmem(8*k,shmid)
      CALL M_bcast_i(COMM_SHMEM,shmid,1)
      CALL attachshmem(shmid,address)
      CALL c_f_pointer(address,IDSHMEM,[k])

      IF (COMM_SHMEM%NODE_ME==1) THEN
         ID=COMM_WORLD%NODE_ME; IDSHMEM=ID
      ENDIF
      CALL M_bcast_i(COMM_SHMEM,ID,1)

      ierror=0
      IF (ID/=IDSHMEM(1)) ierror=1
      CALL M_sum_i(COMM_WORLD,ierror,1)

      CALL detachshmem(address)
      IF (COMM_SHMEM%NODE_ME==1) CALL destroyshmem(shmid)

      IF (ierror>0) THEN
         WRITE(*,*) 'M_divide_shmem: ERROR: not all procs in COMM_SHMEM seem to be on the same physical node.'
         STOP
      ENDIF

      N=COMM_SHMEM%NCPU

!     DO I=1,COMM_WORLD%NCPU
!        IF (COMM_WORLD%NODE_ME==I) THEN
!           WRITE(*,'(A,I4,X,A,I4)') 'global id:',comm_world%node_me,'shmem id:',comm_shmem%node_me
!        ENDIF
!        CALL MPI_barrier(COMM_WORLD%MPI_COMM,ierror)
!     ENDDO
      RETURN
      END SUBROUTINE


!----------------------------------------------------------------------
!
! M_initc: initialise a communicator
!
!----------------------------------------------------------------------

      SUBROUTINE M_initc( COMM)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, ierror
      INTEGER id_in_group

      call MPI_comm_rank( COMM%MPI_COMM, id_in_group, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Initpm: Error in MPI_comm_rank',ierror)

      call MPI_comm_size( COMM%MPI_COMM, COMM%NCPU, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Initpm: Error in MPI_comm_size',ierror)

      COMM%NODE_ME = id_in_group + 1
      CALL init_hard_ids(COMM)

      CALL MPI_barrier( COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Initpm: Error in MPI_barrier',ierror)
      COMM%IONODE =1
      END SUBROUTINE

!----------------------------------------------------------------------
!
! init_hard_ids: map the virtual (1 node_id) to the real node_id
!  (T3E, T3D specific)
!
!----------------------------------------------------------------------

      SUBROUTINE init_hard_ids(COMM)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, tmp, ierror, id_in_group
      INTEGER, ALLOCATABLE :: hid_tmp(:)
# 541

      END SUBROUTINE


      END MODULE
!======================================================================
!
! all other routines are often called with either
! real or complex arrays, vectors or scalars, so I can not put
! them into the F90 module
!
!======================================================================
!----------------------------------------------------------------------
!
! M_exit: exit 1, very simple just exit 1
!
!----------------------------------------------------------------------

      SUBROUTINE M_exit()
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, ierror, id_in_group

!      call MPI_barrier(MPI_comm_world, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('Exitpm: Error in MPI_barrier',ierror)

      call MPI_finalize( ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Exitpm: Error in MPI_finalize',ierror)
      STOP

      RETURN
      END SUBROUTINE


!----------------------------------------------------------------------
!
! M_stop: exits 1 and program because of error, a message is
! printed
!
!----------------------------------------------------------------------

      SUBROUTINE M_stop(message)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      CHARACTER (LEN=*) message
      INTEGER ierror

      WRITE (*,*) message

      call MPI_abort(MPI_comm_world , 1, ierror )
      STOP

      RETURN
      END SUBROUTINE


      SUBROUTINE M_stop_ierr(message, ierror)
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      CHARACTER (LEN=*) message
      INTEGER ierror

      WRITE (*,*) message, ierror

      call MPI_abort(MPI_comm_world , 1, ierror )
      STOP

      RETURN
      END SUBROUTINE

!======================================================================
!
! Send and Receive routines, map directly onto 1
!
!======================================================================

!----------------------------------------------------------------------
!
! M_send_i: send n integers stored in ivec to node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_i (COMM, node, ivec, n)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      INTEGER ivec(n)

      INTEGER status(MPI_status_size), ierror

      call MPI_send( ivec(1), n, MPI_integer, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_send returns',ierror)

      RETURN
      END

!----------------------------------------------------------------------
!
! M_recv_i: receive n integers into array ivec from node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_i(COMM, node, ivec, n )
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      INTEGER ivec(n)
      INTEGER status(MPI_status_size), ierror

      call MPI_recv( ivec(1), n, MPI_integer , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_recv returns',ierror)

      RETURN
      END


!----------------------------------------------------------------------
!
! M_send_z: send n double complex stored in zvec to node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_z (COMM, node, zvec, n)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      COMPLEX(q) :: zvec(n)

      INTEGER status(MPI_status_size), ierror

      call MPI_send( zvec(1), n, MPI_double_complex, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_send returns',ierror)

      RETURN
      END

!----------------------------------------------------------------------
!
! M_recv_z: receive n double complex into array ivec from node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_z(COMM, node, zvec, n )
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      COMPLEX(q) :: zvec(n)
      INTEGER status(MPI_status_size), ierror

      call MPI_recv( zvec(1), n, MPI_double_complex , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_recv returns',ierror)

      RETURN
      END

!----------------------------------------------------------------------
!
! M_send_d: send n double stored in ivec to node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_d (COMM, node, dvec, n)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      REAL(q) :: dvec(n)

      INTEGER status(MPI_status_size), ierror

      call MPI_send( dvec(1), n, MPI_double_precision, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_send returns',ierror)

      RETURN
      END

!----------------------------------------------------------------------
!
! M_recv_d: receive n double  into array ivec from node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_d(COMM, node, dvec, n )
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      REAL(q) :: dvec(n)
      INTEGER status(MPI_status_size), ierror

      call MPI_recv( dvec(1), n, MPI_double_precision , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_recv returns',ierror)

      RETURN
      END


!======================================================================
!
! global sum and maximum routines
!
!======================================================================

# 1119

! split the arrays and copy for ALLREDUCE

!----------------------------------------------------------------------
!
! M_sum_i: performs a global sum on n integers in vector ivec
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_i(COMM, ivec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! check whether n is sensible
      IF (n==0 .OR. COMM%NCPU==1 ) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_sum_i ', n
         STOP
      END IF

!  there is no inplace global sum in 1, thus we have to use
!  a work array

      DO j = 1, n, NITMP
         ichunk = MIN( n-j+1 , NITMP)

         call MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
     &                       MPI_sum, COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_sum, MPI_integer) returns',ierror)

         DO k = 0, ichunk-1
            ivec(j+k) = ITMP_m(k+1)
         ENDDO

      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
!
! M_max_i: performs a global max on n integers in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_i(COMM, ivec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! check whether n is sensible
      IF (n==0 .OR. COMM%NCPU==1 ) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_max_i ', n
         STOP
      END IF

!  there is no inplace global max in 1, thus we have to use
!  a work array

      DO j = 1, n, NITMP
         ichunk = MIN( n-j+1 , NITMP)

         call MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
                             MPI_max, COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_max, MPI_integer) returns',ierror)

         DO k = 0, ichunk-1
            ivec(j+k) = ITMP_m(k+1)
         ENDDO

      ENDDO

      RETURN
      END



!----------------------------------------------------------------------
!
! M_max_d: performs a global max search on n doubles in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_d(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_max_d ', n
         STOP
      END IF

!  there is no inplace global max in 1, thus we have to use
!  a work array
      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         call MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_max, &
                             COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_max, MPI_double_prec) returns',ierror)

         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)

      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
!
! M_sumb_d: performs a global sum on n doubles in vector vec
!  uses MPI_allreduce which is usually very inefficient
!  faster alternative routines can be found below
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_d(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_sumb_d ', n
         STOP
      END IF

!  there is no inplace global sum in 1, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         call MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_sum, &
                             COMM%MPI_COMM, ierror )

         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_sum, MPI_double_prec) returns',ierror)

         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)

      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
!
! M_sumb_d: performs a global sum on n singles in vector vec
!  uses MPI_allreduce which is usually very inefficient
!  faster alternative routines can be found below
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_s(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_sumb_s ', n
         STOP
      END IF

!  there is no inplace global sum in 1, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)
         call MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_real, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_sum, MPI_real) returns',ierror)

         CALL SCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
!
! to make live easier, a global sum for scalars
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_s(COMM, n, v1, v2, v3, v4)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n),v1,v2,v3,v4

      vec=0

      IF (n>0) vec(1)=v1
      IF (n>1) vec(2)=v2
      IF (n>2) vec(3)=v3
      IF (n>3) vec(4)=v4
      IF (n>4) THEN
          WRITE(*,*) 'internal ERROR: invalid n in M_sum_s ', n
          STOP
      END IF

      CALL M_sumb_d(COMM, vec, n)

      IF (n>0) v1=vec(1)
      IF (n>1) v2=vec(2)
      IF (n>2) v3=vec(3)
      IF (n>3) v4=vec(4)
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sum_z: performs a global sum on n complex items in vector vec
!  uses MPI_allreduce which is usually very inefficient
!  faster alternative routines can be found below
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_z(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER j

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_sumb_z ', n
         STOP
      END IF


!  there is no inplace global sum in 1, thus we have to use
!  a work array
      DO j = 1, n, NZTMP
         ichunk = MIN( n-j+1 , NZTMP)

         call MPI_allreduce( vec(j), ZTMP_m(1), ichunk, &
                             MPI_double_complex, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_sum, MPI_double_complex) returns',ierror)

         CALL ZCOPY(ichunk , ZTMP_m(1), 1 ,  vec(j) , 1)

      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
!
! M_prodb_z: performs a global product on n complex numbers in vec
!  uses MPI_allreduce which is usually very inefficient
!
!----------------------------------------------------------------------

      SUBROUTINE M_prodb_z(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_prodb_z ', n
         STOP
      END IF

!  there is no inplace global sum in 1, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         call MPI_allreduce( vec(j), ZTMP_m(1), ichunk, &
                             MPI_double_complex, MPI_prod, &
                             COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_prod, MPI_double_complex) returns',ierror)

         CALL ZCOPY(ichunk , ZTMP_m(1), 1 ,  vec(j) , 1)

      ENDDO

      RETURN
      END



!======================================================================
!
! Global barrier routine
!
!======================================================================

      SUBROUTINE M_barrier(COMM )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      call MPI_barrier( COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_barrier (M_barrier) returns',ierror)

      RETURN
      END


!======================================================================
!
! Global Copy Routines
!
!======================================================================

!----------------------------------------------------------------------
!
! M_bcast_i: copy n integers from root to all nodes
!
!----------------------------------------------------------------------


      SUBROUTINE M_bcast_i(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER vec(n)

      INTEGER ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_i'
         STOP
      END IF

      call MPI_bcast( vec(1), n, MPI_integer, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_i) returns',ierror)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )

 

      RETURN
      END

!----------------------------------------------------------------------
!
! M_bcast_l: copy n logical from root to all nodes
!
!----------------------------------------------------------------------


      SUBROUTINE M_bcast_l(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      LOGICAL vec(n)

      INTEGER ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_l'
         STOP
      END IF

      call MPI_bcast( vec(1), n, MPI_logical, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_l) returns',ierror)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      RETURN
      END

!----------------------------------------------------------------------
!
! M_bcast_i_from: copy n integers from node inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_i_from(COMM, vec, n , inode)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER inode
      INTEGER vec(n)

      INTEGER ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_i'
         STOP
      END IF

      call MPI_bcast( vec(1), n, MPI_integer, inode-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_i) returns',ierror)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      RETURN
      END


!----------------------------------------------------------------------
!
! M_bcast_d: copy n double precision from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_d(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      INTEGER ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_i'
         STOP
      END IF

      call MPI_bcast( vec(1), n,  MPI_double_precision, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_d) returns',ierror)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      RETURN
      END


!----------------------------------------------------------------------
!
! M_bcast_s: copy n single precision from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_s(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)

      INTEGER ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_i'
         STOP
      END IF

      call MPI_bcast( vec(1), n,  MPI_real, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_s) returns',ierror)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )

      RETURN
      END


!----------------------------------------------------------------------
!
! M_bcast_z: copy n double precision complex from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_z(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)

      INTEGER ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_i'
         STOP
      END IF

      call MPI_bcast( vec(1), n,  MPI_double_complex, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_z) returns',ierror)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )

      RETURN
      END

!----------------------------------------------------------------------
!
! M_bcast_z_from: copy n double precision complex from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_z_from(COMM, vec, n, inode )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)

      INTEGER inode,ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_z_from'
         STOP
      END IF

      call MPI_bcast( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_z_from) returns',ierror)

      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      RETURN
      END

!----------------------------------------------------------------------
!
! M_bcast_d_from: copy n double precision from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_d_from(COMM, vec, n, inode )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      INTEGER inode,ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_d_from'
         STOP
      END IF

      call MPI_bcast( vec(1), n,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_d_from) returns',ierror)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      RETURN
      END


!======================================================================
!
! Global Exchange Routine
!
!======================================================================

!----------------------------------------------------------------------
!
! M_alltoallv_z: complex global exchange routine which maps directly onto
!  MPI_alltoallv
!  since MPI_alltoallv is usually very slow, an alternative implementation
!  optimised for clusters exists
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoallv_z(COMM, xsnd, psnd, nsnd, xrcv, prcv, nrcv, rprcv)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      COMPLEX(q) xsnd(*)         ! send buffer
      COMPLEX(q) xrcv(*)         ! receive buffer

      INTEGER psnd(COMM%NCPU+1) ! location of data in send buffer (0 based)
      INTEGER prcv(COMM%NCPU+1) ! location of data in recv buffer (0 based)
      INTEGER rprcv(COMM%NCPU)  ! remote location of data
      INTEGER nsnd(COMM%NCPU+1) ! number of data send to each node
      INTEGER nrcv(COMM%NCPU+1) ! number of data recv from each node

!----------------------------------------------------------------------
# 1902


!----------------------------------------------------------------------
      INTEGER ierror

      call MPI_alltoallv( xsnd(1), nsnd(1), psnd(1), MPI_double_complex, &
                          xrcv(1), nrcv(1), prcv(1), MPI_double_complex, &
                          COMM%MPI_COMM, ierror )

      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: M_alltoallv_z MPI_alltoallv returns',ierror)

!----------------------------------------------------------------------
# 1987


!----------------------------------------------------------------------

      RETURN
      END


!----------------------------------------------------------------------
!
! M_cycle_d: real cyclic exchange routine which maps directly onto
!
!----------------------------------------------------------------------

      SUBROUTINE M_cycle_d(COMM, xsnd, nsnd)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"


      TYPE(communic) COMM
      REAL(q) xsnd(nsnd)         ! send/ receive buffer
      INTEGER nsnd

      INTEGER ierror, i, j, ichunk
      INTEGER :: tag=202
      INTEGER :: request(2)
      INTEGER :: A_OF_STATUSES(2)
      INTEGER :: NDTMP_=3

      DO j = 1, nsnd, NDTMP_
         ichunk = MIN( nsnd-j+1 , NDTMP_)

! initiate the receive from node-1 (note (0._q,0._q) based)
         i = MOD(COMM%NODE_ME-1-1+COMM%NCPU , COMM%NCPU) 
         call MPI_irecv( DTMP_m(1), ichunk, MPI_double_precision, &
              i, tag,  COMM%MPI_COMM, request(1), ierror )
         IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_alltoallv_d irecv returns',ierror)

! initiate the send
         i = MOD(COMM%NODE_ME-1+1 , COMM%NCPU)  ! i (0._q,0._q) based
         call MPI_isend( xsnd(j), ichunk, MPI_double_precision, &
              i, tag,  COMM%MPI_COMM, request(2), ierror )
         IF ( ierror /= MPI_success ) &
               CALL M_stop_ierr('ERROR: M_alltoallv_d irecv returns',ierror)

! wait for send and receive to finish
         call MPI_waitall(2 , request, A_OF_STATUSES, ierror)
         IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_waitall in M_alltoall_d returns',ierror)

! copy result back
         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  xsnd(j) , 1)

      END DO
      END SUBROUTINE M_cycle_d

!----------------------------------------------------------------------
!
! M_alltoall_i: integer routine which maps directly onto
!  MPI_alltoallv
!  this is used only once by VASP
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_i(COMM, xsnd, psnd, nsnd, xrcv, prcv, nrcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER xsnd(*)           ! send buffer
      INTEGER xrcv(*)           ! receive buffer

      INTEGER psnd(COMM%NCPU+1) ! location of data in send buffer (0 based)
      INTEGER prcv(COMM%NCPU+1) ! location of data in recv buffer (0 based)
      INTEGER rprcv(COMM%NCPU)  ! remote location of data
      INTEGER nsnd(COMM%NCPU+1) ! number of data send to each node
      INTEGER nrcv(COMM%NCPU+1) ! number of data recv from each node

! local data
      INTEGER ierror

      call MPI_alltoallv( xsnd(1), nsnd(1), psnd(1), MPI_integer, &
                          xrcv(1), nrcv(1), prcv(1), MPI_integer, &
                          COMM%MPI_COMM, ierror )

      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: imexch MPI_alltoallv returns',ierror)

      RETURN
      END

!----------------------------------------------------------------------
!
! on the T3E M_alltoallv_z can use shmemput instead of 1
! M_alltoallv_z_prepare prepares the additional array rprcv
! which contains the remote locations
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoallv_raddr(COMM, prcv, rprcv)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM

      INTEGER prcv(COMM%NCPU+1) ! location of data in recv buffer (0 based)
      INTEGER rprcv(COMM%NCPU)  ! remote location of data
! local variable
      INTEGER ierror

! only (1._q,0._q) simple MPI_alltoall is required
      CALL MPI_alltoall( prcv(1),  1, MPI_integer, &
                         rprcv(1), 1, MPI_integer, &
                         COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: M_alltoallv_z_prepare',ierror)


      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_alltoallv_simple requires as only input the number of data nsnd send
! from each node to each other node
! it assumes a continous data arrangement on sender and receiver
! and sets up the arrays which are required for MPI_alltoall
! nrcv, psnd, prcv
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoallv_simple(COMM, nsnd, nrcv, psnd, prcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
! input
      INTEGER nsnd(COMM%NCPU)   ! number of data send to each node
! output
      INTEGER psnd(COMM%NCPU+1) ! location of data in send buffer (0 based)
      INTEGER nrcv(COMM%NCPU)   ! number of data recv from each node
      INTEGER prcv(COMM%NCPU+1) ! location of data in recv buffer (0 based)
! local variable
      INTEGER ierror,i

! only (1._q,0._q) simple MPI_alltoall is required
! to find number of received data on each node
      CALL MPI_alltoall( nsnd(1),  1, MPI_integer, &
                         nrcv(1),  1, MPI_integer, &
                         COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: M_alltoallv_z_prepare',ierror)
! now set the locations assuming linear arrangement of data

      psnd(1)=0
      prcv(1)=0

      DO i=1,COMM%NCPU
        psnd(i+1)=psnd(i)+nsnd(i)
        prcv(i+1)=prcv(i)+nrcv(i)
      ENDDO

      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_alltoall_d: complex and real global exchange routine
!     redistributes an array from distribution over bands to
!     distribution over coefficient (or vice versa)
!     original distribution            final distribution
!     |  1  |  2  |  3  |  4  |      |  1  |  1  |  1  |  1  |
!     |  1  |  2  |  3  |  4  |      |  2  |  2  |  2  |  2  |
!     |  1  |  2  |  3  |  4  | <->  |  3  |  3  |  3  |  3  |
!     |  1  |  2  |  3  |  4  |      |  4  |  4  |  4  |  4  |
!
!     xsnd is the array to be redistributed (having n elements)
!     xrcv is the result array with n/NCPU elements received
!     from each processore
!
!     mind that only (n/NCPU) *NCPU data are exchanged
!     it is the responsability of the user to guarantee that n is
!     correct
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_d(COMM, n, xsnd, xrcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st,i 
!----------------------------------------------------------------------
# 2227

!----------------------------------------------------------------------

      INTEGER, SAVE :: tag=201
      INTEGER       :: in
      INTEGER       :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
      INTEGER, PARAMETER :: max_=8000
      INTEGER       :: block, p, sndcount_
      INTEGER       :: actual_proc_group, com_proc_group, &
           proc_group, group_base, i_in_group, irequests

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

!----------------------------------------------------------------------

!----------------------------------------------------------------------

      call MPI_alltoall( xsnd(1), sndcount, MPI_double_precision, &
     &                   xrcv(1), rcvcount, MPI_double_precision, &
     &                   COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: M_alltoall_d returns',ierror)

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('Initpm: Error in MPI_barrier',ierror)

!----------------------------------------------------------------------
# 2440

!----------------------------------------------------------------------


      END SUBROUTINE


!----------------------------------------------------------------------
!
! M_alltoall_s: complex and real global exchange routine
!     redistributes an array from distribution over bands to
!     distribution over coefficient (or vice versa)
!     original distribution            final distribution
!     |  1  |  2  |  3  |  4  |      |  1  |  1  |  1  |  1  |
!     |  1  |  2  |  3  |  4  |      |  2  |  2  |  2  |  2  |
!     |  1  |  2  |  3  |  4  | <->  |  3  |  3  |  3  |  3  |
!     |  1  |  2  |  3  |  4  |      |  4  |  4  |  4  |  4  |
!
!     xsnd is the array to be redistributed (having n elements)
!     xrcv is the result array with n/NCPU elements received
!     from each processore
!
!     mind that only (n/NCPU) *NCPU data are exchanged
!     it is the responsability of the user to guarantee that n is
!     correct
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_s(COMM, n, xsnd, xrcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st,i 

      INTEGER, SAVE :: tag=201
      INTEGER       :: in
      INTEGER       :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
      INTEGER, PARAMETER :: max_=8000
      INTEGER       :: block, p, sndcount_
      INTEGER       :: actual_proc_group, com_proc_group, &
           proc_group, group_base, i_in_group, irequests

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

!----------------------------------------------------------------------

!----------------------------------------------------------------------
      call MPI_alltoall( xsnd(1), sndcount, MPI_real, &
     &                   xrcv(1), rcvcount, MPI_real, &
     &                   COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: M_alltoall_s returns',ierror)

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('Initpm: Error in MPI_barrier',ierror)

!----------------------------------------------------------------------
# 2686

!----------------------------------------------------------------------

      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_alltoall_z: uses M_alltoall_d with twice as many elements
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_z(COMM, n, xsnd, xrcv )
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) xsnd(n), xrcv(n)

      CALL M_alltoall_d(COMM, n*2, xsnd, xrcv)
      END SUBROUTINE

!----------------------------------------------------------------------
!
! z/M_alltoall_d: complex and real global exchange routine
!     redistributes an array from distribution over bands to
!     distribution over coefficient (or vice versa)
!     original distribution            final distribution
!     |  1  |  2  |  3  |  4  |      |  1  |  1  |  1  |  1  |
!     |  1  |  2  |  3  |  4  |      |  2  |  2  |  2  |  2  |
!     |  1  |  2  |  3  |  4  | <->  |  3  |  3  |  3  |  3  |
!     |  1  |  2  |  3  |  4  |      |  4  |  4  |  4  |  4  |
!
!     xsnd is the array to be redistributed (having n elements)
!     xrcv is the result array with n/NCPU elements received
!     from each processore
!
!     mind that only (n/NCPU) *NCPU data are exchanged
!     it is the responsability of the user to guarantee that n is
!     correct
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_d_async(COMM, n, xsnd, xrcv, tag, srequest, rrequest )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER tag
      INTEGER srequest(COMM%NCPU), rrequest(COMM%NCPU)
      REAL(q) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st
      INTEGER i,j,in

      IF (COMM%NCPU == 1) RETURN

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

! local memory copy for data kept on the local node
      DO i = 1,sndcount
         xrcv((COMM%NODE_ME-1)*sndcount + i) = xsnd((COMM%NODE_ME-1)*sndcount + i)
      ENDDO

! initiate send and receive on all nodes
! local copy has already been 1._q each node send NCPU-1 packages
      j=1
      DO in = 0, COMM%NCPU-1
! send to node in + own node id
! such a construct should allow for an efficient use of the network
         
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)
         IF ( COMM%NODE_ME-1 /= i) THEN
            call MPI_isend( xsnd(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, srequest(j), ierror )
            IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_alltoall_d_async MPI_alltoall returns',ierror)
            call MPI_irecv( xrcv(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, rrequest(j), ierror )
            IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_alltoall_d_async MPI_alltoall returns',ierror)
            j=j+1
         ENDIF
      ENDDO
        
      END SUBROUTINE

      SUBROUTINE M_alltoall_wait(COMM, srequest, rrequest )

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER srequest(COMM%NCPU), rrequest(COMM%NCPU)
      INTEGER ierror
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,COMM%NCPU)

! wait for the NCPU-1 outstanding packages

      call MPI_waitall(COMM%NCPU-1, srequest, A_OF_STATUSES, ierror)
            IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_waitall returns',ierror)
      call MPI_waitall(COMM%NCPU-1, rrequest, A_OF_STATUSES, ierror)
            IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_waitall returns',ierror)

      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sumf_d: performs a fast global sum on n doubles in
! vector vec (algorithm by Kresse Georg)
!
! uses complete interchange algorithm (my own invention, but I guess
!  some people must know it)
! exchange data between nodes, sum locally and
! interchange back, this algorithm is faster than typical 1 based
! algorithms (on 8 nodes under MPICH a factor 4)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_d(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n,ncount,nsummed,ndo,i,j, info, n_,mmax
      REAL(q) vec(n)
!----------------------------------------------------------------------
# 2836

!----------------------------------------------------------------------
      REAL(q), ALLOCATABLE :: vec_inter(:)
! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=8000

! quick return if possible
      IF (COMM%NCPU == 1) RETURN
      
      mmax=MIN(n/COMM%NCPU,max_)
      ALLOCATE(vec_inter(mmax*COMM%NCPU))
!----------------------------------------------------------------------

!----------------------------------------------------------------------

      nsummed=0
      n_=n/COMM%NCPU

      DO ndo=0,n_-1,mmax
! forward exchange
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
! sum localy
         DO i=2, COMM%NCPU
           CALL DAXPY(ncount, 1.0_q, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
! replicate data (will be send to each proc)
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
! backward exchange
         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
         CALL DCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
      ENDDO
      
! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         WRITE(0,*) 'internal error in M_sumf_d',n_,nsummed
         STOP
      ENDIF

      IF (n-nsummed /= 0 ) &
        CALL M_sumb_d(COMM, vec(nsummed+1), n-nsummed)

# 2891

      DEALLOCATE(vec_inter)

      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sumf_s: performs a fast global sum on n singles in
! vector vec (algorithm by Kresse Georg)
!
! uses complete interchange algorithm (my own invention, but I guess
!  some people must know it)
! exchange data between nodes, sum locally and
! interchange back, this algorithm is faster than typical 1 based
! algorithms (on 8 nodes under MPICH a factor 4)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_s(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n,ncount,nsummed,ndo,i,j, info, n_,mmax
      REAL(qs) vec(n)
      REAL(qs), ALLOCATABLE :: vec_inter(:)
! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=8000

! quick return if possible
      IF (COMM%NCPU == 1) RETURN
      
      mmax=MIN(n/COMM%NCPU,max_)
      ALLOCATE(vec_inter(mmax*COMM%NCPU))

      nsummed=0
      n_=n/COMM%NCPU

      DO ndo=0,n_-1,mmax
! forward exchange
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_s(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
! sum localy
         DO i=2, COMM%NCPU
           CALL SAXPY(ncount, 1.0, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
! replicate data (will be send to each proc)
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
! backward exchange
         CALL M_alltoall_s(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
         CALL SCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
      ENDDO
      
! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         WRITE(0,*) 'internal error in M_sumf_s',n_,nsummed
         STOP
      ENDIF

      IF (n-nsummed /= 0 ) &
        CALL M_sumb_s(COMM, vec(nsummed+1), n-nsummed)

# 2966

      DEALLOCATE(vec_inter)

      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sumf_z: performs a fast global sum on n complex in
! vector 'vec' (see above)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      CALL M_sumf_d(COMM, vec, 2*n)
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sum_z: performs a sum on n double complex numbers
!  it uses either sumb_d or sumf_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)


      CALL M_sumb_d(COMM, vec, 2*n)
# 3011

            
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sum_d:  performs a sum on n double prec numbers
!  it uses either sumb_d or sumf_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_d(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)


      CALL M_sumb_d(COMM, vec, n)
# 3038

            
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sum_s:  performs a sum on n single prec numbers
!  it uses either sumb_s or sumf_s
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_single(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)


      CALL M_sumb_s(COMM, vec, n)
# 3065

            
      END SUBROUTINE


!----------------------------------------------------------------------
!
! M_sum_master_d:  performs a sum on n double prec numbers to master
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_master_d(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER :: I,thisdata,ierror

      IF (COMM%NCPU == 1 ) RETURN

      DO i=1,n,8000
         thisdata=min(n-i+1,8000)

         call MPI_reduce( vec(i), DTMP_m(1), thisdata, &
                             MPI_double_precision, MPI_sum, &
                             0, COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: M_sum_master_d reduce returns',ierror)

         if (COMM%NODE_ME==1) THEN
            vec(i:i+thisdata-1)=DTMP_m(1:thisdata-1)
                             
         endif
      END DO

      END SUBROUTINE M_sum_master_d

!----------------------------------------------------------------------
!
! M_sum_master_d:  performs a sum on n double complex numbers to master
!
!----------------------------------------------------------------------
      SUBROUTINE M_sum_master_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      CALL M_sum_master_d(COMM,vec,2*n)

      END SUBROUTINE M_sum_master_z

# 3141

