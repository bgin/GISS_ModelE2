#include "rundeck_opts.h"
#ifdef MPI_DEFS_HACK
#include "mpi_defs.h"
#endif

#ifndef CUBED_SPHERE
#define DOMAIN_DECOMP_ATM_IS_1D
#endif

#ifdef NEW_IO
#define USE_DD2D_UTILS
#endif

      MODULE dist_grid_mod
!@sum  This module encapsulates lat-lon decomposition information for the
!@+    message passing implementation.
!@auth Tom Clune GSFC/SSSO/610.3

      use MpiSupport_mod, only: am_i_root
      use Domain_mod
      use Hidden_mod

#ifdef CUBED_SPHERE
#define USE_DD2D_UTILS
#undef CUBED_SPHERE
#endif

#ifdef USE_DD2D_UTILS
      use dd2d_utils, only : dist_grid,init_dist_grid
#endif

      IMPLICIT NONE

#ifdef USE_MPI
! NOTE: USE_MPI is #defined from the "make" command line using option MPI=YES.
#include "mpif.h"
#endif
      SAVE
      PRIVATE ! Except for

!@var DIST_GRID derived type to provide decomposition info
!@+   public components are used to minimize overhead for accessing
!@+   routine components
      PUBLIC :: DIST_GRID
      public :: setMpiCommunicator
!@var  grid Default decomposition; globally accessible for convenience.
      PUBLIC :: grid
!@var INIT_APP set some parameters
      PUBLIC :: INIT_APP
      PUBLIC :: INIT_GRID
      PUBLIC :: DESTROY_GRID
!@var FINISH_APP Cleans up at the end of the run (closes debugging file)
      PUBLIC :: FINISH_APP
!@var GLOBALMIN determine max value across pes
      PUBLIC :: GLOBALMIN
!@var GLOBALMAX determine max value across pes
      PUBLIC :: GLOBALMAX
!@var SUMXPE sum an array over processors without reducing its rank
      PUBLIC :: SUMXPE
!@var GET - extracts bounds information from DIST_GRID object
      PUBLIC :: GET
      PUBLIC :: HERE
      PUBLIC :: LOG_PARALLEL

      PUBLIC :: TRANSP
      PUBLIC :: TRANSPOSE_COLUMN
!@var GLOBALMIN Generic wrapper for Real
      INTERFACE GLOBALMIN
        MODULE PROCEDURE GLOBALMIN_R
      END INTERFACE

!@var GLOBALMAX Generic wrapper for Real/integer
      INTERFACE GLOBALMAX
        MODULE PROCEDURE GLOBALMAX_R
        MODULE PROCEDURE GLOBALMAX_I
        MODULE PROCEDURE GLOBALMAX_I_1D
      END INTERFACE

      INTERFACE TRANSP
        MODULE PROCEDURE TRANSPOSE_ij
        MODULE PROCEDURE TRANSPOSE_ijk
      END INTERFACE

      INTERFACE SUMXPE
        MODULE PROCEDURE SUMXPE_1D
        MODULE PROCEDURE SUMXPE_1D_I
        MODULE PROCEDURE SUMXPE_2D
        MODULE PROCEDURE SUMXPE_3D
        MODULE PROCEDURE SUMXPE_4D
!        MODULE PROCEDURE SUMXPE_5D
      END INTERFACE

!@var ESMF_BCAST Generic routine to broadcast data to all PEs.
      PUBLIC :: ESMF_BCAST
      INTERFACE ESMF_BCAST
        MODULE PROCEDURE ESMF_BCAST_0D
        MODULE PROCEDURE ESMF_BCAST_1D
        MODULE PROCEDURE ESMF_BCAST_2D
        MODULE PROCEDURE ESMF_BCAST_3D
        MODULE PROCEDURE ESMF_BCAST_4D
        MODULE PROCEDURE ESMF_IBCAST_0D
        MODULE PROCEDURE ESMF_IBCAST_1D
        MODULE PROCEDURE ESMF_IBCAST_2D
        MODULE PROCEDURE ESMF_IBCAST_3D
        MODULE PROCEDURE ESMF_IBCAST_4D
      END INTERFACE

!@var BAND_PACK Procedure in which each PE receives data from other PEs
!@+             to fill a contiguous pre-requested range of J indices
      public :: band_pack,band_pack_column
      interface band_pack
        module procedure band_pack_ij
        module procedure band_pack_ijl
      end interface
!@var BAND_PACK_TYPE a data structure needed by BAND_PACK, initialized
!@var via INIT_BAND_PACK_TYPE
      type band_pack_type
!        integer :: im_world
        integer :: j_strt,j_stop
        integer :: j_strt_halo,j_stop_halo
        integer :: jband_strt,jband_stop
        integer, dimension(:), pointer :: scnts,sdspl,sdspl_inplace
        integer, dimension(:), pointer :: rcnts,rdspl,rdspl_inplace
        integer, dimension(:), pointer :: j0_send,j1_send
        integer, dimension(:), pointer :: j0_recv,j1_recv
        integer :: mpi_comm
        integer :: npes_comm
      end type band_pack_type
!@var INIT_BAND_PACK_TYPE initialization routine during which each PE
!@+   requests a range of J indices and sets up the necessary send/receive
!@+   information for the BAND_PACK procedure
      public :: band_pack_type,init_band_pack_type

      PUBLIC SEND_TO_J
      interface SEND_TO_J
         module procedure SEND_TO_J_1D
         module procedure ISEND_TO_J_0D
      end interface

      PUBLIC RECV_FROM_J
      interface RECV_FROM_J
         module procedure RECV_FROM_J_1D
         module procedure IRECV_FROM_J_0D
      end interface


      INTEGER, PARAMETER :: HALO_WIDTH = 1
      integer ::  root


#ifndef USE_DD2D_UTILS
      ! Local grid information
      TYPE DIST_GRID

        type (Hidden_type) :: private
!!$$        type (Domain_type) :: localSubdomain
         INTEGER :: NPES_USED
         ! Parameters for Global domain
         INTEGER :: IM_WORLD        ! Number of Longitudes
         INTEGER :: JM_WORLD        ! Number of latitudes
         ! Parameters for local domain
         LOGICAL :: HAVE_DOMAIN = .false.     ! Whether this PE has any of the domain
         INTEGER :: I_STRT          ! Begin local domain longitude index
         INTEGER :: I_STOP          ! End   local domain longitude index
         INTEGER :: J_STRT          ! Begin local domain latitude  index
         INTEGER :: J_STOP          ! End   local domain latitude  index
         INTEGER :: J_STRT_SKP      ! Begin local domain exclusive of S pole
         INTEGER :: J_STOP_SKP      ! End   local domain exclusive of N pole
         INTEGER :: ni_loc ! for transpose
         ! Parameters for halo of local domain
         INTEGER :: I_STRT_HALO     ! Begin halo longitude index
         INTEGER :: I_STOP_HALO     ! End   halo longitude index
         INTEGER :: J_STRT_HALO     ! Begin halo latitude  index
         INTEGER :: J_STOP_HALO     ! End   halo latitude  index
         ! Parameters for staggered "B" grid
         ! Note that global staggered grid begins at "2".
         INTEGER :: J_STRT_STGR     ! Begin local staggered domain
         INTEGER :: J_STOP_STGR     ! End   local staggered domain

         INTEGER, DIMENSION(:), POINTER :: DJ_MAP
         INTEGER :: DJ
         
      END TYPE DIST_GRID
#endif

      type (DIST_GRID) :: GRID
!not used      TYPE (DIST_GRID) :: GRID_TRANS

      public :: haveLatitude

! Remaining variables are private to the module.

!@var NPES_WORLD number of total processes
      INTEGER :: NPES_WORLD
!@var NP_LON number of azimuthal processes.
      INTEGER :: NP_LON
!@var NP_LAT number of meridional     processes.
      INTEGER :: NP_LAT
!@var MY_PET index of _this_ PET (analagous to MPI rank)
      INTEGER :: my_pet
!@var RANK_LON index of _this_ process in azimuthal set.
      INTEGER :: RANK_LON
!@var RANK_LAT_RANK index of _this_ process in meridional set.
      INTEGER :: RANK_LAT

      INTEGER, PUBLIC :: CHECKSUM_UNIT

      ! temporary public during refactoring
      public :: isPeriodic
      public :: my_pet, npes_world
#ifdef USE_MPI
      public :: esmf_grid_my_pe_loc
      public :: esmf_grid_pe_layout
#endif
      public :: am_i_root

      ! Direction bits
      public :: NORTH, SOUTH
  Integer, Parameter :: NORTH = 2**0
  Integer, Parameter :: SOUTH = 2**1

  public :: getMpiCommunicator
  public :: getNumProcesses
  public :: getNumAllProcesses
  public :: getMpiTag
  public :: incrementMpiTag
  public :: hasSouthPole
  public :: hasNorthPole
  public :: hasPeriodicBC

  public :: getLogUnit


#ifdef USE_MPI
  type axisIndex
     sequence  ! sequence forces the data elements 
                ! to be next to each other in memory 
     integer :: min
     integer :: max
  end type axisIndex
  Public :: axisIndex 
  Public :: getAxisIndex
  interface getAxisIndex
     module procedure getMPIAxisIndex
  end interface
#endif

  integer, parameter :: maxStrLen = 40
  public :: maxStrLen

      CONTAINS

      subroutine init_app
!@sum  This routine initializes the ESMF framework (or MPI)
!@+    The initialization should proceed prior to any grid computations.
!@auth Tom Clune GSFC/SSSO/610.3
      integer :: rc
      character(len=10) :: fileName
      character(len=10) :: hostName
      character(len=80) :: command
      integer :: myUnit
      NPES_WORLD = 1                  ! default NPES = 1 for serial run
      my_pet = 0                      ! default my_pet = root PE for serial run

#ifdef USE_MPI
#ifndef USE_MPP
   call MPI_INIT(rc)
#endif
   call MPI_COMM_SIZE(MPI_COMM_WORLD, NPES_WORLD, rc)
   call MPI_COMM_RANK(MPI_COMM_WORLD, my_pet, rc)
#endif
      if(my_pet == 0) then
         write(*,*)'Num MPI Processes: ', NPES_WORLD
      end if

      NP_LON   = 1
      RANK_LON = 0
      NP_LAT   = NPES_WORLD
      RANK_LAT = my_pet

      
      return
      end subroutine init_app

      SUBROUTINE DESTROY_GRID(grd_dum)
!@sum  This routine releases resources contained within the ESMF grid
!@auth Tom Clune GSFC/SSSO/610.3
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
      END SUBROUTINE DESTROY_GRID

      SUBROUTINE INIT_GRID(grd_dum,IM,JM,LM,width,J_SCM,bc_periodic, &
    &                     CREATE_CAP,npes_max,excess_on_pe0)
!@sum  This routine initializes the ESMF grid as well as data
!@+    structures associated with modelE's domain decomposition
!@auth Tom Clune GSFC/SSSO/610.3
      USE FILEMANAGER, Only : openunit
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
      INTEGER, INTENT(IN) :: IM, JM,LM
      INTEGER, OPTIONAL, INTENT(IN) :: J_SCM ! single column model
      INTEGER, OPTIONAL :: width
      LOGICAL, OPTIONAL, INTENT(IN) :: bc_periodic
      LOGICAL, OPTIONAL, INTENT(IN) :: CREATE_CAP
      INTEGER, OPTIONAL, INTENT(IN) :: npes_max
      LOGICAL, OPTIONAL, INTENT(IN) :: excess_on_pe0
      integer, parameter :: numDims=2
      integer, dimension(numDims) :: grid_size
      integer             :: rc
      INTEGER :: J_EQUATOR
      INTEGER :: I0_DUM, I1_DUM
      INTEGER :: J0_DUM, J1_DUM
      INTEGER :: width_
      INTEGER :: pet
      INTEGER :: NTILES
      INTEGER :: AIbounds(4)
#ifdef USE_MPI
      Type (AxisIndex), Pointer :: AI(:,:)
#endif
      INTEGER :: p
      integer :: npes_used
      integer, dimension(:), allocatable :: pelist
      integer :: group_world,group_used,ierr
      integer :: newCommunicator

      logical :: useCubedSphere

#ifdef CUBED_SPHERE
      useCubedSphere = .true.
#else
      useCubedSphere = .false.
#endif

      if (useCubedSphere) then
        grid_size = (/IM, JM*6/)
!!$$        grd_dum%localSubdomain = newDomain(IM, 6*JM)
      else
        grid_size = (/IM, JM/)
!!$$        grd_dum%localSubdomain = newDomain(IM, JM)
      end if

      grd_dum%IM_WORLD      = IM
      grd_dum%JM_WORLD      = JM

      npes_used = 1
      grd_dum%npes_used = npes_used

#ifdef USE_MPI
      Allocate(grd_dum%dj_map(0:npes_world-1))
      npes_used = min(npes_world, jm-2) ! jm-2 is latlon-specific
      if (present(npes_max)) npes_used = min(npes_max, npes_used)
      grd_dum%npes_used = npes_used

      AIbounds = MPIgridBounds(grd_dum)
      I0_DUM = AIbounds(1)
      I1_DUM = AIbounds(2)
      J0_DUM = AIbounds(3)
      J1_DUM = AIbounds(4)
#else
      I0_DUM = 1
      I1_DUM = IM
      J0_DUM = 1
      J1_DUM = JM

      if (present(J_SCM)) then
         J0_DUM = J_SCM
         J1_DUM = J_SCM
      end if
#endif

      width_ = HALO_WIDTH
      If (Present(width)) width_=width

      ! Wrapped grid
      grd_dum%I_STRT        = I0_DUM
      grd_dum%I_STOP        = I1_DUM
      grd_dum%I_STRT_HALO   = MAX( 1, I0_DUM-width_)
      grd_dum%I_STOP_HALO   = MIN(IM, I1_DUM+width_)
      grd_dum%ni_loc = (RANK_LAT+1)*IM/NPES_used - RANK_LAT*IM/NPES_used

      grd_dum%j_strt        = J0_DUM
      grd_dum%J_STOP        = J1_DUM
#ifdef USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      grd_dum%HAVE_DOMAIN   = J0_DUM <= JM

      grd_dum%J_STRT_SKP = max (   2, J0_DUM)
      grd_dum%J_STOP_SKP = min (JM-1, J1_DUM)

#ifdef USE_MPI
      grd_dum%J_STRT_HALO   = J0_DUM - width_
      grd_dum%J_STOP_HALO   = J1_DUM + width_
      grd_dum%private%numProcesses = npes_used
      grd_dum%private%numAllProcesses = npes_world
      grd_dum%private%mpi_tag = 10  ! initial value

! Create a new MPI communicator including all the PEs with a nonzero
! domain size.  Even when NPES_USED == NPES_WORLD, this is convenient for
! avoiding collisions of MPI tag sequences.

      call mpi_comm_group(MPI_COMM_WORLD,group_world,ierr)
      allocate(pelist(0:npes_used-1))
      do p=0,npes_used-1
        pelist(p) = p
      enddo
      call mpi_group_incl(group_world,npes_used,pelist,group_used,ierr)
      deallocate(pelist)
      call mpi_comm_create(MPI_COMM_WORLD,group_used, newCommunicator, ierr)
      if(.not. grd_dum%HAVE_DOMAIN) newCommunicator = MPI_COMM_NULL
      call setMpiCommunicator(grd_dum, newCommunicator)
#else
      ! I guess we don't need HALO in SCM mode...
      !grd_dum%J_STRT_HALO = MAX(1,  grd_dum % J_STRT - 1)
      !grd_dum%J_STOP_HALO = MIN(JM, grd_dum % J_STOP + 1)
      grd_dum%J_STRT_HALO = MAX(1,  grd_dum % J_STRT)
      grd_dum%J_STOP_HALO = MIN(JM, grd_dum % J_STOP)
#endif

      grd_dum%J_STRT_STGR   = max(2,J0_DUM)
      grd_dum%J_STOP_STGR   = J1_DUM

      grd_dum%private%hasSouthPole = J0_DUM == 1  !(RANK_LAT == 0)
      grd_dum%private%hasNorthPole = J1_DUM == JM .and. J0_DUM <= JM !(RANK_LAT == NP_LAT - 1) &

      J_EQUATOR = JM/2
      grd_dum%private%hasEquator =  (J0_DUM <= J_EQUATOR) .AND. (J1_DUM >= J_EQUATOR)

#ifdef USE_DD2D_UTILS
! need to initialize the dd2d version of dist_grid for I/O
           call init_dist_grid( &
    &     grd_dum%IM_WORLD,grd_dum%JM_WORLD,1,  &
    &     grd_dum%I_STRT,grd_dum%I_STOP, &
    &     grd_dum%j_strt,grd_dum%J_STOP, &
    &     grd_dum%I_STRT_HALO,grd_dum%I_STOP_HALO, &
    &     grd_dum%J_STRT_HALO,grd_dum%J_STOP_HALO,grd_dum)
#endif

      if (present(J_SCM)) then
        ! assume J_SCM is in "general position"
        grd_dum%private%hasSouthPole = .false.
        grd_dum%private%hasNorthPole = .false.
        grd_dum%private%hasEquator    = .false.
      endif

      grd_dum % private%PERIODICBC = isPeriodic(bc_periodic)

     ! assumption: decomposition along "east-west" direction
     ! not used for lat-lon grids
      if(grd_dum % private%PERIODICBC) then
        grd_dum%private%hasSouthPole = .false.
        grd_dum%private%hasNorthPole = .false.
        grd_dum%private%hasEquator    = .false.
      endif

      ! set lookup table PET(J)
      Allocate(grd_dum%private%lookup_pet(1:JM))
      grd_dum%private%lookup_pet(:) = 0

#ifdef USE_MPI
      Allocate(AI(0:NPES_WORLD-1,3))
      call getAxisIndex(grd_dum, AI)

      Do p = 1, npes_world
        grd_dum%private%lookup_pet( AI(p-1,2)%min : AI(p-1,2)%max ) = p-1
      End Do

      do p=0, npes_world-1
         grd_dum%dj_map(p)=AI(p,2)%max - AI(p,2)%min + 1
      end do
      grd_dum%dj=grd_dum%dj_map(my_pet)

      Deallocate(AI)
#endif

      END SUBROUTINE INIT_GRID

      SUBROUTINE GET(grd_dum, I_STRT, I_STOP, &
    &                        I_STRT_HALO, I_STOP_HALO, &
    &                        J_STRT, J_STOP, J_STRT_HALO, J_STOP_HALO, &
    &                        J_STRT_SKP, J_STOP_SKP, &
    &                        J_STRT_STGR, J_STOP_STGR, &
    &                        have_south_pole, have_north_pole)
!@sum  This routine provides general access to various values that characterize the local 
!@+ domain bounds on this process. 
!@auth Tom Clune GSFC/SSSO/610.3
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      INTEGER, OPTIONAL :: I_STRT, I_STOP
      INTEGER, OPTIONAL :: I_STRT_HALO, I_STOP_HALO
      INTEGER, OPTIONAL :: J_STRT, J_STOP
      INTEGER, OPTIONAL :: J_STRT_HALO, J_STOP_HALO
      INTEGER, OPTIONAL :: J_STRT_SKP, J_STOP_SKP
      INTEGER, OPTIONAL :: J_STRT_STGR, J_STOP_STGR
      LOGICAL, OPTIONAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      IF (PRESENT(I_STRT)) I_STRT = grd_dum%I_STRT
      IF (PRESENT(I_STOP)) I_STOP = grd_dum%I_STOP

      IF (PRESENT(I_STRT_HALO)) I_STRT_HALO = grd_dum%I_STRT_HALO
      IF (PRESENT(I_STOP_HALO)) I_STOP_HALO = grd_dum%I_STOP_HALO

      IF (PRESENT(J_STRT)) J_STRT = grd_dum%j_strt
      IF (PRESENT(J_STOP)) J_STOP = grd_dum%J_STOP

      IF (PRESENT(J_STRT_HALO)) J_STRT_HALO = grd_dum%J_STRT_HALO
      IF (PRESENT(J_STOP_HALO)) J_STOP_HALO = grd_dum%J_STOP_HALO

      IF (PRESENT(J_STRT_SKP)) J_STRT_SKP = grd_dum%J_STRT_SKP
      IF (PRESENT(J_STOP_SKP)) J_STOP_SKP = grd_dum%J_STOP_SKP

      IF (PRESENT(J_STRT_STGR)) J_STRT_STGR = grd_dum%J_STRT_STGR
      IF (PRESENT(J_STOP_STGR)) J_STOP_STGR = grd_dum%J_STOP_STGR

      IF (PRESENT(HAVE_SOUTH_POLE)) &
    &             HAVE_SOUTH_POLE= grd_dum%private%hasSouthPole
      IF (PRESENT(HAVE_NORTH_POLE)) &
    &             HAVE_NORTH_POLE= grd_dum%private%hasNorthPole

      END SUBROUTINE GET


      SUBROUTINE FINISH_APP()
!@sum  This routine finalizes the ESMF framework
!@auth Tom Clune GSFC/SSSO/610.3
      USE FILEMANAGER, ONLY : closeunit
      IMPLICIT NONE

      INTEGER :: rc

#ifdef DEBUG_DECOMP
      CALL closeunit(CHECKSUM_UNIT)
      CALL closeunit(grid%private%log_unit)
#endif

#ifdef USE_MPI
   call MPI_Finalize(rc)
#endif

      END SUBROUTINE FINISH_APP

      SUBROUTINE SUMXPE_1D(arr, arr_master, increment)
!@sum  sum a 1-D array over processors without reducing its rank.
!@+ Only the root process receives the final result
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      REAL*8, DIMENSION(:) :: arr
      REAL*8, DIMENSION(:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_MPI
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root()) then
           allocate(arr_tmp(arr_size))
        else
           allocate(arr_tmp(1))
        end if
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + arr_tmp
        endif
        deallocate(arr_tmp)
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr is overwritten by itself after reduction
#ifdef USE_MPI
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_1D

      SUBROUTINE SUMXPE_1D_I(arr, arr_master, increment)
!@sum  sum a 1-D integer array over processors without reducing its rank
!@+ Only the root process receives the final result
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      INTEGER, DIMENSION(:) :: arr
      INTEGER, DIMENSION(:), optional :: arr_master
      logical, intent(in), optional :: increment
      INTEGER, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_MPI
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root()) allocate(arr_tmp(arr_size))
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_INTEGER, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + arr_tmp
          deallocate(arr_tmp)
        endif
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_INTEGER, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr is overwritten by itself after reduction
#ifdef USE_MPI
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size, &
    &        MPI_INTEGER,MPI_SUM,root, &
    &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_1D_I

      SUBROUTINE SUMXPE_2D(arr, arr_master, increment)
!@sum  sum a 2-D array over processors without reducing its rank
!@+ Only the root process receives the final result
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:) :: arr
      REAL*8, DIMENSION(:,:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_MPI
         arr_size = size(arr)
         if(increment_) then
            if(am_i_root()) then
               allocate(arr_tmp(arr_size))
            else
               allocate(arr_tmp(1))
            end if
            call mpi_reduce(arr,arr_tmp,arr_size, &
    &           MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &           MPI_COMM_WORLD, ierr)
            if(am_i_root()) then
              arr_master = arr_master + reshape(arr_tmp,shape(arr))
            endif
            deallocate(arr_tmp)
         else
            call mpi_reduce(arr,arr_master,arr_size, &
    &           MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &           MPI_COMM_WORLD, ierr)
         endif
#else 
         if(increment_) then
            arr_master = arr_master + arr
         else
            arr_master = arr
         endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr is overwritten by itself after reduction
#ifdef USE_MPI
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_2D

      SUBROUTINE SUMXPE_3D(arr, arr_master, increment)
!@sum  sum a 3-D array over processors without reducing its rank
!@+ Only the root process receives the final result
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:,:) :: arr
      REAL*8, DIMENSION(:,:,:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_MPI
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root()) then
           allocate(arr_tmp(arr_size))
        else
           allocate(arr_tmp(1))
        end if
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + reshape(arr_tmp,shape(arr))
        endif
        deallocate(arr_tmp)
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr  is overwritten by itself after reduction
#ifdef USE_MPI
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_3D

      SUBROUTINE SUMXPE_4D(arr, arr_master, increment)
!@sum  sum a 4-D array over processors without reducing its rank
!@+ Only the root process receives the final result
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:,:,:) :: arr
      REAL*8, DIMENSION(:,:,:,:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_MPI
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root()) then
           allocate(arr_tmp(arr_size))
        else
           allocate(arr_tmp(1))
        end if
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + reshape(arr_tmp,shape(arr))
        endif
        deallocate(arr_tmp)
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr  is overwritten by itself after reduction
#ifdef USE_MPI
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_4D




    !---------------------------
#ifdef USE_MPI
      subroutine ESMF_GRID_PE_LAYOUT  (GRID, NX, NY)
!@sum  Get the process topology associated with GRID. This implementation only
!@+ support 1-D decomposition, hence NX=1.
!@auth Tom Clune GSFC/SSSO/610.3
        type (Dist_Grid), intent(IN) :: grid
        integer, intent(OUT)         :: NX, NY

        NX = 1
        NY = getNumProcesses(grid)

      end subroutine ESMF_GRID_PE_LAYOUT

    !---------------------------

      subroutine ESMF_GRID_MY_PE_LOC  (GRID,NX0,NY0)
!@sum  Get the process ranks in lon/lat directions from GRID
!@auth Tom Clune GSFC/SSSO/610.3
        type (Dist_Grid), intent(IN) :: grid
        integer, intent(OUT)          :: NX0, NY0

        NX0 = 0
        NY0 = my_pet

      end subroutine ESMF_GRID_MY_PE_LOC

    !---------------------------
#endif

      SUBROUTINE HERE(file, line)
!@sum a helper routine used to print debugging information
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      Character(Len=*) :: file
      Integer :: line

      INTEGER :: ierr
      Integer, Allocatable :: lines(:)

#ifdef DEBUG_DECOMP
      CALL LOG_PARALLEL(grid, file, line)
      If (AM_I_ROOT()) Then
         WRITE(CHECKSUM_UNIT,*)'HERE: ',file, line
         CALL SYS_FLUSH(CHECKSUM_UNIT)
       End If
#ifdef USE_MPI
       ALLOCATE(lines(npes_world))
       Call MPI_Allgather(line, 1, MPI_INTEGER, lines, 1, MPI_INTEGER, &
    &      MPI_COMM_WORLD, ierr)
       If (Any(lines /= line)) &
    &      call stop_model('HERE: synchronization error -severe.',255)
       Deallocate(lines)
#endif
#endif
#ifdef USE_MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      END SUBROUTINE HERE

      SUBROUTINE LOG_PARALLEL(grd_dum, file, line, i0, i1, x0, x1)
!@sum prints debugging information to a file
!@auth Tom Clune GSFC/SSSO/610.3
      Use FILEMANAGER, only : nameunit
      IMPLICIT NONE

      TYPE(DIST_GRID), INTENT(IN) :: grd_dum
      CHARACTER(Len=*) :: file
      INTEGER          :: line

      INTEGER, OPTIONAL :: i0, i1(:)
      REAL*8, OPTIONAL  :: x0, x1(:)

      INTEGER :: iu
      INTEGER :: n

#ifdef DEBUG_DECOMP
      iu = grd_dum%private%log_unit
      WRITE(iu, *) file, line
      If (PRESENT(i0)) WRITE(iu, *) '   i0=',i0
      If (PRESENT(x0)) WRITE(iu, *) '   x0=',x0
      IF (PRESENT(i1)) THEN
         DO n = 1, Size(i1)
            WRITE(iu, '(10x,i4,1x,a,i6)')n, '   i1=',i1(n)
         END DO
      END IF
      IF (PRESENT(x1)) THEN
         DO n = 1, Size(x1)
            WRITE(iu, '(10x,i4,1x,a,e22.16)')n, '   x1=',x1(n)
         END DO
      END IF
      CALL SYS_FLUSH(iu)
#endif

      END SUBROUTINE LOG_PARALLEL

      SUBROUTINE GLOBALMIN_R(grd_dum, val, val_min)
!@sum determine min value of real variable across processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_min

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_min, 1, MPI_DOUBLE_PRECISION,MPI_MIN, &
    &     getMpiCommunicator(GRD_DUM), ierr)
#else
      val_min = val
#endif

      END SUBROUTINE

      SUBROUTINE GLOBALMAX_R(grd_dum, val, val_max)
!@sum determine max value of real variable across processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_max

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_max, 1, MPI_DOUBLE_PRECISION,MPI_MAX, &
    &     getMpiCommunicator(grd_dum), ierr)
#else
      val_max = val
#endif

      END SUBROUTINE

      SUBROUTINE GLOBALMAX_I(grd_dum, val, val_max)
!@sum  determine max value of integer variable across processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      INTEGER,            INTENT(IN)  :: val
      INTEGER,            INTENT(OUT) :: val_max

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_max, 1, MPI_INTEGER, MPI_MAX, &
    &     getMpiCommunicator(grd_dum), ierr)
#else
      val_max = val
#endif

      END SUBROUTINE

      SUBROUTINE GLOBALMAX_I_1D(grd_dum, val, val_max)
!@sum determine max value of integer 1-D array across processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      INTEGER,            INTENT(IN)  :: val(:)
      INTEGER,            INTENT(OUT) :: val_max(:)

      INTEGER  :: n,ierr

#ifdef USE_MPI
      n = size(val)
      CALL MPI_Allreduce(val, val_max, n, MPI_INTEGER, MPI_MAX, &
    &     getMpiCommunicator(grd_dum), ierr)
#else
      val_max(:) = val(:)
#endif

      END SUBROUTINE


      SUBROUTINE INIT_BAND_PACK_TYPE(grd_src, grd_dst, band_j0,band_j1, &
    &     bandpack)
!@sum initialize the bandpack derived type with the information needed for 
!@+ the band_pack procedure to fill output arrays with data from J indices 
!@+ band_j0 to band_j1
!@auth modelE development team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_src,grd_dst
      INTEGER, INTENT(IN) :: band_j0,band_j1
      TYPE (BAND_PACK_TYPE), intent(OUT) :: bandpack
#ifdef USE_MPI
      integer, dimension(0:npes_world-1) :: &
    &     j0_have,j1_have,j0_requested,j1_requested
      integer :: p, ierr, im,jm, j0send,j1send,j0recv,j1recv,npes

!
! Store some general MPI info (NOTE: for now we assume that the process
! set of grd_src is a subset of that in grd_dst, or vice versa, which
! allows us to take this info from the larger set).
!
      if(getNumProcesses(grd_src) > getNumProcesses(grd_dst)) then
        npes = getNumProcesses(grd_src)
        bandpack%mpi_comm = getMpiCommunicator(grd_src)
      else
        npes = getNumProcesses(grd_dst)
        bandpack%mpi_comm = getMpiCommunicator(grd_dst)
      endif
      bandpack%npes_comm = npes
      allocate(bandpack%j0_send(0:npes-1), &
    &         bandpack%j1_send(0:npes-1), &
    &         bandpack%j0_recv(0:npes-1), &
    &         bandpack%j1_recv(0:npes-1), &
    &         bandpack%scnts(0:npes-1), &
    &         bandpack%sdspl(0:npes-1), &
    &         bandpack%rcnts(0:npes-1), &
    &         bandpack%rdspl(0:npes-1), &
    &         bandpack%sdspl_inplace(0:npes-1), &
    &         bandpack%rdspl_inplace(0:npes-1) &
    &     )
#endif
!      bandpack%im_world = grd_src%im_world
      bandpack%j_strt = grd_src%j_strt
      bandpack%j_stop = grd_src%j_stop
      bandpack%j_strt_halo = grd_src%j_strt_halo
      bandpack%j_stop_halo = grd_src%j_stop_halo
      bandpack%jband_strt = band_j0
      bandpack%jband_stop = band_j1
#ifdef USE_MPI
      im = 1!grd_src%im_world
      jm = grd_src%jm_world
!
! Set up the MPI send/receive information
!
      call mpi_allgather(grd_src%j_strt,1,MPI_INTEGER,j0_have,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      call mpi_allgather(grd_src%j_stop,1,MPI_INTEGER,j1_have,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      call mpi_allgather(band_j0,1,MPI_INTEGER,j0_requested,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      call mpi_allgather(band_j1,1,MPI_INTEGER,j1_requested,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      do p=0,npes-1
        j0send = max(grd_src%j_strt,j0_requested(p))
        j1send = min(grd_src%j_stop,j1_requested(p))
        bandpack%j0_send(p) = j0send
        bandpack%j1_send(p) = j1send
        if(j0send <= j1send) then
          bandpack%scnts(p) = im*(j1send-j0send+1)
          bandpack%sdspl_inplace(p) = im*(j0send-grd_src%j_strt_halo)
        else
          bandpack%scnts(p) = 0
          bandpack%sdspl_inplace(p) = 0
        endif
        j0recv = max(j0_have(p),band_j0)
        j1recv = min(j1_have(p),band_j1)
        bandpack%j0_recv(p) = j0recv
        bandpack%j1_recv(p) = j1recv
        if(j0recv <= j1recv) then
          bandpack%rcnts(p) = im*(j1recv-j0recv+1)
          bandpack%rdspl_inplace(p) = im*(j0recv-band_j0)
        else
          bandpack%rcnts(p) = 0
          bandpack%rdspl_inplace(p) = 0
        endif
      enddo
      bandpack%rdspl(0) = 0
      bandpack%sdspl(0) = 0
      do p=1,npes-1
        bandpack%sdspl(p) = bandpack%sdspl(p-1)+bandpack%scnts(p-1)
        bandpack%rdspl(p) = bandpack%rdspl(p-1)+bandpack%rcnts(p-1)
      enddo
#endif
      RETURN
      END SUBROUTINE INIT_BAND_PACK_TYPE

      SUBROUTINE BAND_PACK_ij(bandpack,ARR,ARR_band)
!@sum use bandpack derived type to fill 2-D output array
!@auth modelE development team

!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,bandpack%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,bandpack%jband_strt:)
#ifdef USE_MPI
      integer, dimension(0:npes_world-1) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr
      integer :: im,npes
      npes = bandpack%npes_comm
      im = size(arr,1)
      scnts(0:npes-1) = im*bandpack%scnts
      sdspl(0:npes-1) = im*bandpack%sdspl_inplace
      rcnts(0:npes-1) = im*bandpack%rcnts
      rdspl(0:npes-1) = im*bandpack%rdspl_inplace
      call mpi_alltoallv(arr, scnts, sdspl, mpi_double_precision, &
    &                   arr_band, rcnts, rdspl, mpi_double_precision, &
    &                   bandpack%mpi_comm, ierr)
#else
      arr_band(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP) = &
    &     arr(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_ij

      SUBROUTINE BAND_PACK_ijl(bandpack,ARR,ARR_band)
!@sum use bandpack derived type to fill 3-D (ijl) output array
!@auth modelE development team

!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,bandpack%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,bandpack%jband_strt:,:)
#ifdef USE_MPI
      integer, dimension(0:npes_world-1) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr,im,lm,i,j,l,n,p
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      integer :: npes
      npes = bandpack%npes_comm
      im = size(arr,1)
      lm = size(arr,3)
      scnts = im*lm*bandpack%scnts
      sdspl = im*lm*bandpack%sdspl
      rcnts = im*lm*bandpack%rcnts
      rdspl = im*lm*bandpack%rdspl
      allocate(bufsend(sum(scnts)),bufrecv(sum(rcnts)))
      n = 0
      do p=0,npes-1
        do l=1,lm
          do j=bandpack%j0_send(p),bandpack%j1_send(p)
            do i=1,im
              n = n + 1
              bufsend(n) = arr(i,j,l)
            enddo
          enddo
        enddo
      enddo
      call mpi_alltoallv(bufsend, scnts, sdspl, mpi_double_precision, &
    &                   bufrecv, rcnts, rdspl, mpi_double_precision, &
    &                   bandpack%mpi_comm, ierr)
      n = 0
      do p=0,npes-1
        do l=1,lm
          do j=bandpack%j0_recv(p),bandpack%j1_recv(p)
            do i=1,im
              n = n + 1
              arr_band(i,j,l) = bufrecv(n)
            enddo
          enddo
        enddo
      enddo
      deallocate(bufsend,bufrecv)
#else
      arr_band(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP,:) = &
    &     arr(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP,:)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_ijl

      SUBROUTINE BAND_PACK_COLUMN(bandpack,ARR,ARR_band)
!@sum use bandpack derived type to fill column output array
!@auth modelE development team

!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,:,bandpack%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,:,bandpack%jband_strt:)
#ifdef USE_MPI
      integer, dimension(0:npes_world-1) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr,im,lm
      integer :: npes
      npes = bandpack%npes_comm
      im = size(arr,2)
      lm = size(arr,1)
      scnts(0:npes-1) = im*lm*bandpack%scnts
      sdspl(0:npes-1) = im*lm*bandpack%sdspl_inplace
      rcnts(0:npes-1) = im*lm*bandpack%rcnts
      rdspl(0:npes-1) = im*lm*bandpack%rdspl_inplace
      call mpi_alltoallv(arr, scnts, sdspl, mpi_double_precision, &
    &                   arr_band, rcnts, rdspl, mpi_double_precision, &
    &                   bandpack%mpi_comm, ierr)
#else
      arr_band(:,:,bandpack%JBAND_STRT:bandpack%JBAND_STOP) = &
    &     arr(:,:,bandpack%JBAND_STRT:bandpack%JBAND_STOP)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_COLUMN

      SUBROUTINE ESMF_BCAST_0D(grd_dum, arr)
!@sum routine to broadcast real value to all processors
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,1,MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(grd_dum), ierr)
#endif

      END SUBROUTINE ESMF_BCAST_0D

      SUBROUTINE ESMF_BCAST_1D(grd_dum, arr)
!@sum routine to broadcast 1-D array data to all processors
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(grd_dum), ierr)
#endif

      END SUBROUTINE ESMF_BCAST_1D

      SUBROUTINE ESMF_BCAST_2D(grd_dum, arr)
!@sum routine to broadcast 2-D array data to all PEs.
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(grd_dum), ierr)
#endif

      END SUBROUTINE ESMF_BCAST_2D

      SUBROUTINE ESMF_BCAST_3D(grd_dum, arr)
!@sum routine to broadcast 3-D array data to all processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:,:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(grd_dum), ierr)
#endif

      END SUBROUTINE ESMF_BCAST_3D

      SUBROUTINE ESMF_BCAST_4D(grd_dum, arr)
!@sum routine to broadcast 4-D array data to all processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:,:,:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(grd_dum), ierr)
#endif

      END SUBROUTINE ESMF_BCAST_4D

      SUBROUTINE ESMF_IBCAST_0D(grd_dum, arr)
!@sum routine to broadcast integer variable data to all processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,1,MPI_INTEGER,root, &
    &     getMpiCommunicator(grd_dum), ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_0D

      SUBROUTINE ESMF_IBCAST_1D(grd_dum, arr)
!@sum routine to broadcast 1-D array integer data to all processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER, root, &
    &     getMpiCommunicator(grd_dum), ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_1D

      SUBROUTINE ESMF_IBCAST_2D(grd_dum, arr)
!@sum routine to broadcast 2-D array integer data to all processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER ,root, &
    &     getMpiCommunicator(grd_dum), ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_2D

      SUBROUTINE ESMF_IBCAST_3D(grd_dum, arr)
!@sum  routine to broadcast 3-D array integer data to all processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:,:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER, root, &
    &     getMpiCommunicator(grd_dum), ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_3D

      SUBROUTINE ESMF_IBCAST_4D(grd_dum, arr)
!@sum  routine to broadcast 4-D array integer data to all processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:,:,:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER ,root, &
    &     getMpiCommunicator(grd_dum), ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_4D


      SUBROUTINE TRANSPOSE_ijk(grid, x_in, x_out, reverse)
!@sum  compute transpose of a decomposed 3-D array
!@auth Tom Clune GSFC/SSSO/610.3
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x_in(:,grid%J_STRT_HALO:,:)
      REAL*8 :: x_out(:,:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER :: I0(0:NPES_WORLD-1), I1(0:NPES_WORLD-1)
      INTEGER :: J0(0:NPES_WORLD-1), J1(0:NPES_WORLD-1)
      REAL*8, ALLOCATABLE :: sbuf(:), rbuf(:)
#ifdef USE_MPI
      TYPE (AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ,nk,k
      INTEGER :: ierr, p, rc
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt, npes
      INTEGER, DIMENSION(0:NPES_WORLD-1) :: scnts, rcnts, sdspl, rdspl
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_MPI
      If (reverse_) Then
         X_IN(:,1:grid%JM_WORLD,:) = X_OUT(:,1:grid%JM_WORLD,:)
      Else
         X_OUT(:,1:grid%JM_WORLD,:) = X_IN(:,1:grid%JM_WORLD,:)
      End If
#else
      npes = getNumProcesses(grid)

      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:npes_world-1,3))
      call getAxisIndex(grid, AI)

      DO p = 0, npes - 1
         J0(p) = AI(p,2)%min
         J1(p) = AI(p,2)%max
      END DO
      DEALLOCATE(AI)

      ni_loc = I1(my_pet) - I0(my_pet) + 1
      nj_loc = J1(my_pet) - J0(my_pet) + 1

      nk = SIZE(X_IN,3)

      ALLOCATE(rbuf(grid%JM_WORLD * ni_loc * nk))
      ALLOCATE(sbuf(grid%IM_WORLD * nj_loc * nk))

      sdspl(0) = 0
      rdspl(0) = 0
      icnt = 0
      DO p = 0, npes -1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  rbuf(icnt) = X_out(i,j,k)
               END DO
            END DO
         ELSE
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  sbuf(icnt) = X_in(i,j,k)
               End Do
            END DO
         END IF
         nip = I1(p) - I0(p) + 1
         njp = J1(p) - J0(p) + 1
         scnts(p) = nj_loc * nip * nk
         rcnts(p) = ni_loc * njp * nk
         If (p > 0) sdspl(p) = sdspl(p-1) + scnts(p-1)
         If (p > 0) rdspl(p) = rdspl(p-1) + rcnts(p-1)
       END DO
      END DO

      If (reverse_) Then
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      End If

      icnt = 0
      DO p = 0, npes - 1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  X_in(i,j,k) = sbuf(icnt)
               End Do
            END DO
         Else
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  X_out(i,j,k) = rbuf(icnt)
               END DO
            END DO
         End If
       END DO
      END DO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)
#endif

      END SUBROUTINE TRANSPOSE_ijk

      SUBROUTINE TRANSPOSE_ij(grid, x_in, x_out, reverse)
!@sum  compute transpose of a decomposed 2-D array
!@auth Tom Clune GSFC/SSSO/610.3
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x_in(:,grid%J_STRT_HALO:)
      REAL*8 :: x_out(:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER :: I0(0:NPES_WORLD-1), I1(0:NPES_WORLD-1)
      INTEGER :: J0(0:NPES_WORLD-1), J1(0:NPES_WORLD-1)
      REAL*8, ALLOCATABLE :: sbuf(:), rbuf(:)
#ifdef USE_MPI
      TYPE (AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ
      INTEGER :: ierr, p, rc
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt, npes
      INTEGER, DIMENSION(0:NPES_WORLD-1) :: scnts, rcnts, sdspl, rdspl
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_MPI
      If (reverse_) Then
         X_IN(:,1:grid%JM_WORLD) = X_OUT(:,1:grid%JM_WORLD)
      Else
         X_OUT(:,1:grid%JM_WORLD) = X_IN(:,1:grid%JM_WORLD)
      End If
#else
      npes = getNumProcesses(grid)

      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:npes_world-1,3))
      call getAxisIndex(grid, AI)

      DO p = 0, npes - 1
         J0(p) = AI(p,2)%min
         J1(p) = AI(p,2)%max
      END DO
      DEALLOCATE(AI)

      ni_loc = I1(my_pet) - I0(my_pet) + 1
      nj_loc = J1(my_pet) - J0(my_pet) + 1


      ALLOCATE(rbuf(grid%JM_WORLD * ni_loc ))
      ALLOCATE(sbuf(grid%IM_WORLD * nj_loc ))

      sdspl(0) = 0
      rdspl(0) = 0
      icnt = 0
      DO p = 0, npes -1
         If (reverse_) Then
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  rbuf(icnt) = X_out(i,j)
               END DO
            END DO
         ELSE
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  sbuf(icnt) = X_in(i,j)
               End Do
            END DO
         END IF
         nip = I1(p) - I0(p) + 1
         njp = J1(p) - J0(p) + 1
         scnts(p) = nj_loc * nip
         rcnts(p) = ni_loc * njp
         If (p > 0) sdspl(p) = sdspl(p-1) + scnts(p-1)
         If (p > 0) rdspl(p) = rdspl(p-1) + rcnts(p-1)
      END DO

      If (reverse_) Then
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      End If

      icnt = 0
      DO p = 0, npes - 1
         If (reverse_) Then
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  X_in(i,j) = sbuf(icnt)
               End Do
            END DO
         Else
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  X_out(i,j) = rbuf(icnt)
               END DO
            END DO
         End If
      END DO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)
#endif

      END SUBROUTINE TRANSPOSE_ij




      SUBROUTINE TRANSPOSE_COLUMN(grid, x, x_tr, reverse)
!@sum  compute transpose of a decomposed column (1-D) array
!@auth Tom Clune GSFC/SSSO/610.3
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x(:,:,grid%J_STRT_HALO:,:)
      REAL*8 :: x_tr(:,:,:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER :: I0(0:NPES_WORLD-1), I1(0:NPES_WORLD-1)
      INTEGER :: J0(0:NPES_WORLD-1), J1(0:NPES_WORLD-1)
      REAL*8, ALLOCATABLE :: sbuf(:,:), rbuf(:,:)
#ifdef USE_MPI
      TYPE (AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ,k
      INTEGER :: ierr, p, rc
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt, npes
      INTEGER, DIMENSION(0:NPES_WORLD-1) :: scnts, rcnts, sdspl, rdspl
      INTEGER :: n, nk
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_MPI
      If (reverse_) Then
         X(:,:,1:grid%JM_WORLD,:) = X_TR(:,:,1:grid%JM_WORLD,:)
      Else
         X_TR(:,:,1:grid%JM_WORLD,:) = X(:,:,1:grid%JM_WORLD,:)
      End If
#else
      npes = getNumProcesses(grid)

      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:npes_world-1,3))
      call getAxisIndex(grid, AI)

      DO p = 0, npes - 1
         J0(p) = AI(p,2)%min
         J1(p) = AI(p,2)%max
      END DO
      DEALLOCATE(AI)

      ni_loc = I1(my_pet) - I0(my_pet) + 1
      nj_loc = J1(my_pet) - J0(my_pet) + 1

      n  = SIZE(X, 1)
      nk = SIZE(X,4)

      ALLOCATE(rbuf(n, grid%JM_WORLD * ni_loc * nk))
      ALLOCATE(sbuf(n, grid%IM_WORLD * nj_loc * nk))

      sdspl(0) = 0
      rdspl(0) = 0
      icnt = 0
      DO p = 0, npes -1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  rbuf(:,icnt) = X_tr(:,i,j,k)
               END DO
            END DO
         ELSE
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  sbuf(:,icnt) = X(:,i,j,k)
               End Do
            END DO
         END IF
         nip = I1(p) - I0(p) + 1
         njp = J1(p) - J0(p) + 1
         scnts(p) = n * nj_loc * nip * nk
         rcnts(p) = n * ni_loc * njp * nk
         If (p > 0) sdspl(p) = sdspl(p-1) + scnts(p-1)
         If (p > 0) rdspl(p) = rdspl(p-1) + rcnts(p-1)
       END DO
      END DO

      If (reverse_) Then
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      End If

      icnt = 0
      DO p = 0, npes - 1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  X(:,i,j,k) = sbuf(:,icnt)
               End Do
            END DO
         Else
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  X_tr(:,i,j,k) = rbuf(:,icnt)
               END DO
            END DO
         End If
       END DO
      END DO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)
#endif

      END SUBROUTINE TRANSPOSE_COLUMN


      logical function haveLatitude(grd_dum, j)
!@sum returns true if latitude j is in the local domain
!@auth Tom Clune GSFC/SSSO/610.3
      type (DIST_GRID), intent(in) :: grd_dum
      integer, intent(in) :: j

      haveLatitude = (j >= grd_dum%j_strt .and. j <= grd_dum%J_STOP)

      end function haveLatitude


      subroutine SEND_TO_J_1D(grd_dum, arr, j_dest, tag)
!@sum Custom routine for ENT component used for distributing initialization data across
!@+ procesess.
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(In) :: arr(:)
      Integer, Intent(In) :: j_dest, tag
      INTEGER :: ierr
#ifdef USE_MPI
      call MPI_Send(arr, Size(arr), MPI_DOUBLE_PRECISION, &
    &     grd_dum%private%lookup_pet(j_dest), tag, MPI_COMM_WORLD, ierr)
#endif
      end subroutine SEND_TO_J_1D

      subroutine ISEND_TO_J_0D(grd_dum, arr, j_dest, tag)
!@sum Custom routine for ENT component used for distributing initialization data across
!@+ procesess.
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(In) :: arr
      Integer, Intent(In) :: j_dest, tag
      INTEGER :: ierr
#ifdef USE_MPI
      call MPI_Send(arr, 1, MPI_INTEGER, &
    &     grd_dum%private%lookup_pet(j_dest), tag, MPI_COMM_WORLD, ierr)
#endif
      end subroutine ISEND_TO_J_0D

      subroutine RECV_FROM_J_1D(grd_dum, arr, j_src, tag)
!@sum Custom routine for ENT component used for distributing initialization data across
!@+ procesess.
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(In) :: arr(:)
      Integer, Intent(In) :: j_src, tag
#ifdef USE_MPI
      INTEGER :: ierr, status(MPI_STATUS_SIZE)
      call MPI_Recv(arr, Size(arr), MPI_DOUBLE_PRECISION, &
    &     grd_dum%private%lookup_pet(j_src), tag, MPI_COMM_WORLD, status, ierr)
#endif
      end subroutine RECV_FROM_J_1D

      subroutine IRECV_FROM_J_0D(grd_dum, arr, j_src, tag)
!@sum Custom routine for ENT component used for distributing initialization data across
!@+ procesess.
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(In) :: arr
      Integer, Intent(In) :: j_src, tag
#ifdef USE_MPI
      INTEGER :: ierr, status(MPI_STATUS_SIZE)
      call MPI_Recv(arr, 1, MPI_INTEGER, &
    &     grd_dum%private%lookup_pet(j_src), tag, MPI_COMM_WORLD, status, ierr)
#endif
      end subroutine IRECV_FROM_J_0D

      logical function isPeriodic(override)
!@sum Helper function to handle optional arguments related to 
!@+ periodic boundaries
!@auth Tom Clune GSFC/SSSO/610.3
        logical, optional, intent(in) :: override

        isPeriodic = .false.
        if (present(override)) isPeriodic = override

      end function isPeriodic

      subroutine setMpiCommunicator(this, comm)
!@sum  Set MPI communicator to grid private data
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(inout) :: this
        integer, intent(in) :: comm
        this%private%MPI_COMM = comm
      end subroutine setMpiCommunicator

      integer function getMpiCommunicator(this) result(comm)
!@sum  Get MPI communicator from grid private data 
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(in) :: this
        comm = this%private%mpi_comm
      end function getMpiCommunicator

      integer function getNumProcesses(this) result(numProcesses)
!@sum  Get number of processes from grid private data 
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(in) :: this
        numProcesses = this%private%numProcesses
      end function getNumProcesses

      integer function getNumAllProcesses(this) result(numAllProcesses)
!@sum  getNumAllProcesses Get number of processes from grid private data 
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(in) :: this
        numAllProcesses = this%private%numAllProcesses
      end function getNumAllProcesses

      subroutine incrementMpiTag(this)
!@sum  increment MPI tag by one
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(inout) :: this
        integer, parameter :: MIN_TAG = 10
        integer, parameter :: MAX_TAG = 128

        integer :: tag
        tag = this%private%mpi_tag
        this%private%mpi_tag = max(mod(tag,MAX_TAG),MIN_TAG) + 1
      end subroutine incrementMpiTag

      integer function getMpiTag(this) result(mpiTag)
!@sum  gets MPI tag from distributed grid private data
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(in) :: this
        mpiTag = this%private%mpi_tag
      end function getMpiTag

      logical function hasSouthPole(this)
!@sum  returns true if GRID has a south pole
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(in) :: this
        hasSouthPole = this%private%hasSouthPole
      end function hasSouthPole

      logical function hasNorthPole(this)
!@sum  returns true if GRID has a north pole
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(in) :: this
        hasNorthPole = this%private%hasNorthPole
      end function hasNorthPole

      logical function hasPeriodicBC(this)
!@sum  returns true if GRID has periodic boundary conditions
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(in) :: this
        hasPeriodicBC = this%private%periodicBC
      end function hasPeriodicBC

      integer function getLogUnit()
!@sum  returns file unit for debugging log files
!@auth Tom Clune GSFC/SSSO/610.3
        use FileManager, only: openUnit
        integer :: unit
        character(len=40) :: logFileName

        integer, parameter :: UNINITIALIZED = -1
        integer, save :: logUnit = UNINITIALIZED

        if (logUnit == UNINITIALIZED) then
          write(logFileName,'(a,i4.4)') 'debug.', my_pet
          call openUnit(logFileName, logUnit, qbin=.false., qold=.false.)
        end if

        getLogUnit = logUnit

      end function getLogUnit

#ifdef USE_MPI

! ----------------------------------------------------------------------
   function MPIgridBounds(grid) result(AIbounds)
!@sum Returns an integer array of the grid's axis indices 
! ----------------------------------------------------------------------
     type (DIST_Grid), intent(in) :: grid
     integer :: AIbounds(4)
     type(AxisIndex), dimension(:,:), pointer :: AI

     allocate(AI(0:npes_world-1,3))
     call getMPIAxisIndex(grid, AI)

     AIbounds(1) = AI(my_pet,1)%min
     AIbounds(2) = AI(my_pet,1)%max
     AIbounds(3) = AI(my_pet,2)%min
     AIbounds(4) = AI(my_pet,2)%max

     deallocate(AI)

   end function MPIgridBounds

! ----------------------------------------------------------------------
   subroutine getMPIAxisIndex(grid, AI)
!@sum Computes and initializes an integer array of the grid's axis indices
!@+ corresponding to the domain decompositions bounds
! ----------------------------------------------------------------------
     type (DIST_Grid), intent(in) :: grid
     type (axisIndex) :: AI(0:,:)

     character(len=maxStrLen) :: IAm='getAxisIndex'     
     integer :: p, npes_end
     integer, allocatable   :: jms(:)

     allocate(jms(0:npes_world-1))

     jms(0:grid%npes_used-1) = getLatDist(grid%jm_world, grid%npes_used)
     if(grid%npes_used<npes_world) jms(grid%npes_used:npes_world-1) = 0

     AI = computeAxisIndex(grid%im_world, jms)

     npes_end = size(AI,1)
     AI(npes_world:npes_end-1,2)%min = AI(npes_world-1,2)%max + 1
     AI(npes_world:npes_end-1,2)%max = AI(npes_world-1,2)%max

     deallocate(jms)

   end subroutine getMPIAxisIndex

! ----------------------------------------------------------------------
   function computeAxisIndex(imGlob, jmsGlob) result(thisAI)
!@sum Computes an integer array of the grid's axis indices corresponding
!@+ to the domain decompositions bounds
! ----------------------------------------------------------------------
     integer, intent(in) :: imGlob, jmsGlob(0:)
     type (axisIndex) :: thisAI(0:npes_world-1,3)
     integer :: p 
     
     do p = 0, npes_world - 1
        thisAI(p,1)%min = 1
        thisAI(p,1)%max = imGlob
        if (p==0) then
           thisAI(p,2)%min = 1
           thisAI(p,2)%max = jmsGlob(p)
        else
           thisAI(p,2)%min = thisAI(p-1,2)%max + 1
           thisAI(p,2)%max = thisAI(p-1,2)%max + jmsGlob(p)
        end if
     end do
     
   end function computeAxisIndex

! ----------------------------------------------------------------------
   function getLatDist(jm, numProcesses) result(latsPerProcess)
!@sum Returns an array of distributed latitudes per process
! ----------------------------------------------------------------------
     ! Contstraint: assumes jm >=4.
     integer, intent(in) :: jm
     integer, intent(in) :: numProcesses
     integer :: latsPerProcess(0:numProcesses-1)
     
     integer :: excess, npes_used, p
     integer :: localAdjustment
     
     latsPerProcess = 0
     
     ! Set minimum requirements per processor
     ! Currently this is 1 lat/proc away from poles
     ! and 2 lat/proc at poles
     select case (npes_world)
     case (1)
        latsPerProcess = jm
        return
     case (2)
        latsPerProcess(0) = jm/2
        latsPerProcess(1) = jm - (jm/2)
        return
     case (3:)
        npes_used = min(numProcesses, jm-2)
        
        ! 1st cut - round down
        latsPerProcess(0:numProcesses-1) = JM/npes_used  ! round down
        
        ! Fix at poles
        latsPerProcess(0) = max(2, latsPerProcess(0))
        latsPerProcess(numProcesses-1) = max(2, latsPerProcess(numProcesses-1))
        
        ! redistribute excess
        excess = JM - sum(latsPerProcess(0:numProcesses-1))
        ! redistribute any remaining excess among interior processors
        do p = 1, numProcesses - 2
           localAdjustment = (p+1)*excess/(numProcesses-2) - (p*excess)/(numProcesses-2)
           latsPerProcess(p) = latsPerProcess(p) + localAdjustment
        end do
     end select
     
   end function getLatDist

#endif

    end module dist_grid_mod

