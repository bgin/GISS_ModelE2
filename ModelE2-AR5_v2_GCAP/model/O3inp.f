! for the moment, this file is included in other files
!#include "rundeck_opts.h"

      module RAD_native_O3
!@sum  Declaring a few new variables necessary for reading ozone
!@+ data for radition code at native GCM horizonal resolution.
!@auth Greg Faluvegi
!@ver  1.0
      use timestream_mod, only : timestream
      implicit none

      SAVE

c      real*8, allocatable, dimension(:,:,:,:) :: delta_O3_max_min
      real*8, allocatable, dimension(:,:,:)::
     &     O3JDAY_native
c#ifdef TRACERS_SPECIAL_Shindell
     &    ,O3JREF_native
c#endif

!@var O3stream interface for reading and time-interpolating O3 files
!@+   See usage notes in timestream_mod
      type(timestream) :: O3stream

      end module RAD_native_O3


      subroutine alloc_RAD_native_O3(grid)
!@SUM  alllocate RAD_native_O3 arrays for current grid
!@auth Greg Faluvegi
!@ver  1.0
      use domain_decomp_atm, only : dist_grid, get
      use RADPAR, only : NLO3
      use RAD_native_O3
      implicit none

      type (dist_grid), intent(in) :: grid
      integer :: J_1H, J_0H, I_0H, I_1H

      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &                 I_STRT_HALO=I_0H, I_STOP_HALO=I_1H )

c      allocate( delta_O3_max_min(I_0H:I_1H,J_0H:J_1H,NLO3,0:12) )
      allocate(    O3JDAY_native(NLO3,I_0H:I_1H,J_0H:J_1H) )
c#ifdef TRACERS_SPECIAL_Shindell
      allocate(    O3JREF_native(NLO3,I_0H:I_1H,J_0H:J_1H) )
c#endif
      
      return
      end subroutine alloc_RAD_native_O3

      subroutine native_UPDO3D(JYEARO,JJDAYO)
!@sum A simplified version of UPDO3D to work at GCM resolution,
!@+ more input files, and with a much simpler time interpolation:
!@+ that is, no DSO3 in the stratosphere of CH4 dependence in
!@+ the troposphere [other than that inherent in the offline files].
!@auth M. Kelley, G. Faluvegi (directly from A. Lacis/R. Ruedy)

      use domain_decomp_atm, only: grid, get
      use param
      use RADPAR, only: plb0,plbO3,S00WM2,RATLS0,NLO3
      use RAD_native_O3
      use timestream_mod, only : init_stream,read_stream
      use pario, only : par_open,par_close,read_dist_data
      implicit none

      integer, intent(in) :: JYEARO,JJDAYO

! Routine expects offline ozone files on nlo3 pressure levels
! and native GCM horizonal resolution and no trend file.

c      real*4,dimension(NLO3):: delta_O3_now
!@dbparam use_sol_Ox_cycle if =1, a cycle of ozone is appled to
!@+ o3year, as a function of the solar constant cycle.
!@var add_sol is [S00WM2(now)-1/2(S00WM2min+S00WM2max)]/
!@+ [S00WM2max-S00WM2min] so that O3(altered) = O3(default) +
!@+ add_sol*delta_O3_max_min
      integer, save :: use_sol_Ox_cycle = 0
c      real*8 :: add_sol

      integer :: i,j,l,jyearx,fid
      logical, save :: init = .false.
      logical :: cyclic
      real*8, allocatable :: o3arr(:,:,:)

      INTEGER :: J_0, J_1, I_0, I_1

      CALL GET(grid, J_STRT    =J_0,  J_STOP    =J_1,
     &               I_STRT    =I_0,  I_STOP    =I_1)

      allocate(o3arr(grid%i_strt_halo:grid%i_stop_halo,
     &               grid%j_strt_halo:grid%j_stop_halo,nlo3))

      jyearx = abs(jyearo)

      if (.not. init) then
        init = .true.

        call sync_param("use_sol_Ox_cycle",use_sol_Ox_cycle)
        if(use_sol_Ox_cycle==1)then
          call stop_model('Greg: address delta_O3 distributed read',255)
        endif

! what is this doing in the original updo3d ??
        if(plbo3(1) < plb0(1)) plbo3(1)=plb0(1)

        cyclic = jyearo < 0

        call init_stream(grid,O3stream,'O3file','O3',
     &       0d0,1d30,'linm2m',jyearx,jjdayo,cyclic=cyclic)

#ifdef TRACERS_SPECIAL_Shindell
! read the 3D field for O3 RCOMPX reference calls
! todo: if the reference O3 corresponds to some year or month
! of the O3 timeseries, retrieve it through the read_stream
! or get_by_index mechanism instead.
        fid = par_open(grid,'Ox_ref','read')
        call read_dist_data(grid,fid,'O3',o3arr)
        call par_close(grid,fid)
        do j=j_0,j_1
        do i=i_0,i_1 
          O3JREF_native(:,I,J)=O3ARR(I,J,:)
        enddo
        enddo
#endif
        
      endif  ! end init

      call read_stream(grid,O3stream,jyearx,jjdayo,o3arr)
      do j=j_0,j_1
      do i=i_0,i_1 
        O3JDAY_native(:,I,J)=O3ARR(I,J,:)
      enddo
      enddo

      deallocate(o3arr)


c      if(use_sol_Ox_cycle==1)then
c        call stop_model('native_UPDO3D: implement sol_Ox_cycle',255)
c        add_sol = (S00WM2*RATLS0-0.5d0*(S0min+S0max))/(S0max-S0min)
c        write(6,661)JJDAYO,S00WM2*RATLS0,S0min,S0max,add_sol
c  661 format('JJDAYO,S00WM2*RATLS0,S0min,S0max,frac=',I4,3F9.2,F7.3)

c      do J=J_0,J_1 ; do  I=I_0,I_1 
c          delta_O3_now(:) = WTMI*delta_O3_max_min(I,J,:,MI) +
c     &                      WTMJ*delta_O3_max_min(I,J,:,MJ)
c          O3JDAY_native(:,I,J) = O3JDAY_native(:,I,J) + 
c     &                           add_sol*delta_O3_now(:)
c      enddo ; enddo
c     endif

      return
      end subroutine native_UPDO3D
