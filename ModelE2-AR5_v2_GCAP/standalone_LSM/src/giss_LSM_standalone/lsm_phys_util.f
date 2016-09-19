!@sum lsm_phys_util.f  Generic time and physics utility subroutines and 
!@sum functions that can be used by different versions of drivers for
!@sum lsm_standalone.
!@auth N.Y.Kiang
!----------------------------------------------------------------------

      module lsm_phys_util

      implicit none

      private

      public UNDEFINT, CDN_max_grid
      public mon,LEAPYR, month, fdiffuse, get_neutral_drag

      character*3, dimension(12),parameter :: mon =
     &  ( / 'JAN','FEB','MAR','APR','MAY','JUN'
     &     ,'JUL','AUG','SEP','OCT','NOV','DEC' /)
      integer, parameter :: UNDEFINT = -2**30
      real*4,save,dimension(:,:),allocatable :: CDN_max_grid !neutral drag coeff.topo

      contains
!======================================================================!


      function LEAPYR(year) 
      logical :: LEAPYR
      integer :: year
      LEAPYR = .false.
      if(mod(year,400) == 0) then
         LEAPYR = .true.
      else
         if(mod(year,4) == 0 .and. mod(year,100) /= 0) then
           LEAPYR = .true.
         endif
      endif
      end function LEAPYR

!----------------------------------------------------------------------

      
!----------------------------------------------------------------------

      function month(yr,yr_sec)
      implicit none
      integer :: month
      integer,intent(in)  :: yr        ! year A.D.
      integer,intent(in)  :: yr_sec    ! seconds from beginning of year
      !-------
      integer, parameter :: daysec = 86400 !24*60*60
      integer, dimension(12) :: monthdays

      integer :: n, countsec

      if (LEAPYR(yr)) then 
         monthdays =
     &     (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
      else
         monthdays =
     &     (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
      endif

      countsec = 0
      n = 0

      do while ((countsec.le.yr_sec).and.(n.lt.12) )
         n = n+1
         countsec=countsec + monthdays(n)*daysec
      end do

      month = n
      !print *,"month",month
      end function month

!----------------------------------------------------------------------

      function fdiffuse(cos_solzen,day,Rg) Result(fdif)
      use CONSTANT
      implicit none
      real*8,intent(in) :: cos_solzen !Cosine solar zenith angle (radians)
      real*8,intent(in) :: day        !Day of year
      real*8,intent(in) :: Rg         !Global incident radiation (W m-2)
      real*8 :: fdif
      !---Local------
      real*8,parameter :: SOLARCONST = 1370.0d0  !Solar constant, (W m-2)
      real*8 :: sbeta !Sine of solar elevation angle
      real*8 :: S0    !Incident rad at top-of-atmosphere at time (W m-2)
      real*8 :: transm !Transmissivity through the atmosphere (fraction)
      real*8 :: R, K   

      sbeta = cos_solzen !cos(zen angle)=sin(elev angle)
      
      if(sbeta.gt.0.0d0)then
        S0=SOLARCONST*(1.0d0+0.033d0*cos(2.0d0*pi*day/365.0d0))*sbeta
        transm=Rg/S0 !Rg=I0/2.3
        R=0.847d0-1.61d0*sbeta+1.04d0*(sbeta**2.d0)
        K=(1.47d0-R)/1.66d0
        if(transm.le.0.22d0)then
          fdif=1.0d0
        elseif(transm.le.0.35)then
          fdif=1.0d0-6.4d0*((transm-0.22d0)**2.d0)
        elseif(transm.le.K)then
          fdif=1.47d0-1.66d0*transm
        else
          fdif=R
        endif
      else
        fdif = 1.d0
      endif
      end function fdiffuse

  !-----------------------------------------------------------------------

      subroutine get_neutral_drag(i0f,i1f,j0f,j1f)
!     Opens the drag coefficent file max(topography,vegetation) from GCM
      use FILEMANAGER, only : openunit, closeunit
      !use PARAM, only : sync_param
      integer,intent(in) :: i0f,i1f,j0f,j1f
      !real*4, dimension(:,:), allocatable :: CDN_max_grid
      !------
      character*80 :: title
      character*80 :: neutral_drag_file
      integer :: iu_nd

      allocate(CDN_max_grid(i0f:i1f,j0f:j1f))
      !call sync_param("CD_coef",neutral_drag_file)
      call openunit("CD_coef",iu_nd,.true.,.true.)
      read(iu_nd) title, CDN_max_grid(i0f:i1f,j0f:j1f)
      !write(0,*) "Set neutral drag file ", iu_nd  !writes to screen
      print *, "Set neutral drag file ", iu_nd  !writes to log file
      
      call closeunit(iu_nd)

      end subroutine get_neutral_drag
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


      end module lsm_phys_util
