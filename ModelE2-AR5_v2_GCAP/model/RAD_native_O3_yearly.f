#include "rundeck_opts.h"

      module RAD_native_O3
!@sum  Declaring a few new variables necessary for reading ozone
!@+ data for radition code at native GCM horizonal resolution.
!@+ This version expects one file per year.
!@auth Greg Faluvegi
!@ver  1.0

      implicit none

      SAVE

      real*8, allocatable, dimension(:,:,:,:) :: O3YEAR
      real*8, allocatable, dimension(:,:,:,:) :: O3ICMA,O3JCMA
      real*8, allocatable, dimension(:,:,:)::O3JDAY_native,O3JREF_native

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

      allocate(           O3YEAR(I_0H:I_1H,J_0H:J_1H,NLO3,0:12) )
      allocate(           O3ICMA(I_0H:I_1H,J_0H:J_1H,NLO3,12) )
      allocate(           O3JCMA(I_0H:I_1H,J_0H:J_1H,NLO3,12) )
      allocate(    O3JDAY_native(NLO3,I_0H:I_1H,J_0H:J_1H) )
      allocate(    O3JREF_native(NLO3,I_0H:I_1H,J_0H:J_1H) )
      
      return
      end subroutine alloc_RAD_native_O3


      subroutine native_UPDO3D(JYEARO,JJDAYO)
!@sum A simplified version of UPDO3D to work at GCM resolution,
!@+ and with no DSO3 in the stratosphere or CH4 dependence in
!@+ the troposphere [other than that inherent in the offline files].
!@+ This version expects one file per year.
!@auth G. Faluvegi (directly from A. Lacis/R. Ruedy)

      use filemanager, only : openunit,closeunit,nameunit
      use domain_decomp_atm, only: grid, get, am_i_root, readt_parallel
      use param
      use RADPAR, only: plb0,plbO3,NLO3
      use RAD_native_O3
      implicit none

! Routine expects offline ozone files on 49 pressure levels
! and native GCM horizonal resolution (IM x JM) and no trend file.

      character*80 :: title
      logical :: qexist
      character*40 :: DDFILE
      integer :: ifile,I,J,L,M,N,JYEARX,JYEARO,JJDAYO,MI,MJ
      real*8 :: WTI,WTJ,WTMJ,WTMI,XMI

      integer, SAVE :: IYR=0, JYRNOW=0, IYRDEC=0, IFIRST=1, JYR
      save ddfile,ifile

      INTEGER :: J_0, J_1, I_0, I_1

      CALL GET(grid, J_STRT    =J_0,  J_STOP    =J_1,
     &               I_STRT    =I_0,  I_STOP    =I_1)

      if(jyearo < 0)then
        jyearo=abs(jyearo)
        print *, 'warning: not allowing ozone cycling.' 
      end if

      if(IFIRST==1) then  ! --------------------- first time through ---
        if(plbo3(1) < plb0(1)) plbo3(1)=plb0(1)                  ! ??
        ! cyclical case no longer allowed so ozone not read here.
!!!!!!!!#ifdef TRACERS_SPECIAL_Shindell
        ! read the 3D field for O3 RCOMPX reference calls
        call openunit ('Ox_ref',ifile,.true.,.true.)
        do L=1,NLO3
          call readt_parallel
     &    (grid,ifile,nameunit(ifile),O3JREF_native(L,:,:),1)
        enddo
        call closeunit(ifile)
!!!!!!!!!#endif
        IFIRST=0
      endif           ! ---------------end of first time through ---


C     To time input data READs, JYEARX is set ahead of JYEARO by 15 days
C     ------------------------------------------------------------------
      JYEARX=JYEARO+(JJDAYO+15)/366

! Next section skipped if O3YEAR already ready for O3JDAY-defining
      IF(JYEARX.ne.JYRNOW)THEN

C****
C**** Get 13 months of O3 data O3YEAR starting with the leading December
C****
      IYR=JYEARX
      JYR=JYEARX+1 ! note 1 years steps assumed

C**** Get first decadal file
      write(ddfile,'(a7,i4.4)') 'O3file_',abs(IYR)
      inquire(file=trim(ddfile),exist=qexist)
      if(.not.qexist)then
        print *,'File '//trim(ddfile)//' nonexistant.'
        call stop_model('nonexistant ozone file',13)
      end if
      call openunit(ddfile,ifile,.true.,.true.)               ! IYR
      do M=1,12 ; do L=1,NLO3
        call readt_parallel
     &  (grid,ifile,nameunit(ifile),O3ICMA(:,:,L,M),1)
      enddo ; enddo
      call closeunit(ifile)

C**** Define December 
!     Read and use prior decadal file to define prior year December
      if(IYRDEC.ne.JYEARX-1) then
        write(ddfile,'(a7,i4.4)') 'O3file_',abs(IYR-1)
        inquire(file=trim(ddfile),exist=qexist)
        if(.not.qexist)then
          print *,'File '//trim(ddfile)//' nonexistant.'
          call stop_model('nonexistant ozone file',13)
        end if
        call openunit(ddfile,ifile,.true.,.true.)          ! KYR
        do M=1,12 ; do L=1,NLO3
          call readt_parallel
     &    (grid,ifile,nameunit(ifile),O3JCMA(:,:,L,M),1)
        enddo ; enddo
        call closeunit(ifile)
        do J=J_0,J_1 ; do L=1,NLO3 ; do I=I_0,I_1     
          O3YEAR(I,J,L,0)=O3JCMA(I,J,L,12)
        enddo; enddo; enddo
        IYRDEC=JYEARX   ! Set flag to indicate December data is current
      endif 

! I do not think this section is needed:
C**** Get next decadal file
!     write(ddfile,'(a7,i4.4)') 'O3file_',abs(JYR)
!     inquire(file=trim(ddfile),exist=qexist)
!     if(.not.qexist)then
!       print *,'File '//trim(ddfile)//' nonexistant.'
!       call stop_model('nonexistant ozone file',13)
!     end if
!     call openunit(ddfile,ifile,.true.,.true.)               ! JYR
!     do M=1,12 ; do L=1,NLO3
!       call readt_parallel
!    &  (grid,ifile,nameunit(ifile),O3JCMA(:,:,L,M),1)
!     enddo ; enddo
!     call closeunit (ifile)

      if(JYEARX.ne.IYRDEC) then ! we are not done with prior December

        if(JYEARX==IYRDEC+1) then  ! copy data from M=12 -> M=0
          O3YEAR(:,:,:,0)=O3YEAR(:,:,:,12) ! DEC from prior year
          IYRDEC=JYEARX  ! Set flag to indicate December data is current
        else 
          call stop_model('GF: not expecting this case.',13)
! was:    Interpolate prior December from the decadal files - start-up
        endif 
      
      end if ! now we are done with prior December


C            Fill in a full year of O3 data 
C            -----------------------------------------------
      do M=1,12 ; do J=J_0,J_1 ; do L=1,NLO3 ; do I=I_0,I_1 
        O3YEAR(I,J,L,M)=O3ICMA(I,J,L,M)
      enddo ; enddo ; enddo ; enddo

      if(JYEARX.ne.IYRDEC)call stop_model('yearly O3 error 2',13)
       ! original code said above if was for cyclical startup case...

      JYRNOW=JYEARX

      ENDIF ! O3YEAR defining was necessary 

C****
C**** O3JDAY_native is interp daily from O3YEAR seasonal data via JJDAYO
C****
C     the formula below yields M near the middle of month M
      XMI=(JJDAYO+JJDAYO+31-(JJDAYO+15)/61+(JJDAYO+14)/61)/61.D0
      MI=XMI
      WTMJ=XMI-MI       !   Intra-year interpolation is linear in JJDAYO
      WTMI=1.D0-WTMJ
      IF(MI > 11) MI=0
      MJ=MI+1
      do J=J_0,J_1 ; do  I=I_0,I_1 
        O3JDAY_native(:,I,J)=WTMI*O3YEAR(I,J,:,MI)+WTMJ*O3YEAR(I,J,:,MJ)
      enddo ; enddo
      return
      end subroutine native_UPDO3D


