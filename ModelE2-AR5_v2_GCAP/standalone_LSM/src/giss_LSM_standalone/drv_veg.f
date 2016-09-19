      module drv_veg
!@sum lsm_force_veg  Module to open, read, close vegetation initialization,
!     and forcing files and interpolate input data for the 
!     GISS GCM land surface module.
!     Need to call deallocate_forcings_veg at end of run.
!     vt0 and vt2=vt0+vdt are time points for data input rows for interpolation.
!     vt1=tyr model time, and vt0<=vt1<=vt2 (vt0<vt2 strictly)
!     Interpolation is linear
!     Created by N.Y. Kiang 3/2010.

      use ent_mod, only : N_PFT, N_CASA_LAYERS
      use lsm_phys_util, only : UNDEFINT

      implicit none

      private

      public allocate_iu_veg, rewind_files_veg
      public open_forcings_veg, close_forcings_veg, deallocate_drv_veg
      public init_interp_forcings_veg, get_interp_forcings_veg

      save

      !* DATASET-SPECIFIC CONSTANTS
      integer, parameter :: VEGDTSEC = 1800 !seconds, veg forcing time step
      integer :: VROWMAX  !Gets set to input parameter at initialization.
      integer :: vrow     !Internal global, row number in veg forcing file
#ifndef MIXED_CANOPY
      integer, allocatable, dimension(:) :: iu_LAI
      integer, allocatable, dimension(:) :: iu_height
#else
      integer :: iu_LAIMIXED
      integer :: iu_heightmixed
#endif
      integer :: vt0, vt2 !Veg input data point times for interpolation.
      real*8, dimension(:,:,:),pointer :: LAI0, height0
      real*8, dimension(:,:,:),pointer :: LAI2, height2


      contains
!------------------------------------------------------------------------
      function get_fname(prefix,n) 
      !Makes character file name with prefix and numbered with integer n.
      !E.g. LAI1, LAI2...LAI16
      character*80 :: get_fname
      character*20,intent(in) :: prefix
      integer, intent(in) :: n
      !-----Local-----
      character*20 :: dstr, nstr
      integer :: digits

         digits = floor(log10(dble(n))) + 1
         write (dstr, '(i1)') digits
         write (nstr,'(i'//trim(dstr)//')') n
         get_fname = trim(prefix)//trim(nstr)
         print *,get_fname
      end function get_fname

!------------------------------------------------------------------------

      subroutine allocate_iu_veg (N_PFT)
      implicit none
      integer, intent(in) :: N_PFT
      !Global ent_mod, only : N_PFT
      !Global integer :: iu_LAI(N_PFT), iu_height(N_PFT)

      allocate(iu_LAI(N_PFT))
      allocate(iu_height(N_PFT))
      
      end subroutine allocate_iu_veg 
!------------------------------------------------------------------------
      subroutine deallocate_iu_veg
      deallocate(iu_LAI)
      deallocate(iu_height)
      end subroutine deallocate_iu_veg
!------------------------------------------------------------------------
      subroutine open_forcings_veg(VROWMAX_IN) 
!@sum open_forcings_veg.  Open LAI and height forcings files if force_VEG.
!vrow is initialized; skip optional
      use filemanager
      implicit none
      !Global ent_mod, only : N_PFT, N_CASA_LAYERS
      !integer,intent(in) :: skip !Number of records to skip at beginning of forcing file.
      !integer,intent(inout) :: vrow !Make global internal to module
      integer,intent(in) :: VROWMAX_IN !Max number of rows in file, read_input.
      !----Local--------
      integer :: p
      character*20 :: pre, fname

#ifdef LSM_DRV_SINGLE
      call openunit("LAI", iu_LAI,.true.,.true.)
      write(*,*) "Opened LAI single."
      call openunit("HEIGHT",iu_height,.true.,.true.)
      write(*,*) "Opened HEIGHT single."
#else   
#ifndef MIXED_CANOPY
      pre = "LAI"
      do p=1,N_PFT
         fname = get_fname(pre,p)
         call openunit(trim(fname),iu_LAI(p),.true.,.true.)
         write(*,*) fname
      enddo

      pre = "HEIGHT"
      do p=1,N_PFT
         fname = get_fname(pre,p)
         call openunit(trim(fname),iu_height(p),.true.,.true.)
         write(*,*) fname
      enddo
#else !MIXED_CANOPY
      call openunit("LAIMIXED",iu_LAIMIXED,.true.,.true.)
      call openunit("HEIGHTMIXED",iu_heightmixed,.true.,.true.)
#endif !MIXED_CANOPY
#endif


      ! skip some records if needed
      !call skip_forcings_veg( skip, vrow )
      vrow = 1 !At first row in forcing file.
      VROWMAX = VROWMAX_IN

      end subroutine open_forcings_veg

!------------------------------------------------------------------------

      subroutine close_forcings_veg()
      use filemanager
      implicit none
      !Global ent_mod, only : N_PFT, N_CASA_LAYERS
      !Global : iu_LAI, iu_height, iu_LAImixed, iu_heightmixed
      !---- Local -----
      integer :: p

#ifndef MIXED_CANOPY
      do p=1,N_PFT
         call closeunit(iu_LAI(p))
      enddo

      do p=1,N_PFT
         call closeunit(iu_height(p))
      enddo
#else !MIXED_CANOPY
      call closeunit(iu_LAIMIXED)
      call closeunit(iu_heightmixed)
#endif !MIXED_CANOPY

      end subroutine close_forcings_veg

!------------------------------------------------------------------
        subroutine rewind_files_veg(vdt,vtstart)
      integer,intent(in) :: vdt, vtstart
      !integer, intent(inout) :: vrow !Global to module
      !----Local----
      integer :: p

#ifndef MIXED_CANOPY
      do p=1,N_PFT
         rewind(iu_LAI(p))
         rewind(iu_height(p))
      enddo
#else
      rewind(iu_LAImixed)
      rewind(iu_heightmixed)
#endif      
      vrow=1
      vt0=vtstart
      vt2=UNDEFINT
      end subroutine rewind_files_veg

!------------------------------------------------------------------------
      subroutine goto_nextfile_veg()
      !integer, intent(inout) :: vrow !Global to module

      call close_forcings_veg()
      call deallocate_drv_veg(.true.)
      !## No naming convention yet for 
      call stop_model("Multiple split-time veg files must be 
     &customized to dataset",255)
      vrow=0
      end subroutine goto_nextfile_veg
!------------------------------------------------------------------------
      subroutine read_forcings_veg(I0,I1,J0,J1,I0f,I1f,J0f,J1f
     &     ,LAI,height,vdt)
      use filemanager
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      integer,intent(in) :: I0f,I1f,J0f,J1f
      real*8, dimension(:,:,:),pointer :: LAI, height
      integer,intent(in) :: vdt
      !integer,intent(inout) :: vrow !Global to module
      !Global ent_mod, only :  N_PFT
      !---- Local -------
      integer :: p, niter, niter1
      real*4, dimension(I0f:I1f,J0f:J1f) :: buf, buf1
      real*4, dimension(N_PFT,I0f:I1f,J0f:J1f) :: vegbuf

      vrow = vrow+1
      if (vrow.gt.VROWMAX) then
         call goto_nextfile_veg()
      endif
            
      !* LAI
      if ( associated( LAI )) then
#ifndef MIXED_CANOPY
         !Get niter check
         read(iu_LAI(1)) niter, buf
         niter1 = niter
         LAI(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
         do p=1,N_PFT        
            read(iu_LAI(p)) niter, buf
            if(niter.ne.niter1) 
     &           call stop_model("read_forcings:bad rec",255)
            LAI(p,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
         enddo
#else !MIXED_CANOPY
        read(iu_LAIMIXED) niter, vegbuf
        LAI(:,I0:I1,J0:J1) = vegbuf(:,I0:I1,J0:J1)
        if(niter.ne.niter1) 
     &       call stop_model("read_forcings:bad rec",255)
#endif !MIXED_CANOPY
      endif

      !* Height
      if ( associated( height )) then
#ifndef MIXED_CANOPY
         do p=1,N_PFT
            read(iu_height(p)) niter, buf
            if(niter.ne.niter1) 
     &           call stop_model("read_forcings:bad rec",255)
            height(p,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
         enddo
#else  !MIXED_CANOPY
         read(iu_heightmixed) niter, vegbuf
         height(:,I0:I1,J0:J1) = vegbuf(:,I0:I1,J0:J1)
         if (niter.ne.niter1)
     &        call stop_model("read_forcings:bad rec",255)
#endif !MIXED_CANOPY
      endif

      end subroutine read_forcings_veg

!------------------------------------------------------------------------
      subroutine skip_forcings_veg( nskip )
      integer, intent(in) :: nskip
      !integer, intent(inout) :: vrow !Global to module
      !-----Local-----
      integer :: i,k

      do k=1,nskip
#ifndef MIXED_CANOPY
         do i=1,N_PFT
            read(iu_LAI(i))
            read(iu_height(i))
            vrow = vrow+1
         enddo
#else !MIXED_CANOPY
!      read(iu_LAIMIXED)
!      read(iu_heightmixed)
         write(*,*) "Code to read global mixed files not done yet" !##
#endif !MIXED_CANOPY 

      enddo
      end subroutine skip_forcings_veg

!------------------------------------------------------------------------

      subroutine interp_veg(I0,I1,J0,J1,vt1,LAI,height)
!@sum interp_veg.  Simple linear interpolation of vegetation inputs
!@+   at time vt1 between vt0 and vt2.
      integer,intent(in) :: I0,I1,J0,J1
      integer,intent(in) :: vt1 !seconds. Model time.
      real*8, dimension(:,:,:),pointer :: LAI, height
      !----Local-------
      integer,dimension(N_PFT,I0:I1,J0:J1) :: m

      if ((vt1.lt.vt0).or.(vt1.gt.vt2).or.(vt1.gt.vt2)) then
         call stop_model("Bad bounds for veg interpolation.",255)
      endif
      if (vt0.eq.vt2) then  
         call stop_model("Unnecessary interpolation: fix code",255)
         write(*,*) "vt0=",vt0,", vt1=",vt1,", vt2=",vt2
      endif

      m = (LAI2(:,:,:)-LAI0(:,:,:))/(vt2-vt0)
      LAI(:,:,:) = m * (vt1-vt0) + LAI0(:,:,:) 

      m = (height2(:,:,:)-height0(:,:,:))/(vt2-vt0)
      height(:,:,:) = m * (vt1-vt0) + height0(:,:,:)

      end subroutine interp_veg
!------------------------------------------------------------------------
      subroutine init_interp_forcings_veg(
     &     tyr,updt,vdt, vtstart, vrowmax_in, N_PFT
     &     ,I0,I1,J0,J1,I0f,I1f,J0f,J1f
     o     ,LAI,height)
      implicit none
      integer,intent(in) :: tyr !(sec) Actual time in run. Multiples of updt.
      integer,intent(in) :: updt !(sec) Veg update interval for model.
      integer,intent(in) :: vdt !(sec) Veg input data interval.
      integer,intent(in) :: vtstart !Seconds. Actual time for first veg input.
      !integer,intent(inout) :: vrow !Row in veg input file. Module global.
      integer,intent(in) :: vrowmax_in
      integer,intent(in) :: N_PFT 
      integer,intent(in) :: I0,I1,J0,J1,I0f,I1f,J0f,J1f
      real*8, dimension(:,:,:),pointer :: LAI, height
      !----Local-----
      integer :: vt0, vt1, vt2 !(sec) vt0=data time, vt1=tyr,vt2=vt0+vdt
      integer :: skip

      !* Allocate data points for interpolation. Global to module.
      allocate( LAI0(N_PFT,I0:I1,J0:J1) )
      allocate( height0(N_PFT,I0:I1,J0:J1) )
      allocate( LAI2(N_PFT,I0:I1,J0:J1) )
      allocate( height2(N_PFT,I0:I1,J0:J1) )

      !* Initialize and zero.
      VROWMAX = vrowmax_in
      vrow = 1
      vt0 = vtstart
      vt1 = tyr
      vt2 = vt0 !Init so defined.
      LAI0(:,:,:) = 0.d0
      height0(:,:,:) = 0.d0
      LAI2(:,:,:) = 0.d0
      height2(:,:,:) = 0.d0

      !* Skip to correct starting row in input data, given model time.
      if (vtstart.lt.vt1) then
         !Start of veg data is before time in run.
         !Skip any rows until just before tyr
         skip = mod(vt1-vt0,vdt) + floor(real (vt1-vt0)/vdt)
         vt0 = vt0 + skip*vdt
         call skip_forcings_veg(skip)
      endif

      !* Read first data point (LAI0, height0)
      call read_forcings_veg(I0,I1,J0,J1,I0f,I1f,J0f,J1f
     &        ,LAI0,height,vdt)         

      !* Determine if interpolation is needed given vdt and updt.
      if ((vdt.le.updt).and.(mod(vdt,updt).eq.0)) then 
         !Easy case: Veg data is divisible into update time step.
         !Never need to interpolate.Could just deallocate LAI2 here.
         LAI(:,:,:) = LAI0(:,:,:) !Initialization
         vt0 = UNDEFINT
         vt2 = UNDEFINT
         nullify(LAI0)
         nullify(LAI)
      else
         !Need to get LAI2 from next +vdt input.
         vt2 = vt0 + vdt
         LAI(:,:,:) = LAI0(:,:,:) !Initialization
         LAI2(:,:,:) = LAI0(:,:,:) !Initialization
         call get_interp_forcings_veg(tyr,updt,vdt
     &        ,I0,I1,J0,J1,I0f,I1f,J0f,J1f
     &        ,LAI,height)
      endif
      end subroutine init_interp_forcings_veg
!------------------------------------------------------------------------

      subroutine get_interp_forcings_veg(tyr,updt,vdt
     &     ,I0,I1,J0,J1,I0f,I1f,J0f,J1f
     &     ,LAI,height)
      !
      implicit none
      integer,intent(in) :: tyr !(sec) Actual time in run. Multiples of updt.
      integer,intent(in) :: updt !(sec) Veg update interval for model.
      integer,intent(in) :: vdt !(sec) Veg input data interval.
      !integer,intent(inout) :: vrow !Row number in veg input file. Global.
      integer,intent(in) :: I0,I1,J0,J1,I0f,I1f,J0f,J1f
!     logical,intent(in) :: init
      real*8, dimension(:,:,:),pointer :: LAI, height
      !----Local-------
      !Global real*8, dimension(:,:,:),pointer :: LAI0, height
      !Global real*8, dimension(:,:,:),pointer :: LAI2, height
      integer:: vt1  !(sec) Gets model time.
      integer,parameter :: daysec= 86400 !60*60*24 seconds per day
      integer :: n,skip

      vt1 = tyr  !tyr = tyr_prev + updt

      if ((vdt.le.updt)) then
         !Interpolation not needed. Assume mod(updt,vdt)=0, divisible.

         !Skip data at smaller time steps than updt,
         skip = (vt1-vt0)/vdt-1 !(vt1-vt0) should = updt
         if (skip.lt.0) 
     &        call stop_model("Oops, didn't update tyr.",255)
         call skip_forcings_veg(skip)
         call read_forcings_veg(I0,I1,J0,J1,I0f,I1f,J0f,J1f
     &        ,LAI,height,vdt)
         vt0=vt1  !Save for next update time.
         vt2=vt1  !Unneeded, really UNDEFINT but for the sake of consistency.

      elseif (vdt.gt.updt) then 
         !Interpolation necessary
         if (vt1.le.vt2) then 
            !Keep current input points. Interpolate at vt1.
            call interp_veg(I0,I1,J0,J1,vt1,LAI,height)
         else           
            !Get next input point. Shift vti, LAIi, heighti.
            vt0=vt2
            vt2=vt2+vdt
            LAI0(:,:,:) = LAI2(:,:,:)
            Height0(:,:,:) = Height2(:,:,:)
            !Read next input data at vt2
            call read_forcings_veg(I0,I1,J0,J1,I0f,I1f,J0f,J1f
     &           ,LAI2,height2,vdt)
            call interp_veg(I0,I1,J0,J1,vt1,LAI,height)
         endif
      endif

      end subroutine get_interp_forcings_veg

!------------------------------------------------------------------------
      subroutine deallocate_drv_veg(force_VEG)
      logical, intent(in) :: force_VEG
      if (force_VEG) then
         deallocate( LAI0 )
         deallocate( height0 )
         deallocate( LAI2 )
         deallocate( height2 )
      endif
      end subroutine deallocate_drv_veg
!------------------------------------------------------------------------

      end module drv_veg
