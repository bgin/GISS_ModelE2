      module ent_forcings
      implicit none

      private

      public open_forcings_file, close_forcings_file,  read_forcings
      public rewind_files

      integer I0_file, I1_file, J0_file, J1_file

      integer iu_SAT  ! - surface air temp
      integer iu_QLAT  ! - latent heat flux
      integer iu_QSEN  ! - sensible heat flux
      integer iu_PREC  ! - precip
      integer iu_MP1  ! - soil matric potential
      integer iu_MP2  ! - soil matric potential
      integer iu_MP3  ! - soil matric potential
      integer iu_MP4  ! - soil matric potential
      integer iu_MP5  ! - soil matric potential
      integer iu_MP6  ! - soil matric potential
      integer iu_SIC1  ! - soil ice
      integer iu_SIC2  ! - soil ice
      integer iu_SIC3  ! - soil ice
      integer iu_SIC4  ! - soil ice
      integer iu_SIC5  ! - soil ice
      integer iu_SIC6  ! - soil ice
      integer iu_VS  ! -  V wind speed
      integer iu_US  ! - U wind speed
      integer iu_CH  ! - heat transf coeff
      integer iu_PS  ! - surface pressure
      integer iu_TP  ! - canopy temperature
      integer iu_DVIS  ! - direct vis
      integer iu_VIS  ! - total vis
      integer iu_QFOL  ! - Qf - foliage surf humidity
      integer iu_ZEN  ! - cos zenith angle
! soil temperature and moisture now from all 6 GISS layers (0-3.5m) -PK 11/29/06
!***note that soil "moisture" is in m, so divide by layer depth to get volum frac*** -PK 
      integer iu_W1   ! - soil water (m), layer 1 (0-0.1 m)
      integer iu_W2   ! - soil water, layer 2 (0.1-0.27 m)
      integer iu_W3   ! - soil water, layer 3 (0.27-0.57 m)
      integer iu_W4   ! - soil water, layer 4 (0.57-1.08 m)
      integer iu_W5   ! - soil water, layer 5 (1.08-1.97 m)
      integer iu_W6   ! - soil water, layer 6 (1.97-3.50 m)
      integer iu_SOILT1   ! - soil temp (C), layer 1
      integer iu_SOILT2   ! - soil temp, layer 2
      integer iu_SOILT3   ! - soil temp, layer 3
      integer iu_SOILT4   ! - soil temp, layer 4
      integer iu_SOILT5   ! - soil temp, layer 5
      integer iu_SOILT6   ! - soil temp, layer 6
!      integer iu_SOILT30cm  !soil temp avg top 30 cm
!      integer iu_SMOIST30cm  !soil volum moist avg top 30 cm

      integer iu_LAI1
      integer iu_LAI2
      integer iu_LAI3
      integer iu_LAI4
      integer iu_LAI5
      integer iu_LAI6
      integer iu_LAI7
      integer iu_LAI8
      integer iu_LAI9
      integer iu_LAI10
      integer iu_LAI11
      integer iu_LAI12
      integer iu_LAI13
      integer iu_LAI14
      integer iu_LAI15
      integer iu_LAI16
      integer iu_LAIMIXED

      integer iu_height1
      integer iu_height2
      integer iu_height3
      integer iu_height4
      integer iu_height5
      integer iu_height6
      integer iu_height7
      integer iu_height8
      integer iu_height9
      integer iu_height10
      integer iu_height11
      integer iu_height12
      integer iu_height13
      integer iu_height14
      integer iu_height15
      integer iu_height16
      integer iu_heightmixed

      contains

      subroutine open_forcings_file( i0f, i1f, j0f, j1f, force_VEG,
     i     N_CASA_LAYERS,skip) !mixed_VEG,skip)
      use filemanager
      integer, intent(in) ::  i0f, i1f, j0f, j1f
      logical, intent(in) :: force_VEG
      integer, intent(in) :: N_CASA_LAYERS
!      logical, intent(in) :: mixed_VEG
      integer skip !#HACK, number of records to skip at beginning of forcing file.
      
      ! These indices could be provided in input file. They describe
      ! the dimensions of arrays in input files and should be compatible
      ! with M, JM, I0, I1, J0, J1 in the main program

      I0_file = i0f
      I1_file = i1f
      J0_file = j0f
      J1_file = j1f


      call openunit("CH",iu_CH,.true.,.true.)
      call openunit("DVIS",iu_DVIS,.true.,.true.)
      call openunit("MP1",iu_MP1,.true.,.true.)
      call openunit("MP2",iu_MP2,.true.,.true.)
      call openunit("MP3",iu_MP3,.true.,.true.)
      call openunit("MP4",iu_MP4,.true.,.true.)
      call openunit("MP5",iu_MP5,.true.,.true.)
      call openunit("MP6",iu_MP6,.true.,.true.)
      call openunit("PREC",iu_PREC,.true.,.true.)
      call openunit("PS",iu_PS,.true.,.true.)
      call openunit("QFOL",iu_QFOL,.true.,.true.)
      call openunit("QLAT",iu_QLAT,.true.,.true.)
      call openunit("QSEN",iu_QSEN,.true.,.true.)
      call openunit("SAT",iu_SAT,.true.,.true.)
      call openunit("SIC1",iu_SIC1,.true.,.true.)
      call openunit("SIC2",iu_SIC2,.true.,.true.)
      call openunit("SIC3",iu_SIC3,.true.,.true.)
      call openunit("SIC4",iu_SIC4,.true.,.true.)
      call openunit("SIC5",iu_SIC5,.true.,.true.)
      call openunit("SIC6",iu_SIC6,.true.,.true.)
      call openunit("TP",iu_TP,.true.,.true.)
      call openunit("US",iu_US,.true.,.true.)
      call openunit("VIS",iu_VIS,.true.,.true.)
      call openunit("VS",iu_VS,.true.,.true.)
      call openunit("ZEN",iu_ZEN,.true.,.true.)
      call openunit("W1",iu_W1,.true.,.true.)
!      if (N_CASA_LAYERS.eq.2) then
      call openunit("W2",iu_W2,.true.,.true.)
      call openunit("W3",iu_W3,.true.,.true.)
      call openunit("W4",iu_W4,.true.,.true.)
      call openunit("W5",iu_W5,.true.,.true.)
      call openunit("W6",iu_W6,.true.,.true.)
!      endif
      call openunit("SOILT1",iu_SOILT1,.true.,.true.)
!      if (N_CASA_LAYERS.eq.2) then
      call openunit("SOILT2",iu_SOILT2,.true.,.true.)
      call openunit("SOILT3",iu_SOILT3,.true.,.true.)
      call openunit("SOILT4",iu_SOILT4,.true.,.true.)
      call openunit("SOILT5",iu_SOILT5,.true.,.true.)
      call openunit("SOILT6",iu_SOILT6,.true.,.true.)
!      endif
!      call openunit("SOILT30cm",iu_SOILT30cm,.true.,.true.)
!      call openunit("SMOIST30cm",iu_SMOIST30cm,.true.,.true.)

      print *,"In ent_forcings before force_VEG, LAI..."  !##
      if ( force_VEG ) then
!         if (.not.mixed_VEG) then
#ifndef MIXED_CANOPY
            call openunit("LAI1",iu_LAI1,.true.,.true.)
            call openunit("LAI2",iu_LAI2,.true.,.true.)
            call openunit("LAI3",iu_LAI3,.true.,.true.)
            call openunit("LAI4",iu_LAI4,.true.,.true.)
            call openunit("LAI5",iu_LAI5,.true.,.true.)
            call openunit("LAI6",iu_LAI6,.true.,.true.)
            call openunit("LAI7",iu_LAI7,.true.,.true.)
            call openunit("LAI8",iu_LAI8,.true.,.true.)
#ifdef PFT_MODEL_ENT
            call openunit("LAI9",iu_LAI9,.true.,.true.)
            call openunit("LAI10",iu_LAI10,.true.,.true.)
            call openunit("LAI11",iu_LAI11,.true.,.true.)
            call openunit("LAI12",iu_LAI12,.true.,.true.)
            call openunit("LAI13",iu_LAI13,.true.,.true.)
            call openunit("LAI14",iu_LAI14,.true.,.true.)
            call openunit("LAI15",iu_LAI15,.true.,.true.)
            call openunit("LAI16",iu_LAI16,.true.,.true.)
#endif !PFT_MODEL_ENT
!         else !mixed_VEG
#else  !MIXED_CANOPY
            call openunit("LAIMIXED",iu_LAIMIXED,.true.,.true.)
!         endif
#endif !MIXED_CANOPY
      endif

      print *,"before force_VEG, HEIGHT..."
      if ( force_VEG ) then
!         if (.not.mixed_VEG) then
#ifndef MIXED_CANOPY
            call openunit("HEIGHT1",iu_height1,.true.,.true.)
            call openunit("HEIGHT2",iu_height2,.true.,.true.)
            call openunit("HEIGHT3",iu_height3,.true.,.true.)
            call openunit("HEIGHT4",iu_height4,.true.,.true.)
            call openunit("HEIGHT5",iu_height5,.true.,.true.)
            call openunit("HEIGHT6",iu_height6,.true.,.true.)
            call openunit("HEIGHT7",iu_height7,.true.,.true.)
            call openunit("HEIGHT8",iu_height8,.true.,.true.)
#ifdef PFT_MODEL_ENT
            call openunit("HEIGHT9",iu_height9,.true.,.true.)
            call openunit("HEIGHT10",iu_height10,.true.,.true.)
            call openunit("HEIGHT11",iu_height11,.true.,.true.)
            call openunit("HEIGHT12",iu_height12,.true.,.true.)
            call openunit("HEIGHT13",iu_height13,.true.,.true.)
            call openunit("HEIGHT14",iu_height14,.true.,.true.)
            call openunit("HEIGHT15",iu_height15,.true.,.true.)
            call openunit("HEIGHT16",iu_height16,.true.,.true.)
#endif !PFT_MODEL_ENT
!         else !mixed_VEG
#else  !MIXED_CANOPY
            call openunit("HEIGHTMIXED",iu_heightmixed,.true.,.true.)
!         endif
#endif !MIXED_CANOPY
      endif

      ! skip some records if needed

      call skip_forcings( skip, force_VEG,N_CASA_LAYERS  )

      end subroutine open_forcings_file


      subroutine close_forcings_file( force_VEG,N_CASA_LAYERS)!,mixed_VEG)
      use filemanager
      logical force_VEG !, mixed_VEG
      integer, intent(in) :: N_CASA_LAYERS

      call closeunit(iu_CH)
      call closeunit(iu_DVIS)
      call closeunit(iu_MP1)
      call closeunit(iu_MP2)
      call closeunit(iu_MP3)
      call closeunit(iu_MP4)
      call closeunit(iu_MP5)
      call closeunit(iu_MP6)
      call closeunit(iu_PREC)
      call closeunit(iu_PS)
      call closeunit(iu_QFOL)
      call closeunit(iu_QLAT)
      call closeunit(iu_QSEN)
      call closeunit(iu_SAT)
      call closeunit(iu_SIC1)
      call closeunit(iu_SIC2)
      call closeunit(iu_SIC3)
      call closeunit(iu_SIC4)
      call closeunit(iu_SIC5)
      call closeunit(iu_SIC6)
      call closeunit(iu_TP)
      call closeunit(iu_US)
      call closeunit(iu_VIS)
      call closeunit(iu_VS)
      call closeunit(iu_ZEN)
      call closeunit(iu_W1)
!      if (N_CASA_LAYERS.eq.2) then
      call closeunit(iu_W2)
      call closeunit(iu_W3)
      call closeunit(iu_W4)
      call closeunit(iu_W5)
      call closeunit(iu_W6)
!      endif
      call closeunit(iu_SOILT1)
!      if (N_CASA_LAYERS.eq.2) then
      call closeunit(iu_SOILT2)
      call closeunit(iu_SOILT3)
      call closeunit(iu_SOILT4)
      call closeunit(iu_SOILT5)
      call closeunit(iu_SOILT6)
!      endif
!      call closeunit(iu_SOILT30cm)
!      call closeunit(iu_SMOIST30cm)

       if ( force_VEG ) then
!          if (.not.mixed_VEG) then
#ifndef MIXED_CANOPY
             call closeunit(iu_LAI1)
             call closeunit(iu_LAI2)
             call closeunit(iu_LAI3)
             call closeunit(iu_LAI4)
             call closeunit(iu_LAI5)
             call closeunit(iu_LAI6)
             call closeunit(iu_LAI7)
             call closeunit(iu_LAI8)
#ifdef PFT_MODEL_ENT
             call closeunit(iu_LAI9)
             call closeunit(iu_LAI10)
             call closeunit(iu_LAI11)
             call closeunit(iu_LAI12)
             call closeunit(iu_LAI13)
             call closeunit(iu_LAI14)
             call closeunit(iu_LAI15)
             call closeunit(iu_LAI16)
#endif !PFT_MODEL_ENT
!          else !mixed_VEG
#else  !MIXED_CANOPY
             call closeunit(iu_LAIMIXED)
!          endif
#endif  !MIXED_CANOPY
       endif

       if ( force_VEG ) then
!         if (.not.mixed_VEG) then
#ifndef MIXED_CANOPY
             call closeunit(iu_height1)
             call closeunit(iu_height2)
             call closeunit(iu_height3)
             call closeunit(iu_height4)
             call closeunit(iu_height5)
             call closeunit(iu_height6)
             call closeunit(iu_height7)
             call closeunit(iu_height8)
#ifdef PFT_MODEL_ENT
             call closeunit(iu_height9)
             call closeunit(iu_height10)
             call closeunit(iu_height11)
             call closeunit(iu_height12)
             call closeunit(iu_height13)
             call closeunit(iu_height14)
             call closeunit(iu_height15)
             call closeunit(iu_height16)
#endif !PFT_MODEL_ENT
!          else
#else  !MIXED_CANOPY
             call closeunit(iu_heightmixed)
!          endif
#endif !MIXED_CANOPY
       endif
     
      end subroutine close_forcings_file


      subroutine skip_forcings( nskip, force_VEG, N_CASA_LAYERS )
      integer, intent(in) :: nskip
      logical, intent(in) :: force_VEG
      integer, intent(in) :: N_CASA_LAYERS
      integer i
      
      if (.not.force_VEG) then
      do i=1,nskip
        read(iu_CH)
        read(iu_DVIS)
        read(iu_MP1)
        read(iu_MP2)
        read(iu_MP3)
        read(iu_MP4)
        read(iu_MP5)
        read(iu_MP6)
        read(iu_PREC)
        read(iu_PS)
        read(iu_QFOL)
        read(iu_QLAT)
        read(iu_QSEN)
        read(iu_SAT)
        read(iu_SIC1)
        read(iu_SIC2)
        read(iu_SIC3)
        read(iu_SIC4)
        read(iu_SIC5)
        read(iu_SIC6)
        read(iu_TP)
        read(iu_US)
        read(iu_VIS)
        read(iu_VS)
        read(iu_ZEN)
        read(iu_W1)
!        if (N_CASA_LAYERS.eq.2) then
        read(iu_W2)
        read(iu_W3)
        read(iu_W4)
        read(iu_W5)
        read(iu_W6)
!        endif
        read(iu_SOILT1)
!        if (N_CASA_LAYERS.eq.2) then
        read(iu_SOILT2)
        read(iu_SOILT3)
        read(iu_SOILT4)
        read(iu_SOILT5)
        read(iu_SOILT6)
!        endif
!        read(iu_SOILT30cm)
!        read(iu_SMOIST30cm)
        if ( force_VEG ) then
          read(iu_LAI1)
          read(iu_LAI2)
          read(iu_LAI3)
          read(iu_LAI4)
          read(iu_LAI5)
          read(iu_LAI6)
          read(iu_LAI7)
          read(iu_LAI8)
#ifdef PFT_MODEL_ENT
          read(iu_LAI9)
          read(iu_LAI10)
          read(iu_LAI11)
          read(iu_LAI12)
          read(iu_LAI13)
          read(iu_LAI14)
          read(iu_LAI15)
          read(iu_LAI16)
#endif
        endif
        if ( force_VEG ) then
          read(iu_height1)
          read(iu_height2)
          read(iu_height3)
          read(iu_height4)
          read(iu_height5)
          read(iu_height6)
          read(iu_height7)
          read(iu_height8)
#ifdef PFT_MODEL_ENT
          read(iu_height9)
          read(iu_height10)
          read(iu_height11)
          read(iu_height12)
          read(iu_height13)
          read(iu_height14)
          read(iu_height15)
          read(iu_height16)
#endif
        endif
      enddo
      endif

      end subroutine skip_forcings


      subroutine read_forcings( I0, I1, J0, J1, 
     &     N_DEPTH, N_CASA_LAYERS,N_PFT,
     &     TairC, TcanopyC, Qf, P_mbar, Ca, Ch, U,
     &     IPARdif, IPARdir, CosZen, 
     &     Soiltemp, Soilmoist, Soilmp, fice, LAI, height, 
     &     force_VEG, do_rewind)!,mixed_VEG
      integer, intent(in) :: I0,I1,J0,J1, N_DEPTH, N_CASA_LAYERS,N_PFT
      real*8, dimension(I0:I1,J0:J1) :: TairC, TcanopyC,Qf,P_mbar,Ca,
     &     Ch, U, IPARdif, IPARdir, CosZen  
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) :: Soilmp, fice
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) :: Soiltemp,Soilmoist
      real*8, dimension(:,:,:),pointer :: LAI, height
      logical :: force_VEG, do_rewind !,mixed_VEG
      !----Local------
      real*4, dimension(I0_file:I1_file,J0_file:J1_file) :: buf, buf1
      real*4, dimension(N_PFT,I0_file:I1_file,J0_file:J1_file) :: vegbuf
      integer niter, niter1, i, j


      if (do_rewind) call rewind_files(force_VEG, N_CASA_LAYERS)

      ! no files for the following data. setting to defaults
      !Qf(I0:I1,J0:J1) = 0.d0
      Ca(I0:I1,J0:J1) = 0.01430229 !(350 ppm @STP in mol m-3) (?old value of 0.0127609 was more like 312 ppm)

!      read(iu_CH) niter, buf
!      print *,'Got here in read_forcings.CH'
!      niter1 = niter
!      Ch(I0:I1,J0:J1) = buf(I0:I1,J0:J1)
      read(iu_CH) niter, buf
      niter1 = niter
      Ch(I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_DVIS) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      IPARdir(I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP1) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP2) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP3) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP4) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP5) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP6) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_PREC) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      != buf(I0:I1,J0:J1)

      read(iu_PS) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      P_mbar(I0:I1,J0:J1) = buf(I0:I1,J0:J1)
!      write(669,*) P_mbar(17,33)

      read(iu_QFOL) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Qf(I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_QLAT) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      != buf(I0:I1,J0:J1)

      read(iu_QSEN) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      != buf(I0:I1,J0:J1)

      read(iu_SAT) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      TairC= buf(I0:I1,J0:J1)

      read(iu_SIC1) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SIC2) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SIC3) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SIC4) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SIC5) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SIC6) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_TP) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      TcanopyC(I0:I1,J0:J1) = buf(I0:I1,J0:J1) 

      read(iu_US) niter, buf1(I0:I1,J0:J1)
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      != buf(I0:I1,J0:J1)

      read(iu_VS) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      != buf(I0:I1,J0:J1)
      !! really should save absolute value
      do j=J0,J1
        do i=I0,I1
          U(i,j) = sqrt( buf1(i,j)**2 + buf(i,j)**2 )
        enddo
      enddo

      read(iu_VIS) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      IPARdif(I0:I1,J0:J1) = buf(I0:I1,J0:J1) - IPARdir(I0:I1,J0:J1)

      read(iu_ZEN) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      CosZen(I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      !***temporary hack for soilmoist and soiltemp*** -PK 7/23/07 
      read(iu_W1) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmoist(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

!      if (N_CASA_LAYERS.eq.2) then
!#ifdef NCASA2
      read(iu_W2) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmoist(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
!#endif
!      endif

      read(iu_W3) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmoist(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_W4) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmoist(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_W5) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmoist(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_W6) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmoist(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SOILT1) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soiltemp(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

!#ifdef NCASA2
      read(iu_SOILT2) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soiltemp(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
!#endif

      read(iu_SOILT3) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soiltemp(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SOILT4) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soiltemp(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SOILT5) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soiltemp(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SOILT6) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soiltemp(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

!* For mixed_VEG, for convenience, we limit the number of prescribed 
!* patches x cohorts not to exceed the size of the LAI and height arrays
!* so that we can just use the same arrays.
!* Instead of LAI1...LAI8 representing veg types 1-8, these are the 
!* cohort LAI values, in order:
!* from oldest to youngest patch, tallest to shortest cohort.
!* So, e.g. for N_PFT=8, there could be 2-3 patches per entcell with
!* 1-3 veg types per patch, not to exceed patches x cohorts = 8.
!* Ditto for N_PFT=16.
      if (force_VEG) then
      if ( associated( LAI )) then !.and.(.not.mixed_VEG) ) then
#ifndef MIXED_CANOPY
        read(iu_LAI1) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI2) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI3) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI4) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI5) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI6) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI7) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(7,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI8) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(8,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
#ifdef PFT_MODEL_ENT
        read(iu_LAI9) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(9,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI10) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(10,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI11) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(11,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI12) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(12,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI13) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(13,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI14) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(14,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI15) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(15,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI16) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(16,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
#endif
!      elseif (associated(LAI).and.(mixed_VEG)) then
#else !MIXED_CANOPY
        read(iu_LAIMIXED) niter, vegbuf
        LAI(:,I0:I1,J0:J1) = vegbuf(:,I0:I1,J0:J1)
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
#endif
      endif
      if ( associated( height )) then !.and.(.not.mixed_VEG) ) then
#ifndef MIXED_CANOPY
        read(iu_height1) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height2) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height3) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height4) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height5) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height6) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height7) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(7,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height8) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(8,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
#ifdef PFT_MODEL_ENT
        read(iu_height9) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(9,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height10) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(10,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height11) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(11,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height12) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(12,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height13) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(13,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height14) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(14,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height15) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(15,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height16) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(16,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
#endif !PFT_MODEL_ENT
!      elseif ( associated( height ).and.mixed_VEG ) then
#else  !MIXED_CANOPY
        read(iu_heightmixed) niter, vegbuf
        height(:,I0:I1,J0:J1) = vegbuf(:,I0:I1,J0:J1)
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
!       endif
#endif !MIXED_CANOPY
      endif !Associated
      endif !force_VEG

      end subroutine read_forcings



      subroutine rewind_files(force_VEG,N_CASA_LAYERS)
      logical,intent(in) ::force_VEG
      integer,intent(in) ::N_CASA_LAYERS
!@sum rewind_files  Reset files back to their beginning to repeat forcings.
      write(*,*) "called rewind_files"
        rewind(iu_CH)
        rewind(iu_DVIS)
        rewind(iu_MP1)
        rewind(iu_MP2)
        rewind(iu_MP3)
        rewind(iu_MP4)
        rewind(iu_MP5)
        rewind(iu_MP6)
        rewind(iu_PREC)
        rewind(iu_PS)
        rewind(iu_QFOL)
        rewind(iu_QLAT)
        rewind(iu_QSEN)
        rewind(iu_SAT)
        rewind(iu_SIC1)
        rewind(iu_SIC2)
        rewind(iu_SIC3)
        rewind(iu_SIC4)
        rewind(iu_SIC5)
        rewind(iu_SIC6)
        rewind(iu_TP)
        rewind(iu_US)
        rewind(iu_VIS)
        rewind(iu_VS)
        rewind(iu_ZEN)
        rewind(iu_W1)
!        if (N_CASA_LAYERS.eq.2) then
        rewind(iu_W2)
        rewind(iu_W3)
        rewind(iu_W4)
        rewind(iu_W5)
        rewind(iu_W6)
!        endif
        rewind(iu_SOILT1)
!        if (N_CASA_LAYERS.eq.2) then
        rewind(iu_SOILT2)
        rewind(iu_SOILT3)
        rewind(iu_SOILT4)
        rewind(iu_SOILT5)
        rewind(iu_SOILT6)
!        endif
!        rewind(iu_SOILT30cm)
!        rewind(iu_SMOIST30cm)
        if ( force_VEG ) then
#ifndef MIXED_CANOPY
          rewind(iu_LAI1)
          rewind(iu_LAI2)
          rewind(iu_LAI3)
          rewind(iu_LAI4)
          rewind(iu_LAI5)
          rewind(iu_LAI6)
          rewind(iu_LAI7)
          rewind(iu_LAI8)
#ifdef PFT_MODEL_ENT
          rewind(iu_LAI9)
          rewind(iu_LAI10)
          rewind(iu_LAI11)
          rewind(iu_LAI12)
          rewind(iu_LAI13)
          rewind(iu_LAI14)
          rewind(iu_LAI15)
          rewind(iu_LAI16)
#endif
#else
!MIXED_CANOPY
          rewind(iu_LAIMIXED)
#endif
!MIXED_CANOPY
        endif
        if ( force_VEG ) then
#ifndef MIXED_CANOPY
          rewind(iu_height1)
          rewind(iu_height2)
          rewind(iu_height3)
          rewind(iu_height4)
          rewind(iu_height5)
          rewind(iu_height6)
          rewind(iu_height7)
          rewind(iu_height8)
#ifdef PFT_MODEL_ENT
          rewind(iu_height9)
          rewind(iu_height10)
          rewind(iu_height11)
          rewind(iu_height12)
          rewind(iu_height13)
          rewind(iu_height14)
          rewind(iu_height15)
          rewind(iu_height16)
#endif
#else
!MIXED_CANOPY
          rewind(iu_heightmixed)
#endif
!MIXED_CANOPY
      endif

      end subroutine rewind_files

      end module ent_forcings
