!@sum Routine to read in meteorological forcings data
!@sum for the NASA GISS Land Surface Model (LSM)
!@sum for a single grid cell with all forcings in a single file.
!@auth M.J. Puma, N.Y. Kiang
!-----------------------------------------------------------------------

      module drv_met

      use CONSTANT, only : undef,pi, 
     &     rhow,lhm,gasc,mair,rgas,tf,sha,grav

      implicit none

      private

      public init_forcings_met, get_forcings_met
      public rewind_files_met, deallocate_drv_met

      !Global variables to module
      integer, parameter :: varnums = 1 !Number of files for met forcings.
      integer,dimension(varnums),save :: iu_vector !Met forcing files.
               
      !ModelE-specific constants.
      real*8,parameter :: CO2ppm = 350.d0! CO2 at land surface [ppm] 
                                         ! 350 ppm = 0.0127609 mol/m3 at STP
      real*8,parameter :: zgs = 10.d0 !height of the surface atmos layer {m]

      !#############################################################!
      !### HARD-CODED SITE PARAMETERS.  SEE ALSO Calc_Ch_site. ###!!!
      !### LATER REPLACE WITH FILE READ ###!!!
      real*8, parameter :: LATD=39.32315
!# Latitude of each site
!#      real*8,parameter :: latd = 39.32315d0  !MMSF
!#      real*8,parameter :: latd = 61.8474d0   !Hyytiala
!#      real*8,parameter :: latd = 38.4067d0   !Vaira
!#      real*8,parameter :: latd = 40.6222d0   !NYC
!#      real*8,parameter :: latd = -2.8567d0   !TNF
!#      real*8,parameter :: latd =  0.3000d0   !Mpala

      contains

!======================================================================!
!NOTE:  get_gswp_forcings uses subroutine COSZT with IM,JM,I0,I1,J0,J1, and 
!       duration of time step, instead of latd input.
!       We may want to replace calc_solarzen with COSZT so that both
!       single and global are consistent.  ##-NK
      subroutine calc_solarzen(td,latdegrees,sbeta1)
      !* Calculate solar zenith angle **in radians**
      !* From Spitters, C. J. T. (1986), AgForMet 38: 231-242.
      use CONSTANT, only : pi
      implicit none
      real*8,intent(in) :: td             ! day(to minute fraction)
      real*8,intent(in) :: latdegrees     ! latitude in degrees
      real*8,parameter :: rad = pi/180.d0 ! Conversion from degrees to radians.
      real*8 :: hour,latrad
      real*8 :: delta                     ! declination angle
      real*8,intent(out) :: sbeta1        ! sbeta1=cos(zen angle)=sin(elev angle)
!      real*8,intent(out) :: solarelev    ! solar elevation angle (rad)
!      real*8,intent(out) :: solarzen     ! solar zenith angle (rad)
      
      hour = (td-floor(td))*24.d0
      latrad = latdegrees*rad
      delta = asin(-sin(rad*23.45d0)*cos(2.d0*pi*(td+10.d0)/365.d0)) 
      sbeta1 = sin(latrad)*sin(delta)+
     &     cos(latrad)*cos(delta)*cos(rad* 15.d0*(hour-12.d0))
      if (sbeta1 < 0.d0) sbeta1 = 0.d0  !**GCM does this too** 

      ! Uncomment to compute the solar elevation and zenth angles below
      !solarelev = asin(sbeta1)  !since asin(cos(zen))=pi/2-zen=elev
      !solarzen = pi/2.d0 - solarelev  

      end subroutine calc_solarzen

!----------------------------------------------------------------------
!##### SITE PARAMS HARD-CODED!!! ALTER BY HAND FOR DIFFERENT SITES#####
      SUBROUTINE calc_Ch_site(CDN_max,Tsurfair,Tgrnd,Psurfair,
     &                   Usurfair,ch)
!@sum  Algorithm based on Mahrt & Ek (1984) and Otles & Gutowski (2005)
!@sum  to compute the atmospheric exchange coefficient (Ch = Cq)
!@sum  with limits based on the values in PBL.f
!@auth MJ Puma, 2008 
      implicit none
      !* Input parameters *!
      real*8,intent(in) :: CDN_max ! neutral drag coef.  [m]
      real*8,intent(in) :: Tsurfair ! atmospheric temperature [K]
      real*8,intent(in) :: Tgrnd ! ground temperature [K]
      real*8,intent(in) :: Psurfair ! air pressure [mb]
      real*8,intent(in) :: Usurfair ! wind speed [m/s]
      !* Local parameters  & variables *!
      real*8,parameter :: kappa = 0.40d0 !von Karman (value  in PBL.f)
      real*8 :: h    !canopy height (m)
      real*8 :: z0m,z0h,z0q ! surface-rough.lngth:momentum,heat,moist[m]
      real*8 :: zd   ! zero-plane displacement height [m]
      real*8 :: Tsp  ! Pot.temperature of the surface atmos.layer
      real*8 :: Tgp  ! Pot.temp of ground@roughness height+ zero plane displ
      real*8 :: Ris  !bulk Richardson number for surface layer
      real*8 :: Dm   !variation of drag coeff for non-neutral stability
      real*8 :: cmn  !Momentum drag coefficient for neutral stability
      real*8 :: chn  !Heat drag coefficient for neutral stability
      real*8 :: cqn  !Moisture drag coefficient for neutral stability
      real*8 :: lzgsbyz0m, lzgsbyz0h, lzgsbyz0q 
      !* Output *!
      real*8,intent(out) :: ch!Atmospheric exchange coeff.(assume ch=cq)

      ! Constants/parameters from PBL.f
      real*8 :: zgs !height of surface atmos layer [m] !!!=10
      real*8,parameter :: smax=0.25d0    ! limit on transfer/drag coef
      real*8,parameter :: smin=0.005d0   ! limit on transfer/drag coef
      real*8,parameter :: cmax=smax*smax ! limit on transfer/drag coef
      real*8,parameter :: cmin=smin*smin ! limit on transfer/drag coef

      ! Compute potential temps at roughness height + zero plane disp)
      Tsp=Tsurfair*(101325.d0/(Psurfair*100.d0))**(gasc/sha)
      Tgp=Tgrnd*(101325.d0/(Psurfair*100.d0))**(gasc/sha) 


      call get_site_params(h,1.5d0,zgs,zd,z0m,z0h,z0q)
      !* For initialization, default LAI=1.5 for Raupach equation.

      ! Atmospheric exchange coefficients for neutral stability
      lzgsbyz0m = log((zgs-zd)/z0m)
      lzgsbyz0h = log((zgs-zd)/z0h)
      lzgsbyz0q = log((zgs-zd)/z0q)
      cmn=kappa*kappa/(lzgsbyz0m**2)
      chn=kappa*kappa/(lzgsbyz0m*lzgsbyz0h)
      cqn=kappa*kappa/(lzgsbyz0m*lzgsbyz0q)

      ! Approximation to the bulk Richardson number for stability
      Ris=((zgs-zd)*grav*(Tsp-Tgp))/(Tgp*(Usurfair**2))

      ! Unstable atmosphere (negative Richardson number)
      if(Ris<0.d0)then
          Dm=1.d0-(15.d0*Ris/(1.d0+75.d0*chn*sqrt(-(zgs/z0m)*Ris)));
          ch=chn*Dm;
      else !Stable atmosphere (positive Richardson number)
          Dm=1.d0/((1.d0+15.d0*Ris)*(sqrt(1.d0+5.d0*Ris)));
          ch=chn*Dm;
      endif
      
!     Limits on ch transfer coefficient from PBL.f
      if (ch > cmax) ch=cmax
      if (ch < cmin) ch=cmin

      return  
      END SUBROUTINE calc_Ch_site

!======================================================================!

      subroutine get_single_data(iu_vector,
     &     I0,I1,J0,J1,
     &         data_srheat          ,
     &         data_trheat          ,
     &         data_ts              ,
     &         data_qs              ,
     &         data_ps              ,
     &         data_ws              ,
     &         data_rainf           ,
     &         data_rainf_c         ,
     &         data_snowf   )

      use filemanager, only : openunit

      integer,intent(in) :: iu_vector(varnums)
      integer,intent(in) :: I0,I1,J0,J1
      real*4,dimension(I0:I1,J0:J1),intent(out) :: 
     &     data_srheat, data_trheat,  data_ts
     &        ,data_qs,     data_ps,      data_ws
     &        ,data_rainf,  data_rainf_c, data_snowf

      read(iu_vector(1),*) data_srheat,data_trheat,data_ts, data_qs,
     &     data_ps,data_ws,data_rainf,data_rainf_c,data_snowf
!      print *, 'in get_single_data'
!      print *, data_srheat,data_trheat,data_ts, data_qs,
!     &     data_ps,data_ws,data_rainf,data_rainf_c,data_snowf


      end subroutine get_single_data

!======================================================================!

      subroutine assign_forcings_single(
     i         I0,I1,J0,J1,
     i         jday,jyear,tyr_sec,dtsec,
     i         tbcs_lsm                           ,
     i         tcanopy_lsm                        ,
     i         data_srheat         ,
     i         data_trheat         ,
     i         data_ts             ,
     i         data_qs             ,
     i         data_ps             ,
     i         data_ws             ,
     i         data_rainf          ,
     i         data_rainf_c        ,
     i         data_snowf          ,
     o         force_Ca             ,
     o         force_cos_zen_angle  ,
     o         force_vis_rad        ,
     o         force_direct_vis_rad ,
     o         force_prec_ms        ,
     o         force_eprec_w        ,
     o         force_sprec_ms       ,
     o         force_seprec_w       ,
     o         force_srheat         ,
     o         force_trheat         ,
     o         force_ts             ,
     o         force_qs             ,
     o         force_ps             ,
     o         force_rhosrf         ,
     o         force_cdh            ,
     o         force_qm1            ,
     o         force_ws             ,
     o         force_pbl_args_ws0   ,
     o         force_pbl_args_tprime,
     o         force_pbl_args_qprime )

      use CONSTANT, only : pi
      use param, only : sync_param
      use lsm_phys_util, only : CDN_max_grid, fdiffuse

      implicit none
      integer, intent(in) :: I0,I1,J0,J1
      integer, intent(in) :: jday      ! day of year
      integer, intent(in) :: jyear     ! year A.D.
      integer, intent(in) :: tyr_sec   ! seconds into year
      real*8, intent(in) ::  dtsec     ! timestep size [sec]
!      real*8,  intent(in) :: latd

      !-----Local variables--------
      integer, parameter :: sec_1day = 86400 ! # of seconds in 1 day
      real*8 :: jtime      ! fraction of day
      real*8:: nday      ! # of timesteps in a day
      real*8 :: td ! time of day of year
      real*8 :: cos_zen_angle ! cos( solar zenith angle)
      real*8 :: rot1,rot2 ! beginning and end of timestep for COSZT
      real*8 :: z0,z0_topo,z0_veg ! surface roughness [m]
      real*8 :: Tgrnd ! ground temperature [deg K] 
      real*8 :: RH ! relative humidity
      real*8 :: disp ! zero plane displacement height [m]
      real*8 :: Ch ! heat transfer coefficient
!     GSWP2 forcings
      real*4, dimension(I0:I1,J0:J1) ::
     &         data_srheat, data_trheat,  data_ts
     &        ,data_qs,     data_ps,      data_ws
     &        ,data_rainf,  data_rainf_c, data_snowf
      real*8, dimension(I0:I1,J0:J1) :: fdif ! frac of radiation that is diffuse
      real*8 :: Rg ! global radiation [W/m2]
      real*8, dimension(I0:I1,J0:J1) ::    ! GISS LSM forcings
     &         tbcs_lsm              ,tcanopy_lsm
     &        ,force_Ca             ,force_cos_zen_angle  
     &        ,force_vis_rad        ,force_direct_vis_rad 
     &        ,force_prec_ms        ,force_eprec_w        
     &        ,force_sprec_ms       ,force_seprec_w
     &        ,force_srheat         ,force_trheat
     &        ,force_ts             ,force_qs
     &        ,force_ps             ,force_rhosrf
     &        ,force_cdh            ,force_qm1
     &        ,force_ws             ,force_pbl_args_ws0   
     &        ,force_pbl_args_tprime,force_pbl_args_qprime
      integer :: i,j
      character*80 :: title
      character*120 :: neutral_drag_file="CD_coef"
      integer :: iu_nd

!*****************************************************************************
!     Calculate mean cosine of zenith angle & fraction of diffuse radiation
!      for the current timestep      
      JTIME = real ((tyr_sec - (jday-1)*86400)/86400.d0) !Fraction of day
      nday = 86400.d0/dtsec
      rot1=(2d0*pi*JTIME)/nday ! beginning of timestep
      rot2=rot1+2d0*pi*dtsec/dble(sec_1day) ! end of timestep
      td = DBLE(jday)+JTIME/nday
      call calc_solarzen(td,LATD,cos_zen_angle)
      Rg = data_srheat(1,1)
      fdif = fdiffuse(cos_zen_angle,td,Rg)
!*****************************************************************************
      loop_latitude_j: do j=J0,J1
         loop_longitude_i: do i = J0,J1
            ! data_ts as landmask
            if (data_ts(i,j) > 0.d0) then
              ! 1) Cosine solar zenith angle
               force_cos_zen_angle(i,j)= cos_zen_angle
              ! 2) Visible radiation [W/m2] 
              !    assume 50% of SW rad.is visible (e.g. Friend 1995 for PAR)
               force_vis_rad(i,j)= 0.5d0*DBLE(data_srheat(i,j))
              ! 3) Direct visible radiation [W/m2]
              force_direct_vis_rad(i,j)=(1-fdif(i,j))*force_vis_rad(i,j)
              ! 4) Total precip.(rain+snow) both convective & large-scale [m/s]
              !    Data [kg/m2/s]
               force_prec_ms(i,j)= ( DBLE(data_rainf(i,j))
     &                              +DBLE(data_snowf(i,j)))/ rhow
              ! 5) Energy of precip.[W/m2]: 0 for rain; lhm*prec for snow
               force_eprec_w(i,j)= -DBLE(data_snowf(i,j))*lhm
              ! 6) Large-scale precipitation [m/s]: currently 0 in modelE
               force_sprec_ms(i,j)= 0.d0 
                  !should be DBLE(data_rainf(i,j)-DBLE(data_rainf_c(i,j)
              ! 7) Energy of large-scale precipitation [W/m2]
               force_seprec_w(i,j)= 0.d0
              ! 8) Incoming shortwave radiation [W/m2]
               force_srheat(i,j)= DBLE(data_srheat(i,j))
              ! 9) Incoming longwave radiation [W/m2]
               force_trheat(i,j)= DBLE(data_trheat(i,j))
              !10) Surface air temperature [K]
               force_ts(i,j)= DBLE(data_ts(i,j))
              !11) Surface air moisture [kg/kg]
               if(data_qs(i,j)<0.0001) data_qs(i,j)=0.0001 
               force_qs(i,j)= DBLE(data_qs(i,j))/
     &                             (1.d0-DBLE(data_qs(i,j)))
              !12) Surface pressure [mb]; data_ps [Pa]
               force_ps(i,j)= DBLE(data_ps(i,j))/100.d0 
              !13) Surface air density [kg/m3]
               force_rhosrf(i,j)= data_ps(i,j)/(rgas*data_ts(i,j))
              !14) Turbulent transfer, Ch (see after 21)
              !15) Amount of water in the 1st atm layer; set to infinity
               force_qm1(i,j)= 1e30 
              !16) Wind speed [m/s]; data_ws [m/s]
               if(data_ws(i,j)<0.01d0) data_ws(i,j)=0.01d0 
               force_ws(i,j)= DBLE(data_ws(i,j))
              !17) Atmospheric CO2 concentration at the land surface [mol/m3]
                !350 ppm = 0.0127609 mol/m3 at STP
               force_Ca(i,j)= CO2ppm*(1.0D-06)*data_ps(i,j)
     &                           /gasc/data_ts(i,j)
              !18) - 21) Planetary boundary layer variables 
               force_pbl_args_ws0 (i,j)    = 0.d0
               force_pbl_args_tprime (i,j) = 0.d0
               force_pbl_args_qprime (i,j) = 0.d0

              !14) Turbulent transfer / humidity transfer coefficient [-]
              !    humidity transfer coef = heat transfer coef (Hansen etal 83)
               !   Surface temperature, Tgrnd, either canopy or ground temp.  
               !   Need Tgrnd to calculate Ch.
               if (tcanopy_lsm(i,j)>-1d30)then
                  Tgrnd = tcanopy_lsm(i,j) + tf ! set to canopy temp; C to K
               elseif (tbcs_lsm(i,j)>-1d30)then
                  Tgrnd = tbcs_lsm(i,j) + tf    ! set to ground temp; C to K
               else ! first timestep
                  Tgrnd = force_ts(i,j)
               endif

               call calc_Ch_site(DBLE(CDN_max_grid(i,j)), 
     &              force_ts(i,j),Tgrnd,force_ps(i,j),force_ws(i,j), Ch)
               force_cdh(i,j)= Ch
!               write(3000,20) i, j, Ch
! 20            format(i10,i10,f15.10)

            else
               force_cos_zen_angle(i,j)= undef
               force_vis_rad(i,j)= undef
               force_direct_vis_rad(i,j)= undef
               force_prec_ms(i,j)= undef
               force_eprec_w(i,j)= undef
               force_sprec_ms(i,j)= undef
               force_seprec_w(i,j)= undef
               force_srheat(i,j)= undef
               force_trheat(i,j)= undef
               force_ts(i,j)= undef
               force_qs(i,j)= undef
               force_ps(i,j)= undef
               force_rhosrf(i,j)= undef
               force_cdh(i,j)= undef
               force_qm1(i,j)= undef
               force_ws(i,j)= undef
               force_Ca(i,j)= undef
               force_pbl_args_ws0 (i,j)    = undef
               force_pbl_args_tprime (i,j) = undef
               force_pbl_args_qprime (i,j) = undef
            endif
         enddo loop_longitude_i
      enddo loop_latitude_j

      end subroutine assign_forcings_single

!-------------------------------------------------------------------
      subroutine open_single_files(iu_vector)
      use filemanager, only : openunit

      implicit none

      character*80 :: filename
      integer,intent(inout) :: iu_vector(varnums)

      filename ='site_forcings'

!     open forcings file
      call openunit(trim(filename),iu_vector(1),.false.,.true.)
      open(iu_vector(1),FILE=filename, STATUS='OLD')

      end subroutine open_single_files

!-------------------------------------------------------------------
      subroutine get_site_params(h,LAI,zgs,zd,z0m,z0h,z0q)
      real*8 :: h               !(m) canopy height
      real*8 :: LAI             !(m2/m2) leaf area index
      real*8 :: zgs             !(m) reference height (eddy flux tower height)
      real*8 :: zd              !(m) displacement height
      real*8 :: z0m             !(m) roughness length for momentum
      real*8 :: z0h             !(m) roughness length for heat
      real*8 :: z0q             !(m) roughness length for humidity
      !----Local-------------

!      Santarem km 83, Raupach (1994) BLM, 71:211-216
      !Constant parameters
      real*8,parameter:: vonKarman = 0.4 !von Karman constant
      real*8,parameter:: CR = 0.3  !roughness-element drag coefficient
      real*8,parameter:: CS = 0.003 !substrate-surface drag coefficient
      real*8,parameter:: Psih = 0.193 !Equation 5, for cw=(zw-d)/(h-d)=2
                  !roughmol/m3ness-sublayer influence function,
                  !departure of velocity profile from the inertial-sublayer
                  !logarithmic relation
      real*8,parameter:: ustarUhmaxratio = 0.3 !(ustar/Uh)max
      real*8,parameter:: cd1 = 7.5  !parameter to fit Eq. 8 for d/h
      !------Local variables------
      real*8 :: ustarUh

!     Santarem km 83, Raupach (1994) BLM, 71:211-216
      h = 35.d0                !35-40 m, S.D.Miller et al 2004
      zgs = 64.d0              !S.D.Miller et al 2004
!     Assumptions.
!     Equation 7
      if (LAI.eq.0.d0) then
         zd = 0.d0              !Better if input is LAI + SAI.
      else
         zd = h*(1.- ((1.-exp(-sqrt(cd1*LAI)))/sqrt(cd1*LAI))) !For LAI>0.
      endif
!     zd =  0.7 * h            !Generic tree canopy, Kaimal & Finnegan 1994
!                              !and Wieringa 1993 review. Equivalent to LAI=1.5.
!     Equation 8
      ustarUh = min( sqrt(CS + CR*LAI/2), ustarUhmaxratio)
!     Equation 4
      z0m = h*(1-zd/h)*exp(-vonKarman/ustarUh - Psih)
      z0h = z0m*exp(-2.d0)     ! exp(-2) = 0.13533528d0
      z0q = z0h

!      Barrow
!       h = 0.25d0 
!       zgs = 4.2d0
!       zd = 0.56d0 * h          !Choudhury et al. 1986 (wheat)
!       z0m = 0.3d0*(0.5d0-zd)   !?
!       z0h = z0m*exp(-2.d0)
!       z0q = z0h 
       
!      Vaira
!      h = 0.5d0 
!      zgs = 4.5d0
!      zd  = 0.56d0 * h         !Choudhury et al. 1986 (wheat)
!      z0m = 0.3d0*(0.5d0-zd)   !?
!      z0h = z0m*exp(-2.d0)     !Similarity theory, Garratt&Hicks(1973)
!      z0q = z0h                !?

!      Tonzi 
!      h = 7.1d0
!      zgs = 23.d0
!      zd  = 0.7 * h        !Miranda et al. 1997 (cerrado) 
!                           !(Kiang 2002 diss., Wiering 1993 review)
!      z0m = 0.9d0          !Kiang diss. 2002, z0m=0.13*h
!      z0h = z0m/7.d0       !?
!      z0h = z0m/25.d0      !Kiang diss. 2002 for sparse canopies.
!      z0q = z0h            !?

!!      Hyytiala
!      h = 13.3d0
!      zgs = 23.3d0
!      zd  = 0.78d0 * h     !Jarvis et al (1976) (coniferous)
!      z0m = 0.075d0* h     !Jarvis et al (1976) (coniferous)
!      z0h = z0m*exp(-2.d0) !Similarity theory, Garratt&Hicks(1973)
!      z0q = z0h

!!     MMSF data
!      h = 26.d0            !(m)
!      zgs = 46.d0
!      zd  = 0.8d0*h        !0.8 from SOURCE?
!      z0m = 2.1d0          !Schmid et al. (2000)
!      z0h = z0m*exp(-2.d0) !Similarity theory, Garratt&Hicks(1973)
!      z0q = z0h            !?


      end subroutine get_site_params

!======================================================================
      subroutine init_forcings_met(
!@sum init_forcings_met  Interface routine with lsm_standalone.      
     &     IM,JM,I0,I1,J0,J1,i0f,i1f,j0f,j1f,
     &     jday,jyear,tyr_sec,dtsec,
     &     tbcs_lsm,tcanopy_lsm,
     &     force_Ca             ,     force_cos_zen_angle  ,
     &     force_vis_rad        ,     force_direct_vis_rad ,
     &     force_prec_ms        ,     force_eprec_w        ,
     &     force_sprec_ms       ,     force_seprec_w       ,
     &     force_srheat         ,     force_trheat         ,
     &     force_ts             ,     force_qs             ,
     &     force_ps             ,     force_rhosrf         ,
     &     force_cdh            ,     force_qm1            ,
     &     force_ws             ,     force_pbl_args_ws0   ,
     &     force_pbl_args_tprime,     force_pbl_args_qprime)

      use filemanager, only : closeunit
      use lsm_phys_util, only : get_neutral_drag

      implicit none

      integer, intent(in) :: IM,JM !Grid resolution
      integer, intent(in) :: I0,I1,J0,J1 !Run boundaries
      integer, intent(in) :: i0f,i1f,j0f,j1f !Input file boundaries.
      integer, intent(in) :: jday      ! day of year
      integer, intent(in) :: jyear       ! year A.D.
      integer, intent(in) :: tyr_sec   ! time from start of year [sec]
      real*8, intent(in) :: dtsec     ! timestep size [sec]
!      real*8,  intent(in) :: latd     !## REPLACE WITH I0,I1,J0,J1
      real*8,dimension(I0:I1,J0:J1), intent(in) :: ! Prognostic vars from LSM
     &     tbcs_lsm,tcanopy_lsm
      real*8,dimension(I0:I1,J0:J1), intent(in) :: ! Forcings to LSM
     &     force_Ca             ,     force_cos_zen_angle  ,
     &     force_vis_rad        ,     force_direct_vis_rad ,
     &     force_prec_ms        ,     force_eprec_w        ,
     &     force_sprec_ms       ,     force_seprec_w       ,
     &     force_srheat         ,     force_trheat         ,
     &     force_ts             ,     force_qs             ,
     &     force_ps             ,     force_rhosrf         ,
     &     force_cdh            ,     force_qm1            ,
     &     force_ws             ,     force_pbl_args_ws0   ,
     &     force_pbl_args_tprime,     force_pbl_args_qprime
      !------- Local ---------

!      integer, parameter :: varnums = 1
!      integer,dimension(varnums),intent(out) :: iu_vector !Global

      real*4, dimension(I0:I1,J0:J1) :: ! Forcings data from files
     &         data_srheat, data_trheat,  data_ts
     &        ,data_qs,     data_ps,      data_ws
     &        ,data_rainf,  data_rainf_c, data_snowf
      !----------------

      call get_neutral_drag(1,IM,1,JM) !Get CDN_max_grid

      call open_single_files(iu_vector)

      call get_single_data(iu_vector
     i      ,I0,I1,J0,J1
     &      ,data_srheat, data_trheat,  data_ts
     &      ,data_qs,     data_ps,      data_ws
     &      ,data_rainf,  data_rainf_c, data_snowf)

!     Obtain forcings for initial timestep -- will become current timestep
      call assign_forcings_single(I0,I1,J0,J1
     i     ,jday,jyear,tyr_sec,dtsec
     i     ,tbcs_lsm               ,tcanopy_lsm
     i     ,data_srheat            , data_trheat
     i     ,data_ts, data_qs, data_ps, data_ws
     i     ,data_rainf , data_rainf_c, data_snowf 
     o     ,force_Ca             ,force_cos_zen_angle  
     o     ,force_vis_rad        ,force_direct_vis_rad 
     o     ,force_prec_ms        ,force_eprec_w        
     o     ,force_sprec_ms       ,force_seprec_w
     o     ,force_srheat         ,force_trheat
     o     ,force_ts             ,force_qs
     o     ,force_ps             ,force_rhosrf
     o     ,force_cdh            ,force_qm1
     o     ,force_ws             ,force_pbl_args_ws0   
     o     ,force_pbl_args_tprime,force_pbl_args_qprime)

      print *, '********END OF INITIAL ASSIGNMENT**************'
      end subroutine init_forcings_met


!======================================================================
      subroutine get_forcings_met(
!@sum get_forcings_met   Interface routine with lsm_standalone.
     i     IM,JM,I0,I1,J0,J1, i0d,i1d,j0d,j1d,
     i     jday, jyear, tyr_sec, dtsec,
     &     tbcs_lsm              ,         tcanopy_lsm           ,   
     &     force_Ca              ,         force_cos_zen_angle   ,
     &     force_vis_rad         ,         force_direct_vis_rad  ,
     &     force_prec_ms         ,         force_eprec_w         ,
     &     force_sprec_ms        ,         force_seprec_w        ,
     &     force_srheat          ,         force_trheat          ,
     &     force_ts              ,         force_qs              ,
     &     force_ps              ,         force_rhosrf          ,
     &     force_cdh             ,         force_qm1             ,
     &     force_ws              ,         force_pbl_args_ws0    ,
     &     force_pbl_args_tprime ,         force_pbl_args_qprime )

      integer, intent(in) :: IM, JM      !Grid resolution
      integer, intent(in) :: I0,I1,J0,J1 !Run grid bounds.
      integer, intent(in) :: i0d,i1d,j0d,j1d !Driver file grid bounds
!      real*8,  intent(in) :: latd
      integer, intent(in) :: jday      ! day of year
      integer, intent(in) :: jyear       ! year A.D.
      integer, intent(in) :: tyr_sec   ! time from start of year [sec]
      real*8,  intent(in) :: dtsec     ! timestep size of simulation[sec]

      !integer,intent(inout) :: iu_vector !Global to module.
      real*8, intent(in),dimension(I0:I1,J0:J1) :: tbcs_lsm, tcanopy_lsm

      real*4, dimension(I0:I1,J0:J1) :: ! GSWP forcings data
     &         data_srheat, data_trheat,  data_ts
     &        ,data_qs,     data_ps,      data_ws
     &        ,data_rainf,  data_rainf_c, data_snowf

      real*8, dimension(I0:I1,J0:J1) :: ! Forcings returned for current time
     &         force_Ca             ,force_cos_zen_angle  
     &        ,force_vis_rad        ,force_direct_vis_rad 
     &        ,force_prec_ms        ,force_eprec_w        
     &        ,force_sprec_ms       ,force_seprec_w
     &        ,force_srheat         ,force_trheat
     &        ,force_ts             ,force_qs
     &        ,force_ps             ,force_rhosrf
     &        ,force_cdh            ,force_qm1
     &        ,force_ws             ,force_pbl_args_ws0   
     &        ,force_pbl_args_tprime,force_pbl_args_qprime

      real*8,dimension(I0:I1,J0:J1) :: ! Forcings at time t0
     &         force_Ca_0             ,force_cos_zen_angle_0  
     &        ,force_vis_rad_0        ,force_direct_vis_rad_0 
     &        ,force_prec_ms_0        ,force_eprec_w_0        
     &        ,force_sprec_ms_0       ,force_seprec_w_0
     &        ,force_srheat_0         ,force_trheat_0
     &        ,force_ts_0             ,force_qs_0
     &        ,force_ps_0             ,force_rhosrf_0
     &        ,force_cdh_0            ,force_qm1_0
     &        ,force_ws_0             ,force_pbl_args_ws0_0   
     &        ,force_pbl_args_tprime_0,force_pbl_args_qprime_0
      
      call get_single_data(iu_vector,
     i     I0,I1,J0,J1,
     o         data_srheat (I0:I1,J0:J1)          ,
     o         data_trheat (I0:I1,J0:J1)          ,
     o         data_ts (I0:I1,J0:J1)              ,
     o         data_qs (I0:I1,J0:J1)              ,
     o         data_ps (I0:I1,J0:J1)              ,
     o         data_ws (I0:I1,J0:J1)              ,
     o         data_rainf( I0:I1,J0:J1)           ,
     o         data_rainf_c (I0:I1,J0:J1)         ,
     o         data_snowf (I0:I1,J0:J1) )
         !print *, 'After get_single_data'

!        Obtain forcings for current timestep

         call assign_forcings_single(
     i         I0,I1,J0,J1,
     i         jday,jyear,tyr_sec,dtsec         ,
     i         tbcs_lsm                           ,
     i         tcanopy_lsm                        ,
     i         data_srheat (I0:I1,J0:J1)         ,
     i         data_trheat (I0:I1,J0:J1)         ,
     i         data_ts (I0:I1,J0:J1)             ,
     i         data_qs (I0:I1,J0:J1)             ,
     i         data_ps (I0:I1,J0:J1)             ,
     i         data_ws (I0:I1,J0:J1)             ,
     i         data_rainf (I0:I1,J0:J1)          ,
     i         data_rainf_c (I0:I1,J0:J1)        ,
     i         data_snowf (I0:I1,J0:J1)          ,
     o         force_Ca (I0:I1,J0:J1)             ,
     o         force_cos_zen_angle (I0:I1,J0:J1)  ,
     o         force_vis_rad (I0:I1,J0:J1)        ,
     o         force_direct_vis_rad (I0:I1,J0:J1) ,
     o         force_prec_ms (I0:I1,J0:J1)        ,
     o         force_eprec_w (I0:I1,J0:J1)        ,
     o         force_sprec_ms (I0:I1,J0:J1)       ,
     o         force_seprec_w (I0:I1,J0:J1)       ,
     o         force_srheat (I0:I1,J0:J1)         ,
     o         force_trheat (I0:I1,J0:J1)         ,
     o         force_ts (I0:I1,J0:J1)             ,
     o         force_qs (I0:I1,J0:J1)             ,
     o         force_ps (I0:I1,J0:J1)             ,
     o         force_rhosrf (I0:I1,J0:J1)         ,
     o         force_cdh (I0:I1,J0:J1)            ,
     o         force_qm1 (I0:I1,J0:J1)            ,
     o         force_ws (I0:I1,J0:J1)             ,
     o         force_pbl_args_ws0 (I0:I1,J0:J1)   ,
     o         force_pbl_args_tprime (I0:I1,J0:J1),
     o         force_pbl_args_qprime (I0:I1,J0:J1) )

      end subroutine get_forcings_met

!------------------------------------------------------------------

      subroutine rewind_files_met
      implicit none
      !integer, intent(in) :: varnums !Global to module
      !integer, intent(in) :: iu_vector(:) !Global to module.
      !-------------
      integer :: n

      do n=1,varnums
         rewind(iu_vector(n))   !Rewind meteorological input files
      enddo

      end subroutine rewind_files_met
!======================================================================
      subroutine deallocate_drv_met
      use filemanager, only : closeunit
      use lsm_phys_util, only : CDN_max_grid
      call closeunit(iu_vector(1))
      deallocate(CDN_max_grid)
      end subroutine deallocate_drv_met
!======================================================================
      end module drv_met
