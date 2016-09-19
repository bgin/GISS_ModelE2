!@sum Routine to read in meteorological forcings data
!@sum from the GSWP2 and from FLUXNET sites for a single grid cell
!@sum for the NASA GISS Land Surface Model (LSM)
!@auth M.J. Puma
!-----------------------------------------------------------------------
!#define RESOL_2x2_5
!#define ECOSYSTEM_SCALE ! Single cell w/FLUXNET met. forcings
!-----------------------------------------------------------------------

      module drv_gswp_force

      use CONSTANT, only : undef,pi, twopi, radian,
     &     rhow,lhm,gasc,mair,rgas,tf,sha,grav

      implicit none
 
      private

      public get_gswp_forcings, init_forcings, get_month

#ifdef ECOSYSTEM_SCALE

      integer, parameter :: im0=1 ! # of long. cells
      integer, parameter :: jm0=1 ! # of lat.  cells
#ifdef RESOL_2x2_5
      real*4, dimension(144,90) :: CDN_temp
#else
      real*4, dimension(72,46) :: CDN_temp
#endif

#else ! GLOBAL RUNS

#ifdef RESOL_2x2_5
      integer, parameter :: im0=144, jm0=90
#else
      integer, parameter :: im0=72, jm0=46
#endif

#endif
      integer, parameter :: num_months = 162 ! # of months of gswp data
      integer, parameter :: varnums = 9      ! # of GSWP2 variables
      integer, parameter :: n_LSM_GSWP = 6   ! # of LSM timesteps per gswp data

      integer, parameter :: sec_1_hr =  3600 ! # of seconds in 1 hr
      integer, parameter :: sec_3_hr = 10800 ! # of seconds in 3 hr
      integer, parameter :: sec_1day = 86400 ! # of seconds in 1 day

      !ModelE-specific constants.
      real*8,parameter :: CO2ppm = 350.d0! CO2 at land surface [ppm]
      real*8,parameter :: zgs = 10.d0 !heigh of the surface atmos layer {m]

      ! Constants/parameters related to heat & humidity transfer coefficient
!USE CONSTANT
!      real*8,parameter :: rho_water = 1d3 ! H2O density @3.98 C (1000kg/m3) CONST.f rhow.
!      real*8,parameter :: lhm = 3.34d5 !latent heat of melting/fusion at 0 C
!      real*8,parameter :: gasc = 8.314510d0 !Ideal gas const(8.314510 J/mol K)
!      real*8,parameter :: mair = 28.9655d0 !molecular weight of dry air [g/mol]
!      real*8,parameter :: rgas = 1d3*gasc/mair !gas constant(287.05 J/K kg)
!      real*8,parameter :: SKIP = -1.E30  !CONST.f undef
!      real*8,parameter :: tfrz = 273.15d0  !Kelvin !##Replace with CONST tf.
      real*8,parameter :: kappa = 0.40d0 !von Karman (value  in PBL.f)
!      real*8,parameter :: cp=1012.d0!Heat capacity of dry air[J kg-1 K-1] ##THIS IS HIGH, for > 100 Celsius. CONST.f sha.
!      real*8,parameter :: grav = 9.8067d0!gravity accel., [m s-2]
!     RGAS = R/M_A = 1000* 8.314510 J/mol K /28.9655 g/mol
!     For values of CO2 much larger than present day (> 4x conc)
!     the molecular weight of dry air M_A could change.
!     Assume that M_O2 = 31.9988 and M_CO2 = 44.00995
!     and current percentages 20.946% and 0.0350% (US Stand. Atm.)
!     Assuming CO2 displaces other gases equally M_A=28.9602 + n*0.00527
!     where n is multiple of present day CO2 conc (350 ppm)
!     For 4xCO2  M_A = 28.9813  => rgas = 286.89
!     For 10xCO2 M_A = 29.0129  => rgas = 286.58

      real*8,save,dimension(im0,jm0) :: ! Forcings at time t1
     &         force_Ca_1             ,force_cos_zen_angle_1  
     &        ,force_vis_rad_1        ,force_direct_vis_rad_1 
     &        ,force_prec_ms_1        ,force_eprec_w_1        
     &        ,force_sprec_ms_1       ,force_seprec_w_1
     &        ,force_srheat_1         ,force_trheat_1
     &        ,force_ts_1             ,force_qs_1
     &        ,force_ps_1             ,force_rhosrf_1
     &        ,force_cdh_1            ,force_qm1_1
     &        ,force_ws_1             ,force_pbl_args_ws0_1   
     &        ,force_pbl_args_tprime_1,force_pbl_args_qprime_1

      real*8,save,dimension(im0,jm0) :: ! Forcings at time t2
     &         force_Ca_2             ,force_cos_zen_angle_2  
     &        ,force_vis_rad_2        ,force_direct_vis_rad_2 
     &        ,force_prec_ms_2        ,force_eprec_w_2        
     &        ,force_sprec_ms_2       ,force_seprec_w_2
     &        ,force_srheat_2         ,force_trheat_2
     &        ,force_ts_2             ,force_qs_2
     &        ,force_ps_2             ,force_rhosrf_2
     &        ,force_cdh_2            ,force_qm1_2
     &        ,force_ws_2             ,force_pbl_args_ws0_2   
     &        ,force_pbl_args_tprime_2,force_pbl_args_qprime_2


      real*8,save,dimension(im0,jm0,n_LSM_GSWP) ::! Interpol. sub-3hr forcings
     &         vector_Ca             ,vector_cos_zen_angle  
     &        ,vector_vis_rad        ,vector_direct_vis_rad 
     &        ,vector_prec_ms        ,vector_eprec_w        
     &        ,vector_sprec_ms       ,vector_seprec_w
     &        ,vector_srheat         ,vector_trheat
     &        ,vector_ts             ,vector_qs
     &        ,vector_ps             ,vector_rhosrf
     &        ,vector_cdh            ,vector_qm1
     &        ,vector_ws             ,vector_pbl_args_ws0   
     &        ,vector_pbl_args_tprime,vector_pbl_args_qprime

      contains
!----------------------------------------------------------------------

      function LEAPYR(year)
      integer :: year
      logical :: LEAPYR
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

      subroutine calc_solarzen(jday,td,latdegrees,sbeta1)
      !* Calculate solar zenith angle **in radians**
      !* From Spitters, C. J. T. (1986), AgForMet 38: 231-242.
      implicit none
      integer,intent(in) :: jday          ! day of year
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

      function fdiffuse(cos_solzen,day,Rg) Result(fdif)
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

!======================================================================!

      subroutine declination(gday,gyr,sind,cosd)
!      Input: ECCEN  = eccentricity of the orbital ellipse
!             OBLIQ  = latitude of Tropic of Cancer
!             OMEGVP = longitude of perihelion (sometimes Pi is added) =
!                    = spatial angle from vernal equinox to perihelion
!                      with Sun as angle vertex
!             DAY    = days measured since 2000 January 1, hour 0
!     
!             EDPY  = Earth days per year
!                     tropical year = 365.2425 (Gregorgian Calendar)
!                     tropical year = 365      (Generic Year)
!             VEDAY = Vernal equinox
!                     79.0 (Generic year Mar 21 hour 0)
!                     79.3125d0 for days from 2000 January 1, hr 0 till vernal
!                          equinox of year 2000 = 31 + 29 + 19 + 7.5/24
!      Output: SIND = sine of declination angle = sin(SUNLAT)
!              COSD = cosine of the declination angle = cos(SUNLAT)

!      Reference for following caculations is:  V.M.Blanco and
!      S.W.McCuskey, 1961, "Basic Physics of the Solar System", pages
!      135 - 151.  Existence of Moon and heavenly bodies other than
!      Earth and Sun are ignored.  Earth is assumed to be spherical.

!@var SIND,COSD orbit related variables computed once a day
!             OBLIQ  = latitude of Tropic of Cancer 
!            TA = true anomaly = spatial angle from perihelion to
!                 current location with Sun as angle vertex
!        TAofVE = TA(VE) = true anomaly of vernal equinox = - OMEGVP 
      integer, intent(in) :: gday      ! gregorian day of year
      integer, intent(in) :: gyr       ! gregorian year

      real*8, parameter :: EDPY=365.2425d0, VEDAY=79d0 
      real*8 :: DOBLIQ,ECCEN,DOMEGVP,DAY
      real*8 :: SIND,COSD

      real*8 :: MA,OMEGVP,OBLIQ,EA,DEA,BSEMI
     &     ,TAofVE,EAofVE,MAofVE,SUNDIS,TA,SUNX,SUNY,SLNORO
     &     ,VEQLON,ROTATE,SUNLON,SUNLAT,SLMEAN
      real*8 ::  EDAYzY,VE2000,year

      DAY  = dble(gday)
      year = dble(gyr)
      call ORBPAR(year,ECCEN,DOBLIQ,DOMEGVP)
      VE2000=VEDAY
      EDAYzY=EDPY
      OMEGVP=DOMEGVP*radian
      OBLIQ=DOBLIQ*radian

!      Determine EAofVE from geometry: tan(EA) = b*sin(TA) / [e+cos(TA)]
!      Determine MAofVE from Kepler's equation: MA = EA - e*sin(EA)
!      Determine MA knowing time from vernal equinox to current day
      BSEMI  = SQRT (1 - ECCEN*ECCEN)
      TAofVE = - OMEGVP
      EAofVE = ATAN2 (BSEMI*SIN(TAofVE), ECCEN+COS(TAofVE))
      MAofVE = EAofVE - ECCEN*SIN(EAofVE)
!      PERIHE = VE2000 - MAofVE*EDAYzY/TWOPI
      MA     = MODULO (TWOPI*(DAY-VE2000)/EDAYzY + MAofVE, TWOPI)
!     
!      Numerically invert Kepler's equation: MA = EA - e*sin(EA)     
      EA  = MA + ECCEN*(SIN(MA) + ECCEN*SIN(2*MA)/2)
   10 dEA = (MA - EA + ECCEN*SIN(EA)) / (1 - ECCEN*COS(EA))
      EA  = EA + dEA
      IF(ABS(dEA).gt.1d-10)  GO TO 10
!     
!      Calculate distance to Sun and true anomaly
      SUNDIS = 1 - ECCEN*COS(EA)
      TA     = ATAN2 (BSEMI*SIN(EA), COS(EA)-ECCEN)
      SIND   = SIN(TA-TAofVE) * SIN(OBLIQ)
      COSD   = SQRT (1 - SIND*SIND)
      
      end subroutine declination

!======================================================================!

      SUBROUTINE ORBPAR (YEAR, ECCEN,OBLIQ,OMEGVP)
!@sum ORBPAR calculates the three orbital parameters as a function of
!@+   YEAR.  The source of these calculations is: Andre L. Berger,
!@+   1978, "Long-Term Variations of Daily Insolation and Quaternary
!@+   Climatic Changes", JAS, v.35, p.2362.  Also useful is: Andre L.
!@+   Berger, May 1978, "A Simple Algorithm to Compute Long Term
!@+   Variations of Daily Insolation", published by Institut
!@+   D'Astronomie de Geophysique, Universite Catholique de Louvain,
!@+   Louvain-la Neuve, No. 18.
!@auth Gary Russell (with extra terms from D. Thresher)
!@auth Michael J. Puma (modified June 2008)
C****
C**** Tables and equations refer to the first reference (JAS).  The
C**** corresponding table or equation in the second reference is
C**** enclosed in parentheses.
C****
      IMPLICIT NONE
C**** Input:
!@var YEAR   = years C.E. are positive, B.C.E are -ve (i.e 4BCE = -3)
      REAL*8, INTENT(IN) :: YEAR
C**** Output:
!@var ECCEN  = eccentricity of orbital ellipse
!@var OBLIQ  = latitude of Tropic of Cancer in degrees
!@var OMEGVP = longitude of perihelion =
!@+          = spatial angle from vernal equinox to perihelion
!@+            in degrees with sun as angle vertex
      REAL*8, INTENT(OUT) :: ECCEN,OBLIQ,OMEGVP
C**** Table 1 (2).  Obliquity relative to mean ecliptic of date: OBLIQD
      REAL*8, PARAMETER, DIMENSION(3,47) :: TABL1 = RESHAPE( (/
     1  -2462.22D0,  31.609970D0,  251.9025D0,
     2   -857.32D0,  32.620499D0,  280.8325D0,
     3   -629.32D0,  24.172195D0,  128.3057D0,
     4   -414.28D0,  31.983780D0,  292.7251D0,
     5   -311.76D0,  44.828339D0,   15.3747D0,
     6    308.94D0,  30.973251D0,  263.7952D0,
     7   -162.55D0,  43.668243D0,  308.4258D0,
     8   -116.11D0,  32.246689D0,  240.0099D0,
     9    101.12D0,  30.599442D0,  222.9725D0,
     O    -67.69D0,  42.681320D0,  268.7810D0,
     1     24.91D0,  43.836456D0,  316.7998D0,
     2     22.58D0,  47.439438D0,  319.6023D0,
     3    -21.16D0,  63.219955D0,  143.8050D0,
     4    -15.65D0,  64.230484D0,  172.7351D0,
     5     15.39D0,   1.010530D0,   28.9300D0,
     6     14.67D0,   7.437771D0,  123.5968D0,
     7    -11.73D0,  55.782181D0,   20.2082D0,
     8     10.27D0,   0.373813D0,   40.8226D0,
     9      6.49D0,  13.218362D0,  123.4722D0,
     O      5.85D0,  62.583237D0,  155.6977D0,
     1     -5.49D0,  63.593765D0,  184.6277D0,
     2     -5.43D0,  76.438309D0,  267.2771D0,
     3      5.16D0,  45.815262D0,   55.0196D0,
     4      5.08D0,   8.448301D0,  152.5268D0,
     5     -4.07D0,  56.792709D0,   49.1382D0,
     6      3.72D0,  49.747849D0,  204.6609D0,
     7      3.40D0,  12.058272D0,   56.5233D0,
     8     -2.83D0,  75.278214D0,  200.3284D0,
     9     -2.66D0,  65.241013D0,  201.6651D0,
     O     -2.57D0,  64.604294D0,  213.5577D0,
     1     -2.47D0,   1.647247D0,   17.0374D0,
     2      2.46D0,   7.811584D0,  164.4194D0,
     3      2.25D0,  12.207832D0,   94.5422D0,
     4     -2.08D0,  63.856659D0,  131.9124D0,
     5     -1.97D0,  56.155991D0,   61.0309D0,
     6     -1.88D0,  77.448837D0,  296.2073D0,
     7     -1.85D0,   6.801054D0,  135.4894D0,
     8      1.82D0,  62.209412D0,  114.8750D0,
     9      1.76D0,  20.656128D0,  247.0691D0,
     O     -1.54D0,  48.344406D0,  256.6113D0,
     1      1.47D0,  55.145462D0,   32.1008D0,
     2     -1.46D0,  69.000534D0,  143.6804D0,
     3      1.42D0,  11.071350D0,   16.8784D0,
     4     -1.18D0,  74.291306D0,  160.6835D0,
     5      1.18D0,  11.047742D0,   27.5932D0,
     6     -1.13D0,   0.636717D0,  348.1074D0,
     7      1.09D0,  12.844549D0,   82.6496D0/), (/3,47/) )
C**** Table 2 (4).  Precessional parameter: ECCEN sin(omega) (unused)
      REAL*8, PARAMETER, DIMENSION(3,46) :: TABL2 = RESHAPE(  (/
     1     .0186080D0,  54.646484D0,   32.012589D0,
     2     .0162752D0,  57.785370D0,  197.181274D0,
     3    -.0130066D0,  68.296539D0,  311.699463D0,
     4     .0098883D0,  67.659821D0,  323.592041D0,
     5    -.0033670D0,  67.286011D0,  282.769531D0,
     6     .0033308D0,  55.638351D0,   90.587509D0,
     7    -.0023540D0,  68.670349D0,  352.522217D0,
     8     .0014002D0,  76.656036D0,  131.835892D0,
     9     .0010070D0,  56.798447D0,  157.536392D0,
     O     .0008570D0,  66.649292D0,  294.662109D0,
     1     .0006499D0,  53.504456D0,  118.253082D0,
     2     .0005990D0,  67.023102D0,  335.484863D0,
     3     .0003780D0,  68.933258D0,  299.806885D0,
     4    -.0003370D0,  56.630219D0,  149.162415D0,
     5     .0003334D0,  86.256454D0,  283.915039D0,
     6     .0003334D0,  23.036499D0,  320.110107D0,
     7     .0002916D0,  89.395340D0,   89.083817D0,
     8     .0002916D0,  26.175385D0,  125.278732D0,
     9     .0002760D0,  69.307068D0,  340.629639D0,
     O    -.0002330D0,  99.906509D0,  203.602081D0,
     1    -.0002330D0,  36.686569D0,  239.796982D0,
     2     .0001820D0,  67.864838D0,  155.484787D0,
     3     .0001772D0,  99.269791D0,  215.494690D0,
     4     .0001772D0,  36.049850D0,  251.689606D0,
     5    -.0001740D0,  56.625275D0,  130.232391D0,
     6    -.0001240D0,  68.856720D0,  214.059708D0,
     7     .0001153D0,  87.266983D0,  312.845215D0,
     8     .0001153D0,  22.025970D0,  291.179932D0,
     9     .0001008D0,  90.405869D0,  118.013870D0,
     O     .0001008D0,  25.164856D0,   96.348694D0,
     1     .0000912D0,  78.818680D0,  160.318298D0,
     2     .0000912D0,  30.474274D0,   83.706894D0,
     3    -.0000806D0, 100.917038D0,  232.532120D0,
     4    -.0000806D0,  35.676025D0,  210.866943D0,
     5     .0000798D0,  81.957565D0,  325.487061D0,
     6     .0000798D0,  33.613159D0,  248.875565D0,
     7    -.0000638D0,  92.468735D0,   80.005234D0,
     8    -.0000638D0,  44.124329D0,    3.393823D0,
     9     .0000612D0, 100.280319D0,  244.424728D0,
     O     .0000612D0,  35.039322D0,  222.759552D0,
     1    -.0000603D0,  98.895981D0,  174.672028D0,
     2    -.0000603D0,  35.676025D0,  210.866943D0,
     3     .0000597D0,  87.248322D0,  342.489990D0,
     4     .0000597D0,  24.028381D0,   18.684967D0,
     5     .0000559D0,  86.630264D0,  324.737793D0,
     6     .0000559D0,  22.662689D0,  279.287354D0/), (/3,46/) )
C**** Table 3 (5).  Eccentricity: ECCEN (unused)
      REAL*8, PARAMETER, DIMENSION(3,42) :: TABL3 = RESHAPE( (/
     1     .01102940D0,   3.138886D0,  165.168686D0,
     2    -.00873296D0,  13.650058D0,  279.687012D0,
     3    -.00749255D0,  10.511172D0,  114.518250D0,
     4     .00672394D0,  13.013341D0,  291.579590D0,
     5     .00581229D0,   9.874455D0,  126.410858D0,
     6    -.00470066D0,   0.636717D0,  348.107422D0,
     7    -.00254464D0,  12.639528D0,  250.756897D0,
     8     .00231485D0,   0.991874D0,   58.574905D0,
     9    -.00221955D0,   9.500642D0,   85.588211D0,
     O     .00201868D0,   2.147012D0,  106.593765D0,
     1    -.00172371D0,   0.373813D0,   40.822647D0,
     2    -.00166112D0,  12.658154D0,  221.112030D0,
     3     .00145096D0,   1.010530D0,   28.930038D0,
     4     .00131342D0,  12.021467D0,  233.004639D0,
     5     .00101442D0,   0.373813D0,   40.822647D0,
     6    -.00088343D0,  14.023871D0,  320.509521D0,
     7    -.00083395D0,   6.277772D0,  330.337402D0,
     8     .00079475D0,   6.277772D0,  330.337402D0,
     9     .00067546D0,  27.300110D0,  199.373871D0,
     O    -.00066447D0,  10.884985D0,  155.340912D0,
     1     .00062591D0,  21.022339D0,  229.036499D0,
     2     .00059751D0,  22.009552D0,   99.823303D0,
     3    -.00053262D0,  27.300110D0,  199.373871D0,
     4    -.00052983D0,   5.641055D0,  342.229980D0,
     5    -.00052983D0,   6.914489D0,  318.444824D0,
     6     .00052836D0,  12.002811D0,  262.649414D0,
     7     .00051457D0,  16.788940D0,   84.855621D0,
     8    -.00050748D0,  11.647654D0,  192.181992D0,
     9    -.00049048D0,  24.535049D0,   75.027847D0,
     O     .00048888D0,  18.870667D0,  294.654541D0,
     1     .00046278D0,  26.026688D0,  223.159103D0,
     2     .00046212D0,   8.863925D0,   97.480820D0,
     3     .00046046D0,  17.162750D0,  125.678268D0,
     4     .00042941D0,   2.151964D0,  125.523788D0,
     5     .00042342D0,  37.174576D0,  325.784668D0,
     6     .00041713D0,  19.748917D0,  252.821732D0,
     7    -.00040745D0,  21.022339D0,  229.036499D0,
     8    -.00040569D0,   3.512699D0,  205.991333D0,
     9    -.00040569D0,   1.765073D0,  124.346024D0,
     O    -.00040385D0,  29.802292D0,   16.435165D0,
     1     .00040274D0,   7.746099D0,  350.172119D0,
     2     .00040068D0,   1.142024D0,  273.759521D0/), (/3,42/) )
C**** Table 4 (1).  Fundamental elements of the ecliptic: ECCEN sin(pi)
      REAL*8, PARAMETER, DIMENSION(3,19) :: TABL4 = RESHAPE( (/
     1     .01860798D0,   4.207205D0,   28.620089D0,
     2     .01627522D0,   7.346091D0,  193.788772D0,
     3    -.01300660D0,  17.857263D0,  308.307024D0,
     4     .00988829D0,  17.220546D0,  320.199637D0,
     5    -.00336700D0,  16.846733D0,  279.376984D0,
     6     .00333077D0,   5.199079D0,   87.195000D0,
     7    -.00235400D0,  18.231076D0,  349.129677D0,
     8     .00140015D0,  26.216758D0,  128.443387D0,
     9     .00100700D0,   6.359169D0,  154.143880D0,
     O     .00085700D0,  16.210016D0,  291.269597D0,
     1     .00064990D0,   3.065181D0,  114.860583D0,
     2     .00059900D0,  16.583829D0,  332.092251D0,
     3     .00037800D0,  18.493980D0,  296.414411D0,
     4    -.00033700D0,   6.190953D0,  145.769910D0,
     5     .00027600D0,  18.867793D0,  337.237063D0,
     6     .00018200D0,  17.425567D0,  152.092288D0,
     7    -.00017400D0,   6.186001D0,  126.839891D0,
     8    -.00012400D0,  18.417441D0,  210.667199D0,
     9     .00001250D0,   0.667863D0,   72.108838D0/), (/3,19/) )
C**** Table 5 (3).  General precession in longitude: psi
      REAL*8, PARAMETER, DIMENSION(3,78) :: TABL5 = RESHAPE( (/
     1    7391.0225890d0,  31.609974d0,   251.9025d0,
     2    2555.1526947d0,  32.620504d0,   280.8325d0,
     3    2022.7629188d0,  24.172203d0,   128.3057d0,
     4   -1973.6517951d0,   0.636717d0,   348.1074d0,
     5    1240.2321818d0,  31.983787d0,   292.7252d0,
     6     953.8679112d0,   3.138886d0,   165.1686d0,
     7    -931.7537108d0,  30.973257d0,   263.7951d0,
     8     872.3795383d0,  44.828336d0,    15.3747d0,
     9     606.3544732d0,   0.991874d0,    58.5749d0,
     O    -496.0274038d0,   0.373813d0,    40.8226d0,
     1     456.9608039d0,  43.668246d0,   308.4258d0,
     2     346.9462320d0,  32.246691d0,   240.0099d0,
     3    -305.8412902d0,  30.599444d0,   222.9725d0,
     4     249.6173246d0,   2.147012d0,   106.5937d0,
     5    -199.1027200d0,  10.511172d0,   114.5182d0,
     6     191.0560889d0,  42.681324d0,   268.7809d0,
     7    -175.2936572d0,  13.650058d0,   279.6869d0,
     8     165.9068833d0,   0.986922d0,    39.6448d0,
     9     161.1285917d0,   9.874455d0,   126.4108d0,
     O     139.7878093d0,  13.013341d0,   291.5795d0,
     1    -133.5228399d0,   0.262904d0,   307.2848d0,
     2     117.0673811d0,   0.004952d0,    18.9300d0,
     3     104.6907281d0,   1.142024d0,   273.7596d0,
     4      95.3227476d0,  63.219948d0,   143.8050d0,
     5      86.7824524d0,   0.205021d0,   191.8927d0,
     6      86.0857729d0,   2.151964d0,   125.5237d0,
     7      70.5893698d0,  64.230478d0,   172.7351d0,
     8     -69.9719343d0,  43.836462d0,   316.7998d0,
     9     -62.5817473d0,  47.439436d0,   319.6024d0,
     O      61.5450059d0,   1.384343d0,    69.7526d0,
     1     -57.9364011d0,   7.437771d0,   123.5968d0,
     2      57.1899832d0,  18.829299d0,   217.6432d0,
     3     -57.0236109d0,   9.500642d0,    85.5882d0,
     4     -54.2119253d0,   0.431696d0,   156.2147d0,
     5      53.2834147d0,   1.160090d0,    66.9489d0,
     6      52.1223575d0,  55.782177d0,    20.2082d0,
     7     -49.0059908d0,  12.639528d0,   250.7568d0,
     8     -48.3118757d0,   1.155138d0,    48.0188d0,
     9     -45.4191685d0,   0.168216d0,     8.3739d0,
     O     -42.2357920d0,   1.647247d0,    17.0374d0,
     1     -34.7971099d0,  10.884985d0,   155.3409d0,
     2      34.4623613d0,   5.610937d0,    94.1709d0,
     3     -33.8356643d0,  12.658184d0,   221.1120d0,
     4      33.6689362d0,   1.010530d0,    28.9300d0,
     5     -31.2521586d0,   1.983748d0,   117.1498d0,
     6     -30.8798701d0,  14.023871d0,   320.5095d0,
     7      28.4640769d0,   0.560178d0,   262.3602d0,
     8     -27.1960802d0,   1.273434d0,   336.2148d0,
     9      27.0860736d0,  12.021467d0,   233.0046d0,
     O     -26.3437456d0,  62.583231d0,   155.6977d0,
     1      24.7253740d0,  63.593761d0,   184.6277d0,
     2      24.6732126d0,  76.438310d0,   267.2772d0,
     3      24.4272733d0,   4.280910d0,    78.9281d0,
     4      24.0127327d0,  13.218362d0,   123.4722d0,
     5      21.7150294d0,  17.818769d0,   188.7132d0,
     6     -21.5375347d0,   8.359495d0,   180.1364d0,
     7      18.1148363d0,  56.792707d0,    49.1382d0,
     8     -16.9603104d0,   8.448301d0,   152.5268d0,
     9     -16.1765215d0,   1.978796d0,    98.2198d0,
     O      15.5567653d0,   8.863925d0,    97.4808d0,
     1      15.4846529d0,   0.186365d0,   221.5376d0,
     2      15.2150632d0,   8.996212d0,   168.2438d0,
     3      14.5047426d0,   6.771027d0,   161.1199d0,
     4     -14.3873316d0,  45.815258d0,    55.0196d0,
     5      13.1351419d0,  12.002811d0,   262.6495d0,
     6      12.8776311d0,  75.278220d0,   200.3284d0,
     7      11.9867234d0,  65.241008d0,   201.6651d0,
     8      11.9385578d0,  18.870667d0,   294.6547d0,
     9      11.7030822d0,  22.009553d0,    99.8233d0,
     O      11.6018181d0,  64.604291d0,   213.5577d0,
     1     -11.2617293d0,  11.498094d0,   154.1631d0,
     2     -10.4664199d0,   0.578834d0,   232.7153d0,
     3      10.4333970d0,   9.237738d0,   138.3034d0,
     4     -10.2377466d0,  49.747842d0,   204.6609d0,
     5      10.1934446d0,   2.147012d0,   106.5938d0,
     6     -10.1280191d0,   1.196895d0,   250.4676d0,
     7      10.0289441d0,   2.133898d0,   332.3345d0,
     8     -10.0034259d0,   0.173168d0,    27.3039d0/), (/3,78/) )
C****
      REAL*8 :: YM1950,SUMC,ARG,ESINPI,ECOSPI,PIE,PSI,FSINFD
      INTEGER :: I
C****
      YM1950 = YEAR-1950.
C****
C**** Obliquity from Table 1 (2):
C****   OBLIQ# = 23.320556 (degrees)             Equation 5.5 (15)
C****   OBLIQ  = OBLIQ# + sum[A cos(ft+delta)]   Equation 1 (5)
C****
      SUMC = 0.
      DO I=1,47
        ARG  = radian*(YM1950*TABL1(2,I)/3600.+TABL1(3,I))
        SUMC = SUMC + TABL1(1,I)*COS(ARG)
      END DO
      OBLIQ = 23.320556D0 + SUMC/3600.
!      OBLIQ  = OBLIQ*radian ! not needed for output in degrees
C****
C**** Eccentricity from Table 4 (1):
C****   ECCEN sin(pi) = sum[M sin(gt+beta)]           Equation 4 (1)
C****   ECCEN cos(pi) = sum[M cos(gt+beta)]           Equation 4 (1)
C****   ECCEN = ECCEN sqrt[sin(pi)^2 + cos(pi)^2]
C****
      ESINPI = 0.
      ECOSPI = 0.
      DO I=1,19
        ARG    = radian*(YM1950*TABL4(2,I)/3600.+TABL4(3,I))
        ESINPI = ESINPI + TABL4(1,I)*SIN(ARG)
        ECOSPI = ECOSPI + TABL4(1,I)*COS(ARG)
      END DO
      ECCEN  = SQRT(ESINPI*ESINPI+ECOSPI*ECOSPI)
C****
C**** Perihelion from Equation 4,6,7 (9) and Table 4,5 (1,3):
C****   PSI# = 50.439273 (seconds of degree)         Equation 7.5 (16)
C****   ZETA =  3.392506 (degrees)                   Equation 7.5 (17)
C****   PSI = PSI# t + ZETA + sum[F sin(ft+delta)]   Equation 7 (9)
C****   PIE = atan[ECCEN sin(pi) / ECCEN cos(pi)]
C****   OMEGVP = PIE + PSI + 3.14159                 Equation 6 (4.5)
C****
      PIE = ATAN2(ESINPI,ECOSPI)
      FSINFD = 0.
      DO I=1,78
        ARG    = radian*(YM1950*TABL5(2,I)/3600.+TABL5(3,I))
        FSINFD = FSINFD + TABL5(1,I)*SIN(ARG)
      END DO
      PSI    = radian*(3.392506D0+(YM1950*50.439273D0+FSINFD)/3600.)
      OMEGVP = MOD(PIE+PSI+.5*TWOPI,TWOPI)
      IF(OMEGVP.lt.0.)  OMEGVP = OMEGVP + TWOPI
      OMEGVP = OMEGVP/radian  ! for output in degrees
C****
      RETURN
      END SUBROUTINE ORBPAR

!======================================================================!
#ifndef ECOSYSTEM_SCALE
      subroutine COSZT (gday,gyr,rot1,rot2,Rg,COSZ,fdif)
!@sum  COSZ0 calculates Earth's zenith angle, weighted by time/sunlight
!@sum  THIS ENTRY COMPUTES THE COSINE OF THE ZENITH ANGLE WEIGHTED BY DAYTIME
!@sum  HOURS FROM rot1 TO rot2, GREENWICH MEAN TIME IN RADIANS.  rot1
!@sum  MUST BE BETWEEN 0 AND 2*PI.  rot2 MUST BE BETWEEN rot1 AND
!@sum  rot1+2*PI.  I=1 MUST LIE ON THE INTERNATIONAL DATE LINE.
!@auth Modified by Michael J. Puma (June 2008)
      integer, intent(in) :: gday                  ! gregorian day of year
      integer, intent(in) :: gyr                   ! gregorian year
      real*8,dimension(im0,jm0) ,intent(in) :: Rg  ! global radiation [W/m2]
!@param FIM,BYIM real values related to number of long. grid boxes
      real*8, parameter :: FIM=im0
      real*8, parameter :: BYIM=1./FIM
!@param  DLON grid spacing in longitude (deg)
      real*8, parameter :: DLON = TWOPI*BYIM

      SAVE
!@var  LON longitude of mid points of primary grid box (radians)
      real*8, dimension(im0) :: LON
!@param  DLAT,DLAT_DG grid spacing in latitude (rad/deg)
      real*8  :: DLAT,DLAT_DG
!@var SINIP,COSIP longitud. sin,cos for pressure grid
      real*8, dimension(im0) :: SINIP,COSIP
!@var SINJ,COSJ sines and cosines for zenith angle calculation
      real*8, dimension(jm0) :: SINJ,COSJ
      real*8 :: PHIS,SPHIS,CPHIS,PHIN,SPHIN,CPHIN,PHIM
!@var SIND = sine of declination angle = sin(SUNLAT)
!@var COSD = cosine of the declination angle = cos(SUNLAT)
      real*8 :: SIND,COSD
      real*8 :: rot1,rot2,drot
      real*8, dimension(im0,jm0) :: COSZ,fdif
      real*8, dimension(im0) :: LT1,LT2,SLT1,SLT2,S2LT1,S2LT2

!     ZERO1 HAS TO EQUAL THE CUT-OFF VALUE FOR COSZ USED IN SOLAR
!     COSZS WORKS CORRECTLY ONLY IF ZERO1 >> 1.D-3
      real*8, parameter :: ZERO1=1.D-2
      real*8 :: S2DAWN,S2DUSK,ECOSZ,ECOSQZ,CLT1,CLT2,ZERO2,CDUSK,DUSK
     &     ,DAWN,SDUSK,SDAWN,CJCD,SJSD,SR1,CR1,SR2,CR2
      integer :: I,J
!!! hack to prevent NaN's
      COSZ(:,:) = 0.d0

!----------------------------------------------------------------------
!Gets dlat,dlon,lon,sinip,cosip       
      DLAT_DG=180./jm0                     ! even spacing (default)
      IF (jm0.eq.46) DLAT_DG=180./(jm0-1)   ! 1/2 box at pole for 4x5
      DLAT=DLAT_DG*radian
!      Set indexes and scalings for the influence of A grid points on
!      adjacent velocity points
!      Calculate relative directions of polar box to nearby U,V points
      DO I=1,im0
        LON(I)=DLON*(I-.5)
        SINIP(I)=SIN(LON(I))
        COSIP(I)=COS(LON(I))
      END DO   
!----------------------------------------------------------------------
!Gets sinj and cosj
!     COMPUTE THE AREA WEIGHTED LATITUDES AND THEIR SINES AND COSINES
!     SINJ,COSJ sines and cosines for zenith angle calculation
      PHIS=-.25*TWOPI
      SPHIS=-1.
      CPHIS=0.
      DO J=1,jm0-1
        PHIN=DLAT*(J-.5*jm0)
        SPHIN=SIN(PHIN)
        CPHIN=COS(PHIN)
        PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
        SINJ(J)=SIN(PHIM)
        COSJ(J)=COS(PHIM)
        PHIS=PHIN
        SPHIS=SPHIN
        CPHIS=CPHIN
      END DO
      PHIN=.25*TWOPI
      SPHIN=1.
      CPHIN=0.
      PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
      SINJ(jm0)=SIN(PHIM)
      COSJ(jm0)=COS(PHIM)
!----------------------------------------------------------------------
!Gets cosd and sind
      call declination(gday,gyr,sind,cosd)
!----------------------------------------------------------------------     
      drot=rot2-rot1
!     COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
      SR1=SIN(rot1)
      CR1=COS(rot1)
      SR2=SIN(rot2)
      CR2=COS(rot2)
!     COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
!        NOON) AND THEIR SINES AND COSINES
      DO I=1,im0
        LT1(I)=rot1+LON(I) ! initial local time
        SLT1(I)=SR1*COSIP(I)+CR1*SINIP(I)
        LT2(I)=rot2+LON(I) ! final local time
        SLT2(I)=SR2*COSIP(I)+CR2*SINIP(I)
      END DO
!     
!     CALCULATION FOR POLAR GRID BOXES   
      DO J=1,jm0,jm0-1
        SJSD=SINJ(J)*SIND
        CJCD=COSJ(J)*COSD
        IF (SJSD+CJCD.GT.ZERO1) THEN
          IF (SJSD-CJCD.LT.0.) THEN
!           AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
            DUSK=ACOS(-SJSD/CJCD)
            SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
            DAWN=-DUSK
            SDAWN=-SDUSK
            COSZ(1,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/TWOPI
          ELSE
!           CONSTANT DAYLIGHT NEAR THE POLES
            COSZ(1,J)=SJSD
          END IF
        ELSE
!         CONSTANT NIGHTIME NEAR THE POLES
          COSZ(1,J)=0.
        END IF
      END DO
!     
!     LOOP OVER NON-POLAR LATITUDES
      DO 500 J=2,jm0-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.GT.ZERO1) THEN
      IF (SJSD-CJCD.LT.0.) THEN
!     COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      DUSK=ACOS(-SJSD/CJCD)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      DAWN=-DUSK
      SDAWN=-SDUSK
!     NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
!     LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 400 I=1,im0
!      FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 220
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
  220 IF (DAWN.LT.LT2(I)-ZERO2) GO TO 240
!     CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      GO TO 400
  240 IF (DAWN.GE.LT1(I)) GO TO 300
      IF (DUSK.LT.LT2(I)) GO TO 260
!     CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/drot
      GO TO 400
  260 IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 280
!     DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I)))/drot
      GO TO 400
!     DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
  280 COSZ(I,J)=(SJSD*(LT2(I)-DAWN-TWOPI+DUSK-LT1(I))+CJCD*
     *  (SLT2(I)-SDAWN+SDUSK-SLT1(I)))/drot
      GO TO 400
  300 IF (DUSK.LT.LT2(I)) GO TO 320
!     NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN))/drot
      GO TO 400
!     NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
  320 COSZ(I,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/drot
  400 CONTINUE
      ELSE
!       CONSTANT DAYLIGHT AT THIS LATITUDE
        DO I=1,im0
          COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/drot
        END DO
      END IF
      ELSE
!       CONSTANT NIGHTIME AT THIS LATITUDE
        COSZ(1:im0,J)=0.
      END IF
  500 CONTINUE

      ! Compute fraction of radiation that is diffuse (maybe add to loop above)
      do i=1,im0
         do j=1,jm0
            fdif(i,j) = fdiffuse(cosz(i,j),dble(gday),Rg(i,j))
         enddo
      enddo

      RETURN
      END subroutine COSZT
#endif      
!======================================================================!

      SUBROUTINE calc_Ch(CDN_max,Tsurfair,Tgrnd,Psurfair,
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


#ifdef ECOSYSTEM_SCALE
!      Santarem km 83
!      zgs = 4.5d0
!      zd  = 0.56d0 *0.5d0 
!      z0m = 0.3d0*(0.5d0-zd)
!      z0h = z0m*exp(-2.d0) ! exp(-2) = 0.13533528d0
!      z0q = z0h

!      Barrow
       zgs = 4.2d0
       zd = 0.56d0 * 0.25d0
       z0m = 0.3d0*(0.5d0-zd)
       z0h = z0m*exp(-2.d0)
       z0q = z0h 
       
!      Vaira
!      zgs = 4.5d0
!      zd  = 0.56d0 *0.5d0 
!      z0m = 0.3d0*(0.5d0-zd)
!      z0h = z0m*exp(-2.d0) ! exp(-2) = 0.13533528d0
!      z0q = z0h

!      Tonzi
!      zgs = 23.d0
!      zd  = 0.56d0 * 7.1d0
!      z0m = 0.9d0
!      z0h = z0m/7.d0 
!      z0q = z0h

!!      Hyytiala
!      zgs = 23.3d0
!      zd  = 0.78d0 *13.3d0 
!      z0m = 0.075d0*13.3d0
!      z0h = z0m*exp(-2.d0) ! exp(-2) = 0.13533528d0
!      z0q = z0h

!!     MMSF data
!      zgs = 46.d0
!      zd  = 0.8d0*26.d0 
!      z0m = 2.1d0
!      z0h = z0m*exp(-2.d0) ! exp(-2) = 0.13533528d0
!      z0q = z0h
#else
      ! Conversions currently in PBL.f ( zgrnd=30./10**roughl(i,j) )
      zgs = 10.d0
      ! Conversions currently in PBL.f ( zgrnd=30./10**roughl(i,j) )
      z0m =30.d0 / (10.d0**CDN_max)
           ! zgs/exp((kappa*kappa)/sqrt(0.001d0*CDN_max))
      z0h = z0m*0.13533528d0 !   exp(-2)
      z0q = z0h
      zd = 0.d0 ! set to 0 because reading in max(ztopo,zveg)
#endif

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
      END SUBROUTINE calc_Ch

!======================================================================!


      SUBROUTINE drv_finterp(ip,trp_flag,val0,val1,val2,val3,
     &                       madtt,valtrp)
!======================================================================
! Description: INTERPOLATION OF ATMOSPHERIC FORCING DATA
!
! Revision 1.2  2002/09/20  20:55:32  dirmeyer
! Introduced a new interpolation type "P" which specifies a temporal
! disaggregation of interpolated precipitation to improve partitioning
! of infiltration/runoff when forcing data interval is large compared to
! the typical length of a convective (or total) rain event.
!
! Revision 1.1  2002/08/16  15:55:18  guo
! Initial revision
!====================================================================== 
!** This routine interpolates atmospheric adtt seconds period forcing
!   data to dtt seconds valtrp time-step for one adtt 
!   interval. Interpolation is performed based on the value of 
!   'trp_flag':
!   "L" or "l" = value represents average over interval ending at current t
!   "N" or "n" = value represents average over interval beginning at current t
!   "C" or "c" = value represents average over interval centered on current t
!   "I" or "i" = instantaneous value at current time (linear interpol.)
!   "P" or "p" = PDF disaggregation applied in time (for precip)
!   "0" (zero) = no interpolation, centered on current time
!   Otherwise  = no interpolation, applied beginning at current time
!************************************************************************
!* NOTE!!!!!!!
!*   For "L", "N", and "C" to conserve the mean over the interval after 
!*   interpolation, 'madtt' MUST BE A MULTIPLE OF 2!
!************************************************************************
!   Algorithm for conserving interpolation scheme - courtesy of 
!   Dag Lohmann, NOAA/NCEP/EMC
!   - P. Dirmeyer, November 2001
!===================================================================
 
      IMPLICIT NONE
  
! Input 
      INTEGER,INTENT(IN)           :: ip       ! Grid point in question
      CHARACTER(len=1),INTENT(IN)  :: trp_flag ! Type of interpolation

      real*8,INTENT(IN)            :: val0   ! last forcing value
                                             ! (for cases N & C)
      real*8,INTENT(IN)            :: val1   ! current forcing value(all cases)
      real*8,INTENT(IN)            :: val2   ! next forcing value 
                                             ! (all cases but default)
      real*8,INTENT(IN)            :: val3   ! ueber-next forcing value 
                                             ! (for cases C & L)
      INTEGER,INTENT(IN)           :: madtt  ! Number of valtrp timesteps in
                                             ! a forcing timestep.
! Output
      real*8, DIMENSION(madtt)       :: valtrp ! Interpolated forcing data vector

! Local
      real*8               :: fac0 ! Weight for last value 
      real*8               :: fac1 ! Weight for current value
      real*8               :: fac2 ! Weight for next value
      real*8               :: rmadtt ! real madtt
      real*8               :: denom  ! denominator of scaling factor for 
                                   ! conserving interpolation
      real*8               :: numer  ! numerator of scaling factor for 
                                   ! conserving interpolation
      INTEGER            :: i
      INTEGER            :: j
      INTEGER            :: ip1  ! i + 1
      INTEGER            :: im1  ! i - 1

      real*8             :: rtdist(120) ! Weights for precip disag PDF
      real*8             :: factor      ! Scaling factor for brevity/intensity
      real*8             :: expo        ! Exponent to fit slope of log-log 
                                        ! relationship
      real*8             :: rtsum       ! Ensure PDF weights add to 1.0
      real*8             :: rtsum0      ! Sum from last time interval
      real*8             :: p0, p1      ! Define edges of each time bin
      real*8             :: rintr1      ! Log-log tail target value

      LOGICAL            :: iniflag
      SAVE iniflag,rtdist
      DATA iniflag/.true./

      rmadtt = float(madtt)

!>>> Generate disaggregation PDF - Just do this once
      IF (iniflag) THEN
         write(6,4030)
 4030    format('Initializing precipitation disaggregation PDF:')
         factor = 66.6666  ! For 33% dry time steps in a rainy forcing interval
!         factor = 29.63  ! For 0% dry time steps in a rainy forcing interval
         expo = -3.0      
         rtsum = 0.0
         rintr1 = factor * float(madtt) ** expo

         DO j = 1, madtt
           rtsum0 = rtsum
           p0 = float(j-1)/float(madtt)
           p1 = float(j)/float(madtt)
           rtdist(j) = rintr1*expo/((expo+1)*(p1-p0)*100)*
     &                  ((100*p1/rintr1)**((expo+1)/expo)-
     &                   (100*p0/rintr1)**((expo+1)/expo))
           rtsum = rtsum + rtdist(j)
           IF (rtsum > 1.0) THEN
             IF (rtsum0 < 1.0) THEN
               rtdist(j) = rtdist(j) + 1.0 - rtsum  ! Stick any remaining 
                                                    ! weight in last interval
             ELSE
               rtdist(j) = 0.0
             ENDIF
           ENDIF
           write(6,4050) j,rtdist(j)
 4050      format('Time disag: ',i3,': ',f11.6)
         ENDDO
         IF (rtsum < 1.0) THEN
             rtdist(1) = rtdist(1) + 1.0 - rtsum  ! Ensure integral 
                                                  ! of PDF = 1.0
         ENDIF
      ENDIF
      iniflag = .false.

      IF (trp_flag == 'I' .or. trp_flag == 'i') THEN
      !** Current value valid at midpoint of interval
        DO j = 1, madtt
          fac1     = float(madtt+1-j)/madtt
          fac2     = 1.0-fac1
          valtrp(j)   = val1*fac1+val2*fac2
        ENDDO

      ELSEIF (trp_flag == 'N' .or. trp_flag == 'n') THEN
      !** Current value is average over next interval
        IF (mod(madtt,2) /= 0) STOP 'drv_interp mod(madtt,2) /= 0'
        DO j = 1, madtt
          fac1 = (2.d0*rmadtt-abs(float(2*j-madtt-1)))/(rmadtt*2.d0)
          fac0 = max(1.0-float(j*2+madtt-1)/(rmadtt*2.0),0.0)
          fac2 = max(1.0-float((madtt+1-j)*2+madtt-1)/(rmadtt*2.0),0.d0)
          denom = 0.5d0*(val0+val2)+3.d0*val1
          numer = 4.d0*val1
          IF (denom > EPSILON(denom)) THEN
             valtrp(j) = (val0*fac0+val1*fac1+val2*fac2) * numer / denom
          ELSE
             valtrp(j) = 0.d0
          ENDIF
        ENDDO

      ELSEIF (trp_flag == 'C' .or. trp_flag == 'c') THEN
      !** Current value is average centered on current time
        IF (mod(madtt,2) /= 0) STOP 'drv_interp mod(madtt,2) /= 0'
        DO j = 1, madtt
          fac1 = (2.d0*rmadtt-(float(2*j-1)))/(rmadtt*2.d0)
          fac0 = 0.d0
          fac2 = 1.d0-fac1
          IF (j > madtt/2) THEN
            denom = 0.5d0*(val1+val3)+3.d0*val2
            numer = 4.d0*val2
          ELSE
            denom = 0.5d0*(val0+val2)+3.d0*val1
            numer = 4.d0*val1
          ENDIF
          IF (denom > EPSILON(denom)) THEN
             valtrp(j) = (val0*fac0+val1*fac1+val2*fac2) * numer / denom
          ELSE
             valtrp(j) = 0.d0
          ENDIF
        ENDDO

      ELSEIF (trp_flag == 'L' .or. trp_flag == 'l') THEN
      !** Current value is average for period ending at current time
        IF (mod(madtt,2) /= 0) STOP 'drv_interp mod(madtt,2) /= 0'
        DO j = 1, madtt
          fac1 = (2.d0*rmadtt-abs(float(2*j-madtt-1)))/(rmadtt*2.d0)
          fac0 = max(1.d0-float(j*2+madtt-1)/(rmadtt*2.d0),0.d0)
          fac2 = max(1.d0-float((madtt+1-j)*2+madtt-1)/(rmadtt*2.d0),
     &         0.d0)
          denom = 0.5d0*(val1+val3)+3.d0*val2
          numer = 4.d0*val2
          IF (denom > EPSILON(denom)) THEN
             valtrp(j) = (val1*fac0+val2*fac1+val3*fac2) * numer / denom
          ELSE
             valtrp(j) = 0.d0
          ENDIF
        ENDDO

      ELSEIF (trp_flag == 'P' .or. trp_flag == 'p') THEN
      !** Current value is from a PDF that disagregates in time 
      ! over the forcing interval
        DO j = 1, madtt
          valtrp(j) = MAX(val2*rtdist(j)*float(madtt),0.d0)
        ENDDO

      ELSEIF (trp_flag == '0') THEN
      !** Current value is applied centered on current time without interpol.
        DO j = 1, madtt
          IF (j > madtt/2) THEN
            valtrp(j) = val2
          ELSE
            valtrp(j) = val1
          ENDIF
        ENDDO

      ELSE
      !** Current value is applied over this interval without interpolation
        DO j = 1, madtt
          valtrp(j) = val1
        ENDDO

      ENDIF

      END SUBROUTINE drv_finterp

!======================================================================!
!     Compute month number from start of year
      subroutine get_month(gyr,itime_3hr,n_month)
      implicit none

      integer,intent(in)  :: gyr         ! year
      integer,intent(in)  :: itime_3hr ! 3 hrly timestep # from yr begin
      integer,intent(out) :: n_month     ! month number

      integer, dimension(12) :: mon_end  ! 3-hr timestep # @ end of month
    

      if (LEAPYR(gyr))then
         mon_end = (/ 248, 480, 728, 968,1216,1456,
     &              1704,1952,2192,2440,2680,2928 /)      
      else
         mon_end = (/ 248, 472, 720, 960,1208,1448,
     &              1696,1944,2184,2432,2672,2920 /)
      endif

      if (itime_3hr <= mon_end(1))then
        n_month = 1
      elseif((itime_3hr>mon_end(1)).AND.(itime_3hr<=mon_end(2)))then
        n_month = 2
      elseif((itime_3hr>mon_end(2)).AND.(itime_3hr<=mon_end(3)))then
        n_month = 3
      elseif((itime_3hr>mon_end(3)).AND.(itime_3hr<=mon_end(4)))then
        n_month = 4
      elseif((itime_3hr>mon_end(4)).AND.(itime_3hr<=mon_end(5)))then
        n_month = 5
      elseif((itime_3hr>mon_end(5)).AND.(itime_3hr<=mon_end(6)))then
        n_month = 6
      elseif((itime_3hr>mon_end(6)).AND.(itime_3hr<=mon_end(7)))then
        n_month = 7
      elseif((itime_3hr>mon_end(7)).AND.(itime_3hr<=mon_end(8)))then
        n_month = 8
      elseif((itime_3hr>mon_end(8)).AND.(itime_3hr<=mon_end(9)))then
        n_month = 9
      elseif((itime_3hr>mon_end(9)).AND.
     &                (itime_3hr<=mon_end(10)))then
        n_month = 10
      elseif((itime_3hr>mon_end(10)).AND.
     &                (itime_3hr<=mon_end(11)))then
        n_month = 11
      elseif((itime_3hr>mon_end(11)).AND.
     &                (itime_3hr<=mon_end(12)))then
        n_month = 12
      endif


      end subroutine get_month

!======================================================================!
#ifndef ECOSYSTEM_SCALE
      subroutine get_gswp_data(n_month_sim,iu_vector,
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
      use domain_decomp, only : mype, array_bcast_r4
      integer,intent(in) :: n_month_sim
      real*4, dimension(im0,jm0),intent(out) :: 
     &         data_srheat, data_trheat,  data_ts
     &        ,data_qs,     data_ps,      data_ws
     &        ,data_rainf,  data_rainf_c, data_snowf

      integer,dimension(varnums),intent(in) :: iu_vector
      integer i
      !print *, 'in get_gswp_data'
      if (mype == 0 ) then
        read(iu_vector(1)) data_srheat
        read(iu_vector(2)) data_trheat
        read(iu_vector(3)) data_ts
        read(iu_vector(4)) data_qs
        read(iu_vector(5)) data_ps
        read(iu_vector(6)) data_ws
        read(iu_vector(7)) data_rainf
        read(iu_vector(8)) data_rainf_c
        read(iu_vector(9)) data_snowf
      endif
      !print *, 'after read statement'
      call array_bcast_r4( data_srheat )
      call array_bcast_r4( data_trheat )
      call array_bcast_r4( data_ts )
      call array_bcast_r4( data_qs )
      call array_bcast_r4( data_ps )
      call array_bcast_r4( data_ws )
      call array_bcast_r4( data_rainf )
      call array_bcast_r4( data_rainf_c )
      call array_bcast_r4( data_snowf )
      !print *, 'after array_bcast'

      end subroutine get_gswp_data

!======================================================================!

      subroutine assign_gswp_forcings(gday,gyr,jtime,nday,dtsec,
     i         tbcs_ij             ,
     i         tcanopy_ij          ,
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

      use param, only : sync_param
      use filemanager, only : openunit,closeunit
      implicit none
      integer, intent(in) :: gday      ! gregorian day of year
      integer, intent(in) :: gyr       ! gregorian year
      real*8, intent(in) :: jtime     ! fraction of day
      integer, intent(in) :: nday      ! # of timesteps in a day
      real*8, intent(in) ::  dtsec     ! timestep size [sec]
      real*8,dimension(im0,jm0), intent(in) :: tbcs_ij,tcanopy_ij

!     local variables
      real*8 :: td ! time of day of year
      real*8 :: rot1,rot2 ! beginning and end of timestep for COSZT
      real*8 :: z0,z0_topo,z0_veg ! surface roughness [m]
      real*8 :: Tgrnd ! ground temperature [deg K] 
      real*8 :: RH ! relative humidity
      real*8 :: disp ! zero plane displacement height [m]
      real*8 :: Ch ! heat transfer coefficient
!     GSWP2 forcings
      real*4, dimension(im0,jm0) ::
     &         data_srheat, data_trheat,  data_ts
     &        ,data_qs,     data_ps,      data_ws
     &        ,data_rainf,  data_rainf_c, data_snowf
      real*8, dimension(im0,jm0) :: cosz1 !mean cosine of zenith angle
      real*8, dimension(im0,jm0) :: fdif ! frac of radiation that is diffuse
      real*4, dimension(im0,jm0), save :: CDN_max ! neutral drag coeff. topo
      real*8, dimension(im0,jm0) :: Rg ! global radiation [W/m2]
      real*8, dimension(im0,jm0) ::    ! GISS LSM forcings
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
      integer :: i,j
      character*80 :: title
      character*120 :: neutral_drag_file="CD_coef"
      integer :: iu_nd
      logical :: first_call=.true.

!*****************************************************************************
!     Calculate mean cosine of zenith angle & fraction of diffuse radiation
!      for the current timestep      
      rot1=(2d0*pi*JTIME)/nday ! beginning of timestep
      rot2=rot1+2d0*pi*dtsec/dble(sec_1day) ! end of timestep
      Rg(:,:) = data_srheat(:,:)
      CALL COSZT (gday,gyr,rot1,rot2,Rg,cosz1,fdif)
!*****************************************************************************
      if (first_call) then
        first_call=.false.
!     Opens the drag coefficent file max(topography,vegetation) from GCM
        call sync_param("neutral_drag_file",neutral_drag_file)
        call openunit(neutral_drag_file,iu_nd,.true.,.true.)
        write(0,*) "reading from ", iu_nd
        read(iu_nd) title, CDN_max
        call closeunit(iu_nd)
      endif
!*****************************************************************************
      loop_latitude_j: do j=1,jm0
         loop_longitude_i: do i = 1,im0
            ! data_ts as landmask
            if (data_ts(i,j) > 0.d0) then
              ! 1) Cosine solar zenith angle
               force_cos_zen_angle(i,j)= cosz1(i,j)
              ! 2) Visible radiation [W/m2] 
              !    assume 50% of SW rad.is visible (e.g. Friend 1995 for PAR)
               force_vis_rad(i,j)= 0.5d0*DBLE(data_srheat(i,j))
              ! 3) Direct visible radiation [W/m2]
              force_direct_vis_rad(i,j)=(1-fdif(i,j))*force_vis_rad(i,j)
              ! 4) Total precip.(rain+snow) both convective & large-scale [m/s]
              !    Data [kg/m2/s]
               force_prec_ms(i,j)= ( DBLE(data_rainf(i,j))
     &                              +DBLE(data_snowf(i,j))) / rhow
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
              !11) Surface mixing(humidity)ratio;data:specific humidity [kg/kg]
               force_qs(i,j)= DBLE(data_qs(i,j))/
     &                             (1.d0-DBLE(data_qs(i,j)))
               if(force_qs(i,j)<0.0001d0) force_qs(i,j)=0.0001d0 
              !12) Surface pressure [mb]; data_ps [Pa]
               force_ps(i,j)= DBLE(data_ps(i,j))/100.d0 
              !13) Surface air density [kg/m3]
               force_rhosrf(i,j)= data_ps(i,j)/(rgas*data_ts(i,j))
              !14) Turbulent transfer, Ch (see after 21)
              !15) Amount of water in the 1st atm layer; set to infinity
               force_qm1(i,j)= 1e30
              !16) Wind speed [m/s]; data_ws [m/s]
               force_ws(i,j)= DBLE(data_ws(i,j))
               if(force_ws(i,j)<0.01d0) force_ws(i,j)=0.01d0 
              !17) Atmospheric CO2 concentration at the land surface [mol/m3]
                !350 ppm = 0.0127609 mol/m3 at STP
               force_Ca(i,j)= CO2ppm*(1.0D-06)*data_ps(i,j)
     &                           /gasc/data_ts(i,j)
              !18) - 21) Planetary boundary layer variables 
               force_pbl_args_ws0 (i,j)    = 0.d0
               force_pbl_args_tprime (i,j) = 0.d0
               force_pbl_args_qprime (i,j) = 0.d0

              !14) Turbulent transfer / humidity transfer coefficient [-]
              !    humidity transfer coef = heat transfer coef
               !   Surface temperature, Tgrnd, either canopy or ground temp.  
               if (tcanopy_ij(i,j)>-1d30)then
                  Tgrnd = tcanopy_ij(i,j) + tf ! set to canopy temp; C to K
               elseif (tbcs_ij(i,j)>-1d30)then
                  Tgrnd = tbcs_ij(i,j) + tf    ! set to ground temp; C to K
               else ! first timestep
                  Tgrnd = force_ts(i,j)
               endif

               call calc_Ch(DBLE(CDN_max(i,j)), force_ts(i,j),Tgrnd,
     &                      force_ps(i,j),force_ws(i,j), Ch)
               force_cdh(i,j)= Ch
!               write(3000,20) i, j, Ch
! 20            format(i10,i10,f15.10)
            else
               force_cos_zen_angle(i,j)= UNDEF
               force_vis_rad(i,j)= UNDEF
               force_direct_vis_rad(i,j)= UNDEF
               force_prec_ms(i,j)= UNDEF
               force_eprec_w(i,j)= UNDEF
               force_sprec_ms(i,j)= UNDEF
               force_seprec_w(i,j)= UNDEF
               force_srheat(i,j)= UNDEF
               force_trheat(i,j)= UNDEF
               force_ts(i,j)= UNDEF
               force_qs(i,j)= UNDEF
               force_ps(i,j)= UNDEF
               force_rhosrf(i,j)= UNDEF
               force_cdh(i,j)= UNDEF
               force_qm1(i,j)= UNDEF
               force_ws(i,j)= UNDEF
               force_Ca(i,j)= UNDEF
               force_pbl_args_ws0 (i,j)    = UNDEF
               force_pbl_args_tprime (i,j) = UNDEF
               force_pbl_args_qprime (i,j) = UNDEF
            endif
         enddo loop_longitude_i
      enddo loop_latitude_j


      end subroutine assign_gswp_forcings

!-------------------------------------------------------------------

      subroutine open_gswp_files(step_num,gyr,itime_3hr,n_month,
     &        n_month_sim,iu_vector,init_flag)
      use param, only : sync_param
      use filemanager, only : openunit

      implicit none

      integer, intent(in) :: gyr         ! gregorian year
      integer,intent(in)  :: itime_3hr ! 3 hrly timestep # from yr begin
      integer,intent(in)  :: n_month     ! month number
      integer,intent(in)  :: n_month_sim ! month # from beginning of sim.
      real*8,intent(in)   :: step_num    ! fraction of 3 hr timestep

      integer, parameter :: num_outfiles = num_months*varnums

      integer :: i,init_flag
      character*120 :: infile
      character*80 :: monthyearBIN(num_months)
      character*80 :: lsm_vars(varnums),basefold_out
      character*80 :: lsm_fold(varnums)

      integer,dimension(varnums),intent(inout) :: iu_vector
      integer,dimension(12) :: mon_begin ! 1st 3hr timestep of month

      DATA lsm_fold/ 'force_srheat/'
     &        ,'force_trheat/', 'force_ts/', 'force_qs/'
     &        ,'force_ps/','force_ws/'
     &        ,'force_rainf/', 'force_rainf_c/', 'force_snowf/'/

      DATA lsm_vars/ 'force_srheat'
     &        ,'force_trheat', 'force_ts', 'force_qs'
     &        ,'force_ps','force_ws'
     &        ,'force_rainf', 'force_rainf_c', 'force_snowf'/

      DATA monthyearBIN/       '198207.bi','198208.bi',
     & '198209.bi','198210.bi','198211.bi','198212.bi',
     & '198301.bi','198302.bi','198303.bi','198304.bi',
     & '198305.bi','198306.bi','198307.bi','198308.bi',
     & '198309.bi','198310.bi','198311.bi','198312.bi',
     & '198401.bi','198402.bi','198403.bi','198404.bi',
     & '198405.bi','198406.bi','198407.bi','198408.bi',
     & '198409.bi','198410.bi','198411.bi','198412.bi',
     & '198501.bi','198502.bi','198503.bi','198504.bi',
     & '198505.bi','198506.bi','198507.bi','198508.bi',
     & '198509.bi','198510.bi','198511.bi','198512.bi',
     & '198601.bi','198602.bi','198603.bi','198604.bi',
     & '198605.bi','198606.bi','198607.bi','198608.bi',
     & '198609.bi','198610.bi','198611.bi','198612.bi',
     & '198701.bi','198702.bi','198703.bi','198704.bi',
     & '198705.bi','198706.bi','198707.bi','198708.bi',
     & '198709.bi','198710.bi','198711.bi','198712.bi',
     & '198801.bi','198802.bi','198803.bi','198804.bi',
     & '198805.bi','198806.bi','198807.bi','198808.bi',
     & '198809.bi','198810.bi','198811.bi','198812.bi',
     & '198901.bi','198902.bi','198903.bi','198904.bi',
     & '198905.bi','198906.bi','198907.bi','198908.bi',
     & '198909.bi','198910.bi','198911.bi','198912.bi',
     & '199001.bi','199002.bi','199003.bi','199004.bi',
     & '199005.bi','199006.bi','199007.bi','199008.bi',
     & '199009.bi','199010.bi','199011.bi','199012.bi',
     & '199101.bi','199102.bi','199103.bi','199104.bi',
     & '199105.bi','199106.bi','199107.bi','199108.bi',
     & '199109.bi','199110.bi','199111.bi','199112.bi',
     & '199201.bi','199202.bi','199203.bi','199204.bi',
     & '199205.bi','199206.bi','199207.bi','199208.bi',
     & '199209.bi','199210.bi','199211.bi','199212.bi',
     & '199301.bi','199302.bi','199303.bi','199304.bi',
     & '199305.bi','199306.bi','199307.bi','199308.bi',
     & '199309.bi','199310.bi','199311.bi','199312.bi',
     & '199401.bi','199402.bi','199403.bi','199404.bi',
     & '199405.bi','199406.bi','199407.bi','199408.bi',
     & '199409.bi','199410.bi','199411.bi','199412.bi',
     & '199501.bi','199502.bi','199503.bi','199504.bi',
     & '199505.bi','199506.bi','199507.bi','199508.bi',
     & '199509.bi','199510.bi','199511.bi','199512.bi' /

      basefold_out = '/Users/nkiang/NancyResearch/GISS/Models/Ent/'//
     &     'Datasets/GSWP/lsm_gswp_4x5'
!      basefold_out = '/discover/nobackup/mpuma/lsm_gswp'
!      basefold_out = '/discover/nobackup/mpuma/lsm_gswp_2x2_5'
!      basefold_out= '/n/Moorcroft_Lab/Users/kim/ent-data/GSWP2'

      call sync_param( "basefold_out", basefold_out )

      if (LEAPYR(gyr))then
         mon_begin = (/   1, 249, 481, 729, 969,1217,
     &                1457,1705,1953,2193,2441,2681 /)      
      else
         mon_begin = (/  1,  249, 473, 721, 961,1209,
     &                1449,1697,1945,2185,2433,2673 /)
      endif
      !print *, 'n_step_3hr = ', itime_3hr 
      !print *, 'mon_begin=', mon_begin(n_month)

!     Check if we are at the beginning of a month      
      if((itime_3hr == mon_begin(n_month)).and.(step_num <1e-6))then

!        Close previously opened files
         if (init_flag==0)then
            call close_outfiles(iu_vector)
         endif

!        open forcings file
         do i = 1,varnums
            infile=basefold_out(1:len_trim(basefold_out))//'/'//
     &         lsm_fold(i)(1:len_trim(lsm_fold(i)))//
     &         lsm_vars(i)(1:len_trim(lsm_vars(i)))//
     &         monthyearBIN(n_month_sim)
     &         (1:len_trim(monthyearBIN(n_month_sim)))
            print *, infile
            print *,  n_month_sim, monthyearBIN(n_month_sim)
     &         (1:len_trim(monthyearBIN(n_month_sim)))
            call openunit(infile,iu_vector(i),.true.,.true.)
!            open(iu_vector(i),FILE=infile, STATUS='UNKNOWN',
!     &              FORM='UNFORMATTED')
         enddo
         !print *, 'iu_vector after openunit', iu_vector

      endif
      end subroutine open_gswp_files

#endif
!======================================================================

      subroutine close_outfiles(iu_vector)
      use filemanager, only : closeunit

      integer, parameter :: varnums = 9 ! number of GSWP2 variables
      integer,dimension(varnums),intent(in) :: iu_vector
      integer i

      do i = 1,varnums
         call closeunit(iu_vector(i))
      enddo
      
      end subroutine close_outfiles
!======================================================================!
#ifdef ECOSYSTEM_SCALE

      subroutine get_single_data(iu_vector,
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

      real*4,dimension(im0,jm0),intent(out) :: data_srheat, data_trheat,
     &                                         data_ts
     &        ,data_qs,     data_ps,      data_ws
     &        ,data_rainf,  data_rainf_c, data_snowf

      integer,intent(in) :: iu_vector

      read(iu_vector,*) data_srheat,data_trheat,data_ts, data_qs,
     &     data_ps,data_ws,data_rainf,data_rainf_c,data_snowf
!      print *, 'in get_single_data'
!      print *, data_srheat,data_trheat,data_ts, data_qs,
!     &     data_ps,data_ws,data_rainf,data_rainf_c,data_snowf


      end subroutine get_single_data

!======================================================================!

      subroutine assign_single_forcings(
     i         i_site, j_site, latd,
     i         gday,gyr,jtime,nday,dtsec,
     i         tbcs_ij             ,
     i         tcanopy_ij          ,
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

      use param, only : sync_param
      use filemanager, only : openunit,closeunit
      implicit none
      integer, intent(in) :: gday      ! gregorian day of year
      integer, intent(in) :: gyr       ! gregorian year
      real*8, intent(in) :: jtime     ! fraction of day
      integer, intent(in) :: nday      ! # of timesteps in a day
      real*8, intent(in) ::  dtsec     ! timestep size [sec]
      real*8, intent(in) :: tbcs_ij,tcanopy_ij
      integer, intent(in) :: i_site, j_site
      real*8,  intent(in) :: latd

!     local variables
      real*8 :: td ! time of day of year
      real*8 :: cos_zen_angle ! cos( solar zenith angle)
      real*8 :: rot1,rot2 ! beginning and end of timestep for COSZT
      real*8 :: z0,z0_topo,z0_veg ! surface roughness [m]
      real*8 :: Tgrnd ! ground temperature [deg K] 
      real*8 :: RH ! relative humidity
      real*8 :: disp ! zero plane displacement height [m]
      real*8 :: Ch ! heat transfer coefficient
!     GSWP2 forcings
      real*4, dimension(im0,jm0) ::
     &         data_srheat, data_trheat,  data_ts
     &        ,data_qs,     data_ps,      data_ws
     &        ,data_rainf,  data_rainf_c, data_snowf
      real*8, dimension(im0,jm0) :: fdif ! frac of radiation that is diffuse
      real*4, dimension(im0,jm0), save :: CDN_max ! neutral drag coeff. topo
      real*8 :: Rg ! global radiation [W/m2]
      real*8, dimension(im0,jm0) ::    ! GISS LSM forcings
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
      integer :: i,j
      character*80 :: title
      character*120 :: neutral_drag_file="CD_coef"
      integer :: iu_nd
      logical :: first_call=.true.

!*****************************************************************************
!     Calculate mean cosine of zenith angle & fraction of diffuse radiation
!      for the current timestep      
      rot1=(2d0*pi*JTIME)/nday ! beginning of timestep
      rot2=rot1+2d0*pi*dtsec/dble(sec_1day) ! end of timestep
      td = DBLE(gday)+JTIME/nday
      call calc_solarzen(gday,td,latd,cos_zen_angle)
      Rg = data_srheat(1,1)
      fdif = fdiffuse(cos_zen_angle,td,Rg)
!*****************************************************************************
      if (first_call) then
        first_call=.false.
!     Opens the drag coefficent file max(topography,vegetation) from GCM
        call sync_param("neutral_drag_file",neutral_drag_file)
        call openunit(neutral_drag_file,iu_nd,.true.,.true.)
        write(0,*) "reading from ", iu_nd
        read(iu_nd) title, CDN_temp
        call closeunit(iu_nd)
        CDN_max=CDN_temp(i_site,j_site)
      endif
!*****************************************************************************
      loop_latitude_j: do j=1,jm0
         loop_longitude_i: do i = 1,im0
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
     &                              +DBLE(data_snowf(i,j)))/ rho_water
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
               if (tcanopy_ij>-1d30)then
                  Tgrnd = tcanopy_ij + tf ! set to canopy temp; C to K
               elseif(tbcs_ij>-1d30)then
                  Tgrnd = tbcs_ij + tf    ! set to ground temp; C to K
               else ! first timestep
                  Tgrnd = force_ts(i,j)     
               endif

               call calc_Ch(DBLE(CDN_max(i,j)), force_ts(i,j),Tgrnd,
     &                      force_ps(i,j),force_ws(i,j), Ch)
               force_cdh(i,j)= Ch
!               write(3000,20) i, j, Ch
! 20            format(i10,i10,f15.10)
            else
               force_cos_zen_angle(i,j)= UNDEF
               force_vis_rad(i,j)= UNDEF
               force_direct_vis_rad(i,j)= UNDEF
               force_prec_ms(i,j)= UNDEF
               force_eprec_w(i,j)= UNDEF
               force_sprec_ms(i,j)= UNDEF
               force_seprec_w(i,j)= UNDEF
               force_srheat(i,j)= UNDEF
               force_trheat(i,j)= UNDEF
               force_ts(i,j)= UNDEF
               force_qs(i,j)= UNDEF
               force_ps(i,j)= UNDEF
               force_rhosrf(i,j)= UNDEF
               force_cdh(i,j)= UNDEF
               force_qm1(i,j)= UNDEF
               force_ws(i,j)= UNDEF
               force_Ca(i,j)= UNDEF
               force_pbl_args_ws0 (i,j)    = UNDEF
               force_pbl_args_tprime (i,j) = UNDEF
               force_pbl_args_qprime (i,j) = UNDEF
            endif
         enddo loop_longitude_i
      enddo loop_latitude_j


      end subroutine assign_single_forcings

!-------------------------------------------------------------------

      subroutine open_single_files(iu_vector)
      use param, only : sync_param
      use filemanager, only : openunit

      implicit none

      character*80 :: lsm_vars(varnums),basefold_out
      character*80 :: filename
      integer,intent(inout) :: iu_vector

      filename ='site_forcings'

!     open forcings file
      call openunit(trim(filename),iu_vector,.false.,.true.)
      open(iu_vector,FILE=filename, STATUS='OLD')

      end subroutine open_single_files

#endif

!======================================================================

      subroutine init_forcings(i_site,j_site,latd,
     &      gday,gyr,tyr_sec,dtsec,iu_vector,tbcs_ij,tcanopy_ij )
      use filemanager, only : closeunit

      implicit none

      integer, intent(in) :: gday      ! gregorian day of year
      integer, intent(in) :: gyr       ! gregorian year
      integer, intent(in) :: tyr_sec   ! time from start of year [sec]
      real*8, intent(in) :: dtsec     ! timestep size [sec]
      integer, intent(in) :: i_site, j_site
      real*8,  intent(in) :: latd

      real*4, dimension(im0,jm0) :: ! GSWP forcings data
     &         data_srheat, data_trheat,  data_ts
     &        ,data_qs,     data_ps,      data_ws
     &        ,data_rainf,  data_rainf_c, data_snowf
#ifdef ECOSYSTEM_SCALE
      integer,intent(inout) :: iu_vector
      real*8, intent(in) :: tbcs_ij, tcanopy_ij
#else
      integer,dimension(varnums),intent(inout) :: iu_vector
      real*8,dimension(im0,jm0), intent(in) :: tbcs_ij, tcanopy_ij
#endif
      integer :: n_month       ! month # from beginning of yr
      integer :: n_month_sim   ! month # from beginning of simulation
      integer :: itime_3hr   ! timestep # (3hr) from start of year
      integer :: n_year        ! year # from beginning of simulation
      integer :: n_month_prev,ip,i,j

      integer :: nday          ! # of timesteps in a day
      real*8 :: itime,jtime   
      real*8  :: step_num_3hr      ! fraction of 3 hr timestep

!     time variables for solar zenith angle calculations
      nday = dble(sec_1day)/dtsec
      itime = dble(tyr_sec)/dtsec
      jtime = mod(itime,dble(nday))

#ifdef ECOSYSTEM_SCALE
      call open_single_files(iu_vector)

      call get_single_data(iu_vector    ,
     o         data_srheat (1:im0,1:jm0)          ,
     o         data_trheat (1:im0,1:jm0)          ,
     o         data_ts (1:im0,1:jm0)              ,
     o         data_qs (1:im0,1:jm0)              ,
     o         data_ps (1:im0,1:jm0)              ,
     o         data_ws (1:im0,1:jm0)              ,
     o         data_rainf( 1:im0,1:jm0)           ,
     o         data_rainf_c (1:im0,1:jm0)         ,
     o         data_snowf (1:im0,1:jm0) )

      call closeunit(iu_vector)

!     Obtain forcings for initial timestep
      call assign_single_forcings(i_site,j_site,latd,
     i         gday,gyr,jtime,nday,dtsec         ,
     i         tbcs_ij                           ,
     i         tcanopy_ij                        ,
     i         data_srheat (1:im0,1:jm0)         ,
     i         data_trheat (1:im0,1:jm0)         ,
     i         data_ts (1:im0,1:jm0)             ,
     i         data_qs (1:im0,1:jm0)             ,
     i         data_ps (1:im0,1:jm0)             ,
     i         data_ws (1:im0,1:jm0)             ,
     i         data_rainf (1:im0,1:jm0)          ,
     i         data_rainf_c (1:im0,1:jm0)        ,
     i         data_snowf (1:im0,1:jm0)          ,
     o         force_Ca_1 (1:im0,1:jm0)             ,
     o         force_cos_zen_angle_1 (1:im0,1:jm0)  ,
     o         force_vis_rad_1 (1:im0,1:jm0)        ,
     o         force_direct_vis_rad_1 (1:im0,1:jm0) ,
     o         force_prec_ms_1 (1:im0,1:jm0)        ,
     o         force_eprec_w_1 (1:im0,1:jm0)        ,
     o         force_sprec_ms_1 (1:im0,1:jm0)       ,
     o         force_seprec_w_1 (1:im0,1:jm0)       ,
     o         force_srheat_1 (1:im0,1:jm0)         ,
     o         force_trheat_1 (1:im0,1:jm0)         ,
     o         force_ts_1 (1:im0,1:jm0)             ,
     o         force_qs_1 (1:im0,1:jm0)             ,
     o         force_ps_1 (1:im0,1:jm0)             ,
     o         force_rhosrf_1 (1:im0,1:jm0)         ,
     o         force_cdh_1 (1:im0,1:jm0)            ,
     o         force_qm1_1 (1:im0,1:jm0)            ,
     o         force_ws_1 (1:im0,1:jm0)             ,
     o         force_pbl_args_ws0_1 (1:im0,1:jm0)   ,
     o         force_pbl_args_tprime_1 (1:im0,1:jm0),
     o         force_pbl_args_qprime_1 (1:im0,1:jm0) )
      print *, 'after assign forcings'
#else
!     3-hourly counters for opening GSWP data (t0 = hour1 not hour0)
      itime_3hr = tyr_sec/sec_3_hr + 1
      step_num_3hr = (dble(tyr_sec+sec_3_hr))/dble(sec_3_hr) - itime_3hr

!     Get month number from start of year
      call get_month(gyr,itime_3hr,n_month)

!     Compute month number counted from start of simulation 
      n_year = gyr - 1982
      if (n_year == 0) then
        n_month_prev = -6
      elseif(n_year == 1) then
        n_month_prev = 6
      else
        n_month_prev = 6 + 12*(n_year-1)
      endif 

      n_month_sim = n_month + n_month_prev
      print *,'gyr,n_year,n_month_prev,n_month,n_month_sim',
     &     gyr,n_year,n_month_prev,n_month,n_month_sim
      call open_gswp_files(step_num_3hr,gyr,itime_3hr
     &           ,n_month,n_month_sim,iu_vector,1)

      call get_gswp_data(n_month_sim,iu_vector    ,
     o         data_srheat (1:im0,1:jm0)          ,
     o         data_trheat (1:im0,1:jm0)          ,
     o         data_ts (1:im0,1:jm0)              ,
     o         data_qs (1:im0,1:jm0)              ,
     o         data_ps (1:im0,1:jm0)              ,
     o         data_ws (1:im0,1:jm0)              ,
     o         data_rainf( 1:im0,1:jm0)           ,
     o         data_rainf_c (1:im0,1:jm0)         ,
     o         data_snowf (1:im0,1:jm0) )

!     Obtain forcings for initial timestep -- will become current timestep
      call assign_gswp_forcings(gday,gyr,jtime,nday,dtsec,
     i         tbcs_ij (1:im0,1:jm0)             ,
     i         tcanopy_ij(1:im0,1:jm0)           ,
     i         data_srheat (1:im0,1:jm0)         ,
     i         data_trheat (1:im0,1:jm0)         ,
     i         data_ts (1:im0,1:jm0)             ,
     i         data_qs (1:im0,1:jm0)             ,
     i         data_ps (1:im0,1:jm0)             ,
     i         data_ws (1:im0,1:jm0)             ,
     i         data_rainf (1:im0,1:jm0)          ,
     i         data_rainf_c (1:im0,1:jm0)        ,
     i         data_snowf (1:im0,1:jm0)          ,
     o         force_Ca_2 (1:im0,1:jm0)             ,
     o         force_cos_zen_angle_2 (1:im0,1:jm0)  ,
     o         force_vis_rad_2 (1:im0,1:jm0)        ,
     o         force_direct_vis_rad_2 (1:im0,1:jm0) ,
     o         force_prec_ms_2 (1:im0,1:jm0)        ,
     o         force_eprec_w_2 (1:im0,1:jm0)        ,
     o         force_sprec_ms_2 (1:im0,1:jm0)       ,
     o         force_seprec_w_2 (1:im0,1:jm0)       ,
     o         force_srheat_2 (1:im0,1:jm0)         ,
     o         force_trheat_2 (1:im0,1:jm0)         ,
     o         force_ts_2 (1:im0,1:jm0)             ,
     o         force_qs_2 (1:im0,1:jm0)             ,
     o         force_ps_2 (1:im0,1:jm0)             ,
     o         force_rhosrf_2 (1:im0,1:jm0)         ,
     o         force_cdh_2 (1:im0,1:jm0)            ,
     o         force_qm1_2 (1:im0,1:jm0)            ,
     o         force_ws_2 (1:im0,1:jm0)             ,
     o         force_pbl_args_ws0_2 (1:im0,1:jm0)   ,
     o         force_pbl_args_tprime_2 (1:im0,1:jm0),
     o         force_pbl_args_qprime_2 (1:im0,1:jm0) )

!        Will become previous timstep (i.e. setting previous = current)
         force_Ca_1(1:im0,1:jm0)=force_Ca_2(1:im0,1:jm0)
         force_cos_zen_angle_1(1:im0,1:jm0)= 
     &                        force_cos_zen_angle_2(1:im0,1:jm0)
         force_vis_rad_1(1:im0,1:jm0) = force_vis_rad_2(1:im0,1:jm0)
         force_direct_vis_rad_1(1:im0,1:jm0) = 
     &                        force_direct_vis_rad_2(1:im0,1:jm0)
         force_prec_ms_1(1:im0,1:jm0) = force_prec_ms_2(1:im0,1:jm0)
         force_eprec_w_1(1:im0,1:jm0) = force_eprec_w_2(1:im0,1:jm0)
         force_sprec_ms_1(1:im0,1:jm0) = force_sprec_ms_2(1:im0,1:jm0)
         force_seprec_w_1(1:im0,1:jm0) = force_seprec_w_2(1:im0,1:jm0)
         force_srheat_1(1:im0,1:jm0) = force_srheat_2(1:im0,1:jm0)
         force_trheat_1(1:im0,1:jm0) = force_trheat_2(1:im0,1:jm0)
         force_ts_1(1:im0,1:jm0) = force_ts_2(1:im0,1:jm0)
         force_qs_1(1:im0,1:jm0) = force_qs_2(1:im0,1:jm0)
         force_ps_1(1:im0,1:jm0) = force_ps_2(1:im0,1:jm0)
         force_rhosrf_1(1:im0,1:jm0) = force_rhosrf_2(1:im0,1:jm0)
         force_cdh_1(1:im0,1:jm0) = force_cdh_2(1:im0,1:jm0)
         force_qm1_1(1:im0,1:jm0) = force_qm1_2(1:im0,1:jm0)
         force_ws_1(1:im0,1:jm0) = force_ws_2(1:im0,1:jm0)
         force_pbl_args_ws0_1(1:im0,1:jm0) =
     &                       force_pbl_args_ws0_2(1:im0,1:jm0)
         force_pbl_args_tprime_1(1:im0,1:jm0) =
     &                       force_pbl_args_tprime_2(1:im0,1:jm0) 
         force_pbl_args_qprime_1(1:im0,1:jm0) = 
     &                       force_pbl_args_qprime_2(1:im0,1:jm0)

#endif
      print *, '********END OF INITIAL ASSIGNMENT**************'
      end subroutine init_forcings

!======================================================================

      subroutine get_gswp_forcings(
     i         i_site, j_site, latd,
     i         gday                  ,
     i         gyr                   ,
     i         tyr_sec               ,
     i         dtsec                 ,
     o         iu_vector             ,
     i         tbcs_ij               ,
     i         tcanopy_ij            ,
     i         end_of_input_flag     ,
     o         force_Ca              ,
     o         force_cos_zen_angle   ,
     o         force_vis_rad         ,
     o         force_direct_vis_rad  ,
     o         force_prec_ms         ,
     o         force_eprec_w         ,
     o         force_sprec_ms        ,
     o         force_seprec_w        ,
     o         force_srheat          ,
     o         force_trheat          ,
     o         force_ts              ,
     o         force_qs              ,
     o         force_ps              ,
     o         force_rhosrf          ,
     o         force_cdh             ,
     o         force_qm1             ,
     o         force_ws              ,
     o         force_pbl_args_ws0    ,
     o         force_pbl_args_tprime ,
     o         force_pbl_args_qprime  )

      implicit none

      integer, intent(in) :: gday      ! gregorian day of year
      integer, intent(in) :: gyr       ! gregorian year
      integer, intent(in) :: tyr_sec   ! time from start of year [sec]
      real*8,  intent(in) :: dtsec     ! timestep size [sec]
      integer, intent(in) :: i_site, j_site
      real*8,  intent(in) :: latd
      logical, intent(in) :: end_of_input_flag

      real*4, dimension(im0,jm0) :: ! GSWP forcings data
     &         data_srheat, data_trheat,  data_ts
     &        ,data_qs,     data_ps,      data_ws
     &        ,data_rainf,  data_rainf_c, data_snowf

      real*8, dimension(im0,jm0) :: ! Forcings returned for current time
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

      real*8,dimension(im0,jm0) :: ! Forcings at time t0
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
#ifdef ECOSYSTEM_SCALE
      integer,intent(inout) :: iu_vector
      real*8, intent(in) :: tbcs_ij, tcanopy_ij
#else
      integer,dimension(varnums),intent(inout) :: iu_vector
      real*8,dimension(im0,jm0), intent(in) :: tbcs_ij, tcanopy_ij
#endif
      real*8, dimension(n_LSM_GSWP) :: valtrp

      integer :: n_month       ! month # from beginning of yr
      integer :: n_month_sim   ! month # from beginning of simulation
      integer :: n_year        ! year # from beginning of simulation
      integer :: n_month_prev,ip,i,j

      integer :: itime_3hr     ! timestep # (3hr) from start of year
      integer :: i_sub3hr      ! counter for sub-3hrly timestep 
      real*8  :: step_num_3hr      ! fraction of 3 hr timestep
      integer :: n3hr          ! # of timesteps in 3hr:( 3hr)/dtsec

      integer :: nday          ! # of timesteps in a day:(1 day)/dtsec
      integer :: itime         ! current time in ITUs (1 ITU = dtsec)
      real*8 :: jtime         ! sub-daily timestep
      logical :: first_read=.true.

!KIM-to update cosz every timestep
      real*8 :: rot1,rot2
      real*8, dimension(im0,jm0) :: cosz1
      real*8, dimension(im0,jm0) :: fdif 
      real*8, dimension(im0,jm0) :: Rg

!     time variables for solar zenith angle calculations
      nday = dble(sec_1day)/dtsec
      itime = dble(tyr_sec)/dtsec
      jtime = mod(itime,nday)

#ifdef ECOSYSTEM_SCALE

         if(first_read)then
            first_read=.false.
            call open_single_files(iu_vector)
         endif

         call get_single_data(iu_vector,
     o         data_srheat (1:im0,1:jm0)          ,
     o         data_trheat (1:im0,1:jm0)          ,
     o         data_ts (1:im0,1:jm0)              ,
     o         data_qs (1:im0,1:jm0)              ,
     o         data_ps (1:im0,1:jm0)              ,
     o         data_ws (1:im0,1:jm0)              ,
     o         data_rainf( 1:im0,1:jm0)           ,
     o         data_rainf_c (1:im0,1:jm0)         ,
     o         data_snowf (1:im0,1:jm0) )
         !print *, 'After get_single_data'

!        Obtain forcings for current timestep
         call assign_single_forcings(
     i         i_site,j_site,latd                ,
     i         gday,gyr,jtime,nday,dtsec         ,
     i         tbcs_ij                           ,
     i         tcanopy_ij                        ,
     i         data_srheat (1:im0,1:jm0)         ,
     i         data_trheat (1:im0,1:jm0)         ,
     i         data_ts (1:im0,1:jm0)             ,
     i         data_qs (1:im0,1:jm0)             ,
     i         data_ps (1:im0,1:jm0)             ,
     i         data_ws (1:im0,1:jm0)             ,
     i         data_rainf (1:im0,1:jm0)          ,
     i         data_rainf_c (1:im0,1:jm0)        ,
     i         data_snowf (1:im0,1:jm0)          ,
     o         force_Ca (1:im0,1:jm0)             ,
     o         force_cos_zen_angle (1:im0,1:jm0)  ,
     o         force_vis_rad (1:im0,1:jm0)        ,
     o         force_direct_vis_rad (1:im0,1:jm0) ,
     o         force_prec_ms (1:im0,1:jm0)        ,
     o         force_eprec_w (1:im0,1:jm0)        ,
     o         force_sprec_ms (1:im0,1:jm0)       ,
     o         force_seprec_w (1:im0,1:jm0)       ,
     o         force_srheat (1:im0,1:jm0)         ,
     o         force_trheat (1:im0,1:jm0)         ,
     o         force_ts (1:im0,1:jm0)             ,
     o         force_qs (1:im0,1:jm0)             ,
     o         force_ps (1:im0,1:jm0)             ,
     o         force_rhosrf (1:im0,1:jm0)         ,
     o         force_cdh (1:im0,1:jm0)            ,
     o         force_qm1 (1:im0,1:jm0)            ,
     o         force_ws (1:im0,1:jm0)             ,
     o         force_pbl_args_ws0 (1:im0,1:jm0)   ,
     o         force_pbl_args_tprime (1:im0,1:jm0),
     o         force_pbl_args_qprime (1:im0,1:jm0) )

#else
!     3-hourly counters for opening GSWP data (t0 = hour1 not hour0)
      itime_3hr = tyr_sec/sec_3_hr + 1
      step_num_3hr = (dble(tyr_sec+sec_3_hr))/dble(sec_3_hr) - itime_3hr
      i_sub3hr = step_num_3hr *(n_LSM_GSWP+1e-06) + 1

!     Get month number from start of year
      call get_month(gyr,itime_3hr,n_month)

!      print *, 'time from start of yr, year =', tyr_sec, gyr
!      print *, 'current timestep number (itime): ',itime 
!      print *, 'subdaily timestep number(jtime): ',jtime 
!      print *, 'number of timesteps in a day (nday):  ',nday  
!      print *, '3hr fraction, step_num_3hr = ', step_num_3hr
!      print *, 'sub-3hrly timestep = ', i_sub3hr
!      print *, 'MONTH #, n_month = ', n_month

!     Compute month number counted from start of simulation 
      n_year = gyr - 1982
      if (n_year == 0) then
        n_month_prev = -6
      elseif(n_year == 1) then
        n_month_prev = 6
      else
        n_month_prev = 6 + 12*(n_year-1)
      endif 

      n_month_sim = n_month + n_month_prev

!     Check to see if current time is a new 3hrly timestep
      if (step_num_3hr <1e-6) then

!        Assign forcings to time t0 (last forcing value)         
         force_Ca_0(1:im0,1:jm0)=force_Ca_1(1:im0,1:jm0)
         force_cos_zen_angle_0(1:im0,1:jm0)= 
     &                        force_cos_zen_angle_1(1:im0,1:jm0)
         force_vis_rad_0(1:im0,1:jm0) = force_vis_rad_1(1:im0,1:jm0)
         force_direct_vis_rad_0(1:im0,1:jm0) = 
     &                        force_direct_vis_rad_1(1:im0,1:jm0)
         force_prec_ms_0(1:im0,1:jm0) = force_prec_ms_1(1:im0,1:jm0)
         force_eprec_w_0(1:im0,1:jm0) = force_eprec_w_1(1:im0,1:jm0)
         force_sprec_ms_0(1:im0,1:jm0) = force_sprec_ms_1(1:im0,1:jm0)
         force_seprec_w_0(1:im0,1:jm0) = force_seprec_w_1(1:im0,1:jm0)
         force_srheat_0(1:im0,1:jm0) = force_srheat_1(1:im0,1:jm0)
         force_trheat_0(1:im0,1:jm0) = force_trheat_1(1:im0,1:jm0)
         force_ts_0(1:im0,1:jm0) = force_ts_1(1:im0,1:jm0)
         force_qs_0(1:im0,1:jm0) = force_qs_1(1:im0,1:jm0)
         force_ps_0(1:im0,1:jm0) = force_ps_1(1:im0,1:jm0)
         force_rhosrf_0(1:im0,1:jm0) = force_rhosrf_1(1:im0,1:jm0)
         force_cdh_0(1:im0,1:jm0) = force_cdh_1(1:im0,1:jm0)
         force_qm1_0(1:im0,1:jm0) = force_qm1_1(1:im0,1:jm0)
         force_ws_0(1:im0,1:jm0) = force_ws_1(1:im0,1:jm0)
         force_pbl_args_ws0_0(1:im0,1:jm0) =
     &                       force_pbl_args_ws0_1(1:im0,1:jm0)
         force_pbl_args_tprime_0(1:im0,1:jm0) =
     &                       force_pbl_args_tprime_1(1:im0,1:jm0) 
         force_pbl_args_qprime_0(1:im0,1:jm0) = 
     &                       force_pbl_args_qprime_1(1:im0,1:jm0)

!        Assign forcings to time t1 (current forcing value)
         force_Ca_1(1:im0,1:jm0)=force_Ca_2(1:im0,1:jm0)
         force_cos_zen_angle_1(1:im0,1:jm0)= 
     &                        force_cos_zen_angle_2(1:im0,1:jm0)
         force_vis_rad_1(1:im0,1:jm0) = force_vis_rad_2(1:im0,1:jm0)
         force_direct_vis_rad_1(1:im0,1:jm0) = 
     &                        force_direct_vis_rad_2(1:im0,1:jm0)
         force_prec_ms_1(1:im0,1:jm0) = force_prec_ms_2(1:im0,1:jm0)
         force_eprec_w_1(1:im0,1:jm0) = force_eprec_w_2(1:im0,1:jm0)
         force_sprec_ms_1(1:im0,1:jm0) = force_sprec_ms_2(1:im0,1:jm0)
         force_seprec_w_1(1:im0,1:jm0) = force_seprec_w_2(1:im0,1:jm0)
         force_srheat_1(1:im0,1:jm0) = force_srheat_2(1:im0,1:jm0)
         force_trheat_1(1:im0,1:jm0) = force_trheat_2(1:im0,1:jm0)
         force_ts_1(1:im0,1:jm0) = force_ts_2(1:im0,1:jm0)
         force_qs_1(1:im0,1:jm0) = force_qs_2(1:im0,1:jm0)
         force_ps_1(1:im0,1:jm0) = force_ps_2(1:im0,1:jm0)
         force_rhosrf_1(1:im0,1:jm0) = force_rhosrf_2(1:im0,1:jm0)
         force_cdh_1(1:im0,1:jm0) = force_cdh_2(1:im0,1:jm0)
         force_qm1_1(1:im0,1:jm0) = force_qm1_2(1:im0,1:jm0)
         force_ws_1(1:im0,1:jm0) = force_ws_2(1:im0,1:jm0)
         force_pbl_args_ws0_1(1:im0,1:jm0) =
     &                       force_pbl_args_ws0_2(1:im0,1:jm0)
         force_pbl_args_tprime_1(1:im0,1:jm0) =
     &                       force_pbl_args_tprime_2(1:im0,1:jm0) 
         force_pbl_args_qprime_1(1:im0,1:jm0) = 
     &                       force_pbl_args_qprime_2(1:im0,1:jm0)


!        Obtain and assign forcings for t2 (next forcing value)
         if(.NOT.end_of_input_flag)then
            call open_gswp_files(step_num_3hr,gyr,itime_3hr
     &           ,n_month,n_month_sim,iu_vector,0)
            call get_gswp_data(n_month_sim,iu_vector,
     o         data_srheat (1:im0,1:jm0)          ,
     o         data_trheat (1:im0,1:jm0)          ,
     o         data_ts (1:im0,1:jm0)              ,
     o         data_qs (1:im0,1:jm0)              ,
     o         data_ps (1:im0,1:jm0)              ,
     o         data_ws (1:im0,1:jm0)              ,
     o         data_rainf( 1:im0,1:jm0)           ,
     o         data_rainf_c (1:im0,1:jm0)         ,
     o         data_snowf (1:im0,1:jm0) )
            call assign_gswp_forcings(gday,gyr,jtime,nday,dtsec,
     i         tbcs_ij(1:im0,1:jm0)              ,
     i         tcanopy_ij(1:im0,1:jm0)           ,
     i         data_srheat (1:im0,1:jm0)         ,
     i         data_trheat (1:im0,1:jm0)         ,
     i         data_ts (1:im0,1:jm0)             ,
     i         data_qs (1:im0,1:jm0)             ,
     i         data_ps (1:im0,1:jm0)             ,
     i         data_ws (1:im0,1:jm0)             ,
     i         data_rainf (1:im0,1:jm0)          ,
     i         data_rainf_c (1:im0,1:jm0)        ,
     i         data_snowf (1:im0,1:jm0)          ,
     o         force_Ca_2 (1:im0,1:jm0)             ,
     o         force_cos_zen_angle_2 (1:im0,1:jm0)  ,
     o         force_vis_rad_2 (1:im0,1:jm0)        ,
     o         force_direct_vis_rad_2 (1:im0,1:jm0) ,
     o         force_prec_ms_2 (1:im0,1:jm0)        ,
     o         force_eprec_w_2 (1:im0,1:jm0)        ,
     o         force_sprec_ms_2 (1:im0,1:jm0)       ,
     o         force_seprec_w_2 (1:im0,1:jm0)       ,
     o         force_srheat_2 (1:im0,1:jm0)         ,
     o         force_trheat_2 (1:im0,1:jm0)         ,
     o         force_ts_2 (1:im0,1:jm0)             ,
     o         force_qs_2 (1:im0,1:jm0)             ,
     o         force_ps_2 (1:im0,1:jm0)             ,
     o         force_rhosrf_2 (1:im0,1:jm0)         ,
     o         force_cdh_2 (1:im0,1:jm0)            ,
     o         force_qm1_2 (1:im0,1:jm0)            ,
     o         force_ws_2 (1:im0,1:jm0)             ,
     o         force_pbl_args_ws0_2 (1:im0,1:jm0)   ,
     o         force_pbl_args_tprime_2 (1:im0,1:jm0),
     o         force_pbl_args_qprime_2 (1:im0,1:jm0) )
         endif
        
         loop_latitude_j: do j=1,jm0
            loop_longitude_i: do i = 1,im0
               
!           drv_finterp(ip,trp_flag,val0,val1,val2,val3,madtt,valtrp)
!           "I" or "i" = instantaneous value at current time (linear interpol.)
!           ip = Grid point in question
!           val0 = last forcing value (for cases N & C)
!           val1 = current forcing value(all cases)
!           val2 = next forcing value (all cases but default)
!           val3 = ueber-next forcing value (for cases C & L)
!           madtt = Number of valtrp timesteps in (MUST BE MULTIPLE OF 2)
!           valtrp = Interpolated forcing data vector
               ip = j + (i-1)*jm0

!              Driver 1
               call drv_finterp(ip,'I',
     &                 force_Ca_0(i,j),
     &                 force_Ca_1(i,j),
     &                 force_Ca_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_Ca(i,j,:) = valtrp(:)
!              Driver 2
               call drv_finterp(ip,'I',
     &                 force_cos_zen_angle_0(i,j),
     &                 force_cos_zen_angle_1(i,j),
     &                 force_cos_zen_angle_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_cos_zen_angle(i,j,:) = valtrp(:)
!              Driver 3
               call drv_finterp(ip,'I',
     &                 force_vis_rad_0(i,j),
     &                 force_vis_rad_1(i,j),
     &                 force_vis_rad_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_vis_rad(i,j,:) = valtrp(:)
!              Driver 4
               call drv_finterp(ip,'I',
     &                 force_direct_vis_rad_0(i,j),
     &                 force_direct_vis_rad_1(i,j),
     &                 force_direct_vis_rad_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_direct_vis_rad(i,j,:) = valtrp(:)
!              Driver 5
               call drv_finterp(ip,'P',
     &                 force_prec_ms_0(i,j),
     &                 force_prec_ms_1(i,j),
     &                 force_prec_ms_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_prec_ms(i,j,:) = valtrp(:)
!              Driver 6
               call drv_finterp(ip,'P',
     &                 force_eprec_w_0(i,j),
     &                 force_eprec_w_1(i,j),
     &                 force_eprec_w_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_eprec_w(i,j,:) = valtrp(:)
!              Driver 7
               call drv_finterp(ip,'P',
     &                 force_sprec_ms_0(i,j),
     &                 force_sprec_ms_1(i,j),
     &                 force_sprec_ms_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_sprec_ms(i,j,:) = valtrp(:)
!              Driver 8
               call drv_finterp(ip,'P',
     &                 force_seprec_w_0(i,j),
     &                 force_seprec_w_1(i,j),
     &                 force_seprec_w_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_seprec_w(i,j,:) = valtrp(:)
!              Driver 9
               call drv_finterp(ip,'I',
     &                 force_srheat_0(i,j),
     &                 force_srheat_1(i,j),
     &                 force_srheat_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_srheat(i,j,:) = valtrp(:)
!              Driver 10
               call drv_finterp(ip,'I',
     &                 force_trheat_0(i,j),
     &                 force_trheat_1(i,j),
     &                 force_trheat_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_trheat(i,j,:) = valtrp(:)
!              Driver 11
               call drv_finterp(ip,'I',
     &                 force_ts_0(i,j),
     &                 force_ts_1(i,j),
     &                 force_ts_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_ts(i,j,:) = valtrp(:)
!              Driver 12
               call drv_finterp(ip,'I',
     &                 force_qs_0(i,j),
     &                 force_qs_1(i,j),
     &                 force_qs_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_qs(i,j,:) = valtrp(:)
!              Driver 13
               call drv_finterp(ip,'I',
     &                 force_ps_0(i,j),
     &                 force_ps_1(i,j),
     &                 force_ps_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_ps(i,j,:) = valtrp(:)
!              Driver 14
               call drv_finterp(ip,'I',
     &                 force_rhosrf_0(i,j),
     &                 force_rhosrf_1(i,j),
     &                 force_rhosrf_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_rhosrf(i,j,:) = valtrp(:)
!              Driver 15
               call drv_finterp(ip,'I',
     &                 force_cdh_0(i,j),
     &                 force_cdh_1(i,j),
     &                 force_cdh_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_cdh(i,j,:) = valtrp(:)
!              Driver 16
               call drv_finterp(ip,'I',
     &                 force_qm1_0(i,j),
     &                 force_qm1_1(i,j),
     &                 force_qm1_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_qm1(i,j,:) = valtrp(:)
!              Driver 17
               call drv_finterp(ip,'I',
     &                 force_ws_0(i,j),
     &                 force_ws_1(i,j),
     &                 force_ws_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_ws(i,j,:) = valtrp(:)
!              Driver 18
               call drv_finterp(ip,'I',
     &                 force_pbl_args_ws0_0(i,j),
     &                 force_pbl_args_ws0_1(i,j),
     &                 force_pbl_args_ws0_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
!              Driver 19
               vector_pbl_args_ws0(i,j,:) = valtrp(:)
              call drv_finterp(ip,'I',
     &                 force_pbl_args_tprime_0(i,j),
     &                 force_pbl_args_tprime_1(i,j),
     &                 force_pbl_args_tprime_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_pbl_args_tprime(i,j,:) = valtrp(:)
!              Driver 20
              call drv_finterp(ip,'I',
     &                 force_pbl_args_qprime_0(i,j),
     &                 force_pbl_args_qprime_1(i,j),
     &                 force_pbl_args_qprime_2(i,j),0.d0,
     &                 n_LSM_GSWP,valtrp)
               vector_pbl_args_qprime(i,j,:) = valtrp(:)

            enddo loop_longitude_i
         enddo loop_latitude_j
         !print *, 'After loop to interpolate forcings'
      endif

!     Assign forcings to LSM based on interpolation of GSWP data
      force_Ca(1:im0,1:jm0)=vector_Ca(1:im0,1:jm0,i_sub3hr)
      force_cos_zen_angle(1:im0,1:jm0)= 
     &                      vector_cos_zen_angle(1:im0,1:jm0,i_sub3hr)
      force_vis_rad(1:im0,1:jm0) = 
     &                      vector_vis_rad(1:im0,1:jm0,i_sub3hr)
      force_direct_vis_rad(1:im0,1:jm0) = 
     &                      vector_direct_vis_rad(1:im0,1:jm0,i_sub3hr)
      force_prec_ms(1:im0,1:jm0) = 
     &                      vector_prec_ms(1:im0,1:jm0,i_sub3hr)
      force_eprec_w(1:im0,1:jm0) = 
     &                      vector_eprec_w(1:im0,1:jm0,i_sub3hr)
      force_sprec_ms(1:im0,1:jm0) = 
     &                      vector_sprec_ms(1:im0,1:jm0,i_sub3hr)
      force_seprec_w(1:im0,1:jm0) = 
     &                      vector_seprec_w(1:im0,1:jm0,i_sub3hr)
      force_srheat(1:im0,1:jm0) = vector_srheat(1:im0,1:jm0,i_sub3hr)
      force_trheat(1:im0,1:jm0) = vector_trheat(1:im0,1:jm0,i_sub3hr)
      force_ts(1:im0,1:jm0) = vector_ts(1:im0,1:jm0,i_sub3hr)
      force_qs(1:im0,1:jm0) = vector_qs(1:im0,1:jm0,i_sub3hr)
      force_ps(1:im0,1:jm0) = vector_ps(1:im0,1:jm0,i_sub3hr)
      force_rhosrf(1:im0,1:jm0) = vector_rhosrf(1:im0,1:jm0,i_sub3hr)
      force_cdh(1:im0,1:jm0) = vector_cdh(1:im0,1:jm0,i_sub3hr)
      force_qm1(1:im0,1:jm0) = vector_qm1(1:im0,1:jm0,i_sub3hr)
      force_ws(1:im0,1:jm0) = vector_ws(1:im0,1:jm0,i_sub3hr)
      force_pbl_args_ws0(1:im0,1:jm0) =
     &                    vector_pbl_args_ws0(1:im0,1:jm0,i_sub3hr)
      force_pbl_args_tprime(1:im0,1:jm0) =
     &                    vector_pbl_args_tprime(1:im0,1:jm0,i_sub3hr) 
      force_pbl_args_qprime(1:im0,1:jm0) = 
     &                    vector_pbl_args_qprime(1:im0,1:jm0,i_sub3hr)

!     Update cosz every timestep rather than use interpolated values (Y.KIM)
      rot1=(2d0*pi*jtime)/nday ! beginning of timestep
      rot2=rot1+2d0*pi*dtsec/dble(sec_1day) ! end of timestep
      Rg(:,:) = force_srheat(:,:)
      CALL COSZT (gday,gyr,rot1,rot2,Rg,cosz1,fdif)
      force_cos_zen_angle = cosz1
      force_direct_vis_rad=(1-fdif)*force_vis_rad
#endif
      end subroutine get_gswp_forcings


      end module drv_gswp_force
