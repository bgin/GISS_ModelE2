E_AR5_CarbonCycle GISS Model E  coupled version          natassa   10/17/2010
generated from E4arobio_g template, AR5 branch 2/6/14

obio rundeck for merged AR5_Ent_devel_2 branch
based on E_AR5_NINT_oR.R GISS Model E  coupled version          larissa   10/08/2010
!! E_AR5_NINT_oR is with NetCDF output;  
               + WMUI_multiplier=2. (to adjust Planetary albedo close to 30%) 
                 (U00a=0.54; U00b=1.0)
E4F40oQ32: E4F40 coupled to 1x1.25deg 32-layer GISS ocean model
E4F40 = modelE as frozen in April 2010:
modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean: coupled to 1x1.25deg 32-layer GISS ocean model (Russell - Schmidt)
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define USE_ENT
#define NEW_IO
#define CHECK_OCEAN                 ! needed to compile aux/file CMPE002
#define TRACERS_ON                  ! include tracers code
#define OBIO_ON_GARYocean           ! obio on Russell ocean
#define TRACERS_OCEAN               ! Gary's Ocean tracers activated
#define TRACERS_OCEAN_INDEP         ! independently defined ocn tracers
#define TRACERS_OceanBiology
#define pCO2_ONLINE
#define OBIO_RAD_coupling
#define TRACERS_GASEXCH_ocean       ! ANY ocean: special tracers to be passed to ocean
#define TRACERS_GASEXCH_ocean_CO2   ! ANY ocean: special tracers to be passed to ocean
#define TRACERS_GASEXCH_land
#define TRACERS_GASEXCH_land_CO2
!!!#define constCO2
!!!#define restoreIRON
!!!#define TRACERS_Alkalinity
!!!#define Jprod_based_on_pp
!!!!#define CHL_from_OBIO               ! ANY ocean: interactive CHL
!!!!#define CHL_from_SeaWIFs
!!!!#define change_PNOICE           ! adjust ice-obio interactions
#define MODIS_LAI
#define NEW_IO_VEG
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
RES_stratF40                        ! horiz/vert resolution, 2x2.5, top at 0.1mb, 40 layers
ORES_1Qx1_L32                       ! ocean horiz res 2x2.5deg, 32 vert layers
OSTRAITS_1QX1_COM                   ! dynamic ocean modules
DIAG_RES_F                          ! diagnostics
FFT144 OFFT288E                     ! Fast Fourier Transform

    ! lat-lon grid specific source codes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
IO_DRV                              ! new i/o 

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

MODEL_COM                           ! model variables 
MODELE                              ! Main and model overhead
ALLOC_DRV                           ! domain decomposition, allocate global distributed arrays
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM     ! land surface and soils + snow model
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + Ent          ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation, albedo, prescribed aerosols
DIAG_COM DIAG DEFACC                ! diagnostics
ODIAG_COM OCEAN_COM OGEOM  ! dynamic ocean modules
OCNDYN OCNDYN2                            ! dynamic ocean routines
OCN_Interp                 ! dynamic ocean routines
OSTRAITS OCNGM OCNKPP                     ! dynamic ocean routines
OCEANR_DIM AFLUXES OFLUXES
ODIAG_PRT                              ! ocean diagnostic print out
OCNFUNTAB                              ! ocean function look up table
SparseCommunicator_mod                 ! sparse gather/scatter module
OCNQUS                     ! routines for advection using the QUS
OCN_Int_LATLON                      ! atm-ocn regrid routines

     ! atmospheric tracers
TRACER_COM  TRACERS_DRV              ! configurable tracer code
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
TRDIAG
TRACER_GASEXCH_CO2                  ! tracer functions needed for gas exch expts
!!!!TRACER_GASEXCH_CFC                 ! tracer functions needed for gas exch expts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!   OCEAN TRACERS       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OCN_TRACER_COM
OCN_TRACER


     ! ocean carbon cycle
obio_dim         |-r8|
obio_incom       |-r8|
obio_com         |-r8|
obio_forc        |-r8|
obio_init        |-r8|
obio_bioinit_g   |-r8|
obio_model       |-r8|
obio_trint       |-r8|
obio_daysetrad   |-r8|
obio_daysetbio   |-r8|
obio_ocalbedo    |-r8|
obio_reflectance |-r8|
obio_sfcirr      |-r8|
obio_edeu        |-r8|
obio_ptend       |-r8|
obio_carbon      |-r8|
obio_update      |-r8|
obio_sinksettl   |-r8|
!!!obio_alkalinity  |-r8|

Components:
shared ESMF_Interface solvers giss_LSM 
dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT
OPTS_giss_LSM = USE_ENT=YES

Data input files:
    ! start from the restart file of an earlier run ...                 ISTART=8
!! AIC=1JAN3001.rsfE114F40oQ32.nc ! initial cond. are the previously spun-up model, 
!!                               no GIC needed, use ISTART=8

!! AIC=1JAN1961.rsfE4F40.MXL65m  ! end of run with KOCEAN=0

    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
AIC=AIC.RES_F40.D771201        ! observed init cond (atm. only)
GIC=GIC.144X90.DEC01.1.ext.nc  ! initial ground conditions
OIC=OIC288X180.D1201                ! Levitus ocean intial conditions
OFTAB=OFTABLE_NEW                   ! ocean function table
AVR=OPF.E1QX1.L32                   ! ocean filter
KBASIN=KB288X180.modelE.BS1         ! ocean basin designations (1 cell Bering Strait)
TOPO_OC=Z1QX1N.BS1                  ! ocean fraction and topography (1 cell Bering Strait)
TOPO=Z2HX2fromZ1QX1N.BS1            ! surface fractions and topography (1 cell Bering Strait)
!!!AIC=/discover/nobackup/projects/giss/prod_input_files/1JAN3051.rsfE119F40oQ32.nc


RVR=RD_Fb.RVR.bin                   ! river direction file (frac. ocean)

CDN=CD144X90.ext
!VEG=V144X90_no_crops.ext
VEG=/discover/nobackup/nkiang/DATA/Vegcov/V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.nc
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp
SOIL=S144X900098M.ext
TOP_INDEX=top_index_144x90_a.ij.ext
ZVAR=ZVAR2X25A             ! topographic variation for gwdrag
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_nmaps_2x2.5bin.dat
GLMELT=GLMELT_144X90_gas.OCN
    ! resolution independent files
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k.1a              ! rad.tables and history files
RADN4=LWCorrTables33k              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2
! updated aerosols need MADAER=3
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust_mass_CakmurMillerJGR06_72x46x20x7x12
RADN7=STRATAER.VOL.SATO.1850-1999.Apr02_hdr
RADN8=cloud.epsilon4.72x46
RADN9=solar.DBbglean.ann850-2000.uvflux_hdr       ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
! ozone files (minimum 1, maximum 9 files + 1 trend file)
GHG=GHG.Mar2009.txt ! use GHG.Jul2009.txt for runs that start before 1850
dH2O=dH2O_by_CH4_monthly
BC_dep=BC.Dry+Wet.depositions.ann
! updated aerosols need MADAER=3
TAero_SUL=SUL_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_SSA=SSA_Koch2008_kg_m2_72x46x20h
TAero_NIT=NIT_Bauer2008_kg_m2_72x46x20_1890-2000h
TAero_OCA=OCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCA=BCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCB=BCB_Koch2008_kg_m2_72x46x20_1890-2000h
! ozone files (minimum 1, maximum 9 files + 1 trend file)
O3file_01=mar2004_o3_shindelltrop_72x46x49x12_1850
O3file_02=mar2004_o3_shindelltrop_72x46x49x12_1890
O3file_03=mar2004_o3_shindelltrop_72x46x49x12_1910
O3file_04=mar2004_o3_shindelltrop_72x46x49x12_1930
O3file_05=mar2004_o3_shindelltrop_72x46x49x12_1950
O3file_06=mar2004_o3_shindelltrop_72x46x49x12_1960
O3file_07=mar2004_o3_shindelltrop_72x46x49x12_1970
O3file_08=mar2005_o3_shindelltrop_72x46x49x12_1980
O3file_09=mar2005_o3_shindelltrop_72x46x49x12_1990
O3trend=mar2005_o3timetrend_46x49x2412_1850_2050

!!!!!!!!!!!!!!!!!!! obio  input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cfle1=abw25b.dat                         ! seawater spectral absorp.
                                         ! and scatt. coefs
cfle2=acbc25b.dat                        ! phytoplankton spectrl absorp.
                                         ! and scatt. coefs
!!!!pco2table=pco2.tbl.asc               ! table to compute pco2 values
                                         ! from sst,sss,dic,alk
                                         ! if not defined pCO2_ONLINE
nitrates_inicond=no3_nodc_annmean.asc    ! initial cond for nitrates (NODC)
silicate_inicond=sio2_nodc_annmean.asc   ! initial cond for silicate (NODC)
!!!dic_inicond=dic_glodap_annmean.asc    ! initial cond for dic (GLODAP)
dic_inicond=preindustrial_ocean_tco2.asc ! initial cond for dic (GLODAP)
alk_inicond=alk_glodap_annmean.asc       ! initial cond/forc for alk(GLODAP)
!!!oasimdirect=oasimdirect_20w_new       ! spectral light components
                                         ! if not def OBIO_RAD_coupling
atmFe_inicond=iron_gocart_1x1mon.asc     ! GOCART iron flux
atmFedirect1=iron_ron_195x180_20w.asc    ! Ron Miller's dust fluxes
facirr=facirr.asc                        ! factors for mean irrad w/in water
eda_esa_ratios=eda_esa_ratios.asc        ! ratios of rad spectrl components
!!!!!!!!!!!!!!!!!!! obio_rad  input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHL_DATA=CHL_WG_2x2.5zavg                !CHL_WG_4x5 or CHL_WG_2x2.5
!! CHL_DATA=CHL_WG_2x2.5                    !CHL_WG_4x5 or CHL_WG_2x2.5
                                         !in Gary'socean grid
                                         !to be used with CHL_from_SeaWIFs


MSU_wts=MSU.RSS.weights.data      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

LAIMAX=/discover/nobackup/nkiang/DATA/Vegcov/V144x90_EntMM16_lai_max_trimmed_scaled_ext1.nc 
HITEent=/discover/nobackup/nkiang/DATA/Vegcov/V144x90_EntMM16_height_trimmed_scaled_ext1.nc
LAI=/discover/nobackup/nkiang/DATA/Vegcov/V144x90_EntMM16_lai_trimmed_scaled_ext1.nc


Label and Namelist:  (next 2 lines)
E_AR5_CarbonCycle (AR5v2 model + obio, cold start, preindustrial CO2)


&&PARAMETERS
! parameters set for coupled ocean runs:
KOCEAN=1            ! ocn is prognostic
OBottom_drag=1      !  Drags at the ocean bottom (NO drags -> OBottom_drag=0)
OCoastal_drag=1     !  Drags at the ocean coasts (NO drags -> OCoastal_drag=0)
ocean_use_qus=1
variable_lk=1
init_flake=1
init_flake=0

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000045 !the default is 15d-6
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.2       ! default is 0.5
CDEF=1.5       ! deformation drag coefficient
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF


! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.54      ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00      ! below 850mb and MC regions; then tune this to get rad.balance

WMUI_multiplier = 2.

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=1.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=1850
!crops_yr=1850  ! if -1, crops in VEG-file is used
!s0_yr=1850
!s0_day=182
!ghg_yr=1850
!ghg_day=182
volc_yr=-1
!volc_day=182
!aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
!albsn_yr=1850
dalbsnX=.024
!o3_yr=-1850
!aer_int_yr=1850    !select desired aerosol emissions year or 0 to use JYEAR
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

!variable_orb_par=0
!orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols

DTO=112.5
DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)

! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time
SUBDD=' '        ! no sub-daily frequency diags
NSUBDD=0         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
KCOPY=1          ! saving acc + rsf
isccp_diags=1    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc

Nssw=48          ! until diurnal diags are fixed, Nssw has to be even
                 ! obio needs that in order to always restart from hour 0
                 ! then we need to do setups for a whole day
Ndisk=720


!parameters that affect CO2 gas exchange
!!! atmCO2=368.6      !uatm for year 2000 
 atmCO2=285.226           !uatm for new preindustrial runs
!!!atmCO2=0.             !prognostic atmCO2
to_volume_MixRat=1    ! for tracer printout
solFe=0.02            ! default iron solubility
!!!solFe=0.05            ! enhanced iron solubility

do_phenology_activegrowth=0

&&END_PARAMETERS

 &INPUTZ
   YEARI=1850,MONTHI=1,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
   YEARE=1860,MONTHE=1,DATEE=1,HOURE=0, KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1850,MONTHE=1,DATEE=2,HOURE=0,IWRITE=1,JWRITE=1,
 &END
### Information below describes your run. Do not delete! ###
Thu Feb  6 13:52:18 EST 2014
ifort version 13.1.3
Git reference: 2d4b950d098ad0a7692fba17258d390d5b14d1b9
