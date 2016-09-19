E_AR5_V2_Dust_NCEP_F40.R GISS Model E  2000 ocn/atm   jan perlwitz 07/16/2014
E_AR5_V2_Dust_NCEP_F40: E_AR5_V2_Dust_F40 + nudged with NCEP winds
(When using, please cite reference for description of dust model:
 Miller, R.L., R.V. Cakmur, J. Perlwitz, I.V. Geogdzhayev, P. Ginoux,
   K.E. Kohfeld, D. Koch, C. Prigent, R. Ruedy, G.A. Schmidt, and I. Tegen,
   2006: Mineral dust aerosols in the NASA Goddard Institute for Space Sciences
   ModelE atmospheric general circulation model. J. Geophys. Res., 111, D06208,
   doi:10.1029/2005JD005796, and others within)

E_AR5_V2_NINT is for AR5 Prime production runs; 
                 U00a=0.65, U00b=1.0, WMUI_multiplier=1.; 
                 grav.wave adjustment: CMTN=0.1, CDEF=1.6 

!! delete lines starting with '!!' unless E4F40 prepares a q-flux ocean run
!! E4qsF40.R GISS Model E  1850 atm, ocn: q-flux 65m             rar 07/15/2009

!! E4qsF40 = E4F40 with 65m q-flux ocean
E4F40 = modelE as frozen in April 2010:
modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 2000
ocean data: prescribed, 1996-2005 climatology
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define NEW_IO                   ! new I/O (netcdf) on
!#define NEW_IO_SUBDD             ! new I/O (netcdf) for subdd on
#define TRAC_ADV_CPU             ! timing index for tracer advection on
#define NUDGE_ON                 ! nudged winds on
#define USE_ENT                  ! include dynamic vegetation model
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! wet deposition and water tracer
#define TRACERS_DUST             ! include dust tracers
#define TRACERS_DUST_Silt4       ! include 4th silt size class of dust
!#define TRACERS_DUST_Silt5       ! include 5th silt size class of dust (experimental, do not use!)
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
#define RAD_O3_GCM_HRES          ! Use GCM horiz resl to input rad code clim Ozone
#define RAD_O3_DECADAL_INPUT
#define NO_HDIURN                ! exclude hdiurn diagnostics
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
RES_stratF40                        ! horiz/vert resolution, 2x2.5, top at 0.1mb, 40 layers
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform

    ! lat-lon grid specific source codes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT_netcdf                ! diagn/post-processing output
IO_DRV                              ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV                             ! advection of T
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

QUS3D                               ! advection of Q/tracers
#include "tracer_dust_source_files"
#include "tracer_shared_source_files"
TRDIAG

#include "modelE4_V2_source_files"
#include "static_ocn_source_files"

NUDGE                               ! code for nudging winds

Components:
#include "E4_components_nc"    /* without "Ent" */
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT    /* needed for "Ent" only */
OPTS_dd2d = NC_IO=PNETCDF

Data input files:
#include "IC_144x90_input_files_nc"
#include "static_ocn_2000_144x90_input_files"

RVR=RD_Fb.RVR.bin          ! river direction file

#include "land144x90_input_files_V2"
#include "nudging_144x89_input_files"
#include "rad_input_files_V2"
#include "TAero2008_input_files_V2"
#include "O3_AR5_V2_144x90_input_files"

#include "dust_tracer_144x90_input_files"

#include "dry_depos_144x90_input_files"

MSU_wts=MSU.RSS.weights.data      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
E_AR5_V2_Dust_NCEP_F40 (E_AR5_V2_Mine_F40 with year 2000 atm., 1996-2005 clim ocn, dust tracers, nudged with NCEP winds)

&&PARAMETERS
#include "static_ocn_params"
#include "sdragF40_params"
#include "gwdragF40_params_V2"

xCDpbl=1.
cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
FT8OPX=1.,1.,1.,1.,1.,1.,1.,1.

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.65      ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00      ! below 850mb and MC regions; then tune this to get rad.balance
WMUI_multiplier = 1.

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=1.      ! activates strat.H2O generated by CH4
KSIALB=0         ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

#include "atmCompos_2000_params"
!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use abs(o3_yr) instead!
!!!!!!!!!!!!!!!!!!!!!!!
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols
#include "aerosol_params"
#include "dust_params_V2"

DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)
#include "diag_params"

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960
&&END_PARAMETERS

 &INPUTZ
 YEARI=1995,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=2006,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1995,MONTHE=12,DATEE=1,HOURE=1,
 &END
