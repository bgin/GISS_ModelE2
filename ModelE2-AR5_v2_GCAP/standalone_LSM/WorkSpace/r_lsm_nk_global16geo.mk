# r_lsm_nk_global.mk
# rundeck for compiling stand-alone version of GISS LSM
# coupled to Ent with FBB photosynthesis for 4x5 resolution
# - Ent 16 PFTs, EntMM cover, Simard heights, MODIS initial LAI, 4x5 resolution for TESTING.
# - GSWP2 meterological forcings CHECK RESOLUTION
# - Run can be global or select grid cells.
# - See lsm_standalone_nk.f:process_input_parameters for additional
#   optional array boundary input parameters and defaults.

# global objects to use:
OBJ_LIST = main_lsm

# components to use:
COMPONENTS = giss_LSM_standalone  ../../model/Ent/  ../../model/giss_LSM    ../../model/shared ../../model/solvers

# specific options for components:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT
OPTS_giss_LSM = USE_ENT=YES OFFLINE_RUN=YES 
OPTS_shared =
OPTS_giss_LSM_standalone = USE_ENT=YES LSM_DRV=GLOBAL
OPTS_solvers = OFFLINE_RUN=YES

INPUT_FILES = \
fort.952=$(GCMSEARCHPATH)/lsm_data_dump_fort.952 \
TOPO=$(GCMSEARCHPATH)/Z72X46N.cor4_nocasp \
SOIL=$(GCMSEARCHPATH)/S4X50093.ext \
TOP_INDEX=$(GCMSEARCHPATH)/top_index_72x46.ij.ext \
VEG=$(VEGPATH)/V72x46_EntMM16_lc_max_trimmed_scaled_nocrops.ij \
CROPS=$(GCMSEARCHPATH)/CROPS_72X46N.cor4 \
soil_textures=$(GCMSEARCHPATH)/soil_textures_top30cm \
SOILCARB_global=$(GCMSEARCHPATH)/soilcarb_top30cm_nmaps_4x5bin.dat \
CD_coef=$(GCMSEARCHPATH)/CD4X500S \
LAImax=$(VEGPATH)/V72x46_EntMM16_lai_max_trimmed_scaled.ij
HITEent=$(VEGPATH)/V72x46_EntMM16_height_trimmed_scaled.ij \
LAIinit=$(VEGPATH)/V72x46_EntMM16_lai_trimmed_scaled_Jan.ij \
# Run settings (one value per line, single quotation marks for strings)
# IM,JM - grid dimensions, e.g. 72x46, 144x90
# I0,I1,J0,J1 - run boundaries; generally global 1,IM,1,JM, or single grid cells.
# year1 - first year of run
# year2 - last year of run
# tyr_sec1 - start of run in seconds from begin of year1
# tyr_sec2 - end of run in seconds from begin of year2, default is end of year.
#    This allows runs less than a year long, or starting in middle of year.
# Mdtsec - time step of input met forcing data set (seconds)
# run settings (one value per line, single quotation marks for strings)
# For spinups, set beginning and end of spinup period, and number of times.
#   If do_spinup==.false., this overrides any spinup time inputs.
# to repeat the period
RUN_PARAMETERS = \
&PARAMETERS \n\
IM = 72 \n\
JM = 46 \n\
I0=1, I1=72, J0=1, J1=46 \n\
year1=1983, year2=1995\n\
tyr_sec1 = 0  \n\
Mdtsec = 10800 \n\
spin_start_yr=1983 \n\
spin_end_yr=1985 \n\
num_times_spin=0 \n\
force_VEG=.false. \n\
do_soilinit=.true. \n\
do_soilresp=.true. \n\
do_phenology_activegrowth=.true. \n\
do_structuralgrowth=.false. \n\
do_frost_hardiness=.true. \n\
do_patchdynamics=.false. \n\
do_init_geo=.true. \n\
do_spinup=.false. \n\
BASEFOLD_OUT='$(GSWP_DATADIR)' \
/

# 2.5 deg long x 2 deg lat Resolution   
#      integer,parameter :: i_site=38,j_site=65 !MMSF
#      integer,parameter :: i_site = 82,j_site =76 !Hyytiala
#      integer,parameter :: i_site=24,j_site=65 !Vaira
#      integer,parameter :: i_site=43,j_site=66 !NYC 
#      integer,parameter :: i_site=51,j_site=44 !TNF
#      integer,parameter :: i_site=87,j_site=46 !Mpala
# 5 deg long x 4 deg lat Resolution   
#      integer,parameter :: i_site=19,j_site=33 !MMSF
#      integer,parameter :: i_site = 40,j_site =39 !Hyytiala
#      integer,parameter :: i_site=12,j_site=33 !Vaira
#      integer,parameter :: i_site=26,j_site=23 !TNF
# Latitude of each site
#      real*8,parameter :: latd = 39.32315d0  !MMSF
#      real*8,parameter :: latd = 61.8474d0   !Hyytiala
#      real*8,parameter :: latd = 38.4067d0   !Vaira
#      real*8,parameter :: latd = 40.6222d0   !NYC
#      real*8,parameter :: latd = -2.8567d0   !TNF
#      real*8,parameter :: latd =  0.3000d0   !Mpala
#!!One year->itime1=17520!spin83_85=52608;86_95=175296

#Locations of BASEFOLD_OUT
#BASEFOLD_OUT=/discover/noback/mpuma/lsm_gswp
#BASEFOLD_OUT=/discover/noback/mpuma/lsm_gswp_2x2_5
#BASEFOLD_OUT=/n/Moorcroft_Lab/Users/kim/ent-data/GSWP2

# use the following parameters for "athena"
#RUN_PARAMETERS = \
#&&PARAMETERS \n\
#basefold_out='/clima1/GSWP2' \n\
#neutral_drag_file='/u/cmrun/CDN5X4.GISS' \n\
#&&END_PARAMETERS

# the following are settings for "quark"
#RUN_PARAMETERS = \
#&&PARAMETERS \n\
#basefold_out='/a11/tmp/GSWP_1983' \n\
#neutral_drag_file='/a11/tmp/GSWP_1983/CDN5X4.GISS' \n\
#&&END_PARAMETERS


