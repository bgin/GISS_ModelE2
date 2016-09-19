# rundeck for compiling stand-alone version of GISS LSM
# coupled to Ent with FBB photosynthesis at 2x2.5 resolution

# global objects to use:
OBJ_LIST = main_lsm

# components to use:
COMPONENTS = giss_LSM_standalone  Ent  giss_LSM    shared solvers

# specific options for components:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES OFFLINE_RUN=YES
OPTS_shared =
OPTS_giss_LSM_standalone = USE_ENT=YES
OPTS_solvers = OFFLINE_RUN=YES

INPUT_FILES = \
fort.951=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.951 \
fort.952=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.952 \
fort.953=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.953 \
SOIL=$(GCMSEARCHPATH)/S144X900098M \
TOP_INDEX=$(GCMSEARCHPATH)/top_index_144x90.ij.ext \
VEG=$(GCMSEARCHPATH)/V144X90_no_crops \
CROPS=$(GCMSEARCHPATH)/CROPS_144X90N.nocasp.ext \
soil_textures=$(GCMSEARCHPATH)/soil_textures_top30cm_2x2.5 \
SOILCARB_global=$(GCMSEARCHPATH)/soilcarb_top30cm_nmaps_2x2.5bin.dat \
CD_coef=$(GCMSEARCHPATH)/CD144X90.ext

# run settings (one value per line, single quotation marks for strings)
# i_site, j_site, & latd are for ecosystem-scale runs; pass dummies otherwise
# For logicals, "0" is ".false." and "1" is ".true."
# run settings (one value per line, single quotation marks for strings)
# i_site, j_site, & latd are for ecosystem-scale runs; pass dummies otherwise
# For logicals, "0" is ".false." and "1" is ".true."
# For spinups, set beginning and end of spinnup period, and number of times
# to repeat the period
RUN_PARAMETERS = \
&&PARAMETERS \n\
i_site=38 \n\
j_site=65 \n\
latd=39.32315d0 \n\
itime1=280512 \n\
jyear=1983 \n\
jday=1 \n\
tyr_sec=0 \n\
dtsec=1800.d0 \n\
spin_start_yr=1983 \n\
spin_end_yr=1985 \n\
num_times_spin=2 \n\
force_VEG_bi=0 \n\
do_soilinit_bi=1 \n\
do_soilresp_bi=1 \n\
do_phenology_activegrowth_bi=0 \n\
do_structuralgrowth_bi=0 \n\
do_frost_hardiness_bi=1 \n\
do_patchdynamics_bi=0 \n\
do_spinup_bi=1 \n\
&&END_PARAMETERS

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