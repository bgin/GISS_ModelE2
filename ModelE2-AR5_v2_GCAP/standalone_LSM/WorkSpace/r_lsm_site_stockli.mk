# rundeck for compiling stand-alone version of GISS LSM
# coupled to Ent with FBB photosynthesis for a single site
# at the ecosystem scale with Stoekli meteorological drivers
# using soil data from 2x2.5 resolution files

# global objects to use:
OBJ_LIST = main_lsm

# components to use:
COMPONENTS = giss_LSM_standalone  Ent  giss_LSM    shared solvers

# specific options for components:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB FLUXNETINIT=YES
OPTS_giss_LSM = USE_ENT=YES OFFLINE_RUN=YES
OPTS_shared =
OPTS_giss_LSM_standalone = USE_ENT=YES
OPTS_solvers = OFFLINE_RUN=YES
#---------------------------------------------------------------------
LSM_DATADIR=/discover/nobackup/mpuma/lsm_single

# Morgan Monroe Input and Parameters
INPUT_FILES = \
fort.951=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.951 \
fort.952=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.952 \
fort.953=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.953 \
SOIL=$(GCMSEARCHPATH)/S144X900098M \
TOP_INDEX=$(GCMSEARCHPATH)/top_index_144x90.ij.ext \
VEG=$(LSM_DATADIR)/V144x90_MMSF.bi \
CROPS=$(GCMSEARCHPATH)/CROPS_144X90N.nocasp.ext \
soil_textures=$(GCMSEARCHPATH)/soil_textures_top30cm_2x2.5 \
SOILCARB_global=$(GCMSEARCHPATH)/soilcarb_top30cm_nmaps_2x2.5bin.dat \
CD_coef=$(GCMSEARCHPATH)/CD144X90.ext \
fluxnet_LAI=$(LSM_DATADIR)/MMSF2005_LAI_hourly.txt \
fluxnet_vheight=$(LSM_DATADIR)/MMSF_vegheight.txt \
site_forcings=$(LSM_DATADIR)/Morgan_Monroe_State_Forest_1999_2006.dat

RUN_PARAMETERS = \
&&PARAMETERS \n\
i_site = 38 \n\
j_site = 65 \n\
latd = 39.32315d0 \n\
itime1=61368 \n\
jyear=1999 \n\
jday = 1 \n\
tyr_sec = 0 \n\
dtsec=3600.d0 \n\
spin_start_yr=1999 \n\
spin_end_yr=1999 \n\
num_times_spin=1 \n\
force_VEG_bi=1 \n\
do_soilinit_bi=1 \n\
do_soilresp_bi=1 \n\
do_phenology_activegrowth_bi=0 \n\
do_structuralgrowth_bi=0 \n\
do_frost_hardiness_bi=1 \n\
do_patchdynamics_bi=0 \n\
do_spinup_bi=0 \n\
&&END_PARAMETERS
#---------------------------------------------------------------------
# Hyytiala Input and Parameters
#INPUT_FILES = \
#fort.951=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.951 \
#fort.952=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.952 \
#fort.953=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.953 \
#SOIL=$(GCMSEARCHPATH)/S144X900098M \
#TOP_INDEX=$(GCMSEARCHPATH)/top_index_144x90.ij.ext \
#VEG=$(LSM_DATADIR)/V144x90_Hyytiala.bi \
#CROPS=$(GCMSEARCHPATH)/CROPS_144X90N.nocasp.ext \
#soil_textures=$(GCMSEARCHPATH)/soil_textures_top30cm_2x2.5 \
#SOILCARB_global=$(GCMSEARCHPATH)/soilcarb_top30cm_nmaps_2x2.5bin.dat \
#CD_coef=$(GCMSEARCHPATH)/CD144X90.ext \
#fluxnet_LAI=$(LSM_DATADIR)/Hyytiala1998_LAI_30min.txt \
#fluxnet_vheight=$(LSM_DATADIR)/Hyytiala_vegheight.txt \
#site_forcings=$(LSM_DATADIR)/Hyytiala_1997_2005.dat
#
#RUN_PARAMETERS = \
#&&PARAMETERS \n\
#i_site = 82 \n\
#j_site = 76 \n\
#latd = 61.8474d0 \n\
#itime1=157776 \n\
#jyear=1997 \n\
#jday = 1 \n\
#tyr_sec = 0 \n\
#dtsec=1800.d0 \n\
#spin_start_yr=1997 \n\
#spin_end_yr=2005 \n\
#num_times_spin=1 \n\
#force_VEG_bi=1 \n\
#do_soilinit_bi=1 \n\
#do_soilresp_bi=1 \n\
#do_phenology_activegrowth_bi=0 \n\
#do_structuralgrowth_bi=0 \n\
#do_frost_hardiness_bi=1 \n\
#do_patchdynamics_bi=0 \n\
#do_spinup_bi=0 \n\
#&&END_PARAMETERS
##---------------------------------------------------------------------
#---------------------------------------------------------------------
## Vaira Input and Parameters
#INPUT_FILES = \
#fort.951=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.951 \
#fort.952=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.952 \
#fort.953=$(GCMSEARCHPATH)/lsm_data_dump_f_fort.953 \
#SOIL=$(GCMSEARCHPATH)/S144X900098M \
#TOP_INDEX=$(GCMSEARCHPATH)/top_index_144x90.ij.ext \
#VEG=$(LSM_DATADIR)/V144x90_Vaira.bi \
#CROPS=$(GCMSEARCHPATH)/CROPS_144X90N.nocasp.ext \
#soil_textures=$(GCMSEARCHPATH)/soil_textures_top30cm_2x2.5 \
#SOILCARB_global=$(GCMSEARCHPATH)/soilcarb_top30cm_nmaps_2x2.5bin.dat \
#CD_coef=$(GCMSEARCHPATH)/CD144X90.ext \
#fluxnet_LAI=$(LSM_DATADIR)/Vaira2002_LAI_30min.txt \
#fluxnet_vheight=$(LSM_DATADIR)/Vaira_vegheight.txt \
#site_forcings=$(LSM_DATADIR)/Vaira_2001_2006.dat
#
#RUN_PARAMETERS = \
#&&PARAMETERS \n\
#i_site = 24 \n\
#j_site = 65 \n\
#latd = 38.4067d0 \n\
#itime1=105168 \n\
#jyear=2001 \n\
#jday = 1 \n\
#tyr_sec = 0 \n\
#dtsec=1800.d0 \n\
#spin_start_yr=2001 \n\
#spin_end_yr=2006 \n\
#num_times_spin=1 \n\
#force_VEG_bi=1 \n\
#do_soilinit_bi=1 \n\
#do_soilresp_bi=1 \n\
#do_phenology_activegrowth_bi=0 \n\
#do_structuralgrowth_bi=0 \n\
#do_frost_hardiness_bi=1 \n\
#do_patchdynamics_bi=0 \n\
#do_spinup_bi=0 \n\
#&&END_PARAMETERS

##---------------------------------------------------------------------
# 2.5 deg long x 2 deg lat Resolution   
#      integer,parameter :: i_site=38,j_site=65 !MMSF
#      integer,parameter :: i_site = 82,j_site =76 !Hyytiala
#      integer,parameter :: i_site=24,j_site=65 !Vaira
#      integer,parameter :: i_site=43,j_site=66 !NYC 
#      integer, parameter :: i_site=51,j_site=44 !TNF
#      integer, parameter :: i_site=87,j_site=46 !Mpala
# 5 deg long x 4 deg lat Resolution   
#      integer,parameter :: i_site=19,j_site=33 !MMSF
#      integer,parameter :: i_site = 40,j_site =39 !Hyytiala
#      integer,parameter :: i_site=12,j_site=33 !Vaira
#      integer, parameter :: i_site=26,j_site=23 !TNF
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
