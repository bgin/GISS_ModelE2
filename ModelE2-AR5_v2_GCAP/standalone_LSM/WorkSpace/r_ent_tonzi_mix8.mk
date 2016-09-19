# rundeck for compiling stand-alone version of GISS LSM
# coupled to Ent with FBB photosynthesis
OBJ_LIST=main_ent

# components to use:  
COMPONENTS = Ent_standalone Ent shared
# specific options for components:  
OPTS_Ent = PS_MODEL=FBB PFT_MODEL=GISS FLUXNET=YES ENT_STANDALONE_DIAG=YES MIXED_CANOPY_OPT=YES
OPTS_Ent_standalone = PFT_MODEL=GISS MIXED_CANOPY_OPT=YES ENT_STANDALONE_DIAG=YES
#OPTS_giss_LSM = USE_ENT=YES
OPTS_shared =  
#OPTS_giss_LSM_standalone = USE_ENT=YES

DATADIR=/raid8/Entdata/Tonzi/Siteforcings_6layer
DATADIR_veg=/raid8/Entdata/Tonzi/Siteforcings_mixedveg8

INPUT_FILES = \
ENTSTRUCT=$(DATADIR_veg)/Entstruct_tonzi8.csv \
VEG=$(DATADIR_veg)/VEGCOV.Tonzi \
CROPS=$(DATADIR)/CROPS.Tonzi \
soil_textures=$(DATADIR)/soil_textures_top30cm.Tonzi \
SOILCARB_site=$(DATADIR)/SOILCARB_site_Tonzi.csv \
CH=$(DATADIR)/CH.Tonzi \
DVIS=$(DATADIR)/DVIS.Tonzi \
MP1=$(DATADIR)/MP1.Tonzi \
MP2=$(DATADIR)/MP2.Tonzi \
MP3=$(DATADIR)/MP3.Tonzi \
MP4=$(DATADIR)/MP4.Tonzi \
MP5=$(DATADIR)/MP5.Tonzi \
MP6=$(DATADIR)/MP6.Tonzi \
PREC=$(DATADIR)/PREC.Tonzi \
PS=$(DATADIR)/PS.Tonzi \
QFOL=$(DATADIR)/QFOL.Tonzi \
QLAT=$(DATADIR)/QLAT.Tonzi \
QSEN=$(DATADIR)/QSEN.Tonzi \
SAT=$(DATADIR)/SurfAirTC.Tonzi \
SIC1=$(DATADIR)/SIC1.Tonzi \
SIC2=$(DATADIR)/SIC2.Tonzi \
SIC3=$(DATADIR)/SIC3.Tonzi \
SIC4=$(DATADIR)/SIC4.Tonzi \
SIC5=$(DATADIR)/SIC5.Tonzi \
SIC6=$(DATADIR)/SIC6.Tonzi \
TP=$(DATADIR)/TcanC.Tonzi \
US=$(DATADIR)/US.Tonzi \
VIS=$(DATADIR)/VIS.Tonzi \
VS=$(DATADIR)/VS.Tonzi \
ZEN=$(DATADIR)/COSZEN.Tonzi \
SOILT1=$(DATADIR)/SOILT1.Tonzi \
SOILT2=$(DATADIR)/SOILT2.Tonzi \
SOILT3=$(DATADIR)/SOILT3.Tonzi \
SOILT4=$(DATADIR)/SOILT4.Tonzi \
SOILT5=$(DATADIR)/SOILT5.Tonzi \
SOILT6=$(DATADIR)/SOILT6.Tonzi \
W1=$(DATADIR)/W1.Tonzi \
W2=$(DATADIR)/W2.Tonzi \
W3=$(DATADIR)/W3.Tonzi \
W4=$(DATADIR)/W4.Tonzi \
W5=$(DATADIR)/W5.Tonzi \
W6=$(DATADIR)/W6.Tonzi \
LAIMIXED=$(DATADIR_veg)/LAIMIXED8.Tonzi \
HEIGHTMIXED=$(DATADIR_veg)/HEIGHTMIXED8.Tonzi \

RUN_PARAMETERS = \
&input_parameters \
year=2002,jday=1,jday2=365,year2=2002,dt=1800.0, \
IM=72,JM=46,I0=13,I1=13,J0=33,J1=33,I0f=13,I1f=13,J0f=33,J1f=33, \
force_VEG=.true., \
do_soilinit=.true., \
do_phenology_activegrowth=.false., \
do_structuralgrowth=.false., \
do_frost_hardiness=.true.,\
do_spinup=.true.\
/
