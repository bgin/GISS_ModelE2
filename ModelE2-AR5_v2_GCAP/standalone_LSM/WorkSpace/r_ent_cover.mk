# Rundeck for compiling program to generate Ent vegetation cover structure.
# Reads in ascii file description, and generates entcells structure for
# a single mixed canopy.  Generates cohort structure only, does not calculate
# Structure can then be read in by an Ent_standalone run.

OBJ_LIST=main_entcover

# components to use:
COMPONENTS = Entcover_standalone Ent shared

# specific options for components:
OPTS_Ent = PS_MODEL=FBB PFT_MODEL=ENT FLUXNET=YES MIXED_CANOPY_OPT=YES
#OPTS_giss_LSM = USE_ENT=YES
OPTS_shared =

DATADIR = /raid8/Entdata/Tonzi/Siteforcings_mixedveg

INPUT_FILES = \
ENTSTRUCT=$(DATADIR)/Entstruct_tonzi16.csv \

RUN_PARAMETERS = \
&input_parameters \
year = 2005,year2=2005,jday=1, jday2=365, \
IM=72, JM=46, I0=13, I1=13, J0=33, J1=33,I0f=13,I1f=13,J0f=33,J1f=33 \
force_VEG=true, do_soilinit=true, do_spinup=true \
/
