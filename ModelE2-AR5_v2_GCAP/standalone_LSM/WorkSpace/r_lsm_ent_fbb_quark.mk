# rundeck for compiling stand-alone version of GISS LSM
# coupled to Ent with FBB photosynthesis

# components to use:
COMPONENTS = giss_LSM_standalone  Ent  giss_LSM    shared

# specific options for components:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES
OPTS_shared =
OPTS_giss_LSM_standalone = USE_ENT=YES

INPUT_FILES = \
fort.951=$(GCMSEARCHPATH)/lsm_data_dump_fort.951 \
fort.952=$(GCMSEARCHPATH)/lsm_data_dump_fort.952 \
fort.953=$(GCMSEARCHPATH)/lsm_data_dump_fort.953 \
SOIL=$(GCMSEARCHPATH)/S4X50093.ext \
TOP_INDEX=$(GCMSEARCHPATH)/top_index_72x46.ij.ext \
VEG=$(GCMSEARCHPATH)/V72X46.1.cor2_no_crops.ext \
CROPS=$(GCMSEARCHPATH)/CROPS_72X46N.cor4.ext \
soil_textures=$(GCMSEARCHPATH)/soil_textures_top30cm \
SOILCARB_global=$(GCMSEARCHPATH)/soilcarb_top30cm_nmaps_4x5bin.dat


# run settings (one value per line, single quotation marks for strings)
RUN_PARAMETERS = \
#&&PARAMETERS \n\
#&&END_PARAMETERS

# the following are settings for "quark"
RUN_PARAMETERS = \
&&PARAMETERS \n\
itime1=48 \n\
basefold_out='/a11/tmp/GSWP_1983' \n\
neutral_drag_file='/a11/tmp/GSWP_1983/CDN5X4.GISS' \n\
&&END_PARAMETERS
