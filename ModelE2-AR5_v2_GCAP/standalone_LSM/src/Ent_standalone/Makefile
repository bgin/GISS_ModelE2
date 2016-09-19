# makefile for component Ent_standalone
.SUFFIXES:

SUPPORTED_OPTIONS = PFT_MODEL MIXED_CANOPY_OPT SOILCARB_SITE_OPT

FSRCS = ent_forcings.f  ent_prog.f

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
#sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk
#CPPFLAGS for ABSOFT and miscellaneous Ent_standalaone features.
#CPPFLAGS += -DCOMPILER_ABSOFT 
CPPFLAGS += -DSOILCARB_SITE  -DENT_STANDALONE_DIAG -DFLUXNETINIT

ifeq ($(PFT_MODEL),ENT)
CPPFLAGS += -DPFT_MODEL_ENT
endif

ifeq ($(MIXED_CANOPY_OPT),YES)
  CPPFLAGS += -DMIXED_CANOPY
endif

ifeq ($(SOILCARB_SITE_OPT),YES)
  CPPFLAGS += -DSOILCARB_SITE
endif

