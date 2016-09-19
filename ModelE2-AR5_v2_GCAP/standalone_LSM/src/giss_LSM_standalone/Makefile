# makefile for the component giss_LSM_standalone
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF LSM_DRV

#FSRCS = lsm_standalone.f drv_gswp_force.f domain_decomp.f
FSRCS = domain_decomp.f
CSRCS = sys_usage.c # mmap_utils.c

# New version of lsm_standalone
ifndef LSM_DRV
	FSRCS += lsm_standalone.f drv_gswp_force.f
else
	ifeq ($(LSM_DRV),GLOBAL)
		FSRCS += lsm_standalone_nk.f drv_met_gswp_nk.f \
		drv_veg.f lsm_phys_util.f
	else 
		ifeq ($(LSM_DRV),SINGLE)
        		FSRCS += lsm_standalone_nk.f drv_met_single_nk.f \
			drv_veg.f lsm_phys_util.f
			CPPFLAGS += -DLSM_DRV_SINGLE
		endif
	endif
endif

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
#sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

CPPFLAGS += -DUSE_ENT
CPPFLAGS += -I$(GISSCLIM_DIR)/include
