# Include source files lists from subdirectories
include $(SRC_DIR)/rp/scft/iterator/sources.mk
include $(SRC_DIR)/rp/scft/sweep/sources.mk

# Combine object and dependency file lists from subdirectories
rp_scft_OBJS= \
    $(rp_scft_iterator_OBJS) \
    $(rp_scft_sweep_OBJS) 

rp_scft_DEPS= \
    $(rp_scft_iterator_DEPS) \
    $(rp_scft_sweep_DEPS) 

# Add files associated with the ScftThermo class template
ifdef PSCF_CPP
  rp_scft_OBJS+= $(BLD_DIR)/rp/scft/ScftThermo.o
  rp_scft_DEPS+= $(BLD_DIR)/rp/scft/ScftThermo.d
endif

ifdef PSCF_CUDA
  rp_scft_OBJS+= $(BLD_DIR)/rp/scft/ScftThermo.ou
  rp_scft_DEPS+= $(BLD_DIR)/rp/scft/ScftThermo.du
endif

