# Include source files lists from subdirectories
include $(SRC_DIR)/rp/scft/iterator/sources.mk
include $(SRC_DIR)/rp/scft/sweep/sources.mk

rp_scft_OBJS=

ifdef PSCF_CPP
rp_scft_CPP_= \
  $(rp_scft_iterator_CPP) \
  $(rp_scft_sweep_CPP) \
  rp/scft/ScftThermo_c.cpp
rp_scft_OBJS+= \
     $(addprefix $(BLD_DIR)/, $(rp_scft_CPP:.cpp=.o))
endif

ifdef PSCF_CUDA
rp_scft_CUDA_= \
  $(rp_scft_iterator_CUDA) \
  $(rp_scft_sweep_CUDA) \
  rp/scft/ScftThermo_u.cu
rp_scft_OBJS += \
     $(addprefix $(BLD_DIR)/, $(rp_scft_CUDA:.cu=.o))
endif

