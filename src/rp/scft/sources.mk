# Include source files lists from subdirectories
include $(SRC_DIR)/rp/scft/iterator/sources.mk
include $(SRC_DIR)/rp/scft/sweep/sources.mk

rp_scft_OBJS=
rp_scft_DEPS=

ifdef PSCF_CPP
  rp_scft_CPP= \
    $(rp_scft_iterator_CPP) \
    $(rp_scft_sweep_CPP) \
    rp/scft/ScftThermo_c.cpp
  rp_scft_OBJS+= \
      $(addprefix $(BLD_DIR)/, $(rp_scft_CPP:.cpp=.o))
  rp_scft_DEPS+= \
      $(addprefix $(BLD_DIR)/, $(rp_scft_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_scft_CU= \
    $(rp_scft_iterator_CU) \
    $(rp_scft_sweep_CU) \
    rp/scft/ScftThermo_u.cu
  rp_scft_OBJS += \
      $(addprefix $(BLD_DIR)/, $(rp_scft_CU:.cu=.o))
  rp_scft_DEPS+= \
      $(addprefix $(BLD_DIR)/, $(rp_scft_CU:.cu=.d))
endif

