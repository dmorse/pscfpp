# Include source files lists from subdirectories
include $(SRC_DIR)/rp/fts/simulator/sources.mk
include $(SRC_DIR)/rp/fts/compressor/sources.mk
include $(SRC_DIR)/rp/fts/brownian/sources.mk
include $(SRC_DIR)/rp/fts/montecarlo/sources.mk
include $(SRC_DIR)/rp/fts/perturbation/sources.mk
include $(SRC_DIR)/rp/fts/ramp/sources.mk

rp_fts_OBJS=
rp_fts_DEPS=

ifdef PSCF_CPP
  rp_fts_CPP= \
    $(rp_fts_simulator_CPP) \
    $(rp_fts_compressor_CPP) \
    $(rp_fts_brownian_CPP) \
    $(rp_fts_montecarlo_CPP) \
    $(rp_fts_perturbation_CPP) \
    $(rp_fts_ramp_CPP)
  rp_fts_OBJS+= \
      $(addprefix $(BLD_DIR)/, $(rp_fts_CPP:.cpp=.o))
  rp_fts_DEPS+= \
      $(addprefix $(BLD_DIR)/, $(rp_fts_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_fts_CU= \
    $(rp_fts_simulator_CU) \
    $(rp_fts_compressor_CU) \
    $(rp_fts_brownian_CU) \
    $(rp_fts_montecarlo_CU) \
    $(rp_fts_perturbation_CU) \
    $(rp_fts_ramp_CU)
  rp_fts_OBJS += \
      $(addprefix $(BLD_DIR)/, $(rp_fts_CU:.cu=.ou))
  rp_fts_DEPS+= \
      $(addprefix $(BLD_DIR)/, $(rp_fts_CU:.cu=.du))
endif

