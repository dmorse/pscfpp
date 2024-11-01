# Include source files lists from subdirectories
include $(SRC_DIR)/rpc/fts/compressor/sources.mk
include $(SRC_DIR)/rpc/fts/montecarlo/sources.mk
include $(SRC_DIR)/rpc/fts/brownian/sources.mk
include $(SRC_DIR)/rpc/fts/perturbation/sources.mk
include $(SRC_DIR)/rpc/fts/ramp/sources.mk
include $(SRC_DIR)/rpc/fts/analyzer/sources.mk
include $(SRC_DIR)/rpc/fts/trajectory/sources.mk

rpc_fts_= \
  $(rpc_fts_compressor_) \
  $(rpc_fts_montecarlo_) \
  $(rpc_fts_brownian_) \
  $(rpc_fts_perturbation_) \
  $(rpc_fts_ramp_) \
  $(rpc_fts_analyzer_) \
  $(rpc_fts_trajectory_) \
  rpc/fts/SimState.cpp \
  rpc/fts/Simulator.cpp \
  rpc/fts/SimulatorFactory.cpp 
  
rpc_fts_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_:.cpp=.o))

