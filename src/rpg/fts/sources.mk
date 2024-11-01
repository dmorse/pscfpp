# Include source files lists from subdirectories
include $(SRC_DIR)/rpg/fts/compressor/sources.mk
include $(SRC_DIR)/rpg/fts/montecarlo/sources.mk
include $(SRC_DIR)/rpg/fts/brownian/sources.mk
include $(SRC_DIR)/rpg/fts/perturbation/sources.mk
include $(SRC_DIR)/rpg/fts/ramp/sources.mk
include $(SRC_DIR)/rpg/fts/trajectory/sources.mk
include $(SRC_DIR)/rpg/fts/analyzer/sources.mk

rpg_fts_= \
  $(rpg_fts_compressor_) \
  $(rpg_fts_mcmove_) \
  $(rpg_fts_bdstep_) \
  $(rpg_fts_perturbation_) \
  $(rpg_fts_ramp_) \
  $(rpg_fts_trajectory_) \
  $(rpg_fts_analyzer_) \
  rpg/fts/SimState.cu \
  rpg/fts/Simulator.cu \
  rpg/fts/SimulatorFactory.cu
  
rpg_fts_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_:.cu=.o))

