# Include source files lists from subdirectories
include $(SRC_DIR)/rpg/simulate/mcmove/sources.mk
include $(SRC_DIR)/rpg/simulate/bdstep/sources.mk
include $(SRC_DIR)/rpg/simulate/trajectory/sources.mk
include $(SRC_DIR)/rpg/simulate/analyzer/sources.mk

rpg_simulate_= \
  $(rpg_simulate_mcmove_) \
  $(rpg_simulate_bdstep_) \
  $(rpg_trajectory_) \
  $(rpg_simulate_analyzer_) \
  rpg/simulate/SimState.cu \
  rpg/simulate/Simulator.cu \
  rpg/simulate/SimulatorFactory.cu
  
rpg_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_simulate_:.cu=.o))

