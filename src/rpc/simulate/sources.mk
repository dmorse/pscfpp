# Include source files lists from subdirectories
include $(SRC_DIR)/rpc/simulate/mcmove/sources.mk
include $(SRC_DIR)/rpc/simulate/bdstep/sources.mk
include $(SRC_DIR)/rpc/simulate/perturbation/sources.mk
include $(SRC_DIR)/rpc/simulate/analyzer/sources.mk
include $(SRC_DIR)/rpc/simulate/trajectory/sources.mk

rpc_simulate_= \
  $(rpc_simulate_mcmove_) \
  $(rpc_simulate_bdstep_) \
  $(rpc_simulate_perturbation_) \
  $(rpc_simulate_analyzer_) \
  $(rpc_trajectory_) \
  rpc/simulate/SimState.cpp \
  rpc/simulate/Simulator.cpp \
  rpc/simulate/SimulatorFactory.cpp 
  
rpc_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_simulate_:.cpp=.o))

