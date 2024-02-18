# Include source files lists from subdirectories
include $(SRC_DIR)/pspc/simulate/mcmove/sources.mk
include $(SRC_DIR)/pspc/simulate/bdstep/sources.mk
include $(SRC_DIR)/pspc/simulate/perturbation/sources.mk
include $(SRC_DIR)/pspc/simulate/analyzer/sources.mk
include $(SRC_DIR)/pspc/simulate/trajectory/sources.mk

pspc_simulate_= \
  $(pspc_simulate_mcmove_) \
  $(pspc_simulate_bdstep_) \
  $(pspc_simulate_perturbation_) \
  $(pspc_simulate_analyzer_) \
  $(pspc_trajectory_) \
  pspc/simulate/Simulator.cpp \
  pspc/simulate/SimulatorFactory.cpp 
  
pspc_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_simulate_:.cpp=.o))

