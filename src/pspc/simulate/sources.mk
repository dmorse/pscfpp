# Include source files lists from subdirectories
include $(SRC_DIR)/pspc/simulate/analyzer/sources.mk
include $(SRC_DIR)/pspc/simulate/mcmove/sources.mk
include $(SRC_DIR)/pspc/simulate/trajectory/sources.mk

pspc_simulate_= \
  $(pspc_simulate_analyzer_) \
  $(pspc_simulate_mcmove_) \
  $(pspc_trajectory_) \
  pspc/simulate/McSimulator.cpp \
  pspc/simulate/McState.cpp \
  pspc/simulate/Simulator.cpp 
  
pspc_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_simulate_:.cpp=.o))

