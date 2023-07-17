# Include source files lists from subdirectories
include $(SRC_DIR)/pspc/simulate/analyzer/sources.mk

pspc_simulate_= \
  $(pspc_simulate_analyzer_) 
  
pspc_simulate_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_simulate_))
pspc_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_simulate_:.cpp=.o))

