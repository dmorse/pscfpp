# Include source files lists from subdirectories
include $(SRC_DIR)/rpc/scft/iterator/sources.mk
include $(SRC_DIR)/rpc/scft/sweep/sources.mk

rpc_scft_= \
  $(rpc_scft_iterator_) \
  $(rpc_scft_sweep_) \
  
rpc_scft_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_scft_:.cpp=.o))

