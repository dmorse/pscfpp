include $(SRC_DIR)/rpc/field/sources.mk
include $(SRC_DIR)/rpc/solvers/sources.mk
include $(SRC_DIR)/rpc/iterator/sources.mk
include $(SRC_DIR)/rpc/sweep/sources.mk
include $(SRC_DIR)/rpc/compressor/sources.mk
include $(SRC_DIR)/rpc/simulate/sources.mk

rpc_= \
  $(rpc_field_) \
  $(rpc_solvers_) \
  $(rpc_iterator_) \
  $(rpc_sweep_) \
  $(rpc_compressor_) \
  $(rpc_simulate_) \
  rpc/System.cpp 

rpc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_:.cpp=.o))

$(rpc_LIB): $(rpc_OBJS)
	$(AR) rcs $(rpc_LIB) $(rpc_OBJS)

