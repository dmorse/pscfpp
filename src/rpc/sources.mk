#-----------------------------------------------------------------------
# Source files in src/rpc and corresponding object files

include $(SRC_DIR)/rpc/field/sources.mk
include $(SRC_DIR)/rpc/solvers/sources.mk
include $(SRC_DIR)/rpc/iterator/sources.mk
include $(SRC_DIR)/rpc/sweep/sources.mk
include $(SRC_DIR)/rpc/compressor/sources.mk
include $(SRC_DIR)/rpc/simulate/sources.mk

# List of source files in src/rpc
rpc_= \
  $(rpc_field_) \
  $(rpc_solvers_) \
  $(rpc_iterator_) \
  $(rpc_sweep_) \
  $(rpc_compressor_) \
  $(rpc_simulate_) \
  rpc/System.cpp 

# List of object file targets
rpc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_:.cpp=.o))

#-----------------------------------------------------------------------
# Path and target rule for the librpc.a library 

rpc_LIBNAME=rpc
rpc_LIB=$(BLD_DIR)/rpc/lib$(rpc_LIBNAME).a

$(rpc_LIB): $(rpc_OBJS)
	$(AR) rcs $(rpc_LIB) $(rpc_OBJS)

