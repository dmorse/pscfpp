#-----------------------------------------------------------------------
# Source and object file lists for src/rpc 

# Include source list files from subdirectories
include $(SRC_DIR)/rpc/environment/sources.mk
include $(SRC_DIR)/rpc/field/sources.mk
include $(SRC_DIR)/rpc/solvers/sources.mk
include $(SRC_DIR)/rpc/scft/sources.mk
include $(SRC_DIR)/rpc/fts/sources.mk
include $(SRC_DIR)/rpc/system/sources.mk

# List of source files in src/rpc
rpc_= \
  $(rpc_environment_) \
  $(rpc_field_) \
  $(rpc_solvers_) \
  $(rpc_scft_) \
  $(rpc_fts_) \
  $(rpc_system_) 

# List of object file targets
rpc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_:.cpp=.o))

#-----------------------------------------------------------------------
# Path and rule for the rpc/librpc.a library 

rpc_LIB=$(BLD_DIR)/rpc/librpc.a

$(rpc_LIB): $(rpc_OBJS)
	$(AR) rcs $(rpc_LIB) $(rpc_OBJS)

