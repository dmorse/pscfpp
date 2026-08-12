#-----------------------------------------------------------------------
# Source and object file lists for src/cpc

# Include source list files from subdirectories
include $(SRC_DIR)/cpc/solvers/sources.mk
include $(SRC_DIR)/cpc/field/sources.mk
include $(SRC_DIR)/cpc/system/sources.mk
include $(SRC_DIR)/cpc/fts/sources.mk

# List of source files in src/cpc
cpc_= \
  $(cpc_solvers_) \
  $(cpc_field_) \
  $(cpc_system_) \
  $(cpc_fts_)

# List of object file targets
cpc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_:.cpp=.o))

#-----------------------------------------------------------------------
# Path and rule for the cpc/libcpc.a library

cpc_LIB=$(BLD_DIR)/cpc/libcpc.a

$(cpc_LIB): $(cpc_OBJS)
	$(AR) rcs $(cpc_LIB) $(cpc_OBJS)

