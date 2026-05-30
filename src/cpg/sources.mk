#-----------------------------------------------------------------------
# Source and object file lists for src/cpg

# Include source list files from subdirectories
include $(SRC_DIR)/cpg/solvers/sources.mk
include $(SRC_DIR)/cpg/field/sources.mk
#include $(SRC_DIR)/cpg/system/sources.mk
#include $(SRC_DIR)/cpg/fts/sources.mk

# List of source files in src/cpg
cpg_= \
  $(cpg_solvers_) \
  $(cpg_field_)

# List of object file targets
cpg_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpg_:.cu=.o))

#-----------------------------------------------------------------------
# Path and rule for the cpg/libcpg.a library

cpg_LIB=$(BLD_DIR)/cpg/libcpg.a

$(cpg_LIB): $(cpg_OBJS)
	$(AR) rcs $(cpg_LIB) $(cpg_OBJS)

