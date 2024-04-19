#-----------------------------------------------------------------------
# Source files in src/rpg and corresponding object file targets

include $(SRC_DIR)/rpg/field/sources.mk
include $(SRC_DIR)/rpg/solvers/sources.mk
include $(SRC_DIR)/rpg/iterator/sources.mk
include $(SRC_DIR)/rpg/sweep/sources.mk
include $(SRC_DIR)/rpg/simulate/sources.mk

# List of source files in src/rpg
rpg_= \
  $(rpg_field_) \
  $(rpg_solvers_) \
  $(rpg_iterator_) \
  $(rpg_sweep_) \
  $(rpg_simulate_) \
  rpg/System.cu

# List of object file targets
rpg_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_:.cu=.o))

#-----------------------------------------------------------------------
# Path and makefile target for the librpg.a library file

rpg_LIBNAME=rpg
rpg_LIB=$(BLD_DIR)/rpg/lib$(rpg_LIBNAME).a

$(rpg_LIB): $(rpg_OBJS)
	$(AR) rcs $(rpg_LIB) $(rpg_OBJS)

