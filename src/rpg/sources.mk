#-----------------------------------------------------------------------
# Source files in src/rpg and corresponding object file targets

include $(SRC_DIR)/rpg/environment/sources.mk
include $(SRC_DIR)/rpg/field/sources.mk
include $(SRC_DIR)/rpg/solvers/sources.mk
include $(SRC_DIR)/rpg/scft/sources.mk
include $(SRC_DIR)/rpg/fts/sources.mk
include $(SRC_DIR)/rpg/system/sources.mk

# List of source files in src/rpg
rpg_= \
  $(rpg_environment_) \
  $(rpg_field_) \
  $(rpg_solvers_) \
  $(rpg_scft_) \
  $(rpg_fts_) \
  $(rpg_system_) 

# List of object file targets
rpg_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_:.cu=.o))

#-----------------------------------------------------------------------
# Path and makefile target for the librpg.a library file

rpg_LIBNAME=rpg
rpg_LIB=$(BLD_DIR)/rpg/lib$(rpg_LIBNAME).a

$(rpg_LIB): $(rpg_OBJS)
	$(AR) rcs $(rpg_LIB) $(rpg_OBJS)

