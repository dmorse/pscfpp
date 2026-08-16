#-----------------------------------------------------------------------
# Include source list files from subdirectories
include $(SRC_DIR)/r1d/field/sources.mk
include $(SRC_DIR)/r1d/solvers/sources.mk
include $(SRC_DIR)/r1d/iterator/sources.mk
include $(SRC_DIR)/r1d/sweep/sources.mk
include $(SRC_DIR)/r1d/system/sources.mk

# List of source files in src/r1d
r1d_OBJS=\
  $(r1d_field_OBJS) \
  $(r1d_solvers_OBJS) \
  $(r1d_iterator_OBJS) \
  $(r1d_sweep_OBJS) \
  $(r1d_system_OBJS) 

# List of dependency files in src/r1d
r1d_DEPS=\
  $(r1d_field_DEPS) \
  $(r1d_solvers_DEPS) \
  $(r1d_iterator_DEPS) \
  $(r1d_sweep_DEPS) \
  $(r1d_system_DEPS) 

#-----------------------------------------------------------------------
# Path and rule for the r1d/libr1d.a library 

r1d_LIB=$(BLD_DIR)/r1d/libr1d.a

$(r1d_LIB): $(r1d_OBJS)
	$(AR) rcs $(r1d_LIB) $(r1d_OBJS)

#-----------------------------------------------------------------------
