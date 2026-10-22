#-----------------------------------------------------------------------
# Include source list files from subdirectories
include $(SRC_DIR)/r1/field/sources.mk
include $(SRC_DIR)/r1/solvers/sources.mk
include $(SRC_DIR)/r1/iterator/sources.mk
include $(SRC_DIR)/r1/sweep/sources.mk
include $(SRC_DIR)/r1/system/sources.mk

# List of source files in src/r1
r1_OBJS=\
  $(r1_field_OBJS) \
  $(r1_solvers_OBJS) \
  $(r1_iterator_OBJS) \
  $(r1_sweep_OBJS) \
  $(r1_system_OBJS) 

# List of dependency files in src/r1
r1_DEPS=\
  $(r1_field_DEPS) \
  $(r1_solvers_DEPS) \
  $(r1_iterator_DEPS) \
  $(r1_sweep_DEPS) \
  $(r1_system_DEPS) 

#-----------------------------------------------------------------------
# Path and rule for the r1/libr1.a library 

r1_LIB=$(BLD_DIR)/r1/libr1.a

$(r1_LIB): $(r1_OBJS)
	$(AR) rcs $(r1_LIB) $(r1_OBJS)

#-----------------------------------------------------------------------
