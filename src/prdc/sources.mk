#-----------------------------------------------------------------------
# Include source list files from subdirectories

include $(SRC_DIR)/prdc/environment/sources.mk

#-----------------------------------------------------------------------
# Object and dependency file lists for src/prdc

prdc_OBJS= \
  $(prdc_environment_OBJS)

prdc_DEPS= \
  $(prdc_environment_DEPS)

#-----------------------------------------------------------------------
# Path and rule for the prdc/libprdc.a library

prdc_LIB=$(BLD_DIR)/prdc/libprdc.a

$(prdc_LIB): $(prdc_OBJS)
	$(AR) rcs $(prdc_LIB) $(prdc_OBJS)

