#-----------------------------------------------------------------------
# Include source list files from subdirectories

include $(SRC_DIR)/prdc/crystal/sources.mk
include $(SRC_DIR)/prdc/field/sources.mk
include $(SRC_DIR)/prdc/fieldIo/sources.mk
include $(SRC_DIR)/prdc/environment/sources.mk

#-----------------------------------------------------------------------
# Object and dependency file lists for src/prdc

prdc_OBJS= \
  $(prdc_crystal_OBJS) \
  $(prdc_field_OBJS) \
  $(prdc_fieldIo_OBJS) \
  $(prdc_environment_OBJS)

prdc_DEPS= \
  $(prdc_crystal_DEPS) \
  $(prdc_field_DEPS) \
  $(prdc_fieldIo_DEPS) \
  $(prdc_environment_DEPS)

#-----------------------------------------------------------------------
# Path and rule for the prdc/libprdc.a library

prdc_LIB=$(BLD_DIR)/prdc/libprdc.a

$(prdc_LIB): $(prdc_OBJS)
	$(AR) rcs $(prdc_LIB) $(prdc_OBJS)

