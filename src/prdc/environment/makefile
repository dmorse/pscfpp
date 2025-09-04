#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(prdc_environment_OBJS) 

clean:
	rm -f $(prdc_environment_OBJS) $(prdc_environment_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_environment_OBJS:.o=.d)
