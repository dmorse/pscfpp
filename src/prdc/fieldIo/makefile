#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(prdc_field_OBJS) 

clean:
	rm -f $(prdc_field_OBJS) $(prdc_field_DEPS) 

.PHONY: all clean
#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_field_DEPS)
