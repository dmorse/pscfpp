#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(prdc_crystal_OBJS) 

clean:
	rm -f $(prdc_crystal_OBJS) $(prdc_crystal_DEPS) 

.PHONY: all clean
#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_crystal_DEPS)
