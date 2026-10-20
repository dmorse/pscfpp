#------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 
 
all: $(prdc_field_OBJS)

clean:
	rm -f $(prdc_field_OBJS)
	rm -f $(prdc_field_DEPS)
	rm -f *.o */*.o
	rm -f *.d */*.d 

.PHONY: all clean
#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_field_DEPS)
