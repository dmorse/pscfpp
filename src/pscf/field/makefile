#------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

#-----------------------------------------------------------------------
# Main targets 
 
all: $(pscf_field_OBJS)

clean:
	rm -f $(pscf_field_OBJS)
	rm -f $(pscf_field_DEPS)
	rm -f *.o */*.o
	rm -f *.d */*.d 

.PHONY: all clean
#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_field_DEPS)
