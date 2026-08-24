#------------------------------------------------------------------------
# src/prdc/makefile
#
# The repository copy of this namespace level makefile is in the src/prdc
# directory. This file is copied verbatim by the configure script to
# create an identical file in the bld/prdc directory.
#------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../config.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#-----------------------------------------------------------------------
# Main targets

all: $(prdc_LIB)

clean:
	rm -f $(prdc_OBJS)
	rm -f $(prdc_DEPS)
	rm -f $(prdc_LIB)
	rm -f *.o */*.o */*/*.o
	rm -f *.d */*.d */*/*.d
	rm -f field/cuda/*.ou field/cuda/*.du 
	cd tests; $(MAKE) clean

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_DEPS)
