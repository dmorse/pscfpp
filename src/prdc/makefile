#------------------------------------------------------------------------
# src/prdc/makefile
#
# The repository copy of this namespace level makefile is in the src/prdc
# directory. This file is copied verbatim by the configure script to 
# create an identical file in the bld/prdc directory.
#------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 
 
all: $(prdc_OBJS) $(prdc_LIB) 

clean:
	rm -f $(prdc_OBJS) $(prdc_OBJS:.o=.d)
	rm -f $(prdc_LIB)
	rm -f *.o */*.o */*/*.o
	rm -f *.d */*.d */*/*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean

.PHONY: all clean veryclean

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
