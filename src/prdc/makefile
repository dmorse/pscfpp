#------------------------------------------------------------------------
# src/prdc/makefile
#
# The repository copy of this namespace level makefile is in the src/prdc
# directory. This file is copied verbatim by the configure script to 
# create the identical copy in the bld/prdc directory.
#------------------------------------------------------------------------
# Include makefile fragments

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/prdc/include.mk

#-----------------------------------------------------------------------
# Main targets 
 
all: $(prdc_OBJS) $(prdc_LIB) 

clean:
	rm -f $(prdc_OBJS) $(prdc_OBJS:.o=.d)
	rm -f $(prdc_LIB)
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean

.PHONY: all clean veryclean

#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
