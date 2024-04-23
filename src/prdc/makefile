#------------------------------------------------------------------------
# src/prdc/makefile
#
# The repository copy of this namespace level makefile is in the src/prdc
# directory. This is copied verbatim by the setup to create the copy in
# in the bld/prdc directory.
#------------------------------------------------------------------------
# Include makefiles

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/prdc/include.mk

#-----------------------------------------------------------------------
# Main targets 
 
.PHONY: all list clean veryclean

all: $(prdc_OBJS) $(prdc_LIB) 

list:
	@echo $(prdc_OBJS)

clean:
	rm -f $(prdc_OBJS) $(prdc_OBJS:.o=.d)
	rm -f $(prdc_LIB)
	rm -f lib*.a
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
