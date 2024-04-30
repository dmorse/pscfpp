#------------------------------------------------------------------------
# src/pscf/makefile
#
# The repository copy of this namespace level makefile is in the src/pscf
# directory. This is copied verbatim by the configure script to create the 
# identical copy in the bld/pscf directory.
#------------------------------------------------------------------------
# Include makefile fragments

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pscf_OBJS) $(pscf_LIB)

clean:
	rm -f $(pscf_OBJS) $(pscf_OBJS:.o=.d) 
	rm -f $(pscf_LIB)
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
ifeq ($(BLD_DIR),$(SRC_DIR))
	-rm -f paths.h
endif

.PHONY: all clean veryclean

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
