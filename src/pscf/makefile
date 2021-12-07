#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the corresponding subdirectory of the bld/ directories by 
# the setup script. Only the copy in the src/ directory is stored in the 
# repository.
#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pscf_OBJS) $(pscf_LIB)

clean:
	rm -f $(pscf_OBJS) $(pscf_OBJS:.o=.d) $(pscf_LIB)
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f *.o */*.o
	-rm -f *.d */*.d
	-rm -f lib*.a
ifeq ($(BLD_DIR),$(SRC_DIR))
	-rm -f paths.h
endif


#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_OBJS:.o=.d)
