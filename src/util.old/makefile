#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(util_OBJS) $(util_LIB)

clean:
	rm -f $(util_OBJS) $(util_OBJS:.o=.d) $(util_LIB)
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	rm -f */*.o */*/*.o */*/*/*.o
	rm -f */*.d */*/*.d */*/*/*.d
	rm -f lib*.a
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f boundary/Boundary.h
endif

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_OBJS:.o=.d)
