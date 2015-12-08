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
include $(BLD_DIR)/chem/config.mk
include $(SRC_DIR)/chem/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/chem/sources.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(chem_OBJS) $(chem_LIB)

clean:
	rm -f $(chem_OBJS) $(chem_OBJS:.o=.d) $(chem_LIB)
	#cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 
	-rm -f lib*.a

#-----------------------------------------------------------------------
# Include dependency files

-include $(chem_OBJS:.o=.d)
