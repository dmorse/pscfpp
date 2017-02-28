#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/pscf/config.mk
include $(BLD_DIR)/pssp/config.mk
include $(SRC_DIR)/pssp/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pssp/sources.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pssp_solvers_OBJS) 

clean:
	rm -f $(pssp_solvers_OBJS) $(pssp_solvers_OBJS:.o=.d)

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 

#-----------------------------------------------------------------------
# Include dependency files

-include $(pssp_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_FD).d 
