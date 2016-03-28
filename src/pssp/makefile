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
include $(BLD_DIR)/pscf/config.mk
include $(BLD_DIR)/pssp/config.mk
include $(SRC_DIR)/pssp/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pssp/sources.mk

#-----------------------------------------------------------------------
# Variable definition

PSCF_PS=$(BLD_DIR)/pssp/pscf_ps
#-----------------------------------------------------------------------
# Main targets 

all: $(pssp_OBJS) $(pssp_LIB)

clean:
	rm -f $(pssp_OBJS) $(pssp_OBJS:.o=.d)
	rm -f $(PSCF_PS).o $(PSCF_PS).d
	rm -f $(pssp_LIB)
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 
	-rm -f lib*.a


# Executable target

$(PSCF_PS_EXE): $(PSCF_PS).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PS_EXE) $(PSCF_PS).o $(LIBS) $(PSCF_GSL_LIB)

# Short name for executable target (for convenience)
pscf_ps:
	$(MAKE) $(PSCF_PS_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pssp_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PS).d 
