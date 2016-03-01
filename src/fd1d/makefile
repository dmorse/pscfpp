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
include $(BLD_DIR)/fd1d/config.mk
include $(SRC_DIR)/fd1d/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/fd1d/sources.mk

#-----------------------------------------------------------------------
# Variable definition

PSCF_FD=$(BLD_DIR)/fd1d/pscf_fd
#-----------------------------------------------------------------------
# Main targets 

all: $(fd1d_OBJS) $(fd1d_LIB) $(PSCF_FD_EXE)

clean:
	rm -f $(fd1d_OBJS) $(fd1d_OBJS:.o=.d) $(fd1d_LIB)
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 
	-rm -f lib*.a


# Executable target

$(PSCF_FD_EXE): $(PSCF_FD).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_FD_EXE) $(PSCF_FD).o $(LIBS) $(PSCF_GSL_LIB)

# Short name for executable target (for convenience)
pscf_fd:
	$(MAKE) $(PSCF_FD_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(fd1d_OBJS:.o=.d)
-include $(PSCF_FD).d 
