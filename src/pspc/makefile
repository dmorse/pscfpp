#------------------------------------------------------------------------
# src/pspc/makefile
#
# The repository copy of this namespace level makefile is in the src/pspc
# directory. This is copied verbatim by the setup script to create the 
# copy in the bld/pspc directory.
#------------------------------------------------------------------------
# Include makefile fragments

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk

#-----------------------------------------------------------------------
# Main pscf_pc program files

# Base path to *.o and *.d files for pscf_pc
PSCF_PC=$(BLD_DIR)/pspc/pscf_pc

# Path to pscf_pc executable file
PSCF_PC_EXE=$(BIN_DIR)/pscf_pc

#-----------------------------------------------------------------------
# Main targets 

.PHONY: all clean veryclean

all: $(pspc_OBJS) $(pspc_LIB) $(PSCF_PC_EXE)

clean:
	rm -f $(pspc_OBJS) $(pspc_OBJS:.o=.d)
	rm -f $(pspc_LIB)
	rm -f $(PSCF_PC).o $(PSCF_PC).d
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f lib*.a

# Target for main pscf_pc executable 
$(PSCF_PC_EXE): $(PSCF_PC).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PC_EXE) $(PSCF_PC).o $(LIBS)

# Note: In the target for pcsf_pc, we include the list $(PSCF_LIBS) of 
# PSCF static libraries as dependencies, but link to the list $(LIBS) 
# of all libraries, which includes relevant external libraries.

# Short name for executable target (for convenience)
pscf_pc:
	$(MAKE) $(PSCF_PC_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pspc_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PC).d 
