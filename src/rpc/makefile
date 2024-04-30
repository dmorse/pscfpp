#------------------------------------------------------------------------
# src/rpc/makefile
#
# The repository copy of this namespace level makefile is in the src/rpc
# directory. This is copied verbatim by the configure script to create 
# the identical copy in the bld/rpc directory.
#------------------------------------------------------------------------
# Include makefile fragments

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpc/include.mk

#-----------------------------------------------------------------------
# Paths to pscf_pc program files

# Base path to *.o and *.d files for pscf_pc
PSCF_PC=$(BLD_DIR)/rpc/pscf_pc

# Path to pscf_pc executable file
PSCF_PC_EXE=$(BIN_DIR)/pscf_pc

#-----------------------------------------------------------------------
# Main targets 

all: $(rpc_OBJS) $(rpc_LIB) $(PSCF_PC_EXE)

clean:
	rm -f $(rpc_OBJS) $(rpc_OBJS:.o=.d)
	rm -f $(rpc_LIB)
	rm -f $(PSCF_PC).o $(PSCF_PC).d
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean

# Rule for pscf_pc executable file
$(PSCF_PC_EXE): $(PSCF_PC).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PC_EXE) $(PSCF_PC).o $(LIBS)

# Note: In the rule for pcsf_pc, we include the list $(PSCF_LIBS) of 
# PSCF static libraries as dependencies, but link to the list $(LIBS) 
# of all libraries, which includes relevant external libraries.

# Short phony target for executable file (for convenience)
pscf_pc:
	$(MAKE) $(PSCF_PC_EXE)

.PHONY: all clean veryclean pscf_pc

#-----------------------------------------------------------------------
# Include dependency files

-include $(rpc_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PC).d 
