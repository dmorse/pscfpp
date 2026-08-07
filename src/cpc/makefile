#------------------------------------------------------------------------
# src/cpc/makefile
#
# The repository copy of this namespace level makefile is in the src/cpc
# directory. This is copied verbatim by the configure script to create 
# an identical file in the bld/cpc directory.
#------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/cpc/sources.mk
include $(SRC_DIR)/cpc/patterns.mk

#-----------------------------------------------------------------------
# Paths to pscf_cpc program files

# Base path to *.o and *.d files for pscf_cpc
PSCF_CPC=$(BLD_DIR)/cpc/pscf_cpc

# Path to pscf_cpc executable file
PSCF_CPC_EXE=$(BIN_DIR)/pscf_cpc

#-----------------------------------------------------------------------
# Main targets 

all: $(cpc_OBJS) $(cpc_LIB)

clean:
	rm -f $(cpc_OBJS) $(cpc_OBJS:.o=.d)
	rm -f $(cpc_LIB)
	rm -f $(PSCF_CPC).o $(PSCF_CPC).d
	rm -f *.o */*.o */*/*.o
	rm -f *.d */*.d */*/*.d
	cd tests; $(MAKE) clean

# Rule for pscf_cpc executable file
$(PSCF_CPC_EXE): $(PSCF_CPC).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_CPC_EXE) $(PSCF_CPC).o $(LIBS)

# Short phony target for executable file (for convenience)
pscf_cpc:
	$(MAKE) $(PSCF_CPC_EXE)

# Note: In the rule for pscf_cpc, we include the list $(PSCF_LIBS) of 
# PSCF static libraries as prerequisites, but link to the list $(LIBS) 
# that includes all relevant external libraries.

.PHONY: all clean pscf_cpc

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(prdc_DEPS)
-include $(cpc_OBJS:.o=.d)
-include $(PSCF_CPC).d 
