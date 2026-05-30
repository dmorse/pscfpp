#------------------------------------------------------------------------
# src/cpg/makefile
#
# The repository copy of this namespace level makefile is in the src/cpg
# directory. This is copied verbatim by the configure script to create 
# the identical copy in in the bld/cpg directory.
#------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/cpg/sources.mk
include $(SRC_DIR)/cpg/patterns.mk

#-----------------------------------------------------------------------
# Paths to pscf_cpg program files

# Base path for *.o and *.d files for pscf_cpg
PSCF_CPG=$(BLD_DIR)/cpg/pscf_cpg

# Path to pscf_cpg executable
PSCF_CPG_EXE=$(BIN_DIR)/pscf_cpg

#-----------------------------------------------------------------------
# Main targets 

all: $(cpg_OBJS) $(cpg_LIB)

clean:
	rm -f $(cpg_OBJS) $(cpg_OBJS:.o=.d)
	rm -f $(cpg_LIB)
	rm -f $(PSCF_CPG).o $(PSCF_CPG).d
	rm -f *.o */*.o */*/*.o
	rm -f *.d */*.d */*/*.d
	rm -f lib*.a 
	cd tests; $(MAKE) clean

# Main pscf_cpg executable file
$(PSCF_CPG_EXE): $(PSCF_CPG).o $(PSCF_LIBS)
	$(NVXX) $(LDFLAGS) -o $(PSCF_CPG_EXE) $(PSCF_CPG).o $(LIBS)

# Short phony target for executable file (for convenience)
pscf_cpg:
	$(MAKE) $(PSCF_CPG_EXE)

# Note: In the linking rule for pscf_cpg, we include the list $(PSCF_LIBS) 
# of PSCF-specific libraries as pre-requisistes but link to the full list 
# that includes relevant external libraries.

.PHONY: all clean pscf_cpg

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(cpg_OBJS:.o=.d)
-include $(PSCF_CPG).d 
