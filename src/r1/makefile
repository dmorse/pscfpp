#------------------------------------------------------------------------
# src/r1/makefile
#
# The repository copy of this namespace level makefile is in the src/r1
# directory. This is copied verbatim by the configure script to create 
# an identical file in the bld/r1 directory.
#------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/r1/sources.mk
include $(SRC_DIR)/r1/patterns.mk

#-----------------------------------------------------------------------
# Paths to pscf_1d program files

# Base path for pscf_1d.o and pscf_1d.d intermediate files 
PSCF_1D=$(BLD_DIR)/r1/pscf_1d

# Path to pscf_1d executable file
PSCF_1D_EXE=$(BIN_DIR)/pscf_1d

#-----------------------------------------------------------------------
# Main targets 

all: $(PSCF_1D_EXE)

clean:
	rm -f $(r1_OBJS)
	rm -f $(r1_DEPS)
	rm -f $(r1_LIB)
	rm -f $(PSCF_1D).o $(PSCF_1D).d
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

# Executable pscf_1d file (linking)
$(PSCF_1D_EXE): $(PSCF_1D).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_1D_EXE) $(PSCF_1D).o $(LIBS) 

# Note: In the rule for pcsf_r1, we include the list $(PSCF_LIBS) of 
# PSCF-specific libraries as prerequisiites, but link to the list 
# $(LIBS) that also includes relevant external libraries.

.PHONY: all clean pscf_1d

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(r1_DEPS)
-include $(PSCF_1D).d
