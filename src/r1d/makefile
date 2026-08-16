#------------------------------------------------------------------------
# src/r1d/makefile
#
# The repository copy of this namespace level makefile is in the src/r1d
# directory. This is copied verbatim by the configure script to create 
# an identical file in the bld/r1d directory.
#------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/r1d/sources.mk
include $(SRC_DIR)/r1d/patterns.mk

#-----------------------------------------------------------------------
# Paths to pscf_r1d program files

# Base path for pscf_r1d.o and pscf_r1d.d intermediate files 
PSCF_R1D=$(BLD_DIR)/r1d/pscf_r1d

# Path to pscf_r1d executable file
PSCF_R1D_EXE=$(BIN_DIR)/pscf_r1d

#-----------------------------------------------------------------------
# Main targets 

all: $(PSCF_R1D_EXE)

clean:
	rm -f $(r1d_OBJS)
	rm -f $(r1d_DEPS)
	rm -f $(r1d_LIB)
	rm -f $(PSCF_R1D).o $(PSCF_R1D).d
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

# Executable pscf_r1d file (linking)
$(PSCF_R1D_EXE): $(PSCF_R1D).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_R1D_EXE) $(PSCF_R1D).o $(LIBS) 

# Note: In the rule for pcsf_r1d, we include the list $(PSCF_LIBS) of 
# PSCF-specific static libraries as prerequisiites, but link to the list 
# $(LIBS) that also includes relevant external libraries.

.PHONY: all clean pscf_r1d

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(r1d_DEPS)
-include $(PSCF_R1D).d
