#------------------------------------------------------------------------
# src/r1d/makefile
#
# The repository copy of this namespace level makefile is in the src/r1d
# directory. This is copied verbatim by the configure script to create 
# the identical copy in the bld/r1d directory.
#------------------------------------------------------------------------
# Include makefile fragments

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/r1d/include.mk

#-----------------------------------------------------------------------
# Paths to pscf_1d program files

# Base path for pscf_1d.o and pscf_1d.d intermediate files 
PSCF_1D=$(BLD_DIR)/r1d/pscf_1d

# Path to pscf_1d executable file
PSCF_1D_EXE=$(BIN_DIR)/pscf_1d

#-----------------------------------------------------------------------
# Main targets 

all: $(r1d_OBJS) $(r1d_LIB) $(PSCF_1D_EXE)

clean:
	rm -f $(r1d_OBJS) $(r1d_OBJS:.o=.d)
	rm -f $(r1d_LIB)
	rm -f $(PSCF_1D).o $(PSCF_1D).d
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean

# Executable pscf_1d file (linking)
$(PSCF_1D_EXE): $(PSCF_1D).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_1D_EXE) $(PSCF_1D).o $(LIBS) 

# Note: In the rule for pcsf_1d, we include the list $(PSCF_LIBS) of 
# PSCF-specific static libraries as dependencies, but link to the list 
# $(LIBS) of all libraries, which includes relevant external libraries.

# Short phony target for pscf_1d executable file (for convenience)
pscf_1d:
	$(MAKE) $(PSCF_1D_EXE)

.PHONY: all clean veryclean pscf_1d

#-----------------------------------------------------------------------
# Include dependency files

-include $(r1d_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_1D).d 
