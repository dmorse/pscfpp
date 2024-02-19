#------------------------------------------------------------------------
# src/r1d/makefile
#
# The repository copy of this namespace level makefile is in the src/r1d
# directory. This is copied verbatim by the setup script to create the 
# copy in the bld/r1d directory.
#------------------------------------------------------------------------
# Include makefile fragments

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/r1d/include.mk

#-----------------------------------------------------------------------
# Paths to main pscf_1d program files

# Path to pscf_1d *.o and *.d files for pscf_1d
PSCF_1D=$(BLD_DIR)/r1d/pscf_1d

# Path to pscf_1d executable
PSCF_1D_EXE=$(BIN_DIR)/pscf_1d

#-----------------------------------------------------------------------
# Main targets 

all: $(r1d_OBJS) $(r1d_LIB) $(PSCF_1D_EXE)

clean:
	rm -f $(r1d_OBJS) $(r1d_OBJS:.o=.d)
	rm -f $(PSCF_1D).o $(PSCF_1D).d
	rm -f $(r1d_LIB)
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f lib*.a

# Executable target (linking)
$(PSCF_1D_EXE): $(PSCF_1D).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_1D_EXE) $(PSCF_1D).o $(LIBS) 

# Short name for executable target (for convenience)
pscf_1d:
	$(MAKE) $(PSCF_1D_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(r1d_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_1D).d 
