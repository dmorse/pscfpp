#------------------------------------------------------------------------
# src/pspc/makefile
#
# The repository copy of this namespace level makefile is in the src/pspc
# directory. This is copied verbatim by the setup script to create the 
# copy in the bld/pspc directory.
#------------------------------------------------------------------------
# Include makefile fragments

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk

#-----------------------------------------------------------------------
# Main program names

# Source file base names 
PSCF_PC1D=$(BLD_DIR)/pspc/pscf_pc1d
PSCF_PC2D=$(BLD_DIR)/pspc/pscf_pc2d
PSCF_PC3D=$(BLD_DIR)/pspc/pscf_pc3d

# List of executable names (single names are defined in pspc/config.mk)
PSCF_PC_EXE = $(PSCF_PC1D_EXE) $(PSCF_PC2D_EXE) $(PSCF_PC3D_EXE)
#-----------------------------------------------------------------------
# Main targets 

all: $(pspc_OBJS) $(pspc_LIB) $(PSCF_PC_EXE)

clean:
	rm -f $(pspc_OBJS) $(pspc_OBJS:.o=.d)
	rm -f $(pspc_LIB)
	rm -f $(PSCF_PC1D).o $(PSCF_PC1D).d
	rm -f $(PSCF_PC2D).o $(PSCF_PC2D).d
	rm -f $(PSCF_PC3D).o $(PSCF_PC3D).d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f *.o */*.o
	-rm -f *.d */*.d
	-rm -f lib*.a

# Executable targets

$(PSCF_PC1D_EXE): $(PSCF_PC1D).o $(PSPC_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PC1D_EXE) $(PSCF_PC1D).o $(LIBS)

$(PSCF_PC2D_EXE): $(PSCF_PC2D).o $(PSPC_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PC2D_EXE) $(PSCF_PC2D).o $(LIBS)

$(PSCF_PC3D_EXE): $(PSCF_PC3D).o $(PSPC_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PC3D_EXE) $(PSCF_PC3D).o $(LIBS)

# Short name for executable target (for convenience)
pscf_pc1d:
	$(MAKE) $(PSCF_PC1D_EXE)

pscf_pc2d:
	$(MAKE) $(PSCF_PC2D_EXE)

pscf_pc3d:
	$(MAKE) $(PSCF_PC3D_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pspc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PC1D).d 
-include $(PSCF_PC2D).d 
-include $(PSCF_PC3D).d 
