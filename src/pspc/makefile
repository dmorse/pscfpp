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
PSCF_PC1=$(BLD_DIR)/pspc/pscf_pc1
PSCF_PC2=$(BLD_DIR)/pspc/pscf_pc2
PSCF_PC3=$(BLD_DIR)/pspc/pscf_pc3
PSCF_PCA=$(BLD_DIR)/pspc/pscf_pc

# List of executable names (single names are defined in pspc/config.mk)
PSCF_PC_EXE = $(PSCF_PC1_EXE) $(PSCF_PC2_EXE) $(PSCF_PC3_EXE) $(PSCF_PCA_EXE)
#-----------------------------------------------------------------------
# Main targets 

all: $(pspc_OBJS) $(pspc_LIB) $(PSCF_PC_EXE)

clean:
	rm -f $(pspc_OBJS) $(pspc_OBJS:.o=.d)
	rm -f $(pspc_LIB)
	rm -f $(PSCF_PC1).o $(PSCF_PC1).d
	rm -f $(PSCF_PC2).o $(PSCF_PC2).d
	rm -f $(PSCF_PC3).o $(PSCF_PC3).d
	rm -f $(PSCF_PCA).o $(PSCF_PCA).d
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f lib*.a

# Executable targets

$(PSCF_PC1_EXE): $(PSCF_PC1).o $(PSPC_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PC1_EXE) $(PSCF_PC1).o $(LIBS)

$(PSCF_PC2_EXE): $(PSCF_PC2).o $(PSPC_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PC2_EXE) $(PSCF_PC2).o $(LIBS)

$(PSCF_PC3_EXE): $(PSCF_PC3).o $(PSPC_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PC3_EXE) $(PSCF_PC3).o $(LIBS)

$(PSCF_PCA_EXE): $(PSCF_PCA).o $(PSPC_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PCA_EXE) $(PSCF_PCA).o $(LIBS)

# Short name for executable target (for convenience)
pscf_pc1:
	$(MAKE) $(PSCF_PC1_EXE)

pscf_pc2:
	$(MAKE) $(PSCF_PC2_EXE)

pscf_pc3:
	$(MAKE) $(PSCF_PC3_EXE)

pscf_pca:
	$(MAKE) $(PSCF_PCA_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pspc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PC1).d 
-include $(PSCF_PC2).d 
-include $(PSCF_PC3).d 
-include $(PSCF_PCA).d 
