#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk

#-----------------------------------------------------------------------
# Main program source file base names 

PSCF_PS1D=$(BLD_DIR)/pspc/pscf_ps1d
PSCF_PS2D=$(BLD_DIR)/pspc/pscf_ps2d
PSCF_PS3D=$(BLD_DIR)/pspc/pscf_ps3d

PSCF_PS_EXE = $(PSCF_PS1D_EXE) $(PSCF_PS2D_EXE) $(PSCF_PS3D_EXE)
#-----------------------------------------------------------------------
# Main targets 

all: $(pspc_OBJS) $(pspc_LIB) $(PSCF_PS_EXE)

clean:
	rm -f $(pspc_OBJS) $(pspc_OBJS:.o=.d)
	rm -f $(pspc_LIB)
	rm -f $(PSCF_PS1D).o $(PSCF_PS1D).d
	rm -f $(PSCF_PS2D).o $(PSCF_PS2D).d
	rm -f $(PSCF_PS3D).o $(PSCF_PS3D).d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 
	-rm -f lib*.a

# Executable targets

$(PSCF_PS1D_EXE): $(PSCF_PS1D).o $(PSPC_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PS1D_EXE) $(PSCF_PS1D).o $(LIBS)

$(PSCF_PS2D_EXE): $(PSCF_PS2D).o $(PSPC_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PS2D_EXE) $(PSCF_PS2D).o $(LIBS)

$(PSCF_PS3D_EXE): $(PSCF_PS3D).o $(PSPC_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PS3D_EXE) $(PSCF_PS3D).o $(LIBS)


# Short name for executable target (for convenience)
pscf_ps1d:
	$(MAKE) $(PSCF_PS1D_EXE)

pscf_ps2d:
	$(MAKE) $(PSCF_PS2D_EXE)

pscf_ps3d:
	$(MAKE) $(PSCF_PS3D_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pspc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PS1D).d 
-include $(PSCF_PS2D).d 
-include $(PSCF_PS3D).d 
