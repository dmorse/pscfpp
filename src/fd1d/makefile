#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the corresponding subdirectory of the bld/ directory by 
# the setup script. Only the copy in the src/ directory is stored in 
# the repository.
#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/fd1d/include.mk

#-----------------------------------------------------------------------
# Variable definition

PSCF_FD=$(BLD_DIR)/fd1d/pscf_fd
#-----------------------------------------------------------------------
# Main targets 

all: $(fd1d_OBJS) $(fd1d_LIB) $(PSCF_FD_EXE)

clean:
	rm -f $(fd1d_OBJS) $(fd1d_OBJS:.o=.d)
	rm -f $(PSCF_FD).o $(PSCF_FD).d
	rm -f $(fd1d_LIB)
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f *.o */*.o
	-rm -f *.d */*.d
	-rm -f lib*.a


# Executable target

$(PSCF_FD_EXE): $(PSCF_FD).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_FD_EXE) $(PSCF_FD).o $(LIBS) 

# Short name for executable target (for convenience)
pscf_fd:
	$(MAKE) $(PSCF_FD_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(fd1d_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_FD).d 
