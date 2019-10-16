#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk

#-----------------------------------------------------------------------
# Variable definition

PSCF_PG1D=$(BLD_DIR)/pspg/pscf_pg1d
PSCF_PG3D=$(BLD_DIR)/pspg/pscf_pg2d
PSCF_PG3D=$(BLD_DIR)/pspg/pscf_pg3d

PSCF_PG_EXE=$(PSCF_PG3D_EXE)
#-----------------------------------------------------------------------
# Main targets 

all: $(pspg_OBJS) $(pspg_LIB) $(PSCF_PG_EXE)

move:
	mv $(BIN_DIR)/pscf_pgsp $(BIN_DIR)/../../pscfpp_test/pscf_pgsp_openmp_for

run1:
	cd ../../../pscfpp_test && make run4

clean:
	rm -f $(pspg_OBJS) $(pspg_OBJS:.o=.d)
	rm -f $(PSCF_PG1D).o $(PSCF_PG1D).d
	rm -f $(PSCF_PG3D).o $(PSCF_PG2D).d
	rm -f $(PSCF_PG3D).o $(PSCF_PG3D).d
	rm -f $(pspg_LIB)
	rm -f $(pspg_LIB)
	rm -f $(pspg_LIB)
	rm -f *.o *.d */*.o */*.d
	#cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f lib*.a

# Executable target

$(PSCF_PG1D_EXE): $(PSCF_PG1D).o $(PSPG_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PG1D_EXE) $(PSCF_PG1D).o $(LIBS)

$(PSCF_PG2D_EXE): $(PSCF_PG2D).o $(PSPG_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PG3D_EXE) $(PSCF_PG2D).o $(LIBS)

$(PSCF_PG3D_EXE): $(PSCF_PG3D).o $(PSPG_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PG3D_EXE) $(PSCF_PG3D).o $(LIBS)

# Short name for executable target (for convenience)
pscf_pg:
	$(MAKE) $(PSCF_PG1D_EXE)
	$(MAKE) $(PSCF_PG3D_EXE)
	$(MAKE) $(PSCF_PG3D_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pspg_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PG1D).d 
-include $(PSCF_PG3D).d 
-include $(PSCF_PG3D).d 
