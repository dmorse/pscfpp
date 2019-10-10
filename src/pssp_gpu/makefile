#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/pssp_gpu/include.mk

#-----------------------------------------------------------------------
# Variable definition

PSCF_PG=$(BLD_DIR)/pssp_gpu/pscf_pg
#-----------------------------------------------------------------------
# Main targets 

all: $(pssp_gpu_OBJS) $(pssp_gpu_LIB) $(PSCF_PG_EXE)

move:
	mv $(BIN_DIR)/pscf_pgsp $(BIN_DIR)/../../pscfpp_test/pscf_pgsp_openmp_for

run1:
	cd ../../../pscfpp_test && make run4

clean:
	rm -f $(pssp_gpu_OBJS) $(pssp_gpu_OBJS:.o=.d)
	rm -f $(PSCF_PG).o $(PSCF_PG).d
	rm -f $(pssp_gpu_LIB)
	rm -f *.o *.d */*.o */*.d
	#cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f lib*.a

# Executable target

$(PSCF_PG_EXE): $(PSCF_PG).o $(PSSP_gpu_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PG_EXE) $(PSCF_PG).o $(LIBS)

# Short name for executable target (for convenience)
pscf_pg:
	$(MAKE) $(PSCF_PG_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pssp_gpu_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PG).d 
