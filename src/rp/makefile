#------------------------------------------------------------------------
# src/rp/makefile
#
# The repository copy of this namespace level makefile is in the src/rp
# directory. This is copied verbatim by the configure script to create 
# the identical copy in in the bld/rp directory.
#------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

#-----------------------------------------------------------------------
# Paths to executable program files

# List of object files for pscf_rp* programs
PSCF_RP_OBJS=

# List of executable files for pscf_rp* programs
PSCF_RP_EXES=

# Paths for pscf_rpc
ifdef PSCF_CPP
PSCF_RPC=$(BLD_DIR)/rp/pscf_rpc
PSCF_RPC_EXE=$(BIN_DIR)/pscf_rpc
PSCF_RP_OBJS+= $(PSCF_RPC).o
PSCF_RP_EXES+= PSCF_RPC_EXE
endif

# Paths for pscf_rpg
ifdef PSCF_CUDA
PSCF_RPG=$(BLD_DIR)/rp/pscf_rpg
PSCF_RPG_EXE=$(BIN_DIR)/pscf_rpg
PSCF_RP_OBJS+= $(PSCF_RPG).o
PSCF_RP_EXES+= PSCF_RPG_EXE
endif

#-----------------------------------------------------------------------
# Main targets 

.PHONY: all clean

all: $(rp_OBJS) $(rp_LIB) $(PSCF_RP_EXES)

clean:
	rm -f $(rp_OBJS) $(rp_OBJS:.o=.d)
	rm -f $(rp_LIB)
	rm -f $(PSCF_RP_OBJS) $(PSCF_RP_OBJS:.o=.d)
	rm -f *.o */*.o */*/*.o
	rm -f *.d */*.d */*/*.d
	rm -f lib*.a 
	cd tests; $(MAKE) clean

# C++ CPU backend
ifdef PSCF_CPP

# Main pscf_rpc executable
$(PSCF_RPC_EXE): $(PSCF_RPC).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_RPC_EXE) $(PSCF_RPC).o $(LIBS)

.PHONY: pscf_rpc
pscf_rpc:
	$(MAKE) $(PSCF_RPC_EXE)

endif

# Cuda GPU backend
ifdef PSCF_CUDA

# Main pscf_rpg executable 
$(PSCF_RPG_EXE): $(PSCF_RPG).o $(PSCF_LIBS)
	$(NVXX) $(LDFLAGS) -o $(PSCF_RPG_EXE) $(PSCF_RPG).o $(LIBS)

.PHONY: pscf_rpg
pscf_rpg:
	$(MAKE) $(PSCF_RPG_EXE)

endif

# Note: In linking rules for executables, the prerequisite list contains
# only the list $(PSCF_LIBS) of relevant PSCF-specific libraries, but the 
# action links the full list $(LIBS) that also includes relevant external
# libraries.

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(rp_OBJS:.o=.d)
-include $(PSCF_RP_OBJS:.o=.d)
