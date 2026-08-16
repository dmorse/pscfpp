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
# Paths associated with main program files

# Lists of object, dependency and executable files for programs
PSCF_RP_OBJS=
PSCF_RP_DEPS=
PSCF_RP_EXES=

# Paths for pscf_rpc CPU program
ifdef PSCF_CPP
  PSCF_RPC=$(BLD_DIR)/rp/pscf_rpc
  PSCF_RPC_EXE=$(BIN_DIR)/pscf_rpc
  PSCF_RP_OBJS+= $(PSCF_RPC).o
  PSCF_RP_DEPS+= $(PSCF_RPC).d
  PSCF_RP_EXES+= $(PSCF_RPC_EXE)
endif

# Paths for pscf_rpg GPU program
ifdef PSCF_CUDA
  PSCF_RPG=$(BLD_DIR)/rp/pscf_rpg
  PSCF_RPG_EXE=$(BIN_DIR)/pscf_rpg
  PSCF_RP_OBJS+= $(PSCF_RPG).ou
  PSCF_RP_DEPS+= $(PSCF_RPG).du
  PSCF_RP_EXES+= $(PSCF_RPG_EXE)
endif

#-----------------------------------------------------------------------
# Main targets 

all: $(PSCF_RP_EXES)

clean:
	rm -f $(rp_OBJS)
	rm -f $(rp_DEPS)
	rm -f $(rp_LIB)
	rm -f $(PSCF_RP_OBJS)
	rm -f $(PSCF_RP_DEPS)
	rm -f *.o */*.o */*/*.o
	rm -f *.d */*.d */*/*.d
	rm -f *.ou */*.ou */*/*.ou
	rm -f *.du */*.du */*/*.du
	rm -f lib*.a 
	cd tests/cpp; $(MAKE) clean
	cd tests/cuda; $(MAKE) clean

ifdef PSCF_CPP

  # Main pscf_rpc executable (linking)
  $(PSCF_RPC_EXE): $(PSCF_RPC).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_RPC_EXE) $(PSCF_RPC).o $(LIBS)

endif

ifdef PSCF_CUDA

  # Main pscf_rpg executable (linking)
  $(PSCF_RPG_EXE): $(PSCF_RPG).ou $(PSCF_LIBS)
	@cp $(PSCF_RPG).ou $(PSCF_RPG).o
	$(NVXX) $(LDFLAGS) -o $(PSCF_RPG_EXE) $(PSCF_RPG).o $(LIBS)
	@rm $(PSCF_RPG).o

endif

# Note: In linking rules for executables, the prerequisite list contains
# only the list $(PSCF_LIBS) of relevant PSCF-specific libraries, but the 
# action links the full list $(LIBS) that also includes relevant external
# libraries.

.PHONY= all clean

#-----------------------------------------------------------------------
# Include all relevant dependency files

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(prdc_DEPS)
-include $(rp_DEPS)
-include $(PSCF_RP_DEPS)
