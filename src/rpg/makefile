#------------------------------------------------------------------------
# src/rpg/makefile
#
# The repository copy of this namespace level makefile is in the src/rpg
# directory. This is copied verbatim by the configure script to create 
# the identical copy in in the bld/rpg directory.
#------------------------------------------------------------------------
# Include *.mk makefile fragments

include ../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/rpg/sources.mk
include $(SRC_DIR)/rpg/patterns.mk

#-----------------------------------------------------------------------
# Paths to pscf_rpg program files

# Base path for *.o and *.d files for pscf_rpg
PSCF_RPG=$(BLD_DIR)/rpg/pscf_rpg

# Path to pscf_rpg executable
PSCF_RPG_EXE=$(BIN_DIR)/pscf_rpg

#-----------------------------------------------------------------------
# Main targets 

all: $(rpg_OBJS) $(rpg_LIB) $(PSCF_RPG_EXE)

clean:
	rm -f $(rpg_OBJS)
	rm -f $(rpg_DEPS)
	rm -f $(rpg_LIB)
	rm -f $(PSCF_RPG).ou $(PSCF_RPG).du
	rm -f *.ou */*.ou */*/*.ou
	rm -f *.du */*.du */*/*.du
	rm -f lib*.a 
	cd tests; $(MAKE) clean

# Main pscf_rpg executable file
$(PSCF_RPG_EXE): $(PSCF_RPG).ou $(PSCF_LIBS)
	@cp $(PSCF_RPG).ou $(PSCF_RPG).o
	$(NVXX) $(LDFLAGS) -o $(PSCF_RPG_EXE) $(PSCF_RPG).o $(LIBS)
	@rm $(PSCF_RPG).o

# Short phony target for executable file (for convenience)
pscf_rpg:
	$(MAKE) $(PSCF_RPG_EXE)

# Note: In the linking rule for pscf_rpg, we include the list $(PSCF_LIBS) 
# of PSCF-specific libraries as pre-requisistes but link to the full list 
# that includes relevant external libraries.

.PHONY: all clean pscf_rpg

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(prdc_DEPS)
-include $(rpg_DEPS)
-include $(PSCF_RPG).du
