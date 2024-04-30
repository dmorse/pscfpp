#------------------------------------------------------------------------
# src/rpg/makefile
#
# The repository copy of this namespace level makefile is in the src/rpg
# directory. This is copied verbatim by the configure script to create 
# the identical copy in in the bld/rpg directory.
#------------------------------------------------------------------------
# Include makefile fragments

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpg/include.mk

#-----------------------------------------------------------------------
# Paths to pscf_pg program files

# Base path for *.o and *.d files for pscf_pg
PSCF_PG=$(BLD_DIR)/rpg/pscf_pg

# Path to pscf_pg executable
PSCF_PG_EXE=$(BIN_DIR)/pscf_pg

#-----------------------------------------------------------------------
# Main targets 

all: $(rpg_OBJS) $(rpg_LIB) $(PSCF_PG_EXE)

clean:
	rm -f $(rpg_OBJS) $(rpg_OBJS:.o=.d)
	rm -f $(rpg_LIB)
	rm -f $(PSCF_PG).o $(PSCF_PG).d
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean

# Main pscf_pg executable file
$(PSCF_PG_EXE): $(PSCF_PG).o $(PSCF_LIBS)
	$(NVXX) $(LDFLAGS) -o $(PSCF_PG_EXE) $(PSCF_PG).o $(LIBS)

# Note: In the linking rule for pscf_pg, we include the list $(PSCF_LIBS) 
# of PSCF-specific libraries as dependencies but link to the list $(LIBS)
# of libraries that includes relevant external libraries.

# Short phony target for executable file (for convenience)
pscf_pg:
	$(MAKE) $(PSCF_PG_EXE)

.PHONY: all clean veryclean pscf_pg

#-----------------------------------------------------------------------
# Include dependency files

-include $(rpg_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PG).d 
