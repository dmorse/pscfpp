#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/rpg/sources.mk
include $(SRC_DIR)/rpg/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(rpg_scft_sweep_OBJS) 

clean:
	rm -f $(rpg_scft_sweep_OBJS) $(rpg_scft_sweep_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(rpg_scft_sweep_OBJS:.o=.d)
