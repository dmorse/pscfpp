#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/rpg/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(rpg_scft_sweep_OBJS) 

clean:
	rm -f $(rpg_scft_sweep_OBJS) $(rpg_scft_sweep_OBJS:.o=.d)

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 

#-----------------------------------------------------------------------
# Include dependency files

-include $(rpg_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
