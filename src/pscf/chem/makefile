#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pscf_chem_OBJS) 

clean:
	rm -f $(pscf_chem_OBJS) $(pscf_chem_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_chem_OBJS:.o=.d)
