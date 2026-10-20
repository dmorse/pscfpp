#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pscf_fieldIo_OBJS) 

clean:
	rm -f $(pscf_fieldIo_OBJS) $(pscf_fieldIo_DEPS) 

.PHONY: all clean
#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_fieldIo_DEPS)
