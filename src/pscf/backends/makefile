#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pscf_backends_OBJS) 

clean:
	rm -f $(pscf_backends_OBJS) 
	rm -f $(pscf_backends_DEPS) 

.PHONY: all clean
#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_backends_DEPS)
