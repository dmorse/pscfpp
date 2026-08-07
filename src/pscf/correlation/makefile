#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

#-----------------------------------------------------------------------
# Main targets

all: $(pscf_correlation_OBJS)

clean:
	rm -f $(pscf_correlation_OBJS)
	rm -f $(pscf_correlation_DEPS)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_correlation_DEPS)
