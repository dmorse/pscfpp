#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

#-----------------------------------------------------------------------
# Main targets

all: $(pscf_cpu_OBJS)

clean:
	rm -f $(pscf_cpu_OBJS)
	rm -f $(pscf_cpu_DEPS)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_cpu_DEPS)
