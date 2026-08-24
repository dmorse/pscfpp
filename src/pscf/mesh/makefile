#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

#-----------------------------------------------------------------------
# Main targets

all: $(pscf_mesh_OBJS)

clean:
	rm -f $(pscf_mesh_OBJS)
	rm -f $(pscf_mesh_DEPS)

.PHONY: all clean
#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_mesh_DEPS)
