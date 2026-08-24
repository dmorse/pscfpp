#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

#-----------------------------------------------------------------------
# Main targets

all: $(pscf_iterator_OBJS)

clean:
	rm -f $(pscf_iterator_OBJS)
	rm -f $(pscf_iterator_DEPS)

.PHONY: all clean
#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_iterator_DEPS)
