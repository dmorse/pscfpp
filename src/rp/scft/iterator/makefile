#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rp_scft_iterator_OBJS) 

clean:
	rm -f $(rp_scft_iterator_OBJS)
	rm -f $(rp_scft_iterator_DEPS)

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(rp_scft_iterator_DEPS)
