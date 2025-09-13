#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/rpg/sources.mk
include $(SRC_DIR)/rpg/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpg_scft_iterator_OBJS) 

clean:
	rm -f $(rpg_scft_iterator_OBJS) $(rpg_scft_iterator_OBJS:.o=.d) 

.PHONY: all clean
#--------------------------------------------------------------------
# Include dependency files

-include $(rpg_scft_iterator_OBJS:.o=.d)
