#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/rpg/sources.mk
include $(SRC_DIR)/rpg/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpg_scft_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(rpg_scft_OBJS) $(rpg_scft_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(rpg_scft_OBJS:.o=.d)
