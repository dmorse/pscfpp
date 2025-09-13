#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/rpg/sources.mk
include $(SRC_DIR)/rpg/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(rpg_environment_OBJS) 

clean:
	rm -f $(rpg_environment_OBJS) $(rpg_environment_OBJS:.o=.d) 

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(rpg_environment_OBJS:.o=.d)
