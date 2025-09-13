#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/rpg/sources.mk
include $(SRC_DIR)/rpg/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpg_fts_trajectory_OBJS) 

clean:
	rm -f $(rpg_fts_trajectory_OBJS) $(rpg_fts_trajectory_OBJS:.o=.d) 

.PHONY: all clean
#--------------------------------------------------------------------
# Include dependency files

-include $(rpg_fts_trajectory_OBJS:.o=.d)
