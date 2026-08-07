#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/rpg/sources.mk
include $(SRC_DIR)/rpg/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpg_fts_brownian_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(rpg_fts_brownian_OBJS) $(rpg_fts_brownian_OBJS:.o=.du) 

#--------------------------------------------------------------------
# Include dependency files

-include $(rpg_fts_brownian_OBJS:.o=.du)
