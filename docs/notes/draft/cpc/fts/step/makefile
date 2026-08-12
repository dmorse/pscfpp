#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/cpc/sources.mk
include $(SRC_DIR)/cpc/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(cpc_fts_step_OBJS) 

clean:
	rm -f $(cpc_fts_step_OBJS) $(cpc_fts_step_OBJS:.o=.d) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(cpc_fts_step_OBJS:.o=.d)
