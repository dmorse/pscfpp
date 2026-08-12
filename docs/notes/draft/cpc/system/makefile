#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/cpc/sources.mk
include $(SRC_DIR)/cpc/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(cpc_system_OBJS) 

clean:
	rm -f $(cpc_system_OBJS) $(cpc_system_OBJS:.o=.d) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(cpc_system_OBJS:.o=.d)
