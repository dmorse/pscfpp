#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpc_system_OBJS) 

clean:
	rm -f $(rpc_system_OBJS) $(rpc_system_OBJS:.o=.d) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(rpc_system_OBJS:.o=.d)
