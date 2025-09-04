#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(rpc_environment_OBJS) 

clean:
	rm -f $(rpc_environment_OBJS) $(rpc_environment_OBJS:.o=.d) 

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(rpc_environment_OBJS:.o=.d)
