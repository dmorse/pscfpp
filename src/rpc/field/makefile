#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpc_field_OBJS) 

clean:
	rm -f $(rpc_field_OBJS) $(rpc_field_OBJS:.o=.d) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(rpc_field_OBJS:.o=.d)
