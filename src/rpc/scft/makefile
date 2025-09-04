#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpc_scft_OBJS) 

clean:
	rm -f $(rpc_scft_OBJS) $(rpc_scft_OBJS:.o=.d) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(rpc_scft_OBJS:.o=.d)
