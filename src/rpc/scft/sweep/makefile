#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(rpc_scft_sweep_OBJS) 

clean:
	rm -f $(rpc_scft_sweep_OBJS) $(rpc_scft_sweep_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(rpc_scft_sweep_OBJS:.o=.d)
