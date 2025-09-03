#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/rpc/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpc_fts_simulator_OBJS) 

clean:
	rm -f $(rpc_fts_simulator_OBJS) $(rpc_fts_simulator_OBJS:.o=.d) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(rpc_fts_simulator_OBJS:.o=.d)
