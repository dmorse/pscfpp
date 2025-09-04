#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(rpc_solvers_OBJS) 

clean:
	rm -f $(rpc_solvers_OBJS) $(rpc_solvers_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(rpc_solvers_OBJS:.o=.d)
