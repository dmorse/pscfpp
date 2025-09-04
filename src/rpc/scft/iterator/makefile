#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/rpc/sources.mk
include $(SRC_DIR)/rpc/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpc_scft_iterator_OBJS) 

clean:
	rm -f $(rpc_scft_iterator_OBJS) $(rpc_scft_iterator_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(rpc_scft_iterator_OBJS:.o=.d)
