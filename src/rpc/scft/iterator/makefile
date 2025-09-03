#--------------------------------------------------------------------
# Include *.mk makefile fragments

SRC_DIR_REL =../../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/rpc/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpc_scft_iterator_OBJS) 

clean:
	rm -f $(rpc_scft_iterator_OBJS) $(rpc_scft_iterator_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(rpc_scft_iterator_OBJS:.o=.d)
