#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/rpc/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpc_fts_perturbation_OBJS) 

clean:
	rm -f $(rpc_fts_perturbation_OBJS) 
	rm -f $(rpc_fts_perturbation_OBJS:.o=.d) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(rpc_fts_perturbation_OBJS:.o=.d)
