#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/rpc/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpc_fts_OBJS) 

clean:
	rm -f $(rpc_fts_OBJS) $(rpc_fts_OBJS:.o=.d) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(rpc_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
