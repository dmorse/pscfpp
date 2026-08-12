#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/cpc/sources.mk
include $(SRC_DIR)/cpc/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(cpc_field_OBJS) 

clean:
	rm -f $(cpc_field_OBJS) $(cpc_field_OBJS:.o=.d) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(cpc_field_OBJS:.o=.d)
