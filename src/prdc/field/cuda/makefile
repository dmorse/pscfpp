#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(prdc_field_cuda_OBJS) 

clean:
	rm -f $(prdc_field_cuda_OBJS) $(prdc_field_cuda_DEPS) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(prdc_field_cuda_DEPS)
