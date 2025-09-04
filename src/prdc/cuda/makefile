#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(prdc_cuda_OBJS) 

clean:
	rm -f $(prdc_cuda_OBJS) $(prdc_cuda_OBJS:.o=.d) 

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(prdc_cuda_OBJS:.o=.d)
