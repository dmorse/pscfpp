#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(prdc_cpu_OBJS) 

clean:
	rm -f $(prdc_cpu_OBJS) $(prdc_cpu_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(prdc_cpu_OBJS:.o=.d)
