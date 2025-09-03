#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/prdc/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(prdc_field_OBJS) 

clean:
	rm -f $(prdc_field_OBJS) $(prdc_field_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_field_OBJS:.o=.d)
