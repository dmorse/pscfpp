#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/prdc/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(prdc_field_OBJS) 

clean:
	rm -f $(prdc_field_OBJS) $(prdc_field_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_OBJS:.o=.d)
