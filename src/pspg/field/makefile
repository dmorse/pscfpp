#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(pspg_field_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(pspg_field_OBJS) $(pspg_field_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(pssp_OBJS:.o=.d)
