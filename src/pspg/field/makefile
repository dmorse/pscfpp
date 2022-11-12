#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk

#--------------------------------------------------------------------
# Main targets 

.PHONY= all clean includes

all: $(pspg_field_OBJS) 

clean:
	rm -f $(pspg_field_OBJS) $(pspg_field_OBJS:.o=.d) 

includes:
	echo $(INCLUDES)

#--------------------------------------------------------------------
# Include dependency files

-include $(pspg_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
