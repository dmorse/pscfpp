#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/rpg/include.mk

#--------------------------------------------------------------------
# Main targets 

.PHONY= all clean includes

all: $(rpg_field_OBJS) 

clean:
	rm -f $(rpg_field_OBJS) $(rpg_field_OBJS:.o=.d) 

includes:
	echo $(INCLUDES)

#--------------------------------------------------------------------
# Include dependency files

-include $(rpg_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
