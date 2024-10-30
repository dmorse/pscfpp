#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/rpg/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpg_scft_iterator_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(rpg_scft_iterator_OBJS) $(rpg_scft_iterator_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(rpg_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
