#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pssp/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(pssp_field_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(pssp_field_OBJS) $(pssp_field_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(pssp_OBJS:.o=.d)
