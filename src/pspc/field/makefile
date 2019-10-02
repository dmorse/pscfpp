#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(pspc_field_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(pspc_field_OBJS) $(pspc_field_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(pspc_OBJS:.o=.d)
