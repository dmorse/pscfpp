#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/cyln/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(cyln_field_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(cyln_field_OBJS) $(cyln_field_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(cyln_OBJS:.o=.d)
