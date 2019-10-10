#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/cyln/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(cyln_misc_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(cyln_misc_OBJS) $(cyln_misc_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(cyln_OBJS:.o=.d)
