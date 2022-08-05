#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(pspc_iterator_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(pspc_iterator_OBJS) $(pspc_iterator_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(pspc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
