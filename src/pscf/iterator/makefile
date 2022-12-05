#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pscf_iterator_OBJS) 

clean:
	rm -f $(pscf_iterator_OBJS) $(pscf_iterator_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_OBJS:.o=.d)
