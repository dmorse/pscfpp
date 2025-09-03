#-----------------------------------------------------------------------
# Include *.mk makefile fragments

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pscf_environment_OBJS) 

clean:
	rm -f $(pscf_environment_OBJS) $(pscf_environment_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_environment_OBJS:.o=.d)
