#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(pscf_field_cpp_OBJS) 

clean:
	rm -f $(pscf_field_cpp_OBJS) $(pscf_field_cpp_DEPS) 

#--------------------------------------------------------------------
# Include dependency files

-include $(pscf_field_cpp_DEPS)
