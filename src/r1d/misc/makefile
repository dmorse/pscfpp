#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/r1d/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(r1d_misc_OBJS) 

clean:
	rm -f $(r1d_misc_OBJS) $(r1d_misc_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(r1d_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
