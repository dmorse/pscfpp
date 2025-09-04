#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/r1d/sources.mk
include $(SRC_DIR)/r1d/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(r1d_domain_OBJS) 

clean:
	rm -f $(r1d_domain_OBJS) $(r1d_domain_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(r1d_domain_OBJS:.o=.d)
