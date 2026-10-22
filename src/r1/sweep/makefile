#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/r1/sources.mk
include $(SRC_DIR)/r1/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(r1_sweep_OBJS) 

clean:
	rm -f $(r1_sweep_OBJS) $(r1_sweep_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(r1_sweep_OBJS:.o=.d)
