#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/r1/sources.mk
include $(SRC_DIR)/r1/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(r1_iterator_OBJS) 

clean:
	rm -f $(r1_iterator_OBJS) $(r1_iterator_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(r1_iterator_OBJS:.o=.d)
