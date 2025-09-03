#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/r1d/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(r1d_iterator_OBJS) 

clean:
	rm -f $(r1d_iterator_OBJS) $(r1d_iterator_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(r1d_iterator_OBJS:.o=.d)
