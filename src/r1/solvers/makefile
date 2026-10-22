#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/r1/sources.mk
include $(SRC_DIR)/r1/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(r1_solvers_OBJS) 

clean:
	rm -f $(r1_solvers_OBJS) $(r1_solvers_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(r1_solvers_OBJS:.o=.d)
