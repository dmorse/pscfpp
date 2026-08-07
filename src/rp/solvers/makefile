#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(rp_solvers_OBJS) 

clean:
	rm -f $(rp_solvers_OBJS)
	rm -f $(rp_solvers_DEPS)
	rm -f *.o* *.d*

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(rp_solvers_DEPS)
