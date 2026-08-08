#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../../config.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rp_fts_simulator_OBJS) 

clean:
	rm -f $(rp_fts_simulator_OBJS)
	rm -f $(rp_fts_simulator_DEPS)
	rm -f *.o *.d
	rm -f *.ou *.du

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(rp_fts_simulator_DEPS)
