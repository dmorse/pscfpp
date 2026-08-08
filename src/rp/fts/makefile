#--------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/rp/sources.mk
include $(SRC_DIR)/rp/patterns.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rp_fts_OBJS) 

clean:
	rm -f $(rp_fts_OBJS)
	rm -f $(rp_fts_DEPS)
	rm -f *.o *.d */*.o */*.d
	rm -f *.ou *.du */*.ou */*.du

.PHONY: all clean

#--------------------------------------------------------------------
# Include dependency files

-include $(rp_fts_DEPS)
