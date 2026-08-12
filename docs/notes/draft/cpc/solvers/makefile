#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/cpc/sources.mk
include $(SRC_DIR)/cpc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(cpc_solvers_OBJS) 

clean:
	rm -f $(cpc_solvers_OBJS) $(cpc_solvers_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(cpc_solvers_OBJS:.o=.d)
