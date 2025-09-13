#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/rpg/sources.mk
include $(SRC_DIR)/rpg/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(rpg_solvers_OBJS) 

clean:
	rm -f $(rpg_solvers_OBJS) $(rpg_solvers_OBJS:.o=.d)

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 

.PHONY: all clean veryclean

#-----------------------------------------------------------------------
# Include dependency files

-include $(rpg_solvers_OBJS:.o=.d)
