#-----------------------------------------------------------------------
# Include *.mk makefile fragments

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/r1d/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(r1d_solvers_OBJS) 

clean:
	rm -f $(r1d_solvers_OBJS) $(r1d_solvers_OBJS:.o=.d)

.PHONY: all clean

#-----------------------------------------------------------------------
# Include dependency files

-include $(r1d_solvers_OBJS:.o=.d)
