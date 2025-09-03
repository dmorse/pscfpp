#-----------------------------------------------------------------------
# Include *.mk makefile fragments

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk

#-----------------------------------------------------------------------
# Main targets 
 
.PHONY: all list clean veryclean

all: $(pscf_cuda_OBJS)

clean:
	rm -f $(pscf_cuda_OBJS) $(pscf_cuda_OBJS:.o=.d)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_cuda_OBJS:.o=.d)
