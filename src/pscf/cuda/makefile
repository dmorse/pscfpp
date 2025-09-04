#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

#-----------------------------------------------------------------------
# Main targets 
 
.PHONY: all list clean veryclean

all: $(pscf_cuda_OBJS)

clean:
	rm -f $(pscf_cuda_OBJS) $(pscf_cuda_OBJS:.o=.d)

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_cuda_OBJS:.o=.d)
