#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk

#-----------------------------------------------------------------------
# Main targets 
#
.PHONY: all list clean veryclean

all: $(pscf_cuda_OBJS)
	#@echo $(pscf_cuda_)
	#@echo $(pscf_cuda_OBJS)

list: 
	@echo $(pscf_cuda_)
	@echo $(pscf_cuda_OBJS)

clean:
	rm -f $(pscf_cuda_OBJS) $(pscf_cuda_OBJS:.o=.d)

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
