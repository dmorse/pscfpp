#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pssp_gpu/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pssp_gpu_inter_OBJS) 

clean:
	rm -f $(pssp_gpu_inter_OBJS) $(pssp_gpu_inter_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(pssp_gpu_OBJS:.o=.d)
