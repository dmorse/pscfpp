#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/prdc/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(prdc_gpu_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(prdc_gpu_OBJS) $(prdc_gpu_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
