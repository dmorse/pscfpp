#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/prdc/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(prdc_cpu_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(prdc_cpu_OBJS) $(prdc_cpu_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
