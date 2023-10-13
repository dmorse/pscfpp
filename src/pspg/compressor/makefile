#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(pspg_compressor_OBJS) 

includes:
	@echo $(INCLUDES)

clean:
	rm -f $(pspg_compressor_OBJS) $(pspg_compressor_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(pspg_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
