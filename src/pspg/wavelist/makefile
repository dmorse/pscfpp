#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(pspg_wavelist_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(pspg_wavelist_OBJS) $(pspg_wavelist_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(pspg_OBJS:.o=.d)
