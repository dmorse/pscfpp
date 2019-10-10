#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pspg_inter_OBJS) 

clean:
	rm -f $(pspg_inter_OBJS) $(pspg_inter_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(pspg_OBJS:.o=.d)
