#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pssp/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pspg_crystal_OBJS) 

clean:
	rm -f $(pspg_crystal_OBJS) $(pspg_crystal_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(pspg_OBJS:.o=.d)
