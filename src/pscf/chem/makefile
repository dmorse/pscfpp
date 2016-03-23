#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/pscf/config.mk
include $(SRC_DIR)/pscf/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pscf_chem_OBJS) 

clean:
	rm -f $(pscf_chem_OBJS) $(pscf_chem_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_OBJS:.o=.d)
