#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspc/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(pspc_simulate_mcmove_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(pspc_simulate_mcmove_OBJS) $(pspc_simulate_mcmove_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(pspc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
