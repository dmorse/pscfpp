#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/rpc/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(rpc_fts_mcmove_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(rpc_fts_mcmove_OBJS) $(rpc_fts_mcmove_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(rpc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
