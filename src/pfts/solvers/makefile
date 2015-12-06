# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/pfts/config.mk
include $(SRC_DIR)/pfts/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pfts/sources.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pfts_solvers_OBJS)

clean:
	rm -f $(pfts_solvers_OBJS) $(pfts_solvers_OBJS:.o=.d)

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 

#-----------------------------------------------------------------------
# Include dependency files

-include $(pfts_solvers_OBJS:.o=.d)
