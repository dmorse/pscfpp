#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/r1d/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(r1d_sweep_OBJS) 

clean:
	rm -f $(r1d_sweep_OBJS) $(r1d_sweep_OBJS:.o=.d)

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 

#-----------------------------------------------------------------------
# Include dependency files

-include $(r1d_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
