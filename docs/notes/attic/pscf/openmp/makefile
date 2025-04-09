#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk
#ifndef PSCF_OPENMP
include $(SRC_DIR)/pscf/openmp/sources.mk
#endif

#-----------------------------------------------------------------------
# Main targets 

all: $(pscf_openmp_OBJS) 

clean:
	rm -f $(pscf_openmp_OBJS) $(pscf_openmp_OBJS:.o=.d) 
	rm -f *.o *.d

#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_OBJS:.o=.d)
