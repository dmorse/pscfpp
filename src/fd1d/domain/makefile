#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/fd1d/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(fd1d_domain_OBJS) 

clean:
	rm -f $(fd1d_domain_OBJS) $(fd1d_domain_OBJS:.o=.d)

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 

#-----------------------------------------------------------------------
# Include dependency files

-include $(fd1d_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
