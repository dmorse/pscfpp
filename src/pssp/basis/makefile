#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pssp/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(pssp_basis_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(pssp_basis_OBJS) $(pssp_basis_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(pssp_OBJS:.o=.d)
