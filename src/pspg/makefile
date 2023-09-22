#------------------------------------------------------------------------
# src/pspg/makefile
#
# The repository copy of this namespace level makefile is in the src/pspg
# directory. This is copied verbatim by the setup to create the copy in
# in the bld/pspg directory.
#------------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk

#-----------------------------------------------------------------------
# Main program source file base name

PSCF_PG=$(BLD_DIR)/pspg/pscf_pg

#-----------------------------------------------------------------------
# Main targets 

all: $(pspg_OBJS) $(pspg_LIB) $(PSCF_PG_EXE)

clean:
	rm -f $(pspg_OBJS) $(pspg_OBJS:.o=.d)
	rm -f $(pspg_LIB)
	rm -f $(PSCF_PG).o $(PSCF_PG).d
	rm -f *.o */*.o
	rm -f *.d */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f lib*.a

# Executable target
$(PSCF_PG_EXE): $(PSCF_PG).o $(PSPG_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PG_EXE) $(PSCF_PG).o $(LIBS)

# Short name for executable target (for convenience)
pscf_pg:
	$(MAKE) $(PSCF_PG_EXE)


#-----------------------------------------------------------------------
# Include dependency files

-include $(pspg_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PG).d 
