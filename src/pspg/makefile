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
# Main program source file base names

PSCF_PG1=$(BLD_DIR)/pspg/pscf_pg1
PSCF_PG2=$(BLD_DIR)/pspg/pscf_pg2
PSCF_PG3=$(BLD_DIR)/pspg/pscf_pg3

PSCF_PG_EXE=$(PSCF_PG3_EXE) $(PSCF_PG2_EXE) $(PSCF_PG1_EXE)
#-----------------------------------------------------------------------
# Main targets 

all: $(pspg_OBJS) $(pspg_LIB) $(PSCF_PG_EXE)

clean:
	rm -f $(pspg_OBJS) $(pspg_OBJS:.o=.d)
	rm -f $(pspg_LIB)
	rm -f $(PSCF_PG1).o $(PSCF_PG1).d
	rm -f $(PSCF_PG3).o $(PSCF_PG2).d
	rm -f $(PSCF_PG3).o $(PSCF_PG3).d
	rm -f *.o *.d */*.o */*.d
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f lib*.a

# Executable target

$(PSCF_PG1_EXE): $(PSCF_PG1).o $(PSPG_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PG1_EXE) $(PSCF_PG1).o $(LIBS)

$(PSCF_PG2_EXE): $(PSCF_PG2).o $(PSPG_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PG2_EXE) $(PSCF_PG2).o $(LIBS)

$(PSCF_PG3_EXE): $(PSCF_PG3).o $(PSPG_LIBS)
	$(CXX) $(LDFLAGS) -o $(PSCF_PG3_EXE) $(PSCF_PG3).o $(LIBS)

# Short name for executable target (for convenience)
pscf_pg1:
	$(MAKE) $(PSCF_PG1_EXE)

pscf_pg2:
	$(MAKE) $(PSCF_PG2_EXE)

pscf_pg3:
	$(MAKE) $(PSCF_PG3_EXE)


#-----------------------------------------------------------------------
# Include dependency files

-include $(pspg_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
-include $(PSCF_PG1).d 
-include $(PSCF_PG3).d 
-include $(PSCF_PG3).d 
