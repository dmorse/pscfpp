pscf_crystal_= \
  pscf/crystal/UnitCell1.cpp \
  pscf/crystal/UnitCell2.cpp \
  pscf/crystal/UnitCell3.cpp \
  pscf/crystal/shiftToMinimum.cpp \
  pscf/crystal/SpaceSymmetry.cpp \
  pscf/crystal/SymmetryGroup.cpp \
  pscf/crystal/Basis.cpp \
  pscf/crystal/groupFile.cpp

pscf_crystal_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_crystal_))
pscf_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_crystal_:.cpp=.o))

# Specialized rule for groupFile.o
$(BLD_DIR)/pscf/crystal/groupFile.o: $(SRC_DIR)/pscf/crystal/groupFile.cpp 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -DDAT_DIR=$(DAT_DIR) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) -DDAT_DIR=$(DAT_DIR) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

