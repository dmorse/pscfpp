pscf_crystal_= \
  pscf/crystal/UnitCell1.cpp \
  pscf/crystal/UnitCell2.cpp \
  pscf/crystal/UnitCell3.cpp \
  pscf/crystal/shiftToMinimum.cpp \
  pscf/crystal/SpaceSymmetry.cpp \
  pscf/crystal/SymmetryGroup.cpp \
  pscf/crystal/SpaceGroup.cpp \
  pscf/crystal/Basis.cpp \
  pscf/crystal/groupFile.cpp

pscf_crystal_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_crystal_))
pscf_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_crystal_:.cpp=.o))

