pscf_crystal_CPP= \
  pscf/crystal/UnitCellBase.cpp \
  pscf/crystal/UnitCell1.cpp \
  pscf/crystal/UnitCell2.cpp \
  pscf/crystal/UnitCell3.cpp \
  pscf/crystal/shiftToMinimum.cpp \
  pscf/crystal/SpaceSymmetry.cpp \
  pscf/crystal/SymmetryGroup.cpp \
  pscf/crystal/SpaceGroup.cpp \
  pscf/crystal/Basis.cpp \
  pscf/crystal/groupFile.cpp \
  pscf/crystal/BFieldComparison.cpp \
  pscf/crystal/getDimension.cpp \
  pscf/crystal/paramIdConversions.cpp \
  pscf/crystal/replicateUnitCell.cpp \
  pscf/crystal/hasVariableAngle.cpp \
  pscf/crystal/sortWaves.cpp

pscf_crystal_OBJS=\
  $(addprefix $(BLD_DIR)/, $(pscf_crystal_CPP:.cpp=.o))

pscf_crystal_DEPS=\
  $(addprefix $(BLD_DIR)/, $(pscf_crystal_CPP:.cpp=.d))

