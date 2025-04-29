prdc_crystal_= \
  prdc/crystal/UnitCell1.cpp \
  prdc/crystal/UnitCell2.cpp \
  prdc/crystal/UnitCell3.cpp \
  prdc/crystal/shiftToMinimum.cpp \
  prdc/crystal/SpaceSymmetry.cpp \
  prdc/crystal/SymmetryGroup.cpp \
  prdc/crystal/SpaceGroup.cpp \
  prdc/crystal/Basis.cpp \
  prdc/crystal/groupFile.cpp \
  prdc/crystal/BFieldComparison.cpp \
  prdc/crystal/getDimension.cpp \
  prdc/crystal/paramIdConversions.cpp \
  prdc/crystal/replicateUnitCell.cpp \
  prdc/crystal/hasVariableAngle.cpp 

prdc_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_crystal_:.cpp=.o))

