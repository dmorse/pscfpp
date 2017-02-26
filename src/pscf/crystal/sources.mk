pscf_crystal_= \
  pscf/crystal/UnitCell.cpp \
  pscf/crystal/shiftToMinimum.cpp 

pscf_crystal_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_crystal_))
pscf_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_crystal_:.cpp=.o))

