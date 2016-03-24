pscf_crystal_= \
  pscf/crystal/UnitCellTmpl.cpp 

pscf_crystal_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_crystal_))
pscf_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_crystal_:.cpp=.o))

