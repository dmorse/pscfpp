pscf_homogeneous_= \
  pscf/homogeneous/Group.cpp \
  pscf/homogeneous/Molecule.cpp 

pscf_homogeneous_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_homogeneous_))
pscf_homogeneous_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_homogeneous_:.cpp=.o))

