pscf_homogeneous_= \
  pscf/homogeneous/Clump.cpp \
  pscf/homogeneous/Molecule.cpp \
  pscf/homogeneous/Mixture.cpp 

pscf_homogeneous_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_homogeneous_:.cpp=.o))

