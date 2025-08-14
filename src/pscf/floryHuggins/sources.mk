pscf_floryHuggins_= \
  pscf/floryHuggins/Clump.cpp \
  pscf/floryHuggins/Molecule.cpp \
  pscf/floryHuggins/Mixture.cpp 

pscf_floryHuggins_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_floryHuggins_:.cpp=.o))

