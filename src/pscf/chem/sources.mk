pscf_chem_= \
  pscf/chem/Monomer.cpp \
  pscf/chem/Vertex.cpp \
  pscf/chem/BlockDescriptor.cpp \
  pscf/chem/Species.cpp 

pscf_chem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_chem_))
pscf_chem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_chem_:.cpp=.o))

