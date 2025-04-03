pscf_chem_= \
  pscf/chem/Monomer.cpp \
  pscf/chem/Vertex.cpp \
  pscf/chem/BlockDescriptor.cpp \
  pscf/chem/SolventDescriptor.cpp \
  pscf/chem/Species.cpp \
  pscf/chem/PolymerType.cpp \
  pscf/chem/PolymerModel.cpp \
  pscf/chem/Debye.cpp 

pscf_chem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_chem_:.cpp=.o))

