pscf_chem_= \
  pscf/chem/Monomer.cpp \
  pscf/chem/Edge.cpp \
  pscf/chem/Vertex.cpp \
  pscf/chem/Species.cpp \
  pscf/chem/SolventSpecies.cpp \
  pscf/chem/PolymerSpecies.cpp \
  pscf/chem/PolymerType.cpp \
  pscf/chem/PolymerModel.cpp \
  pscf/chem/Debye.cpp \
  pscf/chem/VertexIterator.cpp \
  pscf/chem/EdgeIterator.cpp 

pscf_chem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_chem_:.cpp=.o))

