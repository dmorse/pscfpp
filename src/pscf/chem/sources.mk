pscf_chem_CPP= \
  pscf/chem/Monomer.cpp \
  pscf/chem/Edge.cpp \
  pscf/chem/Vertex.cpp \
  pscf/chem/Ensemble.cpp \
  pscf/chem/PolymerType.cpp \
  pscf/chem/PolymerModel.cpp \
  pscf/chem/Species.cpp \
  pscf/chem/PolymerSpecies.cpp \
  pscf/chem/SolventSpecies.cpp \
  pscf/chem/Composition.cpp \
  pscf/chem/VertexIterator.cpp \
  pscf/chem/EdgeIterator.cpp 

pscf_chem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_chem_CPP:.cpp=.o))

pscf_chem_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_chem_CPP:.cpp=.d))
