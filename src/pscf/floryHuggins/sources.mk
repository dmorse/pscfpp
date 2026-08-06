pscf_floryHuggins_CPP= \
  pscf/floryHuggins/FhClump.cpp \
  pscf/floryHuggins/FhMolecule.cpp \
  pscf/floryHuggins/FhMixture.cpp \
  pscf/floryHuggins/FhInteraction.cpp

pscf_floryHuggins_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_floryHuggins_CPP:.cpp=.o))

