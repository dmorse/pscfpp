pscf_interaction_CPP= \
  pscf/interaction/Interaction.cpp 

pscf_interaction_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_interaction_CPP:.cpp=.o))

pscf_interaction_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_interaction_CPP:.cpp=.d))

