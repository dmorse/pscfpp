pscf_inter_= \
  pscf/inter/Interaction.cpp 

pscf_inter_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_inter_:.cpp=.o))

