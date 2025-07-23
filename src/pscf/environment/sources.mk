pscf_environment_= \
  pscf/environment/FieldGeneratorBase.cpp

pscf_environment_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_environment_:.cpp=.o))

