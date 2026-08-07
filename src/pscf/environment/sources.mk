pscf_environment_CPP= \
  pscf/environment/FieldGeneratorBase.cpp

pscf_environment_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_environment_CPP:.cpp=.o))

pscf_environment_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_environment_CPP:.cpp=.d))

