pscf_environment_= \
  pscf/environment/FieldGenerator.cpp \
  pscf/environment/MixAndMatchEnv.cpp

pscf_environment_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_environment_:.cpp=.o))

