pscf_sweep_CPP= \
  pscf/sweep/ParameterModifier.cpp \
  pscf/sweep/ParameterType.cpp

pscf_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_sweep_CPP:.cpp=.o))

pscf_sweep_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_sweep_CPP:.cpp=.d))

