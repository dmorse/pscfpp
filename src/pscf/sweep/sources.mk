pscf_sweep_= \
  pscf/sweep/ParameterModifier.cpp \
  pscf/sweep/ParameterType.cpp

pscf_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_sweep_:.cpp=.o))

