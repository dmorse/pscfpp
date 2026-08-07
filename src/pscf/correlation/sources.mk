pscf_correlation_CPP= \
  pscf/correlation/Debye.cpp \
  pscf/correlation/Polymer.cpp \
  pscf/correlation/Mixture.cpp 

pscf_correlation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_correlation_CPP:.cpp=.o))

pscf_correlation_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_correlation_CPP:.cpp=.d))

