pscf_iterator_CPP= \
  pscf/iterator/AmbdInteraction.cpp \
  pscf/iterator/NanException.cpp

pscf_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_iterator_CPP:.cpp=.o))

pscf_iterator_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_iterator_CPP:.cpp=.d))

