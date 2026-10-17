pscf_fieldIo_CPP= \
  pscf/fieldIo/rFieldIo.cpp

pscf_fieldIo_OBJS=\
  $(addprefix $(BLD_DIR)/, $(pscf_fieldIo_CPP:.cpp=.o))

pscf_fieldIo_DEPS=\
  $(addprefix $(BLD_DIR)/, $(pscf_fieldIo_CPP:.cpp=.d))

