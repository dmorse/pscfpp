prdc_fieldIo_CPP= \
  prdc/fieldIo/rFieldIo.cpp

prdc_fieldIo_OBJS=\
  $(addprefix $(BLD_DIR)/, $(prdc_fieldIo_CPP:.cpp=.o))

prdc_fieldIo_DEPS=\
  $(addprefix $(BLD_DIR)/, $(prdc_fieldIo_CPP:.cpp=.d))

