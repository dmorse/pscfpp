prdc_fieldIo_= \
  prdc/fieldIo/rFieldIo.cpp

prdc_fieldIo_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_fieldIo_:.cpp=.o))

