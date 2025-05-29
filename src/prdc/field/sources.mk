prdc_field_= \
  prdc/field/fieldIoUtil.cpp

prdc_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_field_:.cpp=.o))

