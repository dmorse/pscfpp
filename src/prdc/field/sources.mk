prdc_field_= \
  prdc/field/DomainReal.cpp

prdc_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_field_:.cpp=.o))

