prdc_environment_= \
  prdc/environment/FilmFieldGenMaskBase.cpp

prdc_environment_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_environment_:.cpp=.o))

