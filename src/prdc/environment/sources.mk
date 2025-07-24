prdc_environment_= \
  prdc/environment/Environment.cpp \
  prdc/environment/FieldGenerator.cpp \
  prdc/environment/MixAndMatchEnv.cpp \
  prdc/environment/FilmFieldGenMaskBase.cpp \
  prdc/environment/FilmFieldGenExtBase.cpp

prdc_environment_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_environment_:.cpp=.o))

