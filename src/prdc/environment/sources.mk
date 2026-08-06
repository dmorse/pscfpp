prdc_environment_CPP= \
  prdc/environment/Environment.cpp \
  prdc/environment/FieldGenerator.cpp \
  prdc/environment/MixAndMatchEnv.cpp \
  prdc/environment/FilmFieldGenMaskBase.cpp \
  prdc/environment/FilmFieldGenExtBase.cpp

prdc_environment_OBJS=\
  $(addprefix $(BLD_DIR)/, $(prdc_environment_CPP:.cpp=.o))

