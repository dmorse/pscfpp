rpc_environment_= \
  rpc/environment/FilmFieldGenMask.cpp \
  rpc/environment/FilmFieldGenExt.cpp \
  rpc/environment/EnvironmentFactory.cpp \
  rpc/environment/MixAndMatchEnvs.cpp

rpc_environment_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_environment_:.cpp=.o))

