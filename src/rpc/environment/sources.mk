rpc_environment_= \
  rpc/environment/FilmFieldGenMask.cpp \
  rpc/environment/FilmFieldGenExt.cpp \
  rpc/environment/FilmEnvironment.cpp \
  rpc/environment/EnvironmentFactory.cpp 

rpc_environment_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_environment_:.cpp=.o))

