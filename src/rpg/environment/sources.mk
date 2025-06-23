rpg_environment_= \
  rpg/environment/FilmFieldGenMask.cu \
  rpg/environment/FilmFieldGenExt.cu \
  rpg/environment/EnvironmentFactory.cu \
  rpg/environment/MixAndMatchEnvs.cu

rpg_environment_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_environment_:.cu=.o))

