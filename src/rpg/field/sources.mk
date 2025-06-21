rpg_field_= \
  rpg/field/WFieldContainer.cu \
  rpg/field/CFieldContainer.cu \
  rpg/field/Mask.cu \
  rpg/field/HostDArrayComplex.cu \
  rpg/field/FieldIo.cu \
  rpg/field/Domain.cu \
  rpg/field/FilmFieldGenExt.cu \
  rpg/field/FilmFieldGenMask.cu \
  rpg/field/MixAndMatchEnvs.cu \
  rpg/field/EnvironmentFactory.cu

rpg_field_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(rpg_field_:.cu=.o))

