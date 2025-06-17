rpc_field_= \
  rpc/field/FieldIo.cpp \
  rpc/field/Domain.cpp \
  rpc/field/WFieldContainer.cpp \
  rpc/field/CFieldContainer.cpp \
  rpc/field/Mask.cpp \
  rpc/field/FilmFieldGenExt.cpp \
  rpc/field/FilmFieldGenMask.cpp \
  rpc/field/MixAndMatchEnvs.cpp \
  rpc/field/EnvironmentFactory.cpp
  #rpc/field/BFieldComparison.cpp \

rpc_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_field_:.cpp=.o))

