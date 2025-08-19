rpc_field_= \
  rpc/field/FieldIo.cpp \
  rpc/field/Domain.cpp \
  rpc/field/WFields.cpp \
  rpc/field/CFields.cpp \
  rpc/field/Mask.cpp
  #rpc/field/BFieldComparison.cpp \

rpc_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_field_:.cpp=.o))

