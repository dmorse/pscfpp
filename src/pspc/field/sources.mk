pspc_field_= \
  pspc/field/FieldIo.cpp \
  pspc/field/Domain.cpp \
  pspc/field/WFieldContainer.cpp \
  pspc/field/CFieldContainer.cpp \
  pspc/field/Mask.cpp 
  #pspc/field/BFieldComparison.cpp \

pspc_field_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_field_))
pspc_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_field_:.cpp=.o))

