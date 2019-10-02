pspc_field_= \
  pspc/field/FFT.cpp \
  pspc/field/RField.cpp \
  pspc/field/RFieldDft.cpp \
  pspc/field/FFT.cpp \
  pspc/field/FieldIo.cpp 

pspc_field_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_field_))
pspc_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_field_:.cpp=.o))

