pspc_field_= \
  pspc/field/FFT.cpp \
  pspc/field/RField.cpp \
  pspc/field/RFieldDft.cpp \
  pspc/field/FFT.cpp \
  pspc/field/FieldIo.cpp \
  pspc/field/Domain.cpp \
  pspc/field/BFieldComparison.cpp \
  pspc/field/RFieldComparison.cpp \
  pspc/field/KFieldComparison.cpp \
  pspc/field/WFieldContainer.cpp \
  pspc/field/CFieldContainer.cpp \
  pspc/field/Mask.cpp 


pspc_field_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_field_))
pspc_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_field_:.cpp=.o))

