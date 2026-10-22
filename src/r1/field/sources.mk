r1_field_CPP=\
  r1/field/GeometryMode.cpp \
  r1/field/Domain.cpp \
  r1/field/FieldIo.cpp 

r1_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1_field_CPP:.cpp=.o))

r1_field_DEPS=\
     $(addprefix $(BLD_DIR)/, $(r1_field_CPP:.cpp=.d))

