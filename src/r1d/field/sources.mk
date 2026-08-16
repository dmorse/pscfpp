r1d_field_CPP=\
  r1d/field/GeometryMode.cpp \
  r1d/field/Domain.cpp \
  r1d/field/FieldIo.cpp 

r1d_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_field_CPP:.cpp=.o))

r1d_field_DEPS=\
     $(addprefix $(BLD_DIR)/, $(r1d_field_CPP:.cpp=.d))

