prdc_field_cpp_CPP= \
  prdc/field/cpp/RField.cpp \
  prdc/field/cpp/RFieldDft.cpp \
  prdc/field/cpp/CField.cpp \
  prdc/field/cpp/FFT.cpp \
  prdc/field/cpp/RFieldComparison.cpp \
  prdc/field/cpp/RFieldDftComparison.cpp \
  prdc/field/cpp/CFieldComparison.cpp \
  prdc/field/cpp/FieldBasisConverter.cpp \
  prdc/field/cpp/WaveList.cpp

prdc_field_cpp_OBJS=\
  $(addprefix $(BLD_DIR)/, $(prdc_field_cpp_CPP:.cpp=.o))

prdc_field_cpp_DEPS=\
  $(addprefix $(BLD_DIR)/, $(prdc_field_cpp_CPP:.cpp=.d))

