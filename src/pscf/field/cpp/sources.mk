pscf_field_cpp_CPP= \
  pscf/field/cpp/RField.cpp \
  pscf/field/cpp/RFieldDft.cpp \
  pscf/field/cpp/CField.cpp \
  pscf/field/cpp/FFT.cpp \
  pscf/field/cpp/RFieldComparison.cpp \
  pscf/field/cpp/RFieldDftComparison.cpp \
  pscf/field/cpp/CFieldComparison.cpp \
  pscf/field/cpp/WaveList.cpp

pscf_field_cpp_OBJS=\
  $(addprefix $(BLD_DIR)/, $(pscf_field_cpp_CPP:.cpp=.o))

pscf_field_cpp_DEPS=\
  $(addprefix $(BLD_DIR)/, $(pscf_field_cpp_CPP:.cpp=.d))

