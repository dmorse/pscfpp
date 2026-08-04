prdc_field_cpu_CPP= \
  prdc/field/cpu/RField.cpp \
  prdc/field/cpu/RFieldDft.cpp \
  prdc/field/cpu/CField.cpp \
  prdc/field/cpu/FFT.cpp \
  prdc/field/cpu/RFieldComparison.cpp \
  prdc/field/cpu/RFieldDftComparison.cpp \
  prdc/field/cpu/CFieldComparison.cpp \
  prdc/field/cpu/FieldBasisConverter.cpp \
  prdc/field/cpu/WaveList.cpp

prdc_field_cpu_OBJS=\
  $(addprefix $(BLD_DIR)/, $(prdc_field_cpu_CPP:.cpp=.o))

