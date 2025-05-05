prdc_cpu_= \
  prdc/cpu/complex.cpp \
  prdc/cpu/RField.cpp \
  prdc/cpu/RFieldDft.cpp \
  prdc/cpu/CField.cpp \
  prdc/cpu/FFT.cpp \
  prdc/cpu/RFieldComparison.cpp \
  prdc/cpu/RFieldDftComparison.cpp \
  prdc/cpu/CFieldComparison.cpp \
  prdc/cpu/FieldBasisConverter.cpp \
  prdc/cpu/WaveList.cpp \
  prdc/cpu/Reduce.cpp \
  prdc/cpu/complex.cpp

prdc_cpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_cpu_:.cpp=.o))

