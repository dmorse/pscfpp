prdc_cpu_= \
  prdc/cpu/Field.cpp \
  prdc/cpu/RField.cpp \
  prdc/cpu/RFieldDft.cpp \
  prdc/cpu/FFT.cpp \
  prdc/cpu/RFieldComparison.cpp \
  prdc/cpu/RFieldDftComparison.cpp 

prdc_cpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_cpu_:.cpp=.o))

