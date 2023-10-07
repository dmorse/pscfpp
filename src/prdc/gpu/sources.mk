prdc_gpu_= \
  prdc/gpu/RField.cpp \
  prdc/gpu/RFieldDft.cpp \
  prdc/gpu/FFT.cpp 

prdc_gpu_SRCS=\
     $(addprefix $(SRC_DIR)/, $(prdc_gpu_))

prdc_gpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_gpu_:.cpp=.o))

