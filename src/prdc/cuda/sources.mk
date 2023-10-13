prdc_cuda_= \
  prdc/cuda/RField.cu \
  prdc/cuda/RFieldDft.cu \
  prdc/cuda/FFT.cu 

prdc_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_cuda_:.cu=.o))

