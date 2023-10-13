prdc_gpu_= \
  prdc/gpu/RField.cu \
  prdc/gpu/RFieldDft.cu \
  prdc/gpu/FFT.cu 

prdc_gpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_gpu_:.cu=.o))

