prdc_cuda_= \
  prdc/cuda/RField.cu \
  prdc/cuda/RFieldDft.cu \
  prdc/cuda/CField.cu \
  prdc/cuda/RFieldComparison.cu \
  prdc/cuda/RFieldDftComparison.cu \
  prdc/cuda/CFieldComparison.cu \
  prdc/cuda/FFT.cu \
  prdc/cuda/FFTBatched.cu \
  prdc/cuda/Reduce.cu \
  prdc/cuda/VecOp.cu \
  prdc/cuda/VecOpMisc.cu

prdc_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_cuda_:.cu=.o))

