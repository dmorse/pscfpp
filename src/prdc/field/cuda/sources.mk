prdc_field_cuda_CU= \
  prdc/field/cuda/RField.cu \
  prdc/field/cuda/RFieldDft.cu \
  prdc/field/cuda/CField.cu \
  prdc/field/cuda/HostDArrayComplex.cu \
  prdc/field/cuda/RFieldComparison.cu \
  prdc/field/cuda/RFieldDftComparison.cu \
  prdc/field/cuda/CFieldComparison.cu \
  prdc/field/cuda/FFT.cu \
  prdc/field/cuda/FFTBatched.cu \
  prdc/field/cuda/WaveList.cu 

prdc_field_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_field_cuda_CU:.cu=.o))

prdc_field_cuda_DEPS=\
     $(addprefix $(BLD_DIR)/, $(prdc_field_cuda_CU:.cu=.du))
