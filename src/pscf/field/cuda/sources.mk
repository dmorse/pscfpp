pscf_field_cuda_CU= \
  pscf/field/cuda/RField.cu \
  pscf/field/cuda/RFieldDft.cu \
  pscf/field/cuda/CField.cu \
  pscf/field/cuda/RFieldComparison.cu \
  pscf/field/cuda/RFieldDftComparison.cu \
  pscf/field/cuda/CFieldComparison.cu \
  pscf/field/cuda/FFT.cu \
  pscf/field/cuda/FFTBatched.cu \
  pscf/field/cuda/WaveList.cu 

pscf_field_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_field_cuda_CU:.cu=.ou))

pscf_field_cuda_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_field_cuda_CU:.cu=.du))
