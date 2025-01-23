pscf_cuda_ =\
   pscf/cuda/ThreadGrid.cu \
   pscf/cuda/CudaRandom.cu

pscf_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_cuda_:.cu=.o))

