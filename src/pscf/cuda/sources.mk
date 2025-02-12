pscf_cuda_ =\
   pscf/cuda/ThreadArray.cu \
   pscf/cuda/ThreadMesh.cu \
   pscf/cuda/CudaRandom.cu

pscf_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_cuda_:.cu=.o))

