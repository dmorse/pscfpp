pscf_cuda_ =\
   pscf/cuda/Reduce.cu \
   pscf/cuda/ThreadGrid.cu \
   pscf/cuda/CudaRandom.cu \
   pscf/cuda/VecOp.cu \
   pscf/cuda/VecOpMisc.cu

pscf_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_cuda_:.cu=.o))

