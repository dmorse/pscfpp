pscf_cuda_ =\
   pscf/cuda/KernelWrappers.cu \
   pscf/cuda/LinearAlgebra.cu \
   pscf/cuda/ParallelReductions.cu \
   pscf/cuda/ThreadGrid.cu

pscf_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_cuda_:.cu=.o))

