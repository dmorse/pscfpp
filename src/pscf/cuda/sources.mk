pscf_cuda_ =\
   pscf/cuda/KernelWrappers.cu \
   pscf/cuda/LinearAlgebra.cu \
   pscf/cuda/ParallelReductions.cu \
   pscf/cuda/ThreadGrid.cu \
   pscf/cuda/CudaRandom.cu \
   pscf/cuda/DeviceArray.cu \
   pscf/cuda/HostDArray.cu \
   pscf/cuda/VecOp.cu

pscf_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_cuda_:.cu=.o))

