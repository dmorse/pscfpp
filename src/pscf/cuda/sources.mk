pscf_cuda_ =\
   pscf/cuda/ThreadArray.cu \
   pscf/cuda/ThreadMesh.cu \
   pscf/cuda/DeviceMemory.cu \
   pscf/cuda/CudaVecRandom.cu \
   pscf/cuda/Reduce.cu \
   pscf/cuda/VecOp.cu \
   pscf/cuda/VecOpMisc.cu 

pscf_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_cuda_:.cu=.o))

