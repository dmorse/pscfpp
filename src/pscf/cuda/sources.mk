pscf_cuda_CU =\
   pscf/cuda/ThreadArray.cu \
   pscf/cuda/ThreadMesh.cu \
   pscf/cuda/DeviceMemory.cu \
   pscf/cuda/CudaVecRandom.cu \
   pscf/cuda/Reduce.cu \
   pscf/cuda/VecOp.cu \
   pscf/cuda/VecOpMisc.cu 

pscf_cuda_OBJS=\
   $(addprefix $(BLD_DIR)/, $(pscf_cuda_CU:.cu=.o))

pscf_cuda_DEPS=\
   $(addprefix $(BLD_DIR)/, $(pscf_cuda_CU:.cu=.du))
