pscf_backend_cuda_CU =\
   pscf/backend/cuda/CUT.cu \
   pscf/backend/cuda/ThreadArray.cu \
   pscf/backend/cuda/ThreadMesh.cu \
   pscf/backend/cuda/DeviceMemory.cu \
   pscf/backend/cuda/CudaVecRandom.cu \
   pscf/backend/cuda/Reduce.cu \
   pscf/backend/cuda/VecOp.cu \
   pscf/backend/cuda/VecOpMisc.cu 

pscf_backend_cuda_OBJS=\
   $(addprefix $(BLD_DIR)/, $(pscf_backend_cuda_CU:.cu=.ou))

pscf_backend_cuda_DEPS=\
   $(addprefix $(BLD_DIR)/, $(pscf_backend_cuda_CU:.cu=.du))
