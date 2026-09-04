
include $(SRC_DIR)/pscf/backend/cpp/sources.mk
pscf_backend_OBJS = $(pscf_backend_cpp_OBJS)
pscf_backend_DEPS = $(pscf_backend_cpp_DEPS)

ifdef PSCF_CUDA

  include $(SRC_DIR)/pscf/backend/cuda/sources.mk
  pscf_backend_OBJS += $(pscf_backend_cuda_OBJS)
  pscf_backend_DEPS += $(pscf_backend_cuda_DEPS)
  
  pscf_backend_CU = pscf/backend/CUT.cu 
  pscf_backend_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(pscf_backend_CU:.cu=.ou))
  pscf_backend_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(pscf_backend_CU:.cu=.du))

endif
