
include $(SRC_DIR)/pscf/backend/cpp/sources.mk
pscf_backend_OBJS = $(pscf_backend_cpp_OBJS)
pscf_backend_DEPS = $(pscf_backend_cpp_DEPS)

pscf_backend_CPP = pscf/backend/CPT.cpp
pscf_backend_OBJS+=\
     $(addprefix $(BLD_DIR)/, $(pscf_backend_CPP:.cpp=.o))
pscf_backend_DEPS+=\
     $(addprefix $(BLD_DIR)/, $(pscf_backend_CPP:.cpp=.d))

ifdef PSCF_CUDA
pscf_backend_CU = pscf/backend/CUT.cu 
pscf_backend_OBJS+=\
     $(addprefix $(BLD_DIR)/, $(pscf_backend_CU:.cu=.ou))
pscf_backend_DEPS+=\
     $(addprefix $(BLD_DIR)/, $(pscf_backend_CU:.cu=.du))
endif
