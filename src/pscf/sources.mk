#-----------------------------------------------------------------------
# Source and object file lists for src/pscf 

include $(SRC_DIR)/pscf/backends/sources.mk
include $(SRC_DIR)/pscf/cpu/sources.mk
include $(SRC_DIR)/pscf/math/sources.mk
include $(SRC_DIR)/pscf/mesh/sources.mk
include $(SRC_DIR)/pscf/chem/sources.mk
include $(SRC_DIR)/pscf/interaction/sources.mk
include $(SRC_DIR)/pscf/floryHuggins/sources.mk
include $(SRC_DIR)/pscf/correlation/sources.mk
include $(SRC_DIR)/pscf/environment/sources.mk
include $(SRC_DIR)/pscf/iterator/sources.mk
include $(SRC_DIR)/pscf/sweep/sources.mk

# Serial C++ files
pscf_CPP= \
  $(pscf_backends_CPP) $(pscf_cpu_CPP) \
  $(pscf_math_CPP) $(pscf_mesh_CPP) \
  $(pscf_chem_CPP) $(pscf_interaction_CPP) \
  $(pscf_floryHuggins_CPP) $(pscf_correlation_CPP) \
  $(pscf_environment_CPP) $(pscf_iterator_CPP) \
  $(pscf_sweep_CPP)
pscf_OBJS=\
  $(addprefix $(BLD_DIR)/, $(pscf_CPP:.cpp=.o))
pscf_DEPS=\
  $(addprefix $(BLD_DIR)/, $(pscf_CPP:.cpp=.d))

# CUDA C++ files
ifdef PSCF_CUDA
  include $(SRC_DIR)/pscf/cuda/sources.mk
  pscf_CU= $(pscf_backends_CU) $(pscf_cuda_CU)
  pscf_OBJS+=\
    $(addprefix $(BLD_DIR)/, $(pscf_CU:.cu=.o))
  pscf_DEPS+=\
    $(addprefix $(BLD_DIR)/, $(pscf_CU:.cu=.du))
endif

#-----------------------------------------------------------------------
# Path and rule for libpscf.a library 

pscf_LIB=$(BLD_DIR)/pscf/libpscf.a

$(pscf_LIB): $(pscf_OBJS)
	$(AR) rcs $(pscf_LIB) $(pscf_OBJS) 

#-----------------------------------------------------------------------
