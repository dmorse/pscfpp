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

# Standard C++ files
pscf_OBJS= \
  $(pscf_backends_OBJS) $(pscf_cpu_OBJS) \
  $(pscf_math_OBJS) $(pscf_mesh_OBJS) \
  $(pscf_chem_OBJS) $(pscf_interaction_OBJS) \
  $(pscf_floryHuggins_OBJS) $(pscf_correlation_OBJS) \
  $(pscf_environment_OBJS) $(pscf_iterator_OBJS) \
  $(pscf_sweep_OBJS)

pscf_DEPS= \
  $(pscf_backends_DEPS) $(pscf_cpu_DEPS) \
  $(pscf_math_DEPS) $(pscf_mesh_DEPS) \
  $(pscf_chem_DEPS) $(pscf_interaction_DEPS) \
  $(pscf_floryHuggins_DEPS) $(pscf_correlation_DEPS) \
  $(pscf_environment_DEPS) $(pscf_iterator_DEPS) \
  $(pscf_sweep_DEPS)

# CUDA C++ files
ifdef PSCF_CUDA
  include $(SRC_DIR)/pscf/cuda/sources.mk
  pscf_OBJS+= $(pscf_cuda_OBJS)
  pscf_DEPS+= $(pscf_cuda_DEPS)
endif

## Serial C++ files
#pscf_CPP= \
#  $(pscf_backends_CPP) $(pscf_cpu_CPP) \
#  $(pscf_math_CPP) $(pscf_mesh_CPP) \
#  $(pscf_chem_CPP) $(pscf_interaction_CPP) \
#  $(pscf_floryHuggins_CPP) $(pscf_correlation_CPP) \
#  $(pscf_environment_CPP) $(pscf_iterator_CPP) \
#  $(pscf_sweep_CPP)
#pscf_OBJS=\
#  $(addprefix $(BLD_DIR)/, $(pscf_CPP:.cpp=.o))
#pscf_DEPS=\
#  $(addprefix $(BLD_DIR)/, $(pscf_CPP:.cpp=.d))
#
## CUDA C++ files
#ifdef PSCF_CUDA
#  include $(SRC_DIR)/pscf/cuda/sources.mk
#  pscf_CU= $(pscf_backends_CU) $(pscf_cuda_CU)
#  pscf_OBJS+=\
#    $(addprefix $(BLD_DIR)/, $(pscf_CU:.cu=.ou))
#  pscf_DEPS+=\
#    $(addprefix $(BLD_DIR)/, $(pscf_CU:.cu=.du))
#endif

#-----------------------------------------------------------------------
# Path and rule for libpscf.a library 

pscf_LIB=$(BLD_DIR)/pscf/libpscf.a

$(pscf_LIB): $(pscf_OBJS)
	$(AR) rcs $(pscf_LIB) $(pscf_OBJS) 

#-----------------------------------------------------------------------
