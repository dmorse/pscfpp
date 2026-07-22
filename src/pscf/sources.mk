#-----------------------------------------------------------------------
# Source and object file lists for src/pscf 

# Directories with only *.cpp source files
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

pscf_CPP= \
  $(pscf_cpu_) \
  $(pscf_math_) $(pscf_mesh_) \
  $(pscf_chem_) $(pscf_interaction_) \
  $(pscf_floryHuggins_) $(pscf_correlation_) \
  $(pscf_environment_) $(pscf_iterator_) \
  $(pscf_sweep_)

pscf_OBJS=\
    $(addprefix $(BLD_DIR)/, $(pscf_CPP:.cpp=.o))

# Directories with mixed file types

include $(SRC_DIR)/pscf/backends/sources.mk
pscf_OBJS+= $(pscf_backends_OBJS)

# Directories with only CUDA source files

ifdef PSCF_CUDA
  include $(SRC_DIR)/pscf/cuda/sources.mk
  pscf_OBJS+=$(pscf_cuda_OBJS)
endif

#-----------------------------------------------------------------------
# Path and rule for libpscf.a library 

pscf_LIB=$(BLD_DIR)/pscf/libpscf.a

$(pscf_LIB): $(pscf_OBJS)
	$(AR) rcs $(pscf_LIB) $(pscf_OBJS) 

#-----------------------------------------------------------------------
