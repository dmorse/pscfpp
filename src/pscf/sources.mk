#-----------------------------------------------------------------------
# Source and object file lists for src/pscf 

# Include source list files from subdirectories
include $(SRC_DIR)/pscf/chem/sources.mk
include $(SRC_DIR)/pscf/inter/sources.mk
include $(SRC_DIR)/pscf/math/sources.mk
include $(SRC_DIR)/pscf/mesh/sources.mk
include $(SRC_DIR)/pscf/floryHuggins/sources.mk
include $(SRC_DIR)/pscf/environment/sources.mk
include $(SRC_DIR)/pscf/iterator/sources.mk
include $(SRC_DIR)/pscf/sweep/sources.mk

# C++ source files

pscf_CPP= \
  $(pscf_chem_) $(pscf_inter_) $(pscf_math_) \
  $(pscf_mesh_) $(pscf_crystal_) $(pscf_floryHuggins_) \
  $(pscf_environment_) $(pscf_iterator_) $(pscf_sweep_)

pscf_CPP_OBJS=\
    $(addprefix $(BLD_DIR)/, $(pscf_CPP:.cpp=.o))

pscf_OBJS = $(pscf_CPP_OBJS)

# CUDA source files

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
