#-----------------------------------------------------------------------
# Source files in src/pscf and associated object file targets

# CPP source directories
include $(SRC_DIR)/pscf/chem/sources.mk
include $(SRC_DIR)/pscf/inter/sources.mk
include $(SRC_DIR)/pscf/math/sources.mk
include $(SRC_DIR)/pscf/mesh/sources.mk
include $(SRC_DIR)/pscf/homogeneous/sources.mk
include $(SRC_DIR)/pscf/environment/sources.mk
include $(SRC_DIR)/pscf/iterator/sources.mk
include $(SRC_DIR)/pscf/sweep/sources.mk

# CPP source files

pscf_CPP= \
  $(pscf_chem_) $(pscf_inter_) $(pscf_math_) \
  $(pscf_mesh_) $(pscf_crystal_) $(pscf_homogeneous_) \
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
# Path and makefile target for libpscf.a library 

pscf_LIBNAME=pscf
pscf_LIB=$(BLD_DIR)/pscf/lib$(pscf_LIBNAME).a

$(pscf_LIB): $(pscf_OBJS)
	$(AR) rcs $(pscf_LIB) $(pscf_OBJS) 

