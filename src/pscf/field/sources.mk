#-----------------------------------------------------------------------
# Include source list files from subdirectories

include $(SRC_DIR)/pscf/field/cpp/sources.mk
include $(SRC_DIR)/pscf/field/cuda/sources.mk

#-----------------------------------------------------------------------

# Standard C++ source files from directory pscf/field/cpp
include $(SRC_DIR)/pscf/field/cpp/sources.mk
pscf_field_OBJS=$(pscf_field_cpp_OBJS)
pscf_field_DEPS=$(pscf_field_cpp_DEPS)

# CUDA C++ source files from directory pscf/field/cuda
ifdef PSCF_CUDA
  include $(SRC_DIR)/pscf/field/cuda/sources.mk
  pscf_field_OBJS+= $(pscf_field_cuda_OBJS)
  pscf_field_DEPS+= $(pscf_field_cuda_DEPS)
endif

# - Subdirectory pscf/field/cuda contains only *.cpp  C++ files,
#   which are always compiled.
#
# - Subdirectory pscf/field/cuda contains only *.cu CUDA files,
#   which are compiled only if PSCF_CUDA is defined.
