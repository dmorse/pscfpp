
# Standard C++ source files from directory prdc/field/cpp
include $(SRC_DIR)/prdc/field/cpp/sources.mk
prdc_field_OBJS=$(prdc_field_cpp_OBJS)
prdc_field_DEPS=$(prdc_field_cpp_DEPS)

# CUDA C++ source files from directory prdc/field/cuda
ifdef PSCF_CUDA
  include $(SRC_DIR)/prdc/field/cuda/sources.mk
  prdc_field_OBJS+= $(prdc_field_cuda_OBJS)
  prdc_field_DEPS+= $(prdc_field_cuda_DEPS)
endif

# - Subdirectory prdc/field/cuda contains only *.cpp  C++ files,
#   which are always compiled.
#
# - Subdirectory prdc/field/cuda contains only *.cu CUDA files,
#   which are compiled only if PSCF_CUDA is defined.
