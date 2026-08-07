#-----------------------------------------------------------------------
# Source and object file lists for src/prdc

# Include source list files from subdirectories
include $(SRC_DIR)/prdc/field/cpu/sources.mk
include $(SRC_DIR)/prdc/field/cuda/sources.mk

# Standard C++ source files
prdc_field_CPP= \
  $(prdc_field_cpu_CPP)
prdc_field_OBJS=\
  $(addprefix $(BLD_DIR)/, $(prdc_field_CPP:.cpp=.o))
prdc_field_DEPS=\
  $(addprefix $(BLD_DIR)/, $(prdc_field_CPP:.cpp=.d))

# CUDA C++ source files
ifdef PSCF_CUDA
  prdc_field_CU= $(prdc_field_cuda_CU)
  prdc_field_OBJS +=\
    $(addprefix $(BLD_DIR)/, $(prdc_field_CU:.cu=.ou))
  prdc_field_DEPS +=\
    $(addprefix $(BLD_DIR)/, $(prdc_field_CU:.cu=.du))
endif

