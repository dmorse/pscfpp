#-----------------------------------------------------------------------
# Source and object file lists for src/prdc

# Include source list files from subdirectories
include $(SRC_DIR)/prdc/crystal/sources.mk
include $(SRC_DIR)/prdc/field/cpu/sources.mk
include $(SRC_DIR)/prdc/field/cuda/sources.mk
include $(SRC_DIR)/prdc/fieldIo/sources.mk
include $(SRC_DIR)/prdc/environment/sources.mk

# Standard C++ source files
prdc_CPP= \
  $(prdc_crystal_CPP) \
  $(prdc_field_cpu_CPP) \
  $(prdc_fieldIo_CPP) \
  $(prdc_environment_CPP)
prdc_OBJS=\
  $(addprefix $(BLD_DIR)/, $(prdc_CPP:.cpp=.o))
prdc_DEPSS=\
  $(addprefix $(BLD_DIR)/, $(prdc_CPP:.cpp=.d))

# CUDA C++ source files
ifdef PSCF_CUDA
  prdc_CU= $(prdc_field_cuda_CU)
  prdc_OBJS +=\
    $(addprefix $(BLD_DIR)/, $(prdc_CU:.cu=.o))
  prdc_DEPS +=\
    $(addprefix $(BLD_DIR)/, $(prdc_CU:.cu=.d))
endif

#-----------------------------------------------------------------------
# Path and rule for the prdc/libprdc.a library

prdc_LIB=$(BLD_DIR)/prdc/libprdc.a

$(prdc_LIB): $(prdc_OBJS)
	$(AR) rcs $(prdc_LIB) $(prdc_OBJS)

