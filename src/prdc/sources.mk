#-----------------------------------------------------------------------
# Source files in src/prdc and corresponding object file targets

# Include source list files from subdirectories
include $(SRC_DIR)/prdc/crystal/sources.mk
include $(SRC_DIR)/prdc/cpu/sources.mk
include $(SRC_DIR)/prdc/field/sources.mk
include $(SRC_DIR)/prdc/environment/sources.mk

# C++ source files

prdc_CPP= \
  $(prdc_crystal_) \
  $(prdc_cpu_) \
  $(prdc_field_) \
  $(prdc_environment_) \

prdc_CPP_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_CPP:.cpp=.o))

prdc_OBJS = $(prdc_CPP_OBJS) 

# CUDA C++ source files

ifdef PSCF_CUDA
  include $(SRC_DIR)/prdc/cuda/sources.mk
  prdc_OBJS += $(prdc_cuda_OBJS)
endif

#-----------------------------------------------------------------------
# Path and makefile target for prdc/libprdc.a library 

prdc_LIBNAME=prdc
prdc_LIB=$(BLD_DIR)/prdc/lib$(prdc_LIBNAME).a

$(prdc_LIB): $(prdc_OBJS)
	$(AR) rcs $(prdc_LIB) $(prdc_OBJS)

