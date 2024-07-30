#-----------------------------------------------------------------------
# Source files in src/prdc and corresponding object file targets

# CPP Source Files

include $(SRC_DIR)/prdc/crystal/sources.mk
include $(SRC_DIR)/prdc/cpu/sources.mk
include $(SRC_DIR)/prdc/iterator/sources.mk

prdc_CPP= \
  $(prdc_crystal_) \
  $(prdc_cpu_) \
  $(prdc_iterator_) \

prdc_CPP_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_CPP:.cpp=.o))

prdc_OBJS = $(prdc_CPP_OBJS) 

# CUDA Source Files

ifdef PSCF_CUDA
  include $(SRC_DIR)/prdc/cuda/sources.mk
  prdc_OBJS += $(prdc_cuda_OBJS)
endif

#-----------------------------------------------------------------------
# Path and makefile target for libprdc.a library 

prdc_LIBNAME=prdc
prdc_LIB=$(BLD_DIR)/prdc/lib$(prdc_LIBNAME).a

$(prdc_LIB): $(prdc_OBJS)
	$(AR) rcs $(prdc_LIB) $(prdc_OBJS)

