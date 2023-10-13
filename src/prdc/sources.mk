# CPP Source Files
include $(SRC_DIR)/prdc/crystal/sources.mk
include $(SRC_DIR)/prdc/cpu/sources.mk

prdc_CPP= \
  $(prdc_crystal_) \
  $(prdc_cpu_) 

prdc_CPP_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_CPP:.cpp=.o))

prdc_OBJS=$(pdrc_CPP_OBJS) 

# CUDA Source Files
ifdef PSCF_CUDA

  include $(SRC_DIR)/prdc/gpu/sources.mk

  prdc_CUDA_OBJS=\
	$(addprefix $(BLD_DIR)/, $(prdc_gpu_:.cu=.o))

  prdc_OBJS = $(prdc_CPP_OBJS) $(prdc_CUDA_OBJS)

endif

$(prdc_LIB): $(prdc_OBJS)
	$(AR) rcs $(prdc_LIB) $(prdc_OBJS)

