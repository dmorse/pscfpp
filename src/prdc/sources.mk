include $(SRC_DIR)/prdc/crystal/sources.mk
include $(SRC_DIR)/prdc/cpu/sources.mk
include $(SRC_DIR)/prdc/gpu/sources.mk

prdc_= \
  $(prdc_crystal_) \
  $(prdc_cpu_) 

ifdef PSCF_CUDA
prdc+= $(prdc_gpu_)
endif

prdc_SRCS=\
     $(addprefix $(SRC_DIR)/, $(prdc_))
prdc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_:.cpp=.o))

$(prdc_LIB): $(prdc_OBJS)
	$(AR) rcs $(prdc_LIB) $(prdc_OBJS)

