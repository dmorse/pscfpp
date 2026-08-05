rp_system_OBJS=
rp_system_DEPS=

ifdef PSCF_CPP
  rp_system_CPP= \
    rp/system/System_c.cpp \
    rp/system/SystemConstRef_c.cpp
  rp_system_OBJS += \
       $(addprefix $(BLD_DIR)/, $(rp_system_CPP:.cpp=.o))
  rp_system_DEPS += \
       $(addprefix $(BLD_DIR)/, $(rp_system_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_system_CU= \
    rp/system/System_u.cu \
    rp/system/SystemConstRef_u.cu
  rp_system_OBJS +=\
       $(addprefix $(BLD_DIR)/, $(rp_system_CU:.cu=.o))
  rp_system_DEPS +=\
       $(addprefix $(BLD_DIR)/, $(rp_system_CU:.cu=.d))
endif
