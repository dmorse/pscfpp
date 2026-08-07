rp_system_OBJS=
rp_system_DEPS=

ifdef PSCF_CPP
  rp_system_CPP= \
    rp/system/System.cpp \
    rp/system/SystemConstRef.cpp
  rp_system_OBJS += \
       $(addprefix $(BLD_DIR)/, $(rp_system_CPP:.cpp=.o))
  rp_system_DEPS += \
       $(addprefix $(BLD_DIR)/, $(rp_system_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_system_CU= \
    rp/system/System.cu \
    rp/system/SystemConstRef.cu
  rp_system_OBJS +=\
       $(addprefix $(BLD_DIR)/, $(rp_system_CU:.cu=.ou))
  rp_system_DEPS +=\
       $(addprefix $(BLD_DIR)/, $(rp_system_CU:.cu=.du))
endif
