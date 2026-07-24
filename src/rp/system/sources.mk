rp_system_OBJS=

ifdef PSCF_CPP
rp_system_CPP= \
  rp/system/System_c.cpp 
  #rp/system/SystemConstRef_c.cpp 
rp_system_OBJS +=\
     $(addprefix $(BLD_DIR)/, $(rp_system_CPP:.cpp=.o))
endif

ifdef PSCF_CUDA
rp_system_CUDA= \
  rp/system/System_u.cu 
  #rp/system/SystemConstRef_u.cu 
rp_system_OBJS +=\
     $(addprefix $(BLD_DIR)/, $(rp_system_CUDA:.cu=.o))
endif

