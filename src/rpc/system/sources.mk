rpc_system_= \
  rpc/system/SystemConstRef.cpp

rpc_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_system_:.cpp=.o))

