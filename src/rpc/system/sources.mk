rpc_system_= \
  rpc/system/System.cpp \
  rpc/system/SystemConstRef.cpp \
  rpc/system/Types.cpp

rpc_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_system_:.cpp=.o))

