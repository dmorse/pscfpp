cpc_system_= \
  cpc/system/System.cpp \
  cpc/system/Kernel.cpp 

cpc_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_system_:.cpp=.o))

