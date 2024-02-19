rpc_trajectory_= \
  rpc/simulate/trajectory/TrajectoryReader.cpp \
  rpc/simulate/trajectory/TrajectoryReaderFactory.cpp \
  rpc/simulate/trajectory/FieldConfigReader.cpp 

rpc_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_trajectory_:.cpp=.o))

