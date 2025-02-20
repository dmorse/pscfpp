rpc_fts_trajectory_= \
  rpc/fts/trajectory/TrajectoryReader.cpp \
  rpc/fts/trajectory/TrajectoryReaderFactory.cpp \
  rpc/fts/trajectory/RGridTrajectoryReader.cpp 

rpc_fts_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_trajectory_:.cpp=.o))

