rpg_trajectory_= \
  rpg/simulate/trajectory/TrajectoryReader.cu \
  rpg/simulate/trajectory/TrajectoryReaderFactory.cu \
  rpg/simulate/trajectory/FieldConfigReader.cu 

rpg_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_trajectory_:.cu=.o))

