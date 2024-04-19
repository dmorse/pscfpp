rpg_simulate_trajectory_= \
  rpg/simulate/trajectory/TrajectoryReader.cu \
  rpg/simulate/trajectory/TrajectoryReaderFactory.cu \
  rpg/simulate/trajectory/FieldConfigReader.cu 

rpg_simulate_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_simulate_trajectory_:.cu=.o))

