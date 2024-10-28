rpg_fts_trajectory_= \
  rpg/fts/trajectory/TrajectoryReader.cu \
  rpg/fts/trajectory/TrajectoryReaderFactory.cu \
  rpg/fts/trajectory/FieldConfigReader.cu 

rpg_fts_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_trajectory_:.cu=.o))

