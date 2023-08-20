pspg_trajectory_= \
  pspg/simulate/trajectory/TrajectoryReader.cu \
  pspg/simulate/trajectory/TrajectoryReaderFactory.cu \
  pspg/simulate/trajectory/FieldConfigReader.cu 

pspg_trajectory_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_trajectory_))
pspg_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_trajectory_:.cu=.o))

