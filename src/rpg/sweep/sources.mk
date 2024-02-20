rpg_sweep_= \
  rpg/sweep/FieldState.cu \
  rpg/sweep/BasisFieldState.cu \
  rpg/sweep/Sweep.cu \
  rpg/sweep/LinearSweep.cu \
  rpg/sweep/SweepFactory.cu

rpg_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_sweep_:.cu=.o))

