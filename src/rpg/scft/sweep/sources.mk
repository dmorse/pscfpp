rpg_scft_sweep_= \
  rpg/scft/sweep/FieldState.cu \
  rpg/scft/sweep/BasisFieldState.cu \
  rpg/scft/sweep/Sweep.cu \
  rpg/scft/sweep/LinearSweep.cu \
  rpg/scft/sweep/SweepFactory.cu

rpg_scft_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_scft_sweep_:.cu=.o))

