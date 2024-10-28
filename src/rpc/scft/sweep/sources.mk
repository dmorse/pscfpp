rpc_scft_sweep_= \
  rpc/scft/sweep/FieldState.cpp \
  rpc/scft/sweep/BasisFieldState.cpp \
  rpc/scft/sweep/Sweep.cpp \
  rpc/scft/sweep/LinearSweep.cpp \
  rpc/scft/sweep/SweepFactory.cpp

rpc_scft_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_scft_sweep_:.cpp=.o))

