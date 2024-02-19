rpc_sweep_= \
  rpc/sweep/FieldState.cpp \
  rpc/sweep/BasisFieldState.cpp \
  rpc/sweep/Sweep.cpp \
  rpc/sweep/LinearSweep.cpp \
  rpc/sweep/SweepFactory.cpp

rpc_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_sweep_:.cpp=.o))

