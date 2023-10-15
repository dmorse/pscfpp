pspc_sweep_= \
  pspc/sweep/FieldState.cpp \
  pspc/sweep/BasisFieldState.cpp \
  pspc/sweep/Sweep.cpp \
  pspc/sweep/LinearSweep.cpp \
  pspc/sweep/SweepFactory.cpp

pspc_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_sweep_:.cpp=.o))

