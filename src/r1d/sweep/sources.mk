r1d_sweep_=\
  r1d/sweep/Sweep.cpp \
  r1d/sweep/SweepFactory.cpp \
  r1d/sweep/SweepParameter.cpp \
  r1d/sweep/LinearSweep.cpp

r1d_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_sweep_:.cpp=.o))

