
fd1d_sweep_=\
  fd1d/sweep/Sweep.cpp \
  fd1d/sweep/SweepFactory.cpp \
  fd1d/sweep/SweepParameter.cpp \
  fd1d/sweep/LinearSweep.cpp

fd1d_sweep_SRCS=\
     $(addprefix $(SRC_DIR)/, $(fd1d_sweep_))
fd1d_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_sweep_:.cpp=.o))

