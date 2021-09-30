
fd1d_sweep_=\
  fd1d/sweep/Sweep.cpp \
  fd1d/sweep/CompositionSweep.cpp \
  fd1d/sweep/MuSweep.cpp \
  fd1d/sweep/LengthSweep.cpp \
  fd1d/sweep/SweepFactory.cpp 

fd1d_sweep_SRCS=\
     $(addprefix $(SRC_DIR)/, $(fd1d_sweep_))
fd1d_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_sweep_:.cpp=.o))

