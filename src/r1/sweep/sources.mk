r1_sweep_CPP=\
  r1/sweep/Sweep.cpp \
  r1/sweep/SweepFactory.cpp \
  r1/sweep/SweepParameter.cpp \
  r1/sweep/LinearSweep.cpp

r1_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1_sweep_CPP:.cpp=.o))

r1_sweep_DEPS=\
     $(addprefix $(BLD_DIR)/, $(r1_sweep_CPP:.cpp=.d))

