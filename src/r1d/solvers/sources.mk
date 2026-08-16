r1d_solvers_CPP=\
  r1d/solvers/Propagator.cpp \
  r1d/solvers/Block.cpp \
  r1d/solvers/Polymer.cpp \
  r1d/solvers/Solvent.cpp \
  r1d/solvers/Mixture.cpp

r1d_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_solvers_CPP:.cpp=.o))

r1d_solvers_DEPS=\
     $(addprefix $(BLD_DIR)/, $(r1d_solvers_CPP:.cpp=.d))

