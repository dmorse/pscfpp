r1_solvers_CPP=\
  r1/solvers/Propagator.cpp \
  r1/solvers/Block.cpp \
  r1/solvers/Polymer.cpp \
  r1/solvers/Solvent.cpp \
  r1/solvers/Mixture.cpp

r1_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1_solvers_CPP:.cpp=.o))

r1_solvers_DEPS=\
     $(addprefix $(BLD_DIR)/, $(r1_solvers_CPP:.cpp=.d))

