r1d_solvers_=\
  r1d/solvers/Propagator.cpp \
  r1d/solvers/Block.cpp \
  r1d/solvers/Polymer.cpp \
  r1d/solvers/Mixture.cpp \
  r1d/solvers/Solvent.cpp \

r1d_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_solvers_:.cpp=.o))

