
fd1d_solvers_=\
  fd1d/solvers/Propagator.cpp \
  fd1d/solvers/Block.cpp \
  fd1d/solvers/Polymer.cpp \
  fd1d/solvers/Solvent.cpp \
  fd1d/solvers/Mixture.cpp 

fd1d_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(fd1d_solvers_))
fd1d_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_solvers_:.cpp=.o))

