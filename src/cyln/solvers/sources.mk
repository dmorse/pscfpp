
cyln_solvers_=\
  cyln/solvers/Block.cpp 
  #cyln/solvers/Propagator.cpp \
  #cyln/solvers/Polymer.cpp \
  #cyln/solvers/Solvent.cpp \
  #cyln/solvers/Mixture.cpp 

cyln_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(cyln_solvers_))
cyln_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cyln_solvers_:.cpp=.o))

