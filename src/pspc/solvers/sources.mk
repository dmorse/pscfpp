
pspc_solvers_= \
  pspc/solvers/Propagator.cpp \
  pspc/solvers/Block.cpp \
  pspc/solvers/Polymer.cpp \
  pspc/solvers/Mixture.cpp

  #pspc/solvers/Solvent.cpp \

pspc_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_solvers_))
pspc_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_solvers_:.cpp=.o))

