rpc_solvers_= \
  rpc/solvers/Propagator.cpp \
  rpc/solvers/Block.cpp \
  rpc/solvers/Polymer.cpp \
  rpc/solvers/Solvent.cpp \
  rpc/solvers/Mixture.cpp

rpc_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_solvers_:.cpp=.o))

