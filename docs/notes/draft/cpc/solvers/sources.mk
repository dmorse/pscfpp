cpc_solvers_= \
  cpc/solvers/Propagator.cpp \
  cpc/solvers/Block.cpp \
  cpc/solvers/Polymer.cpp \
  cpc/solvers/Solvent.cpp \
  cpc/solvers/Mixture.cpp

cpc_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_solvers_:.cpp=.o))

