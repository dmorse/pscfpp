
pssp_solvers_=
  #pssp/solvers/Propagator.cpp \
  #pssp/solvers/Block.cpp \
  #pssp/solvers/Polymer.cpp \
  #pssp/solvers/Solvent.cpp \
  #pssp/solvers/Mixture.cpp 

pssp_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_solvers_))
pssp_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_solvers_:.cpp=.o))

