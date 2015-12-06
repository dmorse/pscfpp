pfts_solvers_= pfts/solvers/Species.cpp \
  pfts/solvers/PropagatorStub.cpp \
  pfts/solvers/PolymerStub.cpp \
  pfts/solvers/SolventStub.cpp \
  pfts/solvers/SystemStub.cpp 

pfts_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pfts_solvers_))
pfts_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pfts_solvers_:.cpp=.o))

