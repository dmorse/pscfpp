pfts_solvers_= pfts/solvers/Species.cpp 

pfts_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pfts_solvers_))
pfts_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pfts_solvers_:.cpp=.o))

