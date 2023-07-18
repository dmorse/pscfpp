pspc_simulate_mcmove_= \
  pspc/simulate/mcmove/McMove.cpp \
  pspc/simulate/mcmove/McMoveFactory.cpp \
  pspc/simulate/mcmove/McMoveManager.cpp \
  pspc/simulate/mcmove/RealMove.cpp 
  
pspc_simulate_mcmove_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_simulate_mcmove_))
pspc_simulate_mcmove_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_simulate_mcmove_:.cpp=.o))

