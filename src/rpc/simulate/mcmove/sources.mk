rpc_simulate_mcmove_= \
  rpc/simulate/mcmove/McSimulator.cpp \
  rpc/simulate/mcmove/McState.cpp \
  rpc/simulate/mcmove/McMove.cpp \
  rpc/simulate/mcmove/McMoveFactory.cpp \
  rpc/simulate/mcmove/McMoveManager.cpp \
  rpc/simulate/mcmove/RealMove.cpp \
  rpc/simulate/mcmove/FourierMove.cpp \
  rpc/simulate/mcmove/ForceBiasMove.cpp 
  
rpc_simulate_mcmove_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_simulate_mcmove_:.cpp=.o))

