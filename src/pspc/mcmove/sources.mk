pspc_mcmove_= \
  pspc/mcmove/McMove.cpp \
  pspc/mcmove/McState.cpp \
  pspc/mcmove/McMoveFactory.cpp \
  pspc/mcmove/McMoveManager.cpp \
  pspc/mcmove/McSimulator.cpp \
  pspc/mcmove/RealMove.cpp \
  pspc/mcmove/FourierMove.cpp 
  
pspc_mcmove_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_mcmove_))
pspc_mcmove_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_mcmove_:.cpp=.o))

