rpc_fts_mcmove_= \
  rpc/fts/mcmove/McSimulator.cpp \
  rpc/fts/mcmove/McMove.cpp \
  rpc/fts/mcmove/McMoveFactory.cpp \
  rpc/fts/mcmove/McMoveManager.cpp \
  rpc/fts/mcmove/RealMove.cpp \
  rpc/fts/mcmove/FourierMove.cpp \
  rpc/fts/mcmove/ForceBiasMove.cpp 
  
rpc_fts_mcmove_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_mcmove_:.cpp=.o))

