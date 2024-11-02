rpc_fts_montecarlo_= \
  rpc/fts/montecarlo/McSimulator.cpp \
  rpc/fts/montecarlo/McMove.cpp \
  rpc/fts/montecarlo/McMoveFactory.cpp \
  rpc/fts/montecarlo/McMoveManager.cpp \
  rpc/fts/montecarlo/RealMove.cpp \
  rpc/fts/montecarlo/FourierMove.cpp \
  rpc/fts/montecarlo/ForceBiasMove.cpp 
  
rpc_fts_montecarlo_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_montecarlo_:.cpp=.o))

