rpc_fts_brownian_= \
  rpc/fts/brownian/BdSimulator.cpp \
  rpc/fts/brownian/BdStep.cpp \
  rpc/fts/brownian/BdStepFactory.cpp \
  rpc/fts/brownian/ExplicitBdStep.cpp \
  rpc/fts/brownian/MidstepBdStep.cpp \
  rpc/fts/brownian/PredCorrBdStep.cpp \
  rpc/fts/brownian/LMBdStep.cpp 
  
rpc_fts_brownian_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_brownian_:.cpp=.o))

