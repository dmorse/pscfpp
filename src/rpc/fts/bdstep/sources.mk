rpc_fts_bdstep_= \
  rpc/fts/bdstep/BdSimulator.cpp \
  rpc/fts/bdstep/BdStep.cpp \
  rpc/fts/bdstep/BdStepFactory.cpp \
  rpc/fts/bdstep/ExplicitBdStep.cpp \
  rpc/fts/bdstep/MidstepBdStep.cpp \
  rpc/fts/bdstep/PredCorrBdStep.cpp \
  rpc/fts/bdstep/LMBdStep.cpp 
  
rpc_fts_bdstep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_bdstep_:.cpp=.o))

