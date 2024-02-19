rpc_simulate_bdstep_= \
  rpc/simulate/bdstep/BdSimulator.cpp \
  rpc/simulate/bdstep/BdStep.cpp \
  rpc/simulate/bdstep/BdStepFactory.cpp \
  rpc/simulate/bdstep/ExplicitBdStep.cpp \
  rpc/simulate/bdstep/MidstepBdStep.cpp \
  rpc/simulate/bdstep/PredCorrBdStep.cpp \
  rpc/simulate/bdstep/LMBdStep.cpp 
  
rpc_simulate_bdstep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_simulate_bdstep_:.cpp=.o))

