rpg_simulate_bdstep_= \
  rpg/simulate/bdstep/BdSimulator.cu \
  rpg/simulate/bdstep/BdStep.cu \
  rpg/simulate/bdstep/BdStepFactory.cu \
  rpg/simulate/bdstep/ExplicitBdStep.cu \
  rpg/simulate/bdstep/PredCorrBdStep.cu \
  rpg/simulate/bdstep/LMBdStep.cu 
  
rpg_simulate_bdstep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_simulate_bdstep_:.cu=.o))

