pspg_simulate_bdstep_= \
  pspg/simulate/bdstep/BdSimulator.cu \
  pspg/simulate/bdstep/BdStep.cu \
  pspg/simulate/bdstep/BdStepFactory.cu \
  pspg/simulate/bdstep/ExplicitBdStep.cu \
  pspg/simulate/bdstep/PredCorrBdStep.cu \
  pspg/simulate/bdstep/LMBdStep.cu 
  
pspg_simulate_bdstep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_simulate_bdstep_:.cu=.o))

