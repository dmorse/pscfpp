rpg_fts_bdstep_= \
  rpg/fts/bdstep/BdSimulator.cu \
  rpg/fts/bdstep/BdStep.cu \
  rpg/fts/bdstep/BdStepFactory.cu \
  rpg/fts/bdstep/ExplicitBdStep.cu \
  rpg/fts/bdstep/PredCorrBdStep.cu \
  rpg/fts/bdstep/LMBdStep.cu 
  
rpg_fts_bdstep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_bdstep_:.cu=.o))

