rpg_fts_brownian_= \
  rpg/fts/brownian/BdSimulator.cu \
  rpg/fts/brownian/BdStep.cu \
  rpg/fts/brownian/BdStepFactory.cu \
  rpg/fts/brownian/ExplicitBdStep.cu \
  rpg/fts/brownian/PredCorrBdStep.cu \
  rpg/fts/brownian/LMBdStep.cu 
  
rpg_fts_brownian_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_brownian_:.cu=.ou))

